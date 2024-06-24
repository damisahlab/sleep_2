import os
import numpy as np
import pandas as pd
import datetime

import h5py
import hdf5storage
import collections as cl

import docker

import probeinterface
import spikeinterface as si
import spikeinterface.extractors as se
import spikeinterface.sorters as ss
import spikeinterface.comparison as sc
import spikeinterface.widgets as sw

from spikeinterface.exporters import export_to_phy
from spikeinterface.exporters import export_report

from probeinterface import generate_linear_probe

def mat_to_npy(file_path, start, stop):
    # Import MAT file
    data = hdf5storage.loadmat(file_path) # , variable_names = ['foo', 'bar']

    # Extract sampling frequency
    sampling_rate = data['meta_data']['sampling_rate'].astype(np.int64)[0][0][0]

    # NSx raw timestamps are formatted YYYY-MM-?-DD-HH-MM-SS-MSMSMS
    # (so you do not use the third element in the timestamp; 
    #  also, MNE only accepts timestamps in the UTC timezone)
    raw_time = data['meta_data']['time_stamp'][0][0][0][0].astype(np.int64)
    time_start = datetime.datetime(raw_time[0], raw_time[1], raw_time[3], 
                                   raw_time[4], raw_time[5], raw_time[6], 
                                   raw_time[7], tzinfo = datetime.timezone.utc)

    # Convert tmin/tmax to sample numbers
    tmin = (start - time_start).total_seconds() * sampling_rate
    tmax = (stop - time_start).total_seconds() * sampling_rate

    # Convert time series data from uV to Volts
    data = data['time_series'][0] * 1e-6

    # Since this is only one channel, it loads as only one dimension. We need to add a 
    # dummy second dimension (channels) for se.NumpyRecording() to work.
    data = data[:, np.newaxis]

    # Crop the data to the relevant time
    data = data[tmin.astype(np.int64) : tmax.astype(np.int64) + 1, :]

    return data, sampling_rate

def npy_to_si(data, sampling_rate, ch_num, bandpass_filter = True, re_reference = False):
    # Bring into SpikeInterface
    data = se.NumpyRecording([data], sampling_frequency = sampling_rate) #t_starts, channel_ids
    #sw.plot_timeseries(data, time_range = (0, 5))

    # Some data formats store time series data as values that must be transformed
    # by gain and offset parameters (traces_uV = traces_raw * gains + offset) prior
    # to analysis. The Combinato Docker Image requires these values for analysis,
    # so we set dummy gain/offsets that result in the same values:
    # https://spikeinterface.readthedocs.io/en/latest/modules/extractors/plot_2_working_with_unscaled_traces.html#working-with-unscaled-traces
    #data.set_channel_gains(1)
    #data.set_channel_offsets(0)
    #print(data.get_property('gain_to_uV'))
    #print(data.get_property('offset_to_uV'))

    # Most of the spike-sorters require channel locations in order to 
    # cross-reference putative units between nearby channels. Waveclus
    # is one of the few that does not require this (although it can use
    # channel locations for tetrodes). One solution is to process channels
    # only one at a time, creating a dummy single-channel probe for each.
    # The other option is to put all of the channels on a single probe, but 
    # set their locations quite far away from each other (you could also put 
    # each channel on a different shank or in a different probe group, but 
    # that is more work). https://github.com/SpikeInterface/spikeinterface/issues/313

    # Create a 2D dummy probe object with 1 contact
    if ch_num == 1:

        dummy_probe = probeinterface.Probe()
        dummy_positions = np.array([[0, 0]])
        dummy_probe.set_contacts(positions = dummy_positions)
        #dummy_probe.get_contact_count()
        dummy_probe.set_device_channel_indices(channel_indices = [0])

    # Create a 2d linear probe with 10cm channel spacing
    elif ch_num > 1:
        dummy_probe = generate_linear_probe(num_elec = ch_num,
                                            ypitch = 100000, #inter-probe distance in um
                                            contact_shapes = 'circle',
                                            contact_shape_params = {'radius': 1})
        
        dummy_indices = np.arange(ch_num)
        dummy_probe.set_device_channel_indices(dummy_indices)

    # Specify the data's probe as the dummy probe
    data = data.set_probe(probe = dummy_probe)

    # Apply bandpass to only look at action potentials
    # (or tell spikeinterface that it is already filtered)
    # BlackRock used a 250 Hz - 7.5 KHz hardware bandpass 
    # filter on some of our microelectrode data.
    if bandpass_filter == True:
        data = si.preprocessing.bandpass_filter(data, freq_min = 300, freq_max = 6000)

    elif bandpass_filter == False:
        data.annotate(is_filtered = True)
    
    # Optionally re-reference to common median reference
    if re_reference == True:
        data = si.preprocessing.common_reference(data, reference = 'global', operator = 'median')

    return data

def waveclus(data, sorter_folder):
    # Set custom Waveclus parameters
    waveclus_params = ss.get_default_sorter_params('waveclus')
    waveclus_params['detect_sign'] = 0 # select both positive and negative spikes
    print(waveclus_params)

    # Run Waveclus
    # You need to add the path to the Waveclus package to your system path 
    # (create a new variable called WAVECLUS_PATH)
    sorted_data = ss.run_waveclus(recording = data, 
                                  output_folder = sorter_folder, 
                                  remove_existing_folder = True,
                                  delete_output_folder = True,
                                  docker_image = False,
                                  **waveclus_params)

    print('Units found:', sorted_data.get_unit_ids())

    return sorted_data

def combinato(data, sorter_folder):
    combinato_params = ss.get_default_sorter_params('combinato')
    combinato_params['detect_sign'] = 0 # select both positive and negative spikes
    print(combinato_params)

    # sorted_data = ss.run_sorter(sorter_name = 'combinato',
    #                             recording = data, 
    #                             output_folder = sorter_folder,
    #                             remove_existing_folder = True,
    #                             delete_output_folder = False,
    #                             docker_image = False,
    #                             **combinato_params)

    sorted_data = ss.run_combinato(recording = data, 
                                   output_folder = sorter_folder, 
                                   remove_existing_folder = True,
                                   delete_output_folder = False,
                                   docker_image = False,
                                   **combinato_params)
    
    return sorted_data

def docker_sorter(data, sorter, sorter_folder):
    sorted_data = ss.run_sorter(sorter_name = sorter,
                                recording = data, 
                                output_folder = sorter_folder,
                                remove_existing_folder = True,
                                delete_output_folder = True,
                                docker_image = True)
    
    return sorted_data

def process_spikes(data, sorted_data, waveform_folder, phy_folder, report_folder):
    # Waveform Extractor
    waveforms = si.WaveformExtractor.create(data, sorted_data, folder = waveform_folder, remove_if_exists = True)
    waveforms.set_params(ms_before = 3., ms_after = 4., max_spikes_per_unit = 500)
    waveforms.run_extract_waveforms(n_jobs = -1, chunk_size = 30000, progress_bar = True)

    # Specify desired Quality Control Metrics
    unit_quality_metrics = ['amplitude_cutoff', 'firing_rate', 'isi_violation', 'presence_ratio', 'snr'] 

    # Get default parameters for calculation of quality metrics
    # Set the ISI Violation threshold to 3ms (refractory period)
    qm_parameters = si.qualitymetrics.quality_metric_calculator.get_default_qm_params()
    qm_parameters['isi_violations']['isi_threshold_ms'] = 3
    print(qm_parameters)

    # Computer quality metrics
    #metrics = si.qualitymetrics.compute_quality_metrics(waveform_extractor = waveforms, 
    #                                                    metric_names = unit_quality_metrics,
    #                                                    peak_sign = 'both') # DEPRECATED
    metrics = si.qualitymetrics.compute_quality_metrics(waveform_extractor = waveforms, 
                                                        metric_names = unit_quality_metrics,
                                                        qm_params = qm_parameters)

    # Export for use in Phy
    export_to_phy(waveforms, output_folder = phy_folder, remove_if_exists = True, peak_sign = 'both')

    # Export a report
    export_report(waveforms, output_folder = report_folder)

def wclus_to_pd(file_path):
    '''
    This function takes the 'times_[filename].mat' output file from
    Waveclus and converts it to a Pandas dataframe.
    
    Input:
    - file_path = absolute file path to Waveclus output ('times_[filename].mat')
        (these files contain the unit ID's and spike times in the "cluster_class" variable)

    Output:
    - data = Pandas dataframe (cluster_id, time) 
    '''
    import hdf5storage

    data = hdf5storage.loadmat(file_path)
    data = data['cluster_class']
    data = pd.DataFrame(data, columns = ['cluster_id', 'time_ms'])
    data['cluster_id'] = data['cluster_id'].astype('int64')
    data['seconds'] = data['time_ms'] * 1e-3
    data.drop(columns = ['time_ms'], inplace = True)

    return data