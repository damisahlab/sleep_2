import pandas as pd

def robust_zscore(data):
    '''
    Calculates the robust z-score, which uses the median instead of the
    mean and the median absolute deviation instead of the standard deviation.
    This is a non-parametric version of a z-score and is less vulnerable to
    outlier values.

    Input: 1d numpy series
    Output: 1d numpy series
    '''

    import numpy as np
    from scipy.stats import median_abs_deviation

    rz_median = np.median(data) 
    rz_MAD = median_abs_deviation(data)

    data = ( 0.6745 * (data - rz_median) ) / rz_MAD
    
    return(data)

def ied_event_rejection(event_data, ied_data, rej_int):
    rejected_events = []

    for ch_name in event_data['Channel'].unique():

        # Select IED's for the selected channel
        left_bound = ied_data[ied_data['name'] == ch_name]['start_time']
        right_bound = ied_data[ied_data['name'] == ch_name]['stop_time']
        
        # Create an array of closed intervals for every IED in the channel
        ied_intervals = pd.arrays.IntervalArray.from_arrays(left = left_bound, 
                                                            right = right_bound, 
                                                            closed = 'both')

        # Iterate over every detected SW in the channel
        data_temp = event_data.loc[event_data['Channel'] == ch_name]
        for index, row in data_temp.iterrows():
            
            # Create a closed interval for the SW
            int_temp = pd.Interval(left = row['Start'] - rej_int, 
                                   right = row['End'] + rej_int,
                                   closed = 'both')
            
            # Check if any IED intervals overlap with the SW
            any_overlap = ied_intervals.overlaps(int_temp)
            any_overlap = any_overlap.any()

            if any_overlap == True:
                rejected_events.append(row['ID'])

    print('Number of rejected events: ', len(rejected_events))
    return(rejected_events)

def epoch_rejection(event_data, epoch_data):
    # No need to specify a rejection interval around the epoch
    # since, unlike the IED detection, the rejection interval
    # padding was built into the epoch length in the epoch script.

    rejected_events = []

    # Select IED's for the selected channel
    left_bound = epoch_data['start_time']
    right_bound = epoch_data['stop_time']
    
    # Create an array of closed intervals for every IED in the channel
    epoch_intervals = pd.arrays.IntervalArray.from_arrays(left = left_bound, 
                                                          right = right_bound, 
                                                          closed = 'both')


    for index, row in event_data.iterrows():
        
        # Create a closed interval for the SW
        int_temp = pd.Interval(left = row['Start'], 
                               right = row['End'],
                               closed = 'both')
        
        # Check if any IED intervals overlap with the SW
        any_overlap = epoch_intervals.overlaps(int_temp)
        any_overlap = any_overlap.any()

        if any_overlap == True:
            rejected_events.append(row['ID'])

    print('Number of rejected events: ', len(rejected_events))
    return(rejected_events)

# Mark SW's that are isolated from other SW's (in each channel)
def mark_isolation(event_data, rej_int):
    isolated_events = []

    for ch_name in event_data['Channel'].unique():
        # Select start/end times of SW's in the channel
        left_bound = event_data[event_data['Channel'] == ch_name]['Start']
        right_bound = event_data[event_data['Channel'] == ch_name]['End']
        
        # Create an array of closed intervals for every IED in the channel
        sw_intervals = pd.arrays.IntervalArray.from_arrays(left = left_bound, 
                                                           right = right_bound, 
                                                           closed = 'both')

        # Iterate over every detected SW in the channel
        data_temp = event_data.loc[event_data['Channel'] == ch_name]
        for index, row in data_temp.iterrows():
            
            # Create a closed interval for the SW
            int_temp = pd.Interval(left = row['Start'] - rej_int, 
                                   right = row['End'] + rej_int,
                                   closed = 'both')
            
            # Compute overlaps of this interval with other SW's
            overlap_array = sw_intervals.overlaps(int_temp)

            # Count the number of overlaps with other SW's
            overlaps = overlap_array.tolist().count(True)

            # Consider the SW isolated if there is only
            # one other overlapping SW (the comparison list
            # includes the SW in question, so the minimum is 1):
            if overlaps == 1:
                isolated_events.append(row['ID'])

    print('Number of isolated events: ', len(isolated_events))
    return(isolated_events)

# Mark SW's that have large peak-to-peak amplitudes
def mark_bigpeaks(event_data, threshold, absolute = False):
    from scipy.stats import zscore
    #print('Distribution of positive peaks:\n', event_data['ValPosPeak'].describe())
    #print('Distribution of negative peaks:\n', event_data['ValNegPeak'].describe())

    # If the amplitude threshold is absolute, select SW's with peak-to-peak
    # amplitudes above the specified uV threshold:
    if absolute == True:
        big_peaks = event_data[event_data['PTP'] > threshold, 'ID']

    # If the amplitude threshold is relative, select SW's with peak-to-peak
    # amplitudes greater than the specified number of SD's for PTP amplitudes
    # within its channel; do this for every channel:
    else:
        event_data['zscore'] = event_data.groupby('name')['PTP'].transform(zscore)
        big_peaks = event_data[event_data['zscore'] > threshold, 'ID']

    print('Number of SWs with big peaks: ', len(big_peaks))
    return(big_peaks)

def hilbert_powerphase(data, lower, upper, njobs = 1):
    """
    This function extracts the instantaneous power and phase
    angle of a specified frequency band by using the Hilbert 
    Transform on a time series of voltages.

    Input:
        data = MNE Raw object
        lower = lower bandpass frequency
        upper = upper bandpass frequency
        njobs = number of parallel jobs for bandpass

    Output: 
        data = Pandas long dataframe with Hilbert-derived power 
               and phase angle for the specified frequency band
    """

    import numpy as np

    # Bandpass
    data.filter(l_freq = lower, 
                h_freq = upper,
                n_jobs = njobs)

    # Hilbert Transform returns the
    # analytic signal (complex)
    data.apply_hilbert(envelope = False)

    ### DEPRECATED (do this later once in Pandas)
    # Extract Phase Angle
    #data.apply_function(fun = np.angle, picks = 'all')

    # Convert to Pandas dataframe
    # (this MNE method will automatically scale 'eeg'
    #  electrodes by 1e6 when exporting, so we need to
    #  specify not to do this if we have only phase angles) # deprecated
    data = data.to_data_frame(long_format = True, 
                              time_format = None, # float values in seconds
                              scalings = dict(eeg = 1e6))

    # Extract power and phase angle from Hilbert analytic values
    data['power'] = np.abs(data['value'])**2
    data['phase'] = np.angle(data['value'])

    return(data)

def hilbert_envelope(data, lower, upper, njobs = 1):
    """
    See hilbert_powerphase(). 
    This function extracts the envelope instead.
    """

    import numpy as np

    # Bandpass
    data.filter(l_freq = lower, 
                h_freq = upper,
                n_jobs = njobs)

    # Hilbert Transform returns the envelope
    data.apply_hilbert(envelope = True)

    # Convert to Pandas dataframe
    data = data.to_data_frame(long_format = True, 
                              time_format = None, # float values in seconds
                              scalings = dict(eeg = 1e6))

    data = data.rename(columns={'value' : 'envelope'})

    return(data)

def welch_psd(data, chan_names, sampling_freq, freq_min = 0.3, freq_max = 85, n_jobs = -2):
    """
    Converts an numpy array to a power spectral density using
    the Welch method.

    Input:
    - data = numpy array of (channels, times)

    Output:
    - tfr = long format pandas dataframe with columns...
            channel, frequency, power, log_power
    """

    import numpy as np
    import xarray as xr
    from mne.time_frequency import psd_array_welch
    
    # Compute Power Spectral Density (PSD) using Welch's method
    psd = psd_array_welch(x = data,
                          sfreq = sampling_freq,
                          fmin = freq_min,
                          fmax = freq_max,
                          n_jobs = n_jobs)

    # Extract data from PSD output
    psd_data = psd[0]
    freq_labels = pd.Series(psd[1]).astype('int64').to_list()

    # Convert to Xarray as an intermediate step in
    # getting data into Pandas long (2d) format:
    tfr = xr.DataArray(psd_data,
                       dims = ('channel', 'frequency'),
                       coords = {'channel' : chan_names,
                                 'frequency' : freq_labels})

    tfr = tfr.to_dataframe(name = 'power').reset_index()

    # Calculate log-power
    tfr['log_power'] = 10 * np.log10(tfr['power'])

    return(tfr)