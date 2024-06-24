def epoch_sw(sw_path, bin_width, bin_list, last_bin):
    '''
    Divides slow-wave data into epochs of specified width (in seconds).

    Input:
    - sw_path = path to slow-wave data file
    - bin_width = width of bin (in seconds)
    - bin_list = list of bin cuts (in seconds)
    
    Output:
    - sw = epoched slow-wave data
    '''

    import numpy as np
    import pandas as pd

    sw = pd.read_csv(sw_path)
    sw = sw[['Start', 'End', 'Duration', 'name', 'laterality', 'region']]
    sw['seconds'] = sw['Start']

    # Bin the data with integer bin labels (pandas.cut 
    # by default will create bins open on the left)
    sw['epoch'] = pd.cut(sw['seconds'], bins = bin_list, labels = False, include_lowest = True)

    # Bin the data again, but label using the end time of each bin
    # (using bins[1:] will select every cut point except the first,
    #  which is equivalent to selecting the end of each bin)
    sw['epoch_end'] = pd.cut(sw['seconds'], bins = bin_list, labels = bin_list[1:])
    sw['epoch_end'] = sw['epoch_end'].astype(str).astype(int) 

    ## For Slow-Waves that straddle more than one Epoch, this algorithm
    ## will remove the duration of the slow-wave that extends beyond the 
    ## end of the epoch and then credit that time to the next epoch:

    # Mark SW's that end after the epoch ends
    sw['runover'] = False
    sw.loc[sw['End'] > sw['epoch_end'], 'runover'] = True

    # For run-over Epochs, calculate a new duration with End set to the Epoch End time
    sw.loc[sw['runover'] == True, 'Duration'] = sw.loc[sw['runover'] == True, 'epoch_end'] - sw.loc[sw['runover'] == True, 'Start']

    # For run-over Epochs, credit subsequent Epoch with the run-over time
    # (by creating copies of those SW's, incrementing epoch by 1, and 
    #  then re-setting the duration to the amount of run-over time)
    credit = sw[sw['runover'] == True].copy(deep = True)
    credit['epoch'] = credit['epoch'] + 1
    credit['Duration'] = credit['End'] - credit['epoch_end']

    sw = pd.concat([sw, credit])

    # Sum durations by Epoch/Unit, then divide by epoch length
    sw = sw.groupby(['epoch', 'name', 'laterality', 'region'])['Duration'].sum().reset_index()
    sw['sw_percent'] = sw['Duration'] / bin_width
    sw.sw_percent = sw.sw_percent.round(3)

    # Backfill SW epochs with 0 SW's that are missing
    # since pd.cut() will not have returned any bins
    sw_data = pd.DataFrame()

    for chan_id in sw.name.unique():
        temp = sw[sw.name == chan_id][['epoch', 'Duration', 'sw_percent']]
        temp = temp.set_index('epoch').reindex(index = np.arange(0, last_bin + 1), fill_value = 0).reset_index()

        temp['name'] = chan_id
        temp['laterality'] = sw[sw.name == chan_id]['laterality'].iloc[0]
        temp['region'] = sw[sw.name == chan_id]['region'].iloc[0]

        sw_data = pd.concat([sw_data, temp])

    sw = sw_data.copy(deep = True)

    return(sw)

def epoch_sw_2(sw_path, tmin, tmax, merge_threshold, sampling_freq, bin_list):
    '''
    Same as epoch_sw() except that (1) it will merge slow-waves that
    are within "merge_threshold" of each other and (2) it has the ability
    to use arbitrarily small bins due to the way in which it distributes
    overlapping slow wave durations (the initial epoch_sw() function can
    only distribute overlapping slow-waves over two intervals; this one 
    can distribute overlapping slow-waves over arbitrarily many intervals)

    Input:
    - sw_path = path to slow-wave data file
    - tmin = start crop time
    - tmax = end crop time
    - merge_threshold = how close SW's need to be for merging (in seconds)
    - sampling_freq = sampling frequency in Hz of the FIF file 
                      that the detected SW's were generated from
    - bin_list = list of bin cuts (in seconds)
    
    Output:
    - epoch_data = epoched slow-wave data with percent of SWS

    Inspired by: https://stackoverflow.com/questions/39710296/how-to-separate-time-ranges-intervals-into-bins-if-intervals-occur-over-multiple
    '''

    import numpy as np
    import pandas as pd
    from tqdm import tqdm
    import datetime

    # Load data
    sw = pd.read_csv(sw_path)
    sw = sw[['Start', 'End', 'Channel']]

    # Subtract/add "x/2" seconds to the Start/End times of detected
    # SW's so that you can merge SW's within of "x" seconds each other
    sw['Start'] = sw['Start'] - (merge_threshold / 2)
    sw['End'] = sw['End'] + (merge_threshold / 2)

    # Sort to determine if the SW start time in current row 
    # is greater than the maximum SW end time of all previous rows. 
    sw.sort_values('Start', inplace = True)

    sw_merged = pd.DataFrame()

    for chan_id in sw.Channel.unique():
        
        sw_temp = sw.loc[sw.Channel == chan_id].copy(deep = True)

        # Note that in (<expression>).cumsum(), the expression will be True
        # (and therefore the group number will increment) when Start > End
        # -- that is, when there is no overlap to merge. 
        # The .cumsum() will interpret True = 1 and False = 0
        sw_temp['temp_group'] = (sw_temp['Start'] > sw_temp['End'].shift().cummax()).cumsum()

        # Return minimum start time and maximum finish time from each group
        sw_temp = sw_temp.groupby('temp_group').agg({'Start' : 'min', 
                                                     'End' : 'max'}).reset_index()

        sw_temp = sw_temp[['Start', 'End']]
        sw_temp['channel'] = chan_id

        # Remove merge thresholds from start/end times
        sw_temp['Start'] = sw_temp['Start'] + (merge_threshold / 2)
        sw_temp['End'] = sw_temp['End'] - (merge_threshold / 2)
        
        sw_merged = pd.concat([sw_merged, sw_temp])

    # Convert SW start/stop times in seconds to sample numbers
    sw_merged['start'] = (sw_merged['Start'] * sampling_freq).round(0).astype('int64')
    sw_merged['stop'] = (sw_merged['End'] * sampling_freq).round(0).astype('int64')

    epoch_data = pd.DataFrame()

    for chan_id in tqdm(sw_merged.channel.unique()):

        sw_temp = sw_merged.loc[sw_merged.channel == chan_id].copy(deep = True)

        # Create a series that is as long as the total number 
        # of samples in the recording and fill it with False
        if type(tmin) == datetime.datetime:

            bool_series = pd.Series(False, index = np.arange(((tmax - tmin).seconds * sampling_freq)))
        
        else:

            bool_series = pd.Series(False, index = np.arange(((tmax - tmin) * sampling_freq)))

        for sw_row in sw_temp.itertuples():
            
            # For every SW start/stop time, set the corresponding
            # sample numbers in the boolean series to True
            bool_series[sw_row.start:sw_row.stop] = True

        # Convert the series to a dataframe with the 
        # reset index specifying the sample number
        bool_frame = bool_series.reset_index()
        bool_frame.columns = ['idx', 'sw_ratio']

        # Cut the boolean dataframe into groups (the bin_list 
        # needs to be converted from seconds to sample number)
        bool_frame['epoch'] = pd.cut(bool_frame['idx'], bins = bin_list * sampling_freq, labels = False, include_lowest = True)

        # Average boolean values over each group
        bool_frame = bool_frame.groupby('epoch').mean().reset_index()

        # Formatting
        bool_frame.drop(columns = ['idx'], inplace = True)
        bool_frame['channel'] = chan_id

        epoch_data = pd.concat([epoch_data, bool_frame])
    
    return(epoch_data)

def epoch_spikes(spike_path, bin_width, bin_list, last_bin):
    '''
    Divides spike data into epochs of specified width (in seconds).

    Input:
    - spike_path = path to slow-wave data file
    - bin_width = width of bin (in seconds)
    - bin_list = list of bin cuts (in seconds)
    
    Output:
    - spikes = epoched spike data
    '''

    import numpy as np
    import pandas as pd

    spikes = pd.read_csv(spike_path)
    spikes = spikes[['seconds', 'unit_id', 'unit_laterality', 'unit_region']]

    # Convert each spike time into an assigned Epoch
    spikes['epoch'] = pd.cut(spikes['seconds'], bins = bin_list, labels = False, include_lowest = True)

    # DEPRECATED alternative for SpikeInterface
    #spikes = spikes[(spikes['subject'] == subject) & (spikes['status'] == 'accept')] 
    #spikes = spikes[['time', 'unit', 'laterality', 'region']] 
    #spikes['seconds'] = spikes['time'] / 30000 

    # Average firing rate by unit & epoch
    # Note 1: We use reset_index() in order to use the
    # row names to a unique spike identified that we
    # can use for counting spikes per epoch.
    # Note 2: We don't need to group by laterality 
    # or region, but we do so redundantly in order to 
    # save their values in the final output.
    spikes = spikes.reset_index().groupby(['unit_laterality', 'unit_region', 'unit_id', 'epoch'])['index']
    spikes = spikes.count().reset_index()

    # Calculate average firing rate by epoch (30 seconds) using 
    # 'index', which is the total number of spikes in the epoch:
    spikes['fr'] = spikes['index'] / bin_width
    spikes.fr = spikes.fr.round(3)

    # Backfill spikes epochs with 0 spikes that are missing
    # since pd.cut() will not have returned any bins for them
    spike_data = pd.DataFrame()

    for unit in spikes.unit_id.unique():
        temp = spikes[spikes.unit_id == unit][['epoch', 'fr']]
        temp = temp.set_index('epoch').reindex(index = np.arange(0, last_bin + 1), fill_value = 0).reset_index()

        temp['unit_id'] = unit
        temp['unit_laterality'] = spikes[spikes.unit_id == unit]['unit_laterality'].iloc[0]
        temp['unit_region'] = spikes[spikes.unit_id == unit]['unit_region'].iloc[0]

        spike_data = pd.concat([spike_data, temp])

    spikes = spike_data.copy(deep = True)

    return(spikes)

def make_pca(epoch_data, prefix, num_pc=1, seed=42):
    '''
    Parameters:
        - epoch_data = pandas dataframe with columns [time, group, value]
        - prefix = character string to append to PC names
        - num_pc = number of PC's to return (top num_pc)
        - wide_output = return output as wide output (versus long)
    Returns:
        - pca_data = pandas dataframe with columns [feat_PC, PC, time, value]
    '''

    import pandas as pd
    from scipy.stats import zscore
    from sklearn.decomposition import PCA

    # Pivot wider
    data = epoch_data.pivot(index='group',
                            columns='time',
                            values='value')

    # Transpose matrix
    data = data.T

    # Set PCA parameters
    pca = PCA(n_components=num_pc,
              random_state=seed)

    # Train PCA using your data
    pca.fit(data)

    # Get total variance explained by each PC
    explainer = pca.explained_variance_ratio_ * 100
    print(explainer)

    # Fit PCA to the original data
    pca_data = pca.transform(data)

    # Reformat PCA data
    pca_data = pd.DataFrame(pca_data.T).reset_index()

    pca_data = pd.melt(pca_data,
                       id_vars=['index'],
                       var_name='time',
                       value_name='pc_value')

    pca_data.columns = ['PC', 'time', 'value']
    pca_data['PC'] = pca_data['PC'] + 1
    pca_data['feat_PC'] = prefix + '_' + pca_data['PC'].astype('str')

    # Z-score the values PC-wise
    pca_data['value'] = pca_data.groupby(['PC', 'feat_PC'])[
        'value'].transform(zscore)

    return(pca_data)