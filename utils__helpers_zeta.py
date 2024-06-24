def ZETA_test(spikes, events, regions, sides, ipsilateral_only, time_window):
    '''
    Perform a ZETA test between times of spikes and times of sleep events.

    spikes = Pandas dataframe containing the following columns:
        'seconds' = spike time in seconds from recording start (float)
        'region' = anatomic region of spike (string)
        'laterality' = laterality of spike (string)
        'unit' = unique ID/name of single unit

    events = Pandas dataframe containing the following columns:
        'seconds' = sleep event time in seconds from recording start (float)
        'region' = anatomic region of sleep event (string)
        'laterality' = laterality of sleep event (string)
        'name' = unique ID/name of macroelectrode channel

    regions = list of regions to loop through

    sides = list of lateralities to loop through

    contralateral = boolean of whether or not to compare each unit to the region
                    that is on the contralateral side to it (if True, contralateral;
                    if False, it will compare ipsilateral)

    time_window = window of time around the sleep event to take spike times
                  for calculation of the ZETA score (integer or float)

    Explanation of some ZETA Test parameters:
        - dblUseMaxDur = Time window for which ZETA collects spikes 
          (this window will be centered, so half before and half after event)
        - intResampNum = Number of resamplings (default is 100)
        - intLatencyPeaks = Which latency peaks to return
            1 = ZETA -> Latency of ZETA
            2 = -ZETA -> Latency of largest z-score with inverse sign to ZETA
            3 = Peak -> Peak time of instantaneous firing rate
            4 = Onset time of above peak, defined as the first crossing of peak half-height
        - tplRestrictRange = Range within which to restrict onset/peak latencies (default [-inf inf])
    '''

    import numpy as np
    import pandas as pd
    import itertools as it
    from zetapy import getZeta
    from statsmodels.stats.multitest import fdrcorrection

    zeta_data = pd.DataFrame()
    ifr_data = pd.DataFrame()

    for region, side in it.product(regions, sides):

        # Select list of relevant units
        units = spikes[(spikes['unit_region'] == region) & (spikes['unit_laterality'] == side)]['unit_id'].unique()

        # Select list of relevant channels
        if region == 'CLA':
            if ipsilateral_only == True:
                channels = events[events['laterality'] == side]['name'].unique()
            elif ipsilateral_only == False:
                channels = events['name'].unique()
        else:
            if ipsilateral_only == True:
                channels = events[(events['region'] == region) & (events['laterality'] == side)]['name'].unique()
            elif ipsilateral_only == False:
                channels = events[events['region'] == region]['name'].unique()

        # If there is at least one unit and one macro channel
        # for the specified region-laterality pair:
        if (len(units) > 0) & (len(channels) > 0):

            # Loop through every unit-channel pair
            for unit in units:

                print(side, region, unit)

                for channel in channels:

                    # Get spike/sleep event times for the unit and channel
                    spike_times = spikes[spikes['unit_id'] == unit]['seconds']
                    event_times = events[events['name'] == channel]['seconds']

                    # Set new event times (ZETA only looks at the window after
                    # the event start, so we need to subtract the window time
                    # from the event start to create a window around the event)
                    subtracted_event_times = event_times - time_window
                    double_window = time_window * 2

                    # Run ZETA test for unit-channel pair
                    pval, lat, zdic, ifr = getZeta(arrSpikeTimes = spike_times, 
                                                   arrEventTimes = subtracted_event_times,
                                                   dblUseMaxDur = double_window,
                                                   intResampNum = 100,
                                                   intLatencyPeaks = 4,
                                                   boolReturnZETA = True,
                                                   boolReturnRate = True)
                    
                    # Plot the spike raster plot
                    #plot.eventplot(zdic['vecSpikeT'])
                    #plot.show()

                    # Plot the instantaneous firing rate (IFR)
                    #plot.plot(ifr['vecT'], ifr['vecRate'])
                    #plot.show()

                    # You need to wrap the dictionary in a list [] before passing it
                    # to the Pandas DF constructor or it will error because all the 
                    # values are scalar and therefore do not have an index:
                    zeta_loop = pd.DataFrame([{'unit_id' : unit,
                                               'unit_region' : region,
                                               'unit_side' : side,
                                               'channel' : channel,
                                               'channel_region' : events[events['name'] == channel]['region'].iloc[0],
                                               'channel_side' : events[events['name'] == channel]['laterality'].iloc[0],
                                               'zeta_time' : lat[0],
                                               'peak_time' : lat[2],
                                               'p_value' : pval}])
                    
                    # If you have too few spikes/events to calculate a latency,
                    # then the vector of spike times will be empty and error:
                    if np.isnan(lat[0]) == False:

                        # Return instantaneous firing rate data generated by ZETA
                        # (you can get spike times for raster plots from zdic['vecSpikeT'])
                        ifr_loop = pd.DataFrame({'unit_id' : unit,
                                                 'channel' : channel,
                                                 'seconds' : ifr['vecT'],
                                                 'ifr' : ifr['vecRate']})
                        
                        ifr_data = pd.concat((ifr_data, ifr_loop))

                    zeta_data = pd.concat((zeta_data, zeta_loop))

        else:

            print('No units found for', region, side)

    # FDR correction of p-values (Benjamini-Hochberg)
    zeta_data['fdr_p_value'] = fdrcorrection(zeta_data['p_value'], alpha = 0.05, method = 'indep')[1]

    return(zeta_data, ifr_data)

def spike_epochs(spikes, events, regions, sides, ipsilateral_only, sleep_event, time_window, bin_size):
    '''
    Epoch spike times based on a time window around sleep events. See the documentation
    of the ZETA_test() function for an explanation of common parameters.

    sleep_event = type of sleep event to analyze (this dictates the
                  meta-data that is extracted and saved for sleep events)
    
    bin_size = size of bins (in seconds) into which spike times are categorized
    '''

    import numpy as np
    import pandas as pd
    import itertools as it
    
    data = pd.DataFrame()

    for region, side in it.product(regions, sides):

        # Select list of relevant units
        units = spikes[(spikes['unit_region'] == region) & (spikes['unit_laterality'] == side)]['unit_id'].unique()

        # Select list of relevant channels
        if region == 'CLA':
            if ipsilateral_only == True:
                channels = events[events['laterality'] == side]['name'].unique()
            elif ipsilateral_only == False:
                channels = events['name'].unique()
        else:
            if ipsilateral_only == True:
                channels = events[(events['region'] == region) & (events['laterality'] == side)]['name'].unique()
            elif ipsilateral_only == False:
                channels = events[events['region'] == region]['name'].unique()

        # Loop for every unit-channel pair
        for unit in units:
            for channel in channels:
                
                # Get channel meta-data
                channel_region = events[events['name'] == channel]['region'].iloc[0]
                channel_side = events[events['name'] == channel]['laterality'].iloc[0]

                # Get spike/sleep event times for the unit and channel
                spike_times = spikes[spikes['unit_id'] == unit]['seconds']
                event_times = events[events['name'] == channel]['seconds']

                # Select a time window around every sleep event
                event_windows = pd.arrays.IntervalArray.from_arrays(left = event_times - time_window,
                                                                    right = event_times + time_window,
                                                                    closed = 'both')
                
                # Select spikes contained within each sleep event time window
                for index, event_window in enumerate(event_windows):

                    # Generate booleans for spike times within a window
                    selected_spikes = spike_times.between(left = event_window.left, 
                                                          right = event_window.right, 
                                                          inclusive = 'both')

                    # Use booleans to select the actual spike times
                    event_spikes = spike_times[selected_spikes]

                    # Reset spike times relative to the sleep event time
                    event_spikes = event_spikes - (event_window.left + time_window)

                    # Create an output dictionary to feed to Pandas
                    output_dict = {'unit_id' : unit,
                                   'unit_region' : region,
                                   'unit_side' : side,
                                   'channel' : channel,
                                   'channel_region' : channel_region,
                                   'channel_side' : channel_side,
                                   'event_number' : index,
                                   'spike_time' : event_spikes}
                    
                    # Modify output dictionary depending on sleep event
                    if sleep_event == 'Slow-Wave' or sleep_event == 'K-Complex':

                        event_dict = {'negative_peak' : events[events['name'] == channel]['ValNegPeak'].iloc[index], 
                                      'peak_to_peak' : events[events['name'] == channel]['PTP'].iloc[index], 
                                      'slope' : events[events['name'] == channel]['Slope'].iloc[index]}

                    elif sleep_event == 'Spindle':

                        event_dict = {'amplitude' : events[events['name'] == channel]['Amplitude'].iloc[index], 
                                      'relative_power' : events[events['name'] == channel]['RelPower'].iloc[index], 
                                      'oscillations' : events[events['name'] == channel]['Oscillations'].iloc[index]}
                    
                    output_dict.update(event_dict)
                    
                    # Concatenate data
                    loop_data = pd.DataFrame(output_dict)
                    data = pd.concat((data, loop_data))

    # Add a binned time column
    # (labels need to be one fever than the number of bin edges, hence '1:')
    bin_list = np.arange(-time_window, time_window, bin_size)
    data['time_bin'] = pd.cut(data['spike_time'], bins = bin_list, labels = bin_list[1:])

    return(data)

def response_type(zeta, epochs):
    '''
    This function will define neurons as bursters or suppressors (or none).

    Input: 
        - zeta = Output dataframe from ZETA_test()
        - epochs = Output dataframe from spike_epochs()
    
    Output:
        - output = Modified 'zeta' dataframe with response type appended
    '''

    import numpy as np
    import pandas as pd

    # Set a boolean indicating if spike is before or after the sleep event
    epochs['spike_bin'] = np.where(epochs['spike_time'] <= 0, 'before', 'after')

    # Count the number of spikes before/after, stratified by unit-channel pair
    epochs = epochs.groupby(['unit_id', 'channel', 'spike_bin'])['spike_time'].count().reset_index()

    # Cast dataframe from wide to long so that before/after counts are side-by-side
    epochs = epochs.pivot(index = ['unit_id', 'channel'], columns = 'spike_bin', values = 'spike_time').reset_index()

    # Set unit response type based on how before/after spike counts compare
    epochs['response'] = 'none'
    epochs.loc[epochs['before'] > epochs['after'], 'response'] = 'burst'
    epochs.loc[epochs['before'] < epochs['after'], 'response'] = 'suppress'

    # Merge with ZETA to store response type
    epochs = epochs[['unit_id', 'channel', 'response']]
    output = zeta.merge(epochs, on = ['unit_id', 'channel'])

    return(output)