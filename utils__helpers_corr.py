def correlations(sw, spikes, sw_metric, roi_name):
    '''
    Perform correlation between metrics of slow-wave sleep and unit spiking. 
    This requires SWS metric and spiking data to be binned into identical epochs.
    All units will be correlated with all macroelectrode channels, and p-values 
    of correlations will be FDR-corrected.

    Inputs:
    - sw = Pandas dataframe with columns 'channel', 'laterality', roi, 'epoch', and sw_metric
    - sw_metric = Column in 'sw' that will be used as the SWS metric for correlation
    - spikes = Pandas dataframe with columns 'unit_id', 'unit_laterality', 'unit_region', 'epoch', 'fr'

    Output:
    - corr = Pandas dataframe with correlation information
    '''

    import pandas as pd
    from scipy.stats import spearmanr
    #from dcor import distance_correlation
    from pingouin import distance_corr
    from statsmodels.stats.multitest import fdrcorrection

    corr = pd.DataFrame()

    for unit in spikes['unit_id'].unique():

        # Get laterality and region for the unit
        # (we need to use the .iloc[0] function because the index
        # is not reset to 0 in these slices and we don't know what
        # those index numbers will start with...)
        unit_laterality = spikes[spikes['unit_id'] == unit]['unit_laterality'].iloc[0]
        unit_region = spikes[spikes['unit_id'] == unit]['unit_region'].iloc[0]

        # Compare to all macroelectrode contacts
        ch_list = sw['channel'].unique()

        for chan in ch_list:
            
            # Get channel region/laterality (relevant for claustrum units)
            channel_laterality = sw[sw['channel'] == chan]['laterality'].iloc[0]
            channel_region = sw[sw['channel'] == chan][roi_name].iloc[0]

            # Select only the relevant SW and Spikes data
            sw_temp = sw[sw['channel'] == chan][['epoch', sw_metric]]
            spikes_temp = spikes[spikes['unit_id'] == unit][['epoch', 'fr']]

            # Full outer join (in case any epochs did not have
            # a spike or a SW, in which case one or the other 
            # would be empty and needs to be set to 0)
            corr_temp = sw_temp.merge(spikes_temp, on = 'epoch', how = 'outer')
            corr_temp = corr_temp.fillna(0)

            # Spearman's Rho (non-linear correlation; old and reliable)
            r, r_p = spearmanr(corr_temp[sw_metric], corr_temp['fr'])
            
            # # Distance Correlation (non-linear correlation; new and shiny)
            # # (The implementation from dcor does not return a p-value unless
            # #  you use one of their more complicated functions)
            # d, d_p = distance_corr(x = corr_temp[sw_metric], 
            #                        y = corr_temp['fr'],
            #                        alternative = 'two-sided',
            #                        n_boot = 1000,
            #                        seed = 42)

            # Append to dataframe
            unit_corr = pd.DataFrame({'unit_id' : [unit],
                                      'unit_region' : [unit_region],
                                      'unit_laterality' : [unit_laterality],
                                      'channel' : [chan],
                                      'channel_region' : [channel_region],
                                      'channel_laterality' : [channel_laterality],
                                      'rho' : [r],
                                      'rho_p_value' : [r_p]#,
                                      #'dcor' : [d],
                                      #'dcor_p_value' : [d_p]
                                      })
    
            corr = pd.concat([corr, unit_corr])

    corr['fdr_rho_p_value'] = fdrcorrection(corr['rho_p_value'], alpha = 0.05, method = 'indep')[1]
    #corr['fdr_dcor_p_value'] = fdrcorrection(corr['dcor_p_value'], alpha = 0.05, method = 'indep')[1]

    return(corr)