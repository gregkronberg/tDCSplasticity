'''
'''
import numpy as np
import os
import pandas as pd
from matplotlib import pyplot as plt
from scipy import stats
from collections import OrderedDict 
import copy
import glob
import pickle
from sklearn.cross_decomposition import CCA
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import pdb
import inspect
import math
import itertools
from scipy.optimize import curve_fit

def _initialize_column(df, col, val=None):
    '''
    '''
    if type(col)==str:
        col=[col]
    for col_name in col:
        if col_name not in df:
            df[col_name]=val
    return df

def exp_fit(x, y, p0=[1,-1]):
    '''
    '''
    def exp_func(x, a, b):
        return a*np.exp(b*x)

    popt, pcov = curve_fit(exp_func, x, y, p0)

    return popt

def _build_timeseries_anova_df(df, array_column, columns):
    '''
    '''
    df_new=pd.DataFrame()
    for row_i, row in df.iterrows():
        timeseries = row[array_column]
        print type(timeseries)
        if not type(timeseries)==float or (type(timeseries)==float and not np.isnan(timeseries)):
            if type(timeseries)==float:
                timeseries = np.array([timeseries])
            df_temp = pd.DataFrame(timeseries,)
            df_temp.columns = [array_column]
            print df_temp.columns
            df_temp[array_column+'_i']=df_temp.index.values
            for column in columns:
                df_temp[column]=row[column]
            df_new = df_new.append(df_temp)

    return df_new

def _build_ttest_dict(df, figdf, variable, conditions):
    '''
    '''
    stats_df = _build_timeseries_anova_df(df=df, array_column=variable, columns=conditions)
    stats_df = stats_df.set_index(conditions)

    traces = figdf.index.get_level_values('trace')
    trace_combos = list(itertools.product(traces,traces))
    stats_dict = {}
    print sorted(stats_df.index.unique())
    for combo in trace_combos:
        # print combo[0]
        if combo[0] in stats_df.index.unique() and combo[1] in stats_df.index.unique():
            print combo

            df_0 = stats_df.loc[combo[0]]
            df_1 = stats_df.loc[combo[1]]
            if combo[0] not in stats_dict:

                stats_dict[combo[0]]={
                    'n':df_0.shape[0],
                    'std':np.std(df_0[variable]),
                    'sem':stats.sem(df_0[variable]),
                    'mean':np.mean(df_0[variable]),
                }
            if combo[1] not in stats_dict[combo[0]]:
                stats_dict[combo[0]][combo[1]]= stats.ttest_ind(df_0[variable],df_1[variable]).pvalue

    return stats_dict




# def _to_1sigfig(num, minimum=True):
#     '''
#     '''
#     print 'num',num
#     if type(num)==np.ndarray or type(num)==list or type(num)==tuple:
#         new_val = np.zeros(len(num))
#         sigfigs=[]
#         for i, val in enumerate(num):
#             if val!=0:
# def _get_burst_spiketimes(data, t_window, ):
#     '''
#     '''
#     # zero out artifact
#     #-------------------
#     # copy data
#     data_temp = copy.deepcopy(data)
#     # zero out samples that aren't going to be used
#     data_temp[:t_window[0],:] = np.zeros(data_temp[:t_window[0],:].shape)
#     data_temp[t_window[1]:,:] = np.zeros(data_temp[t_window[1]:,:].shape)
#     # reshape to samples x bursts
#     data_temp = np.reshape(data_temp, (-1, 15), order='F')

#     # create time vector in ms
#     #--------------------------
#     fs_ms = 10. # samples per ms
#     t0 = np.arange(0, data_temp.shape[0]/fs_ms,1/fs_ms)
#     # convert to array with same shape as data (each column is time vector for corresponding burst in data)
#     t_array = np.tile(t0, (data_temp.shape[1],1)).T

#     # center of mass
#     #------------------
#     center_of_mass = np.sum(np.multiply(data_temp, t_array), axis=0)/np.sum(data_temp, axis=0)
#     # weighted_average = np.mean(np.multiply(data_temp, t_array), axis=0).squeeze()/np.mena

#     return center_of_mass  sigfigs.append(int(math.floor(math.log10(abs(val)))))
#             else:
#                 sigfigs.append(0)

#         min_sigfig = min(sigfigs)
#         print min_sigfig
#         for i, val in enumerate(num):
#             if minimum:
#                 new_val[i] = np.round(val, -min_sigfig)
#             else:
#                 new_val[i] = np.round(val, -sigfigs[i])

#     elif type(num)==float or type(num)==np.float64:
#         if num!=0:
#             new_val = np.round(num, -int(math.floor(math.log10(abs(num)))))
#         else:
#             new_val=0

#     print new_val
#     return new_val

def _add_escape(constraints):
    '''add escape quotes to query string based on constraint dictionary
    ==Args==
    -constraints : dictionary specifying constraints to be passed to dataframe as query
            ~can be organized as {constraint column in dataframe}{(condition values from conditions list)}[constraint logical, constraint value] for constraints that are specific to different experimental conditions
            ~can be organized as {constraint column}[constraint logical, constraint value] if the constraint is to be aplied to all conditions
    ==Out==
    -constraints : updated constraints dictionary, any constraint values that are strings have escape quotes added for use in query in string

    ==Updates==
    ==Comments==
    '''
    # add escape quotes to df entries that are strings. needed to search for strings with query function
    # contstraints_all
    for k0, v0 in constraints.iteritems():
        if type(v0)==dict:
            for k1, v1 in v0.iteritems():
                if type(v1[1])==str and '\'' not in v1[1]:
                    constraints[k0][k1][1] = '\'' + v1[1] + '\''
        else:
            if type(v0[1])==str and '\'' not in v0[1]:
                constraints[k0][1] = '\'' + v0[1] + '\''

    return constraints

def _build_constraint_query(constraints_spec, constraints_all, group_key):
    '''
    ==Args==
    ==Out==
    ==Updates==
    ==Comments==
    ''' 
    print 'building constraint query to sort df for:', group_key
    strings = []
    for k0, v0 in constraints_spec.iteritems():
        for k1, v1 in v0.iteritems():
            if type(k1)==str:
                k1=(k1,)

            constraint_vals = set(k1)
            if len(set(constraint_vals) - set(group_key))==0:

                current_string = '{}{}{}'.format(k0, v1[0], v1[1])
                strings.append(current_string)
    constraint_query = '&'.join(strings)


    if len(constraints_all.keys())!=0:
        if any(val not in constraints_spec.keys() for val in constraints_all.keys() ):
            constraint_query_all =  '&'.join([
                '{}{}{}'.format(k,v[0],v[1])
                for k, v in constraints_all.iteritems() if k not in constraints_spec
            ])

            constraint_query +='&' + constraint_query_all 


    return constraint_query

def _align_row_data(self, df, column, align_on='induction_block_0', i_pre=20, i_post=60, include_align=True):
        ''' create columns with normalized slopes that are aligned to the induction block and the final amount of ltp
        ==Args==
        -group_df : group dataframe. must contain columns 'induction_block_0', 'slopes_norm'

        ==Out==
        -group_df : group dataframe. new columns: 'slopes_norm_aligned_0' and 'ltp_final'

        ==Update==
        -group_df.slopes_norm_aligned_0 : normalized slopes are aligned to the first induction block, so that each entry is an array of length 80 (the first 20 samples are pre-indction, the last 60 samples are post-induction)
        -group_df.ltp_final : mean of the last ten normalized slopes (min 51-60) for each trace

        ==Comments==
        '''
        def _get_aligned_data(row):
            '''for each row align slopes data to the induction block
            ==Args==
            -row : row in a pandas dataframe
            ==Out==
            -aligned_data :  a list, where the first element is the aligned slopes
            ==Comments==
            -when fed to the apply function in pandas, a series is returned, where each entry contains the corresponding aligned data
            '''
            # print 'induction block 0:', row['induction_block_0']
            # get indices to align data
            indices = range( int( row[ align_on])-i_pre, int(row[align_on])) + range( int( row[align_on]), int( row[ align_on])+ i_post)
            # get aligned data
            aligned_data = [row[column][indices]]
            return aligned_data

        aligned_column = group_df.apply(lambda row: _get_ltp_final(row), axis=1)
        return aligned_column

def _2array(series, remove_nans=True, remove_nans_axis=0, list_index=0, array_funcs=[], array_func_kws=[]):
    '''
    '''
    
    series = series.dropna()
    series_list = series.tolist()

    if len(series_list)>0 and type(series_list[0])==list:
    #     print 'list to array?'
        series_list = [item[list_index] for item in series_list if type(item)==list and len(item)>list_index]
    #     # print series_list
    #     print np.array(np.array(series_list))

    array = np.array(series_list).squeeze()

    

    # array = np.array(series.tolist()).squeeze()
    # print np.isnan(array)
    # array = array[~np.isnan(array).any(axis=remove_nans_axis)]
    return array

def _sortdf(df, conditions, constraints_spec, constraints_all):
    '''
    '''
    print 'sorting dataframe'
    # add escape quotes to df entries that are strings. needed to search for strings with query function
    # contstraints_all
    constraints_spec = _add_escape(constraints_spec)
    constraints_all = _add_escape(constraints_all)
    # get sorting indices
    # groups is a dictionary with keys that correspond to combinations of conditions (e.g. ('anodal', 'TBS', 'nostim')) and values that are the indices that meet those conditions
    groups = df.groupby(conditions).groups

    # use groupby to sort dataframes by conditions
    df_sorted = {}
    for group_key, group_index in groups.iteritems():
        
        # group_key must be a tuple
        if type(group_key)==str:
            group_key=(group_key,)

        # build constraint query string out of constraint conditions
        constraint_query = _build_constraint_query(constraints_spec, constraints_all, group_key)

        # sort the data that match the experimental conditions and constraints
        df_sorted[group_key] = df.loc[group_index].query(constraint_query)

    return df_sorted

def _load_group_data( directory='', filename='slopes_df.pkl', df=True):
    """ Load group data from folder
    
    ===Args===
    -directory : directory where group data is stored including /
    -filename : file name for group data file, including .pkl
                -file_name cannot contain the string 'data', since this string is used t search for individual data files
    -df : boolean, if true group data is assumed to be a pandas dataframe, otherwise it is assumed to be a nested dictionary

    ===Out===
    -group_data  : if df is True, group_data will be a pandas dataframe.  if no group data is found in the directory an empty dataframe is returned.  if df is False, group_data is returned as a nested dictionary, and if no data is found an empty dictionary is returned

    ===Updates===
    -none

    ===Comments===
    """
    
    # all files in directory
    files = os.listdir(directory)

    # if data file already exists
    if filename in files:
        print 'group data found:', filename

        # if stored as dataframe
        if df:
            # load dataframe
            group_data=pd.read_pickle(directory+filename)
            print 'group data loaded'
        # if stored as dictionary
        else:
            # load dictionary
            with open(directory+filename, 'rb') as pkl_file:
                group_data= pickle.load(pkl_file)
            print 'group data loaded'

    # otherwise create data structure
    else:
        print 'no group data found'
        
        # if dataframe
        if df:
            # create empty dataframe
            group_data = pd.DataFrame()
        # else if dicitonary
        else:
            group_data= {}

    return group_data 

def _remove_outliers(time_series, ind_idx, time_window, std_tol, include_ind=False):
    '''# FIXME docs, comments
    '''
    ts = time_series
    # preallocate
    ts_smooth = copy.copy(ts) 
    ind_idx_array = np.array(ind_idx)
    # iterate over values
    for i, val in enumerate(ts):
        # if the current index is an induction, skip it 
        if i in ind_idx:
            continue

        # get indices of inductions that were before the current index
        ind_cross = np.where(ind_idx_array<=i)[0]
        # if there have not been any inductions yet
        if len(ind_cross)==0:
            # set most recent induction to 0
            ind =[0]
            # otherwise, get most recent induction
        else:
            ind = ind_idx_array[ind_cross]

        # if the current index is within a time_window distanc of the most recent induction
        if include_ind:
            if i-ind[-1]+1 <= time_window:
                vals_i = range(ind[-1]+1,(ind[-1]+1+time_window) )
                current_i = i-ind[-1]+1
                # use forward values to fill out time window
                vals = ts[vals_i]
                # vals = ts[ind[-1]+1:(ind[-1]+1+time_window)]
            # otherwise use the last time window values to fidn outliers
            else:
                vals_i =range(i-time_window,i)
                # current_i=-1
                vals = ts[vals_i]
                current_i =len(vals)-1
        else:
            if i-ind[-1] +1 < time_window:
                if ind[-1]+time_window >= len(ts)+1:
                    vals_i = range(ind[-1],len(ts))
                else:
                    vals_i = range(ind[-1],(ind[-1]+time_window))
                current_i = i-ind[-1]
                # use forward values to fill out time window
                vals=ts[vals_i]
                # vals = ts[ind[-1]:(ind[-1]+time_window)]
            # otherwise use the last time window values to fidn outliers
            else:
                vals_i = range(i-time_window+1,i+1)
                # current_i = -1
                vals=ts[vals_i]
                current_i =len(vals)-1
                # vals = ts[i-time_window+1:i]

        # mean and standard deviation of values
        other_vals = [temp for temp_i, temp in enumerate(vals) if temp_i != current_i]
        vals_std = np.std(other_vals)
        vals_mean = np.mean(other_vals)

        # if value is outside the tolerance set by std_tol, set value to mean of all values in time_window 
        if abs(val-vals_mean) < std_tol*vals_std:
            ts_smooth[i] = val
        else: 
            ts_smooth[i] = vals_mean

    return ts_smooth

def _build_conditions_dict(preprocessed, path):
    ''' FIXME add doc
    '''
    pre = preprocessed
    # dictionary for current data file to be converted to dataframe
    current_dict = {}
    # df = pd.DataFrame()
    # iterate over pathways
    # for path in pre['path_blocks']:

    # FIXME how to handle more than 2 pathways
    # get name of other pathway
    n_paths = len(pre['path_blocks'].keys())
    try:
        path_other = [temp for temp in pre['path_blocks'] if temp != path][0]
    except IndexError:
        path_other = 'None'


    # iterate over inductions
    for ind_i, ind in enumerate(pre['induction_info']):

        # store induction info as entries in current_dict. the induction number is indicated at the end of each dictionary key, e.g. 'induction_0'

        # FIXME indices need to be realigned to each pathway after the pathways are separated.  use ind_idx from _get_slopes_probe
        # store the full induction info dictionary
        current_dict['induction_'+str(ind_i)]=[ind]

        # induction pattern, e.g. TBS, weak5Hz, nostim
        current_dict['induction_pattern_'+str(ind_i)] = ind[path]['protocol']
        try:
            current_dict['induction_pattern_other_'+str(ind_i)] = ind[path_other]['protocol']
        except KeyError:
            current_dict['induction_pattern_other_'+str(ind_i)] = 'None'

        
        # location stimulating electrode, e.g. apical, basal
        current_dict['induction_location_'+str(ind_i)] = ind[path]['location']
        try:
            current_dict['induction_location_other_'+str(ind_i)] = ind[path_other]['location']
        except KeyError:
            current_dict['induction_location_other_'+str(ind_i)] = 'None'


        # induction block
        current_dict['induction_block_'+str(ind_i)] = pre['ind_idx'][path][ind_i]
        # ind['induction_block']

        # electric field magnitude (in V/m)
        current_dict['field_mag_'+str(ind_i)] = ind['field_magnitude']

        # electric field polarity (anodal, cathodal, control)
        current_dict['field_polarity_'+str(ind_i)] = ind['polarity']


    current_dict['comments'] = [pre['comment_dict']]
    current_dict['induction']=[pre['induction_info']]
    current_dict['path']=[path]
    current_dict['filename']=pre['slice_info']['filename']
    current_dict['name']=pre['slice_info']['name']
    current_dict['date']=pre['slice_info']['date']
    current_dict['path_other'] = path_other
    current_dict['age'] = [pre['slice_info']['age']]
    current_dict['height'] = [pre['slice_info']['height']]
    current_dict['hemi'] = [pre['slice_info']['hemi']]
    current_dict['baseline_percent'] = [pre['slice_info']['baseline_percent'][path]]
    current_dict['baseline_max_idx'] = [pre['slice_info']['baseline_max_idx'][path]]

        # current_df = pd.DataFrame(current_dict)
        # df = df.append(current_df, ignore_index=True)

    return current_dict

def _process_new_data_df(group_df, preprocessed_directory='Preprocessed Data/', search_string='.pkl', functions=[], kwlist=[], rerun=[], keep=[], file_limit=[], ):
    ''' process new data and add to group

    ==Args==
    -group_df : pandas dataframe containing group data
    -directory : relative directory containing unprocessed data
    -search_string :  string to identify unprocessed data, typically ".mat"
    -variables : variables to be derived from unprocessed data (e.g. slopes or slopes_norm)

    ==Out==
    -group_df : pandas dataframe containing group data

    ==Update==
    
    ==Comments==
    '''

    # get list of all preprocessed data files
    data_files = glob.glob(preprocessed_directory+'*'+search_string+'*')
    # remove directory and file extension
    data_filenames = [file.split('\\')[-1].split('.')[0] for file in data_files]

    # get list of processed data files (already added to group dataframe)
    if group_df.empty:
        processed_data_filenames = []
    else:
        processed_data_filenames = group_df.name

    # get list of new data files (haven't been added to dataframe)
    new_data_filenames = list(set(data_filenames)-set(processed_data_filenames))
    new_data_files = [data_files[file_i] for file_i, file in enumerate(data_filenames) if file in new_data_filenames]

    # if any functions specified in rerun, set iterfiles to all data files
    if rerun:
        iterfiles = data_files
    # if there are no functions to rerun, only iterate through new files
    else:
        iterfiles = new_data_files

    print 'total data files:', len(data_files) 
    print 'new data files:', len(new_data_files)
    if rerun:
        print 'functions to rerun over all slices:', rerun

    # iterate over new files and update group_data structure
    #`````````````````````````````````````````````````````````````````
    print 'updating group data structure'

    # temporary df for storing rerun function outputs
    df_update = pd.DataFrame()

    # suffixes for merging df_update with group_df
    old_suffix = '_old'
    rerun_suffix = '_rerun'

    # iterate over data files
    for file_i, file in enumerate(iterfiles):

        if file_limit and file_i> file_limit:
            continue
        else:


            # get file name and date
            name = file.split('\\')[-1].split('.')[0]
            date = int(name[:8])

            # # limit number of files processed at a time
            # if file_limit and file_i>file_limit:
            #     continue
            # else:

            # load data file
            with open(file, 'rb') as pkl_file:

                pre = pickle.load(pkl_file)

                print name

            # for each function create df and add to list (either new or rerun)
            new_dfs = []
            rerun_dfs = []
            # iterate over functions
            for func_i, func in enumerate(functions):
                print func


                # if new data file add to new list
                if file in new_data_files:

                    print 'new file found:', name
                    print kwlist[func_i]
                    df_temp = func(pre, df=[], keep=[], **kwlist[func_i])
                    # # print 'test'
                    new_dfs.append(df_temp)

                # if rerun add to old list
                elif func in rerun:
                    # print name
                    # get corresponding row from group_df
                    df_locs = group_df[group_df.name==name].index.values
                    # get df from function
                    df_temp = func(pre, keep=keep, df=group_df.loc[df_locs], **kwlist[func_i])
                    # print 'rerun keys:',(df_temp.keys())
                    # print func
                    rerun_dfs.append(df_temp)

            # if there are new df's, add to bottom of group_df
            if new_dfs:

                # for temp in new_dfs:
                #     print sorted(temp.keys())
                # merge df's across the different functions, drop columns that are common between functions
                # new_df = reduce(lambda left,right: pd.merge(left, right, on=right.columns.intersection(left.columns), how='inner'), new_dfs)
                # print new_dfs
                new_df = reduce(
                    lambda left,right: pd.merge(
                        left, 
                        right[np.append(right.columns.difference(left.columns).values, ['name', 'path'])], 
                        on=['name', 'path'], 
                        how='inner', ), 
                    new_dfs
                    )



                # print new_df
                # new_df = reduce(lambda left,right: pd.merge(left, right[right.columns.difference(left.columns)], on=['name', 'path'], how='inner' , ), new_dfs)# validate='one_to_one'
                # append to bottom of group_df
                # group_df = group_df.append(new_df, ignore_index=True)
                df_update = df_update.append(new_df, ignore_index=True)

                # print df_update.induction_block_0

            if rerun_dfs:
                # for temp in rerun_dfs:
                #     print temp['name'], temp['path']

                # print np.append(rerun_dfs[0].columns.difference(rerun_dfs[1].columns).values, ['name','path'])

                # rerun_df = reduce(lambda left,right: pd.merge(left, right, on=right.columns.intersection(left.columns), how='inner'), rerun_dfs)

                # print rerun_df
                # merge df's across the different functions
                rerun_df = reduce(
                    lambda left,right: pd.merge(
                        left, 
                        right[np.append(right.columns.difference(left.columns).values, ['name', 'path'])], 
                        on=['name', 'path'], 
                        how='inner', ), 
                    rerun_dfs
                    ) # validate='one_to_one'

                df_update = df_update.append(rerun_df, ignore_index=True)

    # update corresponding rows in group_df
    # for key in df_update:
    #     if 'hilbert' in key:
    #         pdb.set_trace()
    #         print df_update[key]

    if group_df.empty:
        group_df=df_update
    elif not df_update.empty:
        group_df = group_df.merge(df_update, on=['name', 'path'], how='outer', suffixes=[old_suffix, rerun_suffix])



    # print group_df[group_df.date_rerun==20181127].induction_block_0_old

    columns = group_df.columns.values

    drop_columns = []
    for column in columns:
        # if 'data_ind_hilbert_smooth'in column:
            # print column, group_df[column]
        if old_suffix in column:
            column_original = column.split(old_suffix)[0]
            column_new = column_original+rerun_suffix
            # print column_new
            # if not a keep column, use new to update old, then drop new
            if not any([temp for temp in keep if temp in column]):
                # FIXME, use overwrite kwarg for update
                group_df[column].update(group_df[column_new])
                drop_columns.append(column_new)
            # if a keep column, use old to update new, drop old
            else:
                group_df[column_new].update(group_df[column])
                drop_columns.append(column)

        # if 'data_ind_hilbert_smooth'in column:
        #     print column, group_df[column]




    # drop_columns = [column for column in columns if rerun_suffix in column]

    # build list of columns to drop
    # drop_columns=[]
    # for column in columns:
    #     if len(keep)==0:
    #         if old_suffix in column:
    #             drop_columns.append(column)
    #     else:
    #         if old_suffix in column and not any([temp for temp in keep if temp in column]):
    #             drop_columns.append(column)

    #         elif rerun_suffix in column and any([temp for temp in keep if temp in column]):
    #             drop_columns.append(column)

    # print drop_columns
    group_df.drop(drop_columns, axis=1, inplace=True)

    # remove suffixes
    rename_dict={}
    for column in group_df.columns.values:

        if old_suffix in column:
            rename_dict[column] = column.split(old_suffix)[0]
        elif rerun_suffix in column:
            rename_dict[column] = column.split(rerun_suffix)[0]

    group_df.rename(columns=rename_dict, inplace=True)
    # pdb.set_trace()
    # print group_df.data_induction_data_filt_iir_band_300_1000_hilbert_sortby_burst_apical[:10]
    # print group_df.data_induction_data_filt_iir_band_300_1000_hilbert_sortby_burst_apical_old[:10]

    return group_df


        # iterate over columns in group_df
        # if column in keep, use old column values and drop suffix
        # otherwise keep new columnn values and drop suffix


        # merge all function dfs on filename and path
        # update appropriate rows and columns of group_df with function dfs, use merge? find indices and use update/combine?

def _postdict2df(postdict, pre):
    ''' convert dictionary of postprocessing variables into dataframe to be merged with group data
    '''
    df = pd.DataFrame()
    for path in pre['path_blocks']:

        current_dict = _build_conditions_dict(pre, path)

        for key in postdict:
            # print key, postdict[key][path]
            if type(postdict[key][path])==float and postdict[key][path] == np.nan:
                current_dict[key]=postdict[key][path]
            else:
                current_dict[key]=[postdict[key][path]]

        for key in current_dict:
            if type(current_dict[key])==list and len(current_dict[key])>1:
                current_dict[key] = [current_dict[key]]

        # convert to dataframe    
        current_df = pd.DataFrame(current_dict)

        # add to group data
        if df.empty:
            df=current_df
        else:
            df = df.append(current_df, ignore_index=True)

    return df

def _cca(df, conditions, constraints_spec, constraints_all, x_variables, y_variables, ):
    '''
    '''
    df_sorted = _sortdf(df=df, conditions=conditions, constraints_spec=constraints_spec, constraints_all=constraints_all)
    cca = {}
    cca_result={}
    for group_key, group_df in df_sorted.iteritems():
        # print any(pd.isna(group_df[x_variables]).any().values)
        if not any(pd.isna(group_df[x_variables]).any().values) and not any(pd.isna(group_df[y_variables]).any().values) :
            cca[group_key] = CCA()

            x = group_df[x_variables].values
            x_norm = x - np.mean(x, axis=0)
            x_norm = x_norm/np.std(x_norm)
            y = group_df[y_variables].values
            y_norm = y - np.mean(y, axis=0)
            y_norm = y_norm/np.std(y_norm)
            polarity = group_df.field_polarity_0_slopes.values
            x_cca, y_cca= cca[group_key].fit(x_norm, y_norm).transform(x_norm, y_norm)
            x_weights = cca[group_key].x_weights_
            y_weights = cca[group_key].y_weights_


            cca_result[group_key]={
            'x':x,
            'y':y,
            'x_norm':x_norm,
            'y_norm':y_norm,
            'x_variables':x_variables,
            'y_variables':y_variables,
            'cca_obj':cca[group_key],
            'x_cca':x_cca,
            'y_cca':y_cca,
            'x_weights':x_weights,
            'y_weights':y_weights,
            'polarity':polarity
            }
    return cca_result

def _plot_scalebar(axes, yscale, xscale, origin, width):
    '''
    '''
    # plot y
    axes.plot([origin[0], origin[0]], [origin[1], origin[1]+yscale], linewidth=width, color='k')
    axes.plot([origin[0], origin[0]+xscale], [origin[1], origin[1]], linewidth=width, color='k')

    return axes

def _plot_vtrace(df_sorted, conditions, figures, variable, colors, markers, titles, mean=True, all_paths=True, **kwargs):
    '''FIXME add docs
    '''
    # FIXME add kwargs to alter figure details
    # create figure groupings (all conditions that will go on the same figure)
    fig = {}
    ax={}
    ylim=[]
    # iterate over figures
    for figure_key in figures:
        if type(figure_key)==str:
            figure_key=(figure_key,)
        # create figure objects
        fig[figure_key], ax[figure_key] = plt.subplots()
        # plot object for each trace
        trace = OrderedDict()
        trace_final = OrderedDict()
        # n for each condition
        n_trial=OrderedDict()
        # pvalue for ttest comparing to corresponding control condition
        pval=OrderedDict()
        # label for plot legend
        label=OrderedDict()
        # iterate over each condition
        for group_key, df in df_sorted.iteritems():
            # print figure_key
            # print group_key
            # print len(set(group_key) - set(figure_key))
            # if type(group_key)==tuple:
            #     group_key = list(group_key)
            # elif type(group_key)==str:
            #     group_key = [group_key]
            # print len(set(group_key) - set(figure_key))
            # if plotting all pathways on the same figure
            if all_paths and not all(temp in group_key for temp in figure_key) and 'figures_any' not in kwargs:
                # print 'huuuuuh'
                # if condition not met continue
                continue
            # if only plotting specific pathways
            elif not all_paths and not group_key[1]==figure_key[0] and group_key[2]==figure_key[1] and 'figures_any' not in kwargs:
                # if condition not met, continue
                # print 'huhhhhhhh'
                continue

            if 'figures_any' in kwargs and len(set(group_key) - set(figure_key)) != 0:

                # print figure_key, group_key,'huh'


            # any(set(group_key)&set(figure_key)):
                continue
            # if conditions are satisfied, plot figure
            else:
                # get color
                color = colors[group_key[0]]
                if 'colors_post' in kwargs:
                    print 'color_found'
                    color_post = kwargs['colors_post'][group_key[0]]
                else:
                    color_post=color    

                color_final = tuple([0.4*val for val in color])
                # list of group key values and list of marker keys,use the first match that you find
                if 'all' in markers.keys():
                    marker=markers['all']
                elif len(list(set(markers.keys()) & set(group_key)))>0:
                    marker_key = list(set(markers.keys()) & set(group_key))[0]
                    marker = markers[marker_key]
                else:
                    marker='.'


                # print set(group_key)
                # print set(markers.keys())
                # print list(set(markers.keys()) & set(group_key))
                # marker_key = list(set(markers.keys()) & set(group_key))[0]
                # # get marker
                # marker = markers[marker_key]
                # convert data to 2d array (slices x time)
                data_array = _2array(df[variable])
                print data_array.shape
                # check that data is 2d (if not, it is due to a lack of data for this condition)
                if len(data_array.shape)>1:
                    # mean across slices
                    data_baseline = np.mean(data_array[:,:,:20], axis=2).squeeze()
                    data_final = np.mean(data_array[:,:,10:], axis=2).squeeze()

                    data_baseline_mean = np.mean(data_baseline, axis=0)
                    data_final_mean = np.mean(data_final, axis=0)

                    #std across slices
                    data_baseline_std = np.std(data_baseline, axis=0)
                    data_final_std = np.std(data_final, axis=0)
                    # sem across slices
                    data_baseline_sem = stats.sem(data_baseline, axis=0)
                    data_final_sem = stats.sem(data_final, axis=0)
                    # time vector
                    t = np.arange(len(data_baseline_mean))*0.1
                    print t.shape, data_baseline_mean.shape
                    # if plotting mean w errorbars
                    # print figure_key
                    # print group_key
                    if mean:
                        if 'shade_error' in kwargs and kwargs['shade_error']:
                            plt.fill_between(t, data_baseline_mean-data_baseline_sem, data_baseline_mean+data_baseline_sem, color=color, alpha=1)
                            
                            # trace[group_key] = ax[figure_key].plot(t, data_baseline_mean, color=color, linewidth=4)
                            plt.fill_between(t, data_final_mean-data_final_sem, data_final_mean+data_final_sem, color=color_post, alpha=1)
                            # trace_final[group_key] = ax[figure_key].plot(t, data_final_mean, color=color_final, linewidth=4)
                            _plot_scalebar(axes=ax[figure_key], xscale=2, yscale=.001, origin=[9,-.0035], width=4)
                            plt.axis('off')
                        else:
                            trace[group_key] = ax[figure_key].errorbar(t, data_baseline_mean, yerr=data_baseline_sem, color=color, linewidth=1)
                            trace_final[group_key] = ax[figure_key].errorbar(t, data_final_mean, yerr=data_final_sem, color=color, linewidth=2)




        ylim = [-.004, .0005]
        # plt.title(titles[figure_key]+' '+variable)
        # plt.legend(trace.values(), label.values())
        if ylim:
            plt.ylim(ylim)
        plt.xlim([0,12])
        plt.show(block=False)

    return fig, ax

def _plot_timeseries(df_sorted, conditions, figures, variable, colors, markers, titles, mean=True, all_paths=True, **kwargs):
    '''FIXME add docs
    '''
    # FIXME add kwargs to alter figure details
    # create figure groupings (all conditions that will go on the same figure)
    fig = {}
    ax={}
    ylim={}
    xlim={}
    # iterate over figures
    peak_value_all=0
    for figure_key in figures:
        if type(figure_key)==str:
            figure_key=(figure_key,)
        # create figure objects
        fig[figure_key], ax[figure_key] = plt.subplots(figsize=(8,6), dpi=100 )
        # plot object for each trace
        trace = OrderedDict()
        # n for each condition
        n_trial=OrderedDict()
        # pvalue for ttest comparing to corresponding control condition
        pval=OrderedDict()
        # label for plot legend
        label=OrderedDict()
        # iterate over each condition
        peak_value_current = 0
        for group_key, df in df_sorted.iteritems():
            # print figure_key
            # print group_key
            # print len(set(group_key) - set(figure_key))
            # if type(group_key)==tuple:
            #     group_key = list(group_key)
            # elif type(group_key)==str:
            #     group_key = [group_key]
            # print len(set(group_key) - set(figure_key))
            # if plotting all pathways on the same figure
            if all_paths and not all(temp in group_key for temp in figure_key) and 'figures_any' not in kwargs:
                # print 'huuuuuh'
                # if condition not met continue
                continue
            # if only plotting specific pathways
            elif not all_paths and not group_key[1]==figure_key[0] and group_key[2]==figure_key[1] and 'figures_any' not in kwargs:
                # if condition not met, continue
                # print 'huhhhhhhh'
                continue

            if 'figures_any' in kwargs and len(set(group_key) - set(figure_key)) != 0:

                # print figure_key, group_key,'huh'


            # any(set(group_key)&set(figure_key)):
                continue
            # if conditions are satisfied, plot figure
            else:
                # get color
                color = colors[group_key[0]]

                marker_match = [key for key in group_key if key in markers.keys()]
                # list of group key values and list of marker keys,use the first match that you find
                if 'all' in markers.keys():
                    marker=markers['all']

                elif len(marker_match)>0:
                    marker = markers[marker_match[0]]


                if 'line_markers' in kwargs:
                    line_markers=kwargs['line_markers']
                else:
                    line_markers={}
                line_marker_match = [key for key in group_key if key in line_markers.keys()]
                # list of group key values and list of marker keys,use the first match that you find
                if 'all' in line_markers.keys():
                    line_marker=line_markers['all']

                elif len(line_marker_match)>0:
                    line_marker = line_markers[line_marker_match[0]]
                # elif len(list(set(markers.keys()) & set(group_key)))==1:
                #     marker_key = list(set(group_key) & set(markers.keys()))[0]
                #     marker = markers[marker_key]

                # elif len(list(set(markers.keys()) & set(group_key)))>1:

                else:
                    line_marker='solid'


                # print set(group_key)
                # print set(markers.keys())
                # print list(set(markers.keys()) & set(group_key))
                # marker_key = list(set(markers.keys()) & set(group_key))[0]
                # # get marker
                # marker = markers[marker_key]
                # convert data to 2d array (slices x time)
                data_array = _2array(df[variable], remove_nans=True, remove_nans_axis=1)
                print data_array.shape
                # check that data is 2d (if not, it is due to a lack of data for this condition)
                if len(data_array.shape)>1:
                    # mean across slices
                    data_mean = np.mean(data_array, axis=0)
                    #std across slices
                    data_std = np.std(data_array, axis=0)
                    # sem across slices
                    data_sem = stats.sem(data_array, axis=0)
                    # time vector
                    t = np.arange(len(data_mean))
                    print t.shape, data_mean.shape
                    # if plotting mean w errorbars
                    # print figure_key
                    # print group_key
                    if mean:
                        print np.max(data_mean+data_sem)
                        if np.max(data_mean+data_sem)>peak_value_current:
                            peak_value_current= np.max(data_mean+data_sem)
                        if peak_value_current>peak_value_all:
                            peak_value_all =  peak_value_current

                        print peak_value_all
                        if 'shade_error' in kwargs and kwargs['shade_error']:
                            plt.fill_between(t, data_mean-data_sem, data_mean+data_sem, color=color, alpha=1)

                            # trace[group_key] = ax[figure_key].plot(t, data_mean, color=color, marker=marker,  markersize=10, linewidth=0, markerfacecolor=color)

                            trace[group_key] = ax[figure_key].plot(t, data_mean, color=color, marker='None',  markersize=10, linewidth=6, markerfacecolor=color, linestyle=line_marker, )
                            # plt.fill_between(t, data_final_mean-data_final_sem, data_final_mean+data_final_sem, color=color_post, alpha=1)
                        else:
                            trace[group_key] = ax[figure_key].errorbar(t, data_mean, yerr=data_sem, color=color, marker=marker,  markersize=15, elinewidth=2, linewidth=3, markerfacecolor=color)
                    # otherwise plot all traces
                    else:
                        trace[group_key] = ax[figure_key].plot(t, data_array.T, color=color, marker=marker,  markersize=10)

                    # get number of slices
                    n_trial[group_key] = data_array.shape[0]

                    # build control group
                    control_group = []
                    for cond_i, cond in enumerate(conditions):
                        if cond == 'field_polarity_0':
                            control_group.append('control')
                        elif cond=='field_mag_0':
                            control_group.append('0')
                        else:
                            control_group.append(group_key[cond_i])
                    control_group= tuple(control_group)

                    if 'slopes_norm' in variable:
                        # final ltp (length n), taken as the average of the last 10 normalized slopes for each slice
                        ltp_final = df_sorted[group_key].ltp_final
                        # corresponding control values
                        ltp_final_control = df_sorted[control_group].ltp_final
                        # ttest compared to control
                        pval[group_key] = stats.ttest_ind(ltp_final, ltp_final_control)[1]
                        ltp_final_mean=np.mean(ltp_final)
                        ltp_final_std=np.std(ltp_final)
                        ltp_final_sem = stats.sem(ltp_final)
                        # label for figure legend including stimulation conditions, n slices, p value compared to the corresponding control
                        label_keys=[]
                        for key in group_key:
                            # print key
                            label_keys.append(key)
                        label_keys = label_keys + ['n={}'.format(n_trial[group_key]), 'p = {0:.3f}'.format(pval[group_key]), 'mean={0:.3f}'.format(ltp_final_mean),'sem={0:.3f}'.format(ltp_final_sem)]
                        # print label_keys
                        label[group_key] = ', '.join(label_keys)
                        # label[group_key] = ', '.join([group_key[0], group_key[1], 'n={}'.format(n_trial[group_key]), 'p = {0:.3f}'.format(pval[group_key])])
            # if not 'common_ylim' in kwargs or not kwargs['common_ylim']:
            # get x and y limits based data
            xlim[figure_key] = ax[figure_key].get_xlim()
            ylim[figure_key] = ax[figure_key].get_ylim()
            # ylim[figure_key] = [0.9, peak_value_current]
            # plt.title(titles[figure_key]+' '+variable)
            if 'show_stats' in kwargs and kwargs['show_stats']:
                plt.legend(trace.values(), label.values())

            

    # if 'common_ylim' in kwargs and kwargs['common_ylim']:
    #     ylim = [0.9, peak_value_all]
    # else:
    #     ylim = [0.9, peak_value_current]
    # print ylim
    # xlim = [0,80]
    for figure_key in ax:
        if 'common_ylim' in kwargs and kwargs['common_ylim'] and 'ylim' not in kwargs:
            ylim[figure_key] = [0.9, peak_value_all]
        elif 'ylim' in kwargs:
            ylim[figure_key] = kwargs['ylim']

        if 'xlim' in kwargs:
            xlim[figure_key] = kwargs['xlim']
        
        # format figure
        ax[figure_key].spines['right'].set_visible(False)
        ax[figure_key].spines['top'].set_visible(False)
        ax[figure_key].spines['left'].set_linewidth(5)
        ax[figure_key].spines['bottom'].set_linewidth(5)
        ax[figure_key].xaxis.set_ticks_position('bottom')
        ax[figure_key].yaxis.set_ticks_position('left')
        ax[figure_key].set_xlabel('Time (min)', fontsize=25, fontweight='heavy')
        ax[figure_key].set_ylabel('Normalized fEPSP slope', fontsize=25, fontweight='heavy')
        xticks = np.arange(0,81, 20)
        yticks = np.round(np.arange(1.,peak_value_all, 0.2), decimals=1)
        ax[figure_key].set_xticks(xticks)
        ax[figure_key].set_xticklabels(xticks, fontsize=20, fontweight='heavy')
        ax[figure_key].set_yticks(yticks)
        ax[figure_key].set_yticklabels(yticks, fontsize=20, fontweight='heavy')
        
        ax[figure_key].set_ylim(ylim[figure_key])
        ax[figure_key].set_xlim(xlim[figure_key])
        plt.figure(fig[figure_key].number)
        plt.tight_layout()
        
        # if ylim:
        #     plt.ylim(ylim)
    plt.show(block=False)
    

    return fig, ax
    
def _plot_bar(df_sorted, figure_params, variable, group_space=1, bar_width=1, bar_spacing=1):
    '''
    '''
    print 'look here:',figure_params.keys()
    fig={}
    ax={}
    n_subgroups={}
    n_traces={}
    xlim={}
    ylim={}
    for figure_key, figure_subgroups in figure_params.iteritems():
        if figure_key!='params':
            fig[figure_key], ax[figure_key] = plt.subplots()
            n_subgroups[figure_key] = len(figure_subgroups.keys()) 
            n_traces[figure_key]={}
            locations = []
            heights=[]
            colors=[]
            fig_args=[]
            xticks=[]
            xticklabels=[]
            cnt=bar_spacing
            for subgroup_key, traces in figure_subgroups.iteritems():
                if subgroup_key!='params':
                    n_traces[figure_key][subgroup_key]=len(traces.keys())
                    cnt+=group_space
                    for trace_key, params in traces.iteritems():
                        if trace_key!='params':
                            trace_args={}
                            # cnt+=bar_spacing
                            if 'location' in params:
                                locations.append(cnt+params['location'])
                            else:
                                locations.append(locations[-1]+1)
                            xticks.append(cnt)
                            xticklabels.append(params['label'])

                            # get data and stats
                            if trace_key in df_sorted:
                                trace_series = df_sorted[trace_key][variable]
                                data_array = (_2array(trace_series, remove_nans=True, remove_nans_axis=1)-1)*100.
                                print type(data_array)
                                print type(data_array)==np.float64
                                if type(data_array)==np.ndarray:
                                    # get stats
                                    # mean across slices
                                    data_mean = np.mean(data_array, axis=0)
                                    #std across slices
                                    data_std = np.std(data_array, axis=0)
                                    # sem across slices
                                    data_sem = stats.sem(data_array, axis=0)

                                    heights.append(data_mean)

                                    trace_args['color'] = params['color']

                                    fig_args.append(trace_args)

                                    colors.append(params['color'])

                                    plt.errorbar(locations[-1], data_mean, data_sem, color=(.5,.5,.5))

            barcontainer = ax[figure_key].bar(locations, heights, width=bar_width, tick_label=xticklabels)
            xlim[figure_key] = ax[figure_key].get_xlim()
            ylim[figure_key] = ax[figure_key].get_ylim()
            # ax[figure_key].set_xticks(xticks, xticklabels,)
            # barcontainer = ax[figure_key].violinplot(locations, heights, width=bar_width, tick_label=xticklabels)
            print 'rotations:', figure_key, figure_params[figure_key]['params']
            ax[figure_key].set_xticklabels(xticklabels, fontsize
                =20, fontweight='heavy', rotation=figure_params[figure_key]['params']['rotation'])
            for bar_i, bar in enumerate(barcontainer):
                bar.set_color(colors[bar_i])

    # get ylim and xlim across all figures
    xlims=[]
    ylims=[]
    for figure_key in ylim:
        xlims.append(xlim[figure_key])
        ylims.append(ylim[figure_key])
    xlim_all = [min([temp[0] for temp in xlims]), max([temp[1] for temp in xlims])]
    ylim_all = [min([temp[0] for temp in ylims]), max([temp[1] for temp in ylims])]

    xlim={}
    ylim={}
    print figure_params.keys()
    for figure_key, axes in ax.iteritems():
        if 'ylim_all' in figure_params['params'] and figure_params['params']['ylim_all']:
            print 'setting ylim'
            ax[figure_key].set_ylim(ylim_all)
            ax[figure_key].set_xlim(xlim_all)

        
        # format figure
        ax[figure_key].spines['right'].set_visible(False)
        ax[figure_key].spines['top'].set_visible(False)
        ax[figure_key].spines['left'].set_linewidth(5)
        ax[figure_key].spines['bottom'].set_linewidth(5)
        ax[figure_key].xaxis.set_ticks_position('bottom')
        ax[figure_key].yaxis.set_ticks_position('left')
        # ax[figure_key].set_xlabel('Time (min)', fontsize=25, fontweight='heavy')
        ax[figure_key].set_ylabel('% LTP', fontsize=25, fontweight='heavy')
        xticks = np.arange(0,81, 20)
        ytickmax = ax[figure_key].get_ylim()[1]
        ytickmin = ax[figure_key].get_ylim()[0]
        yticks = np.round(np.arange(0,ytickmax, 10), decimals=0).astype(int)
        # ax[figure_key].set_xticks(xticks)
        # ax[figure_key].set_xticklabels(xticks, fontsize=20, fontweight='heavy')
        ax[figure_key].set_yticks(yticks)
        ax[figure_key].set_yticklabels(yticks, fontsize=20, fontweight='heavy')
        
        # ax[figure_key].set_ylim(ylim[figure_key])
        # ax[figure_key].set_xlim(xlim)
        plt.figure(fig[figure_key].number)
        plt.tight_layout()

    plt.show(block=False)

    return fig, ax

def _plot_var2var_correlation(df_sorted, figures, variables, colors, markers, titles):
    ''' FIXME adds docs
    '''
    # create figure groupings (all conditions that will go on the same figure)
    fig = {}
    ax={}
    # iterate over figures
    for figure_key in figures:
        # create figure objects
        fig[figure_key], ax[figure_key] = plt.subplots()
        # plot object for each trace
        trace = OrderedDict()
        # n for each condition
        n_trial=OrderedDict()
        # pvalue for ttest comparing to corresponding control condition
        pval=OrderedDict()
        # label for plot legend
        label=OrderedDict()
        # iterate over each condition
        var0_all = []
        var1_all = []
        for group_key, df in df_sorted.iteritems():
            if group_key[1]==figure_key[0] and group_key[2]==figure_key[1]:
            # check that the current conditions belong on the current figure
            # if all(temp in group_key for temp in figure_key):
                color = colors[group_key[0]]
                # print group_key[1]
                marker = markers[group_key[1]]
                var0 = np.array(df[variables[0]])
                var1 = np.array(df[variables[1]])
                var0_all.append(var0)
                var1_all.append(var1)
                trace[group_key] = ax[figure_key].plot(var0, var1, color=color, marker=marker,  markersize=10, linestyle='None')
                # plt.plot(var, ltp, color=color, marker=marker)

        var0_all = np.append(var0_all[0], var0_all[1:])
        var1_all = np.append(var1_all[0], var1_all[1:])
        varmat = np.vstack([var0_all, var1_all])
        # print variables, varmat.shape
        # covariance = np.cov(varmat)
        pearsonsr, pearsonsp = stats.pearsonr(var0_all, var1_all)
        print variables[0], pearsonsr, pearsonsp

        plt.title(titles[figure_key])
        plt.xlabel(variables[0])
        plt.ylabel(variables[1])

        # plot correlation coefficient with label

    plt.show(block=False)

def _plot_trace_brian(df_sorted, figure_params, variable, **kwargs):
    '''FIXME add docs
    '''
    # FIXME add kwargs to alter figure details
    # create figure groupings (all conditions that will go on the same figure)
    fig={}
    ax={}
    n_subgroups={}
    n_traces={}
    xlim={}
    ylim={}
    # iterate over figures
    for figure_key, figure_subgroups in figure_params.iteritems():
        # params kw specifies figure level kwargs
        if figure_key!='params':
            # create figure, passing params as **kwargs
            fig[figure_key], ax[figure_key] = plt.subplots()
            # number of subgroup in figure
            n_subgroups[figure_key] = len(figure_subgroups.keys()) 
            n_traces[figure_key]={}
            # iterate over subgroups of traces
            for subgroup_key, traces in figure_subgroups.iteritems():

                if subgroup_key!='params':
                    # FIXME distribute subgroup padrameters to each trace in the subgroup, with priority to trace parameters
                    n_traces[figure_key][subgroup_key]=len(traces.keys())
                    # get trial id's
                    trace_key_temp = sorted(traces.keys())
                    trace_key_temp = [temp for temp in trace_key_temp if temp!='params'][0]
                    print traces.keys()
                    if'plot_individual_trace' in traces['params'] and len(traces['params']['plot_individual_trace'])>0:
                        trial_ids=[]
                        for iloc in traces['params']['plot_individual_trace']:
                            trial_ids.append(df_sorted[trace_key_temp].trial_id.iloc[iloc])
                    

                    # iterate over individual traces
                    for trace_key, params in traces.iteritems():
                        if trace_key!='params':
                            # get data and stats
                            print df_sorted.keys()
                            trace_series = df_sorted[trace_key][variable]
                            print subgroup_key, trace_key, trace_series
                            data_array = _2array(trace_series, remove_nans=True, remove_nans_axis=1)#*1000
                            # get stats
                            # mean across slices
                            data_mean = np.mean(data_array, axis=0)
                            #std across slices
                            data_std = np.std(data_array, axis=0)
                            # sem across slices
                            data_sem = stats.sem(data_array, axis=0)
                            # time vector
                            t = np.arange(len(data_mean))#/10.
                            if 'plot_mean' in traces['params'] and traces['params']['plot_mean']:
                                # line trace with shaded error
                                #______________________________
                                if 'shade_error' in params and params['shade_error']:

                                    ax[figure_key].plot(t, data_mean, color=params['color'], linewidth=params['linewidth'])
                                    plt.fill_between(t, data_mean-data_sem, data_mean+data_sem, **params['e_params'])

                            elif'plot_individual_trace' in traces['params'] and len(traces['params']['plot_individual_trace'])>0:
                                for trial_id in trial_ids:
                                    trace_data = df_sorted[trace_key][df_sorted[trace_key].trial_id==trial_id][variable]
                                    trace_data = _2array(trace_data)#*1000
                                    t = np.arange(len(trace_data))#/10.
                                    ax[figure_key].plot(t, trace_data, color=params['color'], linewidth=params['linewidth'])

            # get x and y limits based data
            xlim[figure_key] = ax[figure_key].get_xlim()
            ylim[figure_key] = ax[figure_key].get_ylim()
    
    # get ylim and xlim across all figures
    xlims=[]
    ylims=[]
    for figure_key in ylim:
        xlims.append(xlim[figure_key])
        ylims.append(ylim[figure_key])
    xlim_all = [min([temp[0] for temp in xlims]), max([temp[1] for temp in xlims])]
    ylim_all = [min([temp[0] for temp in ylims]), max([temp[1] for temp in ylims])]

    # set common ylim across all figures
    for figure_key, axes in ax.iteritems():
        if 'ylim_all' in figure_params['params'] and figure_params['params']['ylim_all']:
            print 'setting ylim to be the same across all figures'
            if 'ylim' in figure_params['params']:
                ax[figure_key].set_ylim(figure_params['params']['ylim'])
            else:
                ax[figure_key].set_ylim(ylim_all)
        if 'xlim_all' in figure_params['params'] and figure_params['params']['xlim_all']:
            if 'xlim' in figure_params['params']:
                ax[figure_key].set_xlim(figure_params['params']['xlim'])
            else:
                ax[figure_key].set_xlim(xlim_all)

        # format figure
        ax[figure_key].spines['right'].set_visible(False)
        ax[figure_key].spines['top'].set_visible(False)
        ax[figure_key].spines['left'].set_linewidth(5)
        ax[figure_key].spines['bottom'].set_linewidth(5)
        ax[figure_key].xaxis.set_ticks_position('bottom')
        ax[figure_key].yaxis.set_ticks_position('left')
        ax[figure_key].set_xlabel(figure_params['params']['xlabel'], fontsize=25, fontweight='heavy')
        ax[figure_key].set_ylabel(figure_params['params']['ylabel'], fontsize=25, fontweight='heavy')
        for temp in ax[figure_key].get_xticklabels():
            temp.set_fontweight('heavy')
            temp.set_fontsize(15)
        for temp in ax[figure_key].get_yticklabels():
            temp.set_fontweight('heavy')
            temp.set_fontsize(15)
        # print xticks, xticklabels
        # ax[figure_key].set_xticks(xticks)
        # ax[figure_key].set_xticklabels(xticklabels, fontsize=20, fontweight='heavy')
        # ax[figure_key].tick_params(axis='both', labelsize=20, )
        # xticks = figure_params['xticks'] #np.arange(0,81, 20)
        # ytickmax = ax[figure_key].get_ylim()[1]
        # ytickmin = ax[figure_key].get_ylim()[0]
        # yticks = figure_params['yticks']#np.round(np.arange(0,ytickmax, 10), decimals=0).astype(int)
        # ax[figure_key].set_xticks(xticks)
        # ax[figure_key].set_xticklabels(xticks, fontsize=20, fontweight='heavy')
        # ax[figure_key].set_yticks(yticks)
        # ax[figure_key].set_yticklabels(yticks, fontsize=20, fontweight='heavy')
        
        # ax[figure_key].set_ylim(ylim[figure_key])
        # ax[figure_key].set_xlim(xlim)
        plt.figure(fig[figure_key].number)
        plt.tight_layout()
        plt.show(block=False)

    return fig, axes

def _plot_induction_bursts(df_sorted, figure_params, variable, **kwargs):
    '''FIXME add docs
    '''
    # FIXME add kwargs to alter figure details
    # create figure groupings (all conditions that will go on the same figure)
    fig={}
    ax={}
    n_subgroups={}
    n_traces={}
    xlim={}
    ylim={}
    average_bursts = figure_params['params']['average_across_bursts']
    burst_number = figure_params['params']['burst_number']
    induction_number = figure_params['params']['induction_number']
    yscale=figure_params['params']['yscale']
    xscale=figure_params['params']['xscale']
    # iterate over figures
    for figure_key, figure_subgroups in figure_params.iteritems():
        # params kw specifies figure level kwargs
        if figure_key!='params':
            # create figure, passing params as **kwargs
            fig[figure_key], ax[figure_key] = plt.subplots()
            # number of subgroup in figure
            n_subgroups[figure_key] = len(figure_subgroups.keys()) 
            n_traces[figure_key]={}
            # iterate over subgroups of traces
            for subgroup_key, traces in figure_subgroups.iteritems():

                if subgroup_key!='params':
                    # FIXME distribute subgroup padrameters to each trace in the subgroup, with priority to trace parameters
                    n_traces[figure_key][subgroup_key]=len(traces.keys())
                    # get trial id's
                    trace_key_temp = sorted(traces.keys())
                    trace_key_temp = [temp for temp in trace_key_temp if temp!='params'][0]
                    
                    # iterate over individual traces
                    for trace_key, params in traces.iteritems():
                        if trace_key!='params':
                            # get data and stats
                            print df_sorted.keys()
                            trace_series = df_sorted[trace_key][variable]
                            print subgroup_key, trace_key, trace_series
                            # axis0:slices, axis1:time, axis2:bursts
                            data_array = _2array(trace_series, remove_nans=True, remove_nans_axis=1, list_index=induction_number)#*1000
                            if data_array.shape[0]!=0:
                                print data_array.shape
                                # if there is more than 1 slice for the given conditions, 
                                # data array is 3d, slices x samples x bursts, and must be reduced to 2d
                                if len(data_array.shape)>2:
                                    # average over bursts
                                    if average_bursts:
                                        data_array = np.mean(data_array, axis=2).squeeze()
                                    # show a specific burst
                                    else:
                                        data_array= data_array[:,:,burst_number].squeeze()

                                data_array=data_array*yscale

                                # get stats
                                # mean across slices
                                data_mean = np.mean(data_array, axis=0)
                                #std across slices
                                data_std = np.std(data_array, axis=0)
                                # sem across slices
                                data_sem = stats.sem(data_array, axis=0)
                                # time vector
                                t = np.arange(len(data_mean))*xscale
                                if 'plot_mean' in traces['params'] and traces['params']['plot_mean']:
                                    # line trace with shaded error
                                    #______________________________
                                    if 'shade_error' in params and params['shade_error']:

                                        ax[figure_key].plot(t, data_mean, color=params['color'], linewidth=params['linewidth'])
                                        plt.fill_between(t, data_mean-data_sem, data_mean+data_sem, **params['e_params'])

            # get x and y limits based data
            xlim[figure_key] = ax[figure_key].get_xlim()
            ylim[figure_key] = ax[figure_key].get_ylim()
            print 'ylim:', ylim
    
    # get ylim and xlim across all figures
    xlims=[]
    ylims=[]
    for figure_key in ylim:
        xlims.append(xlim[figure_key])
        ylims.append(ylim[figure_key])
    xlim_all = [min([temp[0] for temp in xlims]), max([temp[1] for temp in xlims])]
    ylim_all = [min([temp[0] for temp in ylims]), max([temp[1] for temp in ylims])]

    # set common ylim across all figures
    for figure_key, axes in ax.iteritems():
        if 'ylim_all' in figure_params['params'] and figure_params['params']['ylim_all']:
            print 'setting ylim to be the same across all figures'
            if 'ylim' in figure_params['params']:
                ax[figure_key].set_ylim(figure_params['params']['ylim'])
            else:
                ax[figure_key].set_ylim(ylim_all)
        if 'xlim_all' in figure_params['params'] and figure_params['params']['xlim_all']:
            if 'xlim' in figure_params['params']:
                ax[figure_key].set_xlim(figure_params['params']['xlim'])
            else:
                ax[figure_key].set_xlim(xlim_all)

        # format figure
        ax[figure_key].spines['right'].set_visible(False)
        ax[figure_key].spines['top'].set_visible(False)
        ax[figure_key].spines['left'].set_linewidth(5)
        ax[figure_key].spines['bottom'].set_linewidth(5)
        ax[figure_key].xaxis.set_ticks_position('bottom')
        ax[figure_key].yaxis.set_ticks_position('left')
        ax[figure_key].set_xlabel(figure_params['params']['xlabel'], fontsize=25, fontweight='heavy')
        ax[figure_key].set_ylabel(figure_params['params']['ylabel'], fontsize=25, fontweight='heavy')
        for temp in ax[figure_key].get_xticklabels():
            temp.set_fontweight('heavy')
            temp.set_fontsize(15)
        for temp in ax[figure_key].get_yticklabels():
            temp.set_fontweight('heavy')
            temp.set_fontsize(15)
        # print xticks, xticklabels
        # ax[figure_key].set_xticks(xticks)
        # ax[figure_key].set_xticklabels(xticklabels, fontsize=20, fontweight='heavy')
        # ax[figure_key].tick_params(axis='both', labelsize=20, )
        # xticks = figure_params['xticks'] #np.arange(0,81, 20)
        # ytickmax = ax[figure_key].get_ylim()[1]
        # ytickmin = ax[figure_key].get_ylim()[0]
        # yticks = figure_params['yticks']#np.round(np.arange(0,ytickmax, 10), decimals=0).astype(int)
        # ax[figure_key].set_xticks(xticks)
        # ax[figure_key].set_xticklabels(xticks, fontsize=20, fontweight='heavy')
        # ax[figure_key].set_yticks(yticks)
        # ax[figure_key].set_yticklabels(yticks, fontsize=20, fontweight='heavy')
        
        # ax[figure_key].set_ylim(ylim[figure_key])
        # ax[figure_key].set_xlim(xlim)
        plt.figure(fig[figure_key].number)
        plt.tight_layout()
        plt.show(block=False)

    return fig, axes

def _plot_dose_response_all(df_sorted, figure_params, variable, group_space=1, bar_width=1, bar_spacing=1):
    '''
    '''
    print 'look here:',figure_params['params'].keys()
    if 'markersize' in figure_params['params']:
        markersize=figure_params['params']['markersize']
    fig={}
    ax={}
    n_subgroups={}
    n_traces={}
    xlim={}
    ylim={}
    for figure_key, figure_subgroups in figure_params.iteritems():
        if figure_key!='params':
            fig[figure_key], ax[figure_key] = plt.subplots()
            n_subgroups[figure_key] = len(figure_subgroups.keys()) 
            n_traces[figure_key]={}
            locations = []
            heights=[]
            colors=[]
            fig_args=[]
            xticks=[]
            xticklabels=[]
            cnt=bar_spacing
            for subgroup_key, traces in figure_subgroups.iteritems():
                if subgroup_key!='params':
                    n_traces[figure_key][subgroup_key]=len(traces.keys())
                    cnt+=group_space
                    for trace_key, params in traces.iteritems():
                        if trace_key!='params':
                            trace_args={}
                            # cnt+=bar_spacing
                            if 'location' in params:
                                locations.append(params['location'])
                            else:
                                locations.append(locations[-1]+1)
                            xticks.append(cnt)
                            xticklabels.append(params['label'])

                            current_location = params['location']

                            # FIXME handle single data points
                            # get data and standardts
                            if trace_key in df_sorted:
                                trace_series = df_sorted[trace_key][variable]
                                data_array = (_2array(trace_series, remove_nans=True, remove_nans_axis=1)-1)*100.
                                print data_array
                                print current_location
                                if type(data_array)==np.float64:
                                    current_locations = current_location
                                else:
                                    current_locations = current_location*np.ones(len(data_array))

                                plt.plot(current_locations, data_array, '.', color=params['color'], markersize=markersize)
                                # print trace_key, data_array.shape, type(data_array)
                                # if type(data_array)==np.ndarray:
                                #     # get stats
                                #     # mean across slices
                                #     data_mean = np.mean(data_array, axis=0)
                                #     #std across slices
                                #     data_std = np.std(data_array, axis=0)
                                #     # sem across slices
                                #     data_sem = stats.sem(data_array, axis=0)

                                #     heights.append(data_mean)

                                #     trace_args['color'] = params['color']

                                #     fig_args.append(trace_args)

                                #     colors.append(params['color'])

                                #     plt.errorbar(locations[-1], data_mean, data_sem, color=(.5,.5,.5))

            # barcontainer = ax[figure_key].bar(locations, heights, width=bar_width, tick_label=xticklabels)
            # xlim[figure_key] = ax[figure_key].get_xlim()
            # ylim[figure_key] = ax[figure_key].get_ylim()
            # # ax[figure_key].set_xticks(xticks, xticklabels,)
            # # barcontainer = ax[figure_key].violinplot(locations, heights, width=bar_width, tick_label=xticklabels)
            # print 'rotations:', figure_key, figure_params[figure_key]['params']
            # ax[figure_key].set_xticklabels(xticklabels, fontsize
            #     =20, fontweight='heavy', rotation=figure_params[figure_key]['params']['rotation'])
            # for bar_i, bar in enumerate(barcontainer):
            #     bar.set_color(colors[bar_i])
            plt.show(block=False)

    return fig, ax

def _plot_dose_response_mean(df_sorted, figure_params, variable, group_space=1, bar_width=1, bar_spacing=1):
    '''
    '''
    print 'look here:',figure_params['params'].keys()
    if 'markersize' in figure_params['params']:
        markersize=figure_params['params']['markersize']
    fig={}
    ax={}
    n_subgroups={}
    n_traces={}
    xlim={}
    ylim={}
    for figure_key, figure_subgroups in figure_params.iteritems():
        if figure_key!='params':
            fig[figure_key], ax[figure_key] = plt.subplots()
            n_subgroups[figure_key] = len(figure_subgroups.keys()) 
            n_traces[figure_key]={}
            locations = []
            heights=[]
            colors=[]
            fig_args=[]
            xticks=[]
            xticklabels=[]
            cnt=bar_spacing
            for subgroup_key, traces in figure_subgroups.iteritems():
                if subgroup_key!='params':
                    n_traces[figure_key][subgroup_key]=len(traces.keys())
                    cnt+=group_space
                    for trace_key, params in traces.iteritems():
                        if trace_key!='params':
                            trace_args={}
                            # cnt+=bar_spacing
                            if 'location' in params:
                                locations.append(params['location'])
                            else:
                                locations.append(locations[-1]+1)
                            xticks.append(cnt)
                            xticklabels.append(params['label'])

                            current_location = params['location']

                            # FIXME handle single data points
                            # get data and standardts
                            if trace_key in df_sorted:
                                trace_series = df_sorted[trace_key][variable]
                                data_array = (_2array(trace_series, remove_nans=True, remove_nans_axis=1)-1)*100.
                                # print data_array
                                # print current_location
                                # if type(data_array)==np.float64:
                                #     current_locations = current_location
                                # else:
                                #     current_locations = current_location*np.ones(len(data_array))

                                # plt.plot(current_locations, data_array, '.', color=params['color'], markersize=markersize)
                                # print trace_key, data_array.shape, type(data_array)

                                if type(data_array)==np.float64:
                                    # mean across slices
                                    data_mean = data_array
                                    #std across slices
                                    data_std = 0.
                                    # sem across slices
                                    data_sem = 0.

                                elif type(data_array)==np.ndarray:
                                    # get stats
                                    # mean across slices
                                    data_mean = np.mean(data_array, axis=0)
                                    #std across slices
                                    data_std = np.std(data_array, axis=0)
                                    # sem across slices
                                    data_sem = stats.sem(data_array, axis=0)

                                heights.append(data_mean)

                                trace_args['color'] = params['color']

                                fig_args.append(trace_args)

                                colors.append(params['color'])

                                plt.errorbar(locations[-1], data_mean, data_sem, color=(.5,.5,.5))

            plt.plot(locations, heights, '.', color=params['color'], markersize=markersize)
            # barcontainer = ax[figure_key].bar(locations, heights, width=bar_width, tick_label=xticklabels)
            # xlim[figure_key] = ax[figure_key].get_xlim()
            # ylim[figure_key] = ax[figure_key].get_ylim()
            # # ax[figure_key].set_xticks(xticks, xticklabels,)
            # # barcontainer = ax[figure_key].violinplot(locations, heights, width=bar_width, tick_label=xticklabels)
            # print 'rotations:', figure_key, figure_params[figure_key]['params']
            # ax[figure_key].set_xticklabels(xticklabels, fontsize
            #     =20, fontweight='heavy', rotation=figure_params[figure_key]['params']['rotation'])
            # for bar_i, bar in enumerate(barcontainer):
            #     bar.set_color(colors[bar_i])
            plt.show(block=False)

    return fig, ax

def _build_figure_params(figures_dict):
    '''
    '''
    # build multiindex for trace parameters
    multi_list = []
    level_names = ['figure','subgroup', 'trace']
    for level_1_key, level_1 in figures_dict.iteritems():
        for level_2_key, level_2 in level_1.iteritems():
            for level_3_i, level_3_key in enumerate(level_2):
                multi_list.append((level_1_key, level_2_key, level_3_key))
    multiindex = pd.MultiIndex.from_tuples(multi_list, names=level_names)
    # build dataframe
    figure_df = pd.DataFrame(index=multiindex, dtype='object')

    return figure_df

def _plot_trace_mean(df_sorted, figdf, variable, **kwargs):
    '''FIXME add docs
    '''
    # FIXME add kwargs to alter figure details
    # create figure groupings (all conditions that will go on the same figure)
    fig={}
    ax={}
    n_subgroups={}
    n_traces={}
    xlim={}
    ylim={}
    # set figdf to hierarchical index (figure, subgroup, trace)
    figdf = figdf.reset_index().set_index(['figure','subgroup','trace'])
    # get list of all figure names
    figures = figdf.index.get_level_values('figure').unique().values
    # iterate over figures
    for figkey in figures:
        # create figure, passing params as **kwargs
        fig[figkey], ax[figkey] = plt.subplots()
        # get subgroups list
        subgroups = figdf.loc[figkey].index.get_level_values('subgroup').unique().values
        # subgroups = figdf[figdf.figure==figkey].subgroup.unique()
        # number of subgroup in figure
        n_subgroups[figkey] = len(subgroups) 
        n_traces[figkey]={}
        # iterate over subgroups of traces
        for subkey in subgroups:
            traces = figdf.loc[(figkey, subkey)].index.get_level_values('trace').unique().values
            # FIXME distribute subgroup padrameters to each trace in the subgroup, with priority to trace parameters
            n_traces[figkey][subkey]=len(traces)

            # iterate over traces
            for tracekey in traces:
                params = figdf.loc[(figkey,subkey, tracekey)]
                # get series from df sorted
                trace_series = df_sorted[tracekey][variable]
                # convert to array
                data_array = _2array(trace_series, remove_nans=True, remove_nans_axis=1)#*1000
                # get stats
                # mean across slices
                data_mean = np.mean(data_array, axis=0)
                #std across slices
                data_std = np.std(data_array, axis=0)
                # sem across slices
                data_sem = stats.sem(data_array, axis=0)
                # time vector
                t = np.arange(len(data_mean))#
                # line plot with shaded error
                if figdf.loc[(figkey,subkey,tracekey)].error_style=='shade':
                    ax[figkey].plot(t, data_mean, color=params.trace_color, linewidth=params.trace_linewidth)
                    plt.fill_between(t, data_mean-data_sem, data_mean+data_sem, color=params.error_color, alpha=params.error_alpha)
                # error bar plot
                elif figdf.loc[(figkey,subkey,tracekey)].error_style=='bar':
                    ax[figkey].errorbar(t, data_mean, yerr=data_sem, color=params.trace_color, 
                        marker=params.trace_marker,  
                        markersize=params.markersize, 
                        elinewidth=params.error_linewidth, 
                        linewidth=params.trace_linewidth, 
                        markerfacecolor=params.trace_color, 
                        ecolor=params.error_color)
        # get x and y limits based data
        xlim[figkey] = ax[figkey].get_xlim()
        ylim[figkey] = ax[figkey].get_ylim()

    fig, ax = _plot_format_figures(fig=fig, ax=ax, figdf=figdf, xlim=xlim, ylim=ylim)

    plt.show(block=False)

    return fig, ax

def _plot_bar2(df_sorted, figdf, variable, group_space=1, bar_width=1, bar_spacing=1):
    '''
    '''
    fig={}
    ax={}
    n_subgroups={}
    n_traces={}
    xlim={}
    ylim={}
    # set figdf to hierarchical index (figure, subgroup, trace)
    figdf = figdf.reset_index().set_index(['figure','subgroup','trace'])
    # get list of all figure names
    figures = figdf.index.get_level_values('figure').unique().values
    # iterate over figures
    for figkey in figures:
        fig[figkey], ax[figkey] = plt.subplots()
        # get subgroups list
        subgroups = figdf.loc[figkey].index.get_level_values('subgroup').unique().values
        n_subgroups[figkey] = len(subgroups) 
        n_traces[figkey]={}
        locations = []
        heights=[]
        colors=[]
        fig_args=[]
        xticks=[]
        xticklabels=[]
        cnt=bar_spacing
        # iterate over subgroups of traces
        for subkey in subgroups:
            traces = figdf.loc[(figkey, subkey)].index.get_level_values('trace').unique().values
            # FIXME distribute subgroup padrameters to each trace in the subgroup, with priority to trace parameters
            n_traces[figkey][subkey]=len(traces)
            cnt+=group_space
            # iterate over traces
            for tracekey in traces:
                param = figdf.loc[(figkey,subkey, tracekey)]
                # get series from df sorted
                trace_series = df_sorted[tracekey][variable]
                # convert to array
                data_array = _2array(trace_series, remove_nans=True, remove_nans_axis=1)#*1000
                # print figdf.keys()
                if all(figdf.loc[figkey].fig_topercent):
                    data_array = 100.*(data_array-1)
                # get stats
                if type(data_array)==np.ndarray:
                    # mean across slices
                    data_mean = np.mean(data_array, axis=0)
                    #std across slices
                    data_std = np.std(data_array, axis=0)
                    # sem across slices
                    data_sem = stats.sem(data_array, axis=0)
                    # add plot location
                    print figkey, subkey, tracekey, param.sub_location, param.trace_location
                    plot_location = param.sub_location+param.trace_location
                    locations.append(plot_location)
                    xticks.append(plot_location)
                    xticklabels.append(param.trace_label)
                    colors.append(param.trace_color)
                    heights.append(data_mean)
                    plt.errorbar(locations[-1], data_mean, data_sem, color=param.error_color)

        barcontainer = ax[figkey].bar(locations, heights, width=param.fig_barwidth, tick_label=xticklabels)

        # get x and y lims
        xlim[figkey] = ax[figkey].get_xlim()
        ylim[figkey] = ax[figkey].get_ylim()
        # rotate ticklabels
        ax[figkey].set_xticklabels(xticklabels, fontsize
            =figdf.loc[figkey].fig_xtick_fontsize.unique()[0], fontweight=figdf.loc[figkey].fig_xtick_fontweight.unique()[0], rotation=figdf.loc[figkey].fig_xtick_rotate.unique()[0])
        # set bar color
        for bar_i, bar in enumerate(barcontainer):
            bar.set_color(colors[bar_i])

    fig, ax = _plot_format_figures(fig=fig, ax=ax, figdf=figdf, xlim=xlim, ylim=ylim)

    plt.show(block=False)

    return fig, ax

def _plot_format_figures(fig, ax, figdf, xlim, ylim):
    '''
    '''
    # get ylim and xlim across all figures
    #-------------------------------------
    xlims=[]
    ylims=[]
    for figkey in ylim:
        xlims.append(xlim[figkey])
        ylims.append(ylim[figkey])
    xlim_all = [min([temp[0] for temp in xlims]), max([temp[1] for temp in xlims])]
    ylim_all = [min([temp[0] for temp in ylims]), max([temp[1] for temp in ylims])]

    # iterate over figures and update attributes
    #---------------------------------------------
    for figkey, axes in ax.iteritems():
        ylim_current = ax[figkey].get_ylim()
        xlim_current = ax[figkey].get_xlim()
        # set ylim across all figures
        #----------------------------
        if all(figdf.fig_ylim_all):
            if 'fig_ylim' in figdf.keys():
                ax[figkey].set_ylim(figdf.loc[figkey].fig_ylim.iloc[0])
            elif 'fig_ymin' in figdf.keys():
                ax[figkey].set_ylim([figdf.loc[figkey].fig_ymin.iloc[0], ylim_all[1]])
            elif 'fig_ymax' in figdf.keys():
                ax[figkey].set_ylim([ylim_all[0], figdf.loc[figkey].fig_ymax.iloc[0]])
            else:
                ax[figkey].set_ylim(ylim_all)

        elif 'fig_ymin' in figdf.loc[figkey] and not any(figdf.loc[figkey].fig_ymin.isnull()):
            ax[figkey].set_ylim([figdf.loc[figkey].fig_ymin.unique()[0], ylim_current[1]])
        # set xlim across all figures
        #---------------------------
        if all(figdf.fig_xlim_all):
            if 'fig_xlim' in figdf.keys():
                ax[figkey].set_xlim(figdf.loc[figkey].fig_xlim.iloc[0])
            elif 'fig_xmin' in figdf.keys():
                ax[figkey].set_xlim([figdf.loc[figkey].fig_xmin.iloc[0], xlim_all[1]])
            elif 'fig_xmax' in figdf.keys():
                ax[figkey].set_xlim([xlim_all[0], figdf.loc[figkey].fig_xmax.iloc[0]])
            else:
                ax[figkey].set_xlim(xlim_all)

        elif 'fig_xmin' in figdf.loc[figkey] and not any(figdf.loc[figkey].fig_xmin.isnull()):
            ax[figkey].set_xlim([figdf.loc[figkey].fig_xmin.unique()[0], xlim_current[1]])

        # turn off axes box
        #------------------
        if all(figdf.loc[figkey].fig_boxoff):
            ax[figkey].spines['right'].set_visible(False)
            ax[figkey].spines['top'].set_visible(False)

        # set axes linewidth and tick position
        #----------------------
        ax[figkey].spines['left'].set_linewidth(figdf.loc[figkey].fig_axes_linewidth.unique()[0])
        ax[figkey].spines['bottom'].set_linewidth(figdf.loc[figkey].fig_axes_linewidth.unique()[0])
        ax[figkey].xaxis.set_ticks_position('bottom')
        ax[figkey].yaxis.set_ticks_position('left')

        # set axes labels
        #----------------
        ax[figkey].set_xlabel(figdf.loc[figkey].fig_xlabel.unique()[0], fontsize=figdf.loc[figkey].fig_xlabel_fontsize.unique()[0], fontweight=figdf.loc[figkey].fig_xlabel_fontweight.unique()[0])
        ax[figkey].set_ylabel(figdf.loc[figkey].fig_ylabel.unique()[0], fontsize=figdf.loc[figkey].fig_ylabel_fontsize.unique()[0], fontweight=figdf.loc[figkey].fig_ylabel_fontweight.unique()[0])

        # set x and y ticks
        #--------------------
        xlim[figkey] = ax[figkey].get_xlim()
        ylim[figkey] = ax[figkey].get_ylim()
        xticks = ax[figkey].get_xticks()
        yticks = ax[figkey].get_yticks()
        if 'fig_dyticks' in figdf:
            yticks = np.arange(ylim[figkey][0], ylim[figkey][1], figdf.loc[figkey].fig_dyticks.unique()[0])
        elif 'fig_nyticks' in figdf:
            dstep = float(abs(ylim[figkey][1]-ylim[figkey][0]))/figdf.loc[figkey].fig_nyticks.unique()[0]
            yticks = np.arange(ylim[figkey][0], ylim[figkey][1], dstep)

        if 'fig_dxticks' in figdf:
            xticks = np.arange(xlim[figkey][0], xlim[figkey][1], figdf.loc[figkey].fig_dxticks.unique()[0])
        elif 'fig_nxticks' in figdf:
            dstep = float(abs(xlim[figkey][1]-xlim[figkey][0]))/figdf.loc[figkey].fig_nxticks.unique()[0]
            xticks = np.arange(xlim[figkey][0], xlim[figkey][1], dstep)
        
        ax[figkey].set_yticks(yticks)
        ax[figkey].set_xticks(xticks)

        # set ticklabel attributes
        #--------------------------
        for temp in ax[figkey].get_xticklabels():
            temp.set_fontweight(figdf.loc[figkey].fig_xtick_fontweight.unique()[0])
            temp.set_fontsize(figdf.loc[figkey].fig_xtick_fontsize.unique()[0])
        for temp in ax[figkey].get_yticklabels():
            temp.set_fontweight(figdf.loc[figkey].fig_ytick_fontweight.unique()[0])
            temp.set_fontsize(figdf.loc[figkey].fig_ytick_fontsize.unique()[0])

        # set tight layout
        #-------------------
        if all(figdf.loc[figkey].fig_tight_layout):
            plt.figure(fig[figkey].number)
            plt.tight_layout()

    return fig, ax

class ApplyDF:
    ''' functions to be passed to df to create new columns from existing data
    '''
    def _normalize_column(self, df, colnames=[], colnorm=[],):
        '''
        '''
        # get list of columns to be normalized
        def _norm(row, colnames=colnames, colnorm=colnorm):
            ''' normalize data in colnames columns based on data in colnorm column
            '''
            col_add = '_norm_'+colnorm
            norm = row[colnorm]
            for col in colnames:
                new_data =np.array(row[col])/norm
                row[col+col_add] = new_data
            return row

        kwargs = {'colnames':colnames, 'colnorm':colnorm}
        df = df.apply(_norm, axis=1, **kwargs)
        return df

    def _get_io(self, df, max_idx_col=[], colnames=[], colkeys=[], colkeys_exclude=[]):
        ''' get data for specific pulse from induction data
        '''
        print 'getting io'
        # print colnames
        colnames=copy.copy(colnames)

        # get list of columns to be normalized
        def _get_io(row, max_idx_col=max_idx_col, colnames=colnames, colkeys=colkeys, colkeys_exclude=colkeys_exclude):
            ''' normalize data in colnames columns based on data in colnorm column
            '''
            if len(colkeys)>0:
                columns = row.index.values
                for col in columns:
                    if all([colkey in col for colkey in colkeys]) and not any([temp in col for temp in colkeys_exclude]):
                        
                        colnames.append(col)

            coladd_io_norm = '_io_norm'
            coladd_io_max = '_io_max'
            for col in colnames:
                if not np.isnan(row[max_idx_col]):
                    max_idx = int(row[max_idx_col])
                    # print max_idx
                    max_val = row[col][max_idx]
                    io_vals = row[col][:max_idx+1]/max_val
                    row[col+coladd_io_norm] = io_vals
                    row[col+coladd_io_max] = max_val
                else:
                    row[col+coladd_io_norm]=np.nan
                    row[col+coladd_io_max]=np.nan

            return row

        kwargs = {'colnames':colnames, 'colkeys':colkeys, 'colkeys_exclude':colkeys_exclude, 'max_idx_col':max_idx_col}
        df = df.apply(_get_io, axis=1, **kwargs)
        return df

    def _get_pulse(self, df, colnames=[], colkeys=[], colkeys_exclude=[]):
        ''' get data for specific pulse from induction data
        '''
        # get list of columns to be normalized
        def _get_pulse(row, colnames=colnames, colkeys=colkeys, colkeys_exclude=colkeys_exclude):
            ''' normalize data in colnames columns based on data in colnorm column
            '''
            if len(colkeys)>0:
                columns = row.index.values
                for col in columns:
                    if all([colkey in col for colkey in colkeys]) and not any([colkey in col for colkey in colkeys_exclude]):
                        colnames.append(col)

            coladds = ['_first', '_second', '_third', '_fourth']
            induction_protocol = row.induction_pattern_0
            if induction_protocol=='TBS':
                pulses=4
            elif induction_protocol=='weak5Hz':
                pulses=1
            elif induction_protocol=='nostim':
                pulses=0

            for col in colnames:
                for coladd_i, coladd in enumerate(coladds):
                    # handle nans
                    if type(row[col])==list or type(row[col])==np.ndarray and any(np.isnan(row[col])):
                            new_data=np.nan
                    elif type(row[col])==float and np.isnan(row[col]):
                        new_data=np.nan
                    elif type(row[col])==float:
                        new_data=row[col]
                    else:
                        if pulses==0:
                            mask=[]
                        elif pulses==1:
                            mask = range(0, len(row[col]), pulses)
                        else:
                            mask = range(coladd_i, len(row[col]), pulses)
                    
                        new_data = np.array(row[col])[mask]
                    row[col+coladd] = new_data
            return row

        kwargs = {'colnames':colnames, 'colkeys':colkeys, 'colkeys_exclude':colkeys_exclude}
        df = df.apply(_get_pulse, axis=1, **kwargs)
        return df

    def _mean(self, df, colnames=[], colkeys=[], colkeys_exclude=[]):
        ''' get data for specific pulse from induction data
        '''
        # get list of columns to be normalized
        def _mean(row, colnames=colnames, colkeys=colkeys, colkeys_exclude=colkeys_exclude):
            ''' normalize data in colnames columns based on data in colnorm column
            '''
            if len(colkeys)>0:
                columns = row.index.values
                for col in columns:
                    if all([colkey in col for colkey in colkeys]) and not any([colkey in col for colkey in colkeys_exclude]):
                        colnames.append(col)

            coladd = '_mean'
            for col in colnames:
                row[col+coladd] = np.mean(np.array(row[col]))
            return row

        kwargs = {'colnames':colnames, 'colkeys':colkeys, 'colkeys_exclude':colkeys_exclude}
        df = df.apply(_mean, axis=1, **kwargs)
        return df

    def _last_burst(self, df, colnames=[], colkeys=[], colkeys_exclude=[]):
        ''' get data for specific pulse from induction data
        '''
        print 'getting last burst'
        # get list of columns to be normalized
        def _last_burst(row, colnames=colnames, colkeys=colkeys, colkeys_exclude=colkeys_exclude):
            ''' normalize data in colnames columns based on data in colnorm column
            '''
            colnames=copy.copy(colnames)
            if len(colkeys)>0:
                columns = row.index.values
                for col in columns:
                    if all([colkey in col for colkey in colkeys]) and not any([temp in col for temp in colkeys_exclude]):
                        
                        colnames.append(col)

            coladd = '_lastburst'
            for col in colnames:
                # print row[col]
                if type(row[col])==float and np.isnan(row[col]):
                    row[col+coladd]=np.nan
                elif type(row[col])==list and any(np.isnan(row[col])):
                    row[col+coladd]=np.nan
                else:
                    # print col
                    row[col+coladd] = np.array(row[col])[-1]
            return row

        kwargs = {'colnames':colnames, 'colkeys':colkeys, 'colkeys_exclude':colkeys_exclude}
        df = df.apply(_last_burst, axis=1, **kwargs)
        return df
    
    def _first_burst(self, df, colnames=[], colkeys=[], colkeys_exclude=[]):
        ''' get data for specific pulse from induction data
        '''
        print 'getting first burst'
        # print colnames
        colnames=copy.copy(colnames)

        # get list of columns to be normalized
        def _first_burst(row, colnames=colnames, colkeys=colkeys, colkeys_exclude=colkeys_exclude):
            ''' normalize data in colnames columns based on data in colnorm column
            '''
            if len(colkeys)>0:
                columns = row.index.values
                for col in columns:
                    if all([colkey in col for colkey in colkeys]) and not any([temp in col for temp in colkeys_exclude]):
                        
                        colnames.append(col)
            # print colnames

            coladd = '_firstburst'
            for col in colnames:
                # print row[col]
                if type(row[col])==float and np.isnan(row[col]):
                    row[col+coladd]=np.nan
                elif type(row[col])==list and any(np.isnan(row[col])):
                    row[col+coladd]=np.nan
                else:
                    # print col
                    row[col+coladd] = np.array(row[col])[0]
            return row

        kwargs = {'colnames':colnames, 'colkeys':colkeys, 'colkeys_exclude':colkeys_exclude}
        df = df.apply(_first_burst, axis=1, **kwargs)
        return df

    def _fft(self, df, colnames=[], colkeys=[], colkeys_exclude=[], slice_i=slice(None), unwrap=False, ):
        ''' get data for specific pulse from induction data
        '''
        print 'getting fft'
        # get list of columns to be normalized
        def _fft(row, colnames=colnames, colkeys=colkeys, colkeys_exclude=colkeys_exclude, slice_i=slice_i, unwrap=False):
            ''' normalize data in colnames columns based on data in colnorm column
            '''
            if len(colkeys)>0:
                columns = row.index.values
                for col in columns:
                    if all([colkey in col for colkey in colkeys]) and not any([colkey in col for colkey in colkeys_exclude]):
                        colnames.append(col)

            coladd_abs = '_fft_abs'
            coladd_angle='_fft_angle'
            
            for col in colnames:
                print np.array(row[col]).squeeze()
                # print np.array(row[col]).squeeze()==np.nan
                # print type(np.isnan(np.array(row[col]).squeeze()))
                if type(np.isnan(np.array(row[col]).squeeze()))==np.bool_ and np.isnan(np.array(row[col]).squeeze()):
                    row[col+coladd_abs]=np.nan
                    row[col+coladd_angle]=np.nan
                elif any(np.isnan(np.array(row[col]).squeeze())):
                    row[col+coladd_abs]=np.nan
                    row[col+coladd_angle]=np.nan

                else:
                    slice_is = slice_i.indices(len(np.array(row[col])))
                    if len(range(*slice_is))==0:
                        slicer = slice(None)
                    else:
                        slicer=copy.copy(slice_i)
                    row[col+coladd_abs] = np.abs(np.fft.fft(np.array(row[col])[slicer], n=60))
                    if unwrap:
                        row[col+coladd_angle] = np.unwrap(np.angle(np.fft.fft(np.array(row[col])[slicer], n=60)))
                    else:
                        row[col+coladd_angle] = np.angle(np.fft.fft(np.array(row[col])[slicer], n=60))

            return row

        kwargs = {'colnames':colnames, 'colkeys':colkeys, 'colkeys_exclude':colkeys_exclude, 'slice_i':slice_i, 'unwrap':unwrap}
        df = df.apply(_fft, axis=1, **kwargs)
        return df

    def _get_nth(self, df, colnames=[], colkeys=[], colkeys_exclude=[], n=0):
        ''' get data for specific pulse from induction data
        '''
        print 'getting fft'
        # get list of columns to be normalized
        def _get_nth(row, colnames=colnames, colkeys=colkeys, colkeys_exclude=colkeys_exclude, n=n):
            ''' normalize data in colnames columns based on data in colnorm column
            '''
            if len(colkeys)>0:
                columns = row.index.values
                for col in columns:
                    if all([colkey in col for colkey in colkeys]) and not any([colkey in col for colkey in colkeys_exclude]):
                        colnames.append(col)

            coladd = '_'+str(n)+'th'
            
            for col in colnames:
                print np.array(row[col]).squeeze()
                # print np.array(row[col]).squeeze()==np.nan
                # print type(np.isnan(np.array(row[col]).squeeze()))
                if type(np.isnan(np.array(row[col]).squeeze()))==np.bool_ and np.isnan(np.array(row[col]).squeeze()):
                    row[col+coladd]=np.nan
                    row[col+coladd]=np.nan
                elif any(np.isnan(np.array(row[col]).squeeze())):
                    row[col+coladd]=np.nan
                    row[col+coladd]=np.nan

                else:
                    if len(np.array(row[col]))>n:
                        row[col+coladd] = np.array(row[col])[n]
                    # slice_is = slice_i.indices(len(np.array(row[col])))
                    # if len(range(*slice_is))==0:
                    #     slicer = slice(None)
                    # else:
                    #     slicer=copy.copy(slice_i)
                    # row[col+coladd_abs] = np.abs(np.fft.fft(np.array(row[col])[slicer], n=60))
                    # if unwrap:
                    #     row[col+coladd_angle] = np.unwrap(np.angle(np.fft.fft(np.array(row[col])[slicer], n=60)))
                    # else:
                    #     row[col+coladd_angle] = np.angle(np.fft.fft(np.array(row[col])[slicer], n=60))

            return row

        kwargs = {'colnames':colnames, 'colkeys':colkeys, 'colkeys_exclude':colkeys_exclude, 'n':n}
        df = df.apply(_get_nth, axis=1, **kwargs)
        return df






            
