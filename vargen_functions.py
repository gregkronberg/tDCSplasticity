"""
analysis

"""
from __future__ import division
import numpy as np
import analysis
# import scipy
# import scipy.io
from scipy import stats
from scipy import signal
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import itertools as it
import os
import glob
import cPickle as pickle
import param
import math
import run_control
import copy
import matplotlib.patches as patches
import matplotlib.lines as mlines
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib import cm as colormap
import itertools
import inspect
from matplotlib.ticker import FormatStrFormatter
import sys

kwargs = {'experiment':'exp_1a_dose'}
# specify clopath parameters to be shared across experiments
#------------------------------------------------------------
clopath_param_all= {
    'A_m0':10*100E-5, # depression magnitude parameter (mV^-1)
    'tetam':-70,#-41, # depression threshold (mV)
    'tetap':-67,#-38, # potentiation threshold (mV)
    'tau_x':8,#-38, # time constant for presynaptic filter (ms)
    'tau_m':20,#-38, # time constant for depression lowpass filter
    'tau_p': 3, # time constant (ms) for low pass filter post membrane potential for potentiation
    'A_p':4*100E-5, # amplitude for potentiation (mV^-2)
    'delay':0, # conduction delay (ms)
    'LTD_delay':1, # conduction delay for LTD  
    'lower_bound':0.0,
    'upper_bound':3.0,
    }

clopath_param_all_temp= {
    'A_m0':10*100E-5, # depression magnitude parameter (mV^-1)
    'tetam':-70,#-41, # depression threshold (mV)
    'tetap':-69,#-38, # potentiation threshold (mV)
    'tau_x':8,#-38, # time constant for presynaptic filter (ms)
    'tau_m':20,#-38, # time constant for depression lowpass filter
    'tau_p': 3, # time constant (ms) for low pass filter post membrane potential for potentiation
    'A_p':4*100E-5, # amplitude for potentiation (mV^-2)
    'delay':0, # conduction delay (ms)
    'LTD_delay':1, # conduction delay for LTD  
    'lower_bound':0.0,
    'upper_bound':None,
    }

# clopath_param_all2={
#     'A_m0':3*10E-3, # depression magnitude parameter (mV^-1)
#     'tetam':-70,#-41, # depression threshold (mV)
#     'tetap':-67,#-38, # potentiation threshold (mV)
#     'tau_x':8,#-38, # time constant for presynaptic filter (ms)
#     'tau_m':12,#-38, # time constant for depression lowpass filter
#     'tau_p': 3, # time constant (ms) for low pass filter post membrane potential for potentiation
#     'A_p':6E-3, # amplitude for potentiation (mV^-2)
#     'delay':0, # conduction delay (ms)
#     'LTD_delay':1, # conduction delay for LTD  
#     'lower_bound':0.,
#     }

#########################
# keep these!!!!!!!
###########################
clopath_param_all2={
    'A_m0':4*10E-3, # depression magnitude parameter (mV^-1)
    'tetam':-70,#-41, # depression threshold (mV)
    'tetap':-60,#-38, # potentiation threshold (mV)
    'tau_x':15,#-38, # time constant for presynaptic filter (ms)
    'tau_m':40,#-38, # time constant for depression lowpass filter
    'tau_p': 20, # time constant (ms) for low pass filter post membrane potential for potentiation
    'A_p':40E-3, # amplitude for potentiation (mV^-2)
    'delay':0, # conduction delay (ms)
    'LTD_delay':1, # conduction delay for LTD  
    'lower_bound':0.0,
    'upper_bound':3.0,
    }

clopath_param_all_stdp={
    # stdp parameters
    #------------
    'A_m0':10*100E-5, # depression magnitude parameter (mV^-1)
    'tetam':-71,#-41, # depression threshold (mV)
    'tetap':-50,#-38, # potentiation threshold (mV)
    'tau_x':10,#-38, #  time constant for presynaptic filter (ms)
    'tau_m':20,#-38, # # time constant for depression lowpass filter
    'tau_p': 30, # time constant (ms) for low pass filter post membrane potential for potentiation
    'A_p': 16*100E-5, # amplitude for potentiation (mV^-2)
    'delay':0, # conduction delay (ms)
    'LTD_delay':1,# conduction delay for LTD  
    }

clopath_param_all_stdp2={
    # stdp parameters
    #------------
    'A_m0':2*10*100E-5, # depression magnitude parameter (mV^-1)
    'tetam':-71,#-41, # depression threshold (mV)
    'tetap':-50,#-38, # potentiation threshold (mV)
    'tau_x':10,#-38, #  time constant for presynaptic filter (ms)
    'tau_m':20,#-38, # # time constant for depression lowpass filter
    'tau_p': 30, # time constant (ms) for low pass filter post membrane potential for potentiation
    'A_p': 6*16*100E-5, # amplitude for potentiation (mV^-2)
    'delay':0, # conduction delay (ms)
    'LTD_delay':1,# conduction delay for LTD  
    }

class VarFuncs:
    ''' extract higher level analysis variables from individual cells/experiments.  outputs should be a dataframe for each preprocessed file, which are then combined (see VarGen) into a single dataframe for all experiments.
    '''
    def __init__(self, ):
        '''
        '''
        pass

    def _build_conditions_df(self, pre, var, **kwargs):
        ''' standard dataframe containing common simulation details
        Arguments
        --------------
        ~pre: preprocessed data file
        ~var: variable type that is being converted to df (e.g. 'v')

        Return
        -------------
        ~df_new: dataframe with indices corresponding to traces in pre[var]['data'].  Columns correspond to conditions, which are specified below

        Comments
        --------------
        '''
        # iterate over traces
        df_new = pd.DataFrame(dtype='object')
        df_new['trial_id'] = pre[var]['trial_id']
        df_new['location'] = pre[var]['locations']
        df_new['field'] = pre[var]['field']
        df_new['tree'] = zip(*pre[var]['locations'])[0]
        df_new['sec'] = zip(*pre[var]['locations'])[1]
        df_new['seg'] = zip(*pre[var]['locations'])[2]
        df_new['dt'] = pre[var]['p']['dt']
        df_new['tstop'] = pre[var]['p']['tstop']
        # ac field parameters
        #---------------------
        if 'ac_field' in pre[var]['p']:
            for key, val in pre[var]['p'].iteritems():
                if 'ac_field' in key:
                    df_new[key]=val

        # create columns for parameters in p
        if 'params' in kwargs:
            # Store all parameters
            if kwargs['params']==None or kwargs['params']=='all':
                for param in pre[var]['p'].keys():
                    if param not in df_new:

                        df_new[param]=None
                        for i in df_new.index:
                            value = pre[var]['p'][param]
                            if type(value)==list or type(value)==dict:
                                if len(value)==0:
                                    value=np.nan
                                df_new.at[i,param] = value
                            else:
                                df_new.at[i,param] = value

            else:
                # store specified parameters
                for param in kwargs['params']:
                    if param not in df_new:
                        df_new[param]=None
                        for i in df_new.index:
                            df_new.at[i,param]=pre[var]['p'][param]

        # iterate over traces
        df_new['paths']=None
        df_new['path']=None
        df_new['input_times']=None
        for i, row in df_new.iterrows():
            location = df_new['location'][i]
            # get distance from soma
            #------------------------
            dist = pre[var]['p']['seg_dist'][row.tree][row.sec][row.seg]
            # get morphology entry (overall seg index, location key, x, y, z, dimaeter, parent overall index)
            #----------------------
            morpho = pre[var]['p']['morpho'][row.tree]
            morpho_key = '_'.join([row.tree, str(row.sec), str(row.seg)])
            morpho_val = [seg for sec in morpho for seg in sec if morpho_key in seg][0]
            # get segment distance from soma
            #------------------------
            df_new.at[i, 'seg_dist']=dist
            if 'morpho' not in df_new:
                df_new['morpho']=None
            df_new.at[i, 'morpho']=morpho_val
            # get synaptic input times
            #--------------------------
            # pre['input_times']['data'] is nested list as [burst_i][input_time], locations list contains an entry for each burst (locations are repeated if it receives multiple bursts)
            if 'input_times' in pre:
                if 'input_times' not in df_new:
                    df_new['input_times']=None
                if df_new['location'][i] in pre['input_times']['locations']:
                    # indices of input times for each burst
                    input_time_is = [temp_i for temp_i, temp in enumerate(pre['input_times']['locations']) if temp==location]
                    # print input_time_is
                    # nested list of input times for current location [bursts][input times]
                    input_times = [pre['input_times']['data'][temp_i] for temp_i in input_time_is]
                    # print pre['input_times']['data']
                    # print 'times', input_times
                    input_times = np.array(list((set(np.concatenate(input_times)))))
                    # print input_times

                    # print input_times.shape
                    # flattened array of input times for current location
                    input_times = np.array(input_times).flatten()
                    # print input_times
                    # print input_times.shape
                    # sort input times
                    input_times.sort()
                    # store in df
                    df_new.at[i, 'input_times']=input_times
            # synaptic input path parameters
            #----------------------------------
            paths=[]
            n_paths=len(pre[var]['p']['p_path'].keys())
            df_new.at[i, 'n_paths']=n_paths
            # df_new['paths']=None
            # iterate over input pathways
            for path_key, path in pre[var]['p']['p_path'].iteritems():
                # list of paths that contain the current location
                #------------------------------------------------
                if row.location in path['syn_idx']:
                    temp_i = path['syn_idx'].index(row.location)
                    paths.append(path_key)
                    # synaptic weight at current location
                    #------------------------------------
                    if 'w_idx' in path:
                        df_new.at[i, 'w'] = path['w_idx'][temp_i]
                    df_new.at[i,'path']=path_key
                # iterate over path parameters
                for path_param, param_val in path.iteritems(): 
                    colname = '_'.join(['path', path_key, path_param])
                    if colname not in df_new:
                        df_new[colname]=None
                    # add all path parameters to df
                    #------------------------------
                    df_new.at[i,colname]=param_val
            # add paths that current entry belongs to 
            df_new.at[i,'paths']=paths
        # replace none objects with na
        df_new.fillna(value=pd.np.nan, inplace=True)

        return df_new

    def _get_vtrace(self, pre, df=[], keep=[], **kwargs):
        '''
        '''
        df_new = self._build_conditions_df(pre, var='v',**kwargs)
        for i, val in enumerate(pre['v']['data']):
            if 'data_v' not in df_new:
                df_new['data_v']=None
            # print val.shape
            if len(val.shape)>1:
                print val.shape
            df_new.at[i,'data_v'] = np.array(val)

        return df_new

    def _get_spikes(self, pre, df=[], keep=[], threshold=-30, **kwargs):
        '''
        '''
        print 'threshold', threshold
        # get time step
        dt = pre['v']['p']['dt']
        # build conditions df
        df_new = self._build_conditions_df(pre, var='v', **kwargs)
        # preallocate columns
        df_new['spike_idx']=None
        df_new['spike_train']=None
        df_new['spike_times']=None
        # iterate over voltage traces
        for i, data in enumerate(pre['v']['data']):
            # get spike info
            spike_idx, spike_train, spike_times = analysis.Spikes()._get_spikes(
                data=data, dt=dt, threshold=threshold)
            # add to dataframe
            # spike_idx: list of indices where spike onset occurs, length=number of spikes
            df_new.at[i, 'spike_idx']=spike_idx
            # spike_train: boolean array with same shape as voltage data, where spike onsets are stored as 1
            df_new.at[i, 'spike_train']=spike_train
            # spike_times: list of spike times in ms, length=number of spikes 
            df_new.at[i, 'spike_times']=spike_times

        return df_new

    def _get_w_clopath(self, pre, df=[], keep=[], **kwargs):
        '''
        '''
        df_new = self._build_conditions_df(pre, var='input_times', **kwargs)
        # get corresponding indices in voltage list
        voltage_is = [i for i, loc in enumerate(pre['v']['locations']) if loc in pre['input_times']['locations'] ]
        # convert voltage data to array
        voltage_array = np.array(pre['v']['data'])[voltage_is, :]

        # convert input times to input array
        #--------------------------------------------------------------
        # preallocate
        input_array = np.zeros(voltage_array.shape)
        # time vector
        t = pre['v']['t'] 
        # cut off number of decimals to prevent rounding error
        t = np.around(t, 4)
        # iterate over locations
        for loc_i, input_times in enumerate(pre['input_times']['data']):

            # find input times and to boolean vector
            # iterate over input times
            for t_i, t_t in enumerate(input_times):
                # if input time does not exactly match a sample in time vector    
                if t_t not in t:
                    # set input time to nearest sample
                    input_array[loc_i, int(np.argmin(np.absolute(t-t_t)))]=1
                else:
                    # store input times
                    input_array[loc_i, np.where(t==t_t)]=1
        # get parameters
        param = pre['input_times']['p']
        # update with clopath parameters
        if 'clopath_param' in kwargs:
            param.update(kwargs['clopath_param'])
        # get weight array
        w_clopath = analysis.Clopath()._clopath(x=input_array, u=voltage_array, fs=1/param['dt'], w0=0.5, homeostatic=False, param=param)
        # preallocate column for storing weight time series
        if 'w_clopath' not in df_new:
            df_new['w_clopath']=None
        # iterate over input vectors and store corresponding weight timeseries
        for i, x in enumerate(pre['input_times']['data']):
            df_new.at[i, 'w_clopath'] = copy.copy(w_clopath[i,:])
            df_new.at[i, 'w_final'] = w_clopath[i,-1]
            df_new.at[i, 'w_initial']=w_clopath[i, 0]
            df_new.at[i, 'dw_clopath'] = w_clopath[i,-1]/w_clopath[i,0]

        return df_new

class VarGen:
    ''' apply VarFuncs to generate group variables
    '''
    def __init__(self, ):
        '''
        '''
        pass

    def _process_new_data_df(self, group_df, preprocessed_directory, search_string='data', functions=[], kwlist=[], rerun=[], keep=[], file_limit=[], ):
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
        # iterate over columns in group_df
        # if column in keep, use old column values and drop suffix
        # otherwise keep new columnn values and drop suffix


        # merge all function dfs on filename and path
        # update appropriate rows and columns of group_df with function dfs, use merge? find indices and use update/combine?

        # get list of all preprocessed data files
        print 'search_string', search_string
        data_files = glob.glob(preprocessed_directory+'*'+search_string+'*')
        # remove directory and file extension
        data_filenames = [file.split('\\')[-1].split('.')[0] for file in data_files]
        trial_ids = [file.split('_')[-1] for file in data_filenames]

        # get list of processed data files (already added to group dataframe)
        if group_df.empty:
            processed_trial_ids = []
        else:
            processed_trial_ids = group_df.trial_id

        # get list of new data files (haven't been added to dataframe)
        new_trial_ids = list(set(trial_ids)-set(processed_trial_ids))
        new_data_files = [data_files[trial_id_i] for trial_id_i, trial_id in enumerate(trial_ids) if trial_id in new_trial_ids]

        # if any functions specified in rerun, set iterfiles to all data files
        if rerun:
            iterfiles = data_files
        # if there are no functions to rerun, only iterate through new files
        else:
            iterfiles = new_data_files

        print 'total data files:', len(data_files) 
        print 'new data files:', len(new_data_files)
        if rerun:
            print 'functions to rerun over all slices:', [f.__name__ for f in rerun]

        # iterate over new files and update group_data structure
        #`````````````````````````````````````````````````````````````````
        print 'updating group data structure'

        # temporary df for storing rerun function outputs
        df_update = pd.DataFrame()

        # suffixes for merging df_update with group_df
        old_suffix = '_old'
        rerun_suffix = '_rerun'

        print 'iter files:', len(iterfiles)

        # iterate over data files
        for file_i, file in enumerate(iterfiles):
            # print file_i

            # apply file limit
            #------------------
            if file_limit and file_i> file_limit:
                continue
            else:
                print 'processing file', file_i
                # get trial_id
                #----------------
                trial_id = file.split('\\')[-1].split('.')[0].split('_')[-1]

                # load data file
                #---------------
                with open(file, 'rb') as pkl_file:

                    pre = pickle.load(pkl_file)

                # for each function create df and add to list (either new or rerun)
                new_dfs = []
                rerun_dfs = []
                # iterate over functions
                #-------------------------
                for func_i, func in enumerate(functions):

                    # if new data file add to new list
                    #---------------------------------
                    if file in new_data_files:

                        print 'new file found:', trial_id
                        # print kwlist[func_i]
                        df_temp = func(pre, df=[], keep=[], **kwlist[func_i])
                        # # print 'test'
                        new_dfs.append(df_temp)

                    # if rerunning current func, add to rerun list
                    #-------------------------------------------
                    elif func in rerun:
                        print 'rerunning function', func.__name__, 'on file', file_i
                        # print kwlist[func_i]
                        # get corresponding row from group_df
                        df_locs = group_df[group_df.trial_id==trial_id].index.values
                        # get df from function
                        df_temp = func(pre, keep=keep, df=group_df.loc[df_locs], **kwlist[func_i])
                        # add to rerun list
                        rerun_dfs.append(df_temp)

                # if there are new df's, append to bottom of df_update
                #--------------------------------------------------
                if new_dfs:
                    # combine all new_df's in list
                    new_df = reduce(
                        lambda left,right: 
                        pd.merge(
                            left, 
                            right[np.append(right.columns.difference(left.columns).values, ['trial_id', 'field', 'location'])], 
                            on=['trial_id', 'field', 'location'], 
                            how='inner', ), 
                        new_dfs
                        )
                    # add to update df (will be merged with old df later)
                    df_update = df_update.append(new_df, ignore_index=True)

                # if there are rerun df's, append to bottom of df_update
                #--------------------------------------------------
                if rerun_dfs:
                    # merge df's across the different functions
                    rerun_df = reduce(
                        lambda left,right: 
                        pd.merge(
                            left, 
                            right[np.append(right.columns.difference(left.columns).values, ['trial_id', 'field', 'location'])], 
                            on=['trial_id', 'field','location'], 
                            how='inner', ), 
                        rerun_dfs
                        ) # validate='one_to_one'

                    # add to update df (will be merged with old df later)
                    df_update = df_update.append(rerun_df, ignore_index=True)
        
        # merge old and new df's
        #------------------------------------------------------------            
        # print 'merging old and new data'
        # if no old df is specified, set to group_df to update_df
        if group_df.empty:
            group_df=df_update
        # otherwise merge old df and update df
        elif not df_update.empty:
            group_df = group_df.merge(df_update, on=['trial_id', 'field', 'location'], how='outer', suffixes=[old_suffix, rerun_suffix])
        
        # update columns after merge
        #------------------------------
        # print 'dropping columns after merge'
        # get column values
        columns = group_df.columns.values
        # drop columns, based on keep argument
        drop_columns = []
        # iterate over columns
        for column in columns:
            # if old column
            if old_suffix in column:
                # get original column name
                column_original = column.split(old_suffix)[0]
                # get equivalent rerun column name
                column_new = column_original+rerun_suffix
                # if not a keep column
                if not any([temp for temp in keep if temp in column]):
                    # use new to update old
                    # FIXME, use overwrite kwarg for update
                    group_df[column].update(group_df[column_new])
                    # then drop new
                    drop_columns.append(column_new)
                # if a keep column
                else:
                    # use old to update new, drop old
                    group_df[column_new].update(group_df[column])
                    drop_columns.append(column)
        # drop specified columns
        group_df.drop(drop_columns, axis=1, inplace=True)

        # remove suffixes from merge
        #---------------------------
        # print 'renaming columns after merge'
        # dictionary for updated column names
        rename_dict={}
        # iterate over column names
        for column in group_df.columns.values:
            # check for suffixes and drop 
            if old_suffix in column:
                rename_dict[column] = column.split(old_suffix)[0]
            elif rerun_suffix in column:
                rename_dict[column] = column.split(rerun_suffix)[0]
        # rename columns
        group_df.rename(columns=rename_dict, inplace=True)

        return group_df

    def _vargen(self, variable, directory, functions=[], kwlist=[], rerun=[], keep=[], file_limit=[], write_protocol='.pkl', save=False, size_limit=2E9, **kwargs):
        '''
        '''
        filename = variable+write_protocol
        df = analysis._load_group_data(directory=directory,filename=variable)
        # df = df.reset_index()
        df_copy = copy.deepcopy(df)

        # print sorted(df.keys())
        df = self._process_new_data_df(group_df=df, preprocessed_directory=directory, functions=functions, kwlist=kwlist, rerun=rerun, keep=keep, file_limit=file_limit, **kwargs)
        # print df
        # print df_copy.equals(df)

        df = analysis._save_group_data(df=df, directory=directory, variable=variable, extension='.pkl', check_against=None, **kwargs)

        # # get size of df
        # #---------------------------
        # df_size = sys.getsizeof(df)
        # # if df is too big, break into manageable parts
        # if df_size>size_limit:
        #     print 'df too large to save, breaking into smaller parts'
        #     n_dfs=df_size/size_limit
        #     df_split = np.array_split(df, n_dfs)
        #     for df_i, df_chunk in enumerate(df_split):
        #         # print df_chunk
        #         df_chunk.to_pickle(directory+variable+'_'+str(df_i)+'.pkl')

        # # save group data
        # #----------------

        # elif save:
        #     # check write protocol
        #     if write_protocol=='.pkl':
        #         df.to_pickle(directory+filename)
        #     elif write_protocol=='.h5':
        #         df.to_hdf(directory+filename, key=variable, mode='w', format='fixed')

        # elif not df_copy.equals(df):
        #     print 'saving updated group data'
        #     # check write protocol
        #     if write_protocol=='.pkl':
        #         df.to_pickle(directory+filename)
        #     elif write_protocol=='.h5':
        #         df.to_hdf(directory+filename, key=variable, mode='w', format='fixed')
        return df

    def _vtrace(self, **kwargs):
        '''
        '''
        varfuncs=VarFuncs()
        directory = 'Data/'+kwargs['experiment']+'/'
        variable='vtrace_df'
        functions=[varfuncs._get_vtrace]
        kwlist=[{}]
        rerun=[varfuncs._get_vtrace]
        keep=[]
        file_limit=[]
        df = self._vargen(variable=variable, directory=directory, functions=functions, kwlist=kwlist, rerun=rerun, keep=keep, file_limit=file_limit)
        return df

    def _spikes(self, **kwargs):
        '''
        '''
        varfuncs = VarFuncs()
        directory = 'Data/'+kwargs['experiment']+'/'
        variable='spikes_df'
        functions=[varfuncs._get_spikes]
        kwlist=[{'threshold':-30}]
        rerun=[varfuncs._get_spikes]
        keep=[]
        file_limit=[]
        df = self._vargen(variable=variable, directory=directory, functions=functions, kwlist=kwlist, rerun=rerun, keep=keep, file_limit=file_limit)
        return df

    def _w_clopath(self, **kwargs):
        '''
        '''
        varfuncs=VarFuncs()
        directory = 'Data/'+kwargs['experiment']+'/'
        variable='w_clopath_df'
        functions=[varfuncs._get_w_clopath]
        kwlist=[{'clopath_param':clopath_param_all}]
        # kwlist=[kwargs]
        rerun=[varfuncs._get_w_clopath]
        keep=[]
        file_limit=[]
        df = self._vargen(variable=variable, directory=directory, functions=functions, kwlist=kwlist, rerun=rerun, keep=keep, file_limit=file_limit)
        return df

    def _w_clopath_stdp(self, **kwargs):
        '''
        '''
        varfuncs=VarFuncs()
        directory = 'Data/'+kwargs['experiment']+'/'
        variable='w_clopath_stdp_df'
        functions=[varfuncs._get_w_clopath]
        kwlist=[{'clopath_param':clopath_param_all_stdp}]
        rerun=[varfuncs._get_w_clopath]
        keep=[]
        file_limit=[]
        df = self._vargen(variable=variable, directory=directory, functions=functions, kwlist=kwlist, rerun=rerun, keep=keep, file_limit=file_limit)
        return df

class DfFuncs:
    ''' applied to create new columns after df has been generated, all function should take a dataframe as an argument and return the same dataframe with new columns added

    # FIXME add functionality to check for data that has already been processed to speed up analysis
    '''
    def __init__(self, ):
        '''
        '''
        pass
    
    def _get_area(self, df, index, **kwargs):
        '''
        '''
        # reorganize dataframe
        df = df.set_index(['trial_id','location','field'])
        # column name to add
        coladd = ''
        if coladd not in df or rerun:
            # iterate over rows
            for row_i, row in df.iterrows():
                # time point to measure polarization at
                time_i = int(time/row['dt'])
                # index of corresponding control trace
                control_i = (row_i[0], row_i[1], 0.0)
                # get voltage value from control trace
                control_val = df.loc[control_i].data_v[time_i]
                # get voltage value for current trace
                val = row['data_v'][time_i]
                # set value in 'polarization column'
                df.at[row_i, 'polarization']=val-control_val

        df = df.reset_index()

        return df

    def _get_w_clopath(self, vtrace_df, clopath_param, split_columns=[], rerun=False,**kwargs):
        ''' run clopath algorithm to get synaptic weight evolution for entire df

        Arguments
        ----------
        ~vtrace_df: df must contain a column called input_times (each entry is a list of synaptic input times in ms) and data_v (each entry is a 1d array of voltage data)
        ~clopath_param: clopath rule parameters

        Returns
        ----------
        ~df: new data frame with columns:
                ~w_clopath: timeseriesof synaptic weights
                ~w_final: final value for w_clopath at end of simulation
                ~w_initial: initial value
                ~dw_clopath: final/initial
        '''
        # print progress to terminal
        print 'applying df function:', inspect.stack()[0][3]
        def _get_w(df, clopath_param):
            '''
            '''
            # create voltage array
            v = analysis._2array(df.data_v)
            # drop voltage data
            df = df.drop('data_v', axis=1)
            # fill input times with nan and get values
            input_times = df.input_times.fillna(value=np.nan).values
            # create time vector from dt and tstop
            dt = df.dt.unique()[0]
            tstop = df.tstop.unique()[0]
            # t = np.arange(0., tstop+dt, dt)
            print v.shape[1]
            # t = np.arange(0, dt*v.shape[1]-dt, dt)
            t = np.linspace(0, dt*v.shape[1], num=v.shape[1])
            print t.shape
            # create input array
            input_array = analysis._input_times_2array(input_times, t)
            print v.shape, input_array.shape
            # assert same shape for voltage and input arrays
            assert v.shape==input_array.shape, 'input array and voltage array must be same shape'
            # run clopath algorithm to get weights
            w_clopath = analysis.Clopath()._clopath(x=input_array, u=v, fs=1/dt, w0=0.5, homeostatic=False, param=clopath_param)
            # preallocate column for storing weight time series
            if 'w_clopath' not in df:
                df['w_clopath']=None
            # iterate over dataframe rows and add corresponding clopath data
            for i, row in df.iterrows():
                df.at[i, 'w_clopath'] = copy.copy(w_clopath[i,:])
                df.at[i, 'w_final'] = w_clopath[i,-1]
                df.at[i, 'w_initial']=w_clopath[i, 0]
                df.at[i, 'dw_clopath'] = w_clopath[i,-1]/w_clopath[i,0]
            return df
        
        # reset index and copy df
        df = vtrace_df.reset_index().copy()
        # if voltage traces are different length due to a set of parameters, split dataframe acccording to those parameters (so that voltage data can be made into an array, then recombine dataframes)
        if len(split_columns)>0:
            df = df.set_index(split_columns)
            split_indices = df.index.unique().values
            df_new = pd.DataFrame()
            for index in split_indices:
                print 'index:', index
                print 'df_new:', df_new.shape
                df_temp = df.loc[index].copy()
                df_temp.reset_index(inplace=True)
                df_temp = _get_w(df=df_temp, clopath_param=clopath_param)
                if df_new.empty:
                    print 'df_new empty'
                    df_new = df_temp.copy()
                else:
                    print 'appending'
                    df_new = df_new.append(df_temp, ignore_index=True)
                    print 'df_new:', df_new.shape

        else:
            df_new = _get_w(df=df, clopath_param=clopath_param)

        return df_new

    def _get_spikes(self, vtrace_df, threshold=-40, split_columns=[], **kwargs):
        ''' get spike times for each trace and store in df

        Arguments
        ----------
        ~vtrace_df: df must contain a column called data_v (each entry is a 1d array of voltage data)
        ~threshold: spike detection threshold (mV)
        ~split_columns: list of columns to group the data frame by.  if some voltage traces are different length due to a set of parameters, split dataframe accordingly so that voltage data can be made into arrays, then recombine dataframes after processing the individual arrays

        Returns
        ----------
        ~df: new data frame with columns:
                ~w_clopath: timeseriesof synaptic weights
                ~w_final: final value for w_clopath at end of simulation
                ~w_initial: initial value
                ~dw_clopath: final/initial
        '''
        # print progress to terminal
        print 'applying df function:', inspect.stack()[0][3]
        def _spikes(df, threshold):
            '''
            '''
            # reset dataframe index to range index
            df.reset_index(inplace=True)
            # create voltage array
            v = analysis._2array(df.data_v)
            # drop voltage data
            df = df.drop('data_v', axis=1)
            # create time vector from dt and tstop
            dt = df.dt.unique()[0]
            # run clopath algorithm to get weights
            spike_idx, spike_train, spike_times = analysis.Spikes()._get_spikes(data=v, threshold=threshold, dt=dt,)
            # preallocate column for storing weight time series
            if 'spike_idx' not in df:
                df['spike_idx']=None
            if 'spike_train' not in df:
                df['spike_train']=None
            if 'spike_times' not in df:
                df['spike_times']=None
            # iterate over dataframe rows and add corresponding spike data
            for i, row in df.iterrows():
                df.at[i, 'spike_idx'] = spike_idx[i]
                df.at[i, 'spike_train'] = spike_train[i,:]
                df.at[i, 'spike_times'] = spike_times[i]
            return df
        
        # reset index and copy df
        df = vtrace_df.reset_index().copy()
        # if voltage traces are different length due to a set of parameters, split dataframe acccording to those parameters (so that voltage data can be made into an array, then recombine dataframes)
        if len(split_columns)>0:
            df = df.set_index(split_columns)
            split_indices = df.index.unique().values
            df_new = pd.DataFrame()
            for index in split_indices:
                print 'index:', index
                print 'df_new:', df_new.shape
                df_temp = df.loc[index].copy()
                df_temp.reset_index(inplace=True)
                df_temp = _spikes(df=df_temp, threshold=threshold)
                if df_new.empty:
                    print 'df_new empty'
                    df_new = df_temp.copy()
                else:
                    print 'appending'
                    df_new = df_new.append(df_temp, ignore_index=True)
                    print 'df_new:', df_new.shape

        else:
            df_new = _spikes(df=df, threshold=threshold)

        return df_new
        
    def _get_polarization(self, df, time=10, rerun=False, file_limit=[], **kwargs):
        ''' get change in membrane potential in each recorded compartment in response to DCS
        
        Arguments
        ----
        -df : dataframe with olumn 'data_v' containing voltage data in each compartment, and 'dt' column
        -time : time point to take the difference between stimulated and corresponding control trace (ms)

        Return
        ------
        -df : dataframe with new 'polarization' column

        '''
        # print progress to terminal
        print 'applying df function:', inspect.stack()[0][3]
        # reorganize dataframe
        df = df.set_index(['trial_id','location','field'])
        # column name to add
        coladd = 'polarization'
        # if column doesn't exist or rerun is specified, process all rows
        #----------------------------------------------------------------
        if coladd not in df or rerun:
            # print progress
            print 'processing all rows'
            # iterate over rows
            cnt=-1
            for row_i, row in df.iterrows():
                cnt+=1
                # check for file limit
                if file_limit and cnt> file_limit:
                    print 'file limit reached:', file_limit
                    break
                # time point to measure polarization at
                time_i = int(time/row['dt'])
                # index of corresponding control trace
                control_i = (row_i[0], row_i[1], 0.0,)
                # print df.loc[control_i].data_v.drop_duplicates()
                # get voltage value from control trace
                control_val = df.loc[control_i].data_v.iloc[0][time_i]
                # get voltage value for current trace
                # print row['data_v']
                val = row['data_v'][time_i]
                # set value in 'polarization column'
                df.at[row_i, 'polarization']=val-control_val
        # if column already exists but some values haven't been processed
        #----------------------------------------------------------------
        # FIXME apply this functionality to other df functions!!!!!!!!!!!!!!!
        elif coladd in df and df[coladd].isnull().values.any():
            # get only unprocessed rows (index must match original df)
            df_temp = df[df[coladd].isnull()].copy()
            # print progress
            print 'unprocessed rows found:', df_temp.shape[0]
            # iterate over rows
            cnt=-1
            for row_i, row in df_temp.iterrows():
                cnt+=1
                # check for file limit
                if file_limit and cnt> file_limit:
                    print 'file limit reached:', file_limit
                    break
                # time point to measure polarization at
                time_i = int(time/row['dt'])
                # index of corresponding control trace
                control_i = (row_i[0], row_i[1], 0.0)
                # get voltage value from control trace
                control_val = df.loc[control_i].data_v[time_i]
                # get voltage value for current trace
                val = row['data_v'][time_i]
                # set value in 'polarization column'
                df.at[row_i, 'polarization']=val-control_val
        # reset index
        df = df.reset_index()

        return df

    def _get_xcorr_soma(self, spikes_df, rerun=False, file_limit=[], **kwargs):
        ''' for each recorded compartment get cross correlation of spikes wit hsoma
        
        Arguments
        ---------
        ~spikes_df : dataframe with column 'spike_train' containing boolean vector of spike times, and 'dt' column
        ~variable : column name in df

        Return
        ------
        ~spikes_df : dataframe with new column as 'xcorr_soma'

        Notes
        -----

        '''
        
        # print progress to terminal
        print 'applying df function:', inspect.stack()[0][3]
        # reorganize dataframe
        df = spikes_df.set_index(['trial_id','location','field', 'path'])
        # drop duplicated rows
        df = df[~df.index.duplicated()]
        dt = df.dt.unique()[0]
        # column name to add
        coladd = 'xcorr_soma'
        coladd_t ='xcorr_t'
        # df[coladd]=None
        # if column doesn't exist or rerun is specified, process all rows
        #----------------------------------------------------------------
        if coladd not in df or rerun:
            # print progress
            print 'processing all rows'
            if coladd not in df:
                df[coladd]=None
            if coladd_t not in df:
                df[coladd_t]=None
            # iterate over rows
            cnt=-1
            for row_i, row in df.iterrows():
                print row_i
                cnt+=1
                # check for file limit
                if file_limit and cnt> file_limit:
                    print 'file limit reached:', file_limit
                    break
                # index of corresponding soma trace
                soma_i = (row_i[0], ('soma',0,0), row_i[2], )
                # print soma_i
                # get spike_train from soma trace
                soma_train = df.loc[soma_i]['spike_train'].iloc[0]
                # print soma_train
                # get spike train for current location
                current_train = row['spike_train']
                # print soma_train.shape, current_train.shape
                # get cross correlation
                xcorr_soma = np.correlate(soma_train, current_train, 'full')
                # create time vector 
                xcorr_t = dt*(np.arange(0, len(xcorr_soma))-(len(xcorr_soma)-1)/2)
                # print xcorr_soma.shape
                # set value in 'polarization column'
                print df.loc[row_i].paths
                print df.at[row_i, coladd]
                # print df.at[row_i, coladd].paths
                df.at[row_i, coladd] = xcorr_soma
                df.at[row_i, coladd_t]=xcorr_t
        # if column already exists but some values haven't been processed
        #----------------------------------------------------------------
        # FIXME apply this functionality to other df functions!!!!!!!!!!!!!!!
        elif coladd in df and df[coladd].isnull().values.any():
            # get only unprocessed rows (index must match original df)
            df_temp = df[df[coladd].isnull()].copy()
            # print progress
            print 'unprocessed rows found:', df_temp.shape[0]
            # iterate over rows
            cnt=-1
            for row_i, row in df_temp.iterrows():
                cnt+=1
                # check for file limit
                if file_limit and cnt> file_limit:
                    print 'file limit reached:', file_limit
                    break
                # index of corresponding soma trace
                soma_i = (row_i[0], ('soma',0,0), row_i[2],)
                # print soma_i
                # get spike_train from soma trace
                soma_train = df.loc[soma_i]['spike_train'].iloc[0]
                # print soma_train
                # get spike train for current location
                current_train = row['spike_train']
                # print soma_train.shape, current_train.shape
                # get cross correlation
                xcorr_soma = np.correlate(soma_train, current_train, 'full')
                # create time vector 
                xcorr_t = dt*(np.arange(0, len(xcorr_soma))-(len(xcorr_soma)-1)/2)
                # print xcorr_soma.shape
                # set value in 'polarization column'
                df.at[row_i, coladd] = xcorr_soma
                df.at[row_i, coladd_t]=xcorr_t
        # reset index
        df = df.reset_index()

        return df

    def _get_dcs_effect(self, df, variables, rerun=False, **kwargs):
        ''' get effect of DCS on a given variable at given location
        
        Arguments
        ---------
        ~df : dataframe with olumn 'data_v' containing voltage data in each compartment, and 'dt' column
        ~variable : column name in df

        Return
        ------
        ~df : dataframe with new column as 'original_column_name'+'_dcs_diff'

        Notes
        -----

        '''
        
        # reorganize dataframe
        df = df.set_index(['trial_id','location','field'])
        for variable in variables:
            coladd = variable+'_dcs_diff'
            if coladd not in df or rerun:
                # print progress to terminal
                #-----------------------------
                print 'applying df function:', inspect.stack()[0][3], 'to column:', variable
                # preallocate column with 'object' dtype
                df[coladd]=None
                # iterate over rows
                for row_i, row in df.iterrows():
                    # index of corresponding control trace
                    control_i = (row_i[0], row_i[1], 0.0)
                    # get voltage value from control trace
                    control_val = df.loc[control_i][variable]
                    # get value for current trace
                    val = row[variable]
                    # set value in 'polarization column'
                    df.at[row_i, coladd]=val/control_val

        df = df.reset_index()

        return df

    def _apply_filters(self, df, variables, filters, rerun=False, file_limit=[],**kwargs):
        '''
        Arguments
        ---------
        ~df:
        ~variables:
        ~filters:dictionary of filters [filter_name][b coefficients, a coefficients]
        '''
        # reorganize dataframe
        df = df.set_index(['trial_id','location','field'])
        for variable in variables:
            for filt_key, filt in filters.iteritems():
                coladd = variable+'_filt_'+filt_key
                coladd_hilb = coladd+'_hilb'
                if coladd not in df or rerun:
                    # print progress to terminal
                    #-----------------------------
                    print 'applying filter ', filt_key, ' to ', variable
                    # preallocate column with 'object' dtype
                    df[coladd]=None
                    df[coladd_hilb]=None
                    cnt=0
                    # iterate over rows
                    for row_i, row in df.iterrows():
                        cnt+=1
                        if file_limit and cnt>file_limit:
                            continue
                        # print row_i
                        raw_data = row[variable]
                        # filtered data
                        # filtered_data = signal.filtfilt(filt[0], filt[1], raw_data, method='pad')
                        filtered_data = signal.lfilter(filt[0], filt[1], raw_data,)
                        # 
                        filtered_data_hilb = np.abs(signal.hilbert(filtered_data - np.mean(filtered_data)))
                        df.at[row_i, coladd]=filtered_data
                        df.at[row_i, coladd_hilb]=filtered_data_hilb


        df = df.reset_index()

        return df

    def _apply_filters_cascade(self, df, variables, filters, rerun=False, file_limit=[], **kwargs):
        '''
        Arguments
        ---------
        ~df:
        ~variables:
        ~filters:dictionary of filters in series [filter_name][subfilter_name][b coefficients, a coefficients]
        '''
        # print df.index.names
        # reorganize dataframe
        df = df.set_index(['trial_id','location','field'])
        # iterate over variables
        for variable in variables:
            # iterate over filters
            for filt_key, filt in filters.iteritems():
                # column name
                coladd = variable+'_filt_'+filt_key
                coladd_hilb = coladd+'_hilb'
                # run column?
                if coladd not in df or rerun:
                    # print progress to terminal
                    #-----------------------------
                    print 'applying filter ', filt_key, ' to ', variable
                    # preallocate column with 'object' dtype
                    df[coladd]=None    
                    df[coladd_hilb]=None   
                    cnt=0                 
                    # iterate over rows
                    for row_i, row in df.iterrows():
                        cnt+=1
                        if file_limit and cnt>file_limit:
                            continue
                        # get raw data
                        raw_data = row[variable]
                        # copy to be iteratively filtered
                        filtered_data=copy.deepcopy(raw_data)
                        # iterate over subfilters
                        for subfilt_key, subfilt in filt.iteritems():
                            # print progress to terminal
                            #-----------------------------
                            print 'applying subfilter ', subfilt_key, ' to ', variable
                            # update filtered data
                            # filtered_data = signal.filtfilt(subfilt[0], subfilt[1], filtered_data, method='pad')
                            filtered_data = signal.lfilter(subfilt[0], subfilt[1], filtered_data, )
                        filtered_data_hilb = np.abs(signal.hilbert(filtered_data - np.mean(filtered_data)))
                        # add to dataframe
                        df.at[row_i, coladd]=filtered_data
                        df.at[row_i, coladd_hilb]=filtered_data_hilb

        df = df.reset_index()

        return df