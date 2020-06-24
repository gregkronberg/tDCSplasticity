

import functions
from functions import Clopath
import numpy as np
import pandas as pd
from scipy.optimize import rosen, differential_evolution
from multiprocessing import Pool
import multiprocessing
import uuid
import datetime
import os
import pickle

df_funcs = functions.DfFuncs()

# conditions = ['apical_20hz_anodal', 'apical_20hz_cathodal','apical_20hz_control','apical_tbs_anodal', 'apical_tbs_cathodal','apical_tbs_control', 'basal_20hz_anodal', 'basal_20hz_cathodal','basal_20hz_control','basal_tbs_anodal', 'basal_tbs_cathodal','basal_tbs_control', 'apical_1hz_anodal', 'apical_1hz_cathodal','apical_1hz_control',]
# data = [1.3, 1.15, 1.15, 
#         1.5, 1.25, 1.3,
#         1.8, 1.4, 1.5,
#         1.8, 1.5, 1.55, 
#         0.9, 0.9, 0.8]

def get_experimental_data():
    experimental_data={'apical_20hz':[1.15, 1.15, 1.3],
        'basal_20hz':[1.4, 1.5, 1.8],
        'apical_tbs':[1.25, 1.3, 1.5],
        'basal_tbs':[1.5, 1.55, 1.8],
        'apical_1hz':[.9, .8,.9]
        }
    return experimental_data

def get_data_directories():
    # data_directories = {'apical_tbs':'Data/exp_reduced_neuron_tbs_quick_test/',}
    data_directories = {
    'apical_20hz':'Data/exp_reduced_neuron_20hz/',
    # 'basal_20hz':'Data/exp_reduced_neuron_20hz_basal/',
    'apical_tbs':'Data/exp_reduced_neuron_tbs/',
    # 'basal_tbs':'Data/exp_reduced_neuron_tbs_basal/',
    'apical_1hz':'Data/exp_reduced_neuron_1hz/'
    }
    return data_directories

# load vtrace_dfs
#------------------------------------------
def get_v_df_all(data_directories, ):
    v_df_all = {}
    filename = 'vtrace_df'
    for exp, directory in data_directories.iteritems(): 
        v_df_temp = functions._load_group_data(directory=directory, filename=filename)
        v_df_temp = functions._set_index(v_df_temp, ['path_1_w_mean'])
        v_df_temp = v_df_temp.loc[.006]
        v_df_all[exp] = v_df_temp
    return v_df_all
    # v_df_all[exp] = functions._load_group_data(directory=directory, filename=filename)
# print 'done loading'
# initial guess
#--------------------------------------------
# clopath_param_init = param.ParamClopath().kronberg_2020_reduced()

# clopath_param_template= {
#             'clopath_A_m0':100E-5, # depression magnitude parameter (mV^-1)
#             'clopath_A_p':40E-5, # amplitude for potentiation (mV^-2)
#             'clopath_tetam':-70,#-41, # depression threshold (mV)
#             'clopath_tetap':-61,#-38, # potentiation threshold (mV)
#             'clopath_tau_x':8,#-38, # time constant for presynaptic filter (ms)
#             'clopath_tau_m':20,#-38, # time constant for depression lowpass filter
#             'clopath_tau_p': 3, # time constant (ms) for low pass filter post membrane potential for potentiation
            
#             'clopath_delay':0, # conduction delay (ms)
#             'clopath_LTD_delay':1, # conduction delay for LTD  
#             'clopath_lower_bound':0.0,
#             'clopath_upper_bound':None,#3.0,
#             }

# clopath_param_init = np.array([100E-5, 40E-5, -70, -61, 8, 20, 3])
def get_param_bounds():
    # clopath_param_bounds = [(1E-10, 1E-4), (1E-10, 1E-4), (-73, -60), (-73, -55), (4, 20), (10, 50), (1, 40)]
    clopath_param_bounds = [(1E-6, 1E-2), (1E-6, 1E-2),(-67, -55)]
    return clopath_param_bounds
# function to minimize
#-----------------------------------------------
def f(clopath_param, *args):
    '''
    '''
    # # print 'running'
    # clopath_param_dict= {
    #         'A_m0':clopath_param[0], # depression magnitude parameter (mV^-1)
    #         'A_p':clopath_param[1], # amplitude for potentiation (mV^-2)
    #         'tetam':clopath_param[2],#-70,#-41, # depression threshold (mV)
    #         'tetap':clopath_param[3],#-61,#-38, # potentiation threshold (mV)
    #         'tau_x':clopath_param[4],#8,#-38, # time constant for presynaptic filter (ms)
    #         'tau_m':clopath_param[5],#20,#-38, # time constant for depression lowpass filter
    #         'tau_p':clopath_param[6],# 3, # time constant (ms) for low pass filter post membrane potential for potentiation
            
    #         'delay':0, # conduction delay (ms)
    #         'LTD_delay':1, # conduction delay for LTD  
    #         'lower_bound':0.0,
    #         'upper_bound':None,#3.0,
    #         }

    clopath_param_dict = {
            'clopath_A_m0':clopath_param[0], # depression magnitude parameter (mV^-1)
            'clopath_A_p':clopath_param[1], # amplitude for potentiation (mV^-2)
            'clopath_tetam':-72,#-41, # depression threshold (mV)
            'clopath_tetap':clopath_param[2],#-38, # potentiation threshold (mV)
            'clopath_tau_x':8,#-38, # time constant for presynaptic filter (ms)
            'clopath_tau_m':20,#-38, # time constant for depression lowpass filter
            'clopath_tau_p': 3, # time constant (ms) for low pass filter post membrane potential for potentiation
            
            'clopath_delay':0, # conduction delay (ms)
            'clopath_LTD_delay':1, # conduction delay for LTD  
            'clopath_lower_bound':0.0,
            'clopath_upper_bound':None,#3.0,
            }
    v_df_all=args[0]
    experimental_data=args[1]
    w_df ={}
    output ={}
    error={}
    total_error=0.
    for exp, v_df in v_df_all.iteritems():
        print exp
        if 'apical' in exp:
            tree_key = 'apic'
        elif 'basal' in exp:
            tree_key = 'dend'

        w_df[exp] = df_funcs._get_w_clopath(vtrace_df=v_df, clopath_param=clopath_param_dict, input_times_key='data_input_times')
        w_df[exp] = functions._set_index(w_df[exp], ['field', 'tree_key'])
        # print 'test', w_df[exp].loc[(-5, tree_key),:].dw_clopath
        output[exp] = [w_df[exp].loc[(-5, tree_key),:].dw_clopath.mean(), 
                w_df[exp].loc[(0, tree_key),:].dw_clopath.mean(),
                w_df[exp].loc[(5, tree_key),:].dw_clopath.mean()]
        error[exp] = np.sum((np.array(experimental_data[exp])-np.array(output[exp]))**2)
        print 'output',output[exp]
        total_error += error[exp]
    # total_error = np.sum([_v for _k,_v in error.iteritems()])
    # total_error = error.values()[0]
    print 'total_error',total_error
    # print clopath_param_dict
    # print error
    # total_error = np.abs(np.sum(clopath_param))
    return total_error

# print 'before diff evo call'
# def run_opt():
#     result = differential_evolution(f, bounds=clopath_param_bounds, args=(v_df_all, ), maxiter=2, updating='deferred', popsize=1, workers=2, disp=True)

def _generate_trial_id( ):
    '''
    '''
    # create unique identifier for each trial
    uid = str(uuid.uuid1().int)[-5:]
    now = datetime.datetime.now()
    trial_id = '-'.join(['{:04d}'.format(now.year), '{:02d}'.format(now.month), '{:02d}'.format(now.day), '{:02d}'.format(now.hour), '{:02d}'.format(now.minute), '{:02d}'.format(now.second), '{:02d}'.format(now.microsecond), uid])
    return trial_id

def _save_data(data, file_name=None, data_directory=None, **kwargs): 
    '''
    '''
    if file_name is None:
        file_name = str(_generate_trial_id())
    if data_directory is None:
        data_directory='parameters/'
    # check if folder exists with experiment name
    if os.path.isdir(data_directory) is False:
        print 'making new directory to save data'
        os.mkdir(data_directory)

    # save data as pickle file
    with open(data_directory+file_name+'.pkl', 'wb') as output:
        
        print 'saving data'
        pickle.dump(data, output,protocol=pickle.HIGHEST_PROTOCOL)

if __name__=='__main__':
    multiprocessing.freeze_support()
    experimental_data = get_experimental_data()
    data_directories = get_data_directories()
    v_df_all = get_v_df_all(data_directories=data_directories)
    print 'keys', v_df_all.keys()
    clopath_param_bounds = get_param_bounds()
    n_workers=len(clopath_param_bounds)
    # p = Pool(n_workers)
    result = differential_evolution(f, bounds=clopath_param_bounds, args=(v_df_all, experimental_data), maxiter=2, updating='deferred', popsize=1, workers=n_workers, disp=True,tol=0.1, polish=False)
    result_dict = {
            'A_m0':result.x[0], # depression magnitude parameter (mV^-1)
            'A_p':result.x[1], # amplitude for potentiation (mV^-2)
            'clopath_tetam':-72,#-41, # depression threshold (mV)
            'clopath_tetap':result.x[2],#-38, # potentiation threshold (mV)
            'clopath_tau_x':8,#-38, # time constant for presynaptic filter (ms)
            'clopath_tau_m':20,#-38, # time constant for depression lowpass filter
            'clopath_tau_p': 3, # time constant (ms) for low pass filter post membrane potential for potentiation
            
            'clopath_delay':0, # conduction delay (ms)
            'clopath_LTD_delay':1, # conduction delay for LTD  
            'clopath_lower_bound':0.0,
            'clopath_upper_bound':None,#3.0,
            }
    trial_id = _generate_trial_id()
    filename = 'clopath_param_'+trial_id
    _save_data(data=result_dict, file_name=filename, data_directory='parameters/')
    
    print 'result', result

# if __name__ =='__main__':
#     bounds = [(0,2), (0, 2), (0, 2), (0, 2), (0, 2)]
#     result = differential_evolution(rosen, bounds, updating='deferred',workers=2, disp=True)
    # print result