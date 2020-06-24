"""
analysis

"""
import figsetup_exp_branch_sequence as figsetup
import numpy as np
import analysis
import glob
import pickle
import copy
import matplotlib.pyplot as plt
import itertools
from matplotlib import cm as colormap
from scipy import stats
import os
import functions as fncs

data_local = False
experiment='_'.join(__name__.split('_')[1:]).split('.')[0]
# kwargs = {'experiment':'_'.join(__name__.split('_')[1:]).split('.')[0]}
dpi=350

#############################################################################
# directories and filenames
#############################################################################
# directories
data_directory = 'Data/'
if data_local:
    data_directory='Data local/'
directory = data_directory+experiment+'/'
figure_directory = directory+'Figures/'
# make figure directory if doesnt exist
if os.path.isdir(figure_directory) is False:
    os.mkdir(figure_directory)

# variables to load
variables = ['vtrace_df',]
variables_reload=['vtrace_df']

arrayfuncs = fncs.ArrayFunctions()
dffuncs = fncs.DfFuncs()

#############################################################################
# load variables
#############################################################################
_g = globals()
def _load_variables(variables, global_variables=globals(), variables_reload=[], extension='.pkl'):
    '''
    '''
    for variable in variables:
        # check if exists
        exists = variable in globals()
        if (variable in variables_reload) or not exists :
            filename = variable
            _g[variable]=fncs._load_group_data(directory=directory, filename=filename)
_load_variables(variables=variables, global_variables=_g, variables_reload=variables_reload)

#####################################################################
# average voltage trace, weak path dendrites in paired vs unpaired
#####################################################################
def _vtrace_mean(vtrace_df=vtrace_df, figsetup=figsetup):
    # plot average voltage trace at soma and dendrite
    #------------------------------------------------
    figtype='vtrace_soma_mean_'
    figdf=figsetup.BuildFigDF()._trace_mean()

    # create column for sec_idx string to sort by which section was stimulated
    vtrace_df = dffuncs._to_string(vtrace_df, colnames=['path_1_sec_idx'])

    # truncate arrays to have same length for same sequence delay
    vtrace_df = vtrace_df.set_index(['Cm', 'path_1_syn_limit','path_1_nmda_ampa_ratio','branch_seg_L','path_1_delay', 'path_1_w_mean', 'path_1_sequence_direction', 'field', 'tree_key', 'path_1_sec_idx_str', 'seg_num'])
    # # vtrace_df = dffuncs._truncate_arrays(df=vtrace_df,)

    # vtrace_df['data_v'] = dffuncs._truncate_arrays(vtrace_df['data_v'])
    # figdf['fig_ymin']=-4
    # figdf['fig_dyticks']=1
    # figdf['fig_xmin']=0
    # figdf['fig_xmax']=60
    # vtrace_df['path'] = vtrace_df['path'].fillna('2')
    # vtrace_df = vtrace_df.set_index(['path_1_syn_num','field', 'tree', ])
    kwargs={'dt':1./40}
    array_funcs = [arrayfuncs._truncate_array]#[arrayfuncs._subtract_timeseries]
    array_funcs=[]
    array_func_kws = [{}]#[{'islice':slice(399, 400), 'axis':1}]
    figs, ax = fncs.PlotFuncs()._trace_mean(df=vtrace_df, figdf=figdf, variables='data_v',array_funcs=array_funcs, array_func_kws=array_func_kws, figformat='none', **kwargs)   
    vtrace_df = vtrace_df.reset_index()
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi) 
    return figs, ax, vtrace_df, figdf

figs, ax, vtrace_df, figdf= _vtrace_mean()