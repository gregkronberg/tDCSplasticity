"""
analysis

"""
import figsetup_exp_1a_basal_full_duration as figsetup
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

'''
'''
kwargs = {'experiment':'exp_1a_basal_full_duration'}
dpi=350

#############################################################################
# directories and filenames
#############################################################################
# directories
directory = 'Data local/'+kwargs['experiment']+'/'
figure_directory = 'Data local/'+kwargs['experiment']+'/Figures/'
# make figure directory if doesnt exist
if os.path.isdir(figure_directory) is False:
    os.mkdir(figure_directory)
# filenames    
filename_vtrace = 'vtrace_df.pkl'
filename_w_clopath = 'w_clopath_df.pkl'
filename_spikes = 'spikes_df.pkl'

arrayfuncs = analysis.ArrayFunctions()

#############################################################################
# load variables
#############################################################################
def _load_variables():
    vtrace_df = analysis._load_group_data(directory=directory, filename=filename_vtrace)
    spikes_df = analysis._load_group_data(directory=directory, filename=filename_spikes)
    w_clopath_df = analysis._load_group_data(directory=directory, filename=filename_w_clopath)

    return vtrace_df, w_clopath_df, spikes_df

vtrace_df, w_clopath_df, spikes_df = _load_variables()
#####################################################################
# ltp final bar plot
#####################################################################
def _ltp_final_bar(w_clopath_df=w_clopath_df, figsetup=figsetup):
    # dose response figures
    #-------------------------------------------------------------------
    figtype='ltp_final_'
    figdf = figsetup.BuildFigDF()._ltp_bar()
    w_clopath_df = w_clopath_df.set_index(['path_1_syn_num','field', 'tree', ])
    # figdf.reset_index()
    figs, ax = analysis.PlotFuncs()._bar(df=w_clopath_df, figdf=figdf, variable='dw_clopath')
    w_clopath_df=w_clopath_df.reset_index()
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi)

_ltp_final_bar()

#####################################################################
# cross correlation of dendritic spikes with soma
#####################################################################
def _xcorr_soma(spikes_df=spikes_df, figsetup=figsetup):
    # plot average voltage trace at soma and dendrite
    #------------------------------------------------
    figtype='xcorr_soma_'
    figdf=figsetup.BuildFigDF()._xcorr()
    # figdf['fig_ymin']=-4
    # figdf['fig_dyticks']=1
    # figdf['fig_xmin']=0
    # figdf['fig_xmax']=60
    # vtrace_df['path'] = vtrace_df['path'].fillna('2')
    spikes_df = spikes_df.set_index(['field', 'path', 'path_2_syn_num', 'n_paths', 'tree'])
    kwargs={'dt':1./40, 't':spikes_df.xcorr_t.values[0]}
    array_funcs = []#[arrayfuncs._subtract_timeseries]
    array_func_kws = []#[{'islice':slice(399, 400), 'axis':1}]
    xvals = arrayfuncs._get_trace_xvals
    xvals_kw={'dt':1./40, 'xtype':'xcorr'}
    figs, ax = analysis.PlotFuncs()._trace_mean(
        df=spikes_df, 
        figdf=figdf, 
        variables='xcorr_soma',
        xvals=xvals, 
        xvals_kw=xvals_kw,
        array_funcs=array_funcs, 
        array_func_kws=array_func_kws,
        **kwargs)   
    spikes_df = spikes_df.reset_index()
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi) 

# _xcorr_soma()

#####################################################################
# average voltage trace, weak path dendrites in paired vs unpaired
#####################################################################
def _vtrace_mean(vtrace_df=vtrace_df, figsetup=figsetup):
    # plot average voltage trace at soma and dendrite
    #------------------------------------------------
    figtype='vtrace_soma_mean_'
    figdf=figsetup.BuildFigDF()._trace_mean()
    # figdf['fig_ymin']=-4
    # figdf['fig_dyticks']=1
    # figdf['fig_xmin']=0
    # figdf['fig_xmax']=60
    # vtrace_df['path'] = vtrace_df['path'].fillna('2')
    vtrace_df = vtrace_df.set_index(['path_1_syn_num','field', 'tree', ])
    kwargs={'dt':1./40}
    array_funcs = []#[arrayfuncs._subtract_timeseries]
    array_func_kws = []#[{'islice':slice(399, 400), 'axis':1}]
    figs, ax = analysis.PlotFuncs()._trace_mean(df=vtrace_df, figdf=figdf, variables='data_v',array_funcs=array_funcs, array_func_kws=array_func_kws,**kwargs)   
    vtrace_df = vtrace_df.reset_index()
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi) 

# _vtrace_mean()

#####################################################################
# average weight trace derivative, weak path dendrites in paired vs unpaired
#####################################################################
def _weight_trace_derivative(w_clopath_df=w_clopath_df, figsetup=figsetup):
    # plot average voltage trace at soma and dendrite
    #------------------------------------------------
    figtype='wtrace_diff_dend_mean_'
    figdf=figsetup.BuildFigDF()._wtrace_mean()
    # figdf['fig_ymin']=-4
    # figdf['fig_dyticks']=1
    # figdf['fig_xmin']=0
    # figdf['fig_xmax']=60

    # w_clopath_df['path'] = w_clopath_df['path'].fillna('2')
    w_clopath_df = w_clopath_df.set_index(['path_1_syn_num','field', 'tree', ])
    kwargs={'dt':1./40}
    array_funcs = [np.diff]#[arrayfuncs._subtract_timeseries]
    array_func_kws = [{'axis':1}]#[{'islice':slice(399, 400), 'axis':1}]
    figs, ax = analysis.PlotFuncs()._trace_mean(df=w_clopath_df, figdf=figdf, variables='w_clopath',array_funcs=array_funcs, array_func_kws=array_func_kws, figformat='standard', **kwargs)   
    w_clopath_df = w_clopath_df.reset_index()
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi) 

# _weight_trace_derivative()

def _weight_trace(w_clopath_df=w_clopath_df, figsetup=figsetup):
    # plot average voltage trace at soma and dendrite
    #------------------------------------------------
    figtype='wtrace_dend_mean_'
    figdf=figsetup.BuildFigDF()._wtrace_mean()
    figdf['fig_xscale']=1./1000
    figdf['fig_dxticks']=1
    figdf['fig_ymin']=0
    figdf['fig_ymax']=20.1
    figdf['fig_dyticks']=4
    figdf['fig_ytick_decimals']=0
    # figdf['fig_xmin']=0
    # figdf['fig_xmax']=60

    # w_clopath_df['path'] = w_clopath_df['path'].fillna('2')
    w_clopath_df = w_clopath_df.set_index(['path_1_syn_num','field', 'tree', ])
    kwargs={'dt':1./40}
    array_funcs = []#[arrayfuncs._subtract_timeseries]
    array_func_kws = [{}]#[{'islice':slice(399, 400), 'axis':1}]
    figs, ax = analysis.PlotFuncs()._trace_mean(df=w_clopath_df, figdf=figdf, variables='w_clopath',array_funcs=array_funcs, array_func_kws=array_func_kws, figformat='standard', **kwargs)   
    w_clopath_df = w_clopath_df.reset_index()
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi) 

# _weight_trace()

#####################################################################
# vtrace raw voltage data individual cells
#####################################################################
def _vtrace_individual(vtrace_df=vtrace_df, figsetup=figsetup):
    # plot average voltage trace at soma and dendrite
    #------------------------------------------------
    trial_id = vtrace_df[(vtrace_df.path_1_syn_num==8)].trial_id.unique()[0]
    figdf=figsetup.BuildFigDF()._trace_individual(trial_id=trial_id)
    # voltage trace for first 4 inputs
    #---------------------------------
    figtype='vtrace_individual_beginning_'+trial_id+'_'
    # figdf['fig_ymax']=-40
    # # figdf['fig_dyticks']=1
    figdf['fig_xmin']=0
    figdf['fig_xmax']=100
    # figdf['fig_xmin']=2810
    # figdf['fig_xmax']=3010
    figdf['fig_ymax']=-44.9
    figdf['fig_ymin']=-75.
    figdf['fig_dxticks']=50
    vtrace_df['path'] = vtrace_df['path'].fillna('2')
    vtrace_df = vtrace_df.set_index(['path_1_syn_num','field', 'tree', 'trial_id'])
    kwargs={'dt':1./40}
    array_funcs = []#[arrayfuncs._subtract_timeseries]
    array_func_kws = []#[{'islice':slice(399, 400), 'axis':1}]
    figs, ax = analysis.PlotFuncs()._trace_mean(df=vtrace_df, figdf=figdf, variables='data_v',array_funcs=array_funcs, array_func_kws=array_func_kws,**kwargs)   
    vtrace_df = vtrace_df.reset_index()
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi) 

    # voltage trace for last 4inputs
    #---------------------------------
    figtype='vtrace_individual_end_'+trial_id+'_'
    figdf['fig_xmin']=600
    figdf['fig_xmax']=700
    figdf['fig_ymax']=-44.9
    figdf['fig_ymin']=-75.
    figdf['fig_dxticks']=50
    figdf['fig_dyticks']=10
    vtrace_df['path'] = vtrace_df['path'].fillna('2')
    vtrace_df = vtrace_df.set_index(['path_1_syn_num','field', 'tree', 'trial_id'])
    kwargs={'dt':1./40}
    array_funcs = []#[arrayfuncs._subtract_timeseries]
    array_func_kws = []#[{'islice':slice(399, 400), 'axis':1}]
    figs, ax = analysis.PlotFuncs()._trace_mean(df=vtrace_df, figdf=figdf, variables='data_v',array_funcs=array_funcs, array_func_kws=array_func_kws,**kwargs)   
    vtrace_df = vtrace_df.reset_index()
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi) 

# _vtrace_individual()

#####################################################################
# spike time probability distribution
#####################################################################
def _spike_time_probability(spikes_df=spikes_df, figsetup=figsetup):
    # plot average voltage trace at soma and dendrite
    #------------------------------------------------
    figtype='spike_time_hist_dend_'
    figdf=figsetup.BuildFigDF()._spike_time_hist()
    # figdf['fig_ymin']=-4
    # figdf['fig_dyticks']=1
    # figdf['fig_xmin']=0
    # figdf['fig_xmax']=60
    # vtrace_df['path'] = vtrace_df['path'].fillna('2')
    spikes_df = spikes_df.set_index(['field', 'path', 'path_2_syn_num', 'n_paths', 'tree'])
    kwargs={'dt':1./40, 't':spikes_df.xcorr_t.values[0]}
    array_funcs = [arrayfuncs._flatten]#[arrayfuncs._subtract_timeseries]
    array_func_kws = [{}]#[{'islice':slice(399, 400), 'axis':1}]
    bins = arrayfuncs._get_bins# 'auto'#arrayfuncs._get_bins
    bins_kw={'binmin':20, 'binsize':1}
  
    figs, ax = analysis.PlotFuncs()._hist(
        df=spikes_df, 
        figdf=figdf, 
        variables='spike_times', 
        bins=bins, 
        bins_kw=bins_kw,
        array_funcs=array_funcs, 
        array_func_kws=array_func_kws,
        **kwargs) 
    spikes_df = spikes_df.reset_index()
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi, ) 

# _spike_time_probability()