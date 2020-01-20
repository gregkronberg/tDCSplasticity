"""
analysis

"""
import numpy as np
import analysis
import vargen_functions
import glob
import pickle
import copy
import matplotlib.pyplot as plt
import itertools
from matplotlib import cm as colormap
from scipy import stats
import os
import figsetup_exp_1a_apical as figsetup

'''
'''
kwargs = {'experiment':'exp_1a_apical'}
dpi=350

#############################################################################
# directories and filenames
#############################################################################
# directories
directory = 'Data/'+kwargs['experiment']+'/'
figure_directory = 'Data/'+kwargs['experiment']+'/Figures/'
# make figure directory if doesnt exist
if os.path.isdir(figure_directory) is False:
    os.mkdir(figure_directory)
# filenames    
filename_vtrace = 'vtrace_df.pkl'
filename_w_clopath = 'w_clopath_df.pkl'
# filename_w_clopath_stdp = 'w_clopath_df_stdp.pkl'
filename_spikes = 'spikes_df.pkl'
# functions

arrayfuncs = analysis.ArrayFunctions()

#############################################################################
# load variables
#############################################################################
vtrace_df = analysis._load_group_data(directory=directory, filename=filename_vtrace)
spikes_df = analysis._load_group_data(directory=directory, filename=filename_spikes)
w_clopath_df = analysis._load_group_data(directory=directory, filename=filename_w_clopath)
# w_clopath_stdp_df = analysis._load_group_data(directory=directory, filename=filename_w_clopath_stdp)


#############################################################################
# dw final bar plot
#############################################################################
run_bar=False
if run_bar:
    # setup figdf
    #------------------------------------------------
    figtype='w_bar_'
    figdf=figsetup.BuildFigDF()._w_bar()
    figdf['fig_ymin']=1.4
    figdf['fig_dyticks']=.1
    figdf['fig_ytick_decimals']=1
    # figdf['fig_xmin']=0
    # figdf['fig_xmax']=60
    
    kwargs={'dt':1./40}
    # keyword arguments for array functions
    #--------------------------------------
    fs=40000.
    # filter coefficients
    filt = analysis.Filters(fs=fs).filters_series['iir_band_5_50']
    filt_kws = {'filt_series':filt,'hilbert':False,'islice':slice(400,None)}
    # subtract initial value
    subtract_kws = {'islice':slice(399, 400), 'axis':1}
    # get area indices
    area_kws = {'islice':slice(400, 2000)}
    # array functions 
    #------------------------------------
    # subtract, filter, get area
    array_funcs = []#[arrayfuncs._subtract_timeseries, arrayfuncs._filter_mean, arrayfuncs._get_area]
    array_func_kws = []#[subtract_kws, filt_kws, area_kws]
    # set index in data df
    #--------------------------
    w_clopath_df = w_clopath_df.set_index(['field', 'path_1_syn_num', 'tree'])
    # plots
    #-------
    figs, ax = analysis.PlotFuncs()._bar(df=w_clopath_df, figdf=figdf, variable='dw_clopath',array_funcs=array_funcs, array_func_kws=array_func_kws,**kwargs)
    # reset index   
    vtrace_df = vtrace_df.reset_index()
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi) 

#############################################################################
# epsp area bar plot
#############################################################################
run_bar=False
if run_bar:
    # setup figdf
    #------------------------------------------------
    figtype='epsp_bar_'
    figdf=figsetup.BuildFigDF()._epsp_area()
    figdf['fig_yscale']=1./40.
    figdf['fig_ytick_decimals']=0
    figdf['fig_dyticks']=20
    figdf['fig_ymin']=100
    # figdf['fig_xmax']=60
    vtrace_df = vtrace_df.set_index(['field', 'path_1_syn_num', 'tree'])
    kwargs={'dt':1./40}
    # keyword arguments for array functions
    #--------------------------------------
    fs=40000.
    # filter coefficients
    filt = analysis.Filters(fs=fs).filters_series['iir_band_5_50']
    filt_kws = {'filt_series':filt,'hilbert':False,'islice':slice(400,None)}
    # subtract initial value
    subtract_kws = {'islice':slice(399, 400), 'axis':1}
    # get area indices
    area_kws = {'islice':slice(400, 2000)}
    # array functions 
    #------------------------------------
    # subtract, filter, get area
    array_funcs = [arrayfuncs._subtract_timeseries, arrayfuncs._filter_mean, arrayfuncs._get_area]
    array_func_kws = [subtract_kws, filt_kws, area_kws]
    # plots
    #-------
    figs, ax = analysis.PlotFuncs()._bar(df=vtrace_df, figdf=figdf, variable='data_v',array_funcs=array_funcs, array_func_kws=array_func_kws,**kwargs)
    # reset index   
    vtrace_df = vtrace_df.reset_index()
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi) 

#############################################################################
# hilbert area bar plot
#############################################################################
run_bar=False
if run_bar:
    # setup figdf
    #------------------------------------------------
    figtype='hilbert_bar_'
    figdf=figsetup.BuildFigDF()._hilbert_bar()
    figdf['fig_ymin']=0
    figdf['fig_dyticks']=5
    # figdf['fig_xmin']=0
    # figdf['fig_xmax']=60
    vtrace_df = vtrace_df.set_index(['field', 'path_1_syn_num', 'tree'])
    kwargs={'dt':1./40}
    # keyword arguments for array functions
    #--------------------------------------
    fs=40000.
    # filter coefficients
    filt = analysis.Filters(fs=fs).filters['iir_high_300']
    filt_kws = {'filt':filt,'hilbert':True,'islice':slice(None)}
    # subtract initial value
    subtract_kws = {'islice':slice(399, 400), 'axis':1}
    # get area indices
    area_kws = {'islice':slice(400, 2000), 'axis':0, 'norm':False}
    # array functions 
    #------------------------------------
    # subtract, filter, get area
    array_funcs = [arrayfuncs._subtract_timeseries, arrayfuncs._filter_mean, arrayfuncs._get_area]
    array_func_kws = [subtract_kws, filt_kws, area_kws]
    # plots
    #-------
    figs, ax = analysis.PlotFuncs()._bar(df=vtrace_df, figdf=figdf, variable='data_v',array_funcs=array_funcs, array_func_kws=array_func_kws,**kwargs)
    # reset index   
    vtrace_df = vtrace_df.reset_index()
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi) 

#####################################################################
# dose response
#####################################################################
run_weights=False
if run_weights:
    # dose response figures
    #-------------------------------------------------------------------
    figtype='dose_response_'
    figdf = figsetup.BuildFigDF()._dose_response()
    figs, ax = analysis.PlotFuncs()._dose_response(df=w_clopath_df, figdf=figdf, variable='dw_clopath')
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi)

    # dose response figures
    #-------------------------------------------------------------------
    # figtype='dose_response_stdp_param_'
    # figdf = analysis.BuildFigDF()._dose_response()
    # figs, ax = analysis.PlotFuncs()._dose_response(df=w_clopath_stdp_df, figdf=figdf, variable='dw_clopath')
    # # save figure
    # #------------
    # for fig_key, fig in figs.iteritems():
    #     fname = figure_directory+figtype+str(fig_key)+'.png'
    #     fig.savefig(fname, format='png', dpi=dpi)

#####################################################################
# vtrace raw data
#####################################################################
run_vtrace=False
if run_vtrace:
    # plot average voltage trace at soma and dendrite
    #------------------------------------------------
    figtype='vtrace_'
    figdf=figsetup.BuildFigDF()._trace_mean()
    figdf['fig_ymin']=-4
    figdf['fig_dyticks']=1
    figdf['fig_xmin']=0
    figdf['fig_xmax']=60
    vtrace_df = vtrace_df.set_index(['field', 'path_1_syn_num', 'tree'])
    kwargs={'dt':1./40}
    array_funcs = [arrayfuncs._subtract_timeseries]
    array_func_kws = [{'islice':slice(399, 400), 'axis':1}]
    figs, ax = analysis.PlotFuncs()._trace_mean(df=vtrace_df, figdf=figdf, variable='data_v',array_funcs=array_funcs, array_func_kws=array_func_kws,**kwargs)   
    vtrace_df = vtrace_df.reset_index()
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi) 

#####################################################################
# wtrace raw data
#####################################################################
run_vtrace=True
if run_vtrace:
    # plot average voltage trace at soma and dendrite
    #------------------------------------------------
    figtype='wtrace_'
    figdf=figsetup.BuildFigDF()._wtrace_mean()
    figdf['fig_ymin']=1
    figdf['fig_dyticks']=.2
    figdf['fig_dxticks']=10
    figdf['fig_xmin']=0
    figdf['fig_xmax']=60
    w_clopath_df = w_clopath_df.set_index(['field', 'path_1_syn_num', 'tree'])
    kwargs={'dt':1./40}
    array_funcs = [arrayfuncs._normalize]#[arrayfuncs._subtract_timeseries]
    array_func_kws = [{}]#[{'islice':slice(399, 400), 'axis':1}]
    figs, ax = analysis.PlotFuncs()._trace_mean(df=w_clopath_df, figdf=figdf, variable='w_clopath',array_funcs=array_funcs, array_func_kws=array_func_kws,**kwargs)   
    w_clopath_df = w_clopath_df.reset_index()
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi) 

#####################################################################
# vtrace raw data individual
#####################################################################
run_vtrace=False
if run_vtrace:
    # plot average voltage trace at soma and dendrite
    #------------------------------------------------
    figtype='vtrace_individual_'
    trial_id = vtrace_df.trial_id.iloc[1100]
    figtype='vtrace_individual_'+trial_id+'_'
    figdf=figsetup.BuildFigDF()._trace_individual(trial_id=trial_id)
    figdf['fig_ymax']=-40
    # figdf['fig_dyticks']=1
    figdf['fig_xmin']=0
    figdf['fig_xmax']=60
    vtrace_df = vtrace_df.set_index(['field', 'path_1_syn_num', 'tree', 'trial_id'])
    kwargs={'dt':1./40}
    array_funcs = []#[arrayfuncs._subtract_timeseries]
    array_func_kws = []#[{'islice':slice(399, 400), 'axis':1}]
    figs, ax = analysis.PlotFuncs()._trace_mean(df=vtrace_df, figdf=figdf, variable='data_v',array_funcs=array_funcs, array_func_kws=array_func_kws,**kwargs)   
    vtrace_df = vtrace_df.reset_index()
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi) 

#####################################################################
# vtrace lowpass
#####################################################################
run_vtrace=False
if run_vtrace:
    # plot average voltage trace at soma and dendrite
    #------------------------------------------------
    figtype='vtrace_lowpass_'
    figdf=figsetup.BuildFigDF()._trace_mean()
    # update figdf
    #---------------
    figdf['fig_ymin']=-4
    figdf['fig_dyticks']=1
    figdf['fig_xmin']=0
    figdf['fig_xmax']=60
    vtrace_df = vtrace_df.set_index(['field', 'path_1_syn_num', 'tree'])
    # array functions to pass to plot
    #------------------------------------
    # lowpass filter
    fs=40000.
    filt = analysis.Filters(fs=fs).filters_series['iir_band_5_50']
    filt_kws = {'filt_series':filt, 'hilbert':False, 'islice':slice(400,None)}
    # normalize to onset value
    subtract_kws = {'islice':slice(399, 400), 'axis':1}
    subtract_kws_post = {'islice':slice(399, 400), 'axis':1}
    # array functions and kws
    array_funcs = [arrayfuncs._subtract_timeseries, arrayfuncs._filter_mean, arrayfuncs._subtract_timeseries]
    array_func_kws = [subtract_kws, filt_kws, subtract_kws_post]
    kwargs={'dt':1./40}
    # create plots
    #-------------
    figs, ax = analysis.PlotFuncs()._trace_mean(df=vtrace_df, figdf=figdf, variable='data_v', array_funcs=array_funcs, array_func_kws=array_func_kws, **kwargs)
    # reset df index
    #----------------   
    vtrace_df = vtrace_df.reset_index() 
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi) 

#####################################################################
# vtrace hilbert
#####################################################################
run_vtrace=False
if run_vtrace:
    # plot average voltage trace at soma and dendrite
    #------------------------------------------------
    figtype='vtrace_hilbert_'
    figdf=figsetup.BuildFigDF()._trace_mean()
    # update figdf
    #---------------
    figdf['fig_ymin']=0
    figdf['fig_dyticks']=1
    figdf['fig_xmin']=0
    figdf['fig_xmax']=60
    vtrace_df = vtrace_df.set_index(['field', 'path_1_syn_num', 'tree'])
    # array functions to pass to plot
    #------------------------------------
    # lowpass filter
    fs=40000.
    filt = analysis.Filters(fs=fs).filters['iir_high_300']
    filt_kws = {'filt':filt, 'hilbert':True, 'islice':slice(400,None)}
    # normalize to onset value
    subtract_kws = {'islice':slice(399, 400), 'axis':1}
    subtract_kws_post = {'islice':slice(399, 400), 'axis':0}
    # array functions and kws
    array_funcs = [arrayfuncs._subtract_timeseries, arrayfuncs._filter_mean, ]
    array_func_kws = [subtract_kws, filt_kws, subtract_kws_post]
    kwargs={'dt':1./40}
    # create plots
    #-------------
    figs, ax = analysis.PlotFuncs()._trace_mean(df=vtrace_df, figdf=figdf, variable='data_v', array_funcs=array_funcs, array_func_kws=array_func_kws, **kwargs)
    # reset df index
    #----------------   
    vtrace_df = vtrace_df.reset_index() 
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi) 

#####################################################################
# vtrace high pass
#####################################################################
run_vtrace=False
if run_vtrace:
    # plot average voltage trace at soma and dendrite
    #------------------------------------------------
    filt = 'filt_iir_high_300_hilb'
    variable='data_v_'+filt
    figtype='vtrace_'+filt+'_'
    figdf=figsetup.BuildFigDF()._trace_mean_filtered()
    vtrace_df = vtrace_df.set_index(['field', 'path_1_syn_num', 'tree'])
    kwargs={'dt':1./40}
    array_funcs = []#[arrayfuncs._subtract_timeseries]
    array_func_kws = [{}]#[{'islice':slice(399, 400), 'axis':1}]
    figs, ax = analysis.PlotFuncs()._trace_mean(df=vtrace_df_highpass, figdf=figdf, variable=variable, array_funcs=array_funcs, array_func_kws=array_func_kws,**kwargs)   
    vtrace_df = vtrace_df.reset_index()
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi) 