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
import figsetup_exp_1a_dose_basal as figsetup

'''
'''
kwargs = {'experiment':'exp_1a_dose_basal'}
dpi=350
#############################################################################
# directories and filenames
#############################################################################
# directories
directory = 'Data/'+kwargs['experiment']+'/'
directory_basal = 'Data/exp_1a_dose_basal/'
figure_directory = 'Data/'+kwargs['experiment']+'/Figures/'
# make figure directory if doesnt exist
if os.path.isdir(figure_directory) is False:
    os.mkdir(figure_directory)
# filenames    
filename_vtrace = 'vtrace_df.pkl'
# filename_w_clopath = 'w_clopath_df.pkl'
filename_w_clopath = 'w_clopath_df.pkl'
filename_w_clopath2 = 'w_clopath2_df.pkl'
filename_w_clopath_stdp = 'w_clopath_stdp_df.pkl'
filename_w_clopath_stdp2 = 'w_clopath_stdp2_df.pkl'
filename_spikes = 'spikes_df.pkl'

arrayfuncs = analysis.ArrayFunctions()
buildfigdf = figsetup.BuildFigDF()


#############################################################################
# load variables
#############################################################################
def _load_variables():
    '''
    '''
    vtrace_df = analysis._load_group_data(directory=directory, filename=filename_vtrace)
    w_clopath_df = analysis._load_group_data(directory=directory, filename=filename_w_clopath)
    return vtrace_df, w_clopath_df

# vtrace_df, w_clopath_df = _load_variables()

#####################################################################
# field magnitude vs membrane polarization
#####################################################################
def _polarization_dose_response():
    # dose response figures
    #-------------------------------------------------------------------
    figtype='polarization_dose_response_'
    figdf = figsetup.BuildFigDF()._polarization_dose_response()
    figdf['fig_ytick_assert']=0
    figdf['fig_ytick_decimals']=0
    # figdf['fig_ymax']=0.021
    # figdf['fig_ymin']=1#0.021
    figdf['fig_xmax']=21
    figdf['fig_xmin']=-21
    figdf['fig_xtick_assert']=-20
    figdf['fig_dxticks']=5
    figdf['fig_dyticks']=1
    # figdf.reset_index()
    vtrace_df.set_index(['tree'], inplace=True)
    figs, ax = analysis.PlotFuncs()._var2var_mean(df=vtrace_df, figdf=figdf, variable='polarization', x_variable='field')
    figs, ax = analysis.PlotFuncs()._draw_lines(fig=figs, ax=ax, figdf=figdf)
    vtrace_df.reset_index( inplace=True)
    plt.show(block=False)
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi)

_polarization_dose_response()
#############################################################################
# dose response of mean dcs effect (dw_dcs)
#############################################################################
def _dose_response_dcseffect_mean(w_clopath_df, buildfigdf):

    # dose response figures
    #-------------------------------------------------------------------
    figtype='dose_response_dcseffect_mean_'
    figdf = buildfigdf._dose_response_dcseffect()
    figdf['fig_ymin']=0.79
    # figdf['fig_ymax']=1.61
    # figdf['fig_ymax']=-26
    figdf['fig_xmin']=-21
    figdf['fig_xtick_assert']=-20
    figdf['fig_ytick_assert']=1
    figdf['fig_ytick_decimals']=1
    # figdf['fig_ymax']=1.61
    figdf['fig_dxticks']=5
    figdf['fig_dyticks']=.1
    figdf['fig_ylabel']='Norm. weight mean'

    # array_funcs=[stats.kurtosis]
    # array_func_kws=[{'fisher':False}]
    # array_funcs=[stats.skew]
    # array_func_kws=[{}]
    array_funcs=[np.mean]
    array_func_kws=[{}]
    # array_funcs=[np.std]
    # array_func_kws=[{}]
    # figdf.reset_index()
    w_clopath_df=w_clopath_df.set_index(['field','path_1_syn_num'])
    figs, ax = analysis.PlotFuncs()._var2var_func(df=w_clopath_df, figdf=figdf, x_variable='field', variable='dw_clopath_dcs_diff', figformat='standard', array_funcs=array_funcs, array_func_kws=array_func_kws)
    # figs, ax = analysis.PlotFuncs()._dose_response(df=w_clopath_df, figdf=figdf, variable='dw_clopath_dcs_diff')
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi)

# _dose_response_dcseffect_mean(w_clopath_df, buildfigdf)

#############################################################################
# weight histograms
#############################################################################
def _whist_dcseffect(w_clopath_df, buildfigdf):
    '''
    '''
    # plot histogram of weights
    # ------------------------------------------------
    figtype='whist_'
    figdf=buildfigdf._weight_hist()
    figdf['fig_ymax']=0.031
    # figdf['fig_ymin']=1#0.021
    figdf['fig_xmax']=2.5
    figdf['fig_xmin']=0.5#1.011
    figdf['fig_dyticks']=0.005
    figdf['fig_dxticks']=0.5
    figdf['fig_ytick_decimals']=2
    figdf['fig_xtick_decimals']=1
    # figdf['fig_ylabel']='Weight (a.u.)'
    # figdf.drop('fig_ytick_decimals', axis=1, inplace=True)
    w_clopath_df = w_clopath_df.set_index(['field', 'path_1_syn_num', 'tree'])
    kwargs={
    'dt':1./40,
    'bins':arrayfuncs._get_bins,
    'bins_kw':{'binsize':.01},
    'normalize_weights':True}
    figs, ax = analysis.PlotFuncs()._hist(df=w_clopath_df, figdf=figdf, variables='dw_clopath_dcs_diff',**kwargs)   
    w_clopath_df = w_clopath_df.reset_index()

    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi, bbox_inches='tight') 

# _whist_dcseffect(w_clopath_df, buildfigdf)

#####################################################################
# field magnitude vs membrane polarization
#####################################################################
run_weights=False
if run_weights:

    # dose response figures
    #-------------------------------------------------------------------
    figtype='polarization_dose_response_'
    figdf = figsetup.BuildFigDF()._polarization_dose_response()
    figdf['fig_ytick_assert']=0
    figdf['fig_ytick_decimals']=1
    # figdf['fig_ymax']=0.021
    # figdf['fig_ymin']=1#0.021
    # figdf['fig_xmax']=2.5
    # figdf['fig_xmin']=0.5#1.011
    figdf['fig_dyticks']=1
    # figdf.reset_index()
    vtrace_df.set_index(['tree'], inplace=True)
    figs, ax = analysis.PlotFuncs()._var2var_mean(df=vtrace_df, figdf=figdf, variable='polarization', x_variable='field')
    # figs, ax = analysis.PlotFuncs()._draw_lines(fig=figs, ax=ax, figdf=figdf)
    vtrace_df.reset_index( inplace=True)
    plt.show(block=False)
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

#####################################################################
# weight dose response dcs 2
#####################################################################
run_weights=False
if run_weights:

    # dose response figures
    #-------------------------------------------------------------------
    figtype='dose_response2_'
    figdf = figsetup.BuildFigDF()._dose_response()
    # figdf.reset_index()
    figs, ax = analysis.PlotFuncs()._dose_response(df=w_clopath2_df, figdf=figdf, variable='dw_clopath')
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi)

#####################################################################
# dose response stdp
#####################################################################
run_weights=False
if run_weights:
    # dose response figures
    #-------------------------------------------------------------------
    figtype='dose_response_stdp_'
    figdf = figsetup.BuildFigDF()._dose_response()
    figs, ax = analysis.PlotFuncs()._dose_response(df=w_clopath_stdp_df, figdf=figdf, variable='dw_clopath')
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi)

#####################################################################
# dose response stdp2
#####################################################################
run_weights=False
if run_weights:
    # dose response figures
    #-------------------------------------------------------------------
    figtype='dose_response_stdp2_'
    figdf = figsetup.BuildFigDF()._dose_response()
    figs, ax = analysis.PlotFuncs()._dose_response(df=w_clopath_stdp2_df, figdf=figdf, variable='dw_clopath')
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi)

#####################################################################
# dose response stdp3
#####################################################################
run_weights=False
if run_weights:
    # dose response figures
    #-------------------------------------------------------------------
    figtype='dose_response_stdp3_'
    figdf = figsetup.BuildFigDF()._dose_response()
    figs, ax = analysis.PlotFuncs()._dose_response(df=w_clopath_stdp3_df, figdf=figdf, variable='dw_clopath')
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi)

#####################################################################
# dose response stdp TEMP
#####################################################################
run_weights=False
if run_weights:
    # dose response figures
    #-------------------------------------------------------------------
    figtype='dose_response_stdp_temp_'
    figdf = figsetup.BuildFigDF()._dose_response()
    # w_clopath_stdp_df_temp.reset_index(inplace=True)
    figs, ax = analysis.PlotFuncs()._dose_response(df=w_clopath_stdp_df_temp, figdf=figdf, variable='dw_clopath')
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi)

#############################################################################
# weight histograms
#############################################################################
run_hist=False
if run_hist:
    # plot histogram of weights
    # ------------------------------------------------
    figtype='whist_'
    figdf=figsetup.BuildFigDF()._weight_hist()
    figdf['fig_ymax']=0.021
    # figdf['fig_ymin']=1#0.021
    figdf['fig_xmax']=2.5
    figdf['fig_xmin']=0.5
    figdf['fig_dyticks']=0.005
    # figdf['fig_ytick_decimals']=1
    # figdf['fig_ylabel']='Weight (a.u.)'
    # figdf.drop('fig_ytick_decimals', axis=1, inplace=True)
    w_clopath_df = w_clopath_df.set_index(['field', 'path_1_syn_num', 'tree'])
    kwargs={'dt':1./40}
    figs, ax = analysis.PlotFuncs()._hist(df=w_clopath_df, figdf=figdf, variable='dw_clopath_dcs_diff',**kwargs)   
    w_clopath_df = w_clopath_df.reset_index()

    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi) 

#############################################################################
# weight histograms for stdp parameters
#############################################################################
run_hist=False
if run_hist:
    # plot histogram of weights
    # ------------------------------------------------
    figtype='whist_stdp_'
    figdf=figsetup.BuildFigDF()._weight_hist()
    figdf['fig_ymax']=0.021
    # figdf['fig_ymin']=1#0.021
    figdf['fig_xmax']=2.5
    figdf['fig_xmin']=0.5
    figdf['fig_dyticks']=0.005
    # figdf['fig_ytick_decimals']=1
    # figdf['fig_ylabel']='Weight (a.u.)'
    # figdf.drop('fig_ytick_decimals', axis=1, inplace=True)
    w_clopath_stdp_df = w_clopath_stdp_df.set_index(['field', 'path_1_syn_num', 'tree'])
    kwargs={'dt':1./40}
    figs, ax = analysis.PlotFuncs()._hist(df=w_clopath_stdp_df, figdf=figdf, variable='dw_clopath_dcs_diff',**kwargs)   
    w_clopath_stdp_df = w_clopath_stdp_df.reset_index()

    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi) 

#############################################################################
# weight histograms for stdp3 parameters
#############################################################################
run_hist=False
if run_hist:
    # plot histogram of weights
    # ------------------------------------------------
    figtype='whist_stdp3_'
    figdf=figsetup.BuildFigDF()._weight_hist()
    figdf['fig_ymax']=0.021
    # figdf['fig_ymin']=1#0.021
    figdf['fig_xmax']=2.5
    figdf['fig_xmin']=0.5
    figdf['fig_dyticks']=0.005
    # figdf['fig_ytick_decimals']=1
    # figdf['fig_ylabel']='Weight (a.u.)'
    # figdf.drop('fig_ytick_decimals', axis=1, inplace=True)
    w_clopath_stdp3_df = w_clopath_stdp3_df.set_index(['field', 'path_1_syn_num', 'tree'])
    kwargs={'dt':1./40}
    figs, ax = analysis.PlotFuncs()._hist(df=w_clopath_stdp3_df, figdf=figdf, variable='dw_clopath_dcs_diff',**kwargs)   
    w_clopath_stdp3_df = w_clopath_stdp3_df.reset_index()

    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi) 

#####################################################################
# vtrace raw data
#####################################################################
run_vtrace=False
if run_vtrace:
    # plot average voltage trace at soma and dendrite
    #------------------------------------------------
    figtype='vtrace_'
    figdf=figsetup.BuildFigDF()._trace_mean()
    figdf['fig_ymin']=-3
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
    # figdf['fig_ymin']=
    figdf['fig_dyticks']=1
    vtrace_df = vtrace_df.set_index(['field', 'path_1_syn_num', 'tree'])
    # array functions to pass to plot
    #------------------------------------
    # lowpass filter
    fs=40000.
    filt = analysis.Filters(fs=fs).filters_series['iir_band_5_50']
    filt_kws = {'filt_series':filt, 'hilbert':False, 'islice':slice(None)}
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