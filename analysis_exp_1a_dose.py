"""
analysis

"""
import figsetup_exp_1a_dose as figsetup
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
kwargs = {'experiment':'exp_1a_dose'}
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

def _load_variables_basal():
    '''
    '''
    vtrace_df_basal = analysis._load_group_data(directory=directory_basal, filename=filename_vtrace)
    w_clopath_df_basal = analysis._load_group_data(directory=directory_basal, filename=filename_w_clopath)
    return vtrace_df_basal, w_clopath_df_basal

# vtrace_df_basal, w_clopath_df_basal = _load_variables_basal()

#####################################################################
# field magnitude vs membrane polarization
#####################################################################
def _polarization_dose_response():
    '''
    '''
    # run_weights=True
    # if run_weights:

    # dose response figures
    #-------------------------------------------------------------------
    figtype='polarization_dose_response_combined_'
    figdf = figsetup.BuildFigDF()._polarization_dose_response()
    figdf['fig_ytick_assert']=-1
    figdf['fig_ytick_decimals']=1
    figdf['fig_xtick_assert']=20
    figdf['fig_xtick_decimals']=0
    figdf['fig_ymax']=1.01
    figdf['fig_ymin']=-1.01#0.021
    figdf['fig_xmax']=20.1
    figdf['fig_xmin']=-20.1#1.011
    figdf['fig_dyticks']=0.5
    figdf['fig_dxticks']=5
    # figdf.reset_index()

    # combine apical and basal df
    #-------------------------------
    vtrace_df_crop = vtrace_df[['tree', 'polarization','field']].copy().reset_index()
    vtrace_df_basal_crop = vtrace_df_basal[['tree', 'polarization','field']].copy().reset_index()
    vtrace_df_merge = vtrace_df_crop.append(vtrace_df_basal_crop, ignore_index=True)
    vtrace_df_merge.set_index(['tree'], inplace=True)
    figs, ax = analysis.PlotFuncs()._var2var_mean(df=vtrace_df_merge, figdf=figdf, variable='polarization', x_variable='field')
    figs, ax = analysis.PlotFuncs()._draw_lines(fig=figs, ax=ax, figdf=figdf)
    # vtrace_df.reset_index( inplace=True)
    plt.show(block=False)
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi)

_polarization_dose_response()

#####################################################################
# weight dose response
#####################################################################
run_weights=False
if run_weights:

    # dose response figures
    #-------------------------------------------------------------------
    figtype='dose_response_'
    figdf = figsetup.BuildFigDF()._dose_response()
    figdf['fig_ymin']=1.3
    figdf['fig_ymax']=1.61
    figdf['fig_dyticks']=.1

    # figdf.reset_index()
    figs, ax = analysis.PlotFuncs()._dose_response(df=w_clopath_df, figdf=figdf, variable='dw_clopath')
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi)

#####################################################################
# weight dose response
#####################################################################
# def _dose_response_dcseffect_std(w_clopath_df):

#     # dose response figures
#     #-------------------------------------------------------------------
#     figtype='dose_response_dcseffect_std'
#     figdf = figsetup.BuildFigDF()._dose_response_dcseffect()
#     # figdf['fig_ymin']=1.3
#     # figdf['fig_ymax']=1.61
#     # figdf['fig_dyticks']=.1

#     array_funcs=[np.std]
#     array_func_kws=[{}]
#     # figdf.reset_index()
#     figs, ax = analysis.PlotFuncs()._dose_response(df=w_clopath_df, figdf=figdf, variable='dw_clopath_dcs_diff', figformat=None, array_funcs=array_funcs, array_func_kws=array_func_kws)
#     # figs, ax = analysis.PlotFuncs()._dose_response(df=w_clopath_df, figdf=figdf, variable='dw_clopath_dcs_diff')
#     # save figure
#     #------------
#     for fig_key, fig in figs.iteritems():
#         fname = figure_directory+figtype+str(fig_key)+'.png'
#         fig.savefig(fname, format='png', dpi=dpi)

# _dose_response_dcseffect_std(w_clopath_df)
def _dose_response_dcseffect_std(w_clopath_df):

    # dose response figures
    #-------------------------------------------------------------------
    figtype='dose_response_dcseffect_std_'
    figdf = figsetup.BuildFigDF()._dose_response_dcseffect()
    figdf['fig_ymin']=0
    # figdf['fig_ymax']=1.61
    # figdf['fig_ymax']=-26
    figdf['fig_xmin']=-26
    figdf['fig_xtick_assert']=-25
    # figdf['fig_ymax']=1.61
    figdf['fig_dxticks']=5
    figdf['fig_ylabel']='Norm. weight std'

    # array_funcs=[stats.kurtosis]
    # array_func_kws=[{'fisher':False}]
    # array_funcs=[stats.skew]
    # array_func_kws=[{}]
    array_funcs=[np.std]
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

# _dose_response_dcseffect_std(w_clopath_df)

def _dose_response_dcseffect_mean(w_clopath_df):

    # dose response figures
    #-------------------------------------------------------------------
    figtype='dose_response_dcseffect_mean_'
    figdf = figsetup.BuildFigDF()._dose_response_dcseffect()
    # figdf['fig_ymin']=1.3
    # figdf['fig_ymax']=1.61
    # figdf['fig_ymax']=-26
    figdf['fig_xmin']=-26
    figdf['fig_xtick_assert']=-25
    # figdf['fig_ymax']=1.61
    figdf['fig_dxticks']=5
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

# _dose_response_dcseffect_mean(w_clopath_df)

def _dose_response_dcseffect_kurtosis(w_clopath_df):

    # dose response figures
    #-------------------------------------------------------------------
    figtype='dose_response_dcseffect_kurtosis_'
    figdf = figsetup.BuildFigDF()._dose_response_dcseffect()
    figdf['fig_ymin']=0
    figdf['fig_ymax']=101
    # figdf['fig_ymax']=-26
    figdf['fig_xmin']=-21
    figdf['fig_xtick_assert']=-20
    # figdf['fig_ymax']=1.61
    figdf['fig_dxticks']=10
    figdf['fig_dyticks']=20
    figdf['fig_ylabel']='Norm. weight kurtosis'
    figdf['fig_ytick_decimals']=0

    array_funcs=[stats.kurtosis]
    array_func_kws=[{'fisher':False}]
    # array_funcs=[stats.skew]
    # array_func_kws=[{}]
    # array_funcs=[np.mean]
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

# _dose_response_dcseffect_kurtosis(w_clopath_df)

#####################################################################
# weight dose response low intensity
#####################################################################
run_weights=False
if run_weights:

    # dose response figures
    #-------------------------------------------------------------------
    figtype='dose_response_low_intensity_'
    figdf = figsetup.BuildFigDF()._dose_response()
    figdf['fig_xmin']=-5.
    figdf['fig_xmax']=6.
    # figdf.reset_index()
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
# weight dose response stdp
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
# weight dose response stdp2
#####################################################################
run_weights=False
if run_weights:
    # dose response figures
    #-------------------------------------------------------------------
    figtype='dose_response_stdp_'
    figdf = figsetup.BuildFigDF()._dose_response()
    figs, ax = analysis.PlotFuncs()._dose_response(df=w_clopath_stdp2_df, figdf=figdf, variable='dw_clopath')
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi)

#############################################################################
# weight histograms
#############################################################################
def _whist_dcseffect(w_clopath_df):
    '''
    '''
    # plot histogram of weights
    # ------------------------------------------------
    figtype='whist_'
    figdf=figsetup.BuildFigDF()._weight_hist()
    figdf['fig_ymax']=0.021
    # figdf['fig_ymin']=1#0.021
    figdf['fig_xmax']=2.5
    figdf['fig_xmin']=0.5#1.011
    figdf['fig_dyticks']=0.005
    figdf['fig_dxticks']=0.5
    # figdf['fig_ytick_decimals']=1
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

# _whist_dcseffect(w_clopath_df=w_clopath_df)

#####################################################################
# plot voltage trace for example with large effect
#####################################################################
run_vtrace=False
if run_vtrace:
    # plot average voltage trace at soma and dendrite
    #------------------------------------------------
    figtype='vtrace_individual_large_effect_'
    effect_threshold=2
    effect_range = [1.2,1.3]
    i_1vm = w_clopath_df.field==20
    w_1vm = w_clopath_df[i_1vm]
    large_effect_i = (w_1vm.dw_clopath_dcs_diff>effect_range[0]) & (w_1vm.dw_clopath_dcs_diff<effect_range[1]) 
    w_large_df = w_1vm[large_effect_i]
    trial_id = w_large_df.trial_id.unique()[0]
    # trial_id = vtrace_df.trial_id.iloc[1100]
    figtype='vtrace_individual_'+trial_id+'_'
    figdf=figsetup.BuildFigDF()._trace_individual(trial_id=trial_id)
    figdf['fig_ymax']=-40
    # figdf['fig_dyticks']=1
    figdf['fig_xmin']=0
    figdf['fig_xmax']=60
    vtrace_df = vtrace_df.set_index(['field', 'tree', 'trial_id'])
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
#############################################################################
# weight histograms dcs2
#############################################################################
run_hist=False
if run_hist:
    # plot histogram of weights
    # ------------------------------------------------
    figtype='whist2_'
    figdf=figsetup.BuildFigDF()._weight_hist()
    figdf['fig_ymax']=0.021
    # figdf['fig_ymin']=1#0.021
    figdf['fig_xmax']=2.5
    figdf['fig_xmin']=0.0#1.011
    figdf['fig_dyticks']=0.005
    figdf['fig_dxticks']=0.5
    # figdf['fig_ytick_decimals']=1
    # figdf['fig_ylabel']='Weight (a.u.)'
    # figdf.drop('fig_ytick_decimals', axis=1, inplace=True)
    w_clopath2_df = w_clopath2_df.set_index(['field', 'path_1_syn_num', 'tree'])
    kwargs={'dt':1./40}
    figs, ax = analysis.PlotFuncs()._hist(df=w_clopath2_df, figdf=figdf, variable='dw_clopath_dcs_diff',**kwargs)   
    w_clopath2_df = w_clopath2_df.reset_index()

    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi, bbox_inches='tight') 

#############################################################################
# weight histograms stdp
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
    figdf['fig_xmin']=0.5#1.011
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
# weight histograms stdp2
#############################################################################
run_hist=False
if run_hist:
    # plot histogram of weights
    # ------------------------------------------------
    figtype='whist_stdp2_'
    figdf=figsetup.BuildFigDF()._weight_hist()
    figdf['fig_ymax']=0.021
    # figdf['fig_ymin']=1#0.021
    figdf['fig_xmax']=2.5
    figdf['fig_xmin']=0.5#1.011
    figdf['fig_dyticks']=0.005
    # figdf['fig_ytick_decimals']=1
    # figdf['fig_ylabel']='Weight (a.u.)'
    # figdf.drop('fig_ytick_decimals', axis=1, inplace=True)
    w_clopath_stdp2_df = w_clopath_stdp2_df.set_index(['field', 'path_1_syn_num', 'tree'])
    kwargs={'dt':1./40}
    figs, ax = analysis.PlotFuncs()._hist(df=w_clopath_stdp2_df, figdf=figdf, variable='dw_clopath_dcs_diff',**kwargs)   
    w_clopath_stdp2_df = w_clopath_stdp2_df.reset_index()

    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi) 

#####################################################################
# vtrace
#####################################################################
run_vtrace=False
if run_vtrace:
    # plot average voltage trace at soma and dendrite
    #------------------------------------------------
    figtype='vtrace_'
    figdf=analysis.BuildFigDF()._trace_mean()
    vtrace_df = vtrace_df.reset_index().set_index(['field', 'path_1_syn_num', 'tree'])
    kwargs={'dt':1./40}
    figs, ax = analysis.PlotFuncs()._trace_mean(df=vtrace_df, figdf=figdf, variable='data_v',**kwargs)   
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi) 

#####################################################################
# vtrace highpass filtered hilbert
#####################################################################
run_vtrace=False
if run_vtrace:
    # plot average voltage trace at soma and dendrite
    #------------------------------------------------
    figtype='vtrace_hilbert_'
    figdf=analysis.BuildFigDF()._trace_mean()
    figdf['fig_ymin']=0
    vtrace_df = vtrace_df.reset_index().set_index(['field', 'path_1_syn_num', 'tree'])
    kwargs={'dt':1./40}
    fs=40000.
    filt = analysis.Filters(fs=fs).filters['iir_high_300']
    kws = {'filt':filt, 'hilbert':True}
    array_funcs = [arrayfuncs._filter_mean]
    array_func_kws = [kws]
    figs, ax = analysis.PlotFuncs()._trace_mean(df=vtrace_df, figdf=figdf, variable='data_v', array_funcs=array_funcs, array_func_kws=array_func_kws, **kwargs)   
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
    figdf=analysis.BuildFigDF()._trace_mean()
    figdf['fig_ymin']=-3
    vtrace_df = vtrace_df.reset_index().set_index(['field', 'path_1_syn_num', 'tree'])
    kwargs={'dt':1./40}
    fs=40000.
    filt = analysis.Filters(fs=fs).filters_series['iir_band_5_50']
    filt_kws = {'filt_series':filt, 'hilbert':False}
    subtract_kws = {'islice':slice(399, 400), 'axis':1}
    array_funcs = [arrayfuncs._subtract_timeseries, arrayfuncs._filter_mean]
    array_func_kws = [subtract_kws, filt_kws]
    figs, ax = analysis.PlotFuncs()._trace_mean(df=vtrace_df, figdf=figdf, variable='data_v', array_funcs=array_funcs, array_func_kws=array_func_kws, **kwargs)   
    # save figure
    #------------
    for fig_key, fig in figs.iteritems():
        fname = figure_directory+figtype+str(fig_key)+'.png'
        fig.savefig(fname, format='png', dpi=dpi) 

#####################################################################
# spikes
#####################################################################
run_spikes=False
if run_spikes:
    # load group data
    #-------------------------------------------------------------
    spikes_df = analysis._load_group_data(directory=directory, filename=filename_spikes)

    # update group data vtrace
    #----------------------------------------------------------------
    functions = [funcs._get_spikes]
    kwlist = [{'threshold':-30}]
    rerun=[]
    keep=[]
    spikes_df_temp = copy.deepcopy(spikes_df)
    spikes_df = analysis._process_new_data_df(group_df=spikes_df, preprocessed_directory=directory, functions=functions, kwlist=kwlist, rerun=rerun, keep=keep, file_limit=[])
    # save group data
    #----------------
    if not spikes_temp.equals(spikes_df):
        print 'saving updated group data'
        spikes_df.to_pickle(directory+filename_spikes)