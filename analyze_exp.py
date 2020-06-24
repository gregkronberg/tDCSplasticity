"""
analysis

"""
# import figsetup_exp_reduced_neuron_polarization as figsetup
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
import functions
import inspect
import figsetup

class AnalysisExp(object):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        self.experiment_name = self.__class__.__name__.split('analysis_')[-1]
        self.data_directory = 'Data/'+self.experiment_name+'/'
        self.figure_directory =  'png figures/'+self.experiment_name+'/'
        self.dpi=350
                # make figure directory if doesnt exist
        if os.path.isdir(self.figure_directory) is False:
            os.mkdir(self.figure_directory)

        # plot and df functions
        #---------------------------------
        self.plt_funcs = functions.PlotFuncs()
        self.df_funcs = functions.DfFuncs()
        self.figsetup = getattr(figsetup, 'figsetup_'+self.experiment_name)()

        # generate all figures
        #--------------------------------------------------------------------
        if 'generate_figures' in kwargs and kwargs['generate_figures']:
            self.generate_figures(**kwargs)

    # check globals for variables to avoid reloading large data structures
    def _load_variables(self, variables=[], variables_reload=[], extension='.pkl', **kwargs):
        '''
        '''
        # ensure that variables is iterable
        if type(variables)!=list:
            variables=[variables]
        # iterate over group variables
        for variable in variables:
            # check if variable already exists as an attribute or it is to be reloaded
            if not hasattr(self, variable) or variable in variables_reload:
                filename = variable
                # load group variable
                loaded_variable = functions._load_group_data(directory=self.data_directory, filename=filename)
                # set group variable as an attribute
                setattr(self, variable, loaded_variable)

    def _save_figures(self, figs, figtype, fig_format='png', **kwargs):
        '''
        '''
        # make sure figtype ends in underscore
        if figtype[-1]!='_':
            figtype+='_'

        # save figure
        #------------
        for fig_key, fig in figs.iteritems():
            filename = self.figure_directory+figtype+str(fig_key)+'.'+fig_format
            fig.savefig(filename, format='png', dpi=self.dpi)

    def generate_all_figures(self, **kwargs):
        '''
        '''
        print 'no figure generating methods have been specified'

    def shapeplot(self, plot_variable='polarization', df_variable='vtrace_df', group_by=['field'], figsetup_func=None, fig_descriptor=None, **kwargs):
        '''
        '''
        # get figtype for saving figures
        if fig_descriptor is not None:
            figtype = inspect.stack()[0][3]+'_'+fig_descriptor+'_'
        else:
            figtype = inspect.stack()[0][3]+'_'
        # load voltage trace data
        self._load_variables(variables=df_variable)
        # get data df based on df_variable argument
        data_df = getattr(self, df_variable)
        # check that variable to plot exists in group data
        assert plot_variable in data_df, variable +' not found in group data'
        # set index for group data plotting (index must match figdf)
        data_df = functions._set_index(data_df, group_by)
        # create figdf
        if figsetup_func is None:
            figsetup_func = inspect.stack()[0][3] + '_'+plot_variable
        figdf = getattr(self.figsetup, figsetup_func)()
        # create attribute for the current figdf
        setattr(self, 'figdf_'+inspect.stack()[0][3], figdf)
        # create figures
        figs, ax = self.plt_funcs._shapeplot(df=data_df, figdf=figdf, variable=plot_variable, **kwargs)
        # save figures
        #------------
        self._save_figures(figs=figs, figtype=figtype)

        return figs, ax

    def dose_response(self, y_variable='polarization', x_variable='field', df_variable='vtrace_df', group_by=['tree_key', 'sec_num', 'seg_num'], **kwargs):
        '''
        '''
        # get figtype for saving figures
        figtype = inspect.stack()[0][3] + '_'
        # load voltage trace data
        self._load_variables(variables=df_variable)
        # get data df based on df_variable argument
        data_df = getattr(self, df_variable)
        # check that variable to plot exists in group data
        assert y_variable in data_df, y_variable +' not found in group data'
        assert x_variable in data_df, x_variable +' not found in group data'
        # set index for group data plotting (index must match figdf)
        data_df = functions._set_index(data_df, group_by)
        # create figdf
        figdf = getattr(self.figsetup, 'dose_response_'+y_variable)()
        # create attribute for the current figdf
        setattr(self, 'figdf_'+inspect.stack()[0][3], figdf)
        # create figures
        figs, ax = self.plt_funcs._var2var_mean(df=data_df, figdf=figdf, variable=y_variable, x_variable=x_variable, **kwargs)
        # save figures
        #------------
        self._save_figures(figs=figs, figtype=figtype)

    def var2var(self, y_variable=None, x_variable=None, df_variable='vtrace_df', group_by=['field', 'tree_key', 'sec_num', 'seg_num'], figdf_method=None, **kwargs):
        '''
        '''
        # get figtype for saving figures
        figtype = inspect.stack()[0][3] + '_'
        # load voltage trace data
        if not hasattr(self, df_variable):
            self._load_variables(variables=df_variable)
        # get data df based on df_variable argument
        data_df = getattr(self, df_variable)
        # check that variable to plot exists in group data
        assert y_variable in data_df, y_variable +' not found in group data'
        assert x_variable in data_df, x_variable +' not found in group data'
        # set index for group data plotting (index must match figdf)
        data_df = functions._set_index(data_df, group_by)
        # create figdf
        if figdf_method is None:
            figdf_method = 'var2var_'+x_variable+'_'+y_variable
        figdf = getattr(self.figsetup, figdf_method)()
        # create attribute for the current figdf
        setattr(self, 'figdf_'+inspect.stack()[0][3], figdf)
        # create figures
        figs, ax = self.plt_funcs._var2var_mean(df=data_df, figdf=figdf, variable=y_variable, x_variable=x_variable, **kwargs)
        # save figures
        #------------
        self._save_figures(figs=figs, figtype=figtype)

    def trace_mean(self, variables, df_variable='vtrace_df', xvals=None, xvals_kw={}, figdf_method=None, group_by=None, **kwargs):
        '''
        '''
        plot_type = inspect.stack()[0][3]
        # get figtype for saving figures
        figtype = plot_type + '_'
        # load voltage trace data
        if not hasattr(self, df_variable):
            self._load_variables(variables=df_variable)
        # get data df based on df_variable argument
        data_df = getattr(self, df_variable)
        # check that variable to plot exists in group data
        # assert variables in data_df, y_variable +' not found in group data'
        # assert x_variable in data_df, x_variable +' not found in group data'
        # set index for group data plotting (index must match figdf)
        data_df = functions._set_index(data_df, group_by)
        # create figdf
        if figdf_method is None:
            figdf_method = plot_type
            if type(variables)==list:
                for variable in variables:
                    figdf_method += '_'+ variable
            else:
                figdf_method += '_'+ variable

        figdf = getattr(self.figsetup, figdf_method)()
        # create attribute for the current figdf
        setattr(self, 'figdf_'+inspect.stack()[0][3], figdf)
        # create figures
        figs, ax = self.plt_funcs._trace_mean(df=data_df, figdf=figdf, variables=variables, xvals=xvals, xvals_kw=xvals_kw, **kwargs)
        # figs, ax = self.plt_funcs._var2var_mean(df=data_df, figdf=figdf, variable=y_variable, x_variable=x_variable, **kwargs)
        # save figures
        #------------
        self._save_figures(figs=figs, figtype=figtype)

    # def trace_individual(self, variables, df_variable='vtrace_df', xvals=None, xvals_kw={}, figdf_method=None, group_by=None, group_by_individual=None, individual_column='trial_id',  index=None, random=False, **kwargs):
    #     '''
    #     '''
    #     plot_type = inspect.stack()[0][3]
    #     # get figtype for saving figures
    #     figtype = plot_type + '_'
    #     # load voltage trace data
    #     if not hasattr(self, df_variable):
    #         self._load_variables(variables=df_variable)
    #     # get data df based on df_variable argument
    #     data_df = getattr(self, df_variable)
    #     # check that variable to plot exists in group data
    #     # assert variables in data_df, y_variable +' not found in group data'
    #     # assert x_variable in data_df, x_variable +' not found in group data'
    #     # set index for group data plotting (index must match figdf)
    #     data_df = functions._set_index(data_df, group_by)
    #     # create figdf
    #     if figdf_method is None:
    #         figdf_method = plot_type
    #         if type(variables)==list:
    #             for variable in variables:
    #                 figdf_method += '_'+ variable
    #         else:
    #             figdf_method += '_'+ variable

    #     figdf = getattr(self.figsetup, figdf_method)(**kwargs)
    #     # create attribute for the current figdf
    #     setattr(self, 'figdf_'+inspect.stack()[0][3], figdf)
    #     # create figures
    #     figs, ax = self.plt_funcs._trace_mean(df=data_df, figdf=figdf, variables=variables, xvals=xvals, xvals_kw=xvals_kw, **kwargs)
    #     # figs, ax = self.plt_funcs._var2var_mean(df=data_df, figdf=figdf, variable=y_variable, x_variable=x_variable, **kwargs)
    #     # save figures
    #     #------------
    #     self._save_figures(figs=figs, figtype=figtype)

    def bar(self, variable, df_variable, figdf_method=None, group_by=None, **kwargs):
        '''
        '''
        plot_type = inspect.stack()[0][3]
        # get figtype for saving figures
        figtype = plot_type + '_'
        # load voltage trace data
        if not hasattr(self, df_variable):
            self._load_variables(variables=df_variable)
        # get data df based on df_variable argument
        data_df = getattr(self, df_variable)
        # check that variable to plot exists in group data
        # assert variables in data_df, y_variable +' not found in group data'
        # assert x_variable in data_df, x_variable +' not found in group data'
        # set index for group data plotting (index must match figdf)
        data_df = functions._set_index(data_df, group_by)
        # create figdf
        if figdf_method is None:
            figdf_method = plot_type+'_'+variable
        figdf = getattr(self.figsetup, figdf_method)()
        # create attribute for the current figdf
        setattr(self, 'figdf_'+figdf_method, figdf)
        # create figures
        figs, ax = self.plt_funcs._bar(df=data_df, figdf=figdf, variable=variable, group_space=1, bar_width=1, bar_spacing=1, **kwargs)

        # save figures
        #------------
        self._save_figures(figs=figs, figtype=figtype)

class analysis_exp_reduced_neuron_tbs_conjunctive(AnalysisExp):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        super(analysis_exp_reduced_neuron_tbs_conjunctive, self).__init__(**kwargs)
        if 'variables' in kwargs:
            self._load_variables(**kwargs)

    def generate_figures(self, **kwargs):
        '''
        '''
        # weight plots
        #--------------------------------------------------------------------
        # self._load_variables(variables='w_clopath_df')
        # self.w_clopath_df['path_1_w_mean'] = 10000.*self.w_clopath_df['path_1_w_mean']
        # self.w_clopath_df['path_2_w_mean'] = 10000.*self.w_clopath_df['path_2_w_mean']
        # weight change as a function of initial weight mean
        #-------------------------------------------
        # self.w_clopath_df = functions._set_index(self.w_clopath_df, None)
        # self.var2var(y_variable='dw_clopath', x_variable='path_2_w_mean', df_variable='w_clopath_df', group_by=['field', 'tree_key', 'path_1_w_mean', 'seg_num'], figdf_method='var2var_path_1_w_mean_dw_clopath', **kwargs)
        # weight change as a function of offset delay between pathways
        #-------------------------------------------
        # self.w_clopath_df['path_2_warmup'] = self.w_clopath_df['path_2_warmup']-20
        # self.w_clopath_df = functions._set_index(self.w_clopath_df, None)
        # self.var2var(y_variable='dw_clopath', x_variable='path_2_warmup', df_variable='w_clopath_df', group_by=['field', 'tree_key','path_2_w_mean', 'path_1_w_mean', 'seg_num'], figdf_method='var2var_offset_dw_clopath', **kwargs)
        # weights over time
        #------------------
        # self.trace_mean(variables=['w_clopath'],  df_variable='w_clopath_df', group_by=['field', 'tree_key', 'path_2_w_mean', 'path_1_w_mean', 'seg_num'], figdf_method='trace_mean_w_clopath', dt=.025, figformat='standard', **kwargs)
        #weight bar plots
        #----------------
        # self.bar(variable='dw_clopath',  df_variable='w_clopath_df', group_by=['field', 'tree_key', ], figdf_method='bar_dw_clopath', **kwargs)

        # voltage plots
        #--------------------------------------------------------------------
        self._load_variables(variables='vtrace_df')
        self.vtrace_df['path_1_w_mean'] = 10000.*self.vtrace_df['path_1_w_mean']
        self.vtrace_df['path_2_w_mean'] = 10000.*self.vtrace_df['path_2_w_mean']
        self.vtrace_df['path_2_warmup'] = self.vtrace_df['path_2_warmup']-20
        # # average votlage traces in dendrites
        # #-----------------------------------
        # self.vtrace_df = functions._set_index(self.vtrace_df, None)
        # self.trace_mean(variables=['data_v'],  df_variable='vtrace_df', group_by=['field', 'tree_key', 'path_2_w_mean', 'path_1_w_mean', 'seg_num'], figdf_method='trace_mean_v', dt=.025, figformat='standard', **kwargs)
        # # average votlage traces in dendrites
        # #-----------------------------------
        self.vtrace_df = functions._set_index(self.vtrace_df, None)
        self.trace_mean(variables=['data_v'],  df_variable='vtrace_df', group_by=['field', 'tree_key', 'path_2_warmup','path_2_w_mean', 'path_1_w_mean', 'seg_num'], figdf_method='trace_individual_v', dt=.025, figformat='standard', index=0, **kwargs)

class analysis_exp_reduced_neuron_tbs_conjunctive_high_nmda(AnalysisExp):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        super(analysis_exp_reduced_neuron_tbs_conjunctive_high_nmda, self).__init__(**kwargs)
        if 'variables' in kwargs:
            self._load_variables(**kwargs)

    def generate_figures(self, **kwargs):
        '''
        '''
        # weight plots
        #--------------------------------------------------------------------
        # self._load_variables(variables='w_clopath_df')
        # self.w_clopath_df['path_1_w_mean'] = 10000.*self.w_clopath_df['path_1_w_mean']
        # self.w_clopath_df['path_2_w_mean'] = 10000.*self.w_clopath_df['path_2_w_mean']
        # weight change as a function of initial weight mean
        #-------------------------------------------
        # self.w_clopath_df = functions._set_index(self.w_clopath_df, None)
        # self.var2var(y_variable='dw_clopath', x_variable='path_2_w_mean', df_variable='w_clopath_df', group_by=['field', 'tree_key', 'path_1_w_mean', 'seg_num'], figdf_method='var2var_path_1_w_mean_dw_clopath', **kwargs)
        # weight change as a function of offset delay between pathways
        #-------------------------------------------
        # self.w_clopath_df['path_2_warmup'] = self.w_clopath_df['path_2_warmup']-20
        # self.w_clopath_df = functions._set_index(self.w_clopath_df, None)
        # self.var2var(y_variable='dw_clopath', x_variable='path_2_warmup', df_variable='w_clopath_df', group_by=['field', 'tree_key','path_2_w_mean', 'path_1_w_mean', 'seg_num'], figdf_method='var2var_offset_dw_clopath', **kwargs)
        # weights over time
        #------------------
        # self.trace_mean(variables=['w_clopath'],  df_variable='w_clopath_df', group_by=['field', 'tree_key', 'path_2_w_mean', 'path_1_w_mean', 'seg_num'], figdf_method='trace_mean_w_clopath', dt=.025, figformat='standard', **kwargs)
        #weight bar plots
        #----------------
        # self.bar(variable='dw_clopath',  df_variable='w_clopath_df', group_by=['field', 'tree_key', ], figdf_method='bar_dw_clopath', **kwargs)

        # voltage plots
        ################################################################
        self._load_variables(variables='vtrace_df')
        self.vtrace_df['path_1_w_mean'] = 10000.*self.vtrace_df['path_1_w_mean']
        self.vtrace_df['path_2_w_mean'] = 10000.*self.vtrace_df['path_2_w_mean']
        self.vtrace_df['path_2_warmup'] = self.vtrace_df['path_2_warmup']-20
        # # average votlage traces in dendrites
        # #-------------------------------------------------------------
        # self.vtrace_df = functions._set_index(self.vtrace_df, None)
        # self.trace_mean(variables=['data_v'],  df_variable='vtrace_df', group_by=['field', 'tree_key', 'path_2_w_mean', 'path_1_w_mean', 'seg_num'], figdf_method='trace_mean_v', dt=.025, figformat='standard', **kwargs)
        # individualvotlage traces in dendrites
        #-----------------------------------------------------------------
        self.vtrace_df = functions._set_index(self.vtrace_df, None)
        self.trace_mean(variables=['data_v'],  df_variable='vtrace_df', group_by=['field', 'tree_key', 'path_2_warmup','path_2_w_mean', 'path_1_w_mean', 'seg_num'], figdf_method='trace_individual_v', dt=.025, figformat='standard', index=0, **kwargs)


class analysis_exp_reduced_neuron_1hz(AnalysisExp):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        super(analysis_exp_reduced_neuron_1hz, self).__init__(**kwargs)
        if 'variables' in kwargs:
            self._load_variables(**kwargs)

    def generate_figures(self, **kwargs):
        '''
        '''
        # weight change as a function of initial weight mean
        #--------------------------------------------------------------------
        self._load_variables(variables='w_clopath_df')
        self.w_clopath_df['path_1_w_mean'] = 1000.*self.w_clopath_df['path_1_w_mean']
        # self.w_clopath_df = functions._set_index(self.w_clopath_df, 'path_1_w_mean')
        # self.w_clopath_df = self.w_clopath_df.loc[(5.9, 6., 6.1, 6.2),:]
        # self.w_clopath_df = functions._set_index(self.w_clopath_df, None)
        self.var2var(y_variable='dw_clopath', x_variable='path_1_w_mean', df_variable='w_clopath_df', group_by=['field', 'tree_key', ], figdf_method='var2var_path_1_w_mean_dw_clopath', **kwargs)
        # weights over time
        #--------------------------------------------------------------------
        self.trace_mean(variables=['w_clopath'],  df_variable='w_clopath_df', group_by=['field', 'tree_key', ], figdf_method='trace_mean_w_clopath', dt=.000025, **kwargs)
        # weight bar plots
        #--------------------------------------------------------------------
        self.bar(variable='dw_clopath',  df_variable='w_clopath_df', group_by=['field', 'tree_key', ], figdf_method='bar_dw_clopath', **kwargs)


class analysis_exp_reduced_neuron_tbs_basal(AnalysisExp):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        super(analysis_exp_reduced_neuron_tbs_basal, self).__init__(**kwargs)
        if 'variables' in kwargs:
            self._load_variables(**kwargs)

    def generate_figures(self, **kwargs):
        '''
        '''
        # weight change as a function of initial weight mean
        #--------------------------------------------------------------------
        self._load_variables(variables='w_clopath_df')
        self.w_clopath_df['path_1_w_mean'] = 1000.*self.w_clopath_df['path_1_w_mean']
        # self.w_clopath_df = functions._set_index(self.w_clopath_df, 'path_1_w_mean')
        # self.w_clopath_df = self.w_clopath_df.loc[(5.9, 6., 6.1, 6.2),:]
        self.w_clopath_df = functions._set_index(self.w_clopath_df, None)
        self.var2var(y_variable='dw_clopath', x_variable='path_1_w_mean', df_variable='w_clopath_df', group_by=['field', 'tree_key', ], figdf_method='var2var_path_1_w_mean_dw_clopath', **kwargs)
        # weights over time
        #--------------------------------------------------------------------
        self.trace_mean(variables=['w_clopath'],  df_variable='w_clopath_df', group_by=['field', 'tree_key', ], figdf_method='trace_mean_w_clopath', dt=.025, **kwargs)
        # weight bar plots
        #--------------------------------------------------------------------
        self.bar(variable='dw_clopath',  df_variable='w_clopath_df', group_by=['field', 'tree_key', ], figdf_method='bar_dw_clopath', **kwargs)

class analysis_exp_reduced_neuron_20hz_basal(AnalysisExp):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        super(analysis_exp_reduced_neuron_20hz_basal, self).__init__(**kwargs)
        if 'variables' in kwargs:
            self._load_variables(**kwargs)

    def generate_figures(self, **kwargs):
        '''
        '''
        # weight change as a function of initial weight mean
        #--------------------------------------------------------------------
        self._load_variables(variables='w_clopath_df')
        self.w_clopath_df['path_1_w_mean'] = 1000.*self.w_clopath_df['path_1_w_mean']
        # self.w_clopath_df = functions._set_index(self.w_clopath_df, 'path_1_w_mean')
        # self.w_clopath_df = self.w_clopath_df.loc[(5.9, 6., 6.1, 6.2),:]
        self.w_clopath_df = functions._set_index(self.w_clopath_df, None)
        self.var2var(y_variable='dw_clopath', x_variable='path_1_w_mean', df_variable='w_clopath_df', group_by=['field', 'tree_key', ], figdf_method='var2var_path_1_w_mean_dw_clopath', **kwargs)
        # weights over time
        #--------------------------------------------------------------------
        self.trace_mean(variables=['w_clopath'],  df_variable='w_clopath_df', group_by=['field', 'tree_key', ], figdf_method='trace_mean_w_clopath', dt=.025, **kwargs)
        # weight bar plots
        #--------------------------------------------------------------------
        self.bar(variable='dw_clopath',  df_variable='w_clopath_df', group_by=['field', 'tree_key', ], figdf_method='bar_dw_clopath', **kwargs)

class analysis_exp_reduced_neuron_tbs(AnalysisExp):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        super(analysis_exp_reduced_neuron_tbs, self).__init__(**kwargs)
        if 'variables' in kwargs:
            self._load_variables(**kwargs)

    def generate_figures(self, **kwargs):
        '''
        '''
        # weight change as a function of initial weight mean
        #--------------------------------------------------------------------
        self._load_variables(variables='w_clopath_df')
        self.w_clopath_df['path_1_w_mean'] = 1000.*self.w_clopath_df['path_1_w_mean']
        self.w_clopath_df = functions._set_index(self.w_clopath_df, 'path_1_w_mean')
        self.w_clopath_df = self.w_clopath_df.loc[(5.9, 6., 6.1, 6.2),:]
        self.w_clopath_df = functions._set_index(self.w_clopath_df, None)
        self.var2var(y_variable='dw_clopath', x_variable='path_1_w_mean', df_variable='w_clopath_df', group_by=['field', 'tree_key', ], figdf_method='var2var_path_1_w_mean_dw_clopath', **kwargs)
        # weights over time
        #--------------------------------------------------------------------
        self.trace_mean(variables=['w_clopath'],  df_variable='w_clopath_df', group_by=['field', 'tree_key', ], figdf_method='trace_mean_w_clopath', dt=.025, **kwargs)
        # weight bar plots
        #--------------------------------------------------------------------
        self.bar(variable='dw_clopath',  df_variable='w_clopath_df', group_by=['field', 'tree_key', ], figdf_method='bar_dw_clopath', **kwargs)

class analysis_exp_reduced_neuron_20hz(AnalysisExp):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        super(analysis_exp_reduced_neuron_20hz, self).__init__(**kwargs)
        if 'variables' in kwargs:
            self._load_variables(**kwargs)

    def generate_figures(self, **kwargs):
        '''
        '''
        # weight change as a function of initial weight mean
        #--------------------------------------------------------------------
        self._load_variables(variables='w_clopath_df')
        self.w_clopath_df['path_1_w_mean'] = 1000.*self.w_clopath_df['path_1_w_mean']
        # self.w_clopath_df = functions._set_index(self.w_clopath_df, 'path_1_w_mean')
        # self.w_clopath_df = self.w_clopath_df.loc[(5.9, 6., 6.1, 6.2),:]
        # self.w_clopath_df = functions._set_index(self.w_clopath_df, None)
        self.var2var(y_variable='dw_clopath', x_variable='path_1_w_mean', df_variable='w_clopath_df', group_by=['field', 'tree_key', ], figdf_method='var2var_path_1_w_mean_dw_clopath', **kwargs)
        # weights over time
        #--------------------------------------------------------------------
        self.trace_mean(variables=['w_clopath'],  df_variable='w_clopath_df', group_by=['field', 'tree_key', ], figdf_method='trace_mean_w_clopath', dt=.025, **kwargs)
        # weight bar plots
        #--------------------------------------------------------------------
        self.bar(variable='dw_clopath',  df_variable='w_clopath_df', group_by=['field', 'tree_key', ], figdf_method='bar_dw_clopath', **kwargs)

class analysis_exp_reduced_neuron_polarization_piecewise(AnalysisExp):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        super(analysis_exp_reduced_neuron_polarization_piecewise, self).__init__(**kwargs)
        if 'variables' in kwargs:
            self._load_variables(**kwargs)

    def generate_figures(self, **kwargs):
        '''
        '''
        # shapeplot of polarization mapped from original cell to reduced cell
        #--------------------------------------------------------------------
        # self.shapeplot(plot_variable='polarization_original_cell_mean',df_variable='vtrace_df', group_by=['field'], figsetup_func='shapeplot_polarization', width_scale=50, fig_descriptor='original_polarization', **kwargs)
        # polarization shape plot
        #------------------------------------------------------------------- 
        self.shapeplot(plot_variable='polarization', df_variable='vtrace_df', group_by=['field'], figsetup_func='shapeplot_polarization', width_scale=50, **kwargs)
        # polarization dose response
        #--------------------------------------------------------------------
        # self.dose_response(y_variable='polarization', x_variable='field', df_variable='vtrace_df', group_by=['tree_key', 'sec_num', 'seg_num'], **kwargs)


class analysis_exp_full_neuron_polarization_dose_response(AnalysisExp):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        super(analysis_exp_full_neuron_polarization_dose_response, self).__init__(**kwargs)
        # load variables
        if 'variables' in kwargs:
            self._load_variables(**kwargs)


    def generate_all_figures(self, **kwargs):
        '''
        '''
        # polarization shape plot
        #-------------------------------------------------------------------
        self.shapeplot(plot_variable='polarization', df_variable='vtrace_df', group_by=['field'], **kwargs)
        # polarization dose response
        #--------------------------------------------------------------------
        self.dose_response(y_variable='polarization', x_variable='field', df_variable='vtrace_df', group_by=['tree_key', 'sec_num', 'seg_num'], **kwargs)

class analysis_exp_reduced_neuron_polarization_mirror_estimate(AnalysisExp):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        super(analysis_exp_reduced_neuron_polarization_mirror_estimate, self).__init__(**kwargs)
        if 'variables' in kwargs:
            self._load_variables(**kwargs)

    def generate_figures(self, **kwargs):
        '''
        '''
        # shapeplot of polarization mapped from original cell to reduced cell
        #--------------------------------------------------------------------
        self.shapeplot(plot_variable='polarization_original_cell_mean',df_variable='vtrace_df', group_by=['field'], figsetup_func='shapeplot_polarization', width_scale=50, fig_descriptor='original_polarization', **kwargs)
        # polarization shape plot
        #-------------------------------------------------------------------
        self.shapeplot(plot_variable='polarization', df_variable='vtrace_df', group_by=['field'], figsetup_func='shapeplot_polarization', width_scale=50, **kwargs)
        # polarization dose response
        #--------------------------------------------------------------------
        self.dose_response(y_variable='polarization', x_variable='field', df_variable='vtrace_df', group_by=['tree_key', 'sec_num', 'seg_num'], **kwargs)

class analysis_exp_reduced_neuron_polarization(AnalysisExp):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        super(analysis_exp_reduced_neuron_polarization, self).__init__(**kwargs)
        if 'variables' in kwargs:
            self._load_variables(**kwargs)
            # plot and df functions
        #---------------------------------
        self.plt_funcs = functions.PlotFuncs()
        self.df_funcs = functions.DfFuncs()

    def polarization_dose_response(self, **kwargs):
        '''
        '''
        if not hasattr(self, 'vtrace_df'):
            self._load_variables(variables=['vtrace_df'])
        # plot and df functions
        #---------------------------------
        plt_funcs = functions.PlotFuncs()
        df_funcs = functions.DfFuncs()

        # create column for polarization
        #--------------------------------
        self.vtrace_df = self.df_funcs._get_polarization(df=self.vtrace_df)
        # set index for plotting (index must match figdf)
        #--------------------------------------------------------------
        self.vtrace_df = functions._set_index(self.vtrace_df, ['tree_key', 'sec_num', 'seg_num'])
        # setup figdf
        #--------------------------------------------------------------
        self.figdf = figsetup.BuildFigDF().polarization_dose_response()
        # create plots
        #--------------------------------------------------------------
        self.figs, self.ax = self.plt_funcs._var2var_mean(df=self.vtrace_df, figdf=self.figdf, variable='polarization', x_variable='field')
        # save figure
        #--------------------------------------------------------------
        # get name of the current plotting method
        figtype = inspect.stack()[0][3]
        # iterate figures
        for fig_key, fig in self.figs.iteritems():
            # figure name
            fname = self.figure_directory+figtype+str(fig_key)+'.png'
            # save
            fig.savefig(fname, format='png', dpi=self.dpi)

        return self.figs, self.ax

    def polarization_shapeplot(self, **kwargs):
        '''
        '''
        if not hasattr(self, 'vtrace_df'):
            self._load_variables(variables=['vtrace_df'])
            # create column for polarization
            #--------------------------------
            self.vtrace_df = self.df_funcs._get_polarization(df=self.vtrace_df)

        # set index for plotting (index must match figdf)
        #--------------------------------------------------------------
        self.vtrace_df = functions._set_index(self.vtrace_df, ['field'])
        figtype='polarization_shapeplot_'
        self.figdf = figsetup.BuildFigDF().polarization_shapeplot()
        self.figs, self.ax = self.plt_funcs._shapeplot(df=self.vtrace_df, figdf=self.figdf, variable='polarization', width_scale=10)

        # save figure
        #------------
        for fig_key, fig in self.figs.iteritems():
            fname = self.figure_directory+figtype+str(fig_key)+'.png'
            fig.savefig(fname, format='png', dpi=self.dpi)

        return self.figs, self.ax

        # plt.figure()
        # locations = [('soma',0,0), ('apic',0, 10), ('dend',0, 2)]
        # for location in locations:
        #     plt.plot(self.vtrace_df.loc[[location]].field, self.vtrace_df.loc[[location]].polarization)
        # plt.show(block=False)

class analysis_exp_full_neuron_polarization(AnalysisExp):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        super(analysis_exp_full_neuron_polarization, self).__init__(**kwargs)
        if 'variables' in kwargs:
            self._load_variables(**kwargs)
            # plot and df functions
        #---------------------------------
        self.plt_funcs = functions.PlotFuncs()
        self.df_funcs = functions.DfFuncs()

    def polarization_dose_response(self, **kwargs):
        '''
        '''
        if not hasattr(self, 'vtrace_df'):
            self._load_variables(variables=['vtrace_df'])
        # plot and df functions
        #---------------------------------
        plt_funcs = functions.PlotFuncs()
        df_funcs = functions.DfFuncs()

        # create column for polarization
        #--------------------------------
        self.vtrace_df = self.df_funcs._get_polarization(df=self.vtrace_df)
        # set index for plotting (index must match figdf)
        #--------------------------------------------------------------
        self.vtrace_df = functions._set_index(self.vtrace_df, ['tree_key', 'sec_num', 'seg_num'])
        # setup figdf
        #--------------------------------------------------------------
        self.figdf = figsetup.BuildFigDF().polarization_dose_response()
        # create plots
        #--------------------------------------------------------------
        self.figs, self.ax = self.plt_funcs._var2var_mean(df=self.vtrace_df, figdf=self.figdf, variable='polarization', x_variable='field')
        # save figure
        #--------------------------------------------------------------
        # get name of the current plotting method
        figtype = inspect.stack()[0][3]
        # iterate figures
        for fig_key, fig in self.figs.iteritems():
            # figure name
            fname = self.figure_directory+figtype+str(fig_key)+'.png'
            # save
            fig.savefig(fname, format='png', dpi=self.dpi)

        return self.figs, self.ax

    def polarization_shapeplot(self, **kwargs):
        '''
        '''
        if not hasattr(self, 'vtrace_df'):
            self._load_variables(variables=['vtrace_df'])

        if 'polarization' not in self.vtrace_df.columns:
            # create column for polarization
            #--------------------------------
            self.vtrace_df = self.df_funcs._get_polarization(df=self.vtrace_df)

        # set index for plotting (index must match figdf)
        #--------------------------------------------------------------
        self.vtrace_df = functions._set_index(self.vtrace_df, ['field'])
        figtype='polarization_shapeplot_'
        self.figdf = figsetup.BuildFigDF().polarization_shapeplot()
        self.figs, self.ax = self.plt_funcs._shapeplot(df=self.vtrace_df, figdf=self.figdf, variable='polarization', width_scale=10)

        # save figure
        #------------
        for fig_key, fig in self.figs.iteritems():
            fname = self.figure_directory+figtype+str(fig_key)+'.png'
            fig.savefig(fname, format='png', dpi=self.dpi)

        return self.figs, self.ax



