"""
analysis

"""
from __future__ import division
import numpy as np
import analysis
# import scipy
# import scipy.io
from scipy import stats
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
import functions

class FigSetup(object):
    '''
    '''

    def __init__(self, **kwargs):
        # set dpi for final png image
        #----------------------------
        self.dpi=350

        # colors for plots
        #-------------------
        self.black = (0,0,0)
        self.gray = (0.7,0.7, 0.7)
        self.red = (1,0,0)
        self.red_light = (1,0.7, 0.7)
        self.blue=(0,0,1)
        self.blue_light = (0.7, 0.7, 1)

        # create default df
        self.default_figdf = self._default_figdf()

    def _generate_figdf(self, **kwargs):
        '''
        '''
        self.figdict = self._define_figdict()
        self.figdf = self._build_figdf_from_dict(figdict=self.figdict, default_figdf=self.default_figdf)
        self.figdf = self._update_figure_level_parameters(figdf=self.figdf)
        self.figdf = self._update_trace_level_parameters(figdf=self.figdf)

    def _build_figdf_from_dict(self, figdict, default_figdf=None):
        '''
        '''
        #FIXME handle empy figdict
        # build multiindex for trace parameters
        multi_list = []
        level_names = ['figure','subgroup', 'trace']
        for level_1_key, level_1 in figdict.iteritems():
            for level_2_key, level_2 in level_1.iteritems():
                for level_3_i, level_3_key in enumerate(level_2):
                    multi_list.append((level_1_key, level_2_key, level_3_key))

        multiindex = pd.MultiIndex.from_tuples(multi_list, names=level_names)
        # build dataframe
        figdf = pd.DataFrame(index=multiindex, dtype='object')

        # add default df columns if passed as argument
        #---------------------------------------------
        if default_figdf is not None:
            # set trace level as index
            figdf = figdf.reset_index().set_index('trace')
            # get default parameters
            figdf_default = pd.concat([default_figdf]*len(figdf))
            # set index of defaultdf to match figdf
            figdf_default.index=figdf.index
            # add default df to figdf
            figdf = pd.concat([figdf, figdf_default], axis=1, ignore_index=False)

        return figdf

    def _default_figdf(self):
        '''
        '''
        default_dict={
        # hide the top and right axes boundaries
            'fig_dpi':350,
            'fig_boxoff':True,
            # axes and tick labels
            'fig_axes_linewidth':[4],
            'fig_xlabel_fontsize':[25],
            'fig_ylabel_fontsize':[25],
            'fig_xlabel_fontweight':'heavy',
            'fig_ylabel_fontweight':'heavy',
            'fig_xtick_fontsize':[15],
            'fig_ytick_fontsize':[15],
            'fig_xtick_fontweight':'heavy',
            'fig_ytick_fontweight':'heavy',
            # figure tight layout
            'fig_tight_layout':True,
        }
        default_df = pd.DataFrame(default_dict, dtype='object')

        return default_df

    def _define_fig_dict(self, **kwargs):
        '''
        '''
        figdict = {
                # figure 1
                #-------------------------------------------------
                'figure_1_name':{
                    # subgroup
                    #___________________
                    'subgroup_1_name':[
                        # trace
                        #.................
                        ('trace_1_index_1','trace_1_index_2'),
                        ('trace_2_index_1','trace_2_index_2'),
                    ],
                    'subgroup_2_name':[
                        # trace
                        #...............
                        ('trace_1_index_1','trace_1_index_2'),
                        ('trace_2_index_1','trace_2_index_2'),
                    ],
                },
                # figure 2
                #-------------------------------------------------
                'figure_1_name':{
                    # subgroup
                    #____________________
                    'subgroup_1_name':[
                        # trace
                        #..............
                        ('trace_1_index_1','trace_1_index_2'),
                        ('trace_2_index_1','trace_2_index_2'),
                    ],
                    'subgroup_2_name':[
                        # trace
                        #...............
                        ('trace_1_index_1','trace_1_index_2'),
                        ('trace_2_index_1','trace_2_index_2'),
                    ],
                },
        }

    def _update_figure_level_parameters(self, figdf, **kwargs):
        '''
        '''
        return figdf

    def _update_trace_level_parameters(self, figdf, **kwargs):
        '''
        '''
        return figdf

class figsetup_exp_reduced_neuron_tbs_conjunctive_high_nmda(FigSetup):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        super(figsetup_exp_reduced_neuron_tbs_conjunctive_high_nmda, self).__init__(**kwargs)

    def bar_dw_clopath(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            figdict = {
                # # figure
                # #-------------------------------------------------
                'all':{
                    # subgroup
                    'all':[
                        # trace
                        (5, 'dend',),
                        (0, 'dend',),
                        (-5, 'dend'),
                    ]
                },
            }

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            # figdf['fig_nytickss']=4
            # figdf['fig_nxticks']=10
            figdf['fig_dyticks']=1.
            # figdf['fig_dxticks']=10
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            figdf['fig_ymin']=1.
            figdf['fig_ymax']=7.1
            figdf['fig_xmin']=0.
            # figdf['fig_xmax']=30.
            figdf['fig_ylabel']='Norm. weight'
            # figdf['fig_xlabel']='Time (ms)'
            # # trace level parameters
            # figdf['trace_ealpha']=.7
            # figdf['error_style']='shade'
            # figdf['trace_linewidth']=4
            # figdf['fig_xscale']=1./40
            figdf['fig_barwidth']=1
            # figdf['fig_data_style']='point'
            # figdf['fig_xtick_decimals']=0
            figdf['fig_ytick_decimals']=0
            # figdf['fig_set_xscale']='symlog'

            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # set subgroup parameters
            #-----------------------------------
            figdf = figdf.reset_index().set_index('subgroup')
            figdf.at['all','sub_location']=0

            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]<0:
                    figdf.at[key, 'trace_color']=self.blue
                    figdf.at[key, 'trace_ecolor']=self.blue_light
                    figdf.at[key, 'trace_location']=0
                    figdf.at[key, 'trace_label']='cathodal'
                # apical
                if key[2][0]==0:
                    figdf.at[key, 'trace_color']=self.black
                    figdf.at[key, 'trace_ecolor']=self.gray
                    figdf.at[key, 'trace_location']=1
                    figdf.at[key, 'trace_label']='control'
                # basal
                if key[2][0]>0:
                    figdf.at[key, 'trace_color']=self.red
                    figdf.at[key, 'trace_ecolor']=self.red_light
                    figdf.at[key, 'trace_location']=2
                    figdf.at[key, 'trace_label']='anodal'
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

    def trace_mean_w_clopath(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            segs = {
            'proximal': (0,1,2,3,4,5),
            'distal':(6,7,8,9,10)
            }
            distal_weights = 10.*np.round(np.arange(3., 4., .1),5)
            proximal_weights = 10.*np.round(np.arange(1.4, 2.4, .1), 5)
            figdict={}
            for figloc in segs.keys():
                for proximal_weight in proximal_weights:
                    for distal_weight in distal_weights:
                        figkey = 'w_'+figloc+'_proxweight_'+str(proximal_weight)+'_distweight_'+str(distal_weight)
                        figdict[figkey]={
                        'all':[
                            (5, 'apic', proximal_weight,distal_weight, segs[figloc]),
                            (0, 'apic',proximal_weight, distal_weight, segs[figloc]),
                            (-5, 'apic', proximal_weight, distal_weight, segs[figloc]),
                        ]
                        }

            # figdict = {
            #     # # figure
            #     # #-------------------------------------------------
            #     'proximal_14':{
            #         # subgroup
            #         'all':[
            #             # trace
            #             (5, 'dend',(1.4, 1.5)),
            #             (0, 'dend',(1.4, 1.5)),
            #             (-5, 'dend', (1.4, 1.5)),
            #         ]
            #     },
            # }

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            # figdf['fig_nyticks']=5
            # figdf['fig_nxticks']=10
            # figdf['fig_dyticks']=.5
            # figdf['fig_dxticks']=500.
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            # figdf['fig_ymin']=0.5
            # figdf['fig_ymax']=3.1
            figdf['fig_xmin']=200
            figdf['fig_xmax']=401.
            figdf['fig_ylabel']='Weight'
            figdf['fig_xlabel']='Time (ms)'
            # figdf['fig_dyticks']=.2
            figdf['fig_dxticks']=50
            # # trace level parameters
            figdf['trace_ealpha']=1
            figdf['error_style']='shade'
            figdf['trace_linewidth']=4
            figdf['trace_linestyle']='-'
            # figdf.at[slice(None), 'trace_ecolor'] = gray
            figdf['fig_xscale_log']=False
            figdf['fig_barwidth']=0.8
            figdf['fig_data_style']='point'
            figdf['fig_xtick_decimals']=0
            figdf['fig_ytick_decimals']=1
            figdf['fig_set_xscale']='linear'
            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]<0:
                    figdf.at[key, 'trace_color']=self.blue
                    figdf.at[key, 'trace_ecolor']=self.blue_light
                # apical
                if key[2][0]==0:
                    figdf.at[key, 'trace_color']=self.black
                    figdf.at[key, 'trace_ecolor']=self.gray
                # basal
                if key[2][0]>0:
                    figdf.at[key, 'trace_color']=self.red
                    figdf.at[key, 'trace_ecolor']=self.red_light
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

    def trace_mean_v(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            segs = {
            'proximal': (0,1,2,3,4,5),
            'distal':(6,7,8,9,10)
            }
            distal_weights = np.append(10.*np.round(np.arange(3., 4., .1),5), 0.)
            proximal_weights = np.append(10.*np.round(np.arange(1.4, 2.4, .1), 5), 0.)
            figdict={}
            for figloc in segs.keys():
                for proximal_weight in proximal_weights:
                    for distal_weight in distal_weights:
                        figkey = 'v_'+figloc+'_proxweight_'+str(proximal_weight)+'_distweight_'+str(distal_weight)
                        figdict[figkey]={
                        'all':[
                            (5, 'apic', proximal_weight,distal_weight, segs[figloc]),
                            (0, 'apic',proximal_weight, distal_weight, segs[figloc]),
                            (-5, 'apic', proximal_weight, distal_weight, segs[figloc]),
                        ]
                        }

            # figdict = {
            #     # # figure
            #     # #-------------------------------------------------
            #     'proximal_14':{
            #         # subgroup
            #         'all':[
            #             # trace
            #             (5, 'dend',(1.4, 1.5)),
            #             (0, 'dend',(1.4, 1.5)),
            #             (-5, 'dend', (1.4, 1.5)),
            #         ]
            #     },
            # }

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            # figdf['fig_nyticks']=5
            # figdf['fig_nxticks']=10
            figdf['fig_dyticks']=10
            # figdf['fig_dxticks']=500.
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            figdf['fig_ymin']=-75
            figdf['fig_ymax']=-20
            figdf['fig_xmin']=200
            figdf['fig_xmax']=401.
            figdf['fig_ylabel']='Vm (mV)'
            figdf['fig_xlabel']='Time (ms)'
            # figdf['fig_dyticks']=.2
            figdf['fig_dxticks']=50
            # # trace level parameters
            figdf['trace_ealpha']=1
            figdf['error_style']='shade'
            figdf['trace_linewidth']=4
            figdf['trace_linestyle']='-'
            # figdf.at[slice(None), 'trace_ecolor'] = gray
            figdf['fig_xscale_log']=False
            figdf['fig_barwidth']=0.8
            figdf['fig_data_style']='point'
            figdf['fig_xtick_decimals']=0
            figdf['fig_ytick_decimals']=0
            figdf['fig_set_xscale']='linear'
            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]<0:
                    figdf.at[key, 'trace_color']=self.blue
                    figdf.at[key, 'trace_ecolor']=self.blue_light
                # apical
                if key[2][0]==0:
                    figdf.at[key, 'trace_color']=self.black
                    figdf.at[key, 'trace_ecolor']=self.gray
                # basal
                if key[2][0]>0:
                    figdf.at[key, 'trace_color']=self.red
                    figdf.at[key, 'trace_ecolor']=self.red_light
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

    def trace_individual_v(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            segs = {
            'proximal': (0,1,2,3,4,5),
            'distal':(6,7,8,9,10)
            }
            if 'trial_ids' in kwargs:
                trial_id = kwargs['trial_ids']
            path_offsets_list = [0, 20, 40, 60, 80]
            proximal_weights = [0., 17]
            distal_weights = [0., 32]
            figdict={}
            for figloc in segs.keys():
                for proximal_weight in proximal_weights:
                    for distal_weight in distal_weights:
                        for offset in path_offsets_list:
                            # for trial_id in trial_ids:
                            figkey_list = ['individual','v',figloc,'offset', str(offset), 'proxweight',str(proximal_weight),'distweight',str(distal_weight)]
                            figkey = '_'.join(figkey_list)
                            figdict[figkey]={
                            'all':[
                                (5, 'apic', offset, proximal_weight,distal_weight, segs[figloc]),
                                (0, 'apic',offset, proximal_weight, distal_weight, segs[figloc],),
                                (-5, 'apic', offset, proximal_weight, distal_weight, segs[figloc],),
                            ]
                            }

            # figdict = {
            #     # # figure
            #     # #-------------------------------------------------
            #     'proximal_14':{
            #         # subgroup
            #         'all':[
            #             # trace
            #             (5, 'dend',(1.4, 1.5)),
            #             (0, 'dend',(1.4, 1.5)),
            #             (-5, 'dend', (1.4, 1.5)),
            #         ]
            #     },
            # }

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            # figdf['fig_nyticks']=5
            # figdf['fig_nxticks']=10
            figdf['fig_dyticks']=10
            # figdf['fig_dxticks']=500.
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            figdf['fig_ymin']=-75
            figdf['fig_ymax']=-20
            figdf['fig_xmin']=400
            figdf['fig_xmax']=601.
            figdf['fig_ylabel']='Vm (mV)'
            figdf['fig_xlabel']='Time (ms)'
            # figdf['fig_dyticks']=.2
            figdf['fig_dxticks']=50
            # # trace level parameters
            figdf['trace_ealpha']=1
            figdf['error_style']='shade'
            figdf['trace_linewidth']=4
            figdf['trace_linestyle']='-'
            # figdf.at[slice(None), 'trace_ecolor'] = gray
            figdf['fig_xscale_log']=False
            figdf['fig_barwidth']=0.8
            figdf['fig_data_style']='point'
            figdf['fig_xtick_decimals']=0
            figdf['fig_ytick_decimals']=0
            figdf['fig_set_xscale']='linear'
            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]<0:
                    figdf.at[key, 'trace_color']=self.blue
                    figdf.at[key, 'trace_ecolor']=self.blue_light
                # apical
                if key[2][0]==0:
                    figdf.at[key, 'trace_color']=self.black
                    figdf.at[key, 'trace_ecolor']=self.gray
                # basal
                if key[2][0]>0:
                    figdf.at[key, 'trace_color']=self.red
                    figdf.at[key, 'trace_ecolor']=self.red_light
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

    def var2var_path_1_w_mean_dw_clopath(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            segs = {
            'proximal': (0,1,2,3,4,5),
            'distal':(6,7,8,9,10)
            }
            weights = 10.*np.round(np.arange(3., 4., .1),5)
            figdict={}
            for figloc in segs.keys():
                for weight in weights:
                    figkey = figloc+'_'+str(weight)
                    figdict[figkey]={
                    'all':[
                        (5, 'apic', weight, segs[figloc]),
                        (0, 'apic',weight, segs[figloc]),
                        (-5, 'apic', weight, segs[figloc]),
                    ]
                    }


            # figds

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            # figdf['fig_nyticks']=5
            # figdf['fig_nxticks']=10
            figdf['fig_dyticks']=1.
            figdf['fig_dxticks']=10.*0.5
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            # figdf['fig_ymin']=-2
            # figdf['fig_ymax']=2.1
            figdf['fig_xmin']=10.*1.4
            figdf['fig_xmax']=10*2.41
            figdf['fig_ylabel']='Norm. weight'
            figdf['fig_xlabel']='Initial mean weight (nS)'
            # figdf['fig_dyticks']=.2
            # figdf['fig_dxticks']=20
            # # trace level parameters
            figdf['trace_ealpha']=1
            figdf['error_style']='shade'
            figdf['trace_linewidth']=4
            # figdf.at[slice(None), 'trace_ecolor'] = gray
            figdf['fig_xscale_log']=False
            figdf['fig_barwidth']=0.8
            figdf['fig_data_style']='point'
            figdf['fig_xtick_decimals']=0
            figdf['fig_ytick_decimals']=0
            figdf['fig_set_xscale']='linear'
            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]<0:
                    figdf.at[key, 'trace_color']=self.blue
                    figdf.at[key, 'trace_ecolor']=self.blue_light
                # apical
                if key[2][0]==0:
                    figdf.at[key, 'trace_color']=self.black
                    figdf.at[key, 'trace_ecolor']=self.gray
                # basal
                if key[2][0]>0:
                    figdf.at[key, 'trace_color']=self.red
                    figdf.at[key, 'trace_ecolor']=self.red_light
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

    def var2var_offset_dw_clopath(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            segs = {
            'proximal': (0,1,2,3,4,5),
            'distal':(6,7,8,9,10)
            }
            # weights = 10.*np.round(np.arange(3., 4., .1),5)
            path_offsets_list = [0, 20, 40, 60, 80]
            proximal_weights = [0., 17]
            distal_weights = [0., 32]
            figdict={}
            for figloc in segs.keys():
                for proximal_weight in proximal_weights:
                    for distal_weight in distal_weights:
                        # for offset in path_offsets_list:
                        figkey = 'offset_'+figloc+'_proxweight_'+str(proximal_weight)+'_distweight_'+str(distal_weight)

                        figdict[figkey]={
                        'all':[
                            (5, 'apic',proximal_weight,distal_weight, segs[figloc]),
                            (0, 'apic',proximal_weight,distal_weight, segs[figloc]),
                            (-5, 'apic', proximal_weight,distal_weight, segs[figloc]),
                        ]
                        }


            # figds

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            # figdf['fig_nyticks']=5
            # figdf['fig_nxticks']=10
            figdf['fig_dyticks']=1.
            figdf['fig_dxticks']=20
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            # figdf['fig_ymin']=-2
            # figdf['fig_ymax']=2.1
            figdf['fig_xmin']=0
            figdf['fig_xmax']=81
            figdf['fig_ylabel']='Norm. weight'
            figdf['fig_xlabel']='Delay (ms)'
            # figdf['fig_dyticks']=.2
            # figdf['fig_dxticks']=20
            # # trace level parameters
            figdf['trace_ealpha']=1
            figdf['error_style']='shade'
            figdf['trace_linewidth']=4
            # figdf.at[slice(None), 'trace_ecolor'] = gray
            figdf['fig_xscale_log']=False
            figdf['fig_barwidth']=0.8
            figdf['fig_data_style']='point'
            figdf['fig_xtick_decimals']=0
            figdf['fig_ytick_decimals']=0
            figdf['fig_set_xscale']='linear'
            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]<0:
                    figdf.at[key, 'trace_color']=self.blue
                    figdf.at[key, 'trace_ecolor']=self.blue_light
                # apical
                if key[2][0]==0:
                    figdf.at[key, 'trace_color']=self.black
                    figdf.at[key, 'trace_ecolor']=self.gray
                # basal
                if key[2][0]>0:
                    figdf.at[key, 'trace_color']=self.red
                    figdf.at[key, 'trace_ecolor']=self.red_light
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

class figsetup_exp_reduced_neuron_tbs_conjunctive(FigSetup):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        super(figsetup_exp_reduced_neuron_tbs_conjunctive, self).__init__(**kwargs)

    def bar_dw_clopath(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            figdict = {
                # # figure
                # #-------------------------------------------------
                'all':{
                    # subgroup
                    'all':[
                        # trace
                        (5, 'dend',),
                        (0, 'dend',),
                        (-5, 'dend'),
                    ]
                },
            }

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            # figdf['fig_nytickss']=4
            # figdf['fig_nxticks']=10
            figdf['fig_dyticks']=1.
            # figdf['fig_dxticks']=10
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            figdf['fig_ymin']=1.
            figdf['fig_ymax']=7.1
            figdf['fig_xmin']=0.
            # figdf['fig_xmax']=30.
            figdf['fig_ylabel']='Norm. weight'
            # figdf['fig_xlabel']='Time (ms)'
            # # trace level parameters
            # figdf['trace_ealpha']=.7
            # figdf['error_style']='shade'
            # figdf['trace_linewidth']=4
            # figdf['fig_xscale']=1./40
            figdf['fig_barwidth']=1
            # figdf['fig_data_style']='point'
            # figdf['fig_xtick_decimals']=0
            figdf['fig_ytick_decimals']=0
            # figdf['fig_set_xscale']='symlog'

            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # set subgroup parameters
            #-----------------------------------
            figdf = figdf.reset_index().set_index('subgroup')
            figdf.at['all','sub_location']=0

            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]<0:
                    figdf.at[key, 'trace_color']=self.blue
                    figdf.at[key, 'trace_ecolor']=self.blue_light
                    figdf.at[key, 'trace_location']=0
                    figdf.at[key, 'trace_label']='cathodal'
                # apical
                if key[2][0]==0:
                    figdf.at[key, 'trace_color']=self.black
                    figdf.at[key, 'trace_ecolor']=self.gray
                    figdf.at[key, 'trace_location']=1
                    figdf.at[key, 'trace_label']='control'
                # basal
                if key[2][0]>0:
                    figdf.at[key, 'trace_color']=self.red
                    figdf.at[key, 'trace_ecolor']=self.red_light
                    figdf.at[key, 'trace_location']=2
                    figdf.at[key, 'trace_label']='anodal'
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

    def trace_mean_w_clopath(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            segs = {
            'proximal': (0,1,2,3,4,5),
            'distal':(6,7,8,9,10)
            }
            distal_weights = 10.*np.round(np.arange(3., 4., .1),5)
            proximal_weights = 10.*np.round(np.arange(1.4, 2.4, .1), 5)
            figdict={}
            for figloc in segs.keys():
                for proximal_weight in proximal_weights:
                    for distal_weight in distal_weights:
                        figkey = 'w_'+figloc+'_proxweight_'+str(proximal_weight)+'_distweight_'+str(distal_weight)
                        figdict[figkey]={
                        'all':[
                            (5, 'apic', proximal_weight,distal_weight, segs[figloc]),
                            (0, 'apic',proximal_weight, distal_weight, segs[figloc]),
                            (-5, 'apic', proximal_weight, distal_weight, segs[figloc]),
                        ]
                        }

            # figdict = {
            #     # # figure
            #     # #-------------------------------------------------
            #     'proximal_14':{
            #         # subgroup
            #         'all':[
            #             # trace
            #             (5, 'dend',(1.4, 1.5)),
            #             (0, 'dend',(1.4, 1.5)),
            #             (-5, 'dend', (1.4, 1.5)),
            #         ]
            #     },
            # }

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            # figdf['fig_nyticks']=5
            # figdf['fig_nxticks']=10
            # figdf['fig_dyticks']=.5
            # figdf['fig_dxticks']=500.
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            # figdf['fig_ymin']=0.5
            # figdf['fig_ymax']=3.1
            figdf['fig_xmin']=200
            figdf['fig_xmax']=401.
            figdf['fig_ylabel']='Weight'
            figdf['fig_xlabel']='Time (ms)'
            # figdf['fig_dyticks']=.2
            figdf['fig_dxticks']=50
            # # trace level parameters
            figdf['trace_ealpha']=1
            figdf['error_style']='shade'
            figdf['trace_linewidth']=4
            figdf['trace_linestyle']='-'
            # figdf.at[slice(None), 'trace_ecolor'] = gray
            figdf['fig_xscale_log']=False
            figdf['fig_barwidth']=0.8
            figdf['fig_data_style']='point'
            figdf['fig_xtick_decimals']=0
            figdf['fig_ytick_decimals']=1
            figdf['fig_set_xscale']='linear'
            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]<0:
                    figdf.at[key, 'trace_color']=self.blue
                    figdf.at[key, 'trace_ecolor']=self.blue_light
                # apical
                if key[2][0]==0:
                    figdf.at[key, 'trace_color']=self.black
                    figdf.at[key, 'trace_ecolor']=self.gray
                # basal
                if key[2][0]>0:
                    figdf.at[key, 'trace_color']=self.red
                    figdf.at[key, 'trace_ecolor']=self.red_light
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

    def trace_mean_v(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            segs = {
            'proximal': (0,1,2,3,4,5),
            'distal':(6,7,8,9,10)
            }
            distal_weights = np.append(10.*np.round(np.arange(3., 4., .1),5), 0.)
            proximal_weights = np.append(10.*np.round(np.arange(1.4, 2.4, .1), 5), 0.)
            figdict={}
            for figloc in segs.keys():
                for proximal_weight in proximal_weights:
                    for distal_weight in distal_weights:
                        figkey = 'v_'+figloc+'_proxweight_'+str(proximal_weight)+'_distweight_'+str(distal_weight)
                        figdict[figkey]={
                        'all':[
                            (5, 'apic', proximal_weight,distal_weight, segs[figloc]),
                            (0, 'apic',proximal_weight, distal_weight, segs[figloc]),
                            (-5, 'apic', proximal_weight, distal_weight, segs[figloc]),
                        ]
                        }

            # figdict = {
            #     # # figure
            #     # #-------------------------------------------------
            #     'proximal_14':{
            #         # subgroup
            #         'all':[
            #             # trace
            #             (5, 'dend',(1.4, 1.5)),
            #             (0, 'dend',(1.4, 1.5)),
            #             (-5, 'dend', (1.4, 1.5)),
            #         ]
            #     },
            # }

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            # figdf['fig_nyticks']=5
            # figdf['fig_nxticks']=10
            figdf['fig_dyticks']=10
            # figdf['fig_dxticks']=500.
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            figdf['fig_ymin']=-75
            figdf['fig_ymax']=-20
            figdf['fig_xmin']=200
            figdf['fig_xmax']=401.
            figdf['fig_ylabel']='Vm (mV)'
            figdf['fig_xlabel']='Time (ms)'
            # figdf['fig_dyticks']=.2
            figdf['fig_dxticks']=50
            # # trace level parameters
            figdf['trace_ealpha']=1
            figdf['error_style']='shade'
            figdf['trace_linewidth']=4
            figdf['trace_linestyle']='-'
            # figdf.at[slice(None), 'trace_ecolor'] = gray
            figdf['fig_xscale_log']=False
            figdf['fig_barwidth']=0.8
            figdf['fig_data_style']='point'
            figdf['fig_xtick_decimals']=0
            figdf['fig_ytick_decimals']=0
            figdf['fig_set_xscale']='linear'
            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]<0:
                    figdf.at[key, 'trace_color']=self.blue
                    figdf.at[key, 'trace_ecolor']=self.blue_light
                # apical
                if key[2][0]==0:
                    figdf.at[key, 'trace_color']=self.black
                    figdf.at[key, 'trace_ecolor']=self.gray
                # basal
                if key[2][0]>0:
                    figdf.at[key, 'trace_color']=self.red
                    figdf.at[key, 'trace_ecolor']=self.red_light
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

    def trace_individual_v(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            segs = {
            'proximal': (0,1,2,3,4,5),
            'distal':(6,7,8,9,10)
            }
            if 'trial_ids' in kwargs:
                trial_id = kwargs['trial_ids']
            path_offsets_list = [0, 20, 40, 60, 80]
            proximal_weights = [0., 15]
            distal_weights = [0., 32]
            figdict={}
            for figloc in segs.keys():
                for proximal_weight in proximal_weights:
                    for distal_weight in distal_weights:
                        for offset in path_offsets_list:
                            # for trial_id in trial_ids:
                            figkey_list = ['individual','v',figloc,'offset', str(offset), 'proxweight',str(proximal_weight),'distweight',str(distal_weight)]
                            figkey = '_'.join(figkey_list)
                            figdict[figkey]={
                            'all':[
                                (5, 'apic', offset, proximal_weight,distal_weight, segs[figloc]),
                                (0, 'apic',offset, proximal_weight, distal_weight, segs[figloc],),
                                (-5, 'apic', offset, proximal_weight, distal_weight, segs[figloc],),
                            ]
                            }

            # figdict = {
            #     # # figure
            #     # #-------------------------------------------------
            #     'proximal_14':{
            #         # subgroup
            #         'all':[
            #             # trace
            #             (5, 'dend',(1.4, 1.5)),
            #             (0, 'dend',(1.4, 1.5)),
            #             (-5, 'dend', (1.4, 1.5)),
            #         ]
            #     },
            # }

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            # figdf['fig_nyticks']=5
            # figdf['fig_nxticks']=10
            figdf['fig_dyticks']=10
            # figdf['fig_dxticks']=500.
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            figdf['fig_ymin']=-75
            figdf['fig_ymax']=-20
            figdf['fig_xmin']=400
            figdf['fig_xmax']=601.
            figdf['fig_ylabel']='Vm (mV)'
            figdf['fig_xlabel']='Time (ms)'
            # figdf['fig_dyticks']=.2
            figdf['fig_dxticks']=50
            # # trace level parameters
            figdf['trace_ealpha']=1
            figdf['error_style']='shade'
            figdf['trace_linewidth']=4
            figdf['trace_linestyle']='-'
            # figdf.at[slice(None), 'trace_ecolor'] = gray
            figdf['fig_xscale_log']=False
            figdf['fig_barwidth']=0.8
            figdf['fig_data_style']='point'
            figdf['fig_xtick_decimals']=0
            figdf['fig_ytick_decimals']=0
            figdf['fig_set_xscale']='linear'
            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]<0:
                    figdf.at[key, 'trace_color']=self.blue
                    figdf.at[key, 'trace_ecolor']=self.blue_light
                # apical
                if key[2][0]==0:
                    figdf.at[key, 'trace_color']=self.black
                    figdf.at[key, 'trace_ecolor']=self.gray
                # basal
                if key[2][0]>0:
                    figdf.at[key, 'trace_color']=self.red
                    figdf.at[key, 'trace_ecolor']=self.red_light
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

    def var2var_path_1_w_mean_dw_clopath(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            segs = {
            'proximal': (0,1,2,3,4,5),
            'distal':(6,7,8,9,10)
            }
            weights = 10.*np.round(np.arange(3., 4., .1),5)
            figdict={}
            for figloc in segs.keys():
                for weight in weights:
                    figkey = figloc+'_'+str(weight)
                    figdict[figkey]={
                    'all':[
                        (5, 'apic', weight, segs[figloc]),
                        (0, 'apic',weight, segs[figloc]),
                        (-5, 'apic', weight, segs[figloc]),
                    ]
                    }


            # figds

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            # figdf['fig_nyticks']=5
            # figdf['fig_nxticks']=10
            figdf['fig_dyticks']=1.
            figdf['fig_dxticks']=10.*0.5
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            # figdf['fig_ymin']=-2
            # figdf['fig_ymax']=2.1
            figdf['fig_xmin']=10.*1.4
            figdf['fig_xmax']=10*2.41
            figdf['fig_ylabel']='Norm. weight'
            figdf['fig_xlabel']='Initial mean weight (nS)'
            # figdf['fig_dyticks']=.2
            # figdf['fig_dxticks']=20
            # # trace level parameters
            figdf['trace_ealpha']=1
            figdf['error_style']='shade'
            figdf['trace_linewidth']=4
            # figdf.at[slice(None), 'trace_ecolor'] = gray
            figdf['fig_xscale_log']=False
            figdf['fig_barwidth']=0.8
            figdf['fig_data_style']='point'
            figdf['fig_xtick_decimals']=0
            figdf['fig_ytick_decimals']=0
            figdf['fig_set_xscale']='linear'
            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]<0:
                    figdf.at[key, 'trace_color']=self.blue
                    figdf.at[key, 'trace_ecolor']=self.blue_light
                # apical
                if key[2][0]==0:
                    figdf.at[key, 'trace_color']=self.black
                    figdf.at[key, 'trace_ecolor']=self.gray
                # basal
                if key[2][0]>0:
                    figdf.at[key, 'trace_color']=self.red
                    figdf.at[key, 'trace_ecolor']=self.red_light
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

    def var2var_offset_dw_clopath(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            segs = {
            'proximal': (0,1,2,3,4,5),
            'distal':(6,7,8,9,10)
            }
            # weights = 10.*np.round(np.arange(3., 4., .1),5)
            path_offsets_list = [0, 20, 40, 60, 80]
            proximal_weights = [0., 15]
            distal_weights = [0., 32]
            figdict={}
            for figloc in segs.keys():
                for proximal_weight in proximal_weights:
                    for distal_weight in distal_weights:
                        # for offset in path_offsets_list:
                        figkey = 'offset_'+figloc+'_proxweight_'+str(proximal_weight)+'_distweight_'+str(distal_weight)

                        figdict[figkey]={
                        'all':[
                            (5, 'apic',proximal_weight,distal_weight, segs[figloc]),
                            (0, 'apic',proximal_weight,distal_weight, segs[figloc]),
                            (-5, 'apic', proximal_weight,distal_weight, segs[figloc]),
                        ]
                        }


            # figds

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            # figdf['fig_nyticks']=5
            # figdf['fig_nxticks']=10
            figdf['fig_dyticks']=1.
            figdf['fig_dxticks']=20
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            # figdf['fig_ymin']=-2
            # figdf['fig_ymax']=2.1
            figdf['fig_xmin']=0
            figdf['fig_xmax']=81
            figdf['fig_ylabel']='Norm. weight'
            figdf['fig_xlabel']='Delay (ms)'
            # figdf['fig_dyticks']=.2
            # figdf['fig_dxticks']=20
            # # trace level parameters
            figdf['trace_ealpha']=1
            figdf['error_style']='shade'
            figdf['trace_linewidth']=4
            # figdf.at[slice(None), 'trace_ecolor'] = gray
            figdf['fig_xscale_log']=False
            figdf['fig_barwidth']=0.8
            figdf['fig_data_style']='point'
            figdf['fig_xtick_decimals']=0
            figdf['fig_ytick_decimals']=0
            figdf['fig_set_xscale']='linear'
            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]<0:
                    figdf.at[key, 'trace_color']=self.blue
                    figdf.at[key, 'trace_ecolor']=self.blue_light
                # apical
                if key[2][0]==0:
                    figdf.at[key, 'trace_color']=self.black
                    figdf.at[key, 'trace_ecolor']=self.gray
                # basal
                if key[2][0]>0:
                    figdf.at[key, 'trace_color']=self.red
                    figdf.at[key, 'trace_ecolor']=self.red_light
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

class figsetup_exp_reduced_neuron_1hz(FigSetup):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        super(figsetup_exp_reduced_neuron_1hz, self).__init__(**kwargs)

    def bar_dw_clopath(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            figdict = {
                # # figure
                # #-------------------------------------------------
                'all':{
                    # subgroup
                    'all':[
                        # trace
                        (5, 'apic',),
                        (0, 'apic',),
                        (-5, 'apic'),
                    ]
                },
            }

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            figdf['fig_nyticks']=4
            # figdf['fig_nxticks']=10
            figdf['fig_dyticks']=.2
            # figdf['fig_dxticks']=10
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            figdf['fig_ymin']=.8
            figdf['fig_ymax']=1.01
            figdf['fig_xmin']=0.
            # figdf['fig_xmax']=30.
            figdf['fig_ylabel']='Norm. weight'
            # figdf['fig_xlabel']='Time (ms)'
            # # trace level parameters
            # figdf['trace_ealpha']=.7
            # figdf['error_style']='shade'
            # figdf['trace_linewidth']=4
            # figdf['fig_xscale']=1./40
            figdf['fig_barwidth']=1
            # figdf['fig_data_style']='point'
            # figdf['fig_xtick_decimals']=0
            figdf['fig_ytick_decimals']=2
            # figdf['fig_set_xscale']='symlog'

            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # set subgroup parameters
            #-----------------------------------
            figdf = figdf.reset_index().set_index('subgroup')
            figdf.at['all','sub_location']=0

            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]<0:
                    figdf.at[key, 'trace_color']=self.blue
                    figdf.at[key, 'trace_ecolor']=self.blue_light
                    figdf.at[key, 'trace_location']=0
                    figdf.at[key, 'trace_label']='cathodal'
                # apical
                if key[2][0]==0:
                    figdf.at[key, 'trace_color']=self.black
                    figdf.at[key, 'trace_ecolor']=self.gray
                    figdf.at[key, 'trace_location']=1
                    figdf.at[key, 'trace_label']='control'
                # basal
                if key[2][0]>0:
                    figdf.at[key, 'trace_color']=self.red
                    figdf.at[key, 'trace_ecolor']=self.red_light
                    figdf.at[key, 'trace_location']=2
                    figdf.at[key, 'trace_label']='anodal'
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

    def trace_mean_w_clopath(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            figdict = {
                # # figure
                # #-------------------------------------------------
                'all':{
                    # subgroup
                    'all':[
                        # trace
                        (5, 'apic',),
                        (0, 'apic',),
                        (-5, 'apic'),
                    ]
                },
            }

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            # figdf['fig_nyticks']=5
            # figdf['fig_nxticks']=10
            figdf['fig_dyticks']=.1
            figdf['fig_dxticks']=5
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            figdf['fig_ymin']=0.4
            figdf['fig_ymax']=.61
            figdf['fig_xmin']=0
            figdf['fig_xmax']=13
            figdf['fig_ylabel']='Weight'
            figdf['fig_xlabel']='Time (ss)'
            # figdf['fig_dyticks']=.2
            # figdf['fig_dxticks']=20
            # # trace level parameters
            figdf['trace_ealpha']=1
            figdf['error_style']='shade'
            figdf['trace_linewidth']=4
            figdf['trace_linestyle']='-'
            # figdf.at[slice(None), 'trace_ecolor'] = gray
            figdf['fig_xscale_log']=False
            figdf['fig_barwidth']=0.8
            figdf['fig_data_style']='point'
            figdf['fig_xtick_decimals']=0
            figdf['fig_ytick_decimals']=1
            figdf['fig_set_xscale']='linear'
            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]<0:
                    figdf.at[key, 'trace_color']=self.blue
                    figdf.at[key, 'trace_ecolor']=self.blue_light
                # apical
                if key[2][0]==0:
                    figdf.at[key, 'trace_color']=self.black
                    figdf.at[key, 'trace_ecolor']=self.gray
                # basal
                if key[2][0]>0:
                    figdf.at[key, 'trace_color']=self.red
                    figdf.at[key, 'trace_ecolor']=self.red_light
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

    def var2var_path_1_w_mean_dw_clopath(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            figdict = {
                # # figure
                # #-------------------------------------------------
                'all':{
                    # subgroup
                    'all':[
                        # trace
                        (5, 'apic',),
                        (0, 'apic',),
                        (-5, 'apic'),
                    ]
                },
            }

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            # figdf['fig_nyticks']=5
            # figdf['fig_nxticks']=10
            figdf['fig_dyticks']=.2
            figdf['fig_dxticks']=.1
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            # figdf['fig_ymin']=-2
            # figdf['fig_ymax']=2.1
            figdf['fig_xmin']=5.8
            figdf['fig_xmax']=6.21
            figdf['fig_ylabel']='Normalized weight'
            figdf['fig_xlabel']='Initial mean weight (nS)'
            # figdf['fig_dyticks']=.2
            # figdf['fig_dxticks']=20
            # # trace level parameters
            figdf['trace_ealpha']=1
            figdf['error_style']='shade'
            figdf['trace_linewidth']=4
            # figdf.at[slice(None), 'trace_ecolor'] = gray
            figdf['fig_xscale_log']=False
            figdf['fig_barwidth']=0.8
            figdf['fig_data_style']='point'
            figdf['fig_xtick_decimals']=1
            figdf['fig_ytick_decimals']=1
            figdf['fig_set_xscale']='linear'
            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]<0:
                    figdf.at[key, 'trace_color']=self.blue
                    figdf.at[key, 'trace_ecolor']=self.blue_light
                # apical
                if key[2][0]==0:
                    figdf.at[key, 'trace_color']=self.black
                    figdf.at[key, 'trace_ecolor']=self.gray
                # basal
                if key[2][0]>0:
                    figdf.at[key, 'trace_color']=self.red
                    figdf.at[key, 'trace_ecolor']=self.red_light
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

class figsetup_exp_reduced_neuron_tbs_basal(FigSetup):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        super(figsetup_exp_reduced_neuron_tbs_basal, self).__init__(**kwargs)

    def bar_dw_clopath(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            figdict = {
                # # figure
                # #-------------------------------------------------
                'all':{
                    # subgroup
                    'all':[
                        # trace
                        (5, 'dend',),
                        (0, 'dend',),
                        (-5, 'dend'),
                    ]
                },
            }

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            # figdf['fig_nytickss']=4
            # figdf['fig_nxticks']=10
            figdf['fig_dyticks']=1.
            # figdf['fig_dxticks']=10
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            figdf['fig_ymin']=1.
            figdf['fig_ymax']=7.1
            figdf['fig_xmin']=0.
            # figdf['fig_xmax']=30.
            figdf['fig_ylabel']='Norm. weight'
            # figdf['fig_xlabel']='Time (ms)'
            # # trace level parameters
            # figdf['trace_ealpha']=.7
            # figdf['error_style']='shade'
            # figdf['trace_linewidth']=4
            # figdf['fig_xscale']=1./40
            figdf['fig_barwidth']=1
            # figdf['fig_data_style']='point'
            # figdf['fig_xtick_decimals']=0
            figdf['fig_ytick_decimals']=0
            # figdf['fig_set_xscale']='symlog'

            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # set subgroup parameters
            #-----------------------------------
            figdf = figdf.reset_index().set_index('subgroup')
            figdf.at['all','sub_location']=0

            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]<0:
                    figdf.at[key, 'trace_color']=self.blue
                    figdf.at[key, 'trace_ecolor']=self.blue_light
                    figdf.at[key, 'trace_location']=0
                    figdf.at[key, 'trace_label']='cathodal'
                # apical
                if key[2][0]==0:
                    figdf.at[key, 'trace_color']=self.black
                    figdf.at[key, 'trace_ecolor']=self.gray
                    figdf.at[key, 'trace_location']=1
                    figdf.at[key, 'trace_label']='control'
                # basal
                if key[2][0]>0:
                    figdf.at[key, 'trace_color']=self.red
                    figdf.at[key, 'trace_ecolor']=self.red_light
                    figdf.at[key, 'trace_location']=2
                    figdf.at[key, 'trace_label']='anodal'
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

    def trace_mean_w_clopath(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            figdict = {
                # # figure
                # #-------------------------------------------------
                'all':{
                    # subgroup
                    'all':[
                        # trace
                        (5, 'dend',),
                        (0, 'dend',),
                        (-5, 'dend'),
                    ]
                },
            }

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            # figdf['fig_nyticks']=5
            # figdf['fig_nxticks']=10
            figdf['fig_dyticks']=.5
            figdf['fig_dxticks']=500.
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            figdf['fig_ymin']=0.5
            figdf['fig_ymax']=3.1
            figdf['fig_xmin']=0
            figdf['fig_xmax']=3001.
            figdf['fig_ylabel']='Weight'
            figdf['fig_xlabel']='Time (ms)'
            # figdf['fig_dyticks']=.2
            # figdf['fig_dxticks']=20
            # # trace level parameters
            figdf['trace_ealpha']=1
            figdf['error_style']='shade'
            figdf['trace_linewidth']=4
            figdf['trace_linestyle']='-'
            # figdf.at[slice(None), 'trace_ecolor'] = gray
            figdf['fig_xscale_log']=False
            figdf['fig_barwidth']=0.8
            figdf['fig_data_style']='point'
            figdf['fig_xtick_decimals']=0
            figdf['fig_ytick_decimals']=1
            figdf['fig_set_xscale']='linear'
            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]<0:
                    figdf.at[key, 'trace_color']=self.blue
                    figdf.at[key, 'trace_ecolor']=self.blue_light
                # apical
                if key[2][0]==0:
                    figdf.at[key, 'trace_color']=self.black
                    figdf.at[key, 'trace_ecolor']=self.gray
                # basal
                if key[2][0]>0:
                    figdf.at[key, 'trace_color']=self.red
                    figdf.at[key, 'trace_ecolor']=self.red_light
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

    def var2var_path_1_w_mean_dw_clopath(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            figdict = {
                # # figure
                # #-------------------------------------------------
                'all':{
                    # subgroup
                    'all':[
                        # trace
                        (5, 'dend',),
                        (0, 'dend',),
                        (-5, 'dend'),
                    ]
                },
            }

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            # figdf['fig_nyticks']=5
            # figdf['fig_nxticks']=10
            figdf['fig_dyticks']=1.
            figdf['fig_dxticks']=0.1
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            # figdf['fig_ymin']=-2
            # figdf['fig_ymax']=2.1
            figdf['fig_xmin']=5.8
            figdf['fig_xmax']=6.21
            figdf['fig_ylabel']='Normalized weight'
            figdf['fig_xlabel']='Initial mean weight (nS)'
            # figdf['fig_dyticks']=.2
            # figdf['fig_dxticks']=20
            # # trace level parameters
            figdf['trace_ealpha']=1
            figdf['error_style']='shade'
            figdf['trace_linewidth']=4
            # figdf.at[slice(None), 'trace_ecolor'] = gray
            figdf['fig_xscale_log']=False
            figdf['fig_barwidth']=0.8
            figdf['fig_data_style']='point'
            figdf['fig_xtick_decimals']=1
            figdf['fig_ytick_decimals']=0
            figdf['fig_set_xscale']='linear'
            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]<0:
                    figdf.at[key, 'trace_color']=self.blue
                    figdf.at[key, 'trace_ecolor']=self.blue_light
                # apical
                if key[2][0]==0:
                    figdf.at[key, 'trace_color']=self.black
                    figdf.at[key, 'trace_ecolor']=self.gray
                # basal
                if key[2][0]>0:
                    figdf.at[key, 'trace_color']=self.red
                    figdf.at[key, 'trace_ecolor']=self.red_light
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

class figsetup_exp_reduced_neuron_20hz_basal(FigSetup):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        super(figsetup_exp_reduced_neuron_20hz_basal, self).__init__(**kwargs)

    def bar_dw_clopath(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            figdict = {
                # # figure
                # #-------------------------------------------------
                'all':{
                    # subgroup
                    'all':[
                        # trace
                        (5, 'dend',),
                        (0, 'dend',),
                        (-5, 'dend'),
                    ]
                },
            }

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            # figdf['fig_nyticks']=4
            # figdf['fig_nxticks']=10
            figdf['fig_dyticks']=1.
            # figdf['fig_dxticks']=10
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            figdf['fig_ymin']=1.
            figdf['fig_ymax']=5.1
            figdf['fig_xmin']=0.
            # figdf['fig_xmax']=30.
            figdf['fig_ylabel']='Norm. weight'
            # figdf['fig_xlabel']='Time (ms)'
            # # trace level parameters
            # figdf['trace_ealpha']=.7
            # figdf['error_style']='shade'
            # figdf['trace_linewidth']=4
            # figdf['fig_xscale']=1./40
            figdf['fig_barwidth']=1
            # figdf['fig_data_style']='point'
            # figdf['fig_xtick_decimals']=0
            figdf['fig_ytick_decimals']=0
            # figdf['fig_set_xscale']='symlog'

            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # set subgroup parameters
            #-----------------------------------
            figdf = figdf.reset_index().set_index('subgroup')
            figdf.at['all','sub_location']=0

            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]<0:
                    figdf.at[key, 'trace_color']=self.blue
                    figdf.at[key, 'trace_ecolor']=self.blue_light
                    figdf.at[key, 'trace_location']=0
                    figdf.at[key, 'trace_label']='cathodal'
                # apical
                if key[2][0]==0:
                    figdf.at[key, 'trace_color']=self.black
                    figdf.at[key, 'trace_ecolor']=self.gray
                    figdf.at[key, 'trace_location']=1
                    figdf.at[key, 'trace_label']='control'
                # basal
                if key[2][0]>0:
                    figdf.at[key, 'trace_color']=self.red
                    figdf.at[key, 'trace_ecolor']=self.red_light
                    figdf.at[key, 'trace_location']=2
                    figdf.at[key, 'trace_label']='anodal'
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

    def trace_mean_w_clopath(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            figdict = {
                # # figure
                # #-------------------------------------------------
                'all':{
                    # subgroup
                    'all':[
                        # trace
                        (5, 'dend',),
                        (0, 'dend',),
                        (-5, 'dend'),
                    ]
                },
            }

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            # figdf['fig_nyticks']=5
            # figdf['fig_nxticks']=10
            figdf['fig_dyticks']=.5
            figdf['fig_dxticks']=500.
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            figdf['fig_ymin']=0.5
            figdf['fig_ymax']=2.1
            figdf['fig_xmin']=0
            figdf['fig_xmax']=3001.
            figdf['fig_ylabel']='Weight'
            figdf['fig_xlabel']='Time (ms)'
            # figdf['fig_dyticks']=.2
            # figdf['fig_dxticks']=20
            # # trace level parameters
            figdf['trace_ealpha']=1
            figdf['error_style']='shade'
            figdf['trace_linewidth']=4
            figdf['trace_linestyle']='-'
            # figdf.at[slice(None), 'trace_ecolor'] = gray
            figdf['fig_xscale_log']=False
            figdf['fig_barwidth']=0.8
            figdf['fig_data_style']='point'
            figdf['fig_xtick_decimals']=0
            figdf['fig_ytick_decimals']=1
            figdf['fig_set_xscale']='linear'
            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]<0:
                    figdf.at[key, 'trace_color']=self.blue
                    figdf.at[key, 'trace_ecolor']=self.blue_light
                # apical
                if key[2][0]==0:
                    figdf.at[key, 'trace_color']=self.black
                    figdf.at[key, 'trace_ecolor']=self.gray
                # basal
                if key[2][0]>0:
                    figdf.at[key, 'trace_color']=self.red
                    figdf.at[key, 'trace_ecolor']=self.red_light
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

    def var2var_path_1_w_mean_dw_clopath(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            figdict = {
                # # figure
                # #-------------------------------------------------
                'all':{
                    # subgroup
                    'all':[
                        # trace
                        (5, 'dend',),
                        (0, 'dend',),
                        (-5, 'dend'),
                    ]
                },
            }

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            # figdf['fig_nyticks']=5
            # figdf['fig_nxticks']=10
            figdf['fig_dyticks']=1.
            figdf['fig_dxticks']=0.1
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            # figdf['fig_ymin']=-2
            # figdf['fig_ymax']=2.1
            figdf['fig_xmin']=5.8
            figdf['fig_xmax']=6.21
            figdf['fig_ylabel']='Normalized weight'
            figdf['fig_xlabel']='Initial mean weight (nS)'
            # figdf['fig_dyticks']=.2
            # figdf['fig_dxticks']=20
            # # trace level parameters
            figdf['trace_ealpha']=1
            figdf['error_style']='shade'
            figdf['trace_linewidth']=4
            # figdf.at[slice(None), 'trace_ecolor'] = gray
            figdf['fig_xscale_log']=False
            figdf['fig_barwidth']=0.8
            figdf['fig_data_style']='point'
            figdf['fig_xtick_decimals']=1
            figdf['fig_ytick_decimals']=0
            figdf['fig_set_xscale']='linear'
            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]<0:
                    figdf.at[key, 'trace_color']=self.blue
                    figdf.at[key, 'trace_ecolor']=self.blue_light
                # apical
                if key[2][0]==0:
                    figdf.at[key, 'trace_color']=self.black
                    figdf.at[key, 'trace_ecolor']=self.gray
                # basal
                if key[2][0]>0:
                    figdf.at[key, 'trace_color']=self.red
                    figdf.at[key, 'trace_ecolor']=self.red_light
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

class figsetup_exp_reduced_neuron_20hz(FigSetup):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        super(figsetup_exp_reduced_neuron_20hz, self).__init__(**kwargs)

    def bar_dw_clopath(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            figdict = {
                # # figure
                # #-------------------------------------------------
                'all':{
                    # subgroup
                    'all':[
                        # trace
                        (5, 'apic',),
                        (0, 'apic',),
                        (-5, 'apic'),
                    ]
                },
            }

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            figdf['fig_nyticks']=4
            # figdf['fig_nxticks']=10
            figdf['fig_dyticks']=.05
            # figdf['fig_dxticks']=10
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            figdf['fig_ymin']=1.15
            figdf['fig_ymax']=1.26
            figdf['fig_xmin']=0.
            # figdf['fig_xmax']=30.
            figdf['fig_ylabel']='Norm. weight'
            # figdf['fig_xlabel']='Time (ms)'
            # # trace level parameters
            # figdf['trace_ealpha']=.7
            # figdf['error_style']='shade'
            # figdf['trace_linewidth']=4
            # figdf['fig_xscale']=1./40
            figdf['fig_barwidth']=1
            # figdf['fig_data_style']='point'
            # figdf['fig_xtick_decimals']=0
            figdf['fig_ytick_decimals']=2
            # figdf['fig_set_xscale']='symlog'

            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # set subgroup parameters
            #-----------------------------------
            figdf = figdf.reset_index().set_index('subgroup')
            figdf.at['all','sub_location']=0

            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]<0:
                    figdf.at[key, 'trace_color']=self.blue
                    figdf.at[key, 'trace_ecolor']=self.blue_light
                    figdf.at[key, 'trace_location']=0
                    figdf.at[key, 'trace_label']='cathodal'
                # apical
                if key[2][0]==0:
                    figdf.at[key, 'trace_color']=self.black
                    figdf.at[key, 'trace_ecolor']=self.gray
                    figdf.at[key, 'trace_location']=1
                    figdf.at[key, 'trace_label']='control'
                # basal
                if key[2][0]>0:
                    figdf.at[key, 'trace_color']=self.red
                    figdf.at[key, 'trace_ecolor']=self.red_light
                    figdf.at[key, 'trace_location']=2
                    figdf.at[key, 'trace_label']='anodal'
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

    def trace_mean_w_clopath(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            figdict = {
                # # figure
                # #-------------------------------------------------
                'all':{
                    # subgroup
                    'all':[
                        # trace
                        (5, 'apic',),
                        (0, 'apic',),
                        (-5, 'apic'),
                    ]
                },
            }

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            # figdf['fig_nyticks']=5
            # figdf['fig_nxticks']=10
            figdf['fig_dyticks']=.1
            figdf['fig_dxticks']=500.
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            figdf['fig_ymin']=0.5
            figdf['fig_ymax']=.61
            figdf['fig_xmin']=0
            figdf['fig_xmax']=3001.
            figdf['fig_ylabel']='Weight'
            figdf['fig_xlabel']='Time (ms)'
            # figdf['fig_dyticks']=.2
            # figdf['fig_dxticks']=20
            # # trace level parameters
            figdf['trace_ealpha']=1
            figdf['error_style']='shade'
            figdf['trace_linewidth']=4
            figdf['trace_linestyle']='-'
            # figdf.at[slice(None), 'trace_ecolor'] = gray
            figdf['fig_xscale_log']=False
            figdf['fig_barwidth']=0.8
            figdf['fig_data_style']='point'
            figdf['fig_xtick_decimals']=0
            figdf['fig_ytick_decimals']=1
            figdf['fig_set_xscale']='linear'
            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]<0:
                    figdf.at[key, 'trace_color']=self.blue
                    figdf.at[key, 'trace_ecolor']=self.blue_light
                # apical
                if key[2][0]==0:
                    figdf.at[key, 'trace_color']=self.black
                    figdf.at[key, 'trace_ecolor']=self.gray
                # basal
                if key[2][0]>0:
                    figdf.at[key, 'trace_color']=self.red
                    figdf.at[key, 'trace_ecolor']=self.red_light
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

    def var2var_path_1_w_mean_dw_clopath(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            figdict = {
                # # figure
                # #-------------------------------------------------
                'all':{
                    # subgroup
                    'all':[
                        # trace
                        (5, 'apic',),
                        (0, 'apic',),
                        (-5, 'apic'),
                    ]
                },
            }

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            # figdf['fig_nyticks']=5
            # figdf['fig_nxticks']=10
            figdf['fig_dyticks']=.2
            figdf['fig_dxticks']=.1
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            # figdf['fig_ymin']=-2
            # figdf['fig_ymax']=2.1
            figdf['fig_xmin']=5.8
            figdf['fig_xmax']=6.21
            figdf['fig_ylabel']='Normalized weight'
            figdf['fig_xlabel']='Initial mean weight (nS)'
            # figdf['fig_dyticks']=.2
            # figdf['fig_dxticks']=20
            # # trace level parameters
            figdf['trace_ealpha']=1
            figdf['error_style']='shade'
            figdf['trace_linewidth']=4
            # figdf.at[slice(None), 'trace_ecolor'] = gray
            figdf['fig_xscale_log']=False
            figdf['fig_barwidth']=0.8
            figdf['fig_data_style']='point'
            figdf['fig_xtick_decimals']=1
            figdf['fig_ytick_decimals']=1
            figdf['fig_set_xscale']='linear'
            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]<0:
                    figdf.at[key, 'trace_color']=self.blue
                    figdf.at[key, 'trace_ecolor']=self.blue_light
                # apical
                if key[2][0]==0:
                    figdf.at[key, 'trace_color']=self.black
                    figdf.at[key, 'trace_ecolor']=self.gray
                # basal
                if key[2][0]>0:
                    figdf.at[key, 'trace_color']=self.red
                    figdf.at[key, 'trace_ecolor']=self.red_light
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

class figsetup_exp_reduced_neuron_tbs(FigSetup):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        super(figsetup_exp_reduced_neuron_tbs, self).__init__(**kwargs)

    def bar_dw_clopath(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            figdict = {
                # # figure
                # #-------------------------------------------------
                'all':{
                    # subgroup
                    'all':[
                        # trace
                        (5, 'apic',),
                        (0, 'apic',),
                        (-5, 'apic'),
                    ]
                },
            }

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            figdf['fig_nyticks']=4
            # figdf['fig_nxticks']=10
            figdf['fig_dyticks']=.2
            # figdf['fig_dxticks']=10
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            figdf['fig_ymin']=1.5
            figdf['fig_ymax']=2.51
            figdf['fig_xmin']=0.
            # figdf['fig_xmax']=30.
            figdf['fig_ylabel']='Norm. weight'
            # figdf['fig_xlabel']='Time (ms)'
            # # trace level parameters
            # figdf['trace_ealpha']=.7
            # figdf['error_style']='shade'
            # figdf['trace_linewidth']=4
            # figdf['fig_xscale']=1./40
            figdf['fig_barwidth']=1
            # figdf['fig_data_style']='point'
            # figdf['fig_xtick_decimals']=0
            figdf['fig_ytick_decimals']=2
            # figdf['fig_set_xscale']='symlog'

            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # set subgroup parameters
            #-----------------------------------
            figdf = figdf.reset_index().set_index('subgroup')
            figdf.at['all','sub_location']=0

            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]<0:
                    figdf.at[key, 'trace_color']=self.blue
                    figdf.at[key, 'trace_ecolor']=self.blue_light
                    figdf.at[key, 'trace_location']=0
                    figdf.at[key, 'trace_label']='cathodal'
                # apical
                if key[2][0]==0:
                    figdf.at[key, 'trace_color']=self.black
                    figdf.at[key, 'trace_ecolor']=self.gray
                    figdf.at[key, 'trace_location']=1
                    figdf.at[key, 'trace_label']='control'
                # basal
                if key[2][0]>0:
                    figdf.at[key, 'trace_color']=self.red
                    figdf.at[key, 'trace_ecolor']=self.red_light
                    figdf.at[key, 'trace_location']=2
                    figdf.at[key, 'trace_label']='anodal'
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

    def trace_mean_w_clopath(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            figdict = {
                # # figure
                # #-------------------------------------------------
                'all':{
                    # subgroup
                    'all':[
                        # trace
                        (5, 'apic',),
                        (0, 'apic',),
                        (-5, 'apic'),
                    ]
                },
            }

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            # figdf['fig_nyticks']=5
            # figdf['fig_nxticks']=10
            figdf['fig_dyticks']=.5
            figdf['fig_dxticks']=500.
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            figdf['fig_ymin']=0.5
            figdf['fig_ymax']=1.51
            figdf['fig_xmin']=0
            figdf['fig_xmax']=3001.
            figdf['fig_ylabel']='Weight'
            figdf['fig_xlabel']='Time (ms)'
            # figdf['fig_dyticks']=.2
            # figdf['fig_dxticks']=20
            # # trace level parameters
            figdf['trace_ealpha']=1
            figdf['error_style']='shade'
            figdf['trace_linewidth']=4
            figdf['trace_linestyle']='-'
            # figdf.at[slice(None), 'trace_ecolor'] = gray
            figdf['fig_xscale_log']=False
            figdf['fig_barwidth']=0.8
            figdf['fig_data_style']='point'
            figdf['fig_xtick_decimals']=0
            figdf['fig_ytick_decimals']=1
            figdf['fig_set_xscale']='linear'
            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]<0:
                    figdf.at[key, 'trace_color']=self.blue
                    figdf.at[key, 'trace_ecolor']=self.blue_light
                # apical
                if key[2][0]==0:
                    figdf.at[key, 'trace_color']=self.black
                    figdf.at[key, 'trace_ecolor']=self.gray
                # basal
                if key[2][0]>0:
                    figdf.at[key, 'trace_color']=self.red
                    figdf.at[key, 'trace_ecolor']=self.red_light
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

    def var2var_path_1_w_mean_dw_clopath(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            figdict = {
                # # figure
                # #-------------------------------------------------
                'all':{
                    # subgroup
                    'all':[
                        # trace
                        (5, 'apic',),
                        (0, 'apic',),
                        (-5, 'apic'),
                    ]
                },
            }

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            # figdf['fig_nyticks']=5
            # figdf['fig_nxticks']=10
            figdf['fig_dyticks']=1.
            figdf['fig_dxticks']=0.1
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            # figdf['fig_ymin']=-2
            # figdf['fig_ymax']=2.1
            figdf['fig_xmin']=5.8
            figdf['fig_xmax']=6.21
            figdf['fig_ylabel']='Normalized weight'
            figdf['fig_xlabel']='Initial mean weight (nS)'
            # figdf['fig_dyticks']=.2
            # figdf['fig_dxticks']=20
            # # trace level parameters
            figdf['trace_ealpha']=1
            figdf['error_style']='shade'
            figdf['trace_linewidth']=4
            # figdf.at[slice(None), 'trace_ecolor'] = gray
            figdf['fig_xscale_log']=False
            figdf['fig_barwidth']=0.8
            figdf['fig_data_style']='point'
            figdf['fig_xtick_decimals']=1
            figdf['fig_ytick_decimals']=0
            figdf['fig_set_xscale']='linear'
            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]<0:
                    figdf.at[key, 'trace_color']=self.blue
                    figdf.at[key, 'trace_ecolor']=self.blue_light
                # apical
                if key[2][0]==0:
                    figdf.at[key, 'trace_color']=self.black
                    figdf.at[key, 'trace_ecolor']=self.gray
                # basal
                if key[2][0]>0:
                    figdf.at[key, 'trace_color']=self.red
                    figdf.at[key, 'trace_ecolor']=self.red_light
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

class figsetup_exp_reduced_neuron_polarization_piecewise(FigSetup):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        super(figsetup_exp_reduced_neuron_polarization_piecewise, self).__init__(**kwargs)

    def shapeplot_polarization(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]

            # conditions for each figure
            #----------------------------
            figdict = {
                # figure
                #-------------------------------------------------
                'anodal_5':{
                    # subgroup
                    'anodal':[
                        # trace
                        (5),
                    ]
                },
                'cathodal_5':{
                    # subgroup
                    'cathodal':[
                        # trace
                        (-5),
                    ]
                },
            }

            return figdict
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

    def dose_response_polarization(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            figdict = {
                # # figure
                # #-------------------------------------------------
                'all':{
                    # subgroup
                    'all':[
                        # trace
                        ('soma', 0, 0),
                        ('apic', 0, 10),
                        ('dend',0, 2),
                    ]
                },
            }

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            # figdf['fig_nyticks']=5
            # figdf['fig_nxticks']=10
            figdf['fig_dyticks']=1.
            figdf['fig_dxticks']=1.
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            figdf['fig_ymin']=-2
            figdf['fig_ymax']=2.1
            figdf['fig_xmin']=-20.
            figdf['fig_xmax']=21.
            figdf['fig_ylabel']='Polarization (mV)'
            figdf['fig_xlabel']='Electric field (V/m)'
            # figdf['fig_dyticks']=.2
            # figdf['fig_dxticks']=20
            # # trace level parameters
            figdf['trace_ealpha']=1
            figdf['error_style']='shade'
            figdf['trace_linewidth']=4
            # figdf.at[slice(None), 'trace_ecolor'] = gray
            figdf['fig_xscale_log']=False
            figdf['fig_barwidth']=0.8
            figdf['fig_data_style']='point'
            figdf['fig_xtick_decimals']=0
            figdf['fig_ytick_decimals']=0
            figdf['fig_set_xscale']='linear'
            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]=='soma':
                    figdf.at[key, 'trace_color']=self.blue
                # apical
                if key[2][0]=='apic':
                    figdf.at[key, 'trace_color']=self.black
                # basal
                if key[2][0]=='dend':
                    figdf.at[key, 'trace_color']=self.red
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

class figsetup_exp_full_neuron_polarization_dose_response(FigSetup):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        super(figsetup_exp_full_neuron_polarization_dose_response, self).__init__(**kwargs)

    def shapeplot_polarization(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]

            # conditions for each figure
            #----------------------------
            figdict = {
                # figure
                #-------------------------------------------------
                'anodal_20':{
                    # subgroup
                    'anodal':[
                        # trace
                        (20),
                    ]
                },
                'cathodal_20':{
                    # subgroup
                    'cathodal':[
                        # trace
                        (-20),
                    ]
                },
                'anodal_5':{
                    # subgroup
                    'anodal':[
                        # trace
                        (5),
                    ]
                },
                'cathodal_5':{
                    # subgroup
                    'cathodal':[
                        # trace
                        (-5),
                    ]
                },
            }

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

    def dose_response_polarization(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            figdict = {
                # # figure
                # #-------------------------------------------------
                'all':{
                    # subgroup
                    'all':[
                        # trace
                        ('soma', 0, 0),
                        # ('apic', 0, 10),
                        # ('dend',0, 2),
                    ]
                },
            }

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            # figdf['fig_nyticks']=5
            # figdf['fig_nxticks']=10
            figdf['fig_dyticks']=1.
            figdf['fig_dxticks']=5.
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            figdf['fig_ymin']=-2
            figdf['fig_ymax']=2.1
            figdf['fig_xmin']=-20.
            figdf['fig_xmax']=21.
            figdf['fig_ylabel']='Polarization (mV)'
            figdf['fig_xlabel']='Electric field (V/m)'
            # figdf['fig_dyticks']=.2
            # figdf['fig_dxticks']=20
            # # trace level parameters
            figdf['trace_ealpha']=1
            figdf['error_style']='shade'
            figdf['trace_linewidth']=4
            # figdf.at[slice(None), 'trace_ecolor'] = gray
            figdf['fig_xscale_log']=False
            figdf['fig_barwidth']=0.8
            figdf['fig_data_style']='point'
            figdf['fig_xtick_decimals']=0
            figdf['fig_ytick_decimals']=0
            figdf['fig_set_xscale']='linear'
            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]=='soma':
                    figdf.at[key, 'trace_color']=self.blue
                # apical
                if key[2][0]=='apic':
                    figdf.at[key, 'trace_color']=self.black
                # basal
                if key[2][0]=='dend':
                    figdf.at[key, 'trace_color']=self.red
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

class figsetup_exp_reduced_neuron_polarization_mirror_estimate(FigSetup):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        super(figsetup_exp_reduced_neuron_polarization_mirror_estimate, self).__init__(**kwargs)

    def shapeplot_polarization(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]

            # conditions for each figure
            #----------------------------
            figdict = {
                # figure
                #-------------------------------------------------
                'anodal_20':{
                    # subgroup
                    'anodal':[
                        # trace
                        (20),
                    ]
                },
                'cathodal_20':{
                    # subgroup
                    'cathodal':[
                        # trace
                        (-20),
                    ]
                },
                'anodal_5':{
                    # subgroup
                    'anodal':[
                        # trace
                        (5),
                    ]
                },
                'cathodal_5':{
                    # subgroup
                    'cathodal':[
                        # trace
                        (-5),
                    ]
                },
            }

            return figdict
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

    def dose_response_polarization(self, **kwargs):
        '''
        '''
        def _define_figdict(**kwargs):
            '''
            '''
            # print progress to terminal
            #-----------------------------
            print 'building figdf:', inspect.stack()[0][3]
            # conditions for each figure
            #----------------------------
            figdict = {
                # # figure
                # #-------------------------------------------------
                'all':{
                    # subgroup
                    'all':[
                        # trace
                        ('soma', 0, 0),
                        ('apic', 0, 10),
                        ('dend',0, 2),
                    ]
                },
            }

            return figdict
        
        def _update_figure_level_parameters(figdf, **kwargs):
            '''
            '''
            # fig parameters for all traces
            #---------------------
            # figure level parameters
            # figdf['fig_topercent']=False
            figdf['fig_ylim_all']=False
            figdf['fig_xlim_all']=False
            # figdf['trace_markersize']=10
            # # print figdf.fig_nyticks
            # figdf['fig_nyticks']=5
            # figdf['fig_nxticks']=10
            figdf['fig_dyticks']=1.
            figdf['fig_dxticks']=1.
            # figdf['fig_ylim_all']=True
            # figdf['fig_xlim_all']=True
            figdf['fig_ymin']=-2
            figdf['fig_ymax']=2.1
            figdf['fig_xmin']=-20.
            figdf['fig_xmax']=21.
            figdf['fig_ylabel']='Polarization (mV)'
            figdf['fig_xlabel']='Electric field (V/m)'
            # figdf['fig_dyticks']=.2
            # figdf['fig_dxticks']=20
            # # trace level parameters
            figdf['trace_ealpha']=1
            figdf['error_style']='shade'
            figdf['trace_linewidth']=4
            # figdf.at[slice(None), 'trace_ecolor'] = gray
            figdf['fig_xscale_log']=False
            figdf['fig_barwidth']=0.8
            figdf['fig_data_style']='point'
            figdf['fig_xtick_decimals']=0
            figdf['fig_ytick_decimals']=0
            figdf['fig_set_xscale']='linear'
            return figdf
        
        def _update_trace_level_parameters(figdf, **kwargs):
            '''
            '''
            # preallocate columns as object type
            #---------------------------------------
            init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
            figdf = functions._initialize_column(figdf, init_columns)
            # get all figure, subgroup, trace combinations (each combo specifies a single trace)
            #--------------------------------------------------------------------
            # reset index
            figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
            # get traces as ('figure','subgroup','trace')
            idx_keys = figdf.index.unique().values
            # iterate over traces
            for key in idx_keys:

                # set colors
                #------------------------------------
                # soma
                if key[2][0]=='soma':
                    figdf.at[key, 'trace_color']=self.blue
                # apical
                if key[2][0]=='apic':
                    figdf.at[key, 'trace_color']=self.black
                # basal
                if key[2][0]=='dend':
                    figdf.at[key, 'trace_color']=self.red
            return figdf
        
        def _generate_figdf(**kwargs):
            '''
            '''
            figdict = _define_figdict()
            figdf = self._build_figdf_from_dict(figdict=figdict, default_figdf=self.default_figdf)
            figdf = _update_figure_level_parameters(figdf=figdf)
            figdf = _update_trace_level_parameters(figdf=figdf)
            return figdf

        return _generate_figdf(**kwargs)

class BuildFigDF:
    '''
    '''
    def __init__(self, ):
        '''
        '''
        # set dpi for final png image
        #----------------------------
        self.dpi=350

        # colors for plots
        #-------------------
        self.black = (0,0,0)
        self.gray = (0.7,0.7, 0.7)
        self.red = (1,0,0)
        self.red_light = (1,0.7, 0.7)
        self.blue=(0,0,1)
        self.blue_light = (0.7, 0.7, 1)

    def polarization_shapeplot(self, ):
        '''
        '''
        # print progress to terminal
        #-----------------------------
        print 'building figdf:', inspect.stack()[0][3]

        # conditions for each figure
        #----------------------------
        figdict = {
            # figure
            #-------------------------------------------------
            'anodal_20':{
                # subgroup
                'anodal':[
                    # trace
                    (20),
                ]
            },
            'cathodal_20':{
                # subgroup
                'cathodal':[
                    # trace
                    (-20),
                ]
            },
            'anodal_5':{
                # subgroup
                'anodal':[
                    # trace
                    (5),
                ]
            },
            'cathodal_5':{
                # subgroup
                'cathodal':[
                    # trace
                    (-5),
                ]
            },
        }
        # load default figure parameters and colors
        #------------------------------------------
        default     = self._default()
        black       = self.black
        gray        = self.gray
        red         = self.red
        red_light   = self.red_light
        blue        = self.blue
        blue_light  = self.blue_light

        # create figdf
        #---------------------------------
        figdf = self._build_figdf_from_dict(figdict)
        # set trace level as index
        figdf = figdf.reset_index().set_index('trace')
        # get default parameters
        figdf_default = pd.concat([default]*len(figdf))
        # set index of defaultdf to match figdf
        figdf_default.index=figdf.index
        # add default df to figdf
        figdf = pd.concat([figdf, figdf_default], axis=1, ignore_index=False)

        # fig parameters for all traces
        #---------------------
        # # figure level parameters
        # # figdf['fig_topercent']=False
        # figdf['fig_ylim_all']=False
        # figdf['fig_xlim_all']=False
        # # figdf['trace_markersize']=10
        # # # print figdf.fig_nyticks
        # # figdf['fig_nyticks']=5
        # # figdf['fig_nxticks']=10
        # figdf['fig_dyticks']=4
        # figdf['fig_dxticks']=10
        # # figdf['fig_ylim_all']=True
        # # figdf['fig_xlim_all']=True
        # figdf['fig_ymin']=-74
        # figdf['fig_xmin']=0.
        # # figdf['fig_xmax']=30.
        # figdf['fig_ylabel']='Vm (mV)'
        # figdf['fig_xlabel']='Time (ms)'
        # # figdf['fig_dyticks']=.2
        # # figdf['fig_dxticks']=20
        # # # trace level parameters
        # figdf['trace_ealpha']=.7
        # figdf['error_style']='shade'
        # figdf['trace_linewidth']=4
        # # figdf['fig_xscale']=1./40
        # # figdf['fig_barwidth']=0.8
        # # figdf['fig_data_style']='point'
        # figdf['fig_xtick_decimals']=0
        # figdf['fig_ytick_decimals']=0
        # # figdf['fig_set_xscale']='symlog'



        # # individual trace parameters
        # #----------------------------
        # # preallocate columns as object type
        # figdf['trace_color']=None
        # figdf['trace_ecolor']=None
        # # figdf['fig_xticks']=None
        # # reset index
        # figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])

        # # get all figure, subgroup, trace combinations
        # idx_keys = figdf.index.unique().values
        # # iterate over combinations
        # for key in idx_keys:

        #     # set colors
        #     #------------------------------------
        #     # cathodal
        #     if key[2][0]<0:
        #         figdf.at[key, 'trace_color']=blue
        #         figdf.at[key, 'trace_ecolor']=blue
        #     # control
        #     if key[2][0]==0:
        #         figdf.at[key, 'trace_color']=black
        #         figdf.at[key, 'trace_ecolor']=black
        #     # anodal
        #     if key[2][0]>0:
        #         figdf.at[key, 'trace_color']=red
        #         figdf.at[key, 'trace_ecolor']=red

        return figdf

    def polarization_dose_response(self, ):
        '''
        '''
        # print progress to terminal
        #-----------------------------
        print 'building figdf:', inspect.stack()[0][3]

        # conditions for each figure
        #----------------------------
        figdict = {
            # # figure
            # #-------------------------------------------------
            'all':{
                # subgroup
                'all':[
                    # trace
                    ('soma', 0, 0),
                    ('apic', 0, 10),
                    ('dend',0, 2),
                ]
            },
        }
        # load default figure parameters and colors
        #------------------------------------------
        default   = self._default()
        black       = self.black
        gray        = self.gray
        red         = self.red
        red_light   = self.red_light
        blue        = self.blue
        blue_light  = self.blue_light

        # create figdf
        #---------------------------------
        figdf = self._build_figdf_from_dict(figdict)
        # set trace level as index
        figdf = figdf.reset_index().set_index('trace')
        # get default parameters
        figdf_default = pd.concat([default]*len(figdf))
        # set index of defaultdf to match figdf
        figdf_default.index=figdf.index
        # add default df to figdf
        figdf = pd.concat([figdf, figdf_default], axis=1, ignore_index=False)

        # fig parameters for all traces
        #---------------------
        # figure level parameters
        # figdf['fig_topercent']=False
        figdf['fig_ylim_all']=False
        figdf['fig_xlim_all']=False
        # figdf['trace_markersize']=10
        # # print figdf.fig_nyticks
        # figdf['fig_nyticks']=5
        # figdf['fig_nxticks']=10
        figdf['fig_dyticks']=5.
        figdf['fig_dxticks']=5.
        # figdf['fig_ylim_all']=True
        # figdf['fig_xlim_all']=True
        figdf['fig_ymin']=-15
        figdf['fig_ymax']=15.1
        figdf['fig_xmin']=-20.
        figdf['fig_xmax']=21.
        figdf['fig_ylabel']='Polarization (mV)'
        figdf['fig_xlabel']='Electric field (V/m)'
        # figdf['fig_dyticks']=.2
        # figdf['fig_dxticks']=20
        # # trace level parameters
        figdf['trace_ealpha']=1
        figdf['error_style']='shade'
        figdf['trace_linewidth']=4
        # figdf.at[slice(None), 'trace_ecolor'] = gray
        figdf['fig_xscale_log']=False
        figdf['fig_barwidth']=0.8
        figdf['fig_data_style']='point'
        figdf['fig_xtick_decimals']=0
        figdf['fig_ytick_decimals']=0
        figdf['fig_set_xscale']='linear'


        # individual trace parameters
        #----------------------------
        # preallocate columns as object type
        figdf['trace_color']=None
        figdf['trace_ecolor']=None
        figdf['fig_xticks']=None
        # reset index
        figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])

        # get all figure, subgroup, trace combinations
        idx_keys = figdf.index.unique().values
        # iterate over combinations
        for key in idx_keys:

            # set colors
            #------------------------------------
            # soma
            if key[2][0]=='soma':
                figdf.at[key, 'trace_color']=blue
            # apical
            if key[2][0]=='apic':
                figdf.at[key, 'trace_color']=black
            # basal
            if key[2][0]=='dend':
                figdf.at[key, 'trace_color']=red

        return figdf

    def _polarization_dose_response(self, ):
        '''
        '''
        # print progress to terminal
        #-----------------------------
        print 'building figdf:', inspect.stack()[0][3]

        # conditions for each figure
        #----------------------------
        figdict = {
            'soma_basal':{
                # subgroup
                'soma_basal':[
                    # trace
                    (('soma'),),
                    (('basal'),),
                ]
            },
        }
        # load default figure parameters and colors
        #------------------------------------------
        # default   = self._default()
        default   = analysis._default_figdf()
        black       = self.black
        gray        = self.gray
        red         = self.red
        red_light   = self.red_light
        blue        = self.blue
        blue_light  = self.blue_light

        # create figdf
        #---------------------------------
        figdf = self._build_figdf_from_dict(figdict)
        # set trace level as index
        figdf = figdf.reset_index().set_index('trace')
        # get default parameters
        figdf_default = pd.concat([default]*len(figdf))
        # set index of defaultdf to match figdf
        figdf_default.index=figdf.index
        # add default df to figdf
        figdf = pd.concat([figdf, figdf_default], axis=1, ignore_index=False)

        # fig parameters for all traces
        #---------------------
        # figure level parameters
        # figdf['fig_topercent']=False
        figdf['fig_ylim_all']=False
        figdf['fig_xlim_all']=False
        # figdf['trace_markersize']=10
        # # print figdf.fig_nyticks
        # figdf['fig_nyticks']=5
        # figdf['fig_nxticks']=10
        # figdf['fig_dyticks']=.2
        # figdf['fig_dxticks']=np.exp(10)
        # figdf['fig_dxticks']=5
        # figdf['fig_ylim_all']=True
        # figdf['fig_xlim_all']=True
        # figdf['fig_ymin']=1.2
        # figdf['fig_xmin']=-25.
        # figdf['fig_xmax']=26.
        figdf['fig_ylabel']='Polarization (mV)'
        figdf['fig_xlabel']='Electric field (V/m)'
        # figdf['fig_dyticks']=.2
        # figdf['fig_dxticks']=20
        # # trace level parameters
        figdf['trace_ealpha']=1
        figdf['error_style']='shade'
        figdf['trace_linewidth']=4
        # figdf.at[slice(None), 'trace_ecolor'] = gray
        figdf['fig_xscale_log']=False
        figdf['fig_barwidth']=0.8
        figdf['fig_data_style']='point'
        figdf['fig_xtick_decimals']=0
        figdf['fig_ytick_decimals']=1
        figdf['fig_set_xscale']='linear'
        figdf['fig_hlines']=0
        figdf['fig_vlines']=0


        # individual trace parameters
        #----------------------------
        # preallocate columns as object type
        figdf['trace_color']=None
        figdf['trace_ecolor']=None
        figdf['fig_xticks']=None
        # reset index
        figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])


        # get all figure, subgroup, trace combinations
        idx_keys = figdf.index.unique().values
        # iterate over combinations
        for key in idx_keys:
            if 'soma' in key[2]:
                figdf.at[key, 'trace_color']=black
                figdf.at[key, 'trace_ecolor']=black
            if 'basal' in key[2]:
                figdf.at[key, 'trace_color']=gray
                figdf.at[key, 'trace_ecolor']=gray

        return figdf

    def _trace_mean(self, ):
        '''
        '''
        # print progress to terminal
        #-----------------------------
        print 'building figdf:', inspect.stack()[0][3]

        # conditions for each figure
        #----------------------------
        figdict = {
            # figure
            #-------------------------------------------------
            'all_soma':{
                # subgroup
                'all_soma':[
                    # trace
                    (-20, (4,5,6,7,8), 'soma'),
                    # (-5, ),
                    # (-1, ),
                    # (-0.5, ),
                    (0, (4,5,6,7,8), 'soma'),
                    # (0.5, ),
                    # (1, ),
                    # (5, ),
                    (20, (4,5,6,7,8), 'soma'),
                ]
            },
            'all_axon':{
                # subgroup
                'all_axon':[
                    # trace
                    (-20, (4,5,6,7,8), 'axon'),
                    # (-5, ),
                    # (-1, ),
                    # (-0.5, ),
                    (0, (4,5,6,7,8), 'axon'),
                    # (0.5, ),
                    # (1, ),
                    # (5, ),
                    (20, (4,5,6,7,8), 'axon'),
                ]
            },
            # figure
            #-------------------------------------------------
            'all_dendrite':{
                # subgroup
                'all_dendrite':[
                    # trace
                    (-20, (6,8,10,12), 'basal'),
                    # (-5, ),
                    # (-1, ),
                    # (-0.5, ),
                    (0, (6,8,10,12), 'basal'),
                    # (0.5, ),
                    # (1, ),
                    # (5, ),
                    (20, (6,8,10,12), 'basal'),
                ]
            },
            # figure
            #-------------------------------------------------
            '14_dendrite':{
                # subgroup
                '14_dendrite':[
                    # trace
                    (-20, (14), 'basal'),
                    # (-5, ),
                    # (-1, ),
                    # (-0.5, ),
                    (0, (14), 'basal'),
                    # (0.5, ),
                    # (1, ),
                    # (5, ),
                    (20, (14), 'basal'),
                ]
            },
            # figure
            #-------------------------------------------------
            '8_dendrite':{
                # subgroup
                '8_dendrite':[
                    # trace
                    (-20, (8), 'basal'),
                    # (-5, ),
                    # (-1, ),
                    # (-0.5, ),
                    (0, (8), 'basal'),
                    # (0.5, ),
                    # (1, ),
                    # (5, ),
                    (20, (8), 'basal'),
                ]
            },
            # figure
            #-------------------------------------------------
            '8_soma':{
                # subgroup
                '8_soma':[
                    # trace
                    (-20, (8), 'soma'),
                    # (-5, ),
                    # (-1, ),
                    # (-0.5, ),
                    (0, (8), 'soma'),
                    # (0.5, ),
                    # (1, ),
                    # (5, ),
                    (20, (8), 'soma'),
                ]
            },
            # figure
            #-------------------------------------------------
            '6_dendrite':{
                # subgroup
                '6_dendrite':[
                    # trace
                    (-20, (6), 'basal'),
                    # (-5, ),
                    # (-1, ),
                    # (-0.5, ),
                    (0, (6), 'basal'),
                    # (0.5, ),
                    # (1, ),
                    # (5, ),
                    (20, (6), 'basal'),
                ]
            },

        }
        # load default figure parameters and colors
        #------------------------------------------
        default     = self._default()
        black       = self.black
        gray        = self.gray
        red         = self.red
        red_light   = self.red_light
        blue        = self.blue
        blue_light  = self.blue_light

        # create figdf
        #---------------------------------
        figdf = self._build_figdf_from_dict(figdict)
        # set trace level as index
        figdf = figdf.reset_index().set_index('trace')
        # get default parameters
        figdf_default = pd.concat([default]*len(figdf))
        # set index of defaultdf to match figdf
        figdf_default.index=figdf.index
        # add default df to figdf
        figdf = pd.concat([figdf, figdf_default], axis=1, ignore_index=False)

        # fig parameters for all traces
        #---------------------
        # figure level parameters
        # figdf['fig_topercent']=False
        figdf['fig_ylim_all']=False
        figdf['fig_xlim_all']=False
        # figdf['trace_markersize']=10
        # # print figdf.fig_nyticks
        # figdf['fig_nyticks']=5
        # figdf['fig_nxticks']=10
        figdf['fig_dyticks']=4
        figdf['fig_dxticks']=10
        # figdf['fig_ylim_all']=True
        # figdf['fig_xlim_all']=True
        figdf['fig_ymin']=-74
        figdf['fig_xmin']=0.
        # figdf['fig_xmax']=30.
        figdf['fig_ylabel']='Vm (mV)'
        figdf['fig_xlabel']='Time (ms)'
        # figdf['fig_dyticks']=.2
        # figdf['fig_dxticks']=20
        # # trace level parameters
        figdf['trace_ealpha']=.7
        figdf['error_style']='shade'
        figdf['trace_linewidth']=4
        # figdf['fig_xscale']=1./40
        # figdf['fig_barwidth']=0.8
        # figdf['fig_data_style']='point'
        figdf['fig_xtick_decimals']=0
        figdf['fig_ytick_decimals']=0
        # figdf['fig_set_xscale']='symlog'



        # individual trace parameters
        #----------------------------
        # preallocate columns as object type
        figdf['trace_color']=None
        figdf['trace_ecolor']=None
        # figdf['fig_xticks']=None
        # reset index
        figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])

        # get all figure, subgroup, trace combinations
        idx_keys = figdf.index.unique().values
        # iterate over combinations
        for key in idx_keys:

            # set colors
            #------------------------------------
            # cathodal
            if key[2][0]<0:
                figdf.at[key, 'trace_color']=blue
                figdf.at[key, 'trace_ecolor']=blue
            # control
            if key[2][0]==0:
                figdf.at[key, 'trace_color']=black
                figdf.at[key, 'trace_ecolor']=black
            # anodal
            if key[2][0]>0:
                figdf.at[key, 'trace_color']=red
                figdf.at[key, 'trace_ecolor']=red

        return figdf

    def _trace_mean_filtered(self, ):
        '''
        '''
        # print progress to terminal
        #-----------------------------
        print 'building figdf:', inspect.stack()[0][3]

        # conditions for each figure
        #----------------------------
        figdict = {
            # figure
            #-------------------------------------------------
            'all_soma':{
                # subgroup
                'all_soma':[
                    # trace
                    (-20, (4,5,6,7,8), 'soma'),
                    # (-5, ),
                    # (-1, ),
                    # (-0.5, ),
                    (0, (4,5,6,7,8), 'soma'),
                    # (0.5, ),
                    # (1, ),
                    # (5, ),
                    (20, (4,5,6,7,8), 'soma'),
                ]
            },
            'all_axon':{
                # subgroup
                'all_axon':[
                    # trace
                    (-20, (4,5,6,7,8), 'axon'),
                    # (-5, ),
                    # (-1, ),
                    # (-0.5, ),
                    (0, (4,5,6,7,8), 'axon'),
                    # (0.5, ),
                    # (1, ),
                    # (5, ),
                    (20, (4,5,6,7,8), 'axon'),
                ]
            },
            # figure
            #-------------------------------------------------
            'all_dendrite':{
                # subgroup
                'all_dendrite':[
                    # trace
                    (-20, (6,8,10,12), 'basal'),
                    # (-5, ),
                    # (-1, ),
                    # (-0.5, ),
                    (0, (6,8,10,12), 'basal'),
                    # (0.5, ),
                    # (1, ),
                    # (5, ),
                    (20, (6,8,10,12), 'basal'),
                ]
            },
            # figure
            #-------------------------------------------------
            '14_dendrite':{
                # subgroup
                '14_dendrite':[
                    # trace
                    (-20, (14), 'basal'),
                    # (-5, ),
                    # (-1, ),
                    # (-0.5, ),
                    (0, (14), 'basal'),
                    # (0.5, ),
                    # (1, ),
                    # (5, ),
                    (20, (14), 'basal'),
                ]
            },
            # figure
            #-------------------------------------------------
            '8_dendrite':{
                # subgroup
                '8_dendrite':[
                    # trace
                    (-20, (8), 'basal'),
                    # (-5, ),
                    # (-1, ),
                    # (-0.5, ),
                    (0, (8), 'basal'),
                    # (0.5, ),
                    # (1, ),
                    # (5, ),
                    (20, (8), 'basal'),
                ]
            },
            # figure
            #-------------------------------------------------
            '8_soma':{
                # subgroup
                '8_soma':[
                    # trace
                    (-20, (8), 'soma'),
                    # (-5, ),
                    # (-1, ),
                    # (-0.5, ),
                    (0, (8), 'soma'),
                    # (0.5, ),
                    # (1, ),
                    # (5, ),
                    (20, (8), 'soma'),
                ]
            },
            # figure
            #-------------------------------------------------
            '6_dendrite':{
                # subgroup
                '6_dendrite':[
                    # trace
                    (-20, (6), 'basal'),
                    # (-5, ),
                    # (-1, ),
                    # (-0.5, ),
                    (0, (6), 'basal'),
                    # (0.5, ),
                    # (1, ),
                    # (5, ),
                    (20, (6), 'basal'),
                ]
            },

        }
        # load default figure parameters and colors
        #------------------------------------------
        default     = self._default()
        black       = self.black
        gray        = self.gray
        red         = self.red
        red_light   = self.red_light
        blue        = self.blue
        blue_light  = self.blue_light

        # create figdf
        #---------------------------------
        figdf = self._build_figdf_from_dict(figdict)
        # set trace level as index
        figdf = figdf.reset_index().set_index('trace')
        # get default parameters
        figdf_default = pd.concat([default]*len(figdf))
        # set index of defaultdf to match figdf
        figdf_default.index=figdf.index
        # add default df to figdf
        figdf = pd.concat([figdf, figdf_default], axis=1, ignore_index=False)

        # fig parameters for all traces
        #---------------------
        # figure level parameters
        # figdf['fig_topercent']=False
        figdf['fig_ylim_all']=False
        figdf['fig_xlim_all']=False
        # figdf['trace_markersize']=10
        # # print figdf.fig_nyticks
        # figdf['fig_nyticks']=5
        # figdf['fig_nxticks']=10
        figdf['fig_dyticks']=1
        figdf['fig_dxticks']=10
        # figdf['fig_ylim_all']=True
        # figdf['fig_xlim_all']=True
        figdf['fig_ymin']=-3
        figdf['fig_xmin']=0.
        # figdf['fig_xmax']=30.
        figdf['fig_ylabel']='Vm (mV)'
        figdf['fig_xlabel']='Time (ms)'
        # figdf['fig_dyticks']=.2
        # figdf['fig_dxticks']=20
        # # trace level parameters
        figdf['trace_ealpha']=.7
        figdf['error_style']='shade'
        figdf['trace_linewidth']=4
        # figdf['fig_xscale']=1./40
        # figdf['fig_barwidth']=0.8
        # figdf['fig_data_style']='point'
        figdf['fig_xtick_decimals']=0
        figdf['fig_ytick_decimals']=0
        # figdf['fig_set_xscale']='symlog'



        # individual trace parameters
        #----------------------------
        # preallocate columns as object type
        figdf['trace_color']=None
        figdf['trace_ecolor']=None
        # figdf['fig_xticks']=None
        # reset index
        figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])

        # get all figure, subgroup, trace combinations
        idx_keys = figdf.index.unique().values
        # iterate over combinations
        for key in idx_keys:

            # set colors
            #------------------------------------
            # cathodal
            if key[2][0]<0:
                figdf.at[key, 'trace_color']=blue
                figdf.at[key, 'trace_ecolor']=blue
            # control
            if key[2][0]==0:
                figdf.at[key, 'trace_color']=black
                figdf.at[key, 'trace_ecolor']=black
            # anodal
            if key[2][0]>0:
                figdf.at[key, 'trace_color']=red
                figdf.at[key, 'trace_ecolor']=red

        return figdf

    def _weight_hist(self, ):
        '''
        '''
        # print progress to terminal
        #-----------------------------
        print 'building figdf:', inspect.stack()[0][3]

        # conditions for each figure
        #----------------------------
        figdict = {
            # figure
            #-------------------------------------------------
            'all_dendrite_20':{
                # subgroup
                'all_dendrite_20':[
                    # trace
                    (-20, (5,6,7,8,9), 'basal'),
                    (0, (5,6,7,8,9), 'basal'),
                    (20, (5,6,7,8,9), 'basal'),
                ]
            },
            'all_dendrite_1':{
                # subgroup
                'all_dendrite_1':[
                    # trace
                    (-1, (5,6,7,8,9), 'basal'),
                    (0, (5,6,7,8,9), 'basal'),
                    (1, (5,6,7,8,9), 'basal'),
                ]
            },
            '5-7_20':{
                # subgroup
                '5-7_20':[
                    # trace
                    (-20, (5,6,7), 'basal'),
                    (0, (5,6,7,), 'basal'),
                    (20, (5,6,7,), 'basal'),
                ]
            },
            '5-7_1':{
                # subgroup
                '5-7_1':[
                    # trace
                    (-1, (5,6,7,), 'basal'),
                    (0, (5,6,7,), 'basal'),
                    (1, (5,6,7,), 'basal'),
                ]
            },

            # '5syn_1vm':{
            #     # subgroup
            #     '5syn_1vm':[
            #         # trace
            #         (-1, (5,), 'basal'),
            #         (0, (5,), 'basal'),
            #         (1, (5,), 'basal'),
            #     ]
            # },
            # '6syn_1vm':{
            #     # subgroup
            #     '6syn_1vm':[
            #         # trace
            #         (-1, (6,), 'basal'),
            #         (0, (6,), 'basal'),
            #         (1, (6,), 'basal'),
            #     ]
            # },
            # '7syn_1vm':{
            #     # subgroup
            #     '7syn_1vm':[
            #         # trace
            #         (-1, (7,), 'basal'),
            #         (0, (7,), 'basal'),
            #         (1, (7,), 'basal'),
            #     ]
            # },
            # '8syn_1vm':{
            #     # subgroup
            #     '8syn_1vm':[
            #         # trace
            #         (-1, (8,), 'basal'),
            #         (0, (8,), 'basal'),
            #         (1, (8,), 'basal'),
            #     ]
            # },
            # '9syn_1vm':{
            #     # subgroup
            #     '9syn_1vm':[
            #         # trace
            #         (-1, (9,), 'basal'),
            #         (0, (9,), 'basal'),
            #         (1, (9,), 'basal'),
            #     ]
            # },
            # '5syn_20vm':{
            #     # subgroup
            #     '5syn_20vm':[
            #         # trace
            #         (-20, (5,), 'basal'),
            #         (0, (5,), 'basal'),
            #         (20, (5,), 'basal'),
            #     ]
            # },
            # '6syn_20vm':{
            #     # subgroup
            #     '6syn_20vm':[
            #         # trace
            #         (-20, (6,), 'basal'),
            #         (0, (6,), 'basal'),
            #         (20, (6,), 'basal'),
            #     ]
            # },
            # '7syn_20vm':{
            #     # subgroup
            #     '7syn_20vm':[
            #         # trace
            #         (-20, (7,), 'basal'),
            #         (0, (7,), 'basal'),
            #         (20, (7,), 'basal'),
            #     ]
            # },
            # '8syn_20vm':{
            #     # subgroup
            #     '8syn_20vm':[
            #         # trace
            #         (-20, (8,), 'basal'),
            #         (0, (8,), 'basal'),
            #         (20, (8,), 'basal'),
            #     ]
            # },
            # '9syn_20vm':{
            #     # subgroup
            #     '9syn_20vm':[
            #         # trace
            #         (-20, (9,), 'basal'),
            #         (0, (9,), 'basal'),
            #         (20, (9,), 'basal'),
            #     ]
            # },
        }
        # load default figure parameters and colors
        #------------------------------------------
        default     = self._default()
        black       = self.black
        gray        = self.gray
        red         = self.red
        red_light   = self.red_light
        blue        = self.blue
        blue_light  = self.blue_light

        # create figdf
        #---------------------------------
        figdf = self._build_figdf_from_dict(figdict)
        # set trace level as index
        figdf = figdf.reset_index().set_index('trace')
        # get default parameters
        figdf_default = pd.concat([default]*len(figdf))
        # set index of defaultdf to match figdf
        figdf_default.index=figdf.index
        # add default df to figdf
        figdf = pd.concat([figdf, figdf_default], axis=1, ignore_index=False)

        # fig parameters for all traces
        #---------------------
        # figure level parameters
        # figdf['fig_topercent']=False
        figdf['fig_ylim_all']=False
        figdf['fig_xlim_all']=False
        # figdf['trace_markersize']=10
        # # print figdf.fig_nyticks
        # figdf['fig_nyticks']=5
        # figdf['fig_nxticks']=10
        # figdf['fig_dyticks']=4
        # figdf['fig_dxticks']=10
        # figdf['fig_ylim_all']=True
        # figdf['fig_xlim_all']=True
        # figdf['fig_ymin']=-74
        # figdf['fig_xmin']=0.
        # figdf['fig_xmax']=30.
        figdf['fig_ylabel']='Probability'
        figdf['fig_xlabel']='dW DCS'
        # figdf['fig_dyticks']=.2
        # figdf['fig_dxticks']=20
        # # trace level parameters
        figdf['trace_ealpha']=.7
        figdf['error_style']='shade'
        figdf['trace_linewidth']=3
        # figdf['fig_xscale']=1./40
        # figdf['fig_barwidth']=0.8
        # figdf['fig_data_style']='point'
        # figdf['fig_xtick_decimals']=0
        # figdf['fig_ytick_decimals']=0
        # figdf['fig_set_xscale']='symlog'



        # individual trace parameters
        #----------------------------
        # preallocate columns as object type
        figdf['trace_color']=None
        figdf['trace_ecolor']=None
        # figdf['fig_xticks']=None
        # reset index
        figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])

        # get all figure, subgroup, trace combinations
        idx_keys = figdf.index.unique().values
        # iterate over combinations
        for key in idx_keys:

            # set colors
            #------------------------------------
            # cathodal
            if key[2][0]<0:
                figdf.at[key, 'trace_color']=blue
                figdf.at[key, 'trace_ecolor']=blue
            # control
            if key[2][0]==0:
                figdf.at[key, 'trace_color']=black
                figdf.at[key, 'trace_ecolor']=black
            # anodal
            if key[2][0]>0:
                figdf.at[key, 'trace_color']=red
                figdf.at[key, 'trace_ecolor']=red

        return figdf

    def _default(self):
        '''
        '''
        # set dpi for final png image
        #----------------------------
        self.dpi=350

        # colors for plots
        #-------------------
        self.black = (0,0,0)
        self.gray = (0.7,0.7, 0.7)
        self.red = (1,0,0)
        self.red_light = (1,0.7, 0.7)
        self.blue=(0,0,1)
        self.blue_light = (0.7, 0.7, 1)

        all_dict={
        # hide the top and right axes boundaries
            'fig_dpi':350,
            'fig_boxoff':True,
            # axes and tick labels
            'fig_axes_linewidth':[4],
            'fig_xlabel_fontsize':[25],
            'fig_ylabel_fontsize':[25],
            'fig_xlabel_fontweight':'heavy',
            'fig_ylabel_fontweight':'heavy',
            'fig_xtick_fontsize':[15],
            'fig_ytick_fontsize':[15],
            'fig_xtick_fontweight':'heavy',
            'fig_ytick_fontweight':'heavy',
            # figure tight layout
            'fig_tight_layout':True,
        }
        all_df = pd.DataFrame(all_dict, dtype='object')

        return all_df

    def _build_figdf_from_dict(self, figdict):
        '''
        '''
        # build multiindex for trace parameters
        multi_list = []
        level_names = ['figure','subgroup', 'trace']
        for level_1_key, level_1 in figdict.iteritems():
            for level_2_key, level_2 in level_1.iteritems():
                for level_3_i, level_3_key in enumerate(level_2):
                    multi_list.append((level_1_key, level_2_key, level_3_key))

        multiindex = pd.MultiIndex.from_tuples(multi_list, names=level_names)
        # build dataframe
        figdf = pd.DataFrame(index=multiindex, dtype='object')

        return figdf
        
class PolarizationShapePlot(BuildFigDF):
    '''
    '''
    def __init__(self, **kwargs):

        super(PolarizationShapePlot, self).__init__(**kwargs)

    def _generate_figdf(self, **kwargs):
        '''
        '''
        self.figdict = self._define_figdict()
        self.figdf = self._build_figdf_from_dict(figdict=self.figdict, default_figdf=self.default_figdf)
        
    def _define_figdict(self, **kwargs):
        '''
        '''
        # print progress to terminal
        #-----------------------------
        print 'building figdf:', inspect.stack()[0][3]

        # conditions for each figure
        #----------------------------
        figdict = {
            # figure
            #-------------------------------------------------
            'anodal_20':{
                # subgroup
                'anodal':[
                    # trace
                    (20),
                ]
            },
            'cathodal_20':{
                # subgroup
                'cathodal':[
                    # trace
                    (-20),
                ]
            },
            'anodal_5':{
                # subgroup
                'anodal':[
                    # trace
                    (5),
                ]
            },
            'cathodal_5':{
                # subgroup
                'cathodal':[
                    # trace
                    (-5),
                ]
            },
        }

        return figdict

class PolarizationDoseResponse(BuildFigDF):
    '''
    '''
    def __init__(self, **kwargs):

        super(PolarizationShapePlot, self).__init__(**kwargs)

    def _generate_figdf(self, **kwargs):
        '''
        '''
        self.figdict = self._define_figdict()
        self.figdf = self._build_figdf_from_dict(figdict=self.figdict, default_figdf=self.default_figdf)
        
    def _define_figdict(self, **kwargs):
        '''
        '''
        # print progress to terminal
        #-----------------------------
        print 'building figdf:', inspect.stack()[0][3]
        # conditions for each figure
        #----------------------------
        figdict = {
            # # figure
            # #-------------------------------------------------
            'all':{
                # subgroup
                'all':[
                    # trace
                    ('soma', 0, 0),
                    ('apic', 0, 10),
                    ('dend',0, 2),
                ]
            },
        }

        return figdict

    def _update_figure_level_parameters(self, figdf, **kwargs):
        '''
        '''
        # fig parameters for all traces
        #---------------------
        # figure level parameters
        # figdf['fig_topercent']=False
        figdf['fig_ylim_all']=False
        figdf['fig_xlim_all']=False
        # figdf['trace_markersize']=10
        # # print figdf.fig_nyticks
        # figdf['fig_nyticks']=5
        # figdf['fig_nxticks']=10
        figdf['fig_dyticks']=5.
        figdf['fig_dxticks']=5.
        # figdf['fig_ylim_all']=True
        # figdf['fig_xlim_all']=True
        figdf['fig_ymin']=-15
        figdf['fig_ymax']=15.1
        figdf['fig_xmin']=-20.
        figdf['fig_xmax']=21.
        figdf['fig_ylabel']='Polarization (mV)'
        figdf['fig_xlabel']='Electric field (V/m)'
        # figdf['fig_dyticks']=.2
        # figdf['fig_dxticks']=20
        # # trace level parameters
        figdf['trace_ealpha']=1
        figdf['error_style']='shade'
        figdf['trace_linewidth']=4
        # figdf.at[slice(None), 'trace_ecolor'] = gray
        figdf['fig_xscale_log']=False
        figdf['fig_barwidth']=0.8
        figdf['fig_data_style']='point'
        figdf['fig_xtick_decimals']=0
        figdf['fig_ytick_decimals']=0
        figdf['fig_set_xscale']='linear'

        return figdf

    def _update_trace_level_parameters(self, figdf, **kwargs):
        '''
        '''
        # preallocate columns as object type
        #---------------------------------------
        init_columns = ['trace_color', 'trace_ecolor', 'fig_xticks']
        figdf = functions._initialize_column(figdf, init_columns)
        # get all figure, subgroup, trace combinations (each combo specifies a single trace)
        #--------------------------------------------------------------------
        # reset index
        figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
        # get traces as ('figure','subgroup','trace')
        idx_keys = figdf.index.unique().values
        # iterate over traces
        for key in idx_keys:

            # set colors
            #------------------------------------
            # soma
            if key[2][0]=='soma':
                figdf.at[key, 'trace_color']=self.blue
            # apical
            if key[2][0]=='apic':
                figdf.at[key, 'trace_color']=self.black
            # basal
            if key[2][0]=='dend':
                figdf.at[key, 'trace_color']=self.red

        return figdf

