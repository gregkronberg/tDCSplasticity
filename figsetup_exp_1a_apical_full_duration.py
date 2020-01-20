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

    def _ltp_bar(self, ):
        '''
        '''
        # print progress to terminal
        #-----------------------------
        print 'building figdf:', inspect.stack()[0][3]

        # conditions for each figure
        #----------------------------
        figdict = {

            'TBS_full':{
                # subgroup
                '20Vm':[
                    # trace
                    ((6,8,10,12), -20, ('apical_tuft', 'apical_trunk') ),
                    ((6,8,10,12), 20, ('apical_tuft', 'apical_trunk')),
                    ((6,8,10,12), 0, ('apical_tuft', 'apical_trunk')),
                ],
            },
            'TBS_full_12syn':{
                # subgroup
                '20Vm':[
                    # trace
                    ((12), -20, ('apical_tuft', 'apical_trunk') ),
                    ((12), 20, ('apical_tuft', 'apical_trunk')),
                    ((12), 0, ('apical_tuft', 'apical_trunk')),
                ],
            },
            'TBS_full_8-12syn':{
                # subgroup
                '20Vm':[
                    # trace
                    ((8,10,12), -20, ('apical_tuft', 'apical_trunk') ),
                    ((8,10,12), 20, ('apical_tuft', 'apical_trunk')),
                    ((8,10,12), 0, ('apical_tuft', 'apical_trunk')),
                ],
            },
            'TBS_full_6syn':{
                # subgroup
                '20Vm':[
                    # trace
                    ((6), -20, ('apical_tuft', 'apical_trunk') ),
                    ((6), 20, ('apical_tuft', 'apical_trunk')),
                    ((6), 0, ('apical_tuft', 'apical_trunk')),
                ],
            },
            # 'acs all':{
            #     # subgroup
            #     '20Vm':[
            #         # trace
            #         ((4,5,6,8,9,10,), 20, 0 ),
            #         ((4,5,6,8,9,10,), 20, np.pi),
            #         ((4,5,6,8,9,10,), 0),
            #     ],
            # },
            # 'acs_6_8_10_12':{
            #     # subgroup
            #     '20Vm':[
            #         # trace
            #         ((6,8,10,12), 20, 0 ),
            #         ((6,8,10,12), 20, np.pi),
            #         ((6,8,10,12), 0),
            #     ],
            # },
        }
        # load default figure parameters and colors
        #------------------------------------------
        # default     = self._default()
        default     = analysis._default_figdf()
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
        figdf['fig_dyticks']=.1
        # figdf['fig_dxticks']=10
        # figdf['fig_ylim_all']=True
        # figdf['fig_xlim_all']=True
        figdf['fig_ymin']=1.
        # figdf['fig_ymax']=1.51
        figdf['fig_xmin']=0.
        # figdf['fig_xmax']=30.
        figdf['fig_ylabel']='Norm. weight'
        # figdf['fig_xlabel']='Time (ms)'
        figdf['fig_dyticks']=.1
        figdf['fig_dxticks']=1
        # # trace level parameters
        # figdf['trace_ealpha']=.7
        # figdf['error_style']='shade'
        # figdf['trace_linewidth']=4
        # figdf['fig_xscale']=1./40
        figdf['fig_barwidth']=1
        # figdf['fig_data_style']='point'
        # figdf['fig_xtick_decimals']=0
        figdf['fig_ytick_decimals']=1
        # figdf['fig_set_xscale']='symlog'

        # set subgroup parameters
        #-----------------------------------
        figdf = figdf.reset_index().set_index('subgroup')
        figdf.at['20Vm','sub_location']=0

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
            if key[2][1]==0:
                figdf.at[key, 'trace_location']=1
                figdf.at[key, 'trace_label']='control'
                figdf.at[key, 'trace_color']=black
                figdf.at[key, 'trace_ecolor']=black
            elif key[2][1]>0:
                figdf.at[key, 'trace_location']=2
                figdf.at[key, 'trace_label']='anodal'
                figdf.at[key, 'trace_color']=red
                figdf.at[key, 'trace_ecolor']=red
            elif key[2][1]<0:
                figdf.at[key, 'trace_location']=0
                figdf.at[key, 'trace_label']='cathodal'
                figdf.at[key, 'trace_color']=blue
                figdf.at[key, 'trace_ecolor']=blue

        return figdf

    def _xcorr(self, ):
        '''
        '''
        # print progress to terminal
        #-----------------------------
        print 'building figdf:', inspect.stack()[0][3]

        # conditions for each figure
        #----------------------------
        figdict = {
            'weak_unpaired':{
                # subgroup
                'weak unpaired':[
                    # trace
                    (0, '2', (4,5,6,8,), 1, ('apical_tuft', 'apical_trunk')),
                    (20,'2', (4,5,6,8,), 1, ('apical_tuft', 'apical_trunk')),
                    # (-20,'2', (4,5,6,8,), 1, ('apical_tuft', 'apical_trunk')),
                ],
                # subgroup
                'weak paired':[
                    # trace
                    (0, '2', (4,5,6,8,), 2, ('apical_tuft', 'apical_trunk')),
                    (20,'2', (4,5,6,8,), 2, ('apical_tuft', 'apical_trunk')),
                    # (-20,'2', (4,5,6,8,), 1, ('apical_tuft', 'apical_trunk')),
                ],

                # subgroup
                # 'strong paired':[
                #     # trace
                #     (0, '1', (4,5,6,8,), 2, ('apical_tuft', 'apical_trunk')),
                #     (20,'1', (4,5,6,8,), 2, ('apical_tuft', 'apical_trunk')),
                #     # (-20,'2', (4,5,6,8,), 1, ('apical_tuft', 'apical_trunk')),
                # ],
            },
            'weak_unpaired 4-20':{
                # subgroup
                'weak unpaired 4-20':[
                    # trace
                    (0, '2', (4,5,6,8,10,12,16,20), 1, ('apical_tuft', 'apical_trunk')),
                    (20,'2', (4,5,6,8,10,12,16,20), 1, ('apical_tuft', 'apical_trunk')),
                    # (-20,'2', (4,5,6,8,), 1, ('apical_tuft', 'apical_trunk')),
                ],
                # subgroup
                'weak paired 4-20':[
                    # trace
                    (0, '2', (4,5,6,8,10,12,16,20), 2, ('apical_tuft', 'apical_trunk')),
                    (20,'2', (4,5,6,8,10,12,16,20), 2, ('apical_tuft', 'apical_trunk')),
                    # (-20,'2', (4,5,6,8,), 1, ('apical_tuft', 'apical_trunk')),
                ],
            }
            # 'strong_paired':{
            #     # subgroup
            #     'strong unpaired':[
            #         # trace
            #         (0, '1', (4,5,6,8,), 1, ('apical_tuft', 'apical_trunk')),
            #         (20,'1', (4,5,6,8,), 1, ('apical_tuft', 'apical_trunk')),
            #         # (-20,'2', (4,5,6,8,), 1, ('apical_tuft', 'apical_trunk')),
            #     ],
            #     # subgroup
            #     'strong paired':[
            #         # trace
            #         (0, '1', (4,5,6,8,), 2, ('apical_tuft', 'apical_trunk')),
            #         (20,'1', (4,5,6,8,), 2, ('apical_tuft', 'apical_trunk')),
            #         # (-20,'2', (4,5,6,8,), 1, ('apical_tuft', 'apical_trunk')),
            #     ],
            # },
            # 'weak_paired':{
            #     'strong':[
            #         # trace
            #         (0, '2', (4,5,6,8), 2, ('apical_tuft', 'apical_trunk')),
            #         (20, '2', (4,5,6,8), 2, ('apical_tuft', 'apical_trunk')),
            #         # (-20, '2', (4,5,6,8), 2, ('apical_tuft', 'apical_trunk')),
            #     ],
            #     'strong soma':[
            #         # trace
            #         (0, '2', (4,5,6,8), 2, 'soma'),
            #         (20, '2', (4,5,6,8), 2, 'soma'),
            #         # (-20, '2', (4,5,6,8), 2, ('apical_tuft', 'apical_trunk')),
            #     ]
            # },
        }
        # load default figure parameters and colors
        #------------------------------------------
        # default     = self._default()
        default     = analysis._default_figdf()
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
        # figdf['fig_dyticks']=5
        # figdf['fig_dxticks']=10
        # figdf['fig_ylim_all']=True
        # figdf['fig_xlim_all']=True
        figdf['fig_ymin']=0
        # figdf['fig_ymax']=-40
        figdf['fig_xmin']=-2.
        figdf['fig_xmax']=2.1
        figdf['fig_ylabel']='Probability'
        figdf['fig_xlabel']='dt (soma-dendrite)'
        # figdf['fig_dyticks']=.2
        figdf['fig_dxticks']=1
        # # trace level parameters
        figdf['trace_ealpha']=0
        figdf['error_style']='shade'
        figdf['trace_linewidth']=4
        figdf['trace_linestyle']='-'
        # figdf['fig_xscale']=1./40
        # figdf['fig_barwidth']=0.8
        # figdf['fig_data_style']='point'
        figdf['fig_xtick_decimals']=0
        figdf['fig_ytick_decimals']=2
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
                if key[2][1]=='2':
                    figdf.at[key, 'trace_color']=gray
                    figdf.at[key, 'trace_ecolor']=gray
                else:
                    figdf.at[key, 'trace_color']=black
                    figdf.at[key, 'trace_ecolor']=black
                # figdf.at[key, 'trace_color']=black
                # figdf.at[key, 'trace_ecolor']=black
            # anodal
            if key[2][0]>0:
                if key[2][1]=='2':
                    figdf.at[key, 'trace_color']=red_light
                    figdf.at[key, 'trace_ecolor']=red_light
                else:
                    figdf.at[key, 'trace_color']=red
                    figdf.at[key, 'trace_ecolor']=red
                # figdf.at[key, 'trace_color']=red
                # figdf.at[key, 'trace_ecolor']=red

        return figdf

    def _spike_time_hist(self, ):
        '''
        '''
        # print progress to terminal
        #-----------------------------
        print 'building figdf:', inspect.stack()[0][3]

        # conditions for each figure
        #----------------------------
        figdict = {
            'weak_unpaired':{
                # subgroup
                'weak unpaired':[
                    # trace
                    (0, '2', (4,5,6,8,), 1, ('apical_tuft', 'apical_trunk')),
                    (20,'2', (4,5,6,8,), 1, ('apical_tuft', 'apical_trunk')),
                    # (-20,'2', (4,5,6,8,), 1, ('apical_tuft', 'apical_trunk')),
                ],
                # subgroup
                'weak paired':[
                    # trace
                    (0, '2', (4,5,6,8,), 2, ('apical_tuft', 'apical_trunk')),
                    (20,'2', (4,5,6,8,), 2, ('apical_tuft', 'apical_trunk')),
                    # (-20,'2', (4,5,6,8,), 1, ('apical_tuft', 'apical_trunk')),
                ],
            },
        }
        # load default figure parameters and colors
        #------------------------------------------
        # default     = self._default()
        default     = analysis._default_figdf()
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
        figdf['fig_dyticks']=50
        figdf['fig_dxticks']=5
        # figdf['fig_ylim_all']=True
        # figdf['fig_xlim_all']=True
        figdf['fig_ymin']=0
        figdf['fig_ymax']=300
        figdf['fig_xmin']=20.
        figdf['fig_xmax']=40.1
        figdf['fig_ylabel']='Count'
        figdf['fig_xlabel']='Spike time (ms)'
        # figdf['fig_dyticks']=.2
        # figdf['fig_dxticks']=20
        # # trace level parameters
        figdf['trace_ealpha']=0
        figdf['error_style']='shade'
        figdf['trace_linewidth']=4
        # figdf['fig_xscale']=1./40
        # figdf['fig_barwidth']=0.8
        # figdf['fig_data_style']='point'
        figdf['fig_xtick_decimals']=0
        figdf['fig_ytick_decimals']=0
        # figdf['fig_set_xscale']='symlog'
        figdf['trace_histtype']='step'



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
                if key[2][1]=='2':
                    figdf.at[key, 'trace_color']=gray
                    figdf.at[key, 'trace_ecolor']=gray
                else:
                    figdf.at[key, 'trace_color']=black
                    figdf.at[key, 'trace_ecolor']=black
                # figdf.at[key, 'trace_color']=black
                # figdf.at[key, 'trace_ecolor']=black
            # anodal
            if key[2][0]>0:
                if key[2][1]=='2':
                    figdf.at[key, 'trace_color']=red_light
                    figdf.at[key, 'trace_ecolor']=red_light
                else:
                    figdf.at[key, 'trace_color']=red
                    figdf.at[key, 'trace_ecolor']=red
                # figdf.at[key, 'trace_color']=red
                # figdf.at[key, 'trace_ecolor']=red

        figdf = self._kwargs_from_figdf(figdf=figdf, prefix='trace')

        return figdf

    def _wtrace_mean(self, ):
        '''
        '''
        # print progress to terminal
        #-----------------------------
        print 'building figdf:', #inspect.stack()[0][3]

        # conditions for each figure


        #----------------------------
        figdict = {
            '20Hz':{
                # subgroup
                '20Vm':[
                    # trace
                    ((6,8,10,12), -20, ('apical_tuft', 'apical_trunk') ),
                    ((6,8,10,12), 20, ('apical_tuft', 'apical_trunk')),
                    ((6,8,10,12), 0, ('apical_tuft', 'apical_trunk')),
                ],
            },
            '20Hz_12syn':{
                # subgroup
                '20Vm':[
                    # trace
                    ((12), -20, ('apical_tuft', 'apical_trunk') ),
                    ((12), 20, ('apical_tuft', 'apical_trunk')),
                    ((12), 0, ('apical_tuft', 'apical_trunk')),
                ],
            },
            '20Hz_8-12syn':{
                # subgroup
                '20Vm':[
                    # trace
                    ((8,10,12), -20, ('apical_tuft', 'apical_trunk') ),
                    ((8,10,12), 20, ('apical_tuft', 'apical_trunk')),
                    ((8,10,12), 0, ('apical_tuft', 'apical_trunk')),
                ],
            },
            '20Hz_6syn':{
                # subgroup
                '20Vm':[
                    # trace
                    ((6), -20, ('apical_tuft', 'apical_trunk') ),
                    ((6), 20, ('apical_tuft', 'apical_trunk')),
                    ((6), 0, ('apical_tuft', 'apical_trunk')),
                ],
            },
            '20Hz_8syn':{
                # subgroup
                '20Vm':[
                    # trace
                    ((8), -20, ('apical_tuft', 'apical_trunk') ),
                    ((8), 20, ('apical_tuft', 'apical_trunk')),
                    ((8), 0, ('apical_tuft', 'apical_trunk')),
                ],
            },
        }
        # load default figure parameters and colors
        #------------------------------------------
        # default     = self._default()
        default     = analysis._default_figdf()
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
        figdf['fig_dyticks']=.5
        figdf['fig_dxticks']=500
        # figdf['fig_ylim_all']=True
        # figdf['fig_xlim_all']=True
        # figdf['fig_ymin']=-1
        # figdf['fig_ymax']=1.5
        figdf['fig_xmin']=0.
        # figdf['fig_xmax']=40.1
        figdf['fig_ylabel']='Synaptic weight (a.u.)'
        figdf['fig_xlabel']='Time (s)'
        # figdf['fig_dyticks']=.2
        # figdf['fig_dxticks']=20
        # # trace level parameters
        figdf['trace_ealpha']=.7
        figdf['error_style']='shade'
        figdf['trace_linewidth']=4
        figdf['trace_linestyle']='-'
        # figdf['fig_yscale']=10000.
        # figdf['fig_barwidth']=0.8
        # figdf['fig_data_style']='point'
        figdf['fig_xtick_decimals']=0
        figdf['fig_ytick_decimals']=1
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
            if key[2][1]<0:
                figdf.at[key, 'trace_color']=blue
                figdf.at[key, 'trace_ecolor']=blue
            # control
            if key[2][1]==0:
                if key[2][2]=='soma':
                    figdf.at[key, 'trace_color']=gray
                    figdf.at[key, 'trace_ecolor']=gray
                else:
                    figdf.at[key, 'trace_color']=black
                    figdf.at[key, 'trace_ecolor']=black
                # figdf.at[key, 'trace_color']=black
                # figdf.at[key, 'trace_ecolor']=black
            # anodal
            if key[2][1]>0:
                if key[2][2]=='soma':
                    figdf.at[key, 'trace_color']=red_light
                    figdf.at[key, 'trace_ecolor']=red_light
                else:
                    figdf.at[key, 'trace_color']=red
                    figdf.at[key, 'trace_ecolor']=red
                # figdf.at[key, 'trace_color']=red
                # figdf.at[key, 'trace_ecolor']=red

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
            '20Hz':{
                # subgroup
                '20Vm':[
                    # trace
                    ((6,8,10,12), -20, ('apical_tuft', 'apical_trunk') ),
                    ((6,8,10,12), 20, ('apical_tuft', 'apical_trunk')),
                    ((6,8,10,12), 0, ('apical_tuft', 'apical_trunk')),
                ],
            },
            '20Hz_6syn':{
                # subgroup
                '20Vm':[
                    # trace
                    ((6), -20, ('apical_tuft', 'apical_trunk') ),
                    ((6), 20, ('apical_tuft', 'apical_trunk')),
                    ((6), 0, ('apical_tuft', 'apical_trunk')),
                ],
            },
            '20Hz_8syn':{
                # subgroup
                '20Vm':[
                    # trace
                    ((8), -20, ('apical_tuft', 'apical_trunk') ),
                    ((8), 20, ('apical_tuft', 'apical_trunk')),
                    ((8), 0, ('apical_tuft', 'apical_trunk')),
                ],
            },
            # 'acs':{
            #     # subgroup
            #     '20Vm':[
            #         # trace
            #         ('soma', (4,5,6,8,), 20, 0 ),
            #         ('soma',(4,5,6,8,), 20, np.pi),
            #         ('soma',(4,5,6,8,), 0),
            #     ],
            # },
            # 'acs all':{
            #     # subgroup
            #     '20Vm':[
            #         # trace
            #         ('soma', (4,5,6,8,9,10,11,12), 20, 0 ),
            #         ('soma',(4,5,6,8,9,10,11,12), 20, np.pi),
            #         ('soma',(4,5,6,8,9,10,11,12), 0),
            #     ],
            # },
            # 'acs_4syn':{
            #     # subgroup
            #     '20Vm':[
            #         # trace
            #         ('soma', (6), 20, 0 ),
            #         ('soma',(6), 20, np.pi),
            #         ('soma',(6), 0),
            #     ],
            # },
        }
        # load default figure parameters and colors
        #------------------------------------------
        # default     = self._default()
        default     = analysis._default_figdf()
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
        figdf['fig_dyticks']=2
        figdf['fig_dxticks']=500
        # figdf['fig_ylim_all']=True
        # figdf['fig_xlim_all']=True
        figdf['fig_ymin']=-72.2
        # figdf['fig_ymax']=-67.5
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
        figdf['trace_linestyle']='-'
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

            if key[2][1]==0:
                figdf.at[key, 'trace_location']=1
                figdf.at[key, 'trace_label']='control'
                figdf.at[key, 'trace_color']=black
                figdf.at[key, 'trace_ecolor']=black
            elif key[2][1]>0:
                figdf.at[key, 'trace_location']=2
                figdf.at[key, 'trace_label']='peak'
                figdf.at[key, 'trace_color']=red
                figdf.at[key, 'trace_ecolor']=red
            elif key[2][1]<0:
                figdf.at[key, 'trace_location']=0
                figdf.at[key, 'trace_label']='trough'
                figdf.at[key, 'trace_color']=blue
                figdf.at[key, 'trace_ecolor']=blue

        return figdf

    def _trace_individual(self, trial_id=slice(None)):
        '''
        '''
        # print progress to terminal
        #-----------------------------
        print 'building figdf:', inspect.stack()[0][3]

        # conditions for each figure
        #----------------------------
        figdict = {
            'weak_paired_dend':{
                # subgroup
                'weak':[
                    # trace
                    ((6,8,10,12), 0, ('apical_tuft', 'apical_trunk'), trial_id),
                    ((6,8,10,12), 20,('apical_tuft', 'apical_trunk'), trial_id),
                    ((6,8,10,12), -20,('apical_tuft', 'apical_trunk'), trial_id),
                ],
                # subgroup
                # 'weak soma':[
                #     # trace
                #     (0, '2', 'soma', trial_id),
                #     (20,'2', 'soma', trial_id),
                #     # (-20,'2', (4,5,6,8,), 1, ('apical_tuft', 'apical_trunk')),
                # ],
            },
            'weak_paired_soma':{
                # subgroup
                # 'weak':[
                #     # trace
                #     (0, '2', ('apical_tuft', 'apical_trunk'), trial_id),
                #     (20,'2', ('apical_tuft', 'apical_trunk'), trial_id),
                #     # (-20,'2', (4,5,6,8,), 1, ('apical_tuft', 'apical_trunk')),
                # ],
                # subgroup
                'weak soma':[
                    # trace
                    (0, '2', 'soma', trial_id),
                    (20,'2', 'soma', trial_id),
                    # (-20,'2', (4,5,6,8,), 1, ('apical_tuft', 'apical_trunk')),
                ],
            },
            # 'weak_paired':{
            #     # subgroup
            #     'weak':[
            #         # trace
            #         (0, '2', (4,5,6,8,), 2, ('apical_tuft', 'apical_trunk'), trial_id),
            #         (20,'2', (4,5,6,8,), 2, ('apical_tuft', 'apical_trunk'), trial_id),
            #         # (-20,'2', (4,5,6,8,), 1, ('apical_tuft', 'apical_trunk')),
            #     ],
            #     # subgroup
            #     'weak soma':[
            #         # trace
            #         (0, '2', (4,5,6,8,), 2, 'soma', trial_id),
            #         (20,'2', (4,5,6,8,), 2, 'soma', trial_id),
            #         # (-20,'2', (4,5,6,8,), 1, ('apical_tuft', 'apical_trunk')),
            #     ],
            # },
        }
        # load default figure parameters and colors
        #------------------------------------------
        # default     = self._default()
        default     = analysis._default_figdf()
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
        figdf['fig_dyticks']=10
        figdf['fig_dxticks']=5
        # figdf['fig_ylim_all']=True
        # figdf['fig_xlim_all']=True
        # figdf['fig_ymin']=-70
        # figdf['fig_xmin']=20.
        # figdf['fig_ymax']=-29.9
        # figdf['fig_xmax']=40.1
        figdf['fig_ylabel']='Vm (mV)'
        figdf['fig_xlabel']='Time (ms)'
        # figdf['fig_dyticks']=.2
        # figdf['fig_dxticks']=20
        # # trace level parameters
        figdf['trace_ealpha']=.7
        figdf['error_style']='shade'
        figdf['trace_linewidth']=4
        figdf['trace_linestyle']='-'
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
            # if key[2][2]=='soma':
            #     figdf.at[key, 'trace_linestyle']=(0,(1,1))
            # set colors
            #------------------------------------
            # cathodal
            if key[2][1]<0:
                figdf.at[key, 'trace_color']=blue
                figdf.at[key, 'trace_ecolor']=blue
            # control
            if key[2][1]==0:
                if key[2][2]=='soma':
                    figdf.at[key, 'trace_color']=black
                    figdf.at[key, 'trace_ecolor']=black
                else:
                    figdf.at[key, 'trace_color']=black
                    figdf.at[key, 'trace_ecolor']=black

            # anodal
            if key[2][1]>0:
                if key[2][2]=='soma':
                    figdf.at[key, 'trace_color']=red
                    figdf.at[key, 'trace_ecolor']=red
                else:
                    figdf.at[key, 'trace_color']=red
                    figdf.at[key, 'trace_ecolor']=red

        return figdf
    
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

    def _kwargs_from_figdf(self, figdf, prefix):
        ''' generate a dictionary kwargs to pass to plotting function for each trace to be plotted

        Arguments
        -------------
        ~figdf
        ~prefix: string indicating which columns from figdf to include in the kwarg dictionary.  Assumes that column names are formatted as 'trace_param', where 'param' is a keyword to passed to plotting function for each trace

        Returns
        ----------------
        ~figdf: a new column with name prefix+'_kwargs' is added.  each entry in this column is a dictionary of kwargs for plotting functions

        Comments
        -----------------

        '''
        # get list of all columns
        columns = figdf.columns
        # get columns that contain the specified prefix
        kwarg_columns = [column for column in columns if prefix in column
        ]
        # get suffixes to be keys in kwargs dictionary
        kwarg_keys = [key.split(prefix+'_')[-1] for key in kwarg_columns]
        # preallocate new column
        figdf[prefix+'_kwargs']=None
        # iterate over rows
        for row_i, row in figdf.iterrows():
            # preallocate kwarg dictionary
            kwargs={}
            # iterate over column names
            for kwarg_i, kwarg_column in enumerate(kwarg_columns):
                # get corresponding keys and values
                val = row[kwarg_column]
                key = kwarg_keys[kwarg_i]
                # store in kwarg dictionary
                kwargs[key]=val
            # update new column with kwarg dictionary
            figdf.at[row_i, prefix+'_kwargs'] = copy.deepcopy(kwargs)

        return figdf

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