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
    def _w_bar(self, ):
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
            'all_dendrite':{
                # subgroup
                'all_dendrite':[
                    # trace
                    (-20, (6,8,10,12), 'basal'),
                    (0, (6,8,10,12), 'basal'),
                    (20, (6,8,10,12), 'basal'),
                ]
            },
        }
        # load default figure parameters and colors
        #------------------------------------------
        # default     = self._default()
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
        # figdf['fig_nyticks']=2
        # figdf['fig_nxticks']=3
        # figdf['fig_dyticks']=50
        # figdf['fig_dxticks']=10
        # figdf['fig_ylim_all']=True
        # figdf['fig_xlim_all']=True
        figdf['fig_ymin']=1
        # figdf['fig_xmin']=0.
        # figdf['fig_xmax']=30.
        figdf['fig_ylabel']='Norm. weight (au)'
        figdf['fig_xlabel']=''
        # figdf['fig_dyticks']=.2
        # figdf['fig_dxticks']=20
        # # trace level parameters
        figdf['trace_ealpha']=.7
        figdf['error_style']='shade'
        figdf['trace_linewidth']=4
        # figdf['fig_yscale']=1./40
        figdf['fig_barwidth']=1
        # figdf['fig_data_style']='point'
        # figdf['fig_xtick_decimals']=1
        figdf['fig_ytick_decimals']=1
        # figdf['fig_set_xscale']='symlog'

        # subgroup parameters
        #----------------------
        figdf['sub_location']=0

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
                figdf.at[key, 'trace_label']='Cathodal'
                figdf.at[key, 'trace_color']=blue
                figdf.at[key, 'trace_ecolor']=blue
                figdf.at[key, 'trace_location']=0
            # control
            if key[2][0]==0:
                figdf.at[key, 'trace_label']='Control'
                figdf.at[key, 'trace_color']=black
                figdf.at[key, 'trace_ecolor']=black
                figdf.at[key, 'trace_location']=1
            # anodal
            if key[2][0]>0:
                figdf.at[key, 'trace_label']='Anodal'
                figdf.at[key, 'trace_color']=red
                figdf.at[key, 'trace_ecolor']=red
                figdf.at[key, 'trace_location']=2

        return figdf

    def _shapeplot(self, ):
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
            'anodal':{
                # subgroup
                'anodal':[
                    # trace
                    (20),
                ]
            },
            'cathodal':{
                # subgroup
                'cathodal':[
                    # trace
                    (-20),
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

    def _dose_response(self, ):
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
            'all':{
                # subgroup
                'all':[
                    # trace
                    (-20, ),
                    (-5, ),
                    (-1, ),
                    (-0.5, ),
                    (0, ),
                    (0.5, ),
                    (1, ),
                    (5, ),
                    (20, ),
                ]
            },
            # figure
            # -------------------------------------------------
            '6_8_10':{
                # subgroup
                '6_8_10':[
                    # trace
                    (-20, (6,8,10)),
                    (-5, (6,8,10)),
                    (-1, (6,8,10)),
                    (-0.5, (6,8,10)),
                    (0, (6,8,10)),
                    (0.5, (6,8,10)),
                    (1, (6,8,10)),
                    (5, (6,8,10)),
                    (20, (6,8,10)),
                ]
            },
            # figure
            # -------------------------------------------------
            '6_8_10_12':{
                # subgroup
                '6_8_10_12':[
                    # trace
                    (-20, (6,8,10,12)),
                    (-5, (6,8,10, 12)),
                    (-1, (6,8,10,12)),
                    (-0.5, (6,8,10,12)),
                    (0, (6,8,10,12)),
                    (0.5, (6,8,10,12)),
                    (1, (6,8,10,12)),
                    (5, (6,8,10,12)),
                    (20, (6,8,10,12)),
                ]
            },
            # figure
            # -------------------------------------------------
            '14':{
                # subgroup
                '14':[
                    # trace
                    (-20, 14),
                    (-5, 14),
                    (-1, 14),
                    (-0.5, 14),
                    (0, 14),
                    (0.5, 14),
                    (1, 14),
                    (5, 14),
                    (20, 14),
                ]
            },
            '12':{
                # subgroup
                '12':[
                    # trace
                    (-20, 12),
                    (-5, 12),
                    (-1, 12),
                    (-0.5, 12),
                    (0, 12),
                    (0.5, 12),
                    (1, 12),
                    (5, 12),
                    (20, 12),
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
        # figdf['fig_dyticks']=.2
        figdf['fig_dxticks']=np.exp(10)
        # figdf['fig_ylim_all']=True
        # figdf['fig_xlim_all']=True
        # figdf['fig_ymin']=1.2
        figdf['fig_xmin']=-30.
        figdf['fig_xmax']=30.
        figdf['fig_ylabel']='Normalized weight'
        figdf['fig_xlabel']='Electric field (V/m)'
        # figdf['fig_dyticks']=.2
        # figdf['fig_dxticks']=20
        # # trace level parameters
        # figdf['error_alpha']=1
        # figdf['error_style']='shade'
        # figdf['trace_linewidth']=4
        # figdf.at[slice(None), 'trace_ecolor'] = gray
        figdf['fig_xscale_log']=False
        figdf['fig_barwidth']=0.8
        figdf['fig_data_style']='point'
        figdf['fig_xtick_decimals']=1
        figdf['fig_ytick_decimals']=2
        figdf['fig_set_xscale']='symlog'


        # individual trace parameters
        #----------------------------
        # preallocate columns as object type
        figdf['trace_color']=None
        figdf['trace_ecolor']=None
        figdf['fig_xticks']=None
        # reset index
        figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])

        locations = [-20, -5, -1, 0, 1, 5, 20]
        locations_log = []
        for location in locations:
            if location<0:
                new_loc = -np.log(np.abs(location)+1)
            else:
                new_loc = np.log(np.abs(location)+1)
            locations_log.append(new_loc)



        # get all figure, subgroup, trace combinations
        idx_keys = figdf.index.unique().values
        # iterate over combinations
        for key in idx_keys:

            # set trace location to field magnitude
            #------------------------------------
            # figdf.at[key, 'trace_location'] = locations.index(key[2][0])
            figdf.at[key, 'trace_location'] = key[2][0]
            # figdf.at[key, 'trace_location'] = locations_log[locations.index(key[2][0])]
            # figdf.at[key, 'trace_location'] = np.log(key[2][0]+1)
            figdf.at[key, 'trace_label'] = key[2][0]
            figdf.at[key, 'fig_xticks'] = locations

            # set colors
            #------------------------------------
            # cathodal
            if key[2][0]<0:
                figdf.at[key, 'trace_color']=blue
            # control
            if key[2][0]==0:
                figdf.at[key, 'trace_color']=black
            # anodal
            if key[2][0]>0:
                figdf.at[key, 'trace_color']=red

            figdf.at[key, 'trace_ecolor']=gray





            if 'control' in key[2]:
                if key[2][0]=='weak5Hz':
                    figdf.at[key, 'trace_color']=gray
                    figdf.at[key, 'error_color']=gray
                else:
                    figdf.at[key, 'trace_color']=black
                    figdf.at[key, 'error_color']=black
            elif 'anodal' in key[2]:
                if key[2][0]=='weak5Hz':
                    figdf.at[key, 'trace_color']=red_light
                    figdf.at[key, 'error_color']=red_light
                else:
                    figdf.at[key, 'trace_color']=red
                    figdf.at[key, 'error_color']=red
            elif 'cathodal' in key[2]:
                figdf.at[key, 'trace_color']=blue
                figdf.at[key, 'error_color']=blue
            elif 'trough' in key[2]:
                figdf.at[key, 'trace_color']=red
                figdf.at[key, 'error_color']=red
            elif 'peak' in key[2]:
                figdf.at[key, 'trace_color']=blue
                figdf.at[key, 'error_color']=blue

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

    def _wtrace_mean(self, ):
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
            'all_dendrite':{
                # subgroup
                'all_dendrite':[
                    # trace
                    (-20, (4,5,6,7,8,9,10), 'basal'),
                    (0, (4,5,6,7,8,9,10), 'basal'),
                    (20, (4,5,6,7,8,9,10), 'basal'),
                ]
            },
            'low_dendrite':{
                # subgroup
                'low_dendrite':[
                    # trace
                    (-20, (4,5,6,7,8), 'basal'),
                    (0, (4,5,6,7,8), 'basal'),
                    (20, (4,5,6,7,8), 'basal'),
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
        # figdf['fig_dxticks']=10
        # # figdf['fig_ylim_all']=True
        # # figdf['fig_xlim_all']=True
        # figdf['fig_ymin']=-74
        # figdf['fig_xmin']=0.
        # figdf['fig_xmax']=30.
        figdf['fig_ylabel']='Synapse weight (au)'
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

    def _trace_individual(self, trial_id=slice(None)):
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
                    (-20, (4,5,6,7,8), 'soma', trial_id),
                    # (-5, ),
                    # (-1, ),
                    # (-0.5, ),
                    (0, (4,5,6,7,8), 'soma', trial_id),
                    # (0.5, ),
                    # (1, ),
                    # (5, ),
                    (20, (4,5,6,7,8), 'soma', trial_id),
                ]
            },
            # figure
            #-------------------------------------------------
            'all_dendrite':{
                # subgroup
                'all_dendrite':[
                    # trace
                    (-20, (6,8,10,12), 'basal' , trial_id),
                    # (-5, ),
                    # (-1, ),
                    # (-0.5, ),
                    (0, (6,8,10,12), 'basal', trial_id),
                    # (0.5, ),
                    # (1, ),
                    # (5, ),
                    (20, (6,8,10,12), 'basal', trial_id),
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
        figdf['fig_dyticks']=10
        figdf['fig_dxticks']=10
        # figdf['fig_ylim_all']=True
        # figdf['fig_xlim_all']=True
        figdf['fig_ymin']=-75
        figdf['fig_xmin']=0.
        figdf['fig_ymax']=-25.
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

    def _epsp_area(self, ):
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
            'all_dendrite':{
                # subgroup
                'all_dendrite':[
                    # trace
                    (-20, (4,5,6,7,8), 'basal'),
                    (0, (4,5,6,7,8), 'basal'),
                    (20, (4,5,6,7,8), 'basal'),
                ]
            },
        }
        # load default figure parameters and colors
        #------------------------------------------
        # default     = self._default()
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
        # figdf['fig_nyticks']=2
        # figdf['fig_nxticks']=3
        figdf['fig_dyticks']=50
        # figdf['fig_dxticks']=10
        # figdf['fig_ylim_all']=True
        # figdf['fig_xlim_all']=True
        figdf['fig_ymin']=100
        # figdf['fig_xmin']=0.
        # figdf['fig_xmax']=30.
        figdf['fig_ylabel']='EPSP area (mV*ms)'
        figdf['fig_xlabel']=''
        # figdf['fig_dyticks']=.2
        # figdf['fig_dxticks']=20
        # # trace level parameters
        figdf['trace_ealpha']=.7
        figdf['error_style']='shade'
        figdf['trace_linewidth']=4
        figdf['fig_yscale']=1./40
        figdf['fig_barwidth']=1
        # figdf['fig_data_style']='point'
        figdf['fig_xtick_decimals']=0
        figdf['fig_ytick_decimals']=0
        # figdf['fig_set_xscale']='symlog'

        # subgroup parameters
        #----------------------
        figdf['sub_location']=0

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
                figdf.at[key, 'trace_label']='Cathodal'
                figdf.at[key, 'trace_color']=blue
                figdf.at[key, 'trace_ecolor']=blue
                figdf.at[key, 'trace_location']=0
            # control
            if key[2][0]==0:
                figdf.at[key, 'trace_label']='Control'
                figdf.at[key, 'trace_color']=black
                figdf.at[key, 'trace_ecolor']=black
                figdf.at[key, 'trace_location']=1
            # anodal
            if key[2][0]>0:
                figdf.at[key, 'trace_label']='Anodal'
                figdf.at[key, 'trace_color']=red
                figdf.at[key, 'trace_ecolor']=red
                figdf.at[key, 'trace_location']=2

        return figdf

    def _hilbert_bar(self, ):
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
            'all_dendrite':{
                # subgroup
                'all_dendrite':[
                    # trace
                    (-20, (4,5,6,7,8), 'soma'),
                    (0, (4,5,6,7,8), 'soma'),
                    (20, (4,5,6,7,8), 'soma'),
                ]
            },
        }
        # load default figure parameters and colors
        #------------------------------------------
        # default     = self._default()
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
        # figdf['fig_nyticks']=2
        # figdf['fig_nxticks']=3
        # figdf['fig_dyticks']=50
        # figdf['fig_dxticks']=10
        # figdf['fig_ylim_all']=True
        # figdf['fig_xlim_all']=True
        # figdf['fig_ymin']=100
        # figdf['fig_xmin']=0.
        # figdf['fig_xmax']=30.
        figdf['fig_ylabel']='Pop. spike (mV*ms)'
        figdf['fig_xlabel']=''
        # figdf['fig_dyticks']=.2
        # figdf['fig_dxticks']=20
        # # trace level parameters
        figdf['trace_ealpha']=.7
        figdf['error_style']='shade'
        figdf['trace_linewidth']=4
        figdf['fig_yscale']=1./40
        figdf['fig_barwidth']=1
        # figdf['fig_data_style']='point'
        figdf['fig_xtick_decimals']=0
        figdf['fig_ytick_decimals']=0
        # figdf['fig_set_xscale']='symlog'

        # subgroup parameters
        #----------------------
        figdf['sub_location']=0

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
                figdf.at[key, 'trace_label']='Cathodal'
                figdf.at[key, 'trace_color']=blue
                figdf.at[key, 'trace_ecolor']=blue
                figdf.at[key, 'trace_location']=0
            # control
            if key[2][0]==0:
                figdf.at[key, 'trace_label']='Control'
                figdf.at[key, 'trace_color']=black
                figdf.at[key, 'trace_ecolor']=black
                figdf.at[key, 'trace_location']=1
            # anodal
            if key[2][0]>0:
                figdf.at[key, 'trace_label']='Anodal'
                figdf.at[key, 'trace_color']=red
                figdf.at[key, 'trace_ecolor']=red
                figdf.at[key, 'trace_location']=2

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
            'fig_xlabel_fontsize':[25],#'xx-large',#[25],
            'fig_ylabel_fontsize':[25],#'xx-large',#[25],
            'fig_xlabel_fontweight':[1000],#'extra bold',
            'fig_ylabel_fontweight':[1000],#'extra bold',
            'fig_xtick_fontsize':[20],#'large',#[15],
            'fig_ytick_fontsize':[20],#'large',#[15],
            'fig_xtick_fontweight':[1000], #'extra bold',#'light', #1000, #'heavy',
            'fig_ytick_fontweight':[1000],# 'extra bold',#'light', #1000,#'heavy',
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