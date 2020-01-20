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

    def _dose_response_dcseffect(self, ):
        '''
        '''
        # print progress to terminal
        #-----------------------------
        print 'building figdf:', inspect.stack()[0][3]

        # conditions for each figure
        #----------------------------
        figdict = {
            '6-12':{
                # subgroup
                '6-12':[
                    # trace
                    ((-20,-15,-10,-5,-2,-1 ), (6,7,8,9,10,11,12)),
                    ((20,15,10,5,2,1 ), (6,7,8,9,10,11,12)),
                    ((0,), (6,7,8,9,10,11,12)),

                ]
            },
            # '8-12':{
            #     # subgroup
            #     '8-12':[
            #         # trace
            #         (-20, (8,9,10,11,12)),
            #         (-15, (8,9,10,11,12)),
            #         (-10, (8,9,10,11,12)),
            #         (-5, (8,9,10,11,12)),
            #         (-1,(8,9,10,11,12)),
            #         (-2, (8,9,10,11,12)),
            #         (-0.5, (8,9,10,11,12)),
            #         (0, (8,9,10,11,12)),
            #         (0.5, (8,9,10,11,12)),
            #         (1, (8,9,10,11,12)),
            #         (2, (8,9,10,11,12)),
            #         (5, (8,9,10,11,12)),
            #         (10, (8,9,10,11,12)),
            #         (15, (8,9,10,11,12)),
            #         (20, (8,9,10,11,12)),
            #     ]
            # },
            # '8-12_1-20':{
            #     # subgroup
            #     '8-12':[
            #         # trace
            #         (-20, (8,9,10,11,12)),
            #         (-15, (8,9,10,11,12)),
            #         (-10, (8,9,10,11,12)),
            #         (-5, (8,9,10,11,12)),
            #         (-1,(8,9,10,11,12)),
            #         (-2, (8,9,10,11,12)),
            #         (0, (8,9,10,11,12)),
            #         (1, (8,9,10,11,12)),
            #         (2, (8,9,10,11,12)),
            #         (5, (8,9,10,11,12)),
            #         (10, (8,9,10,11,12)),
            #         (15, (8,9,10,11,12)),
            #         (20, (8,9,10,11,12)),
            #     ]
            # },
            # '6_8_10':{
            #     # subgroup
            #     '6_8_10':[
            #         # trace
            #         (-20, (6,8,10)),
            #         (-15, (6,8,10)),
            #         (-10, (6,8,10)),
            #         (-5, (6,8,10)),
            #         (-1, (6,8,10)),
            #         (-2, (6,8,10)),
            #         (-0.5, (6,8,10)),
            #         (0, (6,8,10)),
            #         (0.5, (6,8,10)),
            #         (1, (6,8,10)),
            #         (2, (6,8,10)),
            #         (5, (6,8,10)),
            #         (10, (6,8,10)),
            #         (15, (6,8,10)),
            #         (20, (6,8,10)),
            #     ]
            # },
            # # figure
            # # -------------------------------------------------
            # '6_8_10_12':{
            #     # subgroup
            #     '6_8_10_12':[
            #         # trace
            #         (-20, (6,8,10,12)),
            #         (-15, (6,8,10, 12)),
            #         (-10, (6,8,10,12)),
            #         (-5, (6,8,10, 12)),
            #         (-2, (6,8,10,12)),
            #         (-1, (6,8,10,12)),
            #         (-0.5, (6,8,10,12)),
            #         (0, (6,8,10,12)),
            #         (0.5, (6,8,10,12)),
            #         (1, (6,8,10,12)),
            #         (2, (6,8,10,12)),
            #         (5, (6,8,10,12)),
            #         (10, (6,8,10,12)),
            #         (15, (6,8,10,12)),
            #         (20, (6,8,10,12)),
            #     ]
            # },
            # # figure
            # # -------------------------------------------------
            # '6':{
            #     # subgroup
            #     '6':[
            #         # trace
            #         (-20, (6,)),
            #         (-15, (6,)),
            #         (-10, (6,)),
            #         (-5, (6,)),
            #         (-2, (6,)),
            #         (-1, (6,)),
            #         (-0.5, (6,)),
            #         (0, (6,)),
            #         (0.5, (6,)),
            #         (1, (6,)),
            #         (2, (6,)),
            #         (5, (6,)),
            #         (10, (6,)),
            #         (15, (6,)),
            #         (20, (6,)),
            #     ]
            # },
            # # figure
            # # -------------------------------------------------
            # '7':{
            #     # subgroup
            #     '7':[
            #         # trace
            #         (-20, (7,)),
            #         (-15, (7,)),
            #         (-10, (7,)),
            #         (-5, (7,)),
            #         (-2, (7,)),
            #         (-1, (7,)),
            #         (-0.5, (7,)),
            #         (0, (7,)),
            #         (0.5, (7,)),
            #         (1, (7,)),
            #         (2, (7,)),
            #         (5, (7,)),
            #         (10, (7,)),
            #         (15, (7,)),
            #         (20, (7,)),
            #     ]
            # },

            # # figure
            # # -------------------------------------------------
            # '8':{
            #     # subgroup
            #     '8':[
            #         # trace
            #         (-20, (8,)),
            #         (-15, (8,)),
            #         (-10, (8,)),
            #         (-5, (8,)),
            #         (-2, (8,)),
            #         (-1, (8,)),
            #         (-0.5, (8,)),
            #         (0, (8,)),
            #         (0.5, (8,)),
            #         (1, (8,)),
            #         (2, (8,)),
            #         (5, (8,)),
            #         (10, (8,)),
            #         (15, (8,)),
            #         (20, (8,)),
            #     ]
            # },
            # # figure
            # # -------------------------------------------------
            # '9':{
            #     # subgroup
            #     '9':[
            #         # trace
            #         (-20, (9,)),
            #         (-15, (9,)),
            #         (-10, (9,)),
            #         (-5, (9,)),
            #         (-2, (9,)),
            #         (-1, (9,)),
            #         (-0.5, (9,)),
            #         (0, (9,)),
            #         (0.5, (9,)),
            #         (1, (9,)),
            #         (2, (9,)),
            #         (5, (9,)),
            #         (10, (9,)),
            #         (15, (9,)),
            #         (20, (9,)),
            #     ]
            # },

            # # figure
            # # -------------------------------------------------
            # '10':{
            #     # subgroup
            #     '10':[
            #         # trace
            #         (-20, (10,)),
            #         (-15, (10,)),
            #         (-10, (10,)),
            #         (-5, (10,)),
            #         (-2, (10,)),
            #         (-1, (10,)),
            #         (-0.5, (10,)),
            #         (0, (10,)),
            #         (0.5, (10,)),
            #         (1, (10,)),
            #         (2, (10,)),
            #         (5, (10,)),
            #         (10, (10,)),
            #         (15, (10,)),
            #         (20, (10,)),
            #     ]
            # },
            # # figure
            # # -------------------------------------------------
            # '11':{
            #     # subgroup
            #     '11':[
            #         # trace
            #         (-20, (11,)),
            #         (-15, (11,)),
            #         (-10, (11,)),
            #         (-5, (11,)),
            #         (-2, (11,)),
            #         (-1, (11,)),
            #         (-0.5, (11,)),
            #         (0, (11,)),
            #         (0.5, (11,)),
            #         (1, (11,)),
            #         (2, (11,)),
            #         (5, (11,)),
            #         (10, (11,)),
            #         (15, (11,)),
            #         (20, (11,)),
            #     ]
            # },
            # '12':{
            #     # subgroup
            #     '12':[
            #         # trace
            #         (-20, 12),
            #         (-15, 12),
            #         (-10, 12),
            #         (-5, 12),
            #         (-1, 12),
            #         (-0.5, 12),
            #         (0, 12),
            #         (0.5, 12),
            #         (1, 12),
            #         (5, 12),
            #         (10, 12),
            #         (15, 12),
            #         (20, 12),
            #     ]
            # },
        }
        # load default figure parameters and colors
        #------------------------------------------
        # default   = self._default()
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
        # figdf['fig_dyticks']=.2
        # figdf['fig_dxticks']=np.exp(10)
        # figdf['fig_dxticks']=5
        # figdf['fig_ylim_all']=True
        # figdf['fig_xlim_all']=True
        # figdf['fig_ymin']=1.2
        # figdf['fig_xmin']=-25.
        # figdf['fig_xmax']=26.
        figdf['fig_ylabel']='Norm. weight '
        figdf['fig_xlabel']='Electric field (V/m)'
        # figdf['fig_dyticks']=.2
        # figdf['fig_dxticks']=20
        # # trace level parameters
        # figdf['error_alpha']=1
        figdf['error_style']='bar'
        figdf['trace_marker']='o'
        figdf['trace_markersize']=10
        figdf['trace_linewidth']=4

        # figdf.at[slice(None), 'trace_ecolor'] = gray
        figdf['fig_xscale_log']=False
        figdf['fig_barwidth']=0.8
        figdf['fig_data_style']='point'
        figdf['fig_xtick_decimals']=0
        figdf['fig_ytick_decimals']=2
        figdf['fig_set_xscale']='linear'


        # individual trace parameters
        #----------------------------
        # preallocate columns as object type
        figdf['trace_color']=None
        figdf['trace_ecolor']=None
        # figdf['fig_xticks']=None
        # reset index
        figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])

        # locations = [-20, -5, -1, 0, 1, 5, 20]
        # locations_log = []
        # for location in locations:
        #     if location<0:
        #         new_loc = -np.log(np.abs(location)+1)
        #     else:
        #         new_loc = np.log(np.abs(location)+1)
        #     locations_log.append(new_loc)



        # get all figure, subgroup, trace combinations
        idx_keys = figdf.index.unique().values
        # iterate over combinations
        for key in idx_keys:

            # set trace location to field magnitude
            #------------------------------------
            # figdf.at[key, 'trace_location'] = locations.index(key[2][0])
            # figdf.at[key, 'trace_location'] = key[2][0]
            # figdf.at[key, 'trace_location'] = locations_log[locations.index(key[2][0])]
            # figdf.at[key, 'trace_location'] = np.log(key[2][0]+1)
            # figdf.at[key, 'trace_label'] = key[2][0]
            # figdf.at[key, 'fig_xticks'] = locations

            # set colors
            #------------------------------------
            print key[2][0]
            # cathodal
            if key[2][0][0]<0:
                figdf.at[key, 'trace_color']=blue
            # control
            if key[2][0][0]==0:
                figdf.at[key, 'trace_color']=black
            # anodal
            if key[2][0][0]>0:
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

    def _dose_response(self, ):
        '''
        '''
        # print progress to terminal
        #-----------------------------
        print 'building figdf:', inspect.stack()[0][3]

        # conditions for each figure
        #----------------------------
        figdict = {
            '6-12':{
                # subgroup
                '6-12':[
                    # trace
                    (-20, (6,7,8,9,10,11,12)),
                    (-15, (6,7,8,9,10,11,12)),
                    (-10, (6,7,8,9,10,11,12)),
                    (-5, (6,7,8,9,10,11,12)),
                    (-1,(6,7,8,9,10,11,12)),
                    (-2, (6,7,8,9,10,11,12)),
                    (-0.5, (6,7,8,9,10,11,12)),
                    (0, (6,7,8,9,10,11,12)),
                    (0.5, (6,7,8,9,10,11,12)),
                    (1, (6,7,8,9,10,11,12)),
                    (2, (6,7,8,9,10,11,12)),
                    (5, (6,7,8,9,10,11,12)),
                    (10, (6,7,8,9,10,11,12)),
                    (15, (6,7,8,9,10,11,12)),
                    (20, (6,7,8,9,10,11,12)),
                ]
            },
            '8-12':{
                # subgroup
                '8-12':[
                    # trace
                    (-20, (8,9,10,11,12)),
                    (-15, (8,9,10,11,12)),
                    (-10, (8,9,10,11,12)),
                    (-5, (8,9,10,11,12)),
                    (-1,(8,9,10,11,12)),
                    (-2, (8,9,10,11,12)),
                    (-0.5, (8,9,10,11,12)),
                    (0, (8,9,10,11,12)),
                    (0.5, (8,9,10,11,12)),
                    (1, (8,9,10,11,12)),
                    (2, (8,9,10,11,12)),
                    (5, (8,9,10,11,12)),
                    (10, (8,9,10,11,12)),
                    (15, (8,9,10,11,12)),
                    (20, (8,9,10,11,12)),
                ]
            },
            '8-12_1-20':{
                # subgroup
                '8-12':[
                    # trace
                    (-20, (8,9,10,11,12)),
                    (-15, (8,9,10,11,12)),
                    (-10, (8,9,10,11,12)),
                    (-5, (8,9,10,11,12)),
                    (-1,(8,9,10,11,12)),
                    (-2, (8,9,10,11,12)),
                    (0, (8,9,10,11,12)),
                    (1, (8,9,10,11,12)),
                    (2, (8,9,10,11,12)),
                    (5, (8,9,10,11,12)),
                    (10, (8,9,10,11,12)),
                    (15, (8,9,10,11,12)),
                    (20, (8,9,10,11,12)),
                ]
            },
            '6_8_10':{
                # subgroup
                '6_8_10':[
                    # trace
                    (-20, (6,8,10)),
                    (-15, (6,8,10)),
                    (-10, (6,8,10)),
                    (-5, (6,8,10)),
                    (-1, (6,8,10)),
                    (-2, (6,8,10)),
                    (-0.5, (6,8,10)),
                    (0, (6,8,10)),
                    (0.5, (6,8,10)),
                    (1, (6,8,10)),
                    (2, (6,8,10)),
                    (5, (6,8,10)),
                    (10, (6,8,10)),
                    (15, (6,8,10)),
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
                    (-15, (6,8,10, 12)),
                    (-10, (6,8,10,12)),
                    (-5, (6,8,10, 12)),
                    (-2, (6,8,10,12)),
                    (-1, (6,8,10,12)),
                    (-0.5, (6,8,10,12)),
                    (0, (6,8,10,12)),
                    (0.5, (6,8,10,12)),
                    (1, (6,8,10,12)),
                    (2, (6,8,10,12)),
                    (5, (6,8,10,12)),
                    (10, (6,8,10,12)),
                    (15, (6,8,10,12)),
                    (20, (6,8,10,12)),
                ]
            },
            # figure
            # -------------------------------------------------
            '6':{
                # subgroup
                '6':[
                    # trace
                    (-20, (6,)),
                    (-15, (6,)),
                    (-10, (6,)),
                    (-5, (6,)),
                    (-2, (6,)),
                    (-1, (6,)),
                    (-0.5, (6,)),
                    (0, (6,)),
                    (0.5, (6,)),
                    (1, (6,)),
                    (2, (6,)),
                    (5, (6,)),
                    (10, (6,)),
                    (15, (6,)),
                    (20, (6,)),
                ]
            },
            # figure
            # -------------------------------------------------
            '7':{
                # subgroup
                '7':[
                    # trace
                    (-20, (7,)),
                    (-15, (7,)),
                    (-10, (7,)),
                    (-5, (7,)),
                    (-2, (7,)),
                    (-1, (7,)),
                    (-0.5, (7,)),
                    (0, (7,)),
                    (0.5, (7,)),
                    (1, (7,)),
                    (2, (7,)),
                    (5, (7,)),
                    (10, (7,)),
                    (15, (7,)),
                    (20, (7,)),
                ]
            },

            # figure
            # -------------------------------------------------
            '8':{
                # subgroup
                '8':[
                    # trace
                    (-20, (8,)),
                    (-15, (8,)),
                    (-10, (8,)),
                    (-5, (8,)),
                    (-2, (8,)),
                    (-1, (8,)),
                    (-0.5, (8,)),
                    (0, (8,)),
                    (0.5, (8,)),
                    (1, (8,)),
                    (2, (8,)),
                    (5, (8,)),
                    (10, (8,)),
                    (15, (8,)),
                    (20, (8,)),
                ]
            },
            # figure
            # -------------------------------------------------
            '9':{
                # subgroup
                '9':[
                    # trace
                    (-20, (9,)),
                    (-15, (9,)),
                    (-10, (9,)),
                    (-5, (9,)),
                    (-2, (9,)),
                    (-1, (9,)),
                    (-0.5, (9,)),
                    (0, (9,)),
                    (0.5, (9,)),
                    (1, (9,)),
                    (2, (9,)),
                    (5, (9,)),
                    (10, (9,)),
                    (15, (9,)),
                    (20, (9,)),
                ]
            },

            # figure
            # -------------------------------------------------
            '10':{
                # subgroup
                '10':[
                    # trace
                    (-20, (10,)),
                    (-15, (10,)),
                    (-10, (10,)),
                    (-5, (10,)),
                    (-2, (10,)),
                    (-1, (10,)),
                    (-0.5, (10,)),
                    (0, (10,)),
                    (0.5, (10,)),
                    (1, (10,)),
                    (2, (10,)),
                    (5, (10,)),
                    (10, (10,)),
                    (15, (10,)),
                    (20, (10,)),
                ]
            },
            # figure
            # -------------------------------------------------
            '11':{
                # subgroup
                '11':[
                    # trace
                    (-20, (11,)),
                    (-15, (11,)),
                    (-10, (11,)),
                    (-5, (11,)),
                    (-2, (11,)),
                    (-1, (11,)),
                    (-0.5, (11,)),
                    (0, (11,)),
                    (0.5, (11,)),
                    (1, (11,)),
                    (2, (11,)),
                    (5, (11,)),
                    (10, (11,)),
                    (15, (11,)),
                    (20, (11,)),
                ]
            },
            '12':{
                # subgroup
                '12':[
                    # trace
                    (-20, 12),
                    (-15, 12),
                    (-10, 12),
                    (-5, 12),
                    (-1, 12),
                    (-0.5, 12),
                    (0, 12),
                    (0.5, 12),
                    (1, 12),
                    (5, 12),
                    (10, 12),
                    (15, 12),
                    (20, 12),
                ]
            },
        }
        # load default figure parameters and colors
        #------------------------------------------
        # default   = self._default()
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
        # figdf['fig_dyticks']=.2
        # figdf['fig_dxticks']=np.exp(10)
        figdf['fig_dxticks']=5
        # figdf['fig_ylim_all']=True
        # figdf['fig_xlim_all']=True
        # figdf['fig_ymin']=1.2
        figdf['fig_xmin']=-25.
        figdf['fig_xmax']=26.
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
        figdf['fig_xtick_decimals']=0
        figdf['fig_ytick_decimals']=2
        figdf['fig_set_xscale']='linear'


        # individual trace parameters
        #----------------------------
        # preallocate columns as object type
        figdf['trace_color']=None
        figdf['trace_ecolor']=None
        # figdf['fig_xticks']=None
        # reset index
        figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])

        # locations = [-20, -5, -1, 0, 1, 5, 20]
        # locations_log = []
        # for location in locations:
        #     if location<0:
        #         new_loc = -np.log(np.abs(location)+1)
        #     else:
        #         new_loc = np.log(np.abs(location)+1)
        #     locations_log.append(new_loc)



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
            # figdf.at[key, 'fig_xticks'] = locations

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

    def _polarization_dose_response(self, ):
        '''
        '''
        # print progress to terminal
        #-----------------------------
        print 'building figdf:', inspect.stack()[0][3]

        # conditions for each figure
        #----------------------------
        figdict = {
            'soma_apical':{
                # subgroup
                'soma_apical':[
                    # trace
                    (('soma'),),
                    (('apical_tuft'),),
                    # (('basal'),),
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
        figdf['trace_markersize']=20
        figdf['trace_ealpha']=1
        figdf['error_style']='shade'
        figdf['trace_linewidth']=6
        figdf['error_linewidth']=4
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
                figdf.at[key, 'trace_marker']='o'
            if 'apical_tuft' in key[2]:
                figdf.at[key, 'trace_color']= gray #(1,0,1)#gray
                figdf.at[key, 'trace_ecolor']=gray#(1,0,1)#gray
                figdf.at[key, 'trace_marker']='x'#gray
            if 'basal' in key[2]:
                figdf.at[key, 'trace_color']=(0,1,1)#gray
                figdf.at[key, 'trace_ecolor']=(0,1,1)#gray
                figdf.at[key, 'trace_marker']='D'

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
                    (-20, (6,8,10,12), 'soma'),
                    # (-5, ),
                    # (-1, ),
                    # (-0.5, ),
                    (0, (6,8,10,12), 'soma'),
                    # (0.5, ),
                    # (1, ),
                    # (5, ),
                    (20, (6,8,10,12), 'soma'),
                ]
            },
            'all_axon':{
                # subgroup
                'all_axon':[
                    # trace
                    (-20, (6,8,10,12), 'axon'),
                    # (-5, ),
                    # (-1, ),
                    # (-0.5, ),
                    (0, (6,8,10,12), 'axon'),
                    # (0.5, ),
                    # (1, ),
                    # (5, ),
                    (20, (6,8,10,12), 'axon'),
                ]
            },
            # figure
            #-------------------------------------------------
            'all_dendrite':{
                # subgroup
                'all_dendrite':[
                    # trace
                    (-20, (6,8,10,12), 'apical_tuft'),
                    # (-5, ),
                    # (-1, ),
                    # (-0.5, ),
                    (0, (6,8,10,12), 'apical_tuft'),
                    # (0.5, ),
                    # (1, ),
                    # (5, ),
                    (20, (6,8,10,12), 'apical_tuft'),
                ]
            },
            # figure
            #-------------------------------------------------
            '14_dendrite':{
                # subgroup
                '14_dendrite':[
                    # trace
                    (-20, (14), 'apical_tuft'),
                    # (-5, ),
                    # (-1, ),
                    # (-0.5, ),
                    (0, (14), 'apical_tuft'),
                    # (0.5, ),
                    # (1, ),
                    # (5, ),
                    (20, (14), 'apical_tuft'),
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
            'all_soma_20':{
                # subgroup
                'all_soma':[
                    # trace
                    (-20, 'soma', trial_id),
                    (0, 'soma', trial_id),
                    (20, 'soma', trial_id),
                ]
            },
            # figure
            #-------------------------------------------------
            'all_dendrite_20':{
                # subgroup
                'all_dendrite':[
                    # trace
                    (-20, 'apical_tuft', trial_id),
                    (0, 'apical_tuft', trial_id),
                    (20, 'apical_tuft', trial_id),
                ]
            },
            'all_soma_1':{
                # subgroup
                'all_soma':[
                    # trace
                    (-1, 'soma', trial_id),
                    (0, 'soma', trial_id),
                    (1, 'soma', trial_id),
                ]
            },
            # figure
            #-------------------------------------------------
            'all_dendrite_1':{
                # subgroup
                'all_dendrite':[
                    # trace
                    (-1, 'apical_tuft', trial_id),
                    (0, 'apical_tuft', trial_id),
                    (1, 'apical_tuft', trial_id),
                ]
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
                    (-20, (6,7,8,9,10,11,12), ('apical_tuft','apical_trunk')),
                    (0, (6,7,8,9,10,11,12), ('apical_tuft','apical_trunk')),
                    (20, (6,7,8,9,10,11,12), ('apical_tuft','apical_trunk')),
                ]
            },
            'all_dendrite_1':{
                # subgroup
                'all_dendrite_1':[
                    # trace
                    (-1, (6,7,8,9,10,11,12), ('apical_tuft','apical_trunk')),
                    (0, (6,7,8,9,10,11,12), ('apical_tuft','apical_trunk')),
                    (1, (6,7,8,9,10,11,12), ('apical_tuft','apical_trunk')),
                ]
            },
            '8-12syn 1vm':{
                # subgroup
                '8-12syn 1vm':[
                    # trace
                    (-1, (8,9,10,11,12), ('apical_tuft','apical_trunk')),
                    (0, (8,9,10,11,12), ('apical_tuft','apical_trunk')),
                    (1, (8,9,10,11,12), ('apical_tuft','apical_trunk')),
                ]
            },
            '8-12syn 2vm':{
                # subgroup
                '8-12syn 2vm':[
                    # trace
                    (-2, (8,9,10,11,12), ('apical_tuft','apical_trunk')),
                    (0, (8,9,10,11,12), ('apical_tuft','apical_trunk')),
                    (2, (8,9,10,11,12), ('apical_tuft','apical_trunk')),
                ]
            },
            '8-12syn 5vm':{
                # subgroup
                '8-12syn 5vm':[
                    # trace
                    (-5, (8,9,10,11,12), ('apical_tuft','apical_trunk')),
                    (0, (8,9,10,11,12), ('apical_tuft','apical_trunk')),
                    (5, (8,9,10,11,12), ('apical_tuft','apical_trunk')),
                ]
            },
            '8-12syn 20vm':{
                # subgroup
                '8-12syn 20vm':[
                    # trace
                    (-20, (8,9,10,11,12), ('apical_tuft','apical_trunk')),
                    (0, (8,9,10,11,12), ('apical_tuft','apical_trunk')),
                    (20, (8,9,10,11,12), ('apical_tuft','apical_trunk')),
                ]
            },
            '12syn 20vm':{
                # subgroup
                '12syn 20vm':[
                    # trace
                    (-20, (12), ('apical_tuft','apical_trunk')),
                    (0, (12), ('apical_tuft','apical_trunk')),
                    (20, (12), ('apical_tuft','apical_trunk')),
                ]
            },
            '12syn 1vm':{
                # subgroup
                '12syn 1vm':[
                    # trace
                    (-1, (12), ('apical_tuft','apical_trunk')),
                    (0, (12), ('apical_tuft','apical_trunk')),
                    (1, (12), ('apical_tuft','apical_trunk')),
                ]
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
        # figdf['fig_dyticks']=4
        # figdf['fig_dxticks']=10
        # figdf['fig_ylim_all']=True
        # figdf['fig_xlim_all']=True
        # figdf['fig_ymin']=-74
        # figdf['fig_xmin']=0.
        # figdf['fig_xmax']=30.
        figdf['fig_ylabel']='Probability'
        figdf['fig_xlabel']='$\Delta W_{DCS}$'
        # figdf['fig_dyticks']=.2
        # figdf['fig_dxticks']=20
        # # trace level parameters
        figdf['trace_ealpha']=.7
        figdf['error_style']='shade'
        figdf['trace_linewidth']=3
        # figdf['fig_xscale']=1./40
        # figdf['fig_barwidth']=0.8
        # figdf['fig_data_style']='point'
        figdf['fig_xtick_decimals']=1
        figdf['fig_ytick_decimals']=3
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


