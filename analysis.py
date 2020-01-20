"""
analysis

"""
from __future__ import division
import numpy as np
from scipy import stats
from scipy import signal
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import itertools 
import os
import glob
import cPickle as pickle
import param
import math
# import run_control
import copy
import matplotlib
import matplotlib.patches as patches
import matplotlib.lines as mlines
from matplotlib.patches import Polygon
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib import cm as colormap
import inspect
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as ticker
import pdb
import sys
from sys import getsizeof, stderr
from collections import deque

# matplotlib.rc('text', usetex=True)
# matplotlib.rcParams['text.latex.preamble']=[r'\boldmath']
# try:
#     from reprlib import repr
# except ImportError:
#     pass


#############################################################################
# common functions
#############################################################################
def _default_figdf():
        '''
        '''
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

def _input_times_2array(input_times, t, t_precision=4, **kwargs):
    '''
    Arguments
    ---------
    ~input_times: nested list of input times [location number][list of input times]
    '''
    # number of locations with input times
    n_locs = len(input_times)
    # number of time steps
    n_steps = len(t)
    # convert input times to input array
    #--------------------------------------------------------------
    # preallocate
    input_array = np.zeros((n_locs, n_steps))
    # cut off number of decimals to prevent rounding error
    t = np.around(t, t_precision)
    # iterate over locations
    for loc_i, times in enumerate(input_times):
        # print type(times)
        # print np.isnan(times)
        if type(times)==float and np.isnan(times):
            continue
        # find input times and to boolean vector
        # iterate over input times
        for t_i, t_t in enumerate(times):
            # if input time does not exactly match a sample in time vector    
            if t_t not in t:
                # set input time to nearest sample
                input_array[loc_i, int(np.argmin(np.absolute(t-t_t)))]=1
            else:
                # store input times
                input_array[loc_i, np.where(t==t_t)]=1

    return input_array

def _format_figures(fig, ax, figdf, tick_max=40):
    '''
    '''
    # get ylim and xlim across all figures
    #-------------------------------------
    xlim={}
    ylim={}
    xlims=[]
    ylims=[]
    # iterate over figures
    for figkey in ax:
        # get ticks
        yticks = ax[figkey].get_yticks()
        xticks = ax[figkey].get_xticks()
        # apply scaling to ticks
        #-------------------------------------
        if 'fig_yscale' in figdf:
            yscale = figdf.loc[figkey].fig_yscale.unique()[0]
            print 'yscale', yscale
            yticks = yscale*yticks
        if 'fig_xscale' in figdf:
            xscale = figdf.loc[figkey].fig_xscale.unique()[0]
            xticks = xscale*xticks
        ax[figkey].set_yticks(yticks)
        ax[figkey].set_xticks(xticks)
        # get updated lims
        xlim[figkey] = _to_1sigfig(ax[figkey].get_xlim())
        ylim[figkey] = _to_1sigfig(ax[figkey].get_ylim())
        # print ylim
        xlims.append(copy.copy(xlim[figkey]))
        ylims.append(copy.copy(ylim[figkey]))
    # find x and y lims across all figures in group
    xlim_all = [min([temp[0] for temp in xlims]), max([temp[1] for temp in xlims])]
    ylim_all = [min([temp[0] for temp in ylims]), max([temp[1] for temp in ylims])]

    # iterate over figures and update attributes
    #---------------------------------------------
    for figkey, axes in ax.iteritems():
        # get x and y lims
        ylim_current = list(ax[figkey].get_ylim())
        xlim_current = list(ax[figkey].get_xlim())
        # set ylim across all figures
        #----------------------------
        if all(figdf.fig_ylim_all):
            xlim_current = xlim_all
            ylim_current = ylim_all
        if 'fig_ymin' in figdf.keys():
            ylim_current[0] = figdf.loc[figkey].fig_ymin.unique()[0]
        if 'fig_xmin' in figdf.keys():
            xlim_current[0] = figdf.loc[figkey].fig_xmin.unique()[0]
        if 'fig_ymax' in figdf.keys():
            ylim_current[1] = figdf.loc[figkey].fig_ymax.unique()[0]
        if 'fig_xmax' in figdf.keys():
            xlim_current[1] = figdf.loc[figkey].fig_xmax.unique()[0]

        ax[figkey].set_ylim(ylim_current)
        ax[figkey].set_xlim(xlim_current)

        # set x and y ticks
        #------------------------------------------------------------------
        print 'setting ticks'
        # get current lims and ticks
        xlim[figkey] = copy.copy(list(ax[figkey].get_xlim()))
        ylim[figkey] = copy.copy(list(ax[figkey].get_ylim()))
        xticks = ax[figkey].get_xticks()
        yticks = ax[figkey].get_yticks()
        nyticks=5
        dyticks = float(abs(ylim[figkey][1]-ylim[figkey][0]))/nyticks
        nxticks=5
        dxticks = float(abs(xlim[figkey][1]-xlim[figkey][0]))/nxticks
        # print 'ylim', ylim[figkey]
        if 'fig_dyticks' in figdf:
            dyticks = figdf.loc[figkey].fig_dyticks.unique()[0]
            nyticks = len(np.arange(ylim[figkey][0], ylim[figkey][1], dyticks))
        elif 'fig_nyticks' in figdf:
            nyticks = figdf.loc[figkey].fig_nyticks.unique()[0]
            dyticks = float(abs(ylim[figkey][1]-ylim[figkey][0]))/nyticks
        # else:
        #     nyticks=5
        #     dyticks = float(abs(ylim[figkey][1]-ylim[figkey][0]))/nyticks

        if nyticks>tick_max:
            nyticks=tick_max
            dyticks = float(abs(ylim[figkey][1]-ylim[figkey][0]))/nyticks
        if 'fig_ytick_assert' in figdf:
            ytick_assert = figdf.loc[figkey].fig_ytick_assert.unique()[0]
        else:
            ytick_assert = ylim[figkey][0]
        # print ylim[figkey]
        # print ytick_assert
        yticks_new_1 = np.flip(np.arange(ytick_assert, ylim[figkey][0], -dyticks))
        # print yticks_new_1
        yticks_new_2 = np.arange(ytick_assert, ylim[figkey][1], dyticks)
        yticks_new = np.append(yticks_new_1, yticks_new_2)

        if 'fig_dxticks' in figdf:
            dxticks = figdf.loc[figkey].fig_dxticks.unique()[0]
            nxticks = len(np.arange(xlim[figkey][0], xlim[figkey][1], dxticks))
        elif 'fig_nxticks' in figdf:
            nxticks = figdf.loc[figkey].fig_nxticks.unique()[0]
            dxticks = float(abs(xlim[figkey][1]-xlim[figkey][0]))/nxticks
        # else:
        #     nxticks=5
        #     dxticks = float(abs(xlim[figkey][1]-xlim[figkey][0]))/nxticks
        if nxticks>tick_max:
            nxticks=tick_max
            dxticks = float(abs(xlim[figkey][1]-xlim[figkey][0]))/nxticks
        if 'fig_xtick_assert' in figdf:
            xtick_assert = figdf.loc[figkey].fig_xtick_assert.unique()[0]
        else:
            xtick_assert = xlim[figkey][0]
            xtick_assert = _to_1sigfig(xtick_assert)
        xticks_new_1 = np.flip(np.arange(xtick_assert, xlim[figkey][0], -dxticks))
        xticks_new_2 = np.arange(xtick_assert, xlim[figkey][1], dxticks)
        xticks_new = np.append(xticks_new_1, xticks_new_2)
        print xticks_new_1, xticks_new_2

        if 'fig_ytick_round' in figdf:
            decimals = figdf.loc[figkey].fig_ytick_round.unique()[0]
            yticks_new = np.round(yticks_new, decimals=decimals)

        if 'fig_xtick_round' in figdf:
            decimals = figdf.loc[figkey].fig_xtick_round.unique()[0]
            xticks_new = np.round(xticks_new, decimals=decimals)

        if 'fig_xticks' in figdf:
            xticks_new = figdf.fig_xticks.values[0]
        if 'fig_yticks' in figdf:
            yticks_new = figdf.fig_yticks.values[0]

        if 'fig_set_xscale' in figdf:
            ax[figkey].set_xscale(figdf.fig_set_xscale.values[0])
        print 'yticks_new',yticks_new
        print 'xticks_new',xticks_new

        # if 'fig_yscale' in figdf:
        #     yscale = figdf.loc[figkey].fig_yscale.unique()[0]
        #     print 'yscale', yscale
        #     yticks_new = yscale*yticks_new
        # if 'fig_xscale' in figdf:
        #     xscale = figdf.loc[figkey].fig_xscale.unique()[0]
        #     xticks_new = xscale*xticks_new
        ax[figkey].set_yticks(yticks_new)
        ax[figkey].set_xticks(xticks_new)
        ax[figkey].set_yticklabels(yticks_new)
        ax[figkey].set_xticklabels(xticks_new)

        # set tick decimal places
        #---------------------------------------------------------------
        if 'fig_xtick_decimals' in figdf:
            dec = figdf.loc[figkey].fig_xtick_decimals.unique()[0]
            ax[figkey].xaxis.set_major_formatter(FormatStrFormatter('%.{}f'.format(dec)))
        if 'fig_ytick_decimals' in figdf:
            dec = figdf.loc[figkey].fig_ytick_decimals.unique()[0]
            ax[figkey].yaxis.set_major_formatter(FormatStrFormatter('%.{}f'.format(dec)))




        #     yticks = np.arange(ylim[figkey][0], ylim[figkey][1], figdf.loc[figkey].fig_dyticks.unique()[0])
        # elif 'fig_nyticks' in figdf:
        #     dstep = float(abs(ylim[figkey][1]-ylim[figkey][0]))/figdf.loc[figkey].fig_nyticks.unique()[0]
        #     yticks = np.arange(ylim[figkey][0], ylim[figkey][1], dstep)

        # if 'fig_dxticks' in figdf:
        #     xticks = np.arange(xlim[figkey][0], xlim[figkey][1], figdf.loc[figkey].fig_dxticks.unique()[0])
        # elif 'fig_nxticks' in figdf:
        #     dstep = float(abs(xlim[figkey][1]-xlim[figkey][0]))/figdf.loc[figkey].fig_nxticks.unique()[0]
        #     xticks = np.arange(xlim[figkey][0], xlim[figkey][1], dstep)
        # print 'resetting ticks'
        # print xticks
        # # print len(yticks)
        # # print xticks
        # ax[figkey].set_yticks(yticks)
        # ax[figkey].set_xticks(xticks)

        

        # set ticklabel attributes
        #--------------------------
        print 'setting ticklabels'
        for temp in ax[figkey].get_xticklabels():
            temp.set_fontweight(figdf.loc[figkey].fig_xtick_fontweight.unique()[0])
            temp.set_fontsize(figdf.loc[figkey].fig_xtick_fontsize.unique()[0])
        for temp in ax[figkey].get_yticklabels():
            temp.set_fontweight(figdf.loc[figkey].fig_ytick_fontweight.unique()[0])
            temp.set_fontsize(figdf.loc[figkey].fig_ytick_fontsize.unique()[0])

        # turn off axes box
        #------------------
        if all(figdf.loc[figkey].fig_boxoff):
            ax[figkey].spines['right'].set_visible(False)
            ax[figkey].spines['top'].set_visible(False)

        if 'fig_axesoff' in figdf.keys() and all(figdf.loc[figkey].fig_axesoff):
            ax[figkey].set_axis_off()

        # set axes linewidth and tick position
        #----------------------
        ax[figkey].spines['left'].set_linewidth(figdf.loc[figkey].fig_axes_linewidth.unique()[0])
        ax[figkey].spines['bottom'].set_linewidth(figdf.loc[figkey].fig_axes_linewidth.unique()[0])
        ax[figkey].xaxis.set_ticks_position('bottom')
        ax[figkey].yaxis.set_ticks_position('left')

        # set axes labels
        #----------------

        # print 'fontsize',figdf['fig_xlabel_fonstize']
        # print 'fontsize',figdf.loc[figkey]['fig_xlabel_fonstize']
        # print 'fontsize',figdf.loc[figkey]['fig_xlabel_fonstize'].unique()[0]
        # print 'fontsize',figdf.loc[figkey].fig_xlabel_fontsize.unique()
        ax[figkey].set_xlabel(figdf.loc[figkey].fig_xlabel.unique()[0], fontsize=figdf.loc[figkey]['fig_xlabel_fontsize'].unique()[0], fontweight=figdf.loc[figkey].fig_xlabel_fontweight.unique()[0])
        ax[figkey].set_ylabel(figdf.loc[figkey].fig_ylabel.unique()[0], fontsize=figdf.loc[figkey]['fig_ylabel_fontsize'].unique()[0], fontweight=figdf.loc[figkey].fig_ylabel_fontweight.unique()[0])

        # set tight layout
        #-------------------
        if all(figdf.loc[figkey].fig_tight_layout):
            plt.figure(fig[figkey].number)
            plt.tight_layout()

    return fig, ax

def _load_group_data( directory='', filename='slopes_df.pkl', df=True):
    """ Load group data from folder
    
    ===Args===
    -directory : directory where group data is stored including /
    -filename : file name for group data file, including .pkl
                -file_name cannot contain the string 'data', since this string is used t search for individual data files
    -df : boolean, if true group data is assumed to be a pandas dataframe, otherwise it is assumed to be a nested dictionary

    ===Out===
    -group_data  : if df is True, group_data will be a pandas dataframe.  if no group data is found in the directory an empty dataframe is returned.  if df is False, group_data is returned as a nested dictionary, and if no data is found an empty dictionary is returned

    ===Updates===
    -none

    ===Comments===
    """
    
    # all files in directory
    # files = os.listdir(directory)
    files = glob.glob(directory+'*'+filename+'*')
    print 'files',files

    if len(files)>0:
        print len(files), 'variable files found'
        group_data=pd.DataFrame()
        for file in files:
            group_data_temp = pd.read_pickle(file)
            print 'group data loaded', file
            # print group_data_temp
            if group_data.empty:
                group_data=group_data_temp
            else:
                group_data=group_data.append(group_data_temp)

        # print 'keys',sorted(group_data.keys())

    # otherwise create data structure
    else:
        print 'no group data found'
        
        # if dataframe
        if df:
            # create empty dataframe
            group_data = pd.DataFrame()
        # else if dicitonary
        else:
            group_data= {}


    # # if data file already exists
    # if filename in files:
    #     print 'group data found:', filename

    #     # if stored as dataframe
    #     if df:
    #         if '.pkl' in filename:
    #             # load dataframe
    #             group_data=pd.read_pickle(directory+filename)
    #             print 'group data loaded'

    #         elif '.h5' in filename:
    #             # load dataframe
    #             group_data=pd.read_hdf(directory+filename)
    #             print 'group data loaded'
    #     # if stored as dictionary
    #     else:
    #         # load dictionary
    #         with open(directory+filename, 'rb') as pkl_file:
    #             group_data= pickle.load(pkl_file)
    #         print 'group data loaded'

    # # otherwise create data structure
    # else:
    #     print 'no group data found'
        
    #     # if dataframe
    #     if df:
    #         # create empty dataframe
    #         group_data = pd.DataFrame()
    #     # else if dicitonary
    #     else:
    #         group_data= {}

    return group_data 

def _save_group_data(df, directory, variable, extension='.pkl', size_limit=5E7, check_against=None, nsplit=10, **kwargs):
    '''
    '''
    if (check_against is not None and not df.equals(check_against)) or check_against is None:
        print 'saving updated group df'
        # get size of df
        #---------------------------
        df_size_temp = sys.getsizeof(df)
        df_size = _total_size(df)
        print 'file size', df_size_temp
        print 'total size', df_size
        # if df is too big, break into manageable parts
        if nsplit is not None:
            print 'splitting group df'
            n_dfs=nsplit
            df_split = np.array_split(df, n_dfs)
            for df_i, df_chunk in enumerate(df_split):
                if extension=='.pkl':
                    # print df_chunk
                    df_chunk.to_pickle(directory+variable+'_'+str(df_i)+extension)
                elif extension=='.h5':
                    df_chunk.to_hdf(directory+variable+'_'+str(df_i)+extension, key=variable, mode='w', format='fixed')
        elif size_limit is not None and df_size>size_limit:
            print 'df too large to save, breaking into smaller parts'
            n_dfs=df_size/size_limit
            df_split = np.array_split(df, n_dfs)
            for df_i, df_chunk in enumerate(df_split):
                if extension=='.pkl':
                    # print df_chunk
                    df_chunk.to_pickle(directory+variable+'_'+str(df_i)+extension)
                elif extension=='.h5':
                    df_chunk.to_hdf(directory+variable+'_'+str(df_i)+extension, key=variable, mode='w', format='fixed')

        else:
            # check write protocol
            if extension=='.pkl':
                df.to_pickle(directory+variable+extension)
            elif extension=='.h5':
                df.to_hdf(directory+variable+extension, key=variable, mode='w', format='fixed')

    return df

def _to_1sigfig(num, minimum=True):
    '''
    '''
    print 'num',num
    if type(num)==np.ndarray or type(num)==list or type(num)==tuple:
        new_val = np.zeros(len(num))
        sigfigs=[]
        for i, val in enumerate(num):
            if val!=0:
                sigfigs.append(int(math.floor(math.log10(abs(val)))))
            else:
                sigfigs.append(0)

        min_sigfig = min(sigfigs)
        print min_sigfig
        for i, val in enumerate(num):
            if minimum:
                new_val[i] = np.round(val, -min_sigfig)
            else:
                new_val[i] = np.round(val, -sigfigs[i])

    elif type(num)==float or type(num)==np.float64:
        if num!=0:
            new_val = np.round(num, -int(math.floor(math.log10(abs(num)))))
        else:
            new_val=0

    print new_val
    return new_val

def _2array(series, remove_nans=True, remove_nans_axis=0, list_index=0, array_funcs=[], array_func_kws=[]):
    '''
    '''
    if remove_nans:
        series = series.dropna()
    series_list = series.tolist()
    # if there are nans in a row, fill whole row with nans
    #-------------------------------------------------------
    if not remove_nans:
        max_size = max([len(row) for row in series_list if type(row)==np.ndarray])
        for row_i, row in enumerate(series_list):
            if row is None or np.isnan(row).any():
                new_row = np.full(max_size, np.nan)
                series_list[row_i]=new_row

    if len(series_list)>0 and type(series_list[0])==list:
    #     print 'list to array?'
        series_list = [item[list_index] for item in series_list if type(item)==list and len(item)>list_index]
    #     # print series_list
    #     print np.array(np.array(series_list))

    array = np.array(series_list).squeeze()

    return array

def _process_new_data_df(group_df, preprocessed_directory, search_string='data', functions=[], kwlist=[], rerun=[], keep=[], file_limit=[], ):
    ''' process new data and add to group

    ==Args==
    -group_df : pandas dataframe containing group data
    -directory : relative directory containing unprocessed data
    -search_string :  string to identify unprocessed data, typically ".mat"
    -variables : variables to be derived from unprocessed data (e.g. slopes or slopes_norm)

    ==Out==
    -group_df : pandas dataframe containing group data

    ==Update==
    
    ==Comments==
    '''
    # iterate over columns in group_df
    # if column in keep, use old column values and drop suffix
    # otherwise keep new columnn values and drop suffix


    # merge all function dfs on filename and path
    # update appropriate rows and columns of group_df with function dfs, use merge? find indices and use update/combine?

    # get list of all preprocessed data files
    data_files = glob.glob(preprocessed_directory+'*'+search_string+'*')
    # remove directory and file extension
    data_filenames = [file.split('\\')[-1].split('.')[0] for file in data_files]
    trial_ids = [file.split('_')[-1] for file in data_filenames]

    # get list of processed data files (already added to group dataframe)
    if group_df.empty:
        processed_trial_ids = []
    else:
        processed_trial_ids = group_df.trial_id

    # get list of new data files (haven't been added to dataframe)
    new_trial_ids = list(set(trial_ids)-set(processed_trial_ids))
    new_data_files = [data_files[trial_id_i] for trial_id_i, trial_id in enumerate(trial_ids) if trial_id in new_trial_ids]

    # if any functions specified in rerun, set iterfiles to all data files
    if rerun:
        iterfiles = data_files
    # if there are no functions to rerun, only iterate through new files
    else:
        iterfiles = new_data_files

    print 'total data files:', len(data_files) 
    print 'new data files:', len(new_data_files)
    if rerun:
        print 'functions to rerun over all slices:', rerun

    # iterate over new files and update group_data structure
    #`````````````````````````````````````````````````````````````````
    print 'updating group data structure'

    # temporary df for storing rerun function outputs
    df_update = pd.DataFrame()

    # suffixes for merging df_update with group_df
    old_suffix = '_old'
    rerun_suffix = '_rerun'

    # iterate over data files
    for file_i, file in enumerate(iterfiles):
        print file_i

        # apply file limit
        #------------------
        if file_limit and file_i> file_limit:
            continue
        else:

            # get trial_id
            #----------------
            trial_id = file.split('\\')[-1].split('.')[0].split('_')[-1]

            # load data file
            #---------------
            with open(file, 'rb') as pkl_file:

                pre = pickle.load(pkl_file)

            # for each function create df and add to list (either new or rerun)
            new_dfs = []
            rerun_dfs = []
            # iterate over functions
            #-------------------------
            for func_i, func in enumerate(functions):

                # if new data file add to new list
                #---------------------------------
                if file in new_data_files:

                    print 'new file found:', trial_id
                    print kwlist[func_i]
                    df_temp = func(pre, df=[], keep=[], **kwlist[func_i])
                    # # print 'test'
                    new_dfs.append(df_temp)

                # if rerunning current func, add to rerun list
                #-------------------------------------------
                elif func in rerun:
                    # get corresponding row from group_df
                    df_locs = group_df[group_df.trial_id==trial_id].index.values
                    # get df from function
                    df_temp = func(pre, keep=keep, df=group_df.loc[df_locs], **kwlist[func_i])
                    # add to rerun list
                    rerun_dfs.append(df_temp)

            # if there are new df's, append to bottom of df_update
            #--------------------------------------------------
            if new_dfs:
                # combine all new_df's in list
                new_df = reduce(
                    lambda left,right: 
                    pd.merge(
                        left, 
                        right[np.append(right.columns.difference(left.columns).values, ['trial_id', 'field', 'location'])], 
                        on=['trial_id', 'field', 'location'], 
                        how='inner', ), 
                    new_dfs
                    )
                # add to update df (will be merged with old df later)
                df_update = df_update.append(new_df, ignore_index=True)

            # if there are rerun df's, append to bottom of df_update
            #--------------------------------------------------
            if rerun_dfs:
                # merge df's across the different functions
                rerun_df = reduce(
                    lambda left,right: 
                    pd.merge(
                        left, 
                        right[np.append(right.columns.difference(left.columns).values, ['trial_id', 'field', 'location'])], 
                        on=['trial_id', 'field','location'], 
                        how='inner', ), 
                    rerun_dfs
                    ) # validate='one_to_one'

                # add to update df (will be merged with old df later)
                df_update = df_update.append(rerun_df, ignore_index=True)
    
    # merge old and new df's
    #------------------------------------------------------------            
    # print 'merging old and new data'
    # if no old df is specified, set to group_df to update_df
    if group_df.empty:
        group_df=df_update
    # otherwise merge old df and update df
    elif not df_update.empty:
        group_df = group_df.merge(df_update, on=['trial_id', 'field', 'location'], how='outer', suffixes=[old_suffix, rerun_suffix])
    
    # update columns after merge
    #------------------------------
    # print 'dropping columns after merge'
    # get column values
    columns = group_df.columns.values
    # drop columns, based on keep argument
    drop_columns = []
    # iterate over columns
    for column in columns:
        # if old column
        if old_suffix in column:
            # get original column name
            column_original = column.split(old_suffix)[0]
            # get equivalent rerun column name
            column_new = column_original+rerun_suffix
            # if not a keep column
            if not any([temp for temp in keep if temp in column]):
                # use new to update old
                # FIXME, use overwrite kwarg for update
                group_df[column].update(group_df[column_new])
                # then drop new
                drop_columns.append(column_new)
            # if a keep column
            else:
                # use old to update new, drop old
                group_df[column_new].update(group_df[column])
                drop_columns.append(column)
    # drop specified columns
    group_df.drop(drop_columns, axis=1, inplace=True)

    # remove suffixes from merge
    #---------------------------
    # print 'renaming columns after merge'
    # dictionary for updated column names
    rename_dict={}
    # iterate over column names
    for column in group_df.columns.values:
        # check for suffixes and drop 
        if old_suffix in column:
            rename_dict[column] = column.split(old_suffix)[0]
        elif rerun_suffix in column:
            rename_dict[column] = column.split(rerun_suffix)[0]
    # rename columns
    group_df.rename(columns=rename_dict, inplace=True)

    return group_df

def _total_size(o, handlers={}, verbose=False):
    """ Returns the approximate memory footprint an object and all of its contents.

    Automatically finds the contents of the following builtin containers and
    their subclasses:  tuple, list, deque, dict, set and frozenset.
    To search other containers, add handlers to iterate over their contents:

        handlers = {SomeContainerClass: iter,
                    OtherContainerClass: OtherContainerClass.get_elements}

    """
    dict_handler = lambda d: itertools.chain.from_iterable(d.items())
    all_handlers = {tuple: iter,
                    list: iter,
                    deque: iter,
                    dict: dict_handler,
                    set: iter,
                    frozenset: iter,
                   }
    all_handlers.update(handlers)     # user handlers take precedence
    seen = set()                      # track which object id's have already been seen
    default_size = getsizeof(0)       # estimate sizeof object without __sizeof__

    def sizeof(o):
        if id(o) in seen:       # do not double count the same object
            return 0
        seen.add(id(o))
        s = getsizeof(o, default_size)

        # if verbose:
        #     print(s, type(o), repr(o), file=stderr)

        for typ, handler in all_handlers.items():
            if isinstance(o, typ):
                s += sum(map(sizeof, handler(o)))
                break
        return s

    return sizeof(o)
#############################################################################
# figure formatting
#############################################################################
class FormatFig:
    '''
    '''
    def __init__(self, ):
        '''
        '''
        pass

    def _standard_figformat(self, figures, axes, figdf, tight=True, xticks=True, yticks=True, xscale=True, yscale=True, xticks_minor=True, **kwargs):
        '''
        '''
        # set fontsizes
        #---------------
        # figures, axes = FormatFig()._set_fontsizes(figures=figures, axes=axes, figdf=figdf)
        # set yticks
        #--------------
        if yticks:
            figures, axes = FormatFig()._set_yticks(figures=figures, axes=axes, figdf=figdf)
        # set xticks
        #--------------
        if xticks:
            figures, axes = FormatFig()._set_xticks(figures=figures, axes=axes, figdf=figdf)

        # set xticks
        #--------------
        if xticks_minor:
            figures, axes = FormatFig()._set_xticks_minor(figures=figures, axes=axes, figdf=figdf)
        # scale yticklabels
        #-------------------
        if yscale:
            figures, axes = FormatFig()._scale_yticklabels(figures=figures, axes=axes, figdf=figdf)
        # scale xticklabels
        #-------------------
        if xscale:
            figures, axes = FormatFig()._scale_xticklabels(figures=figures, axes=axes, figdf=figdf)
        # set yticklabel decimals
        #-----------------------
        figures, axes = FormatFig()._set_yticklabel_decimals(figures=figures, axes=axes, figdf=figdf)
        # set xticklabel decimals
        #-----------------------
        figures, axes = FormatFig()._set_xticklabel_decimals(figures=figures, axes=axes, figdf=figdf)
        # set fontsizes
        #---------------
        figures, axes = FormatFig()._set_fontsizes(figures=figures, axes=axes, figdf=figdf)
        # set tight layout
        #-----------------
        if tight:
            for figkey in figures:
                figures[figkey].tight_layout()
            # plt.tight_layout()

        return figures, axes

    def _set_xticks(self, figures, axes, figdf, tick_max=10,**kwargs):
        '''
        '''

        print 'setting ticks'
        # get xlim across all figures
        #-------------------------------------
        xlim={}
        xlims=[]
        # iterate over figures
        for axkey, ax in axes.iteritems():
            
            # get ticks
            xticks = ax.get_xticks()
            print 'xTICKS',xticks
            # get updated lims (1  significant figure)
            xlim[axkey] = _to_1sigfig(ax.get_xlim())
            xlims.append(copy.copy(xlim[axkey]))
        # find x lims across all figures in group
        xlim_all = [min([temp[0] for temp in xlims]), max([temp[1] for temp in xlims])]

        # iterate over figures
        #--------------------------------------
        for axkey, ax in axes.iteritems():
            # # check if axes are numeric
            # #-----------------------------------------------------------
            # # get current ticklabels
            # xticklabels = [tick.get_text() for tick in ax.get_xticklabels()]
            # while '' in xticklabels:
            #     xticklabels.remove('')
            # print xticklabels
            # # check if tick labels are numeric
            # for label in xticklabels:
            #     try:
            #         label_float = type(float(label))==float
            #     except ValueError:
            #         label_float=False
            # # if ticklabels are numeric, scale them, otherwise leave them alone
            # if not label_float:
            #     continue
            # get ticks and limits
            #----------------------
            xticks = ax.get_xticks()
            # get scale factors
            #--------------------
            if 'fig_xscale' in figdf:
                xscale = figdf.loc[axkey].fig_xscale.values[0]
            else:
                xscale=1
            # get x lims
            xlim_current = list(ax.get_xlim())
            print 'xlim_current', xlim_current

            # set xlim across all figures
            #----------------------------
            if all(figdf.fig_xlim_all):
                xlim_current = xlim_all
            if 'fig_xmin' in figdf.keys():
                xlim_current[0] = figdf.loc[axkey].fig_xmin.unique()[0]/xscale
            if 'fig_xmax' in figdf.keys():
                xlim_current[1] = figdf.loc[axkey].fig_xmax.unique()[0]/xscale
            ax.set_xlim(xlim_current)

            # set x ticks
            #------------------------------------------------------------
            print 'setting ticks'
            # get current lims and ticks
            xlim[axkey] = copy.copy(list(ax.get_xlim()))
            print 'xlim', xlim[axkey]
            xticks = np.array(ax.get_xticks())
            xlim_scaled = np.array(xlim[axkey])

            # get number of ticks and tick spacing
            #------------------------------------------------
            # default nxticks and dxticks
            nxticks=5.
            dxticks = float(abs(xlim_scaled[1]-xlim_scaled[0]))/nxticks
            # specifx dxticks
            if 'fig_dxticks' in figdf:
                dxticks = figdf.loc[axkey].fig_dxticks.unique()[0]/xscale
                nxticks = len(np.arange(xlim_scaled[0], xlim_scaled[1], dxticks))
            # specify nxticks
            elif 'fig_nxticks' in figdf:
                nxticks = figdf.loc[axkey].fig_nxticks.unique()[0]
                dxticks = float(abs(xlim_scaled[1]-xlim_scaled[0]))/nxticks
            # make sure that there isn't a crazy number of ticks
            #-----------------------------------------------------
            if nxticks>tick_max:
                nxticks=tick_max
                dxticks = float(abs(xlim_scaled[1]-xlim_scaled[0]))/nxticks

            # assert a specific tick value
            #------------------------------
            if 'fig_xtick_assert' in figdf:
                xtick_assert = figdf.loc[axkey].fig_xtick_assert.unique()[0]/xscale
            else:
                xtick_assert = xlim_scaled[0]
            # crete new ticks
            #--------------------
            xticks_new_1 = np.flip(np.arange(xtick_assert, xlim_scaled[0], -dxticks))
            # print xticks_new_1
            xticks_new_2 = np.arange(xtick_assert, xlim_scaled[1], dxticks)
            xticks_new = np.append(xticks_new_1, xticks_new_2)
            print 'xticks new',xticks_new_1, xticks_new_2, xticks_new
            print 'xticks new',xtick_assert, xlim_scaled, dxticks, xscale

            # roud tick decimals
            #-------------------------
            # if 'fig_xtick_round' in figdf:
            #     print 'rounding xticks'
            #     decimals = figdf.loc[axkey].fig_xtick_round.unique()[0]
            #     xticks_new = np.round(xticks_new, decimals=decimals)
            # if 'fig_xticks' in figdf:
            #     print 'setting xticks'
            #     xticks_new = figdf.fig_xticks.values[0]
            print 'xticks_new',xticks_new
            print 'xticks_new',xticks_new/xscale
            ax.set_xticks(xticks_new)
            ax.set_xticklabels(xticks_new)

        return figures, axes
    
    def _set_xticks_minor(self, figures, axes, figdf, **kwargs):
        '''
        '''
        print 'setting xticks minor'
        # iterate over figures
        figdf = figdf.reset_index().set_index(['figure', 'subgroup', 'trace'])
        for axkey, ax in axes.iteritems():
            subgroups = figdf.loc[axkey].index.get_level_values('subgroup').unique().values
            xticks_minor=[]
            xticklabels_minor=[]
            if 'xticks_minor_loc' in figdf:
                for subgroup in subgroups:
                    subgroup_loc = figdf.loc[(axkey, subgroup)]['xticks_minor_loc'].unique()[0]
                    xticks_minor.append(subgroup_loc)
                    xticklabels_minor.append(subgroup)

            ax.set_xticks(xticks_minor, minor=True)
            ax.set_xticklabels(xticklabels_minor, minor=True, )

        return figures, axes

    def _set_yticks(self, figures, axes, figdf, tick_max=10,**kwargs):
        '''
        '''
        print 'setting ticks'
        # get ylim across all figures
        #-------------------------------------
        ylim={}
        ylims=[]
        # iterate over figures
        for axkey, ax in axes.iteritems():
            # get ticks
            yticks = ax.get_yticks()
            print 'YTICKS',yticks
            # get updated lims (1  significant figure)
            ylim[axkey] = _to_1sigfig(ax.get_ylim())
            # ylim[axkey] = ax.get_ylim()
            ylims.append(copy.copy(ylim[axkey]))
        # find x and y lims across all figures in group
        ylim_all = [min([temp[0] for temp in ylims]), max([temp[1] for temp in ylims])]

        # iterate over figures
        #--------------------------------------
        for axkey, ax in axes.iteritems():
            # draw figure to create ticklabels
            #----------------------------------
            # figures[axkey].canvas.draw()
            # # check that ticklabels are numeric
            # #----------------------------------------
            # # get current ticklabels
            # yticklabels = [tick.get_text() for tick in ax.get_yticklabels()]
            # while '' in yticklabels:
            #     yticklabels.remove('')
            # print yticklabels
            # # check if tick labels are numeric
            # for label in yticklabels:
            #     try:
            #         label_float = type(float(label))==float
            #     except ValueError:
            #         label_float=False
            # # if ticklabels are numeric, scale them, otherwise leave them alone
            # if not label_float:
            #     print 'yticks not numeric'
            #     continue

            # get ticks and limits
            #----------------------
            yticks = ax.get_yticks()
            # get scale factors
            #--------------------
            if 'fig_yscale' in figdf:
                yscale = figdf.loc[axkey].fig_yscale.values[0]
            else:
                yscale=1
            # get y lims
            ylim_current = list(ax.get_ylim())
            print 'ylim_current', ylim_current

            # set ylim across all figures
            #----------------------------
            if all(figdf.fig_ylim_all):
                ylim_current = ylim_all
            if 'fig_ymin' in figdf.keys():
                ylim_current[0] = figdf.loc[axkey].fig_ymin.unique()[0]/yscale
            if 'fig_ymax' in figdf.keys():
                ylim_current[1] = figdf.loc[axkey].fig_ymax.unique()[0]/yscale
            ax.set_ylim(ylim_current)

            # set y ticks
            #------------------------------------------------------------
            print 'setting ticks'
            # get current lims and ticks
            ylim[axkey] = copy.copy(list(ax.get_ylim()))
            print 'ylim', ylim[axkey]
            yticks = yscale*np.array(ax.get_yticks())
            ylim_scaled = yscale*np.array(ylim[axkey])

            # get number of ticks and tick spacing
            #------------------------------------------------
            # default nyticks and dyticks
            nyticks=5.
            dyticks = float(abs(ylim_scaled[1]-ylim_scaled[0]))/nyticks
            # specify dyticks
            if 'fig_dyticks' in figdf:
                dyticks = figdf.loc[axkey].fig_dyticks.unique()[0]
                nyticks = len(np.arange(ylim_scaled[0], ylim_scaled[1], dyticks))
            # specify nyticks
            elif 'fig_nyticks' in figdf:
                nyticks = figdf.loc[axkey].fig_nyticks.unique()[0]
                dyticks = float(abs(ylim_scaled[1]-ylim_scaled[0]))/nyticks
            # make sure that there isn't a crazy number of ticks
            #-----------------------------------------------------
            if nyticks>tick_max:
                nyticks=tick_max
                dyticks = float(abs(ylim_scaled[1]-ylim_scaled[0]))/nyticks

            # assert a specific tick value
            #------------------------------
            if 'fig_ytick_assert' in figdf:
                ytick_assert = figdf.loc[axkey].fig_ytick_assert.unique()[0]
            else:
                ytick_assert = ylim_scaled[0]
            # crete new ticks
            #--------------------
            yticks_new_1 = np.flip(np.arange(ytick_assert, ylim_scaled[0], -dyticks))
            # print yticks_new_1
            yticks_new_2 = np.arange(ytick_assert, ylim_scaled[1], dyticks)
            yticks_new = np.append(yticks_new_1, yticks_new_2)
            print 'yticks new',yticks_new_1, yticks_new_2
            print 'yticks new',ytick_assert, ylim_scaled, dyticks, yscale

            # roud tick decimals
            #-------------------------
            if 'fig_ytick_round' in figdf:
                decimals = figdf.loc[axkey].fig_ytick_round.unique()[0]
                yticks_new = np.round(yticks_new, decimals=decimals)
            if 'fig_yticks' in figdf:
                yticks_new = figdf.fig_yticks.values[0]
            print 'yticks_new',yticks_new
            ax.set_yticks(yticks_new/yscale)
            ax.set_yticklabels(yticks_new/yscale)

        return figures, axes

    def _set_yticklabel_decimals(self, figures, axes, figdf, **kwargs):
        '''
        '''
        print 'setting tick decimals'
        for axkey, ax in axes.iteritems():
            # get ticks
            #----------------------
            yticks = ax.get_yticks()
            # apply scale factor
            #-------------------------
            if 'fig_yscale' in figdf:
                print 'yscale factor found'
                yscale = figdf.loc[axkey].fig_yscale.values[0]
            else:
                yscale=1
            # apply scaling factor
            yticks = yscale*yticks
            # set tick decimal places
            #---------------------------------------------------------------
            if 'fig_ytick_decimals' in figdf:
                # get decimal places
                dec = figdf.loc[axkey].fig_ytick_decimals.unique()[0]
                # get ticklabels
                yticklabels = ['%.{}f'.format(dec) % float(tick) for tick in yticks]
                print 'yticklabels',yticklabels
                # set ticklabels
                ax.set_yticklabels(yticklabels)
        return figures, axes
    
    def _set_xticklabel_decimals(self, figures, axes, figdf, **kwargs):
        '''
        '''
        print 'setting tick decimals'
        for axkey, ax in axes.iteritems():
            # get ticks
            #-----------------------
            xticks = ax.get_xticks()
            # scale factor
            #-------------------------
            if 'fig_xscale' in figdf:
                print 'xscale factor found'
                xscale = figdf.loc[axkey].fig_xscale.values[0]
            else:
                xscale=1
            # apply scaling factor
            xticks = xscale*xticks
            # set tick decimal places
            #---------------------------------------------------------------
            if 'fig_xtick_decimals' in figdf:
                # get decimal places
                dec = figdf.loc[axkey].fig_xtick_decimals.unique()[0]
                # get ticklabels
                xticklabels = ['%.{}f'.format(dec) % float(tick) for tick in xticks]
                print 'xticklabels',xticklabels
                # set ticklabels
                ax.set_xticklabels(xticklabels)
        return figures, axes

    def _scale_yticklabels(self, figures, axes, figdf, **kwargs):
        '''
        '''
        print 'scaling yticks'
        for axkey, ax in axes.iteritems():
            # get ticks
            yticks = ax.get_yticks()
            #--------------------------------------
            # yticks
            #-------------------------------------
            if 'fig_yscale' in figdf:
                print 'yscale factor found'
                yscale = figdf.loc[axkey].fig_yscale.values[0]
            else:
                yscale=1
            # apply scaling factor
            yticks = yscale*yticks
            # # get current ticklabels
            # yticklabels = [tick.get_text() for tick in ax.get_yticklabels()]
            # print yticklabels
            # # check if tick labels are numeric
            # for label in yticklabels:
            #     try:
            #         label_float = type(float(label))==float
            #     except ValueError:
            #         label_float=False
            # if ticklabels are numeric, scale them, otherwise leave them alone
            # if label_float:
            # update tick decimal places
            #-----------------------------
            # if 'fig_ytick_decimals' in figdf:
            #     # get decimal places
            #     dec = figdf.loc[axkey].fig_ytick_decimals.unique()[0]
            #     # get ticklabels with adjusted decimals
            #     yticklabels = ['%.{}f'.format(dec) % tick for tick in yticks]
            # else:
            #     yticklabels=yticks
            # set ticklabels
            yticklabels=yticks
            ax.set_yticklabels(yticklabels)
        return figures, axes

    def _scale_xticklabels(self, figures, axes, figdf, **kwargs):
        '''
        '''
        print 'scaling xticks'
        for axkey, ax in axes.iteritems():
            # get ticks
            xticks = ax.get_xticks()
            print xticks
            # xticklabels = list(ax.get_xticklabels())
            #--------------------------------------
            # yticks
            #-------------------------------------
            if 'fig_xscale' in figdf:
                print 'xscale factor found'
                xscale = figdf.loc[axkey].fig_xscale.values[0]
            else:
                xscale=1
            # apply scaling factor
            xticks = xscale*xticks
            print xticks
            # # get current ticklabels
            # yticklabels = [tick.get_text() for tick in ax.get_yticklabels()]
            # print yticklabels
            # # check if tick labels are numeric
            # for label in yticklabels:
            #     try:
            #         label_float = type(float(label))==float
            #     except ValueError:
            #         label_float=False
            # if ticklabels are numeric, scale them, otherwise leave them alone
            # if label_float:
            # update tick decimal places
            #-----------------------------
            # if 'fig_xtick_decimals' in figdf:
            #     # get decimal places
            #     dec = figdf.loc[axkey].fig_xtick_decimals.unique()[0]
            #     # get ticklabels with adjusted decimals
            #     xticklabels = ['%.{}f'.format(dec) % tick for tick in xticks]
            # else:
            #     xticklabels=xticks
            # set ticklabels
            xticklabels=xticks
            ax.set_xticklabels(xticklabels)
        return figures, axes

    def _set_fontsizes(self, figures, axes, figdf, **kwargs):
        '''
        '''
        print 'setting axes properties'
        # iterate through figures
        for axkey, ax in axes.iteritems():
            # set ticklabel fontweight and size
            #-----------------------------------
            print 'setting ticklabels'
            # ytick labels
            for label in ax.get_xticklabels():
                # fontweight
                label.set_fontweight(figdf.loc[axkey].fig_xtick_fontweight.values[0])
                # fonstize
                label.set_fontsize(figdf.loc[axkey].fig_xtick_fontsize.values[0])
            # ytick labels
            for label in ax.get_yticklabels():
                # fontweight
                label.set_fontweight(figdf.loc[axkey].fig_ytick_fontweight.values[0])
                # fontsize
                label.set_fontsize(figdf.loc[axkey].fig_ytick_fontsize.values[0])

            for label in ax.get_xticklabels(which='minor'):
                print 'minor label', label
                label.set_fontweight(figdf.loc[axkey].fig_xtick_fontweight.values[0])
                label.set_fontsize(figdf.loc[axkey].fig_xtick_fontsize.values[0])

            # turn off axes box
            #-----------------------------------
            # turn off top and right axes lines
            if all(figdf.loc[axkey].fig_boxoff):
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)
            # turn off all axes lines
            if 'fig_axesoff' in figdf.keys() and all(figdf.loc[axkey].fig_axesoff):
                ax.set_axis_off()

            # set axes linewidth and tick position
            #--------------------------------------
            ax.spines['left'].set_linewidth(figdf.loc[axkey].fig_axes_linewidth.unique()[0])
            ax.spines['bottom'].set_linewidth(figdf.loc[axkey].fig_axes_linewidth.unique()[0])
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')

            # set axes labels
            #----------------
            if 'fig_xlabel' in figdf:
                ax.set_xlabel(figdf.loc[axkey].fig_xlabel.values[0], fontsize=figdf.loc[axkey]['fig_xlabel_fontsize'].values[0], fontweight=figdf.loc[axkey].fig_xlabel_fontweight.values[0])
            if 'fig_ylabel' in figdf:
                ax.set_ylabel(figdf.loc[axkey].fig_ylabel.values[0], fontsize=figdf.loc[axkey]['fig_ylabel_fontsize'].values[0], fontweight=figdf.loc[axkey].fig_ylabel_fontweight.values[0])

            # set tight layout
            #-------------------
            if all(figdf.loc[axkey].fig_tight_layout):
                plt.figure(figures[axkey].number)
                plt.tight_layout()

        return figures, axes

#############################################################################
# plotting functions
#############################################################################
class PlotFuncs:
    ''' each plot function should have coresponding function to generate figdf in BuildFigDF
    '''
    def __init__(self):
        '''
        '''
        pass

    def _figdf_generator(self, df, figdf, variables, remove_nans=True,):
        '''
        '''
        # set figdf to hierarchical index (figure, subgroup, trace)
        figdf = figdf.reset_index().set_index(['figure','subgroup','trace'])
        # get list of all figure names
        figures = figdf.index.get_level_values('figure').unique().values
        # iterate over figures
        for figkey in figures:
            # get subgroups list
            subgroups = figdf.loc[figkey].index.get_level_values('subgroup').unique().values
            # iterate over subgroups of traces
            for subkey in subgroups:
                # get list of tracekeys
                traces = figdf.loc[(figkey, subkey)].index.get_level_values('trace').unique().values
                # iterate over traces
                for tracekey in traces:
                    # get params from figdf
                    params = figdf.loc[(figkey,subkey, tracekey)]
                    # get corresponding series and convert to array
                    #----------------------------------------------
                    # if available, get series, otherwise next iter
                    try:
                        trace_df=df.loc[tracekey,slice(None)]
                    except:
                        print tracekey, 'not in df index'
                        print df.index.names
                        continue
                    try:
                        trace_df=trace_df[variables]
                    except:
                        print variables, 'not in df'
                        continue

                    if remove_nans:
                        # remove nans from trace series
                        na = trace_df.isna()
                        trace_df = trace_df[~na]

                    yield figkey, subkey, tracekey, params, trace_df 

    def _select_kwargs(self, funcs, kwargs):
        ''' remove kwargs from dictionary if the function cannot handle them

            Arguments
            -----------
            ~func: function to check arguments
            ~kwargs: dictionary of keyword arguments with key corresponding to arguments in func

            Returns
            -----------
            ~kwargs_new: new dictionary of keyword arguments with all keys as valid arguments for func

        '''
        kwargs_new={}
        for key, val in kwargs.iteritems():
            for func in funcs:
                if inspect.isclass(func):
                    args, varargs, varkw, defaults=inspect.getargspec(func.__init__)
                elif inspect.isfunction(func):
                    args, varargs, varkw, defaults=inspect.getargspec(func)
                print 'args',args
                print varargs
                print varkw
                print defaults
                if key in args and key not in kwargs_new:
                    kwargs_new[key]=val

        return kwargs_new
    
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
    
    def _draw_lines(self, fig, ax, figdf):
        '''
        '''
        # set figdf to hierarchical index (figure, subgroup, trace)
        figdf = figdf.reset_index().set_index(['figure','subgroup','trace'])
        # get list of all figure names
        figures = figdf.index.get_level_values('figure').unique().values
        # iterate over figures
        for figkey in figures:
            if 'fig_hlines' in figdf:
                yvals = figdf.loc[figkey].fig_hlines.unique()[0]
                xlims = ax[figkey].get_xlim()
                ax[figkey].hlines(yvals, xlims[0], xlims[1])
            if 'fig_vlines' in figdf:
                xvals = figdf.loc[figkey].fig_vlines.unique()[0]
                ylims = ax[figkey].get_ylim()
                ax[figkey].vlines(xvals, ylims[0], ylims[1])
        return fig, ax

    def _draw_rectangle(self, fig, ax, figdf):
        '''
        '''
        print 'plotting rectangle'
        # set figdf to hierarchical index (figure, subgroup, trace)
        figdf = figdf.reset_index().set_index(['figure','subgroup','trace'])
        # get list of all figure names
        figures = figdf.index.get_level_values('figure').unique().values
        # iterate over figures
        for figkey in figures:
            rec = Rectangle((.3, 0), .1, 2, alpha=0.5, facecolor='gray', edgecolor='gray')
            ax[figkey].add_patch(rec)
            plt.draw()
        return fig, ax

    def _dose_response(self, df, figdf, variable, **kwargs):
        '''
        '''
        print 'plotting:', inspect.stack()[0][3]
        # preallocate fig and ax objects
        fig={}
        ax={}
        locations={}
        heights={}
        y_errors={}
        # figdf = figdf.reset_index()
        # set figdf to hierarchical index (figure, subgroup, trace)
        figdf = figdf.reset_index().set_index(['figure','subgroup','trace'])
        print figdf.index.names
        # get list of all figure names
        figures = figdf.index.get_level_values('figure').unique().values
        df = df.set_index(['field', 'path_1_syn_num'])
        # iterate over figures
        for figkey in figures:
            # create figure, passing params as **kwargs
            fig[figkey], ax[figkey] = plt.subplots()
            # get subgroups list
            subgroups = figdf.loc[figkey].index.get_level_values('subgroup').unique().values
            locations[figkey] = []
            heights[figkey]=[]
            y_errors[figkey]=[]
            colors=[]
            ecolors=[]
            xticks=[]
            xticklabels=[]
            # iterate over subgroups of traces
            for subkey in subgroups:
                traces = figdf.loc[(figkey, subkey)].index.get_level_values('trace').unique().values
                # iterate over traces
                for tracekey in traces:
                    # get params from figdf
                    params = figdf.loc[(figkey,subkey, tracekey)]

                    try:
                        # get corresponding series
                        trace_series = df.loc[tracekey, slice(None)][variable]
                    except:
                        print tracekey, 'not in df index'
                        continue

                    # convert to array
                    data_array = _2array(trace_series, remove_nans=True, remove_nans_axis=1)
                    # check that data is correct  shape
                    if type(data_array)==np.ndarray and len(data_array.shape)==1:
                        # print data_array
                        # mean across cells
                        data_mean = np.mean(data_array, axis=0)
                        #std across cells
                        data_std = np.std(data_array, axis=0)
                        # sem across cells
                        data_sem = stats.sem(data_array, axis=0)
                        # check if trace location has been set manually
                        # if 'trace_location' in params:
                        plot_location = params.trace_location
                        locations[figkey].append(plot_location)
                        xticks.append(plot_location)
                        if 'trace_label' in params:
                            xticklabels.append(params.trace_label)
                        else:
                            xticklabels.append(plot_location)

                        colors.append(params.trace_color)
                        ecolors.append(params.trace_ecolor)
                        heights[figkey].append(data_mean)
                        y_errors[figkey].append(data_sem)
                        plt.errorbar(locations[figkey][-1], data_mean, data_sem, ecolor=params.trace_ecolor, elinewidth=2)
                        if 'fig_data_style' in figdf and figdf.fig_data_style.unique()[0]=='point':
                            ax[figkey].plot(locations[figkey][-1], heights[figkey][-1], color=colors[-1], linestyle='None', markersize=30, marker='.')
            if 'fig_barwidth' in figdf:
                width=figdf.fig_barwidth.unique()[0]
            else:
                width=1

            # if 'fig_data_style' in figdf and figdf.fig_data_style.unique()[0]=='point':
            #     ax[figkey].plot(locations[figkey], heights[figkey], linestyle='None', markersize=30, marker='.')

            # barcontainer = ax[figkey].bar(locations[figkey], heights[figkey], width=width, tick_label=xticklabels)
            # ax[figkey].set_xscale('symlog')
            # print xticks
            # ax[figkey].set_xticks(xticks)
            # ax[figkey].set_xticklabels(xticks)


            # set bar color
            # for bar_i, bar in enumerate(barcontainer):
            #     bar.set_color(colors[bar_i])

        # fig, ax = _format_figures(fig=fig, ax=ax, figdf=figdf)

        # fig, ax = FormatFig()._set_fontsizes(figures=fig, axes=ax, figdf=figdf)
        # fig, ax = FormatFig()._set_yticks(figures=fig, axes=ax, figdf=figdf)
        # fig, ax = FormatFig()._scale_ticklabels(figures=fig, axes=ax, figdf=figdf)
        # fig, ax = FormatFig()._set_ticklabel_decimals(figures=fig, axes=ax, figdf=figdf)

        fig, ax = FormatFig()._standard_figformat(figures=fig, axes=ax, figdf=figdf)

        plt.show(block=False)

        # fig, ax = figsetup._format_figures(fig=fig, ax=ax, figdf=figdf, xlim=xlim, ylim=ylim)

        # for figkey in ax:
        # if 'fig_label_pvalues' in figdf and all(figdf.loc[figkey].fig_label_pvalues):
        #     _label_pvalues(ax=ax[figkey], ttests=ttests, traces=tracelist[figkey], x_locations=locations[figkey], y_means=heights[figkey], y_errors=y_errors[figkey])
        # # get x and y lims
        # xlim[figkey] = ax[figkey].get_xlim()
        # ylim[figkey] = ax[figkey].get_ylim()

        # # ylim may change after add significance markers, so reformat
        # fig, ax = figsetup._format_figures(fig=fig, ax=ax, figdf=figdf, xlim=xlim, ylim=ylim)

        return fig, ax

    # def _twinx(self, dfs, figdfs, variables, funcs, func_kws, **kwargs):
    #     '''
    #     '''
    def _var2var_corr(self, df, figdf, variables, df_funcs=[], df_func_kws=[],array_funcs_x=[],array_func_kws_x=[],array_funcs_y=[],array_func_kws_y=[], df_sorted=[], figformat='standard', **kwargs):
        '''FIXME add docs
        '''
        
        # report progress
        print 'plotting:', inspect.stack()[0][3]

        # preallocate figures and axes
        fig={}
        ax={}
        # iteratively apply functions for creating new column to dataframe
        #----------------------------------------------------------------
        if len(df_funcs)>0:
            for df_func_i, df_func in enuemrate(df_funcs):
                df = df_func(df, axis=1, **df_func_kws[df_func_i])
        # iterate over traces to be plotted
        #-----------------------------------
        for figkey, subkey, tracekey, params, trace_df in self._figdf_generator(df, figdf, variables):
            # make new figure if doesnt exist yet
            #-------------------------------------
            if figkey not in ax:
                fig[figkey], ax[figkey] = plt.subplots()
            # get series from df
            #----------------------------------------
            trace_series_x = trace_df[variables[0]]
            trace_series_y = trace_df[variables[1]]
            # check for missing values and remove rows
            #---------------------------------------------
            nax = ~trace_series_x.isna()
            nay = ~trace_series_y.isna()
            # drop any rows that contain nan
            #-----------------------------------------------
            trace_series_x = trace_series_x[nax&nay]
            trace_series_y = trace_series_y[nax&nay]
            # convert to array
            #--------------------------------------------------
            data_array_x = _2array(trace_series_x, remove_nans=True, remove_nans_axis=1).astype(float)
            data_array_y = _2array(trace_series_y, remove_nans=True, remove_nans_axis=1).astype(float)
            # apply array functions
            #----------------------------------------------------
            for i, array_func in enumerate(array_funcs_x):
                data_array_x = array_func(data_array_x, **array_func_kws_x[i])
            for i, array_func in enumerate(array_funcs_y):
                data_array_y = array_func(data_array_y, **array_func_kws_y[i])
            # check array shapes
            #------------------------------
            if type(data_array_x)==np.ndarray and len(data_array_x.shape)==1 and data_array_x.shape[0]>0 and type(data_array_y)==np.ndarray and len(data_array_y.shape)==1 and data_array_y.shape[0]>0:

                # data_x_subgroup[figkey][subkey] = np.append(data_x_subgroup[figkey][subkey], data_array_x)
                # data_y_subgroup[figkey][subkey] = np.append(data_y_subgroup[figkey][subkey], data_array_y)
                # colors[figkey][subkey].append([param.trace_color]*len(data_array_y))
                # markersizes[figkey][subkey].append([param.trace_markersize]*len(data_array_y))

                # if data_array_x.shape[0]>0 and  data_array_x.shape[0]>0:
                # plot
                #------------------------------------------------------
                print data_array_x[0]
                print data_array_y.shape
                ax[figkey].plot(data_array_x, data_array_y, marker='.',linestyle='None', color=params.trace_color, markersize=params.trace_markersize)
        # format figure
        #----------------
        if figformat=='standard':
            fig, ax = FormatFig()._standard_figformat(figures=fig, axes=ax, figdf=figdf)

        plt.show(block=False)

        return fig, ax

    def _var2var_corr_connect(self, df, figdf, x_variable, y_variable, connect_variable, df_funcs=[], df_func_kws=[],array_funcs_x=[],array_func_kws_x=[],array_funcs_y=[],array_func_kws_y=[], df_sorted=[], figformat='standard', **kwargs):
        '''FIXME add docs
        '''
        variables = [x_variable, y_variable, connect_variable]
        index_names = list(df.index.names)
        # get list of all tracekeys
        all_tracekeys = figdf.reset_index().set_index('trace').index.values
        print 'all_tracekeys:', all_tracekeys
        # report progress
        print 'plotting:', inspect.stack()[0][3]

        # preallocate figures and axes
        fig={}
        ax={}
        # iteratively apply functions for creating new column to dataframe
        #----------------------------------------------------------------
        if len(df_funcs)>0:
            for df_func_i, df_func in enuemrate(df_funcs):
                df = df_func(df, axis=1, **df_func_kws[df_func_i])
        # iterate over traces to be plotted
        #-----------------------------------
        for figkey, subkey, tracekey, params, trace_df in self._figdf_generator(df, figdf, variables):
            # make new figure if doesnt exist yet
            #-------------------------------------
            if figkey not in ax:
                fig[figkey], ax[figkey] = plt.subplots()

                # draw connections between any points that share the connect variable
                ###############################################################
                df_connect = df.reset_index().set_index(connect_variable)
                # get list of unique connect values
                #--------------------------------------
                connect_vals = df_connect.index.unique()
                # iterate over connect values
                #--------------------------------------
                for connect_val in connect_vals:
                    # get df for current connect value
                    df_current = df_connect.loc[connect_val]
                    # if there is more than one entry with the same connect variable
                    if len(df_current.shape)>1:
                        # reset index to match figdf tracekey
                        df_current = df_current.reset_index().set_index(index_names)
                        # get list of index values
                        index_values = df_current.index.values
                        # get x values
                        x_values = df_current[x_variable].values
                        # y values
                        y_values = df_current[y_variable].values
                        # pairs
                        xy_pairs = zip(x_values, y_values, index_values)
                        # get combinations of 
                        xy_combos = itertools.combinations(xy_pairs, 2)
                        # iterate over combos
                        for combo in xy_combos:
                            # print 'combo:', (combo[0][-1],)
                            # print (combo[0][-1],) in all_tracekeys
                            # if (combo[0][-1],) in all_tracekeys and (combo[1][-1],) in all_tracekeys:
                            x = [combo[0][0], combo[1][0]]
                            y = [combo[0][1], combo[1][1]]
                            ax[figkey].plot(x, y, 'grey')


            # get series from df
            #----------------------------------------
            trace_series_x = trace_df[variables[0]]
            trace_series_y = trace_df[variables[1]]
            # check for missing values and remove rows
            #---------------------------------------------
            nax = ~trace_series_x.isna()
            nay = ~trace_series_y.isna()
            # drop any rows that contain nan
            #-----------------------------------------------
            trace_series_x = trace_series_x[nax&nay]
            trace_series_y = trace_series_y[nax&nay]
            # convert to array
            #--------------------------------------------------
            data_array_x = _2array(trace_series_x, remove_nans=True, remove_nans_axis=1).astype(float)
            data_array_y = _2array(trace_series_y, remove_nans=True, remove_nans_axis=1).astype(float)
            # apply array functions
            #----------------------------------------------------
            for i, array_func in enumerate(array_funcs_x):
                data_array_x = array_func(data_array_x, **array_func_kws_x[i])
            for i, array_func in enumerate(array_funcs_y):
                data_array_y = array_func(data_array_y, **array_func_kws_y[i])
            # check array shapes
            #------------------------------
            if type(data_array_x)==np.ndarray and len(data_array_x.shape)==1 and data_array_x.shape[0]>0 and type(data_array_y)==np.ndarray and len(data_array_y.shape)==1 and data_array_y.shape[0]>0:

                # data_x_subgroup[figkey][subkey] = np.append(data_x_subgroup[figkey][subkey], data_array_x)
                # data_y_subgroup[figkey][subkey] = np.append(data_y_subgroup[figkey][subkey], data_array_y)
                # colors[figkey][subkey].append([param.trace_color]*len(data_array_y))
                # markersizes[figkey][subkey].append([param.trace_markersize]*len(data_array_y))

                # if data_array_x.shape[0]>0 and  data_array_x.shape[0]>0:
                # plot
                #------------------------------------------------------
                print data_array_x[0]
                print data_array_y.shape
                ax[figkey].plot(data_array_x, data_array_y, marker='.',linestyle='None', color=params.trace_color, markersize=params.trace_markersize)
        # format figure
        #----------------
        if figformat=='standard':
            fig, ax = FormatFig()._standard_figformat(figures=fig, axes=ax, figdf=figdf)

        plt.show(block=False)

        return fig, ax

    def _trace_mean(self, df, figdf, variables, xvals=None, xvals_kw={},array_funcs=[], array_func_kws=[], df_funcs=[], df_func_kws=[], df_sorted=[], figformat='standard', **kwargs):
        '''FIXME add docs
        '''
        
        # report progress
        print 'plotting:', inspect.stack()[0][3]
        # update dt
        #------------------
        if 'dt' in kwargs:
            dt=kwargs['dt']
        else:
            dt=1
        # preallocate figures and axes
        fig={}
        ax={}
        # iteratively apply functions for creating new column to dataframe
        #----------------------------------------------------------------
        if len(df_funcs)>0:
            for df_func_i, df_func in enuemrate(df_funcs):
                df = df_func(df, axis=1, **df_func_kws[df_func_i])
        # iterate over traces to be plotted
        #-----------------------------------
        for figkey, subkey, tracekey, params, trace_df in self._figdf_generator(df, figdf, variables):
            # make new figure if doesnt exist yet
            #-------------------------------------
            if figkey not in ax:
                fig[figkey], ax[figkey] = plt.subplots()
            # convert df to array
            #---------------------
            data_array = _2array(trace_df, remove_nans=True, remove_nans_axis=1)
            print data_array
            print 'data array shape:', data_array.shape

            # iteratively apply array functions
            #-----------------------------------
            if len(array_funcs)>0:
                for array_func_i, array_func in enumerate(array_funcs):
                    data_array=array_func(data_array, **array_func_kws[array_func_i])

            # make sure array is the correct shape
            if type(data_array)==np.ndarray and data_array.shape[0]!=0:
                # if array is 1d, convert to 2d
                if len(data_array.shape)==1:
                    data_array = data_array.reshape((1,-1))
                # mean across slices
                data_mean = np.mean(data_array, axis=0)
                #std across slices
                data_std = np.std(data_array, axis=0)
                # sem across slices
                data_sem = stats.sem(data_array, axis=0)
                # set x values
                #-------------------------------------
                if xvals is None:
                    # default time vector
                    xvals_array = np.linspace(0, len(data_mean)*dt, len(data_mean))
                elif callable(xvals):
                    xvals_array = xvals(data_mean, **xvals_kw)
                else:
                    xvals_array=xvals

                assert 'error_style' in figdf, 'specify error style'
                # line plot with shaded error
                #-----------------------------------
                if figdf.loc[(figkey,subkey,tracekey)].error_style=='shade':
                    # plot mean trace
                    #-------------------------------
                    ax[figkey].plot(xvals_array, data_mean, 
                        color=params.trace_color, 
                        linewidth=params.trace_linewidth,
                        linestyle=params.trace_linestyle)
                    # shade error
                    #-----------------------------------
                    plt.fill_between(xvals_array,data_mean-data_sem,data_mean+data_sem,color=params.trace_ecolor, 
                        alpha=params.trace_ealpha)

                # error bar plot
                #--------------------------------------
                elif figdf.loc[(figkey,subkey,tracekey)].error_style=='bar':
                    ax[figkey].errorbar(xvals_array, data_mean, 
                        yerr=data_sem, 
                        color=params.trace_color, 
                        marker=params.trace_marker,  
                        markersize=params.markersize, 
                        elinewidth=params.error_linewidth, 
                        linewidth=params.trace_linewidth, 
                        markerfacecolor=params.trace_color, 
                        ecolor=params.trace_ecolor)
        # format figure
        #----------------
        if figformat=='standard':
            fig, ax = FormatFig()._standard_figformat(figures=fig, axes=ax, figdf=figdf)

        plt.show(block=False)

        return fig, ax

    def _var2var_mean(self, df, figdf, x_variable, variable, array_funcs=[], array_func_kws=[], df_funcs=[], df_func_kws=[], df_sorted=[], figformat='standard', **kwargs):
        '''FIXME add docs
        '''
        if 'dt' in kwargs:
            dt=kwargs['dt']
        else:
            dt=1
        # report progress
        print 'plotting:', inspect.stack()[0][3]
        # preallocate figures and axes
        fig={}
        ax={}
        # iteratively apply functions for creating new column to dataframe
        #----------------------------------------------------------------
        if len(df_funcs)>0:
            for df_func_i, df_func in enuemrate(df_funcs):
                df = df_func(df, axis=1, **df_func_kws[df_func_i])
        # set figdf to hierarchical index (figure, subgroup, trace)
        figdf = figdf.reset_index().set_index(['figure','subgroup','trace'])
        # get list of all figure names
        figures = figdf.index.get_level_values('figure').unique().values
        # iterate over figures
        for figkey in figures:
            # create figure, passing params as **kwargs
            fig[figkey], ax[figkey] = plt.subplots()
            # get subgroups list
            subgroups = figdf.loc[figkey].index.get_level_values('subgroup').unique().values
            # iterate over subgroups of traces
            for subkey in subgroups:
                traces = figdf.loc[(figkey, subkey)].index.get_level_values('trace').unique().values
                # iterate over traces
                for tracekey in traces:
                    # get params from figdf
                    #-----------------------
                    params = figdf.loc[(figkey,subkey, tracekey)]
                    # get corresponding series and convert to array
                    #----------------------------------------------
                    # if available, get series, otherwise next iter
                    try:
                        trace_series=df.loc[tracekey,slice(None)][[x_variable, variable]]
                    except:
                        print tracekey,'or', variable, 'or', x_variable, 'not in df index'
                        continue

                    # group data by x_variable
                    #-------------------------------------------------
                    trace_series = trace_series.set_index(x_variable)

                    # iterate over x_variable values
                    #------------------------------------------------
                    x_vals = np.sort(trace_series.index.unique().values)
                    y_vals = np.full_like(x_vals, np.nan, dtype=np.double)
                    e_vals = np.full_like(x_vals, np.nan, dtype=np.double)
                    for i, x_val in enumerate(x_vals):
                        print trace_series.loc[x_val, variable]
                        # get corresponding y data and convert to array
                        #----------------------------------------------
                        data_array = _2array(trace_series.loc[x_val, variable], remove_nans=True, remove_nans_axis=1)
                        print 'data array shape:', data_array.shape
                        # iteratively apply array functions
                        #-----------------------------------
                        if len(array_funcs)>0:
                            for array_func_i, array_func in enumerate(array_funcs):
                                data_array=array_func(data_array, **array_func_kws[array_func_i])
                        # make sure array is the correct shape
                        #--------------------------------------
                        if type(data_array)==np.ndarray and data_array.shape[0]!=0:
                            assert len(data_array.shape)==1, 'data array must be 1d'
                            # mean across slices
                            data_mean = np.mean(data_array, axis=0)
                            #std across slices
                            data_std = np.std(data_array, axis=0)
                            # sem across slices
                            data_sem = stats.sem(data_array, axis=0)
                            # update y values and errors for plot
                            y_vals[i]=data_mean
                            e_vals[i]=data_sem

                    print 'x',x_vals
                    print 'y',y_vals
                    print 'e',e_vals
                    # line plot with shaded error
                    #-----------------------------------
                    if figdf.loc[(figkey,subkey,tracekey)].error_style=='shade':
                        ax[figkey].plot(x_vals, y_vals ,color=params.trace_color, linewidth=params.trace_linewidth)
                        # print params.error_color
                        plt.fill_between(x_vals, y_vals-e_vals, y_vals+e_vals, color=params.trace_ecolor, alpha=params.trace_ealpha)
                    # error bar plot
                    #--------------------------------------
                    elif figdf.loc[(figkey,subkey,tracekey)].error_style=='bar':
                        ax[figkey].plot(x_vals, y_vals, 
                            color=params.trace_color, 
                            marker=params.trace_marker,  
                            markersize=params.trace_markersize,  
                            linewidth=params.trace_linewidth, 
                            markerfacecolor=params.trace_color)
                        ax[figkey].errorbar(x_vals, y_vals, yerr=e_vals, color=params.trace_color, 
                            marker=params.trace_marker,  
                            markersize=params.trace_markersize, 
                            elinewidth=params.error_linewidth, 
                            linewidth=params.trace_linewidth, 
                            markerfacecolor=params.trace_color, 
                            ecolor=params.trace_ecolor)

        # fig, ax = _format_figures(fig=fig, ax=ax, figdf=figdf,)
        # plt.axhline(y=1, color='black', linestyle='--', linewidth=4)
        # format figure
        #----------------
        if figformat=='standard':
            fig, ax = FormatFig()._standard_figformat(figures=fig, axes=ax, figdf=figdf)

        # plt.show(block=False)

        return fig, ax

    def _var2var_func(self, df, figdf, x_variable, variable, figformat='standard', array_funcs=[], array_func_kws=[], error_funcs=[], error_func_kws=[],df_funcs=[], df_func_kws=[], df_sorted=[],**kwargs):
        '''FIXME add docs
        '''
        if 'dt' in kwargs:
            dt=kwargs['dt']
        else:
            dt=1
        # report progress
        print 'plotting:', inspect.stack()[0][3]
        # preallocate figures and axes
        fig={}
        ax={}
        # iteratively apply functions for creating new column to dataframe
        #----------------------------------------------------------------
        if len(df_funcs)>0:
            for df_func_i, df_func in enuemrate(df_funcs):
                df = df_func(df, axis=1, **df_func_kws[df_func_i])
        # set figdf to hierarchical index (figure, subgroup, trace)
        figdf = figdf.reset_index().set_index(['figure','subgroup','trace'])
        # get list of all figure names
        figures = figdf.index.get_level_values('figure').unique().values
        # iterate over figures
        for figkey in figures:
            # create figure, passing params as **kwargs
            fig[figkey], ax[figkey] = plt.subplots()
            # get subgroups list
            subgroups = figdf.loc[figkey].index.get_level_values('subgroup').unique().values
            # iterate over subgroups of traces
            for subkey in subgroups:
                traces = figdf.loc[(figkey, subkey)].index.get_level_values('trace').unique().values
                # iterate over traces
                for tracekey in traces:
                    # get params from figdf
                    #-----------------------
                    params = figdf.loc[(figkey,subkey, tracekey)]
                    # get corresponding series and convert to array
                    #----------------------------------------------
                    # if available, get series, otherwise next iter
                    try:
                        # trace_series=df.loc[tracekey,slice(None)][[x_variable, variable]]
                        trace_series=df.loc[tracekey,slice(None)]
                    except:
                        print tracekey,'or', variable, 'or', x_variable, 'not in df index'
                        continue

                    # group data by x_variable
                    #-------------------------------------------------
                    trace_series = trace_series.reset_index().set_index(x_variable)

                    # iterate over x_variable values
                    #------------------------------------------------
                    x_vals = np.sort(trace_series.index.unique().values)
                    # print tracekey
                    # print 'x_vals', x_vals
                    y_vals = np.full_like(x_vals, np.nan, dtype=np.double)
                    e_vals = np.full_like(x_vals, np.nan, dtype=np.double)
                    error_array=0
                    for i, x_val in enumerate(x_vals):
                        # print trace_series
                        # print trace_series.loc[x_val, variable]
                        # get corresponding y data and convert to array
                        #----------------------------------------------
                        data_array = _2array(trace_series.loc[x_val, variable], remove_nans=True, remove_nans_axis=1)
                        
                        # iteratively apply array functions
                        #-----------------------------------
                        if len(error_funcs)>0:
                            for error_func_i, error_func in enumerate(error_funcs):
                                error_array=error_func(data_array, **error_func_kws[error_func_i])

                        print 'data array shape:', data_array.shape
                        # iteratively apply array functions
                        #-----------------------------------
                        if len(array_funcs)>0:
                            for array_func_i, array_func in enumerate(array_funcs):
                                data_array=array_func(data_array, **array_func_kws[array_func_i])

                        # make sure array is the correct shape
                        #--------------------------------------
                        try:
                            float(data_array)
                        except:
                            continue

                        # print type(data_array)
                        # if (type(data_array)==np.ndarray and len(data_array.shape)==0) or type(data_array)==float:
                            # assert len(data_array.shape)==1, 'data array must be 1d'
                            # mean across slices
                            # data_mean = np.mean(data_array, axis=0)
                            # #std across slices
                            # data_std = np.std(data_array, axis=0)
                            # # sem across slices
                            # data_sem = stats.sem(data_array, axis=0)
                            # update y values and errors for plot
                        y_vals[i]=data_array
                        e_vals[i]=error_array

                    x_vals = [float(val) for val in x_vals]
                    print 'x',x_vals
                    print 'y',y_vals

                    # print 'e',e_vals
                    # line plot with shaded error
                    #-----------------------------------
                    if figdf.loc[(figkey,subkey,tracekey)].error_style=='shade':
                        ax[figkey].plot(x_vals, y_vals ,color=params.trace_color, linewidth=params.trace_linewidth)
                        # print params.error_color
                        # plt.fill_between(x_vals, y_vals-e_vals, y_vals+e_vals, color=params.trace_ecolor, alpha=params.trace_ealpha)
                    # error bar plot
                    #--------------------------------------
                    elif figdf.loc[(figkey,subkey,tracekey)].error_style=='bar':
                        ax[figkey].plot(x_vals, y_vals, 
                            color=params.trace_color, 
                            marker=params.trace_marker,  
                            markersize=params.trace_markersize,  
                            linewidth=params.trace_linewidth, 
                            markerfacecolor=params.trace_color)

                        print x_vals, y_vals, e_vals
                        ax[figkey].errorbar(x_vals, y_vals, yerr=e_vals, color=params.trace_color, 
                            marker=params.trace_marker,  
                            markersize=params.trace_markersize, 
                            elinewidth=params.error_linewidth, 
                            linewidth=params.trace_linewidth, 
                            markerfacecolor=params.trace_color, 
                            ecolor=params.trace_ecolor)

        # fig, ax = _format_figures(fig=fig, ax=ax, figdf=figdf,)
        # plt.axhline(y=1, color='black', linestyle='--', linewidth=4)
        # figure formatting
        #--------------------
        # fig, ax = FormatFig()._set_fontsizes(figures=fig, axes=ax, figdf=figdf)
        # fig, ax = FormatFig()._set_yticks(figures=fig, axes=ax, figdf=figdf)
        # fig, ax = FormatFig()._scale_ticklabels(figures=fig, axes=ax, figdf=figdf)
        # format figure
        #----------------
        if figformat=='standard':
            fig, ax = FormatFig()._standard_figformat(figures=fig, axes=ax, figdf=figdf)
        plt.show(block=False)

        return fig, ax

    def _trace_individual(self, df, figdf, variable, index=0, array_funcs=[], array_func_kws=[], df_funcs=[], df_func_kws=[], df_sorted=[],**kwargs):
        '''FIXME add docs
        '''
        # set dt
        #---------
        if 'dt' in kwargs:
            dt=kwargs['dt']
        else:
            dt=1
        # report progress
        #------------------
        print 'plotting:', inspect.stack()[0][3]
        # preallocate figures and axes
        fig={}
        ax={}
        # iteratively apply functions for creating new column to dataframe
        #----------------------------------------------------------------
        if len(df_funcs)>0:
            for df_func_i, df_func in enuemrate(df_funcs):
                df = df_func(df, axis=1, **df_func_kws[df_func_i])
        # set figdf to hierarchical index (figure, subgroup, trace)
        figdf = figdf.reset_index().set_index(['figure','subgroup','trace'])
        # get list of all figure names
        figures = figdf.index.get_level_values('figure').unique().values
        # iterate over figures
        for figkey in figures:
            # create figure, passing params as **kwargs
            fig[figkey], ax[figkey] = plt.subplots()
            # get subgroups list
            subgroups = figdf.loc[figkey].index.get_level_values('subgroup').unique().values
            # iterate over subgroups of traces
            for subkey in subgroups:
                traces = figdf.loc[(figkey, subkey)].index.get_level_values('trace').unique().values
                # iterate over traces
                for tracekey in traces:
                    # get params from figdf
                    #-----------------------
                    params = figdf.loc[(figkey,subkey, tracekey)]
                    # get corresponding series and convert to array
                    #----------------------------------------------
                    # if available, get series, otherwise next iter
                    try:
                        trace_series=df.loc[tracekey,slice(None)][variable]
                    except:
                        print tracekey,'or', variable, 'not in df index'
                        continue
                    # convert to array
                    data_array = _2array(trace_series, remove_nans=True, remove_nans_axis=1)
                    print 'data array shape:', data_array.shape
                    # iteratively apply array functions
                    #-----------------------------------
                    if len(array_funcs)>0:
                        for array_func_i, array_func in enumerate(array_funcs):
                            data_array=array_func(data_array, **array_func_kws[array_func_i])
                    # make sure array is the correct shape
                    if type(data_array)==np.ndarray and data_array.shape[0]!=0:
                        # if array is 1d, convert to 2d
                        if len(data_array.shape)==1:
                            data_array = data_array.reshape((1,-1))
                        # mean across slices
                        data_mean = np.mean(data_array, axis=0)
                        #std across slices
                        data_std = np.std(data_array, axis=0)
                        # sem across slices
                        data_sem = stats.sem(data_array, axis=0)
                        # time vector
                        t = np.arange(0, len(data_mean)*dt, dt)

                        assert 'error_style' in figdf, 'specify error style'
                        # line plot with shaded error
                        #-----------------------------------
                        if figdf.loc[(figkey,subkey,tracekey)].error_style=='shade':
                            ax[figkey].plot(t, data_mean, color=params.trace_color, linewidth=params.trace_linewidth)
                            # print params.error_color
                            plt.fill_between(t, data_mean-data_sem, data_mean+data_sem, color=params.trace_ecolor, alpha=params.trace_ealpha)
                        # error bar plot
                        #--------------------------------------
                        elif figdf.loc[(figkey,subkey,tracekey)].error_style=='bar':
                            ax[figkey].errorbar(t, data_mean, yerr=data_sem, color=params.trace_color, 
                                marker=params.trace_marker,  
                                markersize=params.markersize, 
                                elinewidth=params.error_linewidth, 
                                linewidth=params.trace_linewidth, 
                                markerfacecolor=params.trace_color, 
                                ecolor=params.trace_ecolor)

        fig, ax = _format_figures(fig=fig, ax=ax, figdf=figdf,)

        plt.show(block=False)

        return fig, ax

    def _hist(self, df, figdf, variables, bins='auto', bins_kw={}, normalize_weights=False, cumulative=False, histtype='step', orientation='vertical', array_funcs=[], array_func_kws=[], df_funcs=[], df_func_kws=[],**kwargs):
        '''FIXME add docs
        '''
        print 'plotting:', inspect.stack()[0][3]
        fig={}
        ax={}
        # iteratively apply functions for creating new column to dataframe
        #----------------------------------------------------------------
        if len(df_funcs)>0:
            for df_func_i, df_func in enuemrate(df_funcs):
                df = df_func(df, axis=1, **df_func_kws[df_func_i])
        # iterate over traces to be plotted
        #-----------------------------------
        for figkey, subkey, tracekey, params, trace_df in self._figdf_generator(df, figdf, variables):
            # make new figure if doesnt exist yet
            #-------------------------------------
            if figkey not in ax:
                fig[figkey], ax[figkey] = plt.subplots()
            # convert df to array
            #---------------------
            data_array = _2array(trace_df, remove_nans=True, remove_nans_axis=1)
            print data_array
            print 'data array shape:', data_array.shape

            # iteratively apply array functions
            #-----------------------------------
            if len(array_funcs)>0:
                for array_func_i, array_func in enumerate(array_funcs):
                    data_array=array_func(data_array, **array_func_kws[array_func_i])

            # make sure array is the correct shape
            #-------------------------------------
            if type(data_array)==np.ndarray and data_array.shape[0]!=0:
                if len(data_array.shape)!=1:
                    data_array = data_array.flatten()

                # set kwargs for plotting function
                #-----------------------------------
                # kwargs_hist = self._select_kwargs([plt.hist], params.trace_kwargs)

                # weights
                #--------------------------------
                if normalize_weights:
                    weights=np.ones_like(data_array)/float(len(data_array))
                    # kwargs_hist['weights']=weights
                else:
                    weights=None

                # bins
                #--------------------------------------
                if callable(bins):
                    bin_array = bins(data_array, **bins_kw)
                    # kwargs_hist['bins']=bin_array
                else:
                    bin_array=bins
                    # kwargs_hist['bins']=bins

                # plot
                #--------------------------------------------
                # print kwargs_kosher
                # ax[figkey].hist(data_array, linewidth=params.trace_linewidth, **kwargs_hist)

                # weights=np.ones_like(data_array)
                # weights=np.ones_like(data_array)/float(len(data_array))
                # bins = 100
                # bins='auto'
                # data_range = [min(data_array), max(data_array)]
                # print data_range
                # bins =int(np.floor(len(data_array)/4))
                # bins = np.arange(0, data_range[1], .01)
                # bins = np.arange(0, 60, .5)
                # bins='auto'
                # print data_array
                ax[figkey].hist(data_array, color=params.trace_color, bins=bin_array, cumulative=cumulative, histtype=histtype, weights=weights, orientation=orientation, linewidth=params.trace_linewidth)
        # format figures
        #-------------------------------------------------------------------
        fig, ax = FormatFig()._standard_figformat(figures=fig, axes=ax, figdf=figdf, tight=False)

        plt.show(block=False)

        return fig, ax

    def _shapeplot(self, df, figdf, variable, array_funcs=[], array_func_kws=[], df_funcs=[], df_func_kws=[], df_sorted=[], cmap=colormap.PiYG, **kwargs):
        '''
        '''
        print 'plotting:', inspect.stack()[0][3]
        fig={}
        ax={}
        # set figdf to hierarchical index (figure, subgroup, trace)
        figdf = figdf.reset_index().set_index(['figure','subgroup','trace'])
        # get list of all figure names
        figures = figdf.index.get_level_values('figure').unique().values
        df = df.reset_index().set_index('field')
        # iterate over figures
        for figkey in figures:
            # create figure, passing params as **kwargs
            fig[figkey], ax[figkey] = plt.subplots()
            # get subgroups list
            subgroups = figdf.loc[figkey].index.get_level_values('subgroup').unique().values
            # iterate over subgroups of traces
            for subkey in subgroups:
                traces = figdf.loc[(figkey, subkey)].index.get_level_values('trace').unique().values
                # iterate over traces
                for tracekey in traces:
                    # get params from figdf
                    params = figdf.loc[(figkey,subkey, tracekey)]
                    # get corresponding series and convert to array
                    #----------------------------------------------
                    # if available, get series, otherwise next iter
                    try:
                        data=df.loc[tracekey,slice(None)][variable].values
                        morpho=df.loc[tracekey, slice(None)]['morpho'].values
                    except:
                        print tracekey,'or', variable, 'not in df index'
                        continue

                    patches, colors = ShapePlot().basic(morpho=morpho, data=data, axes=ax[figkey], width_scale=3, colormap=cmap)

                    # plot collection
                    ax[figkey].add_collection(patches)
                    # show colorbar
                    plt.colorbar(patches)
                    # autoscale axes
                    ax[figkey].autoscale()
                    ax[figkey].set_aspect(1.)
                    plt.axis('off')
        plt.show(block=False)

        return fig, ax

    def _bar(self, df, figdf, variable, array_funcs=[], array_func_kws=[], df_funcs=[], df_func_kws=[], group_space=1, bar_width=1, bar_spacing=1, **kwargs):
        '''
        '''
        # preallocate variables to be passed to barplot
        #-----------------------------------------------
        fig={}
        ax={}
        n_subgroups={}
        n_traces={}
        xlim={}
        ylim={}
        data={}
        heights={}
        locations={}
        y_errors={}
        tracelist={}
        # set figdf to hierarchical index (figure, subgroup, trace)
        figdf = figdf.reset_index().set_index(['figure','subgroup','trace'])
        # get list of all figure names
        figures = figdf.index.get_level_values('figure').unique().values
        # get pairwise comparisons between traces
        # ttests  = Stats()._pairwise_ttests(df_sorted, variable, array_funcs=array_funcs, array_func_kws=array_func_kws)
        # iterate over figures
        #--------------------
        for figkey in figures:
            # create figure objects
            fig[figkey], ax[figkey] = plt.subplots()
            # get subgroups list
            subgroups = figdf.loc[figkey].index.get_level_values('subgroup').unique().values
            # preallocate
            #---------------
            locations[figkey] = []
            tracelist[figkey]=[]
            heights[figkey]=[]
            y_errors[figkey]=[]
            colors=[]
            xticks=[]
            xticklabels=[]
            cnt=bar_spacing
            # iterate over subgroups of traces
            #---------------------------------
            for subkey in subgroups:
                # get traces
                traces = figdf.loc[(figkey, subkey)].index.get_level_values('trace').unique().values
                # count subgroups
                cnt+=group_space
                # iterate over traces
                #---------------------
                for tracekey in traces:
                    # trace parameters
                    param = figdf.loc[(figkey,subkey, tracekey)]
                    # get corresponding series and convert to array
                    #----------------------------------------------
                    # if available, get series, otherwise next iter
                    try:
                        trace_series=df.loc[tracekey,slice(None)][variable]
                    except:
                        print tracekey,'or', variable, 'not in df index'
                        continue
                    # convert to array
                    data_array = _2array(trace_series, remove_nans=True, remove_nans_axis=1)
                    print 'data array shape:', data_array.shape
                    # iteratively apply array functions
                    #-----------------------------------
                    if len(array_funcs)>0:
                        for array_func_i, array_func in enumerate(array_funcs):
                            data_array=array_func(data_array, **array_func_kws[array_func_i])
                            print array_func.__name__, data_array.shape
                    print 'data array shape:', data_array.shape
                    print tracekey, data_array
                    # make sure array is the correct shape (1 dimensional with at least one entry in first dimension)
                    #-------------------------------------------
                    if not (type(data_array)==np.ndarray and len(data_array.shape)==1 and data_array.shape[0]!=0):
                        continue
                    # pdb.set_trace()
                    # print figkey
                    # print data_array
                    # mean across slices
                    data_mean = np.mean(data_array, axis=0)
                    #std across slices
                    data_std = np.std(data_array, axis=0)
                    # sem across slices
                    data_sem = stats.sem(data_array, axis=0)
                    # add plot location
                    # print figkey, subkey, tracekey, param.sub_location, param.trace_location
                    print data_mean
                    plot_location = param.sub_location+param.trace_location
                    tracelist[figkey].append(tracekey)
                    locations[figkey].append(plot_location)
                    xticks.append(plot_location)
                    xticklabels.append(param.trace_label)
                    colors.append(param.trace_color)
                    heights[figkey].append(data_mean)
                    y_errors[figkey].append(data_sem)
                    print tracekey, data_mean.shape
                    plt.errorbar(locations[figkey][-1], data_mean, data_sem, color=param.trace_ecolor)
            
            # if len(heights[figkey])>0:
            barcontainer = ax[figkey].bar(locations[figkey], heights[figkey], width=param.fig_barwidth, tick_label=xticklabels)

            # get x and y lims
            xlim[figkey] = ax[figkey].get_xlim()
            ylim[figkey] = ax[figkey].get_ylim()
        
            # set bar color
            for bar_i, bar in enumerate(barcontainer):
                bar.set_color(colors[bar_i])

        # format figure
        #----------------
        fig, ax = FormatFig()._standard_figformat(figures=fig, axes=ax, figdf=figdf, xticks=False, xscale=False, xticks_minor=True)

        plt.show(block=False)

        return fig, ax



    # deprecated plot functions
    #----------------------------------------------------------------------
    def _trace_mean_deprecated(self, df, figdf, variable, array_funcs=[], array_func_kws=[], df_funcs=[], df_func_kws=[], df_sorted=[],**kwargs):
        '''FIXME add docs
        '''
        if 'dt' in kwargs:
            dt=kwargs['dt']
        else:
            dt=1
        # report progress
        print 'plotting:', inspect.stack()[0][3]
        # preallocate figures and axes
        fig={}
        ax={}
        # iteratively apply functions for creating new column to dataframe
        #----------------------------------------------------------------
        if len(df_funcs)>0:
            for df_func_i, df_func in enuemrate(df_funcs):
                df = df_func(df, axis=1, **df_func_kws[df_func_i])

        # set figdf to hierarchical index (figure, subgroup, trace)
        figdf = figdf.reset_index().set_index(['figure','subgroup','trace'])
        # get list of all figure names
        figures = figdf.index.get_level_values('figure').unique().values
        # iterate over figures
        for figkey in figures:
            # create figure, passing params as **kwargs
            fig[figkey], ax[figkey] = plt.subplots()
            # get subgroups list
            subgroups = figdf.loc[figkey].index.get_level_values('subgroup').unique().values
            # iterate over subgroups of traces
            for subkey in subgroups:
                traces = figdf.loc[(figkey, subkey)].index.get_level_values('trace').unique().values
                # iterate over traces
                for tracekey in traces:
                    # get params from figdf
                    #-----------------------
                    params = figdf.loc[(figkey,subkey, tracekey)]
                    # get corresponding series and convert to array
                    #----------------------------------------------
                    # if available, get series, otherwise next iter
                    try:
                        trace_series=df.loc[tracekey,slice(None)][variable]
                    except:
                        print tracekey,'or', variable, 'not in df index'
                        continue
                    # convert to array
                    data_array = _2array(trace_series, remove_nans=True, remove_nans_axis=1)
                    print 'data array shape:', data_array.shape

                    # iteratively apply array functions
                    #-----------------------------------
                    if len(array_funcs)>0:
                        for array_func_i, array_func in enumerate(array_funcs):
                            data_array=array_func(data_array, **array_func_kws[array_func_i])
                    # make sure array is the correct shape
                    if type(data_array)==np.ndarray and data_array.shape[0]!=0:
                        # if array is 1d, convert to 2d
                        if len(data_array.shape)==1:
                            data_array = data_array.reshape((1,-1))
                        # mean across slices
                        data_mean = np.mean(data_array, axis=0)
                        #std across slices
                        data_std = np.std(data_array, axis=0)
                        # sem across slices
                        data_sem = stats.sem(data_array, axis=0)
                        # time vector
                        # if 't' in kwargs:
                        #     t=kwargs['t']
                        # else:
                        # t = np.linspace(0, len(data_mean)*dt, len(data_mean))
                        t = dt*(np.arange(0, len(data_mean))-(len(data_mean)-1)/2)


                        assert 'error_style' in figdf, 'specify error style'
                        # line plot with shaded error
                        #-----------------------------------
                        if figdf.loc[(figkey,subkey,tracekey)].error_style=='shade':
                            ax[figkey].plot(t, data_mean, color=params.trace_color, linewidth=params.trace_linewidth,linestyle=params.trace_linestyle)
                            # print params.error_color
                            plt.fill_between(t, data_mean-data_sem, data_mean+data_sem, color=params.trace_ecolor, alpha=params.trace_ealpha)
                        # error bar plot
                        #--------------------------------------
                        elif figdf.loc[(figkey,subkey,tracekey)].error_style=='bar':
                            ax[figkey].errorbar(t, data_mean, yerr=data_sem, color=params.trace_color, 
                                marker=params.trace_marker,  
                                markersize=params.markersize, 
                                elinewidth=params.error_linewidth, 
                                linewidth=params.trace_linewidth, 
                                markerfacecolor=params.trace_color, 
                                ecolor=params.trace_ecolor)
        # format figure
        #----------------
        fig, ax = FormatFig()._standard_figformat(figures=fig, axes=ax, figdf=figdf)

        # fig, ax = _format_figures(fig=fig, ax=ax, figdf=figdf,)

        plt.show(block=False)

        return fig, ax

    def _hist_deprecated(self, df, figdf, variable, array_funcs=[], array_func_kws=[], df_funcs=[], df_func_kws=[], df_sorted=[],**kwargs):
        '''FIXME add docs
        '''
        print 'plotting:', inspect.stack()[0][3]
        fig={}
        ax={}
        xlim={}
        ylim={}
        # iteratively apply functions for creating new column to dataframe
        #----------------------------------------------------------------
        if len(df_funcs)>0:
            for df_func_i, df_func in enuemrate(df_funcs):
                df = df_func(df, axis=1, **df_func_kws[df_func_i])
        # set figdf to hierarchical index (figure, subgroup, trace)
        figdf = figdf.reset_index().set_index(['figure','subgroup','trace'])
        # get list of all figure names
        figures = figdf.index.get_level_values('figure').unique().values
        # iterate over figures
        for figkey in figures:
            # create figure, passing params as **kwargs
            fig[figkey], ax[figkey] = plt.subplots()
            # get subgroups list
            subgroups = figdf.loc[figkey].index.get_level_values('subgroup').unique().values
            # iterate over subgroups of traces
            for subkey in subgroups:
                traces = figdf.loc[(figkey, subkey)].index.get_level_values('trace').unique().values
                # iterate over traces
                for tracekey in traces:
                    # get params from figdf
                    params = figdf.loc[(figkey,subkey, tracekey)]
                    # get corresponding series and convert to array
                    #----------------------------------------------
                    # if available, get series, otherwise next iter
                    try:
                        trace_series=df.loc[tracekey,slice(None)][variable]
                    except:
                        print tracekey,'or', variable, 'not in df index'
                        continue
                    # remove nans from trace series
                    na = trace_series.isna()
                    trace_series = trace_series[~na]

                    # convert to array
                    data_array = _2array(trace_series, remove_nans=True, remove_nans_axis=1)
                    print data_array
                    print 'data array shape:', data_array.shape
                    # iteratively apply array functions
                    #-----------------------------------
                    if len(array_funcs)>0:
                        for array_func_i, array_func in enumerate(array_funcs):
                            data_array=array_func(data_array, **array_func_kws[array_func_i])

                    print 'data array shape:', data_array.shape
                    # make sure array is the correct shape
                    if type(data_array)==np.ndarray and data_array.shape[0]!=0:
                        if len(data_array.shape)!=1:
                            data_array = data_array.flatten()
                        weights=np.ones_like(data_array)
                        # weights=np.ones_like(data_array)/float(len(data_array))
                        # bins = 100
                        # bins='auto'
                        data_range = [min(data_array), max(data_array)]
                        # print data_range
                        # bins =int(np.floor(len(data_array)/4))
                        # bins = np.arange(0, data_range[1], .01)
                        bins = np.arange(0, 60, .5)
                        # bins='auto'
                        # print data_array
                        ax[figkey].hist(data_array, color=params.trace_color, bins=bins, cumulative=False, histtype='step', weights=weights, orientation='vertical', linewidth=params.trace_linewidth)
                        # ax[figkey].hist(data_array,  bins=bins, cumulative=False, histtype='step', weights=weights, orientation='vertical', params.trace_kwargs)
                        # ax[figkey].invert_xaxis()
                        # plt.title(figkey)

        # fig, ax = FormatFig()._set_fontsizes(figures=fig, axes=ax, figdf=figdf)
        # fig, ax = FormatFig()._set_yticks(figures=fig, axes=ax, figdf=figdf)
        # fig, ax = FormatFig()._scale_ticklabels(figures=fig, axes=ax, figdf=figdf)
        fig, ax = FormatFig()._standard_figformat(figures=fig, axes=ax, figdf=figdf)

        plt.show(block=False)

        return fig, ax
#############################################################################
# array functions
#############################################################################
class ArrayFunctions:
    '''
    '''
    def __init__(self, ):
        '''
        '''
        pass

    def _normalize(self, array, norm_type='divide', norm_slice=slice(1), axis=1, axis_slice=slice(None), **kwargs):
        '''
        '''
        # create slice object with appropriate number of dimensions
        slicer = [slice(None)]*len(array.shape)
        # update slice object
        slicer[axis]=axis_slice
        # slice array
        array = array[slicer]
        # get normalization slicer
        # create slice object
        norm_slicer = copy.deepcopy(slicer)
        norm_slicer[axis] = norm_slice
        if norm_type=='divide':
            print slicer
            print norm_slicer
            array = array/array[norm_slicer]
        elif norm_type=='subtract':
            array = array-array[norm_slicer]

        return array

    def _get_area(self, array, axis=1, islice=slice(None), norm=True, **kwargs):
        ''' get area under curve for timeseries

        Arguments
        ----------
        ~array : array of timeseries, typically cells x samples
        ~axis : axis corresponding to time (to take area under)
        ~islice : python slice object for choosing indices along time dimension
        ~norm : bool indicating whether to normalize values to the first value (after slicing)

        Returns
        --------
        ~area : vector of areas for each cell

        Notes
        --------
        '''
        # create slice object with appropriate number of dimensions
        slicer = [slice(None)]*len(array.shape)
        # update slice object
        slicer[axis]=islice
        # slice array
        array = array[slicer]
        # if  normalizing to first value
        if norm:
            # create slice object
            norm_slicer = copy.deepcopy(slicer)
            norm_slicer[axis] = slice(1)
            # normalize
            array = array-array[norm_slicer]
            # print array

        area = np.sum(array, axis=axis)
        if len(area.shape)==0:
            area=np.reshape(area, (1))


        return area

    def _get_bins(self, data_array, nbins=100, binsize=None, binmin=None, binmax=None):
        ''' generate an array of bins for histograms based on data in data_array

            Arguments
            ------------
            ~data_array: 1d array of values to be binned for histogram
            ~nbins: total number of bins
            ~binsize: step size between bins
            ~binmin: minimum bin value (not size)
            ~binmax: maximum bin value (not size)

            Return
            --------------
            ~bins: array of bin values to be passed to plt.hist

            Comments
            ------------------
            ~data_array must be 1D
            ~binsize overrides nbins if specified
        '''
        print data_array
        data_range = [np.min(data_array), np.max(data_array)]
        print 'data range',data_range
        if binmin is None:
            binmin =data_range[0]
        if binmax is None:
            binmax =data_range[1]
        if binsize is not None:
            if binmin==binmax:
                binmin=binmin-binsize
                binmax=binmax+binsize
            bins = np.arange(binmin, binmax, binsize)
        else:
            bins = np.linspace(binmin, binmax, nbins)

        return bins

    def _get_trace_xvals(self, data_array, dt=1, xtype='time', **kwargs):
        '''
        '''
        print 'getting xvals'
        print data_array
        print data_array.shape
        if xtype=='time':
            xvals_array = np.linspace(0, len(data_array)*dt, len(data_array))
        elif xtype=='xcorr':
            xvals_array = dt*(np.arange(0, len(data_array))-(len(data_array)-1)/2)

        print xvals_array
        print xvals_array.shape

        return xvals_array

    def _filter_mean(self, array, filt=[], filt_series=[], hilbert=False, mean_axis=0, islice=slice(None), slice_axis=1 , **kwargs):
        '''
        '''
        # create slice object with appropriate number of dimensions
        slicer = [slice(None)]*len(array.shape)
        # update slice object
        slicer[slice_axis]=islice
        array_mean = np.mean(array[slicer], axis=mean_axis)

        # single filter
        if filt:
            filtered_array = signal.filtfilt(filt[0], filt[1], array_mean, method='pad', padlen=len(array), padtype='odd')
            # filtered_array = signal.lfilter(filt[0], filt[1], array_mean,)
        # cascade filters
        if filt_series:
            filtered_array=copy.copy(array)
            for subfilt_key, subfilt in filt_series.iteritems():
                filtered_array = signal.filtfilt(subfilt[0], subfilt[1], filtered_array, method='pad', padlen=len(filtered_array),padtype='odd')
                # filtered_array = signal.lfilter(subfilt[0], subfilt[1], filtered_array,)
        # apply hilbert
        if hilbert:
            filtered_array = np.abs(signal.hilbert(filtered_array - np.mean(filtered_array)))

        print 'filtered_mean',np.mean(filtered_array)
        # pdb.set_trace()
        # plt.figure()
        # plt.plot(filtered_array)
        # plt.show(block=False)

        return filtered_array

    def _subtract_timeseries(self, array, islice=slice(1), axis=1, **kwargs):
        '''
        '''
        # create slice object with appropriate number of dimensions
        slicer = [slice(None)]*len(array.shape)
        # update slice object
        slicer[axis]=islice
        # get values to subtract
        sub_vals = copy.copy(array[slicer])
        # subtract values
        array_new = array-sub_vals
        return array_new

    def _flatten(self, array, **kwargs):
        '''
        '''
        flattened = np.array(list(itertools.chain.from_iterable(array)))
        return flattened

    def _reject_outliers(self, array, std_tol=2, axis=0, **kwargs):
        '''
        '''
        # if 'std_tol' in kwargs:
        #     std_tol=kwargs['std_tol']
        # array_mean = np.mean(array, axis=axis)
        # array_std = np.std(array, axis=axis)
        # array_diff = array-array_mean

        array_new = array[abs(array- np.mean(data, axis=axis)) < m * np.std(array, axis=axis)]
        return array_new

    def _ppr(self, array, **kwargs):
        '''
        '''
        ppr = array[:,1]/array[:,0]
        return ppr

    def _last(self, array, **kwargs):
        '''
        '''
        last = array[:,-1]

        return last

    def _first(self, array, **kwargs):
        '''
        '''
        # if data is a vector, convert to n x 1 array
        if len(array.shape)==1:
            data = np.reshape(array, (1,len(array)))
        else:
            data=array
        first = data[:,0]
        return first

    def _slice_mean(self, array, islice, axis=2, **kwargs):
        '''
        '''
        # assert 'axis' in kwargs, 'please specificy axis to operate on'
        # assert '' in kwargs, 'please specificy axis to operate on'
        # axis = kwargs['axis']
        slicer = [slice(None)]*len(array.shape)
        # pdb.set_trace()
        # print slicer
        # print islice
        # print array.shape
        if axis<len(array.shape) and array.shape[axis]>0:
            slicer[axis]=islice
            # print slicer
            # print array[slicer]
            array_new = np.mean(array[slicer], axis=axis)
        else:
            array_new = array
        # array_new = np.mean(array[slicer], axis=axis)
        # else:
        #     array_new=array
        # print array_new.shape
        return array_new

    def _slice(self, array, islice, axis=1, **kwargs):
        '''
        '''
        # assert 'axis' in kwargs, 'please specificy axis to operate on'
        # assert '' in kwargs, 'please specificy axis to operate on'
        # axis = kwargs['axis']
        print 'slicing array',islice
        print array.shape
        slicer = [slice(None)]*len(array.shape)
        # pdb.set_trace()
        # print slicer
        # print islice
        # print array.shape
        if axis<len(array.shape) and array.shape[axis]>0:
            slicer[axis]=islice
            array_new =  copy.copy(array[slicer])
        else:
            array_new=array
        # print slicer
        # print array[slicer]
        # array_new =  copy.copy(array[slicer])
        # array_new = np.mean(array[slicer], axis=axis)
        # print array_new.shape
        return array_new

    def _mean(self, array, **kwargs):
        '''
        '''
        mean = np.mean(array, 2)
        return mean

    def _adapt(self, array, **kwargs):
        '''
        '''
        adapt = array[:,-1]/array[:,0]
        return adapt

    def _nth(self, array, **kwargs):
        '''
        '''
        nth = array[:,kwargs['n']]
        return nth

    def _every_nth(self, array, **kwargs):
        '''
        '''
        pulses = kwargs['pulses']
        pulse_i = kwargs['pulse_i']
        mask = range(pulse_i, array.shape[1], pulses)
        every_nth = array[:,mask]
        return every_nth
#############################################################################
# clopath algorithm
#############################################################################
class Clopath():
    """ Apply Clopath learning rule, given array of synaptic input times and postsynatpic voltage traces
    
        Methods
        ---------
        _clopath: algorithm for updting synaptic weights based on input times and postsynaptic voltage
        
        Attributes
        -----------
        ~p: parameter dictionary
        ~w: weight array (synapses x samples)
        ~u_md: membrane potential for the depression
        ~u_mdbar: homeostatic membrane potential parameter
        ~u_mp: membrane potential for the potentiation
        ~u_sig: thresholded membrane potential
        ~u_md_sig: filtered and thresholded (thetam) membrane potential
        ~u_mp_sig: filtered and thresholded (thetap) membrane potential
        ~x_m0: presynaptic trace
        ~A_m: homeostatic LTD parameter
    """
    def __init__(self, ):
        '''
        '''
        pass

    def _clopath(self, x, u, fs=40, w0=0.5, homeostatic=False, param={}, lower_bound=None, upper_bound=None, **kwargs):
        """ Determine weights from voltage traces and input times with Clopath rule
        
        ===Args===
        -x      :   input array of spikes (1 for spike, 0 for no spike) (compartments x samples)
        -u      :   array of voltage time traces (compartments x samples)
        -fs     :   sampling rate (kHz)
        -w0     :   initial weight 
        -param  :   dictionary of clopath parameters
        -homeostatic : if True, homeostatic term is included in LTD equation
        
        ===Out===
        -w  :    vector of weights
        
        ===Comments===
        -see Clopath et al. 2010 Nature Neuroscience for details of learning rule
        """
        
        # preallocate parameter dictionary
        p = {}
        
        # get number of synapses (to determine shape of preallocated arrays)
        n_syn   = u.shape[0]
        
        # create time vector
        #`````````````````````````````
        # print 'data shape:', u.shape
        dt      = 1./fs
        fs      = int(fs)
        T       = u.shape[1]/fs
        dur_samples = int(T*fs)
        time    = np.arange(0., T, dt)
        
        # plasticity parameters
        #```````````````````````````````````````````````````````````````````
        p['u_ref']      = 9         # reference value for homeostatic process mV
        p['adapt_t']    = 1000*fs   # time to integrate over for homeostatic process (samples)
        p['A_m0']       = 100E-5      # A_m: amplitude for the depression
        p['A_p']        = 38E-5     # A_p: amplitude for the potentiation
        p['tau_p']      = 3         # tau_p: time constant for voltage trace in the potentiation term [ms]
        p['tau_x']      = 8       # tau_x: time constant for presynaptic trace [ms]
        p['tau_m']      = 20         # tau_m: time constant for voltage trace in the depression term [ms]
        p['tetam']      = -70       # tetam: low threshold
        p['tetap']      = -65       # tetap: high threshold in the potentiation term
        p['E_L']        = -70.6     # resting potential [mV], used for LTD adaptation
        p['delay']      = 0        # conduction delay after action potential (ms)
        p['LTD_delay']  = 1         # delay (ms) applied to presynaptic spike before computing LTD term (this allows LTD in response to subthreshold inputs) 
        # update with parameters passed as arguments
        for key, val in param.iteritems():
            p[key] = val


        # print 'clopath parameters:',p
        # preallocate learning rule variables
        #``````````````````````````````````````````````````````````````````````
        self.w           = w0*np.ones((n_syn,dur_samples))       # weights
        self.u_md        = p['E_L']*np.ones((n_syn,dur_samples)) # membrane potential for the depression
        self.u_mdbar     = np.zeros((n_syn,dur_samples))         # homeostatic membrane potential parameter
        self.u_mp        = p['E_L']*np.ones((n_syn,dur_samples)) # membrane potential for the potentiation
        self.u_sig       = p['E_L']*np.ones((n_syn,dur_samples)) # thresholded membrane potential
        self.u_md_sig    = np.zeros((n_syn,dur_samples))         # filtered and thresholded (thetam) membrane potential
        self.u_mp_sig    = np.zeros((n_syn,dur_samples))         # filtered and thresholded (thetap) membrane potential
        self.x_m0        = np.zeros((n_syn,dur_samples))         # presynaptic trace
        self.A_m         = p['A_m0']*np.ones((n_syn,dur_samples))# homeostatic LTD parameter

        # apply learning rule
        #`````````````````````````````````````````````````````````````````````` 
        # main loop
        for i, t in enumerate(time):
            
            # start simulation after specified delay
            if t>p['delay'] and t>p['LTD_delay']:
                             
                # if include homeostatic LTD mechanism
                if homeostatic:
                    
                    # adaptation voltage
                    if t <= p['adapt_t']:
                        self.u_mdbar[:,i]   = np.mean(self.u_md[:,1:i]-p['E_L'],axis=1)
                    else:
                        self.u_mdbar[:,i]   = np.mean(self.u_md[:,i-p['adapt_t']:i-1]-p['E_L'],axis=1)
                
                    # homeostatic modulation of LTD rate based on adaptation voltage
                    self.A_m[:,i]   = p['A_m0']*( self.u_mdbar[:, i-1] **2) /( p['u_ref']**2)   
    
                else:
                    # constant LTD rate
                    self.A_m[:,i]   = p['A_m0']

                # trace of membrane potential (u) with time constant tau_d
                self.u_md[:,i]  = self.u_md[ :, i-1] + dt* ( u[ :, i-1]-self.u_md[ :, i-1])/p['tau_m']
                      
                # trace of membrane potential (u) with time constant tau_p
                self.u_mp[:,i]  = self.u_mp[:, i-1]  +dt*( u[:,i-1]-self.u_mp[:, i-1]) /p['tau_p']
                 
                # trace of input spike train (spikes0)
                self.x_m0[:,i]  = self.x_m0[:,i-1]  +dt*(x[:,i-1]-self.x_m0[:,i-1]) /p['tau_x']
                
                # membrane potential (u) thresholded by thetap
                self.u_sig[:,i] = (u[:,i] > p['tetap']) *( u[:,i] -p['tetap'])
                
                # membrane potential trace (u_md) thresholded by thetam (taken 3ms before since time of a spike is 2ms)
                self.u_md_sig[:,i]  = ( self.u_md[:, i-p['delay']*fs] > p['tetam']) *( self.u_md[:, i-p['delay']*fs] -p['tetam'])                  
                
                # membrane potential trace (u_md) thresholded by thetam (taken 3ms before since time of a spike is 2ms)
                self.u_mp_sig[:,i]  = ( (self.u_mp[:,i-p['delay']*fs] -p['tetam']) > 0) *(self.u_mp[:,i-p['delay']*fs] -p['tetam'])
                
                # update weights
                self.w[:,i] = self.w[:,i-1] - self.A_m[:,i] *self.u_md_sig[:,i] *x[:,i-int(p['LTD_delay']*fs)] + p['A_p']*self.u_sig[:,i] *self.u_mp_sig[:,i] *self.x_m0[:,i]
        # apply upper and lower bounds to weights
        #-----------------------------------------
        if 'lower_bound' in p:
            lower_bound=p['lower_bound']
        if 'upper_bound' in p:
            upper_bound=p['upper_bound']
        print 'weight bounds, ', lower_bound, upper_bound
        if upper_bound is not None or lower_bound is not None:
            self.w = np.clip(self.w, a_min=lower_bound, a_max=upper_bound)

        self.p=p
        # output weight matrix (synapses x samples)
        return self.w

#############################################################################
# ebner algorithm
#############################################################################
class Ebner():
    """ Apply Clopath learning rule, given array of synaptic input times and postsynatpic voltage traces
    
        Methods
        ---------
        _clopath: algorithm for updting synaptic weights based on input times and postsynaptic voltage
        
        Attributes
        -----------
        ~p: parameter dictionary
        ~w: weight array (synapses x samples)
        ~u_md: membrane potential for the depression
        ~u_mdbar: homeostatic membrane potential parameter
        ~u_mp: membrane potential for the potentiation
        ~u_sig: thresholded membrane potential
        ~u_md_sig: filtered and thresholded (thetam) membrane potential
        ~u_mp_sig: filtered and thresholded (thetap) membrane potential
        ~x_m0: presynaptic trace
        ~A_m: homeostatic LTD parameter
    """
    def __init__(self, ):
        '''
        '''
        pass

    def _ebner(self, x, u, fs=40, w0=0.5, homeostatic=False, param={}, lower_bound=None, upper_bound=None, **kwargs):
        """ Determine weights from voltage traces and input times with Clopath rule
        
        ===Args===
        -x      :   input array of spikes (1 for spike, 0 for no spike) (compartments x samples)
        -u      :   array of voltage time traces (compartments x samples)
        -fs     :   sampling rate (kHz)
        -w0     :   initial weight 
        -param  :   dictionary of clopath parameters
        -homeostatic : if True, homeostatic term is included in LTD equation
        
        ===Out===
        -w  :    vector of weights
        
        ===Comments===
        -see Clopath et al. 2010 Nature Neuroscience for details of learning rule
        """
        
        # preallocate parameter dictionary
        p = {}
        
        # get number of synapses (to determine shape of preallocated arrays)
        n_syn   = u.shape[0]
        
        # create time vector
        #`````````````````````````````
        # print 'data shape:', u.shape
        dt      = 1./fs
        fs      = int(fs)
        T       = u.shape[1]/fs
        dur_samples = int(T*fs)
        T_samples = int(T*fs)
        time    = np.arange(0., T, dt)
        
        # plasticity parameters
        #```````````````````````````````````````````````````````````````````
        # ebner parameters from original paper
        #--------------------
        # presynaptic LTD
        #--------------------
        p['tau_T']      = 10 # (ms)
        p['theta_u_T']  = -60 # (mV)
        p['m_T']        = 1.7 # (a.u.)
        p['b_T']        = np.log(p['m_T'])/2 # (a.u.)
        # presynaptic LTP
        #--------------------
        # post component
        p['theta_u_N']  = -30 # (mV) 
        p['tau_N_alpha']  = 7.5# (ms) 
        p['m_N_alpha']  = 2 # (a.u.) 
        p['b_N_alpha']  = np.log(p['m_N_alpha'])/2 # (a.u.)
        p['tau_N_beta']  = 30 # (ms)
        p['m_N_beta']  = 10 # (a.u.)
        p['b_N_beta']  = np.log(p['m_N_beta'])/2 # (a.u.)
        p['theta_N']  = 0.2 # (a.u.)
        # pre component
        p['m_Z']  = 6 # (a.u.)
        p['b_Z']  = np.log(p['m_Z'])/2 # (a.u.)
        p['tau_Z_a'] = 1 # (ms)
        p['tau_Z_b'] = 15 # (ms)
        p['omega_Z'] = np.log(p['tau_Z_b']/p['tau_Z_a'])*p['tau_Z_a']*p['tau_Z_b']/(p['tau_Z_b']-p['tau_Z_a']) # (a.u.)
        p['epsilon_Z'] = 1/(np.exp(-p['omega_Z']/p['tau_Z_b']) - np.exp(-p['omega_Z']/p['tau_Z_a'])) # (a.u.)
        # postsynaptic LTD
        #-------------------------------------
        # pre component
        p['m_G']  = 10 # (a.u.)
        p['b_G']  = np.log(p['m_G'])/2 #(a.u.)
        p['tau_G_a'] = 2 # (ms)
        p['tau_G_b'] = 50 # (ms)
        p['omega_G'] = np.log(p['tau_G_b']/p['tau_G_a'])*p['tau_G_a']*p['tau_G_b']/(p['tau_G_b']-p['tau_G_a']) # (a.u.)
        p['epsilon_G'] = 1/(np.exp(-p['omega_G']/p['tau_G_b']) - np.exp(-p['omega_G']/p['tau_G_a'])) # (a.u.)
        # post component
        p['theta_u_C']  = -68 # (mV)
        p['theta_C_minus']  = 15 # (mV)
        p['theta_C_plus']  = 35 # (mV)

        # postsynaptic LTP
        #-------------------------------------------
        p['tau_K_beta'] = 15 # (ms)
        p['tau_K_gamma'] = 20 # (ms)
        p['m_K_beta']  = 10 # (a.u.)
        p['b_K_beta'] = np.log(p['m_K_beta'])/2 # (a.u.)
        p['s_K_beta'] = 100 # (a.u.)
        p['m_K_alpha']  = 1.5 # (a.u.)
        p['b_K_alpha'] = np.log(p['m_K_alpha'])/2 # (a.u.)

        # whole synapse parameters
        #------------------------------------------------
        p['A_pre_LTD'] = 3E-3 # 2.8E-3, # 1.5E-3 # (a.u.)
        p['A_pre_LTP'] = 33E-4 # 13E-4, # 2.5E-4 # (a.u.)
        p['A_post_LTD'] = 3.6E-4 # 3.6E-4, # 7.5E-4 # (a.u.)
        p['A_post_LTP'] = 20E-2 # 57E-2, # 7.8E-2 # (a.u.)
        p['eta'] = 1/dt # (ms^-1)
        p['w0_pre'] = w0 # (a.u.)
        p['w0_post'] = w0 # (a.u.)
        p['delay'] = 0
        p['LTD_delay']=1


        # update with parameters passed as arguments
        #--------------------------------------------
        for key, val in param.iteritems():
            p[key] = val

        # preallocate learning rule variables
        #``````````````````````````````````````````````````````````````````````
        self.D           = x # presynaptic spike trains
        self.w           = w0*np.ones((n_syn,dur_samples))  # weights
        self.w_pre           = w0*np.ones((n_syn,dur_samples))  # weights
        self.w_post           = w0*np.ones((n_syn,dur_samples))  # weights
        self.w = self.w_pre*self.w_post
        # presynaptic LTD
        #-------------------
        self.u_Tplus     = np.clip((u-p['theta_u_T']), a_min=0, a_max=None)  
        self.T_bar              = 0*np.ones((n_syn,dur_samples))     
        self.T              = 0*np.ones((n_syn,dur_samples))      
        self.E              = 0*np.ones((n_syn,dur_samples)) 
        # presynaptic LTP (postsynaptic component)
        #-----------------------------------------
        self.u_Nplus     = np.clip((u-p['theta_u_N']), a_min=0, a_max=None) 
        self.N_bar_alpha = 0*np.ones((n_syn,dur_samples))   
        self.N_alpha  = 0*np.ones((n_syn,dur_samples)) 
        self.N_bar_beta = 0*np.ones((n_syn,dur_samples))     
        self.N_beta  = 0*np.ones((n_syn,dur_samples))     
        self.N = 0*np.ones((n_syn,dur_samples))    
        # presynaptic LTP (pre component)
        #-----------------------------------------
        self.Z = 0*np.ones((n_syn,dur_samples))   
        self.Z_a = 0*np.ones((n_syn,dur_samples))    
        self.Z_b = 0*np.ones((n_syn,dur_samples))   
        self.X = 0*np.ones((n_syn,dur_samples))    
        # postsynaptic LTD (pre component)
        #-----------------------------------------------------
        self.G = 0*np.ones((n_syn,dur_samples))     
        self.G_a = 0*np.ones((n_syn,dur_samples))     
        self.G_b = 0*np.ones((n_syn,dur_samples))   
        # postsynaptic LTD (post component)
        #----------------
        self.u_Cplus     = np.clip((u-p['theta_u_C']), a_min=0, a_max=None)  
        self.C = 0*np.ones((n_syn,dur_samples))    
        self.P = 0*np.ones((n_syn,dur_samples))  
        # postsynaptic LTP
        #-----------------------------------------------------
        self.K_bar_beta = 0*np.ones((n_syn,dur_samples))     
        self.K_beta = 0*np.ones((n_syn,dur_samples))    
        self.rho = 0*np.ones((n_syn,dur_samples))   
        self.K_alpha = 0*np.ones((n_syn,dur_samples))  
        self.K_gamma = 0*np.ones((n_syn,dur_samples)) 
        self.K = 0*np.ones((n_syn,dur_samples))     


        # apply learning rule
        #`````````````````````````````````````````````````````````````````````` 
        # main loop
        for i, t in enumerate(time):
            
            # start simulation after specified delay
            if t>p['delay'] and t>p['LTD_delay']:
                # presynaptic LTD
                #-----------------------------------------------------------
                # postsynaptic component 
                self.T_bar[:,i] = self.T_bar[:,i-1] + dt*(self.u_Tplus[:,i-1]-self.T_bar[:,i-1])/p['tau_T']
                self.T[:,i] = np.tanh(p['b_T']*self.T_bar[:,i])
                # combine with presynaptic spike train
                self.E[:,i] = self.D[:,i]*self.T[:,i]
                #------------------------------------------------------------

                # presynaptic LTP
                #------------------------------------------------------------
                # post component
                self.N_bar_alpha[:,i] = self.N_bar_alpha[:,i-1] + dt*(self.u_Nplus[:,i-1]-self.N_bar_alpha[:,i-1])/p['tau_N_alpha']
                self.N_alpha[:,i] = np.tanh(p['b_N_alpha']*self.N_bar_alpha[:,i])
                self.N_bar_beta[:,i] = self.N_bar_beta[:,i-1] + dt*(self.N_bar_alpha[:,i-1]-self.N_bar_beta[:,i-1])/p['tau_N_beta']
                self.N_beta[:,i] = np.tanh(p['b_N_beta']*self.N_bar_beta[:,i])
                self.N[:,i] = np.clip(self.N_alpha[:,i]*self.N_beta[:,i]-p['theta_N'], a_min=0, a_max=None)
                # pre component
                self.Z_a[:,i] = self.Z_a[:,i-1] + dt*(p['epsilon_Z']*self.D[:,i-1])-dt*(self.Z_a[:,i-1])/p['tau_Z_a']
                self.Z_b[:,i] = self.Z_b[:,i-1] + dt*(p['epsilon_Z']*self.D[:,i-1])-dt*(self.Z_b[:,i-1])/p['tau_Z_b']
                self.Z[:,i] = np.tanh(p['b_Z']*(self.Z_b[:,i]-self.Z_a[:,i]))
                self.X[:,i] = self.Z[:,i]*self.N[:,i]
                #------------------------------------------------------------

                # postsynaptic LTD
                #----------------------------------------------------------------
                # pre component
                self.G_a[:,i] = self.G_a[:,i-1] + dt*(p['epsilon_G']*self.D[:,i-1])-dt*(self.G_a[:,i-1])/p['tau_G_a']
                self.G_b[:,i] = self.G_b[:,i-1] + dt*(p['epsilon_G']*self.D[:,i-1])-dt*(self.G_b[:,i-1])/p['tau_G_b']
                self.G[:,i] = np.tanh(p['b_G']*(self.G_b[:,i]-self.G_a[:,i]))
                # post component
                self.C[:,i] = self.G[:,i]*self.u_Cplus[:,i]
                self.P[:,i ] = np.clip(self.C[:,i]-p['theta_C_minus'],a_min=0, a_max=None)*np.clip(p['theta_C_plus']-self.C[:,i],a_min=0, a_max=None)*((p['theta_C_plus']-p['theta_C_minus'])/2)**-2
                #----------------------------------------------------------------

                # postsynaptic LTP
                #-----------------------------------------------------------------
                self.N_alpha[:,i] = np.tanh(p['b_N_alpha']*self.N_bar_alpha[:,i])
                self.K_bar_beta[:,i] = self.K_bar_beta[:,i-1] + dt*(self.K_alpha[:,i-1]-self.K_bar_beta[:,i-1])/p['tau_K_beta']
                self.K_beta[:,i] = np.tanh(p['b_K_beta']*p['s_K_beta']*self.K_bar_beta[:,i])
                self.rho[:,i] = 1-self.K_beta[:,i]
                self.K_alpha[:,i] = np.tanh(p['b_K_alpha']*np.clip(self.C[:,i]-p['theta_C_plus'], a_min=0, a_max=None))*self.rho[:,i]
                self.K_gamma[:,i] = self.K_gamma[:,i-1] + dt*(self.K_beta[:,i-1]-self.K_gamma[:,i-1])/p['tau_K_gamma']
                # print i, self.K.shape, self.K_alpha.shape, self.K_beta.shape, self.K_gamma.shape
                self.K[:,i] = self.K_alpha[:,i]*self.K_beta[:,i]*self.K_gamma[:,i]
                #-----------------------------------------------------------------

                # overall weight update
                #----------------------------------------------------------------
                self.w_pre[:,i] = self.w_pre[:,i-1] + dt*(-p['A_pre_LTD']*self.E[:,i-1]+ p['eta']*p['A_pre_LTP']*self.X[:,i-1])
                self.w_post[:,i] = self.w_post[:,i-1] + dt*(-p['eta']*p['A_post_LTD']*self.P[:,i-1]+ p['eta']*p['A_post_LTP']*self.K[:,i-1])
                #----------------------------------------------------------------

        self.w = self.w_pre*self.w_post
        # apply upper and lower bounds to weights
        #-----------------------------------------
        if 'lower_bound' in p:
            lower_bound=p['lower_bound']
        if 'upper_bound' in p:
            upper_bound=p['upper_bound']
        print 'weight bounds, ', lower_bound, upper_bound
        if upper_bound is not None or lower_bound is not None:
            self.w = np.clip(self.w, a_min=lower_bound, a_max=upper_bound)
        # store parameter dictionary
        self.p=p
        # output weight matrix (synapses x samples)
        return self.w
#############################################################################
# get spikes from voltage data
#############################################################################
class Spikes:
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        pass

    def _get_spikes(self, data, threshold,dt):
        '''
        ==Args==
        -data : array of voltage data (mV), traces x samples
        -threshold : threshold voltage for detecting spikes (mV)
        -dt : time step (ms)

        ==Out==
        -spike_idx : nested list containing indices of threshold crossings [trace_i][index]
        -spike_times : nested list containing time (ms) of threshold crossings [trace_i][time]
        -spike_train : boolean numpy array with same dimensions as original voltage data (spike onsets are 1's, all else is zero)
        
        ==Update==
        ==Comments==
        '''
        # if data is 1d array convert to 2d array
        # print data.shape
        if len(data.shape)==1:
            data = data.reshape((1,-1))

        # preallocate variables
        spike_train= np.zeros(data.shape)
        spike_idx=[]
        spike_times=[]


        # iterate over voltage traces
        for row in range(data.shape[0]):
            # print data[row,:]
            # indices of threshold crossings
            spike_i = np.where(np.diff(np.sign(data[row,:]-threshold))>0)[0]
            # time of threshold crossing
            spike_time = spike_i*dt
            # update boolean spike train array
            spike_train[row,spike_i]=1
            # add spike index and times to list
            spike_idx.append(copy.copy(spike_i))
            spike_times.append(copy.copy(spike_time))
        
        return spike_idx, spike_train, spike_times 


    def _get_xcorr_soma(self, spike_train, soma_i, ):
        '''
        '''
        group_data['soma_xcorr']['data']= np.zeros((spike_train.shape[0], spike_train.shape[1]+spike_train.shape[1]-1))
    # def _get_spike_init(self, group_data, bins)
#############################################################################
# filters
#############################################################################
class Filters:
    '''
    '''
    def __init__(self, fs):
        '''
        '''

        # 300-1000 Hz bandpass
        self.filters={}
        nqst = fs/2
        # print nyquist

        self.filters_series= {
        'iir_band_5_50':{},
        'iir_band_300_1000':{}
        }

        self.filters_series['iir_band_5_50']['high_5'] = signal.iirdesign(
            wp = 5./nqst,
            ws= 0.1/nqst,
            gpass=1,
            gstop=20,
            ftype='butter')

        self.filters_series['iir_band_5_50']['low_50'] = signal.iirdesign(
            wp = 50./nqst,
            ws= 100./nqst,
            gpass=1,
            gstop=20,
            ftype='butter')

        self.filters_series['iir_band_300_1000']['high_300'] = signal.iirdesign(
            wp = 300./nqst,
            ws= 200./nqst,
            gpass=1,
            gstop=20,
            ftype='butter')

        self.filters_series['iir_band_300_1000']['low_1000'] = signal.iirdesign(
            wp = 1000./nqst,
            ws= 1100./nqst,
            gpass=1,
            gstop=20,
            ftype='butter')

        # self.filters['iir_high_1000'] = signal.iirdesign(
        #     wp = 1000./nqst,
        #     ws = 600./nqst,
        #     gpass=1.,
        #     gstop=80.,)

        self.filters['iir_high_300'] = signal.iirdesign(
            wp = 300./nqst,
            ws = 200./nqst,
            gpass=1.,
            gstop=80.,)
        
        # self.filters['iir_highpass_1000'] = signal.butter(
        #     5,
        #     [1000./nqst],
        #     btype='high'
        #     )

        self.filters['iir_high_5'] = signal.iirdesign(
            wp = 5./nqst,
            ws= 0.1/nqst,
            gpass=1,
            gstop=20,
            ftype='butter')

    def _filter_raw_data(self, data_probe, data_induction, filters, hilbert_filters=[]):
        '''
        ==Args==
        -data_probe : voltage data during baseline probes. organized as data_probe[channel][data type]
            ~data types: 'data', 'blocks', 'comments', 'fs'
            ~data_probe['apical']['data'] = array of voltage data with dimensions samples x blocks.  blocks where there is no probe (e.g. during induction) have a placeholder column of all zeros
        -data_induction : voltage data during induction blocks. organized as data_induction[channel][data type][induction number]
            ~data types : 'data', 'blocks', 'comments', 'fs'
            ~~data_probe['apical']['data'][induction number] = 1d array of voltage data
        -filters : dictionary as filters[filter key][b coefficients, a coefficients]
        ==Out==
        -data_probe, data_induction
        ==Updates==
        -filtered data are added to data_probe and data_induction as data_probe[channel key]['data_filt_filtername'][samples x blocks]
        data_induction[channel key]['data_filt_filtername'][induction number][samples]
        -filter coefficients are stored as data_probe[channel key]['filters'][b coefficients, a coefficients]
        ==Comments==
        '''
        # iterate over channels in recorded data
        for channel_key, channel in data_probe.iteritems():
            # iterate over paths
            for path, info in channel.iteritems():
                data_probe[channel_key][path]['filters']={}
                data_induction[channel_key][path]['filters']={}
                # iterate over filters
                for filt_key, filt in filters.iteritems():
                    # apply filter to corresponding probe data
                    probe_filtered = signal.filtfilt(filt[0],filt[1],info['data'], axis=0)
                    # store filtered data
                    data_probe[channel_key][path]['data_filt_'+filt_key]=probe_filtered
                    # store filters
                    data_probe[channel_key][path]['filters'][filt_key]=filt

                    # if filt_key in hilbert_filters:
                    #     probe_filtered_hilbert = np.abs(signal.hilbert(probe_filtered, axis=0))
                    #     data_probe[channel_key][path]['data_filt_'+filt_key+'_hilbert'] = probe_filtered_hilbert

                    data_induction[channel_key][path]['filters'][filt_key]=filt
                    data_induction[channel_key][path]['data_filt_'+filt_key]=[]
                    # if filt_key in hilbert_filters:
                    #     data_induction[channel_key][path]['data_filt_'+filt_key+'_hilbert']=[]
                    # iterate over inductions
                    for induction_num, induction in enumerate(data_induction[channel_key][path]['data']):

                        # filter the current induction data
                        induction_filtered = signal.filtfilt(filt[0],filt[1],induction)
                        # store filtered data
                        data_induction[channel_key][path]['data_filt_'+filt_key].append(induction_filtered)

                        # if filt_key in hilbert_filters:
                        #     induction_filtered_hilbert = np.abs(signal.hilbert(induction_filtered, axis=0))
                        #     data_induction[channel_key][path]['data_filt_'+filt_key+'_hilbert'].append(induction_filtered_hilbert)

        return data_probe, data_induction

    def _filter_raw_data_cascade(self, data_probe, data_induction, filters_series,):
        '''
        ==Args==
        -data_probe : voltage data during baseline probes. organized as data_probe[channel][data type]
            ~data types: 'data', 'blocks', 'comments', 'fs'
            ~data_probe['apical']['data'] = array of voltage data with dimensions samples x blocks.  blocks where there is no probe (e.g. during induction) have a placeholder column of all zeros
        -data_induction : voltage data during induction blocks. organized as data_induction[channel][data type][induction number]
            ~data types : 'data', 'blocks', 'comments', 'fs'
            ~~data_probe['apical']['data'][induction number] = 1d array of voltage data
        -filters : dictionary as filters[filter key][b coefficients, a coefficients]
        ==Out==
        -data_probe, data_induction
        ==Updates==
        -filtered data are added to data_probe and data_induction as data_probe[channel key]['data_filt_filtername'][samples x blocks]
        data_induction[channel key]['data_filt_filtername'][induction number][samples]
        -filter coefficients are stored as data_probe[channel key]['filters'][b coefficients, a coefficients]
        ==Comments==
        '''
        # iterate over channels in recorded data
        for channel_key, channel in data_probe.iteritems():
            # iterate over paths
            for path, info in channel.iteritems():
                data_probe[channel_key][path]['filters_cascade']={}
                data_induction[channel_key][path]['filters_cascade']={}
                # iterate over filters
                for filt_key, filt in filters_series.iteritems():
                    data_probe[channel_key][path]['filters_cascade'][filt_key]={}
                    data_induction[channel_key][path]['filters_cascade'][filt_key]={}
                    probe_filtered = copy.deepcopy(info['data'])
                    for subfilt_key, subfilt in filt.iteritems():
                        probe_filtered = signal.filtfilt(subfilt[0], subfilt[1], probe_filtered, axis=0)
                        data_probe[channel_key][path]['filters_cascade'][filt_key][subfilt_key]=subfilt
                        data_induction[channel_key][path]['filters_cascade'][filt_key][subfilt_key]=subfilt

                    # store filtered data
                    data_probe[channel_key][path]['data_filt_'+filt_key]=probe_filtered
                    # hilbert
                    #--------------
                    # if filt_key in hilbert_filters:
                    #     probe_filtered_hilbert = np.abs(signal.hilbert(probe_filtered, axis=0))
                    #     data_probe[channel_key][path]['data_filt_'+filt_key+'_hilbert'] = probe_filtered_hilbert

                    data_induction[channel_key][path]['data_filt_'+filt_key]=[]
                    # hilbert
                    # if filt_key in hilbert_filters:
                    #     data_induction[channel_key][path]['data_filt_'+filt_key+'_hilbert']=[]
                    # iterate over inductions
                    for induction_num, induction in enumerate(data_induction[channel_key][path]['data']):

                        induction_filtered = copy.deepcopy(induction)

                        for subfilt_key, subfilt in filt.iteritems():
                            # filter the current induction data
                            induction_filtered = signal.filtfilt(subfilt[0],subfilt[1],induction_filtered)
                        # store filtered data
                        data_induction[channel_key][path]['data_filt_'+filt_key].append(induction_filtered)
                        # hilbert
                        # if filt_key in hilbert_filters:
                        #     induction_filtered_hilbert = np.abs(signal.hilbert(induction_filtered, axis=0))
                        #     data_induction[channel_key][path]['data_filt_'+filt_key+'_hilbert'].append(induction_filtered_hilbert)

        return data_probe, data_induction

        # load i

    def _get_hilbert(self, data_probe, data_induction, hilbert_filters):
        '''
        '''
        # iterate over channels in recorded data
        for channel_key, channel in data_probe.iteritems():
            # iterate over paths
            for path, info in channel.iteritems():
                # iterate over filters
                for filt_key in hilbert_filters:
                    # apply filter to corresponding probe data
                    probe_hilbert = np.abs(signal.hilbert(info['data_filt_'+filt_key] - np.mean(info['data_filt_'+filt_key], axis=0), axis=0))
                    # store filtered data
                    data_probe[channel_key][path]['data_filt_'+filt_key+'_hilbert']=probe_hilbert
                    data_induction[channel_key][path]['data_filt_'+filt_key+'_hilbert']=[]
                    # iterate over inductions
                    for induction_num, induction in enumerate(data_induction[channel_key][path]['data_filt_'+filt_key]):

                        # filter the current induction data
                        induction_hilbert = np.abs(signal.hilbert(induction-np.mean(induction)))
                        # store filtered data
                        data_induction[channel_key][path]['data_filt_'+filt_key+'_hilbert'].append(induction_hilbert)
        return data_probe, data_induction
#############################################################################
# shapeplot functions
#############################################################################
class ShapePlot():
    """ create plot of neuron morphology with data as a colormap

    similar to ShapePlot in neuron
    """
    
    def __init__(self):
        """
        """
        pass

    #__________________________________________________________________________
    def basic(self, morpho, data, axes, width_scale=1, colormap=colormap.jet):
        """ creates a shapeplot of a given neuron morphology with the given data

        Arguments:
        morpho: morphology structure with dimensions {tree}[section][segment](index, name, x, y, z, diam, parent_index)

        data: single floating point values store in a structure with dimensions {tree}[section][segment]

        axes: pyplot axes object

        colormap: matplotlib.cm colormap
        """
        # create list of points
        if type(morpho)==dict:
            morph_points=[]
            indexes = []
            data_values=[]
            for tree_key, tree in morpho.iteritems():
                for sec_i, sec in enumerate(tree):
                    for seg_i, seg in enumerate(sec):
                        morph_points.append(seg)
                        indexes.append(copy.copy(seg[0]))
                        data_values.append(data[tree_key][sec_i][seg_i])

        else:
            morph_points=morpho
            indexes = [val[0] for val in morph_points]
            data_values = data
        
        # resort lists
        sort_list = np.argsort(indexes)
        morph_points_sorted = [morph_points[i] for i in sort_list]
        data_values_sorted = [data_values[i] for i in sort_list]

        # lists to store segment shapes and their data values (colors)
        patches=[]
        colors=[]

        # iterate through child sections
        for child_i, child in enumerate(morph_points_sorted):
            # if segment is the root segment (parent index is -1)
            if child[-1]==-1:
                # skip, the root will be added as a parent to its children 
                continue
            # if not a root segment
            else:
                # find parent segment
                parent_idx = child[-1]
                parent = [val for i,val in enumerate(morph_points_sorted) if val[0]==parent_idx][0]
            
            # interpolate xyz values between parent and child
            parent_point = (parent[2], parent[3], parent[4])
            child_point = (child[2],child[3],child[4])
            mid_point = self.interpolate_3d(point1=parent_point, point2=child_point, t=0.5)

            # get diameter of parent and child 
            parent_diam = width_scale*parent[5]
            child_diam = width_scale*child[5]
            mid_diam = (parent_diam+child_diam)/2.

            # get data values for parent and child
            parent_color = data_values_sorted[parent[0]]
            child_color = data_values_sorted[child_i]

            # create polygon patch to plot segment
            parent_polygon = self.make_polygon(point1=parent_point, point2=mid_point, d1=parent_diam/2., d2=mid_diam/2.)
            child_polygon =self.make_polygon(point1=child_point, point2=mid_point, d1=child_diam/2., d2=mid_diam/2.)

            # add to list of patches
            patches.append(parent_polygon)
            colors.append(parent_color)
            patches.append(child_polygon)
            colors.append(child_color)
        # print patches

        # create patch collection
        p = PatchCollection(patches, cmap=colormap, alpha=1.)
        # set colors
        p.set_array(np.array(colors))
        # plot collection
        # axes.add_collection(p)
        # # show colorbar
        # plt.colorbar(p)
        # # autoscale axes
        # axes.autoscale()
        # plt.show(block=False)
        return p, colors


    # FIXME show synapse locations on shapeplot
    def show_synapses(self, morpho, data, axes, width_scale=1, colormap=colormap.jet):
        """ creates a shapeplot of a given neuron morphology with the given data

        Arguments:
        morpho: morphology structure with dimensions {tree}[section][segment](index, name, x, y, z, diam, parent_index)

        data: single floating point values store in a structure with dimensions {tree}[section][segment]

        axes: pyplot axes object

        colormap: matplotlib.cm colormap
        """
        # create list of points
        if type(morpho)==dict:
            morph_points=[]
            indexes = []
            data_values=[]
            for tree_key, tree in morpho.iteritems():
                for sec_i, sec in enumerate(tree):
                    for seg_i, seg in enumerate(sec):
                        morph_points.append(seg)
                        indexes.append(copy.copy(seg[0]))
                        data_values.append(data[tree_key][sec_i][seg_i])

        else:
            morph_points=morpho
            indexes = [val[0] for val in morph_points]
            data_values = data
        
        # resort lists
        sort_list = np.argsort(indexes)
        morph_points_sorted = [morph_points[i] for i in sort_list]
        data_values_sorted = [data_values[i] for i in sort_list]

        # lists to store segment shapes and their data values (colors)
        patches=[]
        colors=[]

        # iterate through child sections
        for child_i, child in enumerate(morph_points_sorted):
            # if segment is the root segment (parent index is -1)
            if child[-1]==-1:
                # skip, the root will be added as a parent to its children 
                continue
            # if not a root segment
            else:
                # find parent segment
                parent_idx = child[-1]
                parent = [val for i,val in enumerate(morph_points_sorted) if val[0]==parent_idx][0]
            
            # interpolate xyz values between parent and child
            parent_point = (parent[2], parent[3], parent[4])
            child_point = (child[2],child[3],child[4])
            mid_point = self.interpolate_3d(point1=parent_point, point2=child_point, t=0.5)

            # get diameter of parent and child 
            parent_diam = width_scale*parent[5]
            child_diam = width_scale*child[5]
            mid_diam = (parent_diam+child_diam)/2.

            # get data values for parent and child
            parent_color = data_values_sorted[parent[0]]
            child_color = data_values_sorted[child_i]

            # create polygon patch to plot segment
            parent_polygon = self.make_polygon(point1=parent_point, point2=mid_point, d1=parent_diam/2., d2=mid_diam/2.)
            child_polygon =self.make_polygon(point1=child_point, point2=mid_point, d1=child_diam/2., d2=mid_diam/2.)

            # add to list of patches
            patches.append(parent_polygon)
            colors.append(parent_color)
            patches.append(child_polygon)
            colors.append(child_color)
        # print patches

        # create patch collection
        p = PatchCollection(patches, cmap=colormap, alpha=1.)
        # set colors
        p.set_array(np.array(colors))
        # plot collection
        # axes.add_collection(p)
        # # show colorbar
        # plt.colorbar(p)
        # # autoscale axes
        # axes.autoscale()
        # plt.show(block=False)
        return p, colors

    # create 3d interpolation function
    def interpolate_3d(self, point1, point2, t):
        """
        Arguments:
        point1 and point2 are lists/tuples of 3d points as [x,y,z]

        t is the parameter specifying relative distance between the two points from 0 to 1

        returns a tuple with the 3d coordinates of the requested point along the line (x,y,z)
        """
        x = point1[0] + point2[0]*t  - point1[0]*t
        y = point1[1] + point2[1]*t  - point1[1]*t
        z = point1[2] + point2[2]*t  - point1[2]*t

        return (x,y,z)

    def make_polygon(self, point1, point2, d1, d2):
        """ make a matplotlib.patches.Polygon object for plotting

        Arguments:
        point1, point2: list/tuple containing 3d coordinates of points to be connected

        d1,d2: diameter of each point 

        the function determines the line between the two points, then adds the corresponding diameters to each point along the orthogonal direction to the line, this produces the four points that the polygon

        returns a matplotlib.patches.Polygon object for plotting
        """
        # find slope of connecting line
        dy = point2[1]-point1[1]
        dx = point2[0] - point1[0]
        m = dy/dx

        # find slope of orthogonal line
        m_inv = -1./m

        # find x and y changes to add diameter
        delta_y1 = m_inv*np.sqrt(d1**2/(m_inv**2+1))
        delta_x1  = np.sqrt(d1**2/(m_inv**2+1))
        delta_y2 = m_inv*np.sqrt(d2**2/(m_inv**2+1))
        delta_x2  = np.sqrt(d2**2/(m_inv**2+1))

        # add diameter to first point
        y1_pos = point1[1] + delta_y1
        y1_neg = point1[1] - delta_y1
        x1_pos = point1[0] + delta_x1
        x1_neg = point1[0] - delta_x1

        # add diameter to second point
        y2_pos = point2[1] + delta_y2
        y2_neg = point2[1] - delta_y2
        x2_pos = point2[0] + delta_x2
        x2_neg = point2[0] - delta_x2

        # store points as a 4 x 2 array
        points = np.array([[x1_pos, y1_pos],[x1_neg, y1_neg],[x2_neg, y2_neg],[x2_pos, y2_pos]])

        return Polygon(points)
#######################################################################################################################################################################################################################################
# DEPRECATED
#######################################################################################################################################################################################################################################
class BuildFigDF:
    '''
    '''
    def __init__(self, ):
        '''
        '''
        pass

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

class Stats:
    '''
    '''
    def _pairwise_ttests(self, df_sorted, variable, array_funcs=[], array_func_kws=[]):
        '''
        '''
        # stats

        data={}
        for tracekey in df_sorted.keys():
            # get series from df sorted
            trace_series = df_sorted[tracekey][variable]

            # convert to array
            data[tracekey] = functions._2array(trace_series, remove_nans=True, remove_nans_axis=1, array_funcs=array_funcs, array_func_kws=array_func_kws)#*1000

            # FIXME how to apply array functions here
            if type(data[tracekey])==np.ndarray and data[tracekey].shape[0]>0:
                print tracekey, data[tracekey].shape
                for i, array_func in enumerate(array_funcs):
                    data[tracekey]= array_func(data[tracekey], **array_func_kws[i])

        # get pairwise combinations of traces
        combos = itertools.combinations(data.keys(), 2)
        ttests = {}
        for combo in combos:
            # print combo[0], combo[1]
            # print data[combo[0]]
            # print data[combo[1]]
            ttests[combo] = stats.ttest_ind(data[combo[0]], data[combo[1]])

        return ttests


    def _linregress(self, df_sorted, x_variable, y_variable, array_funcs_x=[], array_func_kws_x=[], array_funcs_y=[], array_func_kws_y=[]):
        '''
        '''
        regressions={}
        data_x={}
        data_y={}
        for tracekey in df_sorted:
            trace_series_x = df_sorted[tracekey][x_variable]
            trace_series_y = df_sorted[tracekey][y_variable]

            data_x[tracekey] = functions._2array(trace_series_x, remove_nans=True, remove_nans_axis=1,)
            data_y[tracekey] = functions._2array(trace_series_y, remove_nans=True, remove_nans_axis=1,)

            # apply array functions
            for i, array_func in enumerate(array_funcs_x):
                data_x[tracekey]= array_func(data_x[tracekey], **array_func_kws[i])
            for i, array_func in enumerate(array_funcs_y):
                data_y[tracekey]= array_func(data_y[tracekey], **array_func_kws[i])

            slope, intercept, r_value, p_value, std_err = stats.linregress(x=data_x[tracekey], y=data_y[tracekey])

            regressions[tracekey]['x_variable'] = x_variable
            regressions[tracekey]['y_variable'] = y_variable
            regressions[tracekey]['slope'] = slope
            regressions[tracekey]['intercept'] = intercept
            regressions[tracekey]['r_value'] = r_value
            regressions[tracekey]['p_value'] = p_value
            regressions[tracekey]['std_err'] = std_err
            
        return regressions

    def _cca(self, df_sorted, x_variables, y_variables):
        '''
        '''
        x_data={}
        y_data={}
        cca={}
        transformed={}
        for tracekey in df_sorted.keys():
            x_data[tracekey] = []
            for x_i, x_var in enumerate(x_variables):

                # get series from df sorted
                trace_series = df_sorted[tracekey][x_var]

                x_data[tracekey].append(functions._2array(trace_series, remove_nans=True, remove_nans_axis=1))

            x_data[tracekey] = np.array(x_data[tracekey]).squeeze().T
            y_data[tracekey] = []
            for y_i, y_var in enumerate(y_variables):

                # get series from df sorted
                trace_series = df_sorted[tracekey][x_var]
                y_data[tracekey].append(functions._2array(trace_series, remove_nans=True, remove_nans_axis=1))
            y_data[tracekey] = np.array(y_data[tracekey]).squeeze().T
            if x_data[tracekey].shape[0]>0:
                cca[tracekey] = CCA()
                cca[tracekey].fit(x_data[tracekey], y_data[tracekey])
                x_c, y_c = cca[tracekey].transform(x_data[tracekey], y_data[tracekey])
                transformed[tracekey]=[x_c, y_c]

        return cca, transformed

    def _anova_ind(self, df_sorted, variable):
        '''
        '''
        data={}
        for tracekey in df_sorted.keys():
            # get series from df sorted
            trace_series = df_sorted[tracekey][variable]

            # convert to array
            data[tracekey] = functions._2array(trace_series, remove_nans=True, remove_nans_axis=1)#*1000

        # get pairwise combinations of traces

        combos = itertools.combinations(data.keys(), 2)
        ttests = {}
        for combo in combos:
            ttests[combo] = stats.ttest_ind(data[combo[0]], data[combo[1]])

class IndexGetter:
    ''' class for slicing data array based on specified conditions

    ==Methods==
    -__init__ : creates variable called conditions to be passed to other methods
                -conditions should be a list of tuples [(condition type, [list of specific conditions])]
                        -e.g. conditions = [('polarity',['anodal','cathodal']),
                        ('syn_num',[[6,8],[10,12]])]


    '''
    def __init__(self, **kwargs):
        '''
        '''
        if 'conditions' in kwargs:
            conditions=kwargs['conditions']
        else:
            conditions=[('polarity', ['anodal','cathodal','control']), ('syn_num',[[6,8,10]]), ('dist', [[[0,200],[0,300]]]), ('path',[['1']])]
        if 'group_data' in kwargs and 'variable' in kwargs:
            self.idx_sets, self.combos = self._get_idx(group_data=kwargs['group_data'], variable=kwargs['variable'], conditions=conditions)

    def _match_idx(self, group_data, variable1, variable2, variable1_idx_sets):
        '''
        '''
        # FIXME
        variable2_idx_sets=[]
        for idx_set_i, idx_set in enumerate(variable1_idx_sets):
            locations1 = [group_data[variable1]['locations'][temp] for temp in idx_set]
            fields1 = [group_data[variable1]['field'][temp] for temp in idx_set]
            trial_ids1 = [group_data[variable1]['trial_id'][temp] for temp in idx_set]



            idx_set2=[]
            for loc_i, loc in enumerate(locations1):
                location2_idx = [ temp_i for temp_i, temp in enumerate(group_data[variable2]['locations']) if temp==loc]
                field2_idx = [ temp_i for temp_i, temp in enumerate(group_data[variable2]['field']) if temp==group_data[variable1]['field'][loc_i]]
                trial_id2_idx = [ temp_i for temp_i, temp in enumerate(group_data[variable2]['trial_id']) if temp==group_data[variable1]['trial_id'][loc_i]]
                idx2 = sorted(list(set.intersection(set(location2_idx), set(field2_idx), set(trial_id2_idx))))[0]
                idx_set2.append(idx2)

            variable2_idx_sets.append(idx_set2)


        return variable2_idx_sets

    def _get_idx(self, group_data, variable, conditions):
        ''' get indices of data entries with specified conditions

        ==Args==
        -group_data : group data structure organize as group_data[variable][info/condition type][recorded segment]
        -conditions : list of tuples to specifify desired conditions [(condition type, [list of specific conditions])]
                        -e.g. conditions = [('polarity',['anodal','cathodal']),
                        ('syn_num',[[6,8],[10,12]])]
                        -all possible combinations of the specified conditions will be returned in idx sets

        ==Out==
        -idx_sets : nested list of indices. first dimension contains an entry for each combination of the conditions specified [condition combination][indices in group_data]
        -combos :  list of conditions for each combination in idx sets (same organization as conditions)
        '''
        # print 'getting indices'
        condition_types = zip(*conditions)[0]
        condition_lists = zip(*conditions)[1]
        condition_combos = [temp for temp in itertools.product(*condition_lists)]
        att_list=dir(self)
        idx={}
        idx_sets=[]
        combos=[]
        for combo_i, combo in enumerate(condition_combos):

            for condition_i, condition in enumerate(combo):
                
                condition_type = condition_types[condition_i]

                att = [temp for temp_i, temp in enumerate(att_list) if condition_type in temp][0]
                
                type_func = getattr(self, att)
                # print type_func

                idx[condition_type] = type_func(group_data=group_data, variable=variable, conditions=condition)
                # if len(idx[condition_type])!=0:
                    # print 'idx found:',condition_type, condition


            condition_sets=[]
            for condition_type_key in idx.keys():
                condition_set = set(idx[condition_type_key])
                condition_sets.append(condition_set)

            # print condition_sets
            idx_sets.append(sorted(list(set.intersection(*condition_sets)))) 
            # print idx_sets
            combo_rezip = zip(condition_types, combo)
            combos.append(combo_rezip)

        return idx_sets, combos

    def _get_dist(self, group_data, variable, conditions=[]):
        '''
        '''

        if len(conditions)==0:

            idx=  [temp_i for temp_i, temp in enumerate(group_data[variable]['path_name']) if temp[0]!='None']
        
        else:# distance info
            dist=[]
            for path_i, path in enumerate(group_data[variable]['path_name']):
                if path[0]=='None':
                    dist.append('None')

                else:
                    dist.append(group_data[variable]['p'][path_i]['p_path'][path[0]]['syn_dist'])


            idx= [temp_i for temp_i, temp in enumerate(dist) if temp in conditions]

        return idx

    def _get_syn_num(self, group_data, variable, conditions=[], **kwargs):
        '''
        '''
        if len(conditions)==0:
            idx = [temp_i for temp_i, temp in enumerate(group_data[variable]['p'])]
        else:
            syn_num = []
            for p_i, p in enumerate(group_data[variable]['p']):
                path_syn_num=[]
                for path_key, path in p['p_path'].iteritems():
                    path_syn_num.append(path['syn_num'])
                syn_num.append(path_syn_num)

            # syn_num = [temp['p_path']['1']['syn_num']for temp_i, temp in enumerate(group_data[variable]['p'])]
            idx = [temp_i for temp_i, temp in enumerate(syn_num) if bool(set(temp)&set(conditions))]
        return idx

    def _get_path(self, group_data, variable, conditions=[]):
        '''
        '''
        if len(conditions)==0:
            idx = [temp_i for temp_i, temp in enumerate(group_data[variable]['path_name'])]

        else:
            idx = [temp_i for temp_i, temp in enumerate(group_data[variable]['path_name']) if bool(set(temp)&set(conditions))]

        return idx

    def _get_polarity(self, group_data, variable, conditions):
        '''
        '''
        idx = [temp_i for temp_i, temp in enumerate(group_data[variable]['polarity']) if temp == conditions]
        return idx

    def _get_magnitude(self, group_data, variable, conditions):
        '''
        '''
        idx = [temp_i for temp_i, temp in enumerate(group_data[variable]['field']) if temp in conditions]
        return idx

    def _get_stdp_dt(self, group_data, variable, conditions):
        '''
        '''
        print conditions
        print group_data[variable]['stdp_dt']
        idx = [temp_i for temp_i, temp in enumerate(group_data[variable]['stdp_dt']) if bool(set(temp)& set(conditions))]

        return idx

    def _get_pulse_freq(self, group_data, variable, conditions):
        '''
        '''
        idx = [[temp_i for temp_i, temp in enumerate(group_data[variable]['pulse_freq']) if temp in conditions]]

    def _get_tree(self, group_data, variable, conditions,):
        '''
        '''
        # print 'getting tree'
        idx = [temp_i for temp_i, temp in enumerate(group_data[variable]['locations']) if bool(set(temp)& set(conditions))]
        return idx
    
    def _find_soma(self, group_data, variable):
        '''
        '''
        soma_idx=[]
        for location_i, location in enumerate(group_data[variable]['locations']):
            # get trial_id
            trial_id = group_data[variable]['trial_id'][location_i]
            field = group_data[variable]['field'][location_i]
            trial_id_idx = [temp_i for temp_i, temp  in enumerate(group_data[variable]['trial_id']) if temp==trial_id]
            field_idx = [temp_i for temp_i, temp  in enumerate(group_data[variable]['field']) if temp==field]

            soma_i = [temp_i for temp_i, temp  in enumerate(group_data[variable]['locations']) if temp==('soma',0,0)]

            soma_idx.append(list(set.intersection(set(trial_id_idx), set(field_idx), set(soma_i)))[0])

        return soma_idx

class Plots():

    def __init__(self):
        '''
        '''
        self.x=1
        pass

    def _range_variable_single_trial(self, data, y_variable, x_variable, **kwargs):
        '''
        '''
        p = data[y_variable]['p']
        t = data[y_variable]['t']
        
        # unique list of recorded locations ()
        locations_unique = list(set(data[y_variable]['locations']))

        nseg = len(locations_unique)
        cols = int(math.ceil(math.sqrt(nseg)))
        print cols
        rows = int(math.ceil(math.sqrt(nseg)))
        fig, ax = plt.subplots(rows, cols)
        # if only one subplot, cast as array
        if type(ax) is not np.ndarray:
            ax = np.array(ax).reshape(1,1)

        print ax
        # iterate over unique locations in data
        row, col = 0,-1

        for loc_i, loc in enumerate(locations_unique):
            if col<cols-1:
                col+=1
            else:
                col=0
                row+=1
            print row, col
            tree, sec_num, seg_num = loc
            
            # get list of indices that match the current
            locs_i = [seg_i for seg_i , seg in enumerate(data[y_variable]['locations']) if seg == locations_unique[loc_i]] 
            print 'locs_i',locs_i
            print data[y_variable]['field']
            data_y = data[y_variable]['data'][locs_i,:]

            # get x_data for the current location
            if x_variable=='t':
                data_x=np.tile(t, (len(locs_i),1))
            else:
                #FIXME x variable indices may be different
                data_x = data[x_variable]['data'][locs_i,:]

            # get list of field magnitudes that match the data for the current location
            fields = [data[y_variable]['field'][i] for i in locs_i ]
            print fields
            
            # get list of plot colors based on field magnitude
            field_colors = [p['field_colors'][p['field'].index(field)] for field_i, field in enumerate(fields)]

            # distance from soma
            seg_dist = p['seg_dist'][tree][sec_num][seg_num]

            # title segment location and distance from soma
            title = loc[0] + ',' + str(loc[1]) + ',' + str(loc[2]) + ',' + ('%.2f'%seg_dist)

            for field_i in range(data_y.shape[0]):
            
                ax[row,col].plot(data_x[field_i, :], data_y[field_i,:], color=field_colors[field_i], linewidth=0.5)

            ax[row,col].set_title(title)
            ax[row,col].set_ylabel(y_variable)
            ax[row,col].set_xlabel(x_variable)    

        plt.show()

        return fig, ax

    def _range_variable_group(self, group_data, x_variable, y_variable, conditions=[('polarity',['anodal','cathodal','control']),('dist',[[[0,200],[0,300]]]), ('syn_num',[[6]]),('path',[['1']])], plot_mean=True, plot_sem=True):
        '''
        '''
        ig_y = IndexGetter(group_data=group_data, variable=y_variable, conditions=conditions)
        self.idx_sets_y=ig_y.idx_sets
        self.combos_y=ig_y.combos

        if x_variable=='t':
            self.idx_sets_x=self.idx_sets_y
            self.combos_y=self.combos_y
        else:
            self.idx_sets_x = ig_y._match_idx(group_data=group_data, variable1=y_variable, variable2=x_variable, variable1_idx_sets=self.idx_sets_y)
            self.combos_x=self.combos_y


        polarities=['cathodal','anodal','control']
        colors=['b','r','k']
        std_colors = ['lightblue','coral','grey']
        fig = plt.figure()
        for idx_set_i, idx_set in enumerate(self.idx_sets_y):
            # print idx_set
            polarity = [cond[1] for cond_i, cond in enumerate(self.combos_y[idx_set_i]) if 'polarity' in cond][0]
            color = colors[polarities.index(polarity)]
            std_color = std_colors[polarities.index(polarity)]
            
            if x_variable=='t':
                x_data = group_data[y_variable]['t']
                x_data_mean=x_data
            else:
                x_data = group_data[x_variable]['data'][self.idx_sets_x[idx_set_i],:]
                x_data_mean=np.mean(x_data, axis=0)
                x_data_std=np.std(x_data, axis=0)
                x_data_sem=stats.sem(x_data,axis=0)

            y_data = group_data[y_variable]['data'][idx_set,:]
            y_data_mean=np.mean(y_data,axis=0)
            y_data_std=np.std(y_data,axis=0)
            y_data_sem=stats.sem(y_data,axis=0)

            if plot_mean and plot_sem:

                plt.plot(x_data_mean, y_data_mean,color=color)
                plt.fill_between(x_data_mean, y_data_mean-y_data_sem, y_data_mean+y_data_sem, color=color, alpha=0.5)

            elif plot_mean:
                plt.plot(x_data_mean, y_data_mean,color=color)

            else:
                plt.plot(x_data.transpose(),y_data.transpose(), color=color)
        # plt.show(block=False)

        return fig

    def _stdp_frequency_weight(self, group_data):
        '''
        '''
        # iterate over stdp delays
        # iterate over stimulation frequencies
            # get corresponding indices in group data
            # get net weight change for each synapse
            # get average and std for current frequency and store in array
            # plot frequency vs average weight change

        stdp_delays = list(set(group_data['clopath']['stdp_dt']))
        stim_freqs = list(set([item for sublist in group_data['clopath']['pulse_freq'] for item in sublist]))

        dw = {}
        markers = ['.','x']
        plt.figure()
        for delay_i, delay  in enumerate(stdp_delays):
            dw[str(delay)]={'freq':[],'dw':[]}
            for stim_freq in stim_freqs:
                dw[str(delay)]['freq'].append(stim_freq)
                conditions = [('path',[['1']]),('stdp_dt',[[delay]]), ('pulse_freq',[[stim_freq]])]
                ig = IndexGetter(group_data=group_data, variable='clopath', conditions=conditions)
                self.idx_sets=ig.idx_sets
                self.combos=ig.combos
                dw[str(delay)]['dw'].append(group_data['clopath']['data'][idx_sets[0],-1]/group_data['clopath']['data'][idx_sets[0],0])

                plt.plot(stim_freq, np.mean(dw[str(delay)]['dw'][-1]), marker=markers[delay_i])
        plt.show(block=False)

    def _dose_response_weight(self, group_data):
        '''
        '''
        # get all intensities
        intensities = np.array(sorted(list(set(group_data['clopath']['field']))))
        x_locations=np.arange(intensities.shape[0])
        dw_mean=np.zeros(intensities.shape[0])
        dw_sem=np.zeros(intensities.shape[0])
        colors=[]
        for intensity_i, intensity in enumerate(intensities):
            intensity_idx = [temp_i for temp_i, temp in enumerate(group_data['clopath']['field']) if temp==intensity]
            dw_mean[intensity_i] = np.mean(group_data['clopath']['data'][intensity_idx,-1]/group_data['clopath']['data'][intensity_idx,0])
            dw_sem[intensity_i] = stats.sem(group_data['clopath']['data'][intensity_idx,-1]/group_data['clopath']['data'][intensity_idx,0])
            if intensity<0.:
                colors.append('blue')
            elif intensity==0.:
                colors.append('black')
            else:
                colors.append('red')

        fig, ax = plt.subplots()
        plt.bar(x_locations, dw_mean, yerr=dw_sem, tick_label=intensities, color=colors)
        plt.ylim((min(dw_mean-dw_sem),None))
        plt.xlabel('Field magnitude (V/m)', fontsize=25, fontweight='heavy')
        plt.ylabel('Norm. weight', fontsize=25, fontweight='heavy')
        for temp in ax.get_xticklabels():
            temp.set_fontweight('heavy')
            temp.set_fontsize(12)
        for temp in ax.get_yticklabels():
            temp.set_fontweight('heavy')
            temp.set_fontsize(12)
        # set axes linewidth and tick position
        #----------------------
        # ax.spines['left'].set_linewidth(10)
        # ax.spines['bottom'].set_linewidth(10)
        # ax.xaxis.set_ticks_position('bottom')
        # ax.yaxis.set_ticks_position('left')
        plt.show(block=False)

        return fig

    def _weight_vs_distance(self, group_data, conditions=[('polarity', ['anodal','cathodal','control']), ('syn_num',[[8,]]), ('dist', [[[0,200],]]), ('path',[['1']])]):
        '''
        '''

        ig = IndexGetter(group_data=group_data, conditions=conditions)
        self.idx_sets=ig.idx_sets
        self.combos=ig.combos

        polarities=['cathodal','anodal','control']
        colors=['b','r','k']
        plt.figure()
        for idx_set_i, idx_set in enumerate(self.idx_sets):
            polarity = [cond[1] for cond_i, cond in enumerate(self.combos[idx_set_i]) if 'polarity' in cond][0]
            color = colors[polarities.index(polarity)]
            plt.plot(group_data['clopath']['t'], np.mean(group_data['clopath']['data'][idx_set,:],axis=0).transpose(),color=color)

        plt.xlabel('time (ms)')
        plt.ylabel('weight (AU)')
        # plt.show(block=False)

    def _final_weight_vs_distance(self, group_data, conditions):
        '''
        '''
        ig = IndexGetter(group_data=group_data, conditions=conditions, variable='clopath')
        self.idx_sets=ig.idx_sets
        self.combos=ig.combos

        polarities=['cathodal','anodal','control']
        colors=['b','r','k']
        plt.figure()
        for idx_set_i, idx_set in enumerate(self.idx_sets):
            polarity = [cond[1] for cond_i, cond in enumerate(self.combos[idx_set_i]) if 'polarity' in cond][0]
            color = colors[polarities.index(polarity)]

    def _xcorr_mean(self, group_data, conditions, dt=.025):
        '''
        '''
        ig = IndexGetter(group_data=group_data, conditions=conditions, variable='soma_xcorr')
        self.idx_sets=ig.idx_sets
        self.combos=ig.combos

        # print self.combos
        polarities=['cathodal','anodal','control']
        colors=['b','r','k']
        fig = plt.figure()
        for idx_set_i, idx_set in enumerate(self.idx_sets):
            polarity = [cond[1] for cond_i, cond in enumerate(self.combos[idx_set_i]) if 'polarity' in cond][0]
            color = colors[polarities.index(polarity)]
            n_samples = group_data['soma_xcorr']['data'].shape[1]
            n_trials=float(len(list(set([group_data['soma_xcorr']['trial_id'][temp] for temp in idx_set]))))
            soma_xcorr_norm = group_data['soma_xcorr']['data'][idx_set,:]/n_trials
            x_vector = dt*(np.arange(n_samples)-(n_samples+1)/2) 
            plt.plot(x_vector, np.mean(soma_xcorr_norm,axis=0).transpose(),color=color)
            plt.xlabel('delay (negative=dendrite first)')
            plt.ylabel('probability density / correlation (spikes^2/trial)')

        return fig

class GroupData:
    '''
    ==Att==
    -group_data :  dictionary for holding data for all simulations in a set of experiments; group_data{variable}{info type}[segment]

    '''

    def __init__(self, directory, file_name, **kwargs):
        '''
        '''
        self.directory=directory
        self.file_name=file_name
        if 'group_data' in kwargs:
            self.group_data=kwargs['group_data']
        else:
            self.group_data=self._load_group_data(directory=directory, file_name=file_name)

        if 'rerun' in kwargs:
            self.rerun=kwargs['rerun']
        else:
            self.rerun=False

    def _standard_run(self, group_data, directory, file_name, search_string, clopath_param, variables, ):
        '''
        '''
        self._add_files_to_group_data(group_data=group_data, directory=directory, search_string=search_string,variables=variables)
        
        self._add_path_to_group_data(group_data=group_data)
        self._add_input_array_to_group_data(group_data=group_data)
        self._add_clopath_to_group_data(group_data=group_data, param=clopath_param, rerun = self.rerun)
        self._add_polarity_to_group_data(group_data=group_data)
        self._add_spikes_to_group_data(group_data=group_data, threshold=-30)
            # self._add_xcorr_to_group_data(group_data=group_data)
        return self.group_data
        # self._save_group_data(group_data=group_data, directory=directory, file_name=file_name)

    def _conditional_run(self, group_data, directory, file_name, search_string, clopath_param, variables, parameter, parameter_value, path_parameter):
        '''
        '''
        self._add_files_to_group_data_conditional_parameter(group_data=group_data, directory=directory, search_string=search_string,variables=variables, parameter=parameter, parameter_value=parameter_value, path_parameter=path_parameter)
        if len(group_data['processed'])>0:
            self._add_path_to_group_data(group_data=group_data)
            self._add_input_array_to_group_data(group_data=group_data)
            self._add_clopath_to_group_data(group_data=group_data, param=clopath_param, rerun = self.rerun)
            self._add_polarity_to_group_data(group_data=group_data)
            self._add_spikes_to_group_data(group_data=group_data, threshold=-30)
        # self._add_xcorr_to_group_data(group_data=group_data)
        return self.group_data
        # self._save_group_data(group_data=group_data, directory=directory, file_name=file_name)
    
    def _load_group_data(self, directory='', file_name=''):
        """ Load group data from folder
        
        ===Args===
        -directory : directory where group data is stored including /
        -file_name : file name for group data file, including .pkl
                    -file_name cannot contain the string 'data', since this string is used t search for individual data files

        ===Out===
        -group_data  : typically a dictionary.  If no file is found with the specified string, an empty dictionary is returned

        ===Updates===
        -none

        ===Comments===
        """
        
        # all files in directory
        files = os.listdir(directory)
        
        # if data file already exists
        if file_name in files:
            # load data
            print 'group data found:', file_name
            with open(directory+file_name, 'rb') as pkl_file:
                group_data= pickle.load(pkl_file)
            print 'group data loaded'
        # otherwise create data structure
        else:
            # data organized as {frequency}{syn distance}{number of synapses}{polarity}[trial]{data type}{tree}[section][segment][spikes]
            print 'no group data found'
            group_data= {}

        return group_data 

        # def _convert_groupdata_to_panel(self, group_data):
        #     '''
        #     '''
        #     group_data_panel= {}
        #     for variable_key, variable in group_data.iteritems():
        #         group_data_panel[variable_key]= pd.DataFrame(variable)

        #     group_panel = pd.Panel(group_data_panel)
        #     return group_panel

        # def _add_files_to_group_panel(self, directory, search_string, group_panel, variables, **kwargs):
        #     '''
        #     '''
        #     # get list of all data files
        #     data_files = glob.glob(directory+search_string)

        #     # get list of processed data files
        #     if group_data:
        #         processed_data_files = group_data['processed']
        #     else:
        #         group_data['processed']=[]
        #         processed_data_files = group_data['processed']

    def _add_files_to_group_data(self, directory, search_string, group_data, variables=['v', 'input_times'], max_iter=1000):
        '''
        ===Args===
        -directory : directory where group data is stored including /
        -search_string : string that is unique individual data files to search directory.  typically '*data*'
        -group_data : group data structure organized as group_data{variable}{data_type}
                -group_data['v']['data_mat'] is a matrix with dimensions segments x samples

        ===Out===
        ===Updates===
        ===Comments===
        -group_data['processed'] = list of data file names that have already been added to the group 
        -variable = variable type that was recorded from neuron simulaitons, e.g. 'v' or 'gbar'
        -data_type for group_data structure:
            -'data_mat' : matrix of time series data for the specified variable
            -'input_mat': matrix of boolean time series of synpatic inputs for the corresponding data_mat (1=active input, 0=inactive input)
            -'conditions' : 
                    -'polarity', 'path', 'tree', 'sec_i', 'sec_num', 'seg_i', 'seg_num', 'input_times'
                    -e.g. group_data['v']['conditions']['polarity'] returns a list of field polarities with indices corresponsponding to rows group_data['v']['data_mat']
            't' : single 0D vector of time values (should be identical for all simulations in a given group data structure)
        '''
        # get list of new data files
        #`````````````````````````````````````````````````````````````
        # get list of all data files
        data_files = glob.glob(directory+search_string)

        # get list of processed data files
        if group_data:
            processed_data_files = group_data['processed']
        else:
            group_data['processed']=[]
            processed_data_files = group_data['processed']

        # get list of new data files
        new_data_files = list(set(data_files)-set(processed_data_files))
        print 'total data files:', len(data_files) 
        print 'new data fies:', len(new_data_files)
        
        # iterate over new files and update group_data structure
        #`````````````````````````````````````````````````````````````````
        print 'updating group data structure'
        # dictionary for temporary storage of individual simulation data

        # iterate over new data files
        for file_i, file in enumerate(new_data_files):

            if file_i<max_iter:
                print 'updating file',file_i,'of',len(new_data_files),'new files'

                # load data file
                with open(file, 'rb') as pkl_file:
                    data = pickle.load(pkl_file)

                # plt.figure()
                # plt.plot(data['v']['data'].T)
                # plt.show(block=False)
                # FIXME add input times to variables for each data files
                # iterate over variables to be updated
                for variable_i, variable in enumerate(variables):

                    if variable not in group_data:
                        group_data[variable] = copy.copy(data[variable])
                        group_data[variable]['p'] = [data[variable]['p'] for temp in  data[variable]['locations']]

                    else:
                        if variable=='input_times':
                            group_data[variable]['data']+=data[variable]['data']
                        else:
                            print group_data[variable]['data'].shape, data[variable]['data'].shape
                            group_data[variable]['data'] = np.append(group_data[variable]['data'], data[variable]['data'], axis=0)

                        group_data[variable]['locations'] += data[variable]['locations']

                        group_data[variable]['trial_id'] += data[variable]['trial_id']

                        group_data[variable]['field'] += data[variable]['field']

                        group_data[variable]['p'] += [data[variable]['p'] for temp in  data[variable]['locations']]

                        if 'syn_types' in group_data[variable]:
                            group_data[variable]['syn_types'] += data[variable]['syn_types']



                    # add file to processed list to keep track of processed files
                    group_data['processed'].append(file)

        print 'finished updating group data structure'

        print 'data shape:',group_data['v']['data'].shape
        # ouput group data structure 
        return group_data

    def _add_files_to_group_data_conditional_parameter(self, directory, search_string, group_data, parameter='pulse_freq', parameter_value=50, path_parameter=False, variables=['v', 'input_times'], max_iter=1000):
        '''
        ===Args===
        -directory : directory where group data is stored including /
        -search_string : string that is unique individual data files to search directory.  typically '*data*'
        -group_data : group data structure organized as group_data{variable}{data_type}
                -group_data['v']['data_mat'] is a matrix with dimensions segments x samples

        ===Out===
        ===Updates===
        ===Comments===
        -group_data['processed'] = list of data file names that have already been added to the group 
        -variable = variable type that was recorded from neuron simulaitons, e.g. 'v' or 'gbar'
        -data_type for group_data structure:
            -'data_mat' : matrix of time series data for the specified variable
            -'input_mat': matrix of boolean time series of synpatic inputs for the corresponding data_mat (1=active input, 0=inactive input)
            -'conditions' : 
                    -'polarity', 'path', 'tree', 'sec_i', 'sec_num', 'seg_i', 'seg_num', 'input_times'
                    -e.g. group_data['v']['conditions']['polarity'] returns a list of field polarities with indices corresponsponding to rows group_data['v']['data_mat']
            't' : single 0D vector of time values (should be identical for all simulations in a given group data structure)
        '''
        # get list of new data files
        #`````````````````````````````````````````````````````````````
        # get list of all data files
        data_files = glob.glob(directory+search_string)

        # get list of processed data files
        if group_data:
            processed_data_files = group_data['processed']
        else:
            group_data['processed']=[]
            processed_data_files = group_data['processed']

        # get list of new data files
        new_data_files = list(set(data_files)-set(processed_data_files))
        print 'total data files:', len(data_files) 
        print 'new data fies:', len(new_data_files)
        
        # iterate over new files and update group_data structure
        #`````````````````````````````````````````````````````````````````
        print 'updating group data structure'
        # dictionary for temporary storage of individual simulation data

        # iterate over new data files
        for file_i, file in enumerate(new_data_files):

            if file_i<max_iter:
                print 'updating file',file_i,'of',len(new_data_files),'new files'

                # load data file
                with open(file, 'rb') as pkl_file:
                    data = pickle.load(pkl_file)

                add_file=False
                if not path_parameter:
                    if parameter_value==data['v']['p'][parameter]:
                        add_file=True
                else: 
                    for path_key, path in data['v']['p']['p_path'].iteritems():
                        # print path_key, parameter
                        if parameter in path:
                            # print path
                            if parameter_value==path[parameter]:
                                print parameter_value, path[parameter]
                                add_file=True
                # print add_file
                if add_file:
                    print 'adding file to group_data'
                    # FIXME add input times to variables for each data files
                    # iterate over variables to be updated
                    for variable_i, variable in enumerate(variables):

                        if variable not in group_data:
                            group_data[variable] = copy.copy(data[variable])
                            group_data[variable]['p'] = [data[variable]['p'] for temp in  data[variable]['locations']]

                        else:
                            if variable=='input_times':
                                group_data[variable]['data']+=data[variable]['data']
                            else:
                                print group_data[variable]['data'].shape, data[variable]['data'].shape
                                group_data[variable]['data'] = np.append(group_data[variable]['data'], data[variable]['data'], axis=0)

                            group_data[variable]['locations'] += data[variable]['locations']

                            group_data[variable]['trial_id'] += data[variable]['trial_id']

                            group_data[variable]['field'] += data[variable]['field']

                            group_data[variable]['p'] += [data[variable]['p'] for temp in  data[variable]['locations']]

                            if 'syn_types' in group_data[variable]:
                                group_data[variable]['syn_types'] += data[variable]['syn_types']



                        # add file to processed list to keep track of processed files
                    group_data['processed'].append(file)

        print 'finished updating group data structure'

        # print 'data shape:',group_data['v']['data'].shape
        # ouput group data structure 
        print group_data.keys()
        return group_data
    
    def _add_clopath_to_group_data(self, group_data, param, rerun=False):
        '''
        ==Args==
        ==Out==
        ==Updates==
        ==Comments==
        '''
        self.clopath = Clopath()
        if 'clopath' not in group_data:

            group_data['clopath']=copy.deepcopy(group_data['v'])
            group_data['clopath']['data']=self.clopath._clopath(x=group_data['input_array']['data'], u=group_data['v']['data'], param=param)

        elif rerun==True:
            group_data['clopath']['data']=self.clopath._clopath(x=group_data['input_array']['data'], u=group_data['v']['data'], param=param)

        elif rerun==False:
            ntotal = group_data['v']['data'].shape[0]
            nclopath =group_data['clopath']['data'].shape[0]
            ndiff = ntotal-nclopath
            group_data['clopath']['data'] = np.append(group_data['clopath']['data'], np.zeros((ndiff, group_data['clopath']['data'].shape[1])), axis=0)
            group_data['clopath']['data'][nclopath:,:] = self.clopath._clopath(x=group_data['input_array']['data'][nclopath:,:], u=group_data['v']['data'][nclopath:,:], param=param)
            for key, val in group_data['v'].iteritems():
                if key != 'data':
                    group_data['clopath'][key]=copy.deepcopy(val)

        return group_data

    def _add_input_array_to_group_data(self, group_data):
        ''' add boolean array of input times to group data

            ==Args==
            -group_data :  organized as group_data{variable}{data_type}
                        -must contain group_data['input_times']['data'] as nested list of input times [segment/burst number][1d array of input times]
                                -note that each location can contain multiple entries, dependending on the number of bursts
            ==Out==
            -group_data : boolean input array (1=input, 0=no input) with same dimensions data['v']['data'] is stored
            ==Updates==
            ==Comments==
        '''
        # check if input array already in group data
        if 'input_array' not in group_data:
            # number of segments already in group_data
            ninput=0

        else:
            # number of segemnts already in group data
            ninput =group_data['input_array']['data'].shape[0]

        # time vector
        t = group_data['v']['t']
        # total number of segments in voltage data (input array to have same dimensions)
        ntotal = group_data['v']['data'].shape[0]
        # difference between final and current number of segments in input array
        ndiff = ntotal-ninput
        
        # print len(group_data['v']['locations'])
        # print len(group_data['v']['trial_id'])
        # preallocate temporary input for new data
        input_array_temp = np.zeros((ndiff,t.shape[0]))
        # iterate over locations in voltage data
        for location_i, location in enumerate(group_data['v']['locations']):
            # if location is not in input array
            if location_i >= ninput:
                # get field magnitude/polarity for the current location
                field = group_data['v']['field'][location_i]
                # get trial id for current location
                trial_id = group_data['v']['trial_id'][location_i]
                # get indices of corresponding locations and field magnitudes in the 'input_times' array
                input_i = [
                temp_i 
                for temp_i,temp in enumerate(group_data['input_times']['locations']) 
                if temp==location 
                and group_data['input_times']['field'][temp_i]==field 
                and group_data['input_times']['trial_id'][temp_i]==trial_id
                ]

                # print input_i
                # print len(group_data['input_times']['locations'])
                # print len(group_data['v']['locations'])
                # list for storing all input times for the current location
                times_all = []
                # iterate over locations indices in the input times array
                for idx in input_i:
                    # get corresponding input times
                    times = group_data['input_times']['data'][idx]
                    # add each time to list
                    for time_i, time in np.ndenumerate(times):
                        times_all.append(time)
                # print 'times:',times_all
                # create boolean vector of same length as time vector (dim: samples)
                x = np.zeros(t.shape)
                # cut off number of decimals to prevent rounding error
                t = np.around(t, 4)
                # find input times and to boolean vector
                for t_i, t_t in enumerate(times_all):
                            
                    if t_t not in t:
                        x[int(np.argmin(np.absolute(t-t_t)))]=1
                        # print np.absolute(t-t_t)
                        # print min
                        print 'test:',int(np.argmin(np.absolute(t-t_t)))
                    else:
                        # store input times
                        x[np.where(t==t_t)]=1

                # add to full group data array
                input_array_temp[location_i-ninput, :]=x
        
        # if input array already in group_data, append the new temporary array
        if 'input_array' in group_data:
            group_data['input_array']['data'] = np.append(group_data['input_array']['data'], input_array_temp, axis=0)

        # if input array not in group data, copy temporary array 
        else:
            group_data['input_array']={'data':input_array_temp}
            group_data['input_array']['data'][ninput:,:]=input_array_temp

        # copy all other condition lists from voltage data (dimensions should be the same)
        for key, val in group_data['v'].iteritems():
            if key!='data':
                group_data['input_array'][key]= copy.copy(val)

        return group_data

    def _add_parameter_to_group_data(self, group_data, parameter, path=False):
        '''
        '''
        for variable_key, variable in group_data.iteritems():
            if variable_key != 'processed':
                if parameter not in variable:
                    variable[parameter]=[]

                nparam = len(variable[parameter])
                ntotal = len(variable['locations'])

                if ntotal > nparam:

                    for location_i, location in enumerate(variable['locations']):
                        if location_i>=nparam:

                            # print variable_key, parameter, variable[parameter]
                            # print type(variable)
                            if not path:
                                variable[parameter].append(variable['p'][location_i][parameter])
                            else:
                                # get all paths that the 
                                variable[parameter].append([])
                                if 'None' not in variable['p_path'][location_i]:
                                    for p_path in variable['p_path'][location_i]:

                                        variable[parameter][-1].append(p_path[parameter])

        return group_data

    def _add_path_to_group_data(self, group_data):
        '''
        ==Args==
        -group_data : as group_data[variable][info type]
        ==Out==

        ==Updates==
        -adds an info type named 'path' organized as [location][path number]. Since each location can belong to multiple paths, each location entry is a list of paths.  If the location doesnt belong to any paths 'None' is put in the list as a place holder
        ==Comments==
        '''
        print 'adding path info'
        # print group_data['processed']
        print 'group_data shape:',group_data['v']['data'].shape
        # iterate over variables
        for variable_key, variable in group_data.iteritems(): 
            
            if variable_key != 'processed':

                # print variable_keae.y, type(variable)
                # add path info list if doesnt exist
                if 'p_path' not in variable:
                    variable['p_path'] =[]
                if 'path_name' not in variable:
                    variable['path_name']=[]

                # print variable_key, type(variable), variable['path_name']
                # number of path entries
                npath = len(variable['p_path'])
                # total number of entries
                ntotal = len(variable['locations'])

                # if there are unprocessed locations
                if ntotal>npath:

                    # iterate over locations
                    for location_i, location in enumerate(variable['locations']):
                        # holder for paths for the current location (can be more than one)
                        p_paths = []
                        path_names=[]
                        # if location hasnt been added to path yet
                        if location_i>=npath:

                            # get parameter dictionary
                            p = variable['p'][location_i]
                            # iterate over paths
                            for path_key, path in p['p_path'].iteritems():
                                # if the current location is in the current path
                                if location in path['syn_idx']:
                                    # add path to list
                                    p_paths.append(path)
                                    path_names.append(path_key)

                            # if location doesnt belong to any paths, add 'None' as placeholder
                            if len(p_paths)==0:
                                p_paths.append('None')
                                path_names.append('None')
                            
                            # add to main list
                            variable['p_path'].append(copy.copy(p_paths))
                            variable['path_name'].append(copy.copy(path_names))

        return group_data

    def _add_polarity_to_group_data(self, group_data):
        '''
        '''
        # iterate over variables
        for variable_key, variable in group_data.iteritems(): 
            
            if variable_key !='processed':

                variable['polarity']=[]
                for field_i, field in enumerate(variable['field']):
                    if float(field)==0:
                        variable['polarity'].append('control')
                    elif float(field)<0:
                        variable['polarity'].append('cathodal')
                    elif float(field)>0:
                        variable['polarity'].append('anodal')
        return group_data

    def _save_group_data(self, group_data, directory, file_name):
        '''
        ===Args===
        ===Out===
        ===Updates===
        ===Comments===
        '''
        print 'saving group data'
        with open(directory+file_name, 'wb') as output:
            pickle.dump(group_data, output,protocol=pickle.HIGHEST_PROTOCOL)

        print 'group data saved as:', file_name 

    def _split_group_data(self, group_data, directory, file_name, conditions):
        '''
        '''
        assert len(conditions)==1, 'multiple conditions given; only one at a time please'

        group_data_split={}
        for variable_key, variable in group_data.iteritems():
            
            if variable_key=='processed':
                group_data_split['processed']=variable

            else:
                ig = IndexGetter(group_data=group_data, variable=variable_key, conditions=conditions)
                idx_set = ig.idx_sets[0]
                combo = ig.combos[0]
                
                group_data_split[variable_key]={}
                for info_type, info in variable.iteritems():

                    if info_type=='t':
                        info_new=copy.deepcopy(info)
                    else:
                        if type(info)==list:

                            info_new = [ info[temp] for temp in idx_set ]

                        elif type(info)==np.ndarray:
                            print info_type, info.shape, len(idx_set)
                            info_new = info[idx_set,:]

                    group_data_split[variable_key][info_type]=info_new

        print conditions
        split_name=''
        for condition in conditions:
            print condition
            add_string = '_' + condition[0] + '_' +str(condition[1])
            split_name+= add_string

        file_name_new = file_name[:file_name.index('.pkl')] + split_name +'.pkl'

        return group_data_split, file_name_new

    def _add_spikes_to_group_data(self, group_data, threshold=-30):
        '''
        '''
        print 'adding spike times to group_data'
        self.Spikes= Spikes()
        if 'spike_times' not in group_data:
            group_data['spike_times']={}
            group_data['spike_idx']={}
            group_data['spike_train']={}
            spike_idx, spike_train, spike_times = self.Spikes._get_spikes(data=group_data['v']['data'], threshold=-30, dt=.025)
            group_data['spike_times']['data']=spike_times
            group_data['spike_idx']['data']=spike_idx
            group_data['spike_train']['data'] = spike_train
            for key, val in group_data['v'].iteritems():
                if key!='data':
                    group_data['spike_times'][key]=copy.deepcopy(val)
                    group_data['spike_idx'][key]=copy.deepcopy(val)
                    group_data['spike_train'][key]=copy.deepcopy(val)

        else:
            ntotal = group_data['v']['data'].shape[0]
            nspikes = len(group_data['spike_times']['data'])
            ndiff = ntotal-nspikes
            if ndiff>0:

                spike_idx, spike_train, spike_times = self.Spikes._get_spikes(data=group_data['v']['data'][-ndiff:,:], threshold=-30, dt=.025)
                group_data['spike_times']['data']+=spike_times
                group_data['spike_idx']['data']+=spike_idx
                group_data['spike_train']['data'] = np.append(group_data['spike_train']['data'], spike_train, axis=0)
                for key, val in group_data['v'].iteritems():
                    if key!='data':
                        group_data['spike_times'][key]=copy.deepcopy(val)
                        group_data['spike_idx'][key]=copy.deepcopy(val)
                        group_data['spike_train'][key]=copy.deepcopy(val)

        return group_data

    def _add_xcorr_to_group_data(self, group_data):
        '''
        '''
        print 'adding cross correlation with soma to group data'
        spike_train = group_data['spike_train']['data']
        if 'soma_xcorr' not in group_data:
            group_data['soma_xcorr']={}
            for key, val in group_data['spike_train'].iteritems():
                if key!= 'data':
                    group_data['soma_xcorr'][key] = copy.deepcopy(val)

            group_data['soma_xcorr']['data']= np.zeros((spike_train.shape[0], spike_train.shape[1]+spike_train.shape[1]-1))
            nxcorr=0

        else:
            ntotal = group_data['spike_train']['data'].shape[0]
            nxcorr =group_data['soma_xcorr']['data'].shape[0]
            ndiff = ntotal-nxcorr
            group_data['soma_xcorr']['data'] = np.append(group_data['soma_xcorr']['data'], np.zeros((ndiff, group_data['soma_xcorr']['data'].shape[1])), axis=0)
            for key, val in group_data['spike_train'].iteritems():
                if key!= 'data':
                    group_data['soma_xcorr'][key] = copy.deepcopy(val)

        soma_idx = IndexGetter()._find_soma(group_data=group_data, variable='spike_train')
        for loc_i, loc in enumerate(group_data['spike_train']['locations']):
            if loc_i >= nxcorr:

                print loc_i, 'of', spike_train.shape[0]
                group_data['soma_xcorr']['data'][loc_i,:] = np.correlate(spike_train[loc_i, :], spike_train[soma_idx[loc_i],:], mode='full')

        group_data['soma_xcorr']['soma_idx']=soma_idx

        return group_data





        # soma_idx = IndexGetter()._find_soma(group_data=group_data, variable='spike_train')
        # spike_train = group_data['spike_train']['data']
        # xcorr = np.zeros((spike_train.shape[0], spike_train.shape[1]+spike_train.shape[1]-1))
        # for loc_i, loc in enumerate(group_data['spike_train']['locations']):
        #     print loc_i, 'of', spike_train.shape[0]
        #     xcorr[loc_i,:] = np.correlate(spike_train[loc_i, :], spike_train[soma_idx[loc_i],:], mode='full')


        # group_data['soma_xcorr']={'data':xcorr}
        # for key, val in group_data['spike_train'].iteritems():
        #     if key!= 'data':
        #         group_data['soma_xcorr'][key] = copy.deepcopy(val)
        # group_data['soma_xcorr']['soma_idx']=soma_idx

        # return group_data, 

