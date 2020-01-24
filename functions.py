'''
'''
import numpy as np
import os
import pandas as pd
from matplotlib import pyplot as plt
from scipy import stats
from collections import OrderedDict 
import copy
import glob
import pickle
from sklearn.cross_decomposition import CCA
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import pdb
import inspect
import math
import itertools
from scipy.optimize import curve_fit
import analysis
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


#############################################################################
# common functions
#############################################################################
# check if dataframe exists in namespace
def _exists(variable):
    '''
    '''
    if variable in globals():
        return True
    else:
        return False

def _initialize_column(df, col, val=None):
    '''
    '''
    if type(col)==str:
        col=[col]
    for col_name in col:
        if col_name not in df:
            df[col_name]=val
    return df
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
# variable generation for group level analysis
#############################################################################
class VarFuncs:
    ''' extract higher level analysis variables from individual cells/experiments.  outputs should be a dataframe for each preprocessed file, which are then combined (see VarGen) into a single dataframe for all experiments.
    '''
    def __init__(self, ):
        '''
        '''
        pass

    def _build_conditions_df(self, pre, var, **kwargs):
        ''' standard dataframe containing common simulation details
        Arguments
        --------------
        ~pre: preprocessed data file
        ~var: variable type that is being converted to df (e.g. 'v')

        Return
        -------------
        ~df_new: dataframe with indices corresponding to traces in pre[var]['data'].  Columns correspond to conditions, which are specified below

        Comments
        --------------
        '''
        # iterate over traces
        df_new = pd.DataFrame(dtype='object')
        df_new['trial_id'] = pre[var]['trial_id']
        df_new['location'] = pre[var]['locations']
        df_new['field'] = pre[var]['field']
        df_new['tree'] = zip(*pre[var]['locations'])[0]
        df_new['sec'] = zip(*pre[var]['locations'])[1]
        df_new['seg'] = zip(*pre[var]['locations'])[2]
        df_new['dt'] = pre[var]['p']['dt']
        df_new['tstop'] = pre[var]['p']['tstop']
        # ac field parameters
        #---------------------
        if 'ac_field' in pre[var]['p']:
            for key, val in pre[var]['p'].iteritems():
                if 'ac_field' in key:
                    df_new[key]=val

        # create columns for parameters in p
        if 'params' in kwargs:
            # Store all parameters
            if kwargs['params']==None or kwargs['params']=='all':
                for param in pre[var]['p'].keys():
                    if param not in df_new:

                        df_new[param]=None
                        for i in df_new.index:
                            value = pre[var]['p'][param]
                            if type(value)==list or type(value)==dict:
                                if len(value)==0:
                                    value=np.nan
                                df_new.at[i,param] = value
                            else:
                                df_new.at[i,param] = value

            else:
                # store specified parameters
                for param in kwargs['params']:
                    if param not in df_new:
                        df_new[param]=None
                        for i in df_new.index:
                            df_new.at[i,param]=pre[var]['p'][param]

        # iterate over traces
        df_new['paths']=None
        df_new['path']=None
        df_new['input_times']=None
        for i, row in df_new.iterrows():
            location = df_new['location'][i]
            # get distance from soma
            #------------------------
            dist = pre[var]['p']['seg_dist'][row.tree][row.sec][row.seg]
            # get morphology entry (overall seg index, location key, x, y, z, dimaeter, parent overall index)
            #----------------------
            morpho = pre[var]['p']['morpho'][row.tree]
            morpho_key = '_'.join([row.tree, str(row.sec), str(row.seg)])
            morpho_val = [seg for sec in morpho for seg in sec if morpho_key in seg][0]
            # get segment distance from soma
            #------------------------
            df_new.at[i, 'seg_dist']=dist
            if 'morpho' not in df_new:
                df_new['morpho']=None
            df_new.at[i, 'morpho']=morpho_val
            # get synaptic input times
            #--------------------------
            # pre['input_times']['data'] is nested list as [burst_i][input_time], locations list contains an entry for each burst (locations are repeated if it receives multiple bursts)
            if 'input_times' in pre:
                if 'input_times' not in df_new:
                    df_new['input_times']=None
                if df_new['location'][i] in pre['input_times']['locations']:
                    # indices of input times for each burst
                    input_time_is = [temp_i for temp_i, temp in enumerate(pre['input_times']['locations']) if temp==location]
                    # print input_time_is
                    # nested list of input times for current location [bursts][input times]
                    input_times = [pre['input_times']['data'][temp_i] for temp_i in input_time_is]
                    # print pre['input_times']['data']
                    # print 'times', input_times
                    input_times = np.array(list((set(np.concatenate(input_times)))))
                    # print input_times

                    # print input_times.shape
                    # flattened array of input times for current location
                    input_times = np.array(input_times).flatten()
                    # print input_times
                    # print input_times.shape
                    # sort input times
                    input_times.sort()
                    # store in df
                    df_new.at[i, 'input_times']=input_times
            # synaptic input path parameters
            #----------------------------------
            paths=[]
            n_paths=len(pre[var]['p']['p_path'].keys())
            df_new.at[i, 'n_paths']=n_paths
            # df_new['paths']=None
            # iterate over input pathways
            for path_key, path in pre[var]['p']['p_path'].iteritems():
                # list of paths that contain the current location
                #------------------------------------------------
                if row.location in path['syn_idx']:
                    temp_i = path['syn_idx'].index(row.location)
                    paths.append(path_key)
                    # synaptic weight at current location
                    #------------------------------------
                    if 'w_idx' in path:
                        df_new.at[i, 'w'] = path['w_idx'][temp_i]
                    df_new.at[i,'path']=path_key
                # iterate over path parameters
                for path_param, param_val in path.iteritems(): 
                    colname = '_'.join(['path', path_key, path_param])
                    if colname not in df_new:
                        df_new[colname]=None
                    # add all path parameters to df
                    #------------------------------
                    df_new.at[i,colname]=param_val
            # add paths that current entry belongs to 
            df_new.at[i,'paths']=paths
        # replace none objects with na
        df_new.fillna(value=pd.np.nan, inplace=True)

        return df_new

    def _pass(self, pre, **kwargs):
        '''
        '''
        return pre

    def _get_vtrace(self, pre, df=[], keep=[], **kwargs):
        '''
        '''
        df_new = self._build_conditions_df(pre, var='v',**kwargs)
        for i, val in enumerate(pre['v']['data']):
            if 'data_v' not in df_new:
                df_new['data_v']=None
            # print val.shape
            if len(val.shape)>1:
                print val.shape
            df_new.at[i,'data_v'] = np.array(val)

        return df_new

    def _get_spikes(self, pre, df=[], keep=[], threshold=-30, **kwargs):
        '''
        '''
        print 'threshold', threshold
        # get time step
        dt = pre['v']['p']['dt']
        # build conditions df
        df_new = self._build_conditions_df(pre, var='v', **kwargs)
        # preallocate columns
        df_new['spike_idx']=None
        df_new['spike_train']=None
        df_new['spike_times']=None
        # iterate over voltage traces
        for i, data in enumerate(pre['v']['data']):
            # get spike info
            spike_idx, spike_train, spike_times = analysis.Spikes()._get_spikes(
                data=data, dt=dt, threshold=threshold)
            # add to dataframe
            # spike_idx: list of indices where spike onset occurs, length=number of spikes
            df_new.at[i, 'spike_idx']=spike_idx
            # spike_train: boolean array with same shape as voltage data, where spike onsets are stored as 1
            df_new.at[i, 'spike_train']=spike_train
            # spike_times: list of spike times in ms, length=number of spikes 
            df_new.at[i, 'spike_times']=spike_times

        return df_new

    def _get_w_clopath(self, pre, df=[], keep=[], **kwargs):
        '''
        '''
        df_new = self._build_conditions_df(pre, var='input_times', **kwargs)
        # get corresponding indices in voltage list
        voltage_is = [i for i, loc in enumerate(pre['v']['locations']) if loc in pre['input_times']['locations'] ]
        # convert voltage data to array
        voltage_array = np.array(pre['v']['data'])[voltage_is, :]

        # convert input times to input array
        #--------------------------------------------------------------
        # preallocate
        input_array = np.zeros(voltage_array.shape)
        # time vector
        t = pre['v']['t'] 
        # cut off number of decimals to prevent rounding error
        t = np.around(t, 4)
        # iterate over locations
        for loc_i, input_times in enumerate(pre['input_times']['data']):

            # find input times and to boolean vector
            # iterate over input times
            for t_i, t_t in enumerate(input_times):
                # if input time does not exactly match a sample in time vector    
                if t_t not in t:
                    # set input time to nearest sample
                    input_array[loc_i, int(np.argmin(np.absolute(t-t_t)))]=1
                else:
                    # store input times
                    input_array[loc_i, np.where(t==t_t)]=1
        # get parameters
        param = pre['input_times']['p']
        # update with clopath parameters
        if 'clopath_param' in kwargs:
            param.update(kwargs['clopath_param'])
        # get weight array
        w_clopath = analysis.Clopath()._clopath(x=input_array, u=voltage_array, fs=1/param['dt'], w0=0.5, homeostatic=False, param=param)
        # preallocate column for storing weight time series
        if 'w_clopath' not in df_new:
            df_new['w_clopath']=None
        # iterate over input vectors and store corresponding weight timeseries
        for i, x in enumerate(pre['input_times']['data']):
            df_new.at[i, 'w_clopath'] = copy.copy(w_clopath[i,:])
            df_new.at[i, 'w_final'] = w_clopath[i,-1]
            df_new.at[i, 'w_initial']=w_clopath[i, 0]
            df_new.at[i, 'dw_clopath'] = w_clopath[i,-1]/w_clopath[i,0]

        return df_new

class VarGen:
    ''' apply VarFuncs to generate group variables
    '''
    def __init__(self, ):
        '''
        '''
        pass

    def _process_new_data_df(self, group_df, preprocessed_directory, search_string='data', functions=[], kwlist=[], rerun=[], keep=[], file_limit=[], ):
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
        print 'search_string', search_string
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
            print 'functions to rerun over all slices:', [f.__name__ for f in rerun]

        # iterate over new files and update group_data structure
        #`````````````````````````````````````````````````````````````````
        print 'updating group data structure'

        # temporary df for storing rerun function outputs
        df_update = pd.DataFrame()

        # suffixes for merging df_update with group_df
        old_suffix = '_old'
        rerun_suffix = '_rerun'

        print 'iter files:', len(iterfiles)

        # iterate over data files
        for file_i, file in enumerate(iterfiles):
            # print file_i

            # apply file limit
            #------------------
            if file_limit and file_i> file_limit:
                continue
            else:
                print 'processing file', file_i
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
                        # print kwlist[func_i]
                        df_temp = func(pre, df=[], keep=[], **kwlist[func_i])
                        # # print 'test'
                        new_dfs.append(df_temp)

                    # if rerunning current func, add to rerun list
                    #-------------------------------------------
                    elif func in rerun:
                        print 'rerunning function', func.__name__, 'on file', file_i
                        # print kwlist[func_i]
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

    def _process_new_data(self, group_df, preprocessed_directory, search_string='data', functions=[], kwlist=[], rerun=[], keep=[], file_limit=[], ):
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
        print 'search_string', search_string
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
            print 'functions to rerun over all slices:', [f.__name__ for f in rerun]

        # iterate over new files and update group_data structure
        #`````````````````````````````````````````````````````````````````
        print 'updating group data structure'

        # temporary df for storing rerun function outputs
        df_update = pd.DataFrame()

        # suffixes for merging df_update with group_df
        old_suffix = '_old'
        rerun_suffix = '_rerun'

        print 'iter files:', len(iterfiles)

        # iterate over data files
        for file_i, file in enumerate(iterfiles):
            # print file_i

            # apply file limit
            #------------------
            if file_limit and file_i> file_limit:
                continue
            else:
                print 'processing file', file_i
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
                        # print kwlist[func_i]
                        df_temp = func(pre, df=[], keep=[], **kwlist[func_i])
                        # # print 'test'
                        new_dfs.append(df_temp)

                    # if rerunning current func, add to rerun list
                    #-------------------------------------------
                    elif func in rerun:
                        print 'rerunning function', func.__name__, 'on file', file_i
                        # print kwlist[func_i]
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

    def _vargen(self, variable, directory, functions=[], kwlist=[], rerun=[], keep=[], file_limit=[], write_protocol='.pkl', save=False, size_limit=2E9, **kwargs):
        '''
        '''
        filename = variable+write_protocol
        df = analysis._load_group_data(directory=directory,filename=variable)
        # df = df.reset_index()
        df_copy = copy.deepcopy(df)

        # print sorted(df.keys())
        df = self._process_new_data_df(group_df=df, preprocessed_directory=directory, functions=functions, kwlist=kwlist, rerun=rerun, keep=keep, file_limit=file_limit, **kwargs)
        # print df
        # print df_copy.equals(df)

        df = analysis._save_group_data(df=df, directory=directory, variable=variable, extension='.pkl', check_against=None, **kwargs)

        # # get size of df
        # #---------------------------
        # df_size = sys.getsizeof(df)
        # # if df is too big, break into manageable parts
        # if df_size>size_limit:
        #     print 'df too large to save, breaking into smaller parts'
        #     n_dfs=df_size/size_limit
        #     df_split = np.array_split(df, n_dfs)
        #     for df_i, df_chunk in enumerate(df_split):
        #         # print df_chunk
        #         df_chunk.to_pickle(directory+variable+'_'+str(df_i)+'.pkl')

        # # save group data
        # #----------------

        # elif save:
        #     # check write protocol
        #     if write_protocol=='.pkl':
        #         df.to_pickle(directory+filename)
        #     elif write_protocol=='.h5':
        #         df.to_hdf(directory+filename, key=variable, mode='w', format='fixed')

        # elif not df_copy.equals(df):
        #     print 'saving updated group data'
        #     # check write protocol
        #     if write_protocol=='.pkl':
        #         df.to_pickle(directory+filename)
        #     elif write_protocol=='.h5':
        #         df.to_hdf(directory+filename, key=variable, mode='w', format='fixed')
        return df

    def _vtrace(self, **kwargs):
        '''
        '''
        varfuncs=VarFuncs()
        directory = 'Data/'+kwargs['experiment']+'/'
        variable='vtrace_df'
        functions=[varfuncs._get_vtrace]
        kwlist=[{}]
        rerun=[varfuncs._get_vtrace]
        keep=[]
        file_limit=[]
        df = self._vargen(variable=variable, directory=directory, functions=functions, kwlist=kwlist, rerun=rerun, keep=keep, file_limit=file_limit)
        return df

    def _spikes(self, **kwargs):
        '''
        '''
        varfuncs = VarFuncs()
        directory = 'Data/'+kwargs['experiment']+'/'
        variable='spikes_df'
        functions=[varfuncs._get_spikes]
        kwlist=[{'threshold':-30}]
        rerun=[varfuncs._get_spikes]
        keep=[]
        file_limit=[]
        df = self._vargen(variable=variable, directory=directory, functions=functions, kwlist=kwlist, rerun=rerun, keep=keep, file_limit=file_limit)
        return df

    def _w_clopath(self, **kwargs):
        '''
        '''
        varfuncs=VarFuncs()
        directory = 'Data/'+kwargs['experiment']+'/'
        variable='w_clopath_df'
        functions=[varfuncs._get_w_clopath]
        kwlist=[{'clopath_param':clopath_param_all}]
        # kwlist=[kwargs]
        rerun=[varfuncs._get_w_clopath]
        keep=[]
        file_limit=[]
        df = self._vargen(variable=variable, directory=directory, functions=functions, kwlist=kwlist, rerun=rerun, keep=keep, file_limit=file_limit)
        return df

    def _w_clopath_stdp(self, **kwargs):
        '''
        '''
        varfuncs=VarFuncs()
        directory = 'Data/'+kwargs['experiment']+'/'
        variable='w_clopath_stdp_df'
        functions=[varfuncs._get_w_clopath]
        kwlist=[{'clopath_param':clopath_param_all_stdp}]
        rerun=[varfuncs._get_w_clopath]
        keep=[]
        file_limit=[]
        df = self._vargen(variable=variable, directory=directory, functions=functions, kwlist=kwlist, rerun=rerun, keep=keep, file_limit=file_limit)
        return df

#############################################################################
# group dataframe manipulation to update before analysis
#############################################################################
class DfFuncs:
    ''' applied to create new columns after df has been generated, all function should take a dataframe as an argument and return the same dataframe with new columns added

    # FIXME add functionality to check for data that has already been processed to speed up analysis
    '''
    def __init__(self, ):
        '''
        '''
        pass
    
    def _to_string(self, df, colnames=[], suffix='_str'):
        ''' convert all cells in column to str
        '''
        for i, colname in enumerate(colnames):
            colname_new = colname+suffix
            df[colname_new] = df[colname].apply(str)
        return df

    def _truncate_arrays(self, df, data_col='data_v', sort_col='path_1_delay', minsize=True):
        '''
        '''
        def _truncate(row, size=None):
                return row[:size]
        # print df
        # df = df.reset_index().set_index(colsort)
        df = df.reset_index().set_index(sort_col)
        df_new=pd.DataFrame()
        for i in df.index.unique():
            df_temp = copy.deepcopy(df.loc[i, data_col])
            sizes = [row.shape[0] for row in df_temp]
            print set(sizes)
            size=min(sizes)
            kwargs = {'size':size}
            df_temp = df_temp.apply(_truncate, **kwargs)
            sizes = [row.shape[0] for row in df_temp]
            print set(sizes)
            df_new=df_new.append(df_temp.reset_index(), ignore_index=True)
            # df.at[i, data_col]=df_temp







            # df_temp = copy.deepcopy(df.loc[i])
            # print df_temp
            # sizes = [row[0].shape[0] for row in df_temp.values]
            # print 'sizes 1:', set(sizes)

            
            # df_temp = df_temp.apply(_truncate)
            # sizes = [row[0].shape[0] for row in df_temp.values]
            # print 'sizes 2:', set(sizes)

            
        # df.reset_index()
        return df_new

    def _get_area(self, df, index, **kwargs):
        '''
        '''
        # reorganize dataframe
        df = df.set_index(['trial_id','location','field'])
        # column name to add
        coladd = ''
        if coladd not in df or rerun:
            # iterate over rows
            for row_i, row in df.iterrows():
                # time point to measure polarization at
                time_i = int(time/row['dt'])
                # index of corresponding control trace
                control_i = (row_i[0], row_i[1], 0.0)
                # get voltage value from control trace
                control_val = df.loc[control_i].data_v[time_i]
                # get voltage value for current trace
                val = row['data_v'][time_i]
                # set value in 'polarization column'
                df.at[row_i, 'polarization']=val-control_val

        df = df.reset_index()

        return df

    def _get_w_clopath(self, vtrace_df, clopath_param, split_columns=[], rerun=False,**kwargs):
        ''' run clopath algorithm to get synaptic weight evolution for entire df

        Arguments
        ----------
        ~vtrace_df: df must contain a column called input_times (each entry is a list of synaptic input times in ms) and data_v (each entry is a 1d array of voltage data)
        ~clopath_param: clopath rule parameters

        Returns
        ----------
        ~df: new data frame with columns:
                ~w_clopath: timeseriesof synaptic weights
                ~w_final: final value for w_clopath at end of simulation
                ~w_initial: initial value
                ~dw_clopath: final/initial
        '''
        # print progress to terminal
        print 'applying df function:', inspect.stack()[0][3]
        def _get_w(df, clopath_param):
            '''
            '''
            # create voltage array
            v = analysis._2array(df.data_v)
            # drop voltage data
            df = df.drop('data_v', axis=1)
            # fill input times with nan and get values
            input_times = df.input_times.fillna(value=np.nan).values
            # create time vector from dt and tstop
            dt = df.dt.unique()[0]
            tstop = df.tstop.unique()[0]
            # t = np.arange(0., tstop+dt, dt)
            print v.shape[1]
            # t = np.arange(0, dt*v.shape[1]-dt, dt)
            t = np.linspace(0, dt*v.shape[1], num=v.shape[1])
            print t.shape
            # create input array
            input_array = analysis._input_times_2array(input_times, t)
            print v.shape, input_array.shape
            # assert same shape for voltage and input arrays
            assert v.shape==input_array.shape, 'input array and voltage array must be same shape'
            # run clopath algorithm to get weights
            w_clopath = analysis.Clopath()._clopath(x=input_array, u=v, fs=1/dt, w0=0.5, homeostatic=False, param=clopath_param)
            # preallocate column for storing weight time series
            if 'w_clopath' not in df:
                df['w_clopath']=None
            # iterate over dataframe rows and add corresponding clopath data
            for i, row in df.iterrows():
                df.at[i, 'w_clopath'] = copy.copy(w_clopath[i,:])
                df.at[i, 'w_final'] = w_clopath[i,-1]
                df.at[i, 'w_initial']=w_clopath[i, 0]
                df.at[i, 'dw_clopath'] = w_clopath[i,-1]/w_clopath[i,0]
            return df
        
        # reset index and copy df
        df = vtrace_df.reset_index().copy()
        # if voltage traces are different length due to a set of parameters, split dataframe acccording to those parameters (so that voltage data can be made into an array, then recombine dataframes)
        if len(split_columns)>0:
            df = df.set_index(split_columns)
            split_indices = df.index.unique().values
            df_new = pd.DataFrame()
            for index in split_indices:
                print 'index:', index
                print 'df_new:', df_new.shape
                df_temp = df.loc[index].copy()
                df_temp.reset_index(inplace=True)
                df_temp = _get_w(df=df_temp, clopath_param=clopath_param)
                if df_new.empty:
                    print 'df_new empty'
                    df_new = df_temp.copy()
                else:
                    print 'appending'
                    df_new = df_new.append(df_temp, ignore_index=True)
                    print 'df_new:', df_new.shape

        else:
            df_new = _get_w(df=df, clopath_param=clopath_param)

        return df_new

    def _get_spikes(self, vtrace_df, threshold=-40, split_columns=[], **kwargs):
        ''' get spike times for each trace and store in df

        Arguments
        ----------
        ~vtrace_df: df must contain a column called data_v (each entry is a 1d array of voltage data)
        ~threshold: spike detection threshold (mV)
        ~split_columns: list of columns to group the data frame by.  if some voltage traces are different length due to a set of parameters, split dataframe accordingly so that voltage data can be made into arrays, then recombine dataframes after processing the individual arrays

        Returns
        ----------
        ~df: new data frame with columns:
                ~w_clopath: timeseriesof synaptic weights
                ~w_final: final value for w_clopath at end of simulation
                ~w_initial: initial value
                ~dw_clopath: final/initial
        '''
        # print progress to terminal
        print 'applying df function:', inspect.stack()[0][3]
        def _spikes(df, threshold):
            '''
            '''
            # reset dataframe index to range index
            df.reset_index(inplace=True)
            # create voltage array
            v = analysis._2array(df.data_v)
            # drop voltage data
            df = df.drop('data_v', axis=1)
            # create time vector from dt and tstop
            dt = df.dt.unique()[0]
            # run clopath algorithm to get weights
            spike_idx, spike_train, spike_times = analysis.Spikes()._get_spikes(data=v, threshold=threshold, dt=dt,)
            # preallocate column for storing weight time series
            if 'spike_idx' not in df:
                df['spike_idx']=None
            if 'spike_train' not in df:
                df['spike_train']=None
            if 'spike_times' not in df:
                df['spike_times']=None
            # iterate over dataframe rows and add corresponding spike data
            for i, row in df.iterrows():
                df.at[i, 'spike_idx'] = spike_idx[i]
                df.at[i, 'spike_train'] = spike_train[i,:]
                df.at[i, 'spike_times'] = spike_times[i]
            return df
        
        # reset index and copy df
        df = vtrace_df.reset_index().copy()
        # if voltage traces are different length due to a set of parameters, split dataframe acccording to those parameters (so that voltage data can be made into an array, then recombine dataframes)
        if len(split_columns)>0:
            df = df.set_index(split_columns)
            split_indices = df.index.unique().values
            df_new = pd.DataFrame()
            for index in split_indices:
                print 'index:', index
                print 'df_new:', df_new.shape
                df_temp = df.loc[index].copy()
                df_temp.reset_index(inplace=True)
                df_temp = _spikes(df=df_temp, threshold=threshold)
                if df_new.empty:
                    print 'df_new empty'
                    df_new = df_temp.copy()
                else:
                    print 'appending'
                    df_new = df_new.append(df_temp, ignore_index=True)
                    print 'df_new:', df_new.shape

        else:
            df_new = _spikes(df=df, threshold=threshold)

        return df_new
        
    def _get_polarization(self, df, time=10, rerun=False, file_limit=[], **kwargs):
        ''' get change in membrane potential in each recorded compartment in response to DCS
        
        Arguments
        ----
        -df : dataframe with olumn 'data_v' containing voltage data in each compartment, and 'dt' column
        -time : time point to take the difference between stimulated and corresponding control trace (ms)

        Return
        ------
        -df : dataframe with new 'polarization' column

        '''
        # print progress to terminal
        print 'applying df function:', inspect.stack()[0][3]
        # reorganize dataframe
        df = df.set_index(['trial_id','location','field'])
        # column name to add
        coladd = 'polarization'
        # if column doesn't exist or rerun is specified, process all rows
        #----------------------------------------------------------------
        if coladd not in df or rerun:
            # print progress
            print 'processing all rows'
            # iterate over rows
            cnt=-1
            for row_i, row in df.iterrows():
                cnt+=1
                # check for file limit
                if file_limit and cnt> file_limit:
                    print 'file limit reached:', file_limit
                    break
                # time point to measure polarization at
                time_i = int(time/row['dt'])
                # index of corresponding control trace
                control_i = (row_i[0], row_i[1], 0.0,)
                # print df.loc[control_i].data_v.drop_duplicates()
                # get voltage value from control trace
                control_val = df.loc[control_i].data_v.iloc[0][time_i]
                # get voltage value for current trace
                # print row['data_v']
                val = row['data_v'][time_i]
                # set value in 'polarization column'
                df.at[row_i, 'polarization']=val-control_val
        # if column already exists but some values haven't been processed
        #----------------------------------------------------------------
        # FIXME apply this functionality to other df functions!!!!!!!!!!!!!!!
        elif coladd in df and df[coladd].isnull().values.any():
            # get only unprocessed rows (index must match original df)
            df_temp = df[df[coladd].isnull()].copy()
            # print progress
            print 'unprocessed rows found:', df_temp.shape[0]
            # iterate over rows
            cnt=-1
            for row_i, row in df_temp.iterrows():
                cnt+=1
                # check for file limit
                if file_limit and cnt> file_limit:
                    print 'file limit reached:', file_limit
                    break
                # time point to measure polarization at
                time_i = int(time/row['dt'])
                # index of corresponding control trace
                control_i = (row_i[0], row_i[1], 0.0)
                # get voltage value from control trace
                control_val = df.loc[control_i].data_v[time_i]
                # get voltage value for current trace
                val = row['data_v'][time_i]
                # set value in 'polarization column'
                df.at[row_i, 'polarization']=val-control_val
        # reset index
        df = df.reset_index()

        return df

    def _get_xcorr_soma(self, spikes_df, rerun=False, file_limit=[], **kwargs):
        ''' for each recorded compartment get cross correlation of spikes wit hsoma
        
        Arguments
        ---------
        ~spikes_df : dataframe with column 'spike_train' containing boolean vector of spike times, and 'dt' column
        ~variable : column name in df

        Return
        ------
        ~spikes_df : dataframe with new column as 'xcorr_soma'

        Notes
        -----

        '''
        
        # print progress to terminal
        print 'applying df function:', inspect.stack()[0][3]
        # reorganize dataframe
        df = spikes_df.set_index(['trial_id','location','field', 'path'])
        # drop duplicated rows
        df = df[~df.index.duplicated()]
        dt = df.dt.unique()[0]
        # column name to add
        coladd = 'xcorr_soma'
        coladd_t ='xcorr_t'
        # df[coladd]=None
        # if column doesn't exist or rerun is specified, process all rows
        #----------------------------------------------------------------
        if coladd not in df or rerun:
            # print progress
            print 'processing all rows'
            if coladd not in df:
                df[coladd]=None
            if coladd_t not in df:
                df[coladd_t]=None
            # iterate over rows
            cnt=-1
            for row_i, row in df.iterrows():
                print row_i
                cnt+=1
                # check for file limit
                if file_limit and cnt> file_limit:
                    print 'file limit reached:', file_limit
                    break
                # index of corresponding soma trace
                soma_i = (row_i[0], ('soma',0,0), row_i[2], )
                # print soma_i
                # get spike_train from soma trace
                soma_train = df.loc[soma_i]['spike_train'].iloc[0]
                # print soma_train
                # get spike train for current location
                current_train = row['spike_train']
                # print soma_train.shape, current_train.shape
                # get cross correlation
                xcorr_soma = np.correlate(soma_train, current_train, 'full')
                # create time vector 
                xcorr_t = dt*(np.arange(0, len(xcorr_soma))-(len(xcorr_soma)-1)/2)
                # print xcorr_soma.shape
                # set value in 'polarization column'
                print df.loc[row_i].paths
                print df.at[row_i, coladd]
                # print df.at[row_i, coladd].paths
                df.at[row_i, coladd] = xcorr_soma
                df.at[row_i, coladd_t]=xcorr_t
        # if column already exists but some values haven't been processed
        #----------------------------------------------------------------
        # FIXME apply this functionality to other df functions!!!!!!!!!!!!!!!
        elif coladd in df and df[coladd].isnull().values.any():
            # get only unprocessed rows (index must match original df)
            df_temp = df[df[coladd].isnull()].copy()
            # print progress
            print 'unprocessed rows found:', df_temp.shape[0]
            # iterate over rows
            cnt=-1
            for row_i, row in df_temp.iterrows():
                cnt+=1
                # check for file limit
                if file_limit and cnt> file_limit:
                    print 'file limit reached:', file_limit
                    break
                # index of corresponding soma trace
                soma_i = (row_i[0], ('soma',0,0), row_i[2],)
                # print soma_i
                # get spike_train from soma trace
                soma_train = df.loc[soma_i]['spike_train'].iloc[0]
                # print soma_train
                # get spike train for current location
                current_train = row['spike_train']
                # print soma_train.shape, current_train.shape
                # get cross correlation
                xcorr_soma = np.correlate(soma_train, current_train, 'full')
                # create time vector 
                xcorr_t = dt*(np.arange(0, len(xcorr_soma))-(len(xcorr_soma)-1)/2)
                # print xcorr_soma.shape
                # set value in 'polarization column'
                df.at[row_i, coladd] = xcorr_soma
                df.at[row_i, coladd_t]=xcorr_t
        # reset index
        df = df.reset_index()

        return df

    def _get_dcs_effect(self, df, variables, rerun=False, **kwargs):
        ''' get effect of DCS on a given variable at given location
        
        Arguments
        ---------
        ~df : dataframe with olumn 'data_v' containing voltage data in each compartment, and 'dt' column
        ~variable : column name in df

        Return
        ------
        ~df : dataframe with new column as 'original_column_name'+'_dcs_diff'

        Notes
        -----

        '''
        
        # reorganize dataframe
        df = df.set_index(['trial_id','location','field'])
        for variable in variables:
            coladd = variable+'_dcs_diff'
            if coladd not in df or rerun:
                # print progress to terminal
                #-----------------------------
                print 'applying df function:', inspect.stack()[0][3], 'to column:', variable
                # preallocate column with 'object' dtype
                df[coladd]=None
                # iterate over rows
                for row_i, row in df.iterrows():
                    # index of corresponding control trace
                    control_i = (row_i[0], row_i[1], 0.0)
                    # get voltage value from control trace
                    control_val = df.loc[control_i][variable]
                    # get value for current trace
                    val = row[variable]
                    # set value in 'polarization column'
                    df.at[row_i, coladd]=val/control_val

        df = df.reset_index()

        return df

    def _apply_filters(self, df, variables, filters, rerun=False, file_limit=[],**kwargs):
        '''
        Arguments
        ---------
        ~df:
        ~variables:
        ~filters:dictionary of filters [filter_name][b coefficients, a coefficients]
        '''
        # reorganize dataframe
        df = df.set_index(['trial_id','location','field'])
        for variable in variables:
            for filt_key, filt in filters.iteritems():
                coladd = variable+'_filt_'+filt_key
                coladd_hilb = coladd+'_hilb'
                if coladd not in df or rerun:
                    # print progress to terminal
                    #-----------------------------
                    print 'applying filter ', filt_key, ' to ', variable
                    # preallocate column with 'object' dtype
                    df[coladd]=None
                    df[coladd_hilb]=None
                    cnt=0
                    # iterate over rows
                    for row_i, row in df.iterrows():
                        cnt+=1
                        if file_limit and cnt>file_limit:
                            continue
                        # print row_i
                        raw_data = row[variable]
                        # filtered data
                        # filtered_data = signal.filtfilt(filt[0], filt[1], raw_data, method='pad')
                        filtered_data = signal.lfilter(filt[0], filt[1], raw_data,)
                        # 
                        filtered_data_hilb = np.abs(signal.hilbert(filtered_data - np.mean(filtered_data)))
                        df.at[row_i, coladd]=filtered_data
                        df.at[row_i, coladd_hilb]=filtered_data_hilb


        df = df.reset_index()

        return df

    def _apply_filters_cascade(self, df, variables, filters, rerun=False, file_limit=[], **kwargs):
        '''
        Arguments
        ---------
        ~df:
        ~variables:
        ~filters:dictionary of filters in series [filter_name][subfilter_name][b coefficients, a coefficients]
        '''
        # print df.index.names
        # reorganize dataframe
        df = df.set_index(['trial_id','location','field'])
        # iterate over variables
        for variable in variables:
            # iterate over filters
            for filt_key, filt in filters.iteritems():
                # column name
                coladd = variable+'_filt_'+filt_key
                coladd_hilb = coladd+'_hilb'
                # run column?
                if coladd not in df or rerun:
                    # print progress to terminal
                    #-----------------------------
                    print 'applying filter ', filt_key, ' to ', variable
                    # preallocate column with 'object' dtype
                    df[coladd]=None    
                    df[coladd_hilb]=None   
                    cnt=0                 
                    # iterate over rows
                    for row_i, row in df.iterrows():
                        cnt+=1
                        if file_limit and cnt>file_limit:
                            continue
                        # get raw data
                        raw_data = row[variable]
                        # copy to be iteratively filtered
                        filtered_data=copy.deepcopy(raw_data)
                        # iterate over subfilters
                        for subfilt_key, subfilt in filt.iteritems():
                            # print progress to terminal
                            #-----------------------------
                            print 'applying subfilter ', subfilt_key, ' to ', variable
                            # update filtered data
                            # filtered_data = signal.filtfilt(subfilt[0], subfilt[1], filtered_data, method='pad')
                            filtered_data = signal.lfilter(subfilt[0], subfilt[1], filtered_data, )
                        filtered_data_hilb = np.abs(signal.hilbert(filtered_data - np.mean(filtered_data)))
                        # add to dataframe
                        df.at[row_i, coladd]=filtered_data
                        df.at[row_i, coladd_hilb]=filtered_data_hilb

        df = df.reset_index()

        return df

class ApplyDF:
    ''' functions to be passed to df to create new columns from existing data
    '''
    def _to_string(self, df, colnames=[], suffix='_str'):
        ''' convert all cells in column to str
        '''
        for i, colname in enumerate(colnames):
            colname_new = colname+suffix
            df[colname_new] = df[colname].apply(str)
        return df

    def _normalize_column(self, df, colnames=[], colnorm=[],):
        '''
        '''
        # get list of columns to be normalized
        def _norm(row, colnames=colnames, colnorm=colnorm):
            ''' normalize data in colnames columns based on data in colnorm column
            '''
            col_add = '_norm_'+colnorm
            norm = row[colnorm]
            for col in colnames:
                new_data =np.array(row[col])/norm
                row[col+col_add] = new_data
            return row

        kwargs = {'colnames':colnames, 'colnorm':colnorm}
        df = df.apply(_norm, axis=1, **kwargs)
        return df

    def _get_io(self, df, max_idx_col=[], colnames=[], colkeys=[], colkeys_exclude=[]):
        ''' get data for specific pulse from induction data
        '''
        print 'getting io'
        # print colnames
        colnames=copy.copy(colnames)

        # get list of columns to be normalized
        def _get_io(row, max_idx_col=max_idx_col, colnames=colnames, colkeys=colkeys, colkeys_exclude=colkeys_exclude):
            ''' normalize data in colnames columns based on data in colnorm column
            '''
            if len(colkeys)>0:
                columns = row.index.values
                for col in columns:
                    if all([colkey in col for colkey in colkeys]) and not any([temp in col for temp in colkeys_exclude]):
                        
                        colnames.append(col)

            coladd_io_norm = '_io_norm'
            coladd_io_max = '_io_max'
            for col in colnames:
                if not np.isnan(row[max_idx_col]):
                    max_idx = int(row[max_idx_col])
                    # print max_idx
                    max_val = row[col][max_idx]
                    io_vals = row[col][:max_idx+1]/max_val
                    row[col+coladd_io_norm] = io_vals
                    row[col+coladd_io_max] = max_val
                else:
                    row[col+coladd_io_norm]=np.nan
                    row[col+coladd_io_max]=np.nan

            return row

        kwargs = {'colnames':colnames, 'colkeys':colkeys, 'colkeys_exclude':colkeys_exclude, 'max_idx_col':max_idx_col}
        df = df.apply(_get_io, axis=1, **kwargs)
        return df

    def _get_pulse(self, df, colnames=[], colkeys=[], colkeys_exclude=[]):
        ''' get data for specific pulse from induction data
        '''
        # get list of columns to be normalized
        def _get_pulse(row, colnames=colnames, colkeys=colkeys, colkeys_exclude=colkeys_exclude):
            ''' normalize data in colnames columns based on data in colnorm column
            '''
            if len(colkeys)>0:
                columns = row.index.values
                for col in columns:
                    if all([colkey in col for colkey in colkeys]) and not any([colkey in col for colkey in colkeys_exclude]):
                        colnames.append(col)

            coladds = ['_first', '_second', '_third', '_fourth']
            induction_protocol = row.induction_pattern_0
            if induction_protocol=='TBS':
                pulses=4
            elif induction_protocol=='weak5Hz':
                pulses=1
            elif induction_protocol=='nostim':
                pulses=0

            for col in colnames:
                for coladd_i, coladd in enumerate(coladds):
                    # handle nans
                    if type(row[col])==list or type(row[col])==np.ndarray and any(np.isnan(row[col])):
                            new_data=np.nan
                    elif type(row[col])==float and np.isnan(row[col]):
                        new_data=np.nan
                    elif type(row[col])==float:
                        new_data=row[col]
                    else:
                        if pulses==0:
                            mask=[]
                        elif pulses==1:
                            mask = range(0, len(row[col]), pulses)
                        else:
                            mask = range(coladd_i, len(row[col]), pulses)
                    
                        new_data = np.array(row[col])[mask]
                    row[col+coladd] = new_data
            return row

        kwargs = {'colnames':colnames, 'colkeys':colkeys, 'colkeys_exclude':colkeys_exclude}
        df = df.apply(_get_pulse, axis=1, **kwargs)
        return df

    def _mean(self, df, colnames=[], colkeys=[], colkeys_exclude=[]):
        ''' get data for specific pulse from induction data
        '''
        # get list of columns to be normalized
        def _mean(row, colnames=colnames, colkeys=colkeys, colkeys_exclude=colkeys_exclude):
            ''' normalize data in colnames columns based on data in colnorm column
            '''
            if len(colkeys)>0:
                columns = row.index.values
                for col in columns:
                    if all([colkey in col for colkey in colkeys]) and not any([colkey in col for colkey in colkeys_exclude]):
                        colnames.append(col)

            coladd = '_mean'
            for col in colnames:
                row[col+coladd] = np.mean(np.array(row[col]))
            return row

        kwargs = {'colnames':colnames, 'colkeys':colkeys, 'colkeys_exclude':colkeys_exclude}
        df = df.apply(_mean, axis=1, **kwargs)
        return df

    def _last_burst(self, df, colnames=[], colkeys=[], colkeys_exclude=[]):
        ''' get data for specific pulse from induction data
        '''
        print 'getting last burst'
        # get list of columns to be normalized
        def _last_burst(row, colnames=colnames, colkeys=colkeys, colkeys_exclude=colkeys_exclude):
            ''' normalize data in colnames columns based on data in colnorm column
            '''
            colnames=copy.copy(colnames)
            if len(colkeys)>0:
                columns = row.index.values
                for col in columns:
                    if all([colkey in col for colkey in colkeys]) and not any([temp in col for temp in colkeys_exclude]):
                        
                        colnames.append(col)

            coladd = '_lastburst'
            for col in colnames:
                # print row[col]
                if type(row[col])==float and np.isnan(row[col]):
                    row[col+coladd]=np.nan
                elif type(row[col])==list and any(np.isnan(row[col])):
                    row[col+coladd]=np.nan
                else:
                    # print col
                    row[col+coladd] = np.array(row[col])[-1]
            return row

        kwargs = {'colnames':colnames, 'colkeys':colkeys, 'colkeys_exclude':colkeys_exclude}
        df = df.apply(_last_burst, axis=1, **kwargs)
        return df
    
    def _first_burst(self, df, colnames=[], colkeys=[], colkeys_exclude=[]):
        ''' get data for specific pulse from induction data
        '''
        print 'getting first burst'
        # print colnames
        colnames=copy.copy(colnames)

        # get list of columns to be normalized
        def _first_burst(row, colnames=colnames, colkeys=colkeys, colkeys_exclude=colkeys_exclude):
            ''' normalize data in colnames columns based on data in colnorm column
            '''
            if len(colkeys)>0:
                columns = row.index.values
                for col in columns:
                    if all([colkey in col for colkey in colkeys]) and not any([temp in col for temp in colkeys_exclude]):
                        
                        colnames.append(col)
            # print colnames

            coladd = '_firstburst'
            for col in colnames:
                # print row[col]
                if type(row[col])==float and np.isnan(row[col]):
                    row[col+coladd]=np.nan
                elif type(row[col])==list and any(np.isnan(row[col])):
                    row[col+coladd]=np.nan
                else:
                    # print col
                    row[col+coladd] = np.array(row[col])[0]
            return row

        kwargs = {'colnames':colnames, 'colkeys':colkeys, 'colkeys_exclude':colkeys_exclude}
        df = df.apply(_first_burst, axis=1, **kwargs)
        return df

    def _fft(self, df, colnames=[], colkeys=[], colkeys_exclude=[], slice_i=slice(None), unwrap=False, ):
        ''' get data for specific pulse from induction data
        '''
        print 'getting fft'
        # get list of columns to be normalized
        def _fft(row, colnames=colnames, colkeys=colkeys, colkeys_exclude=colkeys_exclude, slice_i=slice_i, unwrap=False):
            ''' normalize data in colnames columns based on data in colnorm column
            '''
            if len(colkeys)>0:
                columns = row.index.values
                for col in columns:
                    if all([colkey in col for colkey in colkeys]) and not any([colkey in col for colkey in colkeys_exclude]):
                        colnames.append(col)

            coladd_abs = '_fft_abs'
            coladd_angle='_fft_angle'
            
            for col in colnames:
                print np.array(row[col]).squeeze()
                # print np.array(row[col]).squeeze()==np.nan
                # print type(np.isnan(np.array(row[col]).squeeze()))
                if type(np.isnan(np.array(row[col]).squeeze()))==np.bool_ and np.isnan(np.array(row[col]).squeeze()):
                    row[col+coladd_abs]=np.nan
                    row[col+coladd_angle]=np.nan
                elif any(np.isnan(np.array(row[col]).squeeze())):
                    row[col+coladd_abs]=np.nan
                    row[col+coladd_angle]=np.nan

                else:
                    slice_is = slice_i.indices(len(np.array(row[col])))
                    if len(range(*slice_is))==0:
                        slicer = slice(None)
                    else:
                        slicer=copy.copy(slice_i)
                    row[col+coladd_abs] = np.abs(np.fft.fft(np.array(row[col])[slicer], n=60))
                    if unwrap:
                        row[col+coladd_angle] = np.unwrap(np.angle(np.fft.fft(np.array(row[col])[slicer], n=60)))
                    else:
                        row[col+coladd_angle] = np.angle(np.fft.fft(np.array(row[col])[slicer], n=60))

            return row

        kwargs = {'colnames':colnames, 'colkeys':colkeys, 'colkeys_exclude':colkeys_exclude, 'slice_i':slice_i, 'unwrap':unwrap}
        df = df.apply(_fft, axis=1, **kwargs)
        return df

    def _get_nth(self, df, colnames=[], colkeys=[], colkeys_exclude=[], n=0):
        ''' get data for specific pulse from induction data
        '''
        print 'getting fft'
        # get list of columns to be normalized
        def _get_nth(row, colnames=colnames, colkeys=colkeys, colkeys_exclude=colkeys_exclude, n=n):
            ''' normalize data in colnames columns based on data in colnorm column
            '''
            if len(colkeys)>0:
                columns = row.index.values
                for col in columns:
                    if all([colkey in col for colkey in colkeys]) and not any([colkey in col for colkey in colkeys_exclude]):
                        colnames.append(col)

            coladd = '_'+str(n)+'th'
            
            for col in colnames:
                print np.array(row[col]).squeeze()
                # print np.array(row[col]).squeeze()==np.nan
                # print type(np.isnan(np.array(row[col]).squeeze()))
                if type(np.isnan(np.array(row[col]).squeeze()))==np.bool_ and np.isnan(np.array(row[col]).squeeze()):
                    row[col+coladd]=np.nan
                    row[col+coladd]=np.nan
                elif any(np.isnan(np.array(row[col]).squeeze())):
                    row[col+coladd]=np.nan
                    row[col+coladd]=np.nan

                else:
                    if len(np.array(row[col]))>n:
                        row[col+coladd] = np.array(row[col])[n]
                    # slice_is = slice_i.indices(len(np.array(row[col])))
                    # if len(range(*slice_is))==0:
                    #     slicer = slice(None)
                    # else:
                    #     slicer=copy.copy(slice_i)
                    # row[col+coladd_abs] = np.abs(np.fft.fft(np.array(row[col])[slicer], n=60))
                    # if unwrap:
                    #     row[col+coladd_angle] = np.unwrap(np.angle(np.fft.fft(np.array(row[col])[slicer], n=60)))
                    # else:
                    #     row[col+coladd_angle] = np.angle(np.fft.fft(np.array(row[col])[slicer], n=60))

            return row

        kwargs = {'colnames':colnames, 'colkeys':colkeys, 'colkeys_exclude':colkeys_exclude, 'n':n}
        df = df.apply(_get_nth, axis=1, **kwargs)
        return df
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
                        print df.index.values
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
                    print data_array[0].shape
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

    def _truncate_array(self, array):
        '''
        '''
        #FIXME this is a hack
        if type(array[0])==np.ndarray:
        # if len(array.shape)==1:
            sizes = [row.shape[0] for row in array]
            minsize = min(sizes)
            print minsize, array.shape[0]
            new_array = np.zeros((array.shape[0], minsize))
            for row_i, row in enumerate(array):
                new_array[row_i,:]=  row[:minsize]
        else:
            new_array = array
            
        return new_array


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
#############################################################################
#############################################################################
# deprecated
#############################################################################

    # def _initialize_column(df, col, val=None):
    #     '''
    #     '''
    #     if type(col)==str:
    #         col=[col]
    #     for col_name in col:
    #         if col_name not in df:
    #             df[col_name]=val
    #     return df

    # def exp_fit(x, y, p0=[1,-1]):
    #     '''
    #     '''
    #     def exp_func(x, a, b):
    #         return a*np.exp(b*x)

    #     popt, pcov = curve_fit(exp_func, x, y, p0)

    #     return popt

    # def _build_timeseries_anova_df(df, array_column, columns):
    #     '''
    #     '''
    #     df_new=pd.DataFrame()
    #     for row_i, row in df.iterrows():
    #         timeseries = row[array_column]
    #         print type(timeseries)
    #         if not type(timeseries)==float or (type(timeseries)==float and not np.isnan(timeseries)):
    #             if type(timeseries)==float:
    #                 timeseries = np.array([timeseries])
    #             df_temp = pd.DataFrame(timeseries,)
    #             df_temp.columns = [array_column]
    #             print df_temp.columns
    #             df_temp[array_column+'_i']=df_temp.index.values
    #             for column in columns:
    #                 df_temp[column]=row[column]
    #             df_new = df_new.append(df_temp)

    #     return df_new

    # def _build_ttest_dict(df, figdf, variable, conditions):
    #     '''
    #     '''
    #     stats_df = _build_timeseries_anova_df(df=df, array_column=variable, columns=conditions)
    #     stats_df = stats_df.set_index(conditions)

    #     traces = figdf.index.get_level_values('trace')
    #     trace_combos = list(itertools.product(traces,traces))
    #     stats_dict = {}
    #     print sorted(stats_df.index.unique())
    #     for combo in trace_combos:
    #         # print combo[0]
    #         if combo[0] in stats_df.index.unique() and combo[1] in stats_df.index.unique():
    #             print combo

    #             df_0 = stats_df.loc[combo[0]]
    #             df_1 = stats_df.loc[combo[1]]
    #             if combo[0] not in stats_dict:

    #                 stats_dict[combo[0]]={
    #                     'n':df_0.shape[0],
    #                     'std':np.std(df_0[variable]),
    #                     'sem':stats.sem(df_0[variable]),
    #                     'mean':np.mean(df_0[variable]),
    #                 }
    #             if combo[1] not in stats_dict[combo[0]]:
    #                 stats_dict[combo[0]][combo[1]]= stats.ttest_ind(df_0[variable],df_1[variable]).pvalue

    #     return stats_dict

    # def _add_escape(constraints):
    #     '''add escape quotes to query string based on constraint dictionary
    #     ==Args==
    #     -constraints : dictionary specifying constraints to be passed to dataframe as query
    #             ~can be organized as {constraint column in dataframe}{(condition values from conditions list)}[constraint logical, constraint value] for constraints that are specific to different experimental conditions
    #             ~can be organized as {constraint column}[constraint logical, constraint value] if the constraint is to be aplied to all conditions
    #     ==Out==
    #     -constraints : updated constraints dictionary, any constraint values that are strings have escape quotes added for use in query in string

    #     ==Updates==
    #     ==Comments==
    #     '''
    #     # add escape quotes to df entries that are strings. needed to search for strings with query function
    #     # contstraints_all
    #     for k0, v0 in constraints.iteritems():
    #         if type(v0)==dict:
    #             for k1, v1 in v0.iteritems():
    #                 if type(v1[1])==str and '\'' not in v1[1]:
    #                     constraints[k0][k1][1] = '\'' + v1[1] + '\''
    #         else:
    #             if type(v0[1])==str and '\'' not in v0[1]:
    #                 constraints[k0][1] = '\'' + v0[1] + '\''

    #     return constraints

    # def _build_constraint_query(constraints_spec, constraints_all, group_key):
    #     '''
    #     ==Args==
    #     ==Out==
    #     ==Updates==
    #     ==Comments==
    #     ''' 
    #     print 'building constraint query to sort df for:', group_key
    #     strings = []
    #     for k0, v0 in constraints_spec.iteritems():
    #         for k1, v1 in v0.iteritems():
    #             if type(k1)==str:
    #                 k1=(k1,)

    #             constraint_vals = set(k1)
    #             if len(set(constraint_vals) - set(group_key))==0:

    #                 current_string = '{}{}{}'.format(k0, v1[0], v1[1])
    #                 strings.append(current_string)
    #     constraint_query = '&'.join(strings)


    #     if len(constraints_all.keys())!=0:
    #         if any(val not in constraints_spec.keys() for val in constraints_all.keys() ):
    #             constraint_query_all =  '&'.join([
    #                 '{}{}{}'.format(k,v[0],v[1])
    #                 for k, v in constraints_all.iteritems() if k not in constraints_spec
    #             ])

    #             constraint_query +='&' + constraint_query_all 


    #     return constraint_query

    # def _align_row_data(self, df, column, align_on='induction_block_0', i_pre=20, i_post=60, include_align=True):
    #         ''' create columns with normalized slopes that are aligned to the induction block and the final amount of ltp
    #         ==Args==
    #         -group_df : group dataframe. must contain columns 'induction_block_0', 'slopes_norm'

    #         ==Out==
    #         -group_df : group dataframe. new columns: 'slopes_norm_aligned_0' and 'ltp_final'

    #         ==Update==
    #         -group_df.slopes_norm_aligned_0 : normalized slopes are aligned to the first induction block, so that each entry is an array of length 80 (the first 20 samples are pre-indction, the last 60 samples are post-induction)
    #         -group_df.ltp_final : mean of the last ten normalized slopes (min 51-60) for each trace

    #         ==Comments==
    #         '''
    #         def _get_aligned_data(row):
    #             '''for each row align slopes data to the induction block
    #             ==Args==
    #             -row : row in a pandas dataframe
    #             ==Out==
    #             -aligned_data :  a list, where the first element is the aligned slopes
    #             ==Comments==
    #             -when fed to the apply function in pandas, a series is returned, where each entry contains the corresponding aligned data
    #             '''
    #             # print 'induction block 0:', row['induction_block_0']
    #             # get indices to align data
    #             indices = range( int( row[ align_on])-i_pre, int(row[align_on])) + range( int( row[align_on]), int( row[ align_on])+ i_post)
    #             # get aligned data
    #             aligned_data = [row[column][indices]]
    #             return aligned_data

    #         aligned_column = group_df.apply(lambda row: _get_ltp_final(row), axis=1)
    #         return aligned_column

    # def _2array(series, remove_nans=True, remove_nans_axis=0, list_index=0, array_funcs=[], array_func_kws=[]):
    #     '''
    #     '''
        
    #     series = series.dropna()
    #     series_list = series.tolist()

    #     if len(series_list)>0 and type(series_list[0])==list:
    #     #     print 'list to array?'
    #         series_list = [item[list_index] for item in series_list if type(item)==list and len(item)>list_index]
    #     #     # print series_list
    #     #     print np.array(np.array(series_list))

    #     array = np.array(series_list).squeeze()

        

    #     # array = np.array(series.tolist()).squeeze()
    #     # print np.isnan(array)
    #     # array = array[~np.isnan(array).any(axis=remove_nans_axis)]
    #     return array

    # def _sortdf(df, conditions, constraints_spec, constraints_all):
    #     '''
    #     '''
    #     print 'sorting dataframe'
    #     # add escape quotes to df entries that are strings. needed to search for strings with query function
    #     # contstraints_all
    #     constraints_spec = _add_escape(constraints_spec)
    #     constraints_all = _add_escape(constraints_all)
    #     # get sorting indices
    #     # groups is a dictionary with keys that correspond to combinations of conditions (e.g. ('anodal', 'TBS', 'nostim')) and values that are the indices that meet those conditions
    #     groups = df.groupby(conditions).groups

    #     # use groupby to sort dataframes by conditions
    #     df_sorted = {}
    #     for group_key, group_index in groups.iteritems():
            
    #         # group_key must be a tuple
    #         if type(group_key)==str:
    #             group_key=(group_key,)

    #         # build constraint query string out of constraint conditions
    #         constraint_query = _build_constraint_query(constraints_spec, constraints_all, group_key)

    #         # sort the data that match the experimental conditions and constraints
    #         df_sorted[group_key] = df.loc[group_index].query(constraint_query)

    #     return df_sorted

    # def _load_group_data( directory='', filename='slopes_df.pkl', df=True):
    #     """ Load group data from folder
        
    #     ===Args===
    #     -directory : directory where group data is stored including /
    #     -filename : file name for group data file, including .pkl
    #                 -file_name cannot contain the string 'data', since this string is used t search for individual data files
    #     -df : boolean, if true group data is assumed to be a pandas dataframe, otherwise it is assumed to be a nested dictionary

    #     ===Out===
    #     -group_data  : if df is True, group_data will be a pandas dataframe.  if no group data is found in the directory an empty dataframe is returned.  if df is False, group_data is returned as a nested dictionary, and if no data is found an empty dictionary is returned

    #     ===Updates===
    #     -none

    #     ===Comments===
    #     """
        
    #     # all files in directory
    #     files = os.listdir(directory)

    #     # if data file already exists
    #     if filename in files:
    #         print 'group data found:', filename

    #         # if stored as dataframe
    #         if df:
    #             # load dataframe
    #             group_data=pd.read_pickle(directory+filename)
    #             print 'group data loaded'
    #         # if stored as dictionary
    #         else:
    #             # load dictionary
    #             with open(directory+filename, 'rb') as pkl_file:
    #                 group_data= pickle.load(pkl_file)
    #             print 'group data loaded'

    #     # otherwise create data structure
    #     else:
    #         print 'no group data found'
            
    #         # if dataframe
    #         if df:
    #             # create empty dataframe
    #             group_data = pd.DataFrame()
    #         # else if dicitonary
    #         else:
    #             group_data= {}

    #     return group_data 

    # def _remove_outliers(time_series, ind_idx, time_window, std_tol, include_ind=False):
    #     '''# FIXME docs, comments
    #     '''
    #     ts = time_series
    #     # preallocate
    #     ts_smooth = copy.copy(ts) 
    #     ind_idx_array = np.array(ind_idx)
    #     # iterate over values
    #     for i, val in enumerate(ts):
    #         # if the current index is an induction, skip it 
    #         if i in ind_idx:
    #             continue

    #         # get indices of inductions that were before the current index
    #         ind_cross = np.where(ind_idx_array<=i)[0]
    #         # if there have not been any inductions yet
    #         if len(ind_cross)==0:
    #             # set most recent induction to 0
    #             ind =[0]
    #             # otherwise, get most recent induction
    #         else:
    #             ind = ind_idx_array[ind_cross]

    #         # if the current index is within a time_window distanc of the most recent induction
    #         if include_ind:
    #             if i-ind[-1]+1 <= time_window:
    #                 vals_i = range(ind[-1]+1,(ind[-1]+1+time_window) )
    #                 current_i = i-ind[-1]+1
    #                 # use forward values to fill out time window
    #                 vals = ts[vals_i]
    #                 # vals = ts[ind[-1]+1:(ind[-1]+1+time_window)]
    #             # otherwise use the last time window values to fidn outliers
    #             else:
    #                 vals_i =range(i-time_window,i)
    #                 # current_i=-1
    #                 vals = ts[vals_i]
    #                 current_i =len(vals)-1
    #         else:
    #             if i-ind[-1] +1 < time_window:
    #                 if ind[-1]+time_window >= len(ts)+1:
    #                     vals_i = range(ind[-1],len(ts))
    #                 else:
    #                     vals_i = range(ind[-1],(ind[-1]+time_window))
    #                 current_i = i-ind[-1]
    #                 # use forward values to fill out time window
    #                 vals=ts[vals_i]
    #                 # vals = ts[ind[-1]:(ind[-1]+time_window)]
    #             # otherwise use the last time window values to fidn outliers
    #             else:
    #                 vals_i = range(i-time_window+1,i+1)
    #                 # current_i = -1
    #                 vals=ts[vals_i]
    #                 current_i =len(vals)-1
    #                 # vals = ts[i-time_window+1:i]

    #         # mean and standard deviation of values
    #         other_vals = [temp for temp_i, temp in enumerate(vals) if temp_i != current_i]
    #         vals_std = np.std(other_vals)
    #         vals_mean = np.mean(other_vals)

    #         # if value is outside the tolerance set by std_tol, set value to mean of all values in time_window 
    #         if abs(val-vals_mean) < std_tol*vals_std:
    #             ts_smooth[i] = val
    #         else: 
    #             ts_smooth[i] = vals_mean

    #     return ts_smooth

    # def _build_conditions_dict(preprocessed, path):
    #     ''' FIXME add doc
    #     '''
    #     pre = preprocessed
    #     # dictionary for current data file to be converted to dataframe
    #     current_dict = {}
    #     # df = pd.DataFrame()
    #     # iterate over pathways
    #     # for path in pre['path_blocks']:

    #     # FIXME how to handle more than 2 pathways
    #     # get name of other pathway
    #     n_paths = len(pre['path_blocks'].keys())
    #     try:
    #         path_other = [temp for temp in pre['path_blocks'] if temp != path][0]
    #     except IndexError:
    #         path_other = 'None'


    #     # iterate over inductions
    #     for ind_i, ind in enumerate(pre['induction_info']):

    #         # store induction info as entries in current_dict. the induction number is indicated at the end of each dictionary key, e.g. 'induction_0'

    #         # FIXME indices need to be realigned to each pathway after the pathways are separated.  use ind_idx from _get_slopes_probe
    #         # store the full induction info dictionary
    #         current_dict['induction_'+str(ind_i)]=[ind]

    #         # induction pattern, e.g. TBS, weak5Hz, nostim
    #         current_dict['induction_pattern_'+str(ind_i)] = ind[path]['protocol']
    #         try:
    #             current_dict['induction_pattern_other_'+str(ind_i)] = ind[path_other]['protocol']
    #         except KeyError:
    #             current_dict['induction_pattern_other_'+str(ind_i)] = 'None'

            
    #         # location stimulating electrode, e.g. apical, basal
    #         current_dict['induction_location_'+str(ind_i)] = ind[path]['location']
    #         try:
    #             current_dict['induction_location_other_'+str(ind_i)] = ind[path_other]['location']
    #         except KeyError:
    #             current_dict['induction_location_other_'+str(ind_i)] = 'None'


    #         # induction block
    #         current_dict['induction_block_'+str(ind_i)] = pre['ind_idx'][path][ind_i]
    #         # ind['induction_block']

    #         # electric field magnitude (in V/m)
    #         current_dict['field_mag_'+str(ind_i)] = ind['field_magnitude']

    #         # electric field polarity (anodal, cathodal, control)
    #         current_dict['field_polarity_'+str(ind_i)] = ind['polarity']


    #     current_dict['comments'] = [pre['comment_dict']]
    #     current_dict['induction']=[pre['induction_info']]
    #     current_dict['path']=[path]
    #     current_dict['filename']=pre['slice_info']['filename']
    #     current_dict['name']=pre['slice_info']['name']
    #     current_dict['date']=pre['slice_info']['date']
    #     current_dict['path_other'] = path_other
    #     current_dict['age'] = [pre['slice_info']['age']]
    #     current_dict['height'] = [pre['slice_info']['height']]
    #     current_dict['hemi'] = [pre['slice_info']['hemi']]
    #     current_dict['baseline_percent'] = [pre['slice_info']['baseline_percent'][path]]
    #     current_dict['baseline_max_idx'] = [pre['slice_info']['baseline_max_idx'][path]]

    #         # current_df = pd.DataFrame(current_dict)
    #         # df = df.append(current_df, ignore_index=True)

    #     return current_dict

    # def _process_new_data_df(group_df, preprocessed_directory='Preprocessed Data/', search_string='.pkl', functions=[], kwlist=[], rerun=[], keep=[], file_limit=[], ):
    #     ''' process new data and add to group

    #     ==Args==
    #     -group_df : pandas dataframe containing group data
    #     -directory : relative directory containing unprocessed data
    #     -search_string :  string to identify unprocessed data, typically ".mat"
    #     -variables : variables to be derived from unprocessed data (e.g. slopes or slopes_norm)

    #     ==Out==
    #     -group_df : pandas dataframe containing group data

    #     ==Update==
        
    #     ==Comments==
    #     '''

    #     # get list of all preprocessed data files
    #     data_files = glob.glob(preprocessed_directory+'*'+search_string+'*')
    #     # remove directory and file extension
    #     data_filenames = [file.split('\\')[-1].split('.')[0] for file in data_files]

    #     # get list of processed data files (already added to group dataframe)
    #     if group_df.empty:
    #         processed_data_filenames = []
    #     else:
    #         processed_data_filenames = group_df.name

    #     # get list of new data files (haven't been added to dataframe)
    #     new_data_filenames = list(set(data_filenames)-set(processed_data_filenames))
    #     new_data_files = [data_files[file_i] for file_i, file in enumerate(data_filenames) if file in new_data_filenames]

    #     # if any functions specified in rerun, set iterfiles to all data files
    #     if rerun:
    #         iterfiles = data_files
    #     # if there are no functions to rerun, only iterate through new files
    #     else:
    #         iterfiles = new_data_files

    #     print 'total data files:', len(data_files) 
    #     print 'new data files:', len(new_data_files)
    #     if rerun:
    #         print 'functions to rerun over all slices:', rerun

    #     # iterate over new files and update group_data structure
    #     #`````````````````````````````````````````````````````````````````
    #     print 'updating group data structure'

    #     # temporary df for storing rerun function outputs
    #     df_update = pd.DataFrame()

    #     # suffixes for merging df_update with group_df
    #     old_suffix = '_old'
    #     rerun_suffix = '_rerun'

    #     # iterate over data files
    #     for file_i, file in enumerate(iterfiles):

    #         if file_limit and file_i> file_limit:
    #             continue
    #         else:


    #             # get file name and date
    #             name = file.split('\\')[-1].split('.')[0]
    #             date = int(name[:8])

    #             # # limit number of files processed at a time
    #             # if file_limit and file_i>file_limit:
    #             #     continue
    #             # else:

    #             # load data file
    #             with open(file, 'rb') as pkl_file:

    #                 pre = pickle.load(pkl_file)

    #                 print name

    #             # for each function create df and add to list (either new or rerun)
    #             new_dfs = []
    #             rerun_dfs = []
    #             # iterate over functions
    #             for func_i, func in enumerate(functions):
    #                 print func


    #                 # if new data file add to new list
    #                 if file in new_data_files:

    #                     print 'new file found:', name
    #                     print kwlist[func_i]
    #                     df_temp = func(pre, df=[], keep=[], **kwlist[func_i])
    #                     # # print 'test'
    #                     new_dfs.append(df_temp)

    #                 # if rerun add to old list
    #                 elif func in rerun:
    #                     # print name
    #                     # get corresponding row from group_df
    #                     df_locs = group_df[group_df.name==name].index.values
    #                     # get df from function
    #                     df_temp = func(pre, keep=keep, df=group_df.loc[df_locs], **kwlist[func_i])
    #                     # print 'rerun keys:',(df_temp.keys())
    #                     # print func
    #                     rerun_dfs.append(df_temp)

    #             # if there are new df's, add to bottom of group_df
    #             if new_dfs:

    #                 # for temp in new_dfs:
    #                 #     print sorted(temp.keys())
    #                 # merge df's across the different functions, drop columns that are common between functions
    #                 # new_df = reduce(lambda left,right: pd.merge(left, right, on=right.columns.intersection(left.columns), how='inner'), new_dfs)
    #                 # print new_dfs
    #                 new_df = reduce(
    #                     lambda left,right: pd.merge(
    #                         left, 
    #                         right[np.append(right.columns.difference(left.columns).values, ['name', 'path'])], 
    #                         on=['name', 'path'], 
    #                         how='inner', ), 
    #                     new_dfs
    #                     )



    #                 # print new_df
    #                 # new_df = reduce(lambda left,right: pd.merge(left, right[right.columns.difference(left.columns)], on=['name', 'path'], how='inner' , ), new_dfs)# validate='one_to_one'
    #                 # append to bottom of group_df
    #                 # group_df = group_df.append(new_df, ignore_index=True)
    #                 df_update = df_update.append(new_df, ignore_index=True)

    #                 # print df_update.induction_block_0

    #             if rerun_dfs:
    #                 # for temp in rerun_dfs:
    #                 #     print temp['name'], temp['path']

    #                 # print np.append(rerun_dfs[0].columns.difference(rerun_dfs[1].columns).values, ['name','path'])

    #                 # rerun_df = reduce(lambda left,right: pd.merge(left, right, on=right.columns.intersection(left.columns), how='inner'), rerun_dfs)

    #                 # print rerun_df
    #                 # merge df's across the different functions
    #                 rerun_df = reduce(
    #                     lambda left,right: pd.merge(
    #                         left, 
    #                         right[np.append(right.columns.difference(left.columns).values, ['name', 'path'])], 
    #                         on=['name', 'path'], 
    #                         how='inner', ), 
    #                     rerun_dfs
    #                     ) # validate='one_to_one'

    #                 df_update = df_update.append(rerun_df, ignore_index=True)

    #     # update corresponding rows in group_df
    #     # for key in df_update:
    #     #     if 'hilbert' in key:
    #     #         pdb.set_trace()
    #     #         print df_update[key]

    #     if group_df.empty:
    #         group_df=df_update
    #     elif not df_update.empty:
    #         group_df = group_df.merge(df_update, on=['name', 'path'], how='outer', suffixes=[old_suffix, rerun_suffix])



    #     # print group_df[group_df.date_rerun==20181127].induction_block_0_old

    #     columns = group_df.columns.values

    #     drop_columns = []
    #     for column in columns:
    #         # if 'data_ind_hilbert_smooth'in column:
    #             # print column, group_df[column]
    #         if old_suffix in column:
    #             column_original = column.split(old_suffix)[0]
    #             column_new = column_original+rerun_suffix
    #             # print column_new
    #             # if not a keep column, use new to update old, then drop new
    #             if not any([temp for temp in keep if temp in column]):
    #                 # FIXME, use overwrite kwarg for update
    #                 group_df[column].update(group_df[column_new])
    #                 drop_columns.append(column_new)
    #             # if a keep column, use old to update new, drop old
    #             else:
    #                 group_df[column_new].update(group_df[column])
    #                 drop_columns.append(column)

    #         # if 'data_ind_hilbert_smooth'in column:
    #         #     print column, group_df[column]




    #     # drop_columns = [column for column in columns if rerun_suffix in column]

    #     # build list of columns to drop
    #     # drop_columns=[]
    #     # for column in columns:
    #     #     if len(keep)==0:
    #     #         if old_suffix in column:
    #     #             drop_columns.append(column)
    #     #     else:
    #     #         if old_suffix in column and not any([temp for temp in keep if temp in column]):
    #     #             drop_columns.append(column)

    #     #         elif rerun_suffix in column and any([temp for temp in keep if temp in column]):
    #     #             drop_columns.append(column)

    #     # print drop_columns
    #     group_df.drop(drop_columns, axis=1, inplace=True)

    #     # remove suffixes
    #     rename_dict={}
    #     for column in group_df.columns.values:

    #         if old_suffix in column:
    #             rename_dict[column] = column.split(old_suffix)[0]
    #         elif rerun_suffix in column:
    #             rename_dict[column] = column.split(rerun_suffix)[0]

    #     group_df.rename(columns=rename_dict, inplace=True)
    #     # pdb.set_trace()
    #     # print group_df.data_induction_data_filt_iir_band_300_1000_hilbert_sortby_burst_apical[:10]
    #     # print group_df.data_induction_data_filt_iir_band_300_1000_hilbert_sortby_burst_apical_old[:10]

    #     return group_df


    #         # iterate over columns in group_df
    #         # if column in keep, use old column values and drop suffix
    #         # otherwise keep new columnn values and drop suffix


    #         # merge all function dfs on filename and path
    #         # update appropriate rows and columns of group_df with function dfs, use merge? find indices and use update/combine?

    # def _postdict2df(postdict, pre):
    #     ''' convert dictionary of postprocessing variables into dataframe to be merged with group data
    #     '''
    #     df = pd.DataFrame()
    #     for path in pre['path_blocks']:

    #         current_dict = _build_conditions_dict(pre, path)

    #         for key in postdict:
    #             # print key, postdict[key][path]
    #             if type(postdict[key][path])==float and postdict[key][path] == np.nan:
    #                 current_dict[key]=postdict[key][path]
    #             else:
    #                 current_dict[key]=[postdict[key][path]]

    #         for key in current_dict:
    #             if type(current_dict[key])==list and len(current_dict[key])>1:
    #                 current_dict[key] = [current_dict[key]]

    #         # convert to dataframe    
    #         current_df = pd.DataFrame(current_dict)

    #         # add to group data
    #         if df.empty:
    #             df=current_df
    #         else:
    #             df = df.append(current_df, ignore_index=True)

    #     return df

    # def _cca(df, conditions, constraints_spec, constraints_all, x_variables, y_variables, ):
    #     '''
    #     '''
    #     df_sorted = _sortdf(df=df, conditions=conditions, constraints_spec=constraints_spec, constraints_all=constraints_all)
    #     cca = {}
    #     cca_result={}
    #     for group_key, group_df in df_sorted.iteritems():
    #         # print any(pd.isna(group_df[x_variables]).any().values)
    #         if not any(pd.isna(group_df[x_variables]).any().values) and not any(pd.isna(group_df[y_variables]).any().values) :
    #             cca[group_key] = CCA()

    #             x = group_df[x_variables].values
    #             x_norm = x - np.mean(x, axis=0)
    #             x_norm = x_norm/np.std(x_norm)
    #             y = group_df[y_variables].values
    #             y_norm = y - np.mean(y, axis=0)
    #             y_norm = y_norm/np.std(y_norm)
    #             polarity = group_df.field_polarity_0_slopes.values
    #             x_cca, y_cca= cca[group_key].fit(x_norm, y_norm).transform(x_norm, y_norm)
    #             x_weights = cca[group_key].x_weights_
    #             y_weights = cca[group_key].y_weights_


    #             cca_result[group_key]={
    #             'x':x,
    #             'y':y,
    #             'x_norm':x_norm,
    #             'y_norm':y_norm,
    #             'x_variables':x_variables,
    #             'y_variables':y_variables,
    #             'cca_obj':cca[group_key],
    #             'x_cca':x_cca,
    #             'y_cca':y_cca,
    #             'x_weights':x_weights,
    #             'y_weights':y_weights,
    #             'polarity':polarity
    #             }
    #     return cca_result

    # def _plot_scalebar(axes, yscale, xscale, origin, width):
    #     '''
    #     '''
    #     # plot y
    #     axes.plot([origin[0], origin[0]], [origin[1], origin[1]+yscale], linewidth=width, color='k')
    #     axes.plot([origin[0], origin[0]+xscale], [origin[1], origin[1]], linewidth=width, color='k')

    #     return axes

    # def _plot_vtrace(df_sorted, conditions, figures, variable, colors, markers, titles, mean=True, all_paths=True, **kwargs):
    #     '''FIXME add docs
    #     '''
    #     # FIXME add kwargs to alter figure details
    #     # create figure groupings (all conditions that will go on the same figure)
    #     fig = {}
    #     ax={}
    #     ylim=[]
    #     # iterate over figures
    #     for figure_key in figures:
    #         if type(figure_key)==str:
    #             figure_key=(figure_key,)
    #         # create figure objects
    #         fig[figure_key], ax[figure_key] = plt.subplots()
    #         # plot object for each trace
    #         trace = OrderedDict()
    #         trace_final = OrderedDict()
    #         # n for each condition
    #         n_trial=OrderedDict()
    #         # pvalue for ttest comparing to corresponding control condition
    #         pval=OrderedDict()
    #         # label for plot legend
    #         label=OrderedDict()
    #         # iterate over each condition
    #         for group_key, df in df_sorted.iteritems():
    #             # print figure_key
    #             # print group_key
    #             # print len(set(group_key) - set(figure_key))
    #             # if type(group_key)==tuple:
    #             #     group_key = list(group_key)
    #             # elif type(group_key)==str:
    #             #     group_key = [group_key]
    #             # print len(set(group_key) - set(figure_key))
    #             # if plotting all pathways on the same figure
    #             if all_paths and not all(temp in group_key for temp in figure_key) and 'figures_any' not in kwargs:
    #                 # print 'huuuuuh'
    #                 # if condition not met continue
    #                 continue
    #             # if only plotting specific pathways
    #             elif not all_paths and not group_key[1]==figure_key[0] and group_key[2]==figure_key[1] and 'figures_any' not in kwargs:
    #                 # if condition not met, continue
    #                 # print 'huhhhhhhh'
    #                 continue

    #             if 'figures_any' in kwargs and len(set(group_key) - set(figure_key)) != 0:

    #                 # print figure_key, group_key,'huh'


    #             # any(set(group_key)&set(figure_key)):
    #                 continue
    #             # if conditions are satisfied, plot figure
    #             else:
    #                 # get color
    #                 color = colors[group_key[0]]
    #                 if 'colors_post' in kwargs:
    #                     print 'color_found'
    #                     color_post = kwargs['colors_post'][group_key[0]]
    #                 else:
    #                     color_post=color    

    #                 color_final = tuple([0.4*val for val in color])
    #                 # list of group key values and list of marker keys,use the first match that you find
    #                 if 'all' in markers.keys():
    #                     marker=markers['all']
    #                 elif len(list(set(markers.keys()) & set(group_key)))>0:
    #                     marker_key = list(set(markers.keys()) & set(group_key))[0]
    #                     marker = markers[marker_key]
    #                 else:
    #                     marker='.'


    #                 # print set(group_key)
    #                 # print set(markers.keys())
    #                 # print list(set(markers.keys()) & set(group_key))
    #                 # marker_key = list(set(markers.keys()) & set(group_key))[0]
    #                 # # get marker
    #                 # marker = markers[marker_key]
    #                 # convert data to 2d array (slices x time)
    #                 data_array = _2array(df[variable])
    #                 print data_array.shape
    #                 # check that data is 2d (if not, it is due to a lack of data for this condition)
    #                 if len(data_array.shape)>1:
    #                     # mean across slices
    #                     data_baseline = np.mean(data_array[:,:,:20], axis=2).squeeze()
    #                     data_final = np.mean(data_array[:,:,10:], axis=2).squeeze()

    #                     data_baseline_mean = np.mean(data_baseline, axis=0)
    #                     data_final_mean = np.mean(data_final, axis=0)

    #                     #std across slices
    #                     data_baseline_std = np.std(data_baseline, axis=0)
    #                     data_final_std = np.std(data_final, axis=0)
    #                     # sem across slices
    #                     data_baseline_sem = stats.sem(data_baseline, axis=0)
    #                     data_final_sem = stats.sem(data_final, axis=0)
    #                     # time vector
    #                     t = np.arange(len(data_baseline_mean))*0.1
    #                     print t.shape, data_baseline_mean.shape
    #                     # if plotting mean w errorbars
    #                     # print figure_key
    #                     # print group_key
    #                     if mean:
    #                         if 'shade_error' in kwargs and kwargs['shade_error']:
    #                             plt.fill_between(t, data_baseline_mean-data_baseline_sem, data_baseline_mean+data_baseline_sem, color=color, alpha=1)
                                
    #                             # trace[group_key] = ax[figure_key].plot(t, data_baseline_mean, color=color, linewidth=4)
    #                             plt.fill_between(t, data_final_mean-data_final_sem, data_final_mean+data_final_sem, color=color_post, alpha=1)
    #                             # trace_final[group_key] = ax[figure_key].plot(t, data_final_mean, color=color_final, linewidth=4)
    #                             _plot_scalebar(axes=ax[figure_key], xscale=2, yscale=.001, origin=[9,-.0035], width=4)
    #                             plt.axis('off')
    #                         else:
    #                             trace[group_key] = ax[figure_key].errorbar(t, data_baseline_mean, yerr=data_baseline_sem, color=color, linewidth=1)
    #                             trace_final[group_key] = ax[figure_key].errorbar(t, data_final_mean, yerr=data_final_sem, color=color, linewidth=2)




    #         ylim = [-.004, .0005]
    #         # plt.title(titles[figure_key]+' '+variable)
    #         # plt.legend(trace.values(), label.values())
    #         if ylim:
    #             plt.ylim(ylim)
    #         plt.xlim([0,12])
    #         plt.show(block=False)

    #     return fig, ax

    # def _plot_timeseries(df_sorted, conditions, figures, variable, colors, markers, titles, mean=True, all_paths=True, **kwargs):
    #     '''FIXME add docs
    #     '''
    #     # FIXME add kwargs to alter figure details
    #     # create figure groupings (all conditions that will go on the same figure)
    #     fig = {}
    #     ax={}
    #     ylim={}
    #     xlim={}
    #     # iterate over figures
    #     peak_value_all=0
    #     for figure_key in figures:
    #         if type(figure_key)==str:
    #             figure_key=(figure_key,)
    #         # create figure objects
    #         fig[figure_key], ax[figure_key] = plt.subplots(figsize=(8,6), dpi=100 )
    #         # plot object for each trace
    #         trace = OrderedDict()
    #         # n for each condition
    #         n_trial=OrderedDict()
    #         # pvalue for ttest comparing to corresponding control condition
    #         pval=OrderedDict()
    #         # label for plot legend
    #         label=OrderedDict()
    #         # iterate over each condition
    #         peak_value_current = 0
    #         for group_key, df in df_sorted.iteritems():
    #             # print figure_key
    #             # print group_key
    #             # print len(set(group_key) - set(figure_key))
    #             # if type(group_key)==tuple:
    #             #     group_key = list(group_key)
    #             # elif type(group_key)==str:
    #             #     group_key = [group_key]
    #             # print len(set(group_key) - set(figure_key))
    #             # if plotting all pathways on the same figure
    #             if all_paths and not all(temp in group_key for temp in figure_key) and 'figures_any' not in kwargs:
    #                 # print 'huuuuuh'
    #                 # if condition not met continue
    #                 continue
    #             # if only plotting specific pathways
    #             elif not all_paths and not group_key[1]==figure_key[0] and group_key[2]==figure_key[1] and 'figures_any' not in kwargs:
    #                 # if condition not met, continue
    #                 # print 'huhhhhhhh'
    #                 continue

    #             if 'figures_any' in kwargs and len(set(group_key) - set(figure_key)) != 0:

    #                 # print figure_key, group_key,'huh'


    #             # any(set(group_key)&set(figure_key)):
    #                 continue
    #             # if conditions are satisfied, plot figure
    #             else:
    #                 # get color
    #                 color = colors[group_key[0]]

    #                 marker_match = [key for key in group_key if key in markers.keys()]
    #                 # list of group key values and list of marker keys,use the first match that you find
    #                 if 'all' in markers.keys():
    #                     marker=markers['all']

    #                 elif len(marker_match)>0:
    #                     marker = markers[marker_match[0]]


    #                 if 'line_markers' in kwargs:
    #                     line_markers=kwargs['line_markers']
    #                 else:
    #                     line_markers={}
    #                 line_marker_match = [key for key in group_key if key in line_markers.keys()]
    #                 # list of group key values and list of marker keys,use the first match that you find
    #                 if 'all' in line_markers.keys():
    #                     line_marker=line_markers['all']

    #                 elif len(line_marker_match)>0:
    #                     line_marker = line_markers[line_marker_match[0]]
    #                 # elif len(list(set(markers.keys()) & set(group_key)))==1:
    #                 #     marker_key = list(set(group_key) & set(markers.keys()))[0]
    #                 #     marker = markers[marker_key]

    #                 # elif len(list(set(markers.keys()) & set(group_key)))>1:

    #                 else:
    #                     line_marker='solid'


    #                 # print set(group_key)
    #                 # print set(markers.keys())
    #                 # print list(set(markers.keys()) & set(group_key))
    #                 # marker_key = list(set(markers.keys()) & set(group_key))[0]
    #                 # # get marker
    #                 # marker = markers[marker_key]
    #                 # convert data to 2d array (slices x time)
    #                 data_array = _2array(df[variable], remove_nans=True, remove_nans_axis=1)
    #                 print data_array.shape
    #                 # check that data is 2d (if not, it is due to a lack of data for this condition)
    #                 if len(data_array.shape)>1:
    #                     # mean across slices
    #                     data_mean = np.mean(data_array, axis=0)
    #                     #std across slices
    #                     data_std = np.std(data_array, axis=0)
    #                     # sem across slices
    #                     data_sem = stats.sem(data_array, axis=0)
    #                     # time vector
    #                     t = np.arange(len(data_mean))
    #                     print t.shape, data_mean.shape
    #                     # if plotting mean w errorbars
    #                     # print figure_key
    #                     # print group_key
    #                     if mean:
    #                         print np.max(data_mean+data_sem)
    #                         if np.max(data_mean+data_sem)>peak_value_current:
    #                             peak_value_current= np.max(data_mean+data_sem)
    #                         if peak_value_current>peak_value_all:
    #                             peak_value_all =  peak_value_current

    #                         print peak_value_all
    #                         if 'shade_error' in kwargs and kwargs['shade_error']:
    #                             plt.fill_between(t, data_mean-data_sem, data_mean+data_sem, color=color, alpha=1)

    #                             # trace[group_key] = ax[figure_key].plot(t, data_mean, color=color, marker=marker,  markersize=10, linewidth=0, markerfacecolor=color)

    #                             trace[group_key] = ax[figure_key].plot(t, data_mean, color=color, marker='None',  markersize=10, linewidth=6, markerfacecolor=color, linestyle=line_marker, )
    #                             # plt.fill_between(t, data_final_mean-data_final_sem, data_final_mean+data_final_sem, color=color_post, alpha=1)
    #                         else:
    #                             trace[group_key] = ax[figure_key].errorbar(t, data_mean, yerr=data_sem, color=color, marker=marker,  markersize=15, elinewidth=2, linewidth=3, markerfacecolor=color)
    #                     # otherwise plot all traces
    #                     else:
    #                         trace[group_key] = ax[figure_key].plot(t, data_array.T, color=color, marker=marker,  markersize=10)

    #                     # get number of slices
    #                     n_trial[group_key] = data_array.shape[0]

    #                     # build control group
    #                     control_group = []
    #                     for cond_i, cond in enumerate(conditions):
    #                         if cond == 'field_polarity_0':
    #                             control_group.append('control')
    #                         elif cond=='field_mag_0':
    #                             control_group.append('0')
    #                         else:
    #                             control_group.append(group_key[cond_i])
    #                     control_group= tuple(control_group)

    #                     if 'slopes_norm' in variable:
    #                         # final ltp (length n), taken as the average of the last 10 normalized slopes for each slice
    #                         ltp_final = df_sorted[group_key].ltp_final
    #                         # corresponding control values
    #                         ltp_final_control = df_sorted[control_group].ltp_final
    #                         # ttest compared to control
    #                         pval[group_key] = stats.ttest_ind(ltp_final, ltp_final_control)[1]
    #                         ltp_final_mean=np.mean(ltp_final)
    #                         ltp_final_std=np.std(ltp_final)
    #                         ltp_final_sem = stats.sem(ltp_final)
    #                         # label for figure legend including stimulation conditions, n slices, p value compared to the corresponding control
    #                         label_keys=[]
    #                         for key in group_key:
    #                             # print key
    #                             label_keys.append(key)
    #                         label_keys = label_keys + ['n={}'.format(n_trial[group_key]), 'p = {0:.3f}'.format(pval[group_key]), 'mean={0:.3f}'.format(ltp_final_mean),'sem={0:.3f}'.format(ltp_final_sem)]
    #                         # print label_keys
    #                         label[group_key] = ', '.join(label_keys)
    #                         # label[group_key] = ', '.join([group_key[0], group_key[1], 'n={}'.format(n_trial[group_key]), 'p = {0:.3f}'.format(pval[group_key])])
    #             # if not 'common_ylim' in kwargs or not kwargs['common_ylim']:
    #             # get x and y limits based data
    #             xlim[figure_key] = ax[figure_key].get_xlim()
    #             ylim[figure_key] = ax[figure_key].get_ylim()
    #             # ylim[figure_key] = [0.9, peak_value_current]
    #             # plt.title(titles[figure_key]+' '+variable)
    #             if 'show_stats' in kwargs and kwargs['show_stats']:
    #                 plt.legend(trace.values(), label.values())

                

    #     # if 'common_ylim' in kwargs and kwargs['common_ylim']:
    #     #     ylim = [0.9, peak_value_all]
    #     # else:
    #     #     ylim = [0.9, peak_value_current]
    #     # print ylim
    #     # xlim = [0,80]
    #     for figure_key in ax:
    #         if 'common_ylim' in kwargs and kwargs['common_ylim'] and 'ylim' not in kwargs:
    #             ylim[figure_key] = [0.9, peak_value_all]
    #         elif 'ylim' in kwargs:
    #             ylim[figure_key] = kwargs['ylim']

    #         if 'xlim' in kwargs:
    #             xlim[figure_key] = kwargs['xlim']
            
    #         # format figure
    #         ax[figure_key].spines['right'].set_visible(False)
    #         ax[figure_key].spines['top'].set_visible(False)
    #         ax[figure_key].spines['left'].set_linewidth(5)
    #         ax[figure_key].spines['bottom'].set_linewidth(5)
    #         ax[figure_key].xaxis.set_ticks_position('bottom')
    #         ax[figure_key].yaxis.set_ticks_position('left')
    #         ax[figure_key].set_xlabel('Time (min)', fontsize=25, fontweight='heavy')
    #         ax[figure_key].set_ylabel('Normalized fEPSP slope', fontsize=25, fontweight='heavy')
    #         xticks = np.arange(0,81, 20)
    #         yticks = np.round(np.arange(1.,peak_value_all, 0.2), decimals=1)
    #         ax[figure_key].set_xticks(xticks)
    #         ax[figure_key].set_xticklabels(xticks, fontsize=20, fontweight='heavy')
    #         ax[figure_key].set_yticks(yticks)
    #         ax[figure_key].set_yticklabels(yticks, fontsize=20, fontweight='heavy')
            
    #         ax[figure_key].set_ylim(ylim[figure_key])
    #         ax[figure_key].set_xlim(xlim[figure_key])
    #         plt.figure(fig[figure_key].number)
    #         plt.tight_layout()
            
    #         # if ylim:
    #         #     plt.ylim(ylim)
    #     plt.show(block=False)
        

    #     return fig, ax
        
    # def _plot_bar(df_sorted, figure_params, variable, group_space=1, bar_width=1, bar_spacing=1):
    #     '''
    #     '''
    #     print 'look here:',figure_params.keys()
    #     fig={}
    #     ax={}
    #     n_subgroups={}
    #     n_traces={}
    #     xlim={}
    #     ylim={}
    #     for figure_key, figure_subgroups in figure_params.iteritems():
    #         if figure_key!='params':
    #             fig[figure_key], ax[figure_key] = plt.subplots()
    #             n_subgroups[figure_key] = len(figure_subgroups.keys()) 
    #             n_traces[figure_key]={}
    #             locations = []
    #             heights=[]
    #             colors=[]
    #             fig_args=[]
    #             xticks=[]
    #             xticklabels=[]
    #             cnt=bar_spacing
    #             for subgroup_key, traces in figure_subgroups.iteritems():
    #                 if subgroup_key!='params':
    #                     n_traces[figure_key][subgroup_key]=len(traces.keys())
    #                     cnt+=group_space
    #                     for trace_key, params in traces.iteritems():
    #                         if trace_key!='params':
    #                             trace_args={}
    #                             # cnt+=bar_spacing
    #                             if 'location' in params:
    #                                 locations.append(cnt+params['location'])
    #                             else:
    #                                 locations.append(locations[-1]+1)
    #                             xticks.append(cnt)
    #                             xticklabels.append(params['label'])

    #                             # get data and stats
    #                             if trace_key in df_sorted:
    #                                 trace_series = df_sorted[trace_key][variable]
    #                                 data_array = (_2array(trace_series, remove_nans=True, remove_nans_axis=1)-1)*100.
    #                                 print type(data_array)
    #                                 print type(data_array)==np.float64
    #                                 if type(data_array)==np.ndarray:
    #                                     # get stats
    #                                     # mean across slices
    #                                     data_mean = np.mean(data_array, axis=0)
    #                                     #std across slices
    #                                     data_std = np.std(data_array, axis=0)
    #                                     # sem across slices
    #                                     data_sem = stats.sem(data_array, axis=0)

    #                                     heights.append(data_mean)

    #                                     trace_args['color'] = params['color']

    #                                     fig_args.append(trace_args)

    #                                     colors.append(params['color'])

    #                                     plt.errorbar(locations[-1], data_mean, data_sem, color=(.5,.5,.5))

    #             barcontainer = ax[figure_key].bar(locations, heights, width=bar_width, tick_label=xticklabels)
    #             xlim[figure_key] = ax[figure_key].get_xlim()
    #             ylim[figure_key] = ax[figure_key].get_ylim()
    #             # ax[figure_key].set_xticks(xticks, xticklabels,)
    #             # barcontainer = ax[figure_key].violinplot(locations, heights, width=bar_width, tick_label=xticklabels)
    #             print 'rotations:', figure_key, figure_params[figure_key]['params']
    #             ax[figure_key].set_xticklabels(xticklabels, fontsize
    #                 =20, fontweight='heavy', rotation=figure_params[figure_key]['params']['rotation'])
    #             for bar_i, bar in enumerate(barcontainer):
    #                 bar.set_color(colors[bar_i])

    #     # get ylim and xlim across all figures
    #     xlims=[]
    #     ylims=[]
    #     for figure_key in ylim:
    #         xlims.append(xlim[figure_key])
    #         ylims.append(ylim[figure_key])
    #     xlim_all = [min([temp[0] for temp in xlims]), max([temp[1] for temp in xlims])]
    #     ylim_all = [min([temp[0] for temp in ylims]), max([temp[1] for temp in ylims])]

    #     xlim={}
    #     ylim={}
    #     print figure_params.keys()
    #     for figure_key, axes in ax.iteritems():
    #         if 'ylim_all' in figure_params['params'] and figure_params['params']['ylim_all']:
    #             print 'setting ylim'
    #             ax[figure_key].set_ylim(ylim_all)
    #             ax[figure_key].set_xlim(xlim_all)

            
    #         # format figure
    #         ax[figure_key].spines['right'].set_visible(False)
    #         ax[figure_key].spines['top'].set_visible(False)
    #         ax[figure_key].spines['left'].set_linewidth(5)
    #         ax[figure_key].spines['bottom'].set_linewidth(5)
    #         ax[figure_key].xaxis.set_ticks_position('bottom')
    #         ax[figure_key].yaxis.set_ticks_position('left')
    #         # ax[figure_key].set_xlabel('Time (min)', fontsize=25, fontweight='heavy')
    #         ax[figure_key].set_ylabel('% LTP', fontsize=25, fontweight='heavy')
    #         xticks = np.arange(0,81, 20)
    #         ytickmax = ax[figure_key].get_ylim()[1]
    #         ytickmin = ax[figure_key].get_ylim()[0]
    #         yticks = np.round(np.arange(0,ytickmax, 10), decimals=0).astype(int)
    #         # ax[figure_key].set_xticks(xticks)
    #         # ax[figure_key].set_xticklabels(xticks, fontsize=20, fontweight='heavy')
    #         ax[figure_key].set_yticks(yticks)
    #         ax[figure_key].set_yticklabels(yticks, fontsize=20, fontweight='heavy')
            
    #         # ax[figure_key].set_ylim(ylim[figure_key])
    #         # ax[figure_key].set_xlim(xlim)
    #         plt.figure(fig[figure_key].number)
    #         plt.tight_layout()

    #     plt.show(block=False)

    #     return fig, ax

    # def _plot_var2var_correlation(df_sorted, figures, variables, colors, markers, titles):
    #     ''' FIXME adds docs
    #     '''
    #     # create figure groupings (all conditions that will go on the same figure)
    #     fig = {}
    #     ax={}
    #     # iterate over figures
    #     for figure_key in figures:
    #         # create figure objects
    #         fig[figure_key], ax[figure_key] = plt.subplots()
    #         # plot object for each trace
    #         trace = OrderedDict()
    #         # n for each condition
    #         n_trial=OrderedDict()
    #         # pvalue for ttest comparing to corresponding control condition
    #         pval=OrderedDict()
    #         # label for plot legend
    #         label=OrderedDict()
    #         # iterate over each condition
    #         var0_all = []
    #         var1_all = []
    #         for group_key, df in df_sorted.iteritems():
    #             if group_key[1]==figure_key[0] and group_key[2]==figure_key[1]:
    #             # check that the current conditions belong on the current figure
    #             # if all(temp in group_key for temp in figure_key):
    #                 color = colors[group_key[0]]
    #                 # print group_key[1]
    #                 marker = markers[group_key[1]]
    #                 var0 = np.array(df[variables[0]])
    #                 var1 = np.array(df[variables[1]])
    #                 var0_all.append(var0)
    #                 var1_all.append(var1)
    #                 trace[group_key] = ax[figure_key].plot(var0, var1, color=color, marker=marker,  markersize=10, linestyle='None')
    #                 # plt.plot(var, ltp, color=color, marker=marker)

    #         var0_all = np.append(var0_all[0], var0_all[1:])
    #         var1_all = np.append(var1_all[0], var1_all[1:])
    #         varmat = np.vstack([var0_all, var1_all])
    #         # print variables, varmat.shape
    #         # covariance = np.cov(varmat)
    #         pearsonsr, pearsonsp = stats.pearsonr(var0_all, var1_all)
    #         print variables[0], pearsonsr, pearsonsp

    #         plt.title(titles[figure_key])
    #         plt.xlabel(variables[0])
    #         plt.ylabel(variables[1])

    #         # plot correlation coefficient with label

    #     plt.show(block=False)

    # def _plot_trace_brian(df_sorted, figure_params, variable, **kwargs):
    #     '''FIXME add docs
    #     '''
    #     # FIXME add kwargs to alter figure details
    #     # create figure groupings (all conditions that will go on the same figure)
    #     fig={}
    #     ax={}
    #     n_subgroups={}
    #     n_traces={}
    #     xlim={}
    #     ylim={}
    #     # iterate over figures
    #     for figure_key, figure_subgroups in figure_params.iteritems():
    #         # params kw specifies figure level kwargs
    #         if figure_key!='params':
    #             # create figure, passing params as **kwargs
    #             fig[figure_key], ax[figure_key] = plt.subplots()
    #             # number of subgroup in figure
    #             n_subgroups[figure_key] = len(figure_subgroups.keys()) 
    #             n_traces[figure_key]={}
    #             # iterate over subgroups of traces
    #             for subgroup_key, traces in figure_subgroups.iteritems():

    #                 if subgroup_key!='params':
    #                     # FIXME distribute subgroup padrameters to each trace in the subgroup, with priority to trace parameters
    #                     n_traces[figure_key][subgroup_key]=len(traces.keys())
    #                     # get trial id's
    #                     trace_key_temp = sorted(traces.keys())
    #                     trace_key_temp = [temp for temp in trace_key_temp if temp!='params'][0]
    #                     print traces.keys()
    #                     if'plot_individual_trace' in traces['params'] and len(traces['params']['plot_individual_trace'])>0:
    #                         trial_ids=[]
    #                         for iloc in traces['params']['plot_individual_trace']:
    #                             trial_ids.append(df_sorted[trace_key_temp].trial_id.iloc[iloc])
                        

    #                     # iterate over individual traces
    #                     for trace_key, params in traces.iteritems():
    #                         if trace_key!='params':
    #                             # get data and stats
    #                             print df_sorted.keys()
    #                             trace_series = df_sorted[trace_key][variable]
    #                             print subgroup_key, trace_key, trace_series
    #                             data_array = _2array(trace_series, remove_nans=True, remove_nans_axis=1)#*1000
    #                             # get stats
    #                             # mean across slices
    #                             data_mean = np.mean(data_array, axis=0)
    #                             #std across slices
    #                             data_std = np.std(data_array, axis=0)
    #                             # sem across slices
    #                             data_sem = stats.sem(data_array, axis=0)
    #                             # time vector
    #                             t = np.arange(len(data_mean))#/10.
    #                             if 'plot_mean' in traces['params'] and traces['params']['plot_mean']:
    #                                 # line trace with shaded error
    #                                 #______________________________
    #                                 if 'shade_error' in params and params['shade_error']:

    #                                     ax[figure_key].plot(t, data_mean, color=params['color'], linewidth=params['linewidth'])
    #                                     plt.fill_between(t, data_mean-data_sem, data_mean+data_sem, **params['e_params'])

    #                             elif'plot_individual_trace' in traces['params'] and len(traces['params']['plot_individual_trace'])>0:
    #                                 for trial_id in trial_ids:
    #                                     trace_data = df_sorted[trace_key][df_sorted[trace_key].trial_id==trial_id][variable]
    #                                     trace_data = _2array(trace_data)#*1000
    #                                     t = np.arange(len(trace_data))#/10.
    #                                     ax[figure_key].plot(t, trace_data, color=params['color'], linewidth=params['linewidth'])

    #             # get x and y limits based data
    #             xlim[figure_key] = ax[figure_key].get_xlim()
    #             ylim[figure_key] = ax[figure_key].get_ylim()
        
    #     # get ylim and xlim across all figures
    #     xlims=[]
    #     ylims=[]
    #     for figure_key in ylim:
    #         xlims.append(xlim[figure_key])
    #         ylims.append(ylim[figure_key])
    #     xlim_all = [min([temp[0] for temp in xlims]), max([temp[1] for temp in xlims])]
    #     ylim_all = [min([temp[0] for temp in ylims]), max([temp[1] for temp in ylims])]

    #     # set common ylim across all figures
    #     for figure_key, axes in ax.iteritems():
    #         if 'ylim_all' in figure_params['params'] and figure_params['params']['ylim_all']:
    #             print 'setting ylim to be the same across all figures'
    #             if 'ylim' in figure_params['params']:
    #                 ax[figure_key].set_ylim(figure_params['params']['ylim'])
    #             else:
    #                 ax[figure_key].set_ylim(ylim_all)
    #         if 'xlim_all' in figure_params['params'] and figure_params['params']['xlim_all']:
    #             if 'xlim' in figure_params['params']:
    #                 ax[figure_key].set_xlim(figure_params['params']['xlim'])
    #             else:
    #                 ax[figure_key].set_xlim(xlim_all)

    #         # format figure
    #         ax[figure_key].spines['right'].set_visible(False)
    #         ax[figure_key].spines['top'].set_visible(False)
    #         ax[figure_key].spines['left'].set_linewidth(5)
    #         ax[figure_key].spines['bottom'].set_linewidth(5)
    #         ax[figure_key].xaxis.set_ticks_position('bottom')
    #         ax[figure_key].yaxis.set_ticks_position('left')
    #         ax[figure_key].set_xlabel(figure_params['params']['xlabel'], fontsize=25, fontweight='heavy')
    #         ax[figure_key].set_ylabel(figure_params['params']['ylabel'], fontsize=25, fontweight='heavy')
    #         for temp in ax[figure_key].get_xticklabels():
    #             temp.set_fontweight('heavy')
    #             temp.set_fontsize(15)
    #         for temp in ax[figure_key].get_yticklabels():
    #             temp.set_fontweight('heavy')
    #             temp.set_fontsize(15)
    #         # print xticks, xticklabels
    #         # ax[figure_key].set_xticks(xticks)
    #         # ax[figure_key].set_xticklabels(xticklabels, fontsize=20, fontweight='heavy')
    #         # ax[figure_key].tick_params(axis='both', labelsize=20, )
    #         # xticks = figure_params['xticks'] #np.arange(0,81, 20)
    #         # ytickmax = ax[figure_key].get_ylim()[1]
    #         # ytickmin = ax[figure_key].get_ylim()[0]
    #         # yticks = figure_params['yticks']#np.round(np.arange(0,ytickmax, 10), decimals=0).astype(int)
    #         # ax[figure_key].set_xticks(xticks)
    #         # ax[figure_key].set_xticklabels(xticks, fontsize=20, fontweight='heavy')
    #         # ax[figure_key].set_yticks(yticks)
    #         # ax[figure_key].set_yticklabels(yticks, fontsize=20, fontweight='heavy')
            
    #         # ax[figure_key].set_ylim(ylim[figure_key])
    #         # ax[figure_key].set_xlim(xlim)
    #         plt.figure(fig[figure_key].number)
    #         plt.tight_layout()
    #         plt.show(block=False)

    #     return fig, axes

    # def _plot_induction_bursts(df_sorted, figure_params, variable, **kwargs):
    #     '''FIXME add docs
    #     '''
    #     # FIXME add kwargs to alter figure details
    #     # create figure groupings (all conditions that will go on the same figure)
    #     fig={}
    #     ax={}
    #     n_subgroups={}
    #     n_traces={}
    #     xlim={}
    #     ylim={}
    #     average_bursts = figure_params['params']['average_across_bursts']
    #     burst_number = figure_params['params']['burst_number']
    #     induction_number = figure_params['params']['induction_number']
    #     yscale=figure_params['params']['yscale']
    #     xscale=figure_params['params']['xscale']
    #     # iterate over figures
    #     for figure_key, figure_subgroups in figure_params.iteritems():
    #         # params kw specifies figure level kwargs
    #         if figure_key!='params':
    #             # create figure, passing params as **kwargs
    #             fig[figure_key], ax[figure_key] = plt.subplots()
    #             # number of subgroup in figure
    #             n_subgroups[figure_key] = len(figure_subgroups.keys()) 
    #             n_traces[figure_key]={}
    #             # iterate over subgroups of traces
    #             for subgroup_key, traces in figure_subgroups.iteritems():

    #                 if subgroup_key!='params':
    #                     # FIXME distribute subgroup padrameters to each trace in the subgroup, with priority to trace parameters
    #                     n_traces[figure_key][subgroup_key]=len(traces.keys())
    #                     # get trial id's
    #                     trace_key_temp = sorted(traces.keys())
    #                     trace_key_temp = [temp for temp in trace_key_temp if temp!='params'][0]
                        
    #                     # iterate over individual traces
    #                     for trace_key, params in traces.iteritems():
    #                         if trace_key!='params':
    #                             # get data and stats
    #                             print df_sorted.keys()
    #                             trace_series = df_sorted[trace_key][variable]
    #                             print subgroup_key, trace_key, trace_series
    #                             # axis0:slices, axis1:time, axis2:bursts
    #                             data_array = _2array(trace_series, remove_nans=True, remove_nans_axis=1, list_index=induction_number)#*1000
    #                             if data_array.shape[0]!=0:
    #                                 print data_array.shape
    #                                 # if there is more than 1 slice for the given conditions, 
    #                                 # data array is 3d, slices x samples x bursts, and must be reduced to 2d
    #                                 if len(data_array.shape)>2:
    #                                     # average over bursts
    #                                     if average_bursts:
    #                                         data_array = np.mean(data_array, axis=2).squeeze()
    #                                     # show a specific burst
    #                                     else:
    #                                         data_array= data_array[:,:,burst_number].squeeze()

    #                                 data_array=data_array*yscale

    #                                 # get stats
    #                                 # mean across slices
    #                                 data_mean = np.mean(data_array, axis=0)
    #                                 #std across slices
    #                                 data_std = np.std(data_array, axis=0)
    #                                 # sem across slices
    #                                 data_sem = stats.sem(data_array, axis=0)
    #                                 # time vector
    #                                 t = np.arange(len(data_mean))*xscale
    #                                 if 'plot_mean' in traces['params'] and traces['params']['plot_mean']:
    #                                     # line trace with shaded error
    #                                     #______________________________
    #                                     if 'shade_error' in params and params['shade_error']:

    #                                         ax[figure_key].plot(t, data_mean, color=params['color'], linewidth=params['linewidth'])
    #                                         plt.fill_between(t, data_mean-data_sem, data_mean+data_sem, **params['e_params'])

    #             # get x and y limits based data
    #             xlim[figure_key] = ax[figure_key].get_xlim()
    #             ylim[figure_key] = ax[figure_key].get_ylim()
    #             print 'ylim:', ylim
        
    #     # get ylim and xlim across all figures
    #     xlims=[]
    #     ylims=[]
    #     for figure_key in ylim:
    #         xlims.append(xlim[figure_key])
    #         ylims.append(ylim[figure_key])
    #     xlim_all = [min([temp[0] for temp in xlims]), max([temp[1] for temp in xlims])]
    #     ylim_all = [min([temp[0] for temp in ylims]), max([temp[1] for temp in ylims])]

    #     # set common ylim across all figures
    #     for figure_key, axes in ax.iteritems():
    #         if 'ylim_all' in figure_params['params'] and figure_params['params']['ylim_all']:
    #             print 'setting ylim to be the same across all figures'
    #             if 'ylim' in figure_params['params']:
    #                 ax[figure_key].set_ylim(figure_params['params']['ylim'])
    #             else:
    #                 ax[figure_key].set_ylim(ylim_all)
    #         if 'xlim_all' in figure_params['params'] and figure_params['params']['xlim_all']:
    #             if 'xlim' in figure_params['params']:
    #                 ax[figure_key].set_xlim(figure_params['params']['xlim'])
    #             else:
    #                 ax[figure_key].set_xlim(xlim_all)

    #         # format figure
    #         ax[figure_key].spines['right'].set_visible(False)
    #         ax[figure_key].spines['top'].set_visible(False)
    #         ax[figure_key].spines['left'].set_linewidth(5)
    #         ax[figure_key].spines['bottom'].set_linewidth(5)
    #         ax[figure_key].xaxis.set_ticks_position('bottom')
    #         ax[figure_key].yaxis.set_ticks_position('left')
    #         ax[figure_key].set_xlabel(figure_params['params']['xlabel'], fontsize=25, fontweight='heavy')
    #         ax[figure_key].set_ylabel(figure_params['params']['ylabel'], fontsize=25, fontweight='heavy')
    #         for temp in ax[figure_key].get_xticklabels():
    #             temp.set_fontweight('heavy')
    #             temp.set_fontsize(15)
    #         for temp in ax[figure_key].get_yticklabels():
    #             temp.set_fontweight('heavy')
    #             temp.set_fontsize(15)
    #         # print xticks, xticklabels
    #         # ax[figure_key].set_xticks(xticks)
    #         # ax[figure_key].set_xticklabels(xticklabels, fontsize=20, fontweight='heavy')
    #         # ax[figure_key].tick_params(axis='both', labelsize=20, )
    #         # xticks = figure_params['xticks'] #np.arange(0,81, 20)
    #         # ytickmax = ax[figure_key].get_ylim()[1]
    #         # ytickmin = ax[figure_key].get_ylim()[0]
    #         # yticks = figure_params['yticks']#np.round(np.arange(0,ytickmax, 10), decimals=0).astype(int)
    #         # ax[figure_key].set_xticks(xticks)
    #         # ax[figure_key].set_xticklabels(xticks, fontsize=20, fontweight='heavy')
    #         # ax[figure_key].set_yticks(yticks)
    #         # ax[figure_key].set_yticklabels(yticks, fontsize=20, fontweight='heavy')
            
    #         # ax[figure_key].set_ylim(ylim[figure_key])
    #         # ax[figure_key].set_xlim(xlim)
    #         plt.figure(fig[figure_key].number)
    #         plt.tight_layout()
    #         plt.show(block=False)

    #     return fig, axes

    # def _plot_dose_response_all(df_sorted, figure_params, variable, group_space=1, bar_width=1, bar_spacing=1):
    #     '''
    #     '''
    #     print 'look here:',figure_params['params'].keys()
    #     if 'markersize' in figure_params['params']:
    #         markersize=figure_params['params']['markersize']
    #     fig={}
    #     ax={}
    #     n_subgroups={}
    #     n_traces={}
    #     xlim={}
    #     ylim={}
    #     for figure_key, figure_subgroups in figure_params.iteritems():
    #         if figure_key!='params':
    #             fig[figure_key], ax[figure_key] = plt.subplots()
    #             n_subgroups[figure_key] = len(figure_subgroups.keys()) 
    #             n_traces[figure_key]={}
    #             locations = []
    #             heights=[]
    #             colors=[]
    #             fig_args=[]
    #             xticks=[]
    #             xticklabels=[]
    #             cnt=bar_spacing
    #             for subgroup_key, traces in figure_subgroups.iteritems():
    #                 if subgroup_key!='params':
    #                     n_traces[figure_key][subgroup_key]=len(traces.keys())
    #                     cnt+=group_space
    #                     for trace_key, params in traces.iteritems():
    #                         if trace_key!='params':
    #                             trace_args={}
    #                             # cnt+=bar_spacing
    #                             if 'location' in params:
    #                                 locations.append(params['location'])
    #                             else:
    #                                 locations.append(locations[-1]+1)
    #                             xticks.append(cnt)
    #                             xticklabels.append(params['label'])

    #                             current_location = params['location']

    #                             # FIXME handle single data points
    #                             # get data and standardts
    #                             if trace_key in df_sorted:
    #                                 trace_series = df_sorted[trace_key][variable]
    #                                 data_array = (_2array(trace_series, remove_nans=True, remove_nans_axis=1)-1)*100.
    #                                 print data_array
    #                                 print current_location
    #                                 if type(data_array)==np.float64:
    #                                     current_locations = current_location
    #                                 else:
    #                                     current_locations = current_location*np.ones(len(data_array))

    #                                 plt.plot(current_locations, data_array, '.', color=params['color'], markersize=markersize)
    #                                 # print trace_key, data_array.shape, type(data_array)
    #                                 # if type(data_array)==np.ndarray:
    #                                 #     # get stats
    #                                 #     # mean across slices
    #                                 #     data_mean = np.mean(data_array, axis=0)
    #                                 #     #std across slices
    #                                 #     data_std = np.std(data_array, axis=0)
    #                                 #     # sem across slices
    #                                 #     data_sem = stats.sem(data_array, axis=0)

    #                                 #     heights.append(data_mean)

    #                                 #     trace_args['color'] = params['color']

    #                                 #     fig_args.append(trace_args)

    #                                 #     colors.append(params['color'])

    #                                 #     plt.errorbar(locations[-1], data_mean, data_sem, color=(.5,.5,.5))

    #             # barcontainer = ax[figure_key].bar(locations, heights, width=bar_width, tick_label=xticklabels)
    #             # xlim[figure_key] = ax[figure_key].get_xlim()
    #             # ylim[figure_key] = ax[figure_key].get_ylim()
    #             # # ax[figure_key].set_xticks(xticks, xticklabels,)
    #             # # barcontainer = ax[figure_key].violinplot(locations, heights, width=bar_width, tick_label=xticklabels)
    #             # print 'rotations:', figure_key, figure_params[figure_key]['params']
    #             # ax[figure_key].set_xticklabels(xticklabels, fontsize
    #             #     =20, fontweight='heavy', rotation=figure_params[figure_key]['params']['rotation'])
    #             # for bar_i, bar in enumerate(barcontainer):
    #             #     bar.set_color(colors[bar_i])
    #             plt.show(block=False)

    #     return fig, ax

    # def _plot_dose_response_mean(df_sorted, figure_params, variable, group_space=1, bar_width=1, bar_spacing=1):
    #     '''
    #     '''
    #     print 'look here:',figure_params['params'].keys()
    #     if 'markersize' in figure_params['params']:
    #         markersize=figure_params['params']['markersize']
    #     fig={}
    #     ax={}
    #     n_subgroups={}
    #     n_traces={}
    #     xlim={}
    #     ylim={}
    #     for figure_key, figure_subgroups in figure_params.iteritems():
    #         if figure_key!='params':
    #             fig[figure_key], ax[figure_key] = plt.subplots()
    #             n_subgroups[figure_key] = len(figure_subgroups.keys()) 
    #             n_traces[figure_key]={}
    #             locations = []
    #             heights=[]
    #             colors=[]
    #             fig_args=[]
    #             xticks=[]
    #             xticklabels=[]
    #             cnt=bar_spacing
    #             for subgroup_key, traces in figure_subgroups.iteritems():
    #                 if subgroup_key!='params':
    #                     n_traces[figure_key][subgroup_key]=len(traces.keys())
    #                     cnt+=group_space
    #                     for trace_key, params in traces.iteritems():
    #                         if trace_key!='params':
    #                             trace_args={}
    #                             # cnt+=bar_spacing
    #                             if 'location' in params:
    #                                 locations.append(params['location'])
    #                             else:
    #                                 locations.append(locations[-1]+1)
    #                             xticks.append(cnt)
    #                             xticklabels.append(params['label'])

    #                             current_location = params['location']

    #                             # FIXME handle single data points
    #                             # get data and standardts
    #                             if trace_key in df_sorted:
    #                                 trace_series = df_sorted[trace_key][variable]
    #                                 data_array = (_2array(trace_series, remove_nans=True, remove_nans_axis=1)-1)*100.
    #                                 # print data_array
    #                                 # print current_location
    #                                 # if type(data_array)==np.float64:
    #                                 #     current_locations = current_location
    #                                 # else:
    #                                 #     current_locations = current_location*np.ones(len(data_array))

    #                                 # plt.plot(current_locations, data_array, '.', color=params['color'], markersize=markersize)
    #                                 # print trace_key, data_array.shape, type(data_array)

    #                                 if type(data_array)==np.float64:
    #                                     # mean across slices
    #                                     data_mean = data_array
    #                                     #std across slices
    #                                     data_std = 0.
    #                                     # sem across slices
    #                                     data_sem = 0.

    #                                 elif type(data_array)==np.ndarray:
    #                                     # get stats
    #                                     # mean across slices
    #                                     data_mean = np.mean(data_array, axis=0)
    #                                     #std across slices
    #                                     data_std = np.std(data_array, axis=0)
    #                                     # sem across slices
    #                                     data_sem = stats.sem(data_array, axis=0)

    #                                 heights.append(data_mean)

    #                                 trace_args['color'] = params['color']

    #                                 fig_args.append(trace_args)

    #                                 colors.append(params['color'])

    #                                 plt.errorbar(locations[-1], data_mean, data_sem, color=(.5,.5,.5))

    #             plt.plot(locations, heights, '.', color=params['color'], markersize=markersize)
    #             # barcontainer = ax[figure_key].bar(locations, heights, width=bar_width, tick_label=xticklabels)
    #             # xlim[figure_key] = ax[figure_key].get_xlim()
    #             # ylim[figure_key] = ax[figure_key].get_ylim()
    #             # # ax[figure_key].set_xticks(xticks, xticklabels,)
    #             # # barcontainer = ax[figure_key].violinplot(locations, heights, width=bar_width, tick_label=xticklabels)
    #             # print 'rotations:', figure_key, figure_params[figure_key]['params']
    #             # ax[figure_key].set_xticklabels(xticklabels, fontsize
    #             #     =20, fontweight='heavy', rotation=figure_params[figure_key]['params']['rotation'])
    #             # for bar_i, bar in enumerate(barcontainer):
    #             #     bar.set_color(colors[bar_i])
    #             plt.show(block=False)

    #     return fig, ax

    # def _build_figure_params(figures_dict):
    #     '''
    #     '''
    #     # build multiindex for trace parameters
    #     multi_list = []
    #     level_names = ['figure','subgroup', 'trace']
    #     for level_1_key, level_1 in figures_dict.iteritems():
    #         for level_2_key, level_2 in level_1.iteritems():
    #             for level_3_i, level_3_key in enumerate(level_2):
    #                 multi_list.append((level_1_key, level_2_key, level_3_key))
    #     multiindex = pd.MultiIndex.from_tuples(multi_list, names=level_names)
    #     # build dataframe
    #     figure_df = pd.DataFrame(index=multiindex, dtype='object')

    #     return figure_df

    # def _plot_trace_mean(df_sorted, figdf, variable, **kwargs):
    #     '''FIXME add docs
    #     '''
    #     # FIXME add kwargs to alter figure details
    #     # create figure groupings (all conditions that will go on the same figure)
    #     fig={}
    #     ax={}
    #     n_subgroups={}
    #     n_traces={}
    #     xlim={}
    #     ylim={}
    #     # set figdf to hierarchical index (figure, subgroup, trace)
    #     figdf = figdf.reset_index().set_index(['figure','subgroup','trace'])
    #     # get list of all figure names
    #     figures = figdf.index.get_level_values('figure').unique().values
    #     # iterate over figures
    #     for figkey in figures:
    #         # create figure, passing params as **kwargs
    #         fig[figkey], ax[figkey] = plt.subplots()
    #         # get subgroups list
    #         subgroups = figdf.loc[figkey].index.get_level_values('subgroup').unique().values
    #         # subgroups = figdf[figdf.figure==figkey].subgroup.unique()
    #         # number of subgroup in figure
    #         n_subgroups[figkey] = len(subgroups) 
    #         n_traces[figkey]={}
    #         # iterate over subgroups of traces
    #         for subkey in subgroups:
    #             traces = figdf.loc[(figkey, subkey)].index.get_level_values('trace').unique().values
    #             # FIXME distribute subgroup padrameters to each trace in the subgroup, with priority to trace parameters
    #             n_traces[figkey][subkey]=len(traces)

    #             # iterate over traces
    #             for tracekey in traces:
    #                 params = figdf.loc[(figkey,subkey, tracekey)]
    #                 # get series from df sorted
    #                 trace_series = df_sorted[tracekey][variable]
    #                 # convert to array
    #                 data_array = _2array(trace_series, remove_nans=True, remove_nans_axis=1)#*1000
    #                 # get stats
    #                 # mean across slices
    #                 data_mean = np.mean(data_array, axis=0)
    #                 #std across slices
    #                 data_std = np.std(data_array, axis=0)
    #                 # sem across slices
    #                 data_sem = stats.sem(data_array, axis=0)
    #                 # time vector
    #                 t = np.arange(len(data_mean))#
    #                 # line plot with shaded error
    #                 if figdf.loc[(figkey,subkey,tracekey)].error_style=='shade':
    #                     ax[figkey].plot(t, data_mean, color=params.trace_color, linewidth=params.trace_linewidth)
    #                     plt.fill_between(t, data_mean-data_sem, data_mean+data_sem, color=params.error_color, alpha=params.error_alpha)
    #                 # error bar plot
    #                 elif figdf.loc[(figkey,subkey,tracekey)].error_style=='bar':
    #                     ax[figkey].errorbar(t, data_mean, yerr=data_sem, color=params.trace_color, 
    #                         marker=params.trace_marker,  
    #                         markersize=params.markersize, 
    #                         elinewidth=params.error_linewidth, 
    #                         linewidth=params.trace_linewidth, 
    #                         markerfacecolor=params.trace_color, 
    #                         ecolor=params.error_color)
    #         # get x and y limits based data
    #         xlim[figkey] = ax[figkey].get_xlim()
    #         ylim[figkey] = ax[figkey].get_ylim()

    #     fig, ax = _plot_format_figures(fig=fig, ax=ax, figdf=figdf, xlim=xlim, ylim=ylim)

    #     plt.show(block=False)

    #     return fig, ax

    # def _plot_bar2(df_sorted, figdf, variable, group_space=1, bar_width=1, bar_spacing=1):
    #     '''
    #     '''
    #     fig={}
    #     ax={}
    #     n_subgroups={}
    #     n_traces={}
    #     xlim={}
    #     ylim={}
    #     # set figdf to hierarchical index (figure, subgroup, trace)
    #     figdf = figdf.reset_index().set_index(['figure','subgroup','trace'])
    #     # get list of all figure names
    #     figures = figdf.index.get_level_values('figure').unique().values
    #     # iterate over figures
    #     for figkey in figures:
    #         fig[figkey], ax[figkey] = plt.subplots()
    #         # get subgroups list
    #         subgroups = figdf.loc[figkey].index.get_level_values('subgroup').unique().values
    #         n_subgroups[figkey] = len(subgroups) 
    #         n_traces[figkey]={}
    #         locations = []
    #         heights=[]
    #         colors=[]
    #         fig_args=[]
    #         xticks=[]
    #         xticklabels=[]
    #         cnt=bar_spacing
    #         # iterate over subgroups of traces
    #         for subkey in subgroups:
    #             traces = figdf.loc[(figkey, subkey)].index.get_level_values('trace').unique().values
    #             # FIXME distribute subgroup padrameters to each trace in the subgroup, with priority to trace parameters
    #             n_traces[figkey][subkey]=len(traces)
    #             cnt+=group_space
    #             # iterate over traces
    #             for tracekey in traces:
    #                 param = figdf.loc[(figkey,subkey, tracekey)]
    #                 # get series from df sorted
    #                 trace_series = df_sorted[tracekey][variable]
    #                 # convert to array
    #                 data_array = _2array(trace_series, remove_nans=True, remove_nans_axis=1)#*1000
    #                 # print figdf.keys()
    #                 if all(figdf.loc[figkey].fig_topercent):
    #                     data_array = 100.*(data_array-1)
    #                 # get stats
    #                 if type(data_array)==np.ndarray:
    #                     # mean across slices
    #                     data_mean = np.mean(data_array, axis=0)
    #                     #std across slices
    #                     data_std = np.std(data_array, axis=0)
    #                     # sem across slices
    #                     data_sem = stats.sem(data_array, axis=0)
    #                     # add plot location
    #                     print figkey, subkey, tracekey, param.sub_location, param.trace_location
    #                     plot_location = param.sub_location+param.trace_location
    #                     locations.append(plot_location)
    #                     xticks.append(plot_location)
    #                     xticklabels.append(param.trace_label)
    #                     colors.append(param.trace_color)
    #                     heights.append(data_mean)
    #                     plt.errorbar(locations[-1], data_mean, data_sem, color=param.error_color)

    #         barcontainer = ax[figkey].bar(locations, heights, width=param.fig_barwidth, tick_label=xticklabels)

    #         # get x and y lims
    #         xlim[figkey] = ax[figkey].get_xlim()
    #         ylim[figkey] = ax[figkey].get_ylim()
    #         # rotate ticklabels
    #         ax[figkey].set_xticklabels(xticklabels, fontsize
    #             =figdf.loc[figkey].fig_xtick_fontsize.unique()[0], fontweight=figdf.loc[figkey].fig_xtick_fontweight.unique()[0], rotation=figdf.loc[figkey].fig_xtick_rotate.unique()[0])
    #         # set bar color
    #         for bar_i, bar in enumerate(barcontainer):
    #             bar.set_color(colors[bar_i])

    #     fig, ax = _plot_format_figures(fig=fig, ax=ax, figdf=figdf, xlim=xlim, ylim=ylim)

    #     plt.show(block=False)

    #     return fig, ax

    # def _plot_format_figures(fig, ax, figdf, xlim, ylim):
    #     '''
    #     '''
    #     # get ylim and xlim across all figures
    #     #-------------------------------------
    #     xlims=[]
    #     ylims=[]
    #     for figkey in ylim:
    #         xlims.append(xlim[figkey])
    #         ylims.append(ylim[figkey])
    #     xlim_all = [min([temp[0] for temp in xlims]), max([temp[1] for temp in xlims])]
    #     ylim_all = [min([temp[0] for temp in ylims]), max([temp[1] for temp in ylims])]

    #     # iterate over figures and update attributes
    #     #---------------------------------------------
    #     for figkey, axes in ax.iteritems():
    #         ylim_current = ax[figkey].get_ylim()
    #         xlim_current = ax[figkey].get_xlim()
    #         # set ylim across all figures
    #         #----------------------------
    #         if all(figdf.fig_ylim_all):
    #             if 'fig_ylim' in figdf.keys():
    #                 ax[figkey].set_ylim(figdf.loc[figkey].fig_ylim.iloc[0])
    #             elif 'fig_ymin' in figdf.keys():
    #                 ax[figkey].set_ylim([figdf.loc[figkey].fig_ymin.iloc[0], ylim_all[1]])
    #             elif 'fig_ymax' in figdf.keys():
    #                 ax[figkey].set_ylim([ylim_all[0], figdf.loc[figkey].fig_ymax.iloc[0]])
    #             else:
    #                 ax[figkey].set_ylim(ylim_all)

    #         elif 'fig_ymin' in figdf.loc[figkey] and not any(figdf.loc[figkey].fig_ymin.isnull()):
    #             ax[figkey].set_ylim([figdf.loc[figkey].fig_ymin.unique()[0], ylim_current[1]])
    #         # set xlim across all figures
    #         #---------------------------
    #         if all(figdf.fig_xlim_all):
    #             if 'fig_xlim' in figdf.keys():
    #                 ax[figkey].set_xlim(figdf.loc[figkey].fig_xlim.iloc[0])
    #             elif 'fig_xmin' in figdf.keys():
    #                 ax[figkey].set_xlim([figdf.loc[figkey].fig_xmin.iloc[0], xlim_all[1]])
    #             elif 'fig_xmax' in figdf.keys():
    #                 ax[figkey].set_xlim([xlim_all[0], figdf.loc[figkey].fig_xmax.iloc[0]])
    #             else:
    #                 ax[figkey].set_xlim(xlim_all)

    #         elif 'fig_xmin' in figdf.loc[figkey] and not any(figdf.loc[figkey].fig_xmin.isnull()):
    #             ax[figkey].set_xlim([figdf.loc[figkey].fig_xmin.unique()[0], xlim_current[1]])

    #         # turn off axes box
    #         #------------------
    #         if all(figdf.loc[figkey].fig_boxoff):
    #             ax[figkey].spines['right'].set_visible(False)
    #             ax[figkey].spines['top'].set_visible(False)

    #         # set axes linewidth and tick position
    #         #----------------------
    #         ax[figkey].spines['left'].set_linewidth(figdf.loc[figkey].fig_axes_linewidth.unique()[0])
    #         ax[figkey].spines['bottom'].set_linewidth(figdf.loc[figkey].fig_axes_linewidth.unique()[0])
    #         ax[figkey].xaxis.set_ticks_position('bottom')
    #         ax[figkey].yaxis.set_ticks_position('left')

    #         # set axes labels
    #         #----------------
    #         ax[figkey].set_xlabel(figdf.loc[figkey].fig_xlabel.unique()[0], fontsize=figdf.loc[figkey].fig_xlabel_fontsize.unique()[0], fontweight=figdf.loc[figkey].fig_xlabel_fontweight.unique()[0])
    #         ax[figkey].set_ylabel(figdf.loc[figkey].fig_ylabel.unique()[0], fontsize=figdf.loc[figkey].fig_ylabel_fontsize.unique()[0], fontweight=figdf.loc[figkey].fig_ylabel_fontweight.unique()[0])

    #         # set x and y ticks
    #         #--------------------
    #         xlim[figkey] = ax[figkey].get_xlim()
    #         ylim[figkey] = ax[figkey].get_ylim()
    #         xticks = ax[figkey].get_xticks()
    #         yticks = ax[figkey].get_yticks()
    #         if 'fig_dyticks' in figdf:
    #             yticks = np.arange(ylim[figkey][0], ylim[figkey][1], figdf.loc[figkey].fig_dyticks.unique()[0])
    #         elif 'fig_nyticks' in figdf:
    #             dstep = float(abs(ylim[figkey][1]-ylim[figkey][0]))/figdf.loc[figkey].fig_nyticks.unique()[0]
    #             yticks = np.arange(ylim[figkey][0], ylim[figkey][1], dstep)

    #         if 'fig_dxticks' in figdf:
    #             xticks = np.arange(xlim[figkey][0], xlim[figkey][1], figdf.loc[figkey].fig_dxticks.unique()[0])
    #         elif 'fig_nxticks' in figdf:
    #             dstep = float(abs(xlim[figkey][1]-xlim[figkey][0]))/figdf.loc[figkey].fig_nxticks.unique()[0]
    #             xticks = np.arange(xlim[figkey][0], xlim[figkey][1], dstep)
            
    #         ax[figkey].set_yticks(yticks)
    #         ax[figkey].set_xticks(xticks)

    #         # set ticklabel attributes
    #         #--------------------------
    #         for temp in ax[figkey].get_xticklabels():
    #             temp.set_fontweight(figdf.loc[figkey].fig_xtick_fontweight.unique()[0])
    #             temp.set_fontsize(figdf.loc[figkey].fig_xtick_fontsize.unique()[0])
    #         for temp in ax[figkey].get_yticklabels():
    #             temp.set_fontweight(figdf.loc[figkey].fig_ytick_fontweight.unique()[0])
    #             temp.set_fontsize(figdf.loc[figkey].fig_ytick_fontsize.unique()[0])

    #         # set tight layout
    #         #-------------------
    #         if all(figdf.loc[figkey].fig_tight_layout):
    #             plt.figure(fig[figkey].number)
    #             plt.tight_layout()

    #     return fig, ax


            
