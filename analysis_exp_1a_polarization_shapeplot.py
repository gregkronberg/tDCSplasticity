"""
analysis

"""
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
from matplotlib import cm

'''
'''
kwargs = {'experiment':'exp_1a_polarization_shapeplot'}
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
# functions
funcs = analysis.VarFuncs()

#############################################################################
# load variables
#############################################################################
# load group data
#-------------------------------------------------------------
vtrace_df = analysis._load_group_data(directory=directory, filename=filename_vtrace)

#############################################################################
# shapeplot of membrane polarization
#############################################################################
figtype='polarization_shapeplot_'
figdf = analysis.BuildFigDF()._shapeplot()
figs, ax = analysis.PlotFuncs()._shapeplot(df=vtrace_df, figdf=figdf, variable='polarization')

# save figure
#------------
for fig_key, fig in figs.iteritems():
    fname = figure_directory+figtype+str(fig_key)+'.png'
    fig.savefig(fname, format='png', dpi=dpi)


#############################################################################
# shapeplot of membrane polarization
#############################################################################
figtype='polarization_shapeplot_gray_'
figdf = analysis.BuildFigDF()._shapeplot()
figs, ax = analysis.PlotFuncs()._shapeplot(df=vtrace_df, figdf=figdf, variable='polarization', cmap=cm.binary)

# save figure
#------------
for fig_key, fig in figs.iteritems():
    fname = figure_directory+figtype+str(fig_key)+'.png'
    fig.savefig(fname, format='png', dpi=dpi)



# #####################################################################
# # weights
# #####################################################################
# run_weights=True
# if run_weights:


#     # dose response figures
#     #-------------------------------------------------------------------
#     figtype='dose_response_'
#     figdf = analysis.BuildFigDF()._dose_response()
#     figs, ax = analysis.PlotFuncs()._dose_response(df=w_clopath_df, figdf=figdf, variable='dw_clopath')
#     # save figure
#     #------------
#     for fig_key, fig in figs.iteritems():
#         fname = figure_directory+figtype+str(fig_key)+'.png'
#         fig.savefig(fname, format='png', dpi=dpi)

#     # dose response figures
#     #-------------------------------------------------------------------
#     figtype='dose_response_stdp_param_'
#     figdf = analysis.BuildFigDF()._dose_response()
#     figs, ax = analysis.PlotFuncs()._dose_response(df=w_clopath_stdp_df, figdf=figdf, variable='dw_clopath')
#     # save figure
#     #------------
#     for fig_key, fig in figs.iteritems():
#         fname = figure_directory+figtype+str(fig_key)+'.png'
#         fig.savefig(fname, format='png', dpi=dpi)