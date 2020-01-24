"""
analysis

"""
import vargen_functions
import analysis
import time
import functions

kwargs = {
'experiment': '_'.join(__name__.split('_')[1:]).split('.')[0],#'exp_branch_sequence',
'clopath_param':vargen_functions.clopath_param_all_temp,
'threshold':-30,
}
directory = 'Data/'+kwargs['experiment']+'/'


varfuncs=functions.VarFuncs()
df_funcs=functions.DfFuncs()
vargen_funcs=functions.VarGen()

#############################################################################
# vtrace
#############################################################################
# check if dataframe exists in namespace
try:
    vtrace_df
except NameError:
    var_exists = False
else:
    var_exists = True
def _vtrace(vargen_funcs):
    variable='vtrace_df'
    filename=variable+'.pkl'
    functions=[varfuncs._pass]
    kwlist=[kwargs]
    rerun=[]#[varfuncs._get_vtrace]
    keep=[]
    file_limit=[]
    vtrace_df = vargen_funcs._vargen(variable=variable, directory=directory, functions=functions, kwlist=kwlist, rerun=rerun, keep=keep, file_limit=file_limit)
    return vtrace_df

# if not var_exists:
vtrace_df = _vtrace(vargen_funcs)

# #############################################################################
# # spikes
# #############################################################################
# def _spikes(vtrace_df, df_funcs):
#     '''
#     '''
#     variable='spikes_df'
#     start=time.time()
#     spikes_df = df_funcs._get_spikes(vtrace_df=vtrace_df, threshold=-40)
#     end=time.time()
#     print 'timer:', (end-start)
#     filename=variable+'.pkl'
#     print 'saving', filename
#     spikes_df.to_pickle(directory+filename)
#     return spikes_df

# spikes_df = _spikes(vtrace_df, df_funcs)

# def _spikes_xcorr(spikeS_df, df_funcs):
#     '''
#     '''
#     start=time.time()
#     spikes_df = df_funcs._get_xcorr_soma(spikes_df=spikes_df, threshold=-40, file_limit=[], )
#     end=time.time()
#     print 'timer:', (end-start)
#     filename=variable+'.pkl'
#     print 'saving', filename
#     spikes_df.to_pickle(directory+filename)
#     return spikes_df

# # spikes_df= _spikes_xcorr(spikes_df, df_funcs)

    
# #############################################################################
# # create weight df from vtrace df dcs parameters
# #############################################################################
# def _w_clopath(vtrace_df, vargen_functions):
#     variable='w_clopath_df'
#     start=time.time()
#     w_clopath_df = df_funcs._get_w_clopath(vtrace_df=vtrace_df, clopath_param=vargen_functions.clopath_param_all_temp)
#     end=time.time()
#     print 'timer:', (end-start)
#     filename=variable+'.pkl'
#     print 'saving', filename
#     w_clopath_df.to_pickle(directory+filename)
#     return w_clopath_df

# w_clopath_df=_w_clopath(vtrace_df, vargen_functions)

# # get dcs effect on final weights
# #---------------------------------
# def _get_dcs_effect_on_weights(w_clopath_df):
#     variable='w_clopath_df'
#     variables=['dw_clopath']
#     start=time.time()
#     w_clopath_df = df_funcs._get_dcs_effect(df=w_clopath_df, variables=variables,rerun=True )
#     end=time.time()
#     print 'timer:', (end-start)
#     filename=variable+'.pkl'
#     print 'saving', filename
#     w_clopath_df.to_pickle(directory+filename)

#     return w_clopath_df

# w_clopath_df=_get_dcs_effect_on_weights(w_clopath_df)