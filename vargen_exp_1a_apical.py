"""
vargen
"""
import vargen_functions
import analysis

kwargs = {
'experiment':'exp_1a_apical',
'clopath_param':vargen_functions.clopath_param_all,
'threshold':-30,
}
directory = 'Data/'+kwargs['experiment']+'/'
varfuncs=vargen_functions.VarFuncs()
df_funcs=vargen_functions.DfFuncs()
vargen_funcs=vargen_functions.VarGen()

#############################################################################
# vtrace
#############################################################################
variable='vtrace_df'
filename=variable+'.pkl'
functions=[varfuncs._get_vtrace]
kwlist=[kwargs]
rerun=[]
keep=[]
file_limit=[]
vtrace_df = vargen_funcs._vargen(variable=variable, directory=directory, functions=functions, kwlist=kwlist, rerun=rerun, keep=keep, file_limit=file_limit)

#############################################################################
# spikes
#############################################################################
variable='spikes_df'
functions=[varfuncs._get_spikes]
kwlist=[kwargs]
rerun=[]
keep=[]
file_limit=[]
spikes_df = vargen_funcs._vargen(variable=variable, directory=directory, functions=functions, kwlist=kwlist, rerun=rerun, keep=keep, file_limit=file_limit)

#############################################################################
# weights
#############################################################################
variable='w_clopath_df'
functions=[varfuncs._get_w_clopath]
kwlist=[kwargs]
rerun=[]
keep=[]
file_limit=[]
w_clopath_df = vargen_funcs._vargen(variable=variable, directory=directory, functions=functions, kwlist=kwlist, rerun=rerun, keep=keep, file_limit=file_limit)

# get dcs effect on final weights
#---------------------------------
variables=['dw_clopath']
w_clopath_df = df_funcs._get_dcs_effect(df=w_clopath_df, variables=variables,rerun=True )
filename=variable+'.pkl'
w_clopath_df.to_pickle(directory+filename)