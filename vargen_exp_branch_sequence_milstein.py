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
