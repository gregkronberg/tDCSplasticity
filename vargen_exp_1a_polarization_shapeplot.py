"""
analysis

"""
import vargen_functions

kwargs = {'experiment':'exp_1a_polarization_shapeplot'}
varfuncs=vargen_functions.VarFuncs()
df_funcs=vargen_functions.DfFuncs()
vargen_funcs=vargen_functions.VarGen()
directory = 'Data/'+kwargs['experiment']+'/'
variable='vtrace_df'
filename=variable+'.pkl'
functions=[varfuncs._get_vtrace]
kwlist=[{}]
rerun=[]#[varfuncs._get_vtrace]
keep=[]
file_limit=[]
df = vargen_funcs._vargen(variable=variable, directory=directory, functions=functions, kwlist=kwlist, rerun=rerun, keep=keep, file_limit=file_limit)

df = df_funcs._get_polarization(df)
df.to_pickle(directory+filename)
