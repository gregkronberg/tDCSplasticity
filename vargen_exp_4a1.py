"""
analysis

"""
import vargen_functions
import analysis
import time

kwargs = {
'experiment':'exp_4a1',
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
runme=True
if runme:
    variable='vtrace_df'
    filename=variable+'.pkl'
    functions=[varfuncs._get_vtrace]
    kwlist=[kwargs]
    rerun=[]#[varfuncs._get_vtrace]
    keep=[]
    file_limit=[]
    vtrace_df = vargen_funcs._vargen(variable=variable, directory=directory, functions=functions, kwlist=kwlist, rerun=rerun, keep=keep, file_limit=file_limit)

# get membrane polarization before synaptic inputs arrive
#----------------------------------------------------------
# vtrace_df = df_funcs._get_polarization(df=vtrace_df, time=9.9, rerun=True, file_limit=[])
# print 'saving', filename, 'to', directory
# vtrace_df.to_pickle(directory+filename)
#############################################################################
# spikes
#############################################################################
runme=True
if runme:
    variable='spikes_df'
    start=time.time()
    spikes_df = df_funcs._get_spikes(vtrace_df=vtrace_df, threshold=-40)
    end=time.time()
    print 'timer:', (end-start)
    filename=variable+'.pkl'
    print 'saving', filename
    spikes_df.to_pickle(directory+filename)

    start=time.time()
    spikes_df = df_funcs._get_xcorr_soma(spikes_df=spikes_df, threshold=-40, file_limit=[], )
    end=time.time()
    print 'timer:', (end-start)
    filename=variable+'.pkl'
    print 'saving', filename
    spikes_df.to_pickle(directory+filename)
#############################################################################
# create weight df from vtrace df dcs parameters
#############################################################################
runme=False
if runme:
    variable='w_clopath_df'
    start=time.time()
    w_clopath_df = df_funcs._get_w_clopath(vtrace_df=vtrace_df, clopath_param=vargen_functions.clopath_param_all)
    end=time.time()
    print 'timer:', (end-start)
    filename=variable+'.pkl'
    print 'saving', filename
    w_clopath_df.to_pickle(directory+filename)

# get dcs effect on final weights
#---------------------------------
runme=False
if runme:
    variable='w_clopath_df'
    variables=['dw_clopath']
    start=time.time()
    w_clopath_df = df_funcs._get_dcs_effect(df=w_clopath_df, variables=variables,rerun=True )
    end=time.time()
    print 'timer:', (end-start)
    filename=variable+'.pkl'
    print 'saving', filename
    w_clopath_df.to_pickle(directory+filename)