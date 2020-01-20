"""
analysis

"""
import vargen_functions
import analysis
import time

kwargs = {
'experiment':'exp_5a1',
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
    kwargs['params']=['stdp_dt']
    kwlist=[kwargs]
    rerun=[]#[varfuncs._get_vtrace]
    keep=[]
    file_limit=[]
    vtrace_df = vargen_funcs._vargen(variable=variable, directory=directory, functions=functions, kwlist=kwlist, rerun=rerun, keep=keep, file_limit=file_limit)
#############################################################################
# create weight df from vtrace df stdp params
#############################################################################
runme=False
if runme:
    variable='w_clopath_stdp_df'
    start=time.time()
    w_clopath_stdp_df = df_funcs._get_w_clopath(vtrace_df=vtrace_df, clopath_param=vargen_functions.clopath_param_all_stdp, split_columns=['path_1_pulse_freq'])
    end=time.time()
    print 'timer:', (end-start)
    filename=variable+'.pkl'
    print 'saving', filename
    w_clopath_stdp_df.to_pickle(directory+filename)
#############################################################################
# create weight df from vtrace df stdp2 params
#############################################################################
runme=False
if runme:
    variable='w_clopath_stdp2_df'
    start=time.time()
    w_clopath_stdp2_df = df_funcs._get_w_clopath(vtrace_df=vtrace_df, clopath_param=vargen_functions.clopath_param_all_stdp2, split_columns=['path_1_pulse_freq'])
    end=time.time()
    print 'timer:', (end-start)
    filename=variable+'.pkl'
    print 'saving', filename
    w_clopath_stdp2_df.to_pickle(directory+filename)
#############################################################################
# create weight df from vtrace df dcsparams
#############################################################################
runme=False
if runme:
    variable='w_clopath_df'
    start=time.time()
    w_clopath_df = df_funcs._get_w_clopath(vtrace_df=vtrace_df, clopath_param=vargen_functions.clopath_param_all, split_columns=['path_1_pulse_freq'])
    end=time.time()
    print 'timer:', (end-start)
    filename=variable+'.pkl'
    print 'saving', filename
    w_clopath_df.to_pickle(directory+filename)


#############################################################################
# create weight df from vtrace df dcsparams
#############################################################################
runme=True
if runme:
    variable='w_clopath2_df'
    start=time.time()
    print vargen_functions.clopath_param_all2
    w_clopath2_df = df_funcs._get_w_clopath(vtrace_df=vtrace_df, clopath_param=vargen_functions.clopath_param_all2, split_columns=['path_1_pulse_freq'])
    end=time.time()
    print 'timer:', (end-start)
    filename=variable+'.pkl'
    print 'saving', filename
    w_clopath2_df.to_pickle(directory+filename)