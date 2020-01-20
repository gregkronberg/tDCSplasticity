"""
analysis

"""
import vargen_functions
import analysis
import time

kwargs = {
'experiment':'exp_1a_dose_basal',
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
vtrace_df = df_funcs._get_polarization(df=vtrace_df, time=9, rerun=False, file_limit=[])
variable='vtrace_df'
filename=variable+'.pkl'
print 'saving', filename, 'to', directory
vtrace_df.to_pickle(directory+filename)
#############################################################################
# spikes
#############################################################################
runme=False
if runme:
    variable='spikes_df'
    functions=[varfuncs._get_spikes]
    kwlist=[kwargs]
    rerun=[]
    keep=[]
    file_limit=1
    spikes_df = vargen_funcs._vargen(variable=variable, directory=directory, functions=functions, kwlist=kwlist, rerun=rerun, keep=keep, file_limit=file_limit)
#############################################################################
# weights
#############################################################################
# runme=False
# if runme:
#     variable='w_clopath_df'
#     functions=[varfuncs._get_w_clopath]
#     kwlist=[kwargs]
#     rerun=[]
#     keep=[]
#     file_limit=[]
#     w_clopath_df = vargen_funcs._vargen(variable=variable, directory=directory, functions=functions, kwlist=kwlist, rerun=rerun, keep=keep, file_limit=file_limit)

#     # get dcs effect on final weights
#     #---------------------------------
#     runme=False
#     if runme:
#         variables=['dw_clopath']
#         w_clopath_df = df_funcs._get_dcs_effect(df=w_clopath_df, variables=variables,rerun=True )
#         filename=variable+'.pkl'
#         w_clopath_df.to_pickle(directory+filename)
# #############################################################################
# # weights w stdp params
# #############################################################################
# runme=False
# if runme:
#     variable='w_clopath_stdp_df'
#     functions=[varfuncs._get_w_clopath]
#     kwargs['clopath_param']=vargen_functions.clopath_param_all_stdp
#     kwlist=[kwargs]
#     rerun=[]
#     keep=[]
#     file_limit=[]
#     w_clopath_stdp_df = vargen_funcs._vargen(variable=variable, directory=directory, functions=functions, kwlist=kwlist, rerun=rerun, keep=keep, file_limit=file_limit, write_protocol='.pkl')

#     # get dcs effect on final weights
#     #-------------------------------
#     runme=False
#     if runme:
#         variable='w_clopath_stdp_df'
#         variables=['dw_clopath']
#         w_clopath_stdp_df = df_funcs._get_dcs_effect(df=w_clopath_stdp_df, variables=variables,rerun=True )
#         filename=variable+'.pkl'
#         print 'saving', filename
#         w_clopath_stdp_df.to_pickle(directory+filename)
# #############################################################################
# # create weight df from vtrace df stdp2
# #############################################################################
# runme=False
# if runme:
#     variable='w_clopath_stdp2_df'
#     start=time.time()
#     w_clopath_stdp2_df = df_funcs._get_w_clopath(vtrace_df=vtrace_df, clopath_param=vargen_functions.clopath_param_all_stdp2)
#     end=time.time()
#     print 'timer:', (end-start)
#     filename=variable+'.pkl'
#     print 'saving', filename
#     w_clopath_stdp2_df.to_pickle(directory+filename)

# # get dcs effect on final weights
# #---------------------------------
# runme=False
# if runme:
#     variable='w_clopath_stdp2_df'
#     variables=['dw_clopath']
#     start=time.time()
#     w_clopath_stdp2_df = df_funcs._get_dcs_effect(df=w_clopath_stdp2_df, variables=variables,rerun=True )
#     end=time.time()
#     print 'timer:', (end-start)
#     filename=variable+'.pkl'
#     print 'saving', filename
#     w_clopath_stdp2_df.to_pickle(directory+filename)
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
#############################################################################
# create weight df from vtrace df stdp2
#############################################################################
runme=False
if runme:
    variable='w_clopath_stdp2_df'
    start=time.time()
    w_clopath_stdp2_df = df_funcs._get_w_clopath(vtrace_df=vtrace_df, clopath_param=vargen_functions.clopath_param_all_stdp2)
    end=time.time()
    print 'timer:', (end-start)
    filename=variable+'.pkl'
    print 'saving', filename
    w_clopath_stdp2_df.to_pickle(directory+filename)

# get dcs effect on final weights
#---------------------------------
runme=False
if runme:
    variable='w_clopath_stdp2_df'
    variables=['dw_clopath']
    start=time.time()
    w_clopath_stdp2_df = df_funcs._get_dcs_effect(df=w_clopath_stdp2_df, variables=variables,rerun=True )
    end=time.time()
    print 'timer:', (end-start)
    filename=variable+'.pkl'
    print 'saving', filename
    w_clopath_stdp2_df.to_pickle(directory+filename)
#############################################################################
# create weight df from vtrace df stdp3
#############################################################################
runme=False
if runme:
    variable='w_clopath2_df'
    start=time.time()
    w_clopath2_df = df_funcs._get_w_clopath(vtrace_df=vtrace_df, clopath_param=vargen_functions.clopath_param_all2)
    end=time.time()
    print 'timer:', (end-start)
    filename=variable+'.pkl'
    print 'saving', filename
    w_clopath2_df.to_pickle(directory+filename)

# get dcs effect on final weights
#---------------------------------
runme=False
if runme:
    variable='w_clopath2_df'
    variables=['dw_clopath']
    start=time.time()
    w_clopath2_df = df_funcs._get_dcs_effect(df=w_clopath2_df, variables=variables,rerun=True )
    end=time.time()
    print 'timer:', (end-start)
    filename=variable+'.pkl'
    print 'saving', filename
    w_clopath2_df.to_pickle(directory+filename)