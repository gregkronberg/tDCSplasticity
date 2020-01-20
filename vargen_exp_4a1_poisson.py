"""
analysis

"""
import vargen_functions
import analysis
import time

kwargs = {
'experiment':'exp_4a1_poisson',
'clopath_param':vargen_functions.clopath_param_all_temp,
'threshold':-30,
}
directory = 'Data/'+kwargs['experiment']+'/'
varfuncs=vargen_functions.VarFuncs()
df_funcs=vargen_functions.DfFuncs()
vargen_funcs=vargen_functions.VarGen()

#############################################################################
# vtrace
#############################################################################
def _vtrace(vargen_funcs):
    variable='vtrace_df'
    filename=variable+'.pkl'
    functions=[varfuncs._get_vtrace]
    kwlist=[kwargs]
    rerun=[]#[varfuncs._get_vtrace]
    keep=[]
    file_limit=[]
    vargen_kwargs={'search_string':'data'}
    vtrace_df = vargen_funcs._vargen(variable=variable, directory=directory, functions=functions, kwlist=kwlist, rerun=rerun, keep=keep, file_limit=file_limit, **vargen_kwargs)
    return vtrace_df

vtrace_df = _vtrace(vargen_funcs=vargen_funcs)

#############################################################################
# spikes
#############################################################################
def _spikes(vtrace_df, df_funcs):
    '''
    '''
    variable='spikes_df'
    start=time.time()
    spikes_df = df_funcs._get_spikes(vtrace_df=vtrace_df, threshold=-40)
    end=time.time()
    print 'timer:', (end-start)
    filename=variable+'.pkl'
    print 'saving', filename
    spikes_df=analysis._save_group_data(df=spikes_df, directory=directory, variable='spikes_df', extension='.pkl', )
    # spikes_df.to_pickle(directory+filename)
    return spikes_df

spikes_df = _spikes(vtrace_df=vtrace_df, df_funcs=df_funcs)

def _spikes_xcorr(spikes_df, df_funcs):
    '''
    '''
    start=time.time()
    spikes_df = df_funcs._get_xcorr_soma(spikes_df=spikes_df, threshold=-40, file_limit=[], )
    end=time.time()
    print 'timer:', (end-start)
    filename=variable+'.pkl'
    print 'saving', filename
    # spikes_df.to_pickle(directory+filename)
    spikes_df=analysis._save_group_data(df=spikes_df, directory=directory, variable='spikes_df', extension='.pkl', )
    return spikes_df

# spikes_df= _spikes_xcorr(spikes_df=spikes_df, df_funcs=df_funcs)

    
#############################################################################
# create weight df from vtrace df dcs parameters
#############################################################################
def _w_clopath(vtrace_df, vargen_functions):
    variable='w_clopath_df'
    start=time.time()
    w_clopath_df = df_funcs._get_w_clopath(vtrace_df=vtrace_df, clopath_param=vargen_functions.clopath_param_all_temp)
    end=time.time()
    print 'timer:', (end-start)
    filename=variable+'.pkl'
    print 'saving', filename
    # w_clopath_df.to_pickle(directory+filename)
    w_clopath_df=analysis._save_group_data(df=w_clopath_df, directory=directory, variable='w_clopath_df', extension='.pkl',)
    return w_clopath_df

w_clopath_df=_w_clopath(vtrace_df=vtrace_df, vargen_functions=vargen_functions)

# get dcs effect on final weights
#---------------------------------
def _get_dcs_effect_on_weights(w_clopath_df ):
    variable='w_clopath_df'
    variables=['dw_clopath']
    start=time.time()
    w_clopath_df = df_funcs._get_dcs_effect(df=w_clopath_df, variables=variables,rerun=True )
    end=time.time()
    print 'timer:', (end-start)
    filename=variable+'.pkl'
    print 'saving', filename
    # w_clopath_df.to_pickle(directory+filename)
    w_clopath_df=analysis._save_group_data(df=w_clopath_df, directory=directory, variable='w_clopath_df', extension='.pkl', size_limit=2E9)


    return w_clopath_df

w_clopath_df=_get_dcs_effect_on_weights(w_clopath_df=w_clopath_df,)