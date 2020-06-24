"""
docstring
"""
# imports
from neuron import h
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import itertools as it
import stims
import pickle
import param
import os
import copy
import pandas as pd
import functions as fncs


class Run():
    '''
    Arguments list:
    p - dictionary of parameters

    each experiment appends a list to the appropriate key in the data dictionary
    data is organized as data['type of data'][experiments][sections][time series data vector]
    details of each experiment are tracked via the data['detail'][experiment number], e.g. data['field'][2]
    '''
    def __init__(self):
        '''
        ==Args==
        ==Out==
        ==Update==
        ==Comments==
        '''
        pass

    def _standard_run(self, p, cell, **kwargs):
        '''
        '''
        # update clopath parameters
        self._update_clopath( p=p, syns=cell.syns)
        # activate synapses
        method='_bursts'
        if 'uncage_method' in kwargs:
            method=kwargs['uncage_method']
        self.stim, self.nc = self._activate_synapses_uncage(p_paths=p['p_path'], cell=cell, method=method)
        # set up recording vectors
        self.rec = self._setup_recording_vectors_df(p=p, geo=cell.geo, syns=cell.syns, nc=self.nc)
        # load hoc standard run environment
        h.load_file("stdrun.hoc")
        # run simulation and record data
        self.data = self._run_sims(p=p, rec=self.rec)
        # clear synapses
        self.nc = []
        self.stim=[]

    def _standard_run_df(self, p, cell, **kwargs):
        '''
        '''
        # update clopath parameters
        self._update_clopath( p=p, syns=cell.syns)
        # activate synapses
        method='_bursts'
        if 'uncage_method' in kwargs:
            method=kwargs['uncage_method']
        self.stim, self.nc = self._activate_synapses_uncage(p_paths=p['p_path'], cell=cell, method=method)
        # set up recording vectors
        self.rec = self._setup_recording_vectors_df(p=p, geo=cell.geo, syns=cell.syns, nc=self.nc)
        # load hoc standard run environment
        h.load_file("stdrun.hoc")
        # run simulation and record data
        self.data = self._run_sims_df(p=p, rec=self.rec)
        # clear synapses
        self.nc = []
        self.stim=[]
    
    # update clopath parameters
    def _update_clopath(self, p, syns):
        '''
        ==Args==
        ==Out==
        ==Update==
        ==Comments==
        '''
        # iterate over parameters
        for parameter_key, parameter in p.iteritems():
            # if it is a clopath learning rule parameter
            if 'clopath_' in parameter_key:
                # get parameter name
                p_clopath = parameter_key[len('clopath_'):]

                for tree_key, tree in syns.iteritems():
                    # iterate over sections
                    for sec_i,sec in enumerate(tree):
                        # iterate over segments
                        for seg_i,seg in enumerate(sec):
                            # if segment contains a clopath synapse
                            if 'clopath' in list(seg.keys()): 
                                # set synapse parameter value
                                setattr(seg['clopath'], p_clopath, p['clopath_'+p_clopath])
    
    def _activate_synapses_uncage(self, p_paths,cell, method='_bursts'):
        """ create netstim objects and connect them to synapses

        ==Args==
        -p_paths : parameter dictionary for each synaptic pathway as {path}{parameter:value}
        -method  : string specifing the method for creating NetStim objects.  will call the corresponding method in stims.Uncage

        ==Out==
        -stims : NetStim objects as {path}[synapse number][burst number]
        -nc    : NetCon objects as {path}[synapse number][burst number]

        ==Updates==
        -nc and stims are stored as attributes of the Run instance

        ==Comments==
        """
        # store NetStim objects
        self.stim={}
        # stor NetCon objects
        self.nc={}
        # iterate over paths
        for path_key, path in p_paths.iteritems():
            # create uncage object for activating synapses
            uncage=stims.Uncage()
            # get method for creating NetStim objects
            stim_method = getattr(uncage, method)
            # create NetStims
            self.stim[path_key] = stim_method(path)
            # connect NetStims to synapses with NetCons
            self.nc[path_key] = uncage._connect_synapses(p_path=path, stim=self.stim[path_key], syns=cell.syns)

        return self.stim, self.nc

    def _setup_recording_vectors(self, p, geo, syns, nc):
        ''' set up recording vectors for simulation
        ===Args===
        -p : global parameter dictionary, must contain:
        -----rec_idx : list of tuples [(tree, sec, seg)]
        -----rec_variables :list of variables to record [(variable, variable type, mechanism)]
        -geometry : geo structure

        ===Out===
        -rec : dictionary of recorded variables
                -rec{data{variable}, location, t, p}
                    -if variable is 'input_times', data is a list of 1d arrays containing input times for each burst in each synapse (there can be multiple repeats for each synapse location and each array can be different length depending on the input to each synapse)
                -t and p are single entries, not lists.  They should be the same for all segments
                -location is a list of all recorded segments as [(tree, section, segment)], specified by p['rec_idx']
                -data{variable} is a list of hoc recording vectors 

        ===Updates===
        -hoc recording vectors are created for each segment specified in locations

        ==Comments==
        -for the 'input_times' variable, a list of the corresponding synapse type (eg ampa, nmda) is kept in rec{variable}{'syn_types'}
        -note that each locations list can contain redundant entries due to synapses occuring in multiple pathways
        '''
        print 'setting up recording vectors'
        rec={}
        for variable, variable_type, mechanism in p['rec_variables']:
            if variable=='t':
                continue
            if variable not in rec:
                
                rec[variable]={
                'data':[],
                'locations':[],
                'p':p,
                't':h.Vector()}
                # record time
                rec[variable]['t'].record(h._ref_t)

                if variable == 'input_times':
                    rec[variable]['syn_types']=[]

            #iterate over locations to record from
            for location_i, location in enumerate(p['rec_idx']):
                tree_key,sec_num,seg_num = location
                seg_loc = float(seg_num+1)/(geo[tree_key][sec_num].nseg+1)

                if variable_type == 'range' and  mechanism in dir(geo[tree_key][sec_num](seg_loc)):

                    # point to variable for recording
                    var_rec = getattr(geo[tree_key][sec_num](seg_loc), '_ref_'+variable)
                    
                    # create recording vector
                    rec[variable]['data'].append(h.Vector())
                    rec[variable]['data'][-1].record(var_rec)
                    rec[variable]['locations'].append(location)

                if variable_type == 'syn' and mechanism in syns[tree_key][sec_num][seg_num].keys() and  variable in dir(syns[tree_key][sec_num][seg_num][mechanism]): 
                                
                    # point to variable to record
                    var_rec = getattr(syns[tree_key][sec_num][seg_num][mechanism], '_ref_'+variable)

                    # create recording vector
                    rec[variable]['data'].append(h.Vector())
                    rec[variable]['data'][-1].record(var_rec)
                    rec[variable]['locations'].append(location)

                # synapse input times
                if variable== 'input_times' and variable_type == 'syn':

                    # iterate over synapse pathways
                    for path_key, path in p['p_path'].iteritems():

                        # if the current location is in the current pathway
                        if location in path['syn_idx']:

                            # get index of the location in the pathway
                            seg_i = path['syn_idx'].index(location)

                            # iterate over bursts 
                            for burst_i, burst in enumerate(nc[path_key][seg_i][mechanism]):
                                
                                # add recording vector
                                rec[variable]['data'].append(h.Vector())

                                # set recording from netcon object
                                nc[path_key][seg_i][mechanism][burst_i].record(rec[variable]['data'][-1])
                                
                                # add location
                                rec[variable]['locations'].append(location)
                                
                                # add synapse type
                                rec[variable]['syn_types'].append(mechanism)

        return rec
    
    def _setup_recording_vectors_df(self, p, geo, syns, nc):
        ''' set up recording vectors for simulation
        ===Args===
        -p : global parameter dictionary, must contain:
        -----rec_idx : list of tuples [(tree, sec, seg)]
        -----rec_variables :list of variables to record [(variable, variable type, mechanism)]
        -geometry : geo structure

        ===Out===
        -rec : dictionary of recorded variables
                -rec{data{variable}, location, t, p}
                    -if variable is 'input_times', data is a list of 1d arrays containing input times for each burst in each synapse (there can be multiple repeats for each synapse location and each array can be different length depending on the input to each synapse)
                -t and p are single entries, not lists.  They should be the same for all segments
                -location is a list of all recorded segments as [(tree, section, segment)], specified by p['rec_idx']
                -data{variable} is a list of hoc recording vectors 

        ===Updates===
        -hoc recording vectors are created for each segment specified in locations

        ==Comments==
        -for the 'input_times' variable, a list of the corresponding synapse type (eg ampa, nmda) is kept in rec{variable}{'syn_types'}
        -note that each locations list can contain redundant entries due to synapses occuring in multiple pathways
        '''
        print 'setting up recording vectors'
        rec=pd.DataFrame(dtype=object)
        init_columns = ['location', 'tree_key', 'sec_num', 'seg_num', 'seg_loc','data_t', 'seg_dist', 'morpho', 'n_paths','path', 'paths', 'w', 'mechanism']
        rec = fncs._initialize_column(rec, init_columns)
        # rec['location']=None
        # rec['syn_type']=None
        # rec['tree_key']=None
        # rec['sec_num']=None
        # rec['seg_num']=None
        # rec['p']=None
        # rec['t']=None

        

        #iterate over locations to record from
        for i, location in enumerate(p['rec_idx']):
            # get location info
            #--------------------------------------------------------
            tree_key, sec_num, seg_num = location
            seg_loc = float(seg_num+1)/(geo[tree_key][sec_num].nseg+1)
            rec.at[i, 'location']=location
            rec.at[i, 'tree_key']=tree_key
            rec.at[i, 'sec_num']=sec_num
            rec.at[i, 'seg_num']=seg_num
            rec.at[i, 'seg_loc']=seg_loc

            #update parameters
            #-----------------------------------------------------------
            for param, value in p.iteritems():
                # initialize column
                rec = fncs._initialize_column(rec, param)
                # empty lists to nan
                if type(value)==list or type(value)==dict:
                    if len(value)==0:
                        value=np.nan
                    rec.at[i,param] = value
                else:
                    rec.at[i,param] = value

            # update morphology info
            #--------------------------------------------------------------
            # get distance from soma
            dist = p['seg_dist'][tree_key][sec_num][seg_num]
            rec.at[i, 'seg_dist']=dist
            # get morphology entry (overall seg index, location key, x, y, z, dimaeter, parent overall index)
            morpho = p['morpho'][tree_key]
            morpho_key = '_'.join([tree_key, str(sec_num), str(seg_num)])
            morpho_val = [seg for sec in morpho for seg in sec if morpho_key in seg][0]
            # add morphology info to rec
            rec.at[i, 'morpho']=morpho_val

            # synaptic input path parameters
            #--------------------------------------------------------------
            # keep track of paths that the current segment belongs to
            paths=[]
            # total number of paths
            n_paths=len(p['p_path'].keys())
            rec.at[i, 'n_paths']=n_paths
            # iterate over input pathways
            for path_key, path in p['p_path'].iteritems():
                # list of paths that contain the current location
                if location in path['syn_idx']:
                    temp_i = path['syn_idx'].index(location)
                    paths.append(path_key)
                    # synaptic weight at current location
                    if 'w_idx' in path:
                        rec.at[i, 'w'] = path['w_idx'][temp_i]
                    rec.at[i,'path']=path_key
                # iterate over path parameters
                for path_param, param_val in path.iteritems(): 
                    colname = '_'.join(['path', path_key, path_param])
                    # initaialize column
                    rec = fncs._initialize_column(rec, colname)
                    # add all path parameters to df
                    #------------------------------
                    rec.at[i,colname]=param_val
            # add paths that current entry belongs to 
            rec.at[i,'paths']=paths

            # time vectors
            #----------------------------------------------------------------
            rec.at[i, 'data_t'] = h.Vector()
            rec.loc[i,'data_t'].record(h._ref_t)

            # update recording variables
            #----------------------------------------------------------------
            for variable, variable_type, mechanism in p['rec_variables']:
                # skip time vector
                if variable=='t':
                    continue
                # store mechanism
                rec = fncs._initialize_column(rec, 'mechanism_'+variable)
                rec.at[i, 'mechanism_'+variable]=mechanism
                # initialize columns
                colname = 'data_'+variable
                rec = fncs._initialize_column(rec, colname)

                # range variables
                #-----------------
                if variable_type == 'range' and  mechanism in dir(geo[tree_key][sec_num](seg_loc)):
                    # point to variable for recording
                    var_rec = getattr(geo[tree_key][sec_num](seg_loc), '_ref_'+variable)
                    # create recording vector
                    rec.at[i, colname] = h.Vector()
                    rec.loc[i, colname].record(var_rec)

                # synaptic variables
                #---------------------
                if variable_type == 'syn' and mechanism in syns[tree_key][sec_num][seg_num].keys() and  variable in dir(syns[tree_key][sec_num][seg_num][mechanism]): 
                                
                    # point to variable to record
                    var_rec = getattr(syns[tree_key][sec_num][seg_num][mechanism], '_ref_'+variable)

                    # create recording vector
                    # create recording vector
                    rec.at[i, colname] = h.Vector()
                    rec.loc[i, colname].record(var_rec)

                # synapse input times
                #----------------------
                if variable== 'input_times' and variable_type == 'syn':

                    # iterate over synapse pathways
                    for path_key, path in p['p_path'].iteritems():

                        # if the current location is in the current pathway
                        if location in path['syn_idx']:

                            # get index of the location in the pathway
                            seg_i = path['syn_idx'].index(location)

                            # iterate over bursts 
                            for burst_i, burst in enumerate(nc[path_key][seg_i][mechanism]):
                                
                                # add recording vector
                                rec.at[i, colname]=h.Vector()

                                # set recording from netcon object
                                nc[path_key][seg_i][mechanism][burst_i].record(rec.loc[i,colname])


        return rec

    def _setup_recording_vectors_milstein(self, p, cell, syns, nc):
        ''' set up recording vectors for simulation
        ===Args===
        -p : global parameter dictionary, must contain:
        -----rec_idx : list of tuples [(tree, sec, seg)]
        -----rec_variables :list of variables to record [(variable, variable type, mechanism)]
        -geometry : geo structure

        ===Out===
        -rec : dictionary of recorded variables
                -rec{data{variable}, location, t, p}
                    -if variable is 'input_times', data is a list of 1d arrays containing input times for each burst in each synapse (there can be multiple repeats for each synapse location and each array can be different length depending on the input to each synapse)
                -t and p are single entries, not lists.  They should be the same for all segments
                -location is a list of all recorded segments as [(tree, section, segment)], specified by p['rec_idx']
                -data{variable} is a list of hoc recording vectors 

        ===Updates===
        -hoc recording vectors are created for each segment specified in locations

        ==Comments==
        -for the 'input_times' variable, a list of the corresponding synapse type (eg ampa, nmda) is kept in rec{variable}{'syn_types'}
        -note that each locations list can contain redundant entries due to synapses occuring in multiple pathways
        '''
        print 'setting up recording vectors'
        rec=pd.DataFrame(dtype=object)
        init_columns = ['location', 'tree_key', 'sec_num', 'seg_num', 'seg_loc','data_t', 'seg_dist', 'morpho', 'n_paths','path', 'paths', 'w', 'mechanism']
        rec = fncs._initialize_column(rec, init_columns)
        # rec['location']=None
        # rec['syn_type']=None
        # rec['tree_key']=None
        # rec['sec_num']=None
        # rec['seg_num']=None
        # rec['p']=None
        # rec['t']=None

        #iterate over locations to record from
        for i, location in enumerate(p['rec_idx']):
            # get location info
            #--------------------------------------------------------
            tree_key, sec_num, seg_num, seg_loc = location
            # seg_loc = float(seg_num+1)/(geo[tree_key][sec_num].nseg+1)
            rec.at[i, 'location']=location
            rec.at[i, 'tree_key']=tree_key
            rec.at[i, 'sec_num']=sec_num
            rec.at[i, 'seg_num']=seg_num
            rec.at[i, 'seg_loc']=seg_loc

            #update parameters
            #-----------------------------------------------------------
            for param, value in p.iteritems():
                # initialize column
                rec = fncs._initialize_column(rec, param)
                # empty lists to nan
                if type(value)==list or type(value)==dict:
                    if len(value)==0:
                        value=np.nan
                    rec.at[i,param] = value
                else:
                    rec.at[i,param] = value

            # update morphology info
            #--------------------------------------------------------------
            # get distance from soma
            dist = p['seg_dist'][tree_key][sec_num][seg_num]
            rec.at[i, 'seg_dist']=dist
            # get morphology entry (overall seg index, location key, x, y, z, dimaeter, parent overall index)
            morpho = p['morpho'][tree_key]
            morpho_key = '_'.join([tree_key, str(sec_num), str(seg_num)])
            morpho_val = [seg for sec in morpho for seg in sec if morpho_key in seg][0]
            # add morphology info to rec
            rec.at[i, 'morpho']=morpho_val

            # synaptic input path parameters
            #--------------------------------------------------------------
            # keep track of paths that the current segment belongs to
            paths=[]
            # total number of paths
            n_paths=len(p['p_path'].keys())
            rec.at[i, 'n_paths']=n_paths
            # iterate over input pathways
            for path_key, path in p['p_path'].iteritems():
                # list of paths that contain the current location
                if location in path['syn_idx']:
                    temp_i = path['syn_idx'].index(location)
                    paths.append(path_key)
                    # synaptic weight at current location
                    if 'w_idx' in path:
                        rec.at[i, 'w'] = path['w_idx'][temp_i]
                    rec.at[i,'path']=path_key
                # iterate over path parameters
                for path_param, param_val in path.iteritems(): 
                    colname = '_'.join(['path', path_key, path_param])
                    # initaialize column
                    rec = fncs._initialize_column(rec, colname)
                    # add all path parameters to df
                    #------------------------------
                    rec.at[i,colname]=param_val
            # add paths that current entry belongs to 
            rec.at[i,'paths']=paths

            # time vectors
            #----------------------------------------------------------------
            rec.at[i, 'data_t'] = h.Vector()
            rec.loc[i,'data_t'].record(h._ref_t)

            # update recording variables
            #----------------------------------------------------------------
            for variable, variable_type, mechanism in p['rec_variables']:
                # skip time vector
                if variable=='t':
                    continue
                # store mechanism
                rec = fncs._initialize_column(rec, 'mechanism_'+variable)
                rec.at[i, 'mechanism_'+variable]=mechanism
                # initialize columns
                colname = 'data_'+variable
                rec = fncs._initialize_column(rec, colname)

                # range variables
                #-----------------
                if variable_type == 'range' and  mechanism in dir(geo[tree_key][sec_num](seg_loc)):
                    # point to variable for recording
                    var_rec = getattr(geo[tree_key][sec_num](seg_loc), '_ref_'+variable)
                    # create recording vector
                    rec.at[i, colname] = h.Vector()
                    rec.loc[i, colname].record(var_rec)

                # synaptic variables
                #---------------------
                if variable_type == 'syn' and mechanism in syns[tree_key][sec_num][seg_num].keys() and  variable in dir(syns[tree_key][sec_num][seg_num][mechanism]): 
                                
                    # point to variable to record
                    var_rec = getattr(syns[tree_key][sec_num][seg_num][mechanism], '_ref_'+variable)

                    # create recording vector
                    # create recording vector
                    rec.at[i, colname] = h.Vector()
                    rec.loc[i, colname].record(var_rec)

                # synapse input times
                #----------------------
                if variable== 'input_times' and variable_type == 'syn':

                    # iterate over synapse pathways
                    for path_key, path in p['p_path'].iteritems():

                        # if the current location is in the current pathway
                        if location in path['syn_idx']:

                            # get index of the location in the pathway
                            seg_i = path['syn_idx'].index(location)

                            # iterate over bursts 
                            for burst_i, burst in enumerate(nc[path_key][seg_i][mechanism]):
                                
                                # add recording vector
                                rec.at[i, colname]=h.Vector()

                                # set recording from netcon object
                                nc[path_key][seg_i][mechanism][burst_i].record(rec.loc[i,colname])


        return rec
    
    def _store_recording_vectors_to_data_df(self, data, rec, add_col=['field'], add_col_vals=[], data_key='data_'):
        '''
        '''
        data_temp = copy.deepcopy(rec)

        # add specified columns
        if type(add_col)!=list and type(add_col)!=tuple:
            add_col=[add_col]
        if type(add_col_vals)!=list and type(add_col_vals)!=tuple:
            add_col_vals=[add_col_vals]
        for col_i, col in enumerate(add_col):
            data_temp[col]=add_col_vals[col_i]

        # convert hoc recording vectors to np.array
        for column in data_temp.columns:
            if data_key in column:
                data_temp[column] = data_temp[column].apply(np.array)
        # append to overall data df
        data = data.append(data_temp, ignore_index=True)
        return data

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

    def _run_sims(self, p, rec):
        ''' run simulations, grouping simulations with different applied electric fields
        
        ==Args==
        -p : global parameter dictionary
        -rec : structure containing hoc recording vectors for segments that are specified by p['rec_idx']
                    -rec{variable}{data, location, field, t, p,}[segment number]
                    -t and p are single entries, not lists.  They should be the same for all segments

        ==Out==
        -rec : dictionary of recorded variables
                -data{data{variable}, location, field, t, p}
                    -data{variable} : segments x samples array of time series data for the specified variable.  First dimension matches location and field lists
                            -if variable is input_times, data is a list of 1d arrays containing input times for each burst in each synapse (there can be multiple repeats for each synapse location and each array can be different length depending on the input to each synapse)
                    -location : list of recorded segment locations as [(tree, section, segment)]
                    -field : list of field polarities/intensities (negative:cathodal/soma hyperpolarizing, positive:anodal/soma depolarizing) 
                    -t : 1D array of time vector that is shared for simulations in group
                    -p : single parameter dictionary for all simulations in group

        ==Updates==
        ==Comments==
        -For repeated simulations, hoc recording vectors are overwritten.  To save the data for each simulation, the hoc vectors are copied from the structure rec to the structure data as np arrays
        -for the 'input_times' variable, a list of the corresponding synapse type (eg ampa, nmda) is kept in rec{variable}{'syn_types'}
        -note that each locations list can contain redundant entries due to synapses occuring in multiple pathways
        -
        '''

        nsamples = int(p['tstop']/p['dt'])+1
        data={}
        # preallocate data arrays?
        for variable_key, variable in rec.iteritems():
            if variable_key not in data:

                data[variable_key]={
                'data':np.zeros( (len(p['field'])*len(variable['data']),nsamples)),
                'locations':[],
                'field':[],
                'trial_id':[],
                'p':variable['p'],
                't':variable['t']
                }

        # iterate over field polarities/intensities
        for field_i, field in enumerate(p['field']):

            # insert extracellular field
            dcs = stims.DCS(cell=0, field_angle=p['field_angle'], intensity=field, field_on=p['field_on'], field_off=p['field_off'],)

            if 'ac_field' in p:
                acs = stims.ACS(cell=0, p=p)

            # run time
            h.dt = p['dt']
            h.tstop = p['tstop']
            h.celsius= p['celsius']

            # load standard run environment
            
            print 'running simulation, electric field:', field
            
            # initialize voltage
            h.v_init=p['v_init'] # v_init is an object created in h when stdrun.hoc is called
            
            # run simulation
            h.run()
            print 'simulation finished'

            
            for variable_key, variable in rec.iteritems():

                # handle input times for each synapse. note each list entry can have a different number of input times
                if variable_key == 'input_times':
                    if not isinstance(data[variable_key]['data'],list):
                        data[variable_key]['data']=[]
                    if 'syn_types' not in data[variable_key]:
                        data[variable_key]['syn_types']=[]
                    for loc_i, loc in enumerate(variable['data']):
                        # print loc
                        data[variable_key]['data'].append(np.array(loc))
                    data[variable_key]['field'] += [field for seg_i in variable['locations']] 
                    data[variable_key]['trial_id'] += [p['trial_id'] for seg_i in variable['locations']]
                    data[variable_key]['locations'] += copy.copy(variable['locations'])
                    data[variable_key]['syn_types']+=copy.copy(variable['syn_types'])
                    data[variable_key]['t'] = np.asarray(data[variable_key]['t'])
                else:
                    nseg = len(variable['locations'])
                    print variable_key
                    data[variable_key]['data'][field_i*nseg:(field_i+1)*nseg, :] = np.array(variable['data'])

                    data[variable_key]['field'] += [field for seg_i in variable['locations']] 
                    data[variable_key]['trial_id'] += [p['trial_id'] for seg_i in variable['locations']]
                    data[variable_key]['locations'] += copy.copy(variable['locations'])
                    data[variable_key]['t'] = np.asarray(data[variable_key]['t'])


        return data

    def _run_sims_df(self, p, rec):
        ''' run simulations, grouping simulations with different applied electric fields
        
        ==Args==
        -p : global parameter dictionary
        -rec : structure containing hoc recording vectors for segments that are specified by p['rec_idx']
                    -rec{variable}{data, location, field, t, p,}[segment number]
                    -t and p are single entries, not lists.  They should be the same for all segments

        ==Out==
        -rec : dictionary of recorded variables
                -data{data{variable}, location, field, t, p}
                    -data{variable} : segments x samples array of time series data for the specified variable.  First dimension matches location and field lists
                            -if variable is input_times, data is a list of 1d arrays containing input times for each burst in each synapse (there can be multiple repeats for each synapse location and each array can be different length depending on the input to each synapse)
                    -location : list of recorded segment locations as [(tree, section, segment)]
                    -field : list of field polarities/intensities (negative:cathodal/soma hyperpolarizing, positive:anodal/soma depolarizing) 
                    -t : 1D array of time vector that is shared for simulations in group
                    -p : single parameter dictionary for all simulations in group

        ==Updates==
        ==Comments==
        -For repeated simulations, hoc recording vectors are overwritten.  To save the data for each simulation, the hoc vectors are copied from the structure rec to the structure data as np arrays
        -for the 'input_times' variable, a list of the corresponding synapse type (eg ampa, nmda) is kept in rec{variable}{'syn_types'}
        -note that each locations list can contain redundant entries due to synapses occuring in multiple pathways
        -
        '''

        nsamples = int(p['tstop']/p['dt'])+1
        data={}
        data=pd.DataFrame(dtype=object)
        # preallocate data arrays?
        # for variable_key, variable in rec.iteritems():
        #     if variable_key not in data:

        #         data[variable_key]={
        #         'data':np.zeros( (len(p['field'])*len(variable['data']),nsamples)),
        #         'locations':[],
        #         'field':[],
        #         'trial_id':[],
        #         'p':variable['p'],
        #         't':variable['t']
                # }
        
        
        data['field']=None
        data=pd.DataFrame(dtype=object)
        # iterate over field polarities/intensities
        for field_i, field in enumerate(p['field']):

            # insert extracellular field
            dcs = stims.DCS(cell=0, field_angle=p['field_angle'], intensity=field, field_on=p['field_on'], field_off=p['field_off'],)

            if 'ac_field' in p:
                acs = stims.ACS(cell=0, p=p)

            # run time
            h.dt = p['dt']
            h.tstop = p['tstop']
            h.celsius= p['celsius']

            # load standard run environment
            
            print 'running simulation, electric field:', field
            
            # initialize voltage
            h.v_init=p['v_init'] # v_init is an object created in h when stdrun.hoc is called
            
            # run simulation
            h.run()
            print 'simulation finished'

            data = self._store_recording_vectors_to_data_df(data=data, rec=rec, add_col='field', add_col_vals=field)
            print 
            # for column in rec.columns:
            #     if column not in data:
            #         data[column]=None
            #     for row in rec.index:
            #         if 'data_' in column:
            #             data.at[row, column] = np.array(rec.loc[row, column])
            #         else:
            #             data.at[row, column] = rec.loc[row, column]
            #         data.loc[row, 'field'] = field
                    
            #     data.at[row, column]

            
            # for variable_key, variable in rec.iteritems():

            #     # handle input times for each synapse. note each list entry can have a different number of input times
            #     if variable_key == 'input_times':
            #         if not isinstance(data[variable_key]['data'],list):
            #             data[variable_key]['data']=[]
            #         if 'syn_types' not in data[variable_key]:
            #             data[variable_key]['syn_types']=[]
            #         for loc_i, loc in enumerate(variable['data']):
            #             # print loc
            #             data[variable_key]['data'].append(np.array(loc))
            #         data[variable_key]['field'] += [field for seg_i in variable['locations']] 
            #         data[variable_key]['trial_id'] += [p['trial_id'] for seg_i in variable['locations']]
            #         data[variable_key]['locations'] += copy.copy(variable['locations'])
            #         data[variable_key]['syn_types']+=copy.copy(variable['syn_types'])
            #         data[variable_key]['t'] = np.asarray(data[variable_key]['t'])
            #     else:
            #         nseg = len(variable['locations'])
            #         print variable_key
            #         data[variable_key]['data'][field_i*nseg:(field_i+1)*nseg, :] = np.array(variable['data'])

            #         data[variable_key]['field'] += [field for seg_i in variable['locations']] 
            #         data[variable_key]['trial_id'] += [p['trial_id'] for seg_i in variable['locations']]
            #         data[variable_key]['locations'] += copy.copy(variable['locations'])
            #         data[variable_key]['t'] = np.asarray(data[variable_key]['t'])


        return data
    
    def _save_data_deprecated(self, data, file_name): # save data
        '''
        '''

        p = data.values()[0]['p']
        
        # check if folder exists with experiment name
        if os.path.isdir(p['data_folder']) is False:
            print 'making new directory to save data'
            os.mkdir(p['data_folder'])

        # save data as pickle file
        with open(p['data_folder']+file_name+'.pkl', 'wb') as output:
            
            print 'saving data'
            pickle.dump(data, output,protocol=pickle.HIGHEST_PROTOCOL)

    def _save_data(self, data, file_name, data_directory): # save data
        '''
        '''
        # check if folder exists with experiment name
        if os.path.isdir(data_directory) is False:
            print 'making new directory to save data'
            os.mkdir(data_directory)

        # save data as pickle file
        with open(data_directory+file_name+'.pkl', 'wb') as output:
            
            print 'saving data'
            pickle.dump(data, output,protocol=pickle.HIGHEST_PROTOCOL)

# procedures to be initialized if called as a script
if __name__ =="__main__":
    pass