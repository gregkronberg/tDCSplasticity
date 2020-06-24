"""
run control
"""
# imports
from mpi4py import MPI
# import multiprocessing 
import neuron
h = neuron.h
# from neuron import h
import numpy as np
import matplotlib.pyplot as plt
import param
import neurons
import run
import time
import uuid
import analysis
import sys
import copy
import pickle
import stims
import datetime
import itertools
import inspect
import mech_dicts
import specify_cells as spec
import pandas as pd
import functions as fncs
import functions
import os
import neuron_reduce 
import glob

class Exp(object):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        # get experiment name and directories for saving data and figures
        self.experiment_name = self.__class__.__name__
        self.data_directory = 'Data/'+self.experiment_name+'/'
        self.figure_directory =  'png figures/'+self.experiment_name+'/'
        # run experiment?
        if 'run' in kwargs and kwargs['run']:
            # check for run method and run
            if hasattr(self, 'run'):
                self.run(**kwargs)
            else:
                print 'method "run" does not exist'
    
    def run_and_save(self, cell=None, geo=None, run_method=None, trial=None, trial_id=None, df_funcs=None, df_func_kwargs=None, save=True, **kwargs):
        '''
        '''
        # check for parameter and rec objects
        #-------------------------------------
        assert hasattr(self, 'P'), 'a Parameter object "P" is required to run experiment'
        # update trial number
        #---------------------
        if trial is not None:
            self.P.p['trial'] = trial
        else:
            self.P.p['trial']=0
        # update trial_id
        #---------------------
        if trial_id is not None:
            self.P.p['trial_id'] = trial_id
        else:
            self.P.p['trial_id']=self._generate_trial_id()
        # setup recording
        #----------------------------
        # check that cell has geo attribute
        if cell is not None:
            assert hasattr(cell, 'geo'), 'cell object must have geo attribute to specify hoc sections'
            geo = cell.geo
        elif hasattr(self, 'geo'):
            geo = getattr(self, 'geo')
        # if geo attribute is found
        if geo is not None:
            # check for synapses
            if hasattr(self, 'syns'):
                syns=self.syns
            else:
                syns={}
            # check for netcons
            if hasattr(self, 'nc'):
                nc=self.nc
            else:
                nc={}
            # setup recording vectors
            self.rec = self.setup_recording_vectors_df(p=self.P.p, geo=geo, syns=syns, nc=nc)
        else:
            raise NameError('no hoc geometry specified to record from')
        # start timer
        #--------------------
        start = time.time() 
        # run simulation and record data
        #-------------------------------------
        default_run_method = 'run_sims'
        # check for passed run method
        if run_method is not None:
            self.data = getattr(self, run_method)(**kwargs)
        # otherwise use default run method
        else:
            self.data = getattr(self, default_run_method)(**kwargs)
        # self.data = self.run_sims_milstein(p=P.p, rec=self.rec)
        # end timer
        #------------------
        end = time.time() 
        # print trial and simulation time
        #---------------------------------
        print 'trial_id:'+ str(self.P.p['trial_id']) + ', duration:' + str(end -start) 
        # apply functions to data df before saving
        #---------------------------------------
        if df_funcs is not None:
            # iterate over df functions
            for func_i, df_func in enumerate(df_funcs):
                if df_func_kwargs is not None and len(df_func_kwargs)==len(df_funcs):
                    _kwargs=df_func_kwargs[func_i]
                else:
                    _kwargs ={}
                # apply df_func to data
                self.data = df_func(self.data, **_kwargs)
        # set file name to save data
        #----------------------------
        file_name = str(
            'data_'+
            self.experiment_name+
            '_trial_id_'+str(self.P.p['trial_id'])
            )
        # save data for each trial
        #-----------------------------
        if save:
            self._save_data(data=self.data, file_name=file_name, data_directory=self.data_directory)

    def run_sims(self, **kwargs):
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
        # preallocate data df
        #--------------------
        data=pd.DataFrame(dtype=object)

        # load hoc standard run environment
        #------------------------------------
        h.load_file("stdrun.hoc")

        # set hoc runtime variables
        #-------------
        h.dt = self.P.p['dt']
        h.tstop = self.P.p['tstop']
        h.celsius= self.P.p['celsius']
        h.v_init=self.P.p['v_init'] 
        # v_init is an object created in h when stdrun.hoc is called
        h.v_init=self.P.p['v_init'] 
        
        # run simulation
        h.run()
        print 'simulation finished'

        data = self.store_recording_vectors_to_data_df(data=data, rec=self.rec,**kwargs)

        return data

    def store_recording_vectors_to_data_df(self, data, rec, add_col=[], add_col_vals=[], data_key='data_'):
        ''' note that hoc recording vectors will be overwritten if multiple experiments are run.  hoc vectors need to be converted to numpy arrays and saved in order to not lose the data
        '''
        # temporary copy of rec dataframe
        #--------------------------------
        data_temp = copy.deepcopy(rec)
        # add specified columns to be saved in data but not recorded from hoc
        #--------------------------------------------------------------------
        # make sure column names and values are iterable
        if type(add_col)!=list and type(add_col)!=tuple:
            add_col=[add_col]
        if type(add_col_vals)!=list and type(add_col_vals)!=tuple:
            add_col_vals=[add_col_vals]
        # iterate column names
        for col_i, col in enumerate(add_col):
            # create column
            #--------------
            data_temp[col]=add_col_vals[col_i]

        # convert hoc recording vectors to np.array
        #------------------------------------------
        # iterate columns
        for column in data_temp.columns:

            # data recorded from hoc, should have column name as 'data_variable'
            if data_key in column:
                if column=='data_input_times':
                    print 'converting input_times', np.array(data_temp[column][1])
                # convert to numpy 
                data_temp[column] = data_temp[column].apply(np.array)

        # append to overall data df
        #------------------------------
        data = data.append(data_temp, ignore_index=True)
        return data

    def setup_recording_vectors_df(self, p, geo, syns=None, nc=None):
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
        # initialize rec df
        #-----------------------------
        rec=pd.DataFrame(dtype=object)
        # initialize columns to save information other than recorded data
        #----------------------------------------------------------------
        init_columns = ['location', 'tree_key', 'sec_num', 'seg_num', 'seg_loc','data_t', 'seg_dist', 'morpho', 'n_paths','path', 'paths', 'w', 'mechanism']
        rec = fncs._initialize_column(rec, init_columns)

        #iterate over locations to record from
        #----------------------------------------------
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
                if 'syn_idx' in path and location in path['syn_idx']:
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
                # print variable, variable_type, mechanism
                # print len(syns[tree_key][sec_num])
                # print syns[tree_key][sec_num][seg_num]
                if variable_type == 'syn' and syns is not None and tree_key in syns and mechanism in syns[tree_key][sec_num][seg_num][0].keys() and  variable in dir(syns[tree_key][sec_num][seg_num][0][mechanism]): 
                                
                    print variable_type, variable, var_rec
                    # point to variable to record
                    var_rec = getattr(syns[tree_key][sec_num][seg_num][0][mechanism], '_ref_'+variable)

                    # create recording vector
                    # create recording vector
                    rec.at[i, colname] = h.Vector()
                    rec.loc[i, colname].record(var_rec)

                # synapse input times
                #----------------------
                if variable== 'input_times' and variable_type == 'syn' and nc is not None:
                    # print 'getting input times'

                    # iterate over synapse pathways
                    for path_key, path in p['p_path'].iteritems():
                        # get segment locations from syn_idx
                        syn_idx_locations = zip(*zip(*path['syn_idx'])[:3])

                        # if the current location is in the current pathway
                        if location in syn_idx_locations:

                            # get index of the location in the pathway
                            seg_i = syn_idx_locations.index(location)

                            if type(nc[path_key][seg_i][mechanism])==list:
                                # seg_i = path['syn_idx'].index(location)
                                rec.at[i, colname]=[]
                                # iterate over bursts 
                                for burst_i, burst in enumerate(nc[path_key][seg_i][mechanism]):
                                    
                                    # print 'netcon:', nc[path_key][seg_i][mechanism][burst_i].postloc(), h.cas().hname()

                                    # add recording vector
                                    # rec.at[i, colname]=h.Vector()
                                    rec.at[i, colname].append(h.Vector())

                                    # set recording from netcon object
                                    nc[path_key][seg_i][mechanism][burst_i].record(rec.loc[i,colname][-1])
                            else:
                                rec.at[i, colname] = h.Vector()
                                # set recording from netcon object
                                nc[path_key][seg_i][mechanism].record(rec.loc[i,colname])



        return rec
    
    def setup_recording_vectors_df_network(self, P, network, syns=None, nc=None):
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
        # initialize rec df
        #-----------------------------
        rec=pd.DataFrame(dtype=object)
        p = P.p
        # initialize columns to save information other than recorded data
        #----------------------------------------------------------------
        init_columns = ['location', 'cell','tree_key', 'sec_num', 'seg_num', 'seg_loc','data_t', 'seg_dist', 'morpho', 'n_paths','path', 'paths', 'w', 'mechanism', 'syn_num', 'syn_type', 'syn_source']
        rec = fncs._initialize_column(rec, init_columns)

        #iterate over locations to record from
        #----------------------------------------------
        for i, location in enumerate(p['rec_idx']):
            # get location info
            #--------------------------------------------------------
            cell, tree_key, sec_num, seg_num = location['cell'], location['tree_key'], location['sec_num'], location['seg_num']
            seg_loc = float(seg_num+1)/(geo[tree_key][sec_num].nseg+1)
            rec.at[i, 'location']=location
            rec.at[i, 'cell']=cell
            rec.at[i, 'tree_key']=tree_key
            rec.at[i, 'sec_num']=sec_num
            rec.at[i, 'seg_num']=seg_num
            rec.at[i, 'seg_loc']=seg_loc

            #update parameters
            #-----------------------------------------------------------
            for param, value in p.iteritems():
                # initialize column
                rec = fncs._initialize_column(rec, param)
                # convert empty lists to nan
                if type(value)==list or type(value)==dict:
                    if len(value)==0:
                        value=np.nan
                    rec.at[i,param] = value
                else:
                    rec.at[i,param] = value

            # update morphology info
            #--------------------------------------------------------------
            # get distance from soma
            # FIXME adapt seg_dist for network
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
                if variable_type == 'syn' and syns is not None and tree_key in syns and mechanism in syns[tree_key][sec_num][seg_num].keys() and  variable in dir(syns[tree_key][sec_num][seg_num][mechanism]): 
                                
                    # point to variable to record
                    var_rec = getattr(syns[tree_key][sec_num][seg_num][mechanism], '_ref_'+variable)

                    # create recording vector
                    # create recording vector
                    rec.at[i, colname] = h.Vector()
                    rec.loc[i, colname].record(var_rec)

                # synapse input times
                #----------------------
                if variable== 'input_times' and variable_type == 'syn' and nc is not None:
                    print 'getting input times'

                    # iterate over synapse pathways
                    for path_key, path in p['p_path'].iteritems():
                        print path['syn_idx']
                        print zip(*path['syn_idx'])[:3]
                        syn_idx_locations = zip(*zip(*path['syn_idx'])[:3])
                        print syn_idx_locations
                        # if the current location is in the current pathway
                        if location in syn_idx_locations:
                            print 'location in syn_idx'

                            # get index of the location in the pathway
                            seg_i = syn_idx_locations.index(location)

                            if type(nc[path_key][seg_i][mechanism])==list:
                                # seg_i = path['syn_idx'].index(location)
                                rec.at[i, colname]=[]
                                # iterate over bursts 
                                for burst_i, burst in enumerate(nc[path_key][seg_i][mechanism]):
                                    
                                    # print 'netcon:', nc[path_key][seg_i][mechanism][burst_i].postloc(), h.cas().hname()

                                    print 'location,',location
                                    print path_key, seg_i, mechanism, burst_i
                                    # add recording vector
                                    # rec.at[i, colname]=h.Vector()
                                    rec.at[i, colname].append(h.Vector())

                                    # set recording from netcon object
                                    nc[path_key][seg_i][mechanism][burst_i].record(rec.loc[i,colname][-1])
                            else:
                                print 'recording from vecstim'
                                rec.at[i, colname] = h.Vector()
                                # set recording from netcon object
                                nc[path_key][seg_i][mechanism].record(rec.loc[i,colname])



        return rec
    
    def run_and_save_df_milstein(self, P, cell, trial, trial_id=None, **kwargs):
        '''
        '''
        # store trial number
        P.p['trial']=trial
        # data and figure folder
        P.p['data_directory'] = 'Data/'+P.p['experiment']+'/'
        P.p['fig_folder'] =  'png figures/'+P.p['experiment']+'/'
        
        if trial_id is None:
            trial_id = self._generate_trial_id()
        # # create unique identifier for each trial
        # uid = str(uuid.uuid1().int)[-5:]
        # now = datetime.datetime.now()
        # if trial_id is not None:
        #     trial_id = '-'.join(['{:04d}'.format(now.year), '{:02d}'.format(now.month), '{:02d}'.format(now.day), '{:02d}'.format(now.hour), '{:02d}'.format(now.minute), '{:02d}'.format(now.second), '{:02d}'.format(now.microsecond), uid])
        P.p['trial_id'] = trial_id#str(uuid.uuid4())
                    
        # start timer
        start = time.time() 
        
        # self.run_obj = run.Run()
        # self._standard_run_df(p=P.p, cell=cell, **kwargs)
        # load hoc standard run environment
        h.load_file("stdrun.hoc")
        # run simulation and record data
        self.data = self.run_sims_milstein(p=P.p, rec=self.rec)

        # end timer
        end = time.time() 

        # print trial and simulation time
        print 'trial'+ str(P.p['trial']) + ' duration:' + str(end -start) 
        
        # set file name to save data
        file_name = str(
            'data_'+
            P.p['experiment']+
            '_trial_'+str(P.p['trial'])+
            '_id_'+P.p['trial_id']
            )

        # save data for eahc trial
        self._save_data(data=self.data, file_name=file_name, data_directory=P.p['data_directory'])

    def update_parameters(self, p_update, paths_update, default_p=None, load_fd=True, data_folder='Data/',**kwargs):
        ''' setup parameter dictionary and load cell instance

        ===Args===
        -default_p  : string specifying default parameters.  calls the corresponding method in param.Default
        -p_update   : dictionary of global parameters that are specific to the current experiment. will be used to update the corresponding entries in the default dictionary specified by default_p
        -paths_update   : dictionary path-specific parameters.  Organized as paths_update{path name}{parameters}
        -cell       : string specifying the cell class to be instantiated from cell module as cell.CellClass()

        ===Out===
        -P  : default parameter class containing methods for specifying parameters
        -p  : updated parameter dictionary
                -morpho
                -seg_dist
        -paths  : updated path-specific parameter dictionary as paths{path name}{parameters}
        -cell   : cell class instance

        ===Updates===
        -P.p['fd_parameters', 'data_folder', 'fig_folder', 'seg_dist', 'morpho']

        ===Comments===
        '''
        # instantiate default parameter class
        P = param.Param()

        # load parameters from specified default parameters
        if default_p is not None:
            getattr(P, default_p)()

        # reference to default parameters
        p = P.p
        paths = P.paths

        # load facilitation depression parameters
        if load_fd:
            P._load_fd_parameters(p=p, filename='Data/fd_parameters.pkl')

        # update global parameter dictionary
        p.update(p_update)
        # update path dictionaries
        for key, val in paths_update.iteritems():
            # if path already exists
            if key in paths:
                # update
                paths[key].update(val)
            # if no paths exist
            elif not paths.values():
                # create path from paths_update
                paths[key]=val
            # if other paths exist
            else:
                # copy another path
                paths[key]=copy.copy(paths.values()[0])
                print paths[key]
                print val
                # update path
                paths[key].update(val)
        return P

    def _load_group_variable(self, variable='', directory=None, df=True, **kwargs):
        """ Load group data from folder
        
        ===Args===
        -directory : directory where group data is stored including /
        -filename : file name for group data file, including .pkl
                    -file_name cannot contain the string 'data', since this string is used t search for individual data files
        -df : boolean, if true group data is assumed to be a pandas dataframe, otherwise it is assumed to be a nested dictionary

        ===Out===
        -group_data  : if df is True, group_data will be a pandas dataframe.  if no group data is found in the directory an empty dataframe is returned.  if df is False, group_data is returned as a nested dictionary, and if no data is found an empty dictionary is returned

        ===Updates===
        -none

        ===Comments===
        """
        if directory is None:
            directory = self.data_directory
        # all files in directory
        # files = os.listdir(directory)
        files = glob.glob(directory+'*'+variable+'*')
        print 'files',files

        if len(files)>0:
            print len(files), 'variable files found'
            group_data=pd.DataFrame()
            for file in files:
                group_data_temp = pd.read_pickle(file)
                print 'group data loaded', file
                # print group_data_temp
                if group_data.empty:
                    group_data=group_data_temp
                else:
                    group_data=group_data.append(group_data_temp)

            # print 'keys',sorted(group_data.keys())

        # otherwise create data structure
        else:
            print 'no group data found'
            
            # if dataframe
            if df:
                # create empty dataframe
                group_data = pd.DataFrame()
            # else if dicitonary
            else:
                group_data= {}

        setattr(self, variable, group_data)

        return group_data 
    
    def generate_group_variable_vtrace(self, df_funcs=None, df_func_kwargs=None, **kwargs):
        '''
        '''
        directory = self.data_directory
        # load function for generating variables
        #---------------------------------------
        # varfuncs=fncs.VarFuncs()
        # df_funcs=fncs.DfFuncs()
        # vargen_funcs=fncs.VarGen()
        # check if dataframe exists in namespace
        #---------------------------------------
        if hasattr(self, 'vtrace_df'):
            var_exists = True
        else:
            var_exists = False

        if ('rerun' in kwargs and kwargs['rerun']) or not var_exists:

            variable='vtrace_df'
            filename=variable+'.pkl'
            var_functions=[fncs.VarFuncs()._pass]
            kwlist=[{}]
            rerun=[]#[varfuncs._get_vtrace]
            keep=[]
            file_limit=[]
            self.vtrace_df = fncs.VarGen()._vargen(variable=variable, directory=directory, functions=var_functions, kwlist=kwlist, rerun=rerun, keep=keep, file_limit=file_limit)

            # apply functions to data df before saving
            #---------------------------------------
            if df_funcs is not None:
                # iterate over df functions
                for func_i, df_func in enumerate(df_funcs):
                    if df_func_kwargs is not None and len(df_func_kwargs)==len(df_funcs):
                        _kwargs=df_func_kwargs[func_i]
                    else:
                        _kwargs ={}
                    # apply df_func to data
                    self.vtrace_df = df_func(self.vtrace_df, **_kwargs)
            # save vtrace_df 
            #---------------------------------
            self.vtrace_df = functions._save_group_data(df=self.vtrace_df, directory=directory, variable=variable, extension='.pkl', check_against=None, **kwargs)
    
    def generate_group_variable_w_clopath(self, clopath_param=None, rerun=False, save=True, **kwargs):
        '''
        '''
        directory = self.data_directory
        variable='w_clopath_df'
        filename=variable+'.pkl'
        # load function for generating variables
        #--------------------------------------------------------------------
        varfuncs=fncs.VarFuncs()
        df_funcs=fncs.DfFuncs()
        vargen_funcs=fncs.VarGen()
        # parameters for clopath rule
        #--------------------------------------------------------------------
        # passed as argument
        if clopath_param is not None:
            _clopath_param=clopath_param
        # check for clopath parameters in P for current experiment
        elif hasattr(self, 'P'):
            _clopath_param ={}
            for _key, _val in self.P.p.iteritems():
                if 'clopath' in _key:
                    _clopath_param[_key]=_val
        # otherwise set as default
        else:
            _clopath_param = param.ParamClopath().kronberg_2020()
        # remove clopath tag from parameter names
        clopath_param={}
        for _key, _val in _clopath_param.iteritems():
            if 'clopath' in _key:
                clopath_key = _key.split('clopath_')[-1]
                clopath_param[clopath_key] = _val

        # requires voltage data
        #-----------------------
        assert hasattr(self, 'vtrace_df'), 'voltage data required to create w_clopath_df'
        print 'generating w_clopath_df'
        # rerun
        #-------------------------------------------------------------------
        if rerun:
            start=time.time()
            self.w_clopath_df = df_funcs._get_w_clopath(vtrace_df=self.vtrace_df, clopath_param=clopath_param, w_df=None, **kwargs)
            end=time.time()
            print 'timer:', (end-start)
        # if not rerun, check for previous data and update if it exists
        #--------------------------------------------------------------
        else:
            start=time.time()
            # if w_df does not exist yet, load from directory
            if not hasattr(self, 'w_clopath_df'):
                # load variable
                self.w_clopath_df = functions._load_group_data(directory=directory, filename=variable)
            # if w_df is not empty, pass to df_func to update
            if not self.w_clopath_df.empty:
                self.w_clopath_df = df_funcs._get_w_clopath(vtrace_df=self.vtrace_df, clopath_param=clopath_param, w_df=self.w_clopath_df, **kwargs)
            # if w_df is empty, run df_func on entire v_df
            else:
                self.w_clopath_df = df_funcs._get_w_clopath(vtrace_df=self.vtrace_df, clopath_param=clopath_param, w_df=None, **kwargs)
            end=time.time()
            print 'timer:', (end-start)
        if save:
            # save vtrace_df 
            #---------------------------------
            self.w_clopath_df = functions._save_group_data(df=self.w_clopath_df, directory=directory, variable=variable, extension='.pkl', check_against=None, **kwargs)

    def _save_data(self, data, file_name=None, data_directory=None, **kwargs): # save data
        '''
        '''
        if file_name is None:
            # set file name to save data
            #----------------------------
            file_name = str(
                'data_'+
                self.experiment_name+
                '_trial_id_'+str(self.P.p['trial_id'])
                )
        if data_directory is None:
            data_directory=self.data_directory
        # check if folder exists with experiment name
        if os.path.isdir(data_directory) is False:
            print 'making new directory to save data'
            os.mkdir(data_directory)

        # save data as pickle file
        with open(data_directory+file_name+'.pkl', 'wb') as output:
            
            print 'saving data'
            pickle.dump(data, output,protocol=pickle.HIGHEST_PROTOCOL)

    def _setup_piecewise_field(self, field, P, default=True, testing=False, **kwargs):
        '''
        '''
        def param_default():
            slope_prox=field
            threshold_basal=-10 # border between basal and soma
            threshold_prox=10 # border between soma and prox apical
            threshold_dist=1200 # border between proximal and distal apical
            # slopes
            slope_dist=1.6*slope_prox # distal apical slope
            slope_basal=2.*slope_prox # basal slope
            # offsets
            offset_prox = 200.*slope_prox # proximal apical offset
            offset_dist = -100.*slope_dist + (slope_prox-slope_dist)*threshold_dist+offset_prox # distal apical offset
            offset_basal = 0.*1000.*slope_prox # basal offset
            # store in parameter dictionary in order to get saved
            P.p['field']=slope_prox
            P.p['slope_prox']=slope_prox 
            P.p['slope_dist']=slope_dist
            P.p['slope_basal']=slope_basal
            P.p['offset_prox']=offset_prox
            P.p['offset_dist']=offset_dist
            P.p['offset_basal']=offset_basal
            P.p['threshold_prox']=threshold_prox
            P.p['threshold_dist']=threshold_dist
            P.p['threshold_basal']=threshold_basal
            return P
        
        def param_1():
            slope_prox=field
            threshold_basal=-10 # border between basal and soma
            threshold_prox=10 # border between soma and prox apical
            threshold_dist=1200 # border between proximal and distal apical
            # slopes
            slope_dist=1.6*slope_prox # distal apical slope
            slope_basal=2.4*slope_prox # basal slope
            # offsets
            offset_prox = 200.*slope_prox # proximal apical offset
            offset_dist = -100.*slope_dist + (slope_prox-slope_dist)*threshold_dist+offset_prox # distal apical offset
            offset_basal = 0.*1000.*slope_prox # basal offset
            # store in parameter dictionary in order to get saved
            P.p['field']=slope_prox
            P.p['slope_prox']=slope_prox 
            P.p['slope_dist']=slope_dist
            P.p['slope_basal']=slope_basal
            P.p['offset_prox']=offset_prox
            P.p['offset_dist']=offset_dist
            P.p['offset_basal']=offset_basal
            P.p['threshold_prox']=threshold_prox
            P.p['threshold_dist']=threshold_dist
            P.p['threshold_basal']=threshold_basal
            return P
        
        if default:
            return param_default()
        else:
            return param_1()

    def _generate_trial_id(self, ):
        '''
        '''
        # create unique identifier for each trial
        uid = str(uuid.uuid1().int)[-5:]
        now = datetime.datetime.now()
        trial_id = '-'.join(['{:04d}'.format(now.year), '{:02d}'.format(now.month), '{:02d}'.format(now.day), '{:02d}'.format(now.hour), '{:02d}'.format(now.minute), '{:02d}'.format(now.second), '{:02d}'.format(now.microsecond), uid])
        return trial_id
    #########################################################################
    #########################################################################
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

    def run_sims_milstein(self, p, rec):
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
        
        
        # data['field']=None
        data=pd.DataFrame(dtype=object)

        # run time
        h.dt = p['dt']
        h.tstop = p['tstop']
        h.celsius= p['celsius']

        # load standard run environment
        
        # print 'running simulation, electric field:', field
        
        # initialize voltage
        h.v_init=p['v_init'] # v_init is an object created in h when stdrun.hoc is called
        
        # run simulation
        h.run()
        print 'simulation finished'

        data = self._store_recording_vectors_to_data_df(data=data, rec=rec,)

        return data

    def _standard_parameter_setup(self, default_p, p_update, paths_update, cell_class, load_fd=True, data_folder='Data/',**kwargs):
        ''' setup parameter dictionary and load cell instance

        ===Args===
        -default_p  : string specifying default parameters.  calls the corresponding method in param.Default
        -p_update   : dictionary of global parameters that are specific to the current experiment. will be used to update the corresponding entries in the default dictionary specified by default_p
        -paths_update   : dictionary path-specific parameters.  Organized as paths_update{path name}{parameters}
        -cell       : string specifying the cell class to be instantiated from cell module as cell.CellClass()

        ===Out===
        -P  : default parameter class containing methods for specifying parameters
        -p  : updated parameter dictionary
                -morpho
                -seg_dist
        -paths  : updated path-specific parameter dictionary as paths{path name}{parameters}
        -cell   : cell class instance

        ===Updates===
        -P.p['fd_parameters', 'data_folder', 'fig_folder', 'seg_dist', 'morpho']

        ===Comments===
        '''
        # instantiate default parameter class
        P = param.Param()

        # load parameters from specified default parameters
        getattr(P, default_p)()

        # reference to default parameters
        p = P.p
        paths = P.paths

        # load facilitation depression parameters
        if load_fd:
            P._load_fd_parameters(p=p, filename='Data/fd_parameters.pkl')

        # update global parameter dictionary
        p.update(p_update)
        # update path dictionaries
        for key, val in paths_update.iteritems():
            # if path already exists
            if key in paths:
                # update
                paths[key].update(val)
            # if no paths exist
            elif not paths.values():
                # create path from paths_update
                paths[key]=val
            # if other paths exist
            else:
                # copy another path
                paths[key]=copy.copy(paths.values()[0])
                # update path
                paths[key].update(val)

        # data and figure folder
        p['data_folder'] = data_folder+p['experiment']+'/'
        p['fig_folder'] =  'png figures/'+p['experiment']+'/'

        # load cell and store in parameter dictionary
        cell = getattr(neurons, cell_class)(p)
        cell.geometry(p)
        # insert mechanisms
        cell.mechanisms(p)
        
        # measure distance of each segment from the soma and store in parameter dictionary
        p['seg_dist'] = P._seg_distance(cell)

        # create morphology for shape plots
        p['morpho'] = P._create_morpho(cell.geo)

        return P, cell

    def _update_synapse_parameters(self, P, cell, method='_choose_seg_rand',**kwargs):
        ''' update parameter dictionaries for each pathway before running simulation

        ===Args===
        -p_class    : instance of default parameter class, containing parameter dictionaries and methods for updating parameters for simulations
        -p          : full parameter dictionary
        -paths      : parameter dictionary for separate synaptic pathways, organized as paths{path name}{parameters}
        cell1       : instance of cell class containing geometry and synapse structures (contain all hoc section and synapse objects)

        ===Out===

        ===Updates===
        -p          : updated path dictionaries are added to the main parameter dictionary p (organized as p['p_path']{path name}{parameters})
        -paths      : active synapses and their activity patterns are set and updated for each path

        ===Comments===
        -p  : p should have an entry called 'path_combo', which is a list of paths to be activated during the current simulation
        '''
        p_class=P
        p = P.p
        paths=P.paths



        # add path parameter dictionaries to global p dictionary
        p['p_path']={}

        # list of segments to record from
        p['rec_idx']=[]

        # update pathways and add to global p structure
        for path_key, path in paths.iteritems():
            path['syn_idx']=[]
            path['sequence_delays']=[]
            # if path is included in current combination
            if path_key in p['active_paths']:

                if method=='_choose_seg_rand':

                    # choose active segments for this pathway
                    path['syn_idx'], path['syn_counts'] = p_class._choose_seg_rand(p=p, p_path=path, syns=cell.syns, replace=path['replace_syn'])

                # set weights for each segment in this pathway
                path['w_idx'] = p_class._set_weights_normal(p_path=path)

                # set delays for branch sequences in this pathway 
                if 'delays' in path:
                    path['sequence_delays'] = p_class._set_sequence_delays(syn_idx=path['syn_idx'], delay=path['delay'])
                else:
                    path['sequence_delays'] = p_class._set_sequence_delays(syn_idx=path['syn_idx'], delay=0)

            # record all activated synapses
            p['rec_idx'] += path['syn_idx']
            # update p_path dictionary in global p dictionary
            p['p_path'][path_key]=copy.copy(path)

        # add soma and axon to rec_idx
        if 'soma' in cell.geo:
            p['rec_idx']+=[('soma',0,0)]
        if 'axon' in cell.geo:
            p['rec_idx']+=[('axon',0,0)]

        # p['rec_idx']+=[('soma',0,0),('axon',0,0)]

        # update tstop based on synaptic inputs in each path
        tstops=[]
        warmups=[]
        for path_key, path in paths.iteritems():
            if path_key in p['active_paths']:
                # print path['sequence delays']
                if 'sequence_delays' in path and len(path['sequence_delays'])>0:
                    # print path['sequence_delays']
                    max_delay  = max(path['sequence_delays'])
                else:
                    max_delay=0
                warmups.append(path['warmup'])
                tstops.append(path['warmup'] + max_delay + 1000*(path['bursts']-1)/path['burst_freq'] + 1000*(path['pulses']+1)/path['pulse_freq'] )

            else:
                tstops = [70]
                warmups= [10]

        p['tstop'] = max(tstops)
        p['warmup'] = min(warmups)
        p['field_off'] = p['tstop']
        # FIXME
        p['field_on'] = p['warmup']-10

        return P

    # update tstop based on synaptic inputs in each path
    def update_time_parameters(self, P):
        '''
        '''
        paths=P.paths
        p=P.p
        tstops=[]
        warmups=[]
        for path_key, path in paths.iteritems():
            if path_key in p['active_paths']:
                # print path['sequence delays']
                if 'sequence_delays' in path and len(path['sequence_delays'])>0:
                    # print path['sequence_delays']
                    max_delay  = max(path['sequence_delays'])
                else:
                    max_delay=0
                warmups.append(path['warmup'])
                tstops.append(path['warmup'] + max_delay + 1000*(path['bursts']-1)/path['burst_freq'] + 1000*(path['pulses']+1)/path['pulse_freq'] )

            else:
                tstops = [70]
                warmups= [10]

        p['tstop'] = max(tstops)
        p['warmup'] = min(warmups)
        p['field_off'] = p['tstop']
        # FIXME
        p['field_on'] = p['warmup']-10
        return P

    def _update_synapse_parameters_sequence(self, P, cell, method='_choose_seg_from_branch', reverse=False, **kwargs):
        ''' update parameter dictionaries for each pathway before running simulation specifically for branch sequence 

        ===Args===
        -P : parameter class object
        -cell : cell object
        -method : method for choosing active synapses
        -reverse : reverses sequence order
        ===Return===
        -P : updated parameter object
        ===Updates===
        -p          : updated path dictionaries are added to the main parameter dictionary p (organized as p['p_path']{path name}{parameters})
        -paths      : active synapses and their activity patterns are set and updated for each path
        ===Comments===
        '''
        p_class=P
        p = P.p
        paths=P.paths



        # add path parameter dictionaries to global p dictionary
        p['p_path']={}

        # list of segments to record from
        p['rec_idx']=[]

        # update pathways and add to global p structure
        for path_key, path in paths.iteritems():
            path['syn_idx']=[]
            path['sequence_delays']=[]
            # if path is included in current combination
            if path_key in p['active_paths']:
                # get sec_idx
                if 'sec_idx' in kwargs:
                    sec_idx=kwargs['sec_idx']
                else:
                    terminal_branches = stims._get_terminal_branches(cell.geo)
                    sec_idx={'apical_tuft':[69]}
                # print 'geo:',cell.geo[sec_idx.keys()[0]][sec_idx.values()[0][0]].nseg
                # set branch nseg
                geo = stims._set_branch_nseg(cell.geo, sec_idx, seg_L=p['branch_seg_L'])
                # update synapses after nseg
                # cell.syns = stims._update_synapses_after_nseg(p, cell.geo, cell.syns, sec_idx)
                # print 'syns:',len(cell.syns[sec_idx.keys()[0]][sec_idx.values()[0][0]])
                # print 'geo:',cell.geo[sec_idx.keys()[0]][sec_idx.values()[0][0]].nseg

                cell.mechanisms(p)
                # print 'syns:',len(cell.syns[sec_idx.keys()[0]][sec_idx.values()[0][0]])
                # print 'syns all:', cell.syns
                # print 'geo:',cell.geo[sec_idx.keys()[0]][sec_idx.values()[0][0]].nseg
                # choose segments to activate
                path['syn_idx'], path['syn_counts'] = stims._choose_seg_from_branch(cell.geo, sec_idx)

                # FIXME adjust so that you can choose how far along the branch you want to stimulate
                if 'syn_limit' in path:
                    if len(path['syn_idx'])>path['syn_limit']:
                        # print path['syn_idx']
                        path['syn_idx'] = path['syn_idx'][-path['syn_limit']:]
                        path['syn_counts'] = path['syn_counts'][-path['syn_limit']:]
                # path['syn_idx']=path['syn_idx'][:3]
                # path['syn_counts']=path['syn_counts'][:3]

                if reverse:
                    path['syn_idx'].reverse()
                    path['syn_counts'].reverse()

                # set delays
                path['sequence_delays'] = p_class._set_sequence_delays(syn_idx=path['syn_idx'], delay=path['delay'])
                # set weights for each segment in this pathway
                path['w_idx'] = p_class._set_weights_normal(p_path=path)

            # record all activated synapses
            p['rec_idx'] += path['syn_idx']
            # update p_path dictionary in global p dictionary
            p['p_path'][path_key]=copy.copy(path)

        # add soma and axon to rec_idx
        if 'soma' in cell.geo:
            p['rec_idx']+=[('soma',0,0)]
        if 'axon' in cell.geo:
            p['rec_idx']+=[('axon',0,0)]

        # p['rec_idx']+=[('soma',0,0),('axon',0,0)]

        # update tstop based on synaptic inputs in each path
        tstops=[]
        warmups=[]
        for path_key, path in paths.iteritems():
            if path_key in p['active_paths']:
                # print path['sequence delays']
                if 'sequence_delays' in path and len(path['sequence_delays'])>0:
                    # print path['sequence_delays']
                    max_delay  = max(path['sequence_delays'])
                else:
                    max_delay=0
                warmups.append(path['warmup'])
                tstops.append(path['warmup'] + max_delay + 1000*(path['bursts']-1)/path['burst_freq'] + 1000*(path['pulses']+1)/path['pulse_freq'] )

            else:
                tstops = [70]
                warmups= [10]

        p['tstop'] = max(tstops)
        p['warmup'] = min(warmups)
        p['field_off'] = p['tstop']
        # FIXME
        p['field_on'] = p['warmup']-10
        # cell.mechanisms(p)

        return P

    def _standard_run_and_save(self, P, cell, trial, trial_id=None, **kwargs):
        '''
        '''
        # store trial number
        P.p['trial']=trial
        # data and figure folder
        P.p['data_directory'] = data_folder+p['experiment']+'/'
        P.p['fig_folder'] =  'png figures/'+p['experiment']+'/'
        
        if trial_id is None:
            trial_id = self._generate_trial_id()
        # # create unique identifier for each trial
        # uid = str(uuid.uuid1().int)[-5:]
        # now = datetime.datetime.now()
        # if trial_id is not None:
        #     trial_id = '-'.join(['{:04d}'.format(now.year), '{:02d}'.format(now.month), '{:02d}'.format(now.day), '{:02d}'.format(now.hour), '{:02d}'.format(now.minute), '{:02d}'.format(now.second), '{:02d}'.format(now.microsecond), uid])
        P.p['trial_id'] = trial_id#str(uuid.uuid4())
                    
        # start timer
        start = time.time() 
        
        self.run_obj = run.Run()
        self.run_obj._standard_run(p=P.p, cell=cell, **kwargs)

        # end timer
        end = time.time() 

        # print trial and simulation time
        print 'trial'+ str(P.p['trial']) + ' duration:' + str(end -start) 
        
        # set file name to save data
        file_name = str(
            'data_'+
            P.p['experiment']+
            '_trial_'+str(P.p['trial'])+
            '_id_'+P.p['trial_id']
            )

        # save data for eahc trial
        self.run_obj._save_data(data=self.run_obj.data, file_name=file_name, data_directory=P.p['data_directory'])

        return self.run_obj

    def _standard_run_and_save_df(self, P, cell, trial, trial_id=None, **kwargs):
        '''
        '''
        # store trial number
        P.p['trial']=trial
        # data and figure folder
        P.p['data_directory'] = 'Data/'+P.p['experiment']+'/'
        P.p['fig_folder'] =  'png figures/'+P.p['experiment']+'/'
        
        if trial_id is None:
            trial_id = self._generate_trial_id()
        # # create unique identifier for each trial
        # uid = str(uuid.uuid1().int)[-5:]
        # now = datetime.datetime.now()
        # if trial_id is not None:
        #     trial_id = '-'.join(['{:04d}'.format(now.year), '{:02d}'.format(now.month), '{:02d}'.format(now.day), '{:02d}'.format(now.hour), '{:02d}'.format(now.minute), '{:02d}'.format(now.second), '{:02d}'.format(now.microsecond), uid])
        P.p['trial_id'] = trial_id#str(uuid.uuid4())
                    
        # start timer
        start = time.time() 
        
        self.run_obj = run.Run()
        self.run_obj._standard_run_df(p=P.p, cell=cell, **kwargs)

        # end timer
        end = time.time() 

        # print trial and simulation time
        print 'trial'+ str(P.p['trial']) + ' duration:' + str(end -start) 
        
        # set file name to save data
        file_name = str(
            'data_'+
            P.p['experiment']+
            '_trial_'+str(P.p['trial'])+
            '_id_'+P.p['trial_id']
            )

        # save data for eahc trial
        self.run_obj._save_data(data=self.run_obj.data, file_name=file_name, data_directory=P.p['data_directory'])

        return self.run_obj

    def run_and_save_df_milstein(self, P, cell, trial, trial_id=None, **kwargs):
        '''
        '''
        # store trial number
        P.p['trial']=trial
        # data and figure folder
        P.p['data_directory'] = 'Data/'+P.p['experiment']+'/'
        P.p['fig_folder'] =  'png figures/'+P.p['experiment']+'/'
        
        if trial_id is None:
            trial_id = self._generate_trial_id()
        # # create unique identifier for each trial
        # uid = str(uuid.uuid1().int)[-5:]
        # now = datetime.datetime.now()
        # if trial_id is not None:
        #     trial_id = '-'.join(['{:04d}'.format(now.year), '{:02d}'.format(now.month), '{:02d}'.format(now.day), '{:02d}'.format(now.hour), '{:02d}'.format(now.minute), '{:02d}'.format(now.second), '{:02d}'.format(now.microsecond), uid])
        P.p['trial_id'] = trial_id#str(uuid.uuid4())
                    
        # start timer
        start = time.time() 
        
        # self.run_obj = run.Run()
        # self._standard_run_df(p=P.p, cell=cell, **kwargs)
        # load hoc standard run environment
        h.load_file("stdrun.hoc")
        # run simulation and record data
        self.data = self.run_sims_milstein(p=P.p, rec=self.rec)

        # end timer
        end = time.time() 

        # print trial and simulation time
        print 'trial'+ str(P.p['trial']) + ' duration:' + str(end -start) 
        
        # set file name to save data
        file_name = str(
            'data_'+
            P.p['experiment']+
            '_trial_'+str(P.p['trial'])+
            '_id_'+P.p['trial_id']
            )

        # save data for eahc trial
        self._save_data(data=self.data, file_name=file_name, data_directory=P.p['data_directory'])

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

    def _setup_recording_vectors_milstein(self, P, cell, node_names=True):
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
        p = P.p
        print 'setting up recording vectors'
        rec=pd.DataFrame(dtype=object)
        init_columns = ['location', 'tree_key', 'sec_num', 'seg_num', 'seg_loc','data_t', 'seg_dist', 'morpho', 'n_paths','path', 'paths', 'w', 'mechanism']
        rec = fncs._initialize_column(rec, init_columns)


        var_rec_list = []
        #iterate over locations to record from
        for i, location in enumerate(p['rec_idx']):
            print i, location
            # if node name is specified as a string (assumes a dictionary cell._node_names exists) and locations are specified as [(node name, loc)]
            if node_names:
                name = location[0]
                seg_loc = location[1]
                print name, cell._node_names[name].name
                node = cell._node_names[name]
                tree_key = node.type
                sec_num = [_i for _i, _node in enumerate(cell._node_dict[tree_key]) if _node.name==name][0]
                seg=node.sec(seg_loc)
                # print seg.x
                # print [seg.x for seg in node.sec ]
                # seg_num = [i for i, seg in enumerate(node.sec) if seg==node.sec(seg_loc)]
                # print seg_num
                seg_num=0
            else:
                # get location info
                #--------------------------------------------------------
                tree_key, sec_num, seg_num, seg_loc = location
                node = cell._node_dict[tree_key][sec_num]

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

            # synaptic input path parameters
            #--------------------------------------------------------------
            # keep track of paths that the current segment belongs to
            paths=[]
            # total number of paths
            n_paths=len(P.paths.keys())
            rec.at[i, 'n_paths']=n_paths
            # iterate over input pathways
            for path_key, path in P.paths.iteritems():
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
                print variable, variable_type, mechanism, node.name
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
                if variable_type == 'range' and  mechanism in dir(node.sec(seg_loc)):
                    # point to variable for recording
                    # var_rec = getattr(node.sec(seg_loc), '_ref_'+variable)
                    var_rec_list.append(getattr(node.sec(seg_loc), '_ref_'+variable))
                    # create recording vector
                    rec.at[i, colname] = h.Vector()
                    rec.loc[i, colname].record(var_rec_list[-1])

                # synaptic variables
                #---------------------
                # list synapses in the current node
                current_syns = [(_i, _syn) for _i, _syn in enumerate(node.synapses) if node.sec(_syn.loc)==node.sec(seg_loc)]
                if variable_type == 'syn':
                    # iterate over synapses in the current node
                    for syn_i, syn in current_syns:
                        # get synapse types (e.g. AMPA, NMDA etc)
                        syn_types = syn._syn.keys()
                        # if the requested mechanism to record is in the current synapse and the requested variable exists in the synapse
                        if mechanism in syn_types and variable in dir(syn.target(mechanism)):
                            # point to variable to record
                            var_rec = getattr(node.synapses[syn_i].target(mechanism),'_ref_'+variable)
                            var_rec_list.append(getattr(node.synapses[syn_i].target(mechanism),'_ref_'+variable))
                            # create recording vector
                            rec.at[i, colname] = h.Vector()
                            rec.loc[i, colname].record(var_rec_list[-1])
        return rec
    
    def _store_recording_vectors_to_data_df(self, data, rec, add_col=[], add_col_vals=[], data_key='data_'):
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

    def run_sims_milstein(self, p, rec):
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
        
        
        # data['field']=None
        data=pd.DataFrame(dtype=object)

        # run time
        h.dt = p['dt']
        h.tstop = p['tstop']
        h.celsius= p['celsius']

        # load standard run environment
        
        # print 'running simulation, electric field:', field
        
        # initialize voltage
        h.v_init=p['v_init'] # v_init is an object created in h when stdrun.hoc is called
        
        # run simulation
        h.run()
        print 'simulation finished'

        data = self._store_recording_vectors_to_data_df(data=data, rec=rec,)

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

class exp_branch_sequence_migliore(Exp):
    '''
    '''

    def __init__(self, **kwargs):
        super(exp_branch_sequence_migliore, self).__init__(**kwargs)

    def run(self, **kwargs):

        ''' simulate sequences of inputs on various dendritic branches as in Branco and Hausser 2010

        randomly choose dendritic branch, reduce the size of each segment to be ~10 um, activate segments in a sequence (either towards or away from the terminal), vary the following parameters
        ==Parameters==
        -field : 0, 1, 5, 10, 20 (positive and negative)
        -branches : first five branches in eahc tree
        -sequence delays : 1, 2, 5, 10, ms between activation of neighboring segments (segments are ~10 um apart)
        -synaptic weights : .001, .0015, .002
        '''
        # updates to global parameter dictionary 
        p_update = {
        'experiment' : inspect.stack()[0][3], 
        'trials' : 1,
        'field':[0],
        'rec_variables':[('v','range','v'), ('i','syn','nmda')],
        'active_paths':['1',],
        # set active conductances to zero
        #---------------------------------
        # 'gna':0.,
        # 'ghd':0.,
        # 'gkdr':0., 
        # 'gcalbar': 0.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        # 'KMULT':0.,
        # 'KMULTP':0.,
        # 'Cm':1.,
        'alpha_vspom_nmda':-.062, 
        'v0_block_nmda':0,#0#-5#0,
        'branch_seg_L':5, 
        'dgna':-.0001,
        'AXONM':10
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apical_tuft','apical_trunk'],
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': [],
        'pulses': 1.,
        'pulse_freq': 100.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1.5*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        'delay':.1, 
        'sequence_direction':'in',
        'nmda_ampa_ratio':3.,
        'syn_limit':12
        },}

        # setup cell and updated parameter structures
        self.P, self.cell = self._standard_parameter_setup(
            default_p='migliore_2005',
            cell_class='CellMigliore2005',
            p_update=p_update,
            paths_update=paths_update,
            load_fd=True)

        # get sections with terminal branches
        terminal_branches = stims._get_terminal_branches(self.cell.geo)
        # maximum branches to simulate per tree
        max_branches = 1
        delays = [ 2, ]
        directions = ['in','out']
        weights = .001
        # iterate over branches
        for tree_key, branch_sec in terminal_branches.iteritems():
            if tree_key!='soma' and tree_key!='axon':
                for sec_i, sec_num in enumerate(branch_sec):
                    if sec_i<max_branches:
                        
                        # setup cell and updated parameter structures
                        self.P, self.cell = self._standard_parameter_setup(
                            default_p='migliore_2005',
                            cell_class='CellMigliore2005',
                            p_update=p_update,
                            paths_update=paths_update,
                            load_fd=True)
                        # section to simulate
                        sec_idx = {tree_key:[sec_num]}
                        self.P.paths['1']['sec_idx']=sec_idx
                        # iterate over delays
                        for delay in delays:
                            self.P.paths['1']['delay']=delay
                            # iterate over directions
                            for direction in directions:
                                self.P.paths['1']['sequence_direction']=direction
                                if direction=='in':
                                    reverse=True
                                elif direction=='out':
                                    reverse=False
                                # iterate over trials
                                for weight in weights:
                                    self.P.paths['1']['w_mean']=weight
                                    for trial_i, trial in enumerate(range(self.P.p['trials'])):
                                        self.P = self._update_synapse_parameters_sequence(P=self.P, cell=self.cell, method='_choose_seg_rand', reverse=reverse, sec_idx=sec_idx)
                                        # measure distance of each segment from the soma and store in parameter dictionary
                                        self.P.p['seg_dist'] = self.P._seg_distance(self.cell)

                                        # create morphology for shape plots
                                        self.P.p['morpho'] = self.P._create_morpho(self.cell.geo)

                                        self.run_obj = self._standard_run_and_save_df(P=self.P, cell=self.cell, trial=trial)

class exp_test_place_cells(Exp):
    '''
    '''
    def __init__(self, **kwargs):
        super(exp_test_place_cells, self).__init__(**kwargs)

    def run(self, **kwargs):
        '''
        '''
        # updates to global parameter dictionary 
        p_update = {
        'rec_variables':[('v','range','v')],#,('input_times','syn','ampa'),],
        'active_paths':['1','2'],
        'tstop':1000,
        'gna_inact':1,
        'dgna' : -.000025,          # change in sodium conductance with distance (ohm/cm2/um) from Kim 2015
        'gcalbar': 1.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        'ka_grad': 0.6,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        # 'ghd_grad': 1.5,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apic'],
        'nmda_ampa_ratio':1.,
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': [],
        'pulses': 4.,
        'pulse_freq': 100.,
        'bursts': 15.,
        'burst_freq': 5.,
        'warmup': 20,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        },
        }
        # create param object
        self.P = param.ParamMigliore2005()
        # update parameters
        self.P.p['p_path']={}
        self.P.p.update(p_update)
        self.P.p['p_path'].update(paths_update)

        # create reduced cell
        self.cell = neurons.CellMigliore2005Reduced(P=self.P, cell_number=0)
        self.syns = self.cell.syns

        n_inputs = 15
        self.place_cell = stims.PlaceCell(n_neurons=n_inputs)
        
        
        # record from activated segment
        #------------------------------
        self.P.p['rec_idx']=[('soma', 0, 0),]
        self.P.p['p_path']['1']['syn_idx']=[]


        apical_group = range(0,11)#[6,7,8,9,10]
        field_mags = [-5,0,5]
        field_mags = [0]

        # number of trials 
        #----------------------
        trials=1
        # w_mean
        #------------------------
        # w_mean_default = 5.5E-3
        w_mean_default=2.E-3

        # synaptic weights
        #------------------------------
        self.P.p['p_path']['1']['w_mean'] = w_mean_default
        self.P.p['p_path']['1']['w_std'] = 0.1*self.P.p['p_path']['1']['w_mean']#1.*.001
        self.P.p['p_path']['1']['w_rand'] = False
        # create syn_idx and rec_idx
        #--------------------------------------
        self.P.p['syn_num'] = len(apical_group)
        self.P.p['rec_idx']=[('soma', 0, 0),]
        self.P.p['p_path']['1']['syn_idx']=[]
        for seg_i in apical_group:
            rec_loc = ('apic', 0, seg_i)
            syn_loc = ('apic', 0, seg_i, 0)

            self.P.p['p_path']['1']['syn_idx'].append(syn_loc)
            self.P.p['rec_idx'].append(rec_loc)

        self.syn_map=[]
        for _syn in self.P.p['p_path']['1']['syn_idx']:
            if _syn[2]>5:
                source=0
            else:
                source=1
            self.syn_map.append({
                'cell':self.cell.name,
                'tree_key':_syn[0],
                'sec_num':_syn[1],
                'seg_num':_syn[2],
                'seg_x':[seg for seg in self.cell.geo[_syn[0]][_syn[1]]][_syn[2]].x,
                'syn_num':_syn[3],
                'syn_types':['ampa','nmda','clopath'],
                'weights':[w_mean_default, w_mean_default, w_mean_default],
                'source_i':range(n_inputs),#[_syn[2]], 
                # 'source_i':[source] 
                })
        self.nc = {}
        self.nc['1'] = self.place_cell.connect_synapses(syn_map=self.syn_map, syns=self.syns, sources=self.place_cell.stim)

        # iterate trials
        for trial in range(trials):

            # trial_id
            #-----------------------------------
            trial_id = self._generate_trial_id()
            # iterate over field magnitudes
            #--------------------------------
            for field in field_mags:
                # piecewise liner uniform field over basal, proximal apical, and distal apical dendrites empirically matched to poalrization in full neuron
                self.P = self._setup_piecewise_field(field=field, P=self.P)
                # insert field
                #--------------
                self.dcs = stims.DCS(cell=self.cell, method='piecewise', field_on=5, field_off=self.P.p['tstop'],
                    slope_prox=self.P.p['slope_prox'], 
                    slope_dist=self.P.p['slope_dist'], 
                    slope_basal=self.P.p['slope_basal'], 
                    offset_prox=self.P.p['offset_prox'], 
                    offset_dist=self.P.p['offset_dist'], 
                    offset_basal=self.P.p['offset_basal'], 
                    threshold_prox=self.P.p['threshold_prox'], 
                    threshold_dist=self.P.p['threshold_dist'], 
                    threshold_basal=self.P.p['threshold_basal'],)
                
                # temporary record
                self.temp_record = h.Vector()

                # set recording from netcon object
                self.nc['1'][0][0]['ampa'].record(self.temp_record)

                # runa nd save
                #----------------
                self.run_and_save(cell=self.cell,)

                self.temp_record = np.array(self.temp_record)

        if 'generate_variables' in kwargs and kwargs['generate_variables']:
            self.generate_group_variables(self, **kwargs)

        # # generate group variables
        # #------------------------------------
        # # vtrace
        # self.generate_group_variable_vtrace()
        # # w_clopath
        # clopath_param = param.ParamClopath().kronberg_2020_temp_1()
        # self.generate_group_variable_w_clopath(clopath_param=clopath_param, input_times_key='data_input_times')

    def setup_parallel_parameters(self, **kwargs):
        '''
        '''
        self.parallel_parameters = []
        n_workers = 10
        trials_per_worker=5
        w_means_list=[[5.0E-3, 5.1E-3], [5.2E-3, 5.3E-3],[5.4E-3, 5.5E-3], [5.6E-3, 5.7E-3], [5.8E-3, 5.9E-3], [6.E-3, 6.1E-3], [6.2E-3, 6.3E-3],[6.4E-3, 6.5E-3], [6.6E-3, 6.7E-3], [6.8E-3, 6.9E-3]]
        for w_means in w_means_list:
            self.parallel_parameters.append(
                {'experiment':self.experiment_name,
                'w_means':w_means,
                'trials':trials_per_worker,
                })
        return self.parallel_parameters

    def run_parallel(self, **kwargs):
        '''
        '''
        # print 'HERE', kwargs
        self.run(**kwargs)

    def generate_group_variables(self, **kwargs):
        '''
        '''
        # generate group variables
        #------------------------------------
        # vtrace
        self.generate_group_variable_vtrace()
        # w_clopath
        clopath_param = param.ParamClopath().kronberg_2020_reduced()
        self.generate_group_variable_w_clopath(input_times_key='data_input_times', clopath_param=clopath_param, **kwargs)

class exp_reduced_neuron_tbs_conjunctive_high_excitability(Exp):
    '''
    added R-type calcium channels and reduced A-type potassium channel gradient to 0.95
    NMDA:AMPA ratio = 3
    v0 for NMDA block = 0

    this yields stronger bursting than exp_reduced_neuron_tbs_conjunctive_high_nmda
    '''
    def __init__(self, **kwargs):
        super(exp_reduced_neuron_tbs_conjunctive_high_excitability, self).__init__(**kwargs)

    def run(self, **kwargs):
        '''
        '''
        # updates to global parameter dictionary 
        p_update = {
        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),],
        'active_paths':['1','2'],
        'tstop':1000,
        'gna_inact':1,
        # 'gna':0.,
        'dgna' : 1.*-.000025,          # change in sodium conductance with distance (ohm/cm2/um) from Kim 2015
        'gcalbar': 1.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        'gcabar_r': 2.*.0003 ,          # r-type calcium conductance from Kim et al. 2015 (mho/cm2)
        'gcatbar_t': 0.*.0004 ,          # t-type calcium conductance from Kim et al. 2015 (mho/cm2)
        'ka_grad': .95,#1.,          
        'alpha_vspom_nmda':-.062, 
        'v0_block_nmda':0#0#-5#0,
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apic'],
        'nmda_ampa_ratio':3.,
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': [],
        'pulses': 5.,
        'pulse_freq': 100.,
        'bursts': 15.,
        'burst_freq': 5.,
        'warmup': 20,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        },
        '2':{
        'trees': ['apic'],
        'nmda_ampa_ratio':5.,
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': [],
        'pulses': 5.,
        'pulse_freq': 100.,
        'bursts': 15.,
        'burst_freq': 5.,
        'warmup': 70,
        'w_mean': 1.7*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0, 
        },}

        # create param object
        self.P = param.ParamMigliore2005()
        # update parameters
        self.P.p['p_path']={}
        self.P.p.update(p_update)
        self.P.p['p_path'].update(paths_update)

        # create reduced cell
        self.cell = neurons.CellMigliore2005Reduced(P=self.P, cell_number=0)
        self.syns = self.cell.syns
        
        
        # record from activated segment
        #------------------------------
        self.P.p['rec_idx']=[('soma', 0, 0),]
        self.P.p['p_path']['1']['syn_idx']=[]

        # apical_groups = [[0,1,2,3,4,],[6,7,8,9,10,]]
        # segment numbers on apical tree
        proximal_group = [0,1,2,3,4,5]
        distal_group = [6,7,8,9,10]
        field_mags = [-5,0,5]
        # number of trials 
        #----------------------
        if 'trials' in kwargs:
            trials=kwargs['trials']
        else:
            trials=1
        # w_mean
        #------------------------
        # w_mean_default = 5.5E-3
        w_mean_distal_default=3.
        w_mean_proximal_default=1.4#1.5
        if 'w_means_proximal' in kwargs:
            w_means_proximal = kwargs['w_means_proximal']
        else:
            w_means_proximal=[w_mean_proximal_default]

        if 'w_means_distal' in kwargs:
            w_means_distal = kwargs['w_means_distal']
        else:
            w_means_distal=[w_mean_distal_default]

        if 'path_offsets' in kwargs:
            path_offsets=kwargs['path_offsets']
        else:
            path_offsets=[0]
        for path_offset in path_offsets:
            self.P.p['p_path']['2']['warmup']=self.P.p['p_path']['1']['warmup']+path_offset
            for w_mean_proximal in w_means_proximal:
                
                for w_mean_distal in w_means_distal:
                    # synaptic weights
                    #------------------------------
                    self.P.p['p_path']['1']['w_mean'] = round(w_mean_distal*.001, 5)
                    self.P.p['p_path']['1']['w_std'] = 0.2*self.P.p['p_path']['1']['w_mean']#1.*.001
                    self.P.p['p_path']['1']['w_rand'] = False
                    if self.P.p['p_path']['1']['w_mean']==0:
                        self.P.p['p_path']['1']['w_rand'] = False
                    # create syn_idx and rec_idx
                    #--------------------------------------
                    self.P.p['syn_num'] = len(distal_group)+len(proximal_group)
                    self.P.p['rec_idx']=[('soma', 0, 0),]
                    # distal group
                    #----------------------------------
                    self.P.p['p_path']['1']['syn_idx']=[]
                    for seg_i in distal_group:
                        rec_loc = ('apic', 0, seg_i)
                        syn_loc = ('apic', 0, seg_i, 0)
                        self.P.p['p_path']['1']['syn_idx'].append(syn_loc)
                        self.P.p['rec_idx'].append(rec_loc)
                    # proximal group
                    #-------------------------------------
                    # synaptic weights
                    #------------------------------
                    self.P.p['p_path']['2']['w_mean'] = round(w_mean_proximal*.001, 5)
                    self.P.p['p_path']['2']['w_std'] = 0.2*self.P.p['p_path']['1']['w_mean']#1.*.001
                    self.P.p['p_path']['2']['w_rand'] = False
                    if self.P.p['p_path']['2']['w_mean']==0:
                        self.P.p['p_path']['2']['w_rand'] = False
                    self.P.p['p_path']['2']['syn_idx']=[]
                    for seg_i in proximal_group:
                        rec_loc = ('apic', 0, seg_i)
                        syn_loc = ('apic', 0, seg_i, 0)

                        self.P.p['p_path']['2']['syn_idx'].append(syn_loc)
                        self.P.p['rec_idx'].append(rec_loc)
                    print 'w_mean 1', self.P.p['p_path']['1']['w_mean']
                    print 'w_mean 2', self.P.p['p_path']['2']['w_mean']
                    # iterate trials
                    for trial in range(trials):
                        # path 1
                        #---------
                        # set synapse counts and weights
                        # one synapse per segment
                        self.P.p['p_path']['1']['syn_counts']=np.ones(len(self.P.p['p_path']['1']['syn_idx']))
                        self.P.p['p_path']['2']['syn_counts']=np.ones(len(self.P.p['p_path']['2']['syn_idx']))
                        # weights from normal distribution
                        self.P.p['p_path']['1']['w_idx'] = self.P._set_weights_normal(self.P.p['p_path']['1'])
                        self.P.p['p_path']['2']['w_idx'] = self.P._set_weights_normal(self.P.p['p_path']['2'])
                        # create inputs
                        #----------------------------
                        self.stim1 = stims.Uncage()._bursts_vecstim(self.P.p['p_path']['1'])
                        self.stim2 = stims.Uncage()._bursts_vecstim(self.P.p['p_path']['2'])
                        # connect synapses
                        #--------------------------------
                        self.nc={}
                        self.nc['1'] = stims.Uncage()._connect_synapses_reduced(p_path=self.P.p['p_path']['1'], stim=self.stim1, syns=self.cell.syns, bursts=False)
                        self.nc['2'] = stims.Uncage()._connect_synapses_reduced(p_path=self.P.p['p_path']['2'], stim=self.stim2, syns=self.cell.syns, bursts=False)

                        # trial_id
                        #-----------------------------------
                        trial_id = self._generate_trial_id()
                        # iterate over field magnitudes
                        #--------------------------------
                        for field in field_mags:
                            # piecewise liner uniform field over basal, proximal apical, and distal apical dendrites empirically matched to poalrization in full neuron
                            self.P = self._setup_piecewise_field(field=field, P=self.P)
                            # insert field
                            #--------------
                            self.dcs = stims.DCS(cell=self.cell, method='piecewise', field_on=5, field_off=self.P.p['tstop'],
                                slope_prox=self.P.p['slope_prox'], 
                                slope_dist=self.P.p['slope_dist'], 
                                slope_basal=self.P.p['slope_basal'], 
                                offset_prox=self.P.p['offset_prox'], 
                                offset_dist=self.P.p['offset_dist'], 
                                offset_basal=self.P.p['offset_basal'], 
                                threshold_prox=self.P.p['threshold_prox'], 
                                threshold_dist=self.P.p['threshold_dist'], 
                                threshold_basal=self.P.p['threshold_basal'],)
                            # runa nd save
                            #----------------
                            self.run_and_save(cell=self.cell,)

        if 'generate_variables' in kwargs and kwargs['generate_variables']:
            self.generate_group_variables(self, **kwargs)

        # # generate group variables
        # #------------------------------------
        # # vtrace
        # self.generate_group_variable_vtrace()
        # # w_clopath
        # clopath_param = param.ParamClopath().kronberg_2020_temp_1()
        # self.generate_group_variable_w_clopath(clopath_param=clopath_param, input_times_key='data_input_times')

    def setup_parallel_parameters(self, **kwargs):
        '''
        '''
        self.parallel_parameters = []
        n_workers = 10
        trials_per_worker=10
        # w_means_proximal_list = np.arange(1.4, 2.4, .1)
        # np.append(w_means_proximal_list, [0.])
        # w_means_distal_list = np.arange(3., 4., .1)
        # np.append(w_means_distal_list, [0.])
        # w_means_list=[[5.0E-3, 5.1E-3], [5.2E-3, 5.3E-3],[5.4E-3, 5.5E-3], [5.6E-3, 5.7E-3], [5.8E-3, 5.9E-3], [6.E-3, 6.1E-3], [6.2E-3, 6.3E-3],[6.4E-3, 6.5E-3], [6.6E-3, 6.7E-3], [6.8E-3, 6.9E-3]]
        path_offsets_list = [0, 20, 40, 60, 80]
        w_means_proximal_list = [0., 1.7]
        w_means_distal_list = [0., 3.2]
        for path_offsets in path_offsets_list:
            for w_means_proximal in w_means_proximal_list:
                self.parallel_parameters.append(
                    {'experiment':self.experiment_name,
                    'w_means_proximal':[w_means_proximal],
                    'w_means_distal':w_means_distal_list,
                    'trials':trials_per_worker,
                    'path_offsets':[path_offsets]
                    })
        return self.parallel_parameters

    def run_parallel(self, **kwargs):
        '''
        '''
        # print 'HERE', kwargs
        self.run(**kwargs)

    def generate_group_variables(self, **kwargs):
        '''
        '''
        # generate group variables
        #------------------------------------
        # vtrace
        self.generate_group_variable_vtrace()
        # w_clopath
        clopath_param = param.ParamClopath().kronberg_2020_reduced()
        self.generate_group_variable_w_clopath(input_times_key='data_input_times', clopath_param=clopath_param, **kwargs)

class exp_reduced_neuron_tbs_conjunctive_high_nmda(Exp):
    '''
    pamaters that bring soma near threhsold:
    apical_groups = [[0,1,2,3,4,5,]]
    w_means= np.arange(5.8E-3, 6.3E-3, 1E-3)
    w_std = 0.1*w_mean
    10 simulations for each value of w_mean 
    '''
    def __init__(self, **kwargs):
        super(exp_reduced_neuron_tbs_conjunctive_high_nmda, self).__init__(**kwargs)

    def run(self, **kwargs):
        '''
        '''
        # updates to global parameter dictionary 
        p_update = {
        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),],
        'active_paths':['1','2'],
        'tstop':600,
        'gna_inact':1,
        # 'gna':0.,
        'dgna' : 1.*-.000025,          # change in sodium conductance with distance (ohm/cm2/um) from Kim 2015
        'gcalbar': 1.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        'ka_grad': 1.,#1.,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        # 'ghd_grad': 0.,#1.,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        'alpha_vspom_nmda':-.062, 
        'v0_block_nmda':5#-5#0,
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apic'],
        'nmda_ampa_ratio':5.,
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': [],
        'pulses': 5.,
        'pulse_freq': 100.,
        'bursts': 15.,
        'burst_freq': 5.,
        'warmup': 20,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        },
        '2':{
        'trees': ['apic'],
        'nmda_ampa_ratio':5.,
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': [],
        'pulses': 5.,
        'pulse_freq': 100.,
        'bursts': 15.,
        'burst_freq': 5.,
        'warmup': 70,
        'w_mean': 1.7*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0, 
        },}

        # create param object
        self.P = param.ParamMigliore2005()
        # update parameters
        self.P.p['p_path']={}
        self.P.p.update(p_update)
        self.P.p['p_path'].update(paths_update)

        # create reduced cell
        self.cell = neurons.CellMigliore2005Reduced(P=self.P, cell_number=0)
        self.syns = self.cell.syns
        
        
        # record from activated segment
        #------------------------------
        self.P.p['rec_idx']=[('soma', 0, 0),]
        self.P.p['p_path']['1']['syn_idx']=[]

        # apical_groups = [[0,1,2,3,4,],[6,7,8,9,10,]]
        # segment numbers on apical tree
        proximal_group = [0,1,2,3,4,5]
        distal_group = [6,7,8,9,10]
        field_mags = [-5,0,5]
        # number of trials 
        #----------------------
        if 'trials' in kwargs:
            trials=kwargs['trials']
        else:
            trials=1
        # w_mean
        #------------------------
        # w_mean_default = 5.5E-3
        w_mean_distal_default=3.
        w_mean_proximal_default=1.7#1.5
        if 'w_means_proximal' in kwargs:
            w_means_proximal = kwargs['w_means_proximal']
        else:
            w_means_proximal=[w_mean_proximal_default]

        if 'w_means_distal' in kwargs:
            w_means_distal = kwargs['w_means_distal']
        else:
            w_means_distal=[w_mean_distal_default]

        if 'path_offsets' in kwargs:
            path_offsets=kwargs['path_offsets']
        else:
            path_offsets=[0]
        for path_offset in path_offsets:
            self.P.p['p_path']['2']['warmup']=self.P.p['p_path']['1']['warmup']+path_offset
            for w_mean_proximal in w_means_proximal:
                
                for w_mean_distal in w_means_distal:
                    # synaptic weights
                    #------------------------------
                    self.P.p['p_path']['1']['w_mean'] = round(w_mean_distal*.001, 5)
                    self.P.p['p_path']['1']['w_std'] = 0.2*self.P.p['p_path']['1']['w_mean']#1.*.001
                    self.P.p['p_path']['1']['w_rand'] = True
                    if self.P.p['p_path']['1']['w_mean']==0:
                        self.P.p['p_path']['1']['w_rand'] = False
                    # create syn_idx and rec_idx
                    #--------------------------------------
                    self.P.p['syn_num'] = len(distal_group)+len(proximal_group)
                    self.P.p['rec_idx']=[('soma', 0, 0),]
                    # distal group
                    #----------------------------------
                    self.P.p['p_path']['1']['syn_idx']=[]
                    for seg_i in distal_group:
                        rec_loc = ('apic', 0, seg_i)
                        syn_loc = ('apic', 0, seg_i, 0)
                        self.P.p['p_path']['1']['syn_idx'].append(syn_loc)
                        self.P.p['rec_idx'].append(rec_loc)
                    # proximal group
                    #-------------------------------------
                    # synaptic weights
                    #------------------------------
                    self.P.p['p_path']['2']['w_mean'] = round(w_mean_proximal*.001, 5)
                    self.P.p['p_path']['2']['w_std'] = 0.2*self.P.p['p_path']['1']['w_mean']#1.*.001
                    self.P.p['p_path']['2']['w_rand'] = True
                    if self.P.p['p_path']['2']['w_mean']==0:
                        self.P.p['p_path']['2']['w_rand'] = False
                    self.P.p['p_path']['2']['syn_idx']=[]
                    for seg_i in proximal_group:
                        rec_loc = ('apic', 0, seg_i)
                        syn_loc = ('apic', 0, seg_i, 0)

                        self.P.p['p_path']['2']['syn_idx'].append(syn_loc)
                        self.P.p['rec_idx'].append(rec_loc)
                    print 'w_mean 1', self.P.p['p_path']['1']['w_mean']
                    print 'w_mean 2', self.P.p['p_path']['2']['w_mean']
                    # iterate trials
                    for trial in range(trials):
                        # path 1
                        #---------
                        # set synapse counts and weights
                        # one synapse per segment
                        self.P.p['p_path']['1']['syn_counts']=np.ones(len(self.P.p['p_path']['1']['syn_idx']))
                        self.P.p['p_path']['2']['syn_counts']=np.ones(len(self.P.p['p_path']['2']['syn_idx']))
                        # weights from normal distribution
                        self.P.p['p_path']['1']['w_idx'] = self.P._set_weights_normal(self.P.p['p_path']['1'])
                        self.P.p['p_path']['2']['w_idx'] = self.P._set_weights_normal(self.P.p['p_path']['2'])
                        # create inputs
                        #----------------------------
                        self.stim1 = stims.Uncage()._bursts_vecstim(self.P.p['p_path']['1'])
                        self.stim2 = stims.Uncage()._bursts_vecstim(self.P.p['p_path']['2'])
                        # connect synapses
                        #--------------------------------
                        self.nc={}
                        self.nc['1'] = stims.Uncage()._connect_synapses_reduced(p_path=self.P.p['p_path']['1'], stim=self.stim1, syns=self.cell.syns, bursts=False)
                        self.nc['2'] = stims.Uncage()._connect_synapses_reduced(p_path=self.P.p['p_path']['2'], stim=self.stim2, syns=self.cell.syns, bursts=False)

                        # trial_id
                        #-----------------------------------
                        trial_id = self._generate_trial_id()
                        # iterate over field magnitudes
                        #--------------------------------
                        for field in field_mags:
                            # piecewise liner uniform field over basal, proximal apical, and distal apical dendrites empirically matched to poalrization in full neuron
                            self.P = self._setup_piecewise_field(field=field, P=self.P)
                            # insert field
                            #--------------
                            self.dcs = stims.DCS(cell=self.cell, method='piecewise', field_on=5, field_off=self.P.p['tstop'],
                                slope_prox=self.P.p['slope_prox'], 
                                slope_dist=self.P.p['slope_dist'], 
                                slope_basal=self.P.p['slope_basal'], 
                                offset_prox=self.P.p['offset_prox'], 
                                offset_dist=self.P.p['offset_dist'], 
                                offset_basal=self.P.p['offset_basal'], 
                                threshold_prox=self.P.p['threshold_prox'], 
                                threshold_dist=self.P.p['threshold_dist'], 
                                threshold_basal=self.P.p['threshold_basal'],)
                            # runa nd save
                            #----------------
                            self.run_and_save(cell=self.cell,)

        if 'generate_variables' in kwargs and kwargs['generate_variables']:
            self.generate_group_variables(self, **kwargs)

        # # generate group variables
        # #------------------------------------
        # # vtrace
        # self.generate_group_variable_vtrace()
        # # w_clopath
        # clopath_param = param.ParamClopath().kronberg_2020_temp_1()
        # self.generate_group_variable_w_clopath(clopath_param=clopath_param, input_times_key='data_input_times')

    def setup_parallel_parameters(self, **kwargs):
        '''
        '''
        self.parallel_parameters = []
        n_workers = 10
        trials_per_worker=10
        # w_means_proximal_list = np.arange(1.4, 2.4, .1)
        # np.append(w_means_proximal_list, [0.])
        # w_means_distal_list = np.arange(3., 4., .1)
        # np.append(w_means_distal_list, [0.])
        # w_means_list=[[5.0E-3, 5.1E-3], [5.2E-3, 5.3E-3],[5.4E-3, 5.5E-3], [5.6E-3, 5.7E-3], [5.8E-3, 5.9E-3], [6.E-3, 6.1E-3], [6.2E-3, 6.3E-3],[6.4E-3, 6.5E-3], [6.6E-3, 6.7E-3], [6.8E-3, 6.9E-3]]
        path_offsets_list = [0, 20, 40, 60, 80]
        w_means_proximal_list = [0., 1.7]
        w_means_distal_list = [0., 3.2]
        for path_offsets in path_offsets_list:
            for w_means_proximal in w_means_proximal_list:
                self.parallel_parameters.append(
                    {'experiment':self.experiment_name,
                    'w_means_proximal':[w_means_proximal],
                    'w_means_distal':w_means_distal_list,
                    'trials':trials_per_worker,
                    'path_offsets':[path_offsets]
                    })
        return self.parallel_parameters

    def run_parallel(self, **kwargs):
        '''
        '''
        # print 'HERE', kwargs
        self.run(**kwargs)

    def generate_group_variables(self, **kwargs):
        '''
        '''
        # generate group variables
        #------------------------------------
        # vtrace
        self.generate_group_variable_vtrace()
        # w_clopath
        clopath_param = param.ParamClopath().kronberg_2020_reduced()
        self.generate_group_variable_w_clopath(input_times_key='data_input_times', clopath_param=clopath_param, **kwargs)

class exp_reduced_neuron_tbs_conjunctive(Exp):
    '''
    pamaters that bring soma near threhsold:
    apical_groups = [[0,1,2,3,4,5,]]
    w_means= np.arange(5.8E-3, 6.3E-3, 1E-3)
    w_std = 0.1*w_mean
    10 simulations for each value of w_mean 
    '''
    def __init__(self, **kwargs):
        super(exp_reduced_neuron_tbs_conjunctive, self).__init__(**kwargs)

    def run(self, **kwargs):
        '''
        '''
        # updates to global parameter dictionary 
        p_update = {
        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),],
        'active_paths':['1','2'],
        'tstop':600,
        'gna_inact':1,
        # 'gna':0.,
        'dgna' : 1.*-.000025,          # change in sodium conductance with distance (ohm/cm2/um) from Kim 2015
        'gcalbar': 1.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        'ka_grad': 1.,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        # 'ghd_grad': 0.,#1.,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        'alpha_vspom_nmda':0.7*-.062, 
        'v0_block_nmda':0,
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apic'],
        'nmda_ampa_ratio':2.,
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': [],
        'pulses': 5.,
        'pulse_freq': 100.,
        'bursts': 15.,
        'burst_freq': 5.,
        'warmup': 20,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        },
        '2':{
        'trees': ['apic'],
        'nmda_ampa_ratio':2.,
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': [],
        'pulses': 5.,
        'pulse_freq': 100.,
        'bursts': 15.,
        'burst_freq': 5.,
        'warmup': 70,
        'w_mean': 1.7*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0, 
        },}

        # create param object
        self.P = param.ParamMigliore2005()
        # update parameters
        self.P.p['p_path']={}
        self.P.p.update(p_update)
        self.P.p['p_path'].update(paths_update)

        # create reduced cell
        self.cell = neurons.CellMigliore2005Reduced(P=self.P, cell_number=0)
        self.syns = self.cell.syns
        
        
        # record from activated segment
        #------------------------------
        self.P.p['rec_idx']=[('soma', 0, 0),]
        self.P.p['p_path']['1']['syn_idx']=[]

        # apical_groups = [[0,1,2,3,4,],[6,7,8,9,10,]]
        # segment numbers on apical tree
        proximal_group = [0,1,2,3,4,5]
        distal_group = [6,7,8,9,10]
        field_mags = [-5,0,5]
        # number of trials 
        #----------------------
        if 'trials' in kwargs:
            trials=kwargs['trials']
        else:
            trials=1
        # w_mean
        #------------------------
        # w_mean_default = 5.5E-3
        w_mean_distal_default=3.2
        w_mean_proximal_default=1.5
        if 'w_means_proximal' in kwargs:
            w_means_proximal = kwargs['w_means_proximal']
        else:
            w_means_proximal=[w_mean_proximal_default]

        if 'w_means_distal' in kwargs:
            w_means_distal = kwargs['w_means_distal']
        else:
            w_means_distal=[w_mean_distal_default]

        if 'path_offsets' in kwargs:
            path_offsets=kwargs['path_offsets']
        else:
            path_offsets=[100]
        for path_offset in path_offsets:
            self.P.p['p_path']['2']['warmup']=self.P.p['p_path']['1']['warmup']+path_offset
            for w_mean_proximal in w_means_proximal:
                
                for w_mean_distal in w_means_distal:
                    # synaptic weights
                    #------------------------------
                    self.P.p['p_path']['1']['w_mean'] = round(w_mean_distal*.001, 5)
                    self.P.p['p_path']['1']['w_std'] = 0.2*self.P.p['p_path']['1']['w_mean']#1.*.001
                    self.P.p['p_path']['1']['w_rand'] = True
                    if self.P.p['p_path']['1']['w_mean']==0:
                        self.P.p['p_path']['1']['w_rand'] = False
                    # create syn_idx and rec_idx
                    #--------------------------------------
                    self.P.p['syn_num'] = len(distal_group)+len(proximal_group)
                    self.P.p['rec_idx']=[('soma', 0, 0),]
                    # distal group
                    #----------------------------------
                    self.P.p['p_path']['1']['syn_idx']=[]
                    for seg_i in distal_group:
                        rec_loc = ('apic', 0, seg_i)
                        syn_loc = ('apic', 0, seg_i, 0)
                        self.P.p['p_path']['1']['syn_idx'].append(syn_loc)
                        self.P.p['rec_idx'].append(rec_loc)
                    # proximal group
                    #-------------------------------------
                    # synaptic weights
                    #------------------------------
                    self.P.p['p_path']['2']['w_mean'] = round(w_mean_proximal*.001, 5)
                    self.P.p['p_path']['2']['w_std'] = 0.2*self.P.p['p_path']['1']['w_mean']#1.*.001
                    self.P.p['p_path']['2']['w_rand'] = True
                    if self.P.p['p_path']['2']['w_mean']==0:
                        self.P.p['p_path']['2']['w_rand'] = False
                    self.P.p['p_path']['2']['syn_idx']=[]
                    for seg_i in proximal_group:
                        rec_loc = ('apic', 0, seg_i)
                        syn_loc = ('apic', 0, seg_i, 0)

                        self.P.p['p_path']['2']['syn_idx'].append(syn_loc)
                        self.P.p['rec_idx'].append(rec_loc)
                    print 'w_mean 1', self.P.p['p_path']['1']['w_mean']
                    print 'w_mean 2', self.P.p['p_path']['2']['w_mean']
                    # iterate trials
                    for trial in range(trials):
                        # path 1
                        #---------
                        # set synapse counts and weights
                        # one synapse per segment
                        self.P.p['p_path']['1']['syn_counts']=np.ones(len(self.P.p['p_path']['1']['syn_idx']))
                        self.P.p['p_path']['2']['syn_counts']=np.ones(len(self.P.p['p_path']['2']['syn_idx']))
                        # weights from normal distribution
                        self.P.p['p_path']['1']['w_idx'] = self.P._set_weights_normal(self.P.p['p_path']['1'])
                        self.P.p['p_path']['2']['w_idx'] = self.P._set_weights_normal(self.P.p['p_path']['2'])
                        # create inputs
                        #----------------------------
                        self.stim1 = stims.Uncage()._bursts_vecstim(self.P.p['p_path']['1'])
                        self.stim2 = stims.Uncage()._bursts_vecstim(self.P.p['p_path']['2'])
                        # connect synapses
                        #--------------------------------
                        self.nc={}
                        self.nc['1'] = stims.Uncage()._connect_synapses_reduced(p_path=self.P.p['p_path']['1'], stim=self.stim1, syns=self.cell.syns, bursts=False)
                        self.nc['2'] = stims.Uncage()._connect_synapses_reduced(p_path=self.P.p['p_path']['2'], stim=self.stim2, syns=self.cell.syns, bursts=False)

                        # trial_id
                        #-----------------------------------
                        trial_id = self._generate_trial_id()
                        # iterate over field magnitudes
                        #--------------------------------
                        for field in field_mags:
                            # piecewise liner uniform field over basal, proximal apical, and distal apical dendrites empirically matched to poalrization in full neuron
                            self.P = self._setup_piecewise_field(field=field, P=self.P)
                            # insert field
                            #--------------
                            self.dcs = stims.DCS(cell=self.cell, method='piecewise', field_on=5, field_off=self.P.p['tstop'],
                                slope_prox=self.P.p['slope_prox'], 
                                slope_dist=self.P.p['slope_dist'], 
                                slope_basal=self.P.p['slope_basal'], 
                                offset_prox=self.P.p['offset_prox'], 
                                offset_dist=self.P.p['offset_dist'], 
                                offset_basal=self.P.p['offset_basal'], 
                                threshold_prox=self.P.p['threshold_prox'], 
                                threshold_dist=self.P.p['threshold_dist'], 
                                threshold_basal=self.P.p['threshold_basal'],)
                            # runa nd save
                            #----------------
                            self.run_and_save(cell=self.cell,)

        if 'generate_variables' in kwargs and kwargs['generate_variables']:
            self.generate_group_variables(self, **kwargs)

        # # generate group variables
        # #------------------------------------
        # # vtrace
        # self.generate_group_variable_vtrace()
        # # w_clopath
        # clopath_param = param.ParamClopath().kronberg_2020_temp_1()
        # self.generate_group_variable_w_clopath(clopath_param=clopath_param, input_times_key='data_input_times')

    def setup_parallel_parameters(self, **kwargs):
        '''
        '''
        self.parallel_parameters = []
        n_workers = 10
        trials_per_worker=5
        # w_means_proximal_list = np.arange(1.4, 2.4, .1)
        # np.append(w_means_proximal_list, [0.])
        # w_means_distal_list = np.arange(3., 4., .1)
        # np.append(w_means_distal_list, [0.])
        # w_means_list=[[5.0E-3, 5.1E-3], [5.2E-3, 5.3E-3],[5.4E-3, 5.5E-3], [5.6E-3, 5.7E-3], [5.8E-3, 5.9E-3], [6.E-3, 6.1E-3], [6.2E-3, 6.3E-3],[6.4E-3, 6.5E-3], [6.6E-3, 6.7E-3], [6.8E-3, 6.9E-3]]
        path_offsets_list = [0, 20, 40, 60, 80]
        w_means_proximal_list = [0., 1.5]
        w_means_distal_list = [0., 3.2]
        for path_offsets in path_offsets_list:
            for w_means_proximal in w_means_proximal_list:
                self.parallel_parameters.append(
                    {'experiment':self.experiment_name,
                    'w_means_proximal':[w_means_proximal],
                    'w_means_distal':w_means_distal_list,
                    'trials':trials_per_worker,
                    'path_offsets':[path_offsets]
                    })
        return self.parallel_parameters

    def run_parallel(self, **kwargs):
        '''
        '''
        # print 'HERE', kwargs
        self.run(**kwargs)

    def generate_group_variables(self, **kwargs):
        '''
        '''
        # generate group variables
        #------------------------------------
        # vtrace
        self.generate_group_variable_vtrace()
        # w_clopath
        clopath_param = param.ParamClopath().kronberg_2020_reduced()
        self.generate_group_variable_w_clopath(input_times_key='data_input_times', clopath_param=clopath_param, **kwargs)

class exp_reduced_neuron_tbs_distal(Exp):
    '''
    pamaters that bring soma near threhsold:
    apical_groups = [[0,1,2,3,4,5,]]
    w_means= np.arange(5.8E-3, 6.3E-3, 1E-3)
    w_std = 0.1*w_mean
    10 simulations for each value of w_mean 
    '''
    def __init__(self, **kwargs):
        super(exp_reduced_neuron_tbs_distal, self).__init__(**kwargs)

    def run(self, **kwargs):
        '''
        '''
        # updates to global parameter dictionary 
        p_update = {
        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),],
        'active_paths':['1',],
        'tstop':490,
        'gna_inact':1,
        'dgna' : -.000025,          # change in sodium conductance with distance (ohm/cm2/um) from Kim 2015
        'gcalbar': 1.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apic'],
        'nmda_ampa_ratio':1.,
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': [],
        'pulses': 4.,
        'pulse_freq': 100.,
        'bursts': 15.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        'mod_freq':5,
        'mod_amp':10,
        'mean_rate':10,
        # 'warmup':10,
        'dt':0.025,
        'tstop':100, 
        },}
        # create param object
        self.P = param.ParamMigliore2005()
        # update parameters
        self.P.p['p_path']={}
        self.P.p.update(p_update)
        self.P.p['p_path'].update(paths_update)

        # create reduced cell
        self.cell = neurons.CellMigliore2005Reduced(P=self.P, cell_number=0)
        self.syns = self.cell.syns
        
        
        # record from activated segment
        #------------------------------
        self.P.p['rec_idx']=[('soma', 0, 0),]
        self.P.p['p_path']['1']['syn_idx']=[]

        # apical_groups = [[0,1,2,3,4,],[6,7,8,9,10,]]
        # segment numbers on apical tree
        apical_groups = [[6,7,8,9,10]]
        field_mags = [-5,0,5]
        # number of trials 
        #----------------------
        if 'trials' in kwargs:
            trials=kwargs['trials']
        else:
            trials=1
        # w_mean
        #------------------------
        # w_mean_default = 5.5E-3
        w_mean_default=20.E-3
        if 'w_means' in kwargs:
            w_means = kwargs['w_means']
        else:
            w_means=[w_mean_default]
        for group in apical_groups:
            
            for w_mean in w_means:
                # synaptic weights
                #------------------------------
                self.P.p['p_path']['1']['w_mean'] = w_mean
                self.P.p['p_path']['1']['w_std'] = 0.1*self.P.p['p_path']['1']['w_mean']#1.*.001
                self.P.p['p_path']['1']['w_rand'] = False
                # create syn_idx and rec_idx
                #--------------------------------------
                self.P.p['syn_num'] = len(group)
                self.P.p['rec_idx']=[('soma', 0, 0),]
                self.P.p['p_path']['1']['syn_idx']=[]
                for seg_i in group:
                    rec_loc = ('apic', 0, seg_i)
                    syn_loc = ('apic', 0, seg_i, 0)

                    self.P.p['p_path']['1']['syn_idx'].append(syn_loc)
                    self.P.p['rec_idx'].append(rec_loc)
                # iterate trials
                for trial in range(trials):
                    # set synapse counts and weights
                    # one synapse per segment
                    self.P.p['p_path']['1']['syn_counts']=np.ones(len(self.P.p['p_path']['1']['syn_idx']))
                    # weights from normal distribution
                    self.P.p['p_path']['1']['w_idx'] = self.P._set_weights_normal(self.P.p['p_path']['1'])
                    # create inputs
                    #----------------------------
                    self.stim = stims.Uncage()._bursts_vecstim(self.P.p['p_path']['1'])
                    # connect synapses
                    #--------------------------------
                    self.nc={}
                    self.nc['1'] = stims.Uncage()._connect_synapses_reduced(p_path=self.P.p['p_path']['1'], stim=self.stim, syns=self.cell.syns, bursts=False)
                    # trial_id
                    #-----------------------------------
                    trial_id = self._generate_trial_id()
                    # iterate over field magnitudes
                    #--------------------------------
                    for field in field_mags:
                        # piecewise liner uniform field over basal, proximal apical, and distal apical dendrites empirically matched to poalrization in full neuron
                        self.P = self._setup_piecewise_field(field=field, P=self.P)
                        # insert field
                        #--------------
                        self.dcs = stims.DCS(cell=self.cell, method='piecewise', field_on=5, field_off=self.P.p['tstop'],
                            slope_prox=self.P.p['slope_prox'], 
                            slope_dist=self.P.p['slope_dist'], 
                            slope_basal=self.P.p['slope_basal'], 
                            offset_prox=self.P.p['offset_prox'], 
                            offset_dist=self.P.p['offset_dist'], 
                            offset_basal=self.P.p['offset_basal'], 
                            threshold_prox=self.P.p['threshold_prox'], 
                            threshold_dist=self.P.p['threshold_dist'], 
                            threshold_basal=self.P.p['threshold_basal'],)
                        # runa nd save
                        #----------------
                        self.run_and_save(cell=self.cell,)

        if 'generate_variables' in kwargs and kwargs['generate_variables']:
            self.generate_group_variables(self, **kwargs)

        # # generate group variables
        # #------------------------------------
        # # vtrace
        # self.generate_group_variable_vtrace()
        # # w_clopath
        # clopath_param = param.ParamClopath().kronberg_2020_temp_1()
        # self.generate_group_variable_w_clopath(clopath_param=clopath_param, input_times_key='data_input_times')

    def setup_parallel_parameters(self, **kwargs):
        '''
        '''
        self.parallel_parameters = []
        n_workers = 10
        trials_per_worker=5
        w_means_list=[[5.0E-3, 5.1E-3], [5.2E-3, 5.3E-3],[5.4E-3, 5.5E-3], [5.6E-3, 5.7E-3], [5.8E-3, 5.9E-3], [6.E-3, 6.1E-3], [6.2E-3, 6.3E-3],[6.4E-3, 6.5E-3], [6.6E-3, 6.7E-3], [6.8E-3, 6.9E-3]]
        for w_means in w_means_list:
            self.parallel_parameters.append(
                {'experiment':self.experiment_name,
                'w_means':w_means,
                'trials':trials_per_worker,
                })
        return self.parallel_parameters

    def run_parallel(self, **kwargs):
        '''
        '''
        # print 'HERE', kwargs
        self.run(**kwargs)

    def generate_group_variables(self, **kwargs):
        '''
        '''
        # generate group variables
        #------------------------------------
        # vtrace
        self.generate_group_variable_vtrace()
        # w_clopath
        clopath_param = param.ParamClopath().kronberg_2020_reduced()
        self.generate_group_variable_w_clopath(input_times_key='data_input_times', clopath_param=clopath_param, **kwargs)

class exp_reduced_neuron_1hz(Exp):
    '''
    pamaters that bring soma near threhsold:
    apical_groups = [[0,1,2,3,4,5,]]
    w_means= np.arange(5.8E-3, 6.3E-3, 1E-3)
    w_std = 0.1*w_mean
    10 simulations for each value of w_mean 

    '''
    def __init__(self, **kwargs):
        super(exp_reduced_neuron_1hz, self).__init__(**kwargs)

    def run(self, **kwargs):
        '''
        '''
        # updates to global parameter dictionary 
        p_update = {
        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),],
        'active_paths':['1',],
        'tstop':12090,
        'gna_inact':1,
        'dgna' : -.000025,          # change in sodium conductance with distance (ohm/cm2/um) from Kim 2015
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apic'],
        'nmda_ampa_ratio':1.,
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': [],
        'pulses': 60.,
        'pulse_freq': 1.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        'mod_freq':5,
        'mod_amp':10,
        'mean_rate':10,
        # 'warmup':10,
        'dt':0.025,
        'tstop':100, 
        },}
        # create param object
        self.P = param.ParamMigliore2005()
        # update parameters
        self.P.p['p_path']={}
        self.P.p.update(p_update)
        self.P.p['p_path'].update(paths_update)

        # create reduced cell
        self.cell = neurons.CellMigliore2005Reduced(P=self.P, cell_number=0)
        self.syns = self.cell.syns
        
        
        # record from activated segment
        #------------------------------
        self.P.p['rec_idx']=[('soma', 0, 0),]
        self.P.p['p_path']['1']['syn_idx']=[]

        # apical_groups = [[0,1,2,3,4,],[6,7,8,9,10,]]
        # segment numbers on apical tree
        apical_groups = [[0,1,2,3,4,5,]]
        field_mags = [-5,0,5]
        # number of trials 
        #----------------------
        if 'trials' in kwargs:
            trials=kwargs['trials']
        else:
            trials=1
        # w_mean
        #------------------------
        w_mean_default = 5.5E-3
        if 'w_means' in kwargs:
            w_means = kwargs['w_means']
        else:
            w_means=[w_mean_default]
        for group in apical_groups:
            
            for w_mean in w_means:
                # synaptic weights
                #------------------------------
                self.P.p['p_path']['1']['w_mean'] = w_mean
                self.P.p['p_path']['1']['w_std'] = 0.1*self.P.p['p_path']['1']['w_mean']#1.*.001
                self.P.p['p_path']['1']['w_rand'] = True
                # create syn_idx and rec_idx
                #--------------------------------------
                self.P.p['syn_num'] = len(group)
                self.P.p['rec_idx']=[('soma', 0, 0),]
                self.P.p['p_path']['1']['syn_idx']=[]
                for seg_i in group:
                    rec_loc = ('apic', 0, seg_i)
                    syn_loc = ('apic', 0, seg_i, 0)

                    self.P.p['p_path']['1']['syn_idx'].append(syn_loc)
                    self.P.p['rec_idx'].append(rec_loc)
                # iterate trials
                for trial in range(trials):
                    # set synapse counts and weights
                    # one synapse per segment
                    self.P.p['p_path']['1']['syn_counts']=np.ones(len(self.P.p['p_path']['1']['syn_idx']))
                    # weights from normal distribution
                    self.P.p['p_path']['1']['w_idx'] = self.P._set_weights_normal(self.P.p['p_path']['1'])
                    # create inputs
                    #----------------------------
                    self.stim = stims.Uncage()._bursts_vecstim(self.P.p['p_path']['1'])
                    # connect synapses
                    #--------------------------------
                    self.nc={}
                    self.nc['1'] = stims.Uncage()._connect_synapses_reduced(p_path=self.P.p['p_path']['1'], stim=self.stim, syns=self.cell.syns, bursts=False)
                    # trial_id
                    #-----------------------------------
                    trial_id = self._generate_trial_id()
                    # iterate over field magnitudes
                    #--------------------------------
                    for field in field_mags:
                        # piecewise liner uniform field over basal, proximal apical, and distal apical dendrites empirically matched to poalrization in full neuron
                        self.P = self._setup_piecewise_field(field=field, P=self.P)
                        # insert field
                        #--------------
                        self.dcs = stims.DCS(cell=self.cell, method='piecewise', field_on=5, field_off=self.P.p['tstop'],
                            slope_prox=self.P.p['slope_prox'], 
                            slope_dist=self.P.p['slope_dist'], 
                            slope_basal=self.P.p['slope_basal'], 
                            offset_prox=self.P.p['offset_prox'], 
                            offset_dist=self.P.p['offset_dist'], 
                            offset_basal=self.P.p['offset_basal'], 
                            threshold_prox=self.P.p['threshold_prox'], 
                            threshold_dist=self.P.p['threshold_dist'], 
                            threshold_basal=self.P.p['threshold_basal'],)
                        # runa nd save
                        #----------------
                        self.run_and_save(cell=self.cell,)

        if 'generate_variables' in kwargs and kwargs['generate_variables']:
            self.generate_group_variables(self, **kwargs)

    def setup_parallel_parameters(self, **kwargs):
        '''
        '''
        self.parallel_parameters = []
        n_workers = 5
        trials_per_worker=1
        # w_means_list=[[5.0E-3, 5.1E-3], [5.2E-3, 5.3E-3],[5.4E-3, 5.5E-3], [5.6E-3, 5.7E-3], [5.8E-3, 5.9E-3], [6.E-3, 6.1E-3], [6.2E-3, 6.3E-3],[6.4E-3, 6.5E-3], [6.6E-3, 6.7E-3], [6.8E-3, 6.9E-3]]
        w_means_list=[[5.8E-3], [5.9E-3], [6.E-3], [6.1E-3], [6.2E-3],]
        for worker_i in range(n_workers):
            self.parallel_parameters.append(
                {'experiment':self.experiment_name,
                'w_means':[6E-3],
                'trials':trials_per_worker,
                })
        return self.parallel_parameters

    def run_parallel(self, **kwargs):
        '''
        '''
        # print 'HERE', kwargs
        self.run(**kwargs)

    def generate_group_variables(self, **kwargs):
        '''
        '''
        # generate group variables
        #------------------------------------
        # vtrace
        self.generate_group_variable_vtrace()
        # w_clopath
        clopath_param = param.ParamClopath().kronberg_2020_test_1hz()
        # clopath_param = param.ParamClopath().kronberg_2020_reduced()
        self.generate_group_variable_w_clopath(input_times_key='data_input_times', clopath_param=clopath_param, **kwargs)

class exp_reduced_neuron_20hz_basal(Exp):
    '''
    pamaters that bring soma near threhsold:
    apical_groups = [[0,1,2,3,4,5,]]
    w_means= np.arange(5.8E-3, 6.3E-3, 1E-3)
    w_std = 0.1*w_mean
    10 simulations for each value of w_mean 
    '''
    def __init__(self, **kwargs):
        super(exp_reduced_neuron_20hz_basal, self).__init__(**kwargs)

    def run(self, **kwargs):
        '''
        '''
        # updates to global parameter dictionary 
        p_update = {
        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),],
        'active_paths':['1',],
        'tstop':2900,
        'gna_inact':1,
        'dgna' : -.000025,          # change in sodium conductance with distance (ohm/cm2/um) from Kim 2015
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apic'],
        'nmda_ampa_ratio':1.,
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': [],
        'pulses': 60.,
        'pulse_freq': 20.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        'mod_freq':5,
        'mod_amp':10,
        'mean_rate':10,
        # 'warmup':10,
        'dt':0.025,
        'tstop':100, 
        },}
        # create param object
        self.P = param.ParamMigliore2005()
        # update parameters
        self.P.p['p_path']={}
        self.P.p.update(p_update)
        self.P.p['p_path'].update(paths_update)

        # create reduced cell
        self.cell = neurons.CellMigliore2005Reduced(P=self.P, cell_number=0)
        self.syns = self.cell.syns
        
        
        # record from activated segment
        #------------------------------
        self.P.p['rec_idx']=[('soma', 0, 0),]
        self.P.p['p_path']['1']['syn_idx']=[]

        # apical_groups = [[0,1,2,3,4,],[6,7,8,9,10,]]
        # segment numbers on apical tree
        apical_groups = [[0,1,2,3,4,5,]]
        basal_groups = [[0,1,2,]]
        field_mags = [-5,0,5]
        # number of trials 
        #----------------------
        if 'trials' in kwargs:
            trials=kwargs['trials']
        else:
            trials=1
        # w_mean
        #------------------------
        w_mean_default = 5.5E-3
        if 'w_means' in kwargs:
            w_means = kwargs['w_means']
        else:
            w_means=[w_mean_default]
        # iterate over groups of synapses
        #-----------------------------------
        for group in basal_groups:
            # iterate over mean synaptic weights
            #-----------------------------------
            for w_mean in w_means:
                # synaptic weights
                #------------------------------
                self.P.p['p_path']['1']['w_mean'] = w_mean
                self.P.p['p_path']['1']['w_std'] = 0.1*self.P.p['p_path']['1']['w_mean']#1.*.001
                self.P.p['p_path']['1']['w_rand'] = True
                # create syn_idx and rec_idx
                #--------------------------------------
                self.P.p['syn_num'] = len(group)
                self.P.p['rec_idx']=[('soma', 0, 0),]
                self.P.p['p_path']['1']['syn_idx']=[]
                for seg_i in group:
                    rec_loc = ('dend', 1, seg_i)
                    syn_loc = ('dend', 1, seg_i, 0)

                    self.P.p['p_path']['1']['syn_idx'].append(syn_loc)
                    self.P.p['rec_idx'].append(rec_loc)
                # iterate trials
                #----------------
                for trial in range(trials):
                    # set synapse counts and weights
                    #-------------------------------
                    # one synapse per segment
                    self.P.p['p_path']['1']['syn_counts']=np.ones(len(self.P.p['p_path']['1']['syn_idx']))
                    # weights from normal distribution
                    self.P.p['p_path']['1']['w_idx'] = self.P._set_weights_normal(self.P.p['p_path']['1'])
                    # create inputs from stimulator
                    #----------------------------
                    self.stim = stims.Uncage()._bursts_vecstim(self.P.p['p_path']['1'])
                    # connect synapses
                    #--------------------------------
                    self.nc={}
                    self.nc['1'] = stims.Uncage()._connect_synapses_reduced(p_path=self.P.p['p_path']['1'], stim=self.stim, syns=self.cell.syns, bursts=False)
                    # trial_id
                    #-----------------------------------
                    trial_id = self._generate_trial_id()
                    # iterate over field magnitudes
                    #--------------------------------
                    for field in field_mags:
                        # piecewise liner uniform field over basal, proximal apical, and distal apical dendrites empirically matched to poalrization in full neuron
                        self.P = self._setup_piecewise_field(field=field, P=self.P)
                        # insert field
                        #--------------
                        self.dcs = stims.DCS(cell=self.cell, method='piecewise', field_on=5, field_off=self.P.p['tstop'],
                            slope_prox=self.P.p['slope_prox'], 
                            slope_dist=self.P.p['slope_dist'], 
                            slope_basal=self.P.p['slope_basal'], 
                            offset_prox=self.P.p['offset_prox'], 
                            offset_dist=self.P.p['offset_dist'], 
                            offset_basal=self.P.p['offset_basal'], 
                            threshold_prox=self.P.p['threshold_prox'], 
                            threshold_dist=self.P.p['threshold_dist'], 
                            threshold_basal=self.P.p['threshold_basal'],)
                        # runa nd save
                        #----------------
                        self.run_and_save(cell=self.cell,)

        if 'generate_variables' in kwargs and kwargs['generate_variables']:
            self.generate_group_variables(self, **kwargs)

        # # generate group variables
        # #------------------------------------
        # # vtrace
        # self.generate_group_variable_vtrace()
        # # w_clopath
        # clopath_param = param.ParamClopath().kronberg_2020_temp_1()
        # self.generate_group_variable_w_clopath(clopath_param=clopath_param, input_times_key='data_input_times')

    def setup_parallel_parameters(self, **kwargs):
        '''
        '''
        self.parallel_parameters = []
        n_workers = 5
        trials_per_worker=10
        # w_means_list=[[5.0E-3, 5.1E-3], [5.2E-3, 5.3E-3],[5.4E-3, 5.5E-3], [5.6E-3, 5.7E-3], [5.8E-3, 5.9E-3], [6.E-3, 6.1E-3], [6.2E-3, 6.3E-3],[6.4E-3, 6.5E-3], [6.6E-3, 6.7E-3], [6.8E-3, 6.9E-3]]
        w_means_list=[[5.8E-3], [5.9E-3], [6.E-3], [6.1E-3], [6.2E-3],]
        for w_means in w_means_list:
            self.parallel_parameters.append(
                {'experiment':self.experiment_name,
                'w_means':w_means,
                'trials':trials_per_worker,
                })
        return self.parallel_parameters

    def run_parallel(self, **kwargs):
        '''
        '''
        # print 'HERE', kwargs
        self.run(**kwargs)

    def generate_group_variables(self, **kwargs):
        '''
        '''
        # generate group variables
        #------------------------------------
        # vtrace
        self.generate_group_variable_vtrace()
        # w_clopath
        clopath_param = param.ParamClopath().kronberg_2020_reduced()
        self.generate_group_variable_w_clopath(input_times_key='data_input_times', clopath_param=clopath_param, **kwargs)

class exp_reduced_neuron_tbs_basal(Exp):
    '''
    pamaters that bring soma near threhsold:
    apical_groups = [[0,1,2,3,4,5,]]
    w_means= np.arange(5.8E-3, 6.3E-3, 1E-3)
    w_std = 0.1*w_mean
    10 simulations for each value of w_mean 
    '''
    def __init__(self, **kwargs):
        super(exp_reduced_neuron_tbs_basal, self).__init__(**kwargs)

    def run(self, **kwargs):
        '''
        '''
        # updates to global parameter dictionary 
        p_update = {
        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),],
        'active_paths':['1',],
        'tstop':2900,
        'gna_inact':1,
        'dgna' : -.000025,          # change in sodium conductance with distance (ohm/cm2/um) from Kim 2015
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apic'],
        'nmda_ampa_ratio':1.,
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': [],
        'pulses': 4.,
        'pulse_freq': 100.,
        'bursts': 15.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        'mod_freq':5,
        'mod_amp':10,
        'mean_rate':10,
        # 'warmup':10,
        'dt':0.025,
        'tstop':100, 
        },}
        # create param object
        self.P = param.ParamMigliore2005()
        # update parameters
        self.P.p['p_path']={}
        self.P.p.update(p_update)
        self.P.p['p_path'].update(paths_update)

        # create reduced cell
        self.cell = neurons.CellMigliore2005Reduced(P=self.P, cell_number=0)
        self.syns = self.cell.syns
        
        
        # record from activated segment
        #------------------------------
        self.P.p['rec_idx']=[('soma', 0, 0),]
        self.P.p['p_path']['1']['syn_idx']=[]

        # apical_groups = [[0,1,2,3,4,],[6,7,8,9,10,]]
        # segment numbers on apical tree
        apical_groups = [[0,1,2,3,4,5,]]
        basal_groups = [[0,1,2,]]
        field_mags = [-5,0,5]
        # number of trials 
        #----------------------
        if 'trials' in kwargs:
            trials=kwargs['trials']
        else:
            trials=1
        # w_mean
        #------------------------
        w_mean_default = 5.5E-3
        if 'w_means' in kwargs:
            w_means = kwargs['w_means']
        else:
            w_means=[w_mean_default]
        # iterate over groups of synapses
        #-----------------------------------
        for group in basal_groups:
            # iterate over mean synaptic weights
            #-----------------------------------
            for w_mean in w_means:
                # synaptic weights
                #------------------------------
                self.P.p['p_path']['1']['w_mean'] = w_mean
                self.P.p['p_path']['1']['w_std'] = 0.1*self.P.p['p_path']['1']['w_mean']#1.*.001
                self.P.p['p_path']['1']['w_rand'] = True
                # create syn_idx and rec_idx
                #--------------------------------------
                self.P.p['syn_num'] = len(group)
                self.P.p['rec_idx']=[('soma', 0, 0),]
                self.P.p['p_path']['1']['syn_idx']=[]
                for seg_i in group:
                    rec_loc = ('dend', 1, seg_i)
                    syn_loc = ('dend', 1, seg_i, 0)

                    self.P.p['p_path']['1']['syn_idx'].append(syn_loc)
                    self.P.p['rec_idx'].append(rec_loc)
                # iterate trials
                #----------------
                for trial in range(trials):
                    # set synapse counts and weights
                    #-------------------------------
                    # one synapse per segment
                    self.P.p['p_path']['1']['syn_counts']=np.ones(len(self.P.p['p_path']['1']['syn_idx']))
                    # weights from normal distribution
                    self.P.p['p_path']['1']['w_idx'] = self.P._set_weights_normal(self.P.p['p_path']['1'])
                    # create inputs from stimulator
                    #----------------------------
                    self.stim = stims.Uncage()._bursts_vecstim(self.P.p['p_path']['1'])
                    # connect synapses
                    #--------------------------------
                    self.nc={}
                    self.nc['1'] = stims.Uncage()._connect_synapses_reduced(p_path=self.P.p['p_path']['1'], stim=self.stim, syns=self.cell.syns, bursts=False)
                    # trial_id
                    #-----------------------------------
                    trial_id = self._generate_trial_id()
                    # iterate over field magnitudes
                    #--------------------------------
                    for field in field_mags:
                        # piecewise liner uniform field over basal, proximal apical, and distal apical dendrites empirically matched to poalrization in full neuron
                        self.P = self._setup_piecewise_field(field=field, P=self.P)
                        # insert field
                        #--------------
                        self.dcs = stims.DCS(cell=self.cell, method='piecewise', field_on=5, field_off=self.P.p['tstop'],
                            slope_prox=self.P.p['slope_prox'], 
                            slope_dist=self.P.p['slope_dist'], 
                            slope_basal=self.P.p['slope_basal'], 
                            offset_prox=self.P.p['offset_prox'], 
                            offset_dist=self.P.p['offset_dist'], 
                            offset_basal=self.P.p['offset_basal'], 
                            threshold_prox=self.P.p['threshold_prox'], 
                            threshold_dist=self.P.p['threshold_dist'], 
                            threshold_basal=self.P.p['threshold_basal'],)
                        # runa nd save
                        #----------------
                        self.run_and_save(cell=self.cell,)

        if 'generate_variables' in kwargs and kwargs['generate_variables']:
            self.generate_group_variables(self, **kwargs)

        # # generate group variables
        # #------------------------------------
        # # vtrace
        # self.generate_group_variable_vtrace()
        # # w_clopath
        # clopath_param = param.ParamClopath().kronberg_2020_temp_1()
        # self.generate_group_variable_w_clopath(clopath_param=clopath_param, input_times_key='data_input_times')

    def setup_parallel_parameters(self, **kwargs):
        '''
        '''
        self.parallel_parameters = []
        n_workers = 5
        trials_per_worker=10
        # w_means_list=[[5.0E-3, 5.1E-3], [5.2E-3, 5.3E-3],[5.4E-3, 5.5E-3], [5.6E-3, 5.7E-3], [5.8E-3, 5.9E-3], [6.E-3, 6.1E-3], [6.2E-3, 6.3E-3],[6.4E-3, 6.5E-3], [6.6E-3, 6.7E-3], [6.8E-3, 6.9E-3]]
        w_means_list=[[5.8E-3], [5.9E-3], [6.E-3], [6.1E-3], [6.2E-3],]
        for w_means in w_means_list:
            self.parallel_parameters.append(
                {'experiment':self.experiment_name,
                'w_means':w_means,
                'trials':trials_per_worker,
                })
        return self.parallel_parameters

    def run_parallel(self, **kwargs):
        '''
        '''
        # print 'HERE', kwargs
        self.run(**kwargs)

    def generate_group_variables(self, **kwargs):
        '''
        '''
        # generate group variables
        #------------------------------------
        # vtrace
        self.generate_group_variable_vtrace()
        # w_clopath
        clopath_param = param.ParamClopath().kronberg_2020_reduced()
        self.generate_group_variable_w_clopath(input_times_key='data_input_times', clopath_param=clopath_param, **kwargs)

class exp_reduced_neuron_tbs(Exp):
    '''
    pamaters that bring soma near threhsold:
    apical_groups = [[0,1,2,3,4,5,]]
    w_means= np.arange(5.8E-3, 6.3E-3, 1E-3)
    w_std = 0.1*w_mean
    10 simulations for each value of w_mean 
    '''
    def __init__(self, **kwargs):
        super(exp_reduced_neuron_tbs, self).__init__(**kwargs)

    def run(self, **kwargs):
        '''
        '''
        # updates to global parameter dictionary 
        p_update = {
        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),],
        'active_paths':['1',],
        'tstop':2900,
        'gna_inact':1,
        'dgna' : -.000025,          # change in sodium conductance with distance (ohm/cm2/um) from Kim 2015
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apic'],
        'nmda_ampa_ratio':1.,
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': [],
        'pulses': 4.,
        'pulse_freq': 100.,
        'bursts': 15.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        'mod_freq':5,
        'mod_amp':10,
        'mean_rate':10,
        # 'warmup':10,
        'dt':0.025,
        'tstop':100, 
        },}
        # create param object
        self.P = param.ParamMigliore2005()
        # update parameters
        self.P.p['p_path']={}
        self.P.p.update(p_update)
        self.P.p['p_path'].update(paths_update)

        # create reduced cell
        self.cell = neurons.CellMigliore2005Reduced(P=self.P, cell_number=0)
        self.syns = self.cell.syns
        
        
        # record from activated segment
        #------------------------------
        self.P.p['rec_idx']=[('soma', 0, 0),]
        self.P.p['p_path']['1']['syn_idx']=[]

        # apical_groups = [[0,1,2,3,4,],[6,7,8,9,10,]]
        # segment numbers on apical tree
        apical_groups = [[0,1,2,3,4,5,]]
        field_mags = [-5,0,5]
        # number of trials 
        #----------------------
        if 'trials' in kwargs:
            trials=kwargs['trials']
        else:
            trials=1
        # w_mean
        #------------------------
        w_mean_default = 5.5E-3
        if 'w_means' in kwargs:
            w_means = kwargs['w_means']
        else:
            w_means=[w_mean_default]
        for group in apical_groups:
            
            for w_mean in w_means:
                # synaptic weights
                #------------------------------
                self.P.p['p_path']['1']['w_mean'] = w_mean
                self.P.p['p_path']['1']['w_std'] = 0.1*self.P.p['p_path']['1']['w_mean']#1.*.001
                self.P.p['p_path']['1']['w_rand'] = True
                # create syn_idx and rec_idx
                #--------------------------------------
                self.P.p['syn_num'] = len(group)
                self.P.p['rec_idx']=[('soma', 0, 0),]
                self.P.p['p_path']['1']['syn_idx']=[]
                for seg_i in group:
                    rec_loc = ('apic', 0, seg_i)
                    syn_loc = ('apic', 0, seg_i, 0)

                    self.P.p['p_path']['1']['syn_idx'].append(syn_loc)
                    self.P.p['rec_idx'].append(rec_loc)
                # iterate trials
                for trial in range(trials):
                    # set synapse counts and weights
                    # one synapse per segment
                    self.P.p['p_path']['1']['syn_counts']=np.ones(len(self.P.p['p_path']['1']['syn_idx']))
                    # weights from normal distribution
                    self.P.p['p_path']['1']['w_idx'] = self.P._set_weights_normal(self.P.p['p_path']['1'])
                    # create inputs
                    #----------------------------
                    self.stim = stims.Uncage()._bursts_vecstim(self.P.p['p_path']['1'])
                    # connect synapses
                    #--------------------------------
                    self.nc={}
                    self.nc['1'] = stims.Uncage()._connect_synapses_reduced(p_path=self.P.p['p_path']['1'], stim=self.stim, syns=self.cell.syns, bursts=False)
                    # trial_id
                    #-----------------------------------
                    trial_id = self._generate_trial_id()
                    # iterate over field magnitudes
                    #--------------------------------
                    for field in field_mags:
                        # piecewise liner uniform field over basal, proximal apical, and distal apical dendrites empirically matched to poalrization in full neuron
                        self.P = self._setup_piecewise_field(field=field, P=self.P)
                        # insert field
                        #--------------
                        self.dcs = stims.DCS(cell=self.cell, method='piecewise', field_on=5, field_off=self.P.p['tstop'],
                            slope_prox=self.P.p['slope_prox'], 
                            slope_dist=self.P.p['slope_dist'], 
                            slope_basal=self.P.p['slope_basal'], 
                            offset_prox=self.P.p['offset_prox'], 
                            offset_dist=self.P.p['offset_dist'], 
                            offset_basal=self.P.p['offset_basal'], 
                            threshold_prox=self.P.p['threshold_prox'], 
                            threshold_dist=self.P.p['threshold_dist'], 
                            threshold_basal=self.P.p['threshold_basal'],)
                        # runa nd save
                        #----------------
                        self.run_and_save(cell=self.cell,)

        if 'generate_variables' in kwargs and kwargs['generate_variables']:
            self.generate_group_variables(self, **kwargs)

        # # generate group variables
        # #------------------------------------
        # # vtrace
        # self.generate_group_variable_vtrace()
        # # w_clopath
        # clopath_param = param.ParamClopath().kronberg_2020_temp_1()
        # self.generate_group_variable_w_clopath(clopath_param=clopath_param, input_times_key='data_input_times')

    def setup_parallel_parameters(self, **kwargs):
        '''
        '''
        self.parallel_parameters = []
        n_workers = 10
        trials_per_worker=5
        w_means_list=[[5.0E-3, 5.1E-3], [5.2E-3, 5.3E-3],[5.4E-3, 5.5E-3], [5.6E-3, 5.7E-3], [5.8E-3, 5.9E-3], [6.E-3, 6.1E-3], [6.2E-3, 6.3E-3],[6.4E-3, 6.5E-3], [6.6E-3, 6.7E-3], [6.8E-3, 6.9E-3]]
        for w_means in w_means_list:
            self.parallel_parameters.append(
                {'experiment':self.experiment_name,
                'w_means':w_means,
                'trials':trials_per_worker,
                })
        return self.parallel_parameters

    def run_parallel(self, **kwargs):
        '''
        '''
        # print 'HERE', kwargs
        self.run(**kwargs)

    def generate_group_variables(self, **kwargs):
        '''
        '''
        # generate group variables
        #------------------------------------
        # vtrace
        self.generate_group_variable_vtrace()
        # w_clopath
        clopath_param = param.ParamClopath().kronberg_2020_temp_1()
        self.generate_group_variable_w_clopath(input_times_key='data_input_times', clopath_param=clopath_param, **kwargs)

class exp_reduced_neuron_tbs_quick_test(Exp):
    '''
    pamaters that bring soma near threhsold:
    apical_groups = [[0,1,2,3,4,5,]]
    w_means= np.arange(5.8E-3, 6.3E-3, 1E-3)
    w_std = 0.1*w_mean
    10 simulations for each value of w_mean 
    '''
    def __init__(self, **kwargs):
        super(exp_reduced_neuron_tbs_quick_test, self).__init__(**kwargs)

    def run(self, **kwargs):
        '''
        '''
        # updates to global parameter dictionary 
        p_update = {
        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),],
        'active_paths':['1',],
        'tstop':90,
        'gna_inact':1,
        'dgna' : -.000025,          # change in sodium conductance with distance (ohm/cm2/um) from Kim 2015
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apic'],
        'nmda_ampa_ratio':1.,
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': [],
        'pulses': 4.,
        'pulse_freq': 100.,
        'bursts': 15.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        'mod_freq':5,
        'mod_amp':10,
        'mean_rate':10,
        # 'warmup':10,
        'dt':0.025,
        'tstop':100, 
        },}
        # create param object
        self.P = param.ParamMigliore2005()
        # update parameters
        self.P.p['p_path']={}
        self.P.p.update(p_update)
        self.P.p['p_path'].update(paths_update)

        # create reduced cell
        self.cell = neurons.CellMigliore2005Reduced(P=self.P, cell_number=0)
        self.syns = self.cell.syns
        
        
        # record from activated segment
        #------------------------------
        self.P.p['rec_idx']=[('soma', 0, 0),]
        self.P.p['p_path']['1']['syn_idx']=[]

        # apical_groups = [[0,1,2,3,4,],[6,7,8,9,10,]]
        # segment numbers on apical tree
        apical_groups = [[0,1,2,3,4,5,]]
        field_mags = [-5,0,5]
        # number of trials 
        #----------------------
        if 'trials' in kwargs:
            trials=kwargs['trials']
        else:
            trials=2
        # w_mean
        #------------------------
        # w_mean_default = 5.5E-3
        w_mean_default = 6.E-3
        if 'w_means' in kwargs:
            w_means = kwargs['w_means']
        else:
            w_means=[w_mean_default]
        for group in apical_groups:
            
            for w_mean in w_means:
                # synaptic weights
                #------------------------------
                self.P.p['p_path']['1']['w_mean'] = w_mean
                self.P.p['p_path']['1']['w_std'] = 0.1*self.P.p['p_path']['1']['w_mean']#1.*.001
                self.P.p['p_path']['1']['w_rand'] = True
                # create syn_idx and rec_idx
                #--------------------------------------
                self.P.p['syn_num'] = len(group)
                self.P.p['rec_idx']=[('soma', 0, 0),]
                self.P.p['p_path']['1']['syn_idx']=[]
                for seg_i in group:
                    rec_loc = ('apic', 0, seg_i)
                    syn_loc = ('apic', 0, seg_i, 0)

                    self.P.p['p_path']['1']['syn_idx'].append(syn_loc)
                    self.P.p['rec_idx'].append(rec_loc)
                # iterate trials
                for trial in range(trials):
                    # set synapse counts and weights
                    # one synapse per segment
                    self.P.p['p_path']['1']['syn_counts']=np.ones(len(self.P.p['p_path']['1']['syn_idx']))
                    # weights from normal distribution
                    self.P.p['p_path']['1']['w_idx'] = self.P._set_weights_normal(self.P.p['p_path']['1'])
                    # create inputs
                    #----------------------------
                    self.stim = stims.Uncage()._bursts_vecstim(self.P.p['p_path']['1'])
                    # connect synapses
                    #--------------------------------
                    self.nc={}
                    self.nc['1'] = stims.Uncage()._connect_synapses_reduced(p_path=self.P.p['p_path']['1'], stim=self.stim, syns=self.cell.syns, bursts=False)
                    # trial_id
                    #-----------------------------------
                    trial_id = self._generate_trial_id()
                    # iterate over field magnitudes
                    #--------------------------------
                    for field in field_mags:
                        # piecewise liner uniform field over basal, proximal apical, and distal apical dendrites empirically matched to poalrization in full neuron
                        self.P = self._setup_piecewise_field(field=field, P=self.P)
                        # insert field
                        #--------------
                        self.dcs = stims.DCS(cell=self.cell, method='piecewise', field_on=5, field_off=self.P.p['tstop'],
                            slope_prox=self.P.p['slope_prox'], 
                            slope_dist=self.P.p['slope_dist'], 
                            slope_basal=self.P.p['slope_basal'], 
                            offset_prox=self.P.p['offset_prox'], 
                            offset_dist=self.P.p['offset_dist'], 
                            offset_basal=self.P.p['offset_basal'], 
                            threshold_prox=self.P.p['threshold_prox'], 
                            threshold_dist=self.P.p['threshold_dist'], 
                            threshold_basal=self.P.p['threshold_basal'],)
                        # runa nd save
                        #----------------
                        self.run_and_save(cell=self.cell,)

        if 'generate_variables' in kwargs and kwargs['generate_variables']:
            self.generate_group_variables(self, **kwargs)

        # # generate group variables
        # #------------------------------------
        # # vtrace
        # self.generate_group_variable_vtrace()
        # # w_clopath
        # clopath_param = param.ParamClopath().kronberg_2020_temp_1()
        # self.generate_group_variable_w_clopath(clopath_param=clopath_param, input_times_key='data_input_times')

    def setup_parallel_parameters(self, **kwargs):
        '''
        '''
        self.parallel_parameters = []
        n_workers = 10
        trials_per_worker=5
        w_means_list=[[5.0E-3, 5.1E-3], [5.2E-3, 5.3E-3],[5.4E-3, 5.5E-3], [5.6E-3, 5.7E-3], [5.8E-3, 5.9E-3], [6.E-3, 6.1E-3], [6.2E-3, 6.3E-3],[6.4E-3, 6.5E-3], [6.6E-3, 6.7E-3], [6.8E-3, 6.9E-3]]
        for w_means in w_means_list:
            self.parallel_parameters.append(
                {'experiment':self.experiment_name,
                'w_means':w_means,
                'trials':trials_per_worker,
                })
        return self.parallel_parameters

    def run_parallel(self, **kwargs):
        '''
        '''
        # print 'HERE', kwargs
        self.run(**kwargs)

    def generate_group_variables(self, **kwargs):
        '''
        '''
        # generate group variables
        #------------------------------------
        # vtrace
        self.generate_group_variable_vtrace()
        # w_clopath
        clopath_param = param.ParamClopath().kronberg_2020_temp_1()
        self.generate_group_variable_w_clopath(input_times_key='data_input_times', clopath_param=clopath_param, **kwargs)

class exp_reduced_neuron_20hz(Exp):
    '''
    pamaters that bring soma near threhsold:
    apical_groups = [[0,1,2,3,4,5,]]
    w_means= np.arange(5.8E-3, 6.3E-3, 1E-3)
    w_std = 0.1*w_mean
    10 simulations for each value of w_mean 

    '''
    def __init__(self, **kwargs):
        super(exp_reduced_neuron_20hz, self).__init__(**kwargs)

    def run(self, **kwargs):
        '''
        '''
        # updates to global parameter dictionary 
        p_update = {
        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),],
        'active_paths':['1',],
        'tstop':2900,
        'gna_inact':1,
        'dgna' : -.000025,          # change in sodium conductance with distance (ohm/cm2/um) from Kim 2015
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apic'],
        'nmda_ampa_ratio':1.,
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': [],
        'pulses': 60.,
        'pulse_freq': 20.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        'mod_freq':5,
        'mod_amp':10,
        'mean_rate':10,
        # 'warmup':10,
        'dt':0.025,
        'tstop':100, 
        },}
        # create param object
        self.P = param.ParamMigliore2005()
        # update parameters
        self.P.p['p_path']={}
        self.P.p.update(p_update)
        self.P.p['p_path'].update(paths_update)

        # create reduced cell
        self.cell = neurons.CellMigliore2005Reduced(P=self.P, cell_number=0)
        self.syns = self.cell.syns
        
        
        # record from activated segment
        #------------------------------
        self.P.p['rec_idx']=[('soma', 0, 0),]
        self.P.p['p_path']['1']['syn_idx']=[]

        # apical_groups = [[0,1,2,3,4,],[6,7,8,9,10,]]
        # segment numbers on apical tree
        apical_groups = [[0,1,2,3,4,5,]]
        field_mags = [-5,0,5]
        # number of trials 
        #----------------------
        if 'trials' in kwargs:
            trials=kwargs['trials']
        else:
            trials=1
        # w_mean
        #------------------------
        w_mean_default = 5.5E-3
        if 'w_means' in kwargs:
            w_means = kwargs['w_means']
        else:
            w_means=[w_mean_default]
        for group in apical_groups:
            
            for w_mean in w_means:
                # synaptic weights
                #------------------------------
                self.P.p['p_path']['1']['w_mean'] = w_mean
                self.P.p['p_path']['1']['w_std'] = 0.1*self.P.p['p_path']['1']['w_mean']#1.*.001
                self.P.p['p_path']['1']['w_rand'] = True
                # create syn_idx and rec_idx
                #--------------------------------------
                self.P.p['syn_num'] = len(group)
                self.P.p['rec_idx']=[('soma', 0, 0),]
                self.P.p['p_path']['1']['syn_idx']=[]
                for seg_i in group:
                    rec_loc = ('apic', 0, seg_i)
                    syn_loc = ('apic', 0, seg_i, 0)

                    self.P.p['p_path']['1']['syn_idx'].append(syn_loc)
                    self.P.p['rec_idx'].append(rec_loc)
                # iterate trials
                for trial in range(trials):
                    # set synapse counts and weights
                    # one synapse per segment
                    self.P.p['p_path']['1']['syn_counts']=np.ones(len(self.P.p['p_path']['1']['syn_idx']))
                    # weights from normal distribution
                    self.P.p['p_path']['1']['w_idx'] = self.P._set_weights_normal(self.P.p['p_path']['1'])
                    # create inputs
                    #----------------------------
                    self.stim = stims.Uncage()._bursts_vecstim(self.P.p['p_path']['1'])
                    # connect synapses
                    #--------------------------------
                    self.nc={}
                    self.nc['1'] = stims.Uncage()._connect_synapses_reduced(p_path=self.P.p['p_path']['1'], stim=self.stim, syns=self.cell.syns, bursts=False)
                    # trial_id
                    #-----------------------------------
                    trial_id = self._generate_trial_id()
                    # iterate over field magnitudes
                    #--------------------------------
                    for field in field_mags:
                        # piecewise liner uniform field over basal, proximal apical, and distal apical dendrites empirically matched to poalrization in full neuron
                        self.P = self._setup_piecewise_field(field=field, P=self.P)
                        # insert field
                        #--------------
                        self.dcs = stims.DCS(cell=self.cell, method='piecewise', field_on=5, field_off=self.P.p['tstop'],
                            slope_prox=self.P.p['slope_prox'], 
                            slope_dist=self.P.p['slope_dist'], 
                            slope_basal=self.P.p['slope_basal'], 
                            offset_prox=self.P.p['offset_prox'], 
                            offset_dist=self.P.p['offset_dist'], 
                            offset_basal=self.P.p['offset_basal'], 
                            threshold_prox=self.P.p['threshold_prox'], 
                            threshold_dist=self.P.p['threshold_dist'], 
                            threshold_basal=self.P.p['threshold_basal'],)
                        # runa nd save
                        #----------------
                        self.run_and_save(cell=self.cell,)

        if 'generate_variables' in kwargs and kwargs['generate_variables']:
            self.generate_group_variables(self, **kwargs)

    def setup_parallel_parameters(self, **kwargs):
        '''
        '''
        self.parallel_parameters = []
        n_workers = 10
        trials_per_worker=5
        # w_means_list=[[5.0E-3, 5.1E-3], [5.2E-3, 5.3E-3],[5.4E-3, 5.5E-3], [5.6E-3, 5.7E-3], [5.8E-3, 5.9E-3], [6.E-3, 6.1E-3], [6.2E-3, 6.3E-3],[6.4E-3, 6.5E-3], [6.6E-3, 6.7E-3], [6.8E-3, 6.9E-3]]
        w_means_list=[[5.8E-3], [5.9E-3], [6.E-3], [6.1E-3], [6.2E-3],]
        for w_means in w_means_list:
            self.parallel_parameters.append(
                {'experiment':self.experiment_name,
                'w_means':w_means,
                'trials':trials_per_worker,
                })
        return self.parallel_parameters

    def run_parallel(self, **kwargs):
        '''
        '''
        # print 'HERE', kwargs
        self.run(**kwargs)

    def generate_group_variables(self, **kwargs):
        '''
        '''
        # generate group variables
        #------------------------------------
        # vtrace
        self.generate_group_variable_vtrace()
        # w_clopath
        clopath_param = param.ParamClopath().kronberg_2020_reduced()
        self.generate_group_variable_w_clopath(input_times_key='data_input_times', clopath_param=clopath_param, **kwargs)

class exp_reduced_neuron_polarization_piecewise(Exp):
    '''
    '''
    def __init__(self, **kwargs):
        super(exp_reduced_neuron_polarization_piecewise, self).__init__(**kwargs)

    def run(self, **kwargs):
        '''
        '''
        # updates to global parameter dictionary 
        p_update = {
        'rec_variables':[('v','range','v'),],
        'tstop':70,
        }
        # set up synaptic pathway parameters
        paths_update = {'1':{},}
        # create param object
        self.P = param.ParamMigliore2005()
        # update parameters
        self.P.p['p_path']={}
        self.P.p.update(p_update)
        self.P.p['p_path'].update(paths_update)
        # create reduced cell
        self.cell = neurons.CellMigliore2005Reduced(P=self.P, cell_number=0)
        # rec variables
        self.P.p['rec_variables'] = [('v','range','v')]
        # record from all segments
        self.P.p['rec_idx']=[]
        for tree_key, tree in self.cell.geo.iteritems():
            for sec_i, sec in enumerate(tree):
                for seg_i, seg in enumerate(sec):
                    self.P.p['rec_idx'].append((tree_key, sec_i, seg_i))
        # apply get_polarization  function when generating vtrace_df
        df_funcs = [functions.DfFuncs()._get_polarization]
        df_func_kwargs= [{}]
        # iterate over field magnitudes
        #--------------------------------
        fields=[-5, 0, 5] # proximal apical
        for field in fields: # sets scale of field
            self.P = self._setup_piecewise_field(field=field, P=self.P, )
            print self.P.p['slope_prox']
            # insert field
            #--------------
            self.dcs = stims.DCS(cell=self.cell, method='piecewise', field_on=5, field_off=1000,
                slope_prox=self.P.p['slope_prox'], 
                slope_dist=self.P.p['slope_dist'], 
                slope_basal=self.P.p['slope_basal'], 
                offset_prox=self.P.p['offset_prox'], 
                offset_dist=self.P.p['offset_dist'], 
                offset_basal=self.P.p['offset_basal'], 
                threshold_prox=self.P.p['threshold_prox'], 
                threshold_dist=self.P.p['threshold_dist'], 
                threshold_basal=self.P.p['threshold_basal'],)
            self.run_and_save(cell=self.cell)
        self.generate_group_variable_vtrace(df_funcs=df_funcs, df_func_kwargs=df_func_kwargs)

    def generate_group_variables(self, **kwargs):
        '''
        '''
        # generate group variables
        #------------------------------------
        # vtrace
        #--------
        # df_funcs to ge polarization
        df_funcs = [functions.DfFuncs()._get_polarization]
        df_func_kwargs= [{}]
        # get vtrace
        self.generate_group_variable_vtrace(df_funcs=df_funcs, df_func_kwargs=df_func_kwargs)

class exp_full_neuron_polarization_dose_response(Exp):
    '''
    '''
    def __init__(self, **kwargs):
        super(exp_full_neuron_polarization_dose_response, self).__init__(**kwargs)

    def run(self, **kwargs):
        '''
        '''
        # updates to global parameter dictionary 
        p_update = {
        'rec_variables':[('v','range','v'),],
        'tstop':70,
        }
        # set up synaptic pathway parameters
        paths_update = {'1':{},}
        # create param object
        self.P = param.ParamMigliore2005()
        # update parameters
        self.P.p['p_path']={}
        self.P.p.update(p_update)
        self.P.p['p_path'].update(paths_update)
        # create full cell
        self.cell = neurons.CellMigliore2005(P=self.P, cell_number=0)
        self.syns = self.cell.syns
        # record from all segments
        self.P.p['rec_idx']=[]
        for tree_key, tree in self.cell.geo.iteritems():
            for sec_i, sec in enumerate(tree):
                for seg_i, seg in enumerate(sec):
                    self.P.p['rec_idx'].append((tree_key, sec_i, seg_i))

        trial_id = self._generate_trial_id()
        field_mags = [-20, -10, -5, -2, -1, 0, 1, 2, 5, 10, 20]
        df_funcs = [functions.DfFuncs()._get_polarization]
        df_func_kwargs= [{}]
        for field in field_mags:
            self.P.p['field']=field
            field_components=[0, field, 0]
            # insert field
            #--------------
            self.dcs = stims.DCS(cell=self.cell, cell_type='geo', field_components=field_components, field_on=5, field_off=10000)
            self.run_and_save(cell=self.cell,)

        self.generate_group_variable_vtrace(df_funcs=df_funcs, df_func_kwargs=df_func_kwargs)

class exp_reduced_neuron_polarization(Exp):
    '''
    '''
    def __init__(self, **kwargs):
        super(exp_reduced_neuron_polarization, self).__init__(**kwargs)

    def run(self, **kwargs):
        '''
        '''
        # updates to global parameter dictionary 
        p_update = {
        'rec_variables':[('v','range','v'),],
        'tstop':70,
        }
        # set up synaptic pathway parameters
        paths_update = {'1':{},}
        # create param object
        self.P = param.ParamMigliore2005()
        # update parameters
        self.P.p['p_path']={}
        self.P.p.update(p_update)
        self.P.p['p_path'].update(paths_update)
        # create reduced cell
        self.cell = neurons.CellMigliore2005Reduced(P=self.P, cell_number=0)
        # rec variables
        self.P.p['rec_variables'] = [('v','range','v')]
        # record from all segments
        self.P.p['rec_idx']=[]
        for tree_key, tree in self.cell.geo.iteritems():
            for sec_i, sec in enumerate(tree):
                for seg_i, seg in enumerate(sec):
                    self.P.p['rec_idx'].append((tree_key, sec_i, seg_i))
        # apply get_polarization  function when generating vtrace_df
        df_funcs = [functions.DfFuncs()._get_polarization]
        df_func_kwargs= [{}]
        # iterate over field magnitudes
        #--------------------------------
        field_mags = [-20, -10, -5, -2, -1, 0, 1, 2, 5, 10, 20]
        for field in field_mags:
            self.P.p['field']=field
            field_components=[0, field, 0]
            # insert field
            #--------------
            self.dcs = stims.DCS(cell=self.cell, cell_type='geo', field_components=field_components, field_on=5, field_off=1000)
            self.run_and_save(cell=self.cell)
        self.generate_group_variable_vtrace(df_funcs=df_funcs, df_func_kwargs=df_func_kwargs)

class exp_branch_sequence_milstein(Exp):
    '''
    '''
    def __init__(self, setup=False, run=False, P=None, **kwargs):
        '''
        '''
        # print kwargs
        super(exp_branch_sequence_milstein, self).__init__(**kwargs)
        if run:
            self.run(**kwargs)
        # if setup:
        #     self.setup(**kwargs)
        # if run:
        #     self.run_and_save_df_milstein(P=self.P, cell=self.cell, trial=0)

    def run(self, **kwargs):
        '''
        '''
        # self.P = P
        p_update = {
        'experiment' : 'exp_branch_sequence_milstein', 
        'trials' : 1,
        'dt':0.025,
        'field':[0],
        'rec_variables':[('v','range','v'), ('g', 'syn', 'NMDA_KIN5')],
        'active_paths':['1',],
        'v_init':-68
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apical'],
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': [],
        'pulses': 1.,
        'pulse_freq': 100.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1.,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        'delay':.1, 
        'sequence_direction':'in',
        'nmda_ampa_ratio':3.,
        'syn_limit':12,
        },}

        # stuff to iterate
        #----------------------------------
        directions = ['out',]
        fields=[0,]
        print kwargs
        # print branch_idx
        branch_idx=None#kwargs['branch_idx']
        n_branches=1
        seg_L=10
        sequence_n=5
        sequence_delays=[1,]
        w_means = [20]#kwargs['w_means']#[1, 2, 3, 4, 5]
        nmda_factors = [8]#kwargs['nmda_factors']#[1, 2, 4, 6, 8,]
        # sequence_delays=[4]
        # load cell
        #-------------------------------------------
        # morph_filename = 'EB1-early-bifurcation.swc'
        morph_filename = 'cell26.swc'
        mech_filename='043016 Type A - km2_NMDA_KIN5_Pr'
        mech_dict = mech_dicts.ca1_milstein_high_g_2
        # block sodium channels
        mech_dict['axon']['nax_milstein']['gbar']['value']=0.
        mech_dict['ais']['nax_milstein']['gbar']['value']=0.
        mech_dict['axon_hill']['nax_milstein']['gbar']['value']=0.
        mech_dict['basal']['nas']['gbar']['value']=0.
        mech_dict['soma']['nas']['gbar']['value']=0.
        mech_dict['basal']['nas']['gbar']['slope']=2*-0.0002
        # potassium channels
        mech_dict['basal']['kap']['gkabar']={'value':0.}
        mech_dict['basal']['kad']['gkabar']={'value':0.}
        mech_dict['basal']['kdr']['gkdrbar']['value']=0.
        mech_dict['basal']['h_milstein']['ghbar']['value']=0.
        mech_dict['basal']['calH']['gcalbar']['value']=0.*.00125
        mech_dict['soma']['pas']['g']['value']=.1*3.5e-5
        mech_dict['soma']['cable']['cm']['value']=2

        for w_mean in w_means:
            paths_update['1']['w_mean']=w_mean
            for nmda_factor in nmda_factors:

                mech_dict['basal']['synapse']['NMDA_KIN5']['gmax']['value']=nmda_factor*0.003026
                


                self.cell = spec.CA1_Pyr(morph_filename=morph_filename, mech_filename=mech_filename, mech_dict=mech_dict, full_spines=False)

                # get branch and set synapse locations
                #--------------------------------------------------
                for tree in paths_update['1']['trees']:
                    branches = self.cell.get_terminal_branches(sec_types=tree, max_distance=200)
                    print branch_idx
                    if branch_idx is not None:
                        run_branches= [branch for i, branch in enumerate(branches) if i in branch_idx]
                    else:
                        run_branches = np.random.choice(branches, size=n_branches)

                    print 'run_branches', run_branches
                    # iterate over branches
                    #----------------------------------------------------------
                    for branch in run_branches:
                        self.branch = branch
                        # select branch, adjust nseg, reinit mechanisms
                        #-------------------------------------------------
                        self.cell.set_nseg(branch, seg_L=seg_L)
                        self.cell.update_diameters()
                        self.cell.reinit_mechanisms()
                        # list of locations along section to stimulate
                        #-------------------------------------------------
                        nseg = branch.sec.nseg
                        n_locs=sequence_n
                        if not nseg >=n_locs:
                            continue
                        
                        
                        # iterate over directions
                        #-------------------------------------------------------
                        fig={}
                        for sequence_delay in sequence_delays:
                            plt.close('all')
                            plt.figure(1)
                            plt.figure(2)

                            for direction in directions:
                                paths_update['1']['sequence_direction']=direction
                                locs = np.arange(1, 0, -seg_L/branch.sec.L)[:n_locs]
                                # locs = np.arange(0, 1, seg_L/branch.sec.L)[:n_locs]
                                if direction=='out':
                                    locs = np.flip(locs)
                                syn_locs=[]
                                delays=[]
                                current_delay=0
                                for loc in locs:
                                    delays.append(current_delay)
                                    current_delay+=sequence_delay
                                    syn_locs.append((branch.name, loc))
                                # update parameter dictionaries
                                #--------------------------------------------------------
                                self.P = self.update_parameters(p_update=p_update, paths_update=paths_update, default_p='migliore_2005')
                                self.P.p['delay']=sequence_delay
                                self.P.p['nmda_factor']=nmda_factor
                                self.P.p['leak_conductance'] = mech_dict['soma']['pas']['g']['value']
                                self.P.paths['1']['syn_idx']=syn_locs
                                self.P.paths['1']['syn_counts']= [1 for _loc in syn_locs]
                                self.P.paths['1']['sequence_delays']=delays
                                self.P.paths['1']['w_idx'] = self.P._set_weights_normal(p_path=self.P.paths['1'])
                                self.P.p['branch_name']=self.branch.name
                                self.P.p['branch_tree']=tree
                                syn_types = ['AMPA_KIN', 'NMDA_KIN5']
                                # create stimulators
                                stim = stims.Uncage()._bursts(self.P.paths['1'])

                                # insert synapses
                                #----------------------------------------------------
                                # apical trunk needs to be insert because parameters are inherited from here
                                syns=[]
                                for node in self.cell.trunk:
                                    spec.Synapse(cell=self.cell, node=node, type_list=syn_types,)
                                # insert synapses in the nodes specified by syn_locs
                                for node_i, node_loc in enumerate(syn_locs): 
                                    node_name = node_loc[0]
                                    node = self.cell._node_names[node_name]
                                    loc = node_loc[1]
                                    syns.append(spec.Synapse(cell=self.cell, node=node, type_list=syn_types, source=stim[node_i][0], stochastic=0, loc=loc, weight=self.P.paths['1']['w_mean']))
                                self.cell.init_synaptic_mechanisms()

                                # iterate over field
                                #------------------------------------------------------------
                                for field in fields:
                                    self.P.p['field']=field
                                    dcs=stims.DCS(cell=self.cell, cell_type='milstein', field_components=[0,field,0], field_on=5, field_off=100)
                        

                                    # data and figure folder
                                    #------------------------------------------------------
                                    data_folder = 'Data/'
                                    self.P.p['data_folder'] = data_folder+self.P.p['experiment']+'/'
                                    self.P.p['fig_folder'] =  'png figures/'+self.P.p['experiment']+'/'

                                    # setup recording vectors
                                    #----------------------------------------------------------------
                                    # self.P.p['rec_idx'] = [(self.cell.soma[0].name, 0.),]
                                    print 'record from',branches[0].name
                                    rec_locs = np.arange(0.2,1.1,0.2)
                                    self.P.p['rec_idx']=[]
                                    # for rec_loc in rec_locs:
                                    #     self.P.p['rec_idx'].append((self.branch.name, rec_loc))
                                    self.P.p['rec_idx'].append((self.cell.soma[0].name, 0.))
                                    # self.P.p['rec_idx'] = [(self.branch.name, 1.),]
                                    self.rec = self._setup_recording_vectors_milstein(P=self.P, cell=self.cell)
                                    # update tstop and warmup etc
                                    #-------------------------------------------------------------------
                                    self.update_time_parameters(self.P)
                                    # run and save
                                    #----------------------------------------------------------------
                                    trial=0
                                    self.run_and_save_df_milstein(P=self.P, cell=self.cell, trial=0)

                                    # quick plots
                                    #---------------------
                                    if direction =='out':
                                        color='red'
                                    else:
                                        color='blue'
                                    print color
                                    plt.figure(1)
                                    plt.plot(fncs._2array(self.data.data_v).T, color=color)
                                    plt.figure(2)
                                    plt.plot(fncs._2array(self.data.data_g).T, color=color)
                                    
                                    # clear stimulators and synapses
                                    #----------------------------------------------------
                                    stim=[]
                                    for syn in syns:
                                        syn=[]

                                    self.cell.clear_synapses()
                plt.show(block=False)

    def setup_parallel_parameters(self, **kwargs):
        '''
        '''
        self.parallel_parameters = []
        n_workers = 10
        trials_per_worker=1
        branch_indices = [[0], [1], [2],[3],[4],[5],[6],[7],[8],[9]]
        w_means=[[1,2]]
        nmda_factors = [[1],[2],[3],[4],[6],[8],[10],[12],[14],[16]]
        for w_mean in w_means:
            for nmda_factor in nmda_factors:
                self.parallel_parameters.append(
                    {'experiment':'exp_branch_sequence_milstein',
                    'branch_idx':[7],
                    'w_means':w_mean,
                    'nmda_factors':nmda_factor
                    })
        return self.parallel_parameters

    def run_parallel(self, **kwargs):
        '''
        '''
        # print 'HERE', kwargs
        self.run(**kwargs)
        # self.run_and_save_df_milstein(P=self.P, cell=self.cell, trial=0)

class exp_reduced_neuron_polarization_sigmoid(Exp):
    '''
    '''
    def __init__(self, **kwargs):
        super(exp_reduced_neuron_polarization_sigmoid, self).__init__(**kwargs)

    def run(self, **kwargs):
        '''
        '''
        # updates to global parameter dictionary 
        p_update = {
        'rec_variables':[('v','range','v'),],
        'tstop':70,
        }
        # set up synaptic pathway parameters
        paths_update = {'1':{},}
        # create param object
        self.P = param.ParamMigliore2005()
        # update parameters
        self.P.p['p_path']={}
        self.P.p.update(p_update)
        self.P.p['p_path'].update(paths_update)
        # create reduced cell
        self.cell = neurons.CellMigliore2005Reduced(P=self.P, cell_number=0)
        # rec variables
        self.P.p['rec_variables'] = [('v','range','v')]
        # record from all segments
        self.P.p['rec_idx']=[]
        for tree_key, tree in self.cell.geo.iteritems():
            for sec_i, sec in enumerate(tree):
                for seg_i, seg in enumerate(sec):
                    self.P.p['rec_idx'].append((tree_key, sec_i, seg_i))
        # apply get_polarization  function when generating vtrace_df
        df_funcs = [functions.DfFuncs()._get_polarization]
        df_func_kwargs= [{}]
        # iterate over field magnitudes
        #--------------------------------
        field_mags = [-20, 0, 20]
        field_slopes = [1000.]
        field_thresholds = [-800.]
        field_components=[0,20,0]
        for field in field_mags:
            for field_slope in field_slopes:
                for field_threshold in field_thresholds:
                    self.P.p['field']=field
                    self.P.p['field_slope']=field_slope
                    self.P.p['field_threshold']=field_threshold
                    magnitude_components=[0, field, 0]
                    slope_components=[1, field_slope, 1]
                    threshold_components=[0, field_threshold, 0]

                    # insert field
                    #--------------
                    self.dcs = stims.DCS(cell=self.cell, method='sigmoid', field_components=field_components, field_on=5, field_off=1000, magnitude_components=magnitude_components, slope_components=slope_components, threshold_components=threshold_components, apical_factor=1.5)
                    self.run_and_save(cell=self.cell)
        self.generate_group_variable_vtrace(df_funcs=df_funcs, df_func_kwargs=df_func_kwargs)

class exp_reduced_neuron_polarization_mirror_estimate(Exp):
    ''' map polarization measured in a full neuron model to locations in the corresponding reduced model built by neuron_reduce.  Since multiple compartments in the original model can be mapped to the reduced model, take the average polarization over original compartments.  Then use the mirror estimate (Joucla and Yvert 2009) to calculate the electric field required to achieve this poalrization in the reduced model.  This approach is mainly used to account for proximal apical dendritic branches which are electrotonically close to the soma, but are polarized opposite to the soma.  Simply applying a uniform electric field to the reduced cable leads to proximal compartments that have the same polarization as the soma.
    '''
    def __init__(self, **kwargs):
        super(exp_reduced_neuron_polarization_mirror_estimate, self).__init__(**kwargs)

    def run(self, **kwargs):
        '''
        '''
        # updates to global parameter dictionary 
        p_update = {
        'rec_variables':[('v','range','v'),],
        'tstop':70,
        }
        # set up synaptic pathway parameters
        paths_update = {'1':{},}
        # create param object
        self.P = param.ParamMigliore2005()
        # update parameters
        self.P.p['p_path']={}
        self.P.p.update(p_update)
        self.P.p['p_path'].update(paths_update)
        # create full cell
        self.cell = neurons.CellMigliore2005Reduced(P=self.P, cell_number=0)
        self.full_cell = neurons.CellMigliore2005(cell_number=0)
        # record from all segments
        self.P.p['rec_idx']=[]
        for tree_key, tree in self.cell.geo.iteritems():
            for sec_i, sec in enumerate(tree):
                for seg_i, seg in enumerate(sec):
                    self.P.p['rec_idx'].append((tree_key, sec_i, seg_i))

        ####################################################################
        # check for preloaded mirror estimate
        # otherwise check for full neuron vtrace_df and load
        #####################################################################
        self.vtrace_df_original = functions._load_group_data( directory='Data/exp_full_neuron_polarization_dose_response/', filename='vtrace_df', df=True)
        # get location map between reduced and original cell
        self.loc_map = self.cell._create_loc_map(self.full_cell)
        # delete original cell
        #-------------------------------------------------------
        for _tree_key, _tree in self.full_cell.geo.iteritems():
            for _sec_i, _sec in enumerate(_tree):
                h.delete_section(sec=_sec)
        del self.full_cell
        # apply get_polarization  function when generating vtrace_df
        df_funcs = [functions.DfFuncs()._get_polarization]
        df_func_kwargs= [{}]
        # iterate over field magnitudes
        #--------------------------------
        self.vtrace_df_original = functions._set_index(self.vtrace_df_original, ['location'])
        field_mags = [-20, 0, 20]#[-20, -10, -5, -2, -1, 0, 1, 2, 5, 10, 20]
        for field in field_mags:
            # get polarization map (segments to polarization for a given field magnitude)
            self.polarization_map = self.cell._get_polarization_map(vtrace_df_original=self.vtrace_df_original, loc_map=self.loc_map, field=field)
            # get extracellular voltage map (polarization to extracellular voltage)
            # self.e_map = self.cell._get_mirror_estimate_from_polarization(polarization_map=self.polarization_map, reduced_cell=self.cell, map_method='end_fit', apical_factor=2., basal_factor=1., apical_offset=3., basal_offset=1.)
            self.e_map = self.cell._get_mirror_estimate_from_polarization(polarization_map=self.polarization_map, reduced_cell=self.cell, map_method='end_fit', apical_factor=1., basal_factor=1., apical_offset=1., basal_offset=1.)
            for _k, _v in self.e_map.iteritems():
                if _k[0]=='apic':
                    self.e_map[_k]=1.*_v
                if _k[0]=='dend':
                    self.e_map[_k]=1.*_v

            # print self.e_map

            # set field magnitude 
            self.P.p['field']=field
            # insert dcs from e_map
            self.dcs = stims.DCS(cell=self.cell, e_map=self.e_map, method='from_map', field_on=20, field_off=1000, scale=2.5)
            # run, but don't save until adding columns for original cell data (see below)
            #'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
            self.run_and_save(cell=self.cell, save=False)
            # create original location and polarization columns
            #-------------------------------------------------------------
            # initialize columns
            self.data = functions._initialize_column(self.data, ['location_original_cell', 'polarization_original_cell', 'polarization_original_cell_mean', 'extracellular_voltage'])
            # iterate rows
            for row_i, row in self.data.iterrows():
                # get reduced location
                location = row.location
                # get original locations and mean polarization and add to data df
                #-------------------------------------------------------
                if 'soma' in location:
                    original_locations =  [location]
                elif location in self.loc_map:
                    original_locations = self.loc_map[location]
                else:
                    original_locations=np.nan
                if location in self.polarization_map:
                    original_polarization = self.polarization_map[location]
                    original_polarization_mean = np.mean(original_polarization)
                else:
                    original_polarization = 0.# np.nan
                    original_polarization_mean = 0.#  np.nan
                # get applied extracellular voltage from e_map
                #--------------------------------------------------------
                if type(self.e_map)==dict:
                    extracellular_voltage=self.e_map[location]
                elif type(self.e_map)==list:
                    extracellular_voltage = [_val[1] for _i,_val in enumerate(self.e_map) if _val[0]==location]
                self.data.at[row_i, 'location_original_cell']=original_locations
                self.data.at[row_i, 'polarization_original_cell']=original_polarization
                self.data.at[row_i, 'polarization_original_cell_mean']=original_polarization_mean
                self.data.at[row_i, 'extracellular_voltage']=extracellular_voltage
            # save data 
            #--------------------------------
            self._save_data(data=self.data)


        self.generate_group_variable_vtrace(df_funcs=df_funcs, df_func_kwargs=df_func_kwargs)

class exp_reduced_neuron_polarization_mirror_estimate_test(Exp):
    ''' map polarization measured in a full neuron model to locations in the corresponding reduced model built by neuron_reduce.  Since multiple compartments in the original model can be mapped to the reduced model, take the average polarization over original compartments.  Then use the mirror estimate (Joucla and Yvert 2009) to calculate the electric field required to achieve this poalrization in the reduced model.  This approach is mainly used to account for proximal apical dendritic branches which are electrotonically close to the soma, but are polarized opposite to the soma.  Simply applying a uniform electric field to the reduced cable leads to proximal compartments that have the same polarization as the soma.
    '''
    def __init__(self, **kwargs):
        super(exp_reduced_neuron_polarization_mirror_estimate_test, self).__init__(**kwargs)

    def run(self, **kwargs):
        '''
        '''
        # updates to global parameter dictionary 
        p_update = {
        'rec_variables':[('v','range','v'),],
        'tstop':70,
        }
        # set up synaptic pathway parameters
        paths_update = {'1':{},}
        # create param object
        self.P = param.ParamMigliore2005()
        # update parameters
        self.P.p['p_path']={}
        self.P.p.update(p_update)
        self.P.p['p_path'].update(paths_update)
        # create full cell
        self.cell = neurons.CellMigliore2005Reduced(P=self.P, cell_number=0)
        self.full_cell = neurons.CellMigliore2005(cell_number=0)
        # record from all segments
        self.P.p['rec_idx']=[]
        for tree_key, tree in self.cell.geo.iteritems():
            for sec_i, sec in enumerate(tree):
                for seg_i, seg in enumerate(sec):
                    self.P.p['rec_idx'].append((tree_key, sec_i, seg_i))

        ####################################################################
        # check for preloaded mirror estimate
        # otherwise check for full neuron vtrace_df and load
        #####################################################################
        self.vtrace_df_original = functions._load_group_data( directory='Data/exp_reduced_neuron_polarization/', filename='vtrace_df', df=True)
        # get location map between reduced and original cell
        self.loc_map={}
        for tree_key, tree in self.cell.geo.iteritems():
            for sec_i, sec in enumerate(tree):
                for seg_i, seg in enumerate(sec):
                    _location = (tree_key, sec_i, seg_i)
                    self.loc_map[_location]=_location
        self.loc_map = self.cell._create_loc_map(self.full_cell)

        # apply get_polarization  function when generating vtrace_df
        df_funcs = [functions.DfFuncs()._get_polarization]
        df_func_kwargs= [{}]
        # iterate over field magnitudes
        #--------------------------------
        self.vtrace_df_original = functions._set_index(self.vtrace_df_original, ['field','location',])
        # field_mags = [-20, 0, 20]#[-20, -10, -5, -2, -1, 0, 1, 2, 5, 10, 20]
        # for field in field_mags:
        #     # get polarization map (segments to polarization for a given field magnitude)
        #     self.polarization_map = self.cell._get_polarization_map(vtrace_df_original=self.vtrace_df_original, loc_map=self.loc_map, field=field)
        #     # get extracellular voltage map (polarization to extracellular voltage)
        #     # self.e_map = self.cell._get_mirror_estimate_from_polarization(polarization_map=self.polarization_map, reduced_cell=self.cell, map_method='end_fit', apical_factor=2., basal_factor=1., apical_offset=3., basal_offset=1.)
        #     self.e_map = self.cell._get_mirror_estimate_from_polarization(polarization_map=self.polarization_map, reduced_cell=self.cell, map_method='end_fit', apical_factor=1., basal_factor=1., apical_offset=1., basal_offset=1.)
        #     for _k, _v in self.e_map.iteritems():
        #         if _k[0]=='apic':
        #             self.e_map[_k]=1.*_v
        #         if _k[0]=='dend':
        #             self.e_map[_k]=1.*_v

        #     # print self.e_map

        #     # set field magnitude 
        #     self.P.p['field']=field
        #     # insert dcs from e_map
        #     self.dcs = stims.DCS(cell=self.cell, e_map=self.e_map, method='from_map', field_on=20, field_off=1000, scale=2.5)
        #     # run, but don't save until adding columns for original cell data (see below)
        #     #'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
        #     self.run_and_save(cell=self.cell, save=False)
        #     # create original location and polarization columns
        #     #-------------------------------------------------------------
        #     # initialize columns
        #     self.data = functions._initialize_column(self.data, ['location_original_cell', 'polarization_original_cell', 'polarization_original_cell_mean', 'extracellular_voltage'])
        #     # iterate rows
        #     for row_i, row in self.data.iterrows():
        #         # get reduced location
        #         location = row.location
        #         # get original locations and mean polarization and add to data df
        #         #-------------------------------------------------------
        #         if 'soma' in location:
        #             original_locations =  [location]
        #         elif location in self.loc_map:
        #             original_locations = self.loc_map[location]
        #         else:
        #             original_locations=np.nan
        #         if location in self.polarization_map:
        #             original_polarization = self.polarization_map[location]
        #             original_polarization_mean = np.mean(original_polarization)
        #         else:
        #             original_polarization = 0.# np.nan
        #             original_polarization_mean = 0.#  np.nan
        #         # get applied extracellular voltage from e_map
        #         #--------------------------------------------------------
        #         if type(self.e_map)==dict:
        #             extracellular_voltage=self.e_map[location]
        #         elif type(self.e_map)==list:
        #             extracellular_voltage = [_val[1] for _i,_val in enumerate(self.e_map) if _val[0]==location]
        #         self.data.at[row_i, 'location_original_cell']=original_locations
        #         self.data.at[row_i, 'polarization_original_cell']=original_polarization
        #         self.data.at[row_i, 'polarization_original_cell_mean']=original_polarization_mean
        #         self.data.at[row_i, 'extracellular_voltage']=extracellular_voltage
        #     # save data 
        #     #--------------------------------
        #     self._save_data(data=self.data)


        # self.generate_group_variable_vtrace(df_funcs=df_funcs, df_func_kwargs=df_func_kwargs)

class exp_reduced_neuron_poisson(Exp):
    '''
    '''
    def __init__(self, **kwargs):
        super(exp_reduced_neuron_poisson, self).__init__(**kwargs)

    def run(self, **kwargs):
        '''
        '''
        # updates to global parameter dictionary 
        p_update = {
        'rec_variables':[('v','range','v'),('ica','range','calH')],
        'active_paths':['1',],
        'tstop':1000
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apic'],
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': [],
        'pulses': 4.,
        'pulse_freq': 100.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        'mod_freq':5,
        'mod_amp':10,
        'mean_rate':10,
        'warmup':10,
        'dt':0.025,
        'tstop':1000, 
        },}
        # create reduced cell
        self.cell = neurons.CellMigliore2005Reduced(syns_per_seg=5)
        # parameter object
        self.P = self.cell.P
        self.P.p['p_path']={}
        # update synaptic parameters
        self.P.p.update(p_update)
        self.P.p['p_path'].update(paths_update)
        # choose synapses manually
        #------------------------------
        self.P.p['p_path']['1']['syn_idx'] = [('apic', 0, 10, 0), ('apic', 0, 5, 0),('apic', 0, 2,0), ('dend',0,0,0)]
        self.P.p['p_path']['1']['syn_counts'] = [20, 20, 20,20]
        self.P.p['p_path']['1']['w_idx'] = self.P._set_weights_normal(self.P.p['p_path']['1'])
        # record from activated segment
        #------------------------------
        self.P.p['rec_idx']=[('soma', 0, 0), ('apic', 0, 10), ('apic', 0, 5), ('apic', 0, 1), ('apic', 0, 0)]
        

        # iterate over field magnitudes
        #--------------------------------
        # create inputs
        #----------------------------
        poisson_rates = [50]
        for mean_rate in poisson_rates:
            self.P.p['p_path']['1']['mean_rate']=mean_rate
            self.stim = stims.Uncage()._poisson(self.P.p['p_path']['1'])
            # connect synapses
            self.nc = stims.Uncage()._connect_synapses_reduced(p_path=self.P.p['p_path']['1'], stim=self.stim, syns=self.cell.syns)
            trial_id = self._generate_trial_id()
            field_mags = [-5,]
            for field in field_mags:
                self.P.p['field']=field
                field_components=[0, field, 0]
                # insert field
                #--------------
                self.dcs = stims.DCS(cell=self.cell, cell_type='geo', field_components=field_components, field_on=20, field_off=1000)
                self.run_and_save(cell=self.cell)

class exp_test_neuron_reduce(Exp):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        super(exp_test_neuron_reduce, self).__init__(**kwargs)

    def run(self, **kwargs):
        '''
        '''
        # load mech dictionary
        # specify geometry
        # load cell in milstein framework
        # convert to original geo type
        # insert synapses
        # test neuron reduce 
        # self.P = P
        p_update = {
        'experiment' : 'exp_branch_sequence_milstein', 
        'trials' : 1,
        'dt':0.025,
        'field':[0],
        'rec_variables':[('v','range','v'), ('g', 'syn', 'NMDA_KIN5')],
        'active_paths':['1',],
        'v_init':-68
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apical'],
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': [],
        'pulses': 1.,
        'pulse_freq': 100.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1.,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        'delay':.1, 
        'sequence_direction':'in',
        'nmda_ampa_ratio':3.,
        'syn_limit':12,
        },}
        # load cell
        #-------------------------------------------
        # morph_filename = 'EB1-early-bifurcation.swc'
        morph_filename = 'cell26.swc'
        mech_filename='043016 Type A - km2_NMDA_KIN5_Pr'
        mech_dict = mech_dicts.ca1_milstein_high_g_2
        
        self.cell = spec.CA1_Pyr(morph_filename=morph_filename, mech_filename=mech_filename, mech_dict=mech_dict, full_spines=False)

        self.cell.update_diameters()
        self.cell.reinit_mechanisms()

        # update parameter dictionaries
        #--------------------------------------------------------
        self.P = self.update_parameters(p_update=p_update, paths_update=paths_update, default_p='migliore_2005')

        # insert synapses
        #----------------------------------------------
        syn_types = ['AMPA_KIN', 'NMDA_KIN5']
        
         # insert synapses
        #----------------------------------------------------
        # apical trunk needs to be insert because parameters are inherited from here
        syns=[]
        for node in self.cell.trunk:
            spec.Synapse(cell=self.cell, node=node, type_list=syn_types,)

        # choose synapses
        #-------------------------------------------------
        # list of nodes to insert synapses
        syn_locs = self.cell.choose_syn_locations_random(n=1, sec_types=['basal'], on_spine=True, min_distance=None, max_distance=None, replace=True)
        self.P.paths['1']['syn_idx']=syn_locs
        self.P.paths['1']['syn_counts']= [1 for _loc in syn_locs]
        self.P.paths['1']['sequence_delays']= []
        # create stimulator
        #----------------------------------------------
        # create stimulators
        stim = stims.Uncage()._bursts(self.P.paths['1'])
        self.syns = []
        # insert synapses in the nodes specified by syn_locs
        for node_i, node in enumerate(syn_locs): 
            loc = np.random.rand(1)[0]
            self.syns.append(spec.Synapse(cell=self.cell, node=node, type_list=syn_types, source=stim[0][0], stochastic=0, loc=loc, weight=self.P.paths['1']['w_mean']))
        self.cell.init_synaptic_mechanisms()

        # data and figure folder
        #------------------------------------------------------
        data_folder = 'Data/'
        self.P.p['data_folder'] = data_folder+self.P.p['experiment']+'/'
        self.P.p['fig_folder'] =  'png figures/'+self.P.p['experiment']+'/'

        # setup recording vectors
        #----------------------------------------------------------------
        self.P.p['rec_idx']=[]
        self.P.p['rec_idx'].append((self.cell.soma[0].name, 0.))
        self.rec = self._setup_recording_vectors_milstein(P=self.P, cell=self.cell)

        # update tstop and warmup etc
        #-------------------------------------------------------------------
        self.update_time_parameters(self.P)
        
        # get synapses and netcon objects
        #-------------------------------------------------------------------
        self.cell.get_all_synapses_and_netcons()
        self.all_synapses = self.cell.all_synapses
        self.all_netcons = self.cell.all_netcons

        # get OG hoc style cell
        #-------------------------------------------------------------------
        self.hoccell = spec.HocStyleCell()
        self.hoccell.from_milstein(self.cell)

        # run neuron_reduce
        #--------------------------------------------------------------------
        #apply Neuron_Reduce to simplify the cell
        self.reduced_output = neuron_reduce.subtree_reductor(
            original_cell=self.hoccell,
            synapses_list=self.all_synapses,
            netcons_list=self.all_netcons,
            reduction_frequency=0, 
            total_segments_manual=-1)

        # run and save
        #----------------------------------------------------------------
        # trial=0
        # self.run_and_save_df_milstein(P=self.P, cell=self.cell, trial=0)

    def run_2(self, **kwargs):
        '''
        '''
        # parameter instance
        self.P = param.Param()
        # load migliore parameters
        self.P.migliore_2005()
        self.cell = neurons.CellMigliore2005(p=self.P.p)
        # get OG hoc style cell
        #-------------------------------------------------------------------
        self.hoccell = spec.HocStyleCell()
        self.hoccell.from_geo(self.cell)
        # self.hoccell.soma = self.hoccell.soma[0]
        # list of synapses
        #-------------------------------------------------------------------
        self.synapses = []
        # for tree_key, tree in self.cell.syns.iteritems():
        #     for sec_i, sec in enumerate(tree):
        #         for seg_i, seg in enumerate(sec):
        #             for syn_type, syn in seg.iteritems():
        #                 self.synapses.append(syn)
        # list of netcons
        #--------------------------------------------------------------------
        self.netcons =[]

        # run neuron_reduce
        #--------------------------------------------------------------------
        # apply Neuron_Reduce to simplify the cell
        self.reduced_output = neuron_reduce.subtree_reductor(
            original_cell=self.hoccell,
            synapses_list=self.synapses,
            netcons_list=self.netcons,
            reduction_frequency=0, 
            total_segments_manual=-1)

# parallel functions
#----------------------------------------------------------------------------
# function to pass to parallel context message board
def f_parallel(parameters):
    """ Wrap experiment function so it exists in global namespace

    Arguments: 
    parameters - dictionary with entries 'experiment' and parameters to be passed to Experiment.exp_num.  'experiment' should be of the form 'exp_4a' and specifies which experiment to run
    """
    # get experiment info
    experiment = parameters['experiment']
    
    # get copy of globals dictionary
    g = globals().copy()
    # get experiment class from globals
    experiment_class = g[experiment]
    # create experiment class instance
    experiment_instance = experiment_class()
    # get run function from experiment class
    f = experiment_instance.run_parallel

    # # get specific experiment function
    # f = getattr(exp_instance, experiment)

    # print f
    # print parameters
    # run experiment
    return f(**parameters)

# function for controlling parallel    
def run_parallel(parameters):
    """ Standard run procedure for parallel simulations

    Arguments:
    parameters= must be a list of parameter dictionaries.  Each dictionary in the list is passed to a different worker as the arguments for Experiment.exp_num

    Use ExperimentsParallel.exp_num to design parameters list

    Arguments:
    parameters= must be a list of parameter dictionaries.  Each dictionary in the list is passed to a different worker as the arguments for Experiment.exp_num

    Use ExperimentsParallel.exp_num to design parameters list

    To use multiple workers, python script must be called from the interpreter using syntax:
    'mpiexec -n 10 python script.py'
    the call to mpiexec initializes the mpi with 10 workers

    _run_parallel should be called as:
    if __name__=="__main__":
        _run_parallel(ExperimentsParallel('exp_4a').parameters)
    """

    # make parallel context global
    global pc

    print parameters
    # create parallel context instance
    pc = h.ParallelContext()

    print 'i am', pc.id(), 'of', pc.nhost()
    # start workers, begins an infinitely loop where master workers posts jobs and workers pull jobs until all jobs are finished
    pc.runworker()
    
    print 'length of parameters:',len(parameters)
    # # # distribute experiment and parameters to workers
    for param in parameters:
        # print len(parameters)
        print 'submitted param:',param
        pc.submit(f_parallel, param)
        # print param

    # # continue runnning until all workers are finished
    while pc.working():
        print pc.id(), 'is working'

    # # close parallel context 
    pc.done()

def from_command_line_run_parallel(**kwargs):
    ''' run parallel simulations from command line by passing experiment name

    command line format for parallel: 'mpiexec -n 10 python run_cntrol.py experiment_name parallel'

    command line format for single run: 'python run_cntrol.py experiment_name run'

    experiment_name must point to a class with methods run_parallel and setup_parallel_parameters

    experiment names can be passed under the name==main statement with the kwarg 'experiment'.  this will be overwritten if an experiment name is passed at the command line
    '''
    # set default experiment
    default_experiment = 'exp_reduced_neuron_tbs'
    parallel =True
    if 'experiment' in kwargs:
        default_experiment = kwargs['experiment']
    if 'parallel' in kwargs:
        parallel= kwargs['parallel']
    
    # get command line arguments
    #----------------------------
    n_args = len(sys.argv) # number of arguments
    module_name = str(sys.argv[0]) # module name
    # arbitrary arguments
    if n_args>1:
        command_args = [str(arg) for arg in sys.argv[1:] ]
    else:
        command_args=[]
    # get any arguments that contain experiment names
    #------------------------------------------------
    experiments = [arg for arg in command_args if 'exp_' in arg]
    # if no experiments are passed at command line set to default
    if len(experiments)==0:
        experiments=[default_experiment]
    # iterate over experiments
    for experiment in experiments:
        # get globals
        _g = globals()
        # check that experiment exists
        if experiment in _g:
            # get experiment class from globals
            experiment_class = _g[experiment]
            # just run
            #-------------------------------------
            if 'run' in command_args:
                print 'running experiment', experiment, 'in series'
                experiment_class().run()
            # run parallel?
            #--------------------------------------
            elif 'parallel' in command_args or parallel:
                print 'running experiment', experiment, ' in parallel'
                _parameters = experiment_class().setup_parallel_parameters()
                run_parallel(parameters=_parameters)
            
#############################################################################
#############################################################################
class Experiment:
    """ Impliment experimental procedures.  
    """
    
    def __init__(self, **kwargs):
        """ choose experiment to run

        kwargs must be a dictionary with 'experiment' as a key, the corresponding value points to a given experiment, which is then fed the same kwarg dictionary and run
        """
        if not kwargs:
            pass
        else:
            # retrieve which experiment to run
            experiment = getattr(self, kwargs['experiment'])

            # run experiment
            experiment(**kwargs) 

        """
        EXPERIMENT 1
        Distance and synapse number dependence of DCS effects
        """
    
    def exp_branch_sequence(self, **kwargs):
        ''' simulate sequences of inputs on various dendritic branches as in Branco and Hausser 2010

        randomly choose dendritic branch, reduce the size of each segment to be ~10 um, activate segments in a sequence (either towards or away from the terminal), vary the following parameters
        ==Parameters==
        -field : 0, 1, 5, 10, 20 (positive and negative)
        -branches : first five branches in eahc tree
        -sequence delays : 1, 2, 5, 10, ms between activation of neighboring segments (segments are ~10 um apart)
        -synaptic weights : .001, .0015, .002
        '''
        # updates to global parameter dictionary 
        p_update = {
        'experiment' : inspect.stack()[0][3], 
        'trials' : 1,
        'field':kwargs['field'],
        'rec_variables':[('v','range','v'), ('i','syn','nmda')],
        'active_paths':['1',],
        # set active conductances to zero
        #---------------------------------
        # 'gna':0.,
        # 'ghd':0.,
        # 'gkdr':0., 
        # 'gcalbar': 0.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        # 'KMULT':0.,
        # 'KMULTP':0.,
        # 'Cm':1.,
        'branch_seg_L':10, 
        'dgna':-.0001,
        'AXONM':10
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apical_tuft','apical_trunk'],
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': [],
        'pulses': 1.,
        'pulse_freq': 100.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1.5*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        'delay':.1, 
        'sequence_direction':'in',
        'nmda_ampa_ratio':3.,
        'syn_limit':12
        },}

        # setup cell and updated parameter structures
        self.P, self.cell = self._standard_parameter_setup(
            default_p='migliore_2005',
            cell_class='CellMigliore2005',
            p_update=p_update,
            paths_update=paths_update,
            load_fd=True)

        # get sections with terminal branches
        terminal_branches = stims._get_terminal_branches(self.cell.geo)
        # maximum branches to simulate per tree
        max_branches = 1
        delays = [ 2, ]
        directions = ['in','out']
        weights = kwargs['weight']
        # iterate over branches
        for tree_key, branch_sec in terminal_branches.iteritems():
            if tree_key!='soma' and tree_key!='axon':
                for sec_i, sec_num in enumerate(branch_sec):
                    if sec_i<max_branches:
                        
                        # setup cell and updated parameter structures
                        self.P, self.cell = self._standard_parameter_setup(
                            default_p='migliore_2005',
                            cell_class='CellMigliore2005',
                            p_update=p_update,
                            paths_update=paths_update,
                            load_fd=True)
                        # section to simulate
                        sec_idx = {tree_key:[sec_num]}
                        self.P.paths['1']['sec_idx']=sec_idx
                        # iterate over delays
                        for delay in delays:
                            self.P.paths['1']['delay']=delay
                            # iterate over directions
                            for direction in directions:
                                self.P.paths['1']['sequence_direction']=direction
                                if direction=='in':
                                    reverse=True
                                elif direction=='out':
                                    reverse=False
                                # iterate over trials
                                for weight in weights:
                                    self.P.paths['1']['w_mean']=weight
                                    for trial_i, trial in enumerate(range(self.P.p['trials'])):
                                        self.P = self._update_synapse_parameters_sequence(P=self.P, cell=self.cell, method='_choose_seg_rand', reverse=reverse, sec_idx=sec_idx)
                                        # measure distance of each segment from the soma and store in parameter dictionary
                                        self.P.p['seg_dist'] = self.P._seg_distance(self.cell)

                                        # create morphology for shape plots
                                        self.P.p['morpho'] = self.P._create_morpho(self.cell.geo)

                                        self.run_obj = self._standard_run_and_save_df(P=self.P, cell=self.cell, trial=trial)
    
    def exp_branch_sequence_milstein(self, **kwargs):
        ''' simulate sequences of inputs on various dendritic branches as in Branco and Hausser 2010

        randomly choose dendritic branch, reduce the size of each segment to be ~10 um, activate segments in a sequence (either towards or away from the terminal), vary the following parameters
        ==Parameters==
        -field : 0, 1, 5, 10, 20 (positive and negative)
        -branches : first five branches in eahc tree
        -sequence delays : 1, 2, 5, 10, ms between activation of neighboring segments (segments are ~10 um apart)
        -synaptic weights : .001, .0015, .002
        '''
        # updates to global parameter dictionary 
        p_update = {
        'experiment' : inspect.stack()[0][3], 
        'trials' : 1,
        'field':kwargs['field'],
        'rec_variables':[('v','range','v')],
        'active_paths':['1',],
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apical_tuft','apical_trunk'],
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': [],
        'pulses': 1.,
        'pulse_freq': 100.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1.5*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        'delay':.1, 
        'sequence_direction':'in',
        'nmda_ampa_ratio':3.,
        'syn_limit':12
        },}

        # load cell
        #-------------------------------------------
        morph_filename = 'EB1-early-bifurcation.swc'
        mech_filename='043016 Type A - km2_NMDA_KIN5_Pr'
        mech_dict = mech_dicts.ca1_milstein_1
        cell = spec.CA1_Pyr(morph_filename=morph_filename, mech_filename=mech_filename, mech_dict=mech_dict, full_spines=False)
        # get branch and set synapse locations
        #--------------------------------------------------
        branches = ca1.get_terminal_branches(sec_types=['basal'], max_distance=200)
        # select branch, adjust nseg, reinit mechanisms
        #-------------------------------------------------
        branch=branches[-1]
        seg_L=10.
        ca1.set_nseg(branch, seg_L=seg_L)
        ca1.update_diameters()
        ca1.reinit_mechanisms()
        # list of locations along section to stimulate
        #-------------------------------------------------
        nseg = branch.sec.nseg
        n_locs=8
        locs = np.arange(1, 0, -seg_L/branch.sec.L)[:n_locs]
        # locs = np.arange(0, 1, seg_L/branch.sec.L)[:n_locs]
        locs = np.flip(locs)
        syn_locs=[]
        delays=[]
        current_delay=0
        sequence_delay=6
        for loc in locs:
            delays.append(current_delay)
            current_delay+=sequence_delay
            syn_locs.append((branch.name, loc))
        # update parameter dictionaries
        self.P = self.update_parameters(p_update=p_update, paths_update=paths_update)
        self.P.paths['1']['syn_idx']=syn_locs
        self.P.paths['1']['sequence_delays']=delays
        # create stimulators
        stim = stims.Uncage()._bursts(p_path)
        dcs=stims.DCS(cell=ca1, cell_type='milstein', field_components=[0,-40,0], field_on=5, field_off=100)
        # insert synapses
        #----------------------------------------------------
        # apical trunk needs to be insert because parameters are inherited from here
        for node in ca1.trunk:
            spec.Synapse(cell=ca1, node=node, type_list=syn_types,)
        # insert synapses in the nodes specified by syn_locs
        for node_i, node_loc in enumerate(syn_locs): 
            node = node_loc[0]
            loc = node_loc[1]
            syns.append(spec.Synapse(cell=ca1, node=node, type_list=syn_types, source=stim[node_i][0], stochastic=0, loc=loc, weight=8.))

        #FIXME setup recording vectors, run, and save
        # create standard method when cell is milstein type

    def exp_1a1(self, **kwargs):
        """ Activate a single pathway with a single TBS burst at varying distance from soma: 0-300, 100-400, 200-500, 300-600, 0-200, 100-300, 200-400, 300-500, 400-600

        6, 8, 10, 12, 16,20 synapses are randomly distributed within these distance ranges with conductance 0.001 uS/cm2.  This is chosen to be close to threshold by empirical testing
        
        ==Args==
        -experiment : string containing the name of the experiment, eg 'exp_1a1'
        -trials :  integer number of trials to run
        -syn_dist :  list, or nested list, containing synapse distance requirements for choose_seg_rand. eg [0, 200] chooses synapses between 0 and 200 um from soma
        -syn_num :  list of number of synapses to be activated. note that synapses are chosen randomly with replacement, so the same synapse can be selected multiple times.  In this case the weight is multiplied by the number of times the synapse is selected
        """

        # updates to global parameter dictionary 
        p_update = {
        'experiment' : kwargs['experiment'], 
        'trials' : kwargs['trials'],
        'field':[-20.,0.,20.],
        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),('ica_calH','range','calH'), ('i','syn','nmda')],
        'active_paths':['1',],
        'gcalbar': 1.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apical_tuft','apical_trunk'],
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': kwargs['syn_dist'],
        'pulses': 4.,
        'pulse_freq': 100.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        },}

        # setup cell and updated parameter structures
        self.P, self.cell = self._standard_parameter_setup(
            default_p='migliore_2005',
            cell_class='CellMigliore2005',
            p_update=p_update,
            paths_update=paths_update,
            load_fd=True)

        # iterate over synapse number 
        for syn_num in kwargs['syn_num']:
            # update parameter dictionary
            self.P.paths['1']['syn_num']=syn_num

            # iterate over trials
            for trial_i, trial in enumerate(range(self.P.p['trials'])):

                # create list of active synapses, weights, delays
                # stored in P.p['seg_idx', 'w_list', 'sequence_delays'], 
                self.P = self._update_synapse_parameters(P=self.P, cell=self.cell, method='_choose_seg_rand')

                self.run_obj = self._standard_run_and_save(P=self.P, cell=self.cell, trial=trial)

    def exp_1a1_4compartment(self, **kwargs):
        """ Activate a single pathway with a single TBS burst at varying distance from soma: 0-300, 100-400, 200-500, 300-600, 0-200, 100-300, 200-400, 300-500, 400-600

        6, 8, 10, 12, 16,20 synapses are randomly distributed within these distance ranges with conductance 0.001 uS/cm2.  This is chosen to be close to threshold by empirical testing
        
        ==Args==
        -experiment : string containing the name of the experiment, eg 'exp_1a1'
        -trials :  integer number of trials to run
        -syn_dist :  list, or nested list, containing synapse distance requirements for choose_seg_rand. eg [0, 200] chooses synapses between 0 and 200 um from soma
        -syn_num :  list of number of synapses to be activated. note that synapses are chosen randomly with replacement, so the same synapse can be selected multiple times.  In this case the weight is multiplied by the number of times the synapse is selected
        """

        # updates to global parameter dictionary 
        p_update = {
        'experiment' : inspect.stack()[0][3],#kwargs['experiment'], 
        'trials' : kwargs['trials'],
        'field':[0.],
        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),('ica_calH','range','calH'), ('i','syn','nmda')],
        'active_paths':['1','2','3'],
        'gcalbar': 2.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        # 'gkdr':.1,
        }

        # set up synaptic pathway parameters
        paths_update = {
        '1':{
        'trees': ['apical_prox',],
        'syn_num': 2,
        'nsyns': 1.,
        'syn_dist': [0,10000],
        'pulses': 4.,
        'pulse_freq': 100.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 0*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        },
        '2':{
        'trees': ['apical_dist',],
        'syn_num': 3,
        'nsyns': 1.,
        'syn_dist': [0,10000],
        'pulses': 4.,
        'pulse_freq': 100.,
        'bursts': 2.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        },
        '3':{
        'trees': ['basal',],
        'syn_num': 4,
        'nsyns': 1.,
        'syn_dist': [0,10000],
        'pulses': 4.,
        'pulse_freq': 100.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        }}

        # setup cell and updated parameter structures
        self.P, self.cell = self._standard_parameter_setup(
            default_p='pyramidal_cylinder',
            cell_class='PyramidalCylinder',
            p_update=p_update,
            paths_update=paths_update,
            load_fd=True)

        # update parameter dictionary
        # self.P.paths['1']['syn_num']=syn_num

        # iterate over trials
        for trial_i, trial in enumerate(range(self.P.p['trials'])):

            # create list of active synapses, weights, delays
            # stored in P.p['seg_idx', 'w_list', 'sequence_delays'], 
            self.P = self._update_synapse_parameters(P=self.P, cell=self.cell, method='_choose_seg_rand')

            self.run_obj = self._standard_run_and_save_df(P=self.P, cell=self.cell, trial=trial)

    def exp_1a1_nablock(self, **kwargs):
        """ repeat experiment 1a1 with sodium channels block in soma and axon

        Activate a single pathway with a single TBS burst at varying distance from soma: 0-300, 100-400, 200-500, 300-600, 0-200, 100-300, 200-400, 300-500, 400-600

        6, 8, 10, 12, 16,20 synapses are randomly distributed within these distance ranges with conductance 0.0015 uS/cm2.  This is chosen to be close to threshold by empirical testing
        
        ==Args==
        -experiment : string containing the name of the experiment, eg 'exp_1a1'
        -trials :  integer number of trials to run
        -syn_dist :  list, or nested list, containing synapse distance requirements for choose_seg_rand. eg [0, 200] chooses synapses between 0 and 200 um from soma
        -syn_num :  list of number of synapses to be activated. note that synapses are chosen randomly with replacement, so the same synapse can be selected multiple times.  In this case the weight is multiplied by the number of times the synapse is selected
        """

        # updates to global parameter dictionary 
        p_update = {
        'experiment' : kwargs['experiment'], 
        'trials' : kwargs['trials'],
        'field':[-20.,0.,20.],
        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),('ica_calH','range','calH'), ('i','syn','nmda')],
        'active_paths':['1',],
        'gcalbar': 1.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apical_tuft','apical_trunk'],
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': kwargs['syn_dist'],
        'pulses': 4.,
        'pulse_freq': 100.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        },}

        # setup cell and updated parameter structures
        self.P, self.cell = self._standard_parameter_setup(
            default_p='migliore_2005',
            cell_class='CellMigliore2005',
            p_update=p_update,
            paths_update=paths_update,
            load_fd=True)

        # iterate over synapse number 
        for syn_num in kwargs['syn_num']:
            # update parameter dictionary
            self.P.paths['1']['syn_num']=syn_num

            # iterate over trials
            for trial_i, trial in enumerate(range(self.P.p['trials'])):

                # create list of active synapses, weights, delays
                # stored in P.p['seg_idx', 'w_list', 'sequence_delays'], 
                self.P = self._update_synapse_parameters(P=self.P, cell=self.cell, method='_choose_seg_rand')

                # set sodium conductnace to zero in soma and axon
                for tree_key, tree in self.cell.geo.iteritems():
                    if tree_key=='soma':
                        for sec_i, sec in enumerate(tree):
                            sec.gbar_na3=0.
                    if tree_key=='axon':
                        for sec_i, sec in enumerate(tree):
                            sec.gbar_nax=0.       

                self.run_obj = self._standard_run_and_save(P=self.P, cell=self.cell, trial=trial)

    def exp_1a_dose(self, **kwargs):
        ''' repeat experiment 1a1 with varying electric field intensity

        6,8,10,12,14 synapses at .001 uS per synapse, distributed evenly between 0 and 300 um from soma
        '''
        # updates to global parameter dictionary 
        p_update = {
        'experiment' : kwargs['experiment'], 
        'trials' : kwargs['trials'],
        'field':[-20.,-15, -10, -5.,-2, -1., -0.5, 0., 0.5, 1.,2., 5.,10., 15., 20.],
        # 'field':[-20.,-5., -1., -0.5, 0., 0.5, 1., 5., 20.],
        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),],
        'active_paths':['1',],
        'gcalbar': 1.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apical_tuft','apical_trunk'],
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': kwargs['syn_dist'],
        'pulses': 4.,
        'pulse_freq': 100.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        },}

        # setup cell and updated parameter structures
        self.P, self.cell = self._standard_parameter_setup(
            default_p='migliore_2005',
            cell_class='CellMigliore2005',
            p_update=p_update,
            paths_update=paths_update,
            load_fd=True)

        # iterate over synapse number 
        for syn_num in kwargs['syn_num']:
            # update parameter dictionary
            self.P.paths['1']['syn_num']=syn_num

            # iterate over trials
            for trial_i, trial in enumerate(range(self.P.p['trials'])):

                # create list of active synapses, weights, delays
                # stored in P.p['seg_idx', 'w_list', 'sequence_delays'], 
                self.P = self._update_synapse_parameters(P=self.P, cell=self.cell, method='_choose_seg_rand')

                self.run_obj = self._standard_run_and_save(P=self.P, cell=self.cell, trial=trial)

    def exp_1a_basal(self, **kwargs):
        ''' repeat experiment 1a1 with varying electric field intensity

        6,8,10,12,14 synapses at .001 uS per synapse, distributed evenly between 0 and 300 um from soma
        '''
        # updates to global parameter dictionary 
        p_update = {
        'experiment' : kwargs['experiment'], 
        'trials' : kwargs['trials'],
        'field':[-20., 0., 20.],
        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),],
        'active_paths':['1',],
        'gcalbar': 1.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['basal'],
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': kwargs['syn_dist'],
        'pulses': 4.,
        'pulse_freq': 100.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        },}

        # setup cell and updated parameter structures
        self.P, self.cell = self._standard_parameter_setup(
            default_p='migliore_2005',
            cell_class='CellMigliore2005',
            p_update=p_update,
            paths_update=paths_update,
            load_fd=True)

        # iterate over synapse number 
        for syn_num in kwargs['syn_num']:
            # update parameter dictionary
            self.P.paths['1']['syn_num']=syn_num

            # iterate over trials
            for trial_i, trial in enumerate(range(self.P.p['trials'])):

                # create list of active synapses, weights, delays
                # stored in P.p['seg_idx', 'w_list', 'sequence_delays'], 
                self.P = self._update_synapse_parameters(P=self.P, cell=self.cell, method='_choose_seg_rand')
                # increase simulation duration
                self.P.p['tstop']=100
                self.P.p['field_off']=self.P.p['tstop']

                self.run_obj = self._standard_run_and_save(P=self.P, cell=self.cell, trial=trial)

    def exp_1a_apical(self, **kwargs):
        ''' repeat experiment 1a1 with varying electric field intensity

        6,8,10,12,14 synapses at .001 uS per synapse, distributed evenly between 0 and 300 um from soma
        '''
        # updates to global parameter dictionary 
        p_update = {
        'experiment' : kwargs['experiment'], 
        'trials' : kwargs['trials'],
        'field':[-20., 0., 20.],
        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),],
        'active_paths':['1',],
        'gcalbar': 1.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apical_trunk', 'apical_tuft'],
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': kwargs['syn_dist'],
        'pulses': 4.,
        'pulse_freq': 100.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        },}

        # setup cell and updated parameter structures
        self.P, self.cell = self._standard_parameter_setup(
            default_p='migliore_2005',
            cell_class='CellMigliore2005',
            p_update=p_update,
            paths_update=paths_update,
            load_fd=True)

        # iterate over synapse number 
        for syn_num in kwargs['syn_num']:
            # update parameter dictionary
            self.P.paths['1']['syn_num']=syn_num

            # iterate over trials
            for trial_i, trial in enumerate(range(self.P.p['trials'])):

                # create list of active synapses, weights, delays
                # stored in P.p['seg_idx', 'w_list', 'sequence_delays'], 
                self.P = self._update_synapse_parameters(P=self.P, cell=self.cell, method='_choose_seg_rand')
                # increase simulation duration
                self.P.p['tstop']=100
                self.P.p['field_off']=self.P.p['tstop']

                self.run_obj = self._standard_run_and_save(P=self.P, cell=self.cell, trial=trial)
    
    def exp_1a_dose_basal(self, **kwargs):
        ''' repeat experiment 1a1 with varying electric field intensity

        6,8,10,12,14 synapses at .001 uS per synapse, distributed evenly between 0 and 300 um from soma
        '''
        # updates to global parameter dictionary 
        p_update = {
        'experiment' : kwargs['experiment'], 
        'trials' : kwargs['trials'],
        'field':[-20.,-15, -10, -5.,-2, -1., -0.5, 0., 0.5, 1.,2., 5.,10., 15., 20.],
        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),],
        'active_paths':['1',],
        'gcalbar': 1.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['basal'],
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': kwargs['syn_dist'],
        'pulses': 4.,
        'pulse_freq': 100.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        },}

        # setup cell and updated parameter structures
        self.P, self.cell = self._standard_parameter_setup(
            default_p='migliore_2005',
            cell_class='CellMigliore2005',
            p_update=p_update,
            paths_update=paths_update,
            load_fd=True)

        # iterate over synapse number 
        for syn_num in kwargs['syn_num']:
            # update parameter dictionary
            self.P.paths['1']['syn_num']=syn_num

            # iterate over trials
            for trial_i, trial in enumerate(range(self.P.p['trials'])):

                # create list of active synapses, weights, delays
                # stored in P.p['seg_idx', 'w_list', 'sequence_delays'], 
                self.P = self._update_synapse_parameters(P=self.P, cell=self.cell, method='_choose_seg_rand')
                # increase simulation duration
                self.P.p['tstop']=100

                self.run_obj = self._standard_run_and_save(P=self.P, cell=self.cell, trial=trial)
    
    def exp_1a_polarization_shapeplot(self, **kwargs):
        ''' repeat experiment 1a1 with varying electric field intensity

        6,8,10,12,14 synapses at .001 uS per synapse, distributed evenly between 0 and 300 um from soma
        '''
        # updates to global parameter dictionary 
        p_update = {
        'experiment' : kwargs['experiment'], 
        'trials' : 1,
        'field':[-20., 0.,20.],
        'rec_variables':[('v','range','v'),],
        'active_paths':[],
        'gcalbar': 1.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apical_tuft','apical_trunk'],
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': [0,300],
        'pulses': 4.,
        'pulse_freq': 100.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        },}

        # setup cell and updated parameter structures
        self.P, self.cell = self._standard_parameter_setup(
            default_p='migliore_2005',
            cell_class='CellMigliore2005',
            p_update=p_update,
            paths_update=paths_update,
            load_fd=True)


        # create list of active synapses, weights, delays
        # stored in P.p['seg_idx', 'w_list', 'sequence_delays'], 
        self.P = self._update_synapse_parameters(P=self.P, cell=self.cell, method='_choose_seg_rand')

        self.P.p['rec_idx'] = self.P._create_loc_list(geo=self.cell.geo)

        self.run_obj = self._standard_run_and_save(P=self.P, cell=self.cell, trial=0)

    def exp_1a_ACS(self, **kwargs):
        """ Activate a single pathway with a single TBS burst at varying distance from soma: 0-300, 100-400, 200-500, 300-600, 0-200, 100-300, 200-400, 300-500, 400-600

        6, 8, 10, 12, 16,20 synapses are randomly distributed within these distance ranges with conductance 0.001 uS/cm2.  This is chosen to be close to threshold by empirical testing
        
        ==Args==
        -experiment : string containing the name of the experiment, eg 'exp_1a1'
        -trials :  integer number of trials to run
        -syn_dist :  list, or nested list, containing synapse distance requirements for choose_seg_rand. eg [0, 200] chooses synapses between 0 and 200 um from soma
        -syn_num :  list of number of synapses to be activated. note that synapses are chosen randomly with replacement, so the same synapse can be selected multiple times.  In this case the weight is multiplied by the number of times the synapse is selected
        """

        # updates to global parameter dictionary 
        p_update = {
        'experiment' : kwargs['experiment'], 
        'trials' : kwargs['trials'],
        'field':[0.],
        'ac_field':20.,
        'ac_field_freq':5.,
        'ac_field_phase':np.pi/2,
        'ac_field_on':10.,
        'ac_field_off':300.,
        'ac_field_angle':0.,

        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),],
        'active_paths':['1',],
        'gcalbar': 1.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apical_tuft','apical_trunk'],
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': kwargs['syn_dist'],
        'pulses': 4.,
        'pulse_freq': 100.,
        'bursts': 2.,
        'burst_freq': 5.,
        'warmup': 50,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        },}

        # setup cell and updated parameter structures
        self.P, self.cell = self._standard_parameter_setup(
            default_p='migliore_2005',
            cell_class='CellMigliore2005',
            p_update=p_update,
            paths_update=paths_update,
            load_fd=True)
        
        phases = [0, np.pi, None]
        # iterate over synapse number 
        for syn_num in kwargs['syn_num']:
            # update parameter dictionary
            self.P.paths['1']['syn_num']=syn_num

            # iterate over trials
            for trial_i, trial in enumerate(range(self.P.p['trials'])):
                
                # create list of active synapses, weights, delays
                # stored in P.p['seg_idx', 'w_list', 'sequence_delays'], 
                self.P = self._update_synapse_parameters(P=self.P, cell=self.cell, method='_choose_seg_rand')

                # iterate over peak and trough phases
                for phase in phases:
                    if phase is not None:
                        self.P.p['ac_field_phase']=phase
                        self.P.p['ac_field']=p_update['ac_field']
                    else:
                        self.P.p['ac_field']=0
                        self.P.p['ac_field_phase']=0

                    self.run_obj = self._standard_run_and_save(P=self.P, cell=self.cell, trial=trial,)

    def exp_1a_apical_full_duration(self, **kwargs):
        ''' repeat experiment 1a1 with varying electric field intensity

        6,8,10,12,14 synapses at .001 uS per synapse, distributed evenly between 0 and 300 um from soma
        '''
        # updates to global parameter dictionary 
        p_update = {
        'experiment' : kwargs['experiment'], 
        'trials' : kwargs['trials'],
        'field':[-20., 0., 20.],
        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),],
        'active_paths':['1',],
        'gcalbar': 1.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apical_trunk', 'apical_tuft'],
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': kwargs['syn_dist'],
        'pulses': 4.,
        'pulse_freq': 100.,
        'bursts': 15.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        },}

        # setup cell and updated parameter structures
        self.P, self.cell = self._standard_parameter_setup(
            default_p='migliore_2005',
            cell_class='CellMigliore2005',
            p_update=p_update,
            paths_update=paths_update,
            load_fd=True)

        # iterate over synapse number 
        for syn_num in kwargs['syn_num']:
            # update parameter dictionary
            self.P.paths['1']['syn_num']=syn_num

            # iterate over trials
            for trial_i, trial in enumerate(range(self.P.p['trials'])):

                # create list of active synapses, weights, delays
                # stored in P.p['seg_idx', 'w_list', 'sequence_delays'], 
                self.P = self._update_synapse_parameters(P=self.P, cell=self.cell, method='_choose_seg_rand')
                # increase simulation duration
                # self.P.p['tstop']=100
                self.P.p['field_off']=self.P.p['tstop']

                self.run_obj = self._standard_run_and_save(P=self.P, cell=self.cell, trial=trial)

    def exp_1a_apical_full_duration_20Hz(self, **kwargs):
        ''' repeat experiment 1a1 with varying electric field intensity

        6,8,10,12,14 synapses at .001 uS per synapse, distributed evenly between 0 and 300 um from soma
        '''
        # updates to global parameter dictionary 
        p_update = {
        'experiment' : kwargs['experiment'], 
        'trials' : kwargs['trials'],
        'field':[-20., 0., 20.],
        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),],
        'active_paths':['1',],
        'gcalbar': 1.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apical_trunk', 'apical_tuft'],
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': kwargs['syn_dist'],
        'pulses': 60.,
        'pulse_freq': 20.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        },}

        # setup cell and updated parameter structures
        self.P, self.cell = self._standard_parameter_setup(
            default_p='migliore_2005',
            cell_class='CellMigliore2005',
            p_update=p_update,
            paths_update=paths_update,
            load_fd=True)

        # iterate over synapse number 
        for syn_num in kwargs['syn_num']:
            # update parameter dictionary
            self.P.paths['1']['syn_num']=syn_num

            # iterate over trials
            for trial_i, trial in enumerate(range(self.P.p['trials'])):

                # create list of active synapses, weights, delays
                # stored in P.p['seg_idx', 'w_list', 'sequence_delays'], 
                self.P = self._update_synapse_parameters(P=self.P, cell=self.cell, method='_choose_seg_rand')
                # increase simulation duration
                # self.P.p['tstop']=100
                self.P.p['field_off']=self.P.p['tstop']

                self.run_obj = self._standard_run_and_save(P=self.P, cell=self.cell, trial=trial)

    def exp_1a_basal_full_duration(self, **kwargs):
        ''' repeat experiment 1a1 with varying electric field intensity

        6,8,10,12,14 synapses at .001 uS per synapse, distributed evenly between 0 and 300 um from soma
        '''
        # updates to global parameter dictionary 
        p_update = {
        'experiment' : kwargs['experiment'], 
        'trials' : kwargs['trials'],
        'field':[-20., 0., 20.],
        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),],
        'active_paths':['1',],
        'gcalbar': 1.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['basal'],
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': kwargs['syn_dist'],
        'pulses': 4.,
        'pulse_freq': 100.,
        'bursts': 15.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        },}

        # setup cell and updated parameter structures
        self.P, self.cell = self._standard_parameter_setup(
            default_p='migliore_2005',
            cell_class='CellMigliore2005',
            p_update=p_update,
            paths_update=paths_update,
            load_fd=True)

        # iterate over synapse number 
        for syn_num in kwargs['syn_num']:
            # update parameter dictionary
            self.P.paths['1']['syn_num']=syn_num

            # iterate over trials
            for trial_i, trial in enumerate(range(self.P.p['trials'])):

                # create list of active synapses, weights, delays
                # stored in P.p['seg_idx', 'w_list', 'sequence_delays'], 
                self.P = self._update_synapse_parameters(P=self.P, cell=self.cell, method='_choose_seg_rand')
                # increase simulation duration
                # self.P.p['tstop']=100
                self.P.p['field_off']=self.P.p['tstop']

                self.run_obj = self._standard_run_and_save(P=self.P, cell=self.cell, trial=trial)

    def exp_1a_basal_full_duration_20Hz(self, **kwargs):
        ''' repeat experiment 1a1 with varying electric field intensity

        6,8,10,12,14 synapses at .001 uS per synapse, distributed evenly between 0 and 300 um from soma
        '''
        # updates to global parameter dictionary 
        p_update = {
        'experiment' : kwargs['experiment'], 
        'trials' : kwargs['trials'],
        'field':[-20., 0., 20.],
        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),],
        'active_paths':['1',],
        'gcalbar': 1.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['basal'],
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': kwargs['syn_dist'],
        'pulses': 60.,
        'pulse_freq': 20.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        },}

        # setup cell and updated parameter structures
        self.P, self.cell = self._standard_parameter_setup(
            default_p='migliore_2005',
            cell_class='CellMigliore2005',
            p_update=p_update,
            paths_update=paths_update,
            load_fd=True)

        # iterate over synapse number 
        for syn_num in kwargs['syn_num']:
            # update parameter dictionary
            self.P.paths['1']['syn_num']=syn_num

            # iterate over trials
            for trial_i, trial in enumerate(range(self.P.p['trials'])):

                # create list of active synapses, weights, delays
                # stored in P.p['seg_idx', 'w_list', 'sequence_delays'], 
                self.P = self._update_synapse_parameters(P=self.P, cell=self.cell, method='_choose_seg_rand')
                # increase simulation duration
                # self.P.p['tstop']=100
                self.P.p['field_off']=self.P.p['tstop']

                self.run_obj = self._standard_run_and_save(P=self.P, cell=self.cell, trial=trial)

    def exp_1a_poisson(self, **kwargs):
        ''' repeat experiment 1a1 with varying electric field intensity

        6,8,10,12,14 synapses at .001 uS per synapse, distributed evenly between 0 and 300 um from soma
        '''
        # updates to global parameter dictionary 
        p_update = {
        'experiment' : kwargs['experiment'], 
        'trials' : kwargs['trials'],
        'field':[-20., 0., 20.],
        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),],
        'active_paths':['1',],
        'gcalbar': 1.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apical_trunk', 'apical_tuft'],
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': kwargs['syn_dist'],
        'pulses': 10.,
        'pulse_freq': 20.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0.9,
        'mod_freq':5,
        'mod_amp':38,
        'mean_rate':40
        },}

        # setup cell and updated parameter structures
        self.P, self.cell = self._standard_parameter_setup(
            default_p='migliore_2005',
            cell_class='CellMigliore2005',
            p_update=p_update,
            paths_update=paths_update,
            load_fd=True)

        # iterate over synapse number 
        for syn_num in kwargs['syn_num']:
            # update parameter dictionary
            self.P.paths['1']['syn_num']=syn_num

            # iterate over trials
            for trial_i, trial in enumerate(range(self.P.p['trials'])):

                # create list of active synapses, weights, delays
                # stored in P.p['seg_idx', 'w_list', 'sequence_delays'], 
                self.P = self._update_synapse_parameters(P=self.P, cell=self.cell, method='_choose_seg_rand')
                # increase simulation duration
                self.P.p['tstop']=3000
                self.P.p['field_off']=self.P.p['tstop']
                self.P.p['p_path']['1']['dt']=self.P.p['dt']
                self.P.p['p_path']['1']['tstop']=self.P.p['tstop']


                run_kwargs={'uncage_method':'_poisson'}
                self.run_obj = self._standard_run_and_save(P=self.P, cell=self.cell, trial=trial, **run_kwargs)
        
        '''
        EXPERIMENT 4 
        two pathway associativity and specificity
        '''
    
    def exp_4a1(self, **kwargs):
        """ two pathway experiments

        both pathways are distributed 0-300 um from soma.  pathway 1 receives a theta burst, path 2 receives a single pulse timed to the second pulse of the theta burst

        4,5,6,8,10,12,16,20 synapses at 0-200 or 0-300 um from soma
        
        ==Args==
        -experiment : string containing the name of the experiment, eg 'exp_1a1'
        -trials :  integer number of trials to run
        -syn_dist :  list, or nested list, containing synapse distance requirements for choose_seg_rand. eg [0, 200] chooses synapses between 0 and 200 um from soma
        -syn_num :  list of number of synapses to be activated. note that synapses are chosen randomly with replacement, so the same synapse can be selected multiple times.  In this case the weight is multiplied by the number of times the synapse is selected
        """

        # updates to global parameter dictionary 
        p_update = {
        'experiment' : kwargs['experiment'], 
        'trials' : kwargs['trials'],
        'field':[-20.,0.,20.],
        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),('ica_calH','range','calH'), ('i','syn','nmda')],
        'active_paths':['1','2'],
        'gcalbar': 1.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apical_tuft','apical_trunk'],
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': [],#kwargs['syn_dist'],
        'pulses': 4.,
        'pulse_freq': 100.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        },
        '2':{
        'trees': ['apical_tuft','apical_trunk'],
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': [],#kwargs['syn_dist'],
        'pulses': 1.,
        'pulse_freq': 100.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 20,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        }}

        # setup cell and updated parameter structures
        self.P, self.cell = self._standard_parameter_setup(
            default_p='migliore_2005',
            cell_class='CellMigliore2005',
            p_update=p_update,
            paths_update=paths_update,
            load_fd=True)

        for syn_dist in kwargs['syn_dist']:
            # iterate over synapse number 
            for syn_num in kwargs['syn_num']:
                # update parameter dictionary
                self.P.paths['1']['syn_num']=syn_num
                self.P.paths['2']['syn_num']=syn_num
                # update parameter dictionary
                self.P.paths['1']['syn_dist']=syn_dist
                self.P.paths['2']['syn_dist']=syn_dist

                # iterate over trials
                for trial_i, trial in enumerate(range(self.P.p['trials'])):

                    # create list of active synapses, weights, delays
                    # stored in P.p['seg_idx', 'w_list', 'sequence_delays'], 
                    self.P = self._update_synapse_parameters(P=self.P, cell=self.cell, method='_choose_seg_rand')

                    self.run_obj = self._standard_run_and_save(P=self.P, cell=self.cell, trial=trial)

    def exp_4a2(self, **kwargs):
        """ weak pathway only

        single pulse applied to compare with experiment 4a1
        
        ==Args==
        -experiment : string containing the name of the experiment, eg 'exp_1a1'
        -trials :  integer number of trials to run
        -syn_dist :  list, or nested list, containing synapse distance requirements for choose_seg_rand. eg [0, 200] chooses synapses between 0 and 200 um from soma
        -syn_num :  list of number of synapses to be activated. note that synapses are chosen randomly with replacement, so the same synapse can be selected multiple times.  In this case the weight is multiplied by the number of times the synapse is selected
        """

        # updates to global parameter dictionary 
        p_update = {
        'experiment' : kwargs['experiment'], 
        'trials' : kwargs['trials'],
        'field':[-20.,0.,20.],
        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),('ica_calH','range','calH'), ('i','syn','nmda')],
        'active_paths':['2'],
        'gcalbar': 1.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        }

        # set up synaptic pathway parameters
        paths_update = {
        '2':{
        'trees': ['apical_tuft','apical_trunk'],
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': [],#kwargs['syn_dist'],
        'pulses': 1.,
        'pulse_freq': 100.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 20,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        }}

        # setup cell and updated parameter structures
        self.P, self.cell = self._standard_parameter_setup(
            default_p='migliore_2005',
            cell_class='CellMigliore2005',
            p_update=p_update,
            paths_update=paths_update,
            load_fd=True)

        for syn_dist in kwargs['syn_dist']:
            # iterate over synapse number 
            for syn_num in kwargs['syn_num']:
                # update parameter dictionary
                # self.P.paths['1']['syn_num']=syn_num
                self.P.paths['2']['syn_num']=syn_num
                # update parameter dictionary
                # self.P.paths['1']['syn_dist']=syn_dist
                self.P.paths['2']['syn_dist']=syn_dist

                # iterate over trials
                for trial_i, trial in enumerate(range(self.P.p['trials'])):

                    # create list of active synapses, weights, delays
                    # stored in P.p['seg_idx', 'w_list', 'sequence_delays'], 
                    self.P = self._update_synapse_parameters(P=self.P, cell=self.cell, method='_choose_seg_rand')

                    self.run_obj = self._standard_run_and_save(P=self.P, cell=self.cell, trial=trial)

    def exp_4a1_poisson(self, **kwargs):
        ''' repeat experiment 1a1 with varying electric field intensity

        6,8,10,12,14 synapses at .001 uS per synapse, distributed evenly between 0 and 300 um from soma
        '''
        # updates to global parameter dictionary 
        p_update = {
        'experiment' : kwargs['experiment'], 
        'trials' : kwargs['trials'],
        'field':[20.,-1.,  0.,1.,20.],
        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),],
        'active_paths':['1','2'],
        'gcalbar': 1.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apical_trunk', 'apical_tuft'],
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': kwargs['syn_dist'],
        'pulses': 10.,
        'pulse_freq': 20.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1.5*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':98./100.,
        'mod_freq':5,
        'mod_amp':0,
        'mean_rate':10
        },
        '2':{
        'trees': ['apical_trunk', 'apical_tuft'],
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': kwargs['syn_dist'],
        'pulses': 10.,
        'pulse_freq': 20.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1.5*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':198./200.,
        'mod_freq':5,
        'mod_amp':0,
        'mean_rate':2
        },}

        # setup cell and updated parameter structures
        self.P, self.cell = self._standard_parameter_setup(
            default_p='migliore_2005',
            cell_class='CellMigliore2005',
            p_update=p_update,
            paths_update=paths_update,
            load_fd=True)

        self.P.paths['1']['mean_rate']=kwargs['path_1_mean_rate']
        self.P.paths['2']['mean_rate']=kwargs['path_2_mean_rate']
        # iterate over synapse number 
        for syn_num in kwargs['syn_num']:
            # update parameter dictionary
            self.P.paths['1']['syn_num']=20
            self.P.paths['2']['syn_num']=50
            

            # iterate over trials
            for trial_i, trial in enumerate(range(self.P.p['trials'])):

                # create list of active synapses, weights, delays
                # stored in P.p['seg_idx', 'w_list', 'sequence_delays'], 
                self.P = self._update_synapse_parameters(P=self.P, cell=self.cell, method='_choose_seg_rand')
                # increase simulation duration
                self.P.p['tstop']=1000
                self.P.p['field_off']=self.P.p['tstop']
                self.P.p['p_path']['1']['dt']=self.P.p['dt']
                self.P.p['p_path']['1']['tstop']=self.P.p['tstop']
                self.P.p['p_path']['2']['dt']=self.P.p['dt']
                self.P.p['p_path']['2']['tstop']=self.P.p['tstop']


                run_kwargs={'uncage_method':'_poisson'}
                self.run_obj = self._standard_run_and_save(P=self.P, cell=self.cell, trial=trial, **run_kwargs)
        
        '''
        EXPERIMENT 4 
        two pathway associativity and specificity
        '''

    def exp_4a1_nablock(self, **kwargs):
        """ two pathway experiments

        both pathways are distributed 0-300 um from soma.  pathway 1 receives a theta burst, path 2 receives a single pulse timed to the second pulse of the theta burst

        4,5,6,8,10,12,16,20 synapses at 0-200 or 0-300 um from soma
        
        ==Args==
        -experiment : string containing the name of the experiment, eg 'exp_1a1'
        -trials :  integer number of trials to run
        -syn_dist :  list, or nested list, containing synapse distance requirements for choose_seg_rand. eg [0, 200] chooses synapses between 0 and 200 um from soma
        -syn_num :  list of number of synapses to be activated. note that synapses are chosen randomly with replacement, so the same synapse can be selected multiple times.  In this case the weight is multiplied by the number of times the synapse is selected
        """

        # updates to global parameter dictionary 
        p_update = {
        'experiment' : kwargs['experiment'], 
        'trials' : kwargs['trials'],
        'field':[-20.,0.,20.],
        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),('ica_calH','range','calH'), ('i','syn','nmda')],
        'active_paths':['1','2'],
        'gcalbar': 1.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apical_tuft','apical_trunk'],
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': [],#kwargs['syn_dist'],
        'pulses': 4.,
        'pulse_freq': 100.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        },
        '2':{
        'trees': ['apical_tuft','apical_trunk'],
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': [],#kwargs['syn_dist'],
        'pulses': 1.,
        'pulse_freq': 100.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 20,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        }}

        # setup cell and updated parameter structures
        self.P, self.cell = self._standard_parameter_setup(
            default_p='migliore_2005',
            cell_class='CellMigliore2005',
            p_update=p_update,
            paths_update=paths_update,
            load_fd=True)

        for syn_dist in kwargs['syn_dist']:
            # iterate over synapse number 
            for syn_num in kwargs['syn_num']:
                # update parameter dictionary
                self.P.paths['1']['syn_num']=syn_num
                self.P.paths['2']['syn_num']=syn_num
                # update parameter dictionary
                self.P.paths['1']['syn_dist']=syn_dist
                self.P.paths['2']['syn_dist']=syn_dist

                # iterate over trials
                for trial_i, trial in enumerate(range(self.P.p['trials'])):

                    # create list of active synapses, weights, delays
                    # stored in P.p['seg_idx', 'w_list', 'sequence_delays'], 
                    self.P = self._update_synapse_parameters(P=self.P, cell=self.cell, method='_choose_seg_rand')
                        # set sodium conductnace to zero in soma and axon
                    for tree_key, tree in self.cell.geo.iteritems():
                        if tree_key=='soma':
                            for sec_i, sec in enumerate(tree):
                                sec.gbar_na3=0.
                        if tree_key=='axon':
                            for sec_i, sec in enumerate(tree):
                                sec.gbar_nax=0. 

                    self.run_obj = self._standard_run_and_save(P=self.P, cell=self.cell, trial=trial)

    def exp_4a2_nablock(self, **kwargs):
        """ weak pathway only

        single pulse applied to compare with experiment 4a1
        
        ==Args==
        -experiment : string containing the name of the experiment, eg 'exp_1a1'
        -trials :  integer number of trials to run
        -syn_dist :  list, or nested list, containing synapse distance requirements for choose_seg_rand. eg [0, 200] chooses synapses between 0 and 200 um from soma
        -syn_num :  list of number of synapses to be activated. note that synapses are chosen randomly with replacement, so the same synapse can be selected multiple times.  In this case the weight is multiplied by the number of times the synapse is selected
        """

        # updates to global parameter dictionary 
        p_update = {
        'experiment' : kwargs['experiment'], 
        'trials' : kwargs['trials'],
        'field':[-20.,0.,20.],
        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),('ica_calH','range','calH'), ('i','syn','nmda')],
        'active_paths':['2'],
        'gcalbar': 1.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        }

        # set up synaptic pathway parameters
        paths_update = {
        '2':{
        'trees': ['apical_tuft','apical_trunk'],
        'syn_num': 0,
        'nsyns': 1.,
        'syn_dist': [],#kwargs['syn_dist'],
        'pulses': 1.,
        'pulse_freq': 100.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 20,
        'w_mean': 1*.001,
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        }}

        # setup cell and updated parameter structures
        self.P, self.cell = self._standard_parameter_setup(
            default_p='migliore_2005',
            cell_class='CellMigliore2005',
            p_update=p_update,
            paths_update=paths_update,
            load_fd=True)

        for syn_dist in kwargs['syn_dist']:
            # iterate over synapse number 
            for syn_num in kwargs['syn_num']:
                # update parameter dictionary
                # self.P.paths['1']['syn_num']=syn_num
                self.P.paths['2']['syn_num']=syn_num
                # update parameter dictionary
                # self.P.paths['1']['syn_dist']=syn_dist
                self.P.paths['2']['syn_dist']=syn_dist

                # iterate over trials
                for trial_i, trial in enumerate(range(self.P.p['trials'])):

                    # create list of active synapses, weights, delays
                    # stored in P.p['seg_idx', 'w_list', 'sequence_delays'], 
                    self.P = self._update_synapse_parameters(P=self.P, cell=self.cell, method='_choose_seg_rand')

                    # set sodium conductnace to zero in soma and axon
                    for tree_key, tree in self.cell.geo.iteritems():
                        if tree_key=='soma':
                            for sec_i, sec in enumerate(tree):
                                sec.gbar_na3=0.
                        if tree_key=='axon':
                            for sec_i, sec in enumerate(tree):
                                sec.gbar_nax=0. 

                    self.run_obj = self._standard_run_and_save(P=self.P, cell=self.cell, trial=trial)
        '''
        EXPERIMENT 5
        classic stdp experiments
        '''
    
    def exp_5a1(self, **kwargs):
        ''' reproduce classic stdp experiments in full model

        stdp experiments with delays of 10 or -10 between soma injection and epsp onset repeated 6 times at varying frequency (1,5,10,20,30,50) with 5 synapses at 0-200 um from soma

        adjust clopath parameters to qualitatively reproduce
        '''
        # updates to global parameter dictionary 
        p_update = {
        'experiment' : kwargs['experiment'], 
        'trials' : kwargs['trials'],
        'field':[0.],
        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),('ica_calH','range','calH'), ('i','syn','nmda')],
        'active_paths':['1',],
        'gcalbar': 1.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        'iclamp_amp': 1, # current injection in nA
        'iclamp_dur':5, # current injection duration ms
        'iclamp_loc': ('soma', 0,0), # location of current injection (tree, section, segment)
        'stdp_dt':kwargs['stdp_dt'], # offset between soma current injection and synaptic input (ms). positive=synaptic first, negative=soma first
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apical_tuft','apical_trunk'],
        'syn_num': 6,
        'nsyns': 1.,
        'syn_dist': kwargs['syn_dist'],
        'pulses': 6.,
        'pulse_freq': 100.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1*.001, # mean synaptic weight (microsiemens or micro-ohms)
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        },}

        # setup cell and updated parameter structures
        self.P, self.cell = self._standard_parameter_setup(
            default_p='migliore_2005',
            cell_class='CellMigliore2005',
            p_update=p_update,
            paths_update=paths_update,
            load_fd=True)

        for pulse_freq in kwargs['pulse_freqs']:
            for stdp_dt in kwargs['stdp_dt']:
                self.P.p['stdp_dt']=stdp_dt
                self.P.paths['1']['pulse_freq']=pulse_freq
                self.P.p['iclamp_delays'] = []
                for pulse_i, pulse in enumerate(range(int(self.P.paths['1']['pulses']))):
                    self.P.p['iclamp_delays'].append(self.P.paths['1']['warmup'] + self.P.p['stdp_dt'] + 1000*pulse_i/(self.P.paths['1']['pulse_freq']))
                # iterate over synapse number 
                for syn_num in kwargs['syn_num']:
                    # update parameter dictionary
                    self.P.paths['1']['syn_num']=syn_num

                    # iterate over trials
                    for trial_i, trial in enumerate(range(self.P.p['trials'])):

                        # create list of active synapses, weights, delays
                        # stored in P.p['seg_idx', 'w_list', 'sequence_delays'], 
                        self.P = self._update_synapse_parameters(P=self.P, cell=self.cell, method='_choose_seg_rand')

                        self.iclamp = stims.Intracellular()._insert_IClamp(cell=self.cell, location=self.P.p['iclamp_loc'], p=self.P.p)

                        self.run_obj = self._standard_run_and_save(P=self.P, cell=self.cell, trial=trial)
    
    def exp_5a2(self, **kwargs):
        ''' reproduce classic stdp experiments in full model

        stdp experiments at 20 Hz with varying delays (+/- 1, 2.5,5,7.5, 10,15,20,25) between soma injection and epsp onset repeated 6 times at varying frequency with 5 synapses at 0-200 um from soma

        adjust clopath parameters to qualitatively reproduce
        '''
        # updates to global parameter dictionary 
        p_update = {
        'experiment' : kwargs['experiment'], 
        'trials' : kwargs['trials'],
        'field':[0.],
        'rec_variables':[('v','range','v'),('input_times','syn','ampa'),('ica_calH','range','calH'), ('i','syn','nmda')],
        'active_paths':['1',],
        'gcalbar': 1.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
        'iclamp_amp': 1, # current injection in nA
        'iclamp_dur':5, # current injection duration ms
        'iclamp_loc': ('soma', 0,0), # location of current injection (tree, section, segment)
        'stdp_dt':kwargs['stdp_dt'], # offset between soma current injection and synaptic input (ms). positive=synaptic first, negative=soma first
        }

        # set up synaptic pathway parameters
        paths_update = {'1':{
        'trees': ['apical_tuft','apical_trunk'],
        'syn_num': 6,
        'nsyns': 1.,
        'syn_dist': kwargs['syn_dist'],
        'pulses': 6.,
        'pulse_freq': 100.,
        'bursts': 1.,
        'burst_freq': 5.,
        'warmup': 10,
        'w_mean': 1*.001, # mean synaptic weight (microsiemens or micro-ohms)
        'w_std':0.,
        'w_rand':False,
        'replace_syn':True,
        'syn_frac':0.,
        'noise':0,
        },}

        # setup cell and updated parameter structures
        self.P, self.cell = self._standard_parameter_setup(
            default_p='migliore_2005',
            cell_class='CellMigliore2005',
            p_update=p_update,
            paths_update=paths_update,
            load_fd=True)

        for pulse_freq in kwargs['pulse_freqs']:
            for stdp_dt in kwargs['stdp_dt']:
                self.P.p['stdp_dt']=stdp_dt
                self.P.paths['1']['pulse_freq']=pulse_freq
                self.P.p['iclamp_delays'] = []
                for pulse_i, pulse in enumerate(range(int(self.P.paths['1']['pulses']))):
                    self.P.p['iclamp_delays'].append(self.P.paths['1']['warmup'] + self.P.p['stdp_dt'] + 1000*pulse_i/(self.P.paths['1']['pulse_freq']))
                # iterate over synapse number 
                for syn_num in kwargs['syn_num']:
                    # update parameter dictionary
                    self.P.paths['1']['syn_num']=syn_num

                    # iterate over trials
                    for trial_i, trial in enumerate(range(self.P.p['trials'])):

                        # create list of active synapses, weights, delays
                        # stored in P.p['seg_idx', 'w_list', 'sequence_delays'], 
                        self.P = self._update_synapse_parameters(P=self.P, cell=self.cell, method='_choose_seg_rand')

                        self.iclamp = stims.Intracellular()._insert_IClamp(cell=self.cell, location=self.P.p['iclamp_loc'], p=self.P.p)

                        self.run_obj = self._standard_run_and_save(P=self.P, cell=self.cell, trial=trial)
    
    def update_parameters(self, p_update, paths_update, default_p=None, load_fd=True, data_folder='Data/',**kwargs):
        ''' setup parameter dictionary and load cell instance

        ===Args===
        -default_p  : string specifying default parameters.  calls the corresponding method in param.Default
        -p_update   : dictionary of global parameters that are specific to the current experiment. will be used to update the corresponding entries in the default dictionary specified by default_p
        -paths_update   : dictionary path-specific parameters.  Organized as paths_update{path name}{parameters}
        -cell       : string specifying the cell class to be instantiated from cell module as cell.CellClass()

        ===Out===
        -P  : default parameter class containing methods for specifying parameters
        -p  : updated parameter dictionary
                -morpho
                -seg_dist
        -paths  : updated path-specific parameter dictionary as paths{path name}{parameters}
        -cell   : cell class instance

        ===Updates===
        -P.p['fd_parameters', 'data_folder', 'fig_folder', 'seg_dist', 'morpho']

        ===Comments===
        '''
        # instantiate default parameter class
        P = param.Param()

        # load parameters from specified default parameters
        if default_p is not None:
            getattr(P, default_p)()

        # reference to default parameters
        p = P.p
        paths = P.paths

        # load facilitation depression parameters
        if load_fd:
            P._load_fd_parameters(p=p, filename='Data/fd_parameters.pkl')

        # update global parameter dictionary
        p.update(p_update)
        # update path dictionaries
        for key, val in paths_update.iteritems():
            # if path already exists
            if key in paths:
                # update
                paths[key].update(val)
            # if no paths exist
            elif not paths.values():
                # create path from paths_update
                paths[key]=val
            # if other paths exist
            else:
                # copy another path
                paths[key]=copy.copy(paths.values()[0])
                # update path
                paths[key].update(val)
        return P
    
    def _standard_parameter_setup(self, default_p, p_update, paths_update, cell_class, load_fd=True, data_folder='Data/',**kwargs):
        ''' setup parameter dictionary and load cell instance

        ===Args===
        -default_p  : string specifying default parameters.  calls the corresponding method in param.Default
        -p_update   : dictionary of global parameters that are specific to the current experiment. will be used to update the corresponding entries in the default dictionary specified by default_p
        -paths_update   : dictionary path-specific parameters.  Organized as paths_update{path name}{parameters}
        -cell       : string specifying the cell class to be instantiated from cell module as cell.CellClass()

        ===Out===
        -P  : default parameter class containing methods for specifying parameters
        -p  : updated parameter dictionary
                -morpho
                -seg_dist
        -paths  : updated path-specific parameter dictionary as paths{path name}{parameters}
        -cell   : cell class instance

        ===Updates===
        -P.p['fd_parameters', 'data_folder', 'fig_folder', 'seg_dist', 'morpho']

        ===Comments===
        '''
        # instantiate default parameter class
        P = param.Param()

        # load parameters from specified default parameters
        getattr(P, default_p)()

        # reference to default parameters
        p = P.p
        paths = P.paths

        # load facilitation depression parameters
        if load_fd:
            P._load_fd_parameters(p=p, filename='Data/fd_parameters.pkl')

        # update global parameter dictionary
        p.update(p_update)
        # update path dictionaries
        for key, val in paths_update.iteritems():
            # if path already exists
            if key in paths:
                # update
                paths[key].update(val)
            # if no paths exist
            elif not paths.values():
                # create path from paths_update
                paths[key]=val
            # if other paths exist
            else:
                # copy another path
                paths[key]=copy.copy(paths.values()[0])
                # update path
                paths[key].update(val)

        # data and figure folder
        p['data_folder'] = data_folder+p['experiment']+'/'
        p['fig_folder'] =  'png figures/'+p['experiment']+'/'

        # load cell and store in parameter dictionary
        cell = getattr(neurons, cell_class)(p)
        cell.geometry(p)
        # insert mechanisms
        cell.mechanisms(p)
        
        # measure distance of each segment from the soma and store in parameter dictionary
        p['seg_dist'] = P._seg_distance(cell)

        # create morphology for shape plots
        p['morpho'] = P._create_morpho(cell.geo)

        return P, cell

    def _update_synapse_parameters(self, P, cell, method='_choose_seg_rand',**kwargs):
        ''' update parameter dictionaries for each pathway before running simulation

        ===Args===
        -p_class    : instance of default parameter class, containing parameter dictionaries and methods for updating parameters for simulations
        -p          : full parameter dictionary
        -paths      : parameter dictionary for separate synaptic pathways, organized as paths{path name}{parameters}
        cell1       : instance of cell class containing geometry and synapse structures (contain all hoc section and synapse objects)

        ===Out===

        ===Updates===
        -p          : updated path dictionaries are added to the main parameter dictionary p (organized as p['p_path']{path name}{parameters})
        -paths      : active synapses and their activity patterns are set and updated for each path

        ===Comments===
        -p  : p should have an entry called 'path_combo', which is a list of paths to be activated during the current simulation
        '''
        p_class=P
        p = P.p
        paths=P.paths



        # add path parameter dictionaries to global p dictionary
        p['p_path']={}

        # list of segments to record from
        p['rec_idx']=[]

        # update pathways and add to global p structure
        for path_key, path in paths.iteritems():
            path['syn_idx']=[]
            path['sequence_delays']=[]
            # if path is included in current combination
            if path_key in p['active_paths']:

                if method=='_choose_seg_rand':

                    # choose active segments for this pathway
                    path['syn_idx'], path['syn_counts'] = p_class._choose_seg_rand(p=p, p_path=path, syns=cell.syns, replace=path['replace_syn'])

                # set weights for each segment in this pathway
                path['w_idx'] = p_class._set_weights_normal(p_path=path)

                # set delays for branch sequences in this pathway 
                if 'delays' in path:
                    path['sequence_delays'] = p_class._set_sequence_delays(syn_idx=path['syn_idx'], delay=path['delay'])
                else:
                    path['sequence_delays'] = p_class._set_sequence_delays(syn_idx=path['syn_idx'], delay=0)

            # record all activated synapses
            p['rec_idx'] += path['syn_idx']
            # update p_path dictionary in global p dictionary
            p['p_path'][path_key]=copy.copy(path)

        # add soma and axon to rec_idx
        if 'soma' in cell.geo:
            p['rec_idx']+=[('soma',0,0)]
        if 'axon' in cell.geo:
            p['rec_idx']+=[('axon',0,0)]

        # p['rec_idx']+=[('soma',0,0),('axon',0,0)]

        # update tstop based on synaptic inputs in each path
        tstops=[]
        warmups=[]
        for path_key, path in paths.iteritems():
            if path_key in p['active_paths']:
                # print path['sequence delays']
                if 'sequence_delays' in path and len(path['sequence_delays'])>0:
                    # print path['sequence_delays']
                    max_delay  = max(path['sequence_delays'])
                else:
                    max_delay=0
                warmups.append(path['warmup'])
                tstops.append(path['warmup'] + max_delay + 1000*(path['bursts']-1)/path['burst_freq'] + 1000*(path['pulses']+1)/path['pulse_freq'] )

            else:
                tstops = [70]
                warmups= [10]

        p['tstop'] = max(tstops)
        p['warmup'] = min(warmups)
        p['field_off'] = p['tstop']
        # FIXME
        p['field_on'] = p['warmup']-10

        return P

    def _update_synapse_parameters_sequence(self, P, cell, method='_choose_seg_from_branch', reverse=False, **kwargs):
        ''' update parameter dictionaries for each pathway before running simulation specifically for branch sequence 

        ===Args===
        -P : parameter class object
        -cell : cell object
        -method : method for choosing active synapses
        -reverse : reverses sequence order
        ===Return===
        -P : updated parameter object
        ===Updates===
        -p          : updated path dictionaries are added to the main parameter dictionary p (organized as p['p_path']{path name}{parameters})
        -paths      : active synapses and their activity patterns are set and updated for each path
        ===Comments===
        '''
        p_class=P
        p = P.p
        paths=P.paths



        # add path parameter dictionaries to global p dictionary
        p['p_path']={}

        # list of segments to record from
        p['rec_idx']=[]

        # update pathways and add to global p structure
        for path_key, path in paths.iteritems():
            path['syn_idx']=[]
            path['sequence_delays']=[]
            # if path is included in current combination
            if path_key in p['active_paths']:
                # get sec_idx
                if 'sec_idx' in kwargs:
                    sec_idx=kwargs['sec_idx']
                else:
                    terminal_branches = stims._get_terminal_branches(cell.geo)
                    sec_idx={'apical_tuft':[69]}
                # print 'geo:',cell.geo[sec_idx.keys()[0]][sec_idx.values()[0][0]].nseg
                # set branch nseg
                geo = stims._set_branch_nseg(cell.geo, sec_idx, seg_L=p['branch_seg_L'])
                # update synapses after nseg
                # cell.syns = stims._update_synapses_after_nseg(p, cell.geo, cell.syns, sec_idx)
                # print 'syns:',len(cell.syns[sec_idx.keys()[0]][sec_idx.values()[0][0]])
                # print 'geo:',cell.geo[sec_idx.keys()[0]][sec_idx.values()[0][0]].nseg

                cell.mechanisms(p)
                # print 'syns:',len(cell.syns[sec_idx.keys()[0]][sec_idx.values()[0][0]])
                # print 'syns all:', cell.syns
                # print 'geo:',cell.geo[sec_idx.keys()[0]][sec_idx.values()[0][0]].nseg
                # choose segments to activate
                path['syn_idx'], path['syn_counts'] = stims._choose_seg_from_branch(cell.geo, sec_idx)

                # FIXME adjust so that you can choose how far along the branch you want to stimulate
                if 'syn_limit' in path:
                    if len(path['syn_idx'])>path['syn_limit']:
                        # print path['syn_idx']
                        path['syn_idx'] = path['syn_idx'][-path['syn_limit']:]
                        path['syn_counts'] = path['syn_counts'][-path['syn_limit']:]
                # path['syn_idx']=path['syn_idx'][:3]
                # path['syn_counts']=path['syn_counts'][:3]

                if reverse:
                    path['syn_idx'].reverse()
                    path['syn_counts'].reverse()

                # set delays
                path['sequence_delays'] = p_class._set_sequence_delays(syn_idx=path['syn_idx'], delay=path['delay'])
                # set weights for each segment in this pathway
                path['w_idx'] = p_class._set_weights_normal(p_path=path)

            # record all activated synapses
            p['rec_idx'] += path['syn_idx']
            # update p_path dictionary in global p dictionary
            p['p_path'][path_key]=copy.copy(path)

        # add soma and axon to rec_idx
        if 'soma' in cell.geo:
            p['rec_idx']+=[('soma',0,0)]
        if 'axon' in cell.geo:
            p['rec_idx']+=[('axon',0,0)]

        # p['rec_idx']+=[('soma',0,0),('axon',0,0)]

        # update tstop based on synaptic inputs in each path
        tstops=[]
        warmups=[]
        for path_key, path in paths.iteritems():
            if path_key in p['active_paths']:
                # print path['sequence delays']
                if 'sequence_delays' in path and len(path['sequence_delays'])>0:
                    # print path['sequence_delays']
                    max_delay  = max(path['sequence_delays'])
                else:
                    max_delay=0
                warmups.append(path['warmup'])
                tstops.append(path['warmup'] + max_delay + 1000*(path['bursts']-1)/path['burst_freq'] + 1000*(path['pulses']+1)/path['pulse_freq'] )

            else:
                tstops = [70]
                warmups= [10]

        p['tstop'] = max(tstops)
        p['warmup'] = min(warmups)
        p['field_off'] = p['tstop']
        # FIXME
        p['field_on'] = p['warmup']-10
        # cell.mechanisms(p)

        return P

    def _standard_run_and_save(self, P, cell, trial, trial_id=None, **kwargs):
        '''
        '''
        # store trial number
        P.p['trial']=trial
        # data and figure folder
        P.p['data_directory'] = data_folder+p['experiment']+'/'
        P.p['fig_folder'] =  'png figures/'+p['experiment']+'/'
        
        if trial_id is None:
            trial_id = self._generate_trial_id()
        # # create unique identifier for each trial
        # uid = str(uuid.uuid1().int)[-5:]
        # now = datetime.datetime.now()
        # if trial_id is not None:
        #     trial_id = '-'.join(['{:04d}'.format(now.year), '{:02d}'.format(now.month), '{:02d}'.format(now.day), '{:02d}'.format(now.hour), '{:02d}'.format(now.minute), '{:02d}'.format(now.second), '{:02d}'.format(now.microsecond), uid])
        P.p['trial_id'] = trial_id#str(uuid.uuid4())
                    
        # start timer
        start = time.time() 
        
        self.run_obj = run.Run()
        self.run_obj._standard_run(p=P.p, cell=cell, **kwargs)

        # end timer
        end = time.time() 

        # print trial and simulation time
        print 'trial'+ str(P.p['trial']) + ' duration:' + str(end -start) 
        
        # set file name to save data
        file_name = str(
            'data_'+
            P.p['experiment']+
            '_trial_'+str(P.p['trial'])+
            '_id_'+P.p['trial_id']
            )

        # save data for eahc trial
        self.run_obj._save_data(data=self.run_obj.data, file_name=file_name, data_directory=P.p['data_directory'])

        return self.run_obj

    def _standard_run_and_save_df(self, P, cell, trial, trial_id=None, **kwargs):
        '''
        '''
        # store trial number
        P.p['trial']=trial
        # data and figure folder
        P.p['data_directory'] = 'Data/'+P.p['experiment']+'/'
        P.p['fig_folder'] =  'png figures/'+P.p['experiment']+'/'
        
        if trial_id is None:
            trial_id = self._generate_trial_id()
        # # create unique identifier for each trial
        # uid = str(uuid.uuid1().int)[-5:]
        # now = datetime.datetime.now()
        # if trial_id is not None:
        #     trial_id = '-'.join(['{:04d}'.format(now.year), '{:02d}'.format(now.month), '{:02d}'.format(now.day), '{:02d}'.format(now.hour), '{:02d}'.format(now.minute), '{:02d}'.format(now.second), '{:02d}'.format(now.microsecond), uid])
        P.p['trial_id'] = trial_id#str(uuid.uuid4())
                    
        # start timer
        start = time.time() 
        
        self.run_obj = run.Run()
        self.run_obj._standard_run_df(p=P.p, cell=cell, **kwargs)

        # end timer
        end = time.time() 

        # print trial and simulation time
        print 'trial'+ str(P.p['trial']) + ' duration:' + str(end -start) 
        
        # set file name to save data
        file_name = str(
            'data_'+
            P.p['experiment']+
            '_trial_'+str(P.p['trial'])+
            '_id_'+P.p['trial_id']
            )

        # save data for eahc trial
        self.run_obj._save_data(data=self.run_obj.data, file_name=file_name, data_directory=P.p['data_directory'])

        return self.run_obj

    def _generate_trial_id(self, ):
        '''
        '''
        # create unique identifier for each trial
        uid = str(uuid.uuid1().int)[-5:]
        now = datetime.datetime.now()
        trial_id = '-'.join(['{:04d}'.format(now.year), '{:02d}'.format(now.month), '{:02d}'.format(now.day), '{:02d}'.format(now.hour), '{:02d}'.format(now.minute), '{:02d}'.format(now.second), '{:02d}'.format(now.microsecond), uid])
        return trial_id

##########################################################################################################################################################
# class ExperimentsParallel:
#     """ Organize parameters for distributing to multiple processors
    
#     Contains a function corresponding to each experiment in Experiments class. 

#     Arguments:
#     experiment= experiment number to be run.  should be a string of the form 'exp_num', e.g. 'exp_4a'

#     Each function ExperimentsParallel.exp_num() should output a list of parameter dictionaries. Each element (dictionary) in the list will sent to a different worker/processor.  Parameters should be designed so that simulations on each worker take about the same time (ie load balnced). Specific entries in the dictionary will depend on the details of the experiment. The dictionary will be passed to the corresponding Experiments.exp_num() funtion as **kwargs

#     Once the list of parameter dictionaries is designed, the experiment can be run from command line with the syntax:

#     _run_parallel(ExperimentsParallel('exp_num', **kwargs).parameters)    
#     """
#     def __init__(self, experiment, **kwargs):
#         """ choose experiment to run and pass kwargs
#         """
        
#         # retrieve which experiment to run
#         experiment_function = getattr(self, experiment)

#         if kwargs:
#             kwargs['experiment']=experiment
#             experiment_function(**kwargs)

#         else:
#             kwargs={'experiment':experiment}
#             experiment_function(**kwargs)

#     def quick_run(self, **kwargs):
#         self.parameters=[]
#         n_workers=10
#         for i in range(n_workers):
#             self.parameters.append({'experiment':kwargs['experiment']})

#         return self.parameters

#     def exp_branch_sequence(self, **kwargs):
#         '''
#         '''
#         self.parameters = []
#         n_workers = 9
#         trials_per_worker=1
#         fields = [[-40], [0],  [40]]
#         weights = [[1*.0005],[1*.00075],[1*.001]]
#         for i, field in enumerate(fields):
#             for i, weight in enumerate(weights):
#                 self.parameters.append(
#                     {'experiment':kwargs['experiment'],
#                     'field':field,
#                     'weight':weight,
#                     })
#         return self.parameters

#     def exp_1a1(self, **kwargs):
#         '''
#         '''
#         self.parameters = []
#         n_workers = 9
#         trials_per_worker=10
#         syn_dists = [[0,300],[100,400],[200,500],[300,600],[0,200],[100,300],[200,400],[300,500],[400,600],]
#         for i in range(n_workers):
#             self.parameters.append(
#                 {'experiment':kwargs['experiment'],
#                 'trials':trials_per_worker,
#                 'syn_dist':syn_dists[i],
#                 'syn_num':[4,5,6,8,10,12,16,20]})
#         return self.parameters

#     def exp_1a1_nablock(self, **kwargs):
#         '''
#         '''
#         self.parameters = []
#         n_workers = 9
#         trials_per_worker=10
#         syn_dists = [[0,300],[100,400],[200,500],[300,600],[0,200],[100,300],[200,400],[300,500],[400,600],]
#         for i in range(n_workers):
#             self.parameters.append(
#                 {'experiment':kwargs['experiment'],
#                 'trials':trials_per_worker,
#                 'syn_dist':syn_dists[i],
#                 'syn_num':[4,5,6,8,10,12,16,20]})
#         return self.parameters

#     def exp_1a_dose(self, **kwargs):
#         '''
#         '''
#         self.parameters = []
#         n_workers = 10
#         trials_per_worker=5
#         syn_dists = [0,300]
#         for i in range(n_workers):
#             self.parameters.append(
#                 {'experiment':kwargs['experiment'],
#                 'trials':trials_per_worker,
#                 'syn_dist':syn_dists,
#                 'syn_num':[6,7,8,9,10,11,12]})
#         return self.parameters

#     def exp_1a_basal(self, **kwargs):
#         '''
#         '''
#         self.parameters = []
#         n_workers = 10
#         trials_per_worker=1
#         syn_dists = [0,300]
#         for i in range(n_workers):
#             self.parameters.append(
#                 {'experiment':kwargs['experiment'],
#                 'trials':trials_per_worker,
#                 'syn_dist':syn_dists,
#                 'syn_num':[4,5,6,7,8,9,10]})
#         return self.parameters

#     def exp_1a_apical(self, **kwargs):
#         '''
#         '''
#         self.parameters = []
#         n_workers = 10
#         trials_per_worker=1
#         syn_dists = [0,300]
#         for i in range(n_workers):
#             self.parameters.append(
#                 {'experiment':kwargs['experiment'],
#                 'trials':trials_per_worker,
#                 'syn_dist':syn_dists,
#                 'syn_num':[6,7,8,9,10,11,12]})
#         return self.parameters

#     def exp_1a_apical_full_duration(self, **kwargs):
#         '''
#         '''
#         self.parameters = []
#         n_workers = 10
#         trials_per_worker=1
#         syn_dists = [0,300]
#         for i in range(n_workers):
#             self.parameters.append(
#                 {'experiment':kwargs['experiment'],
#                 'trials':trials_per_worker,
#                 'syn_dist':syn_dists,
#                 'syn_num':[6,8,10,12]})
#         return self.parameters

#     def exp_1a_apical_full_duration_20Hz(self, **kwargs):
#         '''
#         '''
#         self.parameters = []
#         n_workers = 10
#         trials_per_worker=1
#         syn_dists = [0,300]
#         for i in range(n_workers):
#             self.parameters.append(
#                 {'experiment':kwargs['experiment'],
#                 'trials':trials_per_worker,
#                 'syn_dist':syn_dists,
#                 'syn_num':[6,8,10,12]})
#         return self.parameters

#     def exp_1a_basal_full_duration(self, **kwargs):
#         '''
#         '''
#         self.parameters = []
#         n_workers = 10
#         trials_per_worker=1
#         syn_dists = [0,300]
#         for i in range(n_workers):
#             self.parameters.append(
#                 {'experiment':kwargs['experiment'],
#                 'trials':trials_per_worker,
#                 'syn_dist':syn_dists,
#                 'syn_num':[4,6,8,10]})
#         return self.parameters

#     def exp_1a_basal_full_duration_20Hz(self, **kwargs):
#         '''
#         '''
#         self.parameters = []
#         n_workers = 1
#         trials_per_worker=10
#         syn_dists = [0,300]
#         for i in range(n_workers):
#             self.parameters.append(
#                 {'experiment':kwargs['experiment'],
#                 'trials':trials_per_worker,
#                 'syn_dist':syn_dists,
#                 'syn_num':[4,6,8,10]})
#         return self.parameters

#     def exp_1a_dose_basal(self, **kwargs):
#         '''
#         '''
#         self.parameters = []
#         n_workers = 10
#         trials_per_worker=5
#         syn_dists = [0,300]
#         for i in range(n_workers):
#             self.parameters.append(
#                 {'experiment':kwargs['experiment'],
#                 'trials':trials_per_worker,
#                 'syn_dist':syn_dists,
#                 'syn_num':[5,6,7,8,9]})
#         return self.parameters

#     def exp_1a_ACS(self, **kwargs):
#         '''
#         '''
#         self.parameters = []
#         n_workers = 10
#         trials_per_worker=1
#         syn_dists = [[0,300],]
#         for i in range(n_workers):
#             self.parameters.append(
#                 {'experiment':kwargs['experiment'],
#                 'trials':trials_per_worker,
#                 'syn_dist':syn_dists[0],
#                 'syn_num':[6,7,8,9,10,11,12]})
#         return self.parameters

#     def exp_1a_poisson(self, **kwargs):
#         '''
#         '''
#         self.parameters = []
#         n_workers = 10
#         trials_per_worker=1
#         syn_dists = [0,300]
#         for i in range(n_workers):
#             self.parameters.append(
#                 {'experiment':kwargs['experiment'],
#                 'trials':trials_per_worker,
#                 'syn_dist':syn_dists,
#                 'syn_num':[15, 20, 25, 30, 35, 40]})
#         return self.parameters

#     def exp_4a1(self, **kwargs):
#         '''
#         '''
#         self.parameters = []
#         n_workers = 9
#         trials_per_worker=2
#         syn_dists = [[0,200],[0,300]]#[[0,300],[100,400],[200,500],[300,600],[0,200],[100,300],[200,400],[300,500],[400,600],]
#         for i in range(n_workers):
#             self.parameters.append(
#                 {'experiment':kwargs['experiment'],
#                 'trials':trials_per_worker,
#                 'syn_dist':syn_dists,
#                 'syn_num':[4,5,6,8,10,12,16,20]})
#         return self.parameters

#     def exp_4a2(self, **kwargs):
#         '''
#         '''
#         self.parameters = []
#         n_workers = 9
#         trials_per_worker=2
#         syn_dists = [[0,200],[0,300]]#[[0,300],[100,400],[200,500],[300,600],[0,200],[100,300],[200,400],[300,500],[400,600],]
#         for i in range(n_workers):
#             self.parameters.append(
#                 {'experiment':kwargs['experiment'],
#                 'trials':trials_per_worker,
#                 'syn_dist':syn_dists,
#                 'syn_num':[4,5,6,8,10,12,16,20]})
#         return self.parameters

#     def exp_4a1_nablock(self, **kwargs):
#         '''
#         '''
#         self.parameters = []
#         n_workers = 9
#         trials_per_worker=4
#         syn_dists = [[0,300]]#[[0,300],[100,400],[200,500],[300,600],[0,200],[100,300],[200,400],[300,500],[400,600],]
#         for i in range(n_workers):
#             self.parameters.append(
#                 {'experiment':kwargs['experiment'],
#                 'trials':trials_per_worker,
#                 'syn_dist':syn_dists,
#                 'syn_num':[4,5,6,8,10,12,16,20]})
#         return self.parameters

#     def exp_4a2_nablock(self, **kwargs):
#         '''
#         '''
#         self.parameters = []
#         n_workers = 9
#         trials_per_worker=4
#         syn_dists = [[0,300]]#[[0,300],[100,400],[200,500],[300,600],[0,200],[100,300],[200,400],[300,500],[400,600],]
#         for i in range(n_workers):
#             self.parameters.append(
#                 {'experiment':kwargs['experiment'],
#                 'trials':trials_per_worker,
#                 'syn_dist':syn_dists,
#                 'syn_num':[4,5,6,8,10,12,16,20]})
#         return self.parameters

#     def exp_4a1_poisson(self, **kwargs):
#         '''
#         '''
#         self.parameters = []
#         n_workers = 10
#         trials_per_worker=1
#         syn_dists = [0,300]
#         path_1_mean_rate=[10,15]
#         path_2_mean_rate=[1,2,3,4,5]
#         rate_combos = list(itertools.product(path_1_mean_rate, path_2_mean_rate))
#         for i in range(n_workers):
#             self.parameters.append(
#                 {'experiment':kwargs['experiment'],
#                 'trials':trials_per_worker,
#                 'syn_dist':syn_dists,
#                 'path_1_mean_rate':rate_combos[i][0],
#                 'path_2_mean_rate':rate_combos[i][1],
#                 'syn_num':[40]})
#         return self.parameters
#     def exp_5a1(self, **kwargs):
#         '''
#         '''
#         self.parameters = []
#         n_workers = 10
#         trials_per_worker=1
#         syn_dists = [[0,200]]#[[0,300],[100,400],[200,500],[300,600],[0,200],[100,300],[200,400],[300,500],[400,600],]
#         for i in range(n_workers):
#             self.parameters.append(
#                 {'experiment':kwargs['experiment'],
#                 'trials':trials_per_worker,
#                 'syn_dist':syn_dists[0],
#                 'stdp_dt':[-10,10],
#                 'pulse_freqs':[60],
#                 'syn_num':[5]})
#         return self.parameters

#     def exp_5a2(self, **kwargs):
#         '''
#         '''
#         self.parameters = []
#         n_workers = 10
#         trials_per_worker=1
#         syn_dists = [[0,200]]#[[0,300],[100,400],[200,500],[300,600],[0,200],[100,300],[200,400],[300,500],[400,600],]
#         for i in range(n_workers):
#             self.parameters.append(
#                 {'experiment':kwargs['experiment'],
#                 'trials':trials_per_worker,
#                 'syn_dist':syn_dists[0],
#                 'stdp_dt':[-1,-2.5,-5,-7.5,-10,-15,-20,-25,1,2.5,5,7.5,10,15,20,25],
#                 'pulse_freqs':[20],
#                 'syn_num':[5]})
#         return self.parameters

#     def exp_5a3(self, **kwargs):
#         '''
#         '''
#         self.parameters = []
#         n_workers = 10
#         trials_per_worker=1
#         syn_dists = [[0,300]]#[[0,300],[100,400],[200,500],[300,600],[0,200],[100,300],[200,400],[300,500],[400,600],]
#         stdp_dts = [-10,10,-10,10,-10,10,-10,10,-10,10]
#         for i in range(n_workers):
#             self.parameters.append(
#                 {'experiment':kwargs['experiment'],
#                 'trials':trials_per_worker,
#                 'syn_dist':syn_dists[0],
#                 'stdp_dt':stdp_dts[i],
#                 'pulse_freqs':[1,5,10,20,30,40,50],
#                 'syn_num':[6,8,10]})
#         return self.parameters

# # function to pass to parallel context message board
# def _f_parallel(parameters):
#     """ Wrap experiment function so it exists in global namespace

#     Arguments: 
#     parameters - dictionary with entries 'experiment' and parameters to be passed to Experiment.exp_num.  'experiment' should be of the form 'exp_4a' and specifies which experiment to run
#     """
#     # get experiment info
#     experiment = parameters['experiment']
    
#     # create experiment class instance
#     exp_instance = Experiment()

#     # get specific experiment function
#     f = getattr(exp_instance, experiment)

#     print f
#     print parameters
#     # run experiment
#     return f(**parameters)

# # function for controlling parallel    
# def _run_parallel(parameters):
#     """ Standard run procedure for parallel simulations

#     Arguments:
#     parameters= must be a list of parameter dictionaries.  Each dictionary in the list is passed to a different worker as the arguments for Experiment.exp_num

#     Use ExperimentsParallel.exp_num to design parameters list

#     Arguments:
#     parameters= must be a list of parameter dictionaries.  Each dictionary in the list is passed to a different worker as the arguments for Experiment.exp_num

#     Use ExperimentsParallel.exp_num to design parameters list

#     To use multiple workers, python script must be called from the interpreter using syntax:
#     'mpiexec -n 10 python script.py'
#     the call to mpiexec initializes the mpi with 10 workers

#     _run_parallel should be called as:
#     if __name__=="__main__":
#         _run_parallel(ExperimentsParallel('exp_4a').parameters)
#     """

#     # make parallel context global
#     global pc

#     print parameters
#     # create parallel context instance
#     pc = h.ParallelContext()

#     print 'i am', pc.id(), 'of', pc.nhost()
#     # start workers, begins an infinitely loop where master workers posts jobs and workers pull jobs until all jobs are finished
#     pc.runworker()
    
#     print 'length of parameters:',len(parameters)
#     # # # distribute experiment and parameters to workers
#     for param in parameters:
#         # print len(parameters)
#         print 'submitted param:',param
#         pc.submit(_f_parallel, param)
#         # print param

#     # # continue runnning until all workers are finished
#     while pc.working():
#         print pc.id(), 'is working'

#     # # close parallel context 
#     pc.done()

# # function to pass to parallel context message board
# def f_parallel_milstein(parameters):
#     """ Wrap experiment function so it exists in global namespace

#     Arguments: 
#     parameters - dictionary with entries 'experiment' and parameters to be passed to Experiment.exp_num.  'experiment' should be of the form 'exp_4a' and specifies which experiment to run
#     """
#     # get experiment info
#     experiment = parameters['experiment']
    
#     # get copy of globals dictionary
#     g = globals().copy()
#     # get experiment class from globals
#     experiment_class = g[experiment]
#     # create experiment class instance
#     experiment_instance = experiment_class()
#     # get run function from experiment class
#     f = experiment_instance.run_parallel

#     # # get specific experiment function
#     # f = getattr(exp_instance, experiment)

#     # print f
#     # print parameters
#     # run experiment
#     return f(**parameters)

# # function for controlling parallel    
# def run_parallel_milstein(parameters):
#     """ Standard run procedure for parallel simulations

#     Arguments:
#     parameters= must be a list of parameter dictionaries.  Each dictionary in the list is passed to a different worker as the arguments for Experiment.exp_num

#     Use ExperimentsParallel.exp_num to design parameters list

#     Arguments:
#     parameters= must be a list of parameter dictionaries.  Each dictionary in the list is passed to a different worker as the arguments for Experiment.exp_num

#     Use ExperimentsParallel.exp_num to design parameters list

#     To use multiple workers, python script must be called from the interpreter using syntax:
#     'mpiexec -n 10 python script.py'
#     the call to mpiexec initializes the mpi with 10 workers

#     _run_parallel should be called as:
#     if __name__=="__main__":
#         _run_parallel(ExperimentsParallel('exp_4a').parameters)
#     """

#     # make parallel context global
#     global pc

#     print parameters
#     # create parallel context instance
#     pc = h.ParallelContext()

#     print 'i am', pc.id(), 'of', pc.nhost()
#     # start workers, begins an infinitely loop where master workers posts jobs and workers pull jobs until all jobs are finished
#     pc.runworker()
    
#     print 'length of parameters:',len(parameters)
#     # # # distribute experiment and parameters to workers
#     for param in parameters:
#         # print len(parameters)
#         print 'submitted param:',param
#         pc.submit(f_parallel_milstein, param)
#         # print param

#     # # continue runnning until all workers are finished
#     while pc.working():
#         print pc.id(), 'is working'

#     # # close parallel context 
#     pc.done()

if __name__ =="__main__":
    default_experiment = 'exp_reduced_neuron_tbs'
    default_parallel = True
    from_command_line_run_parallel(experiment=default_experiment, parallel=default_parallel)
    # parameters = exp_reduced_neuron_tbs().setup_parallel_parameters()
    # run_parallel(parameters=parameters)

