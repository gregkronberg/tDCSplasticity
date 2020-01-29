"""
run control
"""
# imports
from mpi4py import MPI
# import multiprocessing 
from neuron import h
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

# 
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
        'gcalbar': 1.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
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
        'syn_limit':8
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
        delays = [ 4, ]
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

                # set branch nseg
                geo = stims._set_branch_nseg(cell.geo, sec_idx, seg_L=p['branch_seg_L'])
                # update synapses after nseg
                cell.syns = stims._update_synapses_after_nseg(p, cell.geo, cell.syns, sec_idx)
                # choose segments to activate
                path['syn_idx'], path['syn_counts'] = stims._choose_seg_from_branch(cell.geo, sec_idx)

                if 'syn_limit' in path:
                    if len(path['syn_idx'])>path['syn_limit']:
                        # print path['syn_idx']
                        path['syn_idx'] = path['syn_idx'][:path['syn_limit']]
                        path['syn_counts'] = path['syn_counts'][:path['syn_limit']]
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

class ExperimentsParallel:
    """ Organize parameters for distributing to multiple processors
    
    Contains a function corresponding to each experiment in Experiments class. 

    Arguments:
    experiment= experiment number to be run.  should be a string of the form 'exp_num', e.g. 'exp_4a'

    Each function ExperimentsParallel.exp_num() should output a list of parameter dictionaries. Each element (dictionary) in the list will sent to a different worker/processor.  Parameters should be designed so that simulations on each worker take about the same time (ie load balnced). Specific entries in the dictionary will depend on the details of the experiment. The dictionary will be passed to the corresponding Experiments.exp_num() funtion as **kwargs

    Once the list of parameter dictionaries is designed, the experiment can be run from command line with the syntax:

    _run_parallel(ExperimentsParallel('exp_num', **kwargs).parameters)    
    """
    def __init__(self, experiment, **kwargs):
        """ choose experiment to run and pass kwargs
        """
        
        # retrieve which experiment to run
        experiment_function = getattr(self, experiment)

        if kwargs:
            kwargs['experiment']=experiment
            experiment_function(**kwargs)

        else:
            kwargs={'experiment':experiment}
            experiment_function(**kwargs)

    def quick_run(self, **kwargs):
        self.parameters=[]
        n_workers=10
        for i in range(n_workers):
            self.parameters.append({'experiment':kwargs['experiment']})

        return self.parameters

    def exp_branch_sequence(self, **kwargs):
        '''
        '''
        self.parameters = []
        n_workers = 9
        trials_per_worker=1
        fields = [[-20], [0],  [20]]
        weights = [[.0045],[.005],[.0055]]
        for i, field in enumerate(fields):
            for i, weight in enumerate(weights):
                self.parameters.append(
                    {'experiment':kwargs['experiment'],
                    'field':field,
                    'weight':weight,
                    })
        return self.parameters

    def exp_1a1(self, **kwargs):
        '''
        '''
        self.parameters = []
        n_workers = 9
        trials_per_worker=10
        syn_dists = [[0,300],[100,400],[200,500],[300,600],[0,200],[100,300],[200,400],[300,500],[400,600],]
        for i in range(n_workers):
            self.parameters.append(
                {'experiment':kwargs['experiment'],
                'trials':trials_per_worker,
                'syn_dist':syn_dists[i],
                'syn_num':[4,5,6,8,10,12,16,20]})
        return self.parameters

    def exp_1a1_nablock(self, **kwargs):
        '''
        '''
        self.parameters = []
        n_workers = 9
        trials_per_worker=10
        syn_dists = [[0,300],[100,400],[200,500],[300,600],[0,200],[100,300],[200,400],[300,500],[400,600],]
        for i in range(n_workers):
            self.parameters.append(
                {'experiment':kwargs['experiment'],
                'trials':trials_per_worker,
                'syn_dist':syn_dists[i],
                'syn_num':[4,5,6,8,10,12,16,20]})
        return self.parameters

    def exp_1a_dose(self, **kwargs):
        '''
        '''
        self.parameters = []
        n_workers = 10
        trials_per_worker=5
        syn_dists = [0,300]
        for i in range(n_workers):
            self.parameters.append(
                {'experiment':kwargs['experiment'],
                'trials':trials_per_worker,
                'syn_dist':syn_dists,
                'syn_num':[6,7,8,9,10,11,12]})
        return self.parameters

    def exp_1a_basal(self, **kwargs):
        '''
        '''
        self.parameters = []
        n_workers = 10
        trials_per_worker=1
        syn_dists = [0,300]
        for i in range(n_workers):
            self.parameters.append(
                {'experiment':kwargs['experiment'],
                'trials':trials_per_worker,
                'syn_dist':syn_dists,
                'syn_num':[4,5,6,7,8,9,10]})
        return self.parameters

    def exp_1a_apical(self, **kwargs):
        '''
        '''
        self.parameters = []
        n_workers = 10
        trials_per_worker=1
        syn_dists = [0,300]
        for i in range(n_workers):
            self.parameters.append(
                {'experiment':kwargs['experiment'],
                'trials':trials_per_worker,
                'syn_dist':syn_dists,
                'syn_num':[6,7,8,9,10,11,12]})
        return self.parameters

    def exp_1a_apical_full_duration(self, **kwargs):
        '''
        '''
        self.parameters = []
        n_workers = 10
        trials_per_worker=1
        syn_dists = [0,300]
        for i in range(n_workers):
            self.parameters.append(
                {'experiment':kwargs['experiment'],
                'trials':trials_per_worker,
                'syn_dist':syn_dists,
                'syn_num':[6,8,10,12]})
        return self.parameters

    def exp_1a_apical_full_duration_20Hz(self, **kwargs):
        '''
        '''
        self.parameters = []
        n_workers = 10
        trials_per_worker=1
        syn_dists = [0,300]
        for i in range(n_workers):
            self.parameters.append(
                {'experiment':kwargs['experiment'],
                'trials':trials_per_worker,
                'syn_dist':syn_dists,
                'syn_num':[6,8,10,12]})
        return self.parameters

    def exp_1a_basal_full_duration(self, **kwargs):
        '''
        '''
        self.parameters = []
        n_workers = 10
        trials_per_worker=1
        syn_dists = [0,300]
        for i in range(n_workers):
            self.parameters.append(
                {'experiment':kwargs['experiment'],
                'trials':trials_per_worker,
                'syn_dist':syn_dists,
                'syn_num':[4,6,8,10]})
        return self.parameters

    def exp_1a_basal_full_duration_20Hz(self, **kwargs):
        '''
        '''
        self.parameters = []
        n_workers = 1
        trials_per_worker=10
        syn_dists = [0,300]
        for i in range(n_workers):
            self.parameters.append(
                {'experiment':kwargs['experiment'],
                'trials':trials_per_worker,
                'syn_dist':syn_dists,
                'syn_num':[4,6,8,10]})
        return self.parameters

    def exp_1a_dose_basal(self, **kwargs):
        '''
        '''
        self.parameters = []
        n_workers = 10
        trials_per_worker=5
        syn_dists = [0,300]
        for i in range(n_workers):
            self.parameters.append(
                {'experiment':kwargs['experiment'],
                'trials':trials_per_worker,
                'syn_dist':syn_dists,
                'syn_num':[5,6,7,8,9]})
        return self.parameters

    def exp_1a_ACS(self, **kwargs):
        '''
        '''
        self.parameters = []
        n_workers = 10
        trials_per_worker=1
        syn_dists = [[0,300],]
        for i in range(n_workers):
            self.parameters.append(
                {'experiment':kwargs['experiment'],
                'trials':trials_per_worker,
                'syn_dist':syn_dists[0],
                'syn_num':[6,7,8,9,10,11,12]})
        return self.parameters

    def exp_1a_poisson(self, **kwargs):
        '''
        '''
        self.parameters = []
        n_workers = 10
        trials_per_worker=1
        syn_dists = [0,300]
        for i in range(n_workers):
            self.parameters.append(
                {'experiment':kwargs['experiment'],
                'trials':trials_per_worker,
                'syn_dist':syn_dists,
                'syn_num':[15, 20, 25, 30, 35, 40]})
        return self.parameters

    def exp_4a1(self, **kwargs):
        '''
        '''
        self.parameters = []
        n_workers = 9
        trials_per_worker=2
        syn_dists = [[0,200],[0,300]]#[[0,300],[100,400],[200,500],[300,600],[0,200],[100,300],[200,400],[300,500],[400,600],]
        for i in range(n_workers):
            self.parameters.append(
                {'experiment':kwargs['experiment'],
                'trials':trials_per_worker,
                'syn_dist':syn_dists,
                'syn_num':[4,5,6,8,10,12,16,20]})
        return self.parameters

    def exp_4a2(self, **kwargs):
        '''
        '''
        self.parameters = []
        n_workers = 9
        trials_per_worker=2
        syn_dists = [[0,200],[0,300]]#[[0,300],[100,400],[200,500],[300,600],[0,200],[100,300],[200,400],[300,500],[400,600],]
        for i in range(n_workers):
            self.parameters.append(
                {'experiment':kwargs['experiment'],
                'trials':trials_per_worker,
                'syn_dist':syn_dists,
                'syn_num':[4,5,6,8,10,12,16,20]})
        return self.parameters

    def exp_4a1_nablock(self, **kwargs):
        '''
        '''
        self.parameters = []
        n_workers = 9
        trials_per_worker=4
        syn_dists = [[0,300]]#[[0,300],[100,400],[200,500],[300,600],[0,200],[100,300],[200,400],[300,500],[400,600],]
        for i in range(n_workers):
            self.parameters.append(
                {'experiment':kwargs['experiment'],
                'trials':trials_per_worker,
                'syn_dist':syn_dists,
                'syn_num':[4,5,6,8,10,12,16,20]})
        return self.parameters

    def exp_4a2_nablock(self, **kwargs):
        '''
        '''
        self.parameters = []
        n_workers = 9
        trials_per_worker=4
        syn_dists = [[0,300]]#[[0,300],[100,400],[200,500],[300,600],[0,200],[100,300],[200,400],[300,500],[400,600],]
        for i in range(n_workers):
            self.parameters.append(
                {'experiment':kwargs['experiment'],
                'trials':trials_per_worker,
                'syn_dist':syn_dists,
                'syn_num':[4,5,6,8,10,12,16,20]})
        return self.parameters

    def exp_4a1_poisson(self, **kwargs):
        '''
        '''
        self.parameters = []
        n_workers = 10
        trials_per_worker=1
        syn_dists = [0,300]
        path_1_mean_rate=[10,15]
        path_2_mean_rate=[1,2,3,4,5]
        rate_combos = list(itertools.product(path_1_mean_rate, path_2_mean_rate))
        for i in range(n_workers):
            self.parameters.append(
                {'experiment':kwargs['experiment'],
                'trials':trials_per_worker,
                'syn_dist':syn_dists,
                'path_1_mean_rate':rate_combos[i][0],
                'path_2_mean_rate':rate_combos[i][1],
                'syn_num':[40]})
        return self.parameters
    def exp_5a1(self, **kwargs):
        '''
        '''
        self.parameters = []
        n_workers = 10
        trials_per_worker=1
        syn_dists = [[0,200]]#[[0,300],[100,400],[200,500],[300,600],[0,200],[100,300],[200,400],[300,500],[400,600],]
        for i in range(n_workers):
            self.parameters.append(
                {'experiment':kwargs['experiment'],
                'trials':trials_per_worker,
                'syn_dist':syn_dists[0],
                'stdp_dt':[-10,10],
                'pulse_freqs':[60],
                'syn_num':[5]})
        return self.parameters

    def exp_5a2(self, **kwargs):
        '''
        '''
        self.parameters = []
        n_workers = 10
        trials_per_worker=1
        syn_dists = [[0,200]]#[[0,300],[100,400],[200,500],[300,600],[0,200],[100,300],[200,400],[300,500],[400,600],]
        for i in range(n_workers):
            self.parameters.append(
                {'experiment':kwargs['experiment'],
                'trials':trials_per_worker,
                'syn_dist':syn_dists[0],
                'stdp_dt':[-1,-2.5,-5,-7.5,-10,-15,-20,-25,1,2.5,5,7.5,10,15,20,25],
                'pulse_freqs':[20],
                'syn_num':[5]})
        return self.parameters

    def exp_5a3(self, **kwargs):
        '''
        '''
        self.parameters = []
        n_workers = 10
        trials_per_worker=1
        syn_dists = [[0,300]]#[[0,300],[100,400],[200,500],[300,600],[0,200],[100,300],[200,400],[300,500],[400,600],]
        stdp_dts = [-10,10,-10,10,-10,10,-10,10,-10,10]
        for i in range(n_workers):
            self.parameters.append(
                {'experiment':kwargs['experiment'],
                'trials':trials_per_worker,
                'syn_dist':syn_dists[0],
                'stdp_dt':stdp_dts[i],
                'pulse_freqs':[1,5,10,20,30,40,50],
                'syn_num':[6,8,10]})
        return self.parameters

# function to pass to parallel context message board
def _f_parallel(parameters):
    """ Wrap experiment function so it exists in global namespace

    Arguments: 
    parameters - dictionary with entries 'experiment' and parameters to be passed to Experiment.exp_num.  'experiment' should be of the form 'exp_4a' and specifies which experiment to run
    """
    # get experiment info
    experiment = parameters['experiment']
    
    # create experiment class instance
    exp_instance = Experiment()

    # get specific experiment function
    f = getattr(exp_instance, experiment)

    print f
    print parameters
    # run experiment
    return f(**parameters)

# function for controlling parallel    
def _run_parallel(parameters):
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
        pc.submit(_f_parallel, param)
        # print param

    # # continue runnning until all workers are finished
    while pc.working():
        print pc.id(), 'is working'

    # # close parallel context 
    pc.done()

if __name__ =="__main__":
    # Experiment(experiment='exp_1a1', trials=1, syn_dist=[0,200], syn_num=[8])
    _run_parallel(ExperimentsParallel('exp_branch_sequence').parameters)