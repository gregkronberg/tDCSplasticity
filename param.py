''' parameters
''' 
from neuron import h
import numpy as np
import stims
import copy
import pickle

class BaseParam(object):
    '''
    '''
    def __init__(self):
        '''
        '''
        self.p={}
        self.paths={}

    def _build_mech_df(self, **kwargs):
        '''
        '''
        # self.p = pd.DataFrame()

        # # simulation
        # self.p['cell']=kwargs['cell_list']
        # self.p['trials']=1
        # self.p['dt']=0.05
        # self.p['tstop']=200
        # self.p['celsius']=35
        # # self.p['rec_variables'] = 'v'

        # # self.p['inhibitory_transmitters']=None
        # # self.p['excitatory_transmitters']=None
        # self.p.at[:,'inhibitory_transmitters'] = [[['gaba'] for i in range(len(self.p))]]
        # self.p.at[:,'excitatory_transmitters'] = [[['acetylcholine','glutamate'] for i in range(len(self.p))]]
        # # self.p.at[:,'inhibitory_transmitters']=[['gaba']]
        # # self.p.at[:,'excitatory_transmitters']=[['acetylcholine','glutamate',]]
        # # (value, type, mechanism name, units)
        # df.append({'param':'L','value':5,'type':'cable','mech_name':None,'units':'um'})
        # df.append({'param':'diam','value':5,'type':'cable','mech_name':None,'units':'um'})
        # df.append({'param':'v_init','value':-65.,'type':'hoc','mech_name':None,'units':'mV'})
        # df.append({'param':'cai','value':0.,'type':'ion','mech_name':None,'units':'mM'})
        # df.append({'param':'cao','value':2.,'type':'ion','mech_name':None,'units':'mM'})
        # df.append({'param':'cm','value':1.,'type':'cable','mech_name':None,'units':'uf'})
        # df.append({'param':'cm','value':1.,'type':'cable','mech_name':None,'units':'uf'})

        mech_list =[
        ('L', 5,     'cable',    None,   'um') ,
        ('diam', 5,     'cable',    None,   'um'),
        ('v_init', -65.,  'hoc',      None,   'mV'),
        ('cai', 0.0,   'ion',      None,   'mM'),
        ('cao', 2.0,   'ion',      None,   'mM'),
        ('cm', 1.0,   'cable',    None,   'uf'),
        ('Ra', 100.0, 'cable',    None,   'ohm*cm'),
        ('gmax_Leak', 5.0E-6, 'range',   'Leak', 'S/cm2'), #(S/cm2)
        ('e_Leak', -65.0, 'range',    'Leak', 'mV'),
        ('gmax_k_slow', 2.*0.003, 'range', 'k_slow','S/cm2'),
        ('ek', -60.0, 'ion',      None,   'mV'),
        ('gmax_k_fast', 7.1164395E-5, 'range', 'k_fast', 'S/cm2'),
        ('gmax_ca_boyle', 2*0.003, 'range',  'ca_boyle', 'S/cm2'),
        ('eca', 40.0, 'ion',       None,   'mV'),
        ('Exp2Syn_tau1', 1,     'synapse', 'Exp2Syn', 'ms'),
        ('Exp2Syn_tau2', 5,     'synapse', 'Exp2Syn', 'ms'),
        ('Exp2Syn_e', 0,     'synapse', 'Exp2Syn', 'mV'),
        ('Exp2Syn_i', 20,    'synapse', 'Exp2Syn','nA'),
        ('Exp2Syn_v', 'netcon_presynaptic',    'synapse', 'Exp2Syn','nA'),

        ('GradedSyn_e', 0, 'synapse', 'GradedSyn', 'mV'),
        # self.p['GradedSyn_e_inh']=-60
        ('GradedSyn_vslope', 4., 'synapse', 'GradedSyn', 'mV'),
        ('GradedSyn_vmid', 20., 'synapse', 'GradedSyn', 'mV'),
        ('GradedSyn_gbar', .0001, 'synapse', 'GradedSyn', 'uS'),
        ('GradedSyn_tau', 4., 'synapse', 'GradedSyn', 'ms'),
        ('GradedSyn_vpre', 'pointer_presynaptic', 'synapse', 'GradedSyn', 'ms'),
        ('neuron_to_neuron_elec_syn_conductance', 100E-12, 'synapse', 'neuron_to_neuron_elec_syn', 'uS') ,
        ]

        mech_df_temp = self._mech_list_to_df(mech_list)
        self.mech_df = pd.DataFrame()
        if 'cell_list' in kwargs:
            for cell in kwargs['cell_list']:
                df_temp = copy.deepcopy(mech_df_temp)
                df_temp['cell']=cell
                self.mech_df = self.mech_df.append(df_temp, ignore_index=True)
        else:
            mech_df_temp['cell']=None
            self.mech_df = mech_df_temp
    
    def _mech_list_to_df(self, mech_list, order=['param', 'val', 'mech_type', 'mech_name', 'units', 'origin', 'min_loc', 'max_loc', 'slope'], **kwargs):
        '''
        '''
        mech_dict = {}
        for i, key in enumerate(order):
            mech_dict[key]=[]
            for mech in mech_list:
                val = None
                if i<len(mech):
                    mech_dict[key].append(mech[i])
                else:
                    mech_dict[key].append(None)
                    
        df = pd.DataFrame(mech_dict)
        return df

    def _set_sequence_delays(self, syn_idx, delay):
        ''' set delays for synapses specified by syn_idx
        ==Args==
        -syn_idx : list of activated synapse locations
        -delay : delay between synapses in syn_idx
                -same delay is applied to each synapse
                -FIXME add ability to specify different delays for each synapse
        ==Out==
        -delays : list of delays (from start of simulation) for input to each synapse in syn_idx
        ==Updates==
        ==Comments==
        
        '''
        delays=[]
        for seg_i, seg in enumerate(syn_idx):
            seg_delay = seg_i*delay
            delays.append(seg_i*delay)

        return delays
    
    def _choose_seg_rand(self, p, p_path, syns, replace=False):
        """ choose random segments to activate, given a subtree and distance from soma requirement

        ==Args==
        -p :  global parameter dictionary
                    -object containing distance from soma for each segment as [tree][section][segment]
        -p_path  : synaptic pathway parameter dictionary
                    -trees: subtrees to choose synapses from
                    -syn_frac: fraction of available synapses to choose from
                            -this is overwritten if a specific number of synapses is given
                            -available synapses are determined based on the specified trees and the distance requirement
                    -syn_num: number of synapses to activate (overrides syn_frac)
                    -syn_dist: distance requirement as [min distance, max distance]
                        -can be given as a nested list with multiple distance requirements.  each distance will have the same number of synapses
                        -for more agile control, create multiple paths
        -syns: synapse structure containing hoc synapse mechanism objects as {tree}[section][segment]{synapse type}
        -replace: boolean. True=randomly choose synapses with replacement, e.g. the same synapse can be selected multiple times

        ==Out==
        -syn_idx_unique : unique list of synapse locations to be activated as [(tree, section, segment)]
        -syn_counts_unique :  list of number of occurences of the synapses specified by syn_idx_unique

        ==Updates==
        ==Comments==
        """
        
        print 'selecting synapses to activate'
        
        trees = p_path['trees']
        syn_frac=p_path['syn_frac']
        seg_dist=p['seg_dist']
        syn_num=p_path['syn_num']
        distance=p_path['syn_dist']

        # list of selected synapses as [synapse number](tree, section, segment)
        segs_chosen = []

        # list all segments as [(tree, section, segment)] 
        segs_all = [(tree_key, sec_i, seg_i) for tree_key, tree in syns.iteritems() for sec_i,sec in enumerate(tree) for seg_i,seg in enumerate(tree[sec_i]) if tree_key in trees]

        # there are multiple distance requirements
        if len(distance)>0 and isinstance(distance[0],list):
            print 'distance:',distance

            # for each distance requirement
            for distance_i, distances in enumerate(distance):
                
                # all segments that fit the current distance requirement
                segs_all_dist = [seg for seg_i, seg in enumerate(segs_all) if seg_dist[seg[0]][seg[1]][seg[2]]>distances[0] and seg_dist[seg[0]][seg[1]][seg[2]]<distances[1]] 

                # if different synapse numbers are provided for each distance bin
                if isinstance(syn_num,list) and len(syn_num)>0:
                    
                    # choose segments to activate
                    segs_choose = np.random.choice(len(segs_all_dist), int(syn_num[distance_i]), replace=replace)

                # if a single scalar is given
                elif syn_num:
                    print 'syn_num:', int(syn_num)
                    print 'available segments:',len(segs_all_dist)
                    # choose segments to activate
                    segs_choose = np.random.choice(len(segs_all_dist), int(syn_num), replace=replace)
                
                # if no synapse number is given
                else:
                    # choose segments to activate
                    segs_choose = np.random.choice(len(segs_all_dist), int(syn_frac*len(segs_all_dist)), replace=replace)

                segs_chosen += [segs_all_dist[a] for a in segs_choose]

        # if only one distacne requirement is given
        elif len(distance) > 0:

            # print 'distance requirement'
             # all segments that fit the current distance requirement
            segs_all_dist = [seg for seg_i, seg in enumerate(segs_all) if seg_dist[seg[0]][seg[1]][seg[2]]>distance[0] and seg_dist[seg[0]][seg[1]][seg[2]]<distance[1]]
            print 'available segments:',len(segs_all_dist)

            # if synapse number is given
            if syn_num:
                print 'synapses selected:', int(syn_num)
                # choose segments to activate
                segs_choose = np.random.choice(len(segs_all_dist), int(syn_num), replace=replace)

            else:
                # choose segments to activate
                segs_choose = np.random.choice(len(segs_all_dist), int(syn_frac*len(segs_all_dist)), replace=replace)

            segs_chosen += [segs_all_dist[a] for a in segs_choose]

        # if no distance requirement given
        else:
            if syn_num:
                print 'syn_num:', int(syn_num)
                # choose segments to activate
                segs_choose = np.random.choice(len(segs_all), int(syn_num), replace=replace)

            else:
                # choose segments to activate
                segs_choose = np.random.choice(len(segs_all), int(syn_frac*len(segs_all)), replace=replace)

            segs_chosen += [segs_all_dist[a] for a in segs_choose]

        # list of selected synapse locations (contains repeats)
        syn_idx = segs_chosen

        # list of number of occurences for each synapse in syn_idx
        syn_counts = [syn_idx.count(temp) for temp in syn_idx ]

        # unique list of synapse locations.  note that order may be different from syn_idx
        syn_idx_unique = list(set(syn_idx))

        # list of number of occurences for each synapse in syn_idx_unique
        syn_counts_unique = [syn_counts[syn_idx.index(temp)] for temp in syn_idx_unique]

        # get unique list of synapses and keep track of count

        return syn_idx_unique, syn_counts_unique
    
    def _set_weights_normal(self, p_path):
        """
        sets weights for synapses specified by p_path['syn_idx']

            ==Args==
            -p_path : synaptic pathway parameter dictionary
                    -syn_idx: list of unique synapse locations as[(tree, section, segment)]
                    -syn_counts: list of occurences for each synapse in syn_idx. w_mean will be multiplied by the number of occurances
                    -w_mean: mean synaptic weight in uS
                    -w_std: standard deviation for synaptic weights
                    -w_rand: Boolean. False: all weights are set to w_mean. True: weights are drawn from normal distribution specified by w_mean and w_std

            ==Out==
            -w_idx : list of weights (uS) for synapses specified by syn_idx


            ==Updates==
            ==Comments==
        """
        syn_idx=p_path['syn_idx']
        syn_counts = p_path['syn_counts']
        w_mean=p_path['w_mean']
        w_std=p_path['w_std']
        w_rand=p_path['w_rand']
        w_idx=[]

        # syn_unique = list(set(syn_idx))
        for seg_i, seg in enumerate(syn_idx):


            # tree, sec_num, seg_num = seg
            repeats = syn_counts[seg_i]

            if w_rand:
                val = np.random.normal(repeats*w_mean,w_std)
                w_idx.append(np.clip(a=val, a_min=0, a_max=val))

            else:
                w_idx.append(repeats*w_mean)

        return w_idx

    def _seg_distance(self, cell, cell_type='milstein'):
        """ calculate distance from soma of each segment and store in parameter dictionary

        ==Args==
        -geo  : geometry structure as geo[tree][section][segment]

        ==Out==
        -p['seg_dist']  : structure containing the distance of each segment from the soma, organized in the same way as geo: p['seg_dist'][tree][section][segment distance]

        ==Updates==
        -'seg_dist' is added to p

        ==Comments==
        """

        self.p['seg_dist']={}
        
        # iterate over trees
        for tree_key,tree in cell.geo.iteritems():

            # add dimension for sections
            self.p['seg_dist'][tree_key]=[]
            
            # iterate over sections
            for sec_i,sec in enumerate(tree):
                
                # add dimension for segments
                self.p['seg_dist'][tree_key].append([])
                
                # iterate over segments
                for seg_i,seg in enumerate(sec):
                    
                    # calculate and store distance from soma and store 
                    distance =  h.distance(seg.x, sec=sec)
                    self.p['seg_dist'][tree_key][sec_i].append(distance)

        return self.p['seg_dist']
    
    def _create_morpho(self, geo):
        """ create structure that stores morphology information for plotting with brian2 morphology

        each segment in morpho contains a tuple with seven entries
        (unique_segment_index, name, x, y, z, diam, unique_parent segment_index)

        root segment has index 0, with parent segment index -1
        """

        # initialize morpho structure with same dimensions as geo structure
        morpho = {}
        # iterate over trees
        for tree_key, tree in geo.iteritems():
            morpho[tree_key]=[]
            # iterate over sections
            for sec_i, sec in enumerate(tree):
                morpho[tree_key].append([])
                # iterate over segments
                for seg_i in enumerate(sec):
                    morpho[tree_key][sec_i].append([])

        # find root of cell 
        for tree_key, tree in geo.iteritems():
            for sec_i, sec in enumerate(tree):
                sref = h.SectionRef(sec=sec)
                root = sref.root
                break

        # create new secton list
        nrn_sec_list = h.SectionList()
        # add all seection to list, starting from root
        nrn_sec_list.wholetree()

        # copy nrn section list as a python list
        sec_list = []
        for sec_i_temp, sec_temp in enumerate(nrn_sec_list):
            sec_list.append(sec_temp)

        # nested list for storing segment objects [section_number][segment_number]
        seg_list= []
        # nested list for storing segment indices [section number][segment number]
        seg_list_idx = []
        # nested list for storing index of parent segment, matches seg_list_idx dimesions, [section_number][segment_number]
        parent_list_idx = []
        # keep track of total segment number
        idx = -1
        # iterate through sections in list
        for sec_i, sec in enumerate(sec_list):
            # keep track of the root section
            is_root=False


            # add section dimension to each list
            seg_list.append([])
            seg_list_idx.append([])
            parent_list_idx.append([])

            # reference for current section
            sec_ref =  h.SectionRef(sec=sec)
            
            # find parent section index
            if sec_ref.has_parent():
                parent_sec = sec_ref.parent
                parent_sec_i = [i for i, val in enumerate(sec_list) if parent_sec == val][0]
            else:
                parent_sec_i=-1
                is_root = True

            # iterate through segments in the current section
            for seg_i, seg in enumerate(sec):
                # add to total segments counter and store in lists
                idx+=1
                # copy index count to prevent overwrite during loop
                idx_count = copy.copy(idx)
                # add segment object and index to corresponding list
                seg_list[sec_i].append(seg)
                seg_list_idx[sec_i].append(idx_count)

                # if current segment is not the first in its section 
                if seg_i>0:
                    # set parent to previous segemnt in the section
                    parent_seg_idx = seg_list_idx[sec_i][seg_i-1]
                # else if it is the first segment 
                elif seg_i==0:
                    # if it is the root segment
                    if is_root:
                        parent_seg_idx=-1
                    else:
                        # set parent to the last segment in the parent section
                        parent_seg_idx = seg_list_idx[parent_sec_i][-1]

                # add to list of all parent segments
                parent_list_idx.append(parent_seg_idx)

                # find the current segment in geo structure
                # iterate through geo structure until you find matching segment
                for tree_key_local, tree_local in geo.iteritems():
                    for sec_i_local, sec_local in enumerate(tree_local):
                        for seg_i_local, seg_local in enumerate(sec_local):

                            # if section name and segment index match
                            if (sec.name() == sec_local.name()) and (seg_i == seg_i_local):
                                # segment diameter
                                diam = seg_local.diam
                                # segment xyz coordinates
                                xyz = self._seg_location(sec_local)
                                x = xyz[0][seg_i_local]
                                y = xyz[1][seg_i_local]
                                z = xyz[2][seg_i_local]

                                # segment name
                                name = tree_key_local + '_'+ str(sec_i_local) + '_'  +str(seg_i_local)
                                # create 7-tuple
                                morph_tuple = (idx_count, name, x, y, z, diam, parent_seg_idx)
                                # store in morphology structure
                                morpho[tree_key_local][sec_i_local][seg_i_local] = morph_tuple
            
                                break 
                        else:
                            continue
                        break
                    else:
                        continue
                    break
                else:
                    continue

        return morpho

    def _create_loc_list(self, geo):
        '''
        '''
        locations=[]
        for tree_key, tree in geo.iteritems():
            for sec_i, sec in enumerate(tree):
                for seg_i, seg in enumerate(sec):
                    location = (tree_key, sec_i, seg_i)
                    locations.append(location)

        return locations

    def _seg_location(self, sec):
        """ given a neuron section, output the 3d coordinates of each segment in the section

        ouput is a nested list as [xyz dimension][segment number], with x,y, z dimensions listed in that order

        """
        # number of 3d points in section
        tol =.001
        n3d = int( h.n3d( sec=sec))
        
        # preallocate 3d coordinates
        x = [None]*n3d
        y = [None]*n3d
        z = [None]*n3d
        position_3d =  [None]*n3d
                       
        # loop over 3d coordinates in each section
        for i in range(n3d):
            # retrieve x,y,z
            x[i] = h.x3d(i, sec=sec)
            y[i] = h.y3d(i, sec=sec)
            z[i] = h.z3d(i, sec=sec)

            # calculate total distance of each 3d point from start of section
            if i is 0:
                position_3d[i] = 0
            else:
                position_3d[i] = position_3d[i-1] + np.sqrt((x[i]-x[i-1])**2 + (y[i]-y[i-1])**2 + (z[i]-z[i-1])**2)
        
        seg_x = []
        seg_y = []
        seg_z = []
        for seg_i,seg in enumerate(sec):
                # relative position within section (0-1)
                seg_pos = seg.x            
                
                # segment distance along section in 3D
                seg_dist = seg_pos*position_3d[-1]

                # find first 3D coordinate that contains the segment
                node_i = [dist_i for dist_i,dist in enumerate(position_3d) if dist >= seg_dist]
                
                # if segement occurs exactly at a node set its location to the node location
                if abs(position_3d[node_i[0]] - seg_dist) < tol:
                    seg_x.append( x[ node_i[ 0]])
                    seg_y.append( z[ node_i[ 0]])
                    seg_z.append( z[ node_i[ 0]])

                # otherwise if segment falls between two coordinates, interpolate to get location
                # FIXME clean up
                else:
                    pt1 = position_3d[ node_i[0]-1]
                    pt2 = position_3d[ node_i[0]]
                    scale = (seg_dist-pt1) / (pt2-pt1)
                    interpx = x[ node_i[0]-1] + scale*( x[ node_i[0]] - x[ node_i[0]-1])
                    interpy = y[ node_i[0]-1] + scale*( y[ node_i[0]] - y[ node_i[0]-1])
                    interpz = z[ node_i[0]-1] + scale*( z[ node_i[0]] - z[ node_i[0]-1])
                    seg_x.append( interpx)
                    seg_y.append( interpy)
                    seg_z.append( interpz)
        return [seg_x, seg_y, seg_z]
    
    def _load_fd_parameters(self, P, filename):
        # load parameters
            #````````````````
        with open(filename, 'rb') as pkl_file:
            param_obj = pickle.load(pkl_file)

        params = param_obj.x

        P.p['f_ampa'] = params[0]
        P.p['tau_F_ampa'] = 1E3*params[1]
        P.p['d1_ampa'] = params[2]
        P.p['tau_D1_ampa'] = 1E3*params[3]
        P.p['d2_ampa'] = params[4]
        P.p['tau_D2_ampa'] = 1E3*params[5]
        P.p['d3_ampa'] = params[6]
        P.p['tau_D3_ampa'] = 1E3*params[7]

        return P
    
class ParamMigliore2005(BaseParam):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        super(ParamMigliore2005, self).__init__(**kwargs)
        self.define_p()

    def define_p(self, **kwargs):
        self.p = {
            'experiment' : '',
            'cell_type' : [], 
            'cell_name' : [], 
            'data_folder' : '',
            'fig_folder' : '',
            
            # equivalent cylinder parameters determined by cell.DendriteTransform() of Migliore cell geo5038804.hoc
            'L_basal' : 1600.,
            'L_soma' : 7.5,
            'L_apical_prox' : 1000.,
            'L_apical_dist' : 1000.,
            'diam1_basal' : 1.9,
            'diam1_soma' : 7.5,
            'diam1_apical_prox' : 2.75,
            'diam1_apical_dist' : 2.75,
            'diam2_basal' : 1.9,
            'diam2_soma' : 7.5,
            'diam2_apical_prox' : 2.75,
            'diam2_apical_dist' : 2.75,
            'nsec_basal' : 1,
            'nsec_soma' : 1,
            'nsec_apical_prox' : 1,
            'nsec_apical_dist' : 1,
            'syn_types' : ['ampa', 'nmda', 'clopath'],
            'fixnseg':False,        # determine number of segments in cylinder according to d_lambda rule

            # FIXME, must be set so that variable names are unique
            # set recording variables
                # organized a dictionary of dictionaries [attribute name: [variable type: mechanism]
                # note that if the attribute is part of a synapse object, it will accessed differently than a range variable
                    # range variables can be simply accessed by dot notation directly from a given neuron section
                    # synapse attributes need to be accesed from the synapse object stored in cell.syns
            'rec_variables' : 
            [('v','range','v'),
            ('t','range','t'),
            ('gbar','syn','clopath'),
            ('ica_calH','range','calH'),
            ('i', 'syn','nmda')], 

            # choose y variables to plot [varaibles]
            'plot_variables' : ['v','i','ik_kad','i_hd', 'ica_calH', 'ina_na3', 'gbar'],
            # FIXME, should be a list, where you can choose arbitrary combinations of variables 
            # x variables to plot 
            'x_variables':['t'],
            'group_trees':False,

            # synapse activation
            'syn_frac':[],      # fraction of synapses to activate with choose_seg_rand()
            'trial':0,          # count the current trial number
            'trial_id':0,       # a unique identifier for each trial using uuid64
            'w_rand':[],        # choose synapse weights from a random distribution (Bool)
            'w_std' : [],       # standard deviation of weights distribution, if w_rand is True
            'w_mean': .001,       # mean synaptic weight (microsiemens or micro-ohms)
            'trees': [],        # list of subtrees with active synapses [trees]
            'w_list':[],        # nested list of weights, determined by set_weights().  Weights correspond to segments indexed in seg_idx.  Organized as [tree][section][segment]
            'sec_list':[],      # list of active sections with repeats, each entry corresponds to the section for a given segment in seg_list.  [tree][section number]
            'seg_list':[],      # list of active segments, corresponding to sections in sec_list {tree}[segment number]
            'sec_idx': [],      # list of active sections, without repeats. Indeces in the list correspond to indeces in seg_idx {tree}[section number]
            'seg_idx':[],       # nested list of active segments {tree}[section index][segment number]
            'seg_dist' : {},    # distance of each segment from soma {tree}[section index][segment number]

            # extracellular field stimualation
            'field_angle': 0,   # angle relative to principle cell axis in radians 
            'field':[-20,0,20], # list of stimulation intensities in V/m, negative = cathodal, postivie = anodal
            'field_names':['cathodal','control','anodal'], 
            'field_colors':['b','k','r'],    # plot colors correesponding to entries in field
            'field_on':0,      # stimulation onset time in (ms)
            'field_off': 70,    # stimulation offset time in (ms)
            'dt' : .025,        # integration timestep (ms)
            'warmup': 30,       # simulation warmup time (ms)
            'tstop' : 70,       # simulation duration (ms)

            # bipolar stimulation parameters
            'bursts':1,         # bipolar stimulus bursts
            'pulses':4,         # pulses per bursts 
            'pulse_freq':100,   # pulse frequency within burst (Hz)
            'burst_freq':5,     # burst frequency (Hz)
            'noise' : 0,        # noise in input arrival (see NetCon documentation)

            # branch sequence parameters
            'num_sec':1,
            'seg_L' : 4.,
            'seg_spacing':20,
            'max_seg':[],
            'branch':False,
            'full_path':False,
            'branch_distance':[],
            'branch_seg_distance':[],
            'sequence_delay': 0,
            'sequence_direction':'in',

            # clopath synapse parameters
            'clopath_delay_steps': 1,
            'clopath_A_m':3E-5, # depression magnitude parameter (mV^-1)
            'clopath_tetam':-70,#-41, # depression threshold (mV)
            'clopath_tetap':-65,#-38, # potentiation threshold (mV)
            'clopath_tau_r':8,#-38, # potentiation threshold (mV)
            'clopath_tau_0':30,#-38, # potentiation threshold (mV)
            'clopath_tau_y': 5, # time constant (ms) for low pass filter post membrane potential for potentiation
            'clopath_A_p': 38E-5, # amplitude for potentiation (mV^-2)


            # ampa synapse parameters
            'tau1_ampa' : .1,  # rise time constant (ms)
            'tau2_ampa' : 10.,    # decay time constant   (ms)
            'i_ampa' : 0.18,    # default peak ampa current in uS

            # facilitation depression parameters for AMPA from Varela et al. 1997
            # fit to experimental theta burst and 20 Hz tetanus traces
            'f_ampa':5.,
            'tau_F_ampa':94.,
            'd1_ampa':.45,
            'tau_D1_ampa':540.,
            'd2_ampa':.12,
            'tau_D2_ampa':45.,
            'd3_ampa':.98,
            'tau_D3_ampa':120000.,

            # nmda synapse parameters
            'tau1_nmda' : 1,    # rise time constant (ms)
            'tau2_nmda' : 50,   # decay time constant (ms)
            'v0_block_nmda':0, # shift nmda voltage activation curve (mV)
            'alpha_vspom_nmda':-0.062, # slope of nmda voltage dependence 

            
            # Parameters from Migliore 2005 (signal propogation in oblique dendrites)
            # conductances reported as (nS/um2) in paper, but need to be in (mho/cm2)
            # conversion 10,000*(pS/um2) = 10*(nS/um2) = (mho/cm2) = 1000.*(mS/cm2)
            # *** units in paper are a typo, values are already reported in (mho/cm2) ***
            'Vrest' : -65.,             # leak potential (mV)
            'v_init' :-71.4,  # initialized membrane voltage, based on steady state membrane potential at soma
            'gna' :  1.*0.025,#.025,                # peak sodium conductance (mho/cm2)
            'dgna' : -.000025,          # change in sodium conductance with distance (ohm/cm2/um) from Kim 2015
            'ena' : 55.,                    # sodium reversal potential (mV)
            'gna_inact': 1., # sodium slow inactivation factor (1=none, 0=max inactivation)
            'AXONM' : 100.,              # multiplicative factor for axonal conductance to generate axon potentials in AIS
            'SOMAM':1.5,
            'gkdr' : 1.*0.01,#0.01,             # delayed rectifier potassium peak conductance (mho/cm2)
            'ek' : -80.,                    # potassium reversal potential
            'celsius' : 35.0,               # temperature (degrees C)
            'KMULT' :  1.*0.03,#0.03,           # multiplicative factor for distal A-type potassium conductances
            'KMULTP' : 1.*.03,#0.03,                # multiplicative factor for proximal A-type potassium conductances
            'ghd' : 1.*0.00005,#0.0001,         # peak h-current conductance (mho/cm2)
            'gcalbar': 0.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
            'gcabar_r': 0.*.0003 ,          # r-type calcium conductance from Kim et al. 2015 (mho/cm2)
            'gcatbar_t': 0.*.0004 ,          # t-type calcium conductance from Kim et al. 2015 (mho/cm2)
            'ehd' : -30.,                   # h-current reversal potential (mV)
            'kl_hd' : -6.,#-8.,
            'vhalfl_hd_prox' : -82.,#-73,           # activation threshold for proximal h current (mV)
            'vhalfl_hd_dist' : -90.,#-81,           # activation threshold for distal h-current (mV)
            'vhalfl_kad' : -56.,#-56.,          # inactivation threshold for distal a-type current (mV)
            'vhalfl_kap' : -56.,#-56.,          # inactivation threshold for proximal a-type current (mV)
            'vhalfn_kad' : -1.,#-1.,            # activation threshold for distal a-type urrent (mV)
            'vhalfn_kap' : -1.,#-1.,            # activation threshold for proximal a-type current (mV)
            'RaAll' : 150.,             # axial resistance, all compartments (ohm*cm)
            'RaAx' : 50.,                   # axial resistance, axon (ohm*cm)                   
            'RmAll' : 28000.,           # specific membrane resistance (ohm/cm2)
            'Cm' : 1.,                  # specific membrane capacitance (uf/cm2)
            'ka_grad' : 1.,#1.,#1.,             # slope of a-type potassium channel gradient with distance from soma 
            'ghd_grad' : 3,#1.,#3.,                # slope of h channel gradient with distance from soma
            'ka_cutoff_distance': 350, # distance from soma where ka stops increasing (um)
            'ghd_cutoff_distance': 350, # distance from soma where Ih stops increasing (um)
        }

class ParamClopath(object):
    '''
    
    for experiments comparing 20 Hz and TBS protocols:
    the plasticity rule should have a very slow time constant for LTD (e.g. 20 ms) and a fast time constant for LTP (e.g. 3 ms), but a relatively large amplitude for LTD (e.g. 2.5:1 LTD:LTP).  The threshold for LTP is what needs to be fine tuned, depending on the exact neuron model, as this threshold should be just below the peak of average synaptic inputs.  This allows subthreshold synaptic inputs to make small incremental additions to LTP, while spikes make large contributions to LTP.  For 20 Hz, the small incremental contributions accumuluate, while for TBS spikes dominate
    '''
    def __init__(self, **kwargs):
        '''
        '''
        pass

    def kronberg_2020(self, **kwargs):
        ''' set clopath parameters to values from 2020 brain stim paper
        '''
        self.param= {
            'clopath_A_m0':100E-3, # depression magnitude parameter (mV^-1)
            'clopath_tetam':-70,#-41, # depression threshold (mV)
            'clopath_tetap':-67,#-38, # potentiation threshold (mV)
            'clopath_tau_x':8,#-38, # time constant for presynaptic filter (ms)
            'clopath_tau_m':20,#-38, # time constant for depression lowpass filter
            'clopath_tau_p': 3, # time constant (ms) for low pass filter post membrane potential for potentiation
            'clopath_A_p':40E-3, # amplitude for potentiation (mV^-2)
            'clopath_delay':0, # conduction delay (ms)
            'clopath_LTD_delay':1, # conduction delay for LTD  
            'clopath_lower_bound':0.0,
            'clopath_upper_bound':None,#3.0,
            }
        return self.param

    def kronberg_2020_reduced(self, **kwargs):
        ''' set clopath parameters to values from 2020 brain stim paper
        '''
        self.param= {
            'clopath_A_m0':100E-5, # depression magnitude parameter (mV^-1)
            'clopath_A_p':40E-5, # amplitude for potentiation (mV^-2)
            'clopath_tetam':-70,#-41, # depression threshold (mV)
            'clopath_tetap':-61,#-38, # potentiation threshold (mV)
            'clopath_tau_x':8,#-38, # time constant for presynaptic filter (ms)
            'clopath_tau_m':20,#-38, # time constant for depression lowpass filter
            'clopath_tau_p': 3, # time constant (ms) for low pass filter post membrane potential for potentiation
            
            'clopath_delay':0, # conduction delay (ms)
            'clopath_LTD_delay':1, # conduction delay for LTD  
            'clopath_lower_bound':0.0,
            'clopath_upper_bound':None,#3.0,
            }
        return self.param

    def kronberg_2020_test_1hz(self, **kwargs):
        ''' set clopath parameters to values from 2020 brain stim paper
        '''
        self.param= {
            'clopath_A_m0':20*100E-5, # depression magnitude parameter (mV^-1)
            'clopath_A_p':18*3*40E-5, # amplitude for potentiation (mV^-2)
            'clopath_tetam':-72,#-41, # depression threshold (mV)
            'clopath_tetap':-59,#-38, # potentiation threshold (mV)
            'clopath_tau_x':8,#-38, # time constant for presynaptic filter (ms)
            'clopath_tau_m':20,#-38, # time constant for depression lowpass filter
            'clopath_tau_p': 3, # time constant (ms) for low pass filter post membrane potential for potentiation
            
            'clopath_delay':0, # conduction delay (ms)
            'clopath_LTD_delay':1, # conduction delay for LTD  
            'clopath_lower_bound':0.0,
            'clopath_upper_bound':None,#3.0,
            }
        return self.param

    def kronberg_2020_test_1hz_2(self, **kwargs):
        ''' set clopath parameters to values from 2020 brain stim paper
        '''
        self.param= {
            'clopath_A_m0':100E-5, # depression magnitude parameter (mV^-1)
            'clopath_A_p':40E-5, # amplitude for potentiation (mV^-2)
            'clopath_tetam':-72,#-41, # depression threshold (mV)
            'clopath_tetap':-58,#-38, # potentiation threshold (mV)
            'clopath_tau_x':8,#-38, # time constant for presynaptic filter (ms)
            'clopath_tau_m':20,#-38, # time constant for depression lowpass filter
            'clopath_tau_p': 3, # time constant (ms) for low pass filter post membrane potential for potentiation
            
            'clopath_delay':0, # conduction delay (ms)
            'clopath_LTD_delay':1, # conduction delay for LTD  
            'clopath_lower_bound':0.0,
            'clopath_upper_bound':None,#3.0,
            }
        return self.param


#############################################################################
# deprecated
#############################################################################
class ___Param(object):
    """ base class for experimental parameters
    """
    def __init__(self):
        '''
        '''
        self.p={}
        self.paths={}
        pass

    def default_parameters(self):

      
        exp='default'
        self.p = {
            'experiment' : exp,
            'cell' : [], 
            'data_folder' : 'Data/'+exp+'/',
            'fig_folder' : 'png figures/'+exp+'/',
            
            # equivalent cylinder parameters determined by cell.DendriteTransform() of Migliore cell geo5038804.hoc
            'L_basal' : 1600.,
            'L_soma' : 7.5,
            'L_apical_prox' : 1000.,
            'L_apical_dist' : 1000.,
            'diam1_basal' : 1.9,
            'diam1_soma' : 7.5,
            'diam1_apical_prox' : 2.75,
            'diam1_apical_dist' : 2.75,
            'diam2_basal' : 1.9,
            'diam2_soma' : 7.5,
            'diam2_apical_prox' : 2.75,
            'diam2_apical_dist' : 2.75,
            'nsec_basal' : 1,
            'nsec_soma' : 1,
            'nsec_apical_prox' : 1,
            'nsec_apical_dist' : 1,
            'syn_types' : ['ampa', 'nmda', 'clopath'],
            'fixnseg':False,        # determine number of segments in cylinder according to d_lambda rule
            'nseg':1,

            # FIXME, must be set so that variable names are unique
            # set recording variables
                # organized a dictionary of dictionaries [attribute name: [variable type: mechanism]
                # note that if the attribute is part of a synapse object, it will accessed differently than a range variable
                    # range variables can be simply accessed by dot notation directly from a given neuron section
                    # synapse attributes need to be accesed from the synapse object stored in cell.syns
            'rec_variables' : 
            [('v','range','v'),
            ('gbar','syn','clopath'),
            ('ica_calH','range','calH'),
            ('input_times','syn','ampa')], 

            # choose y variables to plot [varaibles]
            'plot_variables' : ['v','i','ik_kad','i_hd', 'ica_calH', 'ina_na3', 'gbar'],
            # FIXME, should be a list, where you can choose arbitrary combinations of variables 
            # x variables to plot 
            'x_variables':['t'],
            'group_trees':False,

            # synapse activation
            'syn_frac':[],      # fraction of synapses to activate with choose_seg_rand()
            'trial':0,          # count the current trial number
            'trial_id':0,       # a unique identifier for each trial using uuid64
            'w_rand':[],        # choose synapse weights from a random distribution (Bool)
            'w_std' : [],       # standard deviation of weights distribution, if w_rand is True
            'w_mean': [], 		# mean synaptic weight (microsiemens or micro-ohms)
            'trees': [],        # list of subtrees with active synapses [trees]
            'w_list':[],        # nested list of weights, determined by set_weights().  Weights correspond to segments indexed in seg_idx.  Organized as [tree][section][segment]
            'sec_list':[],      # list of active sections with repeats, each entry corresponds to the section for a given segment in seg_list.  [tree][section number]
            'seg_list':[],      # list of active segments, corresponding to sections in sec_list {tree}[segment number]
            'sec_idx': [],      # list of active sections, without repeats. Indeces in the list correspond to indeces in seg_idx {tree}[section number]
            'seg_idx':[],       # nested list of active segments {tree}[section index][segment number]
            'seg_dist' : {},    # distance of each segment from soma {tree}[section index][segment number]

            # extracellular field stimualation
            'field_angle': 0,   # angle relative to principle cell axis in radians 
            'field':[-20,0,20], # list of stimulation intensities in V/m, negative = cathodal, postivie = anodal
            'field_color':['b','k','r'],    # plot colors correesponding to entries in field
            'field_on':20,      # stimulation onset time in (ms)
            'field_off': 70,    # stimulation offset time in (ms)
            'dt' : .025,        # integration timestep (ms)
            'warmup': 30,       # simulation warmup time (ms)
            'tstop' : 70,       # simulation duration (ms)

            # bipolar stimulation parameters
            'bursts':1,         # bipolar stimulus bursts
            'pulses':4,         # pulses per bursts 
            'pulse_freq':100,   # pulse frequency within burst (Hz)
            'burst_freq':5,     # burst frequency (Hz)
            'noise' : 0,        # noise in input arrival (see NetCon documentation)

            # clopath synapse parameters
            'clopath_delay_steps': 1,
            'clopath_tau_0':6, # time constant (ms) for low passed membrane potential for depression
            'clopath_tau_r' : 10, # time constant (ms) for low pass filter presynaptic variable
            'clopath_tau_y': 5, # time constant (ms) for low pass filter post membrane potential for potentiation
            'clopath_A_m':2E-5, # depression magnitude parameter (mV^-1)
            'clopath_A_p': 38E-5, # amplitude for potentiation (mV^-2)
            'clopath_tetam':-60,#-41, # depression threshold (mV)
            'clopath_tetap':-53,#-38, # potentiation threshold (mV)



            # ampa synapse parameters
            'tau1_ampa' : 0.2,  # rise time constant (ms)
            'tau2_ampa' : 2,    # decay time constant   (ms)
            'i_ampa' : 0.18,    # default peak ampa current in uS

            # facilitation depression parameters for AMPA from Varela et al. 1997
            # fit to experimental theta burst and 20 Hz tetanus traces
            'f_ampa':5.,
            'tau_F_ampa':94.,
            'd1_ampa':.45,
            'tau_D1_ampa':540.,
            'd2_ampa':.12,
            'tau_D2_ampa':45.,
            'd3_ampa':.98,
            'tau_D3_ampa':120000.,

            # nmda synapse parameters
            'tau1_nmda' : 1,    # rise time constant (ms)
            'tau2_nmda' : 50,   # decay time constant (ms)

            
            # Parameters from Migliore 2005 (signal propogation in oblique dendrites)
            # conductances reported as (nS/um2) in paper, but need to be in (mho/cm2)
            # conversion 10,000*(pS/um2) = 10*(nS/um2) = (mho/cm2) = .001*(mS/cm2)
            # *** units in paper are a typo, values are already reported in (mho/cm2) ***
            'Vrest' : -65.,             # resting potential (mV)
            'gna' :  1.*0.04,#.025,                # peak sodium conductance (mho/cm2)
            'dgna' : -.000025,          # change in sodium conductance with distance (ohm/cm2/um) from Kim 2015
            'ena' : 55.,                    # sodium reversal potential (mV)
            'gna_inact': 0., # sodium slow inactivation factor (1=none, 0=max inactivation)
            'AXONM' : 50.,              # multiplicative factor for axonal conductance to generate axon potentials in AIS
            'SOMAM':1.5,
            'gkdr' : 1.*0.01,#0.01,             # delayed rectifier potassium peak conductance (mho/cm2)
            'ek' : -90.,                    # potassium reversal potential
            'celsius' : 35.0,               # temperature (degrees C)
            'KMULT' :  1.*0.03,#0.03,           # multiplicative factor for distal A-type potassium conductances
            'KMULTP' : 1.*.03,#0.03,                # multiplicative factor for proximal A-type potassium conductances
            'ghd' : 1.*0.0001,#0.0001,         # peak h-current conductance (mho/cm2)
            'gcalbar': 1.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
            'ehd' : -30.,                   # h-current reversal potential (mV)
            'kl_hd' : -6.,#-8.,
            'vhalfl_hd_prox' : -83.,#-73,           # activation threshold for proximal h current (mV)
            'vhalfl_hd_dist' : -83.,#-81,           # activation threshold for distal h-current (mV)
            'vhalfl_kad' : -56.,#-56.,          # inactivation threshold for distal a-type current (mV)
            'vhalfl_kap' : -56.,#-56.,          # inactivation threshold for proximal a-type current (mV)
            'vhalfn_kad' : -1.,#-1.,            # activation threshold for distal a-type urrent (mV)
            'vhalfn_kap' : -1.,#-1.,            # activation threshold for proximal a-type current (mV)
            'RaAll' : 150.,             # axial resistance, all compartments (ohm*cm)
            'RaAx' : 50.,                   # axial resistance, axon (ohm*cm)                   
            'RmAll' : 28000.,           # specific membrane resistance (ohm/cm2)
            'Cm' : 1.,                  # specific membrane capacitance (uf/cm2)
            'ka_grad' : 1.,#1.,#1.,             # slope of a-type potassium channel gradient with distance from soma 
            'ghd_grad' : 1.5,#1.,#3.,                # slope of h channel gradient with distance from soma 
            }
    
    def migliore_2005(self):
        self.p = {
            'experiment' : '',
            'cell' : [], 
            'data_folder' : '',
            'fig_folder' : '',
            
            # equivalent cylinder parameters determined by cell.DendriteTransform() of Migliore cell geo5038804.hoc
            'L_basal' : 1600.,
            'L_soma' : 7.5,
            'L_apical_prox' : 1000.,
            'L_apical_dist' : 1000.,
            'diam1_basal' : 1.9,
            'diam1_soma' : 7.5,
            'diam1_apical_prox' : 2.75,
            'diam1_apical_dist' : 2.75,
            'diam2_basal' : 1.9,
            'diam2_soma' : 7.5,
            'diam2_apical_prox' : 2.75,
            'diam2_apical_dist' : 2.75,
            'nsec_basal' : 1,
            'nsec_soma' : 1,
            'nsec_apical_prox' : 1,
            'nsec_apical_dist' : 1,
            'syn_types' : ['ampa', 'nmda', 'clopath'],
            'fixnseg':False,        # determine number of segments in cylinder according to d_lambda rule

            # FIXME, must be set so that variable names are unique
            # set recording variables
                # organized a dictionary of dictionaries [attribute name: [variable type: mechanism]
                # note that if the attribute is part of a synapse object, it will accessed differently than a range variable
                    # range variables can be simply accessed by dot notation directly from a given neuron section
                    # synapse attributes need to be accesed from the synapse object stored in cell.syns
            'rec_variables' : 
            [('v','range','v'),
            ('t','range','t'),
            ('gbar','syn','clopath'),
            ('ica_calH','range','calH'),
            ('i', 'syn','nmda')], 

            # choose y variables to plot [varaibles]
            'plot_variables' : ['v','i','ik_kad','i_hd', 'ica_calH', 'ina_na3', 'gbar'],
            # FIXME, should be a list, where you can choose arbitrary combinations of variables 
            # x variables to plot 
            'x_variables':['t'],
            'group_trees':False,

            # synapse activation
            'syn_frac':[],      # fraction of synapses to activate with choose_seg_rand()
            'trial':0,          # count the current trial number
            'trial_id':0,       # a unique identifier for each trial using uuid64
            'w_rand':[],        # choose synapse weights from a random distribution (Bool)
            'w_std' : [],       # standard deviation of weights distribution, if w_rand is True
            'w_mean': .001,       # mean synaptic weight (microsiemens or micro-ohms)
            'trees': [],        # list of subtrees with active synapses [trees]
            'w_list':[],        # nested list of weights, determined by set_weights().  Weights correspond to segments indexed in seg_idx.  Organized as [tree][section][segment]
            'sec_list':[],      # list of active sections with repeats, each entry corresponds to the section for a given segment in seg_list.  [tree][section number]
            'seg_list':[],      # list of active segments, corresponding to sections in sec_list {tree}[segment number]
            'sec_idx': [],      # list of active sections, without repeats. Indeces in the list correspond to indeces in seg_idx {tree}[section number]
            'seg_idx':[],       # nested list of active segments {tree}[section index][segment number]
            'seg_dist' : {},    # distance of each segment from soma {tree}[section index][segment number]

            # extracellular field stimualation
            'field_angle': 0,   # angle relative to principle cell axis in radians 
            'field':[-20,0,20], # list of stimulation intensities in V/m, negative = cathodal, postivie = anodal
            'field_names':['cathodal','control','anodal'], 
            'field_colors':['b','k','r'],    # plot colors correesponding to entries in field
            'field_on':0,      # stimulation onset time in (ms)
            'field_off': 70,    # stimulation offset time in (ms)
            'dt' : .025,        # integration timestep (ms)
            'warmup': 30,       # simulation warmup time (ms)
            'tstop' : 70,       # simulation duration (ms)

            # bipolar stimulation parameters
            'bursts':1,         # bipolar stimulus bursts
            'pulses':4,         # pulses per bursts 
            'pulse_freq':100,   # pulse frequency within burst (Hz)
            'burst_freq':5,     # burst frequency (Hz)
            'noise' : 0,        # noise in input arrival (see NetCon documentation)

            # branch sequence parameters
            'num_sec':1,
            'seg_L' : 4.,
            'seg_spacing':20,
            'max_seg':[],
            'branch':False,
            'full_path':False,
            'branch_distance':[],
            'branch_seg_distance':[],
            'sequence_delay': 0,
            'sequence_direction':'in',

            # clopath synapse parameters
            'clopath_delay_steps': 1,
            'clopath_A_m':3E-5, # depression magnitude parameter (mV^-1)
            'clopath_tetam':-70,#-41, # depression threshold (mV)
            'clopath_tetap':-65,#-38, # potentiation threshold (mV)
            'clopath_tau_r':8,#-38, # potentiation threshold (mV)
            'clopath_tau_0':30,#-38, # potentiation threshold (mV)
            'clopath_tau_y': 5, # time constant (ms) for low pass filter post membrane potential for potentiation
            'clopath_A_p': 38E-5, # amplitude for potentiation (mV^-2)


            # ampa synapse parameters
            'tau1_ampa' : 0.2,  # rise time constant (ms)
            'tau2_ampa' : 2,    # decay time constant   (ms)
            'i_ampa' : 0.18,    # default peak ampa current in uS

            # facilitation depression parameters for AMPA from Varela et al. 1997
            # fit to experimental theta burst and 20 Hz tetanus traces
            'f_ampa':5.,
            'tau_F_ampa':94.,
            'd1_ampa':.45,
            'tau_D1_ampa':540.,
            'd2_ampa':.12,
            'tau_D2_ampa':45.,
            'd3_ampa':.98,
            'tau_D3_ampa':120000.,

            # nmda synapse parameters
            'tau1_nmda' : 1,    # rise time constant (ms)
            'tau2_nmda' : 50,   # decay time constant (ms)

            
            # Parameters from Migliore 2005 (signal propogation in oblique dendrites)
            # conductances reported as (nS/um2) in paper, but need to be in (mho/cm2)
            # conversion 10,000*(pS/um2) = 10*(nS/um2) = (mho/cm2) = 1000.*(mS/cm2)
            # *** units in paper are a typo, values are already reported in (mho/cm2) ***
            'Vrest' : -65.,             # leak potential (mV)
            'v_init' :-71.4,  # initialized membrane voltage, based on steady state membrane potential at soma
            'gna' :  1.*0.025,#.025,                # peak sodium conductance (mho/cm2)
            'dgna' : 0.,#-.000025,          # change in sodium conductance with distance (ohm/cm2/um) from Kim 2015
            'ena' : 55.,                    # sodium reversal potential (mV)
            'gna_inact': 0., # sodium slow inactivation factor (1=none, 0=max inactivation)
            'AXONM' : 100.,              # multiplicative factor for axonal conductance to generate axon potentials in AIS
            'SOMAM':1.5,
            'gkdr' : 1.*0.01,#0.01,             # delayed rectifier potassium peak conductance (mho/cm2)
            'ek' : -90.,                    # potassium reversal potential
            'celsius' : 35.0,               # temperature (degrees C)
            'KMULT' :  1.*0.03,#0.03,           # multiplicative factor for distal A-type potassium conductances
            'KMULTP' : 1.*.03,#0.03,                # multiplicative factor for proximal A-type potassium conductances
            'ghd' : 1.*0.00005,#0.0001,         # peak h-current conductance (mho/cm2)
            'gcalbar': 0.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
            'ehd' : -30.,                   # h-current reversal potential (mV)
            'kl_hd' : -6.,#-8.,
            'vhalfl_hd_prox' : -82.,#-73,           # activation threshold for proximal h current (mV)
            'vhalfl_hd_dist' : -90.,#-81,           # activation threshold for distal h-current (mV)
            'vhalfl_kad' : -56.,#-56.,          # inactivation threshold for distal a-type current (mV)
            'vhalfl_kap' : -56.,#-56.,          # inactivation threshold for proximal a-type current (mV)
            'vhalfn_kad' : -1.,#-1.,            # activation threshold for distal a-type urrent (mV)
            'vhalfn_kap' : -1.,#-1.,            # activation threshold for proximal a-type current (mV)
            'RaAll' : 150.,             # axial resistance, all compartments (ohm*cm)
            'RaAx' : 50.,                   # axial resistance, axon (ohm*cm)                   
            'RmAll' : 28000.,           # specific membrane resistance (ohm/cm2)
            'Cm' : 1.,                  # specific membrane capacitance (uf/cm2)
            'ka_grad' : 1.,#1.,#1.,             # slope of a-type potassium channel gradient with distance from soma 
            'ghd_grad' : 3,#1.,#3.,                # slope of h channel gradient with distance from soma
            'ka_cutoff_distance': 350, # distance from soma where ka stops increasing (um)
            'ghd_cutoff_distance': 350, # distance from soma where Ih stops increasing (um)
        }

    def branco_2010(self):
        self.p = {
            'experiment' : '',
            'cell' : [], 
            'data_folder' : '',
            'fig_folder' : '',
            
            # equivalent cylinder parameters determined by cell.DendriteTransform() of Migliore cell geo5038804.hoc
            'L_basal' : 1600.,
            'L_soma' : 7.5,
            'L_apical_prox' : 1000.,
            'L_apical_dist' : 1000.,
            'diam1_basal' : 1.9,
            'diam1_soma' : 7.5,
            'diam1_apical_prox' : 2.75,
            'diam1_apical_dist' : 2.75,
            'diam2_basal' : 1.9,
            'diam2_soma' : 7.5,
            'diam2_apical_prox' : 2.75,
            'diam2_apical_dist' : 2.75,
            'nsec_basal' : 1,
            'nsec_soma' : 1,
            'nsec_apical_prox' : 1,
            'nsec_apical_dist' : 1,
            'syn_types' : ['ampa', 'nmda', 'clopath'],
            'fixnseg':False,        # determine number of segments in cylinder according to d_lambda rule

            # FIXME, must be set so that variable names are unique
            # set recording variables
                # organized a dictionary of dictionaries [attribute name: [variable type: mechanism]
                # note that if the attribute is part of a synapse object, it will accessed differently than a range variable
                    # range variables can be simply accessed by dot notation directly from a given neuron section
                    # synapse attributes need to be accesed from the synapse object stored in cell.syns
            'rec_variables' : 
            [('v','range','v'),
            ('t','range','t'),
            ('gbar','syn','clopath'),
            ('ica_calH','range','calH'),
            ('i', 'syn','nmda')], 

            # choose y variables to plot [varaibles]
            'plot_variables' : ['v','i','ik_kad','i_hd', 'ica_calH', 'ina_na3', 'gbar'],
            # FIXME, should be a list, where you can choose arbitrary combinations of variables 
            # x variables to plot 
            'x_variables':['t'],
            'group_trees':False,

            # synapse activation
            'syn_frac':[],      # fraction of synapses to activate with choose_seg_rand()
            'trial':0,          # count the current trial number
            'trial_id':0,       # a unique identifier for each trial using uuid64
            'w_rand':[],        # choose synapse weights from a random distribution (Bool)
            'w_std' : [],       # standard deviation of weights distribution, if w_rand is True
            'w_mean': .001,       # mean synaptic weight (microsiemens or micro-ohms)
            'trees': [],        # list of subtrees with active synapses [trees]
            'w_list':[],        # nested list of weights, determined by set_weights().  Weights correspond to segments indexed in seg_idx.  Organized as [tree][section][segment]
            'sec_list':[],      # list of active sections with repeats, each entry corresponds to the section for a given segment in seg_list.  [tree][section number]
            'seg_list':[],      # list of active segments, corresponding to sections in sec_list {tree}[segment number]
            'sec_idx': [],      # list of active sections, without repeats. Indeces in the list correspond to indeces in seg_idx {tree}[section number]
            'seg_idx':[],       # nested list of active segments {tree}[section index][segment number]
            'seg_dist' : {},    # distance of each segment from soma {tree}[section index][segment number]

            # extracellular field stimualation
            'field_angle': 0,   # angle relative to principle cell axis in radians 
            'field':[-20,0,20], # list of stimulation intensities in V/m, negative = cathodal, postivie = anodal
            'field_names':['cathodal','control','anodal'], 
            'field_colors':['b','k','r'],    # plot colors correesponding to entries in field
            'field_on':0,      # stimulation onset time in (ms)
            'field_off': 70,    # stimulation offset time in (ms)
            'dt' : .025,        # integration timestep (ms)
            'warmup': 30,       # simulation warmup time (ms)
            'tstop' : 70,       # simulation duration (ms)

            # bipolar stimulation parameters
            'bursts':1,         # bipolar stimulus bursts
            'pulses':4,         # pulses per bursts 
            'pulse_freq':100,   # pulse frequency within burst (Hz)
            'burst_freq':5,     # burst frequency (Hz)
            'noise' : 0,        # noise in input arrival (see NetCon documentation)

            # branch sequence parameters
            'num_sec':1,
            'seg_L' : 4.,
            'seg_spacing':20,
            'max_seg':[],
            'branch':False,
            'full_path':False,
            'branch_distance':[],
            'branch_seg_distance':[],
            'sequence_delay': 0,
            'sequence_direction':'in',

            # clopath synapse parameters
            'clopath_delay_steps': 1,
            'clopath_A_m':3E-5, # depression magnitude parameter (mV^-1)
            'clopath_tetam':-70,#-41, # depression threshold (mV)
            'clopath_tetap':-65,#-38, # potentiation threshold (mV)
            'clopath_tau_r':8,#-38, # potentiation threshold (mV)
            'clopath_tau_0':30,#-38, # potentiation threshold (mV)
            'clopath_tau_y': 5, # time constant (ms) for low pass filter post membrane potential for potentiation
            'clopath_A_p': 38E-5, # amplitude for potentiation (mV^-2)


            # ampa synapse parameters
            'tau1_ampa' : 0.2,  # rise time constant (ms)
            'tau2_ampa' : 2,    # decay time constant   (ms)
            'i_ampa' : 0.18,    # default peak ampa current in uS

            # facilitation depression parameters for AMPA from Varela et al. 1997
            # fit to experimental theta burst and 20 Hz tetanus traces
            'f_ampa':5.,
            'tau_F_ampa':94.,
            'd1_ampa':.45,
            'tau_D1_ampa':540.,
            'd2_ampa':.12,
            'tau_D2_ampa':45.,
            'd3_ampa':.98,
            'tau_D3_ampa':120000.,

            # nmda synapse parameters
            'tau1_nmda' : 1,    # rise time constant (ms)
            'tau2_nmda' : 50,   # decay time constant (ms)

            
            # Parameters from Migliore 2005 (signal propogation in oblique dendrites)
            # conductances reported as (nS/um2) in paper, but need to be in (mho/cm2)
            # conversion 10,000*(pS/um2) = 10*(nS/um2) = (mho/cm2) = 1000.*(mS/cm2)
            # *** units in paper are a typo, values are already reported in (mho/cm2) ***
            'Vrest' : -65.,             # leak potential (mV)
            'v_init' :-71.4,  # initialized membrane voltage, based on steady state membrane potential at soma
            'gna' :  1.*0.025,#.025,                # peak sodium conductance (mho/cm2)
            'dgna' : 0.,#-.000025,          # change in sodium conductance with distance (ohm/cm2/um) from Kim 2015
            'ena' : 55.,                    # sodium reversal potential (mV)
            'gna_inact': 0., # sodium slow inactivation factor (1=none, 0=max inactivation)
            'AXONM' : 100.,              # multiplicative factor for axonal conductance to generate axon potentials in AIS
            'SOMAM':1.5,
            'gkdr' : 1.*0.01,#0.01,             # delayed rectifier potassium peak conductance (mho/cm2)
            'ek' : -90.,                    # potassium reversal potential
            'celsius' : 35.0,               # temperature (degrees C)
            'KMULT' :  1.*0.03,#0.03,           # multiplicative factor for distal A-type potassium conductances
            'KMULTP' : 1.*.03,#0.03,                # multiplicative factor for proximal A-type potassium conductances
            'ghd' : 1.*0.00005,#0.0001,         # peak h-current conductance (mho/cm2)
            'gcalbar': 0.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
            'ehd' : -30.,                   # h-current reversal potential (mV)
            'kl_hd' : -6.,#-8.,
            'vhalfl_hd_prox' : -82.,#-73,           # activation threshold for proximal h current (mV)
            'vhalfl_hd_dist' : -90.,#-81,           # activation threshold for distal h-current (mV)
            'vhalfl_kad' : -56.,#-56.,          # inactivation threshold for distal a-type current (mV)
            'vhalfl_kap' : -56.,#-56.,          # inactivation threshold for proximal a-type current (mV)
            'vhalfn_kad' : -1.,#-1.,            # activation threshold for distal a-type urrent (mV)
            'vhalfn_kap' : -1.,#-1.,            # activation threshold for proximal a-type current (mV)
            'RaAll' : 150.,             # axial resistance, all compartments (ohm*cm)
            'RaAx' : 50.,                   # axial resistance, axon (ohm*cm)                   
            'RmAll' : 28000.,           # specific membrane resistance (ohm/cm2)
            'Cm' : 1.,                  # specific membrane capacitance (uf/cm2)
            'ka_grad' : 1.,#1.,#1.,             # slope of a-type potassium channel gradient with distance from soma 
            'ghd_grad' : 3,#1.,#3.,                # slope of h channel gradient with distance from soma
            'ka_cutoff_distance': 350, # distance from soma where ka stops increasing (um)
            'ghd_cutoff_distance': 350, # distance from soma where Ih stops increasing (um)
        }

    def pyramidal_cylinder(self):
        self.p = {
            'experiment' : '',
            'cell' : [], 
            'data_folder' : '',
            'fig_folder' : '',
            
            # equivalent cylinder parameters determined by cell.DendriteTransform() of Migliore cell geo5038804.hoc
            'L_basal' : 1600.,
            'L_soma' : 7.5,
            'L_apical_prox' : 900.,
            'L_apical_dist' : 900.,
            'diam1_basal' : 1.9,
            'diam1_soma' : 7.5,
            'diam1_apical_prox' :  2.75,
            'diam1_apical_dist' : 2.75,
            'diam2_basal' : 1.9,
            'diam2_soma' : 7.5,
            'diam2_apical_prox' : 2.75,
            'diam2_apical_dist' : 2.75,
            'nsec_basal' : 1,
            'nsec_soma' : 1,
            'nsec_apical_prox' : 1,
            'nsec_apical_dist' : 1,
            'syn_types' : ['ampa', 'nmda', 'clopath'],
            'fixnseg':False,        # determine number of segments in cylinder according to d_lambda rule

            # FIXME, must be set so that variable names are unique
            # set recording variables
                # organized a dictionary of dictionaries [attribute name: [variable type: mechanism]
                # note that if the attribute is part of a synapse object, it will accessed differently than a range variable
                    # range variables can be simply accessed by dot notation directly from a given neuron section
                    # synapse attributes need to be accesed from the synapse object stored in cell.syns
            'rec_variables' : 
            [('v','range','v'),
            ('t','range','t'),
            ('gbar','syn','clopath'),
            ('ica_calH','range','calH'),
            ('i', 'syn','nmda')], 

            # choose y variables to plot [varaibles]
            'plot_variables' : ['v','i','ik_kad','i_hd', 'ica_calH', 'ina_na3', 'gbar'],
            # FIXME, should be a list, where you can choose arbitrary combinations of variables 
            # x variables to plot 
            'x_variables':['t'],
            'group_trees':False,

            # synapse activation
            'syn_frac':[],      # fraction of synapses to activate with choose_seg_rand()
            'trial':0,          # count the current trial number
            'trial_id':0,       # a unique identifier for each trial using uuid64
            'w_rand':[],        # choose synapse weights from a random distribution (Bool)
            'w_std' : [],       # standard deviation of weights distribution, if w_rand is True
            'w_mean': .001,       # mean synaptic weight (microsiemens or micro-ohms)
            'trees': [],        # list of subtrees with active synapses [trees]
            'w_list':[],        # nested list of weights, determined by set_weights().  Weights correspond to segments indexed in seg_idx.  Organized as [tree][section][segment]
            'sec_list':[],      # list of active sections with repeats, each entry corresponds to the section for a given segment in seg_list.  [tree][section number]
            'seg_list':[],      # list of active segments, corresponding to sections in sec_list {tree}[segment number]
            'sec_idx': [],      # list of active sections, without repeats. Indeces in the list correspond to indeces in seg_idx {tree}[section number]
            'seg_idx':[],       # nested list of active segments {tree}[section index][segment number]
            'seg_dist' : {},    # distance of each segment from soma {tree}[section index][segment number]

            # extracellular field stimualation
            'field_angle': 0,   # angle relative to principle cell axis in radians 
            'field':[-20,0,20], # list of stimulation intensities in V/m, negative = cathodal, postivie = anodal
            'field_names':['cathodal','control','anodal'], 
            'field_colors':['b','k','r'],    # plot colors correesponding to entries in field
            'field_on':0,      # stimulation onset time in (ms)
            'field_off': 70,    # stimulation offset time in (ms)
            'dt' : .025,        # integration timestep (ms)
            'warmup': 30,       # simulation warmup time (ms)
            'tstop' : 70,       # simulation duration (ms)

            # bipolar stimulation parameters
            'bursts':1,         # bipolar stimulus bursts
            'pulses':4,         # pulses per bursts 
            'pulse_freq':100,   # pulse frequency within burst (Hz)
            'burst_freq':5,     # burst frequency (Hz)
            'noise' : 0,        # noise in input arrival (see NetCon documentation)

            # branch sequence parameters
            'num_sec':1,
            'seg_L' : 4.,
            'seg_spacing':20,
            'max_seg':[],
            'branch':False,
            'full_path':False,
            'branch_distance':[],
            'branch_seg_distance':[],
            'sequence_delay': 0,
            'sequence_direction':'in',

            # clopath synapse parameters
            'clopath_delay_steps': 1,
            'clopath_A_m':3E-5, # depression magnitude parameter (mV^-1)
            'clopath_tetam':-70,#-41, # depression threshold (mV)
            'clopath_tetap':-65,#-38, # potentiation threshold (mV)
            'clopath_tau_r':8,#-38, # potentiation threshold (mV)
            'clopath_tau_0':30,#-38, # potentiation threshold (mV)
            'clopath_tau_y': 5, # time constant (ms) for low pass filter post membrane potential for potentiation
            'clopath_A_p': 38E-5, # amplitude for potentiation (mV^-2)


            # ampa synapse parameters
            'tau1_ampa' : 0.2,  # rise time constant (ms)
            'tau2_ampa' : 2,    # decay time constant   (ms)
            'i_ampa' : 0.18,    # default peak ampa current in uS

            # facilitation depression parameters for AMPA from Varela et al. 1997
            # fit to experimental theta burst and 20 Hz tetanus traces
            'f_ampa':5.,
            'tau_F_ampa':94.,
            'd1_ampa':.45,
            'tau_D1_ampa':540.,
            'd2_ampa':.12,
            'tau_D2_ampa':45.,
            'd3_ampa':.98,
            'tau_D3_ampa':120000.,

            # nmda synapse parameters
            'tau1_nmda' : 1,    # rise time constant (ms)
            'tau2_nmda' : 50,   # decay time constant (ms)

            
            # Parameters from Migliore 2005 (signal propogation in oblique dendrites)
            # conductances reported as (nS/um2) in paper, but need to be in (mho/cm2)
            # conversion 10,000*(pS/um2) = 10*(nS/um2) = (mho/cm2) = 1000.*(mS/cm2)
            # *** units in paper are a typo, values are already reported in (mho/cm2) ***
            'Vrest' : -65.,             # leak potential (mV)
            'v_init' :-71.4,  # initialized membrane voltage, based on steady state membrane potential at soma
            'gna' :  1.*0.025,#.025,                # peak sodium conductance (mho/cm2)
            'dgna' : 0.,#-.000025,          # change in sodium conductance with distance (ohm/cm2/um) from Kim 2015
            'ena' : 55.,                    # sodium reversal potential (mV)
            'gna_inact': 0., # sodium slow inactivation factor (1=none, 0=max inactivation)
            'AXONM' : 10.,              # multiplicative factor for axonal conductance to generate axon potentials in AIS
            'SOMAM':1.5,
            'gkdr' : 1.*0.01,#0.01,             # delayed rectifier potassium peak conductance (mho/cm2)
            'ek' : -90.,                    # potassium reversal potential
            'celsius' : 35.0,               # temperature (degrees C)
            'KMULT' :  1.*0.03,#0.03,           # multiplicative factor for distal A-type potassium conductances
            'KMULTP' : 1.*.03,#0.03,                # multiplicative factor for proximal A-type potassium conductances
            'ghd' : 1.*0.00005,#0.0001,         # peak h-current conductance (mho/cm2)
            'gcalbar': 0.*.00125 ,          # L-type calcium conductance from Kim et al. 2015 (mho/cm2)
            'ehd' : -30.,                   # h-current reversal potential (mV)
            'kl_hd' : -6.,#-8.,
            'vhalfl_hd_prox' : -82.,#-73,           # activation threshold for proximal h current (mV)
            'vhalfl_hd_dist' : -90.,#-81,           # activation threshold for distal h-current (mV)
            'vhalfl_kad' : -56.,#-56.,          # inactivation threshold for distal a-type current (mV)
            'vhalfl_kap' : -56.,#-56.,          # inactivation threshold for proximal a-type current (mV)
            'vhalfn_kad' : -1.,#-1.,            # activation threshold for distal a-type urrent (mV)
            'vhalfn_kap' : -1.,#-1.,            # activation threshold for proximal a-type current (mV)
            'RaAll' : 150.,             # axial resistance, all compartments (ohm*cm)
            'RaAx' : 50.,                   # axial resistance, axon (ohm*cm)                   
            'RmAll' : 28000.,           # specific membrane resistance (ohm/cm2)
            'Cm' : 1.,                  # specific membrane capacitance (uf/cm2)
            'ka_grad' : .1,#1.,#1.,#1.,             # slope of a-type potassium channel gradient with distance from soma 
            'ghd_grad' : .1,#1.,#3.,                # slope of h channel gradient with distance from soma
            'ka_cutoff_distance': 350, # distance from soma where ka stops increasing (um)
            'ghd_cutoff_distance': 350, # distance from soma where Ih stops increasing (um)
        }
    
    def _build_mech_df(self, **kwargs):
        '''
        '''
        # self.p = pd.DataFrame()

        # # simulation
        # self.p['cell']=kwargs['cell_list']
        # self.p['trials']=1
        # self.p['dt']=0.05
        # self.p['tstop']=200
        # self.p['celsius']=35
        # # self.p['rec_variables'] = 'v'

        # # self.p['inhibitory_transmitters']=None
        # # self.p['excitatory_transmitters']=None
        # self.p.at[:,'inhibitory_transmitters'] = [[['gaba'] for i in range(len(self.p))]]
        # self.p.at[:,'excitatory_transmitters'] = [[['acetylcholine','glutamate'] for i in range(len(self.p))]]
        # # self.p.at[:,'inhibitory_transmitters']=[['gaba']]
        # # self.p.at[:,'excitatory_transmitters']=[['acetylcholine','glutamate',]]
        # # (value, type, mechanism name, units)
        # df.append({'param':'L','value':5,'type':'cable','mech_name':None,'units':'um'})
        # df.append({'param':'diam','value':5,'type':'cable','mech_name':None,'units':'um'})
        # df.append({'param':'v_init','value':-65.,'type':'hoc','mech_name':None,'units':'mV'})
        # df.append({'param':'cai','value':0.,'type':'ion','mech_name':None,'units':'mM'})
        # df.append({'param':'cao','value':2.,'type':'ion','mech_name':None,'units':'mM'})
        # df.append({'param':'cm','value':1.,'type':'cable','mech_name':None,'units':'uf'})
        # df.append({'param':'cm','value':1.,'type':'cable','mech_name':None,'units':'uf'})

        mech_list =[
        ('L', 5,     'cable',    None,   'um') ,
        ('diam', 5,     'cable',    None,   'um'),
        ('v_init', -65.,  'hoc',      None,   'mV'),
        ('cai', 0.0,   'ion',      None,   'mM'),
        ('cao', 2.0,   'ion',      None,   'mM'),
        ('cm', 1.0,   'cable',    None,   'uf'),
        ('Ra', 100.0, 'cable',    None,   'ohm*cm'),
        ('gmax_Leak', 5.0E-6, 'range',   'Leak', 'S/cm2'), #(S/cm2)
        ('e_Leak', -65.0, 'range',    'Leak', 'mV'),
        ('gmax_k_slow', 2.*0.003, 'range', 'k_slow','S/cm2'),
        ('ek', -60.0, 'ion',      None,   'mV'),
        ('gmax_k_fast', 7.1164395E-5, 'range', 'k_fast', 'S/cm2'),
        ('gmax_ca_boyle', 2*0.003, 'range',  'ca_boyle', 'S/cm2'),
        ('eca', 40.0, 'ion',       None,   'mV'),
        ('Exp2Syn_tau1', 1,     'synapse', 'Exp2Syn', 'ms'),
        ('Exp2Syn_tau2', 5,     'synapse', 'Exp2Syn', 'ms'),
        ('Exp2Syn_e', 0,     'synapse', 'Exp2Syn', 'mV'),
        ('Exp2Syn_i', 20,    'synapse', 'Exp2Syn','nA'),
        ('Exp2Syn_v', 'netcon_presynaptic',    'synapse', 'Exp2Syn','nA'),

        ('GradedSyn_e', 0, 'synapse', 'GradedSyn', 'mV'),
        # self.p['GradedSyn_e_inh']=-60
        ('GradedSyn_vslope', 4., 'synapse', 'GradedSyn', 'mV'),
        ('GradedSyn_vmid', 20., 'synapse', 'GradedSyn', 'mV'),
        ('GradedSyn_gbar', .0001, 'synapse', 'GradedSyn', 'uS'),
        ('GradedSyn_tau', 4., 'synapse', 'GradedSyn', 'ms'),
        ('GradedSyn_vpre', 'pointer_presynaptic', 'synapse', 'GradedSyn', 'ms'),
        ('neuron_to_neuron_elec_syn_conductance', 100E-12, 'synapse', 'neuron_to_neuron_elec_syn', 'uS') ,
        ]

        mech_df_temp = self._mech_list_to_df(mech_list)
        self.mech_df = pd.DataFrame()
        if 'cell_list' in kwargs:
            for cell in kwargs['cell_list']:
                df_temp = copy.deepcopy(mech_df_temp)
                df_temp['cell']=cell
                self.mech_df = self.mech_df.append(df_temp, ignore_index=True)
        else:
            mech_df_temp['cell']=None
            self.mech_df = mech_df_temp
    
    def _mech_list_to_df(self, mech_list, order=['param', 'val', 'mech_type', 'mech_name', 'units', 'origin', 'min_loc', 'max_loc', 'slope'], **kwargs):
        '''
        '''
        mech_dict = {}
        for i, key in enumerate(order):
            mech_dict[key]=[]
            for mech in mech_list:
                val = None
                if i<len(mech):
                    mech_dict[key].append(mech[i])
                else:
                    mech_dict[key].append(None)
                    
        df = pd.DataFrame(mech_dict)
        return df

    def _set_sequence_delays(self, syn_idx, delay):
        ''' set delays for synapses specified by syn_idx
        ==Args==
        -syn_idx : list of activated synapse locations
        -delay : delay between synapses in syn_idx
                -same delay is applied to each synapse
                -FIXME add ability to specify different delays for each synapse
        ==Out==
        -delays : list of delays (from start of simulation) for input to each synapse in syn_idx
        ==Updates==
        ==Comments==
        
        '''
        delays=[]
        for seg_i, seg in enumerate(syn_idx):
            seg_delay = seg_i*delay
            delays.append(seg_i*delay)

        return delays
    
    def _choose_seg_rand(self, p, p_path, syns, replace=False):
        """ choose random segments to activate, given a subtree and distance from soma requirement

            ==Args==
            -p :  global parameter dictionary
                        -object containing distance from soma for each segment as [tree][section][segment]
            -p_path  : synaptic pathway parameter dictionary
                        -trees: subtrees to choose synapses from
                        -syn_frac: fraction of available synapses to choose from
                                -this is overwritten if a specific number of synapses is given
                                -available synapses are determined based on the specified trees and the distance requirement
                        -syn_num: number of synapses to activate (overrides syn_frac)
                        -syn_dist: distance requirement as [min distance, max distance]
                            -can be given as a nested list with multiple distance requirements.  each distance will have the same number of synapses
                            -for more agile control, create multiple paths
            -syns: synapse structure containing hoc synapse mechanism objects as {tree}[section][segment]{synapse type}
            -replace: boolean. True=randomly choose synapses with replacement, e.g. the same synapse can be selected multiple times

            ==Out==
            -syn_idx_unique : unique list of synapse locations to be activated as [(tree, section, segment)]
            -syn_counts_unique :  list of number of occurences of the synapses specified by syn_idx_unique

            ==Updates==
            ==Comments==
        """
        
        print 'selecting synapses to activate'
        
        trees = p_path['trees']
        syn_frac=p_path['syn_frac']
        seg_dist=p['seg_dist']
        syn_num=p_path['syn_num']
        distance=p_path['syn_dist']

        # list of selected synapses as [synapse number](tree, section, segment)
        segs_chosen = []

        # list all segments as [(tree, section, segment)] 
        segs_all = [(tree_key, sec_i, seg_i) for tree_key, tree in syns.iteritems() for sec_i,sec in enumerate(tree) for seg_i,seg in enumerate(tree[sec_i]) if tree_key in trees]

        # there are multiple distance requirements
        if len(distance)>0 and isinstance(distance[0],list):
            print 'distance:',distance

            # for each distance requirement
            for distance_i, distances in enumerate(distance):
                
                # all segments that fit the current distance requirement
                segs_all_dist = [seg for seg_i, seg in enumerate(segs_all) if seg_dist[seg[0]][seg[1]][seg[2]]>distances[0] and seg_dist[seg[0]][seg[1]][seg[2]]<distances[1]] 

                # if different synapse numbers are provided for each distance bin
                if isinstance(syn_num,list) and len(syn_num)>0:
                    
                    # choose segments to activate
                    segs_choose = np.random.choice(len(segs_all_dist), int(syn_num[distance_i]), replace=replace)

                # if a single scalar is given
                elif syn_num:
                    print 'syn_num:', int(syn_num)
                    print 'available segments:',len(segs_all_dist)
                    # choose segments to activate
                    segs_choose = np.random.choice(len(segs_all_dist), int(syn_num), replace=replace)
                
                # if no synapse number is given
                else:
                    # choose segments to activate
                    segs_choose = np.random.choice(len(segs_all_dist), int(syn_frac*len(segs_all_dist)), replace=replace)

                segs_chosen += [segs_all_dist[a] for a in segs_choose]

        # if only one distacne requirement is given
        elif len(distance) > 0:

            # print 'distance requirement'
             # all segments that fit the current distance requirement
            segs_all_dist = [seg for seg_i, seg in enumerate(segs_all) if seg_dist[seg[0]][seg[1]][seg[2]]>distance[0] and seg_dist[seg[0]][seg[1]][seg[2]]<distance[1]]
            print 'available segments:',len(segs_all_dist)

            # if synapse number is given
            if syn_num:
                print 'synapses selected:', int(syn_num)
                # choose segments to activate
                segs_choose = np.random.choice(len(segs_all_dist), int(syn_num), replace=replace)

            else:
                # choose segments to activate
                segs_choose = np.random.choice(len(segs_all_dist), int(syn_frac*len(segs_all_dist)), replace=replace)

            segs_chosen += [segs_all_dist[a] for a in segs_choose]

        # if no distance requirement given
        else:
            if syn_num:
                print 'syn_num:', int(syn_num)
                # choose segments to activate
                segs_choose = np.random.choice(len(segs_all), int(syn_num), replace=replace)

            else:
                # choose segments to activate
                segs_choose = np.random.choice(len(segs_all), int(syn_frac*len(segs_all)), replace=replace)

            segs_chosen += [segs_all_dist[a] for a in segs_choose]

        # list of selected synapse locations (contains repeats)
        syn_idx = segs_chosen

        # list of number of occurences for each synapse in syn_idx
        syn_counts = [syn_idx.count(temp) for temp in syn_idx ]

        # unique list of synapse locations.  note that order may be different from syn_idx
        syn_idx_unique = list(set(syn_idx))

        # list of number of occurences for each synapse in syn_idx_unique
        syn_counts_unique = [syn_counts[syn_idx.index(temp)] for temp in syn_idx_unique]

        # get unique list of synapses and keep track of count

        return syn_idx_unique, syn_counts_unique
    
    def _set_weights_normal(self, p_path):
        """
        sets weights for synapses specified by p_path['syn_idx']

            ==Args==
            -p_path : synaptic pathway parameter dictionary
                    -syn_idx: list of unique synapse locations as[(tree, section, segment)]
                    -syn_counts: list of occurences for each synapse in syn_idx. w_mean will be multiplied by the number of occurances
                    -w_mean: mean synaptic weight in uS
                    -w_std: standard deviation for synaptic weights
                    -w_rand: Boolean. False: all weights are set to w_mean. True: weights are drawn from normal distribution specified by w_mean and w_std

            ==Out==
            -w_idx : list of weights (uS) for synapses specified by syn_idx


            ==Updates==
            ==Comments==
        """
        syn_idx=p_path['syn_idx']
        syn_counts = p_path['syn_counts']
        w_mean=p_path['w_mean']
        w_std=p_path['w_std']
        w_rand=p_path['w_rand']
        w_idx=[]

        # syn_unique = list(set(syn_idx))
        for seg_i, seg in enumerate(syn_idx):


            # tree, sec_num, seg_num = seg
            repeats = syn_counts[seg_i]

            if w_rand:

                w_idx.append(np.random.normal(repeats*w_mean,w_std))

            else:
                w_idx.append(repeats*w_mean)

        return w_idx

    def _seg_distance(self, cell, cell_type='milstein'):
        """ calculate distance from soma of each segment and store in parameter dictionary

        ==Args==
        -geo  : geometry structure as geo[tree][section][segment]

        ==Out==
        -p['seg_dist']  : structure containing the distance of each segment from the soma, organized in the same way as geo: p['seg_dist'][tree][section][segment distance]

        ==Updates==
        -'seg_dist' is added to p

        ==Comments==
        """

        self.p['seg_dist']={}
        
        # iterate over trees
        for tree_key,tree in cell.geo.iteritems():

            # add dimension for sections
            self.p['seg_dist'][tree_key]=[]
            
            # iterate over sections
            for sec_i,sec in enumerate(tree):
                
                # add dimension for segments
                self.p['seg_dist'][tree_key].append([])
                
                # iterate over segments
                for seg_i,seg in enumerate(sec):
                    
                    # calculate and store distance from soma and store 
                    distance =  h.distance(seg.x, sec=sec)
                    self.p['seg_dist'][tree_key][sec_i].append(distance)

        return self.p['seg_dist']
    
    def _create_morpho(self, geo):
        """ create structure that stores morphology information for plotting with brian2 morphology

        each segment in morpho contains a tuple with seven entries
        (unique_segment_index, name, x, y, z, diam, unique_parent segment_index)

        root segment has index 0, with parent segment index -1
        """

        # initialize morpho structure with same dimensions as geo structure
        morpho = {}
        # iterate over trees
        for tree_key, tree in geo.iteritems():
            morpho[tree_key]=[]
            # iterate over sections
            for sec_i, sec in enumerate(tree):
                morpho[tree_key].append([])
                # iterate over segments
                for seg_i in enumerate(sec):
                    morpho[tree_key][sec_i].append([])

        # find root of cell 
        for tree_key, tree in geo.iteritems():
            for sec_i, sec in enumerate(tree):
                sref = h.SectionRef(sec=sec)
                root = sref.root
                break

        # create new secton list
        nrn_sec_list = h.SectionList()
        # add all seection to list, starting from root
        nrn_sec_list.wholetree()

        # copy nrn section list as a python list
        sec_list = []
        for sec_i_temp, sec_temp in enumerate(nrn_sec_list):
            sec_list.append(sec_temp)

        # nested list for storing segment objects [section_number][segment_number]
        seg_list= []
        # nested list for storing segment indices [section number][segment number]
        seg_list_idx = []
        # nested list for storing index of parent segment, matches seg_list_idx dimesions, [section_number][segment_number]
        parent_list_idx = []
        # keep track of total segment number
        idx = -1
        # iterate through sections in list
        for sec_i, sec in enumerate(sec_list):
            # keep track of the root section
            is_root=False


            # add section dimension to each list
            seg_list.append([])
            seg_list_idx.append([])
            parent_list_idx.append([])

            # reference for current section
            sec_ref =  h.SectionRef(sec=sec)
            
            # find parent section index
            if sec_ref.has_parent():
                parent_sec = sec_ref.parent
                parent_sec_i = [i for i, val in enumerate(sec_list) if parent_sec == val][0]
            else:
                parent_sec_i=-1
                is_root = True

            # iterate through segments in the current section
            for seg_i, seg in enumerate(sec):
                # add to total segments counter and store in lists
                idx+=1
                # copy index count to prevent overwrite during loop
                idx_count = copy.copy(idx)
                # add segment object and index to corresponding list
                seg_list[sec_i].append(seg)
                seg_list_idx[sec_i].append(idx_count)

                # if current segment is not the first in its section 
                if seg_i>0:
                    # set parent to previous segemnt in the section
                    parent_seg_idx = seg_list_idx[sec_i][seg_i-1]
                # else if it is the first segment 
                elif seg_i==0:
                    # if it is the root segment
                    if is_root:
                        parent_seg_idx=-1
                    else:
                        # set parent to the last segment in the parent section
                        parent_seg_idx = seg_list_idx[parent_sec_i][-1]

                # add to list of all parent segments
                parent_list_idx.append(parent_seg_idx)

                # find the current segment in geo structure
                # iterate through geo structure until you find matching segment
                for tree_key_local, tree_local in geo.iteritems():
                    for sec_i_local, sec_local in enumerate(tree_local):
                        for seg_i_local, seg_local in enumerate(sec_local):

                            # if section name and segment index match
                            if (sec.name() == sec_local.name()) and (seg_i == seg_i_local):
                                # segment diameter
                                diam = seg_local.diam
                                # segment xyz coordinates
                                xyz = self._seg_location(sec_local)
                                x = xyz[0][seg_i_local]
                                y = xyz[1][seg_i_local]
                                z = xyz[2][seg_i_local]

                                # segment name
                                name = tree_key_local + '_'+ str(sec_i_local) + '_'  +str(seg_i_local)
                                # create 7-tuple
                                morph_tuple = (idx_count, name, x, y, z, diam, parent_seg_idx)
                                # store in morphology structure
                                morpho[tree_key_local][sec_i_local][seg_i_local] = morph_tuple
            
                                break 
                        else:
                            continue
                        break
                    else:
                        continue
                    break
                else:
                    continue

        return morpho

    def _create_loc_list(self, geo):
        '''
        '''
        locations=[]
        for tree_key, tree in geo.iteritems():
            for sec_i, sec in enumerate(tree):
                for seg_i, seg in enumerate(sec):
                    location = (tree_key, sec_i, seg_i)
                    locations.append(location)

        return locations

    def _seg_location(self, sec):
        """ given a neuron section, output the 3d coordinates of each segment in the section

        ouput is a nested list as [xyz dimension][segment number], with x,y, z dimensions listed in that order

        """
        # number of 3d points in section
        tol =.001
        n3d = int( h.n3d( sec=sec))
        
        # preallocate 3d coordinates
        x = [None]*n3d
        y = [None]*n3d
        z = [None]*n3d
        position_3d =  [None]*n3d
                       
        # loop over 3d coordinates in each section
        for i in range(n3d):
            # retrieve x,y,z
            x[i] = h.x3d(i, sec=sec)
            y[i] = h.y3d(i, sec=sec)
            z[i] = h.z3d(i, sec=sec)

            # calculate total distance of each 3d point from start of section
            if i is 0:
                position_3d[i] = 0
            else:
                position_3d[i] = position_3d[i-1] + np.sqrt((x[i]-x[i-1])**2 + (y[i]-y[i-1])**2 + (z[i]-z[i-1])**2)
        
        seg_x = []
        seg_y = []
        seg_z = []
        for seg_i,seg in enumerate(sec):
                # relative position within section (0-1)
                seg_pos = seg.x            
                
                # segment distance along section in 3D
                seg_dist = seg_pos*position_3d[-1]

                # find first 3D coordinate that contains the segment
                node_i = [dist_i for dist_i,dist in enumerate(position_3d) if dist >= seg_dist]
                
                # if segement occurs exactly at a node set its location to the node location
                if abs(position_3d[node_i[0]] - seg_dist) < tol:
                    seg_x.append( x[ node_i[ 0]])
                    seg_y.append( z[ node_i[ 0]])
                    seg_z.append( z[ node_i[ 0]])

                # otherwise if segment falls between two coordinates, interpolate to get location
                # FIXME clean up
                else:
                    pt1 = position_3d[ node_i[0]-1]
                    pt2 = position_3d[ node_i[0]]
                    scale = (seg_dist-pt1) / (pt2-pt1)
                    interpx = x[ node_i[0]-1] + scale*( x[ node_i[0]] - x[ node_i[0]-1])
                    interpy = y[ node_i[0]-1] + scale*( y[ node_i[0]] - y[ node_i[0]-1])
                    interpz = z[ node_i[0]-1] + scale*( z[ node_i[0]] - z[ node_i[0]-1])
                    seg_x.append( interpx)
                    seg_y.append( interpy)
                    seg_z.append( interpz)
        return [seg_x, seg_y, seg_z]
    
    def _load_fd_parameters(self, p, filename):
        # load parameters
            #````````````````
        with open(filename, 'rb') as pkl_file:
            param_obj = pickle.load(pkl_file)

        params = param_obj.x

        p['f_ampa'] = params[0]
        p['tau_F_ampa'] = 1E3*params[1]
        p['d1_ampa'] = params[2]
        p['tau_D1_ampa'] = 1E3*params[3]
        p['d2_ampa'] = params[4]
        p['tau_D2_ampa'] = 1E3*params[5]
        p['d3_ampa'] = params[6]
        p['tau_D3_ampa'] = 1E3*params[7]

        return p