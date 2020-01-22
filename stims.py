"""
implement extracellular stimulation or presynaptic input patterns

"""
import numpy as np
from neuron import h
from sklearn.decomposition import PCA

# FIXME
def _get_principal_axis(geo, include_axon=False):
    '''
    '''
    # arrange segment locations as array with segments x spatial dimensions
    seg_locs, seg_locs_list = _get_seg_locations(geo)
    if not include_axon:
        seg_locs_list = [seg for seg in seg_locs_list if seg[0]!='axon']

    x_vals = np.array([seg[3] for seg in seg_locs_list])
    y_vals = np.array([seg[4] for seg in seg_locs_list])
    z_vals = np.array([seg[5] for seg in seg_locs_list])

    X = np.vstack([x_vals, y_vals, z_vals]).T

    pca = PCA(n_components=1)
    pca.fit(X)
    components = pca.components

    # option to include or ignore axon

def _get_seg_locations(geo):
    '''
    '''
    seg_locs={}
    seg_locs_list=[]
    for tree_key, tree in geo.iteritems():
        seg_locs[tree_key]=[]
        for sec_i, sec in enumerate(tree):
            seg_locs[tree_key].append([])
            sec_xyz = _get_seg_location_from_sec(sec)
            for seg_i, seg in enumerate(sec):
                seg_x = sec_xyz[0][seg_i]
                seg_y = sec_xyz[1][seg_i]
                seg_z = sec_xyz[2][seg_i]
                seg_locs[tree_key][sec_i].append((seg_x, seg_y, seg_z))

                seg_locs_list.append((tree_key, sec_i, seg_i, seg_x, seg_y,seg_z))

    return seg_locs, seg_locs_list

def _get_seg_location_from_sec(self, sec):
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

def _set_branch_nseg(geo, sec_idx, seg_L):
    """ set number of segments for branch that was selected to activate synapses
    
    arguments:

    """

    # iterate over trees in section list
    for tree_key, tree in sec_idx.iteritems():

        for sec_i, sec in enumerate(tree):

            section = geo[tree_key][sec]

            # get section length
            sec_L = section.L
            print 'section length', section.L

            # determine number of segments
            n_seg = int(np.ceil(sec_L/seg_L))

            # # check that number of segments is odd
            if n_seg % 2 != 0:
                n_seg+=1

            # # set number of segments
            section.nseg = n_seg
            print 'nseg', section.nseg
    return geo

def _update_synapses_after_nseg(p, geo, syns, sec_idx):
    ''' if nseg is updated on a branch, synapses need to be reinserted
    ==Args==
    -p : parameter dictionary
    -geo : cell geometry object containing hoc sections and segments
    -syns : cell synapse structure containing synapses as hoc point processes
    -sec_idx : indicates which sections to update synapses on as sec_idx[tree_key][section index]
    ==Return==
    -syns : updated synapse object
    ==Update==
    -syns : reinserts synaptic point processes and updates syns to match the current number of segments in the section specified by sec_idx
    ==Comment==
    '''
    # iterate over trees
    for tree_key, tree in sec_idx.iteritems():
        # iterate over sections
        for sec_i, sec in enumerate(tree):
            # re-initialize the list of synapses
            syns[tree_key][sec]=[]
            # iterate over segments
            for seg_i,seg in enumerate(geo[tree_key][sec]):
                # add dictionary key for each type of synapse
                syns[tree_key][sec].append({'ampa':[],
                'nmda':[],
                'clopath':[]})
                # iterate over types of synapses
                for syn_key,syn in syns[tree_key][sec][seg_i].iteritems():
                    # ampa
                    if syn_key is 'ampa':
                        
                        # adapting exponential synapse based on model in Varela et al. 1997
                        # insert point process and store in syns
                        syns[tree_key][sec][seg_i][syn_key] = h.FDSExp2Syn_D3(geo[tree_key][sec](seg.x))
                        # update parameters
                        syns[tree_key][sec][seg_i][syn_key].f = p['f_ampa']
                        syns[tree_key][sec][seg_i][syn_key].tau_F = p['tau_F_ampa']
                        syns[tree_key][sec][seg_i][syn_key].d1 = p['d1_ampa']
                        syns[tree_key][sec][seg_i][syn_key].tau_D1 = p['tau_D1_ampa']
                        syns[tree_key][sec][seg_i][syn_key].d2 = p['d2_ampa']
                        syns[tree_key][sec][seg_i][syn_key].tau_D2 = p['tau_D2_ampa']
                        syns[tree_key][sec][seg_i][syn_key].d3 = p['d3_ampa']
                        syns[tree_key][sec][seg_i][syn_key].tau_D3 = p['tau_D3_ampa']
                    # nmda
                    elif syn_key is 'nmda':
                        # insert point process and store in syns
                        syns[tree_key][sec][seg_i][syn_key]= h.Exp2SynNMDA(geo[tree_key][sec](seg.x))
                        # update parameters
                        syns[tree_key][sec][seg_i][syn_key].tau1 = p['tau1_nmda']
                        syns[tree_key][sec][seg_i][syn_key].tau2 = p['tau2_nmda']
                        # print syn

                    elif syn_key is 'clopath':
                        # clopath
                        # insert point process and store in syn
                        syns[tree_key][sec][seg_i][syn_key] = h.STDPSynCCNon(geo[tree_key][sec](seg.x))
    return syns

def _get_terminal_branches(geo):
    ''' get sections that correspond to dendritic branch terminals
    ==Args==
    -geo : cell geometry object containing hoc sections and segments
    ==Return==
    -terminal_branches : sections that correspond to terminal branches as terminal_branches[tree_key][section number]
    ==Update==
    ==Comments==
    -terminal branhces are determined based on the section having no children.  This may not be the most general way of finding branches
    '''
    terminal_branches = {}
    # iterate over trees
    for tree_key, tree in geo.iteritems():
        terminal_branches[tree_key]=[]
        # iterate over sections
        for sec_i, sec in enumerate(tree):
            # create hoc sectioref object
            secref = h.SectionRef(sec=sec)
            # if section has no children, it is assumed to be a terminal 
            if secref.nchild()==0:
                # add to list of terminal branches in the current tree
                terminal_branches[tree_key].append(sec_i)
    return terminal_branches

def _choose_seg_from_branch(geo, sec_idx):
    ''' create list of active synapses based on segments in branch
    ==Args==
    -geo : cell geometry object 
    -sec_idx : specifies which sections to activate synapses on. sec_idx[tree_key][section number]
    ==Return==
    -syn_idx_unique : unique list of segment locations to acitvate synapses on as syn_idx_unique[(tree_key, sec_num, seg_num)]
    -syn_count_unique : list with elements corresponding to syn_idx_unique, which specifies the number of synapses to activate at the specified segment
    ==Update==
    ==Comments==
    '''
    syn_idx_unique = []
    syn_count_unique=[]
    for tree_key, tree in sec_idx.iteritems():
        for sec_i, sec in enumerate(tree):
            for seg_i, seg in enumerate(geo[tree_key][sec]):
                location = (tree_key, sec, seg_i)
                syn_idx_unique.append(location)
                syn_count_unique.append(1)

    return syn_idx_unique, syn_count_unique

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

class IPS:
    '''
    Intersectional Pulsed Stimulation from Buzsaki lab
    '''

class Interference:
    '''
    Temporal interference as in Grossman 2017
    '''

class ACS:
    """
    assumes somatodendritic axis is aligned vertically, y is positive for apical dendrites, y is negative for basal dendrites
    """
    def __init__(self, p, cell=0):
        self.insert_e(cell=cell, p=p)

    def insert_e(self, p, cell=0):
        '''
        '''

        intensity=p['ac_field']
        field_angle=p['ac_field_angle']
        field_on=p['ac_field_on'] 
        field_off=p['ac_field_off'] 
        dt=p['dt']
        tstop=p['tstop']
        freq=p['ac_field_freq']
        phase=p['ac_field_phase']

        t = np.arange(0,tstop+dt, dt)
        
        if cell == 0:
            cell=[]
            for sec in h.allsec():
                cell.append(sec)

        # structure to store location and e_extracellular for each segment.  Organized as ['dimension'][section number][segment number]
        location = {'x':[], 'y':[],'z':[],'e':[]}
        
        # vectors for time-dependent control of DCS onset using vector.play
        self.e_vec = []
        self.t_vec = []
        # loop over sections
        for sec_i,sec in enumerate(cell):
            # add dimension for each section to play vectors
            self.e_vec.append([])
            self.t_vec.append([])
            
            # add list for each section to store data
            for dim_key,dim in location.iteritems():
                dim.append([])

            # insert extracellular mechanism
            sec.insert('extracellular')

            # number of 3d points in section
            n3d = int(h.n3d(sec=sec))
            
            # xyz locations of segments
            xyz = self.seg_location(sec)
            # print sec.name(), xyz

            # iterate over segments
            for seg_i,seg in enumerate(sec):
                e=[]
                # xyz location of each segment
                seg_x = xyz[0][seg_i]
                seg_y = xyz[1][seg_i]
                seg_z = xyz[2][seg_i]

                # angle of segment from somato-dendritic axis (neglect z axis)   
                if seg_y == 0:
                    angle = 0
                elif np.isnan(seg_x/float(seg_y)):
                    angle = 0
                else:
                    angle = np.arctan(seg_x/seg_y)
                
                # if y location is negative shift phase by pi
                # FIXME
                if seg_y < -0.001:
                    angle = angle+np.pi
                
                # absolute distance of segment from (0,0) in um
                mag = np.sqrt(seg_x**2 + seg_y**2)
                
                # angle relative to electric field vector, zero angle means along somato-dendritic axis
                angle_field = angle + field_angle
                # convert um to mm
                conversion = .001 

                # calculate extracellular potential
                e = conversion*intensity*mag*np.cos(angle_field)

                e_vec = e*np.sin(phase + 2*np.pi*(freq/1000)*t)

                e_vec[t<field_on]=0.
                e_vec[t>=field_off]=0.

                self.e_vec[sec_i].append(h.Vector(t.shape[0]))
                self.t_vec[sec_i].append(h.Vector(t.shape[0]))

                for e_vec_i, e_vec_val in enumerate(e_vec):
                    self.e_vec[sec_i][seg_i].x[e_vec_i] = e_vec_val
                    self.t_vec[sec_i][seg_i].x[e_vec_i] = t[e_vec_i]

                # apply play method to control DCS field during simulation
                self.e_vec[sec_i][seg_i].play(seg._ref_e_extracellular, self.t_vec[sec_i][seg_i])

    def seg_location(self, sec):
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

class DCS:
    """
    assumes somatodendritic axis is aligned vertically, y is positive for apical dendrites, y is negative for basal dendrites
    """
    def __init__(self, cell=0, intensity=0, field_angle=0, field_on=0, field_off=0):
        '''
        '''
        self.insert_e(cell=cell, intensity=intensity, field_angle=field_angle, field_on=field_on, field_off=field_off)

    def insert_e(self, cell=0, intensity=0, field_angle=0, field_on=0, field_off=1000):
        
        if cell == 0:
            cell=[]
            for sec in h.allsec():
                cell.append(sec)

        # structure to store location and e_extracellular for each segment.  Organized as ['dimension'][section number][segment number]
        location = {'x':[], 'y':[],'z':[],'e':[]}
        
        # vectors for time-dependent control of DCS onset using vector.play
        self.e_vec = []
        self.t_vec = []
        # loop over sections
        for sec_i,sec in enumerate(cell):
            # add dimension for each section to play vectors
            self.e_vec.append([])
            self.t_vec.append([])
            
            # add list for each section to store data
            for dim_key,dim in location.iteritems():
                dim.append([])

            # insert extracellular mechanism
            sec.insert('extracellular')

            # number of 3d points in section
            n3d = int(h.n3d(sec=sec))
            
            # xyz locations of segments
            xyz = self.seg_location(sec)
            # print sec.name(), xyz

            # iterate over segments
            for seg_i,seg in enumerate(sec):
                e=[]
                # xyz location of each segment
                seg_x = xyz[0][seg_i]
                seg_y = xyz[1][seg_i]
                seg_z = xyz[2][seg_i]

                # angle of segment from somato-dendritic axis (neglect z axis)   
                if seg_y == 0:
                    angle = 0
                elif np.isnan(seg_x/float(seg_y)):
                    angle = 0
                else:
                    angle = np.arctan(seg_x/seg_y)
                # if y location is negative shift phase by pi
                # FIXME
                if seg_y < -0.001:
                    angle = angle+np.pi
                
                # absolute distance of segment from (0,0) in um
                mag = np.sqrt(seg_x**2 + seg_y**2)
                
                # angle relative to electric field vector, zero angle means along somato-dendritic axis
                angle_field = angle + field_angle
                # convert um to mm
                conversion = .001 

                # calculate extracellular potential
                e = conversion*intensity*mag*np.cos(angle_field)

                # FIXME
                # create vectors for play mechanism
                self.e_vec[sec_i].append(h.Vector(3))
                self.e_vec[sec_i][seg_i].x[0] = 0
                self.e_vec[sec_i][seg_i].x[1] = e
                self.e_vec[sec_i][seg_i].x[2] = 0

                self.t_vec[sec_i].append(h.Vector(3))
                self.t_vec[sec_i][seg_i].x[0] = 0
                self.t_vec[sec_i][seg_i].x[1] = field_on
                self.t_vec[sec_i][seg_i].x[2] = field_off

                # apply play method to control DCS field during simulation
                self.e_vec[sec_i][seg_i].play(seg._ref_e_extracellular, self.t_vec[sec_i][seg_i])


                # print sec.name(), seg_i, e
                # insert calculated extracellular potential in mV
                # seg.e_extracellular = e

                # print sec.name(), seg_y, seg.x, seg.e_extracellular

    def seg_location(self, sec):
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

class Uncage:
    """ simulate glutamate uncaging experiments
    """
    def __init__(self):
        pass

    def _bursts(self, p_path):
        """ create NetStim objects for activating synapses in bursts

        ==Args==
        entries in p_path dictionary
        -syn_idx  : list of synapses to be activated [synapse number](tree, sectin, segment)
        -delays   : list of onset delays for simulating sequences of inputs.  indices match syn_idx
        -bursts   : number of bursts 
        -pulses   : number of pulses per burst
        -pulse_freq : frequency of pulses within burst (Hz)
        -burst_freq : frwquency of bursts (Hz)
        -warmup     : simulation warmup time, to allow for gate variables to reach steady state.  applied equally to all synapses (ms)
        -noise     : fractional noise for stimulus intervals, ranging from 0 to 1, with 0=no noise, 1=max noise
                    -defined as:    
                        -t0 = (1-noise)*interval
                        -for t<t0, p(t=input)=0
                        -for t>=t0, p(t=input)=(1/(interval - t0))*exp(-(t-t0)/(interval - t0))
                            -this is a negative exponential distribution, which has maximum at t=t0 and decays exponentially afterwards

        ==Out==
        -stim   : list of NetStim objects organized as [synapse number][burst number]

        ==Updates==
        -burst parameters and stim are stored as attributes of the Uncage class

        ==Comments==


        """
        print 'creating NetStim objects'

        warmup=p_path['warmup']
        bursts=p_path['bursts']
        pulses=p_path['pulses']
        pulse_freq=p_path['pulse_freq']
        burst_freq=p_path['burst_freq']
        noise=p_path['noise']
        syn_idx=p_path['syn_idx']
        delays=p_path['sequence_delays']


        fs = 1000. # convert time to ms
        warmup = warmup   # warm up time (ms)
        stim  = [] # list of stim objects [synapse number][burst number]

        print 'syn_idx:', syn_idx
        print delays
        # iterate over synapses in syn_idx
        for seg_i, seg in enumerate(syn_idx):

            # add entry for each synapse
            stim.append([])

            # iterate over bursts
            for burst in range(int(bursts)):

                # create NetStim object for each burst
                stim[seg_i].append(h.NetStim())

                # if a list of onset delays is provided
                if delays:
                    # onset time
                    stim[seg_i][burst].start = warmup + delays[seg_i] + burst*fs/burst_freq
                else:
                    stim[seg_i][burst].start = warmup + burst*fs/burst_freq

                # update interval, noise, number of pulses
                stim[seg_i][burst].interval = fs/pulse_freq
                stim[seg_i][burst].noise  = noise 
                stim[seg_i][burst].number = pulses

        return stim

    def _connect_synapses(self, p_path, stim, syns):
        '''
        ==Args==
        -p : parameter dictionary
                -must contain p['syn_idx'], a zipped list of tuples

        -stim  : nested list of stim objects
                -[segment][burst], with the first index matching entries in p['syn_idx']

        -syns  : structure containing synapse objects
                -syns{tree}[section][segment]{synapse type}
                -synapse types include ampa, nmda, clopath
        ==Out==
        -nc    : structure containing hoc NetCon objects
                -nc[segment]{synapse type}[burst number]
                -segment index matches p['syn_idx']
        ==Updates==
        -NetCon objects are created

        ==Comments==
        '''

        # for storing NetCon objects nc[segment]{synapse type}[burst number]
        nc = []

        # syn_idx_unique=list(set(p_path['syn_idx']))
        # iterate over segments to be activated
        print p_path['syn_idx']
        for seg_i, seg in enumerate(p_path['syn_idx']):

            # add dimension for segments
            nc.append({})

            # get segment location info
            tree, sec_num, seg_num = seg
            print seg
            # iterate over synapse types
            for syntype_key,syntype in syns[tree][sec_num][seg_num].iteritems():

                # create list for storing NetCon objects, one for each burst in stim
                nc[seg_i][syntype_key]=[]

                # iterate over bursts in stim
                for burst_i, burst in enumerate(stim[seg_i]):

                    # create netcon object
                    netcon = h.NetCon(burst, syntype, 0, 0, p_path['w_idx'][seg_i])


                    # add netcon object to list
                    nc[seg_i][syntype_key].append(netcon)

        return nc

    def _poisson_times(self, start, dt, T, mod_freq, mod_amp, mean_rate, refractory=0):
        '''generate poisson spike times with sinusoidal varying rate
        '''
        # time vector
        time=np.arange(start,T,dt)
        mod_freq=mod_freq/1000.
        # sinusoidally varying rate (convert to time to ms)
        rates=.001*(mean_rate + mod_amp*np.sin(2*np.pi*mod_freq*time))
        # set negative rates to zero
        rates=np.clip(rates, a_min=0, a_max=None)
        # maximum rate (peak of the sinusoid)
        rate_max=np.max(rates)
        # list to store spike times
        S=[]
        # keep track times and spike count
        t=0
        k=0
        # generate times until max time is reached
        while t<=T:
            # generate random numbers
            r=np.random.uniform()
            s=np.random.uniform()
            # check that neuron is not refractory
            if np.abs(np.log(r)/rate_max)>refractory:
                # spike time assuming maximum rate
                t = t-np.log(r)/rate_max
                # corresponding index in time and rate vectors
                i = np.abs(time-t).argmin()
                # corresponding rate
                rate=rates[i]
                # apply rate criteria
                if s<=rate/rate_max:
                    # store spike time
                    S.append(t)
        # convert to array
        S=np.array(S)

        return S

    def _poisson(self, p_path):
        '''
        '''
        print 'creating NetStim objects'

        start=p_path['warmup']
        dt=p_path['dt']
        T=p_path['tstop']
        if 'mod_freq' in p_path:
            mod_freq=p_path['mod_freq']
        else:
            mod_freq=0
        if 'mod_amp' in p_path:
            mod_amp=p_path['mod_amp']
        else:
            mod_amp=0
        if 'mean_rate' in p_path:
            mean_rate=p_path['mean_rate']
        else:
            mean_rate=0


        # fs = 1000. # convert time to ms
        # warmup = warmup   # warm up time (ms)
        stim  = [] # list of stim objects [synapse number][burst number]

        # iterate over synapses in syn_idx
        for seg_i, seg in enumerate(p_path['syn_idx']):
            stim.append([])
            # create NetStim object for each burst
            stim[seg_i].append(h.VecStim())
            # get poisson times
            times = self._poisson_times(start=start, dt=dt, T=T, mod_freq=mod_freq, mod_amp=mod_amp, mean_rate=mean_rate)
            # create neuron vector
            vec = h.Vector(times)
            # play times into vecstim
            stim[seg_i][0].play(vec)

        return stim

class Intracellular:
    '''
    '''
    def __init__(self):
        '''
        '''
        pass

    def _insert_IClamp(self, cell, location, p):
        '''
        '''
        tree,sec,seg = location
        seg_loc = float(seg+1)/(cell.geo[tree][sec].nseg+1)
        self.iclamp=[]
        for delay in p['iclamp_delays']:
            self.iclamp.append(h.IClamp(cell.geo[tree][sec](seg_loc)))
            self.iclamp[-1].delay = delay
            self.iclamp[-1].dur = p['iclamp_dur']
            self.iclamp[-1].amp = p['iclamp_amp']
        # self.iclamp = h.IClamp(cell.geo[tree][sec](seg_loc))
        # self.iclamp.delay = p['iclamp_delay']
        # self.iclamp.dur = p['iclamp_dur']
        # self.iclamp.amp = p['iclamp_amp']

        return self.iclamp