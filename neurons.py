"""
create cells and activate subsets of synapses
"""
# imports
import numpy as np
import matplotlib.pyplot as plt
from neuron import h
import stims
import param
import run_control


class CellMigliore2005:
    """ pyramidal neuron based on Migliore et al. 2005

    An instance of this object will creates a cell (hoc objects) at the top level of the hoc interpreter using the hoc files in _init_geometry.  The .geo attribute contains a python mapping to these hoc objects.  The geo object is organized as geo['section tree'][section](segment location)

    the syns attribute creates a container for synapse objects that are added to each segment in the hoc cell.  syns is organized as syns['section tree']['synapse type'][section][segment number]
    
    """
    def __init__(self,p):
        '''
        '''
        pass

    def geometry(self,p):
        '''create cell geometry at hoc top level
        ===Args===
        -p : parameter dictionary
        
        ===Out===
        -geo : geometry structure containing hoc section objects
                -geo{tree}[section][segment]
        -syns : structure containing synapse objects
                -syns{tree}[section][segment]{synapse type}
                -synapse types include ampa, nmda, clopath
        
        ===Updates===
        -hoc objects are created based on the geometry in the loaded hoc file
        
        ===Comments===
        '''
        print 'loading cell geometry:', self.__class__.__name__
        # load cell geometry into hoc interpreter
        h.load_file('geo5038804.hoc')  
        # h.load_file('geoc62564.hoc')
        # set discretization based on dlambda rule (set dlambda in hoc file) 
        h.load_file('fixnseg.hoc')      
        # dictionary for storing geometry ['tree'][sec](seg location)
        self.geo = {}
        # dictionary for storing synapse objects ['tree']['type'][sec][seg]
        self.syns = {}
        # add section trees to geometry dictionary
        self.geo['soma'] = h.soma
        self.geo['axon'] =  h.axon
        self.geo['basal'] = h.dendrite
        self.geo['apical_trunk'] = h.user5
        self.geo['apical_tuft'] = h.apical_dendrite
        
        # set temperature in hoc
        h.celsius = p['celsius']
        # set soma as origin for distance measurements
        h.distance(sec = self.geo['soma'][0])

        return self.geo, self.syns

    def mechanisms(self,p):
        """ insert membrane mechanisms into cell geometry
        
        ==Args==
        -p : parameter dictionary

        ==Out==
        -geo : geometry structure containing hoc section objects
                -geo{tree}[section][segment].mechanism
        -syns : structure of containing synapse mechanisms
                -syns{tree}[section][segment]{synapse type}
        ==Updates==
        -range mechanisms and their parameters are updated according to the parameters in p
        ==Comments==
        self.syns is updated to store an object for each synapse.  It is organized as ['tree']['synapse type'][section][segment].  Note that the last index will depend on how the cell is discretized as the number segments changes in each sections 

        the parameters for each membrane mechanism  are store in a dictionary called p.  See the param module for details.
        """
        print 'loading cell range mechanisms'
        
        # loop over trees
        for tree_key,tree in self.geo.iteritems():
            
            # list to store synapse mechanisms
            self.syns[tree_key] = []

            # loop over sections in tree
            for sec_i,sec in enumerate(tree):
                
                # add dimension for each section
                self.syns[tree_key].append([])

                # common passive biophysics for all sections
                sec.insert('pas')
                # passive conductance (S/cm2)
                sec.g_pas = 1/p['RmAll']            
                # leak reversal potential (mV)  
                sec.e_pas = p['Vrest']              
                # specific capacitance (uf/cm2)
                sec.cm = p['Cm']            
                # axial resistance (ohm cm)         
                sec.Ra = p['RaAll'] 
                                        
                # axon active bipophysics
                if tree_key == 'axon':
                    # voltage gated sodium
                    sec.insert('nax')                       
                    sec.gbar_nax = p['gna']*p['AXONM']
                    # print 'axon sodium conductance:', sec.gbar_nax*10000
                    # delayed rectifier potassium
                    sec.insert('kdr')                       
                    sec.gkdrbar_kdr = p['gkdr']
                    # a-type potassium
                    sec.insert('kap')                       
                    sec.gkabar_kap = p['KMULTP']
                    sec.vhalfl_kap = p['vhalfl_kap']
                    sec.vhalfn_kap = p['vhalfn_kap']
                    # sodium reversal potential 
                    sec.ena = p['ena']      
                    # potassium reversal potential 
                    sec.ek = p['ek']
                    sec.Ra = p['RaAx']


                    for seg_i, seg in enumerate(sec):
                        self.syns[tree_key][sec_i].append({})
                    
                # soma active biophysics
                elif tree_key == 'soma':

                    # voltage gated sodium
                    sec.insert('na3')
                    sec.gbar_na3 = p['gna']*p['SOMAM']
                    sec.ar_na3 = p['gna_inact']
                    # print 'soma sodium conductance:', sec.gbar_na3*10000
                    # h-current         
                    sec.insert('hd')
                    sec.ghdbar_hd = p['ghd']                
                    sec.vhalfl_hd = p['vhalfl_hd_prox']
                    sec.kl_hd = p['kl_hd']
                    sec.ehd_hd = p['ehd']       


                    # delayed rectifier potassium       
                    sec.insert('kdr')
                    sec.gkdrbar_kdr = p['gkdr'] 
                    # a-type potassium      
                    sec.insert('kap')
                    sec.gkabar_kap = p['KMULTP']
                    sec.vhalfl_kap = p['vhalfl_kap']
                    sec.vhalfn_kap = p['vhalfn_kap']

                    sec.insert('calH')
                    sec.gcalbar_calH = p['gcalbar']
                    # sodium reversal potential 
                    sec.ena = p['ena']      
                    # potassium reversal potential 
                    sec.ek = p['ek']    

                    for seg_i,seg in enumerate(sec):
                        self.syns[tree_key][sec_i].append({})

                    
                # dendrites active biophysics
                elif ((tree_key == 'basal') or 
                (tree_key == 'apical_trunk') or 
                (tree_key == 'apical_tuft')):
                    # h-current
                    sec.insert('hd')
                    sec.ghdbar_hd = p['ghd']
                    sec.kl_hd = p['kl_hd']
                    sec.ehd_hd = p['ehd']
                    
                    # voltage gated sodium      
                    sec.insert('na3')
                    sec.gbar_na3 = p['gna']
                    sec.ar_na3 = p['gna_inact']

                    # delayed rectifier potassium   
                    sec.insert('kdr')
                    sec.gkdrbar_kdr = p['gkdr'] 
                    # a-type potassium proximal
                    sec.insert('kap')
                    sec.gkabar_kap = 0  
                    # a-type potassium distal       
                    sec.insert('kad')
                    sec.gkabar_kad = 0  

                    # L-type calcium channel
                    sec.insert('calH')
                    sec.gcalbar_calH = p['gcalbar']

                    # sodium reversal potential 
                    sec.ena = p['ena']
                    # potassium reversal potential
                    sec.ek = p['ek']        

                    # mechanisms that vary with distance from soma
                    # loop over segments
                    for seg_i,seg in enumerate(sec):
                        
                        # print seg_i
                        self.syns[tree_key][sec_i].append({'ampa':[],
                        'nmda':[],
                        'clopath':[]})

                        for syn_key,syn in self.syns[tree_key][sec_i][seg_i].iteritems():
                            
                            if syn_key is 'ampa':
                                
                                # adapting exponential synapse based on model in Varela et al. 1997
                                self.syns[tree_key][sec_i][seg_i][syn_key] = h.FDSExp2Syn_D3(sec(seg.x))
                                self.syns[tree_key][sec_i][seg_i][syn_key].f = p['f_ampa']
                                self.syns[tree_key][sec_i][seg_i][syn_key].tau_F = p['tau_F_ampa']
                                self.syns[tree_key][sec_i][seg_i][syn_key].d1 = p['d1_ampa']
                                self.syns[tree_key][sec_i][seg_i][syn_key].tau_D1 = p['tau_D1_ampa']
                                self.syns[tree_key][sec_i][seg_i][syn_key].d2 = p['d2_ampa']
                                self.syns[tree_key][sec_i][seg_i][syn_key].tau_D2 = p['tau_D2_ampa']
                                self.syns[tree_key][sec_i][seg_i][syn_key].d3 = p['d3_ampa']
                                self.syns[tree_key][sec_i][seg_i][syn_key].tau_D3 = p['tau_D3_ampa']

                                # regular double exponential synapse
                                # self.syns[tree_key][sec_i][seg_i][syn_key] = h.Exp2Syn(sec(seg.x))
                                # self.syns[tree_key][sec_i][seg_i][syn_key].tau1 = p['tau1_ampa']
                                # self.syns[tree_key][sec_i][seg_i][syn_key].tau2 = p['tau2_ampa']
                                # self.syns[tree_key][sec_i][seg_i][syn_key].i = p['i_ampa']
                                # print syn

                            elif syn_key is 'nmda':
                                # print syn_key
                                self.syns[tree_key][sec_i][seg_i][syn_key]= h.Exp2SynNMDA(sec(seg.x))
                                self.syns[tree_key][sec_i][seg_i][syn_key].tau1 = p['tau1_nmda']
                                self.syns[tree_key][sec_i][seg_i][syn_key].tau2 = p['tau2_nmda']
                                # print syn

                            elif syn_key is 'clopath':
                                # print syn_key
                                self.syns[tree_key][sec_i][seg_i][syn_key] = h.STDPSynCCNon(sec(seg.x))

                        # distance from soma
                        seg_dist = h.distance(seg.x,sec=sec)
                        
                        # sodium
                        if abs(p['dgna']*seg_dist)<p['gna']:
                            seg.gbar_na3 = p['gna'] + p['dgna']*seg_dist
                        else:
                            seg.gbar_na3 = 0.
                        
                        # h current
                        if seg_dist < p['ghd_cutoff_distance']:
                            seg.ghdbar_hd = p['ghd']*(1+p['ghd_grad']*seg_dist/100.)
                        else:
                            seg.ghdbar_hd = p['ghd']*(1+p['ghd_grad']*p['ghd_cutoff_distance']/100.)

                        
                        # A-type potassium
                        if seg_dist > 100.: # distal
                            seg.vhalfl_hd = p['vhalfl_hd_dist']
                            seg.vhalfl_kad = p['vhalfl_kad']
                            seg.vhalfn_kad = p['vhalfn_kad']
                            if seg_dist < p['ka_cutoff_distance']:
                                seg.gkabar_kad = p['KMULT']*(1+p['ka_grad']*seg_dist/100.)
                            else:
                                seg.gkabar_kad = p['KMULT']*(1+p['ka_grad']*p['ka_cutoff_distance']/100.)
                        else:   # proximal
                            seg.vhalfl_hd = p['vhalfl_hd_prox']
                            seg.vhalfl_kap = p['vhalfl_kap']
                            seg.vhalfn_kap = p['vhalfn_kap']
                            seg.gkabar_kap = p['KMULTP']*(1+p['ka_grad']*seg_dist/100.)

                        # print tree_key, sec_i, seg_i, dir(seg.calH) 

class PyramidalCylinder:
    """ 4 compartment pyramidal cell with HH dynamics
    """
    def __init__(self, p):
        self.geometry(p)
        self.mechanisms(p)

    def geometry(self, p):
        """
        areas determined from cell geo5038804 from Migliore 2005
        basal: 19966.3598 um2
        apical: 23700.5664916 (tuft) + 11461.4440485 (trunk) = 35162.0105401 um2
        soma: 176.290723263 um2
        axon: 305.021056197 um2
        """
        # list sections
        trees = ['basal', 'soma', 'apical_prox', 'apical_dist']
        areas = [19966.36, 176.29, 35162.01/2., 35162.01/2.]
        
        # store geometry
        self.geo = {}
        # store synapses
        self.syns = {}
        # create sections
        for tree_i, tree in enumerate(trees):
            self.geo[tree] = []
            self.syns[tree] = [] 

            for sec_i in range(p['nsec_'+tree]):
                self.geo[tree].append( h.Section( name=tree))

                # diameter basd on area of full morphology
                diam1 = p['diam1_'+tree]
                diam2 = p['diam2_'+tree] 

                if tree=='soma':    
                    # create 3d specification, with cell arranged vertically
                    h.pt3dadd(0, 0, 0, diam1, sec=self.geo[tree][sec_i])
                    h.pt3dadd(0, p['L_'+tree], 0, diam2, sec=self.geo[tree][sec_i])

                if tree=='basal':
                    h.pt3dadd(0, 0, 0, diam1, sec=self.geo[tree][sec_i])
                    h.pt3dadd(0, -p['L_'+tree], 0, diam2, sec=self.geo[tree][sec_i])

                if tree=='apical_prox':
                    h.pt3dadd(0, p['L_soma'], 0, diam1, sec=self.geo[tree][sec_i])
                    h.pt3dadd(0, p['L_soma']+p['L_'+tree], 0, diam2, sec=self.geo[tree][sec_i])

                if tree=='apical_dist':
                    h.pt3dadd(0, p['L_soma']+p['L_apical_prox'], 0, diam1, sec=self.geo[tree][sec_i])
                    h.pt3dadd(0, p['L_soma']+p['L_apical_prox']+p['L_'+tree], 0, diam2, sec=self.geo[tree][sec_i])

                # add list to store synapses for each section
                self.syns[tree].append([])

                # insert passive mechanism
                self.geo[tree][sec_i].insert('pas')
                # passive conductance (S/cm2)
                self.geo[tree][sec_i].g_pas = 1./p['RmAll']         
                # leak reversal potential (mV)  
                self.geo[tree][sec_i].e_pas = p['Vrest']                
                # specific capacitance (uf/cm2)
                self.geo[tree][sec_i].cm = p['Cm']          
                # axial resistance (ohm cm)         
                self.geo[tree][sec_i].Ra = 1.*p['RaAll'] 

                self.geo[tree][sec_i].L = p['L_'+tree]
                # self.geo[tree][sec_i].diam = p['diam_'+tree]

                self.geo[tree][sec_i].nseg=p['nseg']

        self.geo['basal'][0].connect(self.geo['soma'][0](0),0)
        self.geo['apical_prox'][0].connect(self.geo['soma'][0](1),0)
        self.geo['apical_dist'][0].connect(self.geo['apical_prox'][0](1),0)

        if p['fixnseg']==True:
            h.xopen('fixnseg.hoc')

        # set temperature in hoc
        h.celsius = p['celsius']
        # set soma as origin for distance measurements
        h.distance(sec=self.geo['soma'][0])

    def rotate(self, theta):
        """Rotate the cell about the Z axis.
        """
        for sec in h.allsec():
            for i in range(int(h.n3d(sec=sec))):
                x = h.x3d(i, sec=sec)
                y = h.y3d(i, sec=sec)
                c = np.cos(theta)
                s = np.sin(theta)
                xprime = x * c - y * s
                yprime = x * s + y * c
                h.pt3dchange(i, xprime, yprime, h.z3d(i, sec=sec), h.diam3d(i, sec=sec), sec=sec)

    def mechanisms(self, p):

        for tree_key, tree in self.geo.iteritems():
            for sec_i, sec in enumerate(tree): 

                if tree_key == 'soma':
                    # voltage gated sodium
                    sec.insert('na3')
                    sec.gbar_na3 = p['gna']*p['AXONM']
                    # h-current         
                    sec.insert('hd')
                    sec.ghdbar_hd = p['ghd']                
                    sec.vhalfl_hd = p['vhalfl_hd_prox']
                    sec.kl_hd = p['kl_hd']
                    sec.ehd_hd =  p['ehd']

                    # delayed rectifier potassium       
                    sec.insert('kdr')
                    sec.gkdrbar_kdr = p['gkdr'] 

                    # a-type potassium      
                    sec.insert('kap')
                    sec.gkabar_kap = p['KMULTP']
                    sec.vhalfl_kap = p['vhalfl_kap']
                    sec.vhalfn_kap = p['vhalfn_kap']

                    # L-type calcium channel
                    sec.insert('calH')
                    sec.gcalbar_calH = p['gcalbar']

                    # sodium reversal potential 
                    sec.ena = p['ena']      
                    # potassium reversal potential 
                    sec.ek = p['ek']

                    for seg_i,seg in enumerate(sec):
                        self.syns[tree_key][sec_i].append({})
                
                elif ((tree_key == 'basal') or 
                (tree_key == 'apical_prox') or 
                (tree_key == 'apical_dist')):

                    # h-current
                    sec.insert('hd')
                    sec.ghdbar_hd = p['ghd']
                    sec.kl_hd = p['kl_hd']
                    sec.ehd_hd =  p['ehd']
                    
                    # voltage gated sodium      
                    sec.insert('na3')
                    sec.gbar_na3 = p['gna']

                    # delayed rectifier potassium   
                    sec.insert('kdr')
                    sec.gkdrbar_kdr = p['gkdr'] 
                    # a-type potassium proximal
                    sec.insert('kap')
                    sec.gkabar_kap = 0  
                    # a-type potassium distal       
                    sec.insert('kad')
                    sec.gkabar_kad = 0  

                    # L-type calcium channel
                    sec.insert('calH')
                    sec.gcalbar_calH = p['gcalbar']

                    # sodium reversal potential 
                    sec.ena = p['ena']
                    # potassium reversal potential
                    sec.ek = p['ek']        

                    # mechanisms that vary with distance from soma
                    # loop over segments
                    for seg_i,seg in enumerate(sec):

                        # add segment dimension to syns structure
                        self.syns[tree_key][sec_i].append([])

                        # distance from soma
                        seg_dist = h.distance(seg.x, sec=sec)
                        
                        # sodium
                        seg.gbar_na3 = p['gna'] + p['dgna']*seg_dist
                        # print seg_dist, seg.gbar_na3
                        
                        # h current
                        seg.ghdbar_hd = p['ghd']*(1+p['ghd_grad']*(seg_dist/100.)/(p['L_apical_prox']/200.))

                        # h current
                        if seg_dist < p['ghd_cutoff_distance']*(p['L_apical_prox']/200.):
                            seg.ghdbar_hd = p['ghd']*(1+p['ghd_grad']*(seg_dist/100.)/(p['L_apical_prox']/200.))
                        else:
                            seg.ghdbar_hd = p['ghd']*(1+p['ghd_grad']*(p['ghd_cutoff_distance']/100.)/(p['L_apical_prox']/200.))
                        
                        # A-type potassium
                        if seg_dist > 100.*(p['L_apical_prox']/200.): # distal
                            seg.vhalfl_hd = p['vhalfl_hd_dist']
                            seg.vhalfl_kad = p['vhalfl_kad']
                            seg.vhalfn_kad = p['vhalfn_kad']
                            seg.gkabar_kad = p['KMULT']*(1+p['ka_grad']*(seg_dist/100.)/(p['L_apical_prox']/200.))
                            if seg_dist < p['ka_cutoff_distance']*(p['L_apical_prox']/200.):
                                seg.gkabar_kad = p['KMULT']*(1+p['ka_grad']*(seg_dist/100.)/(p['L_apical_prox']/200.))
                            else:
                                seg.gkabar_kad = p['KMULT']*(1+p['ka_grad']*(p['ka_cutoff_distance']/100.)/(p['L_apical_prox']/200.))
                        else:   # proximal
                            seg.vhalfl_hd = p['vhalfl_hd_prox']
                            seg.vhalfl_kap = p['vhalfl_kap']
                            seg.vhalfn_kap = p['vhalfn_kap']
                            seg.gkabar_kap = p['KMULTP']*(1+p['ka_grad']*(seg_dist/100.)/(p['L_apical_prox']/200.))

                        self.syns[tree_key][sec_i][seg_i] = {
                        'ampa':[],
                        'nmda':[],
                        'clopath':[],
                        }

                        for syn_key,syn in self.syns[tree_key][sec_i][seg_i].iteritems():
                            if syn_key is 'ampa':
                                # Regular ampa synapse
                                # self.syns[tree_key][sec_i][seg_i][syn_key] = h.Exp2Syn(sec(seg.x))
                                # self.syns[tree_key][sec_i][seg_i][syn_key].tau1 = p['tau1_ampa']
                                # self.syns[tree_key][sec_i][seg_i][syn_key].tau2 = p['tau2_ampa']
                                # self.syns[tree_key][sec_i][seg_i][syn_key].i = p['i_ampa']

                                # FD adapting exponential synapse based on model in Varela et al. 1997
                                self.syns[tree_key][sec_i][seg_i][syn_key] = h.FDSExp2Syn_D3(sec(seg.x))
                                self.syns[tree_key][sec_i][seg_i][syn_key].f = p['f_ampa']
                                self.syns[tree_key][sec_i][seg_i][syn_key].tau_F = p['tau_F_ampa']
                                self.syns[tree_key][sec_i][seg_i][syn_key].d1 = p['d1_ampa']
                                self.syns[tree_key][sec_i][seg_i][syn_key].tau_D1 = p['tau_D1_ampa']
                                self.syns[tree_key][sec_i][seg_i][syn_key].d2 = p['d2_ampa']
                                self.syns[tree_key][sec_i][seg_i][syn_key].tau_D2 = p['tau_D2_ampa']
                                self.syns[tree_key][sec_i][seg_i][syn_key].d3 = p['d3_ampa']
                                self.syns[tree_key][sec_i][seg_i][syn_key].tau_D3 = p['tau_D3_ampa']

                            elif syn_key is 'nmda':
                                # print syn_key
                                self.syns[tree_key][sec_i][seg_i][syn_key]= h.Exp2SynNMDA(sec(seg.x))
                                self.syns[tree_key][sec_i][seg_i][syn_key].tau1 = p['tau1_nmda']
                                self.syns[tree_key][sec_i][seg_i][syn_key].tau2 = p['tau2_nmda']
                                # print syn

                            elif syn_key is 'clopath':
                                # print syn_key
                                self.syns[tree_key][sec_i][seg_i][syn_key] = h.STDPSynCCNon(sec(seg.x))                        