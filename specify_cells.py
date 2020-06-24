__author__ = 'Aaron D. Milstein'
# Includes modification of an early version of SWC_neuron.py by Daniele Linaro.
# Includes an extension of BtMorph, created by Ben Torben-Nielsen and modified by Daniele Linaro.
from function_lib import *
from neuron import h  # must be found in system $PYTHONPATH
import pprint
import btmorph  # must be found in system $PYTHONPATH
from itertools import izip
import hoc2swc # converts hoc morphology to swc file based on 3d geometry
import os

# SWC files must use this nonstandard convention to exploit trunk and tuft categorization
swc_types = [soma_type, axon_type, basal_type, apical_type, trunk_type, tuft_type] = [1, 2, 3, 4, 5, 6]
sec_types = ['soma', 'axon_hill', 'ais', 'axon', 'basal', 'trunk', 'apical', 'tuft', 'spine_neck', 'spine_head']

verbose = 0  # Turn on for text reporting during model initialization and simulation
gid_max = 0  # Every new HocCell will receive a global identifier for network simulation, and increment this counter

#-------Wrapper for converting SWC --> BtMorph --> Skeleton morphology in NEURON hoc------------------------
morph_dir='morphologies/'
mech_filename='043016 Type A - km2_NMDA_KIN5_Pr'
morph_filename='EB1-early-bifurcation.swc'
# morph_filename='geo5038804.hoc'
# cell = HocCell(morph_filename, mech_filename)

class HocCell(object):
    def __init__(self, morph_filename=None, mech_filename=None, mech_dict=None):
        """
        :param morph_filename: str : path to .swc file containing morphology
        :param mech_filename: str : path to .pkl file specifying cable parameters and membrane mechanisms
        """
        global gid_max
        self._gid = gid_max
        gid_max += 1
        self.tree = btmorph.STree2()  # Builds a simple tree to store nodes of type 'SHocNode'
        self.index = 0  # Keep track of number of nodes
        self._node_dict = {'soma': [], 'axon': [], 'basal': [], 'trunk': [], 'apical': [], 'tuft': [], 'spine': []}
        # load mech_dict:Refer to function_lib for structure description
        if mech_dict is None:
            self.mech_dict = self.load_mech_dict(mech_filename)  
        else:
            self.mech_dict = mech_dict
            # mechanism dictionary. loads from .pkl or
            # default_mech_dict in function_lib
        if not morph_filename is None:
            self.load_morphology_gk(morph_filename)
            # self.load_morphology_from_swc_gk(morph_filename)
            self.reinit_mechanisms()  # Membrane mechanisms must be reinitialized whenever cable properties (Ra, cm) or
                                      # spatial resolution (nseg) changes.
        self.set_node_names()
        self.update_diameters()
        # self.spike_detector = None
        # if self.axon:
        #     self.init_spike_detector()
        self.random = np.random.RandomState()

    
    def load_morphology_gk(self, morph_filename, morph_dir='morphologies/'):
        '''
        '''
        # check morph file type
        if 'swc' in morph_filename:
            file_type='swc'
        elif 'hoc' in morph_filename:
            file_type='hoc'
        # if hoc morphology
        if file_type=='hoc':
            self.load_morphology_from_hoc_pt3d(morph_filename=morph_filename)
        elif file_type=='swc':
            self.load_morphology_from_swc_gk(morph_filename=morph_filename)
        else:
            raise 'unsupported morphology file: only .hoc or .swc files allowed'

    def load_morphology_from_swc_gk(self, morph_filename):
        """
        This method builds an STree2 comprised of SHocNode nodes associated with hoc sections, connects the hoc
        sections, and initializes various parameters: Ra, cm, L, diam, nseg
        This method implements a standardized soma and axon:
        The soma consists of two cylindrical hoc sections of equal length and diameter, connected (0) to (0).
        The basal dendritic tree is connected to soma[1](1), and the apical tree is connected to soma[0](1).
        The axon is attached to soma[0](0), and consists of three sequential hoc sections:
            1) axon[0] : a tapered cylindrical 'axon hillock' section connected to soma[0](0)
            2) axon[1] : a tapered cylindrical 'axon initial segment' section connected to axon[0](1)
            3) axon[2] : a cylindrical 'axon' section connected to axon[1](1)
        """
        raw_tree = btmorph.STree2()  # import the full tree from an SWC file
        raw_tree.read_SWC_tree_from_file(morph_dir+morph_filename, types=range(10))
        self.raw_tree=raw_tree
        
        # get and set apical trunk nodes if not already specified
        self.trunk_nodes = self.get_apical_trunk_from_raw_tree(raw_tree=self.raw_tree, apical_type=4., tol=1.2)
        self.update_raw_tree_node_types(raw_tree=self.raw_tree, node_list=self.trunk_nodes, new_type=5)

        self.make_soma(raw_tree=self.raw_tree)
        self.make_axon(raw_tree=self.raw_tree)

    def load_morphology_from_hoc_pt3d(self, morph_filename, update_swc_file=True):
        """
        **for btmorph to handle the cell morphology, the soma must have at least 3 nodes**
        This method builds an STree2 comprised of SHocNode nodes associated with hoc sections, connects the hoc
        sections, and initializes various parameters: Ra, cm, L, diam, nseg
        This method implements a standardized soma and axon:
        The soma consists of two cylindrical hoc sections of equal length and diameter, connected (0) to (0).
        The basal dendritic tree is connected to soma[1](1), and the apical tree is connected to soma[0](1).
        The axon is attached to soma[0](0), and consists of three sequential hoc sections:
            1) axon[0] : a tapered cylindrical 'axon hillock' section connected to soma[0](0)
            2) axon[1] : a tapered cylindrical 'axon initial segment' section connected to axon[0](1)
            3) axon[2] : a cylindrical 'axon' section connected to axon[1](1)
        """
        # FIXME update hoc section name to swc section type mapping
        new_map = self.swc_type_from_hoc_section_name
        # hoc2swc.swc_type_from_section_name.__code__ = new_map.__code__ 
        
        # check if corresponding swc file exists
        swc_filename = morph_filename.split('.hoc')[0]+'.swc'
        morph_filenames = os.listdir(morph_dir)
        if swc_filename not in morph_filenames or update_swc_file:
            # load hoc file
            h.load_file(morph_filename) 
            # convert to swc
            hoc2swc.neuron2swc(morph_dir+swc_filename)

        # clear neuron sections
        for sec in h.allsec():
            del sec
        # new btmorph tree
        raw_tree = btmorph.STree2()  # import the full tree from an SWC file
        raw_tree.read_SWC_tree_from_file(morph_dir+swc_filename, types=range(10))
        self.raw_tree=raw_tree
        
        # get and set apical trunk nodes if not already specified
        self.trunk_nodes = self.get_apical_trunk_from_raw_tree(raw_tree=self.raw_tree, apical_type=4., tol=1.2)
        self.update_raw_tree_node_types(raw_tree=self.raw_tree, node_list=self.trunk_nodes, new_type=5)

        self.make_soma(raw_tree=self.raw_tree)
        self.make_axon(raw_tree=self.raw_tree)

    def swc_type_from_hoc_section_name(self, section_name):
        '''
        Returns an integer string of an SWC point type in response to a string name of a NEURON section.

        See column 2 of http://www.neuronland.org/NLMorphologyConverter/MorphologyFormats/SWC/Spec.html

        To map custom section names to parts of SWC cell, override this method. E.g:

        # Create a new name->type map
        def new_map(section_name):
            return "5" if "foo" in section_name else "1"

        # Replace the default map with the new one. Subsequent hoc2swc statements will use the new map.
        from hoc2swc import swc_type_from_section_name
        swc_type_from_section_name.__code__ = new_map.__code__

        :param section_name: name string of a NEURON section
        :return: integer string e.g. "1" or "3" that corresponds to a SWC point type
        '''     
        if 'tuft' in section_name:
            return '5'

        if 'trunk' in section_name:
            return '6'

        if "apic" in section_name:
            return "4"

        if "den" in section_name or 'basal' in section_name:
            return "3"

        if "axon" in section_name or "hillock" in section_name or "initial" in section_name:
            return "2"

        if "soma" in section_name:
            return "1"

        return "5"

    def load_morphology_from_swc(self, morph_filename):
        """
        This method builds an STree2 comprised of SHocNode nodes associated with hoc sections, connects the hoc
        sections, and initializes various parameters: Ra, cm, L, diam, nseg
        This method implements a standardized soma and axon:
        The soma consists of two cylindrical hoc sections of equal length and diameter, connected (0) to (0).
        The basal dendritic tree is connected to soma[1](1), and the apical tree is connected to soma[0](1).
        The axon is attached to soma[0](0), and consists of three sequential hoc sections:
            1) axon[0] : a tapered cylindrical 'axon hillock' section connected to soma[0](0)
            2) axon[1] : a tapered cylindrical 'axon initial segment' section connected to axon[0](1)
            3) axon[2] : a cylindrical 'axon' section connected to axon[1](1)
        """
        raw_tree = btmorph.STree2()  # import the full tree from an SWC file
        raw_tree.read_SWC_tree_from_file(morph_dir+morph_filename, types=range(10))
        self.raw_tree=raw_tree
        soma_length = 14.
        soma_diam = 9.
        for index in range(2):
            node = self.make_section('soma')
            node.sec.L = soma_length/2.
            node.sec.diam = soma_diam
            self._init_cable(node)  # consults the mech_dict to initialize Ra, cm, and nseg
        self.tree.root = self.soma[0]
        self.soma[1].connect(self.soma[0], 0, 0)
        

        for index in range(3):
            self.make_section('axon')
        self.axon[0].type = 'axon_hill'
        self.axon[0].sec.L = 10.
        self.axon[0].set_diam_bounds(3., 2.)  # stores the diameter boundaries for a tapered cylindrical section
        self.axon[1].type = 'ais'
        self.axon[1].sec.L = 15.
        self.axon[1].set_diam_bounds(2., 0.5)  # (2., 1.)
        self.axon[2].sec.L = 500.
        self.axon[2].sec.diam = 0.5  # 1.
        self.axon[0].connect(self.soma[0], 0, 0)
        self.axon[1].connect(self.axon[0], 1, 0)
        self.axon[2].connect(self.axon[1], 1, 0)
        for node in self.axon:
            self._init_cable(node)
        for child in raw_tree.root.children:
            self.make_skeleton(child, self.tree.root)

    def make_axon(self, raw_tree,):
        '''
        '''
        for index in range(3):
            self.make_section('axon')
        self.axon[0].type = 'axon_hill'
        self.axon[0].sec.L = 10.
        self.axon[0].set_diam_bounds(3., 2.)  # stores the diameter boundaries for a tapered cylindrical section
        self.axon[0].p3d = copy.deepcopy(self.soma[1].p3d)
        self.axon[0].p3d.xyz += [0., -self.axon[0].sec.L, 0.]
        self.axon[1].type = 'ais'
        self.axon[1].sec.L = 15.
        self.axon[1].set_diam_bounds(2., 0.5)  # (2., 1.)
        self.axon[1].p3d = copy.deepcopy(self.axon[0].p3d)
        self.axon[1].p3d.xyz += [0., -self.axon[1].sec.L, 0.]
        self.axon[2].sec.L = 500.
        self.axon[2].sec.diam = 0.5  # 1.
        self.axon[2].p3d = copy.deepcopy(self.axon[1].p3d)
        self.axon[2].p3d.xyz += [0., -self.axon[1].sec.L, 0.]
        self.axon[0].connect(self.soma[1], 1, 0)
        self.axon[1].connect(self.axon[0], 1, 0)
        self.axon[2].connect(self.axon[1], 1, 0)

        for node in self.axon:
            self._init_cable(node)
        for child in raw_tree.root.children:
            self.make_skeleton(child, self.tree.root)

    def make_soma(self, raw_tree, soma_L=14., soma_diam=9, ):
        '''
        '''
        # get root SNode
        root = raw_tree.root
        # get SNode2 soma objects
        soma_SNodes = [node for node in raw_tree.get_nodes() if node.content['p3d'].type==soma_type]
        # check if tree root is in the soma
        if root in soma_SNodes:
            # get root p3d
            p3d = root.content['p3d']
        else:
            p3d = soma_SNodes[0].content['p3d']

        # make two soma sections
        for index in range(2):
            node = self.make_section('soma')
            node.sec.L = soma_L/2.
            node.sec.diam = soma_diam
            node.p3d = copy.deepcopy(p3d)
            print node.p3d.xyz
            if index==0:
                node.p3d.xyz +=[0., soma_L/2., 0.] 
                print node.p3d.xyz
            elif index==1:
                node.p3d.xyz += [0., -soma_L/2., 0.] 
                print node.p3d.xyz
            self._init_cable(node)  # consults the mech_dict to initialize Ra, cm, and nseg
        self.tree.root = self.soma[0]
        self.soma[1].connect(self.soma[0], 0, 0)

    def make_section(self, sec_type,**kwargs):
        """
        Create a new hoc section to associate with this node, and this cell, and store information about it in the
        node's content dictionary.
        :param sec_type: str
        :return node: :class:'SHocNode'
        """
        node = SHocNode(self.index)
        if self.index == 0:
            self.tree.root = node
        self.index += 1
        node.type = sec_type
        if sec_type in ['spine_head', 'spine_neck']:
            self._node_dict['spine'].append(node)
        else:
            self._node_dict[sec_type].append(node)
        node.sec = h.Section(name=node.name, cell=self)
        return node

    def test_sec_properties(self, node=None):
        """
        Used for debugging and validating model specification.
        :param node:
        """
        if node is None:
            node = self.tree.root
            #        node.sec.push()
            #        h.psection()
            #        h.pop_section()
        print '{} [nseg: {}, Ra: {}]'.format(node.name, node.sec.nseg, node.sec.Ra)
            #        h('for (x) print (x), diam(x)', sec=node.sec)
        for child in node.children:
            self.test_sec_properties(child)

    def make_skeleton(self, raw_node, parent, length=0, diams=None, ):
        """
        Following construction of soma and axon nodes of type 'SHocNode' in the tree of type 'STree2', this method
        recursively converts dendritic 'SNode2' nodes into 'SHocNode' nodes, and connects them to the appropriate
        somatic nodes. Skeletonized dendritic nodes have only one hoc section for each stretch of unbranched dendrite,
        with length equal to the sum of the lengths of the converted SNode2 nodes.
        Nodes that taper more than 1 um remain tapered, otherwise they are converted into untapered cylinders with
        diameter equal to the mean diameter of the the converted SNode2 nodes.
        Dendrite types that are pre-categorized as basal, apical, trunk, or tuft in the input .swc file are preserved.
        :param raw_node: :class:'SNode2'
        :param parent: :class:'SHocNode'
        :param length: int or float
        :param diams: None or (list: float)
        """
        global verbose
        dend_types=([basal_type, apical_type, trunk_type, tuft_type], ['basal', 'apical', 'trunk', 'tuft'])
        swc = raw_node.get_content()['p3d']
        swc_type = swc.type
        if swc_type in dend_types[0]:
            diam = swc.radius*2
            length += self.get_node_length_swc(raw_node)
            leaves = len(raw_node.children)
            # create a new node when encountering 1) branch points, 2) terminal ends, or 3) change in swc_type
            if leaves > 1 or leaves == 0 or \
                    (leaves == 1 and not raw_node.children[0].get_content()['p3d'].type == swc_type):
                sec_type = dend_types[1][dend_types[0].index(swc_type)]
                new_node = self.make_section(sec_type)
                new_node.sec.L = length
                new_node.p3d = swc
                if (self.tree.is_root(parent)) and (sec_type == 'basal'):
                    parent = self.soma[1]
                new_node.connect(parent)
                if diams is None:
                    new_node.sec.diam = diam
                    self._init_cable(new_node)
                    if verbose:
                        print '{} [nseg: {}, diam: {}, length: {}, parent: {}]'.format(new_node.name, new_node.sec.nseg,diam, length, new_node.parent.name)
                else:
                    diams.append(diam)
                    if len(diams) > 2:
                        mean = np.mean(diams)
                        stdev = np.std(diams)
                        if stdev*2 > 1:  # If 95% of the values are within 1 um, don't taper
                            new_node.set_diam_bounds(mean+stdev, mean-stdev)
                            self._init_cable(new_node)
                            if verbose:
                                print '{} [nseg: {}, diam: ({}:{}), length: {}, parent: {}]'.format(new_node.name,
                                                new_node.sec.nseg, mean+stdev, mean-stdev, length, new_node.parent.name)
                        else:
                            new_node.sec.diam = mean
                            self._init_cable(new_node)
                            if verbose:
                                print '{} [nseg: {}, diam: {}, length: {}, parent: {}]'.format(new_node.name,
                                                                new_node.sec.nseg, mean, length, new_node.parent.name)
                    elif abs(diams[0]-diams[1]) > 1:
                        new_node.set_diam_bounds(diams[0], diams[1])
                        self._init_cable(new_node)
                        if verbose:
                            print '{} [diam: ({}:{}), length: {}, parent: {}]'.format(new_node.name, new_node.sec.nseg,
                                                                    diams[0], diams[1], length, new_node.parent.name)
                    else:
                        mean = np.mean(diams)
                        new_node.sec.diam = mean
                        self._init_cable(new_node)
                        if verbose:
                            print '{} [nseg: {}, diam: {}, length: {}, parent: {}]'.format(new_node.name,
                                                                new_node.sec.nseg, mean, length, new_node.parent.name)
                # Follow all branches from this fork
                for child in raw_node.children:
                    self.make_skeleton(child, new_node)
            else:  # Follow unbranched dendrite
                if diams is None:
                    diams = [diam]
                else:
                    diams.append(diam)
                self.make_skeleton(raw_node.children[0], parent, length, diams)

    def get_nodes_of_subtype(self, sec_type):
        """
        This method searches the node dictionary for nodes of a given sec_type and returns them in a list. Used during
        specification of membrane mechanisms.
        :param sec_type: str
        :return: list of :class:'SHocNode'
        """
        if sec_type in ['axon_hill', 'ais', 'axon']:
            return [node for node in self.axon if node.type == sec_type]
        elif sec_type in ['spine_head', 'spine_neck']:
            return [node for node in self.spine if node.type == sec_type]
        else:
            return self._node_dict[sec_type]

    def load_mech_dict(self, mech_filename=None):
        """
        This method loads the dictionary specifying membrane mechanism parameters. If a .pkl file is not provided, a
        global variable default_mech_dict from function_lib is used.
        :param mech_filename: str
        """
        if not mech_filename is None:
            return read_from_pkl(data_dir+mech_filename+'.pkl')
        else:
            local_mech_dict = copy.deepcopy(default_mech_dict)
            return local_mech_dict

    def _init_cable(self, node):
        """
        If the mechanism dictionary specifies the cable properties 'Ra' or 'cm', then _modify_mechanism() properly sets
        those parameters, and reinitializes the number of segments per section. To avoid redundancy, this
        method passes _modify_mechanism() a copy of the dictionary with the spatial_res parameter removed, since this is
        consulted in setting nseg. However, if spatial_res is the only parameter being specified, it is passed to
        _modify_mechanism()
        :param node: :class:'SHocNode'
        """
        sec_type = node.type
        if sec_type in self.mech_dict and 'cable' in self.mech_dict[sec_type]:
            mech_content = copy.deepcopy(self.mech_dict[sec_type]['cable'])
            if ('Ra' in mech_content) or ('cm' in mech_content):
                if 'spatial_res' in mech_content:
                    del mech_content['spatial_res']
                self._modify_mechanism(node, 'cable', mech_content)
            elif 'spatial_res' in mech_content:
                self._modify_mechanism(node, 'cable', mech_content)
        else:
            node.init_nseg()
            node.reinit_diam()

    def reinit_mechanisms(self, reset_cable=0):
        """
        Once a mechanism dictionary has been loaded, and a morphology has been specified, this method traverses through
        the tree of SHocNode nodes following order of inheritance and properly sets membrane mechanism parameters,
        including gradients and inheritance of parameters from nodes along the path from root. Since cable parameters
        are set during specification of morphology, it is not necessary to immediately reinitialize these parameters
        again. However, they can be manually reinitialized with the reset_cable flag.
        :param reset_cable: boolean
        """
        for sec_type in sec_types:
            if sec_type in self.mech_dict:
                nodes = self.get_nodes_of_subtype(sec_type)
                self._reinit_mech(nodes, reset_cable)

    def init_synaptic_mechanisms(self):
        """
        Spines and synapses are inserted after loading a morphology and specifying membrane mechanisms. This method can
        be executed after inserting synapses. It traverses the dendritic tree in order of inheritance and just sets
        synaptic mechanism parameters specified in the mechanism dictionary.
        """
        for dend_type in ['soma', 'basal', 'trunk', 'apical', 'tuft']:
            if dend_type in self.mech_dict and 'synapse' in self.mech_dict[dend_type] and \
                                                self.sec_type_has_synapses(dend_type):
                for node in self.get_nodes_of_subtype(dend_type):
                    self._modify_mechanism(node, 'synapse', self.mech_dict[dend_type]['synapse'])

    def node_has_synapses(self, node, syn_type=None):
        """
        Checks if a given node contains synapses, or spines with synapses. Can also check for a synaptic point process
        of a specific type.
        :param node: :class:'SHocNode'
        :param syn_type: str
        :return: boolean
        """
        if [syn for syn in node.synapses if syn_type is None or syn_type in syn._syn]:
            return True
        else:
            for spine in node.spines:
                if [syn for syn in spine.synapses if syn_type is None or syn_type in syn._syn]:
                    return True
        return False

    def sec_type_has_synapses(self, sec_type, syn_type=None):
        """
        Checks if any nodes of a given sec_type contain synapses, or spines with synapses. Can also check for a synaptic
        point process of a specific type.
        :param sec_type: str
        :param syn_type: str
        :return: boolean
        """
        for node in self.get_nodes_of_subtype(sec_type):
            if self.node_has_synapses(node, syn_type):
                return True
        return False

    def _reinit_mech(self, nodes, reset_cable=0):
        """
        Given a list of nodes, this method loops through all the mechanisms specified in the mechanism dictionary for
        the hoc section type of each node and updates their associated parameters. If the reset_cable flag is set to 1,
        cable parameters are modified first, then the parameters for all other mechanisms are reinitialized.
        Parameters for synaptic point processes can also be specified in the mechanism dictionary, so one must use the
        method init_synaptic_mechanisms() after inserting synapses.
        :param nodes: list of :class:'SHocNode'
        :param reset_cable: int or boolean
        """
        for node in nodes:
            sec_type = node.type
            if sec_type in self.mech_dict:
                if ('cable' in self.mech_dict[sec_type]) and reset_cable:  
                    # cable properties must be set first, as they
                    # can change nseg, which will affect insertion
                    # of membrane mechanism gradients
                    self._init_cable(node)
                for mech_name in (mech_name for mech_name in self.mech_dict[sec_type]
                                  if not mech_name in ['cable', 'ions']):
                    self._modify_mechanism(node, mech_name, self.mech_dict[sec_type][mech_name])
                # ion-related parameters do not exist until after membrane mechanisms have been inserted
                if 'ions' in self.mech_dict[sec_type]:
                    self._modify_mechanism(node, 'ions', self.mech_dict[sec_type]['ions'])

    def reinitialize_subset_mechanisms(self, sec_type, mech_name):
        """
        During parameter optimization, it is often convenient to reinitialize all the parameters for a single mechanism
        in a subset of compartments. For example, g_pas in basal dendrites that inherit the value from the soma after
        modifying the value in the soma compartment.
        :param sec_type: str
        :param mech_name: str
        :return:
        """
        if sec_type in self.mech_dict and mech_name in self.mech_dict[sec_type]:
            for node in self.get_nodes_of_subtype(sec_type):
                self._modify_mechanism(node, mech_name, self.mech_dict[sec_type][mech_name])

    def _modify_mechanism(self, node, mech_name, mech_content):
        """
        This method loops through all the parameters for a single mechanism specified in the mechanism dictionary and
        calls self._parse_mech_content to interpret the rules and set the values for the given node.
        :param node: :class:'SHocNode'
        :param mech_name: str
        :param mech_content: dict
        """
        if not mech_content is None:
            # only modify synaptic mechanism parameters if synapses have been inserted
            if mech_name == 'synapse':
                if self.node_has_synapses(node):
                    for syn_type in mech_content:
                        for param_name in mech_content[syn_type]:
                            # accommodate multiple dict entries with different location constraints for a single
                            # parameter
                            if type(mech_content[syn_type][param_name]) == dict:
                                self._parse_mech_content(node, mech_name, param_name,
                                                         mech_content[syn_type][param_name], syn_type)
                            else:
                                for mech_content_entry in mech_content[syn_type][param_name]:
                                    self._parse_mech_content(node, mech_name, param_name, mech_content_entry, syn_type)
            else:
                for param_name in mech_content:
                    # accommodate multiple dict entries with different location constraints for a single parameter
                    if type(mech_content[param_name]) == dict:
                        self._parse_mech_content(node, mech_name, param_name, mech_content[param_name])
                    else:
                        for mech_content_entry in mech_content[param_name]:
                            self._parse_mech_content(node, mech_name, param_name, mech_content_entry)
        else:
            node.sec.insert(mech_name)

    def _parse_mech_content(self, node, mech_name, param_name, rules, syn_type=None):
        """
        This method loops through all the segments in a node and sets the value(s) for a single mechanism parameter by
        interpreting the rules specified in the mechanism dictionary. Properly handles ion channel gradients and
        inheritance of values from the closest segment of a specified type of section along the path from root. Also
        handles rules with distance boundaries, and rules to set synaptic (point process) parameters. For gradients,
        specifying a slope implies a linear gradient, while specifying both a slope and a tau implies an exponential
        gradient.
        :param node: :class:'SHocNode'
        :param mech_name: str
        :param param_name: str
        :param rules: dict
        :param syn_type: str
        """
        if mech_name == 'synapse':
            if syn_type is None:
                raise Exception('Cannot set synaptic mechanism parameter: {} without a specified point process'.format(param_name))
            elif not self.node_has_synapses(node, syn_type):
                return None  # ignore mechanism dictionary entries for types of synapses that have not been inserted
        if 'origin' in rules:  # an 'origin' with no 'value' inherits a starting parameter from the origin sec_type
                               # a 'value' with no 'origin' is independent of other sec_types
                               # an 'origin' with a 'value' uses the origin sec_type only as a reference point for
                               # applying a distance-dependent gradient
            if rules['origin'] == 'parent':
                if node.type == 'spine_head':
                    donor = node.parent.parent
                elif node.type == 'spine_neck':
                    donor = node.parent
                else:
                    donor = self.get_dendrite_origin(node)
            elif rules['origin'] in sec_types:
                donor = self._get_node_along_path_to_root(node, rules['origin'])
            else:
                raise Exception('Mechanism: {} parameter: {} cannot reference unknown sec_type: {}'.format(mech_name,param_name, rules['origin']))
        else:
            donor = None
        if 'value' in rules:
            baseline = rules['value']
        elif donor is None:
            raise Exception('Cannot set mechanism: {} parameter: {} without a specified origin or value'.format(
                mech_name, param_name))
        else:
            if (mech_name == 'cable') and (param_name == 'spatial_res'):
                baseline = self._get_spatial_res(donor)
            elif mech_name == 'synapse':
                if self.sec_type_has_synapses(donor.type, syn_type):
                    baseline = self._inherit_mech_param(node, donor, mech_name, param_name, syn_type)
                else:
                    raise Exception('Cannot inherit synaptic mechanism: {} parameter: {} from sec_type: {}'.format(
                                                                            syn_type, param_name, donor.type))
            else:
                baseline = self._inherit_mech_param(node, donor, mech_name, param_name)
        if mech_name == 'cable':  # cable properties can be inherited, but cannot be specified as gradients
            if param_name == 'spatial_res':
                node.init_nseg(baseline)
            else:
                setattr(node.sec, param_name, baseline)
                node.init_nseg(self._get_spatial_res(node))
            node.reinit_diam()
        else:
            if 'min_loc' in rules or 'max_loc' in rules or 'slope' in rules:
                if donor is None:
                    raise Exception('Cannot follow specifications for mechanism: {} parameter: {} without a provided '
                                    'origin'.format(mech_name, param_name))
                if mech_name == 'synapse':
                    self._specify_synaptic_parameter(node, param_name, baseline, rules, syn_type, donor)
                else:
                    if 'min_loc' in rules:
                        min_distance = rules['min_loc']
                    else:
                        min_distance = None
                    if 'max_loc' in rules:
                        max_distance = rules['max_loc']
                    else:
                        max_distance = None
                    min_seg_distance = self.get_distance_to_node(donor, node, 0.5/node.sec.nseg)
                    max_seg_distance = self.get_distance_to_node(donor, node, (0.5 + node.sec.nseg - 1)/node.sec.nseg)
                    # if any part of the section is within the location constraints, insert the mechanism, and specify
                    # the parameter at the segment level
                    if (min_distance is None or max_seg_distance >= min_distance) and \
                            (max_distance is None or min_seg_distance <= max_distance):
                        if not mech_name == 'ions':
                            node.sec.insert(mech_name)
                        if min_distance is None:
                            min_distance = 0.
                        for seg in node.sec:
                            seg_loc = self.get_distance_to_node(donor, node, seg.x)
                            if seg_loc >= min_distance and (max_distance is None or seg_loc <= max_distance):
                                if 'slope' in rules:
                                    seg_loc -= min_distance
                                    if 'tau' in rules:
                                        if 'xhalf' in rules:  # sigmoidal gradient
                                            offset = baseline - rules['slope'] / (1. +
                                                                    np.exp(rules['xhalf'] / rules['tau']))
                                            value = offset + rules['slope'] / (1. +
                                                                    np.exp((rules['xhalf'] - seg_loc) / rules['tau']))
                                        else:  # exponential gradient
                                            offset = baseline - rules['slope']
                                            value = offset + rules['slope'] * np.exp(seg_loc / rules['tau'])
                                    else:  # linear gradient
                                        value = baseline + rules['slope'] * seg_loc
                                    if 'min' in rules and value < rules['min']:
                                        value = rules['min']
                                    #elif value < 0.:
                                    #    value = 0.
                                    elif 'max' in rules and value > rules['max']:
                                        value = rules['max']
                                else:
                                    value = baseline
                            elif 'outside' in rules:        # by default, if only some segments in a section meet the
                                value = rules['outside']    # location constraints, the parameter inherits the
                            else:                           # mechanism's default value. if another value is desired, it
                                value = None                # can be specified via an 'outside' key in the mechanism
                            if not value is None:           # dictionary entry
                                if mech_name == 'ions':
                                    setattr(seg, param_name, value)
                                else:
                                    setattr(getattr(seg, mech_name), param_name, value)
            elif mech_name == 'ions':
                setattr(node.sec, param_name, baseline)
            elif mech_name == 'synapse':
                self._specify_synaptic_parameter(node, param_name, baseline, rules, syn_type)
            else:
                # print node.type, mech_name
                node.sec.insert(mech_name)
                setattr(node.sec, param_name+"_"+mech_name, baseline)

    def _specify_synaptic_parameter(self, node, param_name, baseline, rules, syn_type, donor=None):
        """
        This method interprets an entry from the mechanism dictionary to set parameters associated with a synaptic
        point_process mechanism that has been inserted either into a spine attached to this node, or inserted directly
        in this node. Appropriately implements slopes and inheritances.
        :param node: :class:'SHocNode'
        :param param_name: str
        :param baseline: float
        :param rules: dict
        :param syn_type: str
        :param donor: :class:'SHocNode' or None
        """
        syn_list = []
        syn_list.extend(node.synapses)
        for spine in node.spines:
            syn_list.extend(spine.synapses)
        if 'min_loc' in rules:
            min_distance = rules['min_loc']
        else:
            min_distance = 0.
        if 'max_loc' in rules:
            max_distance = rules['max_loc']
        else:
            max_distance = None
        if 'variance' in rules and rules['variance'] == 'normal':
            normal = True
        else:
            normal = False
        for syn in syn_list:
            if syn_type in syn._syn:  # not all synapses contain every synaptic mechanism
                target = syn.target(syn_type)
                if donor is None:
                    if normal:
                        value = self.random.normal(baseline, baseline / 6.)
                        setattr(target, param_name, value)
                    else:
                        setattr(target, param_name, baseline)
                else:
                    distance = self.get_distance_to_node(donor, node, syn.loc)
                    # note: if only some synapses in a section meet the location constraints, the synaptic parameter
                    # will maintain its default value in all other locations. values for other locations must be
                    # specified with an additional entry in the mechanism dictionary
                    if distance >= min_distance and (max_distance is None or distance <= max_distance):
                        if 'slope' in rules:
                            distance -= min_distance
                            if 'tau' in rules:  # exponential gradient
                                if 'xhalf' in rules:  # sigmoidal gradient
                                    offset = baseline - rules['slope'] / (1. + np.exp(rules['xhalf'] / rules['tau']))
                                    value = offset + rules['slope'] / (1. + np.exp((rules['xhalf'] - distance) /
                                                                                   rules['tau']))
                                else:
                                    offset = baseline - rules['slope']
                                    value = offset + rules['slope'] * np.exp(distance / rules['tau'])
                            else:  # linear gradient
                                value = baseline + rules['slope'] * distance
                            if 'min' in rules and value < rules['min']:
                                value = rules['min']
                            #elif value < 0.:
                            #    value = 0.
                            elif 'max' in rules and value > rules['max']:
                                value = rules['max']
                        else:
                            value = baseline
                        if normal:
                            value = self.random.normal(value, value / 6.)
                        setattr(target, param_name, value)

    def set_special_mech_param_linear_gradient(self, mech_name, param_name, sec_type_list, criterion, end_val):
        """
        This is an admittedly ad-hoc procedure to implement a linearly decreasing gradient of sodium channels in
        terminal branches that is not easily accomplished by the general procedures implementing the mechanism
        dictionary.
        :param mech_name: str
        :param param_name: str
        :param sec_type_list: list of str
        :param criterion: boolean function
        :param end_val: float
        """
        for sec_type in sec_type_list:
            for node in self.get_nodes_of_subtype(sec_type):
                if criterion(node):
                    start_val = getattr(node.sec(0.), param_name+'_'+mech_name)
                    slope = end_val - start_val
                    for seg in node.sec:
                        value = start_val + slope * seg.x
                        setattr(getattr(seg, mech_name), param_name, value)

    def init_spike_detector(self, node=None, loc=1., param='_ref_v', delay=None, weight=None, threshold=None,
                            target=None):
        """
        Converts analog voltage in the specified section to digital spike output. By default, initializes an h.NetCon
        object with voltage as a reference parameter and no target. Can later be connected by h.ParallelContext, or
        re-initialized with a target on a cell contained within the local processing environment.
        :param node: :class:'SHocNode'
        :param loc: float
        :param param: str
        :param delay: float
        :param weight: float
        :param threshold: float
        :return: :class:'h.NetCon'
        """
        if node is None:
            if self.axon:
                node = self.axon[-1]
            else:
                raise Exception('No source node specified for spike detector.')
        if self.spike_detector is not None:
            if delay is None:
                delay = self.spike_detector.delay
            if weight is None:
                weight = self.spike_detector.weight[0]
            if threshold is None:
                threshold = self.spike_detector.threshold
        else:
            if delay is None:
                delay = 0.
            if weight is None:
                weight = 1.
            if threshold is None:
                threshold = -30.
        self.spike_detector = h.NetCon(getattr(node.sec(loc), param), target, sec=node.sec)
        self.spike_detector.delay = delay
        self.spike_detector.weight[0] = weight
        self.spike_detector.threshold = threshold

    def get_dendrite_origin(self, node):
        """
        This method determines the section type of the given node, and returns the node representing the primary branch
        point for the given section type. Basal and trunk sections originate at the soma, and apical and tuft dendrites
        originate at the trunk. For spines, recursively calls with parent node to identify the parent branch first.
        :param node: :class:'SHocNode'
        :return: :class:'SHocNode'
        """
        sec_type = node.type
        if sec_type in ['spine_head', 'spine_neck']:
            return self.get_dendrite_origin(node.parent)
        elif sec_type in ['basal', 'trunk', 'axon_hill', 'ais', 'axon']:
            return self._get_node_along_path_to_root(node, 'soma')
        elif sec_type in ['apical', 'tuft']:
            return self._get_node_along_path_to_root(node, 'trunk')
        elif sec_type == 'soma':
            return node

    def _get_node_along_path_to_root(self, node, sec_type):
        """
        This method follows the path from the given node to the root node, and returns the first node with section type
        sec_type.
        :param node: :class:'SHocNode'
        :param sec_type: str
        :return: :class:'SHocNode'
        """
        parent = node
        while not parent is None:
            # print parent.name
            if parent in self.soma and not sec_type == 'soma':
                parent = None
            elif parent.type == sec_type:
                return parent
            else:
                parent = parent.parent
        raise Exception('The path from node: {} to root does not contain sections of type: {}'.format(node.name,
                                                                                                        sec_type))

    def _get_closest_synapse(self, node, loc, syn_type=None, downstream=True):
        """
        This method finds the closest synapse to the specified location within or downstream of the provided node. Used
        for inheritance of synaptic mechanism parameters. Can also look upstream instead. Can also find the closest
        synapse containing a synaptic point_process of a specific type.
        :param node: :class:'SHocNode'
        :param loc: float
        :param syn_type: str
        :return: :class:'Synapse'
        """

        syn_list = [syn for syn in node.synapses if syn_type is None or syn_type in syn._syn]
        for spine in node.spines:
            syn_list.extend([syn for syn in spine.synapses if syn_type is None or syn_type in syn._syn])
        if not syn_list:
            if downstream:
                for child in [child for child in node.children if child.type == node.type]:
                    target_syn = self._get_closest_synapse(child, 0., syn_type)
                    if target_syn is not None:
                        return target_syn
                return None
            elif node.parent.type == node.type:
                return self._get_closest_synapse(node.parent, 1., syn_type, downstream=False)
        else:
            min_distance = 1.
            target_syn = None
            for syn in syn_list:
                distance = abs(syn.loc - loc)
                if distance < min_distance:
                    min_distance = distance
                    target_syn = syn
            return target_syn

    def _inherit_mech_param(self, node, donor, mech_name, param_name, syn_type=None):
        """
        When the mechanism dictionary specifies that a node inherit a parameter value from a donor node, this method
        returns the value of that parameter found in the section or final segment of the donor node. For synaptic
        mechanism parameters, searches for the closest synapse in the donor node. If the donor node does not contain
        synapses due to location constraints, this method searches first child branches, then parent nodes of the same
        sec_type as the donor node.
        :param node: :class:'SHocNode'
        :param donor: :class:'SHocNode'
        :param mech_name: str
        :param param_name: str
        :param syn_type: str
        :return: float
        """
        try:
            if mech_name in ['cable', 'ions']:
                return getattr(donor.sec, param_name)
            elif mech_name == 'synapse':
                try:  # first look downstream for a nearby synapse, then upstream.
                    return getattr(self._get_closest_synapse(donor, 1., syn_type).target(syn_type), param_name)
                except (AttributeError, KeyError):
                    return getattr(self._get_closest_synapse(donor, 1., syn_type, downstream=False).target(syn_type),
                                   param_name)
            else:
                loc = donor.sec.nseg/(donor.sec.nseg + 1.)
                # accesses the last segment of the section
                return getattr(getattr(donor.sec(loc), mech_name), param_name)
        except (AttributeError, NameError, KeyError):
            if syn_type is None:
                print 'Exception: Mechanism: {} parameter: {} cannot be inherited from sec_type: {}'.format(mech_name,
                                                                                                param_name, donor.type)
            else:
                print 'Exception: Problem inheriting synaptic mechanism: {} parameter {} from sec_type: {}'.format(
                                                                                    syn_type, param_name, donor.type)
            raise KeyError

    def _get_spatial_res(self, node):
        """
        Checks the mechanism dictionary if the section type of this node has a specified spatial resolution factor.
        Used to scale the number of segments per section in the hoc model by a factor of an exponent of 3.
        :param node: :class:'SHocNode
        :return: int
        """
        try:  # if spatial_res has not been specified for the origin type of section, it defaults to 0
            rules = self.mech_dict[node.type]['cable']['spatial_res']
        except KeyError:
            return 0
        if 'value' in rules:
            return rules['value']
        elif 'origin' in rules:
            if rules['origin'] in sec_types:  # if this sec_type also inherits the value, continue following the path
                return self._get_spatial_res(self._get_node_along_path_to_root(node, rules['origin']))
            else:
                print 'Exception: Spatial resolution cannot be inherited from sec_type: {}'.format(rules['origin'])
                raise KeyError
        else:
            print 'Exception: Cannot set spatial resolution without a specified origin or value'
            raise KeyError

    def modify_mech_param(self, sec_type, mech_name, param_name=None, value=None, origin=None, slope=None, tau=None,xhalf=None, min=None, max=None, min_loc=None, max_loc=None, outside=None, syn_type=None, variance=None,replace=True):
        """
        Modifies or inserts new membrane mechanisms into hoc sections of type sec_type. First updates the mechanism
        dictionary, then sets the corresponding hoc parameters. This method is meant to be called manually during
        initial model specification, or during parameter optimization. For modifications to persist across simulations,
        the mechanism dictionary must be saved to a file using self.export_mech_dict() and re-imported during HocCell
        initialization.
        :param sec_type: str
        :param mech_name: str
        :param param_name: str
        :param value: float
        :param origin: str
        :param slope: float
        :param tau: float
        :param xhalf: float
        :param min: float
        :param max: float
        :param min_loc: float
        :param max_loc: float
        :param outside: float
        :param syn_type: str
        :param variance: str
        :param replace: bool
        """
        global verbose
        backup_content = None
        mech_content = None
        if not sec_type in sec_types:
            raise Exception('Cannot specify mechanism: {} parameter: {} for unknown sec_type: {}'.format(mech_name,
                                                                                                param_name, sec_type))
        if mech_name in ['cable', 'ions', 'synapse']:
            if param_name is None:
                raise Exception('No parameter specified for mechanism: {}'.format(mech_name))
        if not param_name is None:
            if value is None and origin is None:
                raise Exception('Cannot set mechanism: {} parameter: {} without a specified origin or value'.format(
                                                                                                mech_name, param_name))
            if mech_name == 'synapse':
                if syn_type is None:
                    raise Exception('Cannot set synaptic mechanism parameter: {} without a specified point '
                                    'process'.format(param_name))
                else:
                    self._modify_synaptic_mech_param(sec_type, param_name, value, origin, slope, tau, xhalf, min, max,
                                                     min_loc, max_loc, outside, syn_type, variance, replace)
                    return
            rules = {}
            if not origin is None:
                if not origin in sec_types+['parent']:
                    raise Exception('Cannot inherit mechanism: {} parameter: {} from unknown sec_type: {}'.format(
                                                                                    mech_name, param_name, origin))
                else:
                    rules['origin'] = origin
            if not value is None:
                rules['value'] = value
            if not slope is None:
                rules['slope'] = slope
            if not tau is None:
                rules['tau'] = tau
            if not xhalf is None:
                rules['xhalf'] = xhalf
            if not min is None:
                rules['min'] = min
            if not max is None:
                rules['max'] = max
            if not min_loc is None:
                rules['min_loc'] = min_loc
            if not max_loc is None:
                rules['max_loc'] = max_loc
            if not outside is None:
                rules['outside'] = outside
            # currently only implemented for synaptic parameters
            if not variance is None:
                rules['variance'] = variance
            mech_content = {param_name: rules}
        if not sec_type in self.mech_dict:  # No mechanisms have been inserted into this type of section yet
            self.mech_dict[sec_type] = {mech_name: mech_content}
        elif not mech_name in self.mech_dict[sec_type]:                 # This mechanism has not yet been inserted into
            backup_content = copy.deepcopy(self.mech_dict[sec_type])    # this type of section
            self.mech_dict[sec_type][mech_name] = mech_content
        elif self.mech_dict[sec_type][mech_name] is None:               # This mechanism has been inserted, but no
            backup_content = copy.deepcopy(self.mech_dict[sec_type])    # parameters have been specified
            self.mech_dict[sec_type][mech_name] = mech_content
        elif param_name in self.mech_dict[sec_type][mech_name]:         # This parameter has already been specified. Now
            backup_content = copy.deepcopy(self.mech_dict[sec_type])    # have to determine whether to replace or extend
            if replace:                                                 # the current dictionary entry.
                self.mech_dict[sec_type][mech_name][param_name] = rules
            elif type(self.mech_dict[sec_type][mech_name][param_name]) == dict:
                self.mech_dict[sec_type][mech_name][param_name] = [self.mech_dict[sec_type][mech_name][param_name],
                                                                   rules]
            elif type(self.mech_dict[sec_type][mech_name][param_name]) == list:
                self.mech_dict[sec_type][mech_name][param_name].append(rules)
        elif not param_name is None:                                    # This mechanism has been inserted, but this
            backup_content = copy.deepcopy(self.mech_dict[sec_type])    # parameter has not yet been specified
            self.mech_dict[sec_type][mech_name][param_name] = rules
        try:
            nodes = self.get_nodes_of_subtype(sec_type)
            if mech_name == 'cable':  # all membrane mechanisms in sections of type sec_type must be reinitialized after
                                      # changing cable properties
                if param_name in ['Ra', 'cm', 'spatial_res']:
                    self._reinit_mech(nodes, reset_cable=1)
                else:
                    print 'Exception: Unknown cable property: {}'.format(param_name)
                    raise KeyError
            else:
                for node in nodes:
                    try:
                        self._modify_mechanism(node, mech_name, mech_content)
                    except (AttributeError, NameError, ValueError, KeyError):
                        if not param_name is None:
                            print 'Exception: Problem modifying mechanism: {} parameter: {}'.format(mech_name,
                                                                                                    param_name)
                        else:
                            print 'Exception: Problem inserting mechanism: {}'.format(mech_name)
                        raise KeyError
        except KeyError:
            if backup_content is None:
                del self.mech_dict[sec_type]
            else:
                self.mech_dict[sec_type] = copy.deepcopy(backup_content)
        finally:
            if verbose:
                pprint.pprint(self.mech_dict)

    def _modify_synaptic_mech_param(self, sec_type, param_name=None, value=None, origin=None, slope=None, tau=None, xhalf=None, min=None, max=None, min_loc=None, max_loc=None, outside=None, syn_type=None,variance=None, replace=True):

        """
        Modifies or inserts new synaptic point process mechanisms into hoc sections of type sec_type. First updates the
        mechanism dictionary, then sets the corresponding hoc parameters. Internal method handles special nested
        dictionary specification for synaptic parameters.
        :param sec_type: str
        :param mech_name: str
        :param param_name: str
        :param value: float
        :param origin: str
        :param slope: float
        :param tau: float
        :param xhalf: float
        :param min: float
        :param max: float
        :param min_loc: float
        :param max_loc: float
        :param outside: float
        :param syn_type: str
        :param variance: str
        :param replace: bool
        """
        global verbose
        backup_content = None
        rules = {}
        if not origin is None:
            if not origin in sec_types + ['parent']:
                raise Exception('Cannot inherit synaptic mechanism: {} parameter: {} from unknown sec_type: {}'.format(
                    syn_type, param_name, origin))
            else:
                rules['origin'] = origin
        if not value is None:
            rules['value'] = value
        if not slope is None:
            rules['slope'] = slope
        if not tau is None:
            rules['tau'] = tau
        if not xhalf is None:
            rules['xhalf'] = xhalf
        if not min is None:
            rules['min'] = min
        if not max is None:
            rules['max'] = max
        if not min_loc is None:
            rules['min_loc'] = min_loc
        if not max_loc is None:
            rules['max_loc'] = max_loc
        if not outside is None:
            rules['outside'] = outside
        if not variance is None:
            rules['variance'] = variance
        mech_content = {param_name: rules}
        if not sec_type in self.mech_dict:  # No mechanisms have been inserted into this type of section yet
            self.mech_dict[sec_type] = {'synapse': {syn_type: mech_content}}
        elif not 'synapse' in self.mech_dict[sec_type]:  # No synaptic mechanisms have been specified in this type of
            backup_content = copy.deepcopy(self.mech_dict[sec_type])  # section yet
            self.mech_dict[sec_type]['synapse'] = {syn_type: mech_content}
        elif not syn_type in self.mech_dict[sec_type]['synapse']:  # This synaptic mechanism has not yet been inserted
            backup_content = copy.deepcopy(self.mech_dict[sec_type])  # into this type of section
            self.mech_dict[sec_type]['synapse'][syn_type] = mech_content
        elif self.mech_dict[sec_type]['synapse'][syn_type] is None:  # This mechanism has been inserted, but no
            backup_content = copy.deepcopy(self.mech_dict[sec_type])  # parameters have been specified
            self.mech_dict[sec_type]['synapse'][syn_type] = mech_content
        elif param_name in self.mech_dict[sec_type]['synapse'][syn_type]:  # This parameter has already been specified.
            backup_content = copy.deepcopy(self.mech_dict[sec_type])  # Now have to determine whether to replace or
            if replace:                                               # extend the current dictionary entry.
                self.mech_dict[sec_type]['synapse'][syn_type][param_name] = rules
            elif type(self.mech_dict[sec_type]['synapse'][syn_type][param_name]) == dict:
                self.mech_dict[sec_type]['synapse'][syn_type][param_name] = \
                    [self.mech_dict[sec_type]['synapse'][syn_type][param_name], rules]
            elif type(self.mech_dict[sec_type]['synapse'][syn_type][param_name]) == list:
                self.mech_dict[sec_type]['synapse'][syn_type][param_name].append(rules)
        elif not param_name is None:  # This mechanism has been inserted, but this
            backup_content = copy.deepcopy(self.mech_dict[sec_type])  # parameter has not yet been specified
            self.mech_dict[sec_type]['synapse'][syn_type][param_name] = rules
        try:
            nodes = self.get_nodes_of_subtype(sec_type)
            for node in nodes:
                try:
                    self._modify_mechanism(node, 'synapse', {syn_type: mech_content})
                except (AttributeError, NameError, ValueError, KeyError):
                    if not param_name is None:
                        print 'Exception: Problem modifying synaptic mechanism: {} parameter: {}'.format(syn_type,
                                                                                                param_name)
                    else:
                        print 'Exception: Problem inserting synaptic mechanism: {}'.format(syn_type)
                    raise KeyError
        except KeyError:
            if backup_content is None:
                del self.mech_dict[sec_type]
            else:
                self.mech_dict[sec_type] = copy.deepcopy(backup_content)
        finally:
            if verbose:
                pprint.pprint(self.mech_dict)

    def export_mech_dict(self, mech_filename=None):
        """
        Following modifications to the mechanism dictionary either during model specification or parameter optimization,
        this method stores the current mech_dict to a pickle file stamped with the date and time. This allows the
        current set of mechanism parameters to be recalled later.
        """
        if mech_filename is None:
            mech_filename = 'mech_dict_'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'.pkl'
        write_to_pkl(data_dir+mech_filename+'.pkl', self.mech_dict)
        print "Exported mechanism dictionary to "+mech_filename+'.pkl'

    def get_node_by_distance_to_soma(self, distance, sec_type):
        """
        Gets the first node of the given section type at least the given distance from a soma node.
        Not particularly useful, since it will always return the same node.
        :param distance: int or float
        :param sec_type: str
        :return: :class:'SHocNode'
        """
        nodes = self._node_dict[sec_type]
        for node in nodes:
            if self.get_distance_to_node(self.tree.root, node) >= distance:
                return node
        raise Exception('No node is {} um from a soma node.'.format(distance))

    def get_distance_to_node(self, root, node, loc=None):
        """
        Returns the distance from the given location on the given node to its connection with a root node.
        :param root: :class:'SHocNode'
        :param node: :class:'SHocNode'
        :param loc: float
        :return: int or float
        """
        length = 0.
        if node in self.soma:
            return length
        if not loc is None:
            length += loc*node.sec.L
        if root in self.soma:
            while not node.parent in self.soma:
                node.sec.push()
                loc = h.parent_connection()
                h.pop_section()
                node = node.parent
                length += loc*node.sec.L
        elif self.node_in_subtree(root, node):
            while not node.parent is root:
                node.sec.push()
                loc = h.parent_connection()
                h.pop_section()
                node = node.parent
                length += loc*node.sec.L
        else:
            return None  # node is not connected to root
        return length

    def get_distance_node_to_node(self, node1, node2, loc=None):
        """
        Returns the distance from the given location on the given node to its connection with a root node.
        :param root: :class:'SHocNode'
        :param node: :class:'SHocNode'
        :param loc: float
        :return: int or float
        """
        # get path between nodes
        # if root is in path, the list stops. create two lists: node1 to root, root to node2
        print 'getting path length between', node1.name, ' and ', node2.name
        path1 = self.tree.path_between_nodes(node1, node2)
        path2=[]
        if self.tree.root not in [node1, node2] and self.tree.root in path1:
            print 'oops, bumped into the root on the way to node2'
            path2 = self.tree.path_between_nodes(node2, self.tree.root)

        L=0.
        for node in path1:
            L+=node.sec.L
        for node in path2:
            L+=node.sec.L
        return L

    # def set_distance_from_soma(self, ):
    #     '''
    #     '''
    #     for tree_key, tree in self._node_dict.iteritems():
    #         for node in tree:
    #             for seg in node.sec:

    def clear_synapses(self, ):
        '''
        '''
        for node_name, node in self._node_names.iteritems():
            node.content['synapses']=[] 

    def set_node_names(self):
        '''
        '''
        self._node_names={}
        for tree_key, tree in self._node_dict.iteritems():
            for node in tree:
                self._node_names[node.name]=node

    def node_in_subtree(self, root, node):
        """
        Checks if a node is contained within a subtree of root.
        :param root: 'class':SNode2 or SHocNode
        :param node: 'class':SNode2 or SHocNode
        :return: boolean
        """
        nodelist = []
        self.tree._gather_nodes(root, nodelist)
        if node in nodelist:
            return True
        else:
            return False

    def get_path_length_swc(self, path):
        """
        Calculates the distance between nodes given a list of SNode2 nodes connected in a path.
        :param path: list : :class:'SNode2'
        :return: int or float
        """
        distance = 0
        for i in range(len(path)-1):
            distance += np.sqrt(np.sum((path[i].get_content()['p3d'].xyz - path[i+1].get_content()['p3d'].xyz)**2))
        return distance

    def get_node_length_swc(self, node):
        """
        Calculates the distance between the center points of an SNode2 node and its parent.
        :param node: :class:'SNode2'
        :return: float
        """
        if not self.tree.is_root(node):
            return np.sqrt(np.sum((node.get_content()['p3d'].xyz - node.parent.get_content()['p3d'].xyz)**2))
        else:
            return np.sqrt(np.sum(node.get_content()['p3d'].xyz**2))

    def get_branch_order(self, node):
        """
        Calculates the branch order of a SHocNode node. The order is defined as 0 for all soma, axon, and apical trunk
        dendrite nodes, but defined as 1 for basal dendrites that branch from the soma, and apical and tuft dendrites
        that branch from the trunk. Increases by 1 after each additional branch point. Makes sure not to count spines.
        :param node: :class:'SHocNode'
        :return: int
        """
        if node.type in ['soma', 'axon_hill', 'ais', 'axon']:
            return 0
        elif node.type == 'trunk':
            children = [child for child in node.parent.children if not child.type == 'spine_neck']
            if len(children) > 1 and children[0].type == 'trunk' and children[1].type == 'trunk':
                return 1
            else:
                return 0
        else:
            order = 0
            path = [branch for branch in self.tree.path_between_nodes(node, self.get_dendrite_origin(node)) if
                    not branch.type in ['soma', 'trunk']]
            for node in path:
                if self.is_terminal(node):
                    order += 1
                elif len([child for child in node.parent.children if not child.type == 'spine_neck']) > 1:
                    order += 1
                elif node.parent.type == 'trunk':
                    order += 1
            return order

    def is_terminal(self, node):
        """
        Calculates if a node is a terminal dendritic branch.
        :param node: :class:'SHocNode'
        :return: bool
        """
        if node.type in ['soma', 'axon_hill', 'ais', 'axon']:
            return False
        else:
            return not bool([child for child in node.children if not child.type == 'spine_neck'])

    def is_bifurcation(self, node, child_type):
        """
        Calculates if a node bifurcates into at least two children of specified type.
        :param node: :class:'SHocNode'
        :param child_type: string
        :return: bool
        """
        return len([child for child in node.children if child.type == child_type]) >= 2

    def get_seg_3d_locations(self, node, loc=None):
        """ given a neuron section, output the 3d coordinates of each segment in the section

        ouput is a nested list as [xyz dimension][segment number], with x,y, z dimensions listed in that order

        """
        # if node.get_parent() is None:
        #     sec_length = np.sqrt(np.sum(node.p3d.xyz**2))
        # else:
        #     sec_length = np.sqrt(np.sum((node.get_content()['p3d'].xyz - node.parent.get_content()['p3d'].xyz)**2))

        node_x = node.p3d.xyz[0]
        node_y = node.p3d.xyz[1]
        node_z = node.p3d.xyz[2]

        if node.type=='soma':
            parent_x =node.p3d.xyz[0]
            parent_z =node.p3d.xyz[2]
            if '0' in node.name:
                parent_y =node.p3d.xyz[1] + node.sec.L
            elif '1' in node.name:
                parent_y =node.p3d.xyz[1] - node.sec.L
        else:
            parent_x = node.parent.p3d.xyz[0]
            parent_y = node.parent.p3d.xyz[1]
            parent_z = node.parent.p3d.xyz[2]

        nseg = node.sec.nseg
        x=[]
        y=[]
        z=[]

        if loc is not None:
            if type(loc)!=list:
                loc = [loc]
            for seg_pos in loc:

                seg_x = parent_x + seg_pos*(node_x-parent_x)
                seg_y = parent_y + seg_pos*(node_y-parent_y)
                seg_z = parent_z + seg_pos*(node_z-parent_z)

                x.append(seg_x)
                y.append(seg_y)
                z.append(seg_z)

            return (x, y, z)
        else: 
            for seg_i, seg in enumerate(node.sec):
                seg_pos = seg.x 

                seg_x = parent_x + seg_pos*(node_x-parent_x)
                seg_y = parent_y + seg_pos*(node_y-parent_y)
                seg_z = parent_z + seg_pos*(node_z-parent_z)

                x.append(seg_x)
                y.append(seg_y)
                z.append(seg_z)

            return (x, y, z)

    def set_spine_p3d(self, ):
        ''' sets p3d attribute for spines to the same xyz location as the parent dendritic segment for both the spine head and spine neck
        '''
        # check for spines 
        if 'spine' in self._node_dict:
            # iterate over spines
            for spine in self._node_dict['spine']:
                # cheeck if spine head or neck, get head, neck, dend nodes
                if spine.type == 'spine_neck':
                    dend=spine.parent
                    neck=spine
                    head=spine.children[0]
                elif spine.type=='spine_head':
                    dend=spine.parent.parent
                    neck=spine.parent
                    head=spine
                # get location (0-1) along parent dendritic segment where spine neck is attached
                dend_loc = h.parent_connection(sec=neck.sec)
                # get 3d coordinates of parent dendritic segment
                dend_x, dend_y, dend_z = self.get_seg_3d_locations(node=dend, loc=dend_loc)
                dend_xyz = np.array([dend_x[0], dend_y[0], dend_z[0]])
                # create p3d in spine head and neck nodes
                neck.p3d = btmorph.P3D2(xyz=dend_xyz, radius=neck.sec.diam/2., type='spine_neck')
                head.p3d = btmorph.P3D2(xyz=dend_xyz, radius=head.sec.diam/2., type='spine_head')

    def get_node_3d_locations(self, node, seg_loc=None, set_content=True):
        """ given a neuron section, output the 3d coordinates of each segment in the section

        ouput is a nested list as [xyz dimension][segment number], with x,y, z dimensions listed in that order

        """
        sec=node.sec
        # number of 3d points in section
        tol =.001
        # specifies the location of each segment so that n3d can be used
        # h.define_shape()
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
        if set_content:
            node.content['seg_x']=seg_x
            node.content['seg_y']=seg_y
            node.content['seg_z']=seg_z
        return [seg_x, seg_y, seg_z]
    
    def get_apical_trunk_from_raw_tree(self, raw_tree, apical_type=4., tol=1.0):
        ''' find nodes that correspond to apical trunk by iterating through the apical tree. at each bifurcation, get the subtrees of each of the children. if one subtrees contains substantially more area (determined by tol), then the apical trunk continues along that subtree until all children has similar subtree sizes
        '''
        # apical_key = 'apical'
        # get apical nodes connected to soma
        root = raw_tree.root
        print root
        print root.children
        trunk = []
        for node in root.children:
            if node.get_content()['p3d'].type==apical_type:
                trunk.append(node)
                while len(trunk[-1].children)>0:
                    if len(trunk[-1].children)==1:
                        trunk.append(trunk[-1].children[0])
                    else:
                        total_subtree_areas=[]
                        for child_i, child in enumerate(trunk[-1].children):
                            subtree = raw_tree.get_sub_tree(child)
                            subtree_nodes = subtree.get_nodes()
                            areas=[]
                            for subtree_node in subtree_nodes:
                                if subtree_node!=subtree.root:
                                    # get length
                                    L = np.sqrt(np.sum((subtree_node.get_content()['p3d'].xyz - subtree_node.parent.get_content()['p3d'].xyz)**2))
                                    r = subtree_node.get_content()['p3d'].radius
                                    area = L*(np.pi*r**2)
                                    areas.append(area)
                            total_area = np.sum(areas)
                            total_subtree_areas.append(total_area)
                        trunk_is = np.argsort(total_subtree_areas)
                        if total_subtree_areas[trunk_is[-1]]/total_subtree_areas[trunk_is[-2]] > tol:
                            trunk_i = np.argmax(total_subtree_areas)
                            trunk.append(trunk[-1].children[trunk_i])
                        else:
                            break
        return trunk

    def update_raw_tree_node_types(self, raw_tree, node_list, new_type):
        '''
        '''
        for node in raw_tree.get_nodes():
            if node in node_list:
                node.content['p3d'].type=new_type

    def set_stochastic_synapses(self, value):
        """
        This method turns stochastic filtering of presynaptic release on or off for all synapses contained in this cell.
        :param value: int in [0, 1]
        """
        for nodelist in self._node_dict.itervalues():
            for node in nodelist:
                for syn in node.synapses:
                    syn.stochastic = value

    def choose_syn_locations_random(self, n, sec_types, on_spine=True, min_distance=None, max_distance=None, replace=True):
        '''
        '''
        if type(sec_types)!=list:
            sec_types = [sec_types]
        # get list of eligible synapses
        potential_nodes=[]
        for sec_type in sec_types:
            for tree_key, tree in self._node_dict.iteritems():
                for node in tree:
                    if node.type==sec_type:
                        # if looking for spines, set the node to its correspondig spine
                        if on_spine and len(node.spines)>0:
                            for spine in node.spines:
                                if spine.type=='spine_head':
                                    node=spine

                        distance_from_soma=self.get_distance_node_to_node(node1=node, node2=self.soma[0])
                        # if no require is given
                        if min_distance is None and max_distance is None:
                            potential_nodes.append(node)
                        # if min and max distances are given
                        elif min_distance is not None and max_distance is not None and distance_from_soma>=min_distance and distance_from_soma<=max_distance:
                                potential_nodes.append(node)
                        # if only min distance is given
                        elif min_distance is not None and distance_from_soma>=min_distance:
                                potential_nodes.append(node)
                        # if only max is given
                        elif max_distance is not None and distance_from_soma<=max_distance:
                            potential_nodes.append(node)


        # choose random number n from list of sections that match criteria
        if len(potential_nodes) >0:
            syns = np.random.choice(potential_nodes, n, replace=replace)
        else:
            raise Exception(n, ' synases requested, but only ', len(potential_nodes), 'synapses found that meet criteria')
        return syns

    def get_terminal_branches(self, sec_types, min_distance=None, max_distance=None):
        '''
        '''
        if type(sec_types)!=list:
            sec_types=[sec_types]
        # get list of eligible synapses
        potential_nodes=[]
        for sec_type in sec_types:
            for tree_key, tree in self._node_dict.iteritems():
                for node in tree:
                    if node.type==sec_type:
                        if self.is_terminal(node):


                            distance_from_soma=self.get_distance_node_to_node(node1=node, node2=self.soma[0])
                            print distance_from_soma
                            # if no require is given
                            if min_distance is None and max_distance is None:
                                potential_nodes.append(node)
                            # if min and max distances are given
                            elif min_distance is not None and max_distance is not None and distance_from_soma>=min_distance and distance_from_soma<=max_distance:
                                    potential_nodes.append(node)
                            # if only min distance is given
                            elif min_distance is not None and distance_from_soma>=min_distance:
                                    potential_nodes.append(node)
                            # if only max is given
                            elif max_distance is not None and distance_from_soma<=max_distance:
                                potential_nodes.append(node)

        return potential_nodes

    def set_nseg(self, node, seg_L=None, nseg=None):
        """ set number of segments for branch that was selected to activate synapses
        
        arguments:

        """


        # get section length
        sec_L = node.sec.L
        print 'section length', node.sec.L

        if seg_L is not None:
            # determine number of segments
            n_seg = int(np.ceil(sec_L/seg_L))
        # # check that number of segments is odd
        if n_seg % 2 == 0:

            n_seg+=1
        print n_seg

        # # set number of segments
        node.sec.nseg = n_seg
        print 'nseg', node.sec.nseg

    def create_geo_syns(self, ):
        ''' FIXME comments
        '''
        self.geo={}
        self.syns={}
        self.hoc_style={}
        for tree_key, tree in self._node_dict.iteritems():
            self.geo[tree_key]=[]
            self.syns[tree_key]=[]
            for node_i, node in enumerate(tree):
                self.geo[tree_key].append(node.sec)
                self.syns[tree_key].append([])
                for seg_i, seg in enumerate(node.sec):
                    loc = seg.x
                    self.syns[tree_key][node_i].append([])
                    for syn in node.synapses:
                        self.syns[tree_key][node_i][seg_i].append(syn)





                        # self.syns[tree_key][node_i][seg_i].append({})
                        # for syn_type in syn._syn.keys():
                        #     syn_object = getattr(h, syn_type)(self.node.sec(loc))
                        #     self.syns[tree_key][node_i][seg_i][syn_type]=syn_object

    def taper_diam(self, sec,zero_bound,one_bound):
        '''FIXME comments
        '''

        for (seg, d) in izip(sec, np.linspace(zero_bound, one_bound, sec.nseg)):
            seg.diam=d

    def update_diameters(self, ):
        '''FIXME comments
        '''
        print 'updating diameters'
        for sec_type in self._node_dict:
            for node in self._node_dict[sec_type]:
                if node != self.tree.root:
                    zero_bound = node.parent.sec(1).diam
                    one_bound = node.sec(1).diam
                    self.taper_diam(sec=node.sec, zero_bound=zero_bound, one_bound=one_bound)

    def get_all_synapses_and_netcons(self, ):
        '''
        '''
        self.all_synapses=[]
        self.all_netcons=[]
        for tree_key, tree in self._node_dict.iteritems():
            for node_i, node in enumerate(tree):
                for syn_object in node.synapses:
                    for syn_type in syn_object._syn.keys():
                        point_process = syn_object._syn[syn_type]['target']
                        netcon = syn_object._syn[syn_type]['netcon']
                        self.all_synapses.append(point_process)
                        self.all_netcons.append(netcon)




    @property
    def gid(self):
        return self._gid

    @property
    def soma(self):
        return self._node_dict['soma']

    @property
    def axon(self):
        return self._node_dict['axon']

    @property
    def basal(self):
        return self._node_dict['basal']

    @property
    def apical(self):
        return self._node_dict['apical']

    @property
    def trunk(self):
        return self._node_dict['trunk']

    @property
    def tuft(self):
        return self._node_dict['tuft']

    @property
    def spine(self):
        return self._node_dict['spine']

#------------------------------Extend SNode2 to interact with NEURON hoc sections------------------------

class HocStyleCell(object):
    ''' create a hoc style cell where section types are attributes and each section type contains a list of hoc sections. used for passing to neuron_reduce
    '''
    def __init__(self, **kwargs):
        '''
        '''
        pass

    def from_milstein(self, cell, **kwargs):
        '''
        '''
        assert hasattr(cell, '_node_dict'),"cell object does not have _node_dict"
        # iterate over sec_types
        for key, val in cell._node_dict.iteritems():
            section_list = [ _node.sec for _node in val ]

            # if sec_type does exist yet
            if not hasattr(self, key):
                if key=='soma':
                    setattr(self, key, section_list[0])
                else:

                    # create list for sec_type
                    setattr(self, key, section_list)
            # for node in val:
            #     getattr(self, key).append(node.sec)

    def from_geo(self, cell, **kwargs):
        '''
        '''
        assert hasattr(cell, 'geo'),'cell object does not have attribute "geo"'

        name_map = {
        'soma': ('soma', 'somatic'),
        'axon': ('axon', 'axonal'),
        'apical_tuft': ('apic', 'apical'),
        'apical_trunk': ('apic', 'apical'),
        'basal': ('dend', 'basal'),
        }
        # iterate over sec_types
        for key, val in cell.geo.iteritems():
            section_list = [ _sec for _sec in val ]
            hoc_section_list = h.SectionList(section_list)

            if key in name_map: 
                sec_name = name_map[key][0]
                seclist_name = name_map[key][1]
                # sections
                if not hasattr(self, sec_name):
                    setattr(self, sec_name, section_list)
                elif hasattr(self, sec_name):
                    updated = getattr(self, sec_name)+section_list
                    setattr(self, sec_name, updated)

                # SectionLists
                if not hasattr(self, seclist_name):
                    setattr(self, seclist_name, hoc_section_list)
                elif hasattr(self, seclist_name):
                    for _sec in section_list:
                        getattr(self, seclist_name).append(_sec)



class SHocNode(btmorph.btstructs2.SNode2):
    """
    Extends SNode2 with some methods for storing and retrieving additional information in the node's content
    dictionary related to running NEURON models specified in the hoc language.
    """

    def __init__(self, index=0):
        """
        :param index: int : unique node identifier
        """
        btmorph.btstructs2.SNode2.__init__(self, index)
        self.content['spines'] = []
        self.content['synapses'] = []

    def get_sec(self):
        """
        Returns the hoc section associated with this node, stored in the node's content dictionary.
        :return: :class:'neuron.h.Section'
        """
        if 'sec' in self.content:
            return self.content['sec']
        else:
            raise Exception('This node does not yet have an associated hoc section.')

    def set_sec(self, sec):
        """
        Stores the hoc section associated with this node in the node's content dictionary.
        :param sec: :class:'neuron.h.Section'
        """
        self.content['sec'] = sec

    sec = property(get_sec, set_sec)

    def get_p3d(self):
        """
        Returns the hoc section associated with this node, stored in the node's content dictionary.
        :return: :class:'neuron.h.Section'
        """
        if 'p3d' in self.content:
            return self.content['p3d']
        else:
            raise Exception('This node does not have p3d information.')

    def set_p3d(self, p3d):
        """
        Stores the hoc section associated with this node in the node's content dictionary.
        :param sec: :class:'neuron.h.Section'
        """
        self.content['p3d'] = p3d

    p3d = property(get_p3d, set_p3d)

    def init_nseg(self, spatial_res=0):
        """
        Initializes the number of hoc segments in this node's hoc section (nseg) based on the AC length constant.
        Must be re-initialized whenever basic cable properties Ra or cm are changed. If the node is a tapered cylinder,
        it should contain at least 3 segments. The spatial resolution parameter increases the number of segments per
        section by a factor of an exponent of 3.
        If a section's nseg has been manually increased beyond the suggestion of the mechanism dictionary, this method
        does not decrease it.
        :param spatial_res: int
        """
        sugg_nseg = d_lambda_nseg(self.sec)
        if not self.get_diam_bounds() is None:
            sugg_nseg = max(sugg_nseg, 3)
        sugg_nseg *= 3**spatial_res
        if self.sec.nseg < sugg_nseg:
            self.sec.nseg = sugg_nseg

    def reinit_diam(self):
        """
        For a node associated with a hoc section that is a tapered cylinder, every time the spatial resolution
        of the section (nseg) is changed, the section diameters must be reinitialized. This method checks the
        node's content dictionary for diameter boundaries and recalibrates the hoc section associated with this node.
        """
        if not self.get_diam_bounds() is None:
            [diam1, diam2] = self.get_diam_bounds()
            h('diam(0:1)={}:{}'.format(diam1, diam2), sec=self.sec)

    def get_diam_bounds(self):
        """
        If the hoc section associated with this node is a tapered cylinder, this method returns a list containing
        the values of the diameters at the 0 and 1 ends of the section, stored in the node's content dictionary.
        Otherwise, it returns None (for non-conical cylinders).
        :return: (list: int) or None
        """
        if 'diam' in self.content:
            return self.content['diam']
        else:
            return None

    def set_diam_bounds(self, diam1, diam2):
        """
        For a node associated with a hoc section that is a tapered cylinder, this stores a list containing the values
        of the diameters at the 0 and 1 ends of the section in the node's content dictionary.
        :param diam1: int
        :param diam2: int
        """
        self.content['diam'] = [diam1, diam2]
        self.reinit_diam()

    def get_type(self):
        """
        NEURON sections are assigned a node type for convenience in order to later specify membrane mechanisms and
        properties for each type of compartment.
        :return: str
        """
        if 'type' in self.content:
            return self.content['type']
        else:
            raise Exception('This node does not yet have a defined type.')

    def set_type(self, type):
        """
        Checks that type is a string in the list of defined section types, and stores the value in the node's content
        dictionary.
        :param type: str
        """
        if type in sec_types:
            self.content['type'] = type
        else:
            raise Exception('That is not a defined type of section.')

    type = property(get_type, set_type)

    def connect(self, parent, pindex=1, cindex=0):
        """
        Connects this SHocNode node to a parent node, and establishes a connection between their associated
        hoc sections.
        :param parent: :class:'SHocNode'
        :param pindex: int in [0,1] Connect to this end of the parent hoc section.
        :param cindex: int in [0,1] Connect this end of the child hoc section
        """
        self.parent=parent
        parent.add_child(self)
        self.sec.connect(parent.sec, pindex, cindex)

    @property
    def name(self):
        """
        Returns a str containing the name of the hoc section associated with this node. Consists of a type descriptor
        and an index identifier.
        :return: str
        """
        if 'type' in self.content:
            return '{0.type}{0.index}'.format(self)
        else:
            raise Exception('This node does not yet have a defined type.')

    @property
    def spines(self):
        """
        Returns a list of the spine head sections attached to the hoc section associated with this node.
        :return: list of :class:'SHocNode' of sec_type == 'spine_head'
        """
        return self.content['spines']

    @property
    def synapses(self):
        """
        Returns a list of the objects of :class:'Synapse' associated with this node.
        :return: list of hoc objects, type depends on .mod file(s) used to implement synapses
        """
        return self.content['synapses']

    # def set_synapses(self)


class CA1_Pyr(HocCell):
    def __init__(self, morph_filename=None, mech_filename=None, mech_dict=None, full_spines=False):
        HocCell.__init__(self, morph_filename=morph_filename, mech_filename=mech_filename, mech_dict=mech_dict)
        self.random.seed(self.gid)  # This cell will always have the same spine and GABA_A synapse locations as long as
                                    # they are inserted in the same order
        if full_spines:
            self.insert_spines_in_subset(['basal', 'trunk', 'apical', 'tuft'])

    def insert_spines_in_subset(self, sec_type_list):
        """
        This method populates the cell tree with spines following spine density information from Erk Bloss &
        Nelson Spruston. Basal dendrites have no spines until the first branch point, and a higher density beyond the
        second branch point. Trunk dendrites have no spines until the first branch point, and an increasing density
        until the tuft branch point(s). Apical dendrites have a density that varies with the distance from the soma of
        their original branch point from the trunk. Terminal tuft branches have a higher density than their parents.
        Should standardize the implementation of the rules for each type of dendrite and import the
        density dictionary from a file, similar to the implementation of the membrane mechanism dictionary.
        :param sec_type_list: list of str
        """
        densities = {'trunk': {'min': 0.2418, 'max': 3.8,
                               'start': min([self.get_distance_to_node(self.tree.root, branch) for branch in self.apical]), 'end': max([self.get_distance_to_node(self.tree.root, branch) for branch in self.trunk])},
                     'basal': {'1': 0., '2': 0.4428, '>2': 1.891},
                     'apical': {'min': 2.273, 'max': 2.688,
                                'start': min([self.get_distance_to_node(self.tree.root, branch) for branch in self.apical]),
                                'end': max([self.get_distance_to_node(self.tree.root, branch) for branch in self.apical if self.get_branch_order(branch) == 1])},
                     'tuft': {'parent': 1.354, 'terminal': 0.7157}
                    }
        if 'basal' in sec_type_list:
            for node in self.basal:
                order = self.get_branch_order(node)
                if order == 2:
                    self.insert_spines_every(node, densities['basal']['2'])
                elif order > 2:
                    self.insert_spines_every(node, densities['basal']['>2'])
        if 'trunk' in sec_type_list:
            for node in self.trunk:
                distance = self.get_distance_to_node(self.tree.root, node)
                if distance >= densities['trunk']['start']:
                    slope = (densities['trunk']['max'] - densities['trunk']['min']) / \
                            (densities['trunk']['end'] - densities['trunk']['start'])
                    density = densities['trunk']['min'] + slope * (distance - densities['trunk']['start'])
                    self.insert_spines_every(node, density)
        if 'apical' in sec_type_list:
            for node in self.apical:
                distance = self.get_distance_to_node(self.tree.root, self.get_dendrite_origin(node), loc=1.)
                slope = (densities['apical']['max'] - densities['apical']['min']) / \
                        (densities['apical']['end'] - densities['apical']['start'])
                density = densities['apical']['min'] + slope * (distance - densities['apical']['start'])
                self.insert_spines_every(node, density)
        if 'tuft' in sec_type_list:
            for node in self.tuft:
                if self.is_terminal(node):
                    self.insert_spines_every(node, densities['tuft']['terminal'])
                else:
                    self.insert_spines_every(node, densities['tuft']['parent'])
        self._reinit_mech(self.spine)

    def insert_spines_every(self, node, density):
        """
        Given a mean spine density in /um, insert spines in the node at the specified density.
        :param node: :class:'SHocNode'
        :param density: float: mean density in /um
        """
        L = node.sec.L
        beta = 1./density
        interval = self.random.exponential(beta)
        while interval < L:
            self.insert_spine(node, interval/L)
            interval += self.random.exponential(beta)

    def insert_spine(self, node, parent_loc, child_loc=0):
        """
        Spines consist of two hoc sections: a cylindrical spine head and a cylindrical spine neck.
        :param node: :class:'SHocNode'
        :param parent_loc: float
        :param child_loc: int
        """
        neck = self.make_section('spine_neck')
        neck.connect(node, parent_loc, child_loc)
        neck.sec.L = 1.58
        neck.sec.diam = 0.077
        self._init_cable(neck)
        head = self.make_section('spine_head')
        head.connect(neck)
        node.spines.append(head)
        head.sec.L = 0.5  # open cylinder, matches surface area of sphere with diam = 0.5
        head.sec.diam = 0.5
        self._init_cable(head)

    def insert_inhibitory_synapses_in_subset(self, sec_type_list=None):
        """

        :param sec_type_list: str
        """
        if sec_type_list is None:
            sec_type_list = ['soma', 'ais', 'basal', 'trunk', 'apical', 'tuft']
        densities = {'soma': 2.857,  # 4.285,
                     'ais': 0.53,
                     'trunk': {'min': 0.3022, 'max': 0.0627,
                               'start': 0.,
                               'end': max([self.get_distance_to_node(self.tree.root, branch) for branch in
                                                                                                self.trunk])},
                     'basal': {'primary': 0.3129, 'intermediate': 0.1728, 'terminal': 0.06543},
                     'apical': {'min': 0.03885, 'max': 0.04512,
                                'start': min([self.get_distance_to_node(self.tree.root, branch) for branch in
                                                                                                self.apical]),
                                'end': max([self.get_distance_to_node(self.tree.root, branch)
                                            for branch in self.apical if self.get_branch_order(branch) == 1])},
                     'tuft': {'parent': 0.2104, 'terminal': 0.1619}
                    }
        if 'soma' in sec_type_list:
            for node in self.soma:
                self.insert_inhibitory_synapse_every(node, densities['soma'])
        if 'ais' in sec_type_list:
            for node in self.get_nodes_of_subtype('ais'):
                self.insert_inhibitory_synapse_every(node, densities['ais'])
        if 'basal' in sec_type_list:
            for node in self.basal:
                if self.is_terminal(node):
                    self.insert_inhibitory_synapse_every(node, densities['basal']['terminal'])
                else:
                    order = self.get_branch_order(node)
                    if order == 1:
                        self.insert_inhibitory_synapse_every(node, densities['basal']['primary'])
                    else:
                        self.insert_inhibitory_synapse_every(node, densities['basal']['intermediate'])
        if 'trunk' in sec_type_list:
            for node in self.trunk:
                distance = self.get_distance_to_node(self.tree.root, node)
                if distance >= densities['trunk']['start']:
                    slope = (densities['trunk']['max'] - densities['trunk']['min']) / \
                            (densities['trunk']['end'] - densities['trunk']['start'])
                    density = densities['trunk']['min'] + slope * (distance - densities['trunk']['start'])
                    self.insert_inhibitory_synapse_every(node, density)
        if 'apical' in sec_type_list:
            for node in self.apical:
                distance = self.get_distance_to_node(self.tree.root, self.get_dendrite_origin(node), loc=1.)
                slope = (densities['apical']['max'] - densities['apical']['min']) / \
                        (densities['apical']['end'] - densities['apical']['start'])
                density = densities['apical']['min'] + slope * (distance - densities['apical']['start'])
                self.insert_inhibitory_synapse_every(node, density)
        if 'tuft' in sec_type_list:
            for node in self.tuft:
                if self.is_terminal(node):
                    self.insert_inhibitory_synapse_every(node, densities['tuft']['terminal'])
                else:
                    self.insert_inhibitory_synapse_every(node, densities['tuft']['parent'])

    def insert_inhibitory_synapse_every(self, node, density, syn_types=['GABA_A_KIN'], stochastic=0):
        """

        :param node: :class:'SHocNode'
        :param density: float: mean density in /um
        :param syn_types: list of str
        :param stochastic: int
        """
        L = node.sec.L
        beta = 1./density
        interval = self.random.exponential(beta)
        while interval < L:
            syn = Synapse(self, node, type_list=syn_types, stochastic=stochastic, loc=interval/L)
            interval += self.random.exponential(beta)

    def zero_na(self):
        """
        Set na channel conductances to zero in all compartments. Used during parameter optimization.
        """
        for sec_type in ['soma', 'axon_hill', 'ais', 'axon', 'basal', 'trunk', 'apical', 'tuft']:
            for na_type in (na_type for na_type in ['nas_kin', 'nat_kin', 'nas', 'nax'] if na_type in
                    self.mech_dict[sec_type]):
                self.modify_mech_param(sec_type, na_type, 'gbar', 0.)

    def zero_h(self):
        """
        Set ih conductances to zero in all compartments. Used during parameter optimization.
        """
        self.modify_mech_param('soma', 'h', 'ghbar', 0.)
        self.mech_dict['trunk']['h']['ghbar']['value'] = 0.
        self.mech_dict['trunk']['h']['ghbar']['slope'] = 0.
        for sec_type in ['basal', 'trunk', 'apical', 'tuft']:
            self.reinitialize_subset_mechanisms(sec_type, 'h')

    def set_terminal_branch_na_gradient(self, gmin=0.):
        """
        This is an admittedly ad-hoc procedure to implement a linearly decreasing gradient of sodium channels in
        terminal branches that is not easily accomplished by the general procedures implementing the mechanism
        dictionary.
        """
        na_type = (na_type for na_type in ['nas_kin', 'nat_kin', 'nas', 'nax']
                   if na_type in self.mech_dict['trunk']).next()
        self.set_special_mech_param_linear_gradient(na_type, 'gbar', ['basal', 'apical', 'tuft'],
                                                    self.is_terminal, gmin)


class Synapse(object):
    """
    The implementation in hoc of synaptic mechanisms that can be triggered is complicated. This container is an attempt
    to wrap all the objects required to deliver synaptic events to a section, and have separable synaptic mechanisms
    (e.g. GluA-Rs and GluN-Rs) respond with individually specifiable weights and kinetics.
    To make model specification and simulation implementation straightforward, synapses are not meant to be moved once
    they are initialized.
    """
    def __init__(self, cell, node, type_list=None, stochastic=1, loc=0.5, delay=0., weight=1., threshold=-30.,
                 stochastic_type=None, source=None, source_node=None, source_param=None, source_loc=0.5):
        """
        Design goals: A source (like a spike detector in a presynaptic neuron) can be specified. If not, a VecStim
        object is used a source, which can be played events at specified times using its .play method. If stochastic,
        all spikes are intercepted by a point process with release probability dynamics and its own unique and
        independent random variable from a uniform distribution. If not, the specified synaptic mechanisms are connected
        directly to the source of spikes.
        :param cell: :class:'HocCell'
        :param node: :class:'SHoCNode'
        :param type_list: list of str
        :param stochastic: int in [0, 1]
        :param loc: float
        :param delay: float
        :param weight: float
        :param threshold: float
        :param stochastic_type: str
        :param source: hoc artificial cell or otherwise hoc object not associated with a section
        :param source_node: :class:'SHocNode' specifies the section containing a range variable to be used as a source
        :param source_param: str
        :param source_loc: str
        """
        self._cell = cell
        self._node = node
        self._stochastic = stochastic
        self._loc = loc
        self._delay = delay
        self._weight = weight
        self._threshold = threshold
        self._syn = {}
        self.randObj = None
        self._node.synapses.append(self)
        if source_param is None:
            source_param = '_ref_v'
        if stochastic_type is None:
            self._stochastic_type = 'Pr'
        else:
            self._stochastic_type = stochastic_type
        if not source is None:
            self._source = {'object': source}
        elif not source_node is None:
            self._source = {'object': getattr(source_node.sec(source_loc), source_param), 'node': source_node}
        else:
            self._source = {'object': h.VecStim()}
        if type_list is None:
            type_list = ['AMPA_KIN']
        elif type(type_list) is not list:
            type_list = [type_list]
        if self.stochastic:
            self._init_stochastic()
        for target in type_list:
            syn = getattr(h, target)(self.node.sec(self._loc))
            self._syn[target] = {'target': syn}
            if self.stochastic:
                self._syn[target]['netcon'] = h.NetCon(self.target(self._stochastic_type), syn)
                self.netcon(target).delay = delay
                self.netcon(target).weight[0] = weight
                self.netcon(target).threshold = threshold
            else:
                self._init_netcon(target)

    def _init_stochastic(self):
        """
        This method constructs and initializes a stochastic filtering mechanism that intercepts spikes delivered to this
        synapse and calculates whether or not to pass a spike to the rest of the specified synaptic mechanisms.
        """
        if self.randObj is None:  # if this synapse has never been stochastic, it needs a new random number generator
            self.randObj = h.Random()
            self.randObj.MCellRan4(self.cell.gid*1e4+1, self.node.index*1e4+self.node.synapses.index(self)+1)
            # a unique sequence for up to ~10,000 spikes per synapse; ~10,000 synapses per node;
            # ~4,290,000 nodes per cell; ~4,290,000 cell in a network
            self.randObj.uniform(0, 1)
        else:  # if this synapse has already been stochastic before, this restarts its random number generator
            self.randObj.seq(self.cell.gid*1e4+1)
        syn = getattr(h, self._stochastic_type)(self.node.sec(self._loc))
        self._syn[self._stochastic_type] = {'target': syn}
        self._init_netcon(self._stochastic_type, delay=0.)
        self.target(self._stochastic_type).setRandObjRef(self.randObj)

    def target(self, target):
        """
        Returns the hoc object for the synaptic mechanism of the specified type
        :param target: str
        :return: :class:'h.HocObject'
        """
        if target in self._syn:
            return self._syn[target]['target']
        else:
            raise KeyError('Synapse type: {} not found at a synapse in {}'.format(target, self._node.name))

    def _init_netcon(self, target, delay=None, weight=None, threshold=None):
        """
        Appropriately initializes new netcon object, depending on whether the current source dictionary specifies a
        hocObject without a section, or a reference variable contained within a section.
        :param target: str
        :param delay = float
        :param weight = float
        :param threshold = float
        """
        if target in self._syn:
            source = self._source['object']
            syn = self._syn[target]
            if weight is None:
                weight = self.weight
            if threshold is None:
                threshold = self.threshold
            if delay is None:
                delay = self.delay
            if 'node' in self._source:
                node = self._source['node']
                syn['netcon'] = h.NetCon(source, syn['target'], sec=node.sec)
            else:
                syn['netcon'] = h.NetCon(source, syn['target'])
            syn['netcon'].delay = delay
            syn['netcon'].weight[0] = weight
            syn['netcon'].threshold = threshold
        else:
            raise KeyError('Synapse type: {} not found at a synapse in {}'.format(target, self._node.name))

    def netcon(self, target):
        """
        Returns the hoc network connection linking the synaptic mechanism of the specified type to a source of spikes.
        :param target: str
        :return: :class:'h.NetCon'
        """
        if target in self._syn:
            return self._syn[target]['netcon']
        else:
            raise KeyError('Synapse type: {} not found at a synapse in {}'.format(target, self._node.name))

    def change_source(self, source=None, node=None, param=None, loc=0.5):
        """
        In order to change the source of a synapse from the default VecStim object to a membrane potential or artificial
        cell, all the netcon objects must be deleted and replaced with new ones. Preserves previously set values for
        delay, weight, and threshold for all synaptic mechanisms associated with this synapse.
        :param source: hoc artificial cell or otherwise hoc object not associated with a section
        :param node: :class:'SHocNode'
        :param param: str corresponding to range variable in section
        :param loc: float
        """
        netcon_dict = {}
        if param is None:
            param = '_ref_v'
        for target in self._syn:
            netcon_dict[target] = {}
            netcon_dict[target]['delay'] = self.netcon(target).delay
            netcon_dict[target]['weight'] = self.netcon(target).weight[0]
            netcon_dict[target]['threshold'] = self.netcon(target).threshold
        if source is None:
            if node is None:
                raise Exception('A source or reference node must be provided to establish a new synaptic connection.')
            else:
                del self._source
                self._source = {'object': getattr(node.sec(loc), param), 'node': node}
        else:
            del self._source
            self._source = {'object': source}
        if self._stochastic:
            del self._syn[self._stochastic_type]['netcon']
            delay = netcon_dict[self._stochastic_type]['delay']
            weight = netcon_dict[self._stochastic_type]['weight']
            threshold = netcon_dict[self._stochastic_type]['threshold']
            self._init_netcon(self._stochastic_type, delay=delay, weight=weight, threshold=threshold)
        else:
            for target in (target for target in self._syn if not target == self._stochastic_type):
                del self._syn[target]['netcon']
                delay = netcon_dict[target]['delay']
                weight = netcon_dict[target]['weight']
                threshold = netcon_dict[target]['threshold']
                self._init_netcon(target, delay=delay, weight=weight, threshold=threshold)

    def get_stochastic(self):
        """
        Returns the value of an internal variable indicating if this synapse has a stochastic filter for spikes.
        :return: int in [0, 1]
        """
        return self._stochastic

    def set_stochastic(self, value):
        """
        Turns on or off stochastic filtering of spikes, preserving delay, weight, and threshold for all synaptic
        mechanisms associated with this synapse.
        :param value: int in [0, 1]
        """
        if not (value == self._stochastic):
            self._stochastic = value
            if value:
                self._init_stochastic()
                for target in (target for target in self._syn if not target == self._stochastic_type):
                    delay = self.netcon(target).delay
                    weight = self.netcon(target).weight[0]
                    threshold = self.netcon(target).threshold
                    del self._syn[target]['netcon']
                    self._syn[target]['netcon'] = h.NetCon(self.target(self._stochastic_type), self.target(target))
                    self.netcon(target).delay = delay
                    self.netcon(target).weight[0] = weight
                    self.netcon(target).threshold = threshold
            else:
                for target in (target for target in self._syn if not target == self._stochastic_type):
                    delay = self.netcon(target).delay
                    weight = self.netcon(target).weight[0]
                    threshold = self.netcon(target).threshold
                    del self._syn[target]['netcon']
                    self._init_netcon(target, delay=delay, weight=weight, threshold=threshold)
                del self._syn[self._stochastic_type]

    stochastic = property(get_stochastic, set_stochastic)

    def get_delay(self):
        """
        Returns the default value of the time delay (ms) between spike and activation for this synapse.
        :return: int or float
        """
        return self._delay

    def set_delay(self, value):
        """
        Changes the value of the time delay (ms) between spike and activation for all synaptic mechanisms associated
        with this synapse, except self._stochastic_type, which retains its current value until set manually.
        :param value: int or float
        """
        self._delay = value
        for target in (target for target in self._syn if not target == self._stochastic_type):
            self.netcon(target).delay = value

    delay = property(get_delay, set_delay)

    def get_weight(self):
        """
        Returns the default value of the activation weight for this synapse.
        :return: float
        """
        return self._weight

    def set_weight(self, value):
        """
        Changes the value of the activation weight for all synaptic mechanisms associated with this synapse.
        :param value: float
        """
        self._weight = value
        for target in (target for target in self._syn):
            self.netcon(target).weight[0] = value

    weight = property(get_weight, set_weight)

    def get_threshold(self):
        """
        Returns the value of the activation threshold for this synapse.
        :return: float
        """
        return self._threshold

    def set_threshold(self, value):
        """
        Changes the value of the activation threshold for this synapse.
        :param value: float
        """
        self._threshold = value
        for target in (target for target in self._syn):
            self.netcon(target).threshold = value

    threshold = property(get_threshold, set_threshold)

    @property
    def source(self):
        """
        Returns the hocObject currently being used as a source.
        :return: :class:'hocObject'
        """
        return self._source['object']

    @property
    def cell(self):
        """
        Returns the cell containing this synapse.
        :return: :class:'HocCell'
        """
        return self._cell

    @property
    def node(self):
        """
        Returns the node containing this synapse.
        :return: :class:'SHocNode'
        """
        return self._node

    @property
    def loc(self):
        """
        Returns the location along the hoc section containing this synapse. For convenience, if the synapse is
        contained in a spine_head, this property method returns the location along the branch section where the
        spine_neck is connected.
        :return: int or float
        """
        if self.node.type == 'spine_head':
            self.node.parent.sec.push()
            loc = h.parent_connection()
            h.pop_section()
            return loc
        else:
            return self._loc
