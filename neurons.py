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
import copy
import specify_cells as spec
import btmorph
import mech_dicts
import neuron_reduce
import math
import functions
import uuid
import datetime

class Cell(object):
    '''
    '''
    def __init__(self,**kwargs):
        '''
        '''
        self.init_parameters_and_cell_name(**kwargs) 
    
    def init_parameters_and_cell_name(self,**kwargs):
        '''
        '''
        # get parameters
        #-----------------------------------------------
        if 'P' in kwargs and kwargs['P'] is not None:
            self.P = kwargs['P']
        else:
            assert hasattr(self, 'P'), 'parameter object required to instantiate cell'
        # get cell name
        #-----------------------------------------------
        # from arguments
        if 'name' in kwargs:
            self.name = kwargs['name']
        # else from parameter dictionary
        elif hasattr(self, 'P') and 'name' in self.P.p:
            self.name=self.P.p['name']
        # else generate cell name
        else:
            # get class name
            class_name = self.__class__.__name__
            # if cell number is specified as argument, use this as id
            if 'cell_number' in kwargs and kwargs['cell_number'] is not None:
                cell_id = str(kwargs['cell_number'])
            else:
                # otherwise generate unique id
                cell_id = self._generate_cell_id()
            # set name
            self.name = '_'.join([class_name, cell_id])
        # update cell name and type in parameter dictionary
        #--------------------------------------------------
        self.P.p['cell_name'] = self.name
        self.P.p['cell_type'] = self.__class__.__name__

    def geometry(self, p):
        '''
        '''
        geo={}
        syns={}
        return geo, syns

    def mechanisms(self, p):
        '''
        '''
        pass

    def _generate_cell_id(self, ):
        '''
        '''
        # create unique identifier for each trial
        uid = str(uuid.uuid1().int)[-5:]
        now = datetime.datetime.now()
        trial_id = '-'.join(['{:04d}'.format(now.year), '{:02d}'.format(now.month), '{:02d}'.format(now.day), '{:02d}'.format(now.hour), '{:02d}'.format(now.minute), '{:02d}'.format(now.second), '{:02d}'.format(now.microsecond), uid])
        return trial_id

    def _seg_distance(self, geo):
        """ calculate distance from soma of each segment and store in parameter dictionary

        ==Args==
        -geo  : geometry structure as geo[tree][section][segment]

        ==Out==
        -p['seg_dist']  : structure containing the distance of each segment from the soma, organized in the same way as geo: p['seg_dist'][tree][section][segment distance]

        ==Updates==
        -'seg_dist' is added to p

        ==Comments==
        """
        # set soma as origin for distance measurements
        h.distance(sec = geo['soma'][0])

        self.p['seg_dist']={}
        
        # iterate over trees
        for tree_key,tree in geo.iteritems():

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
        print 'creating morpho object'
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
        # for tree_key, tree in geo.iteritems():
        #     for sec_i, sec in enumerate(tree):
        #         sref = h.SectionRef(sec=sec)
        #         root = sref.root
        #         break

        # create new secton list
        # nrn_sec_list = h.SectionList()
        # add all seection to list, starting from root
        # nrn_sec_list.wholetree()

        # copy nrn section list as a python list
        # sec_list = []
        # for sec_i_temp, sec_temp in enumerate(nrn_sec_list):
        #     print sec_temp.name()
        #     sec_list.append(sec_temp)
        sec_list = []
        seg_list = []
        seg_list_idx = []
        idx = -1
        for tree_key, tree in geo.iteritems():
            for sec_i, sec in enumerate(tree):
                # nrn_sec_list.append(sec)
                sec_list.append(sec)
                seg_list.append([])
                seg_list_idx.append([])
                for seg_i, seg in enumerate(sec):
                    # add to total segments counter and store in lists
                    idx+=1
                    # copy index count to prevent overwrite during loop
                    idx_count = copy.copy(idx)
                    seg_list[-1].append(seg)
                    seg_list_idx[-1].append(idx_count)


                # print sec.name()
        # for sec_i_temp, sec_temp in enumerate(nrn_sec_list):
        #     sec_list.append(sec_temp)
        #     print sec_temp.name()
        # print sec_list
        # nested list for storing segment objects [section_number][segment_number]
        # seg_list= []
        # nested list for storing segment indices [section number][segment number]
        # seg_list_idx = []
        # nested list for storing index of parent segment, matches seg_list_idx dimesions, [section_number][segment_number]
        parent_list_idx = []
        # keep track of total segment number
        idx = -1
        # print sec_list
        # iterate through sections in list
        for sec_i, sec in enumerate(sec_list):
            # keep track of the root section
            is_root=False


            # add section dimension to each list
            # seg_list.append([])
            # seg_list_idx.append([])
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
                # seg_list[sec_i].append(seg)
                # seg_list_idx[sec_i].append(idx_count)

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

                            # print sec.name(), sec_local.name()
                            # if section name and segment index match
                            if (sec.name() in sec_local.name()) and (seg_i == seg_i_local):
                            # if (sec == sec_local) and (seg_i == seg_i_local):
                            # if (sec.name() == sec_local.name()) and (seg_i == seg_i_local):
                                # print 'true'
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

    def _get_loc(self, geo, segment=None, hoc_name=None, seg_x=None):
        '''
        '''
        if segment is not None:
            hoc_name = segment.sec.hname()
            if '.' in hoc_name:
                hoc_name = hoc_name.split('.')[-1]
            seg_x = segment.x
        else:
            assert hoc_name is not None and seg_x is not None, 'hoc_name and seg_x are both requried if no segment object is passed'
            if '.' in hoc_name:
                hoc_name = hoc_name.split('.')[-1]


        for tree_key, tree in geo.iteritems():
            for sec_i, sec in enumerate(tree):
                if hoc_name in sec.hname():
                    for seg_i, seg in enumerate(sec):
                        if seg==sec(seg_x):
                            location = (tree_key, sec_i, seg_i)

        return location

    def _seg_location(self, sec):
        """ given a neuron section, output the 3d coordinates of each segment in the section

        ouput is a nested list as [xyz dimension][segment number], with x,y, z dimensions listed in that order

        """
        # number of 3d points in section
        tol =.001
        # specifies the location of each segment so that n3d can be used
        # define_shape messes up 3d locations if they arlready specified.  not actually sure what this function does
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
                seg_y.append( y[ node_i[ 0]])
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

    def _get_trees_from_h(self, h):
        ''' HACK
        '''
        trees=[]
        for sec in h.allsec():
            name = sec.name()
            if '[' in name:
                tree_name = name.split('[')[-2].split('.')[-1]
            else:
                tree_name=name
            trees.append(tree_name)
        trees =  list(set(trees))
        return trees

    def _load_geometry(self, geo_filename, geo_filetype='hoc', hoc_name=None, template=False, tree_name_dict=None, **kwargs):
        ''' load cell geometry from hoc file which specifies sections and connections. 
        ==Args==
        -geo_filename : str : hoc file specifying the geometry
        -geo_filetype : str : type of file (for now this only handles hoc # FIXME)
        -hoc_name : str : name at hoc level for the loaded cell template
        -template : bool : if True the hoc file is converted to a template and the geometry is loaded as an instance of the template
        -tree_name_dict : dict : mapping between hoc level tree names (specified in the hoc file) and python level tree names {hoc_name:python name}. e.g. {'apical_dendrite':'apical'}
        ==Return==
        -geo : dict : geometry dictionary as {tree name:hoc list of sections belonging to that tree}
        ==Update==
        ==Comments==
        '''
        # specify hoc level name for the template, otherwise the default is the filename
        if hoc_name is None:
            template_name = geo_filename.split('.')[0].split('_')[0]
        else: 
            template_name= hoc_name
        # dictionary for storing geometry section and segment objects
        
        # load hoc geometry from template
        #--------------------------------------------------------------------
        if geo_filetype == 'hoc' and template:
            # report progress
            print 'loading cell type:', self.__class__.__name__, ', with geometry:', geo_filename
            # convert hoc code to template if it is not already
            self._convert_geo_to_template(geo_filename)
            # load template
            h.load_file(geo_filename)
            # get template at hoc level
            template = getattr(h, template_name)
            # create an instance of the template
            geo = template()
            # instance = template()
        #     print instance.soma[0].name()
        #     # get list of trees in the template
        #     trees = self._get_trees_from_h(h=h)
        #     print trees
        #     # create geo structure 
        #     geo = {}
        #     # iterate over all tree types that current exist in hoc
        #     for tree in trees:
        #         # check that the current cell instance has this tree
        #         if hasattr(instance, tree):
        #             # get the hoc tree
        #             hoc_tree=getattr(instance, tree)
        #             # if the hoc tree is in the dictionary of new tree names, change the tree name at the python level
        #             if tree_name_dict is not None and tree in tree_name_dict:
        #                 tree_key = tree_name_dict[tree]
        #             else:
        #                 tree_key= tree
        #             geo[tree_key] = hoc_tree
        # h.distance(sec = geo['soma'][0])

        return geo

    def _convert_geo_to_template(self, geo_filename, hoc_name=None):
        '''
        '''
        if hoc_name is None:
            name = geo_filename.split('.')[0].split('_')[0]
        else: 
            name= hoc_name
        with open(geo_filename, 'r') as file:
            # list of strings for each line of code
            hocfile = file.readlines()
        
        begin = 'begintemplate '+name+'\n'
        end = 'endtemplate '+name
        all_text = ''.join(hocfile)
        if 'begintemplate' not in all_text:
            hocfile.insert(0, begin)
        if 'endtemplate' not in all_text:
            hocfile.append(end)

        with open(geo_filename, 'w') as file:
            file.writelines(hocfile)

    def _convert_template_to_geo(self, geo_filename):
        '''
        '''
        # name = geo_filename.split('.')[0].split('_')[0]
        with open(geo_filename, 'r') as file:
            # list of strings for each line of code
            hocfile = file.readlines()
        
        begin_key = 'begintemplate '
        end_key = 'endtemplate'
        all_text = ''.join(hocfile)
        remove_lines = [line for i, line in enumerate(hocfile) if (begin_key in line or end_key in line) ]
        for line in remove_lines:
            hocfile.remove(line)

        with open(geo_filename, 'w') as file:
            file.writelines(hocfile)

    def _syns_to_list(self, syns):
        '''
        '''
        syn_list=[]
        for tree_key, tree in syns.iteritems():
            for sec_i, sec in enumerate(tree):
                for seg_i, seg in enumerate(sec):
                    if type(seg)==list:
                        for syn_type_i, syn_type in enumerate(seg):
                            for syn_key, syn in syn_type.iteritems():
                                syn_dict = {
                                'tree_key':tree_key,
                                'sec_num':sec_i,
                                'seg_num':seg_i,
                                'syn_i':syn_type_i,
                                'syn_key':syn_key,
                                'syn':syn,
                                'sec_object':syn.get_segment().sec,
                                'seg_loc':syn.get_segment().x,
                                }
                                syn_list.append(syn_dict)

                    elif type(seg)==dict:
                        for syn_key, syn in seg.iteritems():
                            syn_dict = {
                            'tree_key':tree_key,
                            'sec_num':sec_i,
                            'seg_num':seg_i,
                            # 'syn_i':syn_type_i,
                            'syn_key':syn_key,
                            'syn':syn,
                            'sec_object':syn.get_segment().sec,
                            'seg_loc':syn.get_segment().x
                            }
                            syn_list.append(syn_dict)

        return syn_list

    def _lambda_f(self, sec, f=100.):
        """
        Calculates the AC length constant for the given section at the frequency f
        Used to determine the number of segments per hoc section to achieve the desired spatial and temporal resolution
        :param sec : :class:'h.Section'
        :param f : int
        :return : int
        """
        diam = sec(0.5).diam
        Ra = sec.Ra
        cm = sec.cm
        return 1e5*math.sqrt(diam/(4.*math.pi*f*Ra*cm))

    def _d_lambda_nseg(self, sec, lam=0.1, f=100.):
        """
        The AC length constant for this section and the user-defined fraction is used to determine the maximum size of each
        segment to achieve the desired spatial and temporal resolution. This method returns the number of segments to set
        the nseg parameter for this section. For tapered cylindrical sections, the diam parameter will need to be
        reinitialized after nseg changes.
        :param sec : :class:'h.Section'
        :param lam : int
        :param f : int
        :return : int
        """
        L = sec.L
        return int((L/(lam*self._lambda_f(sec, f))+0.9)/2)*2+1

    def fixnseg(self,):
        '''
        '''
        assert hasattr(self, 'geo'), 'geometry attribute required to find sections and fix nseg'
        for tree_key, tree in self.geo.iteritems():
            for sec_i, sec in enumerate(tree):
                nseg = self._d_lambda_nseg(sec=sec)
                sec.nseg=nseg

    def _get_length_constant(self, sec):
        ''' return length constant (um) of a given section
        '''
        a = sec.diam*1e-4
        L = sec.L*1e-4
        ri = sec.Ra
        rm = 1./sec.g_pas
        area = np.pi*a**2
        circumference = 2.*np.pi*a
        rho_m = rm/circumference
        rho_i = ri/area
        length_constant = np.sqrt(rho_m/rho_i)*1e4
        return length_constant

    def _get_polarization_analytic(self, field, sec, loc):
        '''
        '''
        lam = self._get_length_constant(sec)
        L = sec.L
        L_norm = L/lam
        loc_norm = loc*L/lam
        polarization = .001*field*lam*np.sinh(loc_norm)/np.cosh(L_norm)
        return polarization

    def _get_mirror_estimate_from_polarization(self, vm=None, l_lambda=None, locations=None, polarization_map=None, reduced_cell=None, zero_mean=True, map_method=None, **kwargs):
        ''' mirror estimate is used to 
        '''
        def _from_map(reduced_cell, polarization_map, method=None, apical_factor=2., basal_factor=1.,apical_offset=2., basal_offset=1., **kwargs):
            ''' given a polarization map from segments in a reduced neuron a list of corresponding polarizations in an original neuron, return the voltages that the electric should attempt to reproduce
            ==Args==
            reduced_cell:
            polarization_map: dictionary with keys as location tuples in the reduced model (tree_key, sec_num, seg_num) and values as lists of polarizations recorded from a full model
            method:str:specify how to choose polarization values to replicate with mirror estimate
            ==Return==
            vm:1D array:polarization values to replicate with mirror estimate
            l_lambda:1D array:length constants for each segment in vm
            locs:list:locations for each segment in vm (tree_key, sec_num, seg_num)
            ==Updates==
            ==Comments==
            '''
            if method=='from_mean' or method is None:
                polarization_map_mean = {}
                for key, val in polarization_map.iteritems():
                    polarization_map_mean[key] = np.mean(val)
                vm = []
                l_lambda = []
                locs = []
                for key, val in polarization_map_mean.iteritems():
                    tree_key = key[0]
                    sec_num = key[1]
                    sec = reduced_cell.geo[tree_key][sec_num]
                    # if tree_key == 'dend':
                    #     sec = sec[key[1]]
                    lamb = self._get_length_constant(sec)
                    v = val
                    vm.append(v)
                    l_lambda.append(lamb)
                    locs.append(key)
                vm=np.array(vm)
                l_lambda=np.array(l_lambda)
                return vm, l_lambda, locs

            elif method=='from_median' or method is None:
                polarization_map_median = {}
                for key, val in polarization_map.iteritems():
                    polarization_map_median[key] = np.median(val)
                vm = []
                l_lambda = []
                locs = []
                for key, val in polarization_map_median.iteritems():
                    tree_key = key[0]
                    sec_num = key[1]
                    sec = reduced_cell.geo[tree_key][sec_num]
                    # if tree_key == 'dend':
                    #     sec = sec[key[1]]
                    lamb = self._get_length_constant(sec)
                    v = val
                    vm.append(v)
                    l_lambda.append(lamb)
                    locs.append(key)
                vm=np.array(vm)
                l_lambda=np.array(l_lambda)
                return vm, l_lambda, locs

            elif method == 'linear_fit' or method=='mean_fit' or method=='end_fit':
                x = {}
                y = {}
                for key,val in polarization_map.iteritems():
                    
                    sec_key = (key[0], key[1])
                    if sec_key not in x:
                        x[sec_key]=[]
                        y[sec_key]=[]

                    if method=='mean_fit':
                        x[sec_key].append(key[2])
                        y[sec_key].append(np.mean(val))
                    elif method=='linear_fit' or method=='end_fit':
                        for val_i, val in enumerate(polarization_map[key]):
                            x[sec_key].append(key[2])
                            y[sec_key].append(val)

                vm=[]
                l_lambda=[]
                locs=[]
                for sec_key in x:
                    if method=='end_fit':
                        # get sign of y
                        if len(y[sec_key])>1:
                            max_x = np.max(x[sec_key])
                            min_x = np.min(x[sec_key])
                            max_y_vals = [_val for _i, _val in enumerate(y[sec_key]) if x[sec_key][_i]==max_x]
                            print 'max_y_vals ', max_y_vals 
                            min_y_vals = [_val for _i, _val in enumerate(y[sec_key]) if x[sec_key][_i]==min_x]
                            print 'min_y_vals ', min_y_vals 
                            if np.sign(np.mean(max_y_vals))>0:
                                max_y = np.max(max_y_vals)
                            else:
                                max_y = np.min(max_y_vals)   
                            if np.sign(np.mean(min_y_vals))>0:
                                min_y = np.min(min_y_vals)
                            else:
                                min_y = np.max(min_y_vals)
                            x[sec_key]=[min_x, max_x]
                            y[sec_key]=[min_y, max_y]

                    tree_key = sec_key[0]
                    sec_num = sec_key[1]
                    sec = reduced_cell.geo[tree_key][sec_num]
                    n_seg = sec.nseg
                    lamb = self._get_length_constant(sec)
                    # get new voltage values from linear fit
                    #-------------------------------------------
                    x_array = np.array(x[sec_key])
                    y_array = np.array(y[sec_key])
                    A_mat = np.vstack([x_array, np.ones(len(x_array))]).T
                    
                    m, c = np.linalg.lstsq(A_mat, y_array)[0]
                    if sec_key[0]=='apic':
                        if apical_factor is not None:
                            m = apical_factor*m
                        if apical_offset is not None:
                            c = apical_offset*c
                    if sec_key[0]=='dend':
                        if basal_factor is not None:
                            m = basal_factor*m
                        if basal_offset is not None:
                            c = basal_offset*c
                        m = basal_factor*m

                    x_new = np.arange(n_seg)
                    y_new = x_new*m + c
                    for _i, seg_num in enumerate(x_new):
                        loc = (tree_key, sec_num, int(seg_num))
                        vm.append(y_new[_i])
                        l_lambda.append(lamb)
                        locs.append(loc)
                vm = np.array(vm)
                l_lambda=np.array(l_lambda)

                return vm, l_lambda, locs

        if polarization_map is not None and reduced_cell is not None:
            vm, l_lambda, locations = _from_map(reduced_cell, polarization_map, method=map_method, **kwargs)
        else:
            assert vm is not None and l_lambda is not None, 'if polarization map is not provided, raw vm and l_lambda are required'

        def _mirror_estimate(vm, l_lambda, locations, output_type='dict', **kwargs):
            # calculate mirror estimate based on 
            l_lambda_sum = np.sum(l_lambda**-2)
            A = np.tile(l_lambda**-2, (len(vm), 1))
            I = np.identity(len(vm))
            A = A - I*l_lambda_sum
            A_inv = np.linalg.inv(A)

            # add row for zero mean
            new_row = np.tile(l_lambda**-2, (1,1)) #(1./float(A.shape[1]))*np.ones((1, A.shape[1]))
            A = np.append(A, new_row, axis=0)
            vm=np.append(vm, 0,)

            # v_ext = np.matmul(A_inv, vm*l_lambda_sum)
            # v_ext = np.linalg.solve(A, vm*l_lambda_sum)
            lstsq = np.linalg.lstsq(A, vm*l_lambda_sum)
            residuals = lstsq[1]
            rank = lstsq[2]
            v_ext = lstsq[0]
            v_ext = v_ext-np.mean(v_ext)
            if output_type=='dict':
                assert locations is not None, 'locations list is required to generate map of extracellular voltage at each segment'
                e_map = {}
                for loc_i, loc in enumerate(locations):
                    e_map[loc] = v_ext[loc_i]

                return e_map

            elif locations is not None:
                return zip(locations, v_ext)
            else:
                return v_ext
        e_map= _mirror_estimate(vm=vm, l_lambda=l_lambda, locations=locations, **kwargs)
        return e_map

    def _create_loc_map(self, original_cell, reduced_cell, seg_to_seg, original_to_reduced=False, **kwargs):
        ''' for each location in this reduced cell return a list of corresponding locations in the original cell
        '''
        loc_map={}
        # seg_to_seg maps new locations to original locations, but uses hoc section names, which are mapped to top level location scheme by _get_loc
        for key, val in seg_to_seg.iteritems():
            # top level location in original cell (tree_key, sec_num, seg_num)
            orig_loc = original_cell._get_loc(geo=original_cell.geo, hoc_name=key[0], seg_x=key[1])
            # top level location in reduced cell (tree_key, sec_num, seg_num)
            new_loc = reduced_cell._get_loc(geo=reduced_cell.geo, hoc_name=val[0], seg_x=val[1])
            # specify direction of mapping 
            # original -> reduced
            if original_to_reduced:
                if orig_loc not in loc_map:
                    loc_map[orig_loc]=[]
                loc_map[orig_loc].append(new_loc)
            # reduced -> original
            else:
                if new_loc not in loc_map:
                    loc_map[new_loc]=[]    
                loc_map[new_loc].append(orig_loc)
        return loc_map

class CellMigliore2005(Cell):
    """ pyramidal neuron based on Migliore et al. 2005

    An instance of this object will creates a cell (hoc objects) at the top level of the hoc interpreter using the hoc files in _init_geometry.  The .geo attribute contains a python mapping to these hoc objects.  The geo object is organized as geo['section tree'][section](segment location)

    the syns attribute creates a container for synapse objects that are added to each segment in the hoc cell.  syns is organized as syns['section tree']['synapse type'][section][segment number]
    
    """
    def __init__(self, geo_filename='geo5038804_template.hoc', **kwargs):
        '''
        '''
        # Parameter object hasn't been instantiated, instantiate
        #---------------------------------------------------------------
        if not hasattr(self, 'P'):
            # check kwargs
            #---------------
            if 'P' in kwargs and kwargs['P'] is not None:
                self.P = kwargs['P']
            # otherwise load default
            #-----------------------
            else:
                self.P = param.ParamMigliore2005()
         # initialize from base class
        #----------------------------------------------------------------
        super(CellMigliore2005, self).__init__(**kwargs)

        self.p = self.P.p
        # define geometry and insert synapses
        #-----------------------------------
        self.geo = self.geometry(p=self.P.p, geo_filename=geo_filename, **kwargs)
        self.fixnseg()
        # insert mechanisms 
        #--------------------
        self.mechanisms(p=self.P.p)
        # insert synapses
        #-----------------
        self.syns = self.insert_synapses(p=self.P.p)
        # create morpho object
        #--------------------------
        self.morpho = self._create_morpho(self.geo)
        self.p['morpho'] = self.morpho
        # creata seg_dist entry in p dict
        #---------------------------------
        self.p['seg_dist'] = self._seg_distance(geo=self.geo)

    # geometry geo5038804
    #-------------------------
    def geometry(self, p, geo_filename='geo5038804_template.hoc', **kwargs):
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
        h.load_file(geo_filename)  
        # h.load_file('geoc62564.hoc')
          
        # dictionary for storing geometry ['tree'][sec](seg location)
        geo = {}
        # if geometry file is a template, instantiate cell and get sections from instance
        if 'template' in geo_filename:
            template_name = geo_filename.split('.')[0]
            cell = getattr(h, template_name)()
            geo['soma'] = cell.soma
            geo['axon'] =  cell.axon
            geo['basal'] = cell.dendrite
            geo['apical_trunk'] = cell.user5
            geo['apical_tuft'] = cell.apical_dendrite
        # otherwise get sections directly from hoc
        else:
            geo['soma'] = h.soma
            geo['axon'] =  h.axon
            geo['basal'] = h.dendrite
            geo['apical_trunk'] = h.user5
            geo['apical_tuft'] = h.apical_dendrite

        # set temperature in hoc
        h.celsius = p['celsius']
        # set soma as origin for distance measurements
        h.distance(sec = geo['soma'][0])

        return geo

    def insert_synapses(self, p, **kwargs):
        '''
        '''
        print 'inserting synapses'
        assert hasattr(self, 'geo'), 'cell geometry must be loaded before inserting membrane mechanisms and synapses'
        if 'dendrites' in kwargs:
            dendrites=kwargs['dendrites']
        else:
            dendrites = ['basal', 'apical_trunk', 'apical_tuft']
        syns={}
        # loop over trees
        for tree_key,tree in self.geo.iteritems():
            # list to store synapse mechanisms
            syns[tree_key] = []
            # loop over sections in tree
            for sec_i,sec in enumerate(tree):

                # add dimension for each section
                syns[tree_key].append([])
                # only insert synapses in dendrites
                if tree_key in dendrites:
                    # iterate over segments 
                    for seg_i,seg in enumerate(sec):

                        # preallocate synapse lists
                        syns[tree_key][sec_i].append({'ampa':[],
                        'nmda':[],
                        'clopath':[]})
                        # iterate over synapse types
                        for syn_key,syn in syns[tree_key][sec_i][seg_i].iteritems():
                            # ampa
                            #------------------------------------------
                            if syn_key is 'ampa':
                                # adapting exponential synapse based on model in Varela et al. 1997
                                syns[tree_key][sec_i][seg_i][syn_key] = h.FDSExp2Syn_D3(sec(seg.x))
                                syns[tree_key][sec_i][seg_i][syn_key].f = p['f_ampa']
                                syns[tree_key][sec_i][seg_i][syn_key].tau_F = p['tau_F_ampa']
                                syns[tree_key][sec_i][seg_i][syn_key].d1 = p['d1_ampa']
                                syns[tree_key][sec_i][seg_i][syn_key].tau_D1 = p['tau_D1_ampa']
                                syns[tree_key][sec_i][seg_i][syn_key].d2 = p['d2_ampa']
                                syns[tree_key][sec_i][seg_i][syn_key].tau_D2 = p['tau_D2_ampa']
                                syns[tree_key][sec_i][seg_i][syn_key].d3 = p['d3_ampa']
                                syns[tree_key][sec_i][seg_i][syn_key].tau_D3 = p['tau_D3_ampa']

                                syns[tree_key][sec_i][seg_i][syn_key].tau1 = p['tau1_ampa']
                                syns[tree_key][sec_i][seg_i][syn_key].tau2 = p['tau2_ampa']


                            # nmda
                            #----------------------------
                            elif syn_key is 'nmda':
                                syns[tree_key][sec_i][seg_i][syn_key]= h.Exp2SynNMDA(sec(seg.x))
                                syns[tree_key][sec_i][seg_i][syn_key].tau1 = p['tau1_nmda']
                                syns[tree_key][sec_i][seg_i][syn_key].tau2 = p['tau2_nmda']
                                syns[tree_key][sec_i][seg_i][syn_i][syn_key].v0_block=p['v0_block_nmda']
                                syns[tree_key][sec_i][seg_i][syn_i][syn_key].alpha_vspom=p['alpha_vspom_nmda']
                            # clopath
                            #-------------------------------
                            elif syn_key is 'clopath':
                                # print syn_key
                                syns[tree_key][sec_i][seg_i][syn_key] = h.STDPSynCCNon(sec(seg.x))

        return syns

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
        assert hasattr(self, 'geo'), 'cell geometry must be loaded before inserting membrane mechanisms and synapses'
        dendrites = ['basal', 'apical_trunk', 'apical_tuft']
        # loop over trees
        for tree_key,tree in self.geo.iteritems():
            # loop over sections in tree
            for sec_i,sec in enumerate(tree):
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
                #-----------------------------------------------------
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
                # soma active biophysics
                #----------------------------------------------------------
                elif tree_key == 'soma':
                    # voltage gated sodium
                    sec.insert('na3')
                    sec.gbar_na3 = p['gna']*p['SOMAM']
                    sec.ar_na3 = p['gna_inact']
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
                    # high threshold calcium
                    # sec.insert('calH')
                    # sec.gcalbar_calH = p['gcalbar']
                    # sodium reversal potential 
                    sec.ena = p['ena']      
                    # potassium reversal potential 
                    sec.ek = p['ek']    
                # dendrites active biophysics
                #------------------------------------------------
                elif tree_key in dendrites:
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

                    # # L-type calcium channel
                    sec.insert('calH')
                    sec.gcalbar_calH = p['gcalbar']

                    # # R-type calcium channel
                    sec.insert('car')
                    sec.gcabar_car = p['gcabar_r']

                    # # L-type calcium channel
                    sec.insert('cat')
                    sec.gcatbar_cat = p['gcatbar_t']

                    # sodium reversal potential 
                    sec.ena = p['ena']
                    # potassium reversal potential
                    sec.ek = p['ek']        

                    # mechanisms that vary with distance from soma
                    #--------------------------------------------------------
                    # loop over segments
                    for seg_i,seg in enumerate(sec):
                        # distance from soma
                        #---------------------
                        seg_dist = h.distance(seg.x,sec=sec)
                        # sodium
                        #--------------------
                        if abs(p['dgna']*seg_dist)<p['gna']:
                            seg.gbar_na3 = p['gna'] + p['dgna']*seg_dist
                        else:
                            seg.gbar_na3 = 0.
                        # h current
                        #----------------------
                        if seg_dist < p['ghd_cutoff_distance']:
                            seg.ghdbar_hd = p['ghd']*(1+p['ghd_grad']*seg_dist/100.)
                        else:
                            seg.ghdbar_hd = p['ghd']*(1+p['ghd_grad']*p['ghd_cutoff_distance']/100.)
                        # A-type potassium
                        #-----------------------
                        # distal
                        #---------
                        if seg_dist > 100.: 
                            seg.vhalfl_hd = p['vhalfl_hd_dist']
                            seg.vhalfl_kad = p['vhalfl_kad']
                            seg.vhalfn_kad = p['vhalfn_kad']
                            if seg_dist < p['ka_cutoff_distance']:
                                seg.gkabar_kad = p['KMULT']*(1+p['ka_grad']*seg_dist/100.)
                            else:
                                seg.gkabar_kad = p['KMULT']*(1+p['ka_grad']*p['ka_cutoff_distance']/100.)

                            # T-type calcium
                            #-------------------------
                            seg.gcatbar_cat = p['gcatbar_t']*seg_dist/350.
                        # proximal
                        #----------
                        else:   
                            seg.vhalfl_hd = p['vhalfl_hd_prox']
                            seg.vhalfl_kap = p['vhalfl_kap']
                            seg.vhalfn_kap = p['vhalfn_kap']
                            seg.gkabar_kap = p['KMULTP']*(1+p['ka_grad']*seg_dist/100.)
                            # T-type calcium
                            #-------------------------
                            seg.gcatbar_cat = 0.

                        if seg_dist > 400:
                            seg.gcalbar_calH = p['gcalbar']
                        else:
                            seg.gcalbar_calH =0

    def mechanisms_passive(self, p):
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

                # voltage gated sodium      
                sec.insert('na3')
                sec.gbar_na3 = p['gna']
                sec.ar_na3 = p['gna_inact'] 

                # delayed rectifier potassium   
                sec.insert('kdr')
                sec.gkdrbar_kdr = p['gkdr']

                # dendrites active biophysics
                if ((tree_key == 'basal') or 
                (tree_key == 'apical_trunk') or 
                (tree_key == 'apical_tuft')):
                     
                    # mechanisms that vary with distance from soma
                    # loop over segments
                    # print 'creating synapses'
                    for seg_i,seg in enumerate(sec):
                        # print 'creating synapses'
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

    def insert_synapses_reduced_model(self, p, **kwargs):
        ''' allows for insertion of multiple synapses onto each segment
        '''
        print 'inserting synapses'
        assert hasattr(self, 'geo'), 'cell geometry must be loaded before inserting membrane mechanisms and synapses'
        # names of dendritic trees (in geo dict) to insert synapses onto 
        #---------------------------------------------------------------
        if 'dendrites' in kwargs:
            dendrites=kwargs['dendrites']
        else:
            dendrites = ['basal', 'apical_trunk', 'apical_tuft']
        # number of synapses to insert onto each dendritic segment. default is one.  if short term depression is included, separate synapses are needed so that short term plasticity does not interfere betwixt pathways
        #------------------------------------------------------------------------
        if 'syns_per_seg' in kwargs:
            syns_per_seg=kwargs['syns_per_seg']
        else:
            syns_per_seg=1     
        syns={}
        # loop over trees
        for tree_key,tree in self.geo.iteritems():
            # only insert synapses in dendrites
            if tree_key in dendrites:
                # list to store synapse mechanisms
                syns[tree_key] = []
                # loop over sections in tree
                for sec_i,sec in enumerate(tree):
                    # add dimension for each section
                    syns[tree_key].append([])
                    # iterate over segments 
                    for seg_i,seg in enumerate(sec):
                        syns[tree_key][sec_i].append([])
                        # iterate over synapses per segment
                        for syn_i in range(syns_per_seg):
                            # preallocate synapse lists
                            syns[tree_key][sec_i][seg_i].append({'ampa':[],
                            'nmda':[],
                            'clopath':[]})
                            # iterate over synapse types
                            for syn_key,syn in syns[tree_key][sec_i][seg_i][syn_i].iteritems():
                                # ampa
                                #------------------------------------------
                                if syn_key is 'ampa':
                                    # adapting exponential synapse based on model in Varela et al. 1997
                                    syns[tree_key][sec_i][seg_i][syn_i][syn_key] = h.FDSExp2Syn_D3(sec(seg.x))
                                    syns[tree_key][sec_i][seg_i][syn_i][syn_key].f = p['f_ampa']
                                    syns[tree_key][sec_i][seg_i][syn_i][syn_key].tau_F = p['tau_F_ampa']
                                    syns[tree_key][sec_i][seg_i][syn_i][syn_key].d1 = p['d1_ampa']
                                    syns[tree_key][sec_i][seg_i][syn_i][syn_key].tau_D1 = p['tau_D1_ampa']
                                    syns[tree_key][sec_i][seg_i][syn_i][syn_key].d2 = p['d2_ampa']
                                    syns[tree_key][sec_i][seg_i][syn_i][syn_key].tau_D2 = p['tau_D2_ampa']
                                    syns[tree_key][sec_i][seg_i][syn_i][syn_key].d3 = p['d3_ampa']
                                    syns[tree_key][sec_i][seg_i][syn_i][syn_key].tau_D3 = p['tau_D3_ampa']

                                    syns[tree_key][sec_i][seg_i][syn_i][syn_key].tau1 = p['tau1_ampa']
                                    syns[tree_key][sec_i][seg_i][syn_i][syn_key].tau2 = p['tau2_ampa']

                                # nmda
                                #----------------------------
                                elif syn_key is 'nmda':
                                    syns[tree_key][sec_i][seg_i][syn_i][syn_key]= h.Exp2SynNMDA(sec(seg.x))
                                    syns[tree_key][sec_i][seg_i][syn_i][syn_key].tau1 = p['tau1_nmda']
                                    syns[tree_key][sec_i][seg_i][syn_i][syn_key].tau2 = p['tau2_nmda']
                                    syns[tree_key][sec_i][seg_i][syn_i][syn_key].v0_block=p['v0_block_nmda']
                                    syns[tree_key][sec_i][seg_i][syn_i][syn_key].alpha_vspom=p['alpha_vspom_nmda']
                                # clopath
                                #-------------------------------
                                elif syn_key is 'clopath':
                                    # print syn_key
                                    syns[tree_key][sec_i][seg_i][syn_i][syn_key] = h.STDPSynCCNon(sec(seg.x))

        return syns

class CellMigliore2005Reduced(CellMigliore2005):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        # Parameter object hasn't been instantiated, instantiate
        #---------------------------------------------------------------
        if not hasattr(self, 'P'):
            # check kwargs
            #---------------
            if 'P' in kwargs and kwargs['P'] is not None:
                self.P = kwargs['P']
            # otherwise load default
            #-----------------------
            else:
                self.P = param.ParamMigliore2005()
                kwargs['P']=self.P
        # increase sodium conductance in soma to account for loss of axon inreduced model
        # self.P.p['AXONM']=1000.
        # self.P.p['SOMAM']=1.5
        # self.P.p['gna_inact']=1
        # self.P.p['tau1_ampa']=.2# rise time constant (ms)
        # self.P.p['tau2_ampa']=10.   # decay time constant   (ms)

        # load full cell
        #-----------------
        super(CellMigliore2005Reduced, self).__init__(**kwargs)
        # get syns list
        #------------------
        _synapses = self._syns_to_list(self.syns)
        synapses = [ _syn['syn'] for _syn in _synapses]
        # create temporary netcons list
        #-------------------------------
        # temporary netstim object
        _netstim = h.NetStim()
        netcons=[]
        for syn in synapses:
            netcons.append(h.NetCon(_netstim, syn, ))
        # parameter dictionary
        #-------------------
        self.p=self.P.p
        # reduce cell
        #----------------
        self.reduce(synapses=synapses, netcons=netcons, **kwargs)
        # set 3d geometry of reduced cell
        #-----------------------
        self.set_hoc_pt3d_geometry(**kwargs)
        # create geometry object
        #------------------------
        self.geo = self.create_geo()
        # insert synapses
        #----------------------------
        # load facilitation/depression parameters
        self.P = self.P._load_fd_parameters(P=self.P, filename='Data/fd_parameters.pkl')
        self.P.p['d3_ampa']=1.#0.998
        # self.P.p['d1_ampa']=1#.98
        # self.P.p['d2_ampa']=1#.45*.12#1.
        # self.P.p['tau_F_ampa']=1#50.#1.

        # self.P.p['f_ampa']=0.
        # self.P.p['tau_D1_ampa']=200.
        # self.syns = self.insert_synapses(self.p, dendrites=['apic','dend'])
        self.syns = self.insert_synapses_reduced_model(p=self.P.p, dendrites=['apic','dend'], **kwargs)
        # create morpho object
        #--------------------------
        self.morpho = self._create_morpho(self.geo)
        self.p['morpho'] = self.morpho
        # creata seg_dist entry in p dict
        #---------------------------------
        # print 'h.y3d(2, sec=self.reduced_cell.dend[0])',h.y3d(2, sec=self.reduced_cell.dend[0])
        self.p['seg_dist'] = self._seg_distance(geo=self.geo)
        # print 'h.y3d(2, sec=self.reduced_cell.dend[0])',h.y3d(2, sec=self.reduced_cell.dend[0])

    def reduce(self, **kwargs):
        '''
        '''
        # create hoccell object to pass to neuron_reduce
        #-------------------------------------------------------------------
        self.hoccell = spec.HocStyleCell()
        self.hoccell.from_geo(self)
        # list all synapses
        #-------------------
        if 'synapses' in kwargs and type(kwargs['synapses'])==list:
            self.full_synapses=kwargs['synapses']
        else:
            self.full_synapses=[]
        # list netcons
        #--------------------
        if 'netcons' in kwargs and type(kwargs['netcons'])==list:
            self.full_netcons=kwargs['netcons']
        else:
            self.full_netcons=[]
        # run neuron_reduce
        #--------------------------------------------------------------------
        # apply Neuron_Reduce to simplify the cell
        # self.reduced_cell, self.reduced_synapses, self.reduced_netcons = neuron_reduce.subtree_reductor(
        #     original_cell=self.hoccell,
        #     synapses_list=self.full_synapses,
        #     netcons_list=self.full_netcons,
        #     reduction_frequency=0, 
        #     total_segments_manual=-1)

        self.reduced_cell, self.reduced_synapses, self.reduced_netcons, self.seg_to_seg_text = neuron_reduce.subtree_reductor(
            original_cell=self.hoccell,
            synapses_list=self.full_synapses,
            netcons_list=self.full_netcons,
            reduction_frequency=0, 
            total_segments_manual=-1,
            return_seg_to_seg=True)

    def set_hoc_pt3d_geometry(self, **kwargs):
        ''' given reduced cell, set pt3d points in neuron.  default is that the 0 end of the soma is at (0,0,0).  basal and additional soma segments extend in the negative y direction.  apical dendrites extend in the positive direction.  by default the soma is assumed to be a single segment.
        '''
        # list of potential subtrees to update 3d info
        subtrees = ['soma', 'dend','apic']
        # get length of soma
        soma_L = self.reduced_cell.soma.L
        # create proximal apical branch by copying basal dendritic branch
        #----------------------------------------------------------------
        # print 'getattr(self.reduced_cell, subtree_key)', getattr(self.reduced_cell, 'dend').remove
        move_basal_section_to_apical =False
        if 'move_basal_section_to_apical' in kwargs:
            move_basal_section_to_apical = kwargs['move_basal_section_to_apical']
        # if move_basal_section_to_apical:
        #     assert len(subtrees['dend'])>1, 'two convert a basal section to apical, there must be at least two basal sections'
        #     subtrees['apic'].append(subtrees['dend'].pop())
        # create pt3d geometry
        #-----------------------
        # iterate subtrees
        for subtree_key in subtrees:
            # get list of subtree sections
            subtree = getattr(self.reduced_cell, subtree_key)
            # if the subtree is only a single section, convert to iterable
            if hasattr(subtree, 'nseg'):
                subtree=[subtree]
            # iterate over sections

            for section_i, section in enumerate(subtree):
                # nseg and length
                nseg = section.nseg
                L = section.L
                # clear 3d points
                h.pt3dclear(sec=section)
                if subtree_key=='soma':
                    # add 3d points
                    h.pt3dadd(0.0, 0.0, 0.0, section.diam, sec=section)
                    h.pt3dadd(0.0, -1.*L, 0.0, section.diam, sec=section)
                else:
                    # iterate over segments
                    for segment in section:
                        # x location of current segment
                        seg_x = segment.x
                        # distance of current segment from zero end of the section
                        seg_dist = seg_x*L
                        # xyz values
                        y_val = seg_dist
                        x_val = 0.
                        z_val = 0.
                        # set basal dendrites to negative y values
                        if subtree_key=='dend':
                            # if not moving to apical, y values are negative and added to the length of the soma
                            if not move_basal_section_to_apical or section_i<len(subtree)-1:
                                # adjust for length of the soma
                                y_val += copy.copy(soma_L)
                                y_val *= -1.
                        # soma extends in the negative y direction
                        elif subtree_key =='soma':
                            y_val*=-1
                        # add 3d points
                        h.pt3dadd(x_val, y_val, z_val, segment.diam, sec=section)

    def create_geo(self, **kwargs):
        ''' FIXME comments
        '''
        # list of potential subtrees to update 3d info
        subtrees = ['soma', 'dend','apic']
        geo={}
        for subtree_key in subtrees:
            geo[subtree_key]=[]
            # get list of subtree sections
            subtree = getattr(self.reduced_cell, subtree_key)
            # if the subtree is only a single section, convert to iterable
            if hasattr(subtree, 'nseg'):
                subtree=[subtree]
            # iterate over sections
            for section in subtree:
                geo[subtree_key].append(section)
        return geo

    def _create_loc_map(self, original_cell, **kwargs):
        ''' for each location in this reduced cell return a list of corresponding locations in the original cell
        '''
        loc_map={}
        for key, val in self.seg_to_seg_text.iteritems():
            orig_loc = original_cell._get_loc(geo=original_cell.geo, hoc_name=key[0], seg_x=key[1])
            new_loc = self._get_loc(geo=self.geo, hoc_name=val[0], seg_x=key[1])
            if new_loc not in loc_map:
                loc_map[new_loc]=[]    
            loc_map[new_loc].append(orig_loc)
        return loc_map

    def _get_polarization_map(self, vtrace_df_original, loc_map, field):
        '''
        '''
        v_df = functions._set_index(vtrace_df_original,['field', 'location'])
        polarization_map={}
        for key, val in loc_map.iteritems():
            if key not in polarization_map:
                polarization_map[key]=[]
            for loc in val:
                polarization = v_df.loc[[(field, loc)]].polarization.values[0]
                if np.isnan(polarization):
                    polarization=0.
                polarization_map[key].append(polarization)
        # add soma
        soma_locs = []
        for tree_key, tree in self.geo.iteritems():
            for sec_i, sec in enumerate(tree):
                for seg_i, seg in enumerate(sec):
                    if tree_key=='soma':
                        location = (tree_key, sec_i, seg_i)
                        soma_locs.append(location)
        for loc in soma_locs:
            if loc not in polarization_map:
                polarization_map[loc]=[]
            polarization = v_df.loc[[(field, loc)]].polarization.values[0]
            if np.isnan(polarization):
                polarization=0.
            polarization_map[loc].append(polarization)
        return polarization_map

class CellBranco2010(Cell):
    """ pyramidal neuron based on Migliore et al. 2005

    An instance of this object will creates a cell (hoc objects) at the top level of the hoc interpreter using the hoc files in _init_geometry.  The .geo attribute contains a python mapping to these hoc objects.  The geo object is organized as geo['section tree'][section](segment location)

    the syns attribute creates a container for synapse objects that are added to each segment in the hoc cell.  syns is organized as syns['section tree']['synapse type'][section][segment number]
    
    """
    def __init__(self, p):
        '''
        '''
        self.p=p
        super(CellBranco2010, self).__init__(self.p)
    # def __init__(self,p):
    #     '''
    #     '''
    #     # set geometry and initialize synapses
    #     self.geo, self.syns = self.geometry(p)
    #     # create morphology object for displaying neuron
    #     self.morpho = self._create_morpho(self.geo)
    #     # create
    #     self.mechanisms(p)
        

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
        h.load_file('rc19.hoc')  
        # h.load_file('geoc62564.hoc')
        # set discretization based on dlambda rule (set dlambda in hoc file) 
        h.load_file('fixnseg.hoc')      
        # dictionary for storing geometry ['tree'][sec](seg location)
        self.geo = {}
        # dictionary for storing synapse objects ['tree']['type'][sec][seg]
        self.syns = {}
        # add section trees to geometry dictionary
        self.geo['soma'] = [h.soma]
        self.geo['axon'] =  [h.axon]
        self.geo['basal'] = h.dend
        self.geo['apical'] = h.apic
        
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
                elif (tree_key == 'basal') or (tree_key == 'apical'):
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

class CellGrienberger2017(Cell):
    '''
    '''
    def __init__(self, ):
        '''
        '''
        self.morph_filename='EB2-late-bifurcation.swc'
        self.mech_filename = '043016 Type A - km2_NMDA_KIN5_Pr'
        # self.cell = specify_cells.CA1_Pyr(morph_filename=morph_filename, mech_filename=mech_filename, full_spines=True)

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
        p['nseg']=1
        
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

                # self.geo[tree][sec_i].nseg=p['nseg']

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
                    # sec.insert('calH')
                    # sec.gcalbar_calH = p['gcalbar']

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

class DendriteTransform:
    """"""
    def __init__(self, p):
        cell1 = CellMigliore2005(p)
        apical_transform = self.dendrite_transform(geo=cell1.geo, python_tree=['apical_trunk','apical_tuft'], neuron_tree=['user5', 'apical_dendrite'])
        basal_transform = self.dendrite_transform(geo=cell1.geo, python_tree=['basal'], neuron_tree=['dendrite'])
        print 'apical:', apical_transform['a_cable'], apical_transform['L_cable']
        print 'basal:', basal_transform['a_cable'], basal_transform['L_cable']

    def measure_area(self, tree):
        """
        given a tree measure the total area of all sections in the tree

        tree is a list of sections (hoc objects)
        """
        area_all = []
        for sec_i, sec in enumerate(tree):
            # convert to um to cm (*.0001)
            L = .0001*sec.L
            a = .0001*sec.diam/2.
            rL = sec.Ra
            rm = 1/sec.g_pas
            area = 2*np.pi*a*L
            lam = np.sqrt(a*rm/(2*rL))
            area_all.append(area)

        return sum(area_all)

    def measure_length(self, geo):
        """ measure electrotonic length for each path along a cells dendritic tree
        """
        # keep track of most recent section in each path [paths]
        secs = [geo['soma'][0]]

        # keep track of all sections in paths [paths][sections]
        # does not include soma
        paths=[[]]
        
        # iterate over paths (most recent section)
        for sec_i, sec in enumerate(secs):
            
            # current section
            current_sec_ref = h.SectionRef(sec=sec)
            # current children
            current_children = current_sec_ref.child
            
            # if there are children
            while len(current_children)>0:
                
                # iterate over current children 
                for child_i, child in enumerate(current_children):

                    # add first child to current path
                    if child_i==0:
                        paths[sec_i].append(child)
                        # update current section
                        sec = child
                    
                    # if multiple children
                    if child_i > 0:
                        # create new path to copy previous path when tree splits
                        new_path=[]
                        for section in paths[sec_i]:
                            new_path.append(section)

                        # copy current path in list
                        # if split occurs at soma, do not copy previous list
                        if h.SectionRef(sec=child).parent.name()!='soma':
                            paths.append(new_path)
                        else:
                            # create new list, not including soma
                            paths.append([])
                        
                        # add corresponding child to the current section list
                        secs.append(child)
                        
                        # add corresponding child to new path 
                        paths[sec_i + child_i].append(child)
                
                # update current section and children       
                current_sec_ref = h.SectionRef(sec=sec)
                current_children = current_sec_ref.child


        # calculate electrotonic length of each path
        path_l = [] # [paths][electrotonic section length]
        sec_name = [] # [paths][section names]
        for path_i, path in enumerate(paths):
            path_l.append([])
            sec_name.append([])
            for sec_i, sec in enumerate(path):
                # convert all distances in cm
                # section length
                L = .0001*sec.L
                # section radius
                a = .0001*sec.diam/2
                # membrane resistivity
                rm = 1/sec.g_pas
                # axial resistivity
                rL = sec.Ra
                # space constant lambda
                lam = np.sqrt( (a*rm) / (2*rL) )
                # electrotonic length
                e_length = L/lam
                # electrotonic lengths for all paths and sections [paths][sections]
                path_l[path_i].append(e_length) 
                # keep track of section names [paths][sections]
                sec_name[path_i].append(sec.name())
        # print path_l[0]
        return {'path_l': path_l,
        'sec_name':sec_name}

    def dendrite_transform(self, geo, python_tree, neuron_tree):
        """ equivalent cable transform for dendritic tree
        """
        # FIXME
        rL_cable = 150.
        rm_cable = 28000.
        # area
        A_full=0
        for tree_i, tree in enumerate(python_tree):
            A_full += self.measure_area(geo[tree])

        paths = self.measure_length(geo)
        e_lengths = []
        # iterate over all paths
        for path_i, path in enumerate(paths['path_l']):
            # only keep paths in the neuron_tree argument
            for tree in neuron_tree:
                for sec in paths['sec_name'][path_i]:
                    if tree in sec:
                        e_lengths.append(sum(path))
                        break

                        

        EL_full = np.mean(e_lengths) 

        # convert cm back to um (*10000)
        L_cable = (EL_full**(2./3)) * (A_full*rm_cable/(4.*np.pi*rL_cable))**(1./3)
        a_cable = A_full/(2.*np.pi*L_cable)

        return {'a_cable':10000*a_cable, 'L_cable':10000*L_cable}


class NetworkEIF(object):
    '''
    '''
    def __init__(self, ):
        pass

    def run(self, ):
        '''
        '''
        N_E = 4
        N_I= 2
        
    def syn_double_exp(self, dt, W, g_rise, g_decay,  pre_spikes, t_rise, t_decay, **kwargs):
        '''
        '''
        tp = (t_rise*t_decay)/(t_decay-t_rise)*np.log(t_decay/t_rise)
        factor = -np.exp(-tp/t_rise) + np.exp(-tp/t_decay)
        factor = 1/factor

        g_rise = g_rise -dt*rise/t_rise + factor*W.dot(pre_spikes)
        g_decay = g_decay -dt*decay/t_decay + factor*W.dot(pre_spikes)

        return g_rise, g_decay

    def syn_single_exp(self, dt, W, g, t_decay, pre_spikes, **kwargs):
        '''
        '''
        g = g -dt*g/t_decay + W.dot(pre_spikes)

        return g

    def adex_update(self, i, u, wad, z, V_T, I_syn, I_ext, p, dt=0.1, **kwargs):
        '''
        ==Args==
        :dt: scalar: time step in ms
        :i: integer: index of current time step
        :p: dictionary: parameters (see below for details)

        state variables
        :u: N x T array: membrane voltage
        :wad: N x T array: adaptation current
        :z: N x T array: spike after-potential
        :V_T:N x T array: adaptive threshold
        :spikes: N x T boolean array: binary array of spikes

        input currents
        :I_syn: N x T array: summed synaptic currents
        :I_ext: N x T array: summed external applied current
        
        parameters
        :C: membrane capacitance
        :g_L:  leak conductance
        :E_L: reversal potential
        :V_Tr: adaptive threshold baseline
        :V_spike: spike detection value
        :delta_T:  exponential slope
        :t_wad:  adaptation current time constant
        :a_sub: subthreshold adaptation conductance
        :b_spike: spike triggered adaptation current
        :I_sp: spike after-current
        :t_z: spike after-potential time constant
        :t_VT: threshold adaptation time constant
        :V_Tmax: threshold set point after spike
        :V_reset: membrane voltage reset after spike

        ==Returns==
        :u
        :wad
        :z
        :V_T
        :spikes
        '''
        # voltage update
        #--------------------------------------------------------------------
        du = (1/p['C'])*(
            -p['g_L']*(u[:,i-1]-p['E_L'])  +  
            -p['g_L']*p['delta_T']*np.exp((u[:,i-1]-V_T[:,i-1])/p['delta_T']) - 
            wad[:,i-1] + 
            z[:,i-1] +
            I_ext[:, i-1] + 
            I_syn[:,i-1]) 
        u[:,i] = u[:,i-1] + dt*du

        # reset voltage for cells that cross threshold on previous step
        u[spikes[:,i-1],i] = p['V_reset']

        # get cells that cross on the current step
        spikes[:,i] = u[:,i]>=p['V_spike']

        # set voltage to max if threshold is crossed
        u[spikes[:,i], i] = p['V_spike']

        # adaptation current update
        #-------------------------------------------------------------------
        dwad = (1/p['t_wad'])*( p['a_sub']*( u[:,i-1]-E_L) - wad[:,i-1])
        wad[:,i] = wad[:,i-1] + dt*dwad
        wad[spikes[:,i-1],i] = wad[spikes[:,i-1],i-1] + p['b_spike']

        # spike after-potential
        #-------------------------------------------------------------------
        dz = -(1/p['t_z'])*z[:,i-1]
        z[:,i] = z[:,i-1] + dt*dz
        z[spikes[:,i-1],i] = p['I_sp']

        # threshold
        #------------------------------------------------------------------
        dV_T = (1/p['t_V_T'])*(V_T[:,i-1]-p['V_Tr'])
        V_T[:,i] = V_T[:,i-1] + dt*dV_T
        V_T[spikes[:,i-1],i] = p['V_Tmax']


        return u, spikes, wad, z, V_T


