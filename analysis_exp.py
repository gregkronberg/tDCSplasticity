"""
analysis

"""
import numpy as np
import analysis
import glob
import pickle
import copy
import matplotlib.pyplot as plt
import itertools
from matplotlib import cm as colormap
from scipy import stats
import os

# specify clopath parameters to be shared across experiments
#------------------------------------------------------------
clopath_param_all= {
        'A_m0':10*100E-5, # depression magnitude parameter (mV^-1)
        'tetam':-70,#-41, # depression threshold (mV)
        'tetap':-67,#-38, # potentiation threshold (mV)
        'tau_x':8,#-38, # potentiation threshold (mV)
        'tau_m':20,#-38, # potentiation threshold (mV)
        'tau_p': 3, # time constant (ms) for low pass filter post membrane potential for potentiation
        'A_p':10*40E-5, # amplitude for potentiation (mV^-2)
        'delay':0, # conduction delay (ms)
        'LTD_delay':1,
        }
        # stdp parameters
        #---------------
        # clopath_param= {
        # 'A_m0':10*100E-5, # depression magnitude parameter (mV^-1)
        # 'tetam':-71,#-41, # depression threshold (mV)
        # 'tetap':-50,#-38, # potentiation threshold (mV)
        # 'tau_x':10,#-38, # potentiation threshold (mV)
        # 'tau_m':20,#-38, # potentiation threshold (mV)
        # 'tau_p': 30, # time constant (ms) for low pass filter post membrane potential for potentiation
        # 'A_p': 40*40E-5, # amplitude for potentiation (mV^-2)
        # 'delay':0, # conduction delay (ms)
        # 'LTD_delay':1,
        # }

# dpi value for saving png images
dpi = 350
class Experiment:
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        if not kwargs:
            pass
        else:
            experiment = getattr(self, kwargs['experiment'])

            experiment(**kwargs) 

    def exp_1a1(self, **kwargs):
    	'''
    	'''
        # FIXME data gets too large for pickle to handle. this is fixed in python 3. try breaking into smaller objects or using pandas
    	# get list of data files
        directory = 'Data/'+kwargs['experiment']+'/'
        search_string = '*data*'

        group_file_name = 'group.pkl'

        Group = analysis.GroupData(directory=directory, file_name=group_file_name, **kwargs)
        Plots = analysis.Plots()
        variables = ['v','input_times']

        # updated clopath parameters
        # clopath_param= {
        # 'A_m0':100E-5, # depression magnitude parameter (mV^-1)
        # 'tetam':-70,#-41, # depression threshold (mV)
        # 'tetap':-67,#-38, # potentiation threshold (mV)
        # 'tau_x':8,#-38, # potentiation threshold (mV)
        # 'tau_m':20,#-38, # potentiation threshold (mV)
        # 'tau_p': 3, # time constant (ms) for low pass filter post membrane potential for potentiation
        # 'A_p': 40E-5, # amplitude for potentiation (mV^-2)
        # 'delay':0, # conduction delay (ms)
        # 'LTD_delay':1,
        # }
        clopath_param= {
        'A_m0':10*100E-5, # depression magnitude parameter (mV^-1)
        'tetam':-71,#-41, # depression threshold (mV)
        'tetap':-50,#-38, # potentiation threshold (mV)
        'tau_x':10,#-38, # potentiation threshold (mV)
        'tau_m':20,#-38, # potentiation threshold (mV)
        'tau_p': 30, # time constant (ms) for low pass filter post membrane potential for potentiation
        'A_p': 40*40E-5, # amplitude for potentiation (mV^-2)
        'delay':0, # conduction delay (ms)
        'LTD_delay':1,
        }

        if 'update_group_data' in kwargs and kwargs['update_group_data']:
            self.group_data = Group._standard_run(group_data=Group.group_data, directory=directory, file_name=group_file_name, search_string=search_string, clopath_param=clopath_param, variables=variables)
            Group._save_group_data(group_data=self.group_data, directory=directory, file_name=group_file_name, )

        condition_types=['polarity', 'syn_num', 'dist','path', 'tree']
        polarities=['anodal','cathodal','control']
        paths=['1']
        syn_nums = [4,5,6,8,10,12,16,20]
        dists = [[0,300],[100,400],[200,500],[300,600],[0,200],[100,300],[200,400],[300,500],[400,600],]
        trees = ['apical_tuft','apical_trunk','basal']
        condition_combos = [temp for temp in itertools.product(syn_nums, dists)]
        figs_w=[]
        figs_xcorr=[]
        for combo_i, combo in enumerate(condition_combos):
            conditions_current = [polarities, [[combo[0]]], [[combo[1]]],[paths],[trees]]
            conditions = zip(condition_types, conditions_current)

            figs_w.append(
                Plots._range_variable_group(
                    group_data=self.group_data,
                    x_variable='t',
                    y_variable='clopath',
                    conditions=conditions
                    ))
            plt.xlabel('time (ms)')
            plt.ylabel('weight')
            plt.title(str(combo[0])+'synapses, '+str(combo[1])+'from soma')
            # plt.show(block=False)
            figs_xcorr.append(
                Plots._xcorr_mean(group_data=self.group_data, conditions=conditions ))
            plt.xlabel('delay (ms), negative = dendrite first')
            plt.ylabel('probability')
            plt.title(str(combo[0])+'synapses, '+str(combo[1])+'from soma')
            # plt.show(block=False)
            fig_w_name = directory+'weights_'+str(combo[0])+'syn_'+str(combo[1])+'dist.png'
            fig_xcorr_name = directory+'xcorr_'+str(combo[0])+'syn_'+str(combo[1])+'dist.png'
            # plt.show(block=False)
            figs_w[-1].savefig(fig_w_name, dpi=250)
            figs_xcorr[-1].savefig(fig_xcorr_name, dpi=250)

    def exp_1a1_nablock(self, **kwargs):
        '''
        '''
        # get list of data files
        directory = 'Data/'+kwargs['experiment']+'/'
        search_string = '*data*'

        group_file_name = 'group.pkl'

        Group = analysis.GroupData(directory=directory, file_name=group_file_name, **kwargs)
        Plots = analysis.Plots()
        variables = ['v','input_times']

        # updated clopath parameters
        clopath_param= {
        'A_m0':100E-5, # depression magnitude parameter (mV^-1)
        'tetam':-70,#-41, # depression threshold (mV)
        'tetap':-67,#-38, # potentiation threshold (mV)
        'tau_x':8,#-38, # potentiation threshold (mV)
        'tau_m':20,#-38, # potentiation threshold (mV)
        'tau_p': 3, # time constant (ms) for low pass filter post membrane potential for potentiation
        'A_p': 40E-5, # amplitude for potentiation (mV^-2)
        'delay':0, # conduction delay (ms)
        'LTD_delay':1,
        }

        self.group_data = Group._standard_run(group_data=Group.group_data, directory=directory, file_name=group_file_name, search_string=search_string, clopath_param=clopath_param, variables=variables)
        Group._save_group_data(group_data=self.group_data, directory=directory, file_name=group_file_name, )

        condition_types=['polarity', 'syn_num', 'dist','path', 'tree']
        polarities=['anodal','cathodal','control']
        paths=['1']
        syn_nums = [4,5,6,8,10,12,16,20]
        dists = [[0,300],[100,400],[200,500],[300,600],[0,200],[100,300],[200,400],[300,500],[400,600],]
        trees = ['apical_tuft','apical_trunk','basal']
        condition_combos = [temp for temp in itertools.product(syn_nums, dists)]
        figs_w=[]
        figs_xcorr=[]
        for combo_i, combo in enumerate(condition_combos):
            conditions_current = [polarities, [[combo[0]]], [[combo[1]]],[paths],[trees]]
            conditions = zip(condition_types, conditions_current)

            figs_w.append(
                Plots._range_variable_group(
                    group_data=self.group_data,
                    x_variable='t',
                    y_variable='clopath',
                    conditions=conditions
                    ))
            plt.xlabel('time (ms)')
            plt.ylabel('weight')
            plt.title(str(combo[0])+'synapses, '+str(combo[1])+'from soma')
            figs_xcorr.append(
                Plots._xcorr_mean(group_data=self.group_data, conditions=conditions ))
            plt.xlabel('delay (ms), negative = dendrite first')
            plt.ylabel('probability')
            plt.title(str(combo[0])+'synapses, '+str(combo[1])+'from soma')
            fig_w_name = directory+'weights_'+str(combo[0])+'syn_'+str(combo[1])+'dist.png'
            fig_xcorr_name = directory+'xcorr_'+str(combo[0])+'syn_'+str(combo[1])+'dist.png'
            # plt.show(block=False)
            figs_w[-1].savefig(fig_w_name, dpi=250)
            figs_xcorr[-1].savefig(fig_xcorr_name, dpi=250)

    def exp_1a_dose(self, **kwargs):
        '''
        '''
        # directories and filenames
        #---------------------------------------------------------------
        directory = 'Data/'+kwargs['experiment']+'/'
        figure_directory = 'Data/'+kwargs['experiment']+'/Figures/'
        # make figure directory if doesnt exist
        if os.path.isdir(figure_directory) is False:
            os.mkdir(figure_directory)
        search_string = '*data*'
        filename_vtrace = 'vtrace_df.pkl'
        filename_w_clopath = 'w_clopath_df.pkl'
        filename_spikes = 'spikes_df.pkl'

        funcs = analysis.VarFuncs()

        #####################################################################
        # weights
        #####################################################################
        run_weights=False
        if run_weights:
            # load group data
            #-------------------------------------------------------------
            self.w_clopath_df = analysis._load_group_data(directory=directory, filename=filename_w_clopath)

            # update group data w_clopath
            #----------------------------------------------------------------
            # updated clopath parameters
            clopath_param= {
            'A_m0':10*100E-5, # depression magnitude parameter (mV^-1)
            'tetam':-70,#-41, # depression threshold (mV)
            'tetap':-67,#-38, # potentiation threshold (mV)
            'tau_x':8,#-38, # potentiation threshold (mV)
            'tau_m':20,#-38, # potentiation threshold (mV)
            'tau_p': 3, # time constant (ms) for low pass filter post membrane potential for potentiation
            'A_p':10*40E-5, # amplitude for potentiation (mV^-2)
            'delay':0, # conduction delay (ms)
            'LTD_delay':1,
            }
            
            kwarg = {'clopath_param':clopath_param}
            functions = [funcs._get_w_clopath]
            kwlist = [kwarg]
            rerun=[]
            keep=[]
            w_clopath_df_temp = copy.deepcopy(self.w_clopath_df)
            self.w_clopath_df = analysis._process_new_data_df(group_df=self.w_clopath_df, preprocessed_directory=directory, functions=functions, kwlist=kwlist, rerun=rerun, keep=keep, file_limit=[])

            # save group data
            #----------------
            if not w_clopath_df_temp.equals(self.w_clopath_df):
                print 'saving updated group data'
                self.w_clopath_df.to_pickle(directory+filename_w_clopath)

            # dose response figures
            #-------------------------------------------------------------------
            figtype='dose_response_'
            figdf = analysis.BuildFigDF()._dose_response()
            figs, ax = analysis.PlotFuncs()._dose_response(df=self.w_clopath_df, figdf=figdf, variable='dw_clopath')
            # save figure
            #------------
            for fig_key, fig in figs.iteritems():
                fname = figure_directory+figtype+str(fig_key)+'.png'
                fig.savefig(fname, format='png', dpi=dpi)

        #####################################################################
        # vtrace
        #####################################################################
        run_vtrace=False
        if run_vtrace:
            # load group data
            #-------------------------------------------------------------
            self.vtrace_df = analysis._load_group_data(directory=directory, filename=filename_vtrace)

            # update group data vtrace
            #----------------------------------------------------------------
            functions = [funcs._get_vtrace]
            kwlist = [{}]
            rerun=[]
            keep=[]
            vtrace_df_temp = copy.deepcopy(vtrace_df)
            self.vtrace_df = analysis._process_new_data_df(group_df=self.vtrace_df, preprocessed_directory=directory, functions=functions, kwlist=kwlist, rerun=rerun, keep=keep, file_limit=[])
            # save group data
            #----------------
            if not vtrace_df_temp.equals(self.vtrace_df):
                print 'saving updated group data'
                self.vtrace_df.to_pickle(directory+filename_vtrace)

            # plot average voltage trace at soma and dendrite
            #------------------------------------------------


        #####################################################################
        # spikes
        #####################################################################
        run_spikes=False
        if run_spikes:
            # load group data
            #-------------------------------------------------------------
            self.spikes_df = analysis._load_group_data(directory=directory, filename=filename_spikes)

            # update group data vtrace
            #----------------------------------------------------------------
            functions = [funcs._get_spikes]
            kwlist = [{'threshold':-30}]
            rerun=[]
            keep=[]
            spikes_df_temp = copy.deepcopy(self.spikes_df)
            self.spikes_df = analysis._process_new_data_df(group_df=self.spikes_df, preprocessed_directory=directory, functions=functions, kwlist=kwlist, rerun=rerun, keep=keep, file_limit=[])
            # save group data
            #----------------
            if not spikes_temp.equals(self.spikes_df):
                print 'saving updated group data'
                self.spikes_df.to_pickle(directory+filename_spikes)


    def exp_4a1(self, **kwargs):
        # get list of data files
        directory = 'Data/'+kwargs['experiment']+'/'
        search_string = '*data*'

        group_file_name = 'group.pkl'

        self.Group = analysis.GroupData(directory=directory, file_name=group_file_name, **kwargs)
        self.Plots = analysis.Plots()
        variables = ['v','input_times']

        # updated clopath parameters
        clopath_param= {
        'A_m0':100E-5, # depression magnitude parameter (mV^-1)
        'tetam':-70,#-41, # depression threshold (mV)
        'tetap':-67,#-38, # potentiation threshold (mV)
        'tau_x':8,#-38, # potentiation threshold (mV)
        'tau_m':20,#-38, # potentiation threshold (mV)
        'tau_p': 3, # time constant (ms) for low pass filter post membrane potential for potentiation
        'A_p': 40E-5, # amplitude for potentiation (mV^-2)
        'delay':0, # conduction delay (ms)
        'LTD_delay':1,
        }

        if 'update_group_data' in kwargs and kwargs['update_group_data']:
            self.group_data = self.Group._standard_run(group_data=self.Group.group_data, directory=directory, file_name=group_file_name, search_string=search_string, clopath_param=clopath_param, variables=variables)
            self.Group._save_group_data(group_data=self.group_data, directory=directory, file_name=group_file_name, )

        condition_types=['polarity', 'syn_num', 'dist','path']
        polarities=['anodal','cathodal','control']
        paths=[['1'],['2']]
        syn_nums = [4,5,6,8,10,12,16,20]
        dists = [[0,200],[0,300],]
        condition_combos = [temp for temp in itertools.product(syn_nums, dists, paths)]
        figs=[]
        for combo_i, combo in enumerate(condition_combos):
            conditions_current = [polarities, [[combo[0]]], [[combo[1]]],[[combo[2]]]]
            conditions = zip(condition_types, conditions_current)

            figs.append(
                self.Plots._range_variable_group(
                    group_data=self.group_data,
                    x_variable='t',
                    y_variable='clopath',
                    conditions=conditions
                    ))
            fig_name = directory+'weights_'+str(combo[0])+'syn_'+str(combo[1])+'dist_path'+str(combo[2])+'.png'
            # plt.show(block=False)
            figs[-1].savefig(fig_name, dpi=250)
            plt.xlabel('time (ms)')
            plt.ylabel('weight')
            plt.title('path'+str(combo[2])+','+str(combo[0])+'synapses, '+str(combo[1])+'from soma')
            
    def exp_4a2(self, **kwargs):
        # get list of data files
        directory = 'Data/'+kwargs['experiment']+'/'
        search_string = '*data*'

        group_file_name = 'group.pkl'

        self.Group = analysis.GroupData(directory=directory, file_name=group_file_name, **kwargs)
        self.Plots = analysis.Plots()
        variables = ['v','input_times']

        # updated clopath parameters
        clopath_param= {
        'A_m0':100E-5, # depression magnitude parameter (mV^-1)
        'tetam':-70,#-41, # depression threshold (mV)
        'tetap':-67,#-38, # potentiation threshold (mV)
        'tau_x':8,#-38, # potentiation threshold (mV)
        'tau_m':20,#-38, # potentiation threshold (mV)
        'tau_p': 3, # time constant (ms) for low pass filter post membrane potential for potentiation
        'A_p': 40E-5, # amplitude for potentiation (mV^-2)
        'delay':0, # conduction delay (ms)
        'LTD_delay':1,
        }

        self.group_data = self.Group._standard_run(group_data=self.Group.group_data, directory=directory, file_name=group_file_name, search_string=search_string, clopath_param=clopath_param, variables=variables)
        self.Group._save_group_data(group_data=self.group_data, directory=directory, file_name=group_file_name, )

        condition_types=['polarity', 'syn_num', 'dist','path']
        polarities=['anodal','cathodal','control']
        paths=[['2']]
        syn_nums = [4,5,6,8,10,12,16,20]
        dists = [[0,200],[0,300],]
        condition_combos = [temp for temp in itertools.product(syn_nums, dists)]
        figs=[]
        for combo_i, combo in enumerate(condition_combos):
            conditions_current = [polarities, [[combo[0]]], [[combo[1]]],[paths[0]]]
            conditions = zip(condition_types, conditions_current)

            figs.append(
                self.Plots._range_variable_group(
                    group_data=self.group_data,
                    x_variable='t',
                    y_variable='clopath',
                    conditions=conditions
                    ))
            fig_name = directory+'weights_'+str(combo[0])+'syn_'+str(combo[1])+'dist.png'
            # plt.show(block=False)
            figs[-1].savefig(fig_name, dpi=250)
            plt.xlabel('time (ms)')
            plt.ylabel('weight')
            plt.title(str(combo[0])+'synapses, '+str(combo[1])+'from soma')

    def exp_5a1(self, **kwargs):
        # get list of data files
        directory = 'Data/'+kwargs['experiment']+'/'
        search_string = '*data*'

        # # updated clopath parameters
        # clopath_param= {
        # 'A_m0':10*100E-5, # depression magnitude parameter (mV^-1)
        # 'tetam':-70,#-41, # depression threshold (mV)
        # 'tetap':-50,#-38, # potentiation threshold (mV)
        # 'tau_x':20,#-38, # potentiation threshold (mV)
        # 'tau_m':20,#-38, # potentiation threshold (mV)
        # 'tau_p': 30, # time constant (ms) for low pass filter post membrane potential for potentiation
        # 'A_p': 100*40E-5, # amplitude for potentiation (mV^-2)
        # 'delay':0, # conduction delay (ms)
        # 'LTD_delay':1,
        # }
        # # updated clopath parameters
        clopath_param= {
        'A_m0':10*100E-5, # depression magnitude parameter (mV^-1)
        'tetam':-71,#-41, # depression threshold (mV)
        'tetap':-50,#-38, # potentiation threshold (mV)
        'tau_x':10,#-38, # potentiation threshold (mV)
        'tau_m':20,#-38, # potentiation threshold (mV)
        'tau_p': 30, # time constant (ms) for low pass filter post membrane potential for potentiation
        'A_p': 40*40E-5, # amplitude for potentiation (mV^-2)
        'delay':0, # conduction delay (ms)
        'LTD_delay':1,
        }

        stdp_dts = [-10, 10]
        markers = ['.','x']
        pulse_freqs = [1,5,10,20,30,40,50,75, 100]
        self.Group=[]
        self.group_data= []
        self.dw=[]
        for freq_i, freq in enumerate(pulse_freqs):

            group_file_name = 'group_freq_'+str(freq)+'.pkl'

            self.Group.append(analysis.GroupData(directory=directory, file_name=group_file_name, **kwargs))
            self.Plots = analysis.Plots()
            variables = ['v','input_times']

            self.group_data.append(self.Group[freq_i]._conditional_run(group_data=self.Group[freq_i].group_data, directory=directory, file_name=group_file_name, search_string=search_string, clopath_param=clopath_param, variables=variables, parameter='pulse_freq', parameter_value=freq, path_parameter=True))
            self.group_data[freq_i] = self.Group[freq_i]._add_parameter_to_group_data(group_data=self.group_data[freq_i], parameter='stdp_dt')
            self.group_data[freq_i] = self.Group[freq_i]._add_parameter_to_group_data(group_data=self.group_data[freq_i], parameter='pulse_freq',path=True)
            self.Group[freq_i]._save_group_data(group_data=self.group_data[freq_i], directory=directory, file_name=group_file_name, )

        plt.figure()
        self.dw_mean=[]
        for stdp_dt_i, stdp_dt in enumerate(stdp_dts):
            self.dw.append([])
            self.dw_mean.append([])
            for freq_i, freq in enumerate(pulse_freqs):
                idx = [temp_i for temp_i, temp in enumerate(self.group_data[freq_i]['clopath']['stdp_dt']) if temp==stdp_dt]
                # get weight change
                dw_current = self.group_data[freq_i]['clopath']['data'][idx,-1]/self.group_data[freq_i]['clopath']['data'][idx,0]
                self.dw[stdp_dt_i].append(dw_current)
                self.dw_mean[stdp_dt_i].append(np.mean(dw_current))

                marker = markers[stdp_dt_i]
                plt.plot(freq, np.mean(dw_current), marker=marker, color='black', markersize=20, markeredgewidth=2)
                plt.errorbar(freq, np.mean(dw_current),yerr=stats.sem(dw_current), color='black', markersize=20)
            plt.plot(pulse_freqs, self.dw_mean[stdp_dt_i], color='black', linewidth=4,)
        plt.show(block=False)
        plt.axhline(y=1, linewidth=2, linestyle='--', color='black')
        plt.xlabel('Frequency (Hz)', fontsize=20, fontweight='heavy')
        plt.ylabel('Normalized weight', fontsize=20, fontweight='heavy')
        plt.title('Rate dependence of STDP', fontsize=20, fontweight='heavy')
        # self.Plots._stdp_frequency_weight(group_data=self.group_data)
        # condition_types=['stdp_dt']
        # polarities=['control']
        # paths=[['1']]
        # syn_nums = [5]
        # dists = [[0,200],]
        # stdp_dts = [[-10], [10]]
        # condition_combos = [temp for temp in itertools.product(syn_nums, dists)]
        # figs=[]
        # for syn_num_i, syn_num in enumerate(syn_nums):
        #     conditions_current = [polarities, [[syn_num]], [[dists[0]]],[paths[0]], stdp_dts]
        #     conditions = zip(condition_types, conditions_current)

        #     figs.append(
        #         self.Plots._range_variable_group(
        #             group_data=self.group_data,
        #             x_variable='t',
        #             y_variable='clopath',
        #             conditions=conditions
        #             ))
        #     fig_name = directory+'weights_'+str(syn_num)+'syn_'+str(stdp_dts)+'stdp_dt.png'
        #     # plt.show(block=False)
        #     figs[-1].savefig(fig_name, dpi=250)
        #     plt.xlabel('time (ms)')
        #     plt.ylabel('weight')
        #     plt.title(str(syn_num)+'synapses, '+str(stdp_dts)+'timing')

    def exp_5a2(self, **kwargs):
        # get list of data files
        directory = 'Data/'+kwargs['experiment']+'/'
        search_string = '*data*'

        group_file_name = 'group.pkl'

        self.Group = analysis.GroupData(directory=directory, file_name=group_file_name, **kwargs)
        self.Plots = analysis.Plots()
        variables = ['v','input_times']

        # updated clopath parameters
        clopath_param= {
        'A_m0':10*100E-5, # depression magnitude parameter (mV^-1)
        'tetam':-71,#-41, # depression threshold (mV)
        'tetap':-50,#-38, # potentiation threshold (mV)
        'tau_x':10,#-38, # potentiation threshold (mV)
        'tau_m':20,#-38, # potentiation threshold (mV)
        'tau_p': 30, # time constant (ms) for low pass filter post membrane potential for potentiation
        'A_p': 40*40E-5, # amplitude for potentiation (mV^-2)
        'delay':0, # conduction delay (ms)
        'LTD_delay':1,
        }

        self.group_data = self.Group._standard_run(group_data=self.Group.group_data, directory=directory, file_name=group_file_name, search_string=search_string, clopath_param=clopath_param, variables=variables)
        self.group_data = self.Group._add_parameter_to_group_data(group_data=self.group_data, parameter='stdp_dt')
        self.Group._save_group_data(group_data=self.group_data, directory=directory, file_name=group_file_name, )

        stdp_dts = list(set(self.group_data['clopath']['stdp_dt']))
        input_delay=3
        self.dw=[]
        plt.figure()
        for stdp_dt_i, stdp_dt in enumerate(stdp_dts):
            idx = [temp_i for temp_i, temp in enumerate(self.group_data['clopath']['stdp_dt']) if temp==stdp_dt]
            dw_current = self.group_data['clopath']['data'][idx,-1]/self.group_data['clopath']['data'][idx,0]
            self.dw.append(dw_current)
            plt.plot(stdp_dt+input_delay, np.mean(dw_current), marker='.',markersize=20,  color='black')
            plt.errorbar(stdp_dt+input_delay, np.mean(dw_current), yerr=stats.sem(dw_current), color='black')
        plt.axvline(x=0, linewidth=2, linestyle='--', color='black')
        plt.axhline(y=1, linewidth=2, linestyle='--', color='black')
        plt.xlabel(r'$\Delta$ t (ms)', fontsize=20, fontweight='heavy')
        plt.ylabel('Normalized weight', fontsize=20, fontweight='heavy')
        plt.title('STDP window (20 Hz)', fontsize=20, fontweight='heavy')


        plt.show(block=False)
        # condition_types=['polarity', 'syn_num', 'dist','path', 'stdp_dt']
        # polarities=['control']
        # paths=[['1']]
        # syn_nums = [6,8,10]
        # dists = [[0,300],]
        # stdp_dts = [[-10], [10]]
        # condition_combos = [temp for temp in itertools.product(syn_nums, dists)]
        # figs=[]
        # for syn_num_i, syn_num in enumerate(syn_nums):
        #     conditions_current = [polarities, [[syn_num]], [[dists[0]]],[paths[0]], stdp_dts]
        #     conditions = zip(condition_types, conditions_current)

        #     figs.append(
        #         self.Plots._range_variable_group(
        #             group_data=self.group_data,
        #             x_variable='t',
        #             y_variable='clopath',
        #             conditions=conditions
        #             ))
        #     fig_name = directory+'weights_'+str(syn_num)+'syn_'+str(stdp_dts)+'stdp_dt.png'
        #     # plt.show(block=False)
        #     plt.xlabel('time (ms)')
        #     plt.ylabel('weight')
        #     plt.title(str(syn_num)+'synapses, '+str(stdp_dts)+'timing')
        #     figs[-1].savefig(fig_name, dpi=250)
            
    def exp_5a3(self, **kwargs):
        # get list of data files
        directory = 'Data/'+kwargs['experiment']+'/'
        search_string = '*data*'

        group_file_name = 'group.pkl'

        self.Group = analysis.GroupData(directory=directory, file_name=group_file_name, **kwargs)
        self.Plots = analysis.Plots()
        variables = ['v','input_times']

        # updated clopath parameters
        clopath_param= {
        'A_m0':300E-5, # depression magnitude parameter (mV^-1)
        'tetam':-70,#-41, # depression threshold (mV)
        'tetap':-67,#-38, # potentiation threshold (mV)
        'tau_x':8,#-38, # potentiation threshold (mV)
        'tau_m':20,#-38, # potentiation threshold (mV)
        'tau_p': 3, # time constant (ms) for low pass filter post membrane potential for potentiation
        'A_p': 40E-5, # amplitude for potentiation (mV^-2)
        'delay':0, # conduction delay (ms)
        'LTD_delay':1,
        }

        self.group_data = self.Group._standard_run(group_data=self.Group.group_data, directory=directory, file_name=group_file_name, search_string=search_string, clopath_param=clopath_param, variables=variables,)
        self.group_data = self.Group._add_parameter_to_group_data(group_data=self.group_data, parameter='stdp_dt')
        self.Group._save_group_data(group_data=self.group_data, directory=directory, file_name=group_file_name, )

        condition_types=['polarity', 'syn_num', 'dist','path', 'stdp_dt']
        polarities=['control']
        paths=[['1']]
        syn_nums = [6,8,10]
        dists = [[0,300],]
        stdp_dts = [[-10], [10]]
        condition_combos = [temp for temp in itertools.product(syn_nums, dists)]
        figs=[]
        for syn_num_i, syn_num in enumerate(syn_nums):
            conditions_current = [polarities, [[syn_num]], [[dists[0]]],[paths[0]], stdp_dts]
            conditions = zip(condition_types, conditions_current)

            figs.append(
                self.Plots._range_variable_group(
                    group_data=self.group_data,
                    x_variable='t',
                    y_variable='clopath',
                    conditions=conditions
                    ))
            fig_name = directory+'weights_'+str(syn_num)+'syn_'+str(stdp_dts)+'stdp_dt.png'
            # plt.show(block=False)
            plt.xlabel('time (ms)')
            plt.ylabel('weight')
            plt.title(str(syn_num)+'synapses, '+str(stdp_dts)+'timing')
            figs[-1].savefig(fig_name, dpi=250)

exp = Experiment(experiment='exp_1a_dose', rerun=True)       

if __name__=='__main__':
	Experiment(experiment='exp_1a_dose', rerun=True)






