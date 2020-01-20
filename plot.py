from brian2 import *
import numpy as np
import pandas as pd
from scipy import stats
import copy




def _plot_all(group_dataframe, variable):
    '''
    '''
    gf = group_dataframe[variable]
    plt.figure()
    plt.plot(np.array(gf[gf.field_mag==0].tolist()).T, color='black')
    plt.plot(np.array(gf[gf.field_mag>0].tolist()).T, color='red')
    plt.show(block=False)



class Plot:
    '''
    '''
    def __init__(self,):
        '''
        '''

    def _plot_voltage(self, rec):
        '''
        '''
        self.u  = rec['nrn'].u
        linewidth=4
        colors = ['r','k','b']
        plt.figure()
        for i in range(self.u.shape[0]):
            plt.plot(rec['nrn'].t/ms, self.u[i,:]/mV, color=colors[i], linewidth=linewidth)

        plt.xlabel('Time (ms)', fontsize=20, fontweight='bold')
        plt.ylabel('Membrane potential (mV)', fontsize=20, fontweight='bold')
        plt.ylim([-72,-20])
        plt.show(block=False)


    def _plot_weights(self, rec, npre, npost):
        '''
        '''
        self.w = rec['syn'].w_clopath
        self.w = self.w.reshape(npre,npost,self.w.shape[1])
        markers = ['-','--']
        colors = ['r','k','b']
        linewidth=4
        plt.figure()
        for pre_i in range(npre-1):
            for post_i in range(npost):
                plt.plot(rec['syn'].t/ms, self.w[pre_i,post_i,:], color=colors[post_i], linestyle=markers[pre_i], linewidth=linewidth)
        plt.xlabel('Time (ms)', fontsize=20, fontweight='bold')
        plt.ylabel('Synaptic Weight (AU)', fontsize=20, fontweight='bold')
        plt.ylim([0.4,1.5])
        plt.show(block=False)
    # FIXME
    def _rate_x_location(self, rec):
        '''
        '''
        tol = .01
        integration_window = 1000
        firing_rate_filter = np.ones(integration_window)
        locations=rec['locations']

        anodal_idx = rec['nrn_spikes'].i == 0
        control_idx = rec['nrn_spikes'].i == 1
        cathodal_idx = rec['nrn_spikes'].i == 2
        anodal_spikes = rec['nrn_spikes'].t[anodal_idx]
        control_spikes = rec['nrn_spikes'].t[control_idx]
        cathodal_spikes = rec['nrn_spikes'].t[cathodal_idx]
        time = np.linspace(0, locations.shape[0]*defaultclock.dt, locations.shape[0])
        locations_unique = np.array(list(set(np.round(locations, 0))))
        print locations_unique.shape[0],time.shape[0]
        spike_array_anodal = np.zeros((locations_unique.shape[0],time.shape[0]  ))
        
        for spike_i, spike in enumerate(anodal_spikes):
            # print np.round(spike/ms,2)
            # print np.round(time[10]/ms,2)
            # time_i = np.where(np.round(time/ms,2)==np.round(spike/ms,2))[0]
            # time_i = np.where(abs(time/ms-spike/ms) < tol)
            time_i = np.argmin(abs(time/ms-spike/ms))
            location = locations[time_i]
            # print np.round(spike/ms,2), time_i, location, min(abs(time/ms-spike/ms))

            location_unique_i = np.argmin(abs(locations_unique-location))
            print location_unique_i, time_i
            # print location_unique_i, time_i
            spike_array_anodal[location_unique_i,time_i]=1

        firing_rate_array_anodal = np.zeros(spike_array_anodal.shape)
        for loc_i in range(spike_array_anodal.shape[0]):
            firing_rate_array_anodal[loc_i, :] = np.convolve(spike_array_anodal[loc_i,:], firing_rate_filter, mode='same')


        rec['time'] = time
        rec['locations_unique'] = locations_unique
        rec['anodal_spikes'] = anodal_spikes
        rec['cathodal_spikes'] =  cathodal_spikes
        rec['control_spikes'] = control_spikes

        plt.figure()
        plt.imshow(firing_rate_array_anodal, aspect='auto')
        plt.colorbar()
        plt.show(block=False)

        return spike_array_anodal

    def _time_series(self, rec, variables, indices):
        '''
        ==Args==
        -rec : dictionary of recorded StateMonitor objects
        -variables : list of variables to plot
        ==Out==
        ==Update==
        ==Comments==
        '''
        # object for dimensionless variable
        dimensionless = units.fundamentalunits.DIMENSIONLESS
        # map dimension to specific unit
        dimension_map = [(volt, mV, 'mV'), (amp,pA,'pA'), (siemens,nS,'nS'), (dimensionless,1, 'AU')]
        # iterate over variables
        for variable_i, variable in enumerate(variables):
            # iterate through recorded data types
            for rec_type, rec_data in rec.iteritems():

                if isinstance(rec_data, StateMonitor) and variable in rec_data.recorded_variables.keys():
                    figure()
                    plot_time = getattr(rec_data, 't')
                    plot_data = getattr(rec_data, variable)
                    plot_dimensions = get_dimensions(plot_data)
                    print plot_dimensions
                    plot_units = [temp[1] for temp_i, temp in enumerate(dimension_map) if plot_dimensions==get_dimensions(temp[0])][0]
                    print plot_units
                    plot(plot_time/ms, plot_data[indices[variable_i]].T/plot_units,)
                    show(block=False)

                elif isinstance(rec_data, SpikeMonitor) and variable in list(rec_data.record_variables):
                    figure()
                    plot_time = getattr(rec_data, 't')
                    plot_data = getattr(rec_data, variable)
                    plot_dimensions = get_dimensions(plot_data)
                    plot_units = [temp[1] for temp_i, temp in enumerate(dimension_map) if plot_dimensions==temp[0]][0] 
                    plot(plot_time/ms, plot_data[indices[variable_i]].T/plot_units)
                    show(block=False)