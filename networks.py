import numpy as np
import time as timer

# parameters
#############################################################################
class Param(object):
    '''
    '''
    def __init__(self, **kwargs ):
        '''
        '''
        self.p={}

class ParamOcker2019(Param):
    def __init__(self, **kwargs):
        super(ParamOcker2019, self).__init__(**kwargs)
        self.define_p()

    def define_p(self,):
        '''
        '''
        self.p = {
        'E_L':-70.6, 
        'C':281.,
        'V_Tr':-50.4,
        'V_spike':20.,
        'g_L':30.,
        'delta_T':2.,
        't_wad':144.,
        'a_sub':4.,
        'b_spike':0.805,
        'I_sp':400.,
        't_z':40.,
        't_VT':50.,
        't_ref':2., 
        'V_Tmax':30.4,
        'V_reset':-60.,
        }

class ParamLitwinKumar2014(Param):
    def __init__(self, **kwargs):
        super(ParamLitwinKumar2014, self).__init__(**kwargs)
        self.define_p()

    def define_p(self,):
        '''
        '''
        self.p = {
            # individual cell
            #-------------------
            't_m':20., # membrane time constant: ms
            'E_L':-70., # default leak reversal potential: mV
            'E_L_E':-70., # excitatory leak reversal potential: mV
            'E_L_I':-62., # inhibitory leak reversal potential: mV
            'C':300., # membrane capacitance: pF
            'V_Tr':-52, # resting spike threshold: mV
            'V_spike':20., # peak spike voltage: mV
            'g_L':15., # leak conductance: nS
            'delta_T':2., # EIF slope factor: ms
            't_wad':150., # spike triggered adaptation time scale: ms
            'a_sub':4., # subthreshold adaptation conductance: nS
            'b_spike':0.805, # spike-triggered adaptation current: pA
            # 'I_sp':400., # spike after depolarization current: pA
            'I_sp':0., # spike after depolarization current: pA
            't_z':40., # spiek after depolarization time constant: ms
            't_VT':50., # adaptive threshold time scale: ms
            't_ref':1., # refractory period
            'V_Tmax':-42, # post-spike threshold: mV
            'V_reset':-60., # post-spike reset potential
            'E_rev_E':0., # excitatory synaptic reversal potential: mV
            'E_rev_I':-75., # inhibitory synaptic reversal potential: mV

            # recurrent coupling
            #-------------------
            'N_E': 4000., # number of excitatory neurons: 1
            'N_I':1000., # number of inhibitory neurons: 1
            'p0_EE':0.2, # connection probability for EE connections
            'p0_I':0.2, # connection probability for any I connections
            't_rise_E':1., # rise time constant for E synapses: ms
            't_rise_E':6., # decay time constant for E synapses: ms
            't_rise_E':0.5, # rise time constant for I synapses: ms
            't_rise_E':2., # decay time constant for I synapses: ms 
            'w_EE_min':1.8, # minimum EE synaptic weight: nS
            'w_EE_max':22., # maximum EE synaptic weight: nS
            'w_EE_init':1.8, # initial EE synaptic weight: nS
            'w_EI_min':49., # minimum EI synaptic weight: nS
            'w_EI_max':243., # maximum EI synaptic weight: nS
            'w_EI_init':49., # initial EI synaptic weight: nS
            'w_IE':1.27, # EE synaptic weight, not plastic: nS
            'w_II':16.2, # II synaptic weight, not plastic: nS

            # feedforward input
            #--------------------
            'rate_ff_E': 4.5, # rate of feedforward poisson inputs onto excitatory neurons: kHz
            'rate_ff_I': 2.25, # rate of feedforward poisson inputs onto inhibitory neurons: kHz
        }

        self.p_I = {
            # individual cell
            #-------------------
            't_m':20., # membrane time constant: ms
            'E_L':-62., # default leak reversal potential: mV
            'E_L_E':-70., # excitatory leak reversal potential: mV
            'E_L_I':-62., # inhibitory leak reversal potential: mV
            'C':300., # membrane capacitance: pF
            'V_Tr':-52, # resting spike threshold: mV
            'V_spike':20., # peak spike voltage: mV
            'g_L':15., # leak conductance: nS
            'delta_T':.1, # EIF slope factor: ms
            't_wad':150., # spike triggered adaptation time scale: ms
            'a_sub':0., # subthreshold adaptation conductance: nS
            'b_spike':0., #0.805, # spike-triggered adaptation current: pA
            # 'I_sp':400., # spike after depolarization current: pA
            'I_sp':0., # spike after depolarization current: pA
            't_z':40., # spiek after depolarization time constant: ms
            't_VT':50., # adaptive threshold time scale: ms
            't_ref':1., # refractory period
            'V_Tmax':-52, # post-spike threshold: mV
            'V_reset':-60., # post-spike reset potential
            'E_rev_E':0., # excitatory synaptic reversal potential: mV
            'E_rev_I':-75., # inhibitory synaptic reversal potential: mV
        }

        self.p_clopath= {
            'clopath_A_m0':8E-4, # depression magnitude parameter (pA*mV^-1)
            'clopath_A_p':14E-4, # amplitude for potentiation (pA*mV^-2)
            'clopath_tetam':-70,#-41, # depression threshold (mV)
            'clopath_tetap':-49,#-38, # potentiation threshold (mV)
            'clopath_tau_x':15,#-38, # time constant for presynaptic filter (ms)
            'clopath_tau_m':10,#-38, # time constant for depression lowpass filter
            'clopath_tau_p': 7, # time constant (ms) for low pass filter post membrane potential for potentiation
            
            'clopath_delay':0, # conduction delay (ms)
            'clopath_LTD_delay':0, # conduction delay for LTD  
            'clopath_lower_bound':0.0,
            'clopath_upper_bound':None,#3.0,

            # homeostatic
            #---------------
            'clopath_u_ref' :9,         # reference value for homeostatic process mV
            'clopath_adapt_t':  1000,   # time to integrate over for homeostatic 
            'clopath_homeostatic':False, # apply homeostatic plasticity
            'clopath_w_max':22., # maximum EE synaptic weight: nS
            'clopath_w_min':1.8, # minimum EE synaptic weight: nS
            'clopath_w_init':1.8, # initial EE synaptic weight: nS

            'E_L':-70.,# neuron resting potential, for homeostatic mechanism: mV

        }
        self.p_vogels= {
            't_y': 20,  # time constant for low_pass filtered spike train: ms
            'eta': 1, # learning rate: pA
            'r0': 3, # target firing rate: Hz
            
        }
# networks
#############################################################################
class NetworkLIF(object):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        if 'N_total' in kwargs:
            self.N_total_ph = np.zeros(kwargs['N_total'])

    def build_weight_matrix_random_uniform(self, N_E, N_I, N_total, p0_EE, p0_I, w_max_e, w_max_i):
        '''
        '''
        # binary adjacency matrix
        W0 = np.random.choice([0,1], size=(N_total,N_total), p=[(1-p0_I), p0_I])
        W0_EE = np.random.choice([0,1], size=(N_E,N_E), p=[(1-p0_EE), p0_EE])
        W0[:N_E, :N_E] = W0_EE
        W = np.random.rand(N_total, N_total)
        W[:,:N_E-1] = w_max_e*W[:,:N_E-1]
        W[:,N_E:] = w_max_i*W[:,N_E:]
        W = W*W0
        np.fill_diagonal(W, 0)
        np.fill_diagonal(W0, 0)
        return W, W0

    def neuron_update_exIF(self, i, u, spikes, I_syn, I_ext, p, dt=0.1, **kwargs):
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
        refractory = np.sum(spikes[:,i-int(p['t_ref']/dt):i-1], axis=1, dtype=bool)
        # print refractory.shape
        # voltage update
        #--------------------------------------------------------------------
        du = (1/p['C'])*(

            -p['g_L']*(u[:,i-1]-p['E_L'])  +  
            
            p['g_L']*p['delta_T']*np.exp((u[:,i-1]-p['V_Tr'])/p['delta_T']) - 
            
            I_ext[:, i-1] + 
            
            I_syn[:,i-1]) 
        # print du.shape
        # print u[~refractory,i-1].shape
        # print refractory.shape
        u[~refractory,i] = u[~refractory,i-1] + dt*du[~refractory]

        # reset voltage for cells that cross threshold on previous step
        u[spikes[:,i-1],i] = p['V_reset']

        # get cells that cross on the current step
        spikes[:,i] = u[:,i]>=p['V_spike']

        # set voltage to max if threshold is crossed
        u[spikes[:,i], i] = p['V_spike']

        return u, spikes

    def neuron_update_adexIF(self, i, u, wad, z, V_T, spikes, I_syn, I_ext, p, dt=0.1, **kwargs):
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
            p['g_L']*p['delta_T']*np.exp((u[:,i-1]-V_T[:,i-1])/p['delta_T']) - 
            wad[:,i-1] + 
            z[:,i-1] +
            I_ext[:, i-1] + 
            I_syn[:,i-1]) 
        # print du
        u[:,i] = u[:,i-1] + dt*du

        # reset voltage for cells that cross threshold on previous step
        u[spikes[:,i-1],i] = p['V_reset']

        # get cells that cross on the current step
        spikes[:,i] = u[:,i]>=p['V_spike']

        # set voltage to max if threshold is crossed
        u[spikes[:,i], i] = p['V_spike']

        # adaptation current update
        #-------------------------------------------------------------------
        dwad = (1/p['t_wad'])*( p['a_sub']*( u[:,i-1]-p['E_L']) - wad[:,i-1])
        wad[:,i] = wad[:,i-1] + dt*dwad
        wad[spikes[:,i-1],i] = wad[spikes[:,i-1],i-1] + p['b_spike']

        # spike after-potential
        #-------------------------------------------------------------------
        dz = -(1/p['t_z'])*z[:,i-1]
        z[:,i] = z[:,i-1] + dt*dz
        z[spikes[:,i-1],i] = p['I_sp']

        # threshold
        #------------------------------------------------------------------
        dV_T = (1/p['t_VT'])*(V_T[:,i-1]-p['V_Tr'])
        V_T[:,i] = V_T[:,i-1] - dt*dV_T
        V_T[spikes[:,i-1],i] = p['V_Tmax']


        return u, spikes, wad, z, V_T

    def weight_update_clopath(self, i, x, u, W, W0, u_md, u_mdbar, u_mp, x_m0, A_m, param, dt=0.1, **kwargs):
        """ Determine weights from voltage traces and input times with Clopath rule
        
        ===Args===
        -x      :   input array of spikes (1 for spike, 0 for no spike) (compartments x samples)
        -u      :   array of voltage time traces (compartments x samples)
        -fs     :   sampling rate (kHz)
        -w0     :   initial weight 
        -param  :   dictionary of clopath parameters
        -homeostatic : if True, homeostatic term is included in LTD equation
        
        ===Out===
        -w  :    vector of weights
        
        ===Comments===
        -see Clopath et al. 2010 Nature Neuroscience for details of learning rule
        """
        # remove clopath tag from parameter names
        #-----------------------------------------
        p={}
        for _key, _val in param.iteritems():
            if 'clopath' in _key:
                clopath_key = _key.split('clopath_')[-1]
                p[clopath_key] = _val
            else:
                p[_key]=_val

        # current time in ms
        t = i*dt
            
        # start  after specified delay
        if t>p['delay'] and t>p['LTD_delay']:
                         
            # if include homeostatic LTD mechanism
            if p['homeostatic']:
                
                # adaptation voltage
                if t <= p['adapt_t']:
                    u_mdbar[:,i]   = np.mean(u_md[:,1:i]-p['E_L'],axis=1)
                else:
                    u_mdbar[:,i]   = np.mean(u_md[:,i-int(p['adapt_t']/dt):i-1]-p['E_L'],axis=1)
            
                # homeostatic modulation of LTD rate based on adaptation voltage
                A_m[:,i]   = p['A_m0']*( u_mdbar[:, i-1] **2) /( p['u_ref']**2)   

            else:
                # constant LTD rate
                A_m[:,i]   = p['A_m0']

            # trace of membrane potential (u) with time constant tau_d
            u_md[:,i]  = u_md[ :, i-1] + dt* ( u[ :, i-1]- u_md[ :, i-1])/p['tau_m']
                  
            # trace of membrane potential (u) with time constant tau_p
            u_mp[:,i]  = u_mp[:, i-1]  +dt*( u[:,i-1]-u_mp[:, i-1]) /p['tau_p']
             
            # trace of input spike train (spikes0)
            x_m0[:,i]  = x_m0[:,i-1]  +dt*(x[:,i-1]-x_m0[:,i-1]) /p['tau_x']
            
            # membrane potential (u) thresholded by thetap
            u_sig = (u[:,i] > p['tetap']) *( u[:,i] -p['tetap'])
            
            # membrane potential trace (u_md) thresholded by thetam (taken 3ms before since time of a spike is 2ms)
            u_md_sig  = ( u_md[:, i-int(p['delay']/dt)] > p['tetam']) *( u_md[:, i-int(p['delay']/dt)] -p['tetam'])                  
            
            # membrane potential trace (u_md) thresholded by thetam (taken 3ms before since time of a spike is 2ms)
            u_mp_sig  = ( (u_mp[:,i-int(p['delay']/dt)] -p['tetam']) > 0) *(u_mp[:,i-int(p['delay']/dt)] -p['tetam'])
            
            # update weights
            W[i,:,:] = W0*np.clip(
                W[i-1,:,:] - 
                np.outer(A_m[:,i] *u_md_sig, x[:,i-int(p['LTD_delay']/dt)] ) + 
                np.outer(p['A_p']*u_sig *u_mp_sig, x_m0[:,i]), 
                p['w_min'], 
                p['w_max'])

        return W, W0, u_md, u_mdbar, u_mp, x_m0, A_m

    def weight_update_clopath_temp(self, i, x, u, W, W0, u_md, u_mdbar, u_mp, x_m0, A_m, param, W_LTD_updater, W_LTP_updater, dt=0.1, **kwargs):
        """ Determine weights from voltage traces and input times with Clopath rule
        
        ===Args===
        -x      :   input array of spikes (1 for spike, 0 for no spike) (compartments x samples)
        -u      :   array of voltage time traces (compartments x samples)
        -fs     :   sampling rate (kHz)
        -w0     :   initial weight 
        -param  :   dictionary of clopath parameters
        -homeostatic : if True, homeostatic term is included in LTD equation
        
        ===Out===
        -w  :    vector of weights
        
        ===Comments===
        -see Clopath et al. 2010 Nature Neuroscience for details of learning rule
        """
        # remove clopath tag from parameter names
        #-----------------------------------------
        p={}
        for _key, _val in param.iteritems():
            if 'clopath' in _key:
                clopath_key = _key.split('clopath_')[-1]
                p[clopath_key] = _val
            else:
                p[_key]=_val

        # current time in ms
        t = i*dt
            
        # start  after specified delay
        if t>p['delay'] and t>p['LTD_delay']:
                         
            # if include homeostatic LTD mechanism
            if p['homeostatic']:
                
                # adaptation voltage
                if t <= p['adapt_t']:
                    u_mdbar[:,i]   = np.mean(u_md[:,1:i]-p['E_L'],axis=1)
                else:
                    u_mdbar[:,i]   = np.mean(u_md[:,i-int(p['adapt_t']/dt):i-1]-p['E_L'],axis=1)
            
                # homeostatic modulation of LTD rate based on adaptation voltage
                A_m[:,i]   = p['A_m0']*( u_mdbar[:, i-1] **2) /( p['u_ref']**2)   

            else:
                # constant LTD rate
                A_m[:,i]   = p['A_m0']

            # trace of membrane potential (u) with time constant tau_d
            u_md[:,i]  = u_md[ :, i-1] + dt* ( u[ :, i-1]- u_md[ :, i-1])/p['tau_m']
                  
            # trace of membrane potential (u) with time constant tau_p
            u_mp[:,i]  = u_mp[:, i-1]  +dt*( u[:,i-1]-u_mp[:, i-1]) /p['tau_p']
             
            # trace of input spike train (spikes0)
            x_m0  = x_m0  +dt*(x[:,i-1]-x_m0) /p['tau_x']
            
            # membrane potential (u) thresholded by thetap
            u_sig = (u[:,i] > p['tetap']) *( u[:,i] -p['tetap'])
            
            # membrane potential trace (u_md) thresholded by thetam (taken 3ms before since time of a spike is 2ms)
            u_md_sig  = ( u_md[:, i-int(p['delay']/dt)] > p['tetam']) *( u_md[:, i-int(p['delay']/dt)] -p['tetam'])                  
            
            # membrane potential trace (u_md) thresholded by thetam (taken 3ms before since time of a spike is 2ms)
            u_mp_sig  = ( (u_mp[:,i-int(p['delay']/dt)] -p['tetam']) > 0) *(u_mp[:,i-int(p['delay']/dt)] -p['tetam'])
            
            # update LTD
            np.outer(A_m[:,i] *u_md_sig, x[:,i-int(p['LTD_delay']/dt)], out=W_LTD_updater)
            W  -= W_LTD_updater 
            #update LTP
            np.outer(p['A_p']*u_sig *u_mp_sig, x_m0, out=W_LTP_updater)
            W  += W_LTP_updater
            # enforce bounds
            W = np.clip(W, p['w_min'], p['w_max'])
            # enforce adjacency matrix
            W *= W0
            

        return W, W0, u_md, u_mdbar, u_mp, x_m0, A_m

    def poisson_times(self, rates, dt=.0001, refractory=0., conversion=1000., duration=None, **kwargs):
        '''generate spike times and binary arrays from inhomogenous poisson processes with arbitrary time-varying rates

        ==Args==
        :rates: N x T array: Hz: instantaneous poisson rates.  if 1D, rates are assumed to specify constant rates for N independent processes where T is inferred from the duration argument (see below)
        :dt: scalar: s: time step
        :refractory: scalar: s: refractory period 
        :conversion: scalar: multiplied by spike times to convert to new units 
        :duration: scalar: s: total simulation time in seconds (only used if rates is 1D)

        conversion is used to convert from seconds to miliseconds for compatability with NEURON
        '''
        # get number of neurons (assumed first dimension of rates)
        n_neurons = rates.shape[0]
        # get simulation duration in seconds
        if len(rates.shape)<2:
            if duration is not None:
                T = duration
            else:
                raise 'duration needs to be specified for poisson processes'
        else:
            T = rates.shape[1]*dt
        # create time array in seconds
        time = np.arange(0, T, dt)
        time = np.tile(time, (n_neurons, 1))
        spikes = np.zeros(time.shape)
        # get maximum rate for each process/neuron
        if len(rates.shape)==1:
            rates_max=rates
        else:
            rates_max = np.max(rates, axis=1)
        # nested list to store spike times in 
        spike_times = []
        for nrn in range(n_neurons):
            spike_times.append([])
        # temporary array to store the time of most recent event
        t = np.zeros(n_neurons)
        # if current t is less than total time for any neuron, continue drawing event times
        while np.any(t<=T):
            # generate random numbers
            r=np.random.uniform(size=n_neurons)
            s=np.random.uniform(size=n_neurons)
            # get boolean of neurons that are not refractory
            not_refractory = np.abs(np.log(r)/rates_max)>refractory
            # spike time assuming maximum rate
            t = t-np.log(r)/rates_max
            # tile to compare with full time array 
            t_all = np.tile(t, (time.shape[1], 1)).T
            # corresponding index in time and rate vectors by minimum difference
            i = np.abs(time-t_all).argmin(axis=1)
            # corresponding rate
            idx = [tuple(range(n_neurons)), tuple(i)]
            if len(rates.shape)==1:
                rate=rates
            else:
                rate=rates[idx]
            rate_criteria = s<=rate/rates_max
            for nrn in range(n_neurons):
                if not_refractory[nrn] and rate_criteria[nrn] and t[nrn]<=T:
                    spike_times[nrn].append(conversion*t[nrn])
                    spikes[idx]=1

        return spikes, spike_times

    def poisson_spikes(self, rates, dt=.0001, refractory=0., conversion=1000., duration=None, **kwargs):
        '''generate spike times and binary arrays from inhomogenous poisson processes with arbitrary time-varying rates

        ==Args==
        :rates: N x T array: Hz: instantaneous poisson rates.  if 1D, rates are assumed to specify constant rates for N independent processes where T is inferred from the duration argument (see below)
        :dt: scalar: s: time step
        :refractory: scalar: s: refractory period 
        :conversion: scalar: multiplied by spike times to convert to new units 
        :duration: scalar: s: total simulation time in seconds (only used if rates is 1D)

        conversion is used to convert from seconds to miliseconds for compatability with NEURON
        '''
        # get number of neurons (assumed first dimension of rates)
        n_neurons = rates.shape[0]
        # get simulation duration in seconds
        if len(rates.shape)<2:
            if duration is not None:
                T = duration
            else:
                raise 'duration needs to be specified for poisson processes'
        else:
            T = rates.shape[1]*dt
        # create time array in seconds
        time = np.arange(0, T, dt)
        time = np.tile(time, (n_neurons, 1))
        spikes = np.zeros(time.shape)
        # get maximum rate for each process/neuron
        if len(rates.shape)==1:
            rates_max=rates
        else:
            rates_max = np.max(rates, axis=1)
        # temporary array to store the time of most recent event
        t = np.zeros(n_neurons)
        # if current t is less than total time for any neuron, continue drawing event times
        while np.any(t<=T):
            # generate random numbers
            r=np.random.uniform(size=n_neurons)
            s=np.random.uniform(size=n_neurons)
            # get boolean of neurons that are not refractory
            not_refractory = np.abs(np.log(r)/rates_max)>refractory
            # spike time assuming maximum rate
            t = t-np.log(r)/rates_max
            # tile to compare with full time array 
            t_all = np.tile(t, (time.shape[1], 1)).T
            # corresponding index in time and rate vectors by minimum difference
            i = np.abs(time-t_all).argmin(axis=1)
            # corresponding rate
            idx = [tuple(range(n_neurons)), tuple(i)]
            if len(rates.shape)==1:
                rate=rates
            else:
                rate=rates[idx]
            rate_criteria = s<=rate/rates_max
            for nrn in range(n_neurons):
                if not_refractory[nrn] and rate_criteria[nrn] and t[nrn]<=T:
                    spikes[idx]=1

        return spikes

    def poisson_online(self, i, spikes, next_spike_index, rates, dt=.0001, refractory=0., conversion=1000., duration=None, **kwargs):
        '''generate spike times and binary arrays from inhomogenous poisson processes with arbitrary time-varying rates

        ==Args==
        :rates: N x T array: Hz: instantaneous poisson rates.  if 1D, rates are assumed to specify constant rates for N independent processes where T is inferred from the duration argument (see below)
        :dt: scalar: s: time step
        :refractory: scalar: s: refractory period 
        :conversion: scalar: multiplied by spike times to convert to new units 
        :duration: scalar: s: total simulation time in seconds (only used if rates is 1D)

        conversion is used to convert from seconds to miliseconds for compatability with NEURON
        '''
        # get number of neurons (assumed first dimension of rates)
        n_neurons = rates.shape[0]
        # get simulation duration in seconds
        if len(rates.shape)<2:
            if duration is not None:
                T = duration
            else:
                raise 'duration needs to be specified for poisson processes'
        else:
            T = rates.shape[1]*dt

        t= i*dt
        # if i is zero, generate initial spike times
        if i==0:
            events = np.zeros(next_spike_index.shape)
            while np.any(events==0):
                # find inputs that just spiked
                idx_need_spikes = np.where(events==0)
                # generate random numbers
                r=np.random.uniform(size=len(idx_need_spikes))
                s=np.random.uniform(size=len(idx_need_spikes))
                # get boolean of neurons that are not refractory
                not_refractory = np.abs(np.log(r)/rates_max)>refractory
                # spike time assuming maximum rate
                



                t = t-np.log(r)/rates_max




        # create time array in seconds
        time = np.arange(0, T, dt)
        time = np.tile(time, (n_neurons, 1))
        spikes = np.zeros(time.shape)
        # get maximum rate for each process/neuron
        if len(rates.shape)==1:
            rates_max=rates
        else:
            rates_max = np.max(rates, axis=1)
        # nested list to store spike times in 
        spike_times = []
        for nrn in range(n_neurons):
            spike_times.append([])
        # temporary array to store the time of most recent event
        t = np.zeros(n_neurons)
        # if current t is less than total time for any neuron, continue drawing event times
        while np.any(t<=T):
            # generate random numbers
            r=np.random.uniform(size=n_neurons)
            s=np.random.uniform(size=n_neurons)
            # get boolean of neurons that are not refractory
            not_refractory = np.abs(np.log(r)/rates_max)>refractory
            # spike time assuming maximum rate
            t = t-np.log(r)/rates_max
            # tile to compare with full time array 
            t_all = np.tile(t, (time.shape[1], 1)).T
            # corresponding index in time and rate vectors by minimum difference
            i = np.abs(time-t_all).argmin(axis=1)
            # corresponding rate
            idx = [tuple(range(n_neurons)), tuple(i)]
            if len(rates.shape)==1:
                rate=rates
            else:
                rate=rates[idx]
            rate_criteria = s<=rate/rates_max
            for nrn in range(n_neurons):
                if not_refractory[nrn] and rate_criteria[nrn] and t[nrn]<=T:
                    spike_times[nrn].append(conversion*t[nrn])
                    spikes[idx]=1

        return spikes, spike_times

    def poisson_fixed_rate(self, rate, N, T, dt):
        '''
        '''
        t = np.arange(0, T, dt)
        p0 = rate*T/len(t)
        spikes = np.random.choice([0,1], size=(N,len(t)), p=[(1-p0), p0])
        return spikes

class NetworkLitwinKumar2014(NetworkLIF):
    '''
    '''
    def __init__(self, **kwargs):
        super(NetworkLitwinKumar2014, self).__init__(**kwargs)

    def run(self, ):
        '''
        '''
        P = ParamLitwinKumar2014()
        p = P.p

        # neuron groups
        N_E = 1000
        N_I= 250
        N_total = N_E+N_I

        # generate weight matrix
        p0_EE = 0.2
        p0_I = 0.4
        W = self.build_weight_matrix_random_uniform(N_E=N_E, N_I=N_I, N_total=N_total, p0_EE=p0_EE, p0_I=p0_I)

        # time
        T = 1000
        dt=0.1
        time = np.arange(0, T, dt)
        delay = 1.

        # applied current
        self.I_ext = 10 + 1000.*np.random.rand(N_total, time.shape[0])
        self.I_ext[10:,:]=0.
        # self.I_ext = 5000*np.ones((N_total, time.shape[0]))
        I_syn = np.zeros((N_total, time.shape[0]))
        g_syn = np.zeros((N_total, time.shape[0]))
        E_rev = np.zeros((N_total,))

        # initialize internal neuron state variables
        u = p['E_L']*np.ones((N_total,time.shape[0]))
        spikes = np.zeros((N_total,time.shape[0]), dtype=bool)
        wad = np.zeros((N_total,time.shape[0]))
        z = np.zeros((N_total,time.shape[0]))
        V_T = p['V_Tr']*np.ones((N_total,time.shape[0]))

        for i, t in enumerate(time):
            if t>delay:
                # adaptive exponential neuron voltage
                u, spikes, wad, z, V_T = self.adex_update(i=i, u=u, wad=wad, z=z, V_T=V_T, spikes=spikes, I_syn=I_syn, I_ext=self.I_ext, p=p, dt=dt)
                # update synaptic currents
                g_syn, I_syn  = self.I_syn_single_exp_update(i=i, u=u, g=g_syn, I_syn=I_syn, spikes=spikes, W=W, E_rev=E_rev, tau_g=10, delay=delay, dt=dt)

        self.u = u
    
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

    def I_syn_single_exp_update(self, i, u, g, I_syn, spikes, W, E_rev, tau_g=10, delay=1, dt=0.1):
        '''
        '''
        # print i-1-int(delay/dt)
        # print E_rev.shape
        if len(W.shape)==2:
            # g[:,i] = g[:,i-1] -dt*g[:,i-1]/tau_g + W.dot(spikes[:,i-1-int(delay/dt)])
            g[:, i] = g[:,i-1] -dt*g[:,i-1]/tau_g
            W.dot(spikes[:,i-1-int(delay/dt)], out=self.N_total_ph)
            g[:,i] += self.N_total_ph
            # g[:,i] = g[:,i-1] -dt*g[:,i-1]/tau_g + W.dot(spikes[:,i-1-int(delay/dt)])
        if len(W.shape)==3:
            g[:,i] = g[:,i-1] -dt*g[:,i-1]/tau_g + W[i-1, :,:].dot(spikes[:,i-1-int(delay/dt)])


        I_syn[:,i] = -g[:,i]*(u[:,i]-E_rev)
        return g, I_syn

    def adex_update(self, i, u, wad, z, V_T, spikes, I_syn, I_ext, p, dt=0.1, **kwargs):
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
            p['g_L']*p['delta_T']*np.exp((u[:,i-1]-V_T[:,i-1])/p['delta_T']) - 
            wad[:,i-1] + 
            z[:,i-1] +
            I_ext[:, i-1] + 
            I_syn[:,i-1]) 
        # print du
        u[:,i] = u[:,i-1] + dt*du

        # reset voltage for cells that cross threshold on previous step
        u[spikes[:,i-1],i] = p['V_reset']

        # get cells that cross on the current step
        spikes[:,i] = u[:,i]>=p['V_spike']

        # set voltage to max if threshold is crossed
        u[spikes[:,i], i] = p['V_spike']

        # adaptation current update
        #-------------------------------------------------------------------
        dwad = (1/p['t_wad'])*( p['a_sub']*( u[:,i-1]-p['E_L']) - wad[:,i-1])
        wad[:,i] = wad[:,i-1] + dt*dwad
        wad[spikes[:,i-1],i] = wad[spikes[:,i-1],i-1] + p['b_spike']

        # spike after-potential
        #-------------------------------------------------------------------
        dz = -(1/p['t_z'])*z[:,i-1]
        z[:,i] = z[:,i-1] + dt*dz
        z[spikes[:,i-1],i] = p['I_sp']

        # threshold
        #------------------------------------------------------------------
        dV_T = (1/p['t_VT'])*(V_T[:,i-1]-p['V_Tr'])
        V_T[:,i] = V_T[:,i-1] - dt*dV_T
        V_T[spikes[:,i-1],i] = p['V_Tmax']


        return u, spikes, wad, z, V_T

class NetworkOcker2019(object):
    '''
    '''
    def __init__(self, ):
        '''
        '''
        pass

    def eif_update(self, i, u, spikes, I_syn, I_ext, p, dt=0.1, **kwargs):
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
        refractory = np.sum(spikes[:,i-int(p['t_ref']/dt):i-1], axis=1, dtype=bool)
        # print refractory.shape
        # voltage update
        #--------------------------------------------------------------------
        du = (1/p['C'])*(

            -p['g_L']*(u[:,i-1]-p['E_L'])  +  
            
            p['g_L']*p['delta_T']*np.exp((u[:,i-1]-p['V_Tr'])/p['delta_T']) - 
            
            I_ext[:, i-1] + 
            
            I_syn[:,i-1]) 
        # print du
        u[~refractory,i] = u[~refractory,i-1] + dt*du

        # reset voltage for cells that cross threshold on previous step
        u[spikes[:,i-1],i] = p['V_reset']

        # get cells that cross on the current step
        spikes[:,i] = u[:,i]>=p['V_spike']

        # set voltage to max if threshold is crossed
        u[spikes[:,i], i] = p['V_spike']

        return u, spikes


# experiments
#############################################################################
class Exp(object):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        pass

class exp_litwinkumar_test(Exp):
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        super(exp_litwinkumar_test, self).__init__(**kwargs)

    def run(self, **kwargs):
        '''
        '''
        P = ParamLitwinKumar2014()
        p = P.p

        
        # neuron groups
        #----------------------------------
        N_E = 1000
        N_I= 250
        N_total = N_E+N_I

        # create network object
        net = NetworkLitwinKumar2014(N_total=N_total)

        # generate recurrent weight matrix
        #-----------------------------------
        p0_EE = p['p0_EE']
        p0_I = p['p0_I']
        w_max_e=5#E-1
        w_ff=1.2#1E-1
        w_max_i=50#E-1
        W, W0 = net.build_weight_matrix_random_uniform(N_E=N_E, N_I=N_I, N_total=N_total, p0_EE=p0_EE, p0_I=p0_I, w_max_e=w_max_e, w_max_i=w_max_i)
        W_LTD_updater = np.zeros(W.shape)
        W_LTP_updater = np.zeros(W.shape)

        # time
        #-----------------------------
        T = 50. #(ms)
        dt=0.1 #(ms)
        time = np.arange(0, T, dt)
        nsamples = time.shape[0]
        delay = 1.

        # extend weight matrix in time (samples x post x pre)
        #----------------------------------------------------
        W = np.tile(W, (time.shape[0], 1, 1))

        # generate feedforward spike trains
        #-----------------------------
        # rate of feedforward poisson input (kHz)
        rate_ff_E = 2.5
        rate_ff_I = 2.25
        spikes_input_E = net.poisson_fixed_rate(rate=rate_ff_E, N=N_E, T=T, dt=dt)
        spikes_input_I = net.poisson_fixed_rate(rate=rate_ff_I, N=N_I, T=T, dt=dt)
        self.spikes_input = np.append(spikes_input_E, spikes_input_I, axis=0)
        # self.spikes_input = net.poisson_fixed_rate(rate=rate_ff, N=N_total, T=T, dt=dt)

        # generate feedforward weight matrix
        #-------------------------------------
        # start with overall weight matrix with all zeros and fill diagonal
        self.W_in = np.zeros((N_total, N_total), dtype=float)
        np.fill_diagonal(self.W_in, w_ff)

        # applied external current
        #-------------------------------------------
        # self.I_ext = 10 + 1000.*np.random.rand(N_total, time.shape[0])
        I_ext = np.zeros((N_total, time.shape[0]))

        # initialize synaptic conductance/current
        #-----------------------------------------
        I_syn_E = np.zeros((N_total, time.shape[0]))
        g_syn_E = np.zeros((N_total, time.shape[0]))
        I_syn_I = np.zeros((N_total, time.shape[0]))
        g_syn_I = np.zeros((N_total, time.shape[0]))
        I_syn_ff = np.zeros((N_total, time.shape[0]))
        g_syn_ff = np.zeros((N_total, time.shape[0]))

        # initialize internal neuron state variables
        #----------------------------------------------
        u = p['E_L']*np.ones((N_total,time.shape[0]))
        spikes = np.zeros((N_total,time.shape[0]), dtype=bool)
        wad = np.zeros((N_total,time.shape[0]))
        z = np.zeros((N_total,time.shape[0]))
        V_T = p['V_Tr']*np.ones((N_total,time.shape[0]))

        # initialize clopath state variables for recurrent EE synapses
        #---------------------------------------------------
        # filtered membrane potential for the depression
        u_md        = p['E_L']*np.ones((N_E, nsamples)) 
        # integrated depression potential for homeostatic plasticity
        u_mdbar     = np.zeros((N_E, nsamples))    
        # filtered membrane potential for potentiation     
        u_mp        = p['E_L']*np.ones((N_E, nsamples))
        # filtered presynaptic trace 
        x_m0        = np.zeros((N_E))     
        # homeostatic LTD variable    
        A_m         = P.p_clopath['clopath_A_m0']*np.ones((N_E,nsamples))

        start = timer.time()
        # main loop for numerical integration
        for i, t in enumerate(time):
            if t>delay:
                # combine synaptic inputs from feedforward and recurrent inputs
                I_syn_total = I_syn_E+I_syn_I+I_syn_ff
                # excitatory adaptive exponential neuron state variables
                #-------------------------------------
                if i == 100:
                    start_temp = timer.time()
                u[:N_E,:], spikes[:N_E,:], wad[:N_E,:], z[:N_E,:], V_T[:N_E,:] = net.adex_update(i=i, u=u[:N_E,:], wad=wad[:N_E,:], z=z[:N_E,:], V_T=V_T[:N_E,:], spikes=spikes[:N_E,:], I_syn=I_syn_total[:N_E,:], I_ext=I_ext[:N_E,:], p=p, dt=dt)
                if i == 100:
                    end_temp = timer.time()
                    print 'neuron duration:', str(end_temp-start_temp)
                # inhibitory adaptive exponential neuron state variables
                #-------------------------------------
                # u[N_E:,:], spikes[N_E:,:], wad[N_E:,:], z[N_E:,:], V_T[N_E:,:] = net.adex_update(i=i, u=u[N_E:,:], wad=wad[N_E:,:], z=z[N_E:,:], V_T=V_T[N_E:,:], spikes=spikes[N_E:,:], I_syn=I_syn_total[N_E:,:], I_ext=I_ext[N_E:,:], p=P.p_I, dt=dt)
                if i == 100:
                    start_temp = timer.time()

                u[N_E:,:], spikes[N_E:,:] = net.neuron_update_exIF(i=i, u=u[N_E:,:], spikes=spikes[N_E:,:], I_syn=I_syn_total[N_E:,:], I_ext=I_ext[N_E:,:], p=P.p_I, dt=dt)
                if i == 100:
                    end_temp = timer.time()
                    print 'neuron inhibitory duration:', str(end_temp-start_temp)

                # update recurrent synaptic currents
                #------------------------------------
                # excitatory
                if i == 100:
                    start_temp = timer.time()
                g_syn_E, I_syn_E = net.I_syn_single_exp_update(i=i, u=u, g=g_syn_E, I_syn=I_syn_E, spikes=spikes[:N_E, :], W=W[0, :, :N_E], E_rev=p['E_rev_E'], tau_g=10, delay=delay, dt=dt)
                if i == 100:
                    end_temp = timer.time()
                    print 'excitatory synaptic duration:', str(end_temp-start_temp)
                # inhibitory
                if i == 100:
                    start_temp = timer.time()
                g_syn_I, I_syn_I = net.I_syn_single_exp_update(i=i, u=u, g=g_syn_I, I_syn=I_syn_I, spikes=spikes[N_E:, :], W=W[0, :, N_E:], E_rev=p['E_rev_I'], tau_g=10, delay=delay, dt=dt)
                if i == 100:
                    end_temp = timer.time()
                    print 'inhibitory synaptic duration:', str(end_temp-start_temp)
                # update feedforward excitatory synaptic currents
                #--------------------------------------
                if i == 100:
                    start_temp = timer.time()
                g_syn_ff, I_syn_ff  = net.I_syn_single_exp_update(i=i, u=u, g=g_syn_ff, I_syn=I_syn_ff, spikes=self.spikes_input, W=self.W_in, E_rev=p['E_rev_E'], tau_g=10, delay=delay, dt=dt)
                if i == 100:
                    end_temp = timer.time()
                    print 'feedforawrd synaptic duration:', str(end_temp-start_temp)

                # update recurrent EE weights with clopath rule
                #------------------------------------------------
                
                # W[0,:N_E, :N_E], W0, u_md, u_mdbar, u_mp, x_m0, A_m = net.weight_update_clopath_temp(i=i, x=spikes[:N_E,:], u=u[:N_E,:], W=W[0,:N_E, :N_E], W0=W0[:N_E, :N_E], u_md=u_md, u_mdbar=u_mdbar, u_mp=u_mp, x_m0=x_m0, A_m=A_m, param=P.p_clopath, W_LTP_updater=W_LTP_updater[:N_E, :N_E], W_LTD_updater=W_LTD_updater[:N_E, :N_E], dt=0.1,)
                
        end = timer.time()
        print 'simulation duration:', str(end-start)
        self.u = u
        self.spikes =spikes
        self.W = W
