'''
'''
from brian2 import *
import numpy as np
import pandas as pd
from scipy import stats
import clopath

prefs.codegen.target = 'numpy'

class Param:
    '''
    '''
    def __init__(self, **kwargs):
        '''
        '''
        self.p = self._param(**kwargs)
        self.initial = self._initial(self.p)

    def _param(self, **kwargs):
        '''
        '''

        p = {
            'N_network':1,
            'N_input':1,

            'dt':0.1*ms,
            'run_time':1000*ms,
            # neuron
            #==============================
            'E_L':-70*mV,
            'g_L':40*nS,
            'C':281*pF,
            'delta_T':2*mV,
            't_noise':20*ms,
            't_V_T' : 50*ms,
            'V_Trest' : -55*mV,
            'V_Tmax':-30*mV,
            'threshold_condition':'u>V_T=20*mV',
            'refractory_time':2*ms,
            'u_hold':30*mV,
            'u_reset':-70*mV,
            'I_after' : 400*pA,
            'a_adapt': 4*nS,
            'b_adapt':0.805*pA,
            't_w_adapt':144*ms,
            't_z_after':40*ms,

            'I_input':0*pA,
            'I_field':0*pA,

            't_reset':0.5*ms,
            'spike_hold_time':1*ms,
            'spike_hold_time2': 2*ms - 2*defaultclock.dt,

            # synapse
            #=======================================
            'g_max_ampa': 100*nS,
            't_ampa':2*ms,
            'E_ampa':0*mV,
            'w_ampa':0.5,

            'g_max_nmda':50*nS,
            't_nmda':50*ms,
            'E_nmda':0*mV,
            'w_nmda':0.5,

            # clopath rule
            'tau_lowpass1':20*ms,
            'tau_lowpass2':10*ms,
            'tau_homeo':1000*ms,
            'theta_low' : -65,
            'theta_high':-60,
            'A_LTD': 2E-5,
            'A_LTP':10*40E-5,
            'v_target': 100*mV*mV,
            
            'hold_spike':1,
            'update_ampa_online':0,
            'update_nmda_online':0,
            'homeo_on':1,

            'x_reset':1,
            'w_max_clopath':2,

            'f':5.3,
            't_F':94*ms,
            'd1':0.45,
            't_D1':540*ms,
            'd2':0.12,
            't_D2':45*ms,
            'd3':0.98,
            't_D3':120E3*ms,

            'warmup':10,
            'pulses':4,
            'bursts':1,
            'pulse_freq':100,
            'burst_freq':5
            }

        for kw, val in kwargs.iteritems():
            if kw in p.keys():
                p[kw]=val

        return p

    def _initial(self, p):
        '''
        '''
        initial = {
        'u':p['E_L'],
        'V_T':p['V_Trest'],
        'w_adapt':0*pA,
        'z_after':0*pA,
        'g_ampa':0*nS,
        'g_nmda':0*nS,
        'u_lowpass1':0*mV,
        'u_lowpass2':0*mV,
        'u_homeo':0*mV,
        'x_trace':0,
        'w_clopath':0.5,

        'F':1,
        'D1':1,
        'D2':1,
        'D3':1,
        }
        return initial

class Eq:
    '''
    equation definitions
    '''
    def __init__(self, ):
        '''
        '''
        self.model_nrn_hold_spike = '''
            # membrane voltage (on spike, hold potential at peak)
            #===========================================================
            
            # membrane voltage
            #``````````````````
                du/dt = int(not_refractory)*I_total/C  + I_reset  + ((t-lastspike)>spike_hold_time2)*(1-int(not_refractory))*I_total/C :volt

            # total current 
            #```````````````````
                I_total = (I_L + I_syn + I_exp + I_field + I_input - w_adapt + z_after ) : amp

            # reset current in volt/sec 
            #``````````````````````````
                I_reset = ((t-lastspike)<spike_hold_time2)*((t-lastspike)>spike_hold_time)*(1-int(not_refractory))*((u_reset-u)/t_reset + z_after/C) : volt/second
            
            # leak current
            #``````````````````````````` 
                I_L = -g_L*(u - E_L) : amp

            # exponential current
            #`````````````````````````````
                I_exp = g_L*delta_T*exp((u - V_T)/delta_T) : amp

            # threshold adaptation
            #`````````````````````````````
                dV_T/dt = -(V_T-V_Trest)/t_V_T : volt

            # hyperpolarizing spike frequency adaptation current
            #```````````````````````````````````````````````
                dw_adapt/dt = a_adapt*(u-E_L)/t_w_adapt - w_adapt/t_w_adapt : amp

            # depolarizing spike after potential
            #`````````````````````````````````````
                dz_after/dt = -z_after/t_z_after : amp 

            # synaptic
            #=============================================================
            # ampa
            #``````````````
                dg_ampa/dt = -g_ampa/t_ampa : siemens 
                I_ampa = -g_ampa*(u-E_ampa) : amp

            # nmda
            #`````````````````
                dg_nmda/dt = -g_nmda/t_nmda : siemens 
                B =  1/(1 + exp(-0.062*u/mV)/3.57) : 1 
                I_nmda = -g_nmda*B*(u-E_nmda) : amp

            # total
            #``````````````````````
                I_syn = I_ampa + I_nmda : amp

            # clopath
            #===============================================================
            # low threshold filtered membrane potential
            #```````````````````````````````````````````
                du_lowpass1/dt = (u-u_lowpass1)/tau_lowpass1 : volt 

            # high threshold filtered membrane potential
            #``````````````````````````````````````````````
                du_lowpass2/dt = (u-u_lowpass2)/tau_lowpass2 : volt     

            # homeostatic term
            #````````````````````````````````````````````````
                du_homeo/dt = (u-E_L-u_homeo)/tau_homeo : volt       

            # LTP voltage dependence
            #````````````````````````````````````````````````
                LTP_u = (u_lowpass2/mV - theta_low/mV)*int((u_lowpass2/mV - theta_low/mV) > 0)*(u/mV-theta_high/mV)*int((u/mV-theta_high/mV) >0)  : 1

            # LTD voltage dependence
            #``````````````````````````````````````````````````
                LTD_u = (u_lowpass1/mV - theta_low/mV)*int((u_lowpass1/mV - theta_low/mV) > 0)  : 1

            # homeostatic depression amplitude
            #``````````````````````````````````
                A_LTD_homeo = A_LTD*(u_homeo**2/v_target) : 1  
            
            # parameters
            #===============================================================
            # neuron
            #```````````````````````````````````````````````````````````````
            I_input : amp           # arbitrary current injection
            I_field : amp           # extracellular field current
            I_after : amp           # after spike depolarizing (see z_after)
            C : farad               # membrane capacitance
            g_L : siemens           # leak conductance
            delta_T : volt          # steepness of exponential current (see I_exp)
            t_V_T : second          # treshold adaptation time constant (see V_T)
            a_adapt : siemens       # amplitude of hyperpolarizing adaptation current (see w_adapt)
            t_w_adapt : second      # hyperpolarizing adaptation time constant (see w_adapt)
            t_z_after : second      # depolarizing afterpotential time constant (see z_after)
            u_reset : volt          # reset voltage after spike
            b_adapt : amp           # increment for hyperpolarizing current after spike
            V_Tmax : volt           # spike threshold reset value after spike
            V_Trest: volt           # equillibrium spike threshold
            E_L : volt              # leak potential
            t_reset : second
            spike_hold_time : second
            spike_hold_time2 : second

            # synapse
            #```````````````````````````````````````````````````````````````
            g_max_ampa : siemens
            t_ampa : second
            E_ampa : volt
            g_max_nmda : siemens
            t_nmda : second
            E_nmda : volt
            tau_lowpass1 : second
            tau_lowpass2 : second
            tau_homeo : second
            theta_low : volt
            theta_high : volt
            A_LTD:1
            v_target:volt*volt

            # boolean simulation options
            #```````````````````````````````````````````````````````````````
            hold_spike : 1          # boolean indicating to hold spike at peak potential during spike
            update_ampa_online : 1  # boolean indicating to make ampa conductance plastic during simulation
            update_nmda_online : 1  # boolean indicating to make nmda conductance plastic during simulation
            homeo_on : 1            # boolean indicating to use homeostatic clopath rule
        '''

        self.reset_nrn ='''
            z_after = I_after 
            u = int(hold_spike)*u_hold + int(1-hold_spike)*(u_reset + dt*I_after/C)
            V_T = V_Tmax 
            w_adapt += b_adapt    
        '''

        self.pre_syn_ampa = '''
        
            g_ampa += int(update_ampa_online_post)*w_clopath*g_max_ampa*A + int(1-update_ampa_online_post)*w_ampa*g_max_ampa*A 
        '''

        self.pre_syn_nmda = '''
        
            g_nmda += int(update_nmda_online_post)*w_clopath*g_max_nmda*A + int(1-update_nmda_online_post)*w_nmda*g_max_nmda*A 
        '''

        self.pre_syn_clopath = '''

            w_minus = int(homeo_on_post)*A_LTD_homeo_post*LTD_u_post + int(1-homeo_on_post)*A_LTD_post*LTD_u_post

            w_clopath = clip(w_clopath-w_minus, 0, w_max_clopath)  # apply LTD

            x_trace += dt*x_reset/tau_x  # update presynaptic trace with each input
        '''

        self.model_syn_clopath = '''      
        
            # lowpass presynaptic variable
            #`````````````````````````````````````
                dx_trace/dt = -x_trace/tau_x : 1                          

            # clopath rule for potentiation (depression implented with on_pre)
            #```````````````````````````````````````
                dw_clopath/dt = int(w_clopath<w_max_clopath)*A_LTP*x_trace*LTP_u_post :1
                # dw_clopath/dt = saturated*(w_max_clopath-w_clopath)/dt + (1-saturated)*A_LTP*x_trace*LTP_u_post : 1

                # saturated = int((w_clopath+A_LTP*x_trace*LTP_u_post*dt)>w_max_clopath) : 1 # indicates that next weight update brings synapse to saturation

            # homeostatic depression amplitude
            #`````````````````````````````````````````
                # A_LTD_u = A_LTD*(u_homeo_post**2/v_target) : 1  

            # parameters
            #====================================================================
                w_ampa:1
                tau_x : second
                A_LTP: hertz
                x_reset:1
                w_max_clopath:1
        '''
        # short term plasticity
        #===================================================================
        self.model_syn_stp = '''

            dF/dt = (1-F)/t_F : 1 
            dD1/dt = (1-D1)/t_D1 : 1 
            dD2/dt = (1-D2)/t_D2 : 1 
            dD3/dt = (1-D3)/t_D3 : 1 
            A = F*D1*D2*D3 : 1
            
            # parameters
            #===================================================================
            t_F : second
            f : 1
            t_D1 : second
            d1 : 1
            t_D2 : second
            d2 : 1
            t_D3 : second
            d3 : 1
        '''

        self.pre_syn_stp = '''
            F += f 
            D1 *= d1
            D2 *= d2
            D3 *= d3 
        '''

    def _add_eq(self, *equations):
        '''
        '''
        new_line = '\n'

        joined_equations = new_line.join(equations)

        return joined_equations

class Input:
    '''
    '''
    def __init__(self, ):
        '''
        '''
        pass

    def _poisson(self, N=1, rates=100):
        '''
        '''
        # poisson input group
        #===================================================================
        self.input_poisson = PoissonGroup(N=N, rates=rates*Hz, )

        return self.nrn

    def _tbs(self, p):
        '''
        '''
        warmup = p['warmup']
        pulses = p['pulses']
        bursts = p['bursts']
        pulse_freq=p['pulse_freq']
        burst_freq=p['burst_freq']
        input_times = np.zeros(pulses*bursts)*ms
        indices = np.zeros(pulses*bursts)
        cnt=-1
        for burst in range(bursts):
            for pulse in range(pulses):
                cnt+=1
                time = warmup + 1000*burst/burst_freq + 1000*pulse/pulse_freq
                input_times[cnt] = time*ms
            
        self.input_tbs = SpikeGeneratorGroup(1, indices, input_times )
        
        return self.input_tbs

def _run():

    P = Param()
    Eqs = Eq()
    defaultclock.dt = P.p['dt']
    saturated=0

    eq_nrn = Eqs.model_nrn_hold_spike
    eq_nrn_reset = Eqs.reset_nrn 
    eq_syn = Eqs._add_eq(Eqs.model_syn_clopath, Eqs.model_syn_stp)
    eq_syn_pre = Eqs._add_eq(Eqs.pre_syn_ampa, Eqs.pre_syn_nmda, Eqs.pre_syn_clopath, Eqs.pre_syn_stp)

    nrn = NeuronGroup(N=P.p['N_network'], model= eq_nrn , threshold=P.p['threshold_condition'], reset=Eqs.reset_nrn,   refractory=P.p['refractory_time'],  method='euler')

    # input_nrn = Input()._tbs(P.p)
    input_nrn = PoissonGroup(N=2, rates=[100*Hz, 5*Hz])

    syn = Synapses(input_nrn, nrn, eq_syn, on_pre=eq_syn_pre,)
    syn.connect()

    # rec = {}
    # rec['nrn'] = StateMonitor(nrn, ('u'), record=True)
    # rec['input'] = SpikeMonitor(input_nrn, record=True)
    # rec['syn'] = StateMonitor(syn, ('w_clopath', 'x_trace', ), record=True)

    # set parameters
    #=======================================================================
    for param_key, param_val in P.p.iteritems():
        if param_key in dir(nrn):
            setattr(nrn, param, param_val)
        if param_key in dir(syn):
            setattr(syn, param, param_val)

    # set initial conditions
    #=======================================================================
    for init_key, init_val in P.p.iteritems():
        if init_key in dir(nrn):
            setattr(nrn, init, init_val)
        if param_key in dir(syn):
            setattr(syn, init, init_val)

    run(P.p['run_time'])
    # net = Network(nrn, syn, input_nrn)
    # net.run(P.p['run_time'])

    return rec