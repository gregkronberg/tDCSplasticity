'''
'''
from brian2 import *
import numpy as np
import pandas as pd
from scipy import stats

class Param:
    '''
    '''
    def __init__(self,):
        '''
        '''
        self.p_run = {
        'method':euler,


        }

        self.p_syn = {
        # ampa
        #`````````````````````````
        'g_max_ampa' : 8*100*nS,
        't_ampa' : 2*ms,
        'E_ampa' : 0*mV,
        # nmda
        #`````````````````````````
        'g_max_nmda' : 100*nS,
        't_nmda' : 50*ms,
        'E_nmda' : 0*mV,

        # facilation/depression
        #``````````````````````````
        'f' : 5,
        't_F' : 94*ms,
        'd1' : 0.45,
        't_D1' : 540*ms,
        'd2' : 0.12,
        't_D2' : 45*ms,
        'd3' : 0.98,
        't_D3' : 120E3*ms,

        # clopath plasticity parameters
        #```````````````````````````````````````````````````````````````````
        # reference value for homeostatic process mV
        'v_target' : 100*mV*mV,         
        # amplitude for the depression
        'A_LTD'  : 2E-5 ,   
        # amplitude for the potentiation 
        'A_LTP'   : 5*38E-6/ms,   
        # time constant for voltage trace in the potentiation term 
        'tau_lowpass2'  : 5*ms,  
        # time constant for presynaptic trace        
        'tau_x' : 10*ms,   
        # time constant for voltage trace in the depression term    
        'tau_lowpass1'  : 6*ms,   
        # time constant for homeostatic process      
        'tau_homeo' : 1000*ms,
        # low threshold
        'theta_low'  : -60*mV,  
        # high threshold in the potentiation term    
        'theta_high'    : -53*mV,       
        # maximum weight
        'w_max_clopath':2,
            }

        self.p_nrn = {

        'E_L' : -69*mV, # mV
        'g_L' : 40*nS, # nS
        'delta_T' : 2*mV, # mV
        'C' : 281*pF,
        't_noise' : 20*ms,
        't_V_T' : 50*ms,
        'refractory':2*ms,
        'V_Trest':-55*mV,
        'V_Tmax':-30*mV,

        'refractory_time' : 2.2*ms,
        'spike_hold_time':1*ms, # must be at least 0.2 ms less than refractory time
        'spike_hold_time2':1.8*ms, # must be at least 0.2 ms less than refractory time
        'reset' : -55*mV,
        # time constant for resetting voltage after holding spike (should be equal to dt)
        't_reset' : .1*ms,

        'threshold' : 'u > 20*mV',
        'reset': 'u=-70 mV',


        'I_input':0*pA,
        'I_field':0*pA,
        'I_after':400*pA,
        'a_adapt': 4*nS,
        'b_adapt' :0.805*pA,
        't_w_adapt':144*ms,
        't_z_after':40*ms,
        'u_reset':-69*mV,
        }

    def _set_init(self, p_nrn, p_syn):
        '''
        '''
        self.nrn_init={
        'u':p_nrn['E_L'],
        'V_T':p_nrn['V_Trest'],
        'w_adapt':0,
        'z_after':0,
        }

        self.syn_init={
        'F':1,
        'D1':1,
        'D2':1,
        'D3':1,
        'w_clopath':0.5,
        'x_trace':0,
        'u_homeo':0,
        'u_lowpass1':p_nrn['E_L'],
        'u_lowpass2':p_nrn['E_L'],

        }

class Equations:
    '''
    '''
    def __init__(self,):
        '''
        '''

        # voltage dynamics
        #========================================================================
        self.nrn['adex'] = '''

        du/dt = (I_L + I_exp + I_syn + I_field + I_input - w_adapt + z_after )/C : volt (unless refractory)

        I_L = -g_L*(u - E_L) : amp

        I_exp = g_L*delta_T*exp((u - V_T)/delta_T) : amp

        dV_T/dt = -(V_T-V_Trest)/t_V_T : volt

        dw_adapt/dt = a_adapt*(u-E_L)/t_w_adapt - w_adapt/t_w_adapt : amp

        dz_after/dt = -z_after/t_z_after : amp
        
        # parameters
        #```````````````
        I_input : amp
        I_field : amp
        I_syn : amp
        I_after : amp
        C : farad
        g_L : siemens
        delta_T : 1 
        t_V_T : second
        a_adapt : siemens
        t_w_adapt : second
        t_z_after : second
        u_reset : volt
        b_adapt : amp
        V_Tmax = volt

        '''

        self.nrn['adex_reset']='''
        z_after = I_after
        u = u_reset + dt*I_after/C
        w_adapt += b_adapt
        V_T = V_Tmax 
        '''

        # Synapses
        #========================================================================
        
        # ampa
        #````````````````````````````````````
        self.syn['ampa'] = '''

        dg_ampa/dt = -g_ampa/t_ampa : siemens 
        I_ampa = -g_ampa*(u_post-E_ampa) : amp
        
        g_max_ampa : siemens
        t_ampa : second
        '''

        self.syn['ampa_on_pre'] = '''
        g_ampa += w_ampa*g_max_ampa*A
        '''

        self.syn['ampa_update_online']='''
        w_ampa=w_clopath
        '''

        self.syn['ampa_update_offline']='''
        w_ampa:1
        '''
        # nmda
        #```````````````````````````````````````
        self.syn['nmda'] = '''

        dg_nmda/dt = -g_nmda/t_nmda : siemens 
        B =  1/(1 + exp(-0.062*u_post/mV)/3.57) : 1 
        I_nmda = -g_nmda*B*(u_post-E_nmda) : amp
        
        g_max_nmda : siemens
        t_nmda : second
        '''

        self.syn['nmda_on_pre'] = '''
        g_nmda += w_nmda*g_max_nmda*A
        '''

        # short term plasticity
        #````````````````````````````````````````
        self.syn['stp']='''

        dF/dt = (1-F)/t_F : 1 (clock-driven)
        dD1/dt = (1-D1)/t_D1 : 1 (clock-driven)
        dD2/dt = (1-D2)/t_D2 : 1 (clock-driven)
        dD3/dt = (1-D3)/t_D3 : 1 (clock-driven)
        A = F*D1*D2*D3 : 1
        
        t_F : second
        f : 1
        t_D1 : second
        d1 : 1
        t_D2 : second
        d2 : 1
        t_D3 : second
        d3 : 1
        '''

        self.syn['stp_on_pre']='''
        F += f                    # facilitation/depression
        D1 *= d1
        D2 *= d2
        D3 *= d3 
        '''

        # long term plasticity (clopath rule)
        #''''''''''''''''''''''''''''''''''''''''''''''
        
        self.syn['clopath'] = '''
        
        # low threshold filtered membrane potential
        du_lowpass1/dt = (u_post-u_lowpass1)/tau_lowpass1 : volt (clock-driven)   

        # high threshold filtered membrane potential
            du_lowpass2/dt = (u_post-u_lowpass2)/tau_lowpass2 : volt     

        # homeostatic term
            du_homeo/dt = (u_post-E_L-u_homeo)/tau_homeo : volt       
        
        # lowpass presynaptic variable
            dx_trace/dt = -x_trace/tau_x : 1                          

        # clopath rule for potentiation (depression implented with on_pre)
            dw_clopath/dt = A_LTP*x_trace*(u_lowpass2/mV - theta_low/mV)*int(u_lowpass2/mV - theta_low/mV > 0)*(u_post/mV-theta_high/mV)*int(u_post/mV-theta_high/mV >0)  : 1

        # homeostatic depression amplitude
            A_LTD_u = A_LTD*(u_homeo**2/v_target) : 1  

        tau_lowpass1 : second
        tau_lowpass2 : second
        tau_homeo : second
        tau_x : second
        theta_low : volt
        theta_high : volt
        A_LTP:1
        A_LTD:1
        v_target:volt*volt
        x_reset:1
        w_max_clopath:1

        '''

        self.syn['clopath_on_pre'] = '''

        w_minus = A_LTD_u*(u_lowpass1/mV - theta_low/mV)*int(u_lowpass1/mV - theta_low/mV > 0)   # update LTD on eahc input spike

        w_clopath = clip(w_clopath-w_minus, 0, w_max_clopath)  # apply LTD

        x_trace += x_reset/tau_x  # update presynaptic trace with each input
        '''

        self.p_syn = {
        # ampa
        #`````````````````````````
        'g_max_ampa' : 8*100*nS,
        't_ampa' : 2*ms,
        'E_ampa' : 0*mV,
        # nmda
        #`````````````````````````
        'g_max_nmda' : 100*nS,
        't_nmda' : 50*ms,
        'E_nmda' : 0*mV,

        # facilation/depression
        #``````````````````````````
        'f' : 5,
        't_F' : 94*ms,
        'd1' : 0.45,
        't_D1' : 540*ms,
        'd2' : 0.12,
        't_D2' : 45*ms,
        'd3' : 0.98,
        't_D3' : 120E3*ms,

        # clopath plasticity parameters
        #```````````````````````````````````````````````````````````````````
        # reference value for homeostatic process mV
        'v_target' : 100*mV*mV,         
        # amplitude for the depression
        'A_LTD'  : 2E-5 ,   
        # amplitude for the potentiation 
        'A_LTP'   : 5*38E-6/ms,   
        # time constant for voltage trace in the potentiation term 
        'tau_lowpass2'  : 5*ms,  
        # time constant for presynaptic trace        
        'tau_x' : 10*ms,   
        # time constant for voltage trace in the depression term    
        'tau_lowpass1'  : 6*ms,   
        # time constant for homeostatic process      
        'tau_homeo' : 1000*ms,
        # low threshold
        'theta_low'  : -60*mV,  
        # high threshold in the potentiation term    
        'theta_high'    : -53*mV,       
        # maximum weight
        'w_max':2,
            }

        self.p_nrn = {

        'E_L' : -69*mV, # mV
        'g_L' : 40*nS, # nS
        'delta_T' : 2*mV, # mV
        'C' : 281*pF,
        't_noise' : 20*ms,
        't_V_T' : 50*ms,
        'refractory':2*ms,
        'V_Trest':-55*mV,
        'V_Tmax':-30*mV,

        'refractory_time' : 2.2*ms,
        'spike_hold_time':1*ms, # must be at least 0.2 ms less than refractory time
        'spike_hold_time2':1.8*ms, # must be at least 0.2 ms less than refractory time
        'reset' : -55*mV,
        # time constant for resetting voltage after holding spike (should be equal to dt)
        't_reset' : .1*ms,

        'threshold' : 'u > 20*mV',
        'reset': 'u=-70 mV',


        'I_input':0*pA,
        'I_field':0*pA,
        'I_after':400*pA,
        'a_adapt': 4*nS,
        'b_adapt' :0.805*pA,
        't_w_adapt':144*ms,
        't_z_after':40*ms,
        'u_reset':-69*mV,
        }

class Input:
    '''
    '''
    def __init__(self):
        '''
        '''
        _build_poisson()

    def _build_poisson(self,N=1, rates=100):
        '''
        '''
        self.input_group = PoissonGroup(N=N, rates=rates)

class Nrn:
    '''
    '''
    def __init__(self, ):
        '''
        '''
        self.nrn = _build_adex(N=1)
        self.input_nrn = Input().input_group 
        self.syn = _build_syn(pre=self.input_nrn, post=self.nrn)

    def _build_adex(self, N=1):
        '''
        '''
        self.Eq = Equations()
        self.P = Param()
        self.eq_nrn = self.Eq.nrn['adex']
        self.eq_nrn_reset = self.Eq.nrn['adex_reset']
        self.nrn = NeuronGroup(N, self.eq_nrn, threshold=self.P.p_nrn['threshold'], reset=self.eq_nrn_reset,  )
        for param_key, param_val in self.P.p_nrn.iteritems():
            if param_key in dir(self.nrn):
                setattr(self.nrn, param, param_val)

        return self.nrn

    def _build_syn(self, pre, post, ):
        '''
        '''
        self.eq_syn = self.Eq.syn['ampa'] + self.Eq.syn['nmda'] + self.Eq.syn['stp'] + self.Eq.syn['clopath'] + self.Eq.syn['ampa_update_online']

        self.eq_syn_on_pre = self.Eq.syn['ampa_on_pre'] + self.Eq.syn['nmda_on_pre'] + self.Eq.syn['stp_on_pre'] + self.Eq.syn['clopath_on_pre']

        self.syn = Synapses(pre, post, self.eq_syn, on_pre=self.eq_syn_on_pre)
        for param_key, param_val in self.P.p_syn.iteritems():
            if param_key in dir(self.nrn):
                setattr(self.syn, param, param_val)

        return self.syn

class Run:
    '''
    '''
    def __init__(self,):
        '''
        '''
        self.Nrn=Nrn()
        self.nrn.syn.connect()
        self.rec_nrn={}
        self.rec_nrn['u'] = StateMonitor(self.Nrn.nrn, ('u'), record=True)
        self.rec_syn['w_clopath'] = StateMonitor(self.Nrn.nrn, ('u'), record=True)
        run_time=1*ms
        run(run_time)
        figure()
        plot(rec_nrn['u'].t/ms, rec_nrn['u']/mV)
        figure()
        plot(rec_syn['w_clopath'].t/ms, rec_nrn['w_clopath']/mV)

        show()
class AdexNeuron:
	'''
	'''
	# voltage dynamics for each comparment
    #```````````````````````````````````````````````````````````````````

    self.eqs['v'] = '''

    du/dt = (I_L + I_exp + I_syn + I_field -w)/C : volt (unless refractory)

    I_L = -g_L*(u - E_L) : amp

    I_exp = g_L*delta_T*exp((u - V_T)/delta_T) : amp

    dV_T/dt = -(V_T-V_Trest)/t_V_T : volt

    # ampa and nmda
    #``````````````````````````````````````````````````````````
    dg_nmda/dt = -g_nmda/t_nmda : siemens 
    dg_ampa/dt = -g_ampa/t_ampa : siemens 
    B =  1/(1 + exp(-0.062*u/mV)/3.57) : 1 
    I_nmda = -g_nmda*B*(u-E_nmda) : amp
    I_ampa = -g_ampa*(u-E_ampa) : amp
    I_syn = I_nmda + I_ampa : amp 
    '''

    self.eqs['ampa'] = 
    # synapse equations
    #`````````````````````````````````````````````````````````````````````
    self.eqs_syn = '''
    # # ampa and nmda
    # #``````````````````````````````````````````````````````````
    # dg_nmda/dt = -g_nmda/t_nmda : siemens (clock-driven)
    # dg_ampa/dt = -g_ampa/t_ampa : siemens (clock-driven)
    # B =  1/(1 + exp(-0.062*u_post/mV)/3.57) : 1 
    # I_nmda = -g_nmda*B*(u_post-E_nmda) : amp
    # I_ampa = -g_ampa*(u_post-E_ampa) : amp
    # I_syn_post = I_nmda + I_ampa : amp (summed)

    # facilitation/depression
    #````````````````````````````````````````````````````````
    dF/dt = (1-F)/t_F : 1 (clock-driven)
    dD1/dt = (1-D1)/t_D1 : 1 (clock-driven)
    dD2/dt = (1-D2)/t_D2 : 1 (clock-driven)
    dD3/dt = (1-D3)/t_D3 : 1 (clock-driven)
    A = F*D1*D2*D3 : 1

    # parameters
    #`````````````````````````````````````````````````````````
    w : 1
    # g_max_nmda : siemens
    # g_max_ampa : siemens
    # t_F : second
    # f : 1
    # t_D1 : second
    # d1 : 1
    # t_D2 : second
    # d2 : 1
    # t_D3 : second
    # d3 : 1

    # clopath learning rule
    #`````````````````````````````````````````````````````````
    # low threshold filtered membrane potential
    du_lowpass1/dt = (u_post-u_lowpass1)/tau_lowpass1 : volt (clock-driven)   

    # high threshold filtered membrane potential
        du_lowpass2/dt = (u_post-u_lowpass2)/tau_lowpass2 : volt     

    # homeostatic term
        du_homeo/dt = (u_post-E_L-u_homeo)/tau_homeo : volt       
    
    # lowpass presynaptic variable
        dx_trace/dt = -x_trace/tau_x : 1                          

    # clopath rule for potentiation (depression implented with on_pre)
        dw_clopath/dt = A_LTP*x_trace*(u_lowpass2/mV - theta_low/mV)*int(u_lowpass2/mV - theta_low/mV > 0)*(u_post/mV-theta_high/mV)*int(u_post/mV-theta_high/mV >0)  : 1

    # homeostatic depression amplitude
        A_LTD_u = A_LTD*(u_homeo**2/v_target) : 1   
    '''
    
    self.eqs = """
        dvm/dt = (gL*(EL - vm) + gL*DeltaT*exp((vm - VT)/DeltaT) + I - w)/C : volt
        dw/dt = (a*(vm - EL) - w)/tauw : amp
        I : amp
"""