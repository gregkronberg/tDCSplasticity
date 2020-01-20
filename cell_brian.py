# -*- coding: utf-8 -*-
"""
Simulate 4 compartment neuron with Clopath STDP learning rule using Brian2 simulator

Created on Thu Apr 05 14:29:25 2018

@author: Greg Kronberg
"""
import numpy as np
from brian2 import*
prefs.codegen.target = 'numpy'

def _four_compartment():
    '''
    '''

    # Parameters
    #========================================================================
    n=5 # number of neurons
    E_L = -69*mV # mV
    g_L = 40*nS # nS
    delta_T = 2*mV # mV
    C = 281*pF
    t_noise = 20*ms
    t_V_T = 50*ms
    refractory_time = 2.2*ms
    spike_hold_time=1*ms # must be at least 0.2 ms less than refractory time
    spike_hold_time2=1.8*ms # must be at least 0.2 ms less than refractory time
    reset = -55*mV
    # time constant for resetting voltage after holding spike (should be equal to dt)
    t_reset = .1*ms

    # axial conductance between compartments (g_axial_fromcompartment_tocompartment)
    #````````````````````````````````````````
    g_axial_soma_basal = 1250*nS
    g_axial_basal_soma = 110*nS
    g_axial_soma_proximal = 1250*nS
    g_axial_proximal_soma = 110*nS
    g_axial_proximal_distal = 1500*nS
    g_axial_distal_proximal = 225*nS

    # synapses
    #````````````````````````````````````````````````````````````````````````
    # ampa
    #`````````````````````````
    g_max_ampa = 8*100*nS
    t_ampa = 2*ms
    E_ampa = 0*mV
    # nmda
    #`````````````````````````
    g_max_nmda = 100*nS
    t_nmda = 50*ms
    E_nmda = 0*mV

    # facilation/depression
    #``````````````````````````
    f = 5
    t_F = 94*ms
    d1 = 0.45
    t_D1 = 540*ms
    d2 = 0.12
    t_D2 = 45*ms
    d3 = 0.98
    t_D3 = 120E3*ms

    # clopath plasticity parameters
    #```````````````````````````````````````````````````````````````````
    # reference value for homeostatic process mV
    v_target     = 100*mV*mV         
    # amplitude for the depression
    A_LTD      = 2E-5    
    # amplitude for the potentiation 
    A_LTP        = 5*38E-6/ms   
    # time constant for voltage trace in the potentiation term 
    tau_lowpass2      = 5*ms  
    # time constant for presynaptic trace        
    tau_x      = 10*ms   
    # time constant for voltage trace in the depression term    
    tau_lowpass1      = 6*ms   
    # time constant for homeostatic process      
    tau_homeo = 1000*ms
    # low threshold
    theta_low      = -60*mV  
    # high threshold in the potentiation term    
    theta_high    = -53*mV       
    # maximum weight
    w_max=2

    # compartment-specific thresholds, resets, adaptation
    #```````````````````````````````````````````````````````````````````````
    threshold_soma = '-55*mV'
    threshold_prox = '-30*mV'
    threshold_dist = '-35*mV'
    threshold_basal = '-30*mV'
    reset_soma = '-55*mV'
    reset_prox = '-55*mV'
    reset_dist = '-55*mV'
    reset_basal = '-55*mV'
    V_Trest_soma = -55*mV
    V_Trest_proximal = -30*mV
    V_Trest_distal = -30*mV
    V_Trest_basal = -30*mV
    V_Tmax_soma = -30*mV
    V_Tmax_proximal = -20*mV
    V_Tmax_distal = -20*mV
    V_Tmax_basal = -20*mV
    V_hold_soma = 20*mV
    V_hold_proximal = -20*mV
    V_hold_distal = -20*mV
    V_hold_basal = -20*mV

    # connection parameters
    #``````````````````````````````
    p_connect = 1

    # input stimulus
    #``````````````````````````````
    bursts = 1
    pulses = 4
    burst_freq = 5
    pulse_freq = 100
    warmup = 5

    method = 'euler'
    dcs = .2*mV
    # w=1

    # equations
    #=======================================================================
    # set voltage to holding potential after threshold crossing
    #``````````````````````````````````````````````````````````
    eqs_reset = '''
    u = u_hold
    V_T = V_Tmax
    '''

    # voltage dynamics for each comparment
    #````````````````````````````````````````````````````````````````````
    eqs_compartment = '''
    du/dt = int(not_refractory)*(I_L + I_exp + I_axial + I_syn + I_field)/C  

        + ((t-lastspike)<spike_hold_time2)*((t-lastspike)>spike_hold_time)*(1-int(not_refractory))*(reset-u)/t_reset 

        + ((t-lastspike)>spike_hold_time2)*(1-int(not_refractory))*(I_L + I_exp + I_axial + I_syn + I_field)/C  : volt 

    # du/dt = int(not_refractory)*(I_L + I_exp + I_axial + I_syn +I_field)/C  + ((t-lastspike)>spike_hold_time)*(1-int(not_refractory))*(reset-u)/t_reset : volt

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

    # FIXME
    I_axial = I_axial1 + I_axial2 : amp
    I_axial1 : amp
    I_axial2 :amp
    I_field = I_field1 + I_field2 : amp
    I_field1 : amp
    I_field2 :amp
    # I_syn : amp
    I_ext : amp
    V_Trest : volt
    V_Tmax : volt
    u_hold : volt
    # C : farad
    '''

    # connection between compartments (axial conductance)
    #`````````````````````````````````````````````````````````````````````
    eqs_connect1 = '''
    
    field : volt # axial current due to electric field
    g_axial_in : siemens
    g_axial_out : siemens
    I_axial1_post = g_axial_in*clip(u_pre-u_post, 0*volt, 1000*volt) + g_axial_out*clip(u_pre-u_post, -1000*volt, 0*volt) :  amp (summed)

    I_field1_post =  g_axial_in*clip(field, 0*volt, 1000*volt)  + g_axial_out*clip(field, -1000*volt, 0*volt) : amp (summed) 
    '''

    # connection between compartments (axial conductance), if compartment has a second connection
    #`````````````````````````````````````````````````````````````````````
    eqs_connect2 = '''
    
    field : volt # axial current due to electric field
    g_axial_in : siemens
    g_axial_out : siemens
    I_axial2_post = g_axial_in*clip(u_pre-u_post, 0*volt, 1000*volt) + g_axial_out*clip(u_pre-u_post, -1000*volt, 0*volt) :  amp (summed) 

    I_field2_post =  g_axial_in*clip(field, 0*volt, 1000*volt)  + g_axial_out*clip(field, -1000*volt, 0*volt) : amp (summed) 
    '''

    # synapse equations
    #`````````````````````````````````````````````````````````````````````
    eqs_syn = '''
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

    # to be executed on each presynaptic spike
    #``````````````````````````````````````````````````````````````````````
    pre_syn = '''
    g_nmda += w*g_max_nmda*A  # nmda 
    g_ampa += w*g_max_ampa*A  # ampa
    F += f                    # facilitation/depression
    D1 *= d1
    D2 *= d2
    D3 *= d3 

    w_minus = A_LTD_u*(u_lowpass1/mV - theta_low/mV)*int(u_lowpass1/mV - theta_low/mV > 0)   # update LTD on eahc input spike

    w_clopath = clip(w_clopath-w_minus, 0, w_max)  # apply LTD

    x_trace += 1   # update presynaptic trace with each input
    '''

    # to be executed on each presynaptic spike for online weight updates
    #``````````````````````````````````````````````````````````````````````
    pre_syn_online = '''
    g_nmda += w_clopath*g_max_nmda*A  # nmda 
    g_ampa += w_clopath*g_max_ampa*A  # ampa
    F += f                    # facilitation/depression
    D1 *= d1
    D2 *= d2
    D3 *= d3 

    w_minus = A_LTD_u*(u_lowpass1/mV - theta_low/mV)*int(u_lowpass1/mV - theta_low/mV > 0)   # update LTD on eahc input spike

    w_clopath = clip(w_clopath-w_minus, 0, w_max)  # apply LTD

    x_trace += 1   # update presynaptic trace with each input
    '''

    # create compartments
    #======================================================================
    soma = NeuronGroup(n, eqs_compartment, threshold='u>'+threshold_soma, reset=eqs_reset, refractory=refractory_time,method=method)
    proximal = NeuronGroup(n, eqs_compartment, threshold='u>'+threshold_prox, reset=eqs_reset, refractory=refractory_time, method=method)
    distal = NeuronGroup(n, eqs_compartment, threshold='u>'+threshold_dist, reset=eqs_reset, refractory=refractory_time, method=method)
    basal = NeuronGroup(n, eqs_compartment, threshold='u>'+threshold_basal, reset=eqs_reset, refractory=refractory_time, method=method)

    # update initial conditions and parameters for each compartment
    #======================================================================
    soma.u = E_L
    proximal.u = E_L
    distal.u = E_L
    basal.u = E_L
    soma.V_T= V_Trest_soma
    proximal.V_T = V_Trest_proximal
    distal.V_T = V_Trest_distal
    basal.V_T = V_Trest_basal
    soma.V_Trest= V_Trest_soma
    proximal.V_Trest = V_Trest_proximal
    distal.V_Trest = V_Trest_distal
    basal.V_Trest = V_Trest_basal
    soma.V_Tmax= V_Tmax_soma
    proximal.V_Tmax = V_Tmax_proximal
    distal.V_Tmax = V_Tmax_distal
    basal.V_Tmax = V_Tmax_basal
    soma.u_hold= V_hold_soma
    proximal.u_hold = V_hold_proximal
    distal.u_hold = V_hold_distal
    basal.u_hold = V_hold_basal

    # connect compartments
    #======================================================================
    connect_soma_basal = Synapses(soma, basal, eqs_connect1, method=method)
    connect_soma_proximal = Synapses(soma, proximal, eqs_connect2, method=method)
    connect_proximal_distal = Synapses(proximal, distal, eqs_connect1, method=method)
    connect_basal_soma = Synapses(basal, soma, eqs_connect1, method=method)
    connect_proximal_soma = Synapses(proximal, soma, eqs_connect2, method=method)
    connect_distal_proximal = Synapses(distal, proximal, eqs_connect1, method=method)
    for n_i in range(n):
        connect_soma_basal.connect(i=n_i, j=n_i)
        connect_soma_proximal.connect(i=n_i, j=n_i)
        connect_proximal_distal.connect(i=n_i, j=n_i)
        connect_basal_soma.connect(i=n_i, j=n_i)
        connect_proximal_soma.connect(i=n_i, j=n_i)
        connect_distal_proximal.connect(i=n_i, j=n_i)

    # update axial conductances
    #=======================================================================
    connect_soma_basal.g_axial_in = g_axial_soma_basal
    connect_soma_basal.g_axial_out = g_axial_basal_soma
    connect_soma_proximal.g_axial_in = g_axial_soma_proximal
    connect_soma_proximal.g_axial_out = g_axial_proximal_soma
    connect_proximal_distal.g_axial_in = g_axial_proximal_distal
    connect_proximal_distal.g_axial_out = g_axial_distal_proximal
    connect_basal_soma.g_axial_in = g_axial_basal_soma
    connect_basal_soma.g_axial_out = g_axial_soma_basal
    connect_proximal_soma.g_axial_in = g_axial_proximal_soma
    connect_proximal_soma.g_axial_out = g_axial_soma_proximal
    connect_distal_proximal.g_axial_in = g_axial_distal_proximal
    connect_distal_proximal.g_axial_out = g_axial_proximal_distal

    # apply extracellular field
    #=======================================================================
    connect_soma_basal.field = -dcs
    connect_soma_proximal.field = dcs
    connect_proximal_distal.field = dcs
    connect_basal_soma.field = dcs
    connect_proximal_soma.field = -dcs
    connect_distal_proximal.field = -dcs

    # generate input stimuli
    #======================================================================
    input_times = np.zeros(pulses*bursts)*ms
    indices = np.zeros(pulses*bursts)
    cnt=-1
    for burst in range(bursts):
        for pulse in range(pulses):
            cnt+=1
            time = warmup + 1000*burst/burst_freq + 1000*pulse/pulse_freq
            input_times[cnt] = time*ms
            
    input_spikes = SpikeGeneratorGroup(1, indices, input_times )

    input_syn = Synapses(input_spikes, proximal, eqs_syn, on_pre=pre_syn, method=method)

    input_syn.connect(condition= 'j==i')

    # set FD variables to 1
    #`````````````````````````````
    input_syn.w=1
    input_syn.F=1
    input_syn.D1=1
    input_syn.D2=1
    input_syn.D3=1

    # set initial clopath weights
    #```````````````````````````````````
    input_syn.w_clopath = 0.5
    input_syn.x_trace=0
    input_syn.u_homeo=0
    input_syn.u_lowpass1 = E_L
    input_syn.u_lowpass2 = E_L

            
    # connect neurons
    #========================================================================
    #FIXME
    recurrent_syn = Synapses(soma, proximal, eqs_syn, on_pre=pre_syn, method=method)
    recurrent_syn.connect(condition='j!=i')

    # set FD variables to 1
    #`````````````````````````````
    recurrent_syn.w=1
    recurrent_syn.F=1
    recurrent_syn.D1=1
    recurrent_syn.D2=1
    recurrent_syn.D3=1

    # set initial clopath weights
    #```````````````````````````````````
    recurrent_syn.w_clopath = 0.5
    recurrent_syn.x_trace=0
    recurrent_syn.u_homeo=0
    recurrent_syn.u_lowpass1 = E_L
    recurrent_syn.u_lowpass2 = E_L

    # record variables
    #========================================================================
    # FIXME
    rec_soma = StateMonitor(soma, ('u'), record=True)
    rec_proximal = StateMonitor(proximal, ('I_axial','u'), record=True)
    rec_distal = StateMonitor(distal, ('I_axial','u'), record=True)
    rec_basal = StateMonitor(basal, ('u'), record=True)
    rec_w = StateMonitor(input_syn, ('w_clopath'), record=True)
    rec_w = StateMonitor(input_syn, ('w_clopath'), record=True)

    # run
    #=======================================================================
    run_time = warmup + 1000*(bursts-1)/burst_freq + 1000*(pulses+1)/pulse_freq 
    run(run_time*ms)
        
    # plot
    #=======================================================================
    figure()
    plot(rec_distal.t/ms, rec_distal.I_axial.T/pA)
    figure()
    plot(rec_distal.t/ms, rec_distal.u.T/mV)
    plot(rec_distal.t/ms, rec_proximal.u.T/mV)
    plot(rec_distal.t/ms, rec_soma.u.T/mV)
    plot(rec_distal.t/ms, rec_basal.u.T/mV)
    figure()
    plot(rec_distal.t/ms, rec_soma.u.T/mV)
    figure()
    plot(rec_w.t/ms, rec_w.w_clopath.T)
    
    show()

if __name__ == '__main__':
    _test_run()


    