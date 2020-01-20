'''
'''
from brian2 import *
import numpy as np
import pandas as pd
from scipy import stats
import clopath

prefs.codegen.target = 'numpy'

def _run():
    '''
    '''
# FIWME all synaptic equations need to go in neurongroup
# equations
#========================================================================
    defaultclock.dt=0.025*ms
# voltage dynamics
#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    eq_nrn = '''
        du/dt = (I_L + I_exp + I_syn + I_field + I_input - w_adapt + z_after )/C : volt #(unless refractory)

        I_L = -g_L*(u - E_L) : amp

        I_exp = g_L*delta_T*exp((u - V_T)/delta_T) : amp

        dV_T/dt = -(V_T-V_Trest)/t_V_T : volt

        dw_adapt/dt = a_adapt*(u-E_L)/t_w_adapt - w_adapt/t_w_adapt : amp

        dz_after/dt = -z_after/t_z_after : amp

        # parameters
        #```````````````
        # I_input : amp
        I_field : amp
        I_syn : amp
        # I_after : amp
        # C : farad
        # g_L : siemens
        # delta_T : volt 
        # t_V_T : second
        # a_adapt : siemens
        # t_w_adapt : second
        # t_z_after : second
        # u_reset : volt
        # b_adapt : amp
        # V_Tmax : volt
        # V_Trest:volt
        # E_L : volt
    ''' 

    eq_nrn_hold_spike = '''
        du/dt = int(not_refractory)*(I_L + I_exp + I_syn + I_field + I_input - w_adapt + z_after )/C                                          + ((t-lastspike)<spike_hold_time2)*((t-lastspike)>spike_hold_time)*(1-int(not_refractory))*(u_reset-u)/t_reset     + ((t-lastspike)>spike_hold_time2)*(1-int(not_refractory))*(I_L + I_exp + I_syn + I_field + I_input - w_adapt + z_after )/C   : volt 

        I_L = -g_L*(u - E_L) : amp

        I_exp = g_L*delta_T*exp((u - V_T)/delta_T) : amp

        dV_T/dt = -(V_T-V_Trest)/t_V_T : volt

        dw_adapt/dt = a_adapt*(u-E_L)/t_w_adapt - w_adapt/t_w_adapt : amp

        dz_after/dt = -z_after/t_z_after : amp

        # low threshold filtered membrane potential
        du_lowpass_test/dt = (u-u_lowpass_test)/tau_lowpass1 : volt 

        # parameters
        #```````````````
        # I_input : amp
        I_field : amp
        I_syn : amp
        # I_after : amp
        # C : farad
        # g_L : siemens
        # delta_T : volt 
        # t_V_T : second
        # a_adapt : siemens
        # t_w_adapt : second
        # t_z_after : second
        # u_reset : volt
        # b_adapt : amp
        # V_Tmax : volt
        # V_Trest:volt
        # E_L : volt
    '''

# voltage rest
#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    eq_nrn_reset ='''
        z_after = I_after 
        u = u_reset + dt*I_after/C
        V_T = V_Tmax 
        w_adapt += b_adapt    
    '''

    eq_nrn_reset_hold_spike ='''
        z_after = I_after 
        u = u_hold
        V_T = V_Tmax 
        w_adapt += b_adapt    
    '''
# ampa synapses
#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    eq_syn_ampa = '''
        dg_ampa/dt = -g_ampa/t_ampa : siemens 
        I_ampa = -g_ampa*(u_post-E_ampa) : amp
        
        # w_ampa :1
        # g_max_ampa : siemens
        # t_ampa : second
        # E_ampa : volt
    '''

    eq_syn_ampa_on_pre = '''
        
        g_ampa += w_ampa*g_max_ampa*A 
    '''

    eq_syn_ampa_update_online = '''
        
        w_ampa=w_clopath :1
    '''

    eq_syn_ampa_update_offline = '''
        
        w_ampa:1
    '''

# nmda synapses
#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    eq_syn_nmda = '''

        dg_nmda/dt = -g_nmda/t_nmda : siemens 
        B =  1/(1 + exp(-0.062*u_post/mV)/3.57) : 1 
        I_nmda = -g_nmda*B*(u_post-E_nmda) : amp
        
        # g_max_nmda : siemens
        # t_nmda : second
        # w_nmda:1
        # E_nmda:volt
    '''

    eq_syn_nmda_on_pre = '''
        g_nmda += w_nmda*g_max_nmda*A
    '''

    eq_syn_sum_ampa_nmda = '''
    I_syn_post = I_nmda + I_ampa : amp (summed)
    '''

# short term plasticity
#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    eq_syn_stp = '''

        dF/dt = (1-F)/t_F : 1 (clock-driven)
        dD1/dt = (1-D1)/t_D1 : 1 (clock-driven)
        dD2/dt = (1-D2)/t_D2 : 1 (clock-driven)
        dD3/dt = (1-D3)/t_D3 : 1 (clock-driven)
        A = F*D1*D2*D3 : 1
        
        # t_F : second
        # f : 1
        # t_D1 : second
        # d1 : 1
        # t_D2 : second
        # d2 : 1
        # t_D3 : second
        # d3 : 1
    '''

    eq_syn_stp_on_pre = '''
        F += f 
        D1 *= d1
        D2 *= d2
        D3 *= d3 
    '''

# clopath plasticity rule
#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    eq_syn_clopath = '''
        
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

        # tau_lowpass1 : second
        # tau_lowpass2 : second
        # tau_homeo : second
        # tau_x : second
        # theta_low : volt
        # theta_high : volt
        # A_LTP: hertz
        # A_LTD:1
        # v_target:volt*volt
        # x_reset:1
        # w_max_clopath:1

    '''

    eq_syn_clopath_on_pre = '''

        w_minus = A_LTD_u*(u_lowpass1/mV - theta_low/mV)*int(u_lowpass1/mV - theta_low/mV > 0)   # update LTD on eahc input spike

        w_clopath = clip(w_clopath-w_minus, 0, w_max_clopath)  # apply LTD

        x_trace += dt*x_reset/tau_x  # update presynaptic trace with each input
    '''

# nrn parameters
#============================================================================
    E_L = -69*mV
    g_L = 40*nS
    delta_T = 2*mV
    C = 281*pF
    t_noise = 20*ms
    t_V_T = 50*ms
    refractory = 2*ms
    V_Trest = -55*mV
    V_Tmax = -30*mV
    reset_condition = 'u=-70*mV'
    threshold_condition = 'u>20*mV'
    refractory_time = 1*ms
    u_hold = 10*mV
    I_input = 0*pA
    # I_field = 0*pA
    I_after = 400*pA
    a_adapt = 4*nS
    b_adapt = 0.805*pA
    t_w_adapt = 144*ms
    t_z_after = 40*ms
    u_reset = -60*mV
    I_syn=0*nA

    refractory_time = 2*ms
    spike_hold_time=1*ms # must be at least 0.2 ms less than refractory time
    spike_hold_time2=refractory_time - 2*defaultclock.dt # must be at least 0.2 ms less than refractory time
    # time constant for resetting voltage after holding spike (should be equal to dt)
    t_reset = defaultclock.dt

# synapse parameters
#============================================================================
    # ampa
    #''''''''''''''''''''

    g_max_ampa = 75*nS
    t_ampa = 2*ms
    E_ampa = 0*mV

    # nmda
    #''''''''''''''''''''''
    g_max_nmda = 75*nS
    t_nmda = 50*ms
    E_nmda = 0*mV

    # short term plasticity
    #'''''''''''''''''''''''
    f = 5.3
    t_F = 94*ms
    d1 = 0.45
    t_D1 = 540*ms
    d2 = 0.12
    t_D2 = 45*ms
    d3 = 0.98
    t_D3 = 120E3*ms

    # clopath
    #'''''''''''''''''''''''''
    v_target = 100*mV*mV
    A_LTD = 2E-5
    A_LTP = 5*38E-6/ms
    tau_lowpass2 = 5*ms
    tau_x = 10*ms
    tau_lowpass1 = 6*ms
    tau_homeo = 1000*ms
    theta_low = -60*mV
    theta_high = -40*mV
    w_max_clopath = 2
    x_reset=1
    w_nmda = 0.5
    # w_ampa = 0.5

# input/stimulation parameters
#============================================================================
    pulses=4
    bursts = 4
    pulse_freq = 100
    burst_freq = 5
    warmup = 10

# network parameters
#========================================================================
    N=3


    # create neuron group
    #========================================================================
    # nrn = NeuronGroup(N, eq_nrn , threshold=threshold_condition, reset=eq_nrn_reset, refractory=refractory_time, method='euler')

    nrn = NeuronGroup(N, eq_nrn_hold_spike , threshold=threshold_condition, reset=eq_nrn_reset_hold_spike,   refractory=refractory_time, method='euler')

    # nrn = NeuronGroup(N, eq_nrn, threshold=threshold_condition, reset='u=u_refractory', refractory=refractory_time, method='euler')

    # poisson input group
    #========================================================================
    # input_nrn = PoissonGroup(N=1, rates=10*Hz)


    # theta burst input
    #========================================================================
    input_times = np.zeros(pulses*bursts)*ms
    indices = np.zeros(pulses*bursts)
    cnt=-1
    for burst in range(bursts):
        for pulse in range(pulses):
            cnt+=1
            time = warmup + 1000*burst/burst_freq + 1000*pulse/pulse_freq
            input_times[cnt] = time*ms
            
    input_nrn = SpikeGeneratorGroup(1, indices, input_times )

    # create synapses
    #========================================================================
    eq_syn = eq_syn_ampa + '\n' + eq_syn_nmda + '\n' + eq_syn_stp + '\n' + eq_syn_clopath + '\n' + eq_syn_ampa_update_offline + '\n' + eq_syn_sum_ampa_nmda 

    print eq_syn

    eq_syn_on_pre = eq_syn_ampa_on_pre + '\n' + eq_syn_nmda_on_pre + '\n'+ eq_syn_stp_on_pre + '\n' + eq_syn_clopath_on_pre

    syn = Synapses(input_nrn, nrn, eq_syn, on_pre=eq_syn_on_pre,)
    syn.connect()

    # recording
    #============================================================================
    rec = StateMonitor(nrn, ('u', 'V_T', 'u_lowpass_test'), record=True)
    rec_input = SpikeMonitor(input_nrn, record=True)
    rec_A = StateMonitor(syn, ('A', 'F', 'D1', 'D2', 'D3', 'g_ampa', 'w_clopath', 'x_trace', 'u_lowpass1', 'u_lowpass2', 'u_homeo'), record=True)

    # set initial conditions
    #========================================================================
    # nrn.u['u>20*mV']=20*mV
    # nrn.u['int(not_refractory)'] = u_refractory
    # nrn.u['t-lastspike==refractory_time']=u_reset

    # nrn.u['u>20*mV']=u_hold
    nrn.I_field[0] = 100*pA
    nrn.I_field[1] = 105*pA
    nrn.I_field[2] = -100*pA
    nrn.u = E_L
    nrn.V_T = V_Trest
    nrn.w_adapt = 0*pA
    nrn.z_after = 0*pA
    syn.F=1
    syn.D1=1
    syn.D2=1
    syn.D3=1
    syn.u_lowpass1=E_L
    syn.u_lowpass2=E_L
    # syn.u_homeo=E_L
    syn.w_clopath=0.5
    syn.w_ampa=0.5


    # run
    #============================================================================
    run_time = 300*ms

    run(run_time)

    # figures
    #============================================================================
    figure()
    plot(rec.t/ms, rec.u.T/mV)
    show(block=False)

    return nrn, input_nrn, syn, rec, rec_input, rec_A