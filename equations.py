class AdexBonoClopath:

    def __init__(self, ):
        '''
        '''

        self.neuron = '''

            du/dt = int(not_refractory)*I_total/C  + I_reset  + ((t-lastspike)>spike_hold_time2)*(1-int(not_refractory))*I_total/C :volt

            I_total = (I_L + I_syn + I_exp + I_field + I_input - w_adapt + z_after ) : amp

            I_reset = ((t-lastspike)<spike_hold_time2)*((t-lastspike)>spike_hold_time)*(1-int(not_refractory))*((u_reset-u)/t_reset + z_after/C) : volt/second     

            I_L = -g_L*(u - E_L) : amp

            I_exp = g_L*delta_T*exp((u - V_T)/delta_T) : amp

            dV_T/dt = -(V_T-V_Trest)/t_V_T : volt

            dw_adapt/dt = a_adapt*(u-E_L)/t_w_adapt - w_adapt/t_w_adapt : amp

            dz_after/dt = -z_after/t_z_after : amp 

            # synaptic
            #=======================================
            # ampa
            #``````````````
                dg_ampa/dt = -g_ampa/t_ampa : siemens 
                I_ampa = -g_ampa*(u-E_ampa) : amp

            # nmda
            #`````````````````
                dg_nmda/dt = -g_nmda/t_nmda : siemens 
                B =  1/(1 + exp(-0.062*u/mV)/3.57) : 1 
                I_nmda = -g_nmda*B*(u-E_nmda) : amp

                I_syn = I_ampa + I_nmda : amp

            # clopath
            #```````````````````
            # low threshold filtered membrane potential
                du_lowpass1/dt = (u-u_lowpass1)/tau_lowpass1 : volt 

            # high threshold filtered membrane potential
                du_lowpass2/dt = (u-u_lowpass2)/tau_lowpass2 : volt     

            # homeostatic term
                du_homeo/dt = (u-E_L-u_homeo)/tau_homeo : volt       

            # LTP voltage dependence
                LTP_u = (u_lowpass2/mV - theta_low/mV)*int((u_lowpass2/mV - theta_low/mV) > 0)*(u/mV-theta_high/mV)*int((u/mV-theta_high/mV) >0)  : 1

            # LTD voltage dependence
                LTD_u = (u_lowpass1/mV - theta_low/mV)*int((u_lowpass1/mV - theta_low/mV) > 0)  : 1

            # homeostatic depression amplitude
                #``````````````````````````````````
                A_LTD_homeo = A_LTD*(u_homeo**2/v_target) : 1  
            
            # parameters
            #```````````````
            # I_input : amp
            I_field : amp
            # I_syn : amp
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
        #``````````````````````````````
        self.neuron_reset ='''
            z_after = I_after 
            u = int(hold_spike)*u_hold + int(1-hold_spike)*(u_reset + dt*I_after/C)
            V_T = V_Tmax 
            w_adapt += b_adapt    
        '''
        # ampa synapses
        #'````````````````````````````````
        self.syn_ampa = '''
            # dg_ampa/dt = -g_ampa/t_ampa : siemens 
            # I_ampa = -g_ampa*(u_post-E_ampa) : amp
            
            w_ampa :1
            # g_max_ampa : siemens
            # t_ampa : second
            # E_ampa : volt
            '''

        self.syn_ampa_pre = '''
            g_ampa += int(update_ampa_online)*w_clopath*g_max_ampa*A + int(1-update_ampa_online)*w_ampa*g_max_ampa*A 
            '''

        self.syn_nmda_pre = '''
            g_nmda += int(update_nmda_online)*w_clopath*g_max_nmda*A + int(1-update_nmda_online)*w_nmda*g_max_nmda*A 
            '''
        
        # nmda synapses
        #''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
        self.syn_nmda = '''

            # dg_nmda/dt = -g_nmda/t_nmda : siemens 
            # B =  1/(1 + exp(-0.062*u_post/mV)/3.57) : 1 
            # I_nmda = -g_nmda*B*(u_post-E_nmda) : amp
        
            w_nmda :1
            '''

        # short term plasticity
        #''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
        self.syn_stp = '''

            dF/dt = (1-F)/t_F : 1 
            dD1/dt = (1-D1)/t_D1 : 1 
            dD2/dt = (1-D2)/t_D2 : 1 
            dD3/dt = (1-D3)/t_D3 : 1 
            A = F*D1*D2*D3 : 1
            '''

        self.syn_stp_pre = '''
            F += f 
            D1 *= d1
            D2 *= d2
            D3 *= d3 
            '''

        # clopath plasticity rule
        #''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
        self.syn_clopath = '''      
        
            # lowpass presynaptic variable
                dx_trace/dt = -x_trace/tau_x : 1                          

            # clopath rule for potentiation (depression implented with on_pre)

            # dw_clopath/dt = int(w_clopath<w_max_clopath)*A_LTP*x_trace*LTP_u_post :1
            
            dw_clopath/dt = saturated*(w_max_clopath-w_clopath)/dt + (1-saturated)*A_LTP*x_trace*LTP_u_post : 1

            saturated = int((w_clopath+A_LTP*x_trace*LTP_u_post*dt)>w_max_clopath) : 1 # indicates that next weight update brings synapse to saturation

            w_test:1
            '''

        self.syn_clopath_pre = '''

            w_minus = A_LTD_homeo_post*LTD_u_post
            
            w_clopath = clip(w_clopath-w_minus, 0, w_max_clopath)  # apply LTD

            x_trace += dt*x_reset/tau_x  # update presynaptic trace with each input
            '''

        self.synapse = _add_eq(self.syn_ampa, self.syn_nmda, self.syn_stp, self.syn_clopath)

        self.synapse_pre = _add_eq(self.syn_ampa_pre, self.syn_nmda_pre, self.syn_stp_pre, self.syn_clopath_pre)
    

def _add_eq(*list_args):

        '''
        '''
        equation  = ''
        for eq in list_args:
            equation = equation + '\n' + eq 

        return equation 