COMMENT
Implementation of the STDP rule by Clopath et al., Nat. Neurosci. 13(3):344-352 2010

B. Torben-Nielsen, Hebrew University
C. Clopath, Center for Theoretical Neuroscience, Columbia U. 
ENDCOMMENT

NEURON {
       POINT_PROCESS STDPSynCC
       RANGE tau, e, i : the parameters dealing with the original synapse
       RANGE A_m, A_p, tau_y, tetam,tetap,tau_r,tau_0 : claudia's model
       RANGE t_last_pre : 0 without presynaptic spike, 1 when a presynaptic spike occured
       RANGE u_m1, u_m2,r,g_update,gbar
       RANGE delay_steps,delay_array_u_m1, delay_array_u_m2, pointless_counter 
       RANGE um2s,um1s
       NONSPECIFIC_CURRENT i
}

DEFINE MAX_ARRAY 2000 : delay_steps cannot be higher than this number

UNITS {
      (nA) = (nanoamp)
      (mV) = (millivolt)
      (uS) = (microsiemens)
}

PARAMETER {
	  : parameters for alpha synapse
	  tau = 0.1 (ms) <1e-9,1e9>
	  e = 0	(mV)
	  gbar = 0.05

	  : STDP model parameters
	  tau_0 = 10 : time constant for filtering membrane potential v (called u in the original)
	  tau_r =15 : time constant for low-pass r
	  tau_y = 114 : time constant for post filter for potentiation
	  A_m = 0.00001 : amplitude for depression
	  A_p = 0.00012: amplitude for potentiation
	  tetam = -64.9 : threshold for depression 
	  tetap = -35 : threshold for potentiation 
	  delay_steps = 50 : avoid interference from the AP, set to: AP_width / DT

	  : to make some local variables accessible
	  um1s
	  um2s

	  : internel parameters
	  t_last_pre = -1
	  g_update
	  pointless_counter = 0
}

ASSIGNED {
	 v (mV)
	 i (nA)
	 delay_array_u_m1[MAX_ARRAY]
	 delay_array_u_m2[MAX_ARRAY]
	 delay_array_v[MAX_ARRAY]
}

STATE {
      g (uS)
      u_m1
      u_m2
      r
}

INITIAL {
	g=0
	u_m1 = v 
	u_m2 = v 
	r = 0
	FROM i=0 TO MAX_ARRAY {
	     delay_array_v[i] = v
	     delay_array_u_m1[i] = 0
	     delay_array_u_m2[i] = 0
	}
}

BREAKPOINT {
	   SOLVE state METHOD cnexp
	   i = g*(v - e)
}

DERIVATIVE state { 
	   LOCAL x,u_sig,u_m1_sig,u_m2_sig,set_loc,retrieve_loc 

	  : compute the sigmas
	  if( (v - tetap) >  0) {
	      u_sig = v - tetap 
	      }
	  else { u_sig = 0 }

	  if( (delay_array_u_m1[retrieve_loc] - tetam) > 0) {
	      u_m1_sig = delay_array_u_m1[retrieve_loc] - tetam :
	      um1s = u_m1_sig
	  }
	  else { 
	       u_m1_sig = 0 
	       um1s = 0
	  }
	  if( (delay_array_u_m2[retrieve_loc] - tetam) > 0 ) {
	      u_m2_sig = delay_array_u_m2[retrieve_loc] - tetam :
	      um2s = u_m2_sig
	  }
	  else { 
	       u_m2_sig = 0 
	       um2s = 0
	  }

	  if(t_last_pre == 1) {
	  	x = 1
	    	t_last_pre = 0.5
	    	:printf("x is one")
	  } else {
	    	x = 0
	  }

	  g' = -g/tau
	  u_m1' = (v-u_m1)/tau_0 
	  u_m2' = (v-u_m2)/tau_y
	  r' = (x-r)/ tau_r
	  g_update = - A_m*x*u_m1_sig + A_p*u_sig*r*u_m2_sig 
	  gbar = gbar + g_update

	  : avoid interference of the AP
	  set_loc = fmod(pointless_counter,delay_steps)
	  retrieve_loc = fmod(pointless_counter+delay_steps+1,delay_steps)
	  delay_array_u_m1[set_loc] = u_m1
	  delay_array_u_m2[set_loc] = u_m2
	  pointless_counter = pointless_counter + 1 
}

COMMENT
modified from original AlphaSynapse: when a presynaptic event is received the conductance is always set to the max conductance. Even with multiple successive events, the conductance is not summend and will not go higher than gbar.
ENDCOMMENT
NET_RECEIVE(weight (uS)) {
		   g = gbar : set the max conductance of the synapse
		   t_last_pre = 1
		   printf("received weight=%f, at t=%f\n", weight,t )
}
