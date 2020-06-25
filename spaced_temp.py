'''
'''
import numpy as np
import matplotlib.pyplot as plt



def population_growth_model(**kwargs):
    T = np.arange(10)
    a=2
    b=0.
    c=0
    d=.8
    e=-0.2
    f=1


    activity=1

    plt.figure()
    for activity in [0,1]:
        for b in np.arange(0,0.5, 0.1):
            W=np.zeros(T.shape)
            W[0]=1
            for i, t in enumerate(T):
                if i>0:
                    r = d + e*W[i-1] + f*activity
                    K = a + b*W[i-1] + c*activity
                    dW = r*W[i-1]*(K-W[i-1])/K
                    W[i] += W[i-1]+dW


            plt.plot(W)
    plt.show(block=False)

def statistical_model(**kwargs):
    '''
    '''
    # time vector
    T = np.arange(0, 100)
    # number of synapses
    n_syn = 1000
    # scaling factor for feedback contribution of synaptic strength to input activity
    w_scale=.5/n_syn
    input_scale=1./2
    # parameters for diffusive factors
    tau_rise = 20
    tau_decay= 40
    diffusion_factor= .01/n_syn
    diffusion = np.zeros((T.shape[0]))
    # FIXME option for arbitrary discretization
    # create population of binary synaptic weights
    w = np.zeros((T.shape[0], n_syn ))
    w[0,:] = np.random.randint(2, size=(n_syn))
    # random population of thresholds at each synapse
    threshold = np.zeros(w.shape)
    threshold[0,:] = np.random.random(size=n_syn)
    threshold=np.tile(threshold[0,:], (threshold.shape[0], 1))
    # initialize activity at each synapse
    activity = np.zeros(w.shape)
    # create input vector that specifies the input activity to the network at each time point
    inputs = np.zeros(T.shape)
    # create an input corresponding to each theta burst induction
    inputs[20]=1
    inputs[40]=1
    inputs[60]=1
    inputs[80]=1
    
    d_threshold = np.zeros(w.shape)
    n_crossings = np.zeros(w.shape[0])
    diffusion_rise=np.zeros(w.shape[0])
    diffusion_decay=np.zeros(w.shape[0])
    # step through time
    plt.figure()
    for diff_f in np.arange(0., n_syn/2., n_syn/20.):
        # print(diff_f)
        diffusion_factor = diff_f/n_syn
        def main_loop():
            for t in T:
                # start from t=1
                if t>0:
                    # determine activity (proportional to total synaptic weights)
                    activity[t,:] = w_scale*np.sum(w[t-1,:])*np.random.random(size=n_syn)*inputs[t-1] + input_scale*np.random.random(size=n_syn)*inputs[t-1]
                    # update weights where threshold was crossed
                    w[t,:] = np.where(activity[t-1,:]>threshold[t-1,:], np.ones(w[t-1,:].shape),w[t-1,:] )
                    # update threshold
                    threshold[t,:] = np.clip(threshold[t-1,:]-d_threshold[t-1,:], 0, 1)
                    d_threshold[t, :] = diffusion[t-1]*np.random.random(size=n_syn)/np.sum(w[t-1, :])
                    # update diffusion variable based on threshold crossings
                    #-------------------------------------------------------
                    # get number of crossings in the last time step
                    n_crossings[t] = np.sum(activity[t-1,:]>threshold[t-1,:])
                    # compute rise component of double exponential
                    diffusion_rise[t] = diffusion_rise[t-1]+n_crossings[t-1]-(1/tau_rise)*(diffusion_rise[t-1])
                    # compute decay component of double exponential
                    diffusion_decay[t] = diffusion_decay[t-1]+n_crossings[t-1]-(1/tau_decay)*(diffusion_decay[t-1])
                    # update diffusion variable
                    diffusion[t] = diffusion_factor*(diffusion_decay[t-1]-diffusion_rise[t-1])
                    #FIXME
                    # create a double exponential trace that is proportional to the number of synapses that crossed threshold with some characteristic rise and decay time.  thresholds get lower in proportion to this trace.  this is what models spaced learning effects
                    # # update thresholds based on number of crossing fix me (make this a double exponential with rise and decay time constants)
                    # d_threshold = np.sum(np.where(activity>threshold, np.zeros(threshold.shape)))
            return w, activity, inputs, threshold, diffusion, n_crossings
        w, activity, inputs, threshold, diffusion, n_crossings = main_loop()
        w_norm = np.sum(w, axis=1)/np.sum(w, axis=1)[0]
        
        plt.plot(w_norm)
        plt.ylim([1,2.1])
        plt.show(block=False)
    