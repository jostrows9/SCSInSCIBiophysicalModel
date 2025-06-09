import numpy as np

def spike_times(spikes, dt=1):
    """
    This function trasforms a vector 0 and 1 into a vector
    of spikes times to plot a raster. The units depend on the dt (ms).
    """
    spike_timess = dt*np.where(spikes==1)[0]
    return spike_timess


def inter_spikes_intervals(spikes_times):
    """
    This function trasforms a vector of spikes_times into a vector of spikes time
    intervals.
    """
    isi=np.diff(spikes_times)
    return isi


def inter_spikes_intervals_normalize(spikes,period):
    """
    This function trasforms an array of spikes (NMN x T) to a vector of normalized spike_times.
    """
    n_neurons,_=np.shape(spikes)

    spike_times_aux=[ spike_times(spikes[ineuron]) for ineuron in range(n_neurons)]
    inter_spike_interval=np.concatenate([ inter_spikes_intervals(spike_times_aux[ineuron]) for ineuron in range(n_neurons)] )

    inter_spike_interval_norm=inter_spike_interval/period

    return inter_spike_interval_norm


def bin_fr_hz(data, time, bin_size=100): 
    """
    Convert spike data to firing rate (Hz) by binning spikes. 
    """
    data_flat = [f for fs in data for f in fs] 
    bins = [x for x in range(0, time, bin_size)]
    counts, _ = np.histogram(data_flat, bins)
    return counts/100/len(data)*1000


def firing_rate_to_force(spikes, norm=None):
    """
    Spikes to force folllowing 'Models of Recruiment and Rate Coding Organization
    in Motor Unit Pools' (Fluglevant det al 1993).
    Here we assume that the time bin is 1 ms.
    """
    Nmn, total_time = np.shape(spikes)

    F = np.zeros(total_time)

    # Peak twich force
    RP = 100 # range of twich force
    b = np.log(RP)
    P = np.exp(np.array(range(1,Nmn+1))*b/Nmn)
    np.random.shuffle(P) 

    # Contraction times
    TL = 90 # ms
    RT = 3 # range of contraction times
    c = np.log(RP)/np.log(RT) # this is equivalent to log_RT(RP)
    T = TL/(P**(1/c))
    np.random.shuffle(T)

    # Compute mean ISI. First aproximation of the entire simuation
    sp_times = [spike_times(spikes[ineuron]) for ineuron in range(Nmn)]
    mean_isi = [np.mean( inter_spikes_intervals(sp_times[ineuron]) ) for ineuron in range(Nmn)]

    g = np.zeros(Nmn)
    f = []

    for ineuron in range(Nmn):
        t = np.array(range(int(10*T[ineuron])))
        if len(sp_times[ineuron]==0):
            g[ineuron] = 1
        elif T[ineuron]/mean_isi[ineuron] < 0.4:
            g[ineuron] = 1
        else:
            aux = T[ineuron]/mean_isi[ineuron]
            s = 1-np.exp(-2*(aux)**3) #eq 16
            g[ineuron] = s/aux #eq 17
        f.append((g[ineuron]*P[ineuron]*t/T[ineuron])*np.exp(1-t/T[ineuron]))

    for ineuron in range(Nmn):
        spike_index = np.where(spikes[ineuron]==1)[0]
        for i_initial in spike_index:
            i_final = min(i_initial+len(f[ineuron]), total_time)
            F[i_initial:i_final] = F[i_initial:i_final]+np.array(f[ineuron][:i_final-i_initial])

    if norm is not None:
        F = F/norm

    return F




