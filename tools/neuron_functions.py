from neuron import h
import numpy as np
from scipy import stats
import random as rnd

def create_input_neurons(N, rate, noise, first_spike=0):
    supraspinal_neurons = []
    if type(rate) == np.ndarray:
        for i in range(N):
            if rate[i]==0:
                cell = h.NetStim()
                cell.interval = 1000.0    # Inter-spike interval in ms
                cell.noise = noise
                cell.number = 1e999
                cell.start = 1e999
                supraspinal_neurons.append(cell)
            else:
                cell = h.NetStim()
                cell.interval = 1000.0/ rate[i]   # Inter-spike interval in ms
                cell.noise = noise
                cell.number = 1e999
                cell.start = first_spike
                supraspinal_neurons.append(cell)

    else:
        for _ in range(N):
            cell = h.NetStim()
            cell.interval = 1000.0 / rate  # Inter-spike interval in ms
            cell.noise = noise
            cell.number = 1e999
            cell.start = first_spike
            supraspinal_neurons.append(cell)
    return supraspinal_neurons


def create_spike_recorder_input_neurons(neurons):
    num_neurons = len(neurons)
    spike_times = [h.Vector() for _ in range(num_neurons)]
    spike_detector = [h.NetCon(neurons[i], None) for i in range(num_neurons)]
    for i in range(num_neurons): spike_detector[i].record(spike_times[i])
    return spike_times


def create_spike_recorder_mns(neurons):
    MN_spike_times = [h.Vector() for i in range(len(neurons))]
    MN_spike_detector = []
    for i in range(len(neurons)):
        sp_detector = h.NetCon(neurons[i].soma(0.5)._ref_v, None, sec=neurons[i].soma)
        sp_detector.threshold = -5
        MN_spike_detector.append(sp_detector)
        MN_spike_detector[i].record(MN_spike_times[i])
    return MN_spike_times

def create_exponential_synapses(source,target,W,tau,delay=0,inhibitory=False):
    syn_list=[]
    nc_list=[]

    for itarget in range(len(target)):
        syn_list.append([])
        nc_list.append([])

        for isource in range(len(source)):
            if inhibitory: 
                syn_ = h.Exp2Syn(target[itarget].soma(0.5))
                syn_.tau1 = 1.5
                syn_.tau2 = 2
                syn_.e = -75
            else:
                syn_ = h.ExpSyn(target[itarget].soma(0.5))
                syn_.tau = tau
                
            nc = h.NetCon(source[isource], syn_)
            nc.weight[0] = W[isource,itarget]
            nc_list[-1].append(nc)
            syn_list[-1].append(syn_)
            if type(delay) == np.ndarray: 
                nc.delay = delay[isource,itarget]
            else: 
                nc.delay = delay

    return syn_list,nc_list


def create_inhomongenous_poisson_process(max_fr, tStop, integration_step, offset_x, frequency, offset_y, limit_fr=np.inf): 
    """ Simulate an inhomongenous Poisson process. TODO: full documentation here"""
    
    T = np.arange(0, tStop, integration_step, dtype=float)
    Y = np.maximum(max_fr * np.sin(2 * np.pi * frequency * T - offset_x) + offset_y,0.001)
    Y = np.minimum(Y, limit_fr)
    max_lambda = np.max(Y)
    samples = tStop/integration_step

    total_spikes = stats.poisson.rvs(mu=max_lambda*tStop/1000)
    # create homogenous process
    homogeneous_spikes = np.sort([rnd.uniform(0, samples) for _ in range(total_spikes)])
    # thin process based on inhomogenous rate
    inhomogeneous_spikes = [homogeneous_spikes[i]*integration_step for i in range(total_spikes) if rnd.uniform(0, 1) <= Y[int(homogeneous_spikes[i])]/max_lambda]
    
    return [int(spike) for spike in inhomogeneous_spikes]

def create_inhomogeneous_input_neurons(N, max_rate, tStop=6000, offset_x=0, frequency=.9/2000, offset_y=0, limit_fr=np.inf): 
    input_neurons = []
    for _ in range(N):
        cell = h.VecStim()
        firings = create_inhomongenous_poisson_process(max_rate, tStop=tStop, integration_step=1, offset_x=offset_x, frequency=frequency, offset_y=offset_y, limit_fr=limit_fr)
        vec = h.Vector(firings)
        cell.play(vec)
        input_neurons.append(cell)
    return input_neurons