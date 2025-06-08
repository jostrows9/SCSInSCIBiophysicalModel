import sys
sys.path.append('../SCSinSCIBiophysicalModel')
from neuron import h
import importlib as impl
import numpy as np
import pickle


from tools.general_tools import ensure_dir
import cells as cll
import tools.neuron_functions as nf


def run_mn_pool_supraspinal_scs_simulation(scs_amp: float,
                                   scs_freq: int, 
                                   perc_supra_intact: float = 1, 
                                   supra_inhibit: bool = False, 
                                   simulation_duration: int = 4000,
                                   save_data_folder: str = '', 
                                   plot_sim: bool = False):
    """ 
    Run a simulation of SCS in 
    """

    impl.reload(nf)
    h.load_file('stdrun.hoc')
    np.random.seed(672945) # set the seed so the network is always the same

    # Set SCS parameters
    num_scs_total = 60
    num_scs_effective = int(scs_amp*num_scs_total)

    # Set supraspinal parameters
    num_supraspinal_total = 300
    num_supraspinal = int(num_supraspinal_total*perc_supra_intact)
    rate_supraspinal = 60   # Firing rate (in Hz)
    
    # Set MN parameters
    num_mn = 100
    mn_drug = True
    mn_avg_diameter = 36
    
    # Set synaptic parameters
    synaptic_weight_scs = 0.000148
    synaptic_weight_supra = 0.000148
    if supra_inhibit: 
        synaptic_weight_supra = synaptic_weight_supra*5 # inhibitory synapses are 5x stronger
    shape = 1.2
    tau = 2
    
    print(f"Running simulation for {scs_amp} SCS amp, {scs_freq} SCS frequency, {num_supraspinal} supraspinal inputs") 

    # Create a list to hold the neurons and spike recorders
    supraspinal_neurons = []
    scs_neurons = []
    supraspinal_spikes = []
    scs_pulse_times = []
    
    # Create a population of scs pulses and record their spikes
    if num_scs_effective > 0:
        scs_neurons = nf.create_input_neurons(num_scs_effective,scs_freq,0)
        scs_pulse_times = nf.create_spike_recorder_input_neurons(scs_neurons)

    # Create a population of supraspinal pulses and record their spikes
    if num_supraspinal > 0: 
        supraspinal_neurons = nf.create_input_neurons(num_supraspinal, rate_supraspinal, noise=1, first_spike=0)
        supraspinal_spikes = nf.create_spike_recorder_input_neurons(supraspinal_neurons)

    # Create a population of MNs and record their spikes
    mn_L = mn_avg_diameter + np.random.randn(num_mn)*0.1*mn_avg_diameter
    mns = [cll.MotoneuronNoDendrites("WT", drug=mn_drug, L=mn_L[imn]) for imn in range(num_mn)]
    mn_spike_times = nf.create_spike_recorder_mns(mns)

    # Connect a population of supraspinal fibers to MNs
    W_supraspinal = 0
    if num_supraspinal > 0:
        W_supraspinal = np.random.gamma(shape, scale=synaptic_weight_supra/shape, size=[num_supraspinal,num_mn])
        syn_supraspinal, nc_supraspinal = nf.create_exponential_synapses(supraspinal_neurons, mns, W_supraspinal, tau, inhibitory=supra_inhibit)

    # Connect a population of scs pulses to MNs
    W_scs = 0
    if num_scs_effective>0:
        W_scs = np.random.gamma(shape, scale=synaptic_weight_scs/shape, size=[num_scs_effective, num_mn])
        delay = np.random.lognormal(-0.47, 0.37, size=[num_scs_effective, num_mn]) #Greiner 2021
        syn_scs, nc_scs = nf.create_exponential_synapses(scs_neurons, mns, W_scs, tau, delay)

    h.finitialize()
    h.tstop = simulation_duration
    h.run()

    supraspinal_spikes = [np.array(supraspinal_spikes[i])  if len(supraspinal_spikes[i]) > 0 else [] for i in range(num_supraspinal)]
    scs_pulse_times = [np.array(scs_pulse_times[i])  if len(scs_pulse_times[i]) > 0 else [] for i in range(num_scs_effective)]
    mn_spike_times = [np.array(mn_spike_times[i])  if len(mn_spike_times[i]) > 0 else [] for i in range(num_mn)]

    if save_data_folder != '': 
        ensure_dir(save_data_folder)
        data_filename = f"mnNum_{num_mn}_supraspinalNum_{num_supraspinal}_supraspinalFR_{rate_supraspinal}_supraInhibit_{supra_inhibit}_SCSFreq_{scs_freq}_SCSTotal_{num_scs_total}_SCSAmp_{scs_amp}_SynW_{synaptic_weight_scs}.pickle"

        data={}
        data["mn_spikes"] = mn_spike_times
        data["supraspinal_spikes"] = supraspinal_spikes
        data["scs_frequency"] = scs_freq
        data["scs_amp"] = scs_amp
        data["num_scs_total"] = num_scs_total
        data["W_scs"] = W_scs
        data["supraspinal_rate"] = rate_supraspinal
        data["num_supraspinal"] = num_supraspinal
        data["simulation_duration"] = simulation_duration
        data["num_mn"] = num_mn
        data["synaptic_weight_scs"] = synaptic_weight_scs
        data["synaptic_weight_supra"] = synaptic_weight_supra
        data["W_supraspinal"] = W_supraspinal
        data["mn_L"] = mn_L

        f=open(save_data_folder+data_filename,"wb")
        pickle.dump(data,f)
        f.close()

    import pdb; pdb.set_trace()

    
if __name__ == '__main__': 
    scs_amp = 0.5
    scs_freq = 40
    perc_supra_intact = 0.2

    run_mn_pool_supraspinal_scs_simulation(scs_amp, 
                                           scs_freq, 
                                           perc_supra_intact, 
                                           plot_sim=True)
    
