import numpy as np
import matplotlib.pyplot as plt

def calculate_p2p_amp(emg_traces: np.array, 
                      min_window: np.array=None, 
                      max_window: np.array=None): 
    """
    Calculate the peak-to-peak amplitude (P2P) of EMG responses. 

    Args: 
        - emg_traces: np.array, trace number x time
        - min_window: np.array (Optional), start and stop indices for the minimum value of the response.
        - max_window: np.array (Optional), start and stop indcies for the maximum value of the response.
    """
    if not min_window: 
        min_window = [0, len(emg_traces[0])]
    if not max_window: 
        max_window = [0, len(emg_traces[0])]

    p2p_amps = [np.max(trace[max_window[0]:max_window[1]])
            - np.min(trace[min_window[0]:min_window[1]]) 
            for trace in emg_traces]
    
    return p2p_amps

def compute_isi_distribution(spike_times: np.array, 
                             num_bins: int, 
                             density: bool=True, 
                             plot: bool=False,
                             ax=None):
    """
    Return (and, optionally, plot) ISI distribution. 

    Args: 
        - spike times: np.array, list of sequential spike times 
        - num_bins: int, number of bins 
        - density: bool, normalize the probability density
        - plot: bool, plot histogram
        - ax: matplotlib axes (Optional), if provided, plots the histogram on the axes
    
    Out: 
        - n: array, histogram values
        - bin_edges: array, bin edges
    """
    isi = np.diff(spike_times)
    if (plot) & (ax is not None):
        n, bin_edges, _ = ax.hist(isi, bins=num_bins, density=density)
    else: 
        n, bin_edges, _ = plt.hist(isi, bins=num_bins, density=density)
        if plot: 
            plt.show()

    return n, bin_edges

    

    

