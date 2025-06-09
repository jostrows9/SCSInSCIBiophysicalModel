import numpy as np

def plot_raster_plot(ax, 
                    spike_data, 
                    ylabel: str = '', 
                    marker_color: str = '#4c4d4f'): 
    """
    Plot timeseries spike data as a raster plot. 

    Args: 
        - ax: plot axes to plot on 
        - spike_data: spiking data, neurons x time
        - ylabel: y-axis label 
        - marker_color: marker color
    """
    num_neurons = len(spike_data)
    for i in range(num_neurons): ax.plot(spike_data[i],i*np.ones(len(spike_data[i])), ".",  c=marker_color, markersize=1)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_ylabel(ylabel, fontsize=8)
    ax.yaxis.set_tick_params(labelsize=8)
    ax.xaxis.set_tick_params(labelsize=8)

def plot_time_series_data(ax, 
                    time,
                    time_series_data, 
                    ylabel: str = '', 
                    line_color: str = '#4c4d4f'): 
    """
    Plot data over time. 

    Args: 
        - ax: plot axes to plot on 
        - time: simulation time
        - time_series_data: timeseries data (membrane potential, fr, etc.)
        - ylabel: y-axis label 
        - line_color: line color
    """
    ax.plot(time, time_series_data, color=line_color)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_ylabel(ylabel, fontsize=8)
    ax.yaxis.set_tick_params(labelsize=8)
    ax.yaxis.set_tick_params(labelsize=8)
    