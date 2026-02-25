import matplotlib.pyplot as plt
import numpy as np

def plot_probe_history(time_array, probe_data, fig=None, ax=None, show=True):
    """
    Plots the temperature history of virtual probes.
    
    time_array: NumPy array of time steps [s]
    probe_data: Dictionary mapping probe names to lists/arrays of temperatures [C]
    fig: (Optional) matplotlib figure object
    ax: (Optional) matplotlib axes object
    show: (Optional) Whether to call plt.show()
    """
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    else:
        ax.clear()
    
    for name, temperatures in probe_data.items():
        ax.plot(time_array, temperatures, label=name, linewidth=2)
        
    ax.set_title("Virtual Probe Thermal History", fontsize=14)
    ax.set_xlabel("Time [s]", fontsize=12)
    ax.set_ylabel("Temperature [°C]", fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.legend()
    
    if show and fig is not None:
        plt.savefig("probe_history.png")
        plt.show()

def plot_heatmap(u_matrix, nx_mm=100.0, ny_mm=100.0, fig=None, ax=None, title="PCB Thermal Distribution", clear_fig=False):
    """
    Plots the temperature distribution on the PCB.
    nx_mm, ny_mm: Physical dimensions for the plot axes.
    """
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    
    if clear_fig:
        fig.clear()
        ax = fig.add_subplot(111)
    
    im = ax.imshow(u_matrix, cmap='inferno', origin='lower', extent=[0, nx_mm, 0, ny_mm])
    fig.colorbar(im, ax=ax, label='Temperature [°C]')
    ax.set_title(title)
    ax.set_xlabel('X [mm]')
    ax.set_ylabel('Y [mm]')
    
    return im

def update_heatmap(im, u_matrix, ax, title):
    """
    Updates an existing heatmap image.
    """
    im.set_array(u_matrix)
    im.set_clim(vmin=np.min(u_matrix), vmax=np.max(u_matrix))
    ax.set_title(title)
