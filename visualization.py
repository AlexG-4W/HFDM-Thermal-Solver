"""
Standalone plotting helpers.

Deliberately Qt-free: tests import this module without a QApplication.

plot_heatmap()/update_heatmap() used to live here. They were the last callers
of fig.clear() - plot_heatmap created a colorbar on every call, so the only way
to stop colorbars accumulating was to tear the figure down, which also
destroyed the AxesImage, the axes reference and every probe artist. The board
canvas now keeps its artists alive and retargets a single colorbar; see
board_view.BoardView. Nothing imported these two any more, so they are gone
rather than left as a second, broken way to draw the same thing.
"""

import matplotlib.pyplot as plt

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

