import numpy as np
import matplotlib.pyplot as plt
from visualization import plot_probe_history
import os

def test_plot_probe_history_manual():
    # This test is just to ensure it doesn't crash and generates a file
    time_array = np.linspace(0, 10, 100)
    probe_data = {
        "MCU": 25 + 5 * np.sin(time_array),
        "FET": 30 + 10 * (1 - np.exp(-time_array/2))
    }
    
    # We use a mock to avoid blocking plt.show during automated tests if possible
    # But since we are refactoring it anyway, we will just run it.
    # To avoid blocking, we can patch plt.show
    from unittest.mock import patch
    with patch("matplotlib.pyplot.show"):
        plot_probe_history(time_array, probe_data)
    
    assert os.path.exists("probe_history.png")
    print("Test passed: plot_probe_history ran successfully.")

if __name__ == "__main__":
    test_plot_probe_history_manual()
