import numpy as np
import config
from solver import PCBSolver
import pytest

def test_solver_copper_impact():
    """Verifies that adding copper layers to PCBSolver reduces peak temperature."""
    nx, ny = 50, 50
    Q = np.zeros((ny, nx))
    # Add a central heat source
    Q[20:30, 20:30] = 1e7
    
    # 1. No copper (Baseline FR-4)
    # Using default from config
    config.select_material("FR-4")
    config.select_stackup(layers=0)
    solver_raw = PCBSolver(nx, ny, Q)
    for _ in range(500):
        solver_raw.step()
    max_temp_raw = np.max(solver_raw.u)
    
    # 2. 4-layer 1oz copper
    # We want to test if PCBSolver accepts layers and copper_oz parameters
    try:
        solver_copper = PCBSolver(nx, ny, Q, layers=4, copper_oz=1)
        for _ in range(500):
            solver_copper.step()
        max_temp_copper = np.max(solver_copper.u)
        
        # Copper should significantly reduce peak temperature
        assert max_temp_copper < max_temp_raw
        print(f"Max Temp Raw: {max_temp_raw:.2f}, Copper: {max_temp_copper:.2f}")
    except TypeError:
        pytest.fail("PCBSolver.__init__ failed to accept layers and copper_oz parameters")
