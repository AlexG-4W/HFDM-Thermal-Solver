import numpy as np
import config
from solver import PCBSolver
import pytest

def test_solver_material_parameter():
    """Verifies that passing material_name to solver works correctly."""
    nx, ny = 50, 50
    Q = np.zeros((ny, nx))
    Q[20:30, 20:30] = 1e7
    
    # Run with FR-4 using parameter
    solver_fr4 = PCBSolver(nx, ny, Q, material_name="FR-4")
    for _ in range(100):
        solver_fr4.step()
    max_temp_fr4 = np.max(solver_fr4.u)
    
    # Run with Aluminum using parameter
    solver_alu = PCBSolver(nx, ny, Q, material_name="Aluminum")
    for _ in range(100):
        solver_alu.step()
    max_temp_alu = np.max(solver_alu.u)
    
    assert max_temp_alu < max_temp_fr4
    print(f"Peak Temp FR-4: {max_temp_fr4:.2f}, Aluminum: {max_temp_alu:.2f}")

def test_solver_stability_check():
    """Ensure the solver uses the latest dt from config."""
    config.select_material("Aluminum")
    solver = PCBSolver(10, 10, np.zeros((10, 10)))
    # If the solver didn't use the new dt, it might be unstable or incorrect
    # Here we just check that it runs without NaN
    for _ in range(10):
        u = solver.step()
        assert not np.isnan(u).any()

def test_solver_heatsink_impact():
    """Verifies that a localized h_matrix reduces peak temperature."""
    nx, ny = 50, 50
    Q = np.zeros((ny, nx))
    # Add a central heat source
    Q[20:30, 20:30] = 1e7
    
    # Baseline: Global h = 5.0
    solver_base = PCBSolver(nx, ny, Q, material_name="FR-4")
    for _ in range(200):
        solver_base.step()
    max_temp_base = np.max(solver_base.u)
    
    # Heatsink: h=100 in the center
    h_matrix = np.full((ny, nx), config.h)
    h_matrix[20:30, 20:30] = 100.0
    
    # PCBSolver needs to be updated to accept h_matrix
    try:
        solver_hs = PCBSolver(nx, ny, Q, material_name="FR-4", h_matrix=h_matrix)
        for _ in range(200):
            solver_hs.step()
        max_temp_hs = np.max(solver_hs.u)
        
        assert max_temp_hs < max_temp_base
        print(f"Peak Temp Base: {max_temp_base:.2f}, Heatsink: {max_temp_hs:.2f}")
    except TypeError as e:
        pytest.fail(f"PCBSolver.__init__ failed to accept h_matrix: {e}")

def test_heterogeneous_k():
    """Verifies that PCBSolver handles a spatially variable K matrix."""
    nx, ny = 10, 10
    Q = np.zeros((ny, nx))
    Q[5, 5] = 1e6
    K = np.full((ny, nx), 0.3)
    K[5, :] = 385.0 # High conductivity row
    
    solver = PCBSolver(nx, ny, Q, K_matrix=K)
    
    # Run transient for a bit
    for _ in range(100):
        solver.step()
    
    # Heat should have spread significantly further along the copper row (x-direction at y=5)
    # than in the FR-4 direction (y-direction at x=5)
    t_along_copper = solver.u[5, 7] # 2 cells away in copper
    t_along_fr4 = solver.u[7, 5]    # 2 cells away in FR-4
    
    assert t_along_copper > t_along_fr4
    print(f"Copper direction: {t_along_copper:.2f}, FR-4 direction: {t_along_fr4:.2f}")
