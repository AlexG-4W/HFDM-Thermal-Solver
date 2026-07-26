import numpy as np
import config
from solver import PCBSolver
import pytest


def test_solver_copper_impact():
    """
    Adding copper layers must reduce the peak temperature.

    FIX M-8: the original compared two transient runs of 500 steps each, and
    was broken twice over. select_stackup(layers=0) was a no-op, so the
    "no copper" baseline silently kept config.BOARD_LAYERS = 4 and both solvers
    were the same board - hence the failure `26.805646 < 26.805646`. Even with
    the stackup working, each solver derives its own CFL-limited dt from its
    own k, so 500 steps cover very different spans of physical time and the
    comparison would have been confounded rather than merely degenerate.

    Steady state removes both problems: no dt, no step count, just the
    converged answer for each stackup.
    """
    nx, ny = 50, 50
    Q = np.zeros((ny, nx))
    Q[20:30, 20:30] = 1e7

    config.select_material("FR-4")

    config.select_stackup(layers=0)
    u_raw, _, converged_raw = PCBSolver(nx, ny, Q).solve_steady_state()

    u_copper, _, converged_copper = PCBSolver(
        nx, ny, Q, layers=4, copper_oz=1
    ).solve_steady_state()

    assert converged_raw, "bare-substrate case did not converge"
    assert converged_copper, "4-layer case did not converge"
    assert np.max(u_copper) < np.max(u_raw)


def test_solver_accepts_stackup_parameters():
    """PCBSolver must honour layers / copper_oz rather than ignoring them."""
    nx, ny = 20, 20
    Q = np.zeros((ny, nx))
    bare = PCBSolver(nx, ny, Q, layers=0)
    four = PCBSolver(nx, ny, Q, layers=4, copper_oz=1)
    assert bare.k_eff < four.k_eff
    assert bare.k_eff == pytest.approx(config.K_FR4)


def test_more_copper_spreads_heat_further():
    """Higher k must flatten the profile, not merely lower the peak."""
    nx, ny = 50, 50
    Q = np.zeros((ny, nx))
    Q[20:30, 20:30] = 1e7

    spans = []
    for layers in (0, 4):
        u, _, converged = PCBSolver(
            nx, ny, Q, layers=layers, copper_oz=1
        ).solve_steady_state()
        assert converged
        spans.append(float(np.max(u) - np.min(u)))

    assert spans[1] < spans[0]
