"""
Shared test fixtures.

`config` carries module-level mutable state, and `select_material()` /
`select_stackup()` now genuinely mutate it (FIX M-8 - they used to be `pass`).
Without a restore fixture, a test that switches material or stackup leaks into
every test that runs after it, which makes the suite order-dependent.
"""
import pytest

import config

# Everything select_material / select_stackup / update_derived_properties touch,
# plus the environment knobs the GUI writes at run time.
_MUTABLE_CONFIG = (
    "k", "rho", "cp", "alpha", "dt",
    "BOARD_LAYERS", "COPPER_OZ", "SELECTED_MATERIAL",
    "T_amb", "h", "dx", "d", "t_final",
)


@pytest.fixture(autouse=True)
def restore_config_state():
    """Snapshots config globals before each test and puts them back after."""
    saved = {name: getattr(config, name) for name in _MUTABLE_CONFIG}
    yield
    for name, value in saved.items():
        setattr(config, name, value)
