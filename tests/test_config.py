import config
import pytest

def test_materials_defined():
    # Expect a MATERIALS dictionary to exist
    assert hasattr(config, 'MATERIALS')
    assert "FR-4" in config.MATERIALS
    assert "Aluminum" in config.MATERIALS
    assert "Ceramic" in config.MATERIALS

def test_material_selection():
    """
    Selecting a material must actually change the derived properties.

    FIX M-8: select_material() was a `pass`, so SELECTED_MATERIAL never
    changed and this failed on the first assert. State is restored by the
    autouse fixture in conftest.py.
    """
    initial_alpha = config.alpha
    config.select_material("Aluminum")
    assert config.SELECTED_MATERIAL == "Aluminum"
    assert config.rho == config.MATERIALS["Aluminum"]["rho"]
    assert config.cp == config.MATERIALS["Aluminum"]["cp"]
    assert config.alpha != initial_alpha


def test_material_selection_rejects_unknown():
    """An unknown material must fail loudly, not silently keep the old one."""
    with pytest.raises(ValueError, match="Unknown material"):
        config.select_material("Unobtainium")


def test_stackup_selection():
    """
    Selecting a stackup must change the layer count and the effective k.

    FIX M-8: select_stackup() was also a `pass`, which is why the copper-impact
    test compared a board against itself.
    """
    config.select_stackup(layers=0)
    assert config.BOARD_LAYERS == 0
    k_bare = config.k

    config.select_stackup(layers=4, copper_oz=1)
    assert config.BOARD_LAYERS == 4
    assert config.COPPER_OZ == 1
    assert config.k > k_bare


def test_config_state_is_restored_between_tests():
    """The autouse fixture must undo what the tests above did."""
    assert config.SELECTED_MATERIAL == "FR-4 (Effective Multilayer)"
    assert config.BOARD_LAYERS == 4

def test_stable_dt_calculation():
    # Ensure dt is updated and is stable
    # dt <= dx^2 / (4 * alpha)
    assert config.dt <= 0.9 * (config.dx**2) / (4 * config.alpha)
