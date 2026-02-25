import config
import pytest

def test_materials_defined():
    # Expect a MATERIALS dictionary to exist
    assert hasattr(config, 'MATERIALS')
    assert "FR-4" in config.MATERIALS
    assert "Aluminum" in config.MATERIALS
    assert "Ceramic" in config.MATERIALS

def test_material_selection():
    # Test that changing the material changes alpha
    # Save original state if needed, but here we just test functionality
    initial_alpha = config.alpha
    
    # This is a bit tricky with module-level constants. 
    # Usually, a function like config.select_material("Aluminum") would be better.
    if hasattr(config, 'select_material'):
        config.select_material("Aluminum")
        assert config.SELECTED_MATERIAL == "Aluminum"
        assert config.alpha != initial_alpha
    else:
        pytest.fail("config.select_material not implemented")

def test_stable_dt_calculation():
    # Ensure dt is updated and is stable
    # dt <= dx^2 / (4 * alpha)
    assert config.dt <= 0.9 * (config.dx**2) / (4 * config.alpha)
