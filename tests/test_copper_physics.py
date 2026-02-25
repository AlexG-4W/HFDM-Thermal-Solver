import pytest
import config

def test_copper_constants_defined():
    """Verify that copper material constants are defined in config."""
    assert hasattr(config, 'COPPER_K')
    assert config.COPPER_K == 390.0
    assert hasattr(config, 'COPPER_RHO')
    assert config.COPPER_RHO == 8960.0
    assert hasattr(config, 'COPPER_CP')
    assert config.COPPER_CP == 385.0

def test_calculate_k_eff_exists():
    """Verify that calculate_k_eff function exists."""
    assert hasattr(config, 'calculate_k_eff')

def test_calculate_k_eff_2layer():
    """Verify k_eff calculation for a 2-layer board."""
    # Ensure we use FR-4 substrate
    k_fr4 = config.MATERIALS["FR-4"]["k"]
    # FR-4 (1.6mm) + 2 layers of 1oz copper (35um each)
    expected_k = (390 * 0.070 + k_fr4 * 1.530) / 1.6
    assert config.calculate_k_eff(layers=2, copper_oz=1, substrate_k=k_fr4) == pytest.approx(expected_k, rel=1e-3)

def test_calculate_k_eff_4layer():
    """Verify k_eff calculation for a 4-layer board."""
    k_fr4 = config.MATERIALS["FR-4"]["k"]
    # FR-4 (1.6mm) + 4 layers of 1oz copper (35um each)
    expected_k = (390 * 0.140 + k_fr4 * 1.460) / 1.6
    assert config.calculate_k_eff(layers=4, copper_oz=1, substrate_k=k_fr4) == pytest.approx(expected_k, rel=1e-3)
