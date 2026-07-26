import pytest
import config


def test_copper_constants_defined():
    """
    Verify that copper material constants are defined in config.

    FIX M-8: this asserted on COPPER_K / COPPER_RHO / COPPER_CP, which have
    never existed - config calls them K_CU / RHO_CU / CP_CU. It also had the
    two scalars swapped (COPPER_K == 390, COPPER_CP == 385); copper's thermal
    conductivity is ~385 W/mK and its specific heat ~390 J/kgK, which is what
    config actually holds. The test now uses the real names and values.
    """
    assert hasattr(config, 'K_CU')
    assert config.K_CU == 385.0          # thermal conductivity [W/mK]
    assert hasattr(config, 'RHO_CU')
    assert config.RHO_CU == 8960.0       # density [kg/m^3]
    assert hasattr(config, 'CP_CU')
    assert config.CP_CU == 390.0         # specific heat [J/kgK]


def test_calculate_k_eff_exists():
    """Verify that calculate_k_eff function exists."""
    assert hasattr(config, 'calculate_k_eff')


def test_calculate_k_eff_2layer():
    """Verify k_eff calculation for a 2-layer board."""
    k_fr4 = config.MATERIALS["FR-4"]["k"]
    # FR-4 (1.6mm) + 2 layers of 1oz copper (35um each). Thicknesses in mm;
    # the ratio is unit-free, so it matches the metre-based implementation.
    expected_k = (config.K_CU * 0.070 + k_fr4 * 1.530) / 1.6
    assert config.calculate_k_eff(layers=2, copper_oz=1, substrate_k=k_fr4) \
        == pytest.approx(expected_k, rel=1e-3)


def test_calculate_k_eff_4layer():
    """Verify k_eff calculation for a 4-layer board."""
    k_fr4 = config.MATERIALS["FR-4"]["k"]
    expected_k = (config.K_CU * 0.140 + k_fr4 * 1.460) / 1.6
    assert config.calculate_k_eff(layers=4, copper_oz=1, substrate_k=k_fr4) \
        == pytest.approx(expected_k, rel=1e-3)


def test_calculate_k_eff_monotonic_in_copper():
    """More copper can only raise the effective conductivity."""
    k_fr4 = config.MATERIALS["FR-4"]["k"]
    values = [config.calculate_k_eff(layers=n, copper_oz=1, substrate_k=k_fr4)
              for n in (0, 2, 4, 6)]
    assert values == sorted(values)
    assert values[0] == pytest.approx(k_fr4)      # no copper -> bare substrate


def test_calculate_k_eff_saturates_at_solid_copper():
    """A stackup thicker than the board cannot beat solid copper."""
    assert config.calculate_k_eff(layers=100, copper_oz=1) == config.K_CU
