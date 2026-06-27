import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.KIE_sites.analysis.phase8_biomass_burning import (  # noqa: E402
    aggregate_monthly_bands,
    apply_bb_to_wetland_corrected_decomposition,
    fit_annual_harmonic,
    gfed_ch4_to_tg_month,
    saueressig_cantrell_r_band,
)


def test_gfed_ch4_to_tg_month_converts_g_m2_with_area_to_tg():
    ch4_g_m2 = np.array([[[2.0, 3.0], [4.0, 5.0]]])
    grid_area_m2 = np.array([[10.0, 10.0], [20.0, 20.0]])

    result = gfed_ch4_to_tg_month(ch4_g_m2, grid_area_m2)

    assert result.shape == (1, 2, 2)
    assert np.isclose(result.sum(), (20 + 30 + 80 + 100) * 1e-12)


def test_aggregate_monthly_bands_sums_expected_latitude_bands():
    monthly_tg = np.ones((12, 4, 2))
    lat = np.array([75.0, 45.0, 0.0, -45.0])
    bands = {
        "NH_high": {"lat_min": 60.0, "lat_max": 90.0},
        "NH_mid": {"lat_min": 30.0, "lat_max": 60.0},
        "Tropics": {"lat_min": -30.0, "lat_max": 30.0},
        "SH_extra": {"lat_min": -90.0, "lat_max": -30.0},
    }

    result = aggregate_monthly_bands(monthly_tg, lat, bands)

    for name in bands:
        assert result[name].shape == (12,)
        assert np.allclose(result[name], 2.0)


def test_fit_annual_harmonic_recovers_sine_cosine_coefficients():
    months = np.arange(12)
    values = (
        10.0
        + 2.0 * np.sin(2 * np.pi * months / 12)
        - 3.0 * np.cos(2 * np.pi * months / 12)
    )

    result = fit_annual_harmonic(values)

    assert np.isclose(result["Q_mean_Tg_month"], 10.0, atol=1e-10)
    assert np.isclose(result["B_Q_Tg_month"], 2.0, atol=1e-10)
    assert np.isclose(result["C_Q_Tg_month"], -3.0, atol=1e-10)


def test_apply_bb_subtracts_bb_source_from_existing_wetland_sink():
    phase6_decomp = {
        "Z_sink_13c": [1.0, 2.0],
        "Z_sink_dD": [10.0, 20.0],
    }
    bb_src_13c = complex(-0.1, 0.2)
    bb_src_dD = complex(-1.0, 2.0)

    result = apply_bb_to_wetland_corrected_decomposition(
        phase6_decomp, bb_src_13c, bb_src_dD
    )

    assert np.allclose(result["Z_sink_wetland_plus_bb_13c"], [1.1, 1.8])
    assert np.allclose(result["Z_sink_wetland_plus_bb_dD"], [11.0, 18.0])


def test_saueressig_cantrell_r_band_returns_ordered_positive_values():
    low, high, details = saueressig_cantrell_r_band()

    assert 0 < low < high
    assert np.isclose(details["R_saueressig"], low)
    assert np.isclose(details["R_cantrell"], high)
    assert np.isclose(low, 0.02777140809745613)
    assert np.isclose(high, 0.032285919025439105)
