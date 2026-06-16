import importlib.util
from pathlib import Path
import json

import numpy as np


ANALYSIS_DIR = Path(__file__).resolve().parents[1] / "analysis"
EXPT_DIR = Path(__file__).resolve().parents[1]


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


phase9 = load_module("phase9_sh_source_sensitivity",
                     ANALYSIS_DIR / "phase9_sh_source_sensitivity.py")
phase6 = load_module("phase6_phasor", ANALYSIS_DIR / "phase6_phasor.py")


# ---------------------------------------------------------------------------
# Transport transfer function
# ---------------------------------------------------------------------------
def test_transport_local_band_is_identity():
    """τ = 0 (local / co-located source) must pass through unchanged."""
    H = phase9.transport_transfer(0.0)
    assert np.isclose(H.real, 1.0)
    assert np.isclose(H.imag, 0.0)


def test_transport_attenuates_and_lags_with_tau():
    """Amplitude must decrease monotonically and phase lag grow with τ."""
    taus = [0.0, 0.5, 1.0, 2.0, 3.0]
    mags = [abs(phase9.transport_transfer(t)) for t in taus]
    lags = [-np.angle(phase9.transport_transfer(t)) for t in taus]
    assert all(mags[i] > mags[i + 1] for i in range(len(mags) - 1))
    assert all(lags[i] < lags[i + 1] for i in range(len(lags) - 1))
    assert mags[0] == 1.0


def test_transport_reproduces_phase6_quoted_numbers():
    """At τ=1.3 yr the phase6 comment quotes |H|≈0.11-0.12 and lag≈2.8 months."""
    H = phase9.transport_transfer(1.3)
    abs_H = abs(H)
    lag_months = (-np.angle(H)) / (2 * np.pi) * 12
    assert 0.10 < abs_H < 0.13
    assert 2.5 < lag_months < 3.1


# ---------------------------------------------------------------------------
# Source-phasor construction
# ---------------------------------------------------------------------------
def test_multiband_reduces_to_single_band():
    """A single undamped band must equal phase6's own source-phasor convention."""
    wet = {
        "SH_extra": {"B_Q_Tg_month": 0.0044, "C_Q_Tg_month": 0.0776,
                     "Q_mean_Tg_month": 0.2446},
    }
    z13, zD = phase9.build_total_source([("SH_extra", 0.0)], wet)

    # Replicate phase6's deterministic source construction directly.
    B_mid, C_mid = phase6.convert_wetland_to_phase2_phasor(0.0044, 0.0776)
    z_frac = complex(B_mid, C_mid) / phase9.Q_TOTAL_TG_MONTH
    gap_dD = phase9.BAND_dD["SH_extra"] - phase9.DD_ATM
    gap_13c = phase9.BAND_d13C["SH_extra"] - phase9.D13C_ATM
    z13_ref = gap_13c * z_frac
    zD_ref = gap_dD * z_frac

    assert np.isclose(z13, z13_ref)
    assert np.isclose(zD, zD_ref)


def test_adding_transported_band_grows_source_amplitude():
    """Scenario D (local + remote) must have a larger source phasor than A (local only)."""
    wet = {
        "SH_extra": {"B_Q_Tg_month": 0.0044, "C_Q_Tg_month": 0.0776},
        "Tropics":  {"B_Q_Tg_month": 0.2729, "C_Q_Tg_month": -0.5848},
        "NH_mid":   {"B_Q_Tg_month": -0.2499, "C_Q_Tg_month": -2.4052},
        "NH_high":  {"B_Q_Tg_month": -0.1985, "C_Q_Tg_month": -1.2436},
    }
    _, zD_local = phase9.build_total_source(phase9.scenario_bands("A_local_only", 1.3), wet)
    _, zD_mix = phase9.build_total_source(phase9.scenario_bands("D_transport_mix", 1.3), wet)
    assert abs(zD_mix) > abs(zD_local)


# ---------------------------------------------------------------------------
# Faithful reproduction of phase6 (regression guard)
# ---------------------------------------------------------------------------
def test_scenario_A_reproduces_phase6_for_SH_sites():
    """Scenario A (SH_extra only) must reproduce the published phase6 R_corrected
    for CGO and SPO, confirming the reused code path is faithful."""
    phase6_json = EXPT_DIR / "results" / "phase6_phasor" / "phasor_results.json"
    phase2_json = EXPT_DIR / "results" / "phase2_harmonics" / "harmonic_fits.json"
    wet_json = EXPT_DIR / "data" / "wetland_seasonality.json"
    if not (phase6_json.exists() and phase2_json.exists() and wet_json.exists()):
        import pytest
        pytest.skip("phase2/phase6 outputs not present")

    with open(phase2_json) as f:
        fits = json.load(f)
    with open(wet_json) as f:
        wet = json.load(f)["bands"]
    with open(phase6_json) as f:
        phase6_res = json.load(f)["sites"]

    for site in ("CGO", "SPO"):
        z13, zD = phase9.build_total_source(
            phase9.scenario_bands("A_local_only", 1.3), wet)
        det = phase9.decompose(fits[site], z13, zD)
        # CGO/SPO use the <30° Douglas value (-301) in both phase6 and here.
        assert np.isclose(det["R_corrected"], phase6_res[site]["R_corrected"],
                          atol=2e-3), (
            f"{site}: phase9 R={det['R_corrected']:.4f} vs "
            f"phase6 R={phase6_res[site]['R_corrected']:.4f}")


def test_bulk_R_ordering():
    """Saueressig bulk R must be below Cantrell bulk R (sanity for figure bands)."""
    assert phase9._bulk_R(phase9.ALPHA_13C_SAUERESSIG) < phase9._bulk_R(phase9.ALPHA_13C_CANTRELL)
