#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 14: Southern Hemisphere wetland source-region sensitivity.

This diagnostic keeps Phase 6 unchanged and asks how CGO/SPO wetland source
corrections move if SH_extra-only is perturbed by either additive transport
stress tests or mass-conserving source-region response mixtures.
"""

from pathlib import Path
import json

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

try:
    from . import phase6_phasor as phase6
except ImportError:  # pragma: no cover - direct script execution
    import phase6_phasor as phase6


EXPT_DIR = Path(__file__).resolve().parent.parent
PHASE2_JSON = EXPT_DIR / "results" / "phase2_harmonics" / "harmonic_fits.json"
WETLAND_JSON = EXPT_DIR / "data" / "wetland_seasonality.json"
PHASE6_JSON = EXPT_DIR / "results" / "phase6_phasor" / "phasor_results.json"
OUT_DIR = EXPT_DIR / "results" / "phase14_sh_wetland_sensitivity"
FIG_DIR = EXPT_DIR / "figures"
OUT_JSON = OUT_DIR / "sh_wetland_sensitivity_results.json"
OUT_FIG_COMPONENTS = FIG_DIR / "fig23_sh_wetland_source_components.png"
OUT_FIG_ENVELOPE = FIG_DIR / "fig24_sh_wetland_sensitivity_envelope.png"
OUT_FIG_MASS_CONSERVING = FIG_DIR / "fig25_sh_wetland_mass_conserving_sensitivity.png"

SH_SITES = ("CGO", "SPO")
N_MC = 50000
MASS_CONSERVING_NH_HIGH_FRACTIONS = (0.04, 0.06, 0.08, 0.10)
MASS_CONSERVING_TROPICS_FRACTIONS = (0.00, 0.10, 0.20, 0.30, 0.40, 0.50)
MASS_CONSERVING_TROPICS_ONLY_FRACTIONS = MASS_CONSERVING_TROPICS_FRACTIONS
MASS_CONSERVING_NH_HIGH_ONLY_FRACTIONS = (0.00, *MASS_CONSERVING_NH_HIGH_FRACTIONS)

SOURCE_BAND_DD = {
    "NH_high": {"dD_CH4": -374.0, "sigma": 10.0},
    "NH_mid": {"dD_CH4": -324.0, "sigma": 14.0},
    "Tropics": {"dD_CH4": -301.0, "sigma": 15.0},
    "SH_extra": {"dD_CH4": -301.0, "sigma": 15.0},
}

SCENARIO_ORDER = [
    "sh_only",
    "tropics_low",
    "tropics_nominal",
    "tropics_high",
    "nh_low",
    "nh_nominal",
    "nh_high",
    "full_low",
    "full_nominal",
    "full_high",
]


def rotate_phasor_delay(B, C, lag_months):
    """Delay a phasor peak by lag_months in the Phase 2 B+iC convention."""
    phi = 2 * np.pi * lag_months / 12.0
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)
    b_new = np.asarray(B) * cos_phi + np.asarray(C) * sin_phi
    c_new = -np.asarray(B) * sin_phi + np.asarray(C) * cos_phi
    return b_new, c_new


def _component(name, band, scale=1.0, lag_months=0.0):
    return {
        "name": name,
        "band": band,
        "scale": float(scale),
        "lag_months": float(lag_months),
    }


def _sh_component():
    return _component("SH_extra", "SH_extra", 1.0, 0.0)


def _tropics_component(scale):
    return _component("Tropics", "Tropics", scale, 0.0)


def _nh_components(scale):
    return [
        _component("NH_high_transport", "NH_high", scale, 2.8),
        _component("NH_mid_transport", "NH_mid", scale, 2.8),
    ]


def _mass_conserving_scenario_name(f_nh, f_tropics):
    nh_key = int(round(100 * f_nh))
    tropics_key = int(round(100 * f_tropics))
    return f"mc_nh{nh_key:02d}_tr{tropics_key:02d}"


def _tropics_only_scenario_name(f_tropics):
    tropics_key = int(round(100 * f_tropics))
    return f"tr_only_{tropics_key:02d}"


def _nh_high_only_scenario_name(f_nh):
    nh_key = int(round(100 * f_nh))
    return f"nh_high_only_{nh_key:02d}"


def scenario_definitions():
    """Return additive SH source-region stress-test scenarios."""
    tropics = {"low": 0.10, "nominal": 0.25, "high": 0.50}
    nh = {"low": 0.05, "nominal": 0.12, "high": 0.20}
    scenarios = {
        "sh_only": {
            "label": "SH only",
            "group": "baseline",
            "components": [_sh_component()],
        }
    }
    for key, scale in tropics.items():
        scenarios[f"tropics_{key}"] = {
            "label": f"SH + Tropics {key}",
            "group": "tropics",
            "components": [_sh_component(), _tropics_component(scale)],
        }
    for key, scale in nh.items():
        scenarios[f"nh_{key}"] = {
            "label": f"SH + NH {key}",
            "group": "nh",
            "components": [_sh_component(), *_nh_components(scale)],
        }
    for key in ("low", "nominal", "high"):
        scenarios[f"full_{key}"] = {
            "label": f"SH + Tropics + NH {key}",
            "group": "full",
            "components": [
                _sh_component(),
                _tropics_component(tropics[key]),
                *_nh_components(nh[key]),
            ],
        }
    return {name: scenarios[name] for name in SCENARIO_ORDER}


def mass_conserving_scenario_definitions():
    """Return source-region response mixtures with weights summing to one."""
    scenarios = {}
    for f_nh in MASS_CONSERVING_NH_HIGH_FRACTIONS:
        for f_tropics in MASS_CONSERVING_TROPICS_FRACTIONS:
            f_sh = 1.0 - f_nh - f_tropics
            if f_sh < -1e-12:
                continue
            name = _mass_conserving_scenario_name(f_nh, f_tropics)
            scenarios[name] = {
                "label": f"MC NH {f_nh:.2f}, Tropics {f_tropics:.2f}",
                "group": "mass_conserving",
                "weights": {
                    "SH_extra": float(f_sh),
                    "Tropics": float(f_tropics),
                    "NH_high": float(f_nh),
                },
                "components": [
                    _component("SH_extra", "SH_extra", f_sh, 0.0),
                    _component("Tropics", "Tropics", f_tropics, 0.0),
                    _component("NH_high_transport", "NH_high", f_nh, 2.8),
                ],
            }
    return scenarios


def mass_conserving_tropics_only_scenario_definitions():
    """Return mass-conserving SH/Tropics slices with no NH component."""
    scenarios = {}
    for f_tropics in MASS_CONSERVING_TROPICS_ONLY_FRACTIONS:
        f_sh = 1.0 - f_tropics
        name = _tropics_only_scenario_name(f_tropics)
        scenarios[name] = {
            "label": f"Tropics {f_tropics:.2f}",
            "group": "tropics_only",
            "weights": {
                "SH_extra": float(f_sh),
                "Tropics": float(f_tropics),
            },
            "components": [
                _component("SH_extra", "SH_extra", f_sh, 0.0),
                _component("Tropics", "Tropics", f_tropics, 0.0),
            ],
        }
    return scenarios


def mass_conserving_nh_high_only_scenario_definitions():
    """Return mass-conserving SH/NH_high slices with no tropical component."""
    scenarios = {}
    for f_nh in MASS_CONSERVING_NH_HIGH_ONLY_FRACTIONS:
        f_sh = 1.0 - f_nh
        name = _nh_high_only_scenario_name(f_nh)
        scenarios[name] = {
            "label": f"NH_high {f_nh:.2f}",
            "group": "nh_high_only",
            "weights": {
                "SH_extra": float(f_sh),
                "NH_high": float(f_nh),
            },
            "components": [
                _component("SH_extra", "SH_extra", f_sh, 0.0),
                _component("NH_high_transport", "NH_high", f_nh, 2.8),
            ],
        }
    return scenarios


def load_all():
    """Load Phase 2, wetland, and Phase 6 inputs."""
    phase2 = json.loads(PHASE2_JSON.read_text(encoding="utf-8"))
    wetland = json.loads(WETLAND_JSON.read_text(encoding="utf-8"))
    phase6_results = json.loads(PHASE6_JSON.read_text(encoding="utf-8"))
    return phase2, wetland, phase6_results


def _pair(z):
    return [float(np.real(z)), float(np.imag(z))]


def _month_or_none(z):
    if not np.isfinite(abs(z)) or abs(z) <= 0:
        return None
    return float(phase6.phasor_peak_month(np.real(z), np.imag(z)))


def _component_source_phasor(wetland, component, q_total=phase6.Q_TOTAL_TG_MONTH):
    band = component["band"]
    source_dd = SOURCE_BAND_DD[band]
    band_fit = wetland["bands"][band]

    b_mid, c_mid = phase6.convert_wetland_to_phase2_phasor(
        float(band_fit["B_Q_Tg_month"]),
        float(band_fit["C_Q_Tg_month"]),
    )
    b_lag, c_lag = rotate_phasor_delay(b_mid, c_mid, component["lag_months"])
    b_scaled = float(b_lag) * component["scale"]
    c_scaled = float(c_lag) * component["scale"]
    z_frac = complex(b_scaled, c_scaled) / q_total
    z13 = (phase6.D13C_WETLAND - phase6.D13C_ATM) * z_frac
    zD = (source_dd["dD_CH4"] - phase6.DD_ATM) * z_frac

    return {
        "name": component["name"],
        "band": band,
        "scale": float(component["scale"]),
        "lag_months": float(component["lag_months"]),
        "d13C_source": float(phase6.D13C_WETLAND),
        "dD_source": float(source_dd["dD_CH4"]),
        "dD_sigma": float(source_dd["sigma"]),
        "Z_frac": _pair(z_frac),
        "Z_src_13C": _pair(z13),
        "Z_src_dD": _pair(zD),
        "A_src_13C": float(abs(z13)),
        "A_src_dD": float(abs(zD)),
        "peak_month_13C": _month_or_none(z13),
        "peak_month_dD": _month_or_none(zD),
    }


def build_total_source_phasor(wetland, components, q_total=phase6.Q_TOTAL_TG_MONTH):
    """Build total isotope source phasors from fixed source-region components."""
    source_components = [
        _component_source_phasor(wetland, component, q_total=q_total)
        for component in components
    ]
    z13 = sum((complex(*item["Z_src_13C"]) for item in source_components), 0j)
    zD = sum((complex(*item["Z_src_dD"]) for item in source_components), 0j)
    return {
        "source_components": source_components,
        "Z_src_total_13C": _pair(z13),
        "Z_src_total_dD": _pair(zD),
        "A_src_total_13C": float(abs(z13)),
        "A_src_total_dD": float(abs(zD)),
        "source_peak_month_13C": _month_or_none(z13),
        "source_peak_month_dD": _month_or_none(zD),
    }


def _phase_diff_months(z13, zD):
    p13 = phase6.phasor_peak_month(np.real(z13), np.imag(z13))
    pD = phase6.phasor_peak_month(np.real(zD), np.imag(zD))
    diff = abs(p13 - pD)
    return float(min(diff, 12.0 - diff))


def _analyze_site_deterministic(code, site_fit, phase6_site, source):
    z_obs_13c = complex(site_fit["d13C"]["B"], site_fit["d13C"]["C"])
    z_obs_dD = complex(site_fit["dD"]["B"], site_fit["dD"]["C"])
    z_src_13c = complex(*source["Z_src_total_13C"])
    z_src_dD = complex(*source["Z_src_total_dD"])
    z_sink_13c = z_obs_13c - z_src_13c
    z_sink_dD = z_obs_dD - z_src_dD
    r_corr = abs(z_sink_13c) / abs(z_sink_dD) if abs(z_sink_dD) > 0 else np.nan
    alpha = phase6.ratio_to_alpha_13c(r_corr)

    phase6_r = float(phase6_site["R_corrected"])
    phase6_alpha = float(phase6.ratio_to_alpha_13c(phase6_r))
    return {
        **source,
        "Z_sink_13C": _pair(z_sink_13c),
        "Z_sink_dD": _pair(z_sink_dD),
        "A_sink_13C": float(abs(z_sink_13c)),
        "A_sink_dD": float(abs(z_sink_dD)),
        "R_corrected": float(r_corr),
        "alpha_13C_OH": float(alpha),
        "delta_R_vs_phase6": float(r_corr - phase6_r),
        "delta_alpha_vs_phase6": float(alpha - phase6_alpha),
        "sink_phase_diff_months": _phase_diff_months(z_sink_13c, z_sink_dD),
    }


def _draw_observed_phasors(site_fit, n_mc, rng):
    b13 = float(site_fit["d13C"]["B"])
    c13 = float(site_fit["d13C"]["C"])
    bD = float(site_fit["dD"]["B"])
    cD = float(site_fit["dD"]["C"])
    a13 = np.hypot(b13, c13)
    aD = np.hypot(bD, cD)
    sigma_a13 = phase6.ci_to_sigma(*site_fit["d13C"]["amplitude_ci95"])
    sigma_aD = phase6.ci_to_sigma(*site_fit["dD"]["amplitude_ci95"])

    p13 = np.arctan2(c13, b13)
    pD = np.arctan2(cD, bD)
    peak13 = site_fit["d13C"].get(
        "peak_month_ci95",
        [site_fit["d13C"]["peak_month"] - 1.0, site_fit["d13C"]["peak_month"] + 1.0],
    )
    peakD = site_fit["dD"].get(
        "peak_month_ci95",
        [site_fit["dD"]["peak_month"] - 1.0, site_fit["dD"]["peak_month"] + 1.0],
    )
    sigma_p13 = phase6.ci_to_sigma(*peak13) * (2 * np.pi / 12.0)
    sigma_pD = phase6.ci_to_sigma(*peakD) * (2 * np.pi / 12.0)

    a13_draw = np.maximum(rng.normal(a13, sigma_a13, n_mc), 1e-6)
    aD_draw = np.maximum(rng.normal(aD, sigma_aD, n_mc), 1e-6)
    p13_draw = p13 + rng.normal(0.0, sigma_p13, n_mc)
    pD_draw = pD + rng.normal(0.0, sigma_pD, n_mc)
    return (
        a13_draw * np.cos(p13_draw) + 1j * a13_draw * np.sin(p13_draw),
        aD_draw * np.cos(pD_draw) + 1j * aD_draw * np.sin(pD_draw),
    )


def _draw_source_phasors(wetland, components, n_mc, rng):
    q_total = np.maximum(
        rng.normal(phase6.Q_TOTAL_TG_YR, phase6.Q_TOTAL_TG_YR_SIGMA, n_mc),
        300.0,
    ) / 12.0
    z13 = np.zeros(n_mc, dtype=complex)
    zD = np.zeros(n_mc, dtype=complex)
    for component in components:
        band = component["band"]
        band_fit = wetland["bands"][band]
        source_dd = SOURCE_BAND_DD[band]
        bq = float(band_fit["B_Q_Tg_month"]) * (
            1 + rng.normal(0.0, phase6.WETLAND_BC_FRAC_SIGMA, n_mc)
        )
        cq = float(band_fit["C_Q_Tg_month"]) * (
            1 + rng.normal(0.0, phase6.WETLAND_BC_FRAC_SIGMA, n_mc)
        )
        bq_mid, cq_mid = phase6.convert_wetland_to_phase2_phasor(bq, cq)
        bq_lag, cq_lag = rotate_phasor_delay(
            bq_mid,
            cq_mid,
            component["lag_months"],
        )
        z_frac = (
            component["scale"] * bq_lag + 1j * component["scale"] * cq_lag
        ) / q_total
        d13c_w = rng.normal(phase6.D13C_WETLAND, phase6.D13C_WETLAND_SIGMA, n_mc)
        dD_w = rng.normal(source_dd["dD_CH4"], source_dd["sigma"], n_mc)
        z13 = z13 + (d13c_w - phase6.D13C_ATM) * z_frac
        zD = zD + (dD_w - phase6.DD_ATM) * z_frac
    return z13, zD


def _alpha_samples_from_ratio(r_samples, n_mc, rng):
    f_oh = np.clip(rng.normal(phase6.F_OH, phase6.SIGMA_F_OH, n_mc), 0.5, 0.99)
    f_cl = np.clip(rng.normal(phase6.F_CL, phase6.SIGMA_F_CL, n_mc), 0.0, 0.1)
    f_soil = np.clip(
        rng.normal(phase6.F_SOIL, phase6.SIGMA_F_SOIL, n_mc), 0.0, 0.15
    )
    f_strat = 1.0 - f_oh - f_cl - f_soil
    alpha_d_oh = rng.normal(phase6.ALPHA_D_OH, phase6.SIGMA_ALPHA_D_OH, n_mc)
    alpha_13c_cl = rng.normal(
        phase6.ALPHA_13C_CL, phase6.SIGMA_ALPHA_13C_CL, n_mc
    )
    alpha_d_cl = rng.normal(phase6.ALPHA_D_CL, phase6.SIGMA_ALPHA_D_CL, n_mc)
    return phase6.ratio_to_alpha_13c(
        r_samples,
        alpha_d_oh,
        f_oh,
        f_cl,
        f_soil,
        f_strat,
        alpha_13c_cl,
        alpha_d_cl,
    )


def mc_site_scenario(site_fit, wetland, components, n_mc=N_MC, seed=123):
    """Monte Carlo propagation for one site and one fixed source scenario."""
    rng = np.random.default_rng(seed)
    z_obs_13c, z_obs_dD = _draw_observed_phasors(site_fit, n_mc, rng)
    z_src_13c, z_src_dD = _draw_source_phasors(wetland, components, n_mc, rng)
    z_sink_13c = z_obs_13c - z_src_13c
    z_sink_dD = z_obs_dD - z_src_dD
    r_samples = np.divide(
        np.abs(z_sink_13c),
        np.abs(z_sink_dD),
        out=np.full(n_mc, np.nan),
        where=np.abs(z_sink_dD) > 1e-6,
    )
    alpha_samples = _alpha_samples_from_ratio(r_samples, n_mc, rng)
    mask = np.isfinite(r_samples) & np.isfinite(alpha_samples) & (r_samples > 0) & (
        r_samples < 1
    )
    r_valid = r_samples[mask]
    alpha_valid = alpha_samples[mask]
    if len(r_valid) == 0:
        return {
            "R_corrected_median": np.nan,
            "R_corrected_ci95": [np.nan, np.nan],
            "alpha_13C_OH_median": np.nan,
            "alpha_13C_OH_ci95": [np.nan, np.nan],
            "n_valid": 0,
        }
    return {
        "R_corrected_median": float(np.median(r_valid)),
        "R_corrected_ci95": [
            float(np.percentile(r_valid, 2.5)),
            float(np.percentile(r_valid, 97.5)),
        ],
        "alpha_13C_OH_median": float(np.median(alpha_valid)),
        "alpha_13C_OH_ci95": [
            float(np.percentile(alpha_valid, 2.5)),
            float(np.percentile(alpha_valid, 97.5)),
        ],
        "n_valid": int(len(r_valid)),
    }


def analyze_site_scenarios(
    code,
    site_fit,
    wetland,
    phase6_site,
    scenarios=None,
    include_mc=True,
    n_mc=N_MC,
    seed_base=123,
):
    """Analyze all fixed source scenarios for one SH site."""
    if scenarios is None:
        scenarios = scenario_definitions()
    phase6_r = float(phase6_site["R_corrected"])
    phase6_alpha = float(phase6.ratio_to_alpha_13c(phase6_r))
    output = {
        "phase6_reference": {
            "source_band": phase6_site["source_band"],
            "R_corrected": phase6_r,
            "alpha_13C_OH": phase6_alpha,
        },
        "scenarios": {},
    }
    for i, (name, scenario) in enumerate(scenarios.items()):
        source = build_total_source_phasor(wetland, scenario["components"])
        result = _analyze_site_deterministic(code, site_fit, phase6_site, source)
        result["scenario_label"] = scenario["label"]
        result["scenario_group"] = scenario["group"]
        if include_mc:
            result["mc"] = mc_site_scenario(
                site_fit,
                wetland,
                scenario["components"],
                n_mc=n_mc,
                seed=seed_base + 1000 * SH_SITES.index(code) + i,
            )
        output["scenarios"][name] = result
    return output


def _summarize_by_scenario(sites, scenario_names):
    summary = {}
    for name in scenario_names:
        r_values = np.array(
            [sites[code]["scenarios"][name]["R_corrected"] for code in SH_SITES],
            dtype=float,
        )
        alpha_values = np.array(
            [sites[code]["scenarios"][name]["alpha_13C_OH"] for code in SH_SITES],
            dtype=float,
        )
        summary[name] = {
            "mean_R_corrected": float(np.nanmean(r_values)),
            "min_R_corrected": float(np.nanmin(r_values)),
            "max_R_corrected": float(np.nanmax(r_values)),
            "mean_alpha_13C_OH": float(np.nanmean(alpha_values)),
            "min_alpha_13C_OH": float(np.nanmin(alpha_values)),
            "max_alpha_13C_OH": float(np.nanmax(alpha_values)),
        }
    return summary


def _compare_to_phase6(sites):
    comparison = {}
    for code in SH_SITES:
        scenarios = sites[code]["scenarios"]
        phase6_ref = sites[code]["phase6_reference"]
        deltas_r = np.array(
            [item["delta_R_vs_phase6"] for item in scenarios.values()], dtype=float
        )
        deltas_alpha = np.array(
            [item["delta_alpha_vs_phase6"] for item in scenarios.values()], dtype=float
        )
        comparison[code] = {
            "phase6_R_corrected": phase6_ref["R_corrected"],
            "phase6_alpha_13C_OH": phase6_ref["alpha_13C_OH"],
            "sh_only_delta_R": scenarios["sh_only"]["delta_R_vs_phase6"],
            "sh_only_delta_alpha": scenarios["sh_only"]["delta_alpha_vs_phase6"],
            "scenario_delta_R_range": [
                float(np.nanmin(deltas_r)),
                float(np.nanmax(deltas_r)),
            ],
            "scenario_delta_alpha_range": [
                float(np.nanmin(deltas_alpha)),
                float(np.nanmax(deltas_alpha)),
            ],
        }
    return comparison


def run_analysis(include_mc=True, n_mc=N_MC):
    """Run Phase 14 analysis and return JSON-serializable results."""
    fits, wetland, phase6_results = load_all()
    scenarios = scenario_definitions()
    mass_conserving_scenarios = mass_conserving_scenario_definitions()
    tropics_only_scenarios = mass_conserving_tropics_only_scenario_definitions()
    nh_high_only_scenarios = mass_conserving_nh_high_only_scenario_definitions()
    sites = {
        code: analyze_site_scenarios(
            code,
            fits[code],
            wetland,
            phase6_results["sites"][code],
            scenarios=scenarios,
            include_mc=include_mc,
            n_mc=n_mc,
        )
        for code in SH_SITES
    }
    mass_conserving_sites = {
        code: analyze_site_scenarios(
            code,
            fits[code],
            wetland,
            phase6_results["sites"][code],
            scenarios=mass_conserving_scenarios,
            include_mc=include_mc,
            n_mc=n_mc,
            seed_base=100000,
        )
        for code in SH_SITES
    }
    tropics_only_sites = {
        code: analyze_site_scenarios(
            code,
            fits[code],
            wetland,
            phase6_results["sites"][code],
            scenarios=tropics_only_scenarios,
            include_mc=include_mc,
            n_mc=n_mc,
            seed_base=200000,
        )
        for code in SH_SITES
    }
    nh_high_only_sites = {
        code: analyze_site_scenarios(
            code,
            fits[code],
            wetland,
            phase6_results["sites"][code],
            scenarios=nh_high_only_scenarios,
            include_mc=include_mc,
            n_mc=n_mc,
            seed_base=300000,
        )
        for code in SH_SITES
    }
    return {
        "metadata": {
            "method": "Southern Hemisphere wetland source-region mixture sensitivity",
            "note": "Diagnostic Phase 14 only; Phase 6 SH_extra-only baseline is unchanged.",
            "sites": list(SH_SITES),
            "n_mc": int(n_mc) if include_mc else 0,
        },
        "scenario_definitions": scenarios,
        "source_band_isotopes": {
            band: {
                "d13C_CH4": float(phase6.D13C_WETLAND),
                "d13C_sigma": float(phase6.D13C_WETLAND_SIGMA),
                "dD_CH4": values["dD_CH4"],
                "dD_sigma": values["sigma"],
            }
            for band, values in SOURCE_BAND_DD.items()
        },
        "sites": sites,
        "sh_summary_by_scenario": _summarize_by_scenario(sites, scenarios.keys()),
        "comparison_to_phase6_sh_only": _compare_to_phase6(sites),
        "mass_conserving": {
            "metadata": {
                "method": "Mass-conserving SH/Tropics/NH_high source-region response mixture",
                "note": "Weights sum to one; delayed NH_high is varied over 0.04/0.06/0.08/0.10, and the remaining weight is assigned to SH_extra after the tropical weight is chosen.",
            },
            "nh_high_fractions": list(MASS_CONSERVING_NH_HIGH_FRACTIONS),
            "tropics_fractions": list(MASS_CONSERVING_TROPICS_FRACTIONS),
            "scenario_definitions": mass_conserving_scenarios,
            "sites": mass_conserving_sites,
            "summary_by_scenario": _summarize_by_scenario(
                mass_conserving_sites,
                mass_conserving_scenarios.keys(),
            ),
            "tropics_only": {
                "fractions": list(MASS_CONSERVING_TROPICS_ONLY_FRACTIONS),
                "scenario_definitions": tropics_only_scenarios,
                "sites": tropics_only_sites,
                "summary_by_scenario": _summarize_by_scenario(
                    tropics_only_sites,
                    tropics_only_scenarios.keys(),
                ),
            },
            "nh_high_only": {
                "fractions": list(MASS_CONSERVING_NH_HIGH_ONLY_FRACTIONS),
                "scenario_definitions": nh_high_only_scenarios,
                "sites": nh_high_only_sites,
                "summary_by_scenario": _summarize_by_scenario(
                    nh_high_only_sites,
                    nh_high_only_scenarios.keys(),
                ),
            },
        },
    }


def json_safe(value):
    """Recursively convert NumPy and non-finite values for strict JSON."""
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_safe(item) for item in value]
    if isinstance(value, tuple):
        return [json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, (np.floating, np.integer)):
        item = float(value)
        return item if np.isfinite(item) else None
    if isinstance(value, float):
        return value if np.isfinite(value) else None
    return value


def _plot_vector(ax, z_pair, color, label, linestyle="-"):
    z = complex(*z_pair)
    ax.annotate(
        "",
        xy=(z.real, z.imag),
        xytext=(0, 0),
        arrowprops=dict(
            arrowstyle="-|>",
            color=color,
            lw=2,
            linestyle=linestyle,
            mutation_scale=13,
        ),
    )
    ax.plot(z.real, z.imag, "o", color=color, ms=4, label=label)


def plot_source_components(results):
    """Plot full_nominal source component phasors for CGO and SPO."""
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    colors = {
        "SH_extra": "C0",
        "Tropics": "C2",
        "NH_high": "C3",
        "NH_mid": "C4",
        "total": "black",
    }
    for col, code in enumerate(SH_SITES):
        scenario = results["sites"][code]["scenarios"]["full_nominal"]
        for row, (iso, title) in enumerate((("13C", "d13C"), ("dD", "dD"))):
            ax = axes[row, col]
            x_values = [0.0]
            y_values = [0.0]
            for component in scenario["source_components"]:
                z_pair = component[f"Z_src_{iso}"]
                color = colors.get(component["band"], "0.5")
                _plot_vector(
                    ax,
                    z_pair,
                    color,
                    f"{component['name']} x{component['scale']:.2f}",
                )
                x_values.append(z_pair[0])
                y_values.append(z_pair[1])
            total_pair = scenario[f"Z_src_total_{iso}"]
            _plot_vector(ax, total_pair, colors["total"], "total", linestyle="--")
            x_values.append(total_pair[0])
            y_values.append(total_pair[1])

            pad_x = max(0.02, (max(x_values) - min(x_values)) * 0.25)
            pad_y = max(0.02, (max(y_values) - min(y_values)) * 0.25)
            ax.set_xlim(min(x_values) - pad_x, max(x_values) + pad_x)
            ax.set_ylim(min(y_values) - pad_y, max(y_values) + pad_y)
            ax.axhline(0, color="0.8", lw=0.8)
            ax.axvline(0, color="0.8", lw=0.8)
            ax.set_xlabel("B source coeff")
            ax.set_ylabel("C source coeff")
            ax.set_title(f"{code} source phasors ({title})")
            ax.grid(True, alpha=0.2)
            ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT_FIG_COMPONENTS, dpi=220)
    plt.close(fig)


def plot_sensitivity_envelope(results):
    """Plot R and alpha sensitivity by scenario."""
    scenarios = list(results["scenario_definitions"].keys())
    x = np.arange(len(scenarios))
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))

    for code, marker in zip(SH_SITES, ("o", "s")):
        r_values = [
            results["sites"][code]["scenarios"][name]["R_corrected"]
            for name in scenarios
        ]
        alpha_values = [
            results["sites"][code]["scenarios"][name]["alpha_13C_OH"]
            for name in scenarios
        ]
        axes[0].plot(x, r_values, marker=marker, label=code)
        axes[1].plot(x, alpha_values, marker=marker, label=code)

    axes[0].set_ylabel("R_corrected")
    axes[1].set_ylabel("alpha_13C_OH")
    axes[1].axhline(phase6.ALPHA_13C_SAUERESSIG, color="C2", ls="--", lw=1)
    axes[1].axhline(phase6.ALPHA_13C_CANTRELL, color="C3", ls="--", lw=1)
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(scenarios, rotation=35, ha="right")
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8)
    fig.suptitle("Southern Hemisphere wetland source-region sensitivity")
    fig.tight_layout()
    fig.savefig(OUT_FIG_ENVELOPE, dpi=220)
    plt.close(fig)


def plot_mass_conserving_sensitivity(results):
    """Plot mass-conserving SH source-region mixture sensitivity."""
    mass = results["mass_conserving"]
    nh_values = mass["nh_high_fractions"]
    tropics_values = mass["tropics_fractions"]
    summary = mass["summary_by_scenario"]

    alpha_grid = np.full((len(nh_values), len(tropics_values)), np.nan)
    for i, f_nh in enumerate(nh_values):
        for j, f_tropics in enumerate(tropics_values):
            name = _mass_conserving_scenario_name(f_nh, f_tropics)
            alpha_grid[i, j] = summary[name]["mean_alpha_13C_OH"]

    fig, axes = plt.subplots(1, 3, figsize=(14.8, 4.2), constrained_layout=True)

    ax = axes[0]
    image = ax.imshow(alpha_grid, origin="lower", aspect="auto", cmap="viridis")
    ax.set_xticks(np.arange(len(tropics_values)))
    ax.set_xticklabels([f"{value:.2f}" for value in tropics_values])
    ax.set_yticks(np.arange(len(nh_values)))
    ax.set_yticklabels([f"{value:.2f}" for value in nh_values])
    ax.set_xlabel("Tropics response weight")
    ax.set_ylabel("Delayed NH high response weight")
    ax.set_title("Mass-conserving grid")
    midpoint = 0.5 * (float(np.nanmin(alpha_grid)) + float(np.nanmax(alpha_grid)))
    for i in range(len(nh_values)):
        for j in range(len(tropics_values)):
            value = alpha_grid[i, j]
            in_lab_range = phase6.ALPHA_13C_SAUERESSIG <= value <= phase6.ALPHA_13C_CANTRELL
            ax.text(
                j,
                i,
                f"{value:.4f}" + ("*" if in_lab_range else ""),
                ha="center",
                va="center",
                color="white" if value < midpoint else "black",
                fontsize=7,
            )
            if in_lab_range:
                ax.add_patch(
                    Rectangle(
                        (j - 0.5, i - 0.5),
                        1.0,
                        1.0,
                        fill=False,
                        edgecolor="white",
                        linewidth=1.2,
                    )
                )
    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)

    for ax, slice_key, title, xlabel, color, name_fn in [
        (
            axes[1],
            "tropics_only",
            "Tropics only",
            "Tropics response weight",
            "C2",
            _tropics_only_scenario_name,
        ),
        (
            axes[2],
            "nh_high_only",
            "Delayed NH high only",
            "Delayed NH high response weight",
            "C0",
            _nh_high_only_scenario_name,
        ),
    ]:
        slice_results = mass[slice_key]
        fractions = slice_results["fractions"]
        values = [
            slice_results["summary_by_scenario"][name_fn(fraction)]["mean_alpha_13C_OH"]
            for fraction in fractions
        ]
        ax.axhspan(
            phase6.ALPHA_13C_SAUERESSIG,
            phase6.ALPHA_13C_CANTRELL,
            color="#b7e1b4",
            alpha=0.35,
            label="Lab range",
        )
        ax.plot(fractions, values, marker="o", color=color, lw=1.8)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Mean alpha13C_OH")
        ax.set_title(title)
        ax.set_xticks(fractions)
        ax.grid(alpha=0.25)
        ax.legend(frameon=False, fontsize=7)

    fig.suptitle("Mass-conserving Southern Hemisphere source-region sensitivity")
    fig.savefig(OUT_FIG_MASS_CONSERVING, dpi=220)
    plt.close(fig)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    results = run_analysis(include_mc=True, n_mc=N_MC)
    OUT_JSON.write_text(
        json.dumps(json_safe(results), indent=2, allow_nan=False),
        encoding="utf-8",
    )
    plot_source_components(results)
    plot_sensitivity_envelope(results)
    plot_mass_conserving_sensitivity(results)
    print(f"Wrote {OUT_JSON}")
    print(f"Wrote {OUT_FIG_COMPONENTS}")
    print(f"Wrote {OUT_FIG_ENVELOPE}")
    print(f"Wrote {OUT_FIG_MASS_CONSERVING}")


if __name__ == "__main__":
    main()
