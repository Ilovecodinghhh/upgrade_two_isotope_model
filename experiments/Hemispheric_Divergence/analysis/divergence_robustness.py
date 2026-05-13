#!/usr/bin/env python3
"""
Phase 2.2: Divergence robustness across model configurations.

Tests whether the NH-FF-positive / SH-Mic-positive pattern persists under:
  1. default (sampled KIE + time-varying τ)
  2. fixed_lifetime (τ = 9.0 yr)
  3. cantrell_only (KIE = 1.0054)
  4. saueressig_only (KIE = 1.0039)
  5. low_Cl (Cl = 0.6%)
  6. high_Cl (Cl = 6.5%)
  7. fast_exchange (τ_ex = 0.8 yr, fixed)
  8. slow_exchange (τ_ex = 1.3 yr, fixed)
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]  # repo root
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from scipy import stats as sp_stats
from scipy.optimize import lsq_linear
import json

from common import (
    ModelConfig, load_data,
    sample_KIE, compute_bulk_KIE, compute_lifetime,
    delta_to_fraction_d13C, delta_to_fraction_dD,
    sample_source_signatures_hemi,
    sample_atm_d13C, sample_atm_dD_hemi,
    pad_to_length, KIE_FIXED,
    SINK_FRACTIONS_NH, SINK_FRACTIONS_SH,
    PT_HEMI, LIFETIME_RATIO_NH, LIFETIME_RATIO_SH,
    TAU_EX_MEAN, TAU_EX_STD,
)

EXP_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = EXP_DIR / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

NI = 400  # Fewer iterations for speed across many configs


def compute_trends(arr, years, start=2007):
    mask = (years >= start) & (years <= years[-1])
    yrs = years[mask]
    sub = arr[mask, :]
    slopes = np.zeros(sub.shape[1])
    for k in range(sub.shape[1]):
        col = sub[:, k]
        if np.any(np.isnan(col)):
            slopes[k] = np.nan
            continue
        slopes[k] = sp_stats.linregress(yrs, col).slope
    return slopes


def run_config(data, cfg, tau_ex_override=None, cl_frac_override=None,
               kie_oh13c_override=None, label="default"):
    """Run 2-box model with hemispheric sigs for one configuration."""
    n = data.n_years
    years = data.model_years
    rng = np.random.default_rng(cfg.seed)

    CH4_NH, CH4_SH = data.CH4_NH, data.CH4_SH
    c13_NH, c13_SH = data.c13_NH, data.c13_SH
    c13_glob = data.c13_global

    tau_global = compute_lifetime(years, cfg.lifetime_mode, cfg.tau_fixed)
    tau_NH = tau_global * LIFETIME_RATIO_NH
    tau_SH = tau_global * LIFETIME_RATIO_SH

    FF_NH = np.zeros((n, NI)); FF_SH = np.zeros((n, NI))
    Mic_NH = np.zeros((n, NI)); Mic_SH = np.zeros((n, NI))

    W_NH = np.diag([100.0, 1.0, 0.5])
    W_SH = np.diag([200.0, 1.0, 0.5])

    for k in range(NI):
        if tau_ex_override is not None:
            tau_ex = tau_ex_override
        else:
            tau_ex = max(0.5, rng.normal(TAU_EX_MEAN, TAU_EX_STD))

        kies = sample_KIE(rng, cfg.kie_mode)
        # Override OH_13C KIE if requested (for Cantrell/Saueressig tests)
        if kie_oh13c_override is not None:
            kies['OH_13C'] = kie_oh13c_override

        # Override Cl fraction if requested
        sf_nh = dict(SINK_FRACTIONS_NH)
        sf_sh = dict(SINK_FRACTIONS_SH)
        if cl_frac_override is not None:
            for sf in [sf_nh, sf_sh]:
                old_cl = sf['Cl']
                delta = cl_frac_override - old_cl
                sf['Cl'] = cl_frac_override
                sf['OH'] -= delta  # Redistribute to OH

        K13_NH, KD_NH = compute_bulk_KIE(kies, sf_nh)
        K13_SH, KD_SH = compute_bulk_KIE(kies, sf_sh)
        a13_NH = 1.0 / K13_NH; aD_NH = 1.0 / KD_NH
        a13_SH = 1.0 / K13_SH; aD_SH = 1.0 / KD_SH

        S_NH = np.zeros(n); S_SH = np.zeros(n)
        for i in range(n):
            M_NH = CH4_NH[i] * PT_HEMI; M_NH1 = CH4_NH[i+1] * PT_HEMI
            M_SH = CH4_SH[i] * PT_HEMI; M_SH1 = CH4_SH[i+1] * PT_HEMI
            S_NH[i] = (M_NH1 - M_NH) + M_NH / tau_NH[i] - (M_SH - M_NH) / tau_ex
            S_SH[i] = (M_SH1 - M_SH) + M_SH / tau_SH[i] - (M_NH - M_SH) / tau_ex

        d13C_glob_MC = sample_atm_d13C(data, k, n)
        dD_NH_MC, dD_SH_MC = sample_atm_dD_hemi(data, k, n)

        nc = min(len(c13_glob), n + 1)
        d13C_off = d13C_glob_MC[:nc] - c13_glob[:nc]
        d13C_NH_MC = c13_NH[:nc] + d13C_off
        d13C_SH_MC = c13_SH[:nc] + d13C_off

        f13_NH_atm = delta_to_fraction_d13C(d13C_NH_MC)
        f13_SH_atm = delta_to_fraction_d13C(d13C_SH_MC)
        fD_NH_atm = delta_to_fraction_dD(dD_NH_MC)
        fD_SH_atm = delta_to_fraction_dD(dD_SH_MC)

        d13C_src_NH = np.zeros(n); d13C_src_SH = np.zeros(n)
        dD_src_NH = np.zeros(n); dD_src_SH = np.zeros(n)
        for j in range(n):
            n13_NH = f13_NH_atm[j] * CH4_NH[j] * PT_HEMI
            n13_NH1 = f13_NH_atm[j+1] * CH4_NH[j+1] * PT_HEMI
            n13_SH = f13_SH_atm[j] * CH4_SH[j] * PT_HEMI
            n13_SH1 = f13_SH_atm[j+1] * CH4_SH[j+1] * PT_HEMI
            d13C_src_NH[j] = (n13_NH1 - n13_NH + n13_NH * a13_NH / tau_NH[j] - (n13_SH - n13_NH) / tau_ex) / S_NH[j]
            d13C_src_SH[j] = (n13_SH1 - n13_SH + n13_SH * a13_SH / tau_SH[j] - (n13_NH - n13_SH) / tau_ex) / S_SH[j]

            nD_NH = fD_NH_atm[j] * CH4_NH[j] * PT_HEMI
            nD_NH1 = fD_NH_atm[j+1] * CH4_NH[j+1] * PT_HEMI
            nD_SH = fD_SH_atm[j] * CH4_SH[j] * PT_HEMI
            nD_SH1 = fD_SH_atm[j+1] * CH4_SH[j+1] * PT_HEMI
            dD_src_NH[j] = (nD_NH1 - nD_NH + nD_NH * aD_NH / tau_NH[j] - (nD_SH - nD_NH) / tau_ex) / S_NH[j]
            dD_src_SH[j] = (nD_SH1 - nD_SH + nD_SH * aD_SH / tau_SH[j] - (nD_NH - nD_SH) / tau_ex) / S_SH[j]

        sigs = sample_source_signatures_hemi(rng, data, k, n)

        for j in range(n):
            f13_bb_nh = delta_to_fraction_d13C(sigs['bb_d13C_NH'][j])
            f13_ff_nh = delta_to_fraction_d13C(sigs['ff_d13C_NH'][j])
            f13_mic_nh = delta_to_fraction_d13C(sigs['mic_d13C_NH'][j])
            fD_bb_nh = delta_to_fraction_dD(sigs['bb_dD_NH'][j])
            fD_ff_nh = delta_to_fraction_dD(sigs['ff_dD_NH'][j])
            fD_mic_nh = delta_to_fraction_dD(sigs['mic_dD_NH'][j])

            A_nh = np.array([[1.0, 1.0, 1.0],
                             [f13_bb_nh, f13_ff_nh, f13_mic_nh],
                             [fD_bb_nh, fD_ff_nh, fD_mic_nh]])
            B_nh = np.array([S_NH[j], S_NH[j]*d13C_src_NH[j], S_NH[j]*dD_src_NH[j]])
            try:
                res = lsq_linear(W_NH @ A_nh, W_NH @ B_nh, bounds=(0, S_NH[j]*1.5))
                FF_NH[j,k] = res.x[1]; Mic_NH[j,k] = res.x[2]
            except:
                FF_NH[j,k] = Mic_NH[j,k] = np.nan

            f13_bb_sh = delta_to_fraction_d13C(sigs['bb_d13C_SH'][j])
            f13_ff_sh = delta_to_fraction_d13C(sigs['ff_d13C_SH'][j])
            f13_mic_sh = delta_to_fraction_d13C(sigs['mic_d13C_SH'][j])
            fD_bb_sh = delta_to_fraction_dD(sigs['bb_dD_SH'][j])
            fD_ff_sh = delta_to_fraction_dD(sigs['ff_dD_SH'][j])
            fD_mic_sh = delta_to_fraction_dD(sigs['mic_dD_SH'][j])

            A_sh = np.array([[1.0, 1.0, 1.0],
                             [f13_bb_sh, f13_ff_sh, f13_mic_sh],
                             [fD_bb_sh, fD_ff_sh, fD_mic_sh]])
            B_sh = np.array([S_SH[j], S_SH[j]*d13C_src_SH[j], S_SH[j]*dD_src_SH[j]])
            try:
                res = lsq_linear(W_SH @ A_sh, W_SH @ B_sh, bounds=(0, S_SH[j]*1.5))
                FF_SH[j,k] = res.x[1]; Mic_SH[j,k] = res.x[2]
            except:
                FF_SH[j,k] = Mic_SH[j,k] = np.nan

    return years, FF_NH, FF_SH, Mic_NH, Mic_SH


def main():
    print("=" * 70)
    print("PHASE 2.2: Divergence Robustness")
    print("=" * 70)

    data = load_data(ROOT, two_box=True)

    configs = [
        # (name, cfg, tau_ex_override, cl_frac_override, kie_oh13c_override)
        ("default",       ModelConfig(NI, "sampled", "varying", 9.0, 42), None, None, None),
        ("fixed_tau9",    ModelConfig(NI, "sampled", "fixed", 9.0, 42), None, None, None),
        ("cantrell",      ModelConfig(NI, "sampled", "varying", 9.0, 42), None, None, 1.0054),
        ("saueressig",    ModelConfig(NI, "sampled", "varying", 9.0, 42), None, None, 1.0039),
        ("low_Cl_0.6pct", ModelConfig(NI, "sampled", "varying", 9.0, 42), None, 0.006, None),
        ("high_Cl_6.5pct",ModelConfig(NI, "sampled", "varying", 9.0, 42), None, 0.065, None),
        ("fast_ex_0.8",   ModelConfig(NI, "sampled", "varying", 9.0, 42), 0.8, None, None),
        ("slow_ex_1.3",   ModelConfig(NI, "sampled", "varying", 9.0, 42), 1.3, None, None),
    ]

    rows = []
    for name, cfg, tau_ex_ov, cl_ov, kie_ov in configs:
        print(f"\n  Config: {name}")
        years, FF_NH, FF_SH, Mic_NH, Mic_SH = run_config(
            data, cfg, tau_ex_override=tau_ex_ov, cl_frac_override=cl_ov,
            kie_oh13c_override=kie_ov, label=name)

        nh_ff_slopes = compute_trends(FF_NH, years, start=2007)
        sh_mic_slopes = compute_trends(Mic_SH, years, start=2007)
        glob_ff_slopes = compute_trends(FF_NH + FF_SH, years, start=2007)

        nh_ff_med = np.nanmedian(nh_ff_slopes)
        sh_mic_med = np.nanmedian(sh_mic_slopes)
        glob_ff_med = np.nanmedian(glob_ff_slopes)

        nh_ff_pos = np.nanmean(nh_ff_slopes > 0) * 100
        sh_mic_pos = np.nanmean(sh_mic_slopes > 0) * 100
        glob_ff_neg = np.nanmean(glob_ff_slopes < 0) * 100

        pattern = nh_ff_pos > 50 and sh_mic_pos > 60
        print(f"    NH_FF slope: {nh_ff_med:+.3f} ({nh_ff_pos:.0f}% pos)")
        print(f"    SH_Mic slope: {sh_mic_med:+.3f} ({sh_mic_pos:.0f}% pos)")
        print(f"    Global_FF slope: {glob_ff_med:+.3f} ({glob_ff_neg:.0f}% neg)")
        print(f"    Pattern holds: {'✓' if pattern else '✗'}")

        rows.append({
            "config": name,
            "NH_FF_slope": nh_ff_med,
            "NH_FF_pct_pos": nh_ff_pos,
            "NH_FF_p5": float(np.nanpercentile(nh_ff_slopes, 5)),
            "NH_FF_p95": float(np.nanpercentile(nh_ff_slopes, 95)),
            "SH_Mic_slope": sh_mic_med,
            "SH_Mic_pct_pos": sh_mic_pos,
            "Global_FF_slope": glob_ff_med,
            "Global_FF_pct_neg": glob_ff_neg,
            "pattern_holds": bool(pattern),
        })

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS_DIR / "robustness_table.csv", index=False)
    with open(RESULTS_DIR / "robustness_table.json", 'w') as f:
        json.dump(rows, f, indent=2)
    print(f"\n  Saved to {RESULTS_DIR}/robustness_table.csv")

    n_hold = sum(r['pattern_holds'] for r in rows)
    print(f"\n  Pattern holds in {n_hold}/{len(rows)} configurations")


if __name__ == "__main__":
    main()
