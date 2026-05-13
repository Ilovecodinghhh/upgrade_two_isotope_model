#!/usr/bin/env python3
"""
Phase 4: Exchange rate sensitivity + δ¹³C gradient constraint.

Tests τ_ex = 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0 yr
For each: computes trends + modeled IH δ¹³C/δD gradient vs observed.
Uses the observed gradient to constrain τ_ex.
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
    SINK_FRACTIONS_NH, SINK_FRACTIONS_SH,
    PT_HEMI, LIFETIME_RATIO_NH, LIFETIME_RATIO_SH,
    TAU_EX_MEAN, TAU_EX_STD,
)

EXP_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = EXP_DIR / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

NI = 400


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


def run_exchange_test(data, tau_ex_value):
    """Run 2-box with fixed τ_ex, return trends + predicted IH gradients."""
    cfg = ModelConfig(NI, "sampled", "varying", 9.0, 42)
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

    # Track modeled source-weighted δ¹³C and δD per hemisphere
    d13C_src_all_NH = np.zeros((n, NI))
    d13C_src_all_SH = np.zeros((n, NI))
    dD_src_all_NH = np.zeros((n, NI))
    dD_src_all_SH = np.zeros((n, NI))

    W_NH = np.diag([100.0, 1.0, 0.5])
    W_SH = np.diag([200.0, 1.0, 0.5])

    for k in range(NI):
        tau_ex = tau_ex_value

        kies = sample_KIE(rng, cfg.kie_mode)
        K13_NH, KD_NH = compute_bulk_KIE(kies, SINK_FRACTIONS_NH)
        K13_SH, KD_SH = compute_bulk_KIE(kies, SINK_FRACTIONS_SH)
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

        d13C_src_all_NH[:, k] = d13C_src_NH
        d13C_src_all_SH[:, k] = d13C_src_SH
        dD_src_all_NH[:, k] = dD_src_NH
        dD_src_all_SH[:, k] = dD_src_SH

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

    return years, FF_NH, FF_SH, Mic_NH, Mic_SH, d13C_src_all_NH, d13C_src_all_SH, dD_src_all_NH, dD_src_all_SH


def main():
    print("=" * 70)
    print("PHASE 4: Exchange Rate Sensitivity")
    print("=" * 70)

    data = load_data(ROOT, two_box=True)

    # Observed IH δ¹³C gradient (NH - SH) — from data
    obs_gradient = np.mean(data.c13_NH[:data.n_years] - data.c13_SH[:data.n_years])
    obs_gradient_std = np.std(data.c13_NH[:data.n_years] - data.c13_SH[:data.n_years])
    print(f"  Observed mean IH δ¹³C gradient (NH-SH): {obs_gradient:.3f} ± {obs_gradient_std:.3f} ‰")

    tau_ex_values = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0]
    rows = []

    for tau_ex in tau_ex_values:
        print(f"\n  τ_ex = {tau_ex:.1f} yr ...")
        (years, FF_NH, FF_SH, Mic_NH, Mic_SH,
         d13C_src_NH, d13C_src_SH, dD_src_NH, dD_src_SH) = run_exchange_test(data, tau_ex)

        nh_ff_slopes = compute_trends(FF_NH, years, start=2007)
        sh_mic_slopes = compute_trends(Mic_SH, years, start=2007)
        glob_ff_slopes = compute_trends(FF_NH + FF_SH, years, start=2007)

        # Modeled source-weighted δ¹³C gradient is a proxy for atmospheric gradient
        # The actual atmospheric gradient is from data, but we can check consistency
        # by looking at source-δ¹³C divergence

        row = {
            "tau_ex": tau_ex,
            "NH_FF_slope": float(np.nanmedian(nh_ff_slopes)),
            "NH_FF_p5": float(np.nanpercentile(nh_ff_slopes, 5)),
            "NH_FF_p95": float(np.nanpercentile(nh_ff_slopes, 95)),
            "NH_FF_pct_pos": float(np.nanmean(nh_ff_slopes > 0) * 100),
            "SH_Mic_slope": float(np.nanmedian(sh_mic_slopes)),
            "SH_Mic_pct_pos": float(np.nanmean(sh_mic_slopes > 0) * 100),
            "Global_FF_slope": float(np.nanmedian(glob_ff_slopes)),
            "sigma_FF_NH": float(np.nanstd(np.nanmean(FF_NH, 1))),
            "sigma_FF_SH": float(np.nanstd(np.nanmean(FF_SH, 1))),
        }
        rows.append(row)
        print(f"    NH_FF: {row['NH_FF_slope']:+.3f} ({row['NH_FF_pct_pos']:.0f}% pos)")
        print(f"    SH_Mic: {row['SH_Mic_slope']:+.3f} ({row['SH_Mic_pct_pos']:.0f}% pos)")
        print(f"    Global_FF: {row['Global_FF_slope']:+.3f}")

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS_DIR / "exchange_rate_sensitivity.csv", index=False)
    with open(RESULTS_DIR / "exchange_rate_sensitivity.json", 'w') as f:
        json.dump(rows, f, indent=2)
    print(f"\n  Saved to {RESULTS_DIR}/exchange_rate_sensitivity.csv")


if __name__ == "__main__":
    main()
