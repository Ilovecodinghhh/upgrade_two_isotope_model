#!/usr/bin/env python3
"""
OH KIE Importance — 2×2 One-Box (separate δ¹³C and δD inversions)
==================================================================

Runs the 2×2 one-box model (BB fixed, FF+Mic solved independently
per isotope) under 7 KIE configurations.  This gives σ(FF) from
δ¹³C alone AND σ(FF) from δD alone, letting us see how OH-¹³C KIE
affects the ¹³C-derived FF, and how OH-D KIE affects the δD-derived FF.
"""

import sys, json
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from common import (
    load_data, compute_bulk_KIE, compute_lifetime,
    delta_to_fraction_d13C, delta_to_fraction_dD,
    fraction_to_delta_d13C, fraction_to_delta_dD,
    sample_source_signatures, sample_atm_d13C, sample_atm_dD,
    smooth_5yr, pad_to_length,
    SINK_FRACTIONS_GLOBAL, PT, KIE_FIXED, KIE_DISTRIBUTIONS,
)

OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)

NI = 1000
SEED = 42


def custom_sample_KIE(rng_kie, freeze: dict) -> dict:
    """Sample KIE values, always consuming 8 RNG draws to keep stream aligned.

    Even when a key is frozen, we still draw from rng_kie (and discard)
    so that downstream draws from the data RNG are never shifted.
    """
    kies = {}
    for key, cfg in KIE_DISTRIBUTIONS.items():
        if cfg['dist'] == 'uniform':
            drawn = rng_kie.uniform(cfg['low'], cfg['high'])
        elif cfg['dist'] == 'normal':
            drawn = rng_kie.normal(cfg['mean'], cfg['std'])
        kies[key] = freeze[key] if key in freeze else drawn
    return kies


def run_config(name, freeze, data, seed=SEED):
    print(f"\n{'='*60}")
    print(f"  2×2 Config: {name}")
    print(f"  Frozen: {list(freeze.keys()) if freeze else 'none'}")
    print(f"{'='*60}")

    # Use separate RNG streams: one for KIE, one for data/source signatures
    # This ensures freezing a KIE parameter never shifts the data RNG stream
    rng = np.random.default_rng(seed)
    rng_kie = np.random.default_rng(seed + 1000)
    n = data.n_years
    years = data.model_years
    CH4 = data.CH4_global
    tau = compute_lifetime(years, "varying")

    SumSource = np.zeros(n)
    for i in range(n):
        SumSource[i] = CH4[i+1]*PT - CH4[i]*PT + CH4[i]*PT/tau[i]

    BB = data.BB_annual.copy()

    FF_c = np.zeros((n, NI))
    Mic_c = np.zeros((n, NI))
    FF_d = np.zeros((n, NI))
    Mic_d = np.zeros((n, NI))
    n_neg_c = n_neg_d = 0

    for k in range(NI):
        if (k+1) % 500 == 0:
            print(f"    iter {k+1}/{NI}")

        kies = custom_sample_KIE(rng_kie, freeze)
        KIE_13C, KIE_D = compute_bulk_KIE(kies, SINK_FRACTIONS_GLOBAL)
        alpha_13C = 1.0 / KIE_13C
        alpha_D = 1.0 / KIE_D

        d13C_atm = sample_atm_d13C(data, k, n)
        dD_atm = sample_atm_dD(data, k, n)

        f13 = delta_to_fraction_d13C(d13C_atm)
        fD = delta_to_fraction_dD(dD_atm)
        n13 = f13 * CH4 * PT
        nD = fD * CH4 * PT

        d13C_src_f = np.zeros(n)
        dD_src_f = np.zeros(n)
        for j in range(n):
            d13C_src_f[j] = (n13[j+1] - n13[j] + n13[j]*alpha_13C/tau[j]) / SumSource[j]
            dD_src_f[j] = (nD[j+1] - nD[j] + nD[j]*alpha_D/tau[j]) / SumSource[j]

        sigs = sample_source_signatures(rng, data, k, n)

        for j in range(n):
            S = SumSource[j]
            BB_j = BB[j] if j < len(BB) else data.BB_global_mean

            d13C_src_delta = fraction_to_delta_d13C(d13C_src_f[j])
            dD_src_delta = fraction_to_delta_dD(dD_src_f[j])

            # δ¹³C inversion
            denom_c_val = sigs['ff_d13C'][j] - sigs['mic_d13C'][j]
            if abs(denom_c_val) > 0.1:
                ff_cv = (S * d13C_src_delta - sigs['mic_d13C'][j]*(S - BB_j) -
                         sigs['bb_d13C'][j]*BB_j) / denom_c_val
                mic_cv = S - BB_j - ff_cv
            else:
                ff_cv = mic_cv = np.nan

            if not np.isnan(ff_cv) and (ff_cv < 0 or mic_cv < 0):
                n_neg_c += 1
                ff_cv = max(0, ff_cv); mic_cv = S - BB_j - ff_cv
                if mic_cv < 0: mic_cv = 0; ff_cv = S - BB_j

            FF_c[j, k] = ff_cv
            Mic_c[j, k] = mic_cv

            # δD inversion
            denom_d_val = sigs['ff_dD'][j] - sigs['mic_dD'][j]
            if abs(denom_d_val) > 1.0:
                ff_dv = (S * dD_src_delta - sigs['mic_dD'][j]*(S - BB_j) -
                         sigs['bb_dD'][j]*BB_j) / denom_d_val
                mic_dv = S - BB_j - ff_dv
            else:
                ff_dv = mic_dv = np.nan

            if not np.isnan(ff_dv) and (ff_dv < 0 or mic_dv < 0):
                n_neg_d += 1
                ff_dv = max(0, ff_dv); mic_dv = S - BB_j - ff_dv
                if mic_dv < 0: mic_dv = 0; ff_dv = S - BB_j

            FF_d[j, k] = ff_dv
            Mic_d[j, k] = mic_dv

    FF_c_s = smooth_5yr(FF_c)
    FF_d_s = smooth_5yr(FF_d)

    sigma_ff_c = float(np.nanmean(np.nanstd(FF_c_s, axis=1)))
    sigma_ff_d = float(np.nanmean(np.nanstd(FF_d_s, axis=1)))
    mean_ff_c = float(np.nanmean(FF_c_s))
    mean_ff_d = float(np.nanmean(FF_d_s))

    # Trends
    yr0 = int(years[0])
    i0, i1 = 2005 - yr0, 2007 - yr0 + 1
    base_c = np.nanmean(FF_c_s[i0:i1], axis=0)
    recent_c = np.nanmean(FF_c_s[-3:], axis=0)
    delta_c = recent_c - base_c
    base_d = np.nanmean(FF_d_s[i0:i1], axis=0)
    recent_d = np.nanmean(FF_d_s[-3:], axis=0)
    delta_d = recent_d - base_d

    result = {
        'name': name,
        'frozen_keys': list(freeze.keys()),
        # δ¹³C-derived FF
        'sigma_ff_c13': round(sigma_ff_c, 2),
        'mean_ff_c13': round(mean_ff_c, 1),
        'trend_ff_c13': round(float(np.nanmean(delta_c)), 2),
        'trend_ff_c13_std': round(float(np.nanstd(delta_c)), 2),
        # δD-derived FF
        'sigma_ff_dD': round(sigma_ff_d, 2),
        'mean_ff_dD': round(mean_ff_d, 1),
        'trend_ff_dD': round(float(np.nanmean(delta_d)), 2),
        'trend_ff_dD_std': round(float(np.nanstd(delta_d)), 2),
        # non-physical
        'neg_c13_pct': round(100*n_neg_c / (n*NI), 1),
        'neg_dD_pct': round(100*n_neg_d / (n*NI), 1),
    }

    # Save time series
    ts = pd.DataFrame({
        'year': years,
        'FF_c13_mean': np.nanmean(FF_c_s, 1), 'FF_c13_std': np.nanstd(FF_c_s, 1),
        'FF_dD_mean': np.nanmean(FF_d_s, 1),   'FF_dD_std': np.nanstd(FF_d_s, 1),
        'Mic_c13_mean': np.nanmean(smooth_5yr(Mic_c), 1),
        'Mic_dD_mean': np.nanmean(smooth_5yr(Mic_d), 1),
    })
    ts.to_csv(OUT / f"2x2_{name}_timeseries.csv", index=False)

    print(f"  δ¹³C → σ(FF)={sigma_ff_c:.1f}, mean={mean_ff_c:.0f}, "
          f"ΔFF={np.nanmean(delta_c):+.1f}±{np.nanstd(delta_c):.1f}")
    print(f"  δD   → σ(FF)={sigma_ff_d:.1f}, mean={mean_ff_d:.0f}, "
          f"ΔFF={np.nanmean(delta_d):+.1f}±{np.nanstd(delta_d):.1f}")

    return result


def main():
    print("Loading data...")
    data = load_data(REPO, two_box=False)
    print(f"  n_years={data.n_years}, years={data.model_years[0]}–{data.model_years[-1]}")

    oh13c_mid = KIE_FIXED['OH_13C']
    ohd_mid = KIE_FIXED['OH_D']
    saueressig = KIE_DISTRIBUTIONS['OH_13C']['low']
    cantrell = KIE_DISTRIBUTIONS['OH_13C']['high']

    configs = [
        ("ALL_SAMPLED",       {}),
        ("FIX_OH13C",         {'OH_13C': oh13c_mid}),
        ("FIX_OHD",           {'OH_D': ohd_mid}),
        ("FIX_BOTH_OH",       {'OH_13C': oh13c_mid, 'OH_D': ohd_mid}),
        ("ALL_KIE_FIXED",     dict(KIE_FIXED)),
        ("OH13C_SAUERESSIG",  {'OH_13C': saueressig}),
        ("OH13C_CANTRELL",    {'OH_13C': cantrell}),
    ]

    all_results = {}
    for name, freeze in configs:
        all_results[name] = run_config(name, freeze, data)

    with open(OUT / "2x2_summary.json", "w") as f:
        json.dump(all_results, f, indent=2)

    # Print summary
    print("\n" + "="*80)
    print("2×2 ONE-BOX SUMMARY")
    print("="*80)
    print(f"{'Config':<22} {'σ(FF) δ¹³C':>11} {'σ(FF) δD':>10} "
          f"{'mean FF δ¹³C':>13} {'mean FF δD':>11}")
    print("-"*80)
    for name, r in all_results.items():
        print(f"{name:<22} {r['sigma_ff_c13']:>11.1f} {r['sigma_ff_dD']:>10.1f} "
              f"{r['mean_ff_c13']:>13.0f} {r['mean_ff_dD']:>11.0f}")

    # Variance attribution
    print("\n" + "="*80)
    print("VARIANCE ATTRIBUTION — δ¹³C-derived FF")
    print("="*80)
    s2b_c = all_results['ALL_SAMPLED']['sigma_ff_c13']**2
    for name in ['FIX_OH13C','FIX_OHD','FIX_BOTH_OH','ALL_KIE_FIXED']:
        s2 = all_results[name]['sigma_ff_c13']**2
        red = s2b_c - s2
        pct = 100*red/s2b_c if s2b_c > 0 else 0
        print(f"  {name:<22}: σ²={s2:.0f} → reduction {red:.0f} ({pct:.1f}%)")

    print("\n" + "="*80)
    print("VARIANCE ATTRIBUTION — δD-derived FF")
    print("="*80)
    s2b_d = all_results['ALL_SAMPLED']['sigma_ff_dD']**2
    for name in ['FIX_OH13C','FIX_OHD','FIX_BOTH_OH','ALL_KIE_FIXED']:
        s2 = all_results[name]['sigma_ff_dD']**2
        red = s2b_d - s2
        pct = 100*red/s2b_d if s2b_d > 0 else 0
        print(f"  {name:<22}: σ²={s2:.0f} → reduction {red:.0f} ({pct:.1f}%)")

    print(f"\nDone.  Results → {OUT}/")


if __name__ == "__main__":
    main()
