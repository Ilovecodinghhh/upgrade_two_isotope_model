#!/usr/bin/env python3
"""
OH KIE Importance Experiment
=============================
Quantifies how OH-13C and OH-D KIE uncertainties propagate into
fossil-fuel (FF) emission estimates in the one-box dual-isotope
(3×3) inversion.

Method: Run 5 configurations of the 3×3 one-box model, each with
1000 MC iterations, selectively freezing OH KIE parameters while
sampling everything else:

  1. ALL_SAMPLED    — baseline (all KIE + source sigs + data sampled)
  2. FIX_OH13C      — OH-13C fixed at midpoint, rest sampled
  3. FIX_OHD        — OH-D fixed at midpoint, rest sampled
  4. FIX_BOTH_OH    — both OH-13C and OH-D fixed, rest sampled
  5. ALL_KIE_FIXED  — all 8 KIE parameters fixed, everything else sampled

The FF spread (σ across MC) in each config reveals which KIE
parameter matters most.
"""

import sys, json
from pathlib import Path
import numpy as np
import pandas as pd

# Add repo root to path
REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from common import (
    load_data, sample_KIE, compute_bulk_KIE, compute_lifetime,
    delta_to_fraction_d13C, delta_to_fraction_dD,
    sample_source_signatures, sample_atm_d13C, sample_atm_dD,
    smooth_5yr, pad_to_length,
    SINK_FRACTIONS_GLOBAL, PT, KIE_FIXED, KIE_DISTRIBUTIONS,
)

OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)

NI = 1000
SEED = 42


def custom_sample_KIE(rng, freeze: dict) -> dict:
    """Sample KIE, but override specific keys with fixed values from `freeze`."""
    kies = {}
    for key, cfg in KIE_DISTRIBUTIONS.items():
        if key in freeze:
            kies[key] = freeze[key]
        elif cfg['dist'] == 'uniform':
            kies[key] = rng.uniform(cfg['low'], cfg['high'])
        elif cfg['dist'] == 'normal':
            kies[key] = rng.normal(cfg['mean'], cfg['std'])
    return kies


def run_config(name: str, freeze: dict, data, seed=SEED) -> dict:
    """Run one-box 3x3 model with given frozen KIE parameters."""
    print(f"\n{'='*60}")
    print(f"  Config: {name}")
    print(f"  Frozen: {list(freeze.keys()) if freeze else 'none'}")
    print(f"{'='*60}")

    rng = np.random.default_rng(seed)
    n = data.n_years
    years = data.model_years
    CH4 = data.CH4_global
    tau = compute_lifetime(years, "varying")

    SumSource = np.zeros(n)
    for i in range(n):
        SumSource[i] = CH4[i+1]*PT - CH4[i]*PT + CH4[i]*PT/tau[i]

    FF_comp = np.zeros((n, NI))
    Mic_comp = np.zeros((n, NI))
    BB_comp = np.zeros((n, NI))
    KIE_13C_vals = np.zeros(NI)
    KIE_D_vals = np.zeros(NI)
    n_nonphysical = 0

    for k in range(NI):
        if (k+1) % 500 == 0:
            print(f"    iter {k+1}/{NI}")

        kies = custom_sample_KIE(rng, freeze)
        KIE_13C, KIE_D = compute_bulk_KIE(kies, SINK_FRACTIONS_GLOBAL)
        KIE_13C_vals[k] = KIE_13C
        KIE_D_vals[k] = KIE_D
        alpha_13C = 1.0 / KIE_13C
        alpha_D = 1.0 / KIE_D

        d13C_atm = sample_atm_d13C(data, k, n)
        dD_atm = sample_atm_dD(data, k, n)

        f13 = delta_to_fraction_d13C(d13C_atm)
        fD = delta_to_fraction_dD(dD_atm)
        n13 = f13 * CH4 * PT
        nD = fD * CH4 * PT

        d13C_src = np.zeros(n)
        dD_src = np.zeros(n)
        for j in range(n):
            d13C_src[j] = (n13[j+1] - n13[j] + n13[j]*alpha_13C/tau[j]) / SumSource[j]
            dD_src[j] = (nD[j+1] - nD[j] + nD[j]*alpha_D/tau[j]) / SumSource[j]

        sigs = sample_source_signatures(rng, data, k, n)

        f13_bb = delta_to_fraction_d13C(sigs['bb_d13C'])
        f13_ff = delta_to_fraction_d13C(sigs['ff_d13C'])
        f13_mic = delta_to_fraction_d13C(sigs['mic_d13C'])
        fD_bb = delta_to_fraction_dD(sigs['bb_dD'])
        fD_ff = delta_to_fraction_dD(sigs['ff_dD'])
        fD_mic = delta_to_fraction_dD(sigs['mic_dD'])

        for j in range(n):
            A = np.array([
                [1.0,       1.0,       1.0],
                [f13_bb[j], f13_ff[j], f13_mic[j]],
                [fD_bb[j],  fD_ff[j],  fD_mic[j]],
            ])
            B = np.array([
                SumSource[j],
                SumSource[j] * d13C_src[j],
                SumSource[j] * dD_src[j],
            ])
            try:
                x = np.linalg.solve(A, B)
            except np.linalg.LinAlgError:
                x = np.array([np.nan, np.nan, np.nan])

            if np.any(x < 0) or np.any(~np.isfinite(x)):
                n_nonphysical += 1

            BB_comp[j, k] = x[0]
            FF_comp[j, k] = x[1]
            Mic_comp[j, k] = x[2]

    # Smooth
    FF_s = smooth_5yr(FF_comp)
    Mic_s = smooth_5yr(Mic_comp)
    BB_s = smooth_5yr(BB_comp)

    # Stats
    ff_mean_ts = np.nanmean(FF_s, axis=1)   # (n_years,)
    ff_std_ts = np.nanstd(FF_s, axis=1)

    # Time-mean σ(FF) across all years
    sigma_ff = float(np.nanmean(ff_std_ts))
    # Mean FF level
    mean_ff = float(np.nanmean(ff_mean_ts))
    # Trend: mean of last 3 years minus mean of 2005-2007
    yr0 = int(years[0])
    i0, i1 = 2005 - yr0, 2007 - yr0 + 1
    base = np.nanmean(FF_s[i0:i1], axis=0)
    recent = np.nanmean(FF_s[-3:], axis=0)
    delta_ff = recent - base
    trend_mean = float(np.nanmean(delta_ff))
    trend_std = float(np.nanstd(delta_ff))
    pct_pos = float((delta_ff > 0).sum() / len(delta_ff) * 100)

    # 90% CI of annual FF means
    ff_5 = float(np.nanpercentile(FF_s, 5))
    ff_95 = float(np.nanpercentile(FF_s, 95))

    result = {
        'name': name,
        'frozen_keys': list(freeze.keys()),
        'n_iterations': NI,
        'sigma_ff': round(sigma_ff, 2),
        'mean_ff': round(mean_ff, 1),
        'trend_ff': round(trend_mean, 2),
        'trend_ff_std': round(trend_std, 2),
        'trend_pct_positive': round(pct_pos, 1),
        'ff_5pct': round(ff_5, 1),
        'ff_95pct': round(ff_95, 1),
        'nonphysical_pct': round(100*n_nonphysical/(n*NI), 1),
        'KIE_13C_mean': round(float(KIE_13C_vals.mean()), 6),
        'KIE_13C_std': round(float(KIE_13C_vals.std()), 6),
        'KIE_D_mean': round(float(KIE_D_vals.mean()), 4),
        'KIE_D_std': round(float(KIE_D_vals.std()), 4),
    }

    # Save per-year time series
    ts_df = pd.DataFrame({
        'year': years,
        'FF_mean': ff_mean_ts,
        'FF_std': ff_std_ts,
        'Mic_mean': np.nanmean(Mic_s, axis=1),
        'Mic_std': np.nanstd(Mic_s, axis=1),
        'BB_mean': np.nanmean(BB_s, axis=1),
        'BB_std': np.nanstd(BB_s, axis=1),
    })
    ts_df.to_csv(OUT / f"{name}_timeseries.csv", index=False)

    print(f"  σ(FF) = {sigma_ff:.2f} Tg/yr | mean FF = {mean_ff:.1f} | "
          f"ΔFF = {trend_mean:+.1f}±{trend_std:.1f} | "
          f"{pct_pos:.0f}% positive | nonphys: {result['nonphysical_pct']}%")

    return result


def main():
    print("Loading data...")
    data = load_data(REPO, two_box=False)
    print(f"  n_years = {data.n_years}, years = {data.model_years[0]}–{data.model_years[-1]}")

    # Central values for freezing
    oh13c_mid = KIE_FIXED['OH_13C']  # 1.00465
    ohd_mid = KIE_FIXED['OH_D']      # 1.3105

    configs = [
        ("ALL_SAMPLED",   {}),
        ("FIX_OH13C",     {'OH_13C': oh13c_mid}),
        ("FIX_OHD",       {'OH_D': ohd_mid}),
        ("FIX_BOTH_OH",   {'OH_13C': oh13c_mid, 'OH_D': ohd_mid}),
        ("ALL_KIE_FIXED", dict(KIE_FIXED)),
    ]

    # Also run Saueressig-only and Cantrell-only for OH-13C
    saueressig = KIE_DISTRIBUTIONS['OH_13C']['low']   # 1.0039
    cantrell = KIE_DISTRIBUTIONS['OH_13C']['high']     # 1.0054
    configs.append(("OH13C_SAUERESSIG", {'OH_13C': saueressig}))
    configs.append(("OH13C_CANTRELL",   {'OH_13C': cantrell}))

    all_results = {}
    for name, freeze in configs:
        result = run_config(name, freeze, data)
        all_results[name] = result

    # Save summary
    with open(OUT / "summary.json", "w") as f:
        json.dump(all_results, f, indent=2)

    print("\n" + "="*60)
    print("SUMMARY TABLE")
    print("="*60)
    print(f"{'Config':<22} {'σ(FF)':>8} {'mean FF':>8} {'ΔFF':>8} {'nonphys%':>9}")
    print("-"*60)
    for name, r in all_results.items():
        print(f"{name:<22} {r['sigma_ff']:>8.2f} {r['mean_ff']:>8.1f} "
              f"{r['trend_ff']:>+8.1f} {r['nonphysical_pct']:>8.1f}%")

    # Variance attribution
    print("\n" + "="*60)
    print("VARIANCE ATTRIBUTION (σ² reduction from baseline)")
    print("="*60)
    s2_base = all_results['ALL_SAMPLED']['sigma_ff']**2
    for name in ['FIX_OH13C', 'FIX_OHD', 'FIX_BOTH_OH', 'ALL_KIE_FIXED']:
        s2 = all_results[name]['sigma_ff']**2
        reduction = s2_base - s2
        pct = 100 * reduction / s2_base if s2_base > 0 else 0
        print(f"  {name:<22}: σ²={s2:.1f} → reduction={reduction:.1f} ({pct:.1f}% of baseline)")

    print(f"\n  OH-13C alone: {100*(s2_base - all_results['FIX_OH13C']['sigma_ff']**2)/s2_base:.1f}%")
    print(f"  OH-D alone:   {100*(s2_base - all_results['FIX_OHD']['sigma_ff']**2)/s2_base:.1f}%")
    print(f"  Both OH:      {100*(s2_base - all_results['FIX_BOTH_OH']['sigma_ff']**2)/s2_base:.1f}%")
    print(f"  All KIE:      {100*(s2_base - all_results['ALL_KIE_FIXED']['sigma_ff']**2)/s2_base:.1f}%")

    print(f"\nDone. Results in {OUT}/")


if __name__ == "__main__":
    main()
