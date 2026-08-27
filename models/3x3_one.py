#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
3×3_one.py — Simultaneous δ¹³C + δD, Global One-Box
=====================================================

Approach:
  Solve a 3×3 system using both isotopes SIMULTANEOUSLY to partition
  three sources (BB, FF, Mic):
    S = BB + FF + Mic                     (mass balance)
    S·f¹³C_src = BB·f¹³C_BB + FF·f¹³C_FF + Mic·f¹³C_Mic   (¹³C)
    S·fD_src   = BB·fD_BB   + FF·fD_FF   + Mic·fD_Mic      (D)

  All three sources are free — BB is NOT fixed from CarbonTracker.

Advantage over 2×2:
  - BB is independently constrained by the dual-isotope system
  - Full 3-source partition from data alone

Limitation:
  - The 3×3 matrix can be ill-conditioned (δD row is 100× smaller in
    absolute scale than δ¹³C row), leading to large uncertainties
  - More non-physical solutions than 2×2

Equivalent to v2.0 in the old numbering (but with modular config).

Configurable inputs:
  - KIE: fixed vs sampled
  - Lifetime: fixed τ=9yr vs time-varying τ(t)
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from common import (
    ModelConfig, QualityMonitor, load_data,
    sample_KIE, compute_bulk_KIE, compute_lifetime,
    delta_to_fraction_d13C, delta_to_fraction_dD,
    sample_source_signatures, sample_atm_d13C, sample_atm_dD,
    smooth_5yr, trend_change, pad_to_length,
    SINK_FRACTIONS_GLOBAL, PT,
)

BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "Output_3x3_one"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def run(cfg: ModelConfig):
    print("=" * 70)
    print("3×3 ONE-BOX MODEL (Simultaneous δ¹³C + δD)")
    print(f"  KIE: {cfg.kie_mode} | Lifetime: {cfg.lifetime_mode} | N={cfg.n_iterations}")
    print("=" * 70)

    data = load_data(BASE_DIR, two_box=False)
    n = data.n_years
    years = data.model_years
    CH4 = data.CH4_global
    rng = np.random.default_rng(cfg.seed)

    tau = compute_lifetime(years, cfg.lifetime_mode, cfg.tau_fixed)

    # Total source
    SumSource = np.zeros(n)
    for i in range(n):
        SumSource[i] = CH4[i + 1] * PT - CH4[i] * PT + CH4[i] * PT / tau[i]

    NI = cfg.n_iterations
    BB_comp = np.zeros((n, NI))
    FF_comp = np.zeros((n, NI))
    Mic_comp = np.zeros((n, NI))

    qm = QualityMonitor(n, NI, "Global")

    for k in range(NI):
        if (k + 1) % 200 == 0:
            print(f"  iter {k+1}/{NI}")

        kies = sample_KIE(rng, cfg.kie_mode)
        KIE_13C, KIE_D = compute_bulk_KIE(kies, SINK_FRACTIONS_GLOBAL)
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
            d13C_src[j] = (n13[j+1] - n13[j] + n13[j] * alpha_13C / tau[j]) / SumSource[j]
            dD_src[j] = (nD[j+1] - nD[j] + nD[j] * alpha_D / tau[j]) / SumSource[j]

        sigs = sample_source_signatures(rng, data, k, n)

        f13_bb = delta_to_fraction_d13C(sigs['bb_d13C'])
        f13_ff = delta_to_fraction_d13C(sigs['ff_d13C'])
        f13_mic = delta_to_fraction_d13C(sigs['mic_d13C'])
        fD_bb = delta_to_fraction_dD(sigs['bb_dD'])
        fD_ff = delta_to_fraction_dD(sigs['ff_dD'])
        fD_mic = delta_to_fraction_dD(sigs['mic_dD'])

        for j in range(n):
            A = np.array([
                [1.0,        1.0,        1.0],
                [f13_bb[j],  f13_ff[j],  f13_mic[j]],
                [fD_bb[j],   fD_ff[j],   fD_mic[j]],
            ])
            B = np.array([
                SumSource[j],
                SumSource[j] * d13C_src[j],
                SumSource[j] * dD_src[j],
            ])
            try:
                x = np.linalg.solve(A, B)
                qm.record(j, k, A, x)
            except np.linalg.LinAlgError:
                x = np.array([np.nan, np.nan, np.nan])
                qm.is_nan[j, k] = True
                qm.is_nonphysical[j, k] = True

            BB_comp[j, k] = x[0]
            FF_comp[j, k] = x[1]
            Mic_comp[j, k] = x[2]

    print("\nMC complete!")
    qr = qm.summary()

    # Smooth
    BB_s = smooth_5yr(BB_comp)
    FF_s = smooth_5yr(FF_comp)
    Mic_s = smooth_5yr(Mic_comp)

    print(f"\n{'='*60}")
    print("RESULTS (Smoothed)")
    for lbl, arr in [('BB', BB_s), ('FF', FF_s), ('Mic', Mic_s)]:
        print(f"  {lbl}: {np.nanmean(arr):.1f} ± {np.nanstd(np.nanmean(arr,0)):.1f} Tg/yr")

    print(f"\n  --- TRENDS ---")
    for lbl, arr in [('FF', FF_s), ('Mic', Mic_s), ('BB', BB_s)]:
        d, p = trend_change(arr, years)
        print(f"  {lbl}: Δ={d.mean():+.1f}±{d.std():.1f} ({p:.0f}% pos)")

    # Save
    pd.DataFrame({
        'Year': years,
        'BB_mean': np.nanmean(BB_s, 1), 'BB_std': np.nanstd(BB_s, 1),
        'FF_mean': np.nanmean(FF_s, 1), 'FF_std': np.nanstd(FF_s, 1),
        'Mic_mean': np.nanmean(Mic_s, 1), 'Mic_std': np.nanstd(Mic_s, 1),
        'Lifetime': tau, 'TotalSource': SumSource,
    }).to_csv(OUT_DIR / 'results.csv', index=False)

    with open(OUT_DIR / 'config.json', 'w') as f:
        json.dump({**qr, 'kie_mode': cfg.kie_mode,
                   'lifetime_mode': cfg.lifetime_mode, 'n_iterations': NI}, f, indent=2)

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), dpi=150, sharex=True)
    fig.suptitle(f'3×3 One-Box (KIE={cfg.kie_mode}, τ={cfg.lifetime_mode})', fontsize=13)

    for ax, name, arr, color in [
        (axes[0,0], 'Biomass Burning', BB_s, 'red'),
        (axes[0,1], 'Fossil Fuel', FF_s, 'blue'),
        (axes[1,0], 'Microbial', Mic_s, 'green'),
    ]:
        m = np.nanmean(arr, 1); s = np.nanstd(arr, 1)
        ax.plot(years, m, '-', lw=2.5, color=color)
        ax.fill_between(years, m-s, m+s, alpha=.3, color=color)
        ax.set_title(name); ax.set_ylabel('Tg/yr')
        ax.grid(alpha=.3)

    ax = axes[1,1]
    cond_yr = np.mean(qm.condition_numbers, axis=1)
    ax.bar(years, cond_yr, color='purple', alpha=0.7)
    ax.set_title('Condition Number'); ax.set_ylabel('mean κ(A)')
    ax.set_xlabel('Year'); ax.grid(alpha=.3)
    axes[1,0].set_xlabel('Year')

    plt.tight_layout()
    plt.savefig(OUT_DIR / 'results.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n  Saved to {OUT_DIR}/")
    print("=" * 70)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="3×3 one-box methane model")
    parser.add_argument("--kie", choices=["fixed", "sampled"], default="sampled")
    parser.add_argument("--lifetime", choices=["fixed", "varying"], default="varying")
    parser.add_argument("--tau", type=float, default=9.0)
    parser.add_argument("--iterations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    cfg = ModelConfig(
        n_iterations=args.iterations, kie_mode=args.kie,
        lifetime_mode=args.lifetime, tau_fixed=args.tau, seed=args.seed)
    run(cfg)
