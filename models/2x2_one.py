#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2×2_one.py — BB-Fixed, Separate δ¹³C & δD, Global One-Box
============================================================

Approach:
  Fix BB from CarbonTracker. For each isotope system SEPARATELY,
  solve a 2-equation system for FF and Microbial:
    S_total = FF + Mic + BB_fixed
    S_total × δ_source = FF × δ_FF + Mic × δ_Mic + BB_fixed × δ_BB
  → FF = (S × δ_src − δ_Mic × (S − BB) − δ_BB × BB) / (δ_FF − δ_Mic)
  → Mic = S − BB − FF

  δ¹³C and δD each yield INDEPENDENT FF/Mic estimates.
  Their agreement tests sensitivity to OH/Cl trends and KIE.

Advantages:
  - Well-conditioned (2×2 per isotope, no ill-conditioning)
  - Independent δ¹³C vs δD cross-validation
  - No bounded least squares needed

Limitation:
  - BB not independently constrained by isotopes

Configurable inputs:
  - KIE: fixed vs sampled
  - Lifetime: fixed τ=9yr vs time-varying τ(t)

Based on v3.2 logic minus hemisphere split.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from common import (
    ModelConfig, LoadedData, QualityMonitor,
    load_data, sample_KIE, compute_bulk_KIE, compute_lifetime,
    delta_to_fraction_d13C, delta_to_fraction_dD,
    fraction_to_delta_d13C, fraction_to_delta_dD,
    sample_source_signatures, sample_atm_d13C, sample_atm_dD,
    smooth_5yr, trend_change, pad_to_length,
    SINK_FRACTIONS_GLOBAL, PT, C13_STD, D_STD,
)

BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "Output_2x2_one"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def run(cfg: ModelConfig):
    print("=" * 70)
    print("2×2 ONE-BOX MODEL (BB-Fixed, Separate δ¹³C & δD)")
    print(f"  KIE: {cfg.kie_mode} | Lifetime: {cfg.lifetime_mode} | N={cfg.n_iterations}")
    print("=" * 70)

    data = load_data(BASE_DIR, two_box=False)
    n = data.n_years
    years = data.model_years
    CH4 = data.CH4_global
    rng = np.random.default_rng(cfg.seed)

    # Lifetime
    tau = compute_lifetime(years, cfg.lifetime_mode, cfg.tau_fixed)

    # Total source
    SumSource = np.zeros(n)
    for i in range(n):
        SumSource[i] = CH4[i + 1] * PT - CH4[i] * PT + CH4[i] * PT / tau[i]

    # BB fixed
    BB = data.BB_annual.copy()

    # Results
    FF_d13C = np.zeros((n, cfg.n_iterations))
    Mic_d13C = np.zeros((n, cfg.n_iterations))
    FF_dD = np.zeros((n, cfg.n_iterations))
    Mic_dD = np.zeros((n, cfg.n_iterations))
    n_neg_c, n_neg_d, total = 0, 0, 0

    for k in range(cfg.n_iterations):
        if (k + 1) % 200 == 0:
            print(f"  iter {k + 1}/{cfg.n_iterations}")

        kies = sample_KIE(rng, cfg.kie_mode)
        KIE_13C, KIE_D = compute_bulk_KIE(kies, SINK_FRACTIONS_GLOBAL)
        alpha_13C = 1.0 / KIE_13C
        alpha_D = 1.0 / KIE_D

        # Sample atmospheric observations
        d13C_atm = sample_atm_d13C(data, k, n)
        dD_atm = sample_atm_dD(data, k, n)

        f13 = delta_to_fraction_d13C(d13C_atm)
        fD = delta_to_fraction_dD(dD_atm)

        # Compute isotopic source fractions
        n13 = f13 * CH4 * PT
        nD = fD * CH4 * PT

        d13C_src_f = np.zeros(n)
        dD_src_f = np.zeros(n)
        for j in range(n):
            d13C_src_f[j] = (n13[j + 1] - n13[j] + n13[j] * alpha_13C / tau[j]) / SumSource[j]
            dD_src_f[j] = (nD[j + 1] - nD[j] + nD[j] * alpha_D / tau[j]) / SumSource[j]

        # Sample source signatures
        sigs = sample_source_signatures(rng, data, k, n)

        for j in range(n):
            total += 1
            S = SumSource[j]
            BB_j = BB[j] if j < len(BB) else data.BB_global_mean

            # Convert source fractions back to delta for formula
            d13C_src_delta = fraction_to_delta_d13C(d13C_src_f[j])
            dD_src_delta = fraction_to_delta_dD(dD_src_f[j])

            # --- δ¹³C inversion ---
            denom_c = sigs['ff_d13C'][j] - sigs['mic_d13C'][j]
            if abs(denom_c) > 0.1:
                ff_c = (S * d13C_src_delta - sigs['mic_d13C'][j] * (S - BB_j) -
                        sigs['bb_d13C'][j] * BB_j) / denom_c
                mic_c = S - BB_j - ff_c
            else:
                ff_c, mic_c = np.nan, np.nan

            if not np.isnan(ff_c) and (ff_c < 0 or mic_c < 0):
                n_neg_c += 1
                ff_c = max(0, ff_c)
                mic_c = S - BB_j - ff_c
                if mic_c < 0:
                    mic_c = 0; ff_c = S - BB_j

            FF_d13C[j, k] = ff_c
            Mic_d13C[j, k] = mic_c

            # --- δD inversion ---
            denom_d = sigs['ff_dD'][j] - sigs['mic_dD'][j]
            if abs(denom_d) > 1.0:
                ff_d = (S * dD_src_delta - sigs['mic_dD'][j] * (S - BB_j) -
                        sigs['bb_dD'][j] * BB_j) / denom_d
                mic_d = S - BB_j - ff_d
            else:
                ff_d, mic_d = np.nan, np.nan

            if not np.isnan(ff_d) and (ff_d < 0 or mic_d < 0):
                n_neg_d += 1
                ff_d = max(0, ff_d)
                mic_d = S - BB_j - ff_d
                if mic_d < 0:
                    mic_d = 0; ff_d = S - BB_j

            FF_dD[j, k] = ff_d
            Mic_dD[j, k] = mic_d

    print(f"\n  Negatives: δ¹³C={n_neg_c} ({100*n_neg_c/total:.1f}%), "
          f"δD={n_neg_d} ({100*n_neg_d/total:.1f}%)")

    # Smooth
    FF_d13C_s = smooth_5yr(FF_d13C)
    Mic_d13C_s = smooth_5yr(Mic_d13C)
    FF_dD_s = smooth_5yr(FF_dD)
    Mic_dD_s = smooth_5yr(Mic_dD)

    # Results
    print(f"\n{'='*60}")
    print("RESULTS (Smoothed)")
    print(f"{'='*60}")
    for lbl, arr in [('FF (δ¹³C)', FF_d13C_s), ('Mic (δ¹³C)', Mic_d13C_s),
                     ('FF (δD)', FF_dD_s), ('Mic (δD)', Mic_dD_s)]:
        m = np.nanmean(arr)
        s = np.nanstd(np.nanmean(arr, axis=0))
        print(f"  {lbl}: {m:.1f} ± {s:.1f} Tg/yr")
    print(f"  BB (fixed): {data.BB_global_mean:.1f} Tg/yr")

    # Trends
    print(f"\n  --- TRENDS: Δ(2020–2022 vs 2005–2007) ---")
    for lbl, arr in [('FF δ¹³C', FF_d13C_s), ('Mic δ¹³C', Mic_d13C_s),
                     ('FF δD', FF_dD_s), ('Mic δD', Mic_dD_s)]:
        delta, pct = trend_change(arr, years)
        print(f"  {lbl}: Δ={delta.mean():+.1f}±{delta.std():.1f} ({pct:.0f}% positive)")

    # Save
    import pandas as pd
    pd.DataFrame({
        'Year': years,
        'FF_d13C_mean': np.nanmean(FF_d13C_s, axis=1),
        'FF_d13C_std': np.nanstd(FF_d13C_s, axis=1),
        'Mic_d13C_mean': np.nanmean(Mic_d13C_s, axis=1),
        'Mic_d13C_std': np.nanstd(Mic_d13C_s, axis=1),
        'FF_dD_mean': np.nanmean(FF_dD_s, axis=1),
        'FF_dD_std': np.nanstd(FF_dD_s, axis=1),
        'Mic_dD_mean': np.nanmean(Mic_dD_s, axis=1),
        'Mic_dD_std': np.nanstd(Mic_dD_s, axis=1),
        'BB_fixed': BB[:n],
    }).to_csv(OUT_DIR / 'results.csv', index=False)

    with open(OUT_DIR / 'config.json', 'w') as f:
        json.dump({'kie_mode': cfg.kie_mode, 'lifetime_mode': cfg.lifetime_mode,
                   'n_iterations': cfg.n_iterations,
                   'neg_d13C_pct': round(100*n_neg_c/total, 2),
                   'neg_dD_pct': round(100*n_neg_d/total, 2)}, f, indent=2)

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), dpi=150)
    for ax, lbl, m_c, m_d, s_c, s_d in [
        (axes[0], 'Fossil Fuel',
         np.nanmean(FF_d13C_s, 1), np.nanmean(FF_dD_s, 1),
         2*np.nanstd(FF_d13C_s, 1), 2*np.nanstd(FF_dD_s, 1)),
        (axes[1], 'Microbial',
         np.nanmean(Mic_d13C_s, 1), np.nanmean(Mic_dD_s, 1),
         2*np.nanstd(Mic_d13C_s, 1), 2*np.nanstd(Mic_dD_s, 1)),
    ]:
        ax.plot(years, m_c, 'r-', lw=2, label='δ¹³C')
        ax.fill_between(years, m_c - s_c, m_c + s_c, alpha=0.2, color='red')
        ax.plot(years, m_d, 'b-', lw=2, label='δD')
        ax.fill_between(years, m_d - s_d, m_d + s_d, alpha=0.2, color='blue')
        ax.set_title(f'Global {lbl} (2×2 one-box)')
        ax.set_ylabel('Tg/yr'); ax.set_xlabel('Year')
        ax.legend(); ax.grid(alpha=0.3)

    plt.suptitle(f'2×2 One-Box (KIE={cfg.kie_mode}, τ={cfg.lifetime_mode})', y=1.02)
    plt.tight_layout()
    plt.savefig(OUT_DIR / 'results.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n  Saved to {OUT_DIR}/")
    print("=" * 70)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="2×2 one-box methane model")
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
