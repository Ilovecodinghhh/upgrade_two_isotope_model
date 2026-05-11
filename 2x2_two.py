#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2×2_two.py — BB-Fixed, Separate δ¹³C & δD, Two-Box (NH/SH)
=============================================================

Approach:
  Two-hemisphere version of the 2×2 BB-fixed model.
  For each hemisphere separately:
    - Fix BB from CarbonTracker × hemisphere fraction (GFED4: 55/45)
    - Solve for FF and Mic using each isotope independently
    - Include interhemispheric exchange (τ_ex ~ N(1.0, 0.1) yr)

  Yields 4 time series: FF(δ¹³C), Mic(δ¹³C), FF(δD), Mic(δD) per hemisphere.

Equivalent to v3.2 in the old numbering.

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
    ModelConfig, load_data, sample_KIE, compute_bulk_KIE, compute_lifetime,
    delta_to_fraction_d13C, delta_to_fraction_dD,
    fraction_to_delta_d13C, fraction_to_delta_dD,
    sample_source_signatures, sample_atm_d13C, sample_atm_dD,
    smooth_5yr, trend_change, pad_to_length,
    SINK_FRACTIONS_NH, SINK_FRACTIONS_SH,
    PT, PT_HEMI, LIFETIME_RATIO_NH, LIFETIME_RATIO_SH,
    BB_NH_FRACTION, BB_SH_FRACTION, DD_IH_OFFSET,
    TAU_EX_MEAN, TAU_EX_STD,
)

BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "Output_2x2_two"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def run(cfg: ModelConfig):
    print("=" * 70)
    print("2×2 TWO-BOX MODEL (BB-Fixed, Separate δ¹³C & δD, NH/SH)")
    print(f"  KIE: {cfg.kie_mode} | Lifetime: {cfg.lifetime_mode} | N={cfg.n_iterations}")
    print("=" * 70)

    data = load_data(BASE_DIR, two_box=True)
    n = data.n_years
    years = data.model_years
    rng = np.random.default_rng(cfg.seed)

    CH4_NH, CH4_SH = data.CH4_NH, data.CH4_SH
    c13_NH, c13_SH = data.c13_NH, data.c13_SH
    c13_glob = data.c13_global

    tau_global = compute_lifetime(years, cfg.lifetime_mode, cfg.tau_fixed)
    tau_NH = tau_global * LIFETIME_RATIO_NH
    tau_SH = tau_global * LIFETIME_RATIO_SH

    BB_NH_fixed = data.BB_annual * BB_NH_FRACTION
    BB_SH_fixed = data.BB_annual * BB_SH_FRACTION

    NI = cfg.n_iterations
    FF_NH_c = np.zeros((n, NI)); Mic_NH_c = np.zeros((n, NI))
    FF_SH_c = np.zeros((n, NI)); Mic_SH_c = np.zeros((n, NI))
    FF_NH_d = np.zeros((n, NI)); Mic_NH_d = np.zeros((n, NI))
    FF_SH_d = np.zeros((n, NI)); Mic_SH_d = np.zeros((n, NI))
    n_neg_c, n_neg_d, total = 0, 0, 0

    for k in range(NI):
        if (k + 1) % 200 == 0:
            print(f"  iter {k+1}/{NI}")

        tau_ex = max(0.5, rng.normal(TAU_EX_MEAN, TAU_EX_STD))
        kies = sample_KIE(rng, cfg.kie_mode)
        K13_NH, KD_NH = compute_bulk_KIE(kies, SINK_FRACTIONS_NH)
        K13_SH, KD_SH = compute_bulk_KIE(kies, SINK_FRACTIONS_SH)
        a13_NH = 1.0 / K13_NH; aD_NH = 1.0 / KD_NH
        a13_SH = 1.0 / K13_SH; aD_SH = 1.0 / KD_SH

        # Total source per hemisphere
        S_NH = np.zeros(n); S_SH = np.zeros(n)
        for i in range(n):
            M_NH = CH4_NH[i] * PT_HEMI; M_NH1 = CH4_NH[i+1] * PT_HEMI
            M_SH = CH4_SH[i] * PT_HEMI; M_SH1 = CH4_SH[i+1] * PT_HEMI
            ex_NH = (M_SH - M_NH) / tau_ex
            ex_SH = (M_NH - M_SH) / tau_ex
            S_NH[i] = (M_NH1 - M_NH) + M_NH / tau_NH[i] - ex_NH
            S_SH[i] = (M_SH1 - M_SH) + M_SH / tau_SH[i] - ex_SH

        # Atmospheric observations
        d13C_glob_MC = sample_atm_d13C(data, k, n)
        dD_glob_MC = sample_atm_dD(data, k, n)

        nc = min(len(c13_glob), n + 1)
        d13C_off = d13C_glob_MC[:nc] - c13_glob[:nc]
        d13C_NH_MC = c13_NH[:nc] + d13C_off
        d13C_SH_MC = c13_SH[:nc] + d13C_off
        dD_NH_MC = dD_glob_MC - DD_IH_OFFSET
        dD_SH_MC = dD_glob_MC + DD_IH_OFFSET

        f13_NH = delta_to_fraction_d13C(d13C_NH_MC)
        f13_SH = delta_to_fraction_d13C(d13C_SH_MC)
        fD_NH = delta_to_fraction_dD(dD_NH_MC)
        fD_SH = delta_to_fraction_dD(dD_SH_MC)

        # Isotopic source fractions
        d13C_src_NH = np.zeros(n); d13C_src_SH = np.zeros(n)
        dD_src_NH = np.zeros(n); dD_src_SH = np.zeros(n)
        for j in range(n):
            # ¹³C NH
            n13_NH = f13_NH[j] * CH4_NH[j] * PT_HEMI
            n13_NH1 = f13_NH[j+1] * CH4_NH[j+1] * PT_HEMI
            n13_SH = f13_SH[j] * CH4_SH[j] * PT_HEMI
            ex13_NH = (n13_SH - n13_NH) / tau_ex
            d13C_src_NH[j] = (n13_NH1 - n13_NH + n13_NH * a13_NH / tau_NH[j] - ex13_NH) / S_NH[j]

            # ¹³C SH
            n13_SH1 = f13_SH[j+1] * CH4_SH[j+1] * PT_HEMI
            ex13_SH = (n13_NH - n13_SH) / tau_ex
            d13C_src_SH[j] = (n13_SH1 - n13_SH + n13_SH * a13_SH / tau_SH[j] - ex13_SH) / S_SH[j]

            # D NH
            nD_NH = fD_NH[j] * CH4_NH[j] * PT_HEMI
            nD_NH1 = fD_NH[j+1] * CH4_NH[j+1] * PT_HEMI
            nD_SH = fD_SH[j] * CH4_SH[j] * PT_HEMI
            exD_NH = (nD_SH - nD_NH) / tau_ex
            dD_src_NH[j] = (nD_NH1 - nD_NH + nD_NH * aD_NH / tau_NH[j] - exD_NH) / S_NH[j]

            # D SH
            nD_SH1 = fD_SH[j+1] * CH4_SH[j+1] * PT_HEMI
            exD_SH = (nD_NH - nD_SH) / tau_ex
            dD_src_SH[j] = (nD_SH1 - nD_SH + nD_SH * aD_SH / tau_SH[j] - exD_SH) / S_SH[j]

        sigs = sample_source_signatures(rng, data, k, n)

        for j in range(n):
            total += 1
            d13c_NH_d = fraction_to_delta_d13C(d13C_src_NH[j])
            d13c_SH_d = fraction_to_delta_d13C(d13C_src_SH[j])
            dD_NH_delta = fraction_to_delta_dD(dD_src_NH[j])
            dD_SH_delta = fraction_to_delta_dD(dD_src_SH[j])

            BB_NH_j = BB_NH_fixed[j] if j < len(BB_NH_fixed) else data.BB_global_mean * BB_NH_FRACTION
            BB_SH_j = BB_SH_fixed[j] if j < len(BB_SH_fixed) else data.BB_global_mean * BB_SH_FRACTION

            denom_c = sigs['ff_d13C'][j] - sigs['mic_d13C'][j]
            denom_d = sigs['ff_dD'][j] - sigs['mic_dD'][j]

            def solve_2x2(S, BB_j, delta_src, denom, sig_ff, sig_mic, sig_bb):
                if abs(denom) < (0.1 if 'd13C' in str(sig_ff) else 1.0):
                    return np.nan, np.nan
                ff = (S * delta_src - sig_mic * (S - BB_j) - sig_bb * BB_j) / denom
                mic = S - BB_j - ff
                return ff, mic

            def clamp(ff, mic, S, BB_j):
                neg = False
                if not np.isnan(ff) and (ff < 0 or mic < 0):
                    neg = True
                    ff = max(0, ff); mic = S - BB_j - ff
                    if mic < 0: mic = 0; ff = S - BB_j
                return ff, mic, neg

            # NH δ¹³C
            ff, mic = solve_2x2(S_NH[j], BB_NH_j, d13c_NH_d, denom_c,
                                sigs['ff_d13C'][j], sigs['mic_d13C'][j], sigs['bb_d13C'][j])
            ff, mic, neg = clamp(ff, mic, S_NH[j], BB_NH_j)
            if neg: n_neg_c += 1
            FF_NH_c[j, k] = ff; Mic_NH_c[j, k] = mic

            # SH δ¹³C
            ff, mic = solve_2x2(S_SH[j], BB_SH_j, d13c_SH_d, denom_c,
                                sigs['ff_d13C'][j], sigs['mic_d13C'][j], sigs['bb_d13C'][j])
            ff, mic, neg = clamp(ff, mic, S_SH[j], BB_SH_j)
            if neg: n_neg_c += 1
            FF_SH_c[j, k] = ff; Mic_SH_c[j, k] = mic

            # NH δD
            ff, mic = solve_2x2(S_NH[j], BB_NH_j, dD_NH_delta, denom_d,
                                sigs['ff_dD'][j], sigs['mic_dD'][j], sigs['bb_dD'][j])
            ff, mic, neg = clamp(ff, mic, S_NH[j], BB_NH_j)
            if neg: n_neg_d += 1
            FF_NH_d[j, k] = ff; Mic_NH_d[j, k] = mic

            # SH δD
            ff, mic = solve_2x2(S_SH[j], BB_SH_j, dD_SH_delta, denom_d,
                                sigs['ff_dD'][j], sigs['mic_dD'][j], sigs['bb_dD'][j])
            ff, mic, neg = clamp(ff, mic, S_SH[j], BB_SH_j)
            if neg: n_neg_d += 1
            FF_SH_d[j, k] = ff; Mic_SH_d[j, k] = mic

    print(f"\n  Negatives: δ¹³C={n_neg_c}, δD={n_neg_d}")

    # Global = NH + SH
    FF_G_c = FF_NH_c + FF_SH_c; Mic_G_c = Mic_NH_c + Mic_SH_c
    FF_G_d = FF_NH_d + FF_SH_d; Mic_G_d = Mic_NH_d + Mic_SH_d

    # Smooth
    FF_G_cs = smooth_5yr(FF_G_c); Mic_G_cs = smooth_5yr(Mic_G_c)
    FF_G_ds = smooth_5yr(FF_G_d); Mic_G_ds = smooth_5yr(Mic_G_d)

    print(f"\n{'='*60}")
    print("RESULTS (Global, Smoothed)")
    for lbl, arr in [('FF δ¹³C', FF_G_cs), ('Mic δ¹³C', Mic_G_cs),
                     ('FF δD', FF_G_ds), ('Mic δD', Mic_G_ds)]:
        print(f"  {lbl}: {np.nanmean(arr):.1f} ± {np.nanstd(np.nanmean(arr,0)):.1f} Tg/yr")
    print(f"  BB (fixed): {data.BB_global_mean:.1f} Tg/yr")

    print(f"\n  --- TRENDS ---")
    for lbl, arr in [('FF δ¹³C', FF_G_cs), ('Mic δ¹³C', Mic_G_cs),
                     ('FF δD', FF_G_ds), ('Mic δD', Mic_G_ds)]:
        d, p = trend_change(arr, years)
        print(f"  {lbl}: Δ={d.mean():+.1f}±{d.std():.1f} ({p:.0f}% pos)")

    # Save
    pd.DataFrame({
        'Year': years,
        'FF_G_d13C_mean': np.nanmean(FF_G_cs, 1), 'FF_G_d13C_std': np.nanstd(FF_G_cs, 1),
        'Mic_G_d13C_mean': np.nanmean(Mic_G_cs, 1), 'Mic_G_d13C_std': np.nanstd(Mic_G_cs, 1),
        'FF_G_dD_mean': np.nanmean(FF_G_ds, 1), 'FF_G_dD_std': np.nanstd(FF_G_ds, 1),
        'Mic_G_dD_mean': np.nanmean(Mic_G_ds, 1), 'Mic_G_dD_std': np.nanstd(Mic_G_ds, 1),
        'FF_NH_d13C_mean': np.nanmean(smooth_5yr(FF_NH_c), 1),
        'Mic_NH_d13C_mean': np.nanmean(smooth_5yr(Mic_NH_c), 1),
        'FF_SH_d13C_mean': np.nanmean(smooth_5yr(FF_SH_c), 1),
        'Mic_SH_d13C_mean': np.nanmean(smooth_5yr(Mic_SH_c), 1),
        'BB_fixed': data.BB_annual[:n],
    }).to_csv(OUT_DIR / 'results.csv', index=False)

    with open(OUT_DIR / 'config.json', 'w') as f:
        json.dump({'kie_mode': cfg.kie_mode, 'lifetime_mode': cfg.lifetime_mode,
                   'n_iterations': NI}, f, indent=2)

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), dpi=150, sharex=True)
    fig.suptitle(f'2×2 Two-Box (KIE={cfg.kie_mode}, τ={cfg.lifetime_mode})', fontsize=13)

    for ax, lbl, mc, md in [
        (axes[0,0], 'Global FF', FF_G_cs, FF_G_ds),
        (axes[0,1], 'Global Mic', Mic_G_cs, Mic_G_ds),
    ]:
        m_c = np.nanmean(mc, 1); m_d = np.nanmean(md, 1)
        s_c = 2*np.nanstd(mc, 1); s_d = 2*np.nanstd(md, 1)
        ax.plot(years, m_c, 'r-', lw=2, label='δ¹³C'); ax.fill_between(years, m_c-s_c, m_c+s_c, alpha=.2, color='r')
        ax.plot(years, m_d, 'b-', lw=2, label='δD'); ax.fill_between(years, m_d-s_d, m_d+s_d, alpha=.2, color='b')
        ax.set_title(lbl); ax.set_ylabel('Tg/yr'); ax.legend(); ax.grid(alpha=.3)

    for ax, lbl, ff, mic in [
        (axes[1,0], 'NH (δ¹³C)', smooth_5yr(FF_NH_c), smooth_5yr(Mic_NH_c)),
        (axes[1,1], 'SH (δ¹³C)', smooth_5yr(FF_SH_c), smooth_5yr(Mic_SH_c)),
    ]:
        ax.plot(years, np.nanmean(ff,1), 'b-', lw=2, label='FF')
        ax.plot(years, np.nanmean(mic,1), 'g-', lw=2, label='Mic')
        ax.set_title(lbl); ax.set_ylabel('Tg/yr'); ax.set_xlabel('Year')
        ax.legend(); ax.grid(alpha=.3)

    plt.tight_layout()
    plt.savefig(OUT_DIR / 'results.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n  Saved to {OUT_DIR}/")
    print("=" * 70)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="2×2 two-box methane model")
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
