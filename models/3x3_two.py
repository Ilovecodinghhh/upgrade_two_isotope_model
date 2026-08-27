#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
3×3_two.py — Simultaneous δ¹³C + δD, Two-Box (NH/SH)
=======================================================

Approach:
  Two-hemisphere version of the 3×3 simultaneous model.
  For each hemisphere, solve the full 3×3 system:
    S_hemi = BB + FF + Mic
    S_hemi·f¹³C_src = BB·f¹³C_BB + FF·f¹³C_FF + Mic·f¹³C_Mic
    S_hemi·fD_src   = BB·fD_BB   + FF·fD_FF   + Mic·fD_Mic

  Uses bounded least squares (scipy lsq_linear) to enforce non-negativity,
  since the 3×3 system with δD can be ill-conditioned.

  Includes interhemispheric exchange, hemisphere-specific sink fractions,
  and δD NH/SH offset (±6‰).

Equivalent to v3.1 in the old numbering (optimized 3×3).

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
from scipy.optimize import lsq_linear

from common import (
    ModelConfig, QualityMonitor, load_data,
    sample_KIE, compute_bulk_KIE, compute_lifetime,
    delta_to_fraction_d13C, delta_to_fraction_dD,
    sample_source_signatures, sample_atm_d13C, sample_atm_dD,
    smooth_5yr, trend_change, pad_to_length,
    SINK_FRACTIONS_NH, SINK_FRACTIONS_SH,
    PT, PT_HEMI, LIFETIME_RATIO_NH, LIFETIME_RATIO_SH,
    DD_IH_OFFSET, TAU_EX_MEAN, TAU_EX_STD,
)

BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "Output_3x3_two"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def run(cfg: ModelConfig):
    print("=" * 70)
    print("3×3 TWO-BOX MODEL (Simultaneous δ¹³C + δD, NH/SH)")
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

    NI = cfg.n_iterations
    BB_NH = np.zeros((n, NI)); FF_NH = np.zeros((n, NI)); Mic_NH = np.zeros((n, NI))
    BB_SH = np.zeros((n, NI)); FF_SH = np.zeros((n, NI)); Mic_SH = np.zeros((n, NI))

    qm_NH = QualityMonitor(n, NI, "NH")
    qm_SH = QualityMonitor(n, NI, "SH")

    # Weighting matrices for bounded LS
    W_NH = np.diag([100.0, 1.0, 0.5])
    W_SH = np.diag([200.0, 1.0, 0.5])

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

        f13_NH_atm = delta_to_fraction_d13C(d13C_NH_MC)
        f13_SH_atm = delta_to_fraction_d13C(d13C_SH_MC)
        fD_NH_atm = delta_to_fraction_dD(dD_NH_MC)
        fD_SH_atm = delta_to_fraction_dD(dD_SH_MC)

        # Isotopic source fractions
        d13C_src_NH = np.zeros(n); d13C_src_SH = np.zeros(n)
        dD_src_NH = np.zeros(n); dD_src_SH = np.zeros(n)
        for j in range(n):
            n13_NH = f13_NH_atm[j] * CH4_NH[j] * PT_HEMI
            n13_NH1 = f13_NH_atm[j+1] * CH4_NH[j+1] * PT_HEMI
            n13_SH = f13_SH_atm[j] * CH4_SH[j] * PT_HEMI
            n13_SH1 = f13_SH_atm[j+1] * CH4_SH[j+1] * PT_HEMI
            ex13_NH = (n13_SH - n13_NH) / tau_ex
            ex13_SH = (n13_NH - n13_SH) / tau_ex
            d13C_src_NH[j] = (n13_NH1 - n13_NH + n13_NH * a13_NH / tau_NH[j] - ex13_NH) / S_NH[j]
            d13C_src_SH[j] = (n13_SH1 - n13_SH + n13_SH * a13_SH / tau_SH[j] - ex13_SH) / S_SH[j]

            nD_NH = fD_NH_atm[j] * CH4_NH[j] * PT_HEMI
            nD_NH1 = fD_NH_atm[j+1] * CH4_NH[j+1] * PT_HEMI
            nD_SH = fD_SH_atm[j] * CH4_SH[j] * PT_HEMI
            nD_SH1 = fD_SH_atm[j+1] * CH4_SH[j+1] * PT_HEMI
            exD_NH = (nD_SH - nD_NH) / tau_ex
            exD_SH = (nD_NH - nD_SH) / tau_ex
            dD_src_NH[j] = (nD_NH1 - nD_NH + nD_NH * aD_NH / tau_NH[j] - exD_NH) / S_NH[j]
            dD_src_SH[j] = (nD_SH1 - nD_SH + nD_SH * aD_SH / tau_SH[j] - exD_SH) / S_SH[j]

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

            # NH solve (bounded LS)
            B_nh = np.array([S_NH[j], S_NH[j]*d13C_src_NH[j], S_NH[j]*dD_src_NH[j]])
            A_w = W_NH @ A; B_w = W_NH @ B_nh
            ub = S_NH[j] * 1.5
            try:
                res = lsq_linear(A_w, B_w, bounds=(0, ub))
                x_nh = res.x
                qm_NH.record(j, k, A, x_nh)
            except Exception:
                x_nh = np.array([np.nan, np.nan, np.nan])
                qm_NH.is_nan[j, k] = True; qm_NH.is_nonphysical[j, k] = True
            BB_NH[j,k] = x_nh[0]; FF_NH[j,k] = x_nh[1]; Mic_NH[j,k] = x_nh[2]

            # SH solve
            B_sh = np.array([S_SH[j], S_SH[j]*d13C_src_SH[j], S_SH[j]*dD_src_SH[j]])
            A_w = W_SH @ A; B_w = W_SH @ B_sh
            ub = S_SH[j] * 1.5
            try:
                res = lsq_linear(A_w, B_w, bounds=(0, ub))
                x_sh = res.x
                qm_SH.record(j, k, A, x_sh)
            except Exception:
                x_sh = np.array([np.nan, np.nan, np.nan])
                qm_SH.is_nan[j, k] = True; qm_SH.is_nonphysical[j, k] = True
            BB_SH[j,k] = x_sh[0]; FF_SH[j,k] = x_sh[1]; Mic_SH[j,k] = x_sh[2]

    print("\nMC complete!")
    print("Quality:")
    qr_nh = qm_NH.summary()
    qr_sh = qm_SH.summary()

    # Global
    BB_G = BB_NH + BB_SH; FF_G = FF_NH + FF_SH; Mic_G = Mic_NH + Mic_SH

    # Smooth
    BB_Gs = smooth_5yr(BB_G); FF_Gs = smooth_5yr(FF_G); Mic_Gs = smooth_5yr(Mic_G)
    BB_NHs = smooth_5yr(BB_NH); FF_NHs = smooth_5yr(FF_NH); Mic_NHs = smooth_5yr(Mic_NH)
    BB_SHs = smooth_5yr(BB_SH); FF_SHs = smooth_5yr(FF_SH); Mic_SHs = smooth_5yr(Mic_SH)

    print(f"\n{'='*60}")
    print("RESULTS (Global, Smoothed)")
    for lbl, arr in [('BB', BB_Gs), ('FF', FF_Gs), ('Mic', Mic_Gs)]:
        print(f"  {lbl}: {np.nanmean(arr):.1f} ± {np.nanstd(np.nanmean(arr,0)):.1f} Tg/yr")

    print(f"\n  --- TRENDS ---")
    for lbl, arr in [('FF', FF_Gs), ('Mic', Mic_Gs), ('BB', BB_Gs)]:
        d, p = trend_change(arr, years)
        print(f"  {lbl}: Δ={d.mean():+.1f}±{d.std():.1f} ({p:.0f}% pos)")

    # Save
    pd.DataFrame({
        'Year': years,
        'BB_G_mean': np.nanmean(BB_Gs,1), 'BB_G_std': np.nanstd(BB_Gs,1),
        'FF_G_mean': np.nanmean(FF_Gs,1), 'FF_G_std': np.nanstd(FF_Gs,1),
        'Mic_G_mean': np.nanmean(Mic_Gs,1), 'Mic_G_std': np.nanstd(Mic_Gs,1),
        'BB_NH_mean': np.nanmean(BB_NHs,1), 'FF_NH_mean': np.nanmean(FF_NHs,1), 'Mic_NH_mean': np.nanmean(Mic_NHs,1),
        'BB_SH_mean': np.nanmean(BB_SHs,1), 'FF_SH_mean': np.nanmean(FF_SHs,1), 'Mic_SH_mean': np.nanmean(Mic_SHs,1),
    }).to_csv(OUT_DIR / 'results.csv', index=False)

    with open(OUT_DIR / 'config.json', 'w') as f:
        json.dump({'NH': qr_nh, 'SH': qr_sh,
                   'kie_mode': cfg.kie_mode, 'lifetime_mode': cfg.lifetime_mode,
                   'n_iterations': NI}, f, indent=2)

    # Plot
    fig, axes = plt.subplots(3, 2, figsize=(14, 12), dpi=150, sharex=True)
    fig.suptitle(f'3×3 Two-Box (KIE={cfg.kie_mode}, τ={cfg.lifetime_mode})', fontsize=13)

    for row, (name, color) in enumerate([('BB', 'red'), ('FF', 'blue'), ('Mic', 'green')]):
        for col, (hemi, s_arr) in enumerate([
            ('NH', {'BB': BB_NHs, 'FF': FF_NHs, 'Mic': Mic_NHs}),
            ('SH', {'BB': BB_SHs, 'FF': FF_SHs, 'Mic': Mic_SHs})]):
            ax = axes[row, col]
            m = np.nanmean(s_arr[name], 1); s = np.nanstd(s_arr[name], 1)
            ax.plot(years, m, '-', lw=2.5, color=color)
            ax.fill_between(years, m-s, m+s, alpha=.3, color=color)
            ax.set_ylabel(f'{name} (Tg/yr)')
            if row == 0: ax.set_title(f'{hemi} (5yr smoothed)')
            if row == 2: ax.set_xlabel('Year')
            ax.grid(alpha=.3)

    plt.tight_layout()
    plt.savefig(OUT_DIR / 'results.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n  Saved to {OUT_DIR}/")
    print("=" * 70)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="3×3 two-box methane model")
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
