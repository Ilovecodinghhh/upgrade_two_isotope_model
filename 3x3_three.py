#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
3×3_three.py — Simultaneous δ¹³C + δD, Three-Box (NHext/Trop/SHext)
======================================================================

Three-box extension of 3×3_two.py. For each box, solve the full 3×3 system:
  S_box = BB + FF + Mic
  S_box·f¹³C_src = BB·f¹³C_BB + FF·f¹³C_FF + Mic·f¹³C_Mic
  S_box·fD_src   = BB·fD_BB   + FF·fD_FF   + Mic·fD_Mic

Uses bounded least squares (scipy lsq_linear) to enforce non-negativity.

Inter-box exchange:
  NHext ↔ Trop:  τ_NT ~ N(0.8, 0.1) yr
  Trop ↔ SHext:  τ_TS ~ N(1.2, 0.1) yr

δD atmospheric:  ThreeBox_atm_dD_annual.csv (real 3-box observations)
δD source sigs:  per-box MC matrices (NHext/Trop/SHext)
δ¹³C source sigs: global (shared across boxes — no hemispheric data yet)
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
    sample_source_signatures_three_box, sample_atm_d13C, sample_atm_dD,
    smooth_5yr, trend_change, pad_to_length,
    SINK_FRACTIONS_NHEXT, SINK_FRACTIONS_TROP, SINK_FRACTIONS_SHEXT,
    PT_NHEXT, PT_TROP, PT_SHEXT,
    LIFETIME_RATIO_NHEXT, LIFETIME_RATIO_TROP, LIFETIME_RATIO_SHEXT,
    BB_NHEXT_FRACTION, BB_TROP_FRACTION, BB_SHEXT_FRACTION,
    TAU_EX_NT_MEAN, TAU_EX_NT_STD, TAU_EX_TS_MEAN, TAU_EX_TS_STD,
)

BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "Output_3x3_three"
OUT_DIR.mkdir(parents=True, exist_ok=True)

BOXES = ['NHext', 'Trop', 'SHext']
SINK_FRACS = {
    'NHext': SINK_FRACTIONS_NHEXT,
    'Trop':  SINK_FRACTIONS_TROP,
    'SHext': SINK_FRACTIONS_SHEXT,
}
LT_RATIO = {
    'NHext': LIFETIME_RATIO_NHEXT,
    'Trop':  LIFETIME_RATIO_TROP,
    'SHext': LIFETIME_RATIO_SHEXT,
}
PT_BOX = {'NHext': PT_NHEXT, 'Trop': PT_TROP, 'SHext': PT_SHEXT}
BB_FRAC = {'NHext': BB_NHEXT_FRACTION, 'Trop': BB_TROP_FRACTION, 'SHext': BB_SHEXT_FRACTION}

# Per-box weights for bounded LS
W_BOX = {
    'NHext': np.diag([100.0, 1.0, 0.5]),
    'Trop':  np.diag([150.0, 1.0, 0.5]),
    'SHext': np.diag([200.0, 1.0, 0.5]),
}


def run(cfg: ModelConfig):
    print("=" * 70)
    print("3×3 THREE-BOX MODEL (Simultaneous δ¹³C + δD)")
    print(f"  KIE: {cfg.kie_mode} | Lifetime: {cfg.lifetime_mode} | N={cfg.n_iterations}")
    print("=" * 70)

    data = load_data(BASE_DIR, three_box=True)
    n = data.n_years
    years = data.model_years
    rng = np.random.default_rng(cfg.seed)

    CH4 = {'NHext': data.CH4_NHext, 'Trop': data.CH4_Trop, 'SHext': data.CH4_SHext}
    c13 = {'NHext': data.c13_NHext, 'Trop': data.c13_Trop, 'SHext': data.c13_SHext}
    dD_atm = {'NHext': data.dD_NHext, 'Trop': data.dD_Trop, 'SHext': data.dD_SHext}

    c13_glob = data.c13_global
    tau_global = compute_lifetime(years, cfg.lifetime_mode, cfg.tau_fixed)
    tau = {b: tau_global * LT_RATIO[b] for b in BOXES}

    NI = cfg.n_iterations

    # Storage
    BB = {b: np.zeros((n, NI)) for b in BOXES}
    FF = {b: np.zeros((n, NI)) for b in BOXES}
    Mic = {b: np.zeros((n, NI)) for b in BOXES}
    qm = {b: QualityMonitor(n, NI, b) for b in BOXES}

    for k in range(NI):
        if (k + 1) % 200 == 0:
            print(f"  iter {k+1}/{NI}")

        tau_NT = max(0.3, rng.normal(TAU_EX_NT_MEAN, TAU_EX_NT_STD))
        tau_TS = max(0.3, rng.normal(TAU_EX_TS_MEAN, TAU_EX_TS_STD))

        kies = sample_KIE(rng, cfg.kie_mode)
        KIE = {}
        for b in BOXES:
            K13, KD = compute_bulk_KIE(kies, SINK_FRACS[b])
            KIE[b] = {'a13': 1.0 / K13, 'aD': 1.0 / KD}

        # Total source per box
        S = {b: np.zeros(n) for b in BOXES}
        for i in range(n):
            M = {b: CH4[b][i] * PT_BOX[b] for b in BOXES}
            M1 = {b: CH4[b][i+1] * PT_BOX[b] for b in BOXES}

            # Concentration-driven exchange
            c_nhext = M['NHext'] / PT_NHEXT
            c_trop  = M['Trop'] / PT_TROP
            c_shext = M['SHext'] / PT_SHEXT

            ex_NT_n = (c_trop - c_nhext) * PT_NHEXT / tau_NT
            ex_NT_t = (c_nhext - c_trop) * PT_TROP / tau_NT
            ex_TS_t = (c_shext - c_trop) * PT_TROP / tau_TS
            ex_TS_s = (c_trop - c_shext) * PT_SHEXT / tau_TS

            S['NHext'][i] = (M1['NHext'] - M['NHext']) + M['NHext'] / tau['NHext'][i] - ex_NT_n
            S['Trop'][i]  = (M1['Trop'] - M['Trop'])   + M['Trop'] / tau['Trop'][i]   - ex_NT_t - ex_TS_t
            S['SHext'][i] = (M1['SHext'] - M['SHext']) + M['SHext'] / tau['SHext'][i] - ex_TS_s

        # Atmospheric obs (MC)
        d13C_glob_MC = sample_atm_d13C(data, k, n)
        nc = min(len(c13_glob), n + 1)
        d13C_off = d13C_glob_MC[:nc] - c13_glob[:nc]

        d13C_MC_box = {}
        for b in BOXES:
            base = c13[b][:nc] if len(c13[b]) >= nc else pad_to_length(c13[b], nc)
            d13C_MC_box[b] = base + d13C_off

        dD_glob_MC = sample_atm_dD(data, k, n)
        dD_glob_mean = data.dD_global[:n+1] if len(data.dD_global) >= n+1 else pad_to_length(data.dD_global, n+1)
        dD_offset = dD_glob_MC[:n+1] - dD_glob_mean[:n+1]

        dD_MC_box = {}
        for b in BOXES:
            base = dD_atm[b][:n+1] if dD_atm[b] is not None and len(dD_atm[b]) >= n+1 else pad_to_length(dD_glob_MC, n+1)
            dD_MC_box[b] = base + dD_offset[:len(base)]

        # Isotopic source fractions per box
        d13C_src = {b: np.zeros(n) for b in BOXES}
        dD_src   = {b: np.zeros(n) for b in BOXES}

        for j in range(n):
            for b in BOXES:
                pt = PT_BOX[b]
                tau_b = tau[b][j]
                a13 = KIE[b]['a13']
                aD  = KIE[b]['aD']

                f13 = delta_to_fraction_d13C(d13C_MC_box[b])
                fD  = delta_to_fraction_dD(dD_MC_box[b])

                n13_j  = f13[j]   * CH4[b][j]   * pt
                n13_j1 = f13[j+1] * CH4[b][j+1] * pt
                nD_j  = fD[j]   * CH4[b][j]   * pt
                nD_j1 = fD[j+1] * CH4[b][j+1] * pt

                # Isotopologue exchange
                if b == 'NHext':
                    f13_t = delta_to_fraction_d13C(d13C_MC_box['Trop'])
                    fD_t  = delta_to_fraction_dD(dD_MC_box['Trop'])
                    ex13 = (f13_t[j] * CH4['Trop'][j] - f13[j] * CH4[b][j]) * pt / tau_NT
                    exD  = (fD_t[j]  * CH4['Trop'][j] - fD[j]  * CH4[b][j]) * pt / tau_NT
                elif b == 'SHext':
                    f13_t = delta_to_fraction_d13C(d13C_MC_box['Trop'])
                    fD_t  = delta_to_fraction_dD(dD_MC_box['Trop'])
                    ex13 = (f13_t[j] * CH4['Trop'][j] - f13[j] * CH4[b][j]) * pt / tau_TS
                    exD  = (fD_t[j]  * CH4['Trop'][j] - fD[j]  * CH4[b][j]) * pt / tau_TS
                else:  # Trop — receives from both
                    f13_n = delta_to_fraction_d13C(d13C_MC_box['NHext'])
                    fD_n  = delta_to_fraction_dD(dD_MC_box['NHext'])
                    f13_s = delta_to_fraction_d13C(d13C_MC_box['SHext'])
                    fD_s  = delta_to_fraction_dD(dD_MC_box['SHext'])
                    ex13 = (f13_n[j] * CH4['NHext'][j] - f13[j] * CH4[b][j]) * pt / tau_NT + \
                           (f13_s[j] * CH4['SHext'][j] - f13[j] * CH4[b][j]) * pt / tau_TS
                    exD  = (fD_n[j]  * CH4['NHext'][j] - fD[j]  * CH4[b][j]) * pt / tau_NT + \
                           (fD_s[j]  * CH4['SHext'][j] - fD[j]  * CH4[b][j]) * pt / tau_TS

                d13C_src[b][j] = (n13_j1 - n13_j + n13_j * a13 / tau_b - ex13) / S[b][j]
                dD_src[b][j]   = (nD_j1  - nD_j  + nD_j  * aD  / tau_b - exD)  / S[b][j]

        sigs = sample_source_signatures_three_box(rng, data, k, n)
        f13_bb  = delta_to_fraction_d13C(sigs['bb_d13C'])
        f13_ff  = delta_to_fraction_d13C(sigs['ff_d13C'])
        f13_mic = delta_to_fraction_d13C(sigs['mic_d13C'])

        for j in range(n):
            for b in BOXES:
                fD_bb  = delta_to_fraction_dD(sigs[f'bb_dD_{b}'])
                fD_ff  = delta_to_fraction_dD(sigs[f'ff_dD_{b}'])
                fD_mic = delta_to_fraction_dD(sigs[f'mic_dD_{b}'])

                A = np.array([
                    [1.0,        1.0,        1.0],
                    [f13_bb[j],  f13_ff[j],  f13_mic[j]],
                    [fD_bb[j],   fD_ff[j],   fD_mic[j]],
                ])

                B = np.array([S[b][j],
                              S[b][j] * d13C_src[b][j],
                              S[b][j] * dD_src[b][j]])

                A_w = W_BOX[b] @ A
                B_w = W_BOX[b] @ B
                ub = S[b][j] * 1.5

                try:
                    res = lsq_linear(A_w, B_w, bounds=(0, max(ub, 1.0)))
                    x = res.x
                    qm[b].record(j, k, A, x)
                except Exception:
                    x = np.array([np.nan, np.nan, np.nan])
                    qm[b].is_nan[j, k] = True
                    qm[b].is_nonphysical[j, k] = True

                BB[b][j, k] = x[0]; FF[b][j, k] = x[1]; Mic[b][j, k] = x[2]

    print("\nMC complete!")
    print("Quality:")
    for b in BOXES:
        qm[b].summary()

    # Global
    BB_G = sum(BB[b] for b in BOXES)
    FF_G = sum(FF[b] for b in BOXES)
    Mic_G = sum(Mic[b] for b in BOXES)

    BB_Gs = smooth_5yr(BB_G); FF_Gs = smooth_5yr(FF_G); Mic_Gs = smooth_5yr(Mic_G)

    print(f"\n{'='*60}")
    print("RESULTS (Global, Smoothed)")
    for lbl, arr in [('BB', BB_Gs), ('FF', FF_Gs), ('Mic', Mic_Gs)]:
        print(f"  {lbl}: {np.nanmean(arr):.1f} ± {np.nanstd(np.nanmean(arr,0)):.1f} Tg/yr")

    print(f"\n  --- TRENDS ---")
    for lbl, arr in [('FF', FF_Gs), ('Mic', Mic_Gs), ('BB', BB_Gs)]:
        d, p = trend_change(arr, years)
        print(f"  {lbl}: Δ={d.mean():+.1f}±{d.std():.1f} ({p:.0f}% pos)")

    # Per-box results
    for b in BOXES:
        bb_s = smooth_5yr(BB[b]); ff_s = smooth_5yr(FF[b]); mic_s = smooth_5yr(Mic[b])
        print(f"\n  [{b}] FF: {np.nanmean(ff_s):.1f} ± {np.nanstd(np.nanmean(ff_s,0)):.1f}, "
              f"Mic: {np.nanmean(mic_s):.1f} ± {np.nanstd(np.nanmean(mic_s,0)):.1f}")

    # Save
    result_df = pd.DataFrame({'Year': years})
    for lbl, arr in [('BB_G', BB_Gs), ('FF_G', FF_Gs), ('Mic_G', Mic_Gs)]:
        result_df[f'{lbl}_mean'] = np.nanmean(arr, 1)
        result_df[f'{lbl}_std'] = np.nanstd(arr, 1)
    for b in BOXES:
        for src_name, src_arr in [('BB', BB[b]), ('FF', FF[b]), ('Mic', Mic[b])]:
            s = smooth_5yr(src_arr)
            result_df[f'{src_name}_{b}_mean'] = np.nanmean(s, 1)
            result_df[f'{src_name}_{b}_std'] = np.nanstd(s, 1)
    result_df.to_csv(OUT_DIR / 'results.csv', index=False)

    with open(OUT_DIR / 'config.json', 'w') as f:
        qr = {b: qm[b].summary() for b in BOXES}
        qr['kie_mode'] = cfg.kie_mode
        qr['lifetime_mode'] = cfg.lifetime_mode
        qr['n_iterations'] = NI
        json.dump(qr, f, indent=2)

    # Plot
    fig, axes = plt.subplots(3, 3, figsize=(16, 12), dpi=150, sharex=True)
    fig.suptitle(f'3×3 Three-Box (KIE={cfg.kie_mode}, τ={cfg.lifetime_mode})', fontsize=13)
    colors = {'BB': 'red', 'FF': 'blue', 'Mic': 'green'}

    for col, b in enumerate(BOXES):
        for row, (name, data_dict) in enumerate([('BB', BB), ('FF', FF), ('Mic', Mic)]):
            ax = axes[row, col]
            arr = smooth_5yr(data_dict[b])
            m = np.nanmean(arr, 1); s = np.nanstd(arr, 1)
            ax.plot(years, m, '-', lw=2.5, color=colors[name])
            ax.fill_between(years, m - s, m + s, alpha=0.3, color=colors[name])
            ax.set_ylabel(f'{name} (Tg/yr)')
            if row == 0: ax.set_title(f'{b} (5yr smoothed)')
            if row == 2: ax.set_xlabel('Year')
            ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(OUT_DIR / 'results.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n  Saved to {OUT_DIR}/")
    print("=" * 70)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="3×3 three-box methane model")
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
