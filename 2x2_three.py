#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2×2_three.py — BB-Fixed, Separate δ¹³C & δD, Three-Box (NHext/Trop/SHext)
============================================================================

Three-box extension of 2×2_two.py.  Each box (NHext >30°N, Trop 30°S–30°N,
SHext <30°S) is solved independently for FF and Mic using each isotope,
with BB fixed from CarbonTracker × latitude fraction.

Inter-box exchange:
  NHext ↔ Trop:  τ_NT ~ N(0.8, 0.1) yr
  Trop ↔ SHext:  τ_TS ~ N(1.2, 0.1) yr
  No direct NHext ↔ SHext exchange.

δD atmospheric:  ThreeBox_atm_dD_annual.csv (station-level, 2005–2024)
δD source sigs:  per-box MC matrices (1998–2021, 1000 iterations)
δ¹³C:            real NH/SH split → interpolated to 3 boxes
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
    sample_source_signatures_three_box, sample_atm_d13C, sample_atm_dD,
    smooth_5yr, trend_change, pad_to_length,
    SINK_FRACTIONS_NHEXT, SINK_FRACTIONS_TROP, SINK_FRACTIONS_SHEXT,
    PT_NHEXT, PT_TROP, PT_SHEXT,
    LIFETIME_RATIO_NHEXT, LIFETIME_RATIO_TROP, LIFETIME_RATIO_SHEXT,
    BB_NHEXT_FRACTION, BB_TROP_FRACTION, BB_SHEXT_FRACTION,
    TAU_EX_NT_MEAN, TAU_EX_NT_STD, TAU_EX_TS_MEAN, TAU_EX_TS_STD,
    DD_IH_OFFSET,
)

BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "Output_2x2_three"
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


def run(cfg: ModelConfig):
    print("=" * 70)
    print("2×2 THREE-BOX MODEL (BB-Fixed, Separate δ¹³C & δD)")
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
    BB_fixed = {b: data.BB_annual * BB_FRAC[b] for b in BOXES}

    NI = cfg.n_iterations
    # Storage: FF and Mic per box per isotope
    FF = {b: {'c': np.zeros((n, NI)), 'd': np.zeros((n, NI))} for b in BOXES}
    Mic = {b: {'c': np.zeros((n, NI)), 'd': np.zeros((n, NI))} for b in BOXES}
    n_neg = {'c': 0, 'd': 0}
    total = 0

    for k in range(NI):
        if (k + 1) % 200 == 0:
            print(f"  iter {k+1}/{NI}")

        # Sample exchange times
        tau_NT = max(0.3, rng.normal(TAU_EX_NT_MEAN, TAU_EX_NT_STD))
        tau_TS = max(0.3, rng.normal(TAU_EX_TS_MEAN, TAU_EX_TS_STD))

        # KIE per box
        kies = sample_KIE(rng, cfg.kie_mode)
        KIE = {}
        for b in BOXES:
            K13, KD = compute_bulk_KIE(kies, SINK_FRACS[b])
            KIE[b] = {'a13': 1.0 / K13, 'aD': 1.0 / KD}

        # Total source per box (mass balance with inter-box exchange)
        S = {b: np.zeros(n) for b in BOXES}
        for i in range(n):
            M = {b: CH4[b][i] * PT_BOX[b] for b in BOXES}
            M1 = {b: CH4[b][i+1] * PT_BOX[b] for b in BOXES}

            # Exchange fluxes (positive = into box)
            ex_NT_n = (M['Trop'] / PT_TROP - M['NHext'] / PT_NHEXT) * PT_NHEXT / tau_NT
            ex_NT_t = -ex_NT_n * PT_NHEXT / PT_TROP  # conservation
            ex_TS_t = (M['SHext'] / PT_SHEXT - M['Trop'] / PT_TROP) * PT_TROP / tau_TS
            ex_TS_s = -ex_TS_t * PT_TROP / PT_SHEXT

            S['NHext'][i] = (M1['NHext'] - M['NHext']) + M['NHext'] / tau['NHext'][i] - ex_NT_n
            S['Trop'][i]  = (M1['Trop'] - M['Trop'])   + M['Trop'] / tau['Trop'][i]   - ex_NT_t - ex_TS_t
            S['SHext'][i] = (M1['SHext'] - M['SHext']) + M['SHext'] / tau['SHext'][i] - ex_TS_s

        # Atmospheric obs (MC-sampled)
        d13C_glob_MC = sample_atm_d13C(data, k, n)
        nc = min(len(c13_glob), n + 1)
        d13C_off = d13C_glob_MC[:nc] - c13_glob[:nc]

        # δ¹³C per box: real mean + MC offset
        d13C_MC_box = {}
        for b in BOXES:
            base = c13[b][:nc] if len(c13[b]) >= nc else pad_to_length(c13[b], nc)
            d13C_MC_box[b] = base + d13C_off

        # δD per box: use annual observations + small MC perturbation from global
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
            for bi, b in enumerate(BOXES):
                pt = PT_BOX[b]
                tau_b = tau[b][j]
                a13 = KIE[b]['a13']
                aD  = KIE[b]['aD']

                f13 = delta_to_fraction_d13C(d13C_MC_box[b])
                fD  = delta_to_fraction_dD(dD_MC_box[b])

                # ¹³C
                n13_j  = f13[j]   * CH4[b][j]   * pt
                n13_j1 = f13[j+1] * CH4[b][j+1] * pt

                # D
                nD_j  = fD[j]   * CH4[b][j]   * pt
                nD_j1 = fD[j+1] * CH4[b][j+1] * pt

                # Exchange: need neighbor isotope ratios
                if b == 'NHext':
                    f13_nb = delta_to_fraction_d13C(d13C_MC_box['Trop'])
                    fD_nb  = delta_to_fraction_dD(dD_MC_box['Trop'])
                    n13_nb = f13_nb[j] * CH4['Trop'][j] * PT_TROP
                    nD_nb  = fD_nb[j]  * CH4['Trop'][j] * PT_TROP
                    # Exchange flux of isotopologue
                    ex13 = (n13_nb / PT_TROP - n13_j / pt) * pt / tau_NT
                    exD  = (nD_nb / PT_TROP  - nD_j / pt)  * pt / tau_NT
                elif b == 'Trop':
                    # NH side
                    f13_nh = delta_to_fraction_d13C(d13C_MC_box['NHext'])
                    fD_nh  = delta_to_fraction_dD(dD_MC_box['NHext'])
                    n13_nh = f13_nh[j] * CH4['NHext'][j] * PT_NHEXT
                    nD_nh  = fD_nh[j]  * CH4['NHext'][j] * PT_NHEXT
                    # SH side
                    f13_sh = delta_to_fraction_d13C(d13C_MC_box['SHext'])
                    fD_sh  = delta_to_fraction_dD(dD_MC_box['SHext'])
                    n13_sh = f13_sh[j] * CH4['SHext'][j] * PT_SHEXT
                    nD_sh  = fD_sh[j]  * CH4['SHext'][j] * PT_SHEXT
                    # NH exchange (into Trop)
                    ex13_nt = -(n13_nb_dummy := n13_j / pt - n13_nh / PT_NHEXT) * PT_NHEXT / tau_NT * (-1)
                    # Simpler: just compute concentration-driven exchange
                    ex13 = (n13_nh / PT_NHEXT - n13_j / pt) * pt / tau_NT + \
                           (n13_sh / PT_SHEXT - n13_j / pt) * pt / tau_TS
                    exD  = (nD_nh / PT_NHEXT - nD_j / pt) * pt / tau_NT + \
                           (nD_sh / PT_SHEXT - nD_j / pt) * pt / tau_TS
                else:  # SHext
                    f13_nb = delta_to_fraction_d13C(d13C_MC_box['Trop'])
                    fD_nb  = delta_to_fraction_dD(dD_MC_box['Trop'])
                    n13_nb = f13_nb[j] * CH4['Trop'][j] * PT_TROP
                    nD_nb  = fD_nb[j]  * CH4['Trop'][j] * PT_TROP
                    ex13 = (n13_nb / PT_TROP - n13_j / pt) * pt / tau_TS
                    exD  = (nD_nb / PT_TROP  - nD_j / pt)  * pt / tau_TS

                d13C_src[b][j] = (n13_j1 - n13_j + n13_j * a13 / tau_b - ex13) / S[b][j]
                dD_src[b][j]   = (nD_j1  - nD_j  + nD_j  * aD  / tau_b - exD)  / S[b][j]

        sigs = sample_source_signatures_three_box(rng, data, k, n)

        for j in range(n):
            total += 1
            for b in BOXES:
                d13c_delta = fraction_to_delta_d13C(d13C_src[b][j])
                dD_delta   = fraction_to_delta_dD(dD_src[b][j])

                BB_j = BB_fixed[b][j] if j < len(BB_fixed[b]) else data.BB_global_mean * BB_FRAC[b]

                denom_c = sigs['ff_d13C'][j] - sigs['mic_d13C'][j]
                denom_d = sigs[f'ff_dD_{b}'][j] - sigs[f'mic_dD_{b}'][j]

                # δ¹³C solve
                if abs(denom_c) > 0.1:
                    ff = (S[b][j] * d13c_delta - sigs['mic_d13C'][j] * (S[b][j] - BB_j) - sigs['bb_d13C'][j] * BB_j) / denom_c
                    mic = S[b][j] - BB_j - ff
                    if ff < 0 or mic < 0:
                        n_neg['c'] += 1
                        ff = max(0, ff); mic = S[b][j] - BB_j - ff
                        if mic < 0: mic = 0; ff = S[b][j] - BB_j
                else:
                    ff = np.nan; mic = np.nan
                FF[b]['c'][j, k] = ff; Mic[b]['c'][j, k] = mic

                # δD solve
                if abs(denom_d) > 1.0:
                    ff = (S[b][j] * dD_delta - sigs[f'mic_dD_{b}'][j] * (S[b][j] - BB_j) - sigs[f'bb_dD_{b}'][j] * BB_j) / denom_d
                    mic = S[b][j] - BB_j - ff
                    if ff < 0 or mic < 0:
                        n_neg['d'] += 1
                        ff = max(0, ff); mic = S[b][j] - BB_j - ff
                        if mic < 0: mic = 0; ff = S[b][j] - BB_j
                else:
                    ff = np.nan; mic = np.nan
                FF[b]['d'][j, k] = ff; Mic[b]['d'][j, k] = mic

    print(f"\n  Negatives: δ¹³C={n_neg['c']}, δD={n_neg['d']}")

    # Global = sum of 3 boxes
    FF_G_c = sum(FF[b]['c'] for b in BOXES)
    Mic_G_c = sum(Mic[b]['c'] for b in BOXES)
    FF_G_d = sum(FF[b]['d'] for b in BOXES)
    Mic_G_d = sum(Mic[b]['d'] for b in BOXES)

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
    result_df = pd.DataFrame({'Year': years})
    for lbl, arr in [('FF_G_d13C', FF_G_cs), ('Mic_G_d13C', Mic_G_cs),
                     ('FF_G_dD', FF_G_ds), ('Mic_G_dD', Mic_G_ds)]:
        result_df[f'{lbl}_mean'] = np.nanmean(arr, 1)
        result_df[f'{lbl}_std'] = np.nanstd(arr, 1)
    for b in BOXES:
        for iso in ['c', 'd']:
            iso_label = 'd13C' if iso == 'c' else 'dD'
            result_df[f'FF_{b}_{iso_label}_mean'] = np.nanmean(smooth_5yr(FF[b][iso]), 1)
            result_df[f'Mic_{b}_{iso_label}_mean'] = np.nanmean(smooth_5yr(Mic[b][iso]), 1)
    result_df['BB_fixed'] = data.BB_annual[:n]
    result_df.to_csv(OUT_DIR / 'results.csv', index=False)

    with open(OUT_DIR / 'config.json', 'w') as f:
        json.dump({'kie_mode': cfg.kie_mode, 'lifetime_mode': cfg.lifetime_mode,
                   'n_iterations': NI}, f, indent=2)

    # Plot
    fig, axes = plt.subplots(2, 3, figsize=(16, 8), dpi=150, sharex=True)
    fig.suptitle(f'2×2 Three-Box (KIE={cfg.kie_mode}, τ={cfg.lifetime_mode})', fontsize=13)
    colors = {'d13C': '#2166ac', 'dD': '#b2182b'}

    for col, b in enumerate(BOXES):
        for row, (lbl, data_dict) in enumerate([('FF', FF), ('Mic', Mic)]):
            ax = axes[row, col]
            for iso, iso_lbl in [('c', 'δ¹³C'), ('d', 'δD')]:
                arr = smooth_5yr(data_dict[b][iso])
                m = np.nanmean(arr, 1); s = 2 * np.nanstd(arr, 1)
                clr = colors['d13C'] if iso == 'c' else colors['dD']
                ax.plot(years, m, '-', lw=2, color=clr, label=iso_lbl)
                ax.fill_between(years, m - s, m + s, alpha=0.2, color=clr)
            ax.set_title(f'{lbl} — {b}')
            ax.set_ylabel('Tg/yr')
            if row == 1: ax.set_xlabel('Year')
            ax.legend(fontsize=7); ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(OUT_DIR / 'results.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n  Saved to {OUT_DIR}/")
    print("=" * 70)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="2×2 three-box methane model")
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
