#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 2: Dual-Isotope Joint WLS Inversion — KIE Sensitivity
=============================================================
Uses both δ¹³C and δD simultaneously in a weighted least-squares
3×2 over-determined system to partition FF and Mic.
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import lsq_linear
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from common import (
    load_data, sample_KIE, compute_bulk_KIE, compute_lifetime,
    delta_to_fraction_d13C, delta_to_fraction_dD,
    fraction_to_delta_d13C, fraction_to_delta_dD,
    sample_source_signatures, sample_atm_d13C, sample_atm_dD,
    smooth_5yr, trend_change, pad_to_length,
    SINK_FRACTIONS_GLOBAL, PT,
)

OUT_DIR = Path(__file__).resolve().parent / "results" / "phase2_dual_isotope"
FIG_DIR = Path(__file__).resolve().parent / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

N_ITER = 1000
SEED = 42

# Weights for WLS (inverse variance)
W_MASS = 1.0
W_C13 = 1.0 / (1.0**2)   # σ(δ¹³C) ≈ 1‰
W_D = 1.0 / (5.0**2)     # σ(δD) ≈ 5‰


def run_dual_isotope(oh13c_mode: str, label: str):
    """Run joint dual-isotope WLS with a specific OH-¹³C KIE setting."""
    print(f"\n{'='*60}")
    print(f"  Phase 2 — Dual isotope | OH_13C = {oh13c_mode} ({label})")
    print(f"{'='*60}")

    data = load_data(REPO_ROOT, two_box=False)
    n = data.n_years
    years = data.model_years
    CH4 = data.CH4_global
    rng = np.random.default_rng(SEED)

    tau = compute_lifetime(years, 'varying')

    SumSource = np.zeros(n)
    for i in range(n):
        SumSource[i] = CH4[i + 1] * PT - CH4[i] * PT + CH4[i] * PT / tau[i]

    BB = data.BB_annual.copy()

    FF_all = np.zeros((n, N_ITER))
    Mic_all = np.zeros((n, N_ITER))
    residuals = np.zeros((n, N_ITER))
    n_neg = 0
    total = 0

    for k in range(N_ITER):
        if (k + 1) % 200 == 0:
            print(f"    iter {k+1}/{N_ITER}")

        kies = sample_KIE(rng, 'sampled')
        if oh13c_mode == 'saueressig':
            kies['OH_13C'] = 1.0039
        elif oh13c_mode == 'cantrell':
            kies['OH_13C'] = 1.0054

        KIE_13C, KIE_D = compute_bulk_KIE(kies, SINK_FRACTIONS_GLOBAL)
        alpha_13C = 1.0 / KIE_13C
        alpha_D = 1.0 / KIE_D

        d13C_atm = sample_atm_d13C(data, k, n)
        dD_atm = sample_atm_dD(data, k, n)

        f13 = delta_to_fraction_d13C(d13C_atm)
        fD = delta_to_fraction_dD(dD_atm)
        n13 = f13 * CH4[:n+1] * PT
        nD = fD * CH4[:n+1] * PT

        sigs = sample_source_signatures(rng, data, k, n)

        for j in range(n):
            total += 1
            S = SumSource[j]
            BB_j = BB[j] if j < len(BB) else data.BB_global_mean

            # Compute isotopic source deltas
            d13C_src_f = (n13[j+1] - n13[j] + n13[j] * alpha_13C / tau[j]) / S
            dD_src_f = (nD[j+1] - nD[j] + nD[j] * alpha_D / tau[j]) / S
            d13C_src = fraction_to_delta_d13C(d13C_src_f)
            dD_src = fraction_to_delta_dD(dD_src_f)

            # Build 3×2 WLS system
            # A·x = b where x = [FF, Mic]
            A = np.array([
                [1.0, 1.0],
                [sigs['ff_d13C'][j], sigs['mic_d13C'][j]],
                [sigs['ff_dD'][j], sigs['mic_dD'][j]],
            ])
            b = np.array([
                S - BB_j,
                S * d13C_src - BB_j * sigs['bb_d13C'][j],
                S * dD_src - BB_j * sigs['bb_dD'][j],
            ])

            # Apply weights
            w = np.sqrt(np.array([W_MASS, W_C13, W_D]))
            Aw = A * w[:, None]
            bw = b * w

            # Solve bounded least squares
            try:
                result = lsq_linear(Aw, bw, bounds=(0, np.inf))
                ff, mic = result.x
                residuals[j, k] = result.cost
            except Exception:
                ff, mic = np.nan, np.nan

            if not np.isnan(ff) and (ff < 0 or mic < 0):
                n_neg += 1

            FF_all[j, k] = ff
            Mic_all[j, k] = mic

    pct_neg = 100 * n_neg / max(1, total)
    print(f"  Negatives: {n_neg} ({pct_neg:.1f}%)")

    FF_s = smooth_5yr(FF_all)
    Mic_s = smooth_5yr(Mic_all)

    delta_ff, pct_ff = trend_change(FF_s, years)
    delta_mic, pct_mic = trend_change(Mic_s, years)

    print(f"  FF trend:  Δ={delta_ff.mean():+.1f} ± {delta_ff.std():.1f} Tg/yr ({pct_ff:.0f}% positive)")
    print(f"  Mic trend: Δ={delta_mic.mean():+.1f} ± {delta_mic.std():.1f} Tg/yr ({pct_mic:.0f}% positive)")

    np.savez(OUT_DIR / f"run_{label}.npz",
             FF=FF_s, Mic=Mic_s, years=years,
             delta_ff=delta_ff, delta_mic=delta_mic,
             mean_residual=float(np.nanmean(residuals)))

    return {
        'label': label,
        'oh13c_mode': oh13c_mode,
        'FF_trend_mean': float(delta_ff.mean()),
        'FF_trend_std': float(delta_ff.std()),
        'FF_pct_positive': float(pct_ff),
        'Mic_trend_mean': float(delta_mic.mean()),
        'Mic_trend_std': float(delta_mic.std()),
        'Mic_pct_positive': float(pct_mic),
        'neg_pct': round(pct_neg, 2),
        'mean_residual': float(np.nanmean(residuals)),
    }


def make_figure(summaries):
    """Histogram of FF/Mic trends for the 3 KIE settings (dual isotope)."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=300)
    colors = {'A_saueressig': 'tab:blue', 'B_cantrell': 'tab:red', 'C_sampled': 'gray'}
    labels_map = {'A_saueressig': 'Saueressig (1.0039)',
                  'B_cantrell': 'Cantrell (1.0054)',
                  'C_sampled': 'Sampled [1.0039–1.0054]'}

    for label in ['A_saueressig', 'B_cantrell', 'C_sampled']:
        dat = np.load(OUT_DIR / f"run_{label}.npz")
        for ax_idx, (key, title) in enumerate([('delta_ff', 'Fossil Fuel'), ('delta_mic', 'Microbial')]):
            arr = dat[key]
            m, s = arr.mean(), arr.std()
            axes[ax_idx].hist(arr, bins=50, alpha=0.5, color=colors[label],
                              label=f"{labels_map[label]}: {m:+.1f}±{s:.1f}", density=True)

    for ax, title in zip(axes, ['Fossil Fuel', 'Microbial']):
        ax.axvline(0, color='black', ls='--', lw=1)
        ax.set_xlabel('Δ Emissions (Tg/yr): 2020–2022 vs 2005–2007')
        ax.set_ylabel('Density')
        ax.set_title(f'{title} — Dual Isotope (δ¹³C + δD)')
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(FIG_DIR / "phase2_dual_isotope_trends.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"\n  Figure saved: {FIG_DIR / 'phase2_dual_isotope_trends.png'}")


if __name__ == "__main__":
    results = []
    results.append(run_dual_isotope('saueressig', 'A_saueressig'))
    results.append(run_dual_isotope('cantrell', 'B_cantrell'))
    results.append(run_dual_isotope('sampled', 'C_sampled'))

    with open(OUT_DIR / "summary.json", 'w') as f:
        json.dump(results, f, indent=2)

    make_figure(results)
    print("\n✓ Phase 2 complete.")
