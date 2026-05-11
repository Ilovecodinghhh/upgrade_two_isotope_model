#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 5: Weight & Cl Fraction Sensitivity Sweep
=================================================
Tests:
1. How does the optimal δD weight in WLS affect KSR?
2. How does the Cl sink fraction (0.6% vs 3.5%) affect KIE sensitivity?

From MASTER_DATA_INVENTORY critical uncertainties:
- Cl fraction ranges 0.6% (Thanwerdas) to 3.5% (default) — factor of 6
- Cl has the largest δD KIE (α=1.52), so this massively affects δD budget
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
    PT,
)

BASE = Path(__file__).resolve().parent
OUT_DIR = BASE / "results" / "phase5_weight_Cl_sweep"
FIG_DIR = BASE / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

N_ITER = 500  # smaller for sweep (many configs)
SEED = 42

# Cl fraction configurations from MASTER_DATA_INVENTORY
CL_CONFIGS = {
    'thanwerdas_low': {'OH': 0.899, 'Cl': 0.006, 'Strat': 0.030, 'Soil': 0.065},
    'default':        {'OH': 0.835, 'Cl': 0.035, 'Strat': 0.070, 'Soil': 0.060},
    'high_Cl':        {'OH': 0.805, 'Cl': 0.065, 'Strat': 0.070, 'Soil': 0.060},
}

# δD weight sweep: from 0 (pure δ¹³C) to 1.0 (equal weight with δ¹³C)
DD_WEIGHTS = [0.0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0]


def run_single_config(sink_fracs: dict, w_dD: float, oh13c_mode: str):
    """Run 1-box dual-isotope with specified Cl fraction and δD weight."""
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

    for k in range(N_ITER):
        kies = sample_KIE(rng, 'sampled')
        if oh13c_mode == 'saueressig':
            kies['OH_13C'] = 1.0039
        elif oh13c_mode == 'cantrell':
            kies['OH_13C'] = 1.0054

        KIE_13C, KIE_D = compute_bulk_KIE(kies, sink_fracs)
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
            S = SumSource[j]
            BB_j = BB[j] if j < len(BB) else data.BB_global_mean

            d13C_src_f = (n13[j+1] - n13[j] + n13[j] * alpha_13C / tau[j]) / S
            dD_src_f = (nD[j+1] - nD[j] + nD[j] * alpha_D / tau[j]) / S
            d13C_src = fraction_to_delta_d13C(d13C_src_f)
            dD_src = fraction_to_delta_dD(dD_src_f)

            if w_dD == 0:
                # Pure δ¹³C analytic
                denom = sigs['ff_d13C'][j] - sigs['mic_d13C'][j]
                if abs(denom) > 0.1:
                    ff = (S * d13C_src - sigs['mic_d13C'][j] * (S - BB_j) -
                          sigs['bb_d13C'][j] * BB_j) / denom
                    mic = S - BB_j - ff
                else:
                    ff, mic = S / 3, S / 3
                ff = max(0, ff)
                mic = max(0, S - BB_j - ff)
            else:
                # WLS with specified δD weight
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
                W_C13 = 1.0 / (1.0**2)
                W_D_eff = w_dD / (5.0**2)
                w_vec = np.sqrt(np.array([1.0, W_C13, W_D_eff]))
                Aw = A * w_vec[:, None]
                bw = b * w_vec
                try:
                    result = lsq_linear(Aw, bw, bounds=(0, np.inf))
                    ff, mic = result.x
                except Exception:
                    ff = max(0, (S - BB_j) / 2)
                    mic = max(0, (S - BB_j) / 2)

            FF_all[j, k] = ff
            Mic_all[j, k] = mic

    FF_s = smooth_5yr(FF_all)
    Mic_s = smooth_5yr(Mic_all)
    delta_ff, _ = trend_change(FF_s, years)
    delta_mic, _ = trend_change(Mic_s, years)
    return delta_ff, delta_mic, np.std(FF_s, axis=1).mean(), np.std(Mic_s, axis=1).mean()


def main():
    print("=" * 60)
    print("  Phase 5 — Weight & Cl Fraction Sensitivity Sweep")
    print("=" * 60)

    all_results = {}

    # === Test 1: δD weight sweep (default Cl fraction) ===
    print("\n--- Test 1: δD Weight Sweep (Cl=3.5%) ---")
    weight_results = []
    for w_dD in DD_WEIGHTS:
        ksr_data = {}
        for oh13c in ['saueressig', 'cantrell']:
            dff, dmic, unc_ff, unc_mic = run_single_config(
                CL_CONFIGS['default'], w_dD, oh13c)
            ksr_data[oh13c] = {
                'ff_mean': float(dff.mean()),
                'mic_mean': float(dmic.mean()),
                'ff_std': float(dff.std()),
                'mic_std': float(dmic.std()),
                'unc_ff': float(unc_ff),
                'unc_mic': float(unc_mic),
            }

        spread_ff = abs(ksr_data['cantrell']['ff_mean'] - ksr_data['saueressig']['ff_mean'])
        spread_mic = abs(ksr_data['cantrell']['mic_mean'] - ksr_data['saueressig']['mic_mean'])
        avg_unc_ff = (ksr_data['cantrell']['ff_std'] + ksr_data['saueressig']['ff_std']) / 2
        avg_unc_mic = (ksr_data['cantrell']['mic_std'] + ksr_data['saueressig']['mic_std']) / 2

        entry = {
            'w_dD': w_dD,
            'spread_ff': spread_ff,
            'spread_mic': spread_mic,
            'avg_unc_ff': avg_unc_ff,
            'avg_unc_mic': avg_unc_mic,
            'ksr_data': ksr_data,
        }
        weight_results.append(entry)
        print(f"  w_dD={w_dD:.2f}: spread_FF={spread_ff:.2f}, unc_FF={avg_unc_ff:.1f} | "
              f"spread_Mic={spread_mic:.2f}, unc_Mic={avg_unc_mic:.1f}")

    # === Test 2: Cl fraction sweep (optimal weight from Test 1) ===
    print("\n--- Test 2: Cl Fraction Sweep ---")
    cl_results = []
    for cl_name, cl_fracs in CL_CONFIGS.items():
        for w_dD in [0.0, 0.1, 0.5]:  # key weights
            ksr_data = {}
            for oh13c in ['saueressig', 'cantrell']:
                dff, dmic, unc_ff, unc_mic = run_single_config(cl_fracs, w_dD, oh13c)
                ksr_data[oh13c] = {
                    'ff_mean': float(dff.mean()),
                    'mic_mean': float(dmic.mean()),
                    'ff_std': float(dff.std()),
                    'mic_std': float(dmic.std()),
                }

            spread_ff = abs(ksr_data['cantrell']['ff_mean'] - ksr_data['saueressig']['ff_mean'])
            spread_mic = abs(ksr_data['cantrell']['mic_mean'] - ksr_data['saueressig']['mic_mean'])

            entry = {
                'cl_config': cl_name,
                'cl_frac': cl_fracs['Cl'],
                'w_dD': w_dD,
                'spread_ff': spread_ff,
                'spread_mic': spread_mic,
            }
            cl_results.append(entry)
            print(f"  Cl={cl_name} ({cl_fracs['Cl']:.3f}), w_dD={w_dD:.1f}: "
                  f"spread_FF={spread_ff:.2f}, spread_Mic={spread_mic:.2f}")

    # Save
    all_results = {
        'weight_sweep': weight_results,
        'cl_sweep': cl_results,
    }
    with open(OUT_DIR / "summary.json", 'w') as f:
        json.dump(all_results, f, indent=2, default=str)

    # === Figures ===

    # Figure 6a: δD weight vs KIE spread and total uncertainty
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=300)

    weights = [r['w_dD'] for r in weight_results]
    spreads_ff = [r['spread_ff'] for r in weight_results]
    spreads_mic = [r['spread_mic'] for r in weight_results]
    uncs_ff = [r['avg_unc_ff'] for r in weight_results]
    uncs_mic = [r['avg_unc_mic'] for r in weight_results]

    ax = axes[0]
    ax.plot(weights, spreads_ff, 'ro-', lw=2, markersize=8, label='KIE Spread (FF)')
    ax.plot(weights, spreads_mic, 'bo-', lw=2, markersize=8, label='KIE Spread (Mic)')
    ax.set_xlabel('δD Weight in WLS')
    ax.set_ylabel('|Cantrell − Saueressig| (Tg/yr)')
    ax.set_title('KIE Sensitivity vs δD Weight')
    ax.legend()
    ax.grid(alpha=0.3)
    ax.set_xscale('symlog', linthresh=0.01)

    ax = axes[1]
    ax.plot(weights, uncs_ff, 'ro-', lw=2, markersize=8, label='Total σ (FF)')
    ax.plot(weights, uncs_mic, 'bo-', lw=2, markersize=8, label='Total σ (Mic)')
    ax.set_xlabel('δD Weight in WLS')
    ax.set_ylabel('Average σ across years (Tg/yr)')
    ax.set_title('Total Uncertainty vs δD Weight')
    ax.legend()
    ax.grid(alpha=0.3)
    ax.set_xscale('symlog', linthresh=0.01)

    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig6_weight_sweep.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"\n  Figure 6 saved: {FIG_DIR / 'fig6_weight_sweep.png'}")

    # Figure 7: Cl fraction × weight interaction
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=300)

    for cl_name, marker, color in [('thanwerdas_low', 's', 'green'),
                                    ('default', 'o', 'blue'),
                                    ('high_Cl', '^', 'red')]:
        subset = [r for r in cl_results if r['cl_config'] == cl_name]
        ws = [r['w_dD'] for r in subset]
        sff = [r['spread_ff'] for r in subset]
        smic = [r['spread_mic'] for r in subset]

        axes[0].plot(ws, sff, f'{marker}-', color=color, lw=2, markersize=8,
                     label=f"Cl={subset[0]['cl_frac']:.3f} ({cl_name})")
        axes[1].plot(ws, smic, f'{marker}-', color=color, lw=2, markersize=8,
                     label=f"Cl={subset[0]['cl_frac']:.3f} ({cl_name})")

    for ax, title in zip(axes, ['Fossil Fuel', 'Microbial']):
        ax.set_xlabel('δD Weight in WLS')
        ax.set_ylabel('KIE Spread (Tg/yr)')
        ax.set_title(f'{title} — Cl Fraction × δD Weight')
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig7_Cl_weight_interaction.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Figure 7 saved: {FIG_DIR / 'fig7_Cl_weight_interaction.png'}")

    print("\n✓ Phase 5 complete.")


if __name__ == "__main__":
    main()
