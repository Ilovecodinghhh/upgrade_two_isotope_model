#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 3: Comparison — KSR Calculation, Statistical Tests, Publication Figures
==============================================================================
Loads Phase 1 and Phase 2 results and produces the key comparison metrics.
"""

import json
from pathlib import Path

import numpy as np
from scipy.stats import ks_2samp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent
P1_DIR = BASE / "results" / "phase1_d13C_only"
P2_DIR = BASE / "results" / "phase2_dual_isotope"
OUT_DIR = BASE / "results" / "phase3_comparison"
FIG_DIR = BASE / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

N_BOOT = 1000
RNG = np.random.default_rng(123)


def load_run(directory, label):
    d = np.load(directory / f"run_{label}.npz")
    return d


def compute_ksr_bootstrap(p1_A, p1_B, p2_A, p2_B, key, n_boot=N_BOOT):
    """Compute KIE Sensitivity Ratio with bootstrap CI."""
    arr_s_A = p1_A[key]
    arr_s_B = p1_B[key]
    arr_d_A = p2_A[key]
    arr_d_B = p2_B[key]

    n = len(arr_s_A)

    # Point estimate
    spread_single = abs(arr_s_B.mean() - arr_s_A.mean())
    spread_dual = abs(arr_d_B.mean() - arr_d_A.mean())
    ksr = spread_single / max(spread_dual, 1e-6)

    # Bootstrap
    ksr_boot = []
    for _ in range(n_boot):
        idx = RNG.choice(n, size=n, replace=True)
        ss = abs(arr_s_B[idx].mean() - arr_s_A[idx].mean())
        sd = abs(arr_d_B[idx].mean() - arr_d_A[idx].mean())
        ksr_boot.append(ss / max(sd, 1e-6))

    ksr_boot = np.array(ksr_boot)
    ci_low = float(np.percentile(ksr_boot, 2.5))
    ci_high = float(np.percentile(ksr_boot, 97.5))

    return float(ksr), ci_low, ci_high


def main():
    print("=" * 60)
    print("  Phase 3 — Comparison & KSR Calculation")
    print("=" * 60)

    # Load runs
    p1_A = load_run(P1_DIR, 'A_saueressig')
    p1_B = load_run(P1_DIR, 'B_cantrell')
    p1_C = load_run(P1_DIR, 'C_sampled')
    p2_A = load_run(P2_DIR, 'A_saueressig')
    p2_B = load_run(P2_DIR, 'B_cantrell')
    p2_C = load_run(P2_DIR, 'C_sampled')

    years = p1_A['years']

    # --- KSR ---
    ksr_ff, ci_ff_lo, ci_ff_hi = compute_ksr_bootstrap(p1_A, p1_B, p2_A, p2_B, 'delta_ff')
    ksr_mic, ci_mic_lo, ci_mic_hi = compute_ksr_bootstrap(p1_A, p1_B, p2_A, p2_B, 'delta_mic')

    print(f"\n  KSR (Fossil Fuel):  {ksr_ff:.2f}  [95% CI: {ci_ff_lo:.2f} – {ci_ff_hi:.2f}]")
    print(f"  KSR (Microbial):    {ksr_mic:.2f}  [95% CI: {ci_mic_lo:.2f} – {ci_mic_hi:.2f}]")

    # --- KS tests ---
    ks_ff_single = ks_2samp(p1_A['delta_ff'], p1_B['delta_ff'])
    ks_ff_dual = ks_2samp(p2_A['delta_ff'], p2_B['delta_ff'])
    ks_mic_single = ks_2samp(p1_A['delta_mic'], p1_B['delta_mic'])
    ks_mic_dual = ks_2samp(p2_A['delta_mic'], p2_B['delta_mic'])

    print(f"\n  KS test (FF,  δ¹³C-only):     D={ks_ff_single.statistic:.3f}, p={ks_ff_single.pvalue:.2e}")
    print(f"  KS test (FF,  dual-isotope):  D={ks_ff_dual.statistic:.3f}, p={ks_ff_dual.pvalue:.2e}")
    print(f"  KS test (Mic, δ¹³C-only):     D={ks_mic_single.statistic:.3f}, p={ks_mic_single.pvalue:.2e}")
    print(f"  KS test (Mic, dual-isotope):  D={ks_mic_dual.statistic:.3f}, p={ks_mic_dual.pvalue:.2e}")

    # --- Uncertainty reduction ---
    unc_red_ff = []
    unc_red_mic = []
    for j in range(len(years)):
        s_single_ff = np.nanstd(p1_C['FF'][:, :][j])  # sampled-KIE run
        s_dual_ff = np.nanstd(p2_C['FF'][:, :][j])
        red_ff = (1 - s_dual_ff / max(s_single_ff, 1e-6)) * 100
        unc_red_ff.append(float(red_ff))

        s_single_mic = np.nanstd(p1_C['Mic'][:, :][j])
        s_dual_mic = np.nanstd(p2_C['Mic'][:, :][j])
        red_mic = (1 - s_dual_mic / max(s_single_mic, 1e-6)) * 100
        unc_red_mic.append(float(red_mic))

    print(f"\n  Mean uncertainty reduction (FF):  {np.mean(unc_red_ff):.1f}%")
    print(f"  Mean uncertainty reduction (Mic): {np.mean(unc_red_mic):.1f}%")

    # --- Save summary ---
    summary = {
        'KSR_FF': {'value': ksr_ff, 'CI_95': [ci_ff_lo, ci_ff_hi]},
        'KSR_Mic': {'value': ksr_mic, 'CI_95': [ci_mic_lo, ci_mic_hi]},
        'uncertainty_reduction_FF_pct': {'mean': float(np.mean(unc_red_ff)), 'by_year': unc_red_ff},
        'uncertainty_reduction_Mic_pct': {'mean': float(np.mean(unc_red_mic)), 'by_year': unc_red_mic},
        'KS_test_FF': {
            'single_D': float(ks_ff_single.statistic), 'single_p': float(ks_ff_single.pvalue),
            'dual_D': float(ks_ff_dual.statistic), 'dual_p': float(ks_ff_dual.pvalue),
        },
        'KS_test_Mic': {
            'single_D': float(ks_mic_single.statistic), 'single_p': float(ks_mic_single.pvalue),
            'dual_D': float(ks_mic_dual.statistic), 'dual_p': float(ks_mic_dual.pvalue),
        },
    }
    with open(OUT_DIR / "summary.json", 'w') as f:
        json.dump(summary, f, indent=2)

    # === FIGURES ===

    # --- Figure 1: 2×2 histogram comparison ---
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), dpi=300)

    panels = [
        (axes[0, 0], 'delta_ff', 'Fossil Fuel — δ¹³C-only', p1_A, p1_B, ks_ff_single),
        (axes[0, 1], 'delta_ff', 'Fossil Fuel — Dual Isotope', p2_A, p2_B, ks_ff_dual),
        (axes[1, 0], 'delta_mic', 'Microbial — δ¹³C-only', p1_A, p1_B, ks_mic_single),
        (axes[1, 1], 'delta_mic', 'Microbial — Dual Isotope', p2_A, p2_B, ks_mic_dual),
    ]
    panel_labels = ['(a)', '(b)', '(c)', '(d)']
    ksr_vals = [ksr_ff, ksr_ff, ksr_mic, ksr_mic]

    for idx, (ax, key, title, dA, dB, ks_res) in enumerate(panels):
        ax.hist(dA[key], bins=50, alpha=0.5, color='tab:blue', label='Saueressig', density=True)
        ax.hist(dB[key], bins=50, alpha=0.5, color='tab:red', label='Cantrell', density=True)
        ax.axvline(0, color='black', ls='--', lw=0.8)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel('Δ Emissions (Tg/yr)', fontsize=9)
        ax.set_ylabel('Density', fontsize=9)
        ax.legend(fontsize=7, loc='upper right')
        ax.grid(alpha=0.2)

        # Annotate
        ksr_label = ksr_vals[idx] if idx % 2 == 1 else None
        p_str = f"p={ks_res.pvalue:.2e}"
        annot = p_str
        if idx % 2 == 1:
            ksr_name = 'FF' if idx < 2 else 'Mic'
            annot += f"\nKSR={ksr_vals[idx]:.2f}"
        ax.text(0.03, 0.95, f"{panel_labels[idx]}\n{annot}",
                transform=ax.transAxes, fontsize=8, va='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.suptitle('KIE Sensitivity: δ¹³C-only vs Dual-Isotope (δ¹³C + δD)', fontsize=12, y=1.02)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig1_KSR_summary.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"\n  Figure 1 saved: {FIG_DIR / 'fig1_KSR_summary.png'}")

    # --- Figure 2: Uncertainty time series ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=300)

    for ax, (key, title, unc_arr) in zip(axes, [
        ('FF', 'Fossil Fuel', unc_red_ff), ('Mic', 'Microbial', unc_red_mic)
    ]):
        sig_single = 2 * np.nanstd(p1_C[key], axis=1)
        sig_dual = 2 * np.nanstd(p2_C[key], axis=1)
        ax.plot(years, sig_single, 'r-', lw=2, label='δ¹³C-only (2σ)')
        ax.plot(years, sig_dual, 'b-', lw=2, label='Dual isotope (2σ)')
        ax.fill_between(years, sig_dual, sig_single, alpha=0.2, color='green',
                         label=f'Reduction (mean {np.mean(unc_arr):.1f}%)')
        ax.set_xlabel('Year')
        ax.set_ylabel('2σ Uncertainty (Tg/yr)')
        ax.set_title(f'{title} — Uncertainty Band Width')
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig2_uncertainty_timeseries.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Figure 2 saved: {FIG_DIR / 'fig2_uncertainty_timeseries.png'}")

    # --- Figure 3: Emission time series ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=300)

    for ax, (key, title) in zip(axes, [('FF', 'Fossil Fuel'), ('Mic', 'Microbial')]):
        med_s = np.nanmedian(p1_C[key], axis=1)
        lo_s = np.nanpercentile(p1_C[key], 2.5, axis=1)
        hi_s = np.nanpercentile(p1_C[key], 97.5, axis=1)

        med_d = np.nanmedian(p2_C[key], axis=1)
        lo_d = np.nanpercentile(p2_C[key], 2.5, axis=1)
        hi_d = np.nanpercentile(p2_C[key], 97.5, axis=1)

        ax.plot(years, med_s, 'r-', lw=2, label='δ¹³C-only median')
        ax.fill_between(years, lo_s, hi_s, alpha=0.15, color='red')
        ax.plot(years, med_d, 'b-', lw=2, label='Dual isotope median')
        ax.fill_between(years, lo_d, hi_d, alpha=0.15, color='blue')

        ax.set_xlabel('Year')
        ax.set_ylabel('Emissions (Tg/yr)')
        ax.set_title(f'{title} — Sampled KIE (Run C)')
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig3_emission_timeseries.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Figure 3 saved: {FIG_DIR / 'fig3_emission_timeseries.png'}")

    print(f"\n✓ Phase 3 complete. Summary at {OUT_DIR / 'summary.json'}")


if __name__ == "__main__":
    main()
