#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 8: Fine-Resolution Threshold Sweep + Temporal Stability
==============================================================
Phase 6b sampled 7 thresholds (25, 50, 75, 100, 150, 200, 300 Tg/yr) with
a coarse step. Phase 8 does:

(1) Fine sweep: thresholds 30..220 in 10-Tg steps
    - Pinpoint the threshold that maximises (i) KSR and (ii) discriminant power
    - Compute bootstrap-CI on the discriminant difference (Cantrell − Saueressig)
      at each threshold to find where it first becomes statistically significant.

(2) Temporal stability: split the record into 3 sub-periods (1999–2006, 2007–2014,
    2015–2022). For each sub-period we compute:
      - agreement rate (Cantrell, Saueressig)
      - discriminant Δ
    This tests whether the agreement-rate signal is stable across the renewed
    growth (2007), the plateau, and the post-2014 acceleration — i.e. is the
    discriminant a robust property of the atmospheric record or is it driven
    by one specific epoch?

Output:
  results/phase8_fine_thresholds/summary_N5000.json
  figures/fig13_fine_threshold_N5000.png
  figures/fig14_temporal_stability_N5000.png
"""

import json
import sys
from pathlib import Path

import numpy as np
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
    smooth_5yr, trend_change,
    SINK_FRACTIONS_GLOBAL, PT,
)

BASE = Path(__file__).resolve().parent
OUT_DIR = BASE / "results" / "phase8_fine_thresholds"
FIG_DIR = BASE / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

N_ITER = 5000
SEED = 42
FINE_THRESHOLDS = list(range(30, 221, 10))  # 30..220 step 10  (20 values)


def run_full(oh13c_mode: str):
    """Same as Phase 6b run_full_inversions but trimmed."""
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

    FF_c13 = np.full((n, N_ITER), np.nan)
    FF_dD = np.full((n, N_ITER), np.nan)
    Mic_c13 = np.full((n, N_ITER), np.nan)

    for k in range(N_ITER):
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
        n13 = f13 * CH4[:n + 1] * PT
        nD = fD * CH4[:n + 1] * PT
        sigs = sample_source_signatures(rng, data, k, n)

        for j in range(n):
            S = SumSource[j]
            BB_j = BB[j] if j < len(BB) else data.BB_global_mean
            d13C_src_f = (n13[j+1] - n13[j] + n13[j] * alpha_13C / tau[j]) / S
            d13C_src = fraction_to_delta_d13C(d13C_src_f)
            denom_c = sigs['ff_d13C'][j] - sigs['mic_d13C'][j]
            if abs(denom_c) > 0.1:
                ff_c = (S * d13C_src - sigs['mic_d13C'][j] * (S - BB_j) -
                        sigs['bb_d13C'][j] * BB_j) / denom_c
                FF_c13[j, k] = ff_c
                Mic_c13[j, k] = S - BB_j - ff_c

            dD_src_f = (nD[j+1] - nD[j] + nD[j] * alpha_D / tau[j]) / S
            dD_src = fraction_to_delta_dD(dD_src_f)
            denom_d = sigs['ff_dD'][j] - sigs['mic_dD'][j]
            if abs(denom_d) > 0.5:
                FF_dD[j, k] = (S * dD_src - sigs['mic_dD'][j] * (S - BB_j) -
                               sigs['bb_dD'][j] * BB_j) / denom_d

    return FF_c13, FF_dD, Mic_c13, years


def trend_for_filtered(FF_c13, Mic_c13, agree_mask, years):
    n, N = FF_c13.shape
    good = (agree_mask.sum(axis=0) >= n * 0.8)
    n_good = int(good.sum())
    if n_good < 30:
        return float('nan'), float('nan'), n_good
    FF_s = smooth_5yr(np.clip(FF_c13[:, good], 0, None))
    Mic_s = smooth_5yr(np.clip(Mic_c13[:, good], 0, None))
    delta_ff = trend_change(FF_s, years)[0]
    delta_mic = trend_change(Mic_s, years)[0]
    return float(delta_ff.mean()), float(delta_mic.mean()), n_good


def bootstrap_rate(diffs, valid, threshold, n_boot=2000):
    """Bootstrap a single rate from a 2D (years × iters) flattened diff matrix."""
    rng = np.random.default_rng(101)
    flat_valid = valid.flatten()
    a = ((diffs < threshold) & valid).flatten()[flat_valid]
    nv = len(a)
    out = np.empty(n_boot)
    for b in range(n_boot):
        out[b] = a[rng.choice(nv, nv, replace=True)].mean()
    return out


def main():
    print("=" * 60)
    print("  Phase 8 — Fine Threshold Sweep + Temporal Stability")
    print("=" * 60)

    print("\n  Running full inversions (Saueressig, Cantrell, sampled)...")
    runs = {}
    for oh in ['saueressig', 'cantrell', 'sampled']:
        FF_c, FF_d, Mic_c, yrs = run_full(oh)
        runs[oh] = {'FF_c13': FF_c, 'FF_dD': FF_d, 'Mic_c13': Mic_c, 'years': yrs}
        print(f"    {oh}: done")

    years = runs['saueressig']['years']
    n = len(years)

    # === Part 1: fine threshold sweep ===
    print(f"\n  Fine threshold sweep ({len(FINE_THRESHOLDS)} thresholds)...")
    sweep = {}
    for t in FINE_THRESHOLDS:
        # Compute rate, discriminant, KSR
        s_diffs = np.abs(runs['saueressig']['FF_c13'] - runs['saueressig']['FF_dD'])
        c_diffs = np.abs(runs['cantrell']['FF_c13'] - runs['cantrell']['FF_dD'])
        s_valid = ~(np.isnan(runs['saueressig']['FF_c13']) | np.isnan(runs['saueressig']['FF_dD']))
        c_valid = ~(np.isnan(runs['cantrell']['FF_c13']) | np.isnan(runs['cantrell']['FF_dD']))

        s_agree = (s_diffs < t) & s_valid
        c_agree = (c_diffs < t) & c_valid
        s_rate = s_agree.sum() / max(s_valid.sum(), 1)
        c_rate = c_agree.sum() / max(c_valid.sum(), 1)

        # Bootstrap CI on the difference
        s_boot = bootstrap_rate(s_diffs, s_valid, t)
        c_boot = bootstrap_rate(c_diffs, c_valid, t)
        diff_boot = c_boot - s_boot
        diff_ci = (float(np.percentile(diff_boot, 2.5)),
                   float(np.percentile(diff_boot, 97.5)))
        diff_significant = diff_ci[0] > 0  # CI fully above zero

        # KSR
        s_ff_mean, _, s_ng = trend_for_filtered(
            runs['saueressig']['FF_c13'], runs['saueressig']['Mic_c13'], s_agree, years)
        c_ff_mean, _, c_ng = trend_for_filtered(
            runs['cantrell']['FF_c13'], runs['cantrell']['Mic_c13'], c_agree, years)
        spread_c13 = 1.9823  # from Phase 1 baseline (constant)
        if not (np.isnan(s_ff_mean) or np.isnan(c_ff_mean)):
            spread_agree = abs(c_ff_mean - s_ff_mean)
            ksr = spread_c13 / max(spread_agree, 1e-6)
        else:
            ksr = float('nan')

        sweep[t] = {
            'rate_S': float(s_rate),
            'rate_C': float(c_rate),
            'diff_pp': float((c_rate - s_rate) * 100),
            'diff_CI_pp': [diff_ci[0]*100, diff_ci[1]*100],
            'significant': bool(diff_significant),
            'KSR_FF': float(ksr) if not np.isnan(ksr) else None,
            'n_good_S': s_ng,
            'n_good_C': c_ng,
        }

    # Find optima
    valid_ksrs = [(t, sweep[t]['KSR_FF']) for t in FINE_THRESHOLDS if sweep[t]['KSR_FF'] is not None]
    best_ksr_t, best_ksr = max(valid_ksrs, key=lambda x: x[1])
    best_disc_t, best_disc = max(((t, sweep[t]['diff_pp']) for t in FINE_THRESHOLDS),
                                  key=lambda x: x[1])
    sig_thresholds = [t for t in FINE_THRESHOLDS if sweep[t]['significant']]
    print(f"    Best KSR: t={best_ksr_t} (KSR={best_ksr:.2f})")
    print(f"    Max discriminant: t={best_disc_t} ({best_disc:.1f} pp)")
    print(f"    Significant range: t ∈ [{min(sig_thresholds)}, {max(sig_thresholds)}]")

    # === Part 2: temporal stability ===
    print("\n  Temporal stability across 3 sub-periods...")
    epochs = {
        'epoch1_1999_2006': (1999, 2006),
        'epoch2_2007_2014': (2007, 2014),
        'epoch3_2015_2022': (2015, 2022),
    }
    epoch_results = {}
    for ename, (lo, hi) in epochs.items():
        mask = (years >= lo) & (years <= hi)
        ep = {}
        for oh in ['saueressig', 'cantrell']:
            d = np.abs(runs[oh]['FF_c13'][mask] - runs[oh]['FF_dD'][mask])
            v = ~(np.isnan(runs[oh]['FF_c13'][mask]) | np.isnan(runs[oh]['FF_dD'][mask]))
            agree = (d < 100) & v
            rate = agree.sum() / max(v.sum(), 1)
            # Bootstrap
            flat_v = v.flatten()
            a = ((d < 100) & v).flatten()[flat_v]
            rng = np.random.default_rng(202)
            nv = len(a)
            boot = np.array([a[rng.choice(nv, nv, replace=True)].mean() for _ in range(2000)])
            ci = (float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5)))
            ep[oh] = {'rate': float(rate), 'CI': list(ci)}
        diff = ep['cantrell']['rate'] - ep['saueressig']['rate']
        sig = ep['saueressig']['CI'][1] < ep['cantrell']['CI'][0]
        ep['discriminant_pp'] = float(diff * 100)
        ep['significant'] = bool(sig)
        epoch_results[ename] = ep
        print(f"    {ename}: S={ep['saueressig']['rate']:.1%}, "
              f"C={ep['cantrell']['rate']:.1%}, Δ={diff*100:+.1f}pp, sig={sig}")

    # === Save ===
    summary = {
        'fine_threshold_sweep': {str(t): sweep[t] for t in FINE_THRESHOLDS},
        'optimal': {
            'best_KSR_threshold': best_ksr_t,
            'best_KSR_value': float(best_ksr),
            'max_discriminant_threshold': best_disc_t,
            'max_discriminant_pp': float(best_disc),
            'significant_threshold_min': int(min(sig_thresholds)) if sig_thresholds else None,
            'significant_threshold_max': int(max(sig_thresholds)) if sig_thresholds else None,
        },
        'temporal_stability': epoch_results,
    }
    with open(OUT_DIR / "summary_N5000.json", 'w') as f:
        json.dump(summary, f, indent=2)

    # === Figure 13: fine threshold ===
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), dpi=300)

    ts = np.array(FINE_THRESHOLDS)
    rates_S = np.array([sweep[t]['rate_S'] for t in ts]) * 100
    rates_C = np.array([sweep[t]['rate_C'] for t in ts]) * 100
    diffs = np.array([sweep[t]['diff_pp'] for t in ts])
    ci_lo = np.array([sweep[t]['diff_CI_pp'][0] for t in ts])
    ci_hi = np.array([sweep[t]['diff_CI_pp'][1] for t in ts])
    ksrs = np.array([sweep[t]['KSR_FF'] if sweep[t]['KSR_FF'] is not None else np.nan for t in ts])

    ax = axes[0]
    ax.plot(ts, rates_S, 'b-o', lw=2, markersize=4, label='Saueressig')
    ax.plot(ts, rates_C, 'r-o', lw=2, markersize=4, label='Cantrell')
    ax.set_xlabel('Threshold (Tg/yr)')
    ax.set_ylabel('Agreement Rate (%)')
    ax.set_title('(a) Fine-Resolution Agreement Rate')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    ax = axes[1]
    ax.plot(ts, diffs, 'g-s', lw=2, markersize=5)
    ax.fill_between(ts, ci_lo, ci_hi, color='g', alpha=0.2, label='95% CI')
    ax.axhline(0, color='black', lw=0.5)
    ax.axvline(best_disc_t, color='purple', ls='--', lw=1,
               label=f'max Δ at {best_disc_t}')
    ax.set_xlabel('Threshold (Tg/yr)')
    ax.set_ylabel('Δ Agreement Rate (pp)')
    ax.set_title('(b) Discriminant Power (Cantrell − Saueressig)')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    ax = axes[2]
    ax.plot(ts, ksrs, 'm-^', lw=2, markersize=5)
    ax.axhline(1.0, color='gray', ls='--', lw=1, label='KSR=1')
    ax.axvline(best_ksr_t, color='orange', ls='--', lw=1,
               label=f'max KSR at {best_ksr_t}')
    ax.set_xlabel('Threshold (Tg/yr)')
    ax.set_ylabel('KSR (FF)')
    ax.set_title('(c) KSR(FF) — Higher = more KIE damping')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    plt.suptitle('Phase 8a: Fine-Resolution Threshold Sweep', fontsize=12)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig13_fine_threshold_N5000.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"\n  Figure 13 saved: {FIG_DIR / 'fig13_fine_threshold_N5000.png'}")

    # === Figure 14: temporal stability ===
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=300)

    enames = list(epochs.keys())
    labels = ['1999–2006\n(plateau)', '2007–2014\n(renewed growth)', '2015–2022\n(acceleration)']
    rates_S_ep = [epoch_results[e]['saueressig']['rate'] * 100 for e in enames]
    rates_C_ep = [epoch_results[e]['cantrell']['rate'] * 100 for e in enames]
    ci_S = [epoch_results[e]['saueressig']['CI'] for e in enames]
    ci_C = [epoch_results[e]['cantrell']['CI'] for e in enames]
    err_S = [[r - c[0]*100 for r, c in zip(rates_S_ep, ci_S)],
             [c[1]*100 - r for r, c in zip(rates_S_ep, ci_S)]]
    err_C = [[r - c[0]*100 for r, c in zip(rates_C_ep, ci_C)],
             [c[1]*100 - r for r, c in zip(rates_C_ep, ci_C)]]

    ax = axes[0]
    x = np.arange(len(enames))
    w = 0.35
    ax.bar(x - w/2, rates_S_ep, w, yerr=err_S, color='tab:blue', alpha=0.7,
           label='Saueressig', capsize=4, edgecolor='black')
    ax.bar(x + w/2, rates_C_ep, w, yerr=err_C, color='tab:red', alpha=0.7,
           label='Cantrell', capsize=4, edgecolor='black')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel('Agreement Rate (%)')
    ax.set_title('(a) Agreement Rate by Epoch')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, axis='y')

    ax = axes[1]
    discs = [epoch_results[e]['discriminant_pp'] for e in enames]
    sigs = [epoch_results[e]['significant'] for e in enames]
    bar_colors = ['green' if s else 'gray' for s in sigs]
    ax.bar(x, discs, color=bar_colors, alpha=0.75, edgecolor='black')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel('Discriminant Δ (pp)')
    ax.set_title('(b) Cantrell − Saueressig by Epoch\n(green = significant)')
    ax.grid(alpha=0.3, axis='y')
    ax.axhline(0, color='black', lw=0.5)
    for xi, d in zip(x, discs):
        ax.text(xi, d + 0.5, f'{d:+.1f}', ha='center', fontsize=9, fontweight='bold')

    plt.suptitle('Phase 8b: Temporal Stability of the KIE Discriminant', fontsize=12)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig14_temporal_stability_N5000.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Figure 14 saved: {FIG_DIR / 'fig14_temporal_stability_N5000.png'}")

    print("\n  === SUMMARY ===")
    print(f"  Best KSR threshold: {best_ksr_t} Tg/yr (KSR={best_ksr:.2f})")
    print(f"  Max discriminant threshold: {best_disc_t} Tg/yr ({best_disc:.1f} pp)")
    print(f"  Statistically significant range: {min(sig_thresholds)}–{max(sig_thresholds)} Tg/yr")
    print(f"  Discriminant stable across all 3 epochs: "
          f"{all(epoch_results[e]['significant'] for e in enames)}")
    print("\n✓ Phase 8 complete.")


if __name__ == "__main__":
    main()
