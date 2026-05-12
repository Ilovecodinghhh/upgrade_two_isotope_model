#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 6b: Agreement Threshold Sensitivity + KIE Discriminant
==============================================================
Extends Phase 6 with:
1. Sweep agreement threshold (25, 50, 75, 100, 150, 200 Tg/yr)
2. Test agreement rate as OH-¹³C KIE discriminant
3. Bootstrap CI on agreement rates to test if Cantrell vs Saueressig
   difference is statistically significant
4. Time-varying vs fixed lifetime effect on agreement
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import bootstrap
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

BASE = Path(__file__).resolve().parent
OUT_DIR = BASE / "results" / "phase6b_threshold_sweep"
FIG_DIR = BASE / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

N_ITER = 1000
SEED = 42

# Threshold sweep
THRESHOLDS = [25, 50, 75, 100, 150, 200, 300]


def run_full_inversions(oh13c_mode: str, lifetime_mode: str = 'varying'):
    """Run both δ¹³C and δD inversions, return per-year per-iteration FF values."""
    data = load_data(REPO_ROOT, two_box=False)
    n = data.n_years
    years = data.model_years
    CH4 = data.CH4_global
    rng = np.random.default_rng(SEED)

    tau = compute_lifetime(years, lifetime_mode)

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
        n13 = f13 * CH4[:n+1] * PT
        nD = fD * CH4[:n+1] * PT

        sigs = sample_source_signatures(rng, data, k, n)

        for j in range(n):
            S = SumSource[j]
            BB_j = BB[j] if j < len(BB) else data.BB_global_mean

            # δ¹³C
            d13C_src_f = (n13[j+1] - n13[j] + n13[j] * alpha_13C / tau[j]) / S
            d13C_src = fraction_to_delta_d13C(d13C_src_f)
            denom_c = sigs['ff_d13C'][j] - sigs['mic_d13C'][j]
            if abs(denom_c) > 0.1:
                ff_c = (S * d13C_src - sigs['mic_d13C'][j] * (S - BB_j) -
                        sigs['bb_d13C'][j] * BB_j) / denom_c
                FF_c13[j, k] = ff_c
                Mic_c13[j, k] = S - BB_j - ff_c

            # δD
            dD_src_f = (nD[j+1] - nD[j] + nD[j] * alpha_D / tau[j]) / S
            dD_src = fraction_to_delta_dD(dD_src_f)
            denom_d = sigs['ff_dD'][j] - sigs['mic_dD'][j]
            if abs(denom_d) > 0.5:
                ff_d = (S * dD_src - sigs['mic_dD'][j] * (S - BB_j) -
                        sigs['bb_dD'][j] * BB_j) / denom_d
                FF_dD[j, k] = ff_d

    return FF_c13, FF_dD, Mic_c13, years


def compute_agreement_metrics(FF_c13, FF_dD, Mic_c13, years, threshold):
    """Compute agreement rate and filtered trend for a given threshold."""
    n, N = FF_c13.shape
    agreement = np.abs(FF_c13 - FF_dD) < threshold
    # Ignore NaNs
    valid = ~(np.isnan(FF_c13) | np.isnan(FF_dD))
    agreement = agreement & valid

    # Per-year agreement rate
    agree_rate = np.zeros(n)
    for j in range(n):
        v = valid[j].sum()
        agree_rate[j] = agreement[j].sum() / max(v, 1)

    # Overall agreement
    overall_rate = agreement.sum() / max(valid.sum(), 1)

    # Filter: keep iterations where ≥80% of years agree
    good_iters = (agreement.sum(axis=0) >= n * 0.8)
    n_good = good_iters.sum()

    if n_good > 30:
        FF_filt = np.clip(FF_c13[:, good_iters], 0, None)
        Mic_filt = np.clip(Mic_c13[:, good_iters], 0, None)
        FF_s = smooth_5yr(FF_filt)
        Mic_s = smooth_5yr(Mic_filt)
        delta_ff, _ = trend_change(FF_s, years)
        delta_mic, _ = trend_change(Mic_s, years)
    else:
        delta_ff = np.array([np.nan])
        delta_mic = np.array([np.nan])

    return {
        'overall_rate': float(overall_rate),
        'agree_rate_by_year': agree_rate,
        'n_good_iters': int(n_good),
        'delta_ff_mean': float(np.nanmean(delta_ff)),
        'delta_ff_std': float(np.nanstd(delta_ff)),
        'delta_mic_mean': float(np.nanmean(delta_mic)),
        'delta_mic_std': float(np.nanstd(delta_mic)),
    }


def bootstrap_agreement_rate(FF_c13, FF_dD, threshold, n_boot=2000):
    """Bootstrap CI for overall agreement rate."""
    valid = ~(np.isnan(FF_c13) | np.isnan(FF_dD))
    diffs = np.abs(FF_c13 - FF_dD)
    agree_flat = ((diffs < threshold) & valid).flatten()
    valid_flat = valid.flatten()

    # Only use valid entries
    agree_valid = agree_flat[valid_flat]
    rng = np.random.default_rng(99)

    rates = np.zeros(n_boot)
    n_valid = len(agree_valid)
    for i in range(n_boot):
        idx = rng.choice(n_valid, size=n_valid, replace=True)
        rates[i] = agree_valid[idx].mean()

    return float(np.percentile(rates, 2.5)), float(np.percentile(rates, 97.5))


def main():
    print("=" * 60)
    print("  Phase 6b — Threshold Sweep + KIE Discriminant")
    print("=" * 60)

    # === Part 1: Run inversions for both KIE values ===
    print("\n  Running inversions...")
    results = {}
    for oh13c in ['saueressig', 'cantrell', 'sampled']:
        print(f"    OH_13C = {oh13c}")
        FF_c, FF_d, Mic_c, years = run_full_inversions(oh13c, 'varying')
        results[oh13c] = {'FF_c13': FF_c, 'FF_dD': FF_d, 'Mic_c13': Mic_c, 'years': years}

    # === Part 2: Threshold sweep ===
    print("\n  Threshold sweep...")
    sweep_results = {}
    for thresh in THRESHOLDS:
        sweep_results[thresh] = {}
        for oh13c in ['saueressig', 'cantrell', 'sampled']:
            r = results[oh13c]
            m = compute_agreement_metrics(r['FF_c13'], r['FF_dD'], r['Mic_c13'], r['years'], thresh)
            sweep_results[thresh][oh13c] = m

        # KSR at this threshold
        spread_c13 = abs(
            compute_agreement_metrics(results['cantrell']['FF_c13'], results['cantrell']['FF_dD'],
                                      results['cantrell']['Mic_c13'], years, 9999)['delta_ff_mean'] -
            compute_agreement_metrics(results['saueressig']['FF_c13'], results['saueressig']['FF_dD'],
                                      results['saueressig']['Mic_c13'], years, 9999)['delta_ff_mean']
        )
        if sweep_results[thresh]['cantrell']['n_good_iters'] > 30 and \
           sweep_results[thresh]['saueressig']['n_good_iters'] > 30:
            spread_agree = abs(sweep_results[thresh]['cantrell']['delta_ff_mean'] -
                               sweep_results[thresh]['saueressig']['delta_ff_mean'])
            ksr = spread_c13 / max(spread_agree, 1e-6)
        else:
            ksr = np.nan
            spread_agree = np.nan

        sweep_results[thresh]['KSR_FF'] = float(ksr)
        sweep_results[thresh]['spread_c13'] = float(spread_c13)
        sweep_results[thresh]['spread_agree'] = float(spread_agree)

        rate_s = sweep_results[thresh]['saueressig']['overall_rate']
        rate_c = sweep_results[thresh]['cantrell']['overall_rate']
        print(f"    threshold={thresh:3d}: agree_S={rate_s:.1%}, agree_C={rate_c:.1%}, "
              f"KSR={ksr:.2f}, n_good(S)={sweep_results[thresh]['saueressig']['n_good_iters']}, "
              f"n_good(C)={sweep_results[thresh]['cantrell']['n_good_iters']}")

    # === Part 3: Bootstrap CIs on agreement rate difference ===
    print("\n  Bootstrap CIs on agreement rate difference (threshold=100)...")
    ci_s = bootstrap_agreement_rate(results['saueressig']['FF_c13'],
                                    results['saueressig']['FF_dD'], 100)
    ci_c = bootstrap_agreement_rate(results['cantrell']['FF_c13'],
                                    results['cantrell']['FF_dD'], 100)

    rate_s_100 = sweep_results[100]['saueressig']['overall_rate']
    rate_c_100 = sweep_results[100]['cantrell']['overall_rate']
    diff = rate_c_100 - rate_s_100

    print(f"    Saueressig: {rate_s_100:.1%} [{ci_s[0]:.1%}, {ci_s[1]:.1%}]")
    print(f"    Cantrell:   {rate_c_100:.1%} [{ci_c[0]:.1%}, {ci_c[1]:.1%}]")
    print(f"    Difference: {diff:.1%}")
    # If CIs don't overlap, difference is significant
    significant = ci_s[1] < ci_c[0]
    print(f"    Significant (non-overlapping CIs): {'YES' if significant else 'NO'}")

    # === Part 4: Fixed vs varying lifetime ===
    print("\n  Lifetime comparison (threshold=100)...")
    lifetime_results = {}
    for lt_mode in ['varying', 'fixed']:
        lifetime_results[lt_mode] = {}
        for oh13c in ['saueressig', 'cantrell']:
            FF_c, FF_d, Mic_c, yrs = run_full_inversions(oh13c, lt_mode)
            m = compute_agreement_metrics(FF_c, FF_d, Mic_c, yrs, 100)
            lifetime_results[lt_mode][oh13c] = m
            print(f"    τ={lt_mode}, OH_13C={oh13c}: agree={m['overall_rate']:.1%}")

    # === Save ===
    summary = {
        'threshold_sweep': {str(k): {
            'KSR_FF': v['KSR_FF'],
            'spread_c13': v['spread_c13'],
            'spread_agree': v['spread_agree'],
            'saueressig_rate': v['saueressig']['overall_rate'],
            'cantrell_rate': v['cantrell']['overall_rate'],
            'saueressig_n_good': v['saueressig']['n_good_iters'],
            'cantrell_n_good': v['cantrell']['n_good_iters'],
        } for k, v in sweep_results.items()},
        'bootstrap_CI': {
            'saueressig': {'rate': rate_s_100, 'CI_95': list(ci_s)},
            'cantrell': {'rate': rate_c_100, 'CI_95': list(ci_c)},
            'difference': diff,
            'significant': bool(significant),
        },
        'lifetime_effect': {
            lt: {oh: {'rate': m['overall_rate'], 'n_good': m['n_good_iters']}
                 for oh, m in v.items()}
            for lt, v in lifetime_results.items()
        },
    }
    with open(OUT_DIR / "summary.json", 'w') as f:
        json.dump(summary, f, indent=2)

    # === Figures ===

    # Figure 9: Threshold sweep — KSR and agreement rates
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), dpi=300)

    # (a) Agreement rate vs threshold
    ax = axes[0]
    rates_s = [sweep_results[t]['saueressig']['overall_rate'] for t in THRESHOLDS]
    rates_c = [sweep_results[t]['cantrell']['overall_rate'] for t in THRESHOLDS]
    rates_m = [sweep_results[t]['sampled']['overall_rate'] for t in THRESHOLDS]
    ax.plot(THRESHOLDS, [r*100 for r in rates_s], 'b-o', lw=2, label='Saueressig (1.0039)')
    ax.plot(THRESHOLDS, [r*100 for r in rates_c], 'r-o', lw=2, label='Cantrell (1.0054)')
    ax.plot(THRESHOLDS, [r*100 for r in rates_m], 'k--o', lw=1, label='Sampled')
    ax.set_xlabel('Agreement Threshold (Tg/yr)')
    ax.set_ylabel('Overall Agreement Rate (%)')
    ax.set_title('(a) Agreement Rate vs Threshold')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # (b) KSR vs threshold
    ax = axes[1]
    ksrs = [sweep_results[t]['KSR_FF'] for t in THRESHOLDS]
    ax.plot(THRESHOLDS, ksrs, 'g-s', lw=2, markersize=8)
    ax.axhline(1.0, color='gray', ls='--', lw=1, label='KSR=1 (no improvement)')
    ax.set_xlabel('Agreement Threshold (Tg/yr)')
    ax.set_ylabel('KSR (FF)')
    ax.set_title('(b) KIE Sensitivity Ratio vs Threshold')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    valid_ksrs_plot = [k for k in ksrs if not np.isnan(k)]
    ax.set_ylim(0, max(valid_ksrs_plot) * 1.2 if valid_ksrs_plot else 5)

    # (c) Agreement rate difference (KIE discriminant)
    ax = axes[2]
    diffs = [sweep_results[t]['cantrell']['overall_rate'] -
             sweep_results[t]['saueressig']['overall_rate'] for t in THRESHOLDS]
    ax.bar(range(len(THRESHOLDS)), [d*100 for d in diffs], color='purple', alpha=0.7)
    ax.set_xticks(range(len(THRESHOLDS)))
    ax.set_xticklabels([str(t) for t in THRESHOLDS])
    ax.set_xlabel('Agreement Threshold (Tg/yr)')
    ax.set_ylabel('Rate(Cantrell) − Rate(Saueressig) (pp)')
    ax.set_title('(c) KIE Discriminant Power\n(Cantrell gives higher agreement)')
    ax.grid(alpha=0.3)
    ax.axhline(0, color='black', lw=0.5)

    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig9_threshold_sweep.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"\n  Figure 9 saved: {FIG_DIR / 'fig9_threshold_sweep.png'}")

    # Figure 10: Agreement rate time series by KIE (at optimal threshold)
    # Find threshold with max KSR
    valid_ksrs = [(t, k) for t, k in zip(THRESHOLDS, ksrs) if not np.isnan(k)]
    if valid_ksrs:
        best_thresh = max(valid_ksrs, key=lambda x: x[1])[0]
    else:
        best_thresh = 100

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=300)

    ax = axes[0]
    for oh13c, color, name in [('saueressig', 'tab:blue', 'Saueressig'),
                                ('cantrell', 'tab:red', 'Cantrell'),
                                ('sampled', 'gray', 'Sampled')]:
        rate_by_year = sweep_results[best_thresh][oh13c]['agree_rate_by_year']
        ax.plot(years, rate_by_year * 100, '-', color=color, lw=2, label=name)
    ax.set_xlabel('Year')
    ax.set_ylabel('Agreement Rate (%)')
    ax.set_title(f'(a) Per-Year Agreement Rate (threshold={best_thresh} Tg/yr)')
    ax.legend()
    ax.grid(alpha=0.3)

    # (b) Lifetime effect
    ax = axes[1]
    lt_data = [
        ('varying', 'saueressig', 'Varying τ + Saueressig', 'tab:blue', '-'),
        ('varying', 'cantrell', 'Varying τ + Cantrell', 'tab:red', '-'),
        ('fixed', 'saueressig', 'Fixed τ + Saueressig', 'tab:blue', '--'),
        ('fixed', 'cantrell', 'Fixed τ + Cantrell', 'tab:red', '--'),
    ]
    bar_labels = []
    bar_vals = []
    bar_colors = []
    for lt, oh, name, color, ls in lt_data:
        bar_labels.append(name)
        bar_vals.append(lifetime_results[lt][oh]['overall_rate'] * 100)
        bar_colors.append(color if lt == 'varying' else 'lightblue' if 'Saueressig' in name else 'lightsalmon')

    bars = ax.bar(range(len(bar_labels)), bar_vals, color=bar_colors, alpha=0.7, edgecolor='black')
    ax.set_xticks(range(len(bar_labels)))
    ax.set_xticklabels(bar_labels, rotation=20, ha='right', fontsize=8)
    ax.set_ylabel('Agreement Rate (%)')
    ax.set_title('(b) Lifetime Mode × KIE Effect on Agreement')
    ax.grid(alpha=0.3, axis='y')
    for bar, val in zip(bars, bar_vals):
        ax.text(bar.get_x() + bar.get_width()/2, val + 0.5, f'{val:.1f}%',
                ha='center', fontsize=8)

    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig10_agreement_timeseries.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Figure 10 saved: {FIG_DIR / 'fig10_agreement_timeseries.png'}")

    print(f"\n  === SUMMARY ===")
    print(f"  Best threshold for KSR: {best_thresh} Tg/yr (KSR={max(valid_ksrs, key=lambda x: x[1])[1]:.2f})")
    print(f"  Agreement rate difference (Cantrell − Saueressig): {diff:.1%}")
    print(f"  Statistically significant: {'YES' if significant else 'NO'}")
    print(f"\n✓ Phase 6b complete.")


if __name__ == "__main__":
    main()
