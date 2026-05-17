#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 9: Editorial Assessment Fixes
====================================
Addresses the three quantitative gaps identified in EditorialAssessment_v1.md:

  [HIGH]   Item 3 — Increase MC iterations to N=5000 and add KSR confidence
                     intervals (bootstrap). Also report CI on the agreement-rate
                     discriminant at the headline threshold (90 Tg/yr).

  [MEDIUM] Item 5 — Test agreement-filter discriminant under Thanwerdas low-Cl
                     (0.6%) and high-Cl (6.5%) scenarios. Phase 5 only tested
                     Cl sensitivity for the WLS approach, not for the agreement
                     filter.

  [MEDIUM] Item 6 — Sweep the 80% year-agreement parameter (currently hardcoded).
                     Test 60%, 70%, 80%, 90%, 95% and show whether the
                     discriminant is insensitive.

Output:
  results/phase9_editorial_fixes/
    high_n_summary.json        — N=5000 headline numbers + KSR bootstrap CI
    cl_sensitivity.json        — discriminant at Cl = 0.6%, 3.5%, 6.5%
    year_agreement_sweep.json  — discriminant at year-agreement ∈ {60..95%}

  figures/
    fig15_high_n.png           — KSR + discriminant with N=5000 and CIs
    fig16_cl_sensitivity.png   — discriminant under 3 Cl scenarios
    fig17_year_agree_sweep.png — discriminant vs year-agreement threshold
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
OUT_DIR = BASE / "results" / "phase9_editorial_fixes"
FIG_DIR = BASE / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Core inversion — same as Phase 8 but parameterised for Cl fraction and N
# ---------------------------------------------------------------------------

def run_inversions(oh13c_mode: str, n_iter: int = 1000, seed: int = 42,
                   cl_fraction: float | None = None):
    """Run independent δ¹³C and δD inversions.

    Parameters
    ----------
    oh13c_mode : 'saueressig' | 'cantrell' | 'sampled'
    n_iter : number of MC iterations
    seed : RNG seed
    cl_fraction : if not None, override the Cl sink fraction (renormalise OH
                  to keep OH + Cl + Strat + Soil = 1).

    Returns (FF_c13, FF_dD, Mic_c13, years)  — each array is (n_years, n_iter)
    """
    data = load_data(REPO_ROOT, two_box=False)
    n = data.n_years
    years = data.model_years
    CH4 = data.CH4_global
    rng = np.random.default_rng(seed)
    tau = compute_lifetime(years, 'varying')

    # Build sink fractions (optionally override Cl)
    sf = dict(SINK_FRACTIONS_GLOBAL)  # copy
    if cl_fraction is not None:
        old_cl = sf['Cl']
        sf['Cl'] = cl_fraction
        # Redistribute the difference into OH to keep sum = 1
        sf['OH'] += (old_cl - cl_fraction)

    SumSource = np.zeros(n)
    for i in range(n):
        SumSource[i] = CH4[i + 1] * PT - CH4[i] * PT + CH4[i] * PT / tau[i]
    BB = data.BB_annual.copy()

    FF_c13 = np.full((n, n_iter), np.nan)
    FF_dD  = np.full((n, n_iter), np.nan)
    Mic_c13 = np.full((n, n_iter), np.nan)

    for k in range(n_iter):
        if k % 1000 == 0 and k > 0:
            print(f"        iter {k}/{n_iter}")

        kies = sample_KIE(rng, 'sampled')
        if oh13c_mode == 'saueressig':
            kies['OH_13C'] = 1.0039
        elif oh13c_mode == 'cantrell':
            kies['OH_13C'] = 1.0054

        KIE_13C, KIE_D = compute_bulk_KIE(kies, sf)
        alpha_13C = 1.0 / KIE_13C
        alpha_D   = 1.0 / KIE_D

        d13C_atm = sample_atm_d13C(data, k, n)
        dD_atm   = sample_atm_dD(data, k, n)
        f13 = delta_to_fraction_d13C(d13C_atm)
        fD  = delta_to_fraction_dD(dD_atm)
        n13 = f13 * CH4[:n + 1] * PT
        nD  = fD  * CH4[:n + 1] * PT
        sigs = sample_source_signatures(rng, data, k, n)

        for j in range(n):
            S = SumSource[j]
            BB_j = BB[j] if j < len(BB) else data.BB_global_mean

            # δ¹³C branch
            d13C_src_f = (n13[j+1] - n13[j] + n13[j] * alpha_13C / tau[j]) / S
            d13C_src = fraction_to_delta_d13C(d13C_src_f)
            denom_c = sigs['ff_d13C'][j] - sigs['mic_d13C'][j]
            if abs(denom_c) > 0.1:
                ff_c = (S * d13C_src - sigs['mic_d13C'][j] * (S - BB_j)
                        - sigs['bb_d13C'][j] * BB_j) / denom_c
                FF_c13[j, k] = ff_c
                Mic_c13[j, k] = S - BB_j - ff_c

            # δD branch
            dD_src_f = (nD[j+1] - nD[j] + nD[j] * alpha_D / tau[j]) / S
            dD_src = fraction_to_delta_dD(dD_src_f)
            denom_d = sigs['ff_dD'][j] - sigs['mic_dD'][j]
            if abs(denom_d) > 0.5:
                FF_dD[j, k] = (S * dD_src - sigs['mic_dD'][j] * (S - BB_j)
                               - sigs['bb_dD'][j] * BB_j) / denom_d

    return FF_c13, FF_dD, Mic_c13, years


# ---------------------------------------------------------------------------
# Agreement metrics — parameterised year-agreement fraction
# ---------------------------------------------------------------------------

def agreement_metrics(FF_c13, FF_dD, Mic_c13, years, threshold,
                      year_agree_frac: float = 0.80):
    """Compute agreement-rate and filtered trend.

    year_agree_frac: minimum fraction of years that must agree per iteration
    to keep that iteration (the '80% parameter').
    """
    n, N = FF_c13.shape
    valid = ~(np.isnan(FF_c13) | np.isnan(FF_dD))
    agree = (np.abs(FF_c13 - FF_dD) < threshold) & valid

    overall_rate = agree.sum() / max(valid.sum(), 1)

    # Keep iterations where ≥year_agree_frac of years agree
    good = (agree.sum(axis=0) >= n * year_agree_frac)
    n_good = int(good.sum())

    if n_good >= 20:
        FF_s = smooth_5yr(np.clip(FF_c13[:, good], 0, None))
        Mic_s = smooth_5yr(np.clip(Mic_c13[:, good], 0, None))
        delta_ff = trend_change(FF_s, years)[0]
        delta_mic = trend_change(Mic_s, years)[0]
        ff_trend_mean = float(delta_ff.mean())
        mic_trend_mean = float(delta_mic.mean())
    else:
        ff_trend_mean = float('nan')
        mic_trend_mean = float('nan')

    return {
        'rate': float(overall_rate),
        'n_good': n_good,
        'ff_trend': ff_trend_mean,
        'mic_trend': mic_trend_mean,
    }


def bootstrap_rate(FF_c13, FF_dD, threshold, n_boot=2000, seed=101):
    """Bootstrap CI for overall agreement rate."""
    valid = ~(np.isnan(FF_c13) | np.isnan(FF_dD))
    agree = ((np.abs(FF_c13 - FF_dD) < threshold) & valid).flatten()
    valid_flat = valid.flatten()
    a = agree[valid_flat]
    rng = np.random.default_rng(seed)
    nv = len(a)
    out = np.empty(n_boot)
    for b in range(n_boot):
        out[b] = a[rng.choice(nv, nv, replace=True)].mean()
    return out


def bootstrap_ksr(runs_s, runs_c, threshold, baseline_spread,
                  n_boot=2000, seed=303, year_agree_frac=0.80):
    """Bootstrap CI for KSR at a given threshold.

    Resamples MC iteration indices (with replacement) and recomputes
    the Cantrell–Saueressig FF trend spread, then KSR = baseline / spread.
    """
    n, N = runs_s['FF_c13'].shape
    years = runs_s['years']
    rng = np.random.default_rng(seed)

    ksr_samples = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.choice(N, N, replace=True)

        # Saueressig filtered
        s_fc = runs_s['FF_c13'][:, idx]
        s_fd = runs_s['FF_dD'][:, idx]
        s_mc = runs_s['Mic_c13'][:, idx]
        s_valid = ~(np.isnan(s_fc) | np.isnan(s_fd))
        s_agree = (np.abs(s_fc - s_fd) < threshold) & s_valid
        s_good = (s_agree.sum(axis=0) >= n * year_agree_frac)
        # Cantrell filtered
        c_fc = runs_c['FF_c13'][:, idx]
        c_fd = runs_c['FF_dD'][:, idx]
        c_mc = runs_c['Mic_c13'][:, idx]
        c_valid = ~(np.isnan(c_fc) | np.isnan(c_fd))
        c_agree = (np.abs(c_fc - c_fd) < threshold) & c_valid
        c_good = (c_agree.sum(axis=0) >= n * year_agree_frac)

        if s_good.sum() < 20 or c_good.sum() < 20:
            ksr_samples[b] = np.nan
            continue

        s_ff_s = smooth_5yr(np.clip(s_fc[:, s_good], 0, None))
        c_ff_s = smooth_5yr(np.clip(c_fc[:, c_good], 0, None))
        s_trend = trend_change(s_ff_s, years)[0].mean()
        c_trend = trend_change(c_ff_s, years)[0].mean()
        spread = abs(c_trend - s_trend)
        ksr_samples[b] = baseline_spread / max(spread, 1e-6)

    return ksr_samples


# ===================================================================
# TASK 1 — HIGH-N run (N=5000) + KSR bootstrap CIs
# ===================================================================

def task_high_n():
    N_ITER = 5000
    SEED = 42
    THRESHOLDS = [50, 90, 100, 150]  # key thresholds only

    print("=" * 60)
    print("  Task 1: High-N run (N=5000) + KSR bootstrap")
    print("=" * 60)

    runs = {}
    for oh in ['saueressig', 'cantrell']:
        print(f"    Running {oh} (N={N_ITER})...")
        FF_c, FF_d, Mic_c, yrs = run_inversions(oh, n_iter=N_ITER, seed=SEED)
        runs[oh] = {'FF_c13': FF_c, 'FF_dD': FF_d, 'Mic_c13': Mic_c, 'years': yrs}

    years = runs['saueressig']['years']
    n = len(years)

    # Baseline δ¹³C-only spread (unfiltered = threshold=infinity effectively)
    s_unf = agreement_metrics(runs['saueressig']['FF_c13'], runs['saueressig']['FF_dD'],
                              runs['saueressig']['Mic_c13'], years, 1e6)
    c_unf = agreement_metrics(runs['cantrell']['FF_c13'], runs['cantrell']['FF_dD'],
                              runs['cantrell']['Mic_c13'], years, 1e6)
    baseline_spread = abs(c_unf['ff_trend'] - s_unf['ff_trend'])
    print(f"    Baseline δ¹³C-only spread: {baseline_spread:.3f} Tg/yr")

    results = {}
    for t in THRESHOLDS:
        print(f"    Threshold {t}...")
        s_m = agreement_metrics(runs['saueressig']['FF_c13'], runs['saueressig']['FF_dD'],
                                runs['saueressig']['Mic_c13'], years, t)
        c_m = agreement_metrics(runs['cantrell']['FF_c13'], runs['cantrell']['FF_dD'],
                                runs['cantrell']['Mic_c13'], years, t)

        spread = abs(c_m['ff_trend'] - s_m['ff_trend'])
        ksr = baseline_spread / max(spread, 1e-6)

        # Bootstrap agreement rate CIs
        s_boot = bootstrap_rate(runs['saueressig']['FF_c13'], runs['saueressig']['FF_dD'], t)
        c_boot = bootstrap_rate(runs['cantrell']['FF_c13'], runs['cantrell']['FF_dD'], t)
        diff_boot = c_boot - s_boot
        disc_ci = (float(np.percentile(diff_boot, 2.5)) * 100,
                   float(np.percentile(diff_boot, 97.5)) * 100)

        # Bootstrap KSR CI
        ksr_boot = bootstrap_ksr(runs['saueressig'], runs['cantrell'], t, baseline_spread)
        ksr_valid = ksr_boot[~np.isnan(ksr_boot)]
        ksr_ci = (float(np.percentile(ksr_valid, 2.5)),
                  float(np.percentile(ksr_valid, 97.5))) if len(ksr_valid) > 100 else (np.nan, np.nan)

        results[t] = {
            'rate_S': s_m['rate'],
            'rate_C': c_m['rate'],
            'n_good_S': s_m['n_good'],
            'n_good_C': c_m['n_good'],
            'discriminant_pp': float((c_m['rate'] - s_m['rate']) * 100),
            'discriminant_CI_pp': list(disc_ci),
            'KSR': float(ksr),
            'KSR_CI_95': list(ksr_ci),
            'ff_trend_S': s_m['ff_trend'],
            'ff_trend_C': c_m['ff_trend'],
        }
        print(f"      S={s_m['rate']:.1%} (n={s_m['n_good']}), "
              f"C={c_m['rate']:.1%} (n={c_m['n_good']}), "
              f"Δ={results[t]['discriminant_pp']:.1f}pp, "
              f"KSR={ksr:.2f} [{ksr_ci[0]:.2f}, {ksr_ci[1]:.2f}]")

    summary = {
        'N_iterations': N_ITER,
        'baseline_spread_Tg': baseline_spread,
        'thresholds': {str(t): results[t] for t in THRESHOLDS},
    }

    with open(OUT_DIR / "high_n_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)

    # --- Figure 15: High-N results ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), dpi=300)

    ts = THRESHOLDS
    # (a) Agreement rates with CIs
    ax = axes[0]
    s_rates = [results[t]['rate_S'] * 100 for t in ts]
    c_rates = [results[t]['rate_C'] * 100 for t in ts]
    x = np.arange(len(ts))
    w = 0.35
    ax.bar(x - w/2, s_rates, w, color='tab:blue', alpha=0.7,
           label='Saueressig', edgecolor='black')
    ax.bar(x + w/2, c_rates, w, color='tab:red', alpha=0.7,
           label='Cantrell', edgecolor='black')
    ax.set_xticks(x)
    ax.set_xticklabels([str(t) for t in ts])
    ax.set_xlabel('Threshold (Tg/yr)')
    ax.set_ylabel('Agreement Rate (%)')
    ax.set_title(f'(a) Agreement Rate (N={N_ITER})')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, axis='y')

    # (b) Discriminant with CI
    ax = axes[1]
    discs = [results[t]['discriminant_pp'] for t in ts]
    ci_lo = [results[t]['discriminant_CI_pp'][0] for t in ts]
    ci_hi = [results[t]['discriminant_CI_pp'][1] for t in ts]
    err_lo = [d - lo for d, lo in zip(discs, ci_lo)]
    err_hi = [hi - d for d, hi in zip(discs, ci_hi)]
    ax.bar(x, discs, color='green', alpha=0.7, edgecolor='black',
           yerr=[err_lo, err_hi], capsize=5)
    ax.axhline(0, color='black', lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([str(t) for t in ts])
    ax.set_xlabel('Threshold (Tg/yr)')
    ax.set_ylabel('Discriminant Δ (pp)')
    ax.set_title('(b) Cantrell − Saueressig + 95% CI')
    ax.grid(alpha=0.3, axis='y')
    for xi, d in zip(x, discs):
        ax.text(xi, d + 1, f'{d:.1f}', ha='center', fontsize=9, fontweight='bold')

    # (c) KSR with CI
    ax = axes[2]
    ksrs = [results[t]['KSR'] for t in ts]
    ksr_lo = [results[t]['KSR_CI_95'][0] for t in ts]
    ksr_hi = [results[t]['KSR_CI_95'][1] for t in ts]
    err_lo = [k - lo for k, lo in zip(ksrs, ksr_lo)]
    err_hi = [hi - k for k, hi in zip(ksrs, ksr_hi)]
    ax.bar(x, ksrs, color='purple', alpha=0.7, edgecolor='black',
           yerr=[err_lo, err_hi], capsize=5)
    ax.axhline(1.0, color='gray', ls='--', lw=1, label='KSR=1')
    ax.set_xticks(x)
    ax.set_xticklabels([str(t) for t in ts])
    ax.set_xlabel('Threshold (Tg/yr)')
    ax.set_ylabel('KSR (FF)')
    ax.set_title('(c) KSR + 95% bootstrap CI')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, axis='y')
    for xi, k in zip(x, ksrs):
        ax.text(xi, k + 0.1, f'{k:.2f}', ha='center', fontsize=9, fontweight='bold')

    plt.suptitle(f'Phase 9a: High-N (N={N_ITER}) Results with Bootstrap CIs',
                 fontsize=12)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig15_high_n.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"\n  Figure 15 saved.")

    return results


# ===================================================================
# TASK 2 — Cl fraction sensitivity for the agreement filter
# ===================================================================

def task_cl_sensitivity():
    N_ITER = 5000
    SEED = 42
    THRESHOLD = 90  # headline threshold
    CL_SCENARIOS = {
        'Thanwerdas_low': 0.006,
        'Default':        0.035,
        'High_Cl':        0.065,
    }

    print("\n" + "=" * 60)
    print("  Task 2: Cl-fraction sensitivity for agreement filter")
    print("=" * 60)

    results = {}
    for cl_name, cl_frac in CL_SCENARIOS.items():
        print(f"\n  Cl = {cl_frac:.1%} ({cl_name})")
        runs = {}
        for oh in ['saueressig', 'cantrell']:
            print(f"    Running {oh}...")
            FF_c, FF_d, Mic_c, yrs = run_inversions(
                oh, n_iter=N_ITER, seed=SEED, cl_fraction=cl_frac)
            runs[oh] = {'FF_c13': FF_c, 'FF_dD': FF_d, 'Mic_c13': Mic_c, 'years': yrs}

        years = runs['saueressig']['years']
        s_m = agreement_metrics(runs['saueressig']['FF_c13'], runs['saueressig']['FF_dD'],
                                runs['saueressig']['Mic_c13'], years, THRESHOLD)
        c_m = agreement_metrics(runs['cantrell']['FF_c13'], runs['cantrell']['FF_dD'],
                                runs['cantrell']['Mic_c13'], years, THRESHOLD)

        # Bootstrap discriminant CI
        s_boot = bootstrap_rate(runs['saueressig']['FF_c13'], runs['saueressig']['FF_dD'], THRESHOLD)
        c_boot = bootstrap_rate(runs['cantrell']['FF_c13'], runs['cantrell']['FF_dD'], THRESHOLD)
        diff_boot = c_boot - s_boot
        disc_ci = (float(np.percentile(diff_boot, 2.5)) * 100,
                   float(np.percentile(diff_boot, 97.5)) * 100)
        significant = disc_ci[0] > 0

        results[cl_name] = {
            'cl_fraction': cl_frac,
            'rate_S': s_m['rate'],
            'rate_C': c_m['rate'],
            'n_good_S': s_m['n_good'],
            'n_good_C': c_m['n_good'],
            'discriminant_pp': float((c_m['rate'] - s_m['rate']) * 100),
            'discriminant_CI_pp': list(disc_ci),
            'significant': significant,
        }
        print(f"    S={s_m['rate']:.1%}, C={c_m['rate']:.1%}, "
              f"Δ={results[cl_name]['discriminant_pp']:.1f}pp, sig={significant}")

    with open(OUT_DIR / "cl_sensitivity.json", 'w') as f:
        json.dump(results, f, indent=2)

    # --- Figure 16: Cl sensitivity ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=300)

    names = list(CL_SCENARIOS.keys())
    cl_vals = [f"Cl={CL_SCENARIOS[n]:.1%}" for n in names]

    # (a) Paired bars
    ax = axes[0]
    x = np.arange(len(names))
    w = 0.35
    s_rates = [results[n]['rate_S'] * 100 for n in names]
    c_rates = [results[n]['rate_C'] * 100 for n in names]
    ax.bar(x - w/2, s_rates, w, color='tab:blue', alpha=0.7,
           label='Saueressig', edgecolor='black')
    ax.bar(x + w/2, c_rates, w, color='tab:red', alpha=0.7,
           label='Cantrell', edgecolor='black')
    ax.set_xticks(x)
    ax.set_xticklabels(cl_vals)
    ax.set_ylabel('Agreement Rate (%)')
    ax.set_title(f'(a) Agreement Rate at T={THRESHOLD} Tg/yr')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, axis='y')

    # (b) Discriminant with CI
    ax = axes[1]
    discs = [results[n]['discriminant_pp'] for n in names]
    ci_lo = [results[n]['discriminant_CI_pp'][0] for n in names]
    ci_hi = [results[n]['discriminant_CI_pp'][1] for n in names]
    err_lo = [d - lo for d, lo in zip(discs, ci_lo)]
    err_hi = [hi - d for d, hi in zip(discs, ci_hi)]
    sigs = [results[n]['significant'] for n in names]
    colors = ['green' if s else 'gray' for s in sigs]
    ax.bar(x, discs, color=colors, alpha=0.7, edgecolor='black',
           yerr=[err_lo, err_hi], capsize=5)
    ax.axhline(0, color='black', lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(cl_vals)
    ax.set_ylabel('Discriminant Δ (pp)')
    ax.set_title('(b) Discriminant + 95% CI (green = sig)')
    ax.grid(alpha=0.3, axis='y')
    for xi, d, s in zip(x, discs, sigs):
        marker = '***' if s else 'n.s.'
        ax.text(xi, d + 1.5, f'{d:.1f}\n{marker}', ha='center', fontsize=9, fontweight='bold')

    plt.suptitle('Phase 9b: Cl Fraction Sensitivity of Agreement-Filter Discriminant',
                 fontsize=12)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig16_cl_sensitivity.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"\n  Figure 16 saved.")

    return results


# ===================================================================
# TASK 3 — Year-agreement fraction sweep
# ===================================================================

def task_year_agree_sweep():
    N_ITER = 5000
    SEED = 42
    THRESHOLD = 90
    YEAR_FRACS = [0.60, 0.70, 0.80, 0.90, 0.95]

    print("\n" + "=" * 60)
    print("  Task 3: Year-agreement fraction sweep")
    print("=" * 60)

    # Run inversions once (re-use across year_agree_frac values)
    runs = {}
    for oh in ['saueressig', 'cantrell']:
        print(f"    Running {oh} (N={N_ITER})...")
        FF_c, FF_d, Mic_c, yrs = run_inversions(oh, n_iter=N_ITER, seed=SEED)
        runs[oh] = {'FF_c13': FF_c, 'FF_dD': FF_d, 'Mic_c13': Mic_c, 'years': yrs}

    years = runs['saueressig']['years']

    results = {}
    for frac in YEAR_FRACS:
        s_m = agreement_metrics(runs['saueressig']['FF_c13'], runs['saueressig']['FF_dD'],
                                runs['saueressig']['Mic_c13'], years, THRESHOLD,
                                year_agree_frac=frac)
        c_m = agreement_metrics(runs['cantrell']['FF_c13'], runs['cantrell']['FF_dD'],
                                runs['cantrell']['Mic_c13'], years, THRESHOLD,
                                year_agree_frac=frac)

        disc = (c_m['rate'] - s_m['rate']) * 100  # always same overall rate (doesn't depend on year_agree)

        results[str(frac)] = {
            'year_agree_frac': frac,
            'n_good_S': s_m['n_good'],
            'n_good_C': c_m['n_good'],
            'ff_trend_S': s_m['ff_trend'],
            'ff_trend_C': c_m['ff_trend'],
            'discriminant_pp': float(disc),
        }
        # The overall agreement rate doesn't change with year_agree_frac
        # (it only affects which iterations are *kept* for trend estimation).
        # So discriminant_pp is constant — what changes is n_good and trend estimates.
        # We also compute KSR at each frac.
        if not (np.isnan(s_m['ff_trend']) or np.isnan(c_m['ff_trend'])):
            spread = abs(c_m['ff_trend'] - s_m['ff_trend'])
            # Need baseline spread (unfiltered)
            s_unf = agreement_metrics(runs['saueressig']['FF_c13'], runs['saueressig']['FF_dD'],
                                      runs['saueressig']['Mic_c13'], years, 1e6, frac)
            c_unf = agreement_metrics(runs['cantrell']['FF_c13'], runs['cantrell']['FF_dD'],
                                      runs['cantrell']['Mic_c13'], years, 1e6, frac)
            bl = abs(c_unf['ff_trend'] - s_unf['ff_trend'])
            ksr = bl / max(spread, 1e-6)
            results[str(frac)]['KSR'] = float(ksr)
        else:
            results[str(frac)]['KSR'] = None

        print(f"    frac={frac:.0%}: n_good(S)={s_m['n_good']}, n_good(C)={c_m['n_good']}, "
              f"KSR={results[str(frac)].get('KSR', 'N/A')}")

    with open(OUT_DIR / "year_agreement_sweep.json", 'w') as f:
        json.dump(results, f, indent=2)

    # --- Figure 17: Year-agreement sweep ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=300)

    fracs = YEAR_FRACS
    labels = [f'{f:.0%}' for f in fracs]

    # (a) n_good iterations retained
    ax = axes[0]
    x = np.arange(len(fracs))
    w = 0.35
    s_ng = [results[str(f)]['n_good_S'] for f in fracs]
    c_ng = [results[str(f)]['n_good_C'] for f in fracs]
    ax.bar(x - w/2, s_ng, w, color='tab:blue', alpha=0.7,
           label='Saueressig', edgecolor='black')
    ax.bar(x + w/2, c_ng, w, color='tab:red', alpha=0.7,
           label='Cantrell', edgecolor='black')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlabel('Year-Agreement Fraction')
    ax.set_ylabel('Retained Iterations (of 5000)')
    ax.set_title(f'(a) Sample Size at T={THRESHOLD} Tg/yr')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, axis='y')

    # (b) KSR across year-agreement fracs
    ax = axes[1]
    ksrs = [results[str(f)]['KSR'] for f in fracs]
    valid_ksrs = [(f, k) for f, k in zip(fracs, ksrs) if k is not None]
    if valid_ksrs:
        vf, vk = zip(*valid_ksrs)
        ax.plot([f'{f:.0%}' for f in vf], vk, 'go-', lw=2, markersize=8)
        ax.axhline(1.0, color='gray', ls='--', lw=1, label='KSR=1')
        ax.set_ylabel('KSR (FF)')
        ax.set_xlabel('Year-Agreement Fraction')
        ax.set_title('(b) KSR vs Year-Agreement Threshold')
        ax.legend(fontsize=9)
        ax.grid(alpha=0.3)
        for fi, ki in zip(vf, vk):
            ax.text(f'{fi:.0%}', ki + 0.05, f'{ki:.2f}', ha='center', fontsize=9)
    else:
        ax.text(0.5, 0.5, 'Insufficient data', ha='center', va='center',
                transform=ax.transAxes)

    plt.suptitle('Phase 9c: Sensitivity to Year-Agreement Fraction',
                 fontsize=12)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig17_year_agree_sweep.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"\n  Figure 17 saved.")

    return results


# ===================================================================
# MAIN
# ===================================================================

def main():
    print("=" * 60)
    print("  Phase 9: Editorial Assessment Fixes")
    print("  (N=5000, Cl sensitivity, year-agree sweep)")
    print("=" * 60)

    t1 = task_high_n()
    t2 = task_cl_sensitivity()
    t3 = task_year_agree_sweep()

    print("\n" + "=" * 60)
    print("  Phase 9 complete. Summary:")
    print("=" * 60)
    print(f"\n  Task 1 (High-N): See {OUT_DIR / 'high_n_summary.json'}")
    print(f"  Task 2 (Cl):     See {OUT_DIR / 'cl_sensitivity.json'}")
    print(f"  Task 3 (Year%):  See {OUT_DIR / 'year_agreement_sweep.json'}")
    print(f"\n  Figures: fig15, fig16, fig17 in {FIG_DIR}")
    print("\n✓ Phase 9 complete.")


if __name__ == "__main__":
    main()
