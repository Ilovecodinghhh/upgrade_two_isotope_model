#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 6: Bayesian Agreement Framework
=======================================
Instead of WLS (hard coupling), uses δD as independent validation:
1. Solve δ¹³C-only → FF/Mic posterior
2. Solve δD-only → FF/Mic posterior
3. Compute intersection (agreement zone)
4. Test whether agreement zone reduces KIE sensitivity

This follows the approach of Riddell-Young (2025): solve independently,
then check consistency.
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
    smooth_5yr, trend_change, pad_to_length,
    SINK_FRACTIONS_GLOBAL, PT,
)

BASE = Path(__file__).resolve().parent
OUT_DIR = BASE / "results" / "phase6_bayesian"
FIG_DIR = BASE / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

N_ITER = 1000
SEED = 42


def solve_d13C_only(S, BB_j, d13C_src, sigs, j):
    """Analytic δ¹³C-only solution."""
    denom = sigs['ff_d13C'][j] - sigs['mic_d13C'][j]
    if abs(denom) > 0.1:
        ff = (S * d13C_src - sigs['mic_d13C'][j] * (S - BB_j) -
              sigs['bb_d13C'][j] * BB_j) / denom
        mic = S - BB_j - ff
        return ff, mic
    return np.nan, np.nan


def solve_dD_only(S, BB_j, dD_src, sigs, j):
    """Analytic δD-only solution."""
    denom = sigs['ff_dD'][j] - sigs['mic_dD'][j]
    if abs(denom) > 0.5:  # δD has larger separations, so larger threshold
        ff = (S * dD_src - sigs['mic_dD'][j] * (S - BB_j) -
              sigs['bb_dD'][j] * BB_j) / denom
        mic = S - BB_j - ff
        return ff, mic
    return np.nan, np.nan


def run_agreement(oh13c_mode: str, label: str):
    """Run independent δ¹³C and δD inversions + agreement filter."""
    print(f"\n  Phase 6 | OH_13C={oh13c_mode} ({label})")

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

    # Store all solutions
    FF_c13 = np.full((n, N_ITER), np.nan)
    Mic_c13 = np.full((n, N_ITER), np.nan)
    FF_dD = np.full((n, N_ITER), np.nan)
    Mic_dD = np.full((n, N_ITER), np.nan)
    FF_agree = np.full((n, N_ITER), np.nan)
    Mic_agree = np.full((n, N_ITER), np.nan)
    agreement_mask = np.zeros((n, N_ITER), dtype=bool)

    for k in range(N_ITER):
        if (k + 1) % 500 == 0:
            print(f"      iter {k+1}/{N_ITER}")

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

            # δ¹³C source delta
            d13C_src_f = (n13[j+1] - n13[j] + n13[j] * alpha_13C / tau[j]) / S
            d13C_src = fraction_to_delta_d13C(d13C_src_f)

            # δD source delta
            dD_src_f = (nD[j+1] - nD[j] + nD[j] * alpha_D / tau[j]) / S
            dD_src = fraction_to_delta_dD(dD_src_f)

            # Solve independently
            ff_c, mic_c = solve_d13C_only(S, BB_j, d13C_src, sigs, j)
            ff_d, mic_d = solve_dD_only(S, BB_j, dD_src, sigs, j)

            FF_c13[j, k] = ff_c
            Mic_c13[j, k] = mic_c
            FF_dD[j, k] = ff_d
            Mic_dD[j, k] = mic_d

            # Agreement check: do the two independent solutions agree?
            # "Agree" = FF estimates within 1σ of each other
            # Use a generous threshold: |FF_c13 - FF_dD| < tolerance
            # Tolerance = 50% of total source (generous, ~275 Tg)
            # OR more strict: both give same sign for trend direction
            if not np.isnan(ff_c) and not np.isnan(ff_d):
                # Agreement: FF from both methods within 100 Tg/yr of each other
                if abs(ff_c - ff_d) < 100:
                    agreement_mask[j, k] = True
                    # Use δ¹³C value (better constrained) but validated by δD
                    FF_agree[j, k] = ff_c
                    Mic_agree[j, k] = mic_c
                else:
                    # Disagreement: reject this iteration for this year
                    FF_agree[j, k] = np.nan
                    Mic_agree[j, k] = np.nan
            else:
                FF_agree[j, k] = ff_c  # fall back to δ¹³C when δD unavailable
                Mic_agree[j, k] = mic_c

    # Compute statistics for agreed iterations
    # Clamp negatives
    FF_c13 = np.clip(FF_c13, 0, None)
    Mic_c13 = np.clip(Mic_c13, 0, None)
    FF_agree_clean = np.where(agreement_mask | np.isnan(FF_agree), FF_agree, FF_c13)

    # For trend, use only δ¹³C (baseline) and agreement-filtered
    FF_c13_s = smooth_5yr(np.nan_to_num(FF_c13, nan=0))
    Mic_c13_s = smooth_5yr(np.nan_to_num(Mic_c13, nan=0))

    # Agreement-filtered: replace rejected iterations with NaN, then compute stats
    # For smoothing, we need to handle NaN — use nanmean approach
    agree_rate = agreement_mask.mean(axis=1)
    print(f"    Agreement rate: {agree_rate.mean()*100:.1f}% (range {agree_rate.min()*100:.0f}–{agree_rate.max()*100:.0f}%)")

    # For trend comparison, compute per-iteration:
    # Only use iterations where ALL years agree (strict)
    all_agree = agreement_mask.all(axis=0)  # iterations where every year agrees
    n_all_agree = all_agree.sum()
    print(f"    Iterations with ALL years agreeing: {n_all_agree}/{N_ITER} ({100*n_all_agree/N_ITER:.0f}%)")

    # Less strict: use per-year agreement, compute stats on available iterations
    # For each year, compute trend using only agreeing iterations
    # Simplified: just use the agreed FF values, replace NaN with δ¹³C solution
    FF_filtered = np.where(agreement_mask, FF_c13, np.nan)
    FF_filtered = np.clip(FF_filtered, 0, None)

    # Fill NaN with per-year median for smoothing
    FF_for_smooth = np.copy(FF_c13)
    for j in range(n):
        valid = agreement_mask[j, :]
        if valid.sum() > 50:
            # Replace non-agreeing with NaN to exclude from trend
            FF_for_smooth[j, ~valid] = np.nan

    # Compute trend only on agreeing iterations
    delta_ff_c13, pct_ff_c13 = trend_change(FF_c13_s, years)

    # For agreement-filtered, compute trend per iteration only if enough years agree
    # Use a relaxed approach: for each iteration, check if at least 80% of years agree
    good_iters = (agreement_mask.sum(axis=0) >= n * 0.8)
    n_good = good_iters.sum()
    print(f"    Iterations with ≥80% agreement: {n_good}/{N_ITER}")

    if n_good > 50:
        FF_agree_s = smooth_5yr(np.clip(FF_c13[:, good_iters], 0, None))
        Mic_agree_s = smooth_5yr(np.clip(Mic_c13[:, good_iters], 0, None))
        delta_ff_agree, pct_ff_agree = trend_change(FF_agree_s, years)
        delta_mic_agree, pct_mic_agree = trend_change(Mic_agree_s, years)
    else:
        delta_ff_agree = delta_ff_c13
        delta_mic_agree = trend_change(Mic_c13_s, years)[0]

    delta_mic_c13, pct_mic_c13 = trend_change(Mic_c13_s, years)

    print(f"    δ¹³C-only: FF Δ={delta_ff_c13.mean():+.1f}±{delta_ff_c13.std():.1f}, "
          f"Mic Δ={delta_mic_c13.mean():+.1f}±{delta_mic_c13.std():.1f}")
    if n_good > 50:
        print(f"    Agreement-filtered: FF Δ={delta_ff_agree.mean():+.1f}±{delta_ff_agree.std():.1f}, "
              f"Mic Δ={delta_mic_agree.mean():+.1f}±{delta_mic_agree.std():.1f}")

    # Save
    np.savez(OUT_DIR / f"run_{label}.npz",
             FF_c13=FF_c13_s, Mic_c13=Mic_c13_s,
             years=years,
             delta_ff_c13=delta_ff_c13, delta_mic_c13=delta_mic_c13,
             delta_ff_agree=delta_ff_agree, delta_mic_agree=delta_mic_agree,
             agreement_rate=agree_rate,
             n_good_iters=n_good)

    return {
        'label': label,
        'oh13c_mode': oh13c_mode,
        'delta_ff_c13': delta_ff_c13,
        'delta_mic_c13': delta_mic_c13,
        'delta_ff_agree': delta_ff_agree,
        'delta_mic_agree': delta_mic_agree,
        'agree_rate': float(agree_rate.mean()),
        'n_good_iters': int(n_good),
        'ff_std_c13': float(delta_ff_c13.std()),
        'ff_std_agree': float(delta_ff_agree.std()) if n_good > 50 else None,
    }


def main():
    print("=" * 60)
    print("  Phase 6 — Bayesian Agreement Framework")
    print("=" * 60)

    results = {}
    for oh13c, lbl in [('saueressig', 'A_saueressig'),
                        ('cantrell', 'B_cantrell'),
                        ('sampled', 'C_sampled')]:
        results[lbl] = run_agreement(oh13c, lbl)

    # KSR for agreement-filtered
    spread_c13_ff = abs(results['B_cantrell']['delta_ff_c13'].mean() -
                        results['A_saueressig']['delta_ff_c13'].mean())
    spread_agree_ff = abs(results['B_cantrell']['delta_ff_agree'].mean() -
                          results['A_saueressig']['delta_ff_agree'].mean())
    ksr_agree_ff = spread_c13_ff / max(spread_agree_ff, 1e-6)

    spread_c13_mic = abs(results['B_cantrell']['delta_mic_c13'].mean() -
                         results['A_saueressig']['delta_mic_c13'].mean())
    spread_agree_mic = abs(results['B_cantrell']['delta_mic_agree'].mean() -
                           results['A_saueressig']['delta_mic_agree'].mean())
    ksr_agree_mic = spread_c13_mic / max(spread_agree_mic, 1e-6)

    print(f"\n  === AGREEMENT KSR ===")
    print(f"  KSR_agree (FF):  {ksr_agree_ff:.2f}  (spreads: c13={spread_c13_ff:.2f}, agree={spread_agree_ff:.2f})")
    print(f"  KSR_agree (Mic): {ksr_agree_mic:.2f}  (spreads: c13={spread_c13_mic:.2f}, agree={spread_agree_mic:.2f})")

    # Uncertainty reduction from filtering
    std_c13_ff = results['C_sampled']['ff_std_c13']
    std_agree_ff = results['C_sampled']['ff_std_agree']
    if std_agree_ff:
        red_ff = (1 - std_agree_ff / std_c13_ff) * 100
        print(f"\n  Uncertainty reduction (FF): {red_ff:.1f}%")

    # Save summary
    summary = {
        'KSR_agree_FF': float(ksr_agree_ff),
        'KSR_agree_Mic': float(ksr_agree_mic),
        'spread_c13_FF': float(spread_c13_ff),
        'spread_agree_FF': float(spread_agree_ff),
        'spread_c13_Mic': float(spread_c13_mic),
        'spread_agree_Mic': float(spread_agree_mic),
        'agreement_rates': {k: v['agree_rate'] for k, v in results.items()},
        'n_good_iters': {k: v['n_good_iters'] for k, v in results.items()},
    }
    with open(OUT_DIR / "summary.json", 'w') as f:
        json.dump(summary, f, indent=2)

    # === Figure 8: Agreement framework results ===
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), dpi=300)

    # (a) FF trends: δ¹³C-only for Saueressig vs Cantrell
    ax = axes[0, 0]
    ax.hist(results['A_saueressig']['delta_ff_c13'], bins=40, alpha=0.5,
            color='tab:blue', density=True, label='Saueressig (δ¹³C)')
    ax.hist(results['B_cantrell']['delta_ff_c13'], bins=40, alpha=0.5,
            color='tab:red', density=True, label='Cantrell (δ¹³C)')
    ax.axvline(0, color='black', ls='--', lw=0.8)
    ax.set_title(f'(a) FF Trends — δ¹³C-only\nspread = {spread_c13_ff:.2f} Tg/yr')
    ax.set_xlabel('Δ FF Emissions (Tg/yr)')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # (b) FF trends: agreement-filtered
    ax = axes[0, 1]
    ax.hist(results['A_saueressig']['delta_ff_agree'], bins=40, alpha=0.5,
            color='tab:blue', density=True, label='Saueressig (filtered)')
    ax.hist(results['B_cantrell']['delta_ff_agree'], bins=40, alpha=0.5,
            color='tab:red', density=True, label='Cantrell (filtered)')
    ax.axvline(0, color='black', ls='--', lw=0.8)
    ax.set_title(f'(b) FF Trends — Agreement-Filtered\nspread = {spread_agree_ff:.2f} Tg/yr, KSR = {ksr_agree_ff:.2f}')
    ax.set_xlabel('Δ FF Emissions (Tg/yr)')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # (c) Mic trends: δ¹³C-only
    ax = axes[1, 0]
    ax.hist(results['A_saueressig']['delta_mic_c13'], bins=40, alpha=0.5,
            color='tab:blue', density=True, label='Saueressig (δ¹³C)')
    ax.hist(results['B_cantrell']['delta_mic_c13'], bins=40, alpha=0.5,
            color='tab:red', density=True, label='Cantrell (δ¹³C)')
    ax.axvline(0, color='black', ls='--', lw=0.8)
    ax.set_title(f'(c) Mic Trends — δ¹³C-only\nspread = {spread_c13_mic:.2f} Tg/yr')
    ax.set_xlabel('Δ Mic Emissions (Tg/yr)')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # (d) Agreement rate over time
    ax = axes[1, 1]
    for lbl, color, name in [('A_saueressig', 'tab:blue', 'Saueressig'),
                              ('B_cantrell', 'tab:red', 'Cantrell'),
                              ('C_sampled', 'gray', 'Sampled')]:
        dat = np.load(OUT_DIR / f"run_{lbl}.npz")
        ax.plot(dat['years'], dat['agreement_rate'] * 100, '-', color=color,
                lw=2, label=name)
    ax.set_xlabel('Year')
    ax.set_ylabel('Agreement Rate (%)')
    ax.set_title('(d) δ¹³C–δD Agreement Rate Over Time')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax.set_ylim(0, 100)

    plt.suptitle(f'Phase 6: Bayesian Agreement Framework\nKSR(FF)={ksr_agree_ff:.2f}, KSR(Mic)={ksr_agree_mic:.2f}',
                 fontsize=12)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig8_agreement_framework.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"\n  Figure 8 saved: {FIG_DIR / 'fig8_agreement_framework.png'}")
    print("\n✓ Phase 6 complete.")


if __name__ == "__main__":
    main()
