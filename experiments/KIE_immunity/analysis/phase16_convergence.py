#!/usr/bin/env python3
"""
phase16_convergence.py — MC convergence analysis (Review B5)
=============================================================
Tests whether 400 iterations is sufficient by running at multiple N.
"""
import json
import numpy as np
from core import run_2box_flex, trend_stats, sigma_ff, load_data, REPO_ROOT, OUT_DIR
from common import smooth_5yr

def main():
    print("=" * 60)
    print("PHASE 16: MC convergence (Review issue B5)")
    print("=" * 60)

    data = load_data(REPO_ROOT, two_box=True)

    n_values = [50, 100, 200, 400, 600, 800, 1000]

    results = {}
    print(f"\n{'N_iter':>7} {'σ(FF)':>7} {'ΔFF med':>8} {'[5%':>7} {'95%]':>7}")
    print("-" * 42)

    for N in n_values:
        FF = run_2box_flex(data, N, 42)
        sig = sigma_ff(FF)
        med, lo, hi = trend_stats(FF, data.model_years)

        # Also variance decomposition at this N
        var_total = np.nanmean(np.nanvar(smooth_5yr(FF)[8:], axis=1))
        FF_fk = run_2box_flex(data, N, 42, fix_kie=True)
        var_fk = np.nanmean(np.nanvar(smooth_5yr(FF_fk)[8:], axis=1))
        kie_pct = max(0, (var_total - var_fk) / var_total * 100) if var_total > 0 else 0

        print(f"{N:>7} {sig:>6.1f}  {med:>+7.1f}  {lo:>+6.1f}  {hi:>+6.1f}")

        results[str(N)] = {
            'n_iter': N,
            'sigma_ff': sig,
            'trend_median': med, 'trend_5pct': lo, 'trend_95pct': hi,
            'kie_pct': kie_pct,
        }

    # Check convergence: is 400 within 5% of 1000?
    sig_400 = results['400']['sigma_ff']
    sig_1000 = results['1000']['sigma_ff']
    pct_diff = abs(sig_400 - sig_1000) / sig_1000 * 100
    print(f"\n  σ(FF) at N=400 vs N=1000: {sig_400:.1f} vs {sig_1000:.1f} ({pct_diff:.1f}% diff)")
    if pct_diff < 5:
        print("  ✓ N=400 converged (within 5% of N=1000)")
    else:
        print("  ⚠ N=400 may not be converged — consider using N=1000")

    with open(OUT_DIR / "phase16_convergence.json", 'w') as f:
        json.dump(results, f, indent=2, default=float)
    print(f"\n  Saved: {OUT_DIR / 'phase16_convergence.json'}")


if __name__ == "__main__":
    main()
