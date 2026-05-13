#!/usr/bin/env python3
"""
phase18_diagnostics.py — Solver diagnostics & regression trends (Review B4/B7/C4)
===================================================================================
Reports solver failure rate, bound-hit rate, and linear regression trends.
"""
import json
import numpy as np
from core import (run_2box_flex, trend_stats, sigma_ff, compute_trend_regression,
                  load_data, REPO_ROOT, OUT_DIR)

def main():
    print("=" * 60)
    print("PHASE 18: Solver diagnostics + regression trends")
    print("=" * 60)

    data = load_data(REPO_ROOT, two_box=True)
    N = 400; SEED = 42

    # Run with diagnostics
    FF, diag = run_2box_flex(data, N, SEED, track_diagnostics=True)

    print(f"\n  Solver diagnostics:")
    print(f"    Total solves:    {diag['total_solves']}")
    print(f"    Failures:        {diag['solver_failures']} ({diag['failure_rate_pct']:.2f}%)")
    print(f"    Bound hits:      {diag['bound_hits']} ({diag['bound_hit_rate_pct']:.1f}%)")

    # Step-change trend
    med, lo, hi = trend_stats(FF, data.model_years)
    print(f"\n  Step-change trend (mean(2010-2018) - mean(2000-2006)):")
    print(f"    ΔFF = {med:+.1f} [{lo:+.1f}, {hi:+.1f}] Tg/yr")

    # Linear regression trend
    reg = compute_trend_regression(FF, data.model_years)
    print(f"\n  Linear regression trend (2000–2020):")
    print(f"    Slope = {reg['slope_median']:+.2f} [{reg['slope_5pct']:+.2f}, {reg['slope_95pct']:+.2f}] Tg/yr²")
    print(f"    Median p-value = {reg['pvalue_median']:.3f}")
    print(f"    % iterations significant (p<0.05): {reg['pct_significant']:.1f}%")

    results = {
        'diagnostics': diag,
        'trend_step': {'median': med, '5pct': lo, '95pct': hi},
        'trend_regression': reg,
    }

    with open(OUT_DIR / "phase18_diagnostics.json", 'w') as f:
        json.dump(results, f, indent=2, default=float)
    print(f"\n  Saved: {OUT_DIR / 'phase18_diagnostics.json'}")


if __name__ == "__main__":
    main()
