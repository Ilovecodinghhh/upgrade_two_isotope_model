#!/usr/bin/env python3
"""
phase17_seed_sensitivity.py — Seed sensitivity (Review B5/2.2.4)
=================================================================
Tests robustness to RNG seed choice.
"""
import json
import numpy as np
from core import run_2box_flex, trend_stats, sigma_ff, load_data, REPO_ROOT, OUT_DIR

def main():
    print("=" * 60)
    print("PHASE 17: Seed sensitivity (Review issue B5)")
    print("=" * 60)

    data = load_data(REPO_ROOT, two_box=True)
    N = 400

    seeds = [42, 123, 314, 777, 2024]

    results = {}
    print(f"\n{'Seed':>6} {'σ(FF)':>7} {'ΔFF med':>8} {'[5%':>7} {'95%]':>7}")
    print("-" * 42)

    for seed in seeds:
        FF = run_2box_flex(data, N, seed)
        sig = sigma_ff(FF)
        med, lo, hi = trend_stats(FF, data.model_years)
        print(f"{seed:>6} {sig:>6.1f}  {med:>+7.1f}  {lo:>+6.1f}  {hi:>+6.1f}")

        results[str(seed)] = {
            'seed': seed,
            'sigma_ff': sig,
            'trend_median': med, 'trend_5pct': lo, 'trend_95pct': hi,
        }

    sigmas = [r['sigma_ff'] for r in results.values()]
    trends = [r['trend_median'] for r in results.values()]
    print(f"\n  σ(FF) range: {min(sigmas):.1f} – {max(sigmas):.1f} (spread {max(sigmas)-min(sigmas):.1f})")
    print(f"  ΔFF range:   {min(trends):+.1f} – {max(trends):+.1f} (spread {max(trends)-min(trends):.1f})")

    with open(OUT_DIR / "phase17_seed_sensitivity.json", 'w') as f:
        json.dump(results, f, indent=2, default=float)
    print(f"\n  Saved: {OUT_DIR / 'phase17_seed_sensitivity.json'}")


if __name__ == "__main__":
    main()
