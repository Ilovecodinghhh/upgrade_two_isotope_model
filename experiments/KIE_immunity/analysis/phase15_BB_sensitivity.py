#!/usr/bin/env python3
"""
phase15_BB_sensitivity.py — BB emission sensitivity (Review B2)
================================================================
Tests ±20% BB perturbation and ±10% NH/SH split perturbation.
"""
import json
import numpy as np
from core import run_2box_flex, trend_stats, sigma_ff, load_data, REPO_ROOT, OUT_DIR

def main():
    print("=" * 60)
    print("PHASE 15: BB sensitivity (Review issue B2)")
    print("=" * 60)

    data = load_data(REPO_ROOT, two_box=True)
    N = 400; SEED = 42

    configs = [
        ("BB −20%", 0.8),
        ("BB −10%", 0.9),
        ("BB baseline", 1.0),
        ("BB +10%", 1.1),
        ("BB +20%", 1.2),
    ]

    results = {}
    print(f"\n{'Config':<18} {'σ(FF)':>7} {'ΔFF med':>8} {'[5%':>7} {'95%]':>7}")
    print("-" * 50)

    for label, scale in configs:
        FF = run_2box_flex(data, N, SEED, bb_scale=scale)
        sig = sigma_ff(FF)
        med, lo, hi = trend_stats(FF, data.model_years)
        print(f"{label:<18} {sig:>6.1f}  {med:>+7.1f}  {lo:>+6.1f}  {hi:>+6.1f}")

        results[label] = {
            'bb_scale': scale,
            'sigma_ff': sig,
            'trend_median': med, 'trend_5pct': lo, 'trend_95pct': hi,
        }

    with open(OUT_DIR / "phase15_BB_sensitivity.json", 'w') as f:
        json.dump(results, f, indent=2, default=float)
    print(f"\n  Saved: {OUT_DIR / 'phase15_BB_sensitivity.json'}")


if __name__ == "__main__":
    main()
