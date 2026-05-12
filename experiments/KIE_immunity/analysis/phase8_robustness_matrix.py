#!/usr/bin/env python3
"""Phase 8 — Combined robustness matrix: τ × OH_D × Cl."""
import json
import itertools
import numpy as np
from core import run_2box_flex, trend_stats, sigma_ff, load_data, REPO_ROOT, OUT_DIR

def main():
    print("=" * 70)
    print("PHASE 8: Combined robustness matrix (τ × OH_D × Cl)")
    print("=" * 70)

    data = load_data(REPO_ROOT, two_box=True)
    N = 400

    tau_values = [8.5, 9.0, 9.5]
    oh_d_values = [1.294, 1.327]
    cl_values = [0.006, 0.035, 0.065]

    results = {}
    n_total = len(tau_values) * len(oh_d_values) * len(cl_values)
    n_negative = 0
    n_robust = 0  # CI doesn't include zero

    print(f"\n{'τ':>5} {'OH_D':>6} {'Cl%':>5} │ {'ΔFF med':>8} {'[5%':>7} {'95%]':>7} {'σ(FF)':>7} │ Sign")
    print("─" * 70)

    for tau, oh_d, cl in itertools.product(tau_values, oh_d_values, cl_values):
        FF = run_2box_flex(data, N, 42,
                          tau_mode="fixed", tau_fixed=tau,
                          oh_d_fixed=oh_d, cl_frac=cl)
        med, lo, hi = trend_stats(FF, data.model_years)
        sig = sigma_ff(FF)

        key = f"τ={tau}_OHD={oh_d}_Cl={cl}"
        is_neg = med < 0
        ci_zero = lo < 0 < hi
        sign_str = "−" if is_neg else "+"
        robust_str = "✓" if (is_neg and not ci_zero) else ("≈0" if ci_zero else "⚠+")

        if is_neg:
            n_negative += 1
        if is_neg and not ci_zero:
            n_robust += 1

        print(f"{tau:>5.1f} {oh_d:>6.3f} {cl*100:>4.1f}% │ {med:>+7.1f}  {lo:>+6.1f}  {hi:>+6.1f}  {sig:>6.1f} │ {sign_str} {robust_str}")

        results[key] = {
            'tau': tau, 'oh_d': oh_d, 'cl_frac': cl,
            'trend_median': med, 'trend_5pct': lo, 'trend_95pct': hi,
            'sigma_ff': sig,
        }

    print(f"\n{'=' * 70}")
    print("ROBUSTNESS SUMMARY")
    print(f"{'=' * 70}")
    print(f"  Total cells:         {n_total}")
    print(f"  ΔFF negative:        {n_negative}/{n_total} ({n_negative/n_total*100:.0f}%)")
    print(f"  Robustly negative:   {n_robust}/{n_total} ({n_robust/n_total*100:.0f}%)")
    print(f"  CI includes zero:    {n_negative - n_robust}/{n_total}")
    print(f"  ΔFF positive:        {n_total - n_negative}/{n_total}")

    if n_negative == n_total:
        print("\n  ✓✓ FF TREND REVERSAL IS IRON-CLAD across all 18 cells")
    elif n_negative > n_total * 0.8:
        print(f"\n  ✓ FF trend reversal holds in {n_negative}/{n_total} cells — ROBUST with caveats")
    elif n_negative > n_total * 0.5:
        print(f"\n  ⚠ FF trend reversal holds in {n_negative}/{n_total} cells — CONDITIONAL")
    else:
        print(f"\n  ❌ FF trend reversal fails in majority of cells — NOT ROBUST")

    results['_summary'] = {
        'n_total': n_total,
        'n_negative': n_negative,
        'n_robust': n_robust,
        'pct_negative': n_negative/n_total*100,
        'pct_robust': n_robust/n_total*100,
    }

    with open(OUT_DIR / "phase8_robustness_matrix.json", 'w') as f:
        json.dump(results, f, indent=2, default=float)
    print(f"\n  Saved: {OUT_DIR / 'phase8_robustness_matrix.json'}")


if __name__ == "__main__":
    main()
