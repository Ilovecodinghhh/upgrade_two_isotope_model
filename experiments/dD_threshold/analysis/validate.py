#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
validate.py — Smoke test for the refactored 2-box runner.
"""

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import load_data
from core import run_twobox, ci_width

N = 50
SEED = 42


def main():
    data = load_data(REPO_ROOT, two_box=True)
    passed = True

    # 1. Dual mode smoke test
    print("Running dual mode (N=50)...")
    res_dual = run_twobox(data, 1.0, N, SEED, mode="dual")
    ff_dual = res_dual['FF_G']

    # 2. d13C_only mode smoke test
    print("Running d13C_only mode (N=50)...")
    res_c13 = run_twobox(data, 1.0, N, SEED, mode="d13C_only")
    ff_c13 = res_c13['FF_G']

    # Check finite and positive (>0 for at least 90% of iterations)
    for label, ff in [("dual", ff_dual), ("d13C_only", ff_c13)]:
        finite_frac = np.isfinite(ff).mean()
        pos_frac = (ff > 0).mean()
        print(f"  {label}: finite={finite_frac:.2%}, positive={pos_frac:.2%}")
        if finite_frac < 0.90:
            print(f"  FAIL: {label} has <90% finite values")
            passed = False
        if pos_frac < 0.90:
            print(f"  FAIL: {label} has <90% positive values")
            passed = False

    # Check dual CI < d13C_only CI at multiplier=1
    ci_dual = ci_width(ff_dual)
    ci_c13 = ci_width(ff_c13)
    print(f"  CI widths: dual={ci_dual:.1f}, d13C_only={ci_c13:.1f}")
    if ci_dual >= ci_c13:
        print("  FAIL: dual CI should be < d13C_only CI at multiplier=1")
        passed = False

    if passed:
        print("\nPASS")
    else:
        print("\nFAIL")
        sys.exit(1)


if __name__ == "__main__":
    main()
