#!/usr/bin/env python3
"""
run_all.py — Master rerun script for KIE_immunity experiment
==============================================================
Runs all phases sequentially to produce a single consistent dataset.
Addresses Review A1: all results from one definitive run.
"""
import subprocess
import sys
import json
import time
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parent
OUT_DIR = ANALYSIS_DIR.parent / "results"

# Order matters for dependencies
PHASES = [
    "variance_decomposition.py",  # Primary variance decomp (bug-fixed)
    "compare_basu2022.py",        # Basu KIE spread comparison (bug-fixed)
    "phase5_tau_sensitivity.py",
    "phase6_OHD_sensitivity.py",
    "phase7_Cl_sensitivity.py",
    "phase8_robustness_matrix.py",
    "phase9_bootstrap_variance.py",  # Depends on variance_decomposition
    "phase11_tau_ex.py",
    "phase12_edgar_validation.py",
    "phase14_W_sensitivity.py",      # NEW: Review A2
    "phase15_BB_sensitivity.py",     # NEW: Review B2
    "phase16_convergence.py",        # NEW: Review B5
    "phase17_seed_sensitivity.py",   # NEW: Review B5
    "phase18_diagnostics.py",        # NEW: Review B4/B7/C4
    "phase13_summary_table.py",      # Last: depends on phase9, basu
]


def run_phase(script):
    print(f"\n{'='*70}")
    print(f"  RUNNING: {script}")
    print(f"{'='*70}")
    t0 = time.time()
    result = subprocess.run(
        [sys.executable, str(ANALYSIS_DIR / script)],
        cwd=str(ANALYSIS_DIR),
        capture_output=False,
    )
    elapsed = time.time() - t0
    if result.returncode != 0:
        print(f"  ❌ FAILED: {script} (exit code {result.returncode})")
        return False
    print(f"  ✓ Done in {elapsed:.1f}s")
    return True


def main():
    print("="*70)
    print("  KIE_IMMUNITY: FULL RERUN (post-review v4)")
    print(f"  Output: {OUT_DIR}")
    print("="*70)

    t_start = time.time()
    failed = []

    for script in PHASES:
        if not run_phase(script):
            failed.append(script)

    elapsed = time.time() - t_start
    print(f"\n{'='*70}")
    print(f"  COMPLETE in {elapsed/60:.1f} minutes")
    if failed:
        print(f"  ❌ Failed phases: {', '.join(failed)}")
    else:
        print(f"  ✓ All {len(PHASES)} phases succeeded")

        # Write version marker
        version = {
            'version': 'v4-post-review',
            'timestamp': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
            'phases': len(PHASES),
            'changes': [
                'A3/B6: Fixed Strat_D (1.050→1.179) and Soil_D (1.103→1.083) in variance_decomposition.py and compare_basu2022.py',
                'A2: Added phase14_W_sensitivity.py',
                'B2: Added phase15_BB_sensitivity.py',
                'B4: Added compute_trend_regression() to core.py',
                'B5: Added phase16_convergence.py and phase17_seed_sensitivity.py',
                'B7/C4: Added solver diagnostics tracking to core.py + phase18_diagnostics.py',
                'core.py: W matrix now parameterized, bb_scale added',
            ],
        }
        with open(OUT_DIR / "version.json", 'w') as f:
            json.dump(version, f, indent=2)
        print(f"  Saved: {OUT_DIR / 'version.json'}")

    print("="*70)


if __name__ == "__main__":
    main()
