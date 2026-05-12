#!/usr/bin/env python3
"""Phase 9 — Bootstrap confidence intervals on variance decomposition."""
import json
import numpy as np
from core import run_2box_flex, load_data, REPO_ROOT, OUT_DIR
from common import smooth_5yr

def bootstrap_variance(FF_total, FF_fix_kie, FF_fix_sigs, FF_fix_tau,
                       n_boot=1000, seed=123):
    """Bootstrap-resample MC iterations to get CIs on variance components."""
    rng = np.random.default_rng(seed)
    n_iter = FF_total.shape[1]

    sigmas = []
    kie_pcts = []
    sig_pcts = []
    tau_pcts = []
    resid_pcts = []

    for _ in range(n_boot):
        idx = rng.choice(n_iter, size=n_iter, replace=True)

        v_tot = np.nanmean(np.nanvar(smooth_5yr(FF_total[:, idx])[8:], axis=1))
        v_fk  = np.nanmean(np.nanvar(smooth_5yr(FF_fix_kie[:, idx])[8:], axis=1))
        v_fs  = np.nanmean(np.nanvar(smooth_5yr(FF_fix_sigs[:, idx])[8:], axis=1))
        v_ft  = np.nanmean(np.nanvar(smooth_5yr(FF_fix_tau[:, idx])[8:], axis=1))

        if v_tot < 1e-9:
            continue

        kie = max(0, v_tot - v_fk)
        sig = max(0, v_tot - v_fs)
        tau = max(0, v_tot - v_ft)
        total = kie + sig + tau + max(0, v_tot - kie - sig - tau)

        sigmas.append(np.sqrt(v_tot))
        kie_pcts.append(kie/total*100 if total > 0 else 0)
        sig_pcts.append(sig/total*100 if total > 0 else 0)
        tau_pcts.append(tau/total*100 if total > 0 else 0)
        resid_pcts.append(max(0, v_tot - kie - sig - tau)/total*100 if total > 0 else 0)

    def ci(arr):
        return float(np.median(arr)), float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))

    return {
        'sigma_ff': ci(sigmas),
        'kie_pct': ci(kie_pcts),
        'sig_pct': ci(sig_pcts),
        'tau_pct': ci(tau_pcts),
        'resid_pct': ci(resid_pcts),
    }


def main():
    print("=" * 60)
    print("PHASE 9: Bootstrap CIs on variance decomposition")
    print("=" * 60)

    data = load_data(REPO_ROOT, two_box=True)
    N = 400

    configs = {
        'd13C_only': dict(tau_mode="varying"),
        'dual_offset': dict(tau_mode="varying"),
        'dual_real_hemi': dict(tau_mode="varying"),
    }

    results = {}

    for label in ['d13C_only', 'dual_offset', 'dual_real_hemi']:
        print(f"\n  Running {label}...")

        if label == 'd13C_only':
            # Need d13C-only runs — import from variance_decomposition
            from variance_decomposition import run_2box
            FF_tot = run_2box(data, "d13C_only", N, 42, use_real_hemi_dD=False)
            FF_fk = run_2box(data, "d13C_only", N, 42, fix_kie=True, use_real_hemi_dD=False)
            FF_fs = run_2box(data, "d13C_only", N, 42, fix_sigs=True, use_real_hemi_dD=False)
            FF_ft = run_2box(data, "d13C_only", N, 42, fix_tau=True, use_real_hemi_dD=False)
        elif label == 'dual_offset':
            from variance_decomposition import run_2box
            FF_tot = run_2box(data, "dual", N, 42, use_real_hemi_dD=False)
            FF_fk = run_2box(data, "dual", N, 42, fix_kie=True, use_real_hemi_dD=False)
            FF_fs = run_2box(data, "dual", N, 42, fix_sigs=True, use_real_hemi_dD=False)
            FF_ft = run_2box(data, "dual", N, 42, fix_tau=True, use_real_hemi_dD=False)
        else:
            FF_tot = run_2box_flex(data, N, 42)
            FF_fk = run_2box_flex(data, N, 42, fix_kie=True)
            FF_fs = run_2box_flex(data, N, 42, fix_sigs=True)
            FF_ft = run_2box_flex(data, N, 42, tau_mode="fixed", tau_fixed=9.0)

        boot = bootstrap_variance(FF_tot, FF_fk, FF_fs, FF_ft)

        print(f"    σ(FF): {boot['sigma_ff'][0]:.1f} [{boot['sigma_ff'][1]:.1f}, {boot['sigma_ff'][2]:.1f}]")
        print(f"    KIE%:  {boot['kie_pct'][0]:.1f} [{boot['kie_pct'][1]:.1f}, {boot['kie_pct'][2]:.1f}]")
        print(f"    Sig%:  {boot['sig_pct'][0]:.1f} [{boot['sig_pct'][1]:.1f}, {boot['sig_pct'][2]:.1f}]")
        print(f"    τ%:    {boot['tau_pct'][0]:.1f} [{boot['tau_pct'][1]:.1f}, {boot['tau_pct'][2]:.1f}]")
        print(f"    Resid: {boot['resid_pct'][0]:.1f} [{boot['resid_pct'][1]:.1f}, {boot['resid_pct'][2]:.1f}]")

        results[label] = boot

    # Key comparison
    print(f"\n{'=' * 60}")
    print("KEY TEST: Is real-hemi Sig% significantly > offset Sig%?")
    print(f"{'=' * 60}")
    off_sig_hi = results['dual_offset']['sig_pct'][2]
    real_sig_lo = results['dual_real_hemi']['sig_pct'][1]
    if real_sig_lo > off_sig_hi:
        print(f"  ✓ YES: real-hemi Sig% [{results['dual_real_hemi']['sig_pct'][1]:.1f}, {results['dual_real_hemi']['sig_pct'][2]:.1f}]")
        print(f"         offset   Sig% [{results['dual_offset']['sig_pct'][1]:.1f}, {results['dual_offset']['sig_pct'][2]:.1f}]")
        print(f"         CIs don't overlap → statistically significant difference")
    else:
        print(f"  ⚠ CIs overlap: real [{results['dual_real_hemi']['sig_pct'][1]:.1f}, {results['dual_real_hemi']['sig_pct'][2]:.1f}] vs offset [{results['dual_offset']['sig_pct'][1]:.1f}, {results['dual_offset']['sig_pct'][2]:.1f}]")

    with open(OUT_DIR / "phase9_bootstrap.json", 'w') as f:
        json.dump(results, f, indent=2, default=float)
    print(f"\n  Saved: {OUT_DIR / 'phase9_bootstrap.json'}")


if __name__ == "__main__":
    main()
