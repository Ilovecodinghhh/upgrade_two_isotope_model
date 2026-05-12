#!/usr/bin/env python3
"""Phase 13 — Summary table with bootstrap CIs for all headline numbers."""
import sys
import json
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
ANALYSIS_DIR = Path(__file__).resolve().parent.parent / "analysis"
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(ANALYSIS_DIR))

from common import load_data, smooth_5yr
from core import run_2box_flex, trend_stats, sigma_ff, OUT_DIR
from variance_decomposition import run_2box


def bootstrap_metric(metric_fn, FF, n_boot=1000, seed=123):
    """Bootstrap a metric function over MC iterations (axis=1)."""
    rng = np.random.default_rng(seed)
    n_iter = FF.shape[1]
    vals = []
    for _ in range(n_boot):
        idx = rng.choice(n_iter, size=n_iter, replace=True)
        vals.append(metric_fn(FF[:, idx]))
    return float(np.median(vals)), float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def main():
    print("=" * 60)
    print("PHASE 13: Summary table with bootstrap CIs")
    print("=" * 60)

    data = load_data(REPO_ROOT, two_box=True)
    N = 400; SEED = 42

    # Load pre-computed bootstrap results
    with open(OUT_DIR / "phase9_bootstrap.json") as f:
        boot = json.load(f)
    with open(OUT_DIR / "basu_comparison_v2.json") as f:
        basu = json.load(f)

    configs = [
        ('δ¹³C-only', 'd13C_only'),
        ('Dual (offset)', 'dual_offset'),
        ('Dual (real hemi)', 'dual_real_hemi'),
    ]

    # Run models for trend bootstrap
    print("  Running models for trend CIs...")
    FF_runs = {}
    FF_runs['d13C_only'] = run_2box(data, "d13C_only", N, SEED, use_real_hemi_dD=False)
    FF_runs['dual_offset'] = run_2box(data, "dual", N, SEED, use_real_hemi_dD=False)
    FF_runs['dual_real_hemi'] = run_2box_flex(data, N, SEED)

    years = data.model_years

    def trend_fn(FF):
        FF_s = smooth_5yr(FF)
        yr0 = int(years[0])
        pre = FF_s[max(0,2000-yr0):2007-yr0, :].mean(axis=0)
        post = FF_s[max(0,2010-yr0):2019-yr0, :].mean(axis=0)
        return float(np.nanmedian(post - pre))

    rows = []
    print(f"\n{'Config':<22} {'σ(FF)':>14} {'ΔFF trend':>18} {'KIE%':>14} {'Sig%':>14} {'KIE spread':>14}")
    print("─" * 100)

    for label, key in configs:
        # σ(FF) from bootstrap
        sig_med, sig_lo, sig_hi = boot[key]['sigma_ff']

        # ΔFF trend with bootstrap
        t_med, t_lo, t_hi = bootstrap_metric(trend_fn, FF_runs[key])

        # KIE% from bootstrap
        kie_med, kie_lo, kie_hi = boot[key]['kie_pct']

        # Sig% from bootstrap
        sig_pct_med, sig_pct_lo, sig_pct_hi = boot[key]['sig_pct']

        # KIE spread from Basu comparison
        basu_key = 'our_d13C_only' if key == 'd13C_only' else ('our_dual' if key == 'dual_offset' else 'our_dual_real_hemi')
        if basu_key in basu:
            kie_spread = basu[basu_key].get('kie_spread', None)
        else:
            kie_spread = None

        row = {
            'config': label,
            'sigma_ff': f"{sig_med:.1f} [{sig_lo:.1f}, {sig_hi:.1f}]",
            'dff_trend': f"{t_med:+.1f} [{t_lo:+.1f}, {t_hi:+.1f}]",
            'kie_pct': f"{kie_med:.1f} [{kie_lo:.1f}, {kie_hi:.1f}]",
            'sig_pct': f"{sig_pct_med:.1f} [{sig_pct_lo:.1f}, {sig_pct_hi:.1f}]",
            'kie_spread': f"{kie_spread:.1f}" if kie_spread else "—",
        }
        rows.append(row)

        print(f"{label:<22} {row['sigma_ff']:>14} {row['dff_trend']:>18} {row['kie_pct']:>14} {row['sig_pct']:>14} {row['kie_spread']:>14}")

    # Basu reference row
    print(f"{'Basu 2022 (3D)':<22} {'—':>14} {'—':>18} {'—':>14} {'—':>14} {'13.0':>14}")

    # LaTeX table
    latex = r"""\begin{table}[ht]
\centering
\caption{Headline results: fossil-fuel emission uncertainty across model configurations. 
Values show median [95\% bootstrap CI] from 1000 resamples of 400 MC iterations.}
\label{tab:headline}
\begin{tabular}{lccccc}
\toprule
Configuration & $\sigma(\mathrm{FF})$ & $\Delta\mathrm{FF}$ trend & KIE\% & Sig\% & KIE spread \\
 & (Tg/yr) & (Tg/yr) & & & (Tg/yr) \\
\midrule
"""
    for r in rows:
        latex += f"  {r['config']} & {r['sigma_ff']} & {r['dff_trend']} & {r['kie_pct']} & {r['sig_pct']} & {r['kie_spread']} \\\\\n"
    latex += r"""  \midrule
  Basu et al.\ (2022, 3D) & --- & --- & --- & --- & 13.0 \\
\bottomrule
\end{tabular}
\end{table}
"""

    with open(OUT_DIR / "table1.tex", 'w') as f:
        f.write(latex)
    print(f"\n  Saved: {OUT_DIR / 'table1.tex'}")

    # CSV
    import csv
    with open(OUT_DIR / "table1.csv", 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['config', 'sigma_ff', 'dff_trend', 'kie_pct', 'sig_pct', 'kie_spread'])
        w.writeheader()
        w.writerows(rows)
    print(f"  Saved: {OUT_DIR / 'table1.csv'}")

    # JSON
    with open(OUT_DIR / "phase13_summary.json", 'w') as f:
        json.dump(rows, f, indent=2)
    print(f"  Saved: {OUT_DIR / 'phase13_summary.json'}")


if __name__ == "__main__":
    main()
