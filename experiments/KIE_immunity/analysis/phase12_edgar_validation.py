#!/usr/bin/env python3
"""Phase 12 — Compare FF estimates with bottom-up inventories."""
import sys
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
ANALYSIS_DIR = Path(__file__).resolve().parent.parent / "analysis"
FIG_DIR = Path(__file__).resolve().parent.parent / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(ANALYSIS_DIR))

from common import load_data, smooth_5yr
from core import run_2box_flex, OUT_DIR


def load_edgar_ff(repo_root):
    """Load EDGAR 8.0 Coal + ONG totals (Tg/yr)."""
    base = repo_root / "ImportantReferences/Riddell-Young2025PNAS_DS/Riddell-Young_2025_MassBalancePackage/Riddell-Young_2025_MassBalancePackage/data"
    
    coal = pd.read_csv(base / "EDGAR8_Coal.csv")
    ong = pd.read_csv(base / "EDGAR8_ONG.csv")
    
    # Sum all countries per year, convert kt → Tg
    year_cols = [c for c in coal.columns if c.isdigit() and int(c) >= 1997]
    
    coal_total = coal[year_cols].sum(axis=0) / 1000  # kt → Tg
    ong_total = ong[year_cols].sum(axis=0) / 1000
    
    years = np.array([int(y) for y in year_cols])
    ff_total = coal_total.values + ong_total.values
    
    return years, ff_total, coal_total.values, ong_total.values


def load_carbontracker_ff(repo_root):
    """Load CarbonTracker CH4 FF posterior."""
    df = pd.read_excel(repo_root / "rel/data/CarbonTracker_CH4.xlsx")
    return df['Date'].values, df['FF Post '].values, df['unc.1'].values


def main():
    print("=" * 60)
    print("PHASE 12: EDGAR / CarbonTracker validation")
    print("=" * 60)
    
    data = load_data(REPO_ROOT, two_box=True)
    N = 400; SEED = 42
    
    # Our model
    print("  Running dual real-hemi model...")
    FF = run_2box_flex(data, N, SEED)
    FF_s = smooth_5yr(FF)
    years_model = data.model_years
    ff_med = np.nanmedian(FF_s, axis=1)
    ff_lo = np.nanpercentile(FF_s, 5, axis=1)
    ff_hi = np.nanpercentile(FF_s, 95, axis=1)
    
    # EDGAR
    print("  Loading EDGAR 8.0...")
    years_edgar, edgar_ff, edgar_coal, edgar_ong = load_edgar_ff(REPO_ROOT)
    
    # CarbonTracker
    print("  Loading CarbonTracker CH4...")
    years_ct, ct_ff, ct_unc = load_carbontracker_ff(REPO_ROOT)
    
    # Trend analysis: post-2007 change
    def post2007_trend(years, values):
        mask_pre = (years >= 2000) & (years <= 2006)
        mask_post = (years >= 2010) & (years <= 2018)
        if mask_pre.sum() == 0 or mask_post.sum() == 0:
            return np.nan
        return np.mean(values[mask_post]) - np.mean(values[mask_pre])
    
    edgar_trend = post2007_trend(years_edgar, edgar_ff)
    ct_trend = post2007_trend(years_ct, ct_ff)
    
    print(f"\n  Post-2007 FF trends:")
    print(f"    Our model (real hemi):  {np.nanmedian(ff_med[np.isin(years_model, range(2010,2019))]) - np.nanmedian(ff_med[np.isin(years_model, range(2000,2007))]):.1f} Tg/yr")
    print(f"    EDGAR 8.0 (Coal+ONG):  {edgar_trend:+.1f} Tg/yr")
    print(f"    CarbonTracker FF post: {ct_trend:+.1f} Tg/yr")
    
    # ========== FIGURE ==========
    fig, ax = plt.subplots(figsize=(120/25.4, 80/25.4), dpi=300)
    plt.rcParams.update({'font.size': 8, 'axes.labelsize': 8})
    
    # Our model
    ax.plot(years_model, ff_med, color='#1b7837', linewidth=1.5, label='This study (dual real-hemi)', zorder=5)
    ax.fill_between(years_model, ff_lo, ff_hi, color='#1b7837', alpha=0.15, zorder=4)
    
    # EDGAR
    mask_e = (years_edgar >= 1997) & (years_edgar <= 2022)
    ax.plot(years_edgar[mask_e], edgar_ff[mask_e], color='#d95f02', linewidth=1.2,
            linestyle='--', marker='s', markersize=3, label='EDGAR 8.0 (Coal+ONG)', zorder=3)
    
    # CarbonTracker
    ax.plot(years_ct, ct_ff, color='#7570b3', linewidth=1.2,
            linestyle=':', marker='^', markersize=3, label='CarbonTracker CH₄ (posterior)', zorder=3)
    
    ax.axvline(2007, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.set_xlabel('Year')
    ax.set_ylabel('FF CH₄ emissions (Tg/yr)')
    ax.set_xlim(1998, 2022)
    ax.legend(loc='best', fontsize=7, framealpha=0.9)
    ax.set_title('Fossil-fuel methane: top-down vs bottom-up', fontsize=9, fontweight='bold')
    
    plt.tight_layout()
    fig.savefig(FIG_DIR / 'fig_edgar_validation.png', dpi=300, bbox_inches='tight')
    fig.savefig(FIG_DIR / 'fig_edgar_validation.pdf', bbox_inches='tight')
    print(f"\n  Saved: {FIG_DIR / 'fig_edgar_validation.png'}")
    plt.close()
    
    # Save results
    results = {
        'edgar_trend': edgar_trend,
        'ct_trend': ct_trend,
        'edgar_years': years_edgar.tolist(),
        'edgar_ff': edgar_ff.tolist(),
        'ct_years': years_ct.tolist(),
        'ct_ff': ct_ff.tolist(),
    }
    with open(OUT_DIR / "phase12_edgar.json", 'w') as f:
        json.dump(results, f, indent=2, default=float)
    print(f"  Saved: {OUT_DIR / 'phase12_edgar.json'}")


if __name__ == "__main__":
    main()
