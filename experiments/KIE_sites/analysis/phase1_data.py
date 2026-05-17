#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
phase1_data.py — Data extraction and pairing for the KIE_sites experiment
==========================================================================

Loads event-level δ¹³C (NOAA/INSTAAR) and δD (Riddell-Young 2025 compilation)
data for all co-located sites, pairs them by date, computes monthly means,
and produces summary plots and statistics.

Output:
  results/phase1_data/
    site_monthly_{CODE}.csv    — monthly means for each site
    site_summary.json          — per-site data quality summary
  figures/
    fig1_timeseries_grid.png   — δ¹³C and δD time series for all sites
    fig1_data_coverage.png     — temporal coverage heatmap
"""

import sys
from pathlib import Path
import json

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================================================================
# PATHS
# ============================================================================
# Script is at experiments/KIE_sites/analysis/phase1_data.py
# Repo root is 4 levels up: analysis → KIE_sites → experiments → repo_root
REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
SITES_DIR = REPO_ROOT / "sitesdata"
D13C_EVENT_DIR = SITES_DIR / "isotope_d13C" / "noaa_instaar_2023_event"
DD_RAW_DIR = SITES_DIR / "isotope_dD" / "raw_observations"
RESULTS_DIR = Path(__file__).resolve().parent.parent / "results" / "phase1_data"
FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================================
# CO-LOCATED SITE DEFINITIONS
# ============================================================================
# Each entry: site code for display, δ¹³C file code (INSTAAR), δD file codes
# (one or more labs), latitude, and a "clean MBL" flag for later filtering.
#
# δD file IDs follow the Riddell-Young naming: lowercase site + optional lab suffix.
# Multiple δD files at the same location (e.g. INSTAAR + MPI) are listed
# separately so we can analyse each lab independently or merge later.

COLOCATED_SITES = [
    {"code": "ALT", "lat":  82.45, "d13c_file": "alt", "dD_files": ["alt", "altMPI"],       "mbl": False},
    {"code": "ZEP", "lat":  78.91, "d13c_file": "zep", "dD_files": ["zep", "zepIMAU", "nyaNIPR"], "mbl": True},
    {"code": "BRW", "lat":  71.32, "d13c_file": "brw", "dD_files": ["brw", "brwIMAU"],       "mbl": True},
    {"code": "CBA", "lat":  55.21, "d13c_file": "cba", "dD_files": ["cba"],                  "mbl": True},
    {"code": "MHD", "lat":  53.33, "d13c_file": "mhd", "dD_files": ["mhd"],                  "mbl": True},
    {"code": "AZR", "lat":  38.77, "d13c_file": "azr", "dD_files": ["azr"],                  "mbl": False},
    {"code": "MLO", "lat":  19.54, "d13c_file": "mlo", "dD_files": ["mlo", "mloIMAU"],       "mbl": False},
    {"code": "KUM", "lat":  19.56, "d13c_file": "kum", "dD_files": ["kum"],                  "mbl": True},
    {"code": "ASC", "lat":  -7.97, "d13c_file": "asc", "dD_files": ["asc"],                  "mbl": True},
    {"code": "SMO", "lat": -14.25, "d13c_file": "smo", "dD_files": ["smo", "smoIMAU"],       "mbl": True},
    {"code": "CGO", "lat": -40.68, "d13c_file": "cgo", "dD_files": ["cgo", "cgoIMAU"],       "mbl": True},
    {"code": "SPO", "lat": -89.98, "d13c_file": "spo", "dD_files": ["spo"],                  "mbl": True},
]


# ============================================================================
# DATA LOADING
# ============================================================================

def load_d13C_event(site_code_lower: str) -> pd.DataFrame:
    """Load NOAA/INSTAAR δ¹³C event data for one site.

    Returns DataFrame with columns: decimal_year, year, month, d13C, d13C_unc, qcflag.
    Only keeps measurements with clean QC flags (first char = '.').
    """
    # Find the file — naming convention: ch4c13_{code}_surface-flask_7_sil_event.txt
    pattern = f"ch4c13_{site_code_lower}_*_event.txt"
    matches = list(D13C_EVENT_DIR.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"No δ¹³C file matching {pattern} in {D13C_EVENT_DIR}")
    fpath = matches[0]

    # Read header_lines count from first line
    with open(fpath) as f:
        first_line = f.readline()
    n_header = int(first_line.split(":")[-1].strip())

    # Read data, skipping the header block
    df = pd.read_csv(
        fpath, sep=r"\s+", skiprows=n_header,
        names=["site", "year", "month", "day", "hour", "minute", "second",
               "datetime", "decimal_year", "flask_id", "value", "value_unc",
               "lat", "lon", "alt", "elev", "intake_ht", "method",
               "event_num", "instrument", "analysis_dt", "qcflag"],
        dtype={"qcflag": str},
        na_values=["-999.999", "-999.99", "-9.999"],
    )

    # QC filter: keep only rows where the rejection flag (1st char) is '.'
    df = df[df["qcflag"].str[0] == "."].copy()

    # Keep relevant columns
    df = df.rename(columns={"value": "d13C", "value_unc": "d13C_unc"})
    df = df[["decimal_year", "year", "month", "d13C", "d13C_unc", "qcflag"]].copy()
    df = df.dropna(subset=["d13C"])
    df = df.sort_values("decimal_year").reset_index(drop=True)
    return df


def load_dD_event(site_id: str) -> pd.DataFrame:
    """Load δD event data from the Riddell-Young compilation for one site/lab.

    Returns DataFrame with columns: decimal_year, year, month, dD.
    """
    fpath = DD_RAW_DIR / f"{site_id}_01D0_dat.txt"
    if not fpath.exists():
        raise FileNotFoundError(f"No δD file: {fpath}")

    # Simple 2-column format: decimal_year  dD_permil
    df = pd.read_csv(fpath, sep=r"\s+", header=None, names=["decimal_year", "dD"])
    df = df.dropna(subset=["dD"])
    df["year"] = np.floor(df["decimal_year"]).astype(int)
    df["month"] = np.clip(
        ((df["decimal_year"] - df["year"]) * 12 + 1).astype(int), 1, 12
    )
    df = df.sort_values("decimal_year").reset_index(drop=True)
    return df


# ============================================================================
# MONTHLY AGGREGATION
# ============================================================================

def compute_monthly_means(d13c_df: pd.DataFrame,
                          dD_df: pd.DataFrame) -> pd.DataFrame:
    """Compute monthly means of δ¹³C and δD, then merge on (year, month).

    Returns a DataFrame with columns:
        year, month, decimal_year, d13C_mean, d13C_std, d13C_n,
        dD_mean, dD_std, dD_n
    """
    # δ¹³C monthly stats
    c13_monthly = (
        d13c_df.groupby(["year", "month"])
        .agg(d13C_mean=("d13C", "mean"),
             d13C_std=("d13C", "std"),
             d13C_n=("d13C", "count"))
        .reset_index()
    )

    # δD monthly stats
    dD_monthly = (
        dD_df.groupby(["year", "month"])
        .agg(dD_mean=("dD", "mean"),
             dD_std=("dD", "std"),
             dD_n=("dD", "count"))
        .reset_index()
    )

    # Merge on (year, month) — inner join keeps only months with both isotopes
    merged = pd.merge(c13_monthly, dD_monthly, on=["year", "month"], how="inner")
    merged["decimal_year"] = merged["year"] + (merged["month"] - 0.5) / 12.0
    merged = merged.sort_values("decimal_year").reset_index(drop=True)
    return merged


# ============================================================================
# MAIN: LOAD ALL SITES, AGGREGATE, SUMMARISE
# ============================================================================

def main():
    print("=" * 70)
    print("Phase 1: Data extraction and pairing for KIE_sites experiment")
    print("=" * 70)

    summary = {}

    for site in COLOCATED_SITES:
        code = site["code"]
        print(f"\n--- {code} (lat {site['lat']:+.1f}°) ---")

        # Load δ¹³C
        try:
            d13c = load_d13C_event(site["d13c_file"])
            print(f"  δ¹³C: {len(d13c)} obs, {d13c['year'].min()}–{d13c['year'].max()}")
        except FileNotFoundError as e:
            print(f"  δ¹³C: MISSING — {e}")
            continue

        # Load δD — try each lab file; use the first available from INSTAAR
        # (same-flask pairing).  Store all labs for reference.
        dD_all_labs = {}
        for dD_id in site["dD_files"]:
            try:
                dD_all_labs[dD_id] = load_dD_event(dD_id)
                n = len(dD_all_labs[dD_id])
                yr0 = dD_all_labs[dD_id]["year"].min()
                yr1 = dD_all_labs[dD_id]["year"].max()
                print(f"  δD [{dD_id}]: {n} obs, {yr0}–{yr1}")
            except FileNotFoundError:
                print(f"  δD [{dD_id}]: file not found, skipping")

        if not dD_all_labs:
            print(f"  ⚠ No δD data for {code}, skipping")
            continue

        # Use the INSTAAR δD file (same flasks as δ¹³C) as the primary pairing.
        # INSTAAR files have no lab suffix (e.g., "mlo" not "mloIMAU").
        primary_dD_id = site["d13c_file"]  # INSTAAR δD uses same site code
        if primary_dD_id in dD_all_labs:
            dD_primary = dD_all_labs[primary_dD_id]
            dD_lab = "INSTAAR"
        else:
            # Fallback: use whichever lab has the most observations
            primary_dD_id = max(dD_all_labs, key=lambda k: len(dD_all_labs[k]))
            dD_primary = dD_all_labs[primary_dD_id]
            dD_lab = primary_dD_id
            print(f"  ⚠ No INSTAAR δD for {code}, using {dD_lab} as primary")

        # Compute monthly means and pair
        monthly = compute_monthly_means(d13c, dD_primary)
        n_paired_months = len(monthly)
        print(f"  Paired months (both isotopes): {n_paired_months}")

        if n_paired_months == 0:
            print(f"  ⚠ No overlapping months for {code}, skipping")
            continue

        # Save monthly CSV
        out_csv = RESULTS_DIR / f"site_monthly_{code}.csv"
        monthly.to_csv(out_csv, index=False, float_format="%.4f")

        # Per-site summary
        overlap_start = monthly["decimal_year"].min()
        overlap_end = monthly["decimal_year"].max()
        n_years_overlap = overlap_end - overlap_start

        summary[code] = {
            "latitude": site["lat"],
            "mbl_site": site["mbl"],
            "dD_primary_lab": dD_lab,
            "d13C_n_obs": int(len(d13c)),
            "dD_n_obs": int(len(dD_primary)),
            "paired_months": int(n_paired_months),
            "overlap_start": round(float(overlap_start), 2),
            "overlap_end": round(float(overlap_end), 2),
            "overlap_years": round(float(n_years_overlap), 1),
            "d13C_mean": round(float(monthly["d13C_mean"].mean()), 3),
            "d13C_seasonal_range": round(float(
                monthly["d13C_mean"].max() - monthly["d13C_mean"].min()
            ), 3),
            "dD_mean": round(float(monthly["dD_mean"].mean()), 2),
            "dD_seasonal_range": round(float(
                monthly["dD_mean"].max() - monthly["dD_mean"].min()
            ), 2),
            "dD_other_labs": [k for k in dD_all_labs if k != primary_dD_id],
        }

    # Save summary JSON
    with open(RESULTS_DIR / "site_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\n✓ Summary saved to {RESULTS_DIR / 'site_summary.json'}")

    # ==========================================================================
    # FIGURE 1: Time series grid — δ¹³C and δD for all sites
    # ==========================================================================
    sites_with_data = [s for s in COLOCATED_SITES if s["code"] in summary]
    n_sites = len(sites_with_data)

    fig, axes = plt.subplots(n_sites, 2, figsize=(14, 2.2 * n_sites),
                             sharex=True, squeeze=False)
    fig.suptitle("Co-located δ¹³C-CH₄ and δD-CH₄ Monthly Means", fontsize=13, y=0.995)

    for i, site in enumerate(sites_with_data):
        code = site["code"]
        monthly = pd.read_csv(RESULTS_DIR / f"site_monthly_{code}.csv")

        # Left panel: δ¹³C
        ax_c = axes[i, 0]
        ax_c.plot(monthly["decimal_year"], monthly["d13C_mean"],
                  "o-", ms=2.5, lw=0.8, color="C0")
        ax_c.set_ylabel("δ¹³C (‰)", fontsize=8)
        ax_c.set_title(f"{code} ({site['lat']:+.1f}°)", fontsize=9, loc="left")
        ax_c.tick_params(labelsize=7)

        # Right panel: δD
        ax_d = axes[i, 1]
        ax_d.plot(monthly["decimal_year"], monthly["dD_mean"],
                  "o-", ms=2.5, lw=0.8, color="C1")
        ax_d.set_ylabel("δD (‰)", fontsize=8)
        ax_d.set_title(f"{code} ({site['lat']:+.1f}°)", fontsize=9, loc="left")
        ax_d.tick_params(labelsize=7)

    axes[-1, 0].set_xlabel("Year", fontsize=9)
    axes[-1, 1].set_xlabel("Year", fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.99])
    fig.savefig(FIGURES_DIR / "fig1_timeseries_grid.png", dpi=300)
    plt.close(fig)
    print(f"✓ Figure saved: {FIGURES_DIR / 'fig1_timeseries_grid.png'}")

    # ==========================================================================
    # FIGURE 2: Temporal coverage heatmap (months with paired data per site)
    # ==========================================================================
    fig2, ax2 = plt.subplots(figsize=(12, 4))

    # Build a matrix: rows = sites (sorted by latitude, N→S), cols = year-month
    all_year_months = set()
    site_months_dict = {}
    for site in sites_with_data:
        code = site["code"]
        monthly = pd.read_csv(RESULTS_DIR / f"site_monthly_{code}.csv")
        ym_set = set(zip(monthly["year"].astype(int), monthly["month"].astype(int)))
        site_months_dict[code] = ym_set
        all_year_months |= ym_set

    if all_year_months:
        ym_sorted = sorted(all_year_months)
        ym_to_idx = {ym: j for j, ym in enumerate(ym_sorted)}
        n_cols = len(ym_sorted)

        # Sort sites by latitude (N → S)
        sorted_sites = sorted(sites_with_data, key=lambda s: -s["lat"])
        labels = [f"{s['code']} ({s['lat']:+.0f}°)" for s in sorted_sites]
        matrix = np.zeros((len(sorted_sites), n_cols))

        for i, site in enumerate(sorted_sites):
            for ym in site_months_dict[site["code"]]:
                matrix[i, ym_to_idx[ym]] = 1.0

        ax2.imshow(matrix, aspect="auto", cmap="YlGn", interpolation="nearest",
                   vmin=0, vmax=1)
        ax2.set_yticks(range(len(labels)))
        ax2.set_yticklabels(labels, fontsize=8)

        # X-axis: show January of each year
        jan_indices = [ym_to_idx[(y, 1)] for y, m in ym_sorted
                       if m == 1 and (y, 1) in ym_to_idx]
        jan_labels = [str(y) for y, m in ym_sorted
                      if m == 1 and (y, 1) in ym_to_idx]
        # Deduplicate (ym_sorted may have duplicates from set iteration)
        seen = set()
        jan_idx_clean, jan_lbl_clean = [], []
        for idx, lbl in zip(jan_indices, jan_labels):
            if lbl not in seen:
                seen.add(lbl)
                jan_idx_clean.append(idx)
                jan_lbl_clean.append(lbl)
        ax2.set_xticks(jan_idx_clean)
        ax2.set_xticklabels(jan_lbl_clean, fontsize=7, rotation=45)
        ax2.set_xlabel("Year", fontsize=9)
        ax2.set_title("Paired δ¹³C + δD monthly data coverage (green = data present)",
                       fontsize=10)

    fig2.tight_layout()
    fig2.savefig(FIGURES_DIR / "fig1_data_coverage.png", dpi=300)
    plt.close(fig2)
    print(f"✓ Figure saved: {FIGURES_DIR / 'fig1_data_coverage.png'}")

    # ==========================================================================
    # Print summary table
    # ==========================================================================
    print("\n" + "=" * 70)
    print("SITE SUMMARY")
    print("=" * 70)
    print(f"{'Site':<6} {'Lat':>6} {'Lab':<9} {'Months':>6} {'Overlap':>12} "
          f"{'Δδ¹³C':>7} {'ΔδD':>7}")
    print("-" * 70)
    for code, info in sorted(summary.items(), key=lambda x: -x[1]["latitude"]):
        print(f"{code:<6} {info['latitude']:>+6.1f} {info['dD_primary_lab']:<9} "
              f"{info['paired_months']:>6} "
              f"{info['overlap_start']:.0f}–{info['overlap_end']:.0f}  "
              f"{info['d13C_seasonal_range']:>6.3f}‰ "
              f"{info['dD_seasonal_range']:>6.2f}‰")
    print(f"\nTotal sites with paired data: {len(summary)}")


if __name__ == "__main__":
    main()
