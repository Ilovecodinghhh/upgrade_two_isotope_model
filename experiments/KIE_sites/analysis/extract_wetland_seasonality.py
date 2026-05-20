#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_wetland_seasonality.py — Wetland CH₄ emission seasonality per site
==========================================================================

Reads the Li et al. (2026) ESSD monthly 1°×1° gridded wetland CH₄ emission
dataset and extracts the seasonal cycle for each latitude band relevant to
the 12 KIE_sites monitoring stations.

Workflow
--------
1.  Open Li2026ESSD_DS.nc (monthly, 1°×1°, kg CH₄ cell⁻¹ month⁻¹).
2.  For each latitude band, sum over longitude → zonal total (kg/month).
3.  Subset to 2005–2010 (the INSTAAR co-located δ¹³C + δD period).
4.  Compute monthly climatology (12-month mean over 6 years).
5.  Fit annual harmonic:  Q(t) = Q̄ + B·sin(2πm/12) + C·cos(2πm/12),
    where m = 0 … 11 (Jan … Dec).
6.  Assign each KIE site to its source-region band.
7.  Save structured JSON + summary figure.

Latitude bands
--------------
  NH_high   : 60–90 °N   (Arctic/boreal wetlands)
  NH_mid    : 30–60 °N   (temperate/boreal)
  Tropics   : 30 °S–30 °N
  SH_extra  : 90–30 °S   (austral extratropical)
  Global    : all latitudes (flux-weighted; used for CGO, SPO)

Harmonic convention
-------------------
  Q(m) = Q_mean + B_Q·sin(2π·m/12) + C_Q·cos(2π·m/12)

  where m = month index (0 = Jan, 11 = Dec).
  This matches Phase 2 convention (B·sin + C·cos) but with m in months
  rather than fractional years.  The phase6_phasor.py script converts
  between the two.

  Peak month = atan2(−C_Q, B_Q) × 12/(2π)    (0-indexed from Jan)
  Amplitude  = √(B_Q² + C_Q²)

Output
------
  data/wetland_seasonality.json
  figures/fig7_wetland_seasonality.png

References
----------
  Li et al. (2026) ESSD — doi:10.5281/zenodo.18870108
"""

from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================================================================
# PATHS
# ============================================================================
EXPT_DIR  = Path(__file__).resolve().parent.parent          # experiments/KIE_sites/
ROOT_DIR  = EXPT_DIR.parent.parent                          # repo root
NC_PATH   = ROOT_DIR / "ImportantReferences" / "Li2026ESSD_DS.nc"
OUT_JSON  = EXPT_DIR / "data"    / "wetland_seasonality.json"
OUT_FIG   = EXPT_DIR / "figures" / "fig7_wetland_seasonality.png"

(EXPT_DIR / "data").mkdir(parents=True, exist_ok=True)
(EXPT_DIR / "figures").mkdir(parents=True, exist_ok=True)

# ============================================================================
# CONSTANTS
# ============================================================================
KG_TO_TG = 1e-9   # 1 Tg = 1e9 kg

# Latitude bands  (lat_min inclusive, lat_max exclusive for slicing)
BANDS = {
    "NH_high":  {"lat_min": 60.0, "lat_max": 90.0,  "label": "NH high (60–90°N)"},
    "NH_mid":   {"lat_min": 30.0, "lat_max": 60.0,  "label": "NH mid (30–60°N)"},
    "Tropics":  {"lat_min":-30.0, "lat_max": 30.0,  "label": "Tropics (30°S–30°N)"},
    "SH_extra": {"lat_min":-90.0, "lat_max":-30.0,  "label": "SH extra (90–30°S)"},
}

# Site → source-region assignment
# Background marine/Arctic stations sample upstream emissions, not local.
SITE_BAND = {
    "ALT": "NH_high",   # 82.5°N  Arctic
    "ZEP": "NH_high",   # 78.9°N  Svalbard
    "BRW": "NH_high",   # 71.3°N  Barrow
    "CBA": "NH_mid",    # 55.2°N  Cold Bay
    "MHD": "NH_mid",    # 53.3°N  Mace Head
    "AZR": "NH_mid",    # 38.8°N  Azores (excluded)
    "MLO": "Tropics",   # 19.5°N  Mauna Loa (excluded)
    "KUM": "Tropics",   # 19.6°N  Kumukahi
    "ASC": "Tropics",   # −8.0°S  Ascension (excluded)
    "SMO": "Tropics",   # −14.2°S Samoa (excluded)
    "CGO": "SH_extra",  # −40.7°S Cape Grim — local SH wetlands; NH signal attenuated by transport
    "SPO": "SH_extra",  # −90.0°S South Pole — local SH wetlands; NH signal maximally attenuated
}

MONTH_NAMES = ["Jan","Feb","Mar","Apr","May","Jun",
               "Jul","Aug","Sep","Oct","Nov","Dec"]


# ============================================================================
# LOAD NETCDF
# ============================================================================
def load_li2026():
    """Load Li2026 ESSD NetCDF and return lat, lon, time arrays + wetch4.

    Returns
    -------
    lat   : 1-D array (180,)  in degrees N (−89.5 … 89.5)
    lon   : 1-D array (360,)  in degrees E (−179.5 … 179.5)
    time  : 1-D array (312,)  months-since-1999-12-01  → Jan 2000 = 1
    data  : 3-D array (312, 180, 360)  kg CH₄ cell⁻¹ month⁻¹, NaN-filled
    """
    import netCDF4 as nc
    ds = nc.Dataset(str(NC_PATH), "r")
    lat  = ds.variables["lat"][:]
    lon  = ds.variables["lon"][:]
    time = ds.variables["time"][:]
    data = ds.variables["wetch4"][:]              # masked array
    data = np.where(data.mask, np.nan, data.data) if hasattr(data, "mask") else np.array(data)
    ds.close()
    return np.asarray(lat), np.asarray(lon), np.asarray(time, dtype=int), data


# ============================================================================
# ZONAL SUMMATION
# ============================================================================
def zonal_sum(data, lat, lat_min, lat_max):
    """Sum emissions over longitude within a latitude band.

    Parameters
    ----------
    data : (ntime, nlat, nlon)  kg/month per cell, may contain NaN
    lat  : (nlat,)
    lat_min, lat_max : band boundaries (degrees N)

    Returns
    -------
    timeseries : (ntime,) in Tg/month
    """
    mask = (lat >= lat_min) & (lat < lat_max)
    band = data[:, mask, :]                       # (ntime, nlat_band, nlon)
    total_kg = np.nansum(band, axis=(1, 2))       # (ntime,)
    return total_kg * KG_TO_TG


def global_sum(data):
    """Sum emissions over all lat/lon → global total (Tg/month)."""
    return np.nansum(data, axis=(1, 2)) * KG_TO_TG


# ============================================================================
# SUBSET 2005–2010 & MONTHLY CLIMATOLOGY
# ============================================================================
def subset_2005_2010(timeseries, time_months):
    """Extract Jan 2005 – Dec 2010 (72 months).

    time_months is months-since-1999-12-01.  Jan 2000 = 1.
    Jan 2005 = 61, Dec 2010 = 132.
    """
    idx = (time_months >= 61) & (time_months <= 132)
    return timeseries[idx]


def monthly_climatology(ts_72):
    """Reshape 72 months (6 full years) into (6, 12), return mean (12,)."""
    assert len(ts_72) == 72, f"Expected 72 months, got {len(ts_72)}"
    return ts_72.reshape(6, 12).mean(axis=0)


# ============================================================================
# HARMONIC FIT
# ============================================================================
def fit_annual_harmonic(clim_12):
    """Fit  Q(m) = Q̄ + B·sin(2πm/12) + C·cos(2πm/12)  to 12-month climatology.

    Parameters
    ----------
    clim_12 : (12,) monthly means (Tg/month)

    Returns
    -------
    dict with Q_mean, B_Q, C_Q, amplitude, phase_month, frac_amplitude
    """
    m = np.arange(12)
    omega = 2 * np.pi / 12

    # Design matrix  [1, sin, cos]
    A = np.column_stack([np.ones(12), np.sin(omega * m), np.cos(omega * m)])
    coeffs, _, _, _ = np.linalg.lstsq(A, clim_12, rcond=None)

    Q_mean = coeffs[0]
    B_Q    = coeffs[1]
    C_Q    = coeffs[2]
    amp    = np.sqrt(B_Q**2 + C_Q**2)
    # Peak month: Q is maximised when sin(ωm+φ) = 1
    # Q = Q̄ + amp·sin(ωm + φ)  where  B = amp·cos(φ), C = amp·sin(φ)
    # → peak at ωm = π/2 − φ → m_peak = (π/2 − atan2(C,B)) / ω
    phase_rad  = np.arctan2(C_Q, B_Q)
    peak_month = ((np.pi / 2 - phase_rad) / omega) % 12

    frac_amp = amp / Q_mean if Q_mean > 0 else np.nan

    # Residual variance check
    fitted = A @ coeffs
    ss_res = np.sum((clim_12 - fitted)**2)
    ss_tot = np.sum((clim_12 - Q_mean)**2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan

    return {
        "Q_mean_Tg_month": round(float(Q_mean), 4),
        "B_Q_Tg_month":    round(float(B_Q), 4),
        "C_Q_Tg_month":    round(float(C_Q), 4),
        "amplitude_Tg_month": round(float(amp), 4),
        "phase_rad":        round(float(phase_rad), 4),
        "peak_month_index": round(float(peak_month), 2),   # 0=Jan … 6=Jul
        "peak_month_name":  MONTH_NAMES[int(round(peak_month)) % 12],
        "frac_amplitude":   round(float(frac_amp), 4),
        "R2":               round(float(R2), 4),
        "annual_Tg_yr":     round(float(Q_mean * 12), 1),
        "monthly_clim_Tg":  [round(float(v), 4) for v in clim_12],
    }


# ============================================================================
# MAIN
# ============================================================================
def main():
    print("=" * 70)
    print("Extracting wetland CH₄ emission seasonality from Li et al. (2026)")
    print("=" * 70)

    # ── Load data ──
    lat, lon, time, data = load_li2026()
    print(f"Grid: {data.shape}  (time={len(time)}, lat={len(lat)}, lon={len(lon)})")

    # ── Process each latitude band ──
    results = {"metadata": {}, "bands": {}, "site_assignment": {}}
    results["metadata"] = {
        "source": "Li et al. (2026) ESSD, doi:10.5281/zenodo.18870108",
        "variable": "wetch4 (kg CH4 / cell / month)",
        "period": "2005-2010 monthly climatology",
        "harmonic": "Q(m) = Q_mean + B_Q*sin(2*pi*m/12) + C_Q*cos(2*pi*m/12), m=0..11 (Jan..Dec)",
        "date_built": "2026-05-20",
    }

    for band_key, band_info in BANDS.items():
        ts = zonal_sum(data, lat, band_info["lat_min"], band_info["lat_max"])
        ts_sub = subset_2005_2010(ts, time)
        clim = monthly_climatology(ts_sub)
        fit = fit_annual_harmonic(clim)
        fit["label"] = band_info["label"]
        results["bands"][band_key] = fit

        print(f"\n{band_info['label']}:")
        print(f"  Annual: {fit['annual_Tg_yr']:.1f} Tg/yr")
        print(f"  Peak: {fit['peak_month_name']} (index {fit['peak_month_index']:.1f})")
        print(f"  Frac amplitude: {fit['frac_amplitude']:.3f}")
        print(f"  R²: {fit['R2']:.3f}")
        clim_str = "  Monthly: " + " ".join(f"{v:.2f}" for v in clim)
        print(clim_str)

    # ── Global (for CGO/SPO) ──
    ts_global = global_sum(data)
    ts_global_sub = subset_2005_2010(ts_global, time)
    clim_global = monthly_climatology(ts_global_sub)
    fit_global = fit_annual_harmonic(clim_global)
    fit_global["label"] = "Global (all latitudes)"
    results["bands"]["Global"] = fit_global

    print(f"\nGlobal:")
    print(f"  Annual: {fit_global['annual_Tg_yr']:.1f} Tg/yr")
    print(f"  Peak: {fit_global['peak_month_name']} (index {fit_global['peak_month_index']:.1f})")
    print(f"  Frac amplitude: {fit_global['frac_amplitude']:.3f}")
    print(f"  R²: {fit_global['R2']:.3f}")

    # ── Assign sites ──
    for site, band in SITE_BAND.items():
        if band in results["bands"]:
            results["site_assignment"][site] = {
                "source_band": band,
                "label": results["bands"][band]["label"],
                "B_Q_Tg_month": results["bands"][band]["B_Q_Tg_month"],
                "C_Q_Tg_month": results["bands"][band]["C_Q_Tg_month"],
                "Q_mean_Tg_month": results["bands"][band]["Q_mean_Tg_month"],
                "amplitude_Tg_month": results["bands"][band]["amplitude_Tg_month"],
                "peak_month_name": results["bands"][band]["peak_month_name"],
                "frac_amplitude": results["bands"][band]["frac_amplitude"],
            }

    # ── Save JSON ──
    with open(OUT_JSON, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n✓ Saved: {OUT_JSON}")

    # ── Figure ──
    plot_seasonality(results)

    # ── Summary table ──
    print("\n" + "=" * 70)
    print("SITE ASSIGNMENTS")
    print("=" * 70)
    print(f"{'Site':<5} {'Band':<12} {'Q_mean':>8} {'B_Q':>8} {'C_Q':>8} "
          f"{'Amp':>8} {'Peak':>5} {'fAmp':>6}")
    print("-" * 70)
    for site in ["ALT","ZEP","BRW","CBA","MHD","AZR","MLO","KUM","ASC","SMO","CGO","SPO"]:
        s = results["site_assignment"].get(site, {})
        if s:
            print(f"{site:<5} {s['source_band']:<12} "
                  f"{s['Q_mean_Tg_month']:>8.3f} {s['B_Q_Tg_month']:>8.4f} "
                  f"{s['C_Q_Tg_month']:>8.4f} {s['amplitude_Tg_month']:>8.4f} "
                  f"{s['peak_month_name']:>5} {s['frac_amplitude']:>6.3f}")


# ============================================================================
# FIGURE
# ============================================================================
def plot_seasonality(results):
    """Fig 7: Monthly wetland emission climatology + harmonic fit per band."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharey=False)
    axes = axes.flatten()

    band_keys = ["NH_high", "NH_mid", "Tropics", "SH_extra", "Global"]
    colours   = ["#d62728", "#ff7f0e", "#2ca02c", "#1f77b4", "#9467bd"]
    m = np.arange(12)
    m_fine = np.linspace(0, 11, 200)
    omega = 2 * np.pi / 12

    for i, (bk, col) in enumerate(zip(band_keys, colours)):
        ax = axes[i]
        bd = results["bands"][bk]
        clim = np.array(bd["monthly_clim_Tg"])
        fitted = (bd["Q_mean_Tg_month"]
                  + bd["B_Q_Tg_month"] * np.sin(omega * m_fine)
                  + bd["C_Q_Tg_month"] * np.cos(omega * m_fine))

        ax.bar(m, clim, color=col, alpha=0.5, width=0.7, label="Li2026 climatology")
        ax.plot(m_fine, fitted, "-", color=col, lw=2, label="Annual harmonic fit")
        ax.axhline(bd["Q_mean_Tg_month"], color="gray", ls=":", lw=1, label=f"Mean = {bd['Q_mean_Tg_month']:.2f}")

        ax.set_title(f"{bd['label']}\n{bd['annual_Tg_yr']:.0f} Tg/yr, peak {bd['peak_month_name']}, "
                     f"frac amp = {bd['frac_amplitude']:.2f}, R² = {bd['R2']:.3f}",
                     fontsize=9)
        ax.set_xticks(m)
        ax.set_xticklabels([mn[0] for mn in MONTH_NAMES], fontsize=8)
        ax.set_ylabel("Tg CH₄ / month", fontsize=9)
        ax.legend(fontsize=7, loc="upper left")

    # Panel (f): all bands normalised
    ax = axes[5]
    for bk, col in zip(band_keys, colours):
        bd = results["bands"][bk]
        clim = np.array(bd["monthly_clim_Tg"])
        if bd["Q_mean_Tg_month"] > 0:
            ax.plot(m, clim / bd["Q_mean_Tg_month"], "o-", color=col, ms=4,
                    label=f"{bd['label']} (fA={bd['frac_amplitude']:.2f})")
    ax.axhline(1.0, color="gray", ls=":", lw=1)
    ax.set_title("(f) Normalised seasonal cycle (all bands)", fontsize=9)
    ax.set_xticks(m)
    ax.set_xticklabels([mn[0] for mn in MONTH_NAMES], fontsize=8)
    ax.set_ylabel("Q / Q̄", fontsize=9)
    ax.legend(fontsize=7)

    fig.suptitle("Li et al. (2026) Wetland CH₄ Emission Seasonality — 2005–2010 Climatology",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT_FIG, dpi=300)
    plt.close(fig)
    print(f"✓ Figure saved: {OUT_FIG}")


if __name__ == "__main__":
    main()
