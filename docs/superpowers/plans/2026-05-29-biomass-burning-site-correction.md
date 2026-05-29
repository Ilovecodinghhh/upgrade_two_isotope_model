# Biomass Burning Site Correction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a standalone KIE_sites biomass-burning source-correction phase that uses GFED5 monthly CH4 emissions, compares against the existing wetland-only correction, and visualizes the incremental BB effect.

**Architecture:** Do not modify existing Phase 1-7 scripts. Create a new Phase 8 script that builds a seasonal zonal BB product from GFED5, applies BB phasor subtraction on top of existing Phase 6 wetland-corrected results, and writes new JSON/CSV/figure outputs under `results/phase8_biomass_burning/` and `figures/`.

**Tech Stack:** Python, numpy, xarray/netCDF4, matplotlib, json, pathlib, pytest.

---

## File Structure

- Create: `experiments/KIE_sites/analysis/phase8_biomass_burning.py`
  - Owns the full new phase: GFED5 loading, zonal seasonal harmonics, BB isotope signature lookup, BB phasor correction, comparison tables, and figures.
- Create: `experiments/KIE_sites/tests/test_phase8_biomass_burning.py`
  - Unit tests for GFED unit conversion, band aggregation, harmonic fitting, and additive BB phasor subtraction using small synthetic arrays.
- Create output at runtime: `experiments/KIE_sites/data/biomass_burning_seasonality.json`
  - GFED5 seasonal zonal BB harmonics and site assignments.
- Create output at runtime: `experiments/KIE_sites/results/phase8_biomass_burning/bb_correction_results.json`
  - Per-site observed, wetland-only, and wetland+BB comparison results.
- Create output at runtime: `experiments/KIE_sites/results/phase8_biomass_burning/bb_correction_summary.csv`
  - Compact table for manuscripts/review.
- Create figures at runtime:
  - `experiments/KIE_sites/figures/fig15_bb_seasonality_by_band.png`
  - `experiments/KIE_sites/figures/fig16_bb_source_phasor_comparison.png`
  - `experiments/KIE_sites/figures/fig17_bb_correction_comparison.png`

Do not edit `phase6_phasor.py`; Phase 8 may read `phase6_phasor.py` constants/helpers by import or duplicate tiny pure helpers if import proves brittle.

## Task 1: Synthetic Tests for BB Seasonal Utilities

**Files:**
- Create: `experiments/KIE_sites/tests/test_phase8_biomass_burning.py`
- Create later: `experiments/KIE_sites/analysis/phase8_biomass_burning.py`

- [ ] **Step 1: Add tests for unit conversion and band sums**

```python
import numpy as np

from experiments.KIE_sites.analysis.phase8_biomass_burning import (
    BANDS,
    aggregate_monthly_bands,
    gfed_ch4_to_tg_month,
)


def test_gfed_ch4_to_tg_month_converts_g_m2_with_area_to_tg():
    ch4_g_m2 = np.array([[[2.0, 3.0], [4.0, 5.0]]])
    grid_area_m2 = np.array([[10.0, 10.0], [20.0, 20.0]])

    result = gfed_ch4_to_tg_month(ch4_g_m2, grid_area_m2)

    assert result.shape == (1, 2, 2)
    assert np.isclose(result.sum(), (20 + 30 + 80 + 100) * 1e-12)


def test_aggregate_monthly_bands_sums_expected_latitude_bands():
    monthly_tg = np.ones((12, 4, 2))
    lat = np.array([75.0, 45.0, 0.0, -45.0])
    bands = {
        "NH_high": {"lat_min": 60.0, "lat_max": 90.0},
        "NH_mid": {"lat_min": 30.0, "lat_max": 60.0},
        "Tropics": {"lat_min": -30.0, "lat_max": 30.0},
        "SH_extra": {"lat_min": -90.0, "lat_max": -30.0},
    }

    result = aggregate_monthly_bands(monthly_tg, lat, bands)

    for name in bands:
        assert result[name].shape == (12,)
        assert np.allclose(result[name], 2.0)
```

- [ ] **Step 2: Add tests for harmonic fitting and BB source addition**

```python
import numpy as np

from experiments.KIE_sites.analysis.phase8_biomass_burning import (
    fit_annual_harmonic,
    source_phasor,
    apply_bb_to_wetland_corrected_decomposition,
)


def test_fit_annual_harmonic_recovers_sine_cosine_coefficients():
    months = np.arange(12)
    values = 10.0 + 2.0 * np.sin(2 * np.pi * months / 12) - 3.0 * np.cos(2 * np.pi * months / 12)

    result = fit_annual_harmonic(values)

    assert np.isclose(result["Q_mean_Tg_month"], 10.0, atol=1e-10)
    assert np.isclose(result["B_Q_Tg_month"], 2.0, atol=1e-10)
    assert np.isclose(result["C_Q_Tg_month"], -3.0, atol=1e-10)


def test_apply_bb_subtracts_bb_source_from_existing_wetland_sink():
    phase6_decomp = {
        "Z_sink_13c": [1.0, 2.0],
        "Z_sink_dD": [10.0, 20.0],
    }
    bb_src_13c = complex(-0.1, 0.2)
    bb_src_dD = complex(-1.0, 2.0)

    result = apply_bb_to_wetland_corrected_decomposition(phase6_decomp, bb_src_13c, bb_src_dD)

    assert np.allclose(result["Z_sink_wetland_plus_bb_13c"], [1.1, 1.8])
    assert np.allclose(result["Z_sink_wetland_plus_bb_dD"], [11.0, 18.0])
```

- [ ] **Step 3: Run tests to verify they fail before implementation**

Run:

```bash
pytest experiments/KIE_sites/tests/test_phase8_biomass_burning.py -v
```

Expected: FAIL because `phase8_biomass_burning.py` does not exist yet.

## Task 2: Implement GFED5 Seasonal Zonal BB Product

**Files:**
- Create: `experiments/KIE_sites/analysis/phase8_biomass_burning.py`
- Output: `experiments/KIE_sites/data/biomass_burning_seasonality.json`

- [ ] **Step 1: Add constants and paths**

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 8: biomass-burning source phasor correction for KIE_sites.

This phase is standalone. It does not modify Phase 1-7; it reads Phase 2
harmonics and Phase 6 wetland-only correction outputs, then writes new BB
comparison outputs.
"""

from pathlib import Path
import csv
import json

import numpy as np
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

EXPT_DIR = Path(__file__).resolve().parent.parent
MODEL_DIR = EXPT_DIR.parent.parent
GFED_DIR = (
    MODEL_DIR
    / "ImportantReferences"
    / "Riddell-Young2025PNAS_DS"
    / "Riddell-Young_2025_MassBalancePackage"
    / "Riddell-Young_2025_MassBalancePackage"
    / "data"
    / "GFED5_Beta"
)
PHASE2_JSON = EXPT_DIR / "results" / "phase2_harmonics" / "harmonic_fits.json"
PHASE6_JSON = EXPT_DIR / "results" / "phase6_phasor" / "phasor_results.json"
OUT_DATA_JSON = EXPT_DIR / "data" / "biomass_burning_seasonality.json"
OUT_DIR = EXPT_DIR / "results" / "phase8_biomass_burning"
FIG_DIR = EXPT_DIR / "figures"

OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

KG_TO_TG = 1e-9
G_TO_TG = 1e-12
CLIMATOLOGY_YEARS = tuple(range(2005, 2011))
CLEAN_SITES = ["ALT", "ZEP", "BRW", "CBA", "MHD", "KUM", "CGO", "SPO"]

BANDS = {
    "NH_high": {"lat_min": 60.0, "lat_max": 90.0, "label": "NH high (60-90N)"},
    "NH_mid": {"lat_min": 30.0, "lat_max": 60.0, "label": "NH mid (30-60N)"},
    "Tropics": {"lat_min": -30.0, "lat_max": 30.0, "label": "Tropics (30S-30N)"},
    "SH_extra": {"lat_min": -90.0, "lat_max": -30.0, "label": "SH extra (90-30S)"},
}

SITE_BAND = {
    "ALT": "NH_high",
    "ZEP": "NH_high",
    "BRW": "NH_high",
    "CBA": "NH_mid",
    "MHD": "NH_mid",
    "AZR": "NH_mid",
    "MLO": "Tropics",
    "KUM": "Tropics",
    "ASC": "Tropics",
    "SMO": "Tropics",
    "CGO": "SH_extra",
    "SPO": "SH_extra",
}

SIGNATURE_REGION_BY_BAND = {
    "NH_high": "NHext",
    "NH_mid": "NHext",
    "Tropics": "Trop",
    "SH_extra": "SHext",
}

D13C_ATM = -47.3
DD_ATM = -86.0
Q_TOTAL_TG_MONTH = 580.0 / 12.0
```

- [ ] **Step 2: Implement GFED conversion and band aggregation**

```python
def gfed_ch4_to_tg_month(ch4_g_m2_month, grid_area_m2):
    """Convert GFED CH4 from g CH4 m^-2 month^-1 to Tg CH4 month^-1 per grid cell."""
    return np.asarray(ch4_g_m2_month, dtype=np.float64) * np.asarray(grid_area_m2, dtype=np.float64) * G_TO_TG


def aggregate_monthly_bands(monthly_tg, lat, bands=BANDS):
    """Sum monthly grid-cell Tg emissions into latitude bands."""
    lat = np.asarray(lat, dtype=np.float64)
    monthly_tg = np.asarray(monthly_tg, dtype=np.float64)
    out = {}
    for name, cfg in bands.items():
        mask = (lat >= cfg["lat_min"]) & (lat < cfg["lat_max"])
        out[name] = np.nansum(monthly_tg[:, mask, :], axis=(1, 2))
    out["Global"] = np.nansum(monthly_tg, axis=(1, 2))
    return out
```

- [ ] **Step 3: Implement GFED yearly loader**

```python
def load_gfed_year(year):
    """Load one GFED5 annual file and return monthly Tg grid plus latitudes."""
    path = GFED_DIR / f"GFED5_Beta_monthly_{year}.nc"
    if not path.exists():
        raise FileNotFoundError(f"Missing GFED5 file: {path}")
    with xr.open_dataset(path) as ds:
        ch4 = ds["CH4"].values
        area = ds["grid_area"].values
        lat = ds["lat"].values
    monthly_tg = gfed_ch4_to_tg_month(ch4, area)
    return monthly_tg, lat
```

- [ ] **Step 4: Implement harmonic fitting and JSON builder**

```python
MONTH_NAMES = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]


def fit_annual_harmonic(clim_12):
    """Fit Q(m) = mean + B*sin(2*pi*m/12) + C*cos(2*pi*m/12)."""
    clim_12 = np.asarray(clim_12, dtype=np.float64)
    m = np.arange(12)
    omega = 2 * np.pi / 12.0
    design = np.column_stack([np.ones(12), np.sin(omega * m), np.cos(omega * m)])
    coeffs, _, _, _ = np.linalg.lstsq(design, clim_12, rcond=None)
    fitted = design @ coeffs
    mean, bq, cq = coeffs
    amp = float(np.hypot(bq, cq))
    phase_rad = float(np.arctan2(cq, bq))
    peak_month = float(((np.pi / 2 - phase_rad) / omega) % 12)
    ss_res = float(np.sum((clim_12 - fitted) ** 2))
    ss_tot = float(np.sum((clim_12 - mean) ** 2))
    return {
        "Q_mean_Tg_month": float(mean),
        "B_Q_Tg_month": float(bq),
        "C_Q_Tg_month": float(cq),
        "amplitude_Tg_month": amp,
        "phase_rad": phase_rad,
        "peak_month_index": peak_month,
        "peak_month_name": MONTH_NAMES[int(round(peak_month)) % 12],
        "frac_amplitude": float(amp / mean) if mean > 0 else float("nan"),
        "R2": float(1 - ss_res / ss_tot) if ss_tot > 0 else float("nan"),
        "annual_Tg_yr": float(mean * 12.0),
        "monthly_clim_Tg": [float(x) for x in clim_12],
    }
```

- [ ] **Step 5: Build and save `biomass_burning_seasonality.json`**

```python
def build_bb_seasonality():
    """Build GFED5 BB monthly climatology and harmonic coefficients by latitude band."""
    monthly_by_band = {name: [] for name in [*BANDS.keys(), "Global"]}
    for year in CLIMATOLOGY_YEARS:
        monthly_tg, lat = load_gfed_year(year)
        band_totals = aggregate_monthly_bands(monthly_tg, lat, BANDS)
        for name, values in band_totals.items():
            monthly_by_band[name].append(values)

    bands = {}
    for name, year_arrays in monthly_by_band.items():
        climatology = np.vstack(year_arrays).mean(axis=0)
        bands[name] = fit_annual_harmonic(climatology)
        if name in BANDS:
            bands[name]["label"] = BANDS[name]["label"]
        else:
            bands[name]["label"] = "Global"

    site_assignment = {}
    for site, band in SITE_BAND.items():
        site_assignment[site] = {"band": band, **bands[band]}

    payload = {
        "metadata": {
            "source": "GFED5 Beta monthly CH4 emissions",
            "gfed_dir": str(GFED_DIR),
            "variable": "CH4 (g CH4 m^-2 month^-1) multiplied by grid_area (m2)",
            "climatology_years": list(CLIMATOLOGY_YEARS),
            "harmonic": "Q(m) = Q_mean + B_Q*sin(2*pi*m/12) + C_Q*cos(2*pi*m/12), m=0..11",
            "note": "Standalone Phase 8 input; existing wetland phase files are not modified.",
        },
        "bands": bands,
        "site_assignment": site_assignment,
    }
    OUT_DATA_JSON.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_DATA_JSON, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, allow_nan=False)
    return payload
```

- [ ] **Step 6: Run unit tests**

Run:

```bash
pytest experiments/KIE_sites/tests/test_phase8_biomass_burning.py -v
```

Expected: PASS for conversion, aggregation, and harmonic tests.

## Task 3: Implement BB Source Signatures and Phasor Correction

**Files:**
- Modify: `experiments/KIE_sites/analysis/phase8_biomass_burning.py`
- Read: `rel/data/BB_d13C_{NHext,Trop,SHext}_MC.csv`
- Read: `rel/data/BB_dD_{NHext,Trop,SHext}_MC.csv`
- Read: `experiments/KIE_sites/results/phase6_phasor/phasor_results.json`

- [ ] **Step 1: Add CSV signature loader**

```python
def load_bb_signature(region, isotope, years=CLIMATOLOGY_YEARS):
    """Return mean and std for BB isotope signature over selected years and MC draws."""
    if isotope not in {"d13C", "dD"}:
        raise ValueError(f"Unsupported isotope: {isotope}")
    path = MODEL_DIR / "rel" / "data" / f"BB_{isotope}_{region}_MC.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing BB signature file: {path}")

    values = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        for row in reader:
            year = int(float(row[0]))
            if year in years:
                values.extend(float(x) for x in row[1:] if x != "")
    arr = np.asarray(values, dtype=np.float64)
    return {"mean": float(np.nanmean(arr)), "sigma": float(np.nanstd(arr)), "file": str(path)}


def load_bb_signatures_by_band():
    """Map each BB seasonal band to matching existing three-box BB isotope signatures."""
    out = {}
    for band, region in SIGNATURE_REGION_BY_BAND.items():
        out[band] = {
            "region": region,
            "d13C": load_bb_signature(region, "d13C"),
            "dD": load_bb_signature(region, "dD"),
        }
    return out
```

- [ ] **Step 2: Add generic source phasor helpers**

```python
def rotate_month_index_to_midpoint(B_Q, C_Q):
    """Rotate month-index harmonic coefficients into Phase 2 month-midpoint basis."""
    delta = 2 * np.pi * 0.5 / 12.0
    cos_d = np.cos(delta)
    sin_d = np.sin(delta)
    return B_Q * cos_d + C_Q * sin_d, -B_Q * sin_d + C_Q * cos_d


def source_phasor(B_Q, C_Q, q_total, delta_source, delta_atm):
    """Compute isotope source phasor in permil from source seasonality and isotope gap."""
    B_mid, C_mid = rotate_month_index_to_midpoint(B_Q, C_Q)
    z_frac = complex(B_mid, C_mid) / q_total
    return (delta_source - delta_atm) * z_frac


def apply_bb_to_wetland_corrected_decomposition(phase6_decomp, bb_src_13c, bb_src_dD):
    """Subtract BB source phasor from existing Phase 6 wetland-corrected sink phasor."""
    z_wet_sink_13c = complex(*phase6_decomp["Z_sink_13c"])
    z_wet_sink_dD = complex(*phase6_decomp["Z_sink_dD"])
    z_wet_bb_sink_13c = z_wet_sink_13c - bb_src_13c
    z_wet_bb_sink_dD = z_wet_sink_dD - bb_src_dD
    r = abs(z_wet_bb_sink_13c) / abs(z_wet_bb_sink_dD) if abs(z_wet_bb_sink_dD) > 0 else np.nan
    return {
        "Z_src_bb_13c": [float(bb_src_13c.real), float(bb_src_13c.imag)],
        "Z_src_bb_dD": [float(bb_src_dD.real), float(bb_src_dD.imag)],
        "A_src_bb_13c": float(abs(bb_src_13c)),
        "A_src_bb_dD": float(abs(bb_src_dD)),
        "Z_sink_wetland_plus_bb_13c": [float(z_wet_bb_sink_13c.real), float(z_wet_bb_sink_13c.imag)],
        "Z_sink_wetland_plus_bb_dD": [float(z_wet_bb_sink_dD.real), float(z_wet_bb_sink_dD.imag)],
        "A_sink_wetland_plus_bb_13c": float(abs(z_wet_bb_sink_13c)),
        "A_sink_wetland_plus_bb_dD": float(abs(z_wet_bb_sink_dD)),
        "R_wetland_plus_bb": float(r),
    }
```

- [ ] **Step 3: Implement per-site correction loop**

```python
def load_json(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def run_bb_correction(bb_seasonality):
    """Apply BB correction on top of Phase 6 wetland-only results."""
    phase6 = load_json(PHASE6_JSON)
    signatures = load_bb_signatures_by_band()
    results = {
        "metadata": {
            "method": "Standalone Phase 8 BB source phasor correction applied after Phase 6 wetland correction",
            "phase6_input": str(PHASE6_JSON),
            "bb_seasonality_input": str(OUT_DATA_JSON),
            "existing_phase_files_modified": False,
        },
        "bb_signatures": signatures,
        "sites": {},
    }

    for code in CLEAN_SITES:
        if code not in phase6["sites"]:
            continue
        site_bb = bb_seasonality["site_assignment"][code]
        band = site_bb["band"]
        sig = signatures[band]
        bq = site_bb["B_Q_Tg_month"]
        cq = site_bb["C_Q_Tg_month"]
        bb_src_13c = source_phasor(bq, cq, Q_TOTAL_TG_MONTH, sig["d13C"]["mean"], D13C_ATM)
        bb_src_dD = source_phasor(bq, cq, Q_TOTAL_TG_MONTH, sig["dD"]["mean"], DD_ATM)
        phase6_decomp = phase6["sites"][code]["deterministic"]
        corrected = apply_bb_to_wetland_corrected_decomposition(phase6_decomp, bb_src_13c, bb_src_dD)
        results["sites"][code] = {
            "band": band,
            "bb_signature_region": sig["region"],
            "R_obs": phase6_decomp["R_obs"],
            "R_wetland_only": phase6_decomp["R_corrected"],
            "A_src_wetland_13c": phase6_decomp["A_src_13c"],
            "A_src_wetland_dD": phase6_decomp["A_src_dD"],
            **corrected,
            "delta_R_vs_wetland_only": corrected["R_wetland_plus_bb"] - phase6_decomp["R_corrected"],
        }
    return results
```

- [ ] **Step 4: Save JSON and CSV summary**

```python
def save_results(results):
    out_json = OUT_DIR / "bb_correction_results.json"
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, allow_nan=False)

    out_csv = OUT_DIR / "bb_correction_summary.csv"
    fields = [
        "site", "band", "bb_signature_region", "R_obs", "R_wetland_only",
        "R_wetland_plus_bb", "delta_R_vs_wetland_only",
        "A_src_wetland_13c", "A_src_bb_13c", "A_src_wetland_dD", "A_src_bb_dD",
    ]
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for site, row in results["sites"].items():
            writer.writerow({key: site if key == "site" else row.get(key, "") for key in fields})
    return out_json, out_csv
```

- [ ] **Step 5: Re-run tests**

Run:

```bash
pytest experiments/KIE_sites/tests/test_phase8_biomass_burning.py -v
```

Expected: PASS.

## Task 4: Add Comparison Visualizations

**Files:**
- Modify: `experiments/KIE_sites/analysis/phase8_biomass_burning.py`
- Output figures under `experiments/KIE_sites/figures/`

- [ ] **Step 1: Plot BB seasonality by band**

```python
def plot_bb_seasonality(bb_seasonality):
    fig, ax = plt.subplots(figsize=(8, 4.8))
    months = np.arange(1, 13)
    for name, band in bb_seasonality["bands"].items():
        if name == "Global":
            continue
        ax.plot(months, band["monthly_clim_Tg"], marker="o", label=band["label"])
    ax.set_xlabel("Month")
    ax.set_ylabel("GFED5 BB CH4 (Tg/month)")
    ax.set_title("Fig 15: Biomass burning seasonal emissions by latitude band")
    ax.set_xticks(months)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    out = FIG_DIR / "fig15_bb_seasonality_by_band.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    return out
```

- [ ] **Step 2: Plot wetland vs BB source phasor amplitudes**

```python
def plot_source_phasor_comparison(results):
    sites = list(results["sites"].keys())
    x = np.arange(len(sites))
    wet13 = [results["sites"][s]["A_src_wetland_13c"] for s in sites]
    bb13 = [results["sites"][s]["A_src_bb_13c"] for s in sites]
    wetd = [results["sites"][s]["A_src_wetland_dD"] for s in sites]
    bbd = [results["sites"][s]["A_src_bb_dD"] for s in sites]

    fig, axes = plt.subplots(2, 1, figsize=(8.5, 6.2), sharex=True)
    width = 0.38
    axes[0].bar(x - width / 2, wet13, width, label="Wetland")
    axes[0].bar(x + width / 2, bb13, width, label="BB")
    axes[0].set_ylabel("A source d13C (per mil)")
    axes[0].legend(frameon=False)
    axes[1].bar(x - width / 2, wetd, width, label="Wetland")
    axes[1].bar(x + width / 2, bbd, width, label="BB")
    axes[1].set_ylabel("A source dD (per mil)")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(sites)
    fig.suptitle("Fig 16: Wetland and biomass-burning source phasor amplitudes")
    fig.tight_layout()
    out = FIG_DIR / "fig16_bb_source_phasor_comparison.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    return out
```

- [ ] **Step 3: Plot R comparison**

```python
def plot_correction_comparison(results):
    sites = list(results["sites"].keys())
    x = np.arange(len(sites))
    r_obs = [results["sites"][s]["R_obs"] for s in sites]
    r_wet = [results["sites"][s]["R_wetland_only"] for s in sites]
    r_bb = [results["sites"][s]["R_wetland_plus_bb"] for s in sites]

    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    ax.plot(x, r_obs, "o-", label="Observed")
    ax.plot(x, r_wet, "s-", label="Wetland corrected")
    ax.plot(x, r_bb, "^-", label="Wetland + BB corrected")
    ax.set_xticks(x)
    ax.set_xticklabels(sites)
    ax.set_ylabel("R = A(d13C) / A(dD)")
    ax.set_title("Fig 17: Incremental effect of biomass-burning source correction")
    ax.legend(frameon=False)
    fig.tight_layout()
    out = FIG_DIR / "fig17_bb_correction_comparison.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    return out
```

- [ ] **Step 4: Call plots from `main()`**

```python
def main():
    print("=" * 70)
    print("Phase 8: Biomass-burning source correction")
    print("=" * 70)
    bb_seasonality = build_bb_seasonality()
    results = run_bb_correction(bb_seasonality)
    out_json, out_csv = save_results(results)
    fig15 = plot_bb_seasonality(bb_seasonality)
    fig16 = plot_source_phasor_comparison(results)
    fig17 = plot_correction_comparison(results)
    print(f"Saved: {OUT_DATA_JSON}")
    print(f"Saved: {out_json}")
    print(f"Saved: {out_csv}")
    print(f"Saved: {fig15}")
    print(f"Saved: {fig16}")
    print(f"Saved: {fig17}")


if __name__ == "__main__":
    main()
```

## Task 5: Run Phase 8 and Verify Outputs

**Files:**
- Run: `experiments/KIE_sites/analysis/phase8_biomass_burning.py`
- Verify outputs in `experiments/KIE_sites/data/`, `experiments/KIE_sites/results/phase8_biomass_burning/`, and `experiments/KIE_sites/figures/`

- [ ] **Step 1: Run tests**

Run:

```bash
pytest experiments/KIE_sites/tests/test_phase8_biomass_burning.py -v
```

Expected: all tests pass.

- [ ] **Step 2: Run the new standalone phase**

Run:

```bash
python experiments/KIE_sites/analysis/phase8_biomass_burning.py
```

Expected: creates `biomass_burning_seasonality.json`, `bb_correction_results.json`, `bb_correction_summary.csv`, and three figure files.

- [ ] **Step 3: Inspect result summary**

Run:

```bash
python -c "import json; p='experiments/KIE_sites/results/phase8_biomass_burning/bb_correction_results.json'; d=json.load(open(p, encoding='utf-8')); print(d['metadata']['existing_phase_files_modified']); print(sorted(d['sites'])); print({k: round(v['delta_R_vs_wetland_only'], 5) for k,v in d['sites'].items()})"
```

Expected:

```text
False
['ALT', 'BRW', 'CBA', 'CGO', 'KUM', 'MHD', 'SPO', 'ZEP']
```

The final dictionary values should be finite numbers.

- [ ] **Step 4: Verify no existing phase scripts changed**

Run:

```bash
git diff -- experiments/KIE_sites/analysis/phase1_data.py experiments/KIE_sites/analysis/phase2_harmonics.py experiments/KIE_sites/analysis/phase3_synthesis.py experiments/KIE_sites/analysis/phase4_deconv.py experiments/KIE_sites/analysis/phase5_kie.py experiments/KIE_sites/analysis/phase6_phasor.py experiments/KIE_sites/analysis/phase7_yearly_stability.py
```

Expected: no diff output.

- [ ] **Step 5: Commit**

```bash
git add experiments/KIE_sites/analysis/phase8_biomass_burning.py \
        experiments/KIE_sites/tests/test_phase8_biomass_burning.py \
        experiments/KIE_sites/data/biomass_burning_seasonality.json \
        experiments/KIE_sites/results/phase8_biomass_burning \
        experiments/KIE_sites/figures/fig15_bb_seasonality_by_band.png \
        experiments/KIE_sites/figures/fig16_bb_source_phasor_comparison.png \
        experiments/KIE_sites/figures/fig17_bb_correction_comparison.png
git commit -m "feat: add standalone biomass burning site correction phase"
```

## Self-Review

- Spec coverage: standalone new phase, no edits to existing phases, GFED5 raw data use, BB correction, and comparison visualization are covered.
- Placeholder scan: no TBD/TODO placeholders remain.
- Type consistency: function names used in tests match implementation tasks.
- Scope: focused on Phase 8 only; existing Phase 1-7 remain read-only inputs.
