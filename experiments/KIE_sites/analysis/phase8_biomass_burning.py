#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 8: biomass-burning source phasor correction for KIE_sites.

This phase is standalone. It does not modify Phase 1-7; it reads existing
Phase 6 wetland-only correction outputs, then writes new BB comparison outputs.
"""

from pathlib import Path
import csv
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr


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
PHASE6_JSON = EXPT_DIR / "results" / "phase6_phasor" / "phasor_results.json"
OUT_DATA_JSON = EXPT_DIR / "data" / "biomass_burning_seasonality.json"
OUT_DIR = EXPT_DIR / "results" / "phase8_biomass_burning"
FIG_DIR = EXPT_DIR / "figures"

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

ALPHA_13C_SAUERESSIG = 1.0039
ALPHA_13C_CANTRELL = 1.0054
ALPHA_D_OH = 1.294
ALPHA_13C_CL = 1.066
ALPHA_D_CL = 1.508
ALPHA_13C_SOIL = 1.022
ALPHA_D_SOIL = 1.066
ALPHA_13C_STRAT = 1.013
ALPHA_D_STRAT = 1.16
F_OH = 0.84
F_CL = 0.035
F_SOIL = 0.06
F_STRAT = 0.065

MONTH_NAMES = [
    "Jan",
    "Feb",
    "Mar",
    "Apr",
    "May",
    "Jun",
    "Jul",
    "Aug",
    "Sep",
    "Oct",
    "Nov",
    "Dec",
]


def json_safe(value):
    """Convert numpy values recursively and replace non-finite floats with None."""
    if isinstance(value, dict):
        return {k: json_safe(v) for k, v in value.items()}
    if isinstance(value, list):
        return [json_safe(v) for v in value]
    if isinstance(value, tuple):
        return [json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        value = float(value)
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def write_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(json_safe(payload), f, indent=2, allow_nan=False)


def load_json(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def gfed_ch4_to_tg_month(ch4_g_m2_month, grid_area_m2):
    """Convert GFED CH4 from g CH4 m^-2 month^-1 to Tg CH4 month^-1 per cell."""
    return (
        np.asarray(ch4_g_m2_month, dtype=np.float64)
        * np.asarray(grid_area_m2, dtype=np.float64)
        * G_TO_TG
    )


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


def fit_annual_harmonic(clim_12):
    """Fit Q(m) = mean + B*sin(2*pi*m/12) + C*cos(2*pi*m/12)."""
    clim_12 = np.asarray(clim_12, dtype=np.float64)
    if clim_12.shape != (12,):
        raise ValueError(f"Expected 12 monthly values, got shape {clim_12.shape}")
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


def build_bb_seasonality():
    """Build GFED5 BB monthly climatology and harmonic coefficients by band."""
    monthly_by_band = {name: [] for name in [*BANDS.keys(), "Global"]}
    yearly_totals = {}

    for year in CLIMATOLOGY_YEARS:
        print(f"  Loading GFED5 BB {year}...")
        monthly_tg, lat = load_gfed_year(year)
        band_totals = aggregate_monthly_bands(monthly_tg, lat, BANDS)
        yearly_totals[str(year)] = {
            name: float(np.sum(values)) for name, values in band_totals.items()
        }
        for name, values in band_totals.items():
            monthly_by_band[name].append(values)

    bands = {}
    for name, year_arrays in monthly_by_band.items():
        climatology = np.vstack(year_arrays).mean(axis=0)
        bands[name] = fit_annual_harmonic(climatology)
        bands[name]["label"] = BANDS.get(name, {"label": "Global"})["label"]

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
            "yearly_total_Tg": yearly_totals,
        },
        "bands": bands,
        "site_assignment": site_assignment,
    }
    write_json(OUT_DATA_JSON, payload)
    return payload


def load_bb_signature(region, isotope, years=CLIMATOLOGY_YEARS):
    """Return mean and std for a BB isotope signature over selected years."""
    if isotope not in {"d13C", "dD"}:
        raise ValueError(f"Unsupported isotope: {isotope}")
    path = MODEL_DIR / "rel" / "data" / f"BB_{isotope}_{region}_MC.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing BB signature file: {path}")

    values = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            year = int(float(row[0]))
            if year in years:
                values.extend(float(x) for x in row[1:] if x != "")
    arr = np.asarray(values, dtype=np.float64)
    return {
        "mean": float(np.nanmean(arr)),
        "sigma": float(np.nanstd(arr)),
        "file": str(path),
    }


def load_bb_signatures_by_band():
    """Map each BB seasonal band to matching three-box BB isotope signatures."""
    out = {}
    for band, region in SIGNATURE_REGION_BY_BAND.items():
        out[band] = {
            "region": region,
            "d13C": load_bb_signature(region, "d13C"),
            "dD": load_bb_signature(region, "dD"),
        }
    return out


def rotate_month_index_to_midpoint(B_Q, C_Q):
    """Rotate month-index harmonic coefficients into Phase 2 midpoint basis."""
    delta = 2 * np.pi * 0.5 / 12.0
    cos_d = np.cos(delta)
    sin_d = np.sin(delta)
    return B_Q * cos_d + C_Q * sin_d, -B_Q * sin_d + C_Q * cos_d


def source_phasor(B_Q, C_Q, q_total, delta_source, delta_atm):
    """Compute isotope source phasor in permil from source seasonality."""
    B_mid, C_mid = rotate_month_index_to_midpoint(B_Q, C_Q)
    z_frac = complex(B_mid, C_mid) / q_total
    return (delta_source - delta_atm) * z_frac


def bulk_sink_ratio_for_alpha(alpha_13c_oh):
    """Return bulk sink amplitude ratio R for a candidate OH 13C KIE."""
    eps_13c = (
        F_OH * (alpha_13c_oh - 1)
        + F_CL * (ALPHA_13C_CL - 1)
        + F_SOIL * (ALPHA_13C_SOIL - 1)
        + F_STRAT * (ALPHA_13C_STRAT - 1)
    ) * 1000.0
    eps_d = (
        F_OH * (ALPHA_D_OH - 1)
        + F_CL * (ALPHA_D_CL - 1)
        + F_SOIL * (ALPHA_D_SOIL - 1)
        + F_STRAT * (ALPHA_D_STRAT - 1)
    ) * 1000.0
    return float(eps_13c / eps_d)


def saueressig_cantrell_r_band():
    """Return ordered R bounds corresponding to Saueressig and Cantrell alpha values."""
    r_saueressig = bulk_sink_ratio_for_alpha(ALPHA_13C_SAUERESSIG)
    r_cantrell = bulk_sink_ratio_for_alpha(ALPHA_13C_CANTRELL)
    lo = min(r_saueressig, r_cantrell)
    hi = max(r_saueressig, r_cantrell)
    return lo, hi, {
        "R_saueressig": r_saueressig,
        "R_cantrell": r_cantrell,
        "alpha_saueressig": ALPHA_13C_SAUERESSIG,
        "alpha_cantrell": ALPHA_13C_CANTRELL,
    }


def phasor_peak_month(z):
    """Return 0-indexed peak month for a B+iC phasor."""
    phase_rad = np.arctan2(z.imag, z.real)
    return float(((np.pi / 2 - phase_rad) / (2 * np.pi / 12.0)) % 12)


def apply_bb_to_wetland_corrected_decomposition(phase6_decomp, bb_src_13c, bb_src_dD):
    """Subtract BB source phasor from an existing wetland-corrected sink phasor."""
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
        "peak_src_bb_13c": phasor_peak_month(bb_src_13c),
        "peak_src_bb_dD": phasor_peak_month(bb_src_dD),
        "Z_sink_wetland_plus_bb_13c": [
            float(z_wet_bb_sink_13c.real),
            float(z_wet_bb_sink_13c.imag),
        ],
        "Z_sink_wetland_plus_bb_dD": [
            float(z_wet_bb_sink_dD.real),
            float(z_wet_bb_sink_dD.imag),
        ],
        "A_sink_wetland_plus_bb_13c": float(abs(z_wet_bb_sink_13c)),
        "A_sink_wetland_plus_bb_dD": float(abs(z_wet_bb_sink_dD)),
        "peak_sink_wetland_plus_bb_13c": phasor_peak_month(z_wet_bb_sink_13c),
        "peak_sink_wetland_plus_bb_dD": phasor_peak_month(z_wet_bb_sink_dD),
        "R_wetland_plus_bb": float(r),
    }


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
            "q_total_Tg_month": Q_TOTAL_TG_MONTH,
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
        phase6_decomp = phase6["sites"][code]
        corrected = apply_bb_to_wetland_corrected_decomposition(
            phase6_decomp, bb_src_13c, bb_src_dD
        )
        results["sites"][code] = {
            "band": band,
            "bb_signature_region": sig["region"],
            "bb_d13C_mean": sig["d13C"]["mean"],
            "bb_d13C_sigma": sig["d13C"]["sigma"],
            "bb_dD_mean": sig["dD"]["mean"],
            "bb_dD_sigma": sig["dD"]["sigma"],
            "bb_Q_mean_Tg_month": site_bb["Q_mean_Tg_month"],
            "bb_B_Q_Tg_month": bq,
            "bb_C_Q_Tg_month": cq,
            "bb_amplitude_Tg_month": site_bb["amplitude_Tg_month"],
            "bb_frac_amplitude": site_bb["frac_amplitude"],
            "bb_peak_month": site_bb["peak_month_index"],
            "R_obs": phase6_decomp["R_obs"],
            "R_wetland_only": phase6_decomp["R_corrected"],
            "A_src_wetland_13c": phase6_decomp["A_src_13c"],
            "A_src_wetland_dD": phase6_decomp["A_src_dD"],
            **corrected,
            "delta_R_vs_wetland_only": corrected["R_wetland_plus_bb"]
            - phase6_decomp["R_corrected"],
            "delta_A_src_13c_bb_minus_wetland": corrected["A_src_bb_13c"]
            - phase6_decomp["A_src_13c"],
            "delta_A_src_dD_bb_minus_wetland": corrected["A_src_bb_dD"]
            - phase6_decomp["A_src_dD"],
        }
    return results


def save_results(results):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_json = OUT_DIR / "bb_correction_results.json"
    write_json(out_json, results)

    out_csv = OUT_DIR / "bb_correction_summary.csv"
    fields = [
        "site",
        "band",
        "bb_signature_region",
        "R_obs",
        "R_wetland_only",
        "R_wetland_plus_bb",
        "delta_R_vs_wetland_only",
        "A_src_wetland_13c",
        "A_src_bb_13c",
        "A_src_wetland_dD",
        "A_src_bb_dD",
        "bb_Q_mean_Tg_month",
        "bb_amplitude_Tg_month",
        "bb_peak_month",
    ]
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for site, row in results["sites"].items():
            csv_row = {"site": site}
            csv_row.update({key: row.get(key, "") for key in fields if key != "site"})
            writer.writerow(csv_row)
    return out_json, out_csv


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


def plot_correction_comparison(results):
    sites = list(results["sites"].keys())
    x = np.arange(len(sites))
    r_obs = [results["sites"][s]["R_obs"] for s in sites]
    r_wet = [results["sites"][s]["R_wetland_only"] for s in sites]
    r_bb = [results["sites"][s]["R_wetland_plus_bb"] for s in sites]

    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    band_lo, band_hi, band_details = saueressig_cantrell_r_band()
    ax.axhspan(
        band_lo,
        band_hi,
        color="0.85",
        alpha=0.45,
        zorder=0,
        label=(
            "Saueressig-Cantrell band "
            f"(R={band_lo:.4f}-{band_hi:.4f})"
        ),
    )
    ax.axhline(
        band_details["R_saueressig"],
        color="C2",
        ls="--",
        lw=1.2,
        zorder=1,
        label=f"Saueressig (R={band_details['R_saueressig']:.4f})",
    )
    ax.axhline(
        band_details["R_cantrell"],
        color="C3",
        ls="--",
        lw=1.2,
        zorder=1,
        label=f"Cantrell (R={band_details['R_cantrell']:.4f})",
    )
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
