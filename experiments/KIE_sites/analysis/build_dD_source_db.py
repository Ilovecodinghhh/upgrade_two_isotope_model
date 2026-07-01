#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_dD_source_db.py — Site-specific δD-CH₄ source signature database
========================================================================

Replaces the single global δD_source = −310 ‰ assumption with latitude-
appropriate values for each of the 12 KIE_sites stations.

Method
------
1.  Read monthly precipitation δ²H from OIPC (already queried and saved).
2.  Compute annual-mean and growing-season-mean δ²H_precip for each site.
3.  Apply Douglas et al. (2021) regressions (Table S2) to predict the
    wetland CH₄ δ²H emitted at each latitude:
        δ²H-CH₄ = slope × δ²H-H₂O + intercept
    Three regression variants are computed; the primary is the
    *wetland, growing-season precipitation* line (highest R² for wetlands).
4.  Assign per-site uncertainty from regression SE + RMSE.
5.  Cross-validate against Douglas (2021) Table 1 zonal-mean values.
6.  Save structured JSON database + summary figure.

References
----------
- Douglas et al. (2021) Biogeosciences 18, 3505–3527, Table S2.
- Bowen & Revenaugh (2003) WRR 39, 1299 (OIPC).
- Waldron et al. (1999) GCA 63, 2237–2245.
- Chanton et al. (2006) JGR 111, G04004.

Output
------
    data/dD_source_database.json
    figures/fig6_dD_source_vs_latitude.png
"""

from pathlib import Path
import json, csv, math

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================================================================
# PATHS
# ============================================================================
EXPT_DIR  = Path(__file__).resolve().parent.parent          # experiments/KIE_sites/
DATA_DIR  = EXPT_DIR / "data"
FIG_DIR   = EXPT_DIR / "figures"
OIPC_CSV  = DATA_DIR / "oipc_precipitation_dD.csv"
OUT_JSON  = DATA_DIR / "dD_source_database.json"
OUT_FIG   = FIG_DIR  / "fig6_dD_source_vs_latitude.png"

FIG_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================================
# DOUGLAS 2021 REGRESSIONS  (Table S2)
# δ²H-CH₄ = slope × δ²H-H₂O + intercept
# ============================================================================
REGRESSIONS = {
    # ── Primary ──
    "wetland_gs": {
        "label": "Wetlands, growing-season δ²Hp",
        "slope": 0.705, "slope_se": 0.074,
        "intercept": -284.5, "intercept_se": 6.1,
        "R2": 0.633, "RMSE": 22.4, "n": 55,
    },
    # ── Cross-checks ──
    "wetland_ann": {
        "label": "Wetlands, annual δ²Hp",
        "slope": 0.575, "slope_se": 0.058,
        "intercept": -286.9, "intercept_se": 5.7,
        "R2": 0.651, "RMSE": 21.9, "n": 55,
    },
    "wetland_best": {
        "label": "Wetlands, best-estimate δ²H-H₂O",
        "slope": 0.51, "slope_se": 0.058,
        "intercept": -299.5, "intercept_se": 5.2,
        "R2": 0.59, "RMSE": 23.7, "n": 55,
    },
    "all_best": {
        "label": "All freshwater, best-estimate δ²H-H₂O",
        "slope": 0.44, "slope_se": 0.05,
        "intercept": -297.8, "intercept_se": 4.5,
        "R2": 0.418, "RMSE": 27.4, "n": 129,
    },
    # ── Waldron 1999 for comparison ──
    "waldron1999": {
        "label": "Waldron et al. (1999) global",
        "slope": 0.675, "slope_se": 0.10,
        "intercept": -284.0, "intercept_se": 6.0,
        "R2": 0.50, "RMSE": 26.3, "n": 51,
    },
}

# ============================================================================
# DOUGLAS 2021 TABLE 1 — ZONAL MEANS (for cross-validation)
# ============================================================================
DOUGLAS_ZONAL = {
    "wetland": [
        {"band": "<30°N",   "lat_range": (-90, 30), "dD_CH4": -301, "sigma": 15, "flux_Tg": 115},
        {"band": "30–60°N", "lat_range": (30, 60),  "dD_CH4": -324, "sigma": 14, "flux_Tg": 25},
        {"band": ">60°N",   "lat_range": (60, 90),  "dD_CH4": -374, "sigma": 10, "flux_Tg": 9},
    ],
    "inland_water": [
        {"band": "<30°N",   "lat_range": (-90, 30), "dD_CH4": -301, "sigma": 12, "flux_Tg": 80},
        {"band": "30–60°N", "lat_range": (30, 60),  "dD_CH4": -308, "sigma": 18, "flux_Tg": 64},
        {"band": ">60°N",   "lat_range": (60, 90),  "dD_CH4": -347, "sigma": 9,  "flux_Tg": 16},
    ],
    "rice": [
        {"band": "<30°N",   "lat_range": (-90, 30), "dD_CH4": -324, "sigma": 8, "flux_Tg": 19},
        {"band": "30–60°N", "lat_range": (30, 60),  "dD_CH4": -325, "sigma": 8, "flux_Tg": 12},
    ],
}

# Non-wetland source δD from Sherwood et al. (2017) and Douglas (2021) Table 1
NON_WETLAND_SOURCES = {
    "enteric_fermentation": {"dD_CH4": -308, "sigma": 28, "flux_Tg": 111},
    "landfills":            {"dD_CH4": -297, "sigma":  6, "flux_Tg": 65},
    "coal_mining":          {"dD_CH4": -232, "sigma":  5, "flux_Tg": 42},
    "oil_and_gas":          {"dD_CH4": -189, "sigma":  2, "flux_Tg": 79},
    "biomass_burning":      {"dD_CH4": -211, "sigma": 15, "flux_Tg": 17},
    "termites":             {"dD_CH4": -343, "sigma": 50, "flux_Tg":  9},
    "permafrost":           {"dD_CH4": -374, "sigma": 15, "flux_Tg":  1},
}


# ============================================================================
# GROWING-SEASON DEFINITION
# ============================================================================
def growing_season_months(lat: float) -> list[int]:
    """Return month indices (1-12) for the growing season at a given latitude.

    Approximation: months where mean temperature > 0 °C.
    - Tropics (|lat| < 30°): all months
    - NH mid/high (lat ≥ 30°):  Apr–Oct  (months 4–10)
    - SH mid/high (lat ≤ −30°): Oct–Apr  (months 10,11,12,1,2,3,4)
    """
    if abs(lat) < 30:
        return list(range(1, 13))
    elif lat >= 30:
        return list(range(4, 11))     # Apr–Oct
    else:
        return [10, 11, 12, 1, 2, 3, 4]  # Oct–Apr


# ============================================================================
# LOAD OIPC DATA
# ============================================================================
def load_oipc() -> list[dict]:
    """Read OIPC precipitation δ²H from CSV, return list of site dicts."""
    sites = []
    with open(OIPC_CSV) as f:
        reader = csv.DictReader(
            (row for row in f if not row.startswith("#")),
        )
        for row in reader:
            months = {m: float(row[m]) for m in
                      ["Jan","Feb","Mar","Apr","May","Jun",
                       "Jul","Aug","Sep","Oct","Nov","Dec"]}
            sites.append({
                "code": row["site"].strip(),
                "lat":  float(row["lat"]),
                "lon":  float(row["lon"]),
                "elev_m": float(row["elev_m"]),
                "monthly_dD_precip": months,
            })
    return sites


# ============================================================================
# COMPUTE MEANS
# ============================================================================
MONTH_NAMES = ["Jan","Feb","Mar","Apr","May","Jun",
               "Jul","Aug","Sep","Oct","Nov","Dec"]

def compute_means(site: dict) -> dict:
    """Compute annual and growing-season mean δ²H-precip for one site."""
    vals = [site["monthly_dD_precip"][m] for m in MONTH_NAMES]
    annual_mean = np.mean(vals)

    gs_months = growing_season_months(site["lat"])
    gs_vals = [vals[m - 1] for m in gs_months]
    gs_mean = np.mean(gs_vals)

    return {
        "annual_mean_dD_precip": round(float(annual_mean), 1),
        "growing_season_mean_dD_precip": round(float(gs_mean), 1),
        "growing_season_months": gs_months,
    }


# ============================================================================
# PREDICT δ²H-CH₄ VIA REGRESSION
# ============================================================================
def predict_dD_CH4(dD_H2O: float, reg: dict) -> dict:
    """Apply one regression to predict δ²H-CH₄ and its uncertainty.

    Prediction uncertainty combines:
      - regression slope × δ²H-H₂O uncertainty  (small for OIPC, ignore)
      - RMSE of the regression (dominant term)
      - intercept SE (systematic)
    Total 1σ ≈ sqrt(RMSE² + (slope_se × dD_H2O)² + intercept_se²)
    """
    pred = reg["slope"] * dD_H2O + reg["intercept"]
    # Prediction interval (1σ): dominated by scatter (RMSE)
    sigma = math.sqrt(
        reg["RMSE"]**2 +
        (reg["slope_se"] * dD_H2O)**2 +
        reg["intercept_se"]**2
    )
    return {
        "predicted_dD_CH4": round(pred, 1),
        "sigma": round(sigma, 1),
        "regression_used": reg["label"],
        "dD_H2O_input": round(dD_H2O, 1),
    }


# ============================================================================
# ZONAL CROSS-CHECK
# ============================================================================
def get_douglas_zonal(lat: float, source_type: str = "wetland") -> dict:
    """Look up the Douglas 2021 zonal mean for a given latitude."""
    zones = DOUGLAS_ZONAL.get(source_type, DOUGLAS_ZONAL["wetland"])
    for z in zones:
        lo, hi = z["lat_range"]
        if lo <= lat < hi:
            return {"band": z["band"], "dD_CH4": z["dD_CH4"], "sigma": z["sigma"]}
    # Fall through: SH sites map to <30°N band (Douglas only resolves NH bands)
    return {"band": "<30°N (SH default)", "dD_CH4": zones[0]["dD_CH4"],
            "sigma": zones[0]["sigma"]}


# ============================================================================
# MAIN DATABASE BUILDER
# ============================================================================
def build_database() -> dict:
    """Build the full site-specific δD source database."""
    sites = load_oipc()
    db = {"metadata": {}, "sites": {}, "non_wetland_sources": NON_WETLAND_SOURCES,
          "regressions": REGRESSIONS, "douglas_zonal": DOUGLAS_ZONAL}

    db["metadata"] = {
        "description": "Site-specific δ²H-CH₄ source signatures for KIE_sites experiment",
        "primary_regression": "wetland_gs (Douglas 2021, wetlands, growing-season δ²Hp)",
        "oipc_source": "OIPC v3.1, Bowen & Revenaugh (2003)",
        "date_built": "2026-05-20",
    }

    print(f"{'Site':<5} {'Lat':>6} {'δ²Hp_ann':>9} {'δ²Hp_gs':>8} "
          f"{'Pred_gs':>8} {'±1σ':>5} {'Zonal':>6} {'Δ':>5}")
    print("-" * 65)

    for site in sites:
        code = site["code"]
        lat  = site["lat"]
        means = compute_means(site)

        # Predictions from each regression
        predictions = {}
        for key, reg in REGRESSIONS.items():
            if "gs" in key:
                dD_input = means["growing_season_mean_dD_precip"]
            else:
                dD_input = means["annual_mean_dD_precip"]
            predictions[key] = predict_dD_CH4(dD_input, reg)

        # Primary: wetland growing-season
        primary = predictions["wetland_gs"]

        # Zonal cross-check
        zonal = get_douglas_zonal(lat)

        # Recommended value: use primary regression prediction
        # but note zonal mean as independent check
        delta = primary["predicted_dD_CH4"] - zonal["dD_CH4"]

        print(f"{code:<5} {lat:>+6.1f} "
              f"{means['annual_mean_dD_precip']:>+9.1f} "
              f"{means['growing_season_mean_dD_precip']:>+8.1f} "
              f"{primary['predicted_dD_CH4']:>+8.1f} "
              f"{primary['sigma']:>5.1f} "
              f"{zonal['dD_CH4']:>+6.0f} "
              f"{delta:>+5.0f}")

        # ── Recommendation logic ──
        # OIPC is preferred only at sites near real wetlands where
        # the regression is in-domain.  At remote/ocean/extreme-lat
        # stations the Douglas zonal mean better represents the
        # emission-weighted upstream source.
        use_oipc = code in ("BRW", "CBA")   # near extensive wetlands
        if use_oipc:
            rec_val   = primary["predicted_dD_CH4"]
            rec_sigma = primary["sigma"]
            rec_src   = "OIPC_regression"
        else:
            rec_val   = float(zonal["dD_CH4"])
            rec_sigma = float(zonal["sigma"])
            rec_src   = "Douglas_zonal"

        db["sites"][code] = {
            "lat": lat,
            "lon": site["lon"],
            "elev_m": site["elev_m"],
            "monthly_dD_precip": site["monthly_dD_precip"],
            "annual_mean_dD_precip": means["annual_mean_dD_precip"],
            "growing_season_mean_dD_precip": means["growing_season_mean_dD_precip"],
            "growing_season_months": means["growing_season_months"],
            "predictions": predictions,
            "oipc_prediction": {
                "dD_CH4": primary["predicted_dD_CH4"],
                "sigma": primary["sigma"],
                "method": "wetland_gs",
            },
            "douglas_zonal": {
                "dD_CH4": float(zonal["dD_CH4"]),
                "sigma": float(zonal["sigma"]),
                "band": zonal["band"],
            },
            "recommended": {
                "dD_CH4": rec_val,
                "sigma":  rec_sigma,
                "source": rec_src,
            },
        }

    return db


# ============================================================================
# FIGURE
# ============================================================================
def plot_database(db: dict) -> None:
    """Figure 6: Predicted δD-CH₄ source vs latitude, with Douglas zonal bands."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

    sites = db["sites"]
    codes = sorted(sites.keys(), key=lambda c: -sites[c]["lat"])
    lats = [sites[c]["lat"] for c in codes]

    # ── Panel (a): OIPC vs Zonal vs Recommended ──
    oipc_vals  = [sites[c]["oipc_prediction"]["dD_CH4"] for c in codes]
    oipc_sigs  = [sites[c]["oipc_prediction"]["sigma"] for c in codes]
    zonal_vals = [sites[c]["douglas_zonal"]["dD_CH4"] for c in codes]
    zonal_sigs = [sites[c]["douglas_zonal"]["sigma"] for c in codes]
    rec_vals   = [sites[c]["recommended"]["dD_CH4"] for c in codes]
    rec_sigs   = [sites[c]["recommended"]["sigma"] for c in codes]

    ax1.errorbar(lats, oipc_vals, yerr=oipc_sigs, fmt="o", color="C0",
                 ms=5, capsize=3, label="OIPC regression", alpha=0.5)
    ax1.errorbar(lats, zonal_vals, yerr=zonal_sigs, fmt="s", color="C1",
                 ms=5, capsize=3, label="Douglas zonal mean", alpha=0.5)
    ax1.scatter(lats, rec_vals, s=120, c="black", marker="*", zorder=5,
                label="Recommended")

    # Douglas zonal bands (wetlands)
    for z in DOUGLAS_ZONAL["wetland"]:
        lo, hi = z["lat_range"]
        ax1.fill_between([lo, hi],
                         z["dD_CH4"] - z["sigma"], z["dD_CH4"] + z["sigma"],
                         alpha=0.12, color="C1")
        ax1.hlines(z["dD_CH4"], lo, hi, colors="C1", linestyles="--", lw=1)

    # Global constant reference
    ax1.axhline(-310, color="red", ls=":", lw=1.5,
                label="Old global constant (−310 ‰)")

    for c, lat, rv in zip(codes, lats, rec_vals):
        src = sites[c]["recommended"]["source"]
        tag = "O" if src == "OIPC_regression" else "Z"
        ax1.annotate(f"{c} [{tag}]", (lat, rv),
                     textcoords="offset points", xytext=(6, -12), fontsize=6.5)

    ax1.set_xlabel("Latitude (°)", fontsize=10)
    ax1.set_ylabel("δ²H-CH₄ source (‰ VSMOW)", fontsize=10)
    ax1.set_title("(a) OIPC vs Douglas zonal — recommended marked ★", fontsize=10)
    ax1.legend(fontsize=7, loc="lower left")
    ax1.set_xlim(-100, 100)
    ax1.set_ylim(-600, -240)

    # ── Panel (b): Precipitation δ²H vs latitude ──
    ann_vals = [sites[c]["annual_mean_dD_precip"] for c in codes]
    gs_vals = [sites[c]["growing_season_mean_dD_precip"] for c in codes]

    ax2.plot(lats, ann_vals, "s-", color="C1", label="Annual mean δ²Hp", ms=6)
    ax2.plot(lats, gs_vals, "o-", color="C0", label="Growing-season mean δ²Hp", ms=6)

    for c, lat in zip(codes, lats):
        ax2.annotate(c, (lat, sites[c]["growing_season_mean_dD_precip"]),
                     textcoords="offset points", xytext=(6, 6), fontsize=7)

    ax2.set_xlabel("Latitude (°)", fontsize=10)
    ax2.set_ylabel("Precipitation δ²H (‰ VSMOW)", fontsize=10)
    ax2.set_title("(b) OIPC precipitation δ²H at KIE sites", fontsize=10)
    ax2.legend(fontsize=8)
    ax2.set_xlim(-100, 100)

    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=300)
    plt.close(fig)
    print(f"\n✓ Figure saved: {OUT_FIG}")


# ============================================================================
# MAIN
# ============================================================================
def main():
    print("=" * 65)
    print("Building site-specific δD-CH₄ source database")
    print("=" * 65)

    db = build_database()

    # Save JSON
    with open(OUT_JSON, "w") as f:
        json.dump(db, f, indent=2)
    print(f"\n✓ Database saved: {OUT_JSON}")

    # Figure
    plot_database(db)

    # ── Summary ──
    print("\n" + "=" * 75)
    print("RECOMMENDED δD SOURCE SIGNATURES")
    print("=" * 75)
    print(f"\n{'Site':<5} {'Lat':>6} {'OIPC':>7} {'Zonal':>7} "
          f"{'Rec':>7} {'±1σ':>5} {'Source':<18} {'Old':>6} {'Shift':>6}")
    print("-" * 75)
    for code in sorted(db["sites"], key=lambda c: -db["sites"][c]["lat"]):
        s = db["sites"][code]
        oipc  = s["oipc_prediction"]["dD_CH4"]
        zonal = s["douglas_zonal"]["dD_CH4"]
        rec   = s["recommended"]["dD_CH4"]
        sig   = s["recommended"]["sigma"]
        src   = s["recommended"]["source"]
        old   = -310.0
        print(f"{code:<5} {s['lat']:>+6.1f} {oipc:>+7.0f} {zonal:>+7.0f} "
              f"{rec:>+7.0f} {sig:>5.0f} {src:<18} {old:>+6.0f} {rec-old:>+6.0f}")


if __name__ == "__main__":
    main()
