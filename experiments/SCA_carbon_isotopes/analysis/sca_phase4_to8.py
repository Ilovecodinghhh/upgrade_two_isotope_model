#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 4-8 analysis outputs for the SCA carbon-isotope experiment."""

from __future__ import annotations

from io import StringIO
from pathlib import Path
import argparse

import numpy as np
import pandas as pd
import requests
from scipy import stats

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


PROJECT_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"
FIGURES_DIR = PROJECT_ROOT / "figures"
RAW_DATA_DIR = PROJECT_ROOT / "SCA_noaa_ch4_isotopes"
NOAA_SITE_METADATA_URL = "https://gml.noaa.gov/dv/site/"
DEFAULT_VALUE_COLUMNS = ("sca_harmonic", "sca_detrended_range", "sca_raw_range")


def fit_linear_trend(years: np.ndarray | pd.Series, values: np.ndarray | pd.Series) -> dict[str, float]:
    """Fit an ordinary least-squares line and return slope diagnostics."""
    x = np.asarray(years, dtype=float)
    y = np.asarray(values, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]
    n = int(len(x))
    if n < 2:
        return {
            "n_years": n,
            "slope_per_year": np.nan,
            "intercept": np.nan,
            "stderr": np.nan,
            "ci95_low": np.nan,
            "ci95_high": np.nan,
            "p_value": np.nan,
            "r2": np.nan,
        }

    slope, intercept = np.polyfit(x, y, deg=1)
    fitted = intercept + slope * x
    residuals = y - fitted
    ss_res = float(np.sum(residuals**2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    sxx = float(np.sum((x - np.mean(x)) ** 2))

    if n > 2 and sxx > 0:
        residual_variance = ss_res / (n - 2)
        stderr = float(np.sqrt(residual_variance / sxx))
        if stderr > 0:
            t_stat = float(slope / stderr)
            p_value = float(2.0 * stats.t.sf(abs(t_stat), df=n - 2))
        else:
            p_value = 0.0 if abs(slope) > 0 else 1.0
        t_crit = float(stats.t.ppf(0.975, df=n - 2))
        ci95_low = float(slope - t_crit * stderr)
        ci95_high = float(slope + t_crit * stderr)
    else:
        stderr = np.nan
        p_value = np.nan
        ci95_low = np.nan
        ci95_high = np.nan

    return {
        "n_years": n,
        "slope_per_year": float(slope),
        "intercept": float(intercept),
        "stderr": stderr,
        "ci95_low": ci95_low,
        "ci95_high": ci95_high,
        "p_value": p_value,
        "r2": float(r2) if np.isfinite(r2) else np.nan,
    }


def _pairwise_slopes(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    slopes = []
    for i in range(len(x) - 1):
        dx = x[i + 1 :] - x[i]
        dy = y[i + 1 :] - y[i]
        valid = dx != 0
        if np.any(valid):
            slopes.extend((dy[valid] / dx[valid]).tolist())
    return np.asarray(slopes, dtype=float)


def fit_sen_trend(
    years: np.ndarray | pd.Series,
    values: np.ndarray | pd.Series,
    bootstrap_iterations: int = 1000,
    random_seed: int = 42,
) -> dict[str, float]:
    """Fit a Theil-Sen median slope and bootstrap confidence interval."""
    x = np.asarray(years, dtype=float)
    y = np.asarray(values, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]
    n = int(len(x))
    slopes = _pairwise_slopes(x, y)
    if n < 2 or len(slopes) == 0:
        return {
            "n_years": n,
            "sen_slope_per_year": np.nan,
            "sen_intercept": np.nan,
            "sen_ci95_low": np.nan,
            "sen_ci95_high": np.nan,
        }

    sen_slope = float(np.median(slopes))
    sen_intercept = float(np.median(y - sen_slope * x))
    bootstrap_slopes = []
    rng = np.random.default_rng(random_seed)
    for _ in range(int(bootstrap_iterations)):
        indices = rng.integers(0, n, size=n)
        sample_slopes = _pairwise_slopes(x[indices], y[indices])
        if len(sample_slopes):
            bootstrap_slopes.append(float(np.median(sample_slopes)))
    if bootstrap_slopes:
        ci95_low, ci95_high = np.percentile(bootstrap_slopes, [2.5, 97.5])
    else:
        ci95_low = np.nan
        ci95_high = np.nan

    return {
        "n_years": n,
        "sen_slope_per_year": sen_slope,
        "sen_intercept": sen_intercept,
        "sen_ci95_low": float(ci95_low),
        "sen_ci95_high": float(ci95_high),
    }


def compute_period_trends(
    period_inputs: pd.DataFrame,
    value_columns: tuple[str, ...] = DEFAULT_VALUE_COLUMNS,
    min_years: int = 3,
    sen_bootstrap_iterations: int = 1000,
) -> pd.DataFrame:
    """Compute site-period SCA trends for each requested SCA metric."""
    rows = []
    df = period_inputs.copy()
    df["usable"] = df["usable"].astype(bool)
    df = df.loc[df["usable"]].copy()
    for (period, site), group in df.groupby(["period", "site"], sort=True):
        group = group.sort_values("year")
        for metric in value_columns:
            if metric not in group.columns:
                continue
            metric_group = group.dropna(subset=[metric])
            if len(metric_group) < min_years:
                continue
            trend = fit_linear_trend(metric_group["year"], metric_group[metric])
            sen_trend = fit_sen_trend(
                metric_group["year"],
                metric_group[metric],
                bootstrap_iterations=sen_bootstrap_iterations,
            )
            rows.append(
                {
                    "site": site,
                    "period": period,
                    "metric": metric,
                    "start_year": int(metric_group["year"].min()),
                    "end_year": int(metric_group["year"].max()),
                    "n_years": int(trend["n_years"]),
                    "mean_sca": float(metric_group[metric].mean()),
                    "slope_per_year": trend["slope_per_year"],
                    "slope_per_decade": trend["slope_per_year"] * 10.0,
                    "ci95_low": trend["ci95_low"],
                    "ci95_high": trend["ci95_high"],
                    "ci95_low_per_decade": trend["ci95_low"] * 10.0,
                    "ci95_high_per_decade": trend["ci95_high"] * 10.0,
                    "p_value": trend["p_value"],
                    "r2": trend["r2"],
                    "sen_slope_per_year": sen_trend["sen_slope_per_year"],
                    "sen_slope_per_decade": sen_trend["sen_slope_per_year"] * 10.0,
                    "sen_intercept": sen_trend["sen_intercept"],
                    "sen_ci95_low": sen_trend["sen_ci95_low"],
                    "sen_ci95_high": sen_trend["sen_ci95_high"],
                    "sen_ci95_low_per_decade": sen_trend["sen_ci95_low"] * 10.0,
                    "sen_ci95_high_per_decade": sen_trend["sen_ci95_high"] * 10.0,
                }
            )
    return pd.DataFrame(rows).sort_values(["period", "metric", "site"]).reset_index(drop=True)


def normalize_noaa_site_metadata(
    noaa_table: pd.DataFrame,
    site_codes: list[str] | tuple[str, ...] | set[str],
    metadata_source_url: str = NOAA_SITE_METADATA_URL,
) -> pd.DataFrame:
    """Normalize the NOAA GML site table to the station codes in this experiment."""
    site_order = [str(site).upper() for site in site_codes]
    site_set = set(site_order)
    df = noaa_table.copy()
    df["site"] = df["Code"].astype(str).str.upper()
    df = df.loc[df["site"].isin(site_set)].copy()
    df = df.rename(
        columns={
            "Name": "site_name",
            "Country": "country",
            "Latitude": "latitude",
            "Longitude": "longitude",
            "Elevation (meters)": "elevation_m",
            "Project": "project",
        }
    )
    df["latitude"] = pd.to_numeric(df["latitude"], errors="coerce")
    df["longitude"] = pd.to_numeric(df["longitude"], errors="coerce")
    df["elevation_m"] = pd.to_numeric(df["elevation_m"], errors="coerce")
    df["hemisphere"] = np.where(df["latitude"] >= 0, "NH", "SH")
    df["latitude_band"] = np.select(
        [
            df["latitude"].abs() < 23.5,
            df["latitude"] >= 23.5,
            df["latitude"] <= -23.5,
        ],
        ["tropical", "northern_mid_high", "southern_mid_high"],
        default="unknown",
    )
    df["metadata_source_url"] = metadata_source_url
    df["site_order"] = df["site"].map({site: index for index, site in enumerate(site_order)})
    columns = [
        "site",
        "site_name",
        "country",
        "latitude",
        "longitude",
        "elevation_m",
        "hemisphere",
        "latitude_band",
        "project",
        "metadata_source_url",
    ]
    return df.sort_values("site_order")[columns].reset_index(drop=True)


def compare_with_legacy_trends(
    new_trends: pd.DataFrame,
    legacy_trends: pd.DataFrame,
    period: str = "2002_2022",
    metric: str = "sca_harmonic",
) -> pd.DataFrame:
    """Align new fixed-period trends with the legacy exploratory trend table."""
    new_trends = new_trends.copy()
    if "slope_per_decade" not in new_trends.columns:
        new_trends["slope_per_decade"] = new_trends["slope_per_year"] * 10.0
    if "p_value" not in new_trends.columns:
        new_trends["p_value"] = np.nan
    if "r2" not in new_trends.columns:
        new_trends["r2"] = np.nan
    new_subset = new_trends.loc[
        (new_trends["period"] == period) & (new_trends["metric"] == metric),
        ["site", "period", "metric", "slope_per_year", "slope_per_decade", "p_value", "r2"],
    ].copy()
    new_subset = new_subset.rename(
        columns={
            "slope_per_year": "new_slope_per_year",
            "slope_per_decade": "new_slope_per_decade",
            "p_value": "new_p_value",
            "r2": "new_r2",
        }
    )
    legacy_subset = legacy_trends.copy()
    legacy_subset["site"] = legacy_subset["Site"].astype(str).str.upper()
    legacy_subset = legacy_subset.rename(
        columns={
            "Slope_per_year": "legacy_slope_per_year",
            "P_value": "legacy_p_value",
            "R2": "legacy_r2",
        }
    )
    legacy_subset["legacy_slope_per_decade"] = legacy_subset["legacy_slope_per_year"] * 10.0
    legacy_subset["legacy_method"] = "legacy_ssa_reconstruction_range"
    comparison = new_subset.merge(
        legacy_subset[
            [
                "site",
                "legacy_slope_per_year",
                "legacy_slope_per_decade",
                "legacy_p_value",
                "legacy_r2",
                "legacy_method",
            ]
        ],
        on="site",
        how="inner",
    )
    comparison["slope_delta_new_minus_legacy"] = (
        comparison["new_slope_per_year"] - comparison["legacy_slope_per_year"]
    )
    comparison["slope_delta_decade_new_minus_legacy"] = (
        comparison["new_slope_per_decade"] - comparison["legacy_slope_per_decade"]
    )
    return comparison.sort_values("site").reset_index(drop=True)


def _load_noaa_site_metadata_table(metadata_source_url: str) -> pd.DataFrame:
    response = requests.get(metadata_source_url, timeout=30)
    response.raise_for_status()
    tables = pd.read_html(StringIO(response.text))
    if not tables:
        raise ValueError(f"No site metadata tables found at {metadata_source_url}")
    return tables[0]


def _load_period_inputs(results_dir: Path) -> pd.DataFrame:
    files = sorted(results_dir.glob("period_inputs_*.csv"))
    if not files:
        raise FileNotFoundError(f"No period input tables found in {results_dir}")
    frames = [pd.read_csv(path) for path in files]
    return pd.concat(frames, ignore_index=True)


def _summarize_trends(trends: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (period, metric), group in trends.groupby(["period", "metric"], sort=True):
        slopes = group["slope_per_decade"].dropna()
        if slopes.empty:
            continue
        rows.append(
            {
                "period": period,
                "metric": metric,
                "n_sites": int(slopes.size),
                "median_slope_per_decade": float(slopes.median()),
                "q25_slope_per_decade": float(slopes.quantile(0.25)),
                "q75_slope_per_decade": float(slopes.quantile(0.75)),
                "n_positive": int((slopes > 0).sum()),
                "n_negative": int((slopes < 0).sum()),
                "n_p_lt_0_05": int((group["p_value"] < 0.05).sum()),
            }
        )
    return pd.DataFrame(rows).sort_values(["period", "metric"]).reset_index(drop=True)


def _summarize_trend_robustness(trends: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (period, metric), group in trends.groupby(["period", "metric"], sort=True):
        valid = group.dropna(subset=["slope_per_decade", "sen_slope_per_decade"]).copy()
        if valid.empty:
            continue
        ols = valid["slope_per_decade"]
        sen = valid["sen_slope_per_decade"]
        rows.append(
            {
                "period": period,
                "metric": metric,
                "n_sites": int(len(valid)),
                "median_ols_slope_per_decade": float(ols.median()),
                "median_sen_slope_per_decade": float(sen.median()),
                "median_sen_minus_ols_per_decade": float((sen - ols).median()),
                "n_same_sign": int((np.sign(ols) == np.sign(sen)).sum()),
                "n_ols_positive": int((ols > 0).sum()),
                "n_sen_positive": int((sen > 0).sum()),
                "n_sen_ci_excludes_zero": int(
                    (
                        (valid["sen_ci95_low_per_decade"] > 0)
                        | (valid["sen_ci95_high_per_decade"] < 0)
                    ).sum()
                ),
            }
        )
    return pd.DataFrame(rows).sort_values(["period", "metric"]).reset_index(drop=True)


def _plot_trends_by_latitude(trends: pd.DataFrame, figures_dir: Path) -> Path:
    data = trends.loc[trends["metric"] == "sca_harmonic"].copy()
    periods = [period for period in ["2002_2022", "2016_2022"] if period in set(data["period"])]
    fig, axes = plt.subplots(1, len(periods), figsize=(5.2 * len(periods), 4.3), sharey=True)
    if len(periods) == 1:
        axes = [axes]
    colors = {"NH": "#315b8a", "SH": "#b85c38"}
    for ax, period in zip(axes, periods):
        group = data.loc[data["period"] == period]
        for hemisphere, sub in group.groupby("hemisphere"):
            ax.scatter(
                sub["latitude"],
                sub["slope_per_decade"],
                s=45,
                color=colors.get(hemisphere, "#555555"),
                label=hemisphere,
                alpha=0.85,
            )
        ax.axhline(0, color="#777777", linewidth=0.9)
        ax.axvline(0, color="#cccccc", linewidth=0.8)
        ax.set_title(period.replace("_", "-"))
        ax.set_xlabel("Latitude")
        ax.grid(True, linewidth=0.4, alpha=0.35)
    axes[0].set_ylabel("Harmonic SCA trend (per mil per decade)")
    axes[-1].legend(frameon=False, loc="best")
    fig.tight_layout()
    out = figures_dir / "phase6_sca_harmonic_trends_by_latitude.png"
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def _plot_period_comparison(trends: pd.DataFrame, figures_dir: Path) -> Path:
    data = trends.loc[trends["metric"] == "sca_harmonic"].copy()
    pivot = data.pivot_table(index="site", columns="period", values="slope_per_decade", aggfunc="first")
    fig, ax = plt.subplots(figsize=(5.2, 4.8))
    if {"2002_2022", "2016_2022"}.issubset(pivot.columns):
        sub = pivot.dropna(subset=["2002_2022", "2016_2022"])
        ax.scatter(sub["2002_2022"], sub["2016_2022"], s=48, color="#36624a", alpha=0.85)
        for site, row in sub.iterrows():
            ax.text(row["2002_2022"], row["2016_2022"], site, fontsize=7, ha="left", va="bottom")
        low = float(np.nanmin(sub[["2002_2022", "2016_2022"]].to_numpy()))
        high = float(np.nanmax(sub[["2002_2022", "2016_2022"]].to_numpy()))
        pad = max((high - low) * 0.08, 0.01)
        ax.plot([low - pad, high + pad], [low - pad, high + pad], color="#999999", linewidth=0.9)
    ax.axhline(0, color="#cccccc", linewidth=0.8)
    ax.axvline(0, color="#cccccc", linewidth=0.8)
    ax.set_xlabel("2002-2022 trend (per mil per decade)")
    ax.set_ylabel("2016-2022 trend (per mil per decade)")
    ax.set_title("Period comparison for harmonic SCA trends")
    ax.grid(True, linewidth=0.4, alpha=0.35)
    fig.tight_layout()
    out = figures_dir / "phase6_sca_harmonic_period_comparison.png"
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def _plot_method_sensitivity(trends: pd.DataFrame, figures_dir: Path) -> Path:
    metrics = list(DEFAULT_VALUE_COLUMNS)
    periods = [period for period in ["2002_2022", "2016_2022", "2020_2022"] if period in set(trends["period"])]
    fig, axes = plt.subplots(1, len(periods), figsize=(5.0 * len(periods), 4.3), sharey=True)
    if len(periods) == 1:
        axes = [axes]
    for ax, period in zip(axes, periods):
        data = trends.loc[trends["period"] == period]
        positions = np.arange(len(metrics))
        box_data = [data.loc[data["metric"] == metric, "slope_per_decade"].dropna().to_numpy() for metric in metrics]
        ax.boxplot(box_data, positions=positions, widths=0.5, showfliers=False)
        for index, metric in enumerate(metrics):
            vals = box_data[index]
            if len(vals):
                jitter = np.linspace(-0.12, 0.12, len(vals)) if len(vals) > 1 else np.array([0.0])
                ax.scatter(np.full(len(vals), positions[index]) + jitter, vals, s=22, alpha=0.6, color="#444444")
        ax.axhline(0, color="#999999", linewidth=0.9)
        ax.set_xticks(positions)
        ax.set_xticklabels(["harmonic", "detrended range", "raw range"], rotation=25, ha="right")
        ax.set_title(period.replace("_", "-"))
        ax.grid(True, axis="y", linewidth=0.4, alpha=0.35)
    axes[0].set_ylabel("SCA trend (per mil per decade)")
    fig.tight_layout()
    out = figures_dir / "phase6_sca_method_sensitivity.png"
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def _plot_coverage_heatmap(eligibility: pd.DataFrame, metadata: pd.DataFrame, figures_dir: Path) -> Path:
    merged = eligibility.merge(metadata[["site", "latitude"]], on="site", how="left")
    site_order = (
        merged[["site", "latitude"]]
        .drop_duplicates()
        .sort_values("latitude", ascending=False)["site"]
        .tolist()
    )
    periods = [period for period in ["2002_2022", "2016_2022", "2020_2022"] if period in set(merged["period"])]
    matrix = np.full((len(site_order), len(periods)), np.nan)
    for row_index, site in enumerate(site_order):
        for col_index, period in enumerate(periods):
            match = merged.loc[(merged["site"] == site) & (merged["period"] == period), "usable_fraction"]
            if not match.empty:
                matrix[row_index, col_index] = float(match.iloc[0])

    fig_height = max(5.0, 0.26 * len(site_order))
    fig, ax = plt.subplots(figsize=(5.4, fig_height))
    image = ax.imshow(matrix, aspect="auto", vmin=0, vmax=1, cmap="viridis")
    ax.set_xticks(np.arange(len(periods)))
    ax.set_xticklabels([period.replace("_", "-") for period in periods], rotation=25, ha="right")
    ax.set_yticks(np.arange(len(site_order)))
    ax.set_yticklabels(site_order)
    for row_index in range(matrix.shape[0]):
        for col_index in range(matrix.shape[1]):
            value = matrix[row_index, col_index]
            if np.isfinite(value):
                text_color = "white" if value < 0.55 else "black"
                ax.text(col_index, row_index, f"{value:.2f}", ha="center", va="center", fontsize=7, color=text_color)
    ax.set_title("Usable-year fraction by fixed period")
    cbar = fig.colorbar(image, ax=ax, shrink=0.85)
    cbar.set_label("Usable fraction")
    fig.tight_layout()
    out = figures_dir / "phase6_period_coverage_heatmap.png"
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def run_phase4_to8_outputs(
    results_dir: Path | str = RESULTS_DIR,
    figures_dir: Path | str = FIGURES_DIR,
    metadata_table: pd.DataFrame | None = None,
    make_figures: bool = True,
    legacy_trends: pd.DataFrame | None | str = "auto",
) -> dict[str, Path]:
    """Generate Phase 4-8 trend, metadata, comparison, and figure outputs."""
    results_path = Path(results_dir)
    figures_path = Path(figures_dir)
    results_path.mkdir(parents=True, exist_ok=True)
    figures_path.mkdir(parents=True, exist_ok=True)

    period_inputs = _load_period_inputs(results_path)
    trends = compute_period_trends(period_inputs)
    site_codes = sorted(period_inputs["site"].astype(str).str.upper().unique())
    if metadata_table is None:
        metadata_table = _load_noaa_site_metadata_table(NOAA_SITE_METADATA_URL)
    metadata = normalize_noaa_site_metadata(metadata_table, site_codes=site_codes)
    trends_with_metadata = trends.merge(metadata, on="site", how="left")
    trend_summary = _summarize_trends(trends)
    trend_robustness_summary = _summarize_trend_robustness(trends)

    outputs = {
        "site_period_trends": results_path / "site_period_trends.csv",
        "site_metadata": results_path / "site_metadata.csv",
        "site_period_trends_with_metadata": results_path / "site_period_trends_with_metadata.csv",
        "period_metric_summary": results_path / "period_metric_summary.csv",
        "trend_robustness_summary": results_path / "trend_robustness_summary.csv",
    }
    trends.to_csv(outputs["site_period_trends"], index=False)
    metadata.to_csv(outputs["site_metadata"], index=False)
    trends_with_metadata.to_csv(outputs["site_period_trends_with_metadata"], index=False)
    trend_summary.to_csv(outputs["period_metric_summary"], index=False)
    trend_robustness_summary.to_csv(outputs["trend_robustness_summary"], index=False)

    if isinstance(legacy_trends, str) and legacy_trends == "auto":
        legacy_path = RAW_DATA_DIR / "sca_trend_summary.csv"
        legacy_trends = pd.read_csv(legacy_path) if legacy_path.exists() else None
    if legacy_trends is not None:
        comparison = compare_with_legacy_trends(trends, legacy_trends)
        outputs["legacy_trend_comparison"] = results_path / "legacy_sca_trend_comparison.csv"
        comparison.to_csv(outputs["legacy_trend_comparison"], index=False)

    if make_figures:
        outputs["figure_trends_by_latitude"] = _plot_trends_by_latitude(trends_with_metadata, figures_path)
        outputs["figure_period_comparison"] = _plot_period_comparison(trends_with_metadata, figures_path)
        outputs["figure_method_sensitivity"] = _plot_method_sensitivity(trends_with_metadata, figures_path)
        eligibility_path = results_path / "site_period_eligibility.csv"
        if eligibility_path.exists():
            eligibility = pd.read_csv(eligibility_path)
            outputs["figure_coverage_heatmap"] = _plot_coverage_heatmap(eligibility, metadata, figures_path)

    return outputs


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run SCA carbon-isotope Phase 4-8 analysis outputs.")
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--figures-dir", type=Path, default=FIGURES_DIR)
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args(argv)

    outputs = run_phase4_to8_outputs(
        results_dir=args.results_dir,
        figures_dir=args.figures_dir,
        make_figures=not args.no_figures,
    )
    for name, path in outputs.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
