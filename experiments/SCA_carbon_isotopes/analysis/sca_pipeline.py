#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Reproducible SCA pipeline for NOAA CH4 carbon isotope monthly data."""

from __future__ import annotations

from pathlib import Path
import argparse

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent
RAW_DATA_DIR = PROJECT_ROOT / "SCA_noaa_ch4_isotopes"
RESULTS_DIR = PROJECT_ROOT / "results"
STANDARD_DIRS = ("data", "results", "figures", "docs")
DEFAULT_PERIODS = {
    "2002_2022": (2002, 2022),
    "2016_2022": (2016, 2022),
    "2020_2022": (2020, 2022),
}


MONTHLY_COLUMNS = ["site", "year", "month", "date", "value", "source_file"]


def ensure_project_structure(root: Path | str = PROJECT_ROOT) -> list[Path]:
    """Create the Phase 1 reproducible project directories."""
    base = Path(root)
    created = []
    for name in STANDARD_DIRS:
        path = base / name
        path.mkdir(parents=True, exist_ok=True)
        created.append(path)
    return created


def load_noaa_ch4c13_file(path: Path | str) -> pd.DataFrame:
    """Load one NOAA/INSTAAR CH4C13 monthly text file."""
    source = Path(path)
    df = pd.read_csv(
        source,
        comment="#",
        sep=r"\s+",
        names=["site", "year", "month", "value"],
        engine="python",
    )
    df = df.dropna(subset=["site", "year", "month", "value"]).copy()
    df["site"] = df["site"].astype(str).str.upper()
    df["year"] = pd.to_numeric(df["year"], errors="coerce").astype("Int64")
    df["month"] = pd.to_numeric(df["month"], errors="coerce").astype("Int64")
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df = df.dropna(subset=["year", "month", "value"]).copy()
    df["year"] = df["year"].astype(int)
    df["month"] = df["month"].astype(int)
    df["date"] = pd.to_datetime(
        {
            "year": df["year"],
            "month": df["month"],
            "day": np.full(len(df), 15, dtype=int),
        }
    )
    df["source_file"] = source.name
    return df[MONTHLY_COLUMNS].sort_values(["site", "date"]).reset_index(drop=True)


def load_all_noaa_ch4c13(input_dir: Path | str = RAW_DATA_DIR) -> pd.DataFrame:
    """Load all site-level NOAA CH4C13 monthly files in the input directory."""
    directory = Path(input_dir)
    files = sorted(directory.glob("ch4c13_*_surface-flask_7_sil_month.txt"))
    if not files:
        raise FileNotFoundError(f"No NOAA CH4C13 monthly files found in {directory}")
    frames = [load_noaa_ch4c13_file(path) for path in files]
    return pd.concat(frames, ignore_index=True).sort_values(["site", "date"]).reset_index(drop=True)


def _quarter(month: pd.Series) -> pd.Series:
    return ((month.astype(int) - 1) // 3) + 1


def compute_year_coverage(
    monthly: pd.DataFrame,
    min_months: int = 8,
    min_quarters: int = 3,
) -> pd.DataFrame:
    """Summarize whether each site-year has enough monthly coverage."""
    df = monthly.copy()
    df["quarter"] = _quarter(df["month"])
    coverage = (
        df.groupby(["site", "year"], as_index=False)
        .agg(
            n_months=("month", "nunique"),
            n_quarters=("quarter", "nunique"),
            first_month=("month", "min"),
            last_month=("month", "max"),
        )
        .sort_values(["site", "year"])
        .reset_index(drop=True)
    )
    coverage["usable"] = (
        (coverage["n_months"] >= min_months)
        & (coverage["n_quarters"] >= min_quarters)
    )
    return coverage


def _decimal_year(dates: pd.Series) -> np.ndarray:
    return dates.dt.year.to_numpy(dtype=float) + (dates.dt.month.to_numpy(dtype=float) - 0.5) / 12.0


def add_site_linear_detrended_value(monthly: pd.DataFrame) -> pd.DataFrame:
    """Add `value_detrended` after removing each site's linear trend."""
    frames = []
    for _, site_df in monthly.groupby("site", sort=True):
        site_df = site_df.sort_values("date").copy()
        if len(site_df) < 2:
            site_df["value_detrended"] = site_df["value"]
        else:
            x = _decimal_year(site_df["date"])
            y = site_df["value"].to_numpy(dtype=float)
            slope, intercept = np.polyfit(x, y, deg=1)
            fitted = intercept + slope * x
            site_df["value_detrended"] = y - fitted + float(np.nanmean(y))
        frames.append(site_df)
    return pd.concat(frames, ignore_index=True).sort_values(["site", "date"]).reset_index(drop=True)


def _harmonic_sca(months: pd.Series, values: pd.Series) -> float:
    """Return annual harmonic peak-to-trough amplitude for one site-year."""
    if len(values) < 3:
        return np.nan
    x = 2.0 * np.pi * (months.to_numpy(dtype=float) - 0.5) / 12.0
    design = np.column_stack([np.ones_like(x), np.sin(x), np.cos(x)])
    coeffs, _, _, _ = np.linalg.lstsq(design, values.to_numpy(dtype=float), rcond=None)
    _, b_sin, c_cos = coeffs
    return float(2.0 * np.hypot(b_sin, c_cos))


def compute_site_year_sca(
    monthly: pd.DataFrame,
    min_months: int = 8,
    min_quarters: int = 3,
) -> pd.DataFrame:
    """Compute raw, detrended, and harmonic SCA for every site-year."""
    df = add_site_linear_detrended_value(monthly)
    coverage = compute_year_coverage(df, min_months=min_months, min_quarters=min_quarters)

    rows = []
    for (site, year), group in df.groupby(["site", "year"], sort=True):
        raw = float(group["value"].max() - group["value"].min())
        detrended = float(group["value_detrended"].max() - group["value_detrended"].min())
        harmonic = _harmonic_sca(group["month"], group["value_detrended"])
        rows.append(
            {
                "site": site,
                "year": int(year),
                "sca_raw_range": raw,
                "sca_detrended_range": detrended,
                "sca_harmonic": harmonic,
            }
        )

    sca = pd.DataFrame(rows)
    out = coverage.merge(sca, on=["site", "year"], how="left")
    return out.sort_values(["site", "year"]).reset_index(drop=True)


def _normalize_periods(periods: dict[str, tuple[int, int]] | None) -> dict[str, tuple[int, int]]:
    if periods is None:
        periods = DEFAULT_PERIODS
    normalized = {}
    for name, bounds in periods.items():
        start_year, end_year = bounds
        start_year = int(start_year)
        end_year = int(end_year)
        if end_year < start_year:
            raise ValueError(f"Period {name!r} ends before it starts: {start_year}-{end_year}")
        normalized[str(name)] = (start_year, end_year)
    return normalized


def _safe_period_name(period: str) -> str:
    return "".join(char if char.isalnum() or char in ("-", "_") else "_" for char in str(period))


def compute_period_eligibility(
    coverage: pd.DataFrame,
    periods: dict[str, tuple[int, int]] | None = None,
    min_usable_fraction: float = 0.70,
    diagnostic_year_threshold: int = 5,
) -> pd.DataFrame:
    """Summarize whether each site has enough usable years in fixed periods."""
    normalized_periods = _normalize_periods(periods)
    df = coverage.copy()
    df["usable"] = df["usable"].astype(bool)
    rows = []
    for site, site_df in df.groupby("site", sort=True):
        for period, (start_year, end_year) in normalized_periods.items():
            period_mask = site_df["year"].between(start_year, end_year)
            period_df = site_df.loc[period_mask]
            expected_years = end_year - start_year + 1
            observed_years = int(period_df["year"].nunique())
            usable_years = int(period_df.loc[period_df["usable"], "year"].nunique())
            usable_fraction = usable_years / expected_years
            rows.append(
                {
                    "site": site,
                    "period": period,
                    "start_year": start_year,
                    "end_year": end_year,
                    "expected_years": expected_years,
                    "observed_years": observed_years,
                    "usable_years": usable_years,
                    "usable_fraction": usable_fraction,
                    "eligible": usable_fraction >= min_usable_fraction,
                    "diagnostic_only": expected_years < diagnostic_year_threshold,
                }
            )
    return pd.DataFrame(rows).sort_values(["period", "site"]).reset_index(drop=True)


def filter_sca_for_period(
    sca: pd.DataFrame,
    eligibility: pd.DataFrame,
    period: str,
    require_eligible: bool = True,
) -> pd.DataFrame:
    """Return usable site-year SCA rows within a fixed period."""
    period_eligibility = eligibility.loc[eligibility["period"] == period].copy()
    if period_eligibility.empty:
        raise ValueError(f"Period {period!r} is not present in eligibility table")

    start_years = period_eligibility["start_year"].unique()
    end_years = period_eligibility["end_year"].unique()
    if len(start_years) != 1 or len(end_years) != 1:
        raise ValueError(f"Period {period!r} has inconsistent year bounds")

    selected_sites = set(period_eligibility["site"])
    if require_eligible:
        selected_sites = set(period_eligibility.loc[period_eligibility["eligible"].astype(bool), "site"])

    out = sca.loc[
        sca["site"].isin(selected_sites)
        & sca["year"].between(int(start_years[0]), int(end_years[0]))
        & sca["usable"].astype(bool)
    ].copy()
    out.insert(0, "period", period)
    return out.sort_values(["site", "year"]).reset_index(drop=True)


def run_pipeline(
    input_dir: Path | str = RAW_DATA_DIR,
    results_dir: Path | str = RESULTS_DIR,
    min_months: int = 8,
    min_quarters: int = 3,
    periods: dict[str, tuple[int, int]] | None = None,
    min_usable_fraction: float = 0.70,
) -> dict[str, Path]:
    """Run Phase 1-3 outputs from raw monthly data."""
    ensure_project_structure(PROJECT_ROOT)
    results_path = Path(results_dir)
    results_path.mkdir(parents=True, exist_ok=True)
    normalized_periods = _normalize_periods(periods)
    monthly = load_all_noaa_ch4c13(input_dir)
    monthly_clean = add_site_linear_detrended_value(monthly)
    coverage = compute_year_coverage(monthly_clean, min_months=min_months, min_quarters=min_quarters)
    sca = compute_site_year_sca(monthly_clean, min_months=min_months, min_quarters=min_quarters)
    eligibility = compute_period_eligibility(
        coverage,
        periods=normalized_periods,
        min_usable_fraction=min_usable_fraction,
    )
    usable_sca = sca.loc[sca["usable"].astype(bool)].copy()

    outputs = {
        "monthly_clean": results_path / "site_monthly_clean.csv",
        "year_coverage": results_path / "site_year_coverage.csv",
        "site_yearly_sca": results_path / "site_yearly_sca.csv",
        "site_yearly_sca_usable": results_path / "site_yearly_sca_usable.csv",
        "period_eligibility": results_path / "site_period_eligibility.csv",
    }
    monthly_clean.to_csv(outputs["monthly_clean"], index=False)
    coverage.to_csv(outputs["year_coverage"], index=False)
    sca.to_csv(outputs["site_yearly_sca"], index=False)
    usable_sca.to_csv(outputs["site_yearly_sca_usable"], index=False)
    eligibility.to_csv(outputs["period_eligibility"], index=False)
    for period in normalized_periods:
        key = f"period_inputs_{_safe_period_name(period)}"
        outputs[key] = results_path / f"{key}.csv"
        filter_sca_for_period(sca, eligibility, period).to_csv(outputs[key], index=False)
    return outputs


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run the SCA carbon-isotope Phase 1-3 pipeline.")
    parser.add_argument("--input-dir", type=Path, default=RAW_DATA_DIR)
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--min-months", type=int, default=8)
    parser.add_argument("--min-quarters", type=int, default=3)
    parser.add_argument("--min-usable-fraction", type=float, default=0.70)
    args = parser.parse_args(argv)

    outputs = run_pipeline(
        input_dir=args.input_dir,
        results_dir=args.results_dir,
        min_months=args.min_months,
        min_quarters=args.min_quarters,
        min_usable_fraction=args.min_usable_fraction,
    )
    for name, path in outputs.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
