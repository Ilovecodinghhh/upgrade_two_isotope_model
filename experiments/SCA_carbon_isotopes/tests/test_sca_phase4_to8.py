from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.SCA_carbon_isotopes.analysis.sca_phase4_to8 import (  # noqa: E402
    compare_with_legacy_trends,
    compute_period_trends,
    fit_linear_trend,
    fit_sen_trend,
    normalize_noaa_site_metadata,
    run_phase4_to8_outputs,
)


def test_fit_linear_trend_reports_slope_and_uncertainty_fields():
    years = np.array([2000, 2001, 2002, 2003, 2004])
    values = 1.0 + 0.2 * (years - 2000)

    trend = fit_linear_trend(years, values)

    assert trend["n_years"] == 5
    assert np.isclose(trend["slope_per_year"], 0.2)
    assert np.isclose(trend["intercept"], -399.0)
    assert trend["ci95_low"] <= trend["slope_per_year"] <= trend["ci95_high"]
    assert np.isclose(trend["r2"], 1.0)


def test_fit_sen_trend_reports_robust_slope_and_bootstrap_interval():
    years = np.array([2000, 2001, 2002, 2003, 2004])
    values = 1.0 + 0.2 * (years - 2000)

    trend = fit_sen_trend(years, values, bootstrap_iterations=200, random_seed=7)

    assert trend["n_years"] == 5
    assert np.isclose(trend["sen_slope_per_year"], 0.2)
    assert np.isclose(trend["sen_intercept"], -399.0)
    assert trend["sen_ci95_low"] <= trend["sen_slope_per_year"] <= trend["sen_ci95_high"]


def test_compute_period_trends_returns_one_row_per_site_period_metric():
    years = np.arange(2000, 2005)
    period_inputs = pd.DataFrame(
        {
            "period": "test",
            "site": "ABC",
            "year": years,
            "usable": True,
            "sca_harmonic": 0.5 + 0.10 * (years - 2000),
            "sca_detrended_range": 0.8 + 0.05 * (years - 2000),
        }
    )

    trends = compute_period_trends(
        period_inputs,
        value_columns=("sca_harmonic", "sca_detrended_range"),
    )

    by_metric = trends.set_index("metric")
    assert set(by_metric.index) == {"sca_harmonic", "sca_detrended_range"}
    assert np.isclose(by_metric.loc["sca_harmonic", "slope_per_year"], 0.10)
    assert np.isclose(by_metric.loc["sca_detrended_range", "slope_per_decade"], 0.50)
    assert np.isclose(by_metric.loc["sca_harmonic", "sen_slope_per_year"], 0.10)
    assert np.isclose(by_metric.loc["sca_detrended_range", "sen_slope_per_decade"], 0.50)
    assert "sen_ci95_low_per_decade" in trends.columns
    assert "sen_ci95_high_per_decade" in trends.columns
    assert by_metric.loc["sca_harmonic", "start_year"] == 2000
    assert by_metric.loc["sca_harmonic", "end_year"] == 2004


def test_normalize_noaa_site_metadata_filters_sites_and_assigns_bands():
    noaa_table = pd.DataFrame(
        {
            "Code": ["abc", "DEF", "GHI"],
            "Name": ["Equatorial", "Northern", "Southern"],
            "Country": ["X", "Y", "Z"],
            "Latitude": [1.0, 45.0, -60.0],
            "Longitude": [100.0, -120.0, 30.0],
            "Elevation (meters)": [10.0, 20.0, 30.0],
            "Project": ["Surface Flasks", "Surface Flasks", "Other"],
        }
    )

    metadata = normalize_noaa_site_metadata(noaa_table, site_codes=["ABC", "GHI"])

    assert metadata["site"].tolist() == ["ABC", "GHI"]
    by_site = metadata.set_index("site")
    assert by_site.loc["ABC", "hemisphere"] == "NH"
    assert by_site.loc["ABC", "latitude_band"] == "tropical"
    assert by_site.loc["GHI", "hemisphere"] == "SH"
    assert by_site.loc["GHI", "latitude_band"] == "southern_mid_high"
    assert by_site.loc["ABC", "metadata_source_url"].startswith("https://gml.noaa.gov/")


def test_compare_with_legacy_trends_aligns_sites_and_reports_slope_delta():
    new_trends = pd.DataFrame(
        {
            "site": ["ABC", "DEF"],
            "period": ["2002_2022", "2002_2022"],
            "metric": ["sca_harmonic", "sca_harmonic"],
            "slope_per_year": [0.02, -0.01],
        }
    )
    legacy_trends = pd.DataFrame(
        {
            "Site": ["ABC", "DEF"],
            "Slope_per_year": [0.01, -0.03],
            "P_value": [0.2, 0.1],
            "R2": [0.5, 0.4],
        }
    )

    comparison = compare_with_legacy_trends(
        new_trends,
        legacy_trends,
        period="2002_2022",
        metric="sca_harmonic",
    )

    by_site = comparison.set_index("site")
    assert np.isclose(by_site.loc["ABC", "slope_delta_new_minus_legacy"], 0.01)
    assert np.isclose(by_site.loc["DEF", "legacy_slope_per_decade"], -0.30)
    assert by_site.loc["ABC", "legacy_method"] == "legacy_ssa_reconstruction_range"


def test_run_phase4_to8_outputs_writes_trends_and_metadata_tables(tmp_path):
    results_dir = tmp_path / "results"
    figures_dir = tmp_path / "figures"
    results_dir.mkdir()
    rows = []
    for period in ["2002_2022", "2016_2022", "2020_2022"]:
        for year in range(2020, 2023):
            rows.append(
                {
                    "period": period,
                    "site": "ABC",
                    "year": year,
                    "usable": True,
                    "sca_raw_range": 0.7 + 0.02 * (year - 2020),
                    "sca_detrended_range": 0.6 + 0.03 * (year - 2020),
                    "sca_harmonic": 0.5 + 0.04 * (year - 2020),
                }
            )
        pd.DataFrame(rows[-3:]).to_csv(results_dir / f"period_inputs_{period}.csv", index=False)

    metadata_table = pd.DataFrame(
        {
            "Code": ["ABC"],
            "Name": ["Synthetic Station"],
            "Country": ["Testland"],
            "Latitude": [-40.0],
            "Longitude": [150.0],
            "Elevation (meters)": [12.0],
            "Project": ["Surface Flasks"],
        }
    )

    outputs = run_phase4_to8_outputs(
        results_dir=results_dir,
        figures_dir=figures_dir,
        metadata_table=metadata_table,
        make_figures=False,
        legacy_trends=None,
    )

    assert outputs["site_period_trends"].exists()
    assert outputs["site_metadata"].exists()
    assert outputs["site_period_trends_with_metadata"].exists()
    assert outputs["trend_robustness_summary"].exists()
    trends = pd.read_csv(outputs["site_period_trends"])
    assert set(trends["metric"]) == {"sca_raw_range", "sca_detrended_range", "sca_harmonic"}
    assert "sen_slope_per_decade" in trends.columns
    robustness = pd.read_csv(outputs["trend_robustness_summary"])
    assert set(["median_ols_slope_per_decade", "median_sen_slope_per_decade"]).issubset(
        robustness.columns
    )
    merged = pd.read_csv(outputs["site_period_trends_with_metadata"])
    assert merged.loc[0, "latitude"] == -40.0
    assert merged.loc[0, "hemisphere"] == "SH"
