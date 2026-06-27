from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.SCA_carbon_isotopes.analysis.sca_pipeline import (  # noqa: E402
    compute_period_eligibility,
    compute_site_year_sca,
    compute_year_coverage,
    ensure_project_structure,
    filter_sca_for_period,
    load_noaa_ch4c13_file,
    run_pipeline,
)


def test_load_noaa_ch4c13_file_parses_comments_and_monthly_values(tmp_path):
    path = tmp_path / "ch4c13_abc_surface-flask_7_sil_month.txt"
    path.write_text(
        "\n".join(
            [
                "# data_fields: site year month value",
                "abc 2020 1 -47.1",
                "abc 2020 2 -47.2",
            ]
        ),
        encoding="utf-8",
    )

    df = load_noaa_ch4c13_file(path)

    assert list(df.columns) == ["site", "year", "month", "date", "value", "source_file"]
    assert df["site"].tolist() == ["ABC", "ABC"]
    assert df["date"].dt.day.tolist() == [15, 15]
    assert np.allclose(df["value"], [-47.1, -47.2])


def test_compute_year_coverage_requires_months_and_quarters():
    dates = pd.to_datetime(
        [
            "2020-01-15",
            "2020-02-15",
            "2020-03-15",
            "2020-04-15",
            "2020-05-15",
            "2020-06-15",
            "2020-07-15",
            "2020-10-15",
            "2021-01-15",
            "2021-02-15",
            "2021-03-15",
            "2021-04-15",
            "2021-05-15",
            "2021-06-15",
            "2021-07-15",
        ]
    )
    df = pd.DataFrame(
        {
            "site": "ABC",
            "date": dates,
            "year": dates.year,
            "month": dates.month,
            "value": np.arange(len(dates), dtype=float),
            "source_file": "synthetic.txt",
        }
    )

    coverage = compute_year_coverage(df, min_months=8, min_quarters=3)

    by_year = coverage.set_index("year")
    assert bool(by_year.loc[2020, "usable"])
    assert not bool(by_year.loc[2021, "usable"])
    assert by_year.loc[2020, "n_months"] == 8
    assert by_year.loc[2020, "n_quarters"] == 4


def test_compute_site_year_sca_reports_raw_detrended_and_harmonic_amplitudes():
    months = np.tile(np.arange(1, 13), 3)
    years = np.repeat([2020, 2021, 2022], 12)
    dates = pd.to_datetime([f"{year}-{month:02d}-15" for year, month in zip(years, months)])
    amplitude = 0.25
    seasonal = amplitude * np.cos(2 * np.pi * (months - 7) / 12.0)
    trend = 0.01 * np.arange(len(months))
    values = -47.0 + trend + seasonal
    df = pd.DataFrame(
        {
            "site": "ABC",
            "date": dates,
            "year": years,
            "month": months,
            "value": values,
            "source_file": "synthetic.txt",
        }
    )

    sca = compute_site_year_sca(df, min_months=8, min_quarters=3)

    assert set(["sca_raw_range", "sca_detrended_range", "sca_harmonic"]).issubset(sca.columns)
    assert (sca["usable"]).all()
    assert (sca["sca_raw_range"] > sca["sca_detrended_range"]).all()
    assert np.allclose(sca["sca_harmonic"], 2 * amplitude, atol=0.04)


def test_compute_period_eligibility_uses_expected_years_not_observed_years():
    coverage = pd.DataFrame(
        {
            "site": ["ABC"] * 5 + ["DEF"] * 4,
            "year": [2000, 2001, 2002, 2003, 2004, 2000, 2001, 2003, 2004],
            "n_months": [12, 12, 12, 12, 7, 12, 12, 12, 12],
            "n_quarters": [4, 4, 4, 4, 3, 4, 4, 4, 4],
            "usable": [True, True, True, True, False, True, True, True, True],
        }
    )

    eligibility = compute_period_eligibility(
        coverage,
        periods={"test": (2000, 2004)},
        min_usable_fraction=0.70,
    )

    by_site = eligibility.set_index("site")
    assert by_site.loc["ABC", "expected_years"] == 5
    assert by_site.loc["ABC", "observed_years"] == 5
    assert by_site.loc["ABC", "usable_years"] == 4
    assert bool(by_site.loc["ABC", "eligible"])
    assert by_site.loc["DEF", "expected_years"] == 5
    assert by_site.loc["DEF", "observed_years"] == 4
    assert by_site.loc["DEF", "usable_years"] == 4
    assert bool(by_site.loc["DEF", "eligible"])
    assert np.isclose(by_site.loc["DEF", "usable_fraction"], 0.8)


def test_filter_sca_for_period_keeps_only_usable_years_from_eligible_sites():
    sca = pd.DataFrame(
        {
            "site": ["ABC", "ABC", "ABC", "DEF", "DEF"],
            "year": [2000, 2001, 2002, 2000, 2001],
            "usable": [True, False, True, True, True],
            "sca_harmonic": [0.1, 0.2, 0.3, 0.4, 0.5],
        }
    )
    eligibility = pd.DataFrame(
        {
            "site": ["ABC", "DEF"],
            "period": ["test", "test"],
            "start_year": [2000, 2000],
            "end_year": [2002, 2002],
            "eligible": [True, False],
        }
    )

    filtered = filter_sca_for_period(sca, eligibility, "test")

    assert filtered["site"].tolist() == ["ABC", "ABC"]
    assert filtered["year"].tolist() == [2000, 2002]
    assert filtered["period"].tolist() == ["test", "test"]


def test_run_pipeline_writes_reproducible_phase1_phase3_outputs(tmp_path):
    input_dir = tmp_path / "raw"
    input_dir.mkdir()
    lines = ["# data_fields: site year month value"]
    for year in [2020, 2021]:
        for month in range(1, 13):
            value = -47.0 + 0.2 * np.cos(2 * np.pi * (month - 7) / 12.0)
            lines.append(f"abc {year} {month} {value:.6f}")
    (input_dir / "ch4c13_abc_surface-flask_7_sil_month.txt").write_text(
        "\n".join(lines),
        encoding="utf-8",
    )

    created = ensure_project_structure(tmp_path)
    outputs = run_pipeline(
        input_dir=input_dir,
        results_dir=tmp_path / "results",
        periods={"two_years": (2020, 2021)},
    )

    assert {path.name for path in created} == {"data", "results", "figures", "docs"}
    assert all(path.exists() for path in created)
    assert all(path.exists() for path in outputs.values())
    yearly = pd.read_csv(outputs["site_yearly_sca"])
    assert yearly["usable"].all()
    assert np.allclose(yearly["sca_harmonic"], 0.4, atol=0.02)
    eligibility = pd.read_csv(outputs["period_eligibility"])
    assert eligibility.loc[0, "period"] == "two_years"
    assert bool(eligibility.loc[0, "eligible"])
    usable = pd.read_csv(outputs["site_yearly_sca_usable"])
    assert len(usable) == len(yearly)
    period_inputs = pd.read_csv(outputs["period_inputs_two_years"])
    assert period_inputs["period"].tolist() == ["two_years", "two_years"]


def test_legacy_scripts_do_not_use_absolute_desktop_data_path():
    analysis_root = ROOT / "experiments" / "SCA_carbon_isotopes" / "SCA_noaa_ch4_isotopes"
    for script_name in ["Amplitute_TimeSeries.py", "Monthly_variation.py"]:
        text = (analysis_root / script_name).read_text(encoding="utf-8")
        assert "D:\\0desktop" not in text
