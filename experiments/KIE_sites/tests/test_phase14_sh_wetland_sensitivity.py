import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.KIE_sites.analysis import phase6_phasor  # noqa: E402
from experiments.KIE_sites.analysis.phase14_sh_wetland_sensitivity import (  # noqa: E402
    MASS_CONSERVING_NH_HIGH_FRACTIONS,
    MASS_CONSERVING_NH_HIGH_ONLY_FRACTIONS,
    MASS_CONSERVING_TROPICS_FRACTIONS,
    MASS_CONSERVING_TROPICS_ONLY_FRACTIONS,
    analyze_site_scenarios,
    build_total_source_phasor,
    load_all,
    mass_conserving_nh_high_only_scenario_definitions,
    mass_conserving_scenario_definitions,
    mass_conserving_tropics_only_scenario_definitions,
    rotate_phasor_delay,
    run_analysis,
    scenario_definitions,
)


def test_rotate_phasor_delay_moves_july_peak_three_months_later():
    b_lagged, c_lagged = rotate_phasor_delay(0.0, -1.0, lag_months=3.0)

    peak = phase6_phasor.phasor_peak_month(b_lagged, c_lagged)

    assert np.isclose(peak, 9.0)


def test_sh_only_reproduces_phase6_southern_hemisphere_deterministic_ratios():
    fits, wetland, phase6 = load_all()
    scenarios = scenario_definitions()

    for code in ("CGO", "SPO"):
        site_results = analyze_site_scenarios(
            code,
            fits[code],
            wetland,
            phase6["sites"][code],
            scenarios={"sh_only": scenarios["sh_only"]},
            include_mc=False,
        )

        assert np.isclose(
            site_results["scenarios"]["sh_only"]["R_corrected"],
            phase6["sites"][code]["R_corrected"],
            atol=1e-12,
        )


def test_full_nominal_source_amplitude_is_larger_than_sh_only():
    _, wetland, _ = load_all()
    scenarios = scenario_definitions()

    sh_only = build_total_source_phasor(wetland, scenarios["sh_only"]["components"])
    full_nominal = build_total_source_phasor(
        wetland, scenarios["full_nominal"]["components"]
    )

    assert full_nominal["A_src_total_13C"] > sh_only["A_src_total_13C"]
    assert full_nominal["A_src_total_dD"] > sh_only["A_src_total_dD"]
    assert len(full_nominal["source_components"]) == 4


def test_mass_conserving_scenarios_keep_source_weights_normalized():
    scenarios = mass_conserving_scenario_definitions()

    assert tuple(MASS_CONSERVING_NH_HIGH_FRACTIONS) == (0.04, 0.06, 0.08, 0.10)
    assert tuple(MASS_CONSERVING_TROPICS_FRACTIONS) == (
        0.00,
        0.10,
        0.20,
        0.30,
        0.40,
        0.50,
    )
    assert len(scenarios) == 24

    for name, scenario in scenarios.items():
        components = scenario["components"]
        weights = [component["scale"] for component in components]
        bands = [component["band"] for component in components]

        assert name.startswith("mc_nh")
        assert np.isclose(sum(weights), 1.0)
        assert all(weight >= 0.0 for weight in weights)
        assert bands == ["SH_extra", "Tropics", "NH_high"]

        nh_component = components[-1]
        assert nh_component["name"] == "NH_high_transport"
        assert np.isclose(nh_component["lag_months"], 2.8)


def test_mass_conserving_one_dimensional_slices_keep_weights_normalized():
    tropics_only = mass_conserving_tropics_only_scenario_definitions()
    nh_high_only = mass_conserving_nh_high_only_scenario_definitions()

    assert tuple(MASS_CONSERVING_TROPICS_ONLY_FRACTIONS) == (
        0.00,
        0.10,
        0.20,
        0.30,
        0.40,
        0.50,
    )
    assert tuple(MASS_CONSERVING_NH_HIGH_ONLY_FRACTIONS) == (
        0.00,
        0.04,
        0.06,
        0.08,
        0.10,
    )

    for scenario in tropics_only.values():
        components = scenario["components"]
        assert [component["band"] for component in components] == ["SH_extra", "Tropics"]
        assert np.isclose(sum(component["scale"] for component in components), 1.0)

    for scenario in nh_high_only.values():
        components = scenario["components"]
        assert [component["band"] for component in components] == ["SH_extra", "NH_high"]
        assert np.isclose(sum(component["scale"] for component in components), 1.0)
        assert np.isclose(components[-1]["lag_months"], 2.8)


def test_run_analysis_output_is_strict_json_without_mc():
    output = run_analysis(include_mc=False)

    encoded = json.dumps(output, allow_nan=False)
    decoded = json.loads(encoded)

    assert set(decoded["sites"]) == {"CGO", "SPO"}
    assert "full_nominal" in decoded["sh_summary_by_scenario"]
    assert "mc_nh04_tr50" in decoded["mass_conserving"]["summary_by_scenario"]
    assert "tr_only_50" in decoded["mass_conserving"]["tropics_only"]["summary_by_scenario"]
    assert "nh_high_only_10" in decoded["mass_conserving"]["nh_high_only"]["summary_by_scenario"]
    assert "CGO" in decoded["comparison_to_phase6_sh_only"]
