import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.KIE_sites.analysis.make_manuscript_figures import (  # noqa: E402
    ganesan_banded_shift,
    plot_ganesan_banded_shifts,
    plot_mass_conserving_slice,
)


def test_mass_conserving_slice_marks_main_mc_reference_line():
    phase14 = {
        "mass_conserving": {
            "tropics_only": {
                "fractions": [0.0, 0.1],
                "summary_by_scenario": {
                    "tr_only_00": {"mean_alpha_13C_OH": 1.0052},
                    "tr_only_10": {"mean_alpha_13C_OH": 1.0039},
                },
            }
        }
    }
    fig, ax = plt.subplots()

    plot_mass_conserving_slice(
        ax,
        phase14,
        "tropics_only",
        "Tropics only",
        "Tropics response weight",
        "#2ca02c",
        reference_alpha=1.0046,
    )

    matching = [
        line
        for line in ax.lines
        if line.get_label() == "Main SH MC central"
        and np.allclose(line.get_ydata(), [1.0046, 1.0046])
    ]
    plt.close(fig)

    assert matching


def test_ganesan_banded_shift_is_zero_when_source_signature_is_unchanged():
    path = (
        ROOT
        / "experiments"
        / "KIE_sites"
        / "results"
        / "phase6_phasor"
        / "phasor_results.json"
    )
    phase6 = json.loads(path.read_text(encoding="utf-8"))

    assert ganesan_banded_shift(phase6["sites"]["CGO"]) == 0.0
    assert ganesan_banded_shift(phase6["sites"]["SPO"]) == 0.0


def test_ganesan_shift_panel_uses_delta_axis_without_lab_range():
    fig, ax = plt.subplots()
    plot_ganesan_banded_shifts(
        ax,
        ["CGO", "SPO"],
        [0.0, 0.0],
    )

    assert ax.get_xlabel() == "Change in corrected ratio, ΔR"
    assert not ax.patches
    assert {line.get_label() for line in ax.lines} >= {
        "Uniform baseline (ΔR = 0)",
        "Banded d13C shift",
    }
    plt.close(fig)
