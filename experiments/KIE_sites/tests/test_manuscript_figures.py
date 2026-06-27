import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.KIE_sites.analysis.make_manuscript_figures import (  # noqa: E402
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
