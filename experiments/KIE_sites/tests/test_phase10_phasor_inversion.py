import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.KIE_sites.analysis.phase10_phasor_inversion import (  # noqa: E402
    invert_sink_phasors,
    phasor_residual,
)


def test_invert_sink_phasors_recovers_known_ratio_and_alpha():
    sink_13c = complex(0.2, -0.1)
    sink_dD = complex(4.0, -2.0)

    result = invert_sink_phasors(sink_13c, sink_dD)

    assert np.isclose(result["R_sink"], 0.05)
    assert result["alpha_13c_oh"] > 1.0


def test_phasor_residual_subtracts_source_and_sink_from_observed():
    observed = complex(2.0, 3.0)
    source = complex(-1.0, 0.5)
    sink = complex(3.0, 2.5)

    residual = phasor_residual(observed, source, sink)

    assert residual == complex(0.0, 0.0)
