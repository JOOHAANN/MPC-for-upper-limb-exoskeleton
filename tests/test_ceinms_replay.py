"""Independent polynomial benchmarks for reference-replay numerics."""

from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(
    0, str(Path(__file__).resolve().parents[1] / "experiments/multi_subject")
)
from ceinms_force_replay import MuscleCurve, causal_velocity


class ReplayNumerics(unittest.TestCase):
    def test_cubic_curve_matches_analytic_polynomial(self):
        x = np.array([0, 0.12, 0.31, 0.58, 0.9, 1.4])
        polynomial = lambda q: 2 * q**3 - 3 * q**2 + 0.4 * q + 2
        curve = MuscleCurve(x, polynomial(x))
        q = np.linspace(x[0], x[-1], 101)
        np.testing.assert_allclose(curve(q), polynomial(q), atol=1e-12)
        np.testing.assert_allclose(
            curve.derivative(q), 6 * q**2 - 6 * q + 0.4, atol=1e-11
        )

    def test_rolling_velocity_retains_cubic_derivative_on_irregular_clock(self):
        t = np.cumsum(np.linspace(0.03, 0.12, 31))
        position = np.column_stack(
            [0.04 * t**3 - 0.02 * t**2 + 0.3, 0.1 * t**2 + 0.02 * t]
        )
        expected = np.column_stack([0.12 * t**2 - 0.04 * t, 0.2 * t + 0.02])
        actual = causal_velocity(t, position)
        np.testing.assert_allclose(actual[3:], expected[3:], atol=1e-10)


if __name__ == "__main__":
    unittest.main()
