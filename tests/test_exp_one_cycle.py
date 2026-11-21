"""Check that one complete M pair survives unrelated later data problems."""

from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments/multi_subject"))
from review_exp_one_cycle import find_m_bouts


class OneMCycle(unittest.TestCase):
    def setUp(self):
        self.t = np.arange(0, 12.01, 0.01)
        # One complete double peak followed by a repetition that ends early.
        self.y = np.interp(
            self.t,
            [0, 1, 3, 4.2, 5.5, 7, 9, 11, 12.01],
            [0.1, 0.1, 0.7, 0.28, 0.6, 0.1, 0.1, 0.65, 0.5],
        )
        self.six = np.tile(self.y[:, None], (1, 6))

    def test_keep_complete_pair_before_a_truncated_repetition(self):
        candidates = find_m_bouts(self.t, self.y, self.six, [self.t], {})
        clean = [c for c in candidates if c["cleanSourceWindow"]]
        self.assertTrue(clean)
        self.assertTrue(any(c["peakEpochSeconds"] == [3.0, 5.5] for c in clean))
        self.assertTrue(all(c["endEpochSeconds"] < 8 for c in clean))

    def test_gap_inside_selected_pair_is_rejected(self):
        clock = self.t[(self.t < 3.9) | (self.t > 4.4)]
        candidates = find_m_bouts(self.t, self.y, self.six, [clock], {})
        self.assertTrue(candidates)
        self.assertFalse(any(c["cleanSourceWindow"] for c in candidates))

    def test_conflicting_timestamp_is_excluded_locally(self):
        outside = find_m_bouts(
            self.t, self.y, self.six, [self.t], {}, excluded_intervals=[(9.7, 10.1)]
        )
        self.assertTrue(any(c["cleanSourceWindow"] for c in outside))
        inside = find_m_bouts(
            self.t, self.y, self.six, [self.t], {}, excluded_intervals=[(3.8, 4.2)]
        )
        self.assertFalse(any(c["cleanSourceWindow"] for c in inside))


if __name__ == "__main__":
    unittest.main()
