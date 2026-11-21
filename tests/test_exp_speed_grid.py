"""Guard grid priorities and reporting when CEINMS support is absent."""

import sys
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments/multi_subject"))
from prepare_subjects import combine
import fit_subject as fitting


class SpeedGrid(unittest.TestCase):
    def fixture(self, weight=None):
        n = 301
        cycle = dict(
            id="grid_case",
            trialPath="fixture",
            loadKg=10,
            speedCode=35,
            cadenceBPM=20,
            split="train",
            angleUsable=True,
            angleTorqueCorrelation=0.0,
            q=1.1 - 0.3 * np.sin(np.linspace(0, 2 * np.pi, n)),
            torque=np.full(n, np.nan),
            torqueReferenceValid=np.zeros(n, dtype=bool),
            back=np.full(n, 0.2),
            mappedBack=np.full(n, 0.2),
            front=np.full(n, 0.03),
        )
        if weight is not None:
            cycle["fitWeight"] = weight
        seq = combine([cycle], "fixture", "fixture")
        generator = dict(mu=np.zeros(10), scale=np.ones(10), beta=np.zeros(11))
        return fitting.record(seq, generator)

    @staticmethod
    def fake_rollout(p, neural, cfg, rec, source, steps, intent_decoder=None):
        n = len(rec["seq"]["q"])
        out = np.zeros((n, 24))
        out[:, 0] = 1.5
        out[:, 3] = 80
        out[:, 4] = 4
        out[:, 2] = 76
        return out, False, np.full((n, 2), 0.03)

    def test_priority_weights_change_loss_without_switching_model_or_inputs(self):
        old = self.fixture()
        priority = self.fixture(4.0)
        self.assertNotIn("fitWeight", old["segments"][0])
        self.assertEqual(priority["segments"][0]["fitWeight"], 4.0)
        np.testing.assert_array_equal(old["base"], priority["base"])
        p = (fitting.LOW + fitting.HIGH) / 2
        neural = (fitting.NLOW + fitting.NHIGH) / 2
        with patch.object(fitting, "rollout", self.fake_rollout):
            a = fitting.residual(p, neural, np.zeros(1), [old], 0, p, neural, 10)
            b = fitting.residual(p, neural, np.zeros(1), [priority], 0, p, neural, 10)
        # The final 44 entries are parameter priors, which are not case-weighted.
        np.testing.assert_allclose(b[:-44], 2 * a[:-44])
        np.testing.assert_array_equal(b[-44:], a[-44:])

    def test_absent_torque_is_unavailable_not_zero_error(self):
        rec = self.fixture()
        p = (fitting.LOW + fitting.HIGH) / 2
        models = [dict(parameters=p, neuralParameters=np.ones(10))]
        with patch.object(fitting, "rollout", self.fake_rollout):
            rows = fitting.score_rows(models, np.zeros(1), [rec])
        self.assertIsNone(rows[0]["torqueRMSENm"])
        self.assertEqual(rows[0]["torqueReferenceCoverageFraction"], 0.0)


if __name__ == "__main__":
    unittest.main()
