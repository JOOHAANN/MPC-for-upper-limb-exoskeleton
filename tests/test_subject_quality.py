"""Critical quality-gate and bounded-controller regression checks.

Run: python -m unittest discover -s tests -p test_subject_quality.py
"""

import sys
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments/multi_subject"))
from quality_gate import native_tables

sys.path.insert(0, str(ROOT / "experiments/legacy_mpc"))
from run_hill_mpc import mpc_command


def fixture(folder, shifted=None, contradictory_duplicate=False):
    t = 1000.0 + np.arange(50) * 0.01
    force = 10.0 + np.arange(50) * 0.2
    arm = 0.04 + np.arange(50) * 0.0001
    values = {
        "ik.sto": np.sin(t),
        "Activations.sto": np.full(50, 0.3),
        "emg.sto": np.full(50, 0.2),
        "Torque.sto": force * arm,
        "MuscleForces.sto": force,
        "ma_L5_S1_Flex_Ext.sto": arm,
        "lmt.sto": np.full(50, 0.25),
    }
    for name, y in values.items():
        times = t + 0.02 if name == shifted else t
        data = np.c_[times, y]
        if name == "ma_L5_S1_Flex_Ext.sto":
            duplicate = data[25].copy()
            if contradictory_duplicate:
                duplicate[1] += 0.01
            data = np.insert(data, 26, duplicate, axis=0)
        with (folder / name).open("w", encoding="utf-8") as f:
            f.write("Fixture\nnRows=0\nnColumns=2\nendheader\ntime muscle\n")
            np.savetxt(f, data, fmt="%.12f")


class SubjectQuality(unittest.TestCase):
    def test_aligned_exports_and_identical_duplicate(self):
        with TemporaryDirectory() as tmp:
            fixture(Path(tmp))
            _, _, identity, exact, issues, repairs = native_tables(Path(tmp))
            self.assertEqual(issues, [])
            self.assertEqual(exact, 1.0)
            self.assertLess(identity["maximumErrorNm"], 1e-9)
            self.assertIn(
                "remove_identical_duplicate_rows:ma_L5_S1_Flex_Ext.sto", repairs
            )

    def test_shifted_torque_fails_even_on_same_regular_clock_grid(self):
        with TemporaryDirectory() as tmp:
            fixture(Path(tmp), shifted="Torque.sto")
            *_, issues, _ = native_tables(Path(tmp))
            self.assertIn("force_arm_torque_identity_failure", issues)

    def test_activation_emg_clock_mismatch(self):
        with TemporaryDirectory() as tmp:
            fixture(Path(tmp), shifted="Activations.sto")
            *_, issues, _ = native_tables(Path(tmp))
            self.assertIn("activation_emg_clock_mismatch", issues)

    def test_conflicting_duplicate_is_not_silently_averaged(self):
        with TemporaryDirectory() as tmp:
            fixture(Path(tmp), contradictory_duplicate=True)
            *_, issues, _ = native_tables(Path(tmp))
            self.assertIn(
                "conflicting_duplicate_timestamps:ma_L5_S1_Flex_Ext.sto", issues
            )

    def test_unattainable_torque_target_preserves_actuator_rate(self):
        class ToyPlant:
            def preview(self, state, controls, options):
                out = np.zeros((16, 27))
                out[:, 0] = 1.5
                out[:, 3] = 180 - 1.2 * np.repeat(controls, 4)
                return out

        state = np.zeros(34)
        state[28] = 40.0
        options = dict(
            points=16,
            blocks=4,
            predictionDt=0.025,
            predictionSubsteps=10,
            rateNmPerSec=200.0,
            controlPeriod=0.02,
            capNm=80.0,
            angleWeight=10.0,
            backTorqueWeight=64.0,
            velocityWeight=100.0,
            torqueConstraintPriority=True,
            predictionErrorGuardNm=[2.0] * 16,
        )
        command, flag, _ = mpc_command(
            ToyPlant(), state, state.copy(), np.full(4, 40.0), 100.0, options
        )
        self.assertEqual(flag, 1)
        self.assertAlmostEqual(command, 44.0, places=5)
        self.assertLess(180 - 1.2 * command, 180 - 1.2 * 40)
        self.assertGreater(180 - 1.2 * command, 100.0)

    def test_threshold_tracking_releases_unnecessary_assistance(self):
        class ToyPlant:
            def preview(self, state, controls, options):
                out = np.zeros((16, 27))
                out[:, 0] = 1.5
                out[:, 3] = 150.0 - np.repeat(controls, 4)
                return out

        state = np.zeros(34)
        state[28] = 60.0
        options = dict(
            points=16,
            blocks=4,
            predictionDt=0.025,
            predictionSubsteps=10,
            rateNmPerSec=200.0,
            controlPeriod=0.02,
            capNm=80.0,
            angleWeight=10.0,
            backTorqueWeight=64.0,
            velocityWeight=100.0,
            torqueConstraintPriority=True,
            predictionErrorGuardNm=[2.0] * 16,
            backTorqueUnderWeight=5.0,
            trackingHeadroomNm=3.0,
            assistanceEffortWeight=4.0,
        )
        command, flag, _ = mpc_command(
            ToyPlant(), state, state.copy(), np.full(4, 60.0), 100.0, options
        )
        self.assertEqual(flag, 1)
        self.assertGreaterEqual(command, 56.0 - 1e-7)
        self.assertLess(command, 60.0)

    def test_tracking_does_not_raise_an_already_low_natural_torque(self):
        class ToyPlant:
            def preview(self, state, controls, options):
                out = np.zeros((16, 27))
                out[:, 0] = 1.5
                out[:, 3] = 70.0 - np.repeat(controls, 4)
                return out

        state = np.zeros(34)
        options = dict(
            points=16,
            blocks=4,
            predictionDt=0.025,
            controlPeriod=0.02,
            rateNmPerSec=200.0,
            capNm=80.0,
            torqueConstraintPriority=True,
            predictionErrorGuardNm=[3.0] * 16,
            backTorqueUnderWeight=5.0,
            trackingHeadroomNm=3.0,
            assistanceEffortWeight=4.0,
        )
        command, flag, _ = mpc_command(
            ToyPlant(), state, state.copy(), np.zeros(4), 100.0, options
        )
        self.assertEqual(flag, 1)
        self.assertEqual(command, 0.0)


if __name__ == "__main__":
    unittest.main()
