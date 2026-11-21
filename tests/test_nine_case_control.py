"""Guard contact physics and causality of the new activation forecaster."""

import sys
from pathlib import Path
import unittest
from types import SimpleNamespace
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments/multi_subject"))
from prepare_subjects import combine
from fit_subject import many, absolute_jacobian, rollout
from scipy.io import loadmat
from run_hill_mpc import configure_map_preview


class NineCaseControl(unittest.TestCase):
    def test_decoder_modes_share_identical_mechanics_at_zero_tonic_gain(self):
        package = loadmat(
            ROOT / "data/processed/hill_compliant/neural_models.mat",
            simplify_cells=True,
        )
        model = many(package["models"])[0]
        r = model["neuralParameters"].copy()
        r[9] = 0.0
        n = 60
        inputs = np.tile([0.2, 0.03], (n, 1))
        seq = dict(
            q=np.full(n, 1.5),
            m=np.full(n, 5.0),
            waiting=np.zeros(n),
            dt=0.01,
            standingAngle=1.5,
            v0=0.0,
        )
        rec = dict(
            seq=seq, base=[inputs, inputs], difference=[np.zeros(n), np.zeros(n)]
        )
        subtractive, bad_a, _ = rollout(
            model["parameters"], r, package["controller"], rec, 0, 6, "subtractive"
        )
        self.assertFalse(bad_a)
        for mode in ("divisive", "soft_tonic"):
            with self.subTest(mode=mode):
                result, failed, _ = rollout(
                    model["parameters"], r, package["controller"], rec, 0, 6, mode
                )
                self.assertFalse(failed)
                np.testing.assert_array_equal(subtractive, result)

    def test_normalized_derivative_remains_resolvable_near_both_bounds(self):
        # Model solves have finite resolution: tiny relative steps can vanish.
        def resolved_model(z):
            return np.round(np.array([3 * z[0], 0.7 * z[1]]), 7)

        point = np.array([1e-5, 1 - 1e-5])
        jac = absolute_jacobian(resolved_model, 2e-5)(point)
        np.testing.assert_allclose(jac, np.diag([3.0, 0.7]), atol=1e-7)

    def test_contact_schedule_keeps_unloaded_bending_outside_pick_and_place(self):
        n = 301
        r = dict(
            id="fixture",
            trialPath="fixture",
            loadKg=10,
            speedCode=45,
            cadenceBPM=20,
            q=np.ones(n),
            torque=np.ones(n),
            back=np.zeros(n),
            mappedBack=np.zeros(n),
            front=np.zeros(n),
            split="train",
            angleUsable=True,
            angleTorqueCorrelation=-1,
            loadContactStartSeconds=0.6,
            loadContactEndSeconds=2.2,
        )
        s = combine([r], "fixture", "fixture")
        seg = many(s["segments"])[0]
        start = int(seg["startIndex"]) - 1
        self.assertEqual(s["m"][start + 30], 0)
        self.assertGreater(s["m"][start + 120], 9.9)
        self.assertLess(s["m"][start + 290], 0.01)

    def fixture(self):
        n = 200
        t = np.arange(n) * 0.01
        planned = np.c_[0.2 + 0.1 * np.sin(t), 0.03 + 0.01 * np.cos(t)]
        seq = dict(
            dt=0.01, plannedMAPInputs=planned, m=np.ones(n) * 10, waiting=np.zeros(n)
        )
        rng = np.random.default_rng(3)
        predictor = dict(
            lags=[0, 2, 5, 10],
            scale=np.ones(8),
            coefficients=rng.normal(size=(8, 12)) * 0.05,
        )
        opt = dict(
            inputForecast="arx",
            activationResidualPredictor=predictor,
            points=6,
            predictionDt=0.025,
        )
        return seq, opt, planned

    def test_arx_forecast_cannot_read_future_experimental_activation(self):
        seq, opt, planned = self.fixture()
        observed = planned + 0.02
        plant = SimpleNamespace()
        configure_map_preview(plant, seq, observed, 80, opt)
        baseline = plant.map_activation.copy()
        changed = observed.copy()
        changed[81:] = 0.9
        configure_map_preview(plant, seq, changed, 80, opt)
        np.testing.assert_array_equal(plant.map_activation, baseline)

    def test_exact_generated_map_needs_no_residual_prediction(self):
        seq, opt, planned = self.fixture()
        plant = SimpleNamespace()
        configure_map_preview(plant, seq, planned, 80, opt)
        query = 80 + np.arange(1, 7) * 0.025 / 0.01
        expected = np.column_stack(
            [np.interp(query, np.arange(len(planned)), planned[:, j]) for j in (0, 1)]
        )
        np.testing.assert_allclose(plant.map_activation, expected, atol=1e-14)


if __name__ == "__main__":
    unittest.main()
