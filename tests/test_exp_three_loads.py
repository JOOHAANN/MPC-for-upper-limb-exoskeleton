"""Excitation-only MAP input and missing-torque calibration regressions."""

import sys
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments/multi_subject"))
from excitation_only_map import excitation_signal, CHANNELS
import excitation_only_map as excitation_module
import train_exp_map as map_training
from prepare_exp_three_loads import supported
from prepare_subjects import combine
import fit_subject as fitting


class ExpThreeLoads(unittest.TestCase):
    def test_fresh_five_load_map_fits_with_excitation_files_only(self):
        with TemporaryDirectory() as temp:
            base = Path(temp)
            subject = "synthetic_A01"
            time = np.arange(0, 9, 0.01)
            for load in (5, 7.5, 10, 12.5, 15):
                for speed in (25, 45):
                    folder = (
                        base
                        / subject
                        / "CalibrationData"
                        / f"dynamics_stoop_{load:g}kg_{speed}"
                    )
                    folder.mkdir(parents=True)
                    signal = 0.04 + (0.2 + 0.004 * load) * (
                        np.exp(-(((time - 2.5) / 1.1) ** 2))
                        + np.exp(-(((time - 5.5) / 1.1) ** 2))
                    )
                    with (folder / "emgFilt.sto").open("w") as f:
                        f.write(
                            "Excitation\nendheader\ntime " + " ".join(CHANNELS) + "\n"
                        )
                        np.savetxt(
                            f,
                            np.c_[1000 + time, np.tile(signal[:, None], (1, 8))],
                            fmt="%.12f",
                        )
            with patch.object(excitation_module, "DATA", base), patch.object(
                map_training, "DATA", base
            ):
                gate, review = excitation_module.build_review(subject, base / "review")
                self.assertFalse(review["requiresIK"])
                self.assertEqual(len(review["conditions"]), 10)
                map_training.fit_subject_map(subject, gate, review, base / "model")
            self.assertTrue((base / "model" / (subject + "_map.mat")).exists())
            self.assertEqual(len(list(base.rglob("ik.sto"))), 0)
            self.assertEqual(len(list(base.rglob("Torque.sto"))), 0)

    def test_excitation_input_needs_no_ik_or_force_exports(self):
        with TemporaryDirectory() as temp:
            p = Path(temp) / "emgFilt.sto"
            times = 1000 + np.arange(1501) * 0.002
            with p.open("w") as f:
                f.write("Excitation\nendheader\ntime " + " ".join(CHANNELS) + "\n")
                np.savetxt(f, np.c_[times, np.full((len(times), 8), 0.2)], fmt="%.12f")
            out = excitation_signal(p)
            self.assertEqual(len(list(Path(temp).iterdir())), 1)
            self.assertLess(
                np.max(abs(out["activation"] - np.expm1(0.02) / np.expm1(0.1))), 1e-12
            )
            self.assertEqual(out["timeEpochSeconds"][0], 1000.0)

    def test_force_support_rejects_gaps_and_endpoint_extrapolation(self):
        clock = np.array([0.0, 0.01, 0.02, 0.5, 0.51])
        q = np.array([-0.01, 0.015, 0.1, 0.499, 0.5, 0.505, 0.52])
        np.testing.assert_array_equal(
            supported(clock, q), [False, True, False, False, True, True, False]
        )

    def test_missing_torque_cannot_enter_fit_loss(self):
        n = 201
        valid = np.ones(n, dtype=bool)
        valid[100:] = False
        r = dict(
            id="anchor",
            trialPath="trial",
            q=1.3 - 0.3 * np.sin(np.linspace(0, np.pi, n)),
            torque=np.where(valid, 100.0, np.nan),
            torqueReferenceValid=valid,
            back=np.full(n, 0.2),
            mappedBack=np.full(n, 0.2),
            front=np.full(n, 0.03),
            loadKg=5,
            speedCode=45,
            cadenceBPM=30,
            split="train",
            angleUsable=True,
            angleTorqueCorrelation=-0.5,
        )
        seq = combine([r], "fixture", "fixture")
        generator = dict(mu=np.zeros(10), scale=np.ones(10), beta=np.zeros(11))
        rec = fitting.record(seq, generator)
        p = (fitting.LOW + fitting.HIGH) / 2
        neural = (fitting.NLOW + fitting.NHIGH) / 2

        def fake_rollout(p, r, cfg, rec, source, steps):
            size = len(rec["seq"]["q"])
            out = np.zeros((size, 24))
            out[:, 0] = 1.4
            out[:, 2] = 90
            out[:, 3] = 100
            return out, False, np.full((size, 2), 0.03)

        with patch.object(fitting, "rollout", fake_rollout):
            a = fitting.residual(p, neural, np.zeros(1), [rec], 0, p, neural, 10)
            mask = ~np.asarray(rec["seq"]["torqueReferenceValid"], dtype=bool)
            rec["seq"]["torque"][mask] = 1e8
            b = fitting.residual(p, neural, np.zeros(1), [rec], 0, p, neural, 10)
            np.testing.assert_array_equal(a, b)
            self.assertTrue(np.isfinite(b).all())
            rec["seq"]["torque"][~mask] += 10
            c = fitting.residual(p, neural, np.zeros(1), [rec], 0, p, neural, 10)
            self.assertFalse(np.array_equal(b, c))


if __name__ == "__main__":
    unittest.main()
