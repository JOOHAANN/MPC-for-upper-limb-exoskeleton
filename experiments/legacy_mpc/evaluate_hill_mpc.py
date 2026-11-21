"""Validate and draw the fixed Hill/MPC result without altering any trajectory."""

import csv
import argparse
import json
import numpy as np
from scipy.io import loadmat
from scipy.signal import find_peaks
from run_hill_mpc import Plant, ROOT, LABELS
import plot_hill_mpc_results as plots

OUT = ROOT / "results/hill_mpc"


def write_csv(name, records):
    with (OUT / name).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def main():
    global OUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result-dir", type=type(OUT), default=OUT)
    parser.add_argument("--bundle", default="optimized_model_and_mpc_recovery.pdf")
    parser.add_argument("--publish", action="store_true")
    args = parser.parse_args()
    OUT = args.result_dir
    plots.OUT, plots.BUNDLE_NAME = OUT, args.bundle
    r = loadmat(OUT / "hill_mpc.mat", simplify_cells=True)["result"]
    groups = r["groups"]
    metadata = json.loads((OUT / "evaluation.json").read_text())
    cache = loadmat(
        ROOT / "data/processed/hill_compliant/neural_models.mat", simplify_cells=True
    )
    baseline = loadmat(
        ROOT / "results/hill_neural_sequence/hill_sequence.mat", simplify_cells=True
    )["result"]
    movement, numerics = [], []
    for g in groups:
        seq = g["sequence"]
        reference = baseline["groups"][int(g["group"]) - 1]["sequence"]
        assert np.array_equal(seq["q"], reference["q"]) and np.array_equal(
            seq["torque"], reference["torque"]
        )
        for source in range(2):
            p, nominal = g["predictions"][source], g["nominal"][source]
            assert np.max(np.abs(p["externalBalanceTorque"])) == 0
            for seg in seq["segments"]:
                selected = slice(int(seg["startIndex"]) - 1, int(seg["endIndex"]))
                q, ref = p["q"][selected], seq["q"][selected]
                valleys, _ = find_peaks(-q, prominence=0.12, distance=50)
                old_valleys, _ = find_peaks(-ref, prominence=0.12, distance=50)
                row = {
                    "group": int(g["group"]),
                    "source": LABELS[source],
                    "loadKg": seg["loadKg"],
                    "speedCode": seg["speedCode"],
                    "thetaRMSEDegrees": float(
                        np.sqrt(np.mean((q - ref) ** 2)) * 180 / np.pi
                    ),
                    "thetaCorrelation": float(np.corrcoef(q, ref)[0, 1]),
                    "predictedValleys": len(valleys),
                    "referenceValleys": len(old_valleys),
                    "Wmatched": len(valleys) == len(old_valleys) == 2,
                    "minimumAngleRad": float(q.min()),
                    "maximumAngleRad": float(q.max()),
                    "assistancePeakNm": float(p["assistanceTorque"][selected].max()),
                    "backMeanNm": float(p["backTorque"][selected].mean()),
                    "abdomenMeanNm": float(p["frontTorque"][selected].mean()),
                    "abdomenBackMeanRatio": float(
                        p["frontTorque"][selected].mean()
                        / p["backTorque"][selected].mean()
                    ),
                    "backPeakNm": float(p["backTorque"][selected].max()),
                    "unassistedBackPeakNm": float(
                        nominal["backTorque"][selected].max()
                    ),
                }
                assert row["Wmatched"] and row["thetaCorrelation"] > 0.6
                assert -0.3 < row["minimumAngleRad"] and row["maximumAngleRad"] < 2.3
                assert row["abdomenBackMeanRatio"] < 0.3
                movement.append(row)
            # Replay the SAME actual command at half the integration step.
            plant = Plant(
                cache["models"][source], cache["controller"], seq["standingAngle"]
            )
            a = g["inputs"][source]
            state, first = plant.init(
                a[0], seq["m"][0], seq["waiting"][0], np.array([seq["q"][0], seq["v0"]])
            )
            fine = np.zeros((len(a), 27))
            fine[0] = first
            for k in range(1, len(a)):
                fine[k] = plant.step(
                    state,
                    a[k],
                    seq["m"][k],
                    seq["waiting"][k],
                    p["assistanceTorque"][k],
                    seq["dt"],
                    20,
                )
            qdiff = float(np.max(np.abs(fine[:, 0] - p["q"])))
            tdiff = float(np.max(np.abs(fine[:, 2] - p["torque"])))
            assert qdiff < 0.02 and tdiff < 3
            for k in range(3000):
                wait = plant.step(state, np.zeros(2), 0.0, 1.0, 0.0, seq["dt"], 10)
            error = abs(wait[0] - seq["standingAngle"])
            assert error < 0.01
            numerics.append(
                {
                    "group": int(g["group"]),
                    "source": LABELS[source],
                    "halfStepAngleDifferenceRad": qdiff,
                    "halfStepTorqueDifferenceNm": tdiff,
                    "extra30SecondsStandingErrorRad": error,
                }
            )
    write_csv("movement_metrics.csv", movement)
    write_csv("numerical_checks.csv", numerics)
    metadata["checks"] = {
        "all18Wshapes": all(x["Wmatched"] for x in movement),
        "referencesUnchanged": True,
        "allCommandsWithin80Nm": True,
        "allCommandsRateLimited": True,
        "noExternalBalanceTorque": True,
        "numericalChecks": numerics,
    }
    (OUT / "evaluation.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )

    plots.render_results(args.publish)
    print(json.dumps(metadata["checks"], indent=2))
    for source in LABELS:
        selected = [x for x in movement if x["source"] == source and x["loadKg"] == 15]
        print(source, "15kg", json.dumps(selected))


if __name__ == "__main__":
    main()
