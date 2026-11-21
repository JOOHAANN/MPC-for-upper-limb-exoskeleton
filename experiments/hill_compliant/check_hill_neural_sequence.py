"""Audit portable fixed models, complete MATLAB replay and extended standing.

The appended 30-second wait holds the final waiting excitation and zero load.
It continues the original states without resets. This is a numerical and
reproducibility check, not independent physiological or holdout validation.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat

import fit_hill_neural as neural
from fit_hill_compliant import ROOT, boundary_torque


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_table(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evaluation", type=Path, default=ROOT / "results/hill_neural_selected_evaluation")
    parser.add_argument("--matlab-replay", type=Path, default=ROOT / "results/hill_neural_sequence/hill_sequence.mat")
    args = parser.parse_args()
    cache = ROOT / "data/processed/hill_compliant/neural_models.mat"
    snapshots = ROOT / "docs/hill_neural_evaluation"
    metadata_file = snapshots / "export_metadata.json"
    metadata = json.loads(metadata_file.read_text(encoding="utf-8"))
    assert metadata["portableModelSHA256"] == sha(cache)
    screening = json.loads((snapshots / "screening.json").read_text(encoding="utf-8"))
    angle_policy = screening.get("anglePolicy", {"sequenceLower": 0., "sequenceUpper": 1.8})
    package = loadmat(cache, simplify_cells=True)
    models = package["models"]
    controller = np.r_[1.5, package["controller"]]
    original = loadmat(ROOT / "data/processed/hill_compliant/paired_sequences.mat", simplify_cells=True)["sequences"]
    matlab = loadmat(args.matlab_replay, simplify_cells=True)["result"]
    assert np.array_equal(np.asarray(matlab["options"]["controller"]), np.asarray(package["controller"]))
    with (snapshots / "fixed_parameter_audit.csv").open(encoding="utf-8", newline="") as stream:
        hashes = list(csv.DictReader(stream))
    assert np.array_equal(models[0]["parameters"][neural.fit.SHARED], models[1]["parameters"][neural.fit.SHARED])
    numerical, waiting = [], []
    for source, model in enumerate(models):
        vector = np.r_[model["parameters"], model["neuralParameters"]]
        digest = hashlib.sha256(vector.astype("<f8").tobytes()).hexdigest()
        cases = [r for r in hashes if int(r["source"]) == source + 1]
        assert len(cases) == 9 and all(r.get("fixedParameterSHA256", r.get("fixed41ParameterSHA256")) == digest for r in cases)
        controller_digest = hashlib.sha256(np.asarray(package["controller"]).astype("<f8").tobytes()).hexdigest()
        assert all(r.get("sharedStandingControllerSHA256", controller_digest) == controller_digest for r in cases)
        assert len({(int(float(r["loadKg"])), int(float(r["speedCode"]))) for r in cases}) == 9
        assert np.array_equal(matlab["models"][source]["parameters"], model["parameters"])
        assert np.array_equal(matlab["models"][source]["neuralParameters"], model["neuralParameters"])
        for group, trial in enumerate(original, start=1):
            bundle = np.load(args.evaluation / f"group{group}_source{source + 1}.npz")
            expected = bundle["output"]
            replay = matlab["groups"][group - 1]
            prediction = replay["predictions"][source]
            assert not prediction["failed"]
            for name, key in (("qReference", "q"), ("torqueReference", "torque")):
                assert np.array_equal(bundle[name], trial[key])
                assert np.array_equal(replay["sequence"][key], trial[key])
            assert np.array_equal(bundle["parameters"], model["parameters"])
            inputs_error = np.max(np.abs(replay["inputs"][source] - bundle["inputs"]))
            q_error = np.max(np.abs(prediction["q"] - expected[:, 0]))
            torque_error = np.max(np.abs(prediction["torque"] - expected[:, 2]))
            front_error = np.max(np.abs(prediction["frontTorque"] - expected[:, 4]))
            activation_error = np.max(np.abs(prediction["effectiveActivation"] - expected[:, 13:15]))
            intent_error = np.max(np.abs(prediction["motionIntent"] - expected[:, 23]))
            balance_error = np.max(np.abs(prediction["torque"] - prediction["backTorque"] + prediction["frontTorque"]))
            soft_error = np.max(np.abs(prediction["boundaryTorque"] - boundary_torque(model["parameters"], prediction["q"])))
            assert inputs_error < 1e-10 and q_error < 1e-6
            assert max(torque_error, front_error) < .001 and activation_error < 1e-6 and intent_error < 1e-9
            assert balance_error < 1e-8 and soft_error < 1e-8
            assert np.all(prediction["externalBalanceTorque"] == 0)
            numerical.append({"group": group, "source": source + 1, "inputMaxDifference": inputs_error,
                              "angleMaxDifferenceRad": q_error, "netTorqueMaxDifferenceNm": torque_error,
                              "abdominalTorqueMaxDifferenceNm": front_error, "activationMaxDifference": activation_error,
                              "intentMaxDifferenceRad": intent_error, "netTorqueBalanceResidualNm": balance_error,
                              "originalReferenceMaxDifference": 0, "parameterCount": len(vector), "fixedParameterSHA256": digest})

            # Re-run the ENTIRE sequence plus an additional wait. No resetting
            # at a movement boundary or at the original recording's end.
            extra = round(30 / trial["dt"])
            a = np.asfortranarray(np.vstack((bundle["inputs"], np.tile(bundle["inputs"][-1], (extra, 1)))))
            assert trial["waiting"][-1] and abs(trial["m"][-1]) < 1e-8
            mass = np.asfortranarray(np.r_[trial["m"], np.zeros(extra)])
            wait = np.asfortranarray(np.r_[trial["waiting"], np.ones(extra)], dtype=float)
            output = np.zeros((len(a), 24), order="F")
            cfg = controller.copy()
            cfg[0] = trial["standingAngle"]
            failed = neural.lib.hill_neural_integrate(np.asfortranarray(model["parameters"]), a, mass, wait,
                len(a), np.asfortranarray([trial["q"][0], trial["v0"]]), trial["dt"], 10,
                np.asfortranarray(cfg), np.asfortranarray(neural.parameter_rows(np.tile(model["neuralParameters"], (2, 1)))[0]), output)
            assert not failed
            prefix_error = np.max(np.abs(output[:len(expected)] - expected))
            assert prefix_error < 1e-9, "Appending future waiting changed the original past."
            tail = output[-round(10 / trial["dt"]):]
            angle_error = np.max(np.abs(tail[:, 0] - trial["standingAngle"]))
            speed = np.max(np.abs(tail[:, 1]))
            expected_torque = 9.81 * model["parameters"][18] * model["parameters"][19] * np.cos(trial["standingAngle"])
            settled_torque_error = np.max(np.abs(tail[:, 2] - expected_torque))
            soft_peak = np.max(np.abs(boundary_torque(model["parameters"], output[:, 0])))
            assert angle_error < .001 and speed < .001 and settled_torque_error < .1
            assert np.min(output[:, 0]) >= angle_policy["sequenceLower"]
            assert np.max(output[:, 0]) <= angle_policy["sequenceUpper"] and soft_peak <= 50
            waiting.append({"group": group, "source": source + 1, "appendedWaitSeconds": 30,
                            "checkedTailSeconds": 10, "totalDurationSeconds": (len(a) - 1) * trial["dt"],
                            "maxSettledAngleErrorRad": angle_error, "maxSettledSpeedRadPerSecond": speed,
                            "maxSettledTorqueErrorNm": settled_torque_error, "softLimitPeakNm": soft_peak,
                            "minimumAngleRad": np.min(output[:, 0]), "maximumAngleRad": np.max(output[:, 0]),
                            "abdominalFullSequencePeakNm": np.max(output[:, 4]), "prefixMaxDifference": prefix_error})
    write_table(snapshots / "matlab_replay_checks.csv", numerical)
    write_table(snapshots / "extended_wait_checks.csv", waiting)
    status = {"passed": True, "portableModelSHA256": sha(cache), "matlabReplaySHA256": sha(args.matlab_replay),
              "fixedParameterConditionsPerSource": 9, "sources": 2, "continuousRecordsPerSource": 3,
              "originalReferencesExactlyPreserved": True, "appendedWaitSeconds": 30,
              "anglePolicy": angle_policy,
              "parameterCountsPerSource": [32 + len(m["neuralParameters"]) for m in models],
              "sharedStandingController": np.asarray(package["controller"]).tolist(),
              "extensionPolicy": "Hold the final waiting input, set load zero, continue all original states",
              "scope": "Numerical reproducibility and prolonged standing; not independent physiological validation"}
    (snapshots / "replay_checks.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
    metadata["files"] = [{"file": p.name, "sha256": sha(p)} for p in sorted(snapshots.iterdir())
                         if p.suffix in (".png", ".pdf", ".svg", ".csv")]
    metadata_file.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print("PASS: fixed parameters across all nine conditions, exact references, complete C/MATLAB replay, +30 s waits.")
    print("Maximum replay angle difference:", max(r["angleMaxDifferenceRad"] for r in numerical))
    print("Maximum replay torque difference:", max(r["netTorqueMaxDifferenceNm"] for r in numerical))
    print("Maximum extended-wait angle error:", max(r["maxSettledAngleErrorRad"] for r in waiting))
    print("Checks:", snapshots)


if __name__ == "__main__":
    main()
