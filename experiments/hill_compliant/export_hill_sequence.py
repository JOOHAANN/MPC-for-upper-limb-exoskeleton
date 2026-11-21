"""Export audited continuous-sequence figures and tables for the repository.

Run after run_hill_sequence and test_hill_sequence. References are checked
against the original full MAT records. Abdominal torque is independently
reconstructed from tendon strain and moment arm before snapshots are copied.
"""
import csv
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results/hill_sequence"
OUT = ROOT / "docs/hill_sequence_evaluation"
GROUPS = ("5_35_10_55_15_45", "5_45_10_35_15_55", "5_55_10_45_15_35")


def read_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(name, records):
    with (OUT / name).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def values(records, key):
    return np.array([float(row[key]) for row in records])


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    result = loadmat(RESULTS / "hill_sequence.mat", simplify_cells=True)["result"]
    assert len(result["groups"]) == 3
    audit, protocol, force_audit, torque_audit = [], [], [], []
    for group, (name, run) in enumerate(zip(GROUPS, result["groups"]), 1):
        assert isinstance(run, dict), "Export requires all three groups."
        trial = run["sequence"]
        comparison = read_csv(RESULTS / f"abdominal_torque_comparison_group{group}.csv")
        assert np.allclose(values(comparison, "time"), trial["t"], atol=1e-12, rtol=0)
        original_path = ROOT / f"data/processed/training_map__muscle_tendon__{name}_done__build_data/Combination_raw.mat"
        original = loadmat(original_path)["Combination"][::10]
        assert np.array_equal(trial["q"], original[:, 2])
        assert np.array_equal(trial["torque"], original[:, 0])
        assert np.isclose(float(trial["standingAngle"]), 1.5, atol=1e-12, rtol=0)
        audit.append({
            "group": group, "originalDataset": name, "samples": len(original),
            "durationSeconds": float(trial["t"][-1]), "standingAngleRad": float(trial["standingAngle"]),
            "maximumReferenceAngleDifferenceRad": 0.,
            "maximumReferenceTorqueDifferenceNm": 0.,
            "originalFilledWaitingTorqueNm": 20.,
            "filledReferenceSamples": int(np.sum(trial["referenceIsFilled"])),
            "originalMatSHA256": digest(original_path),
        })
        for segment in trial["segments"]:
            protocol.append({"group": group, **segment})
        for source, prediction in enumerate(run["predictions"], 1):
            par = result["models"][source - 1]["parameters"]
            assert not prediction["failed"]
            assert np.all(prediction["externalBalanceTorque"] == 0)
            radius = par[21]
            lo = radius / par[10]
            width = np.sin(par[27])
            slack = par[23] - lo * np.sqrt(par[11]**2 - width**2)
            lmt = par[23] + radius * (prediction["q"] - 1.4)
            lt = lmt - lo * np.sqrt(np.maximum(0, prediction["frontLengthNormalized"]**2 - width**2))
            tendon_force = (par[8] / radius) * np.maximum(0, (lt - slack) / (par[25] * slack))**2
            force = prediction["muscleForces"]
            prefix = "experimental" if source == 1 else "map"
            for column, key in (("AbdominalTorqueNm", "frontTorque"),
                                ("ActiveCETorqueNm", "frontActive"),
                                ("PassivePETorqueNm", "frontPassive"),
                                ("FiberDampingTorqueNm", "frontFiberDamping")):
                assert np.allclose(values(comparison, prefix + column), prediction[key], atol=1e-9, rtol=0)
            assert np.allclose(values(comparison, prefix + "AbdominalJointContributionNm"),
                               -prediction["frontTorque"], atol=1e-9, rtol=0)
            assert np.allclose(values(comparison, prefix + "HillActivation"),
                               prediction["effectiveActivation"][:, 1], atol=1e-12, rtol=0)
            base = values(comparison, prefix + "BaseFrontInput")
            assert np.all((base >= 0) & (base <= 1))
            if source == 1:
                assert np.allclose(base, trial["front"], atol=1e-12, rtol=0)
            geometry_residual = np.max(np.abs(tendon_force - force["frontTendonN"]))
            component_residual = np.max(np.abs(force["frontTendonN"] - force["frontActiveN"]
                                               - force["frontPassiveN"] - force["frontFiberDampingN"]))
            trace = read_csv(RESULTS / f"group{group}_source{source}.csv")
            assert np.max(np.abs(values(trace, "frontTendonForceN") - tendon_force)) < 1e-5
            assert geometry_residual < 1e-5 and component_residual < 1e-5
            force_audit.append({"group": group, "source": source,
                                "tendonGeometryResidualN": float(geometry_residual),
                                "forceComponentResidualN": float(component_residual)})
            torque_residual = np.max(np.abs(radius * tendon_force - prediction["frontTorque"]))
            components = prediction["frontActive"] + prediction["frontPassive"] + prediction["frontFiberDamping"]
            component_torque_residual = np.max(np.abs(prediction["frontTorque"] - components))
            joint_residual = np.max(np.abs(prediction["torque"] - prediction["backTorque"] + prediction["frontTorque"]))
            assert torque_residual < 1e-7 and component_torque_residual < 1e-7 and joint_residual < 1e-9
            assert np.allclose(values(trace, "abdominalJointContributionNm"),
                               -prediction["frontTorque"], atol=1e-9, rtol=0)
            torque_audit.append({"group": group, "source": source,
                                 "tendonMomentArmResidualNm": float(torque_residual),
                                 "torqueComponentResidualNm": float(component_torque_residual),
                                 "signedJointBalanceResidualNm": float(joint_residual)})
    convergence = read_csv(RESULTS / "convergence.csv")
    stability = read_csv(RESULTS / "stability.csv")
    movements = read_csv(RESULTS / "movement_metrics.csv")
    waiting = read_csv(RESULTS / "waiting_metrics.csv")
    forces = read_csv(RESULTS / "muscle_force_metrics.csv")
    torques = read_csv(RESULTS / "muscle_torque_metrics.csv")
    assert len(convergence) == len(stability) == len(waiting) == len(forces) == len(torques) == 6
    assert len(movements) == 18 and np.all(values(movements, "Wmatched") == 1)
    assert max(values(convergence, "angleDifferenceRad")) < .01
    assert max(values(convergence, "torqueDifferenceNm")) < 3
    assert max(values(stability, "finalAngleErrorRad")) < .001
    assert max(values(waiting, "maxSettledAngleErrorRad")) < .01
    OUT.mkdir(parents=True, exist_ok=True)
    write_csv("reference_audit.csv", audit)
    write_csv("protocol.csv", protocol)
    write_csv("force_audit.csv", force_audit)
    write_csv("torque_audit.csv", torque_audit)
    summary = []
    for source in ("Experimental activation", "MAP activation"):
        selected = [row for row in movements if row["source"] == source]
        waits = [row for row in waiting if row["source"] == source]
        muscle = [row for row in forces if row["source"] == source]
        muscle_torque = [row for row in torques if row["source"] == source]
        assert len(selected) == 9
        summary.append({
            "source": source, "movementConditions": len(selected),
            "meanMovementThetaRMSEDegrees": float(values(selected, "thetaRMSEDegrees").mean()),
            "meanMovementTorqueRMSENm": float(values(selected, "torqueRMSENm").mean()),
            "twoValleyMatches": int(values(selected, "Wmatched").sum()),
            "maximumSettledWaitingErrorRad": float(values(waits, "maxSettledAngleErrorRad").max()),
            "minimumPredictedAngleRad": float(values(waits, "minAngleRad").min()),
            "maximumPredictedAngleRad": float(values(waits, "maxAngleRad").max()),
            "abdominalTendonPeakN": float(values(muscle, "frontTendonPeakN").max()),
            "backTendonPeakN": float(values(muscle, "backTendonPeakN").max()),
            "abdominalFlexionTorquePeakNm": float(values(muscle_torque, "abdominalFlexionTorquePeakNm").max()),
            "backHillTorquePeakNm": float(values(muscle_torque, "backHillTorquePeakNm").max()),
            "maximumSoftLimitTorqueNm": float(values(waits, "softLimitPeakNm").max()),
            "externalBalanceTorqueNm": 0.,
        })
    write_csv("summary.csv", summary)
    for filename in ("movement_metrics.csv", "waiting_metrics.csv", "muscle_force_metrics.csv", "muscle_torque_metrics.csv",
                     "convergence.csv", "stability.csv"):
        shutil.copy2(RESULTS / filename, OUT / filename)
    files = []
    for group in (1, 2, 3):
        csv_name = f"abdominal_torque_comparison_group{group}.csv"
        shutil.copy2(RESULTS / csv_name, OUT / csv_name)
        files.append({"file": csv_name, "sha256": digest(RESULTS / csv_name)})
        for stem in ("torque_theta", "abdominal_torque", "abdominal_torque_comparison", "standing_balance", "sequence"):
            for suffix in (".png", ".pdf"):
                name = f"{stem}_group{group}{suffix}"
                path = RESULTS / name
                shutil.copy2(path, OUT / name)
                files.append({"file": name, "sha256": digest(path)})
    metadata = {
        "resultMatSHA256": digest(RESULTS / "hill_sequence.mat"),
        "calibratedModelsSHA256": digest(ROOT / "data/processed/hill_compliant/calibrated_models.mat"),
        "controller": np.asarray(result["options"]["controller"]).tolist(),
        "dtSeconds": .01, "substepsPerSample": int(result["options"]["steps"]),
        "angleReferencePolicy": "Exact original full Scale records, including their 1.5-rad waiting seams",
        "waitingTorquePolicy": "Original filled 20-Nm value displayed; excluded from movement fit metrics",
        "forcePolicy": "Tendon tension in N; CE, PE and fiber damping include pennation projection",
        "abdominalTorquePolicy": "Primary abdominal plots use Nm from the Hill solver; positive flexion magnitude, negative signed joint contribution",
        "evaluationPolicy": "Continuous calibrated replay, not an additional independent holdout",
        "files": files,
    }
    (OUT / "export_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    for record in summary:
        print(record)
    print("PASS: original references, abdominal torque units, Hill decomposition and signed joint contribution.")
    print("Exported figures and tables:", OUT)


if __name__ == "__main__":
    main()
