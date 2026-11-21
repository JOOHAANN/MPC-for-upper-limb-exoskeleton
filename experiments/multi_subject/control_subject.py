"""Tune one causal MPC policy per subject, then replay every available case.

Controller selection uses training cycles only. The Hill models and MAP are
frozen throughout tuning and testing. Measured angle/torque never enter the
online MPC; a separate unassisted Hill state supplies its motion target.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.io import loadmat, savemat
from scipy.signal import find_peaks

from audit_subjects import ROOT, OUT
from prepare_subjects import combine, structs
from fit_subject import many, record, rollout

sys.path.insert(0, str(ROOT / "experiments/legacy_mpc"))
from run_hill_mpc import Plant, simulate_case, prediction_dict, LABELS, OUTPUT_NAMES


def options(weight=128.0, margin=35.0, velocity=100.0, recovery=12.0):
    return {
        "points": 16,
        "blocks": 4,
        "predictionDt": 0.025,
        "predictionSubsteps": 10,
        "controlPeriod": 0.02,
        "actualSubsteps": 10,
        "rateNmPerSec": 200.0,
        "capNm": 80.0,
        "angleWeight": 10.0,
        "backTorqueWeight": weight,
        "backTorqueMarginNm": margin,
        "inputForecast": "hold",
        "inputForecastTime": 0.10,
        "inputSlopeWindow": 0.05,
        "inputSlopeLimit": 3.0,
        "velocityWeight": velocity,
        "recoveryBackTargetNm": recovery,
        "recoveryBackWeight": 512.0,
        "recoveryVelocityWeight": 10.0,
        "earlyRecovery": True,
        "recoveryContactToleranceKg": 0.5,
        "recoveryIntentToleranceRad": 0.10,
        "recoveryAngleToleranceRad": 0.025,
        "recoveryVelocityToleranceRadPerSec": 0.04,
    }


def metrics(actual, nominal, seq, flags, source):
    rows = []
    for seg in many(seq["segments"]):
        ids = slice(int(seg["startIndex"]) - 1, int(seg["endIndex"]))
        limit = float(np.interp(seg["loadKg"], [5, 10, 15], [100, 140, 180]))
        excess = np.maximum(0, actual[ids, 3] - limit)
        q = actual[ids, 0]
        q0 = nominal[ids, 0]
        ref = seq["q"][ids]
        rows.append(
            {
                "recordId": seg["recordId"],
                "source": LABELS[source],
                "loadKg": float(seg["loadKg"]),
                "speedCode": int(seg["speedCode"]),
                "split": seg["split"],
                "angleReferenceUsable": bool(seg["angleUsable"]),
                "backThresholdNm": limit,
                "backPeakNm": float(actual[ids, 3].max()),
                "nominalBackPeakNm": float(nominal[ids, 3].max()),
                "backMeanNm": float(actual[ids, 3].mean()),
                "abdomenMeanNm": float(actual[ids, 4].mean()),
                "abdomenPeakNm": float(actual[ids, 4].max()),
                "abdomenBackRatio": float(
                    actual[ids, 4].sum() / max(actual[ids, 3].sum(), 1e-9)
                ),
                "excessPeakNm": float(excess.max()),
                "excessMeanNm": float(excess.mean()),
                "excessOver1NmSeconds": float(np.sum(excess > 1) * seq["dt"]),
                "thetaRMSEDegrees": float(
                    np.sqrt(np.mean((q - ref) ** 2)) * 180 / np.pi
                ),
                "shapeCorrelation": (
                    float(np.corrcoef(q, ref)[0, 1]) if np.std(q) > 1e-8 else 0.0
                ),
                "nominalAngleDeviationRad": float(np.sqrt(np.mean((q - q0) ** 2))),
                "nominalValleys": len(find_peaks(-q0, prominence=0.12, distance=50)[0]),
                "predictedValleys": len(
                    find_peaks(-q, prominence=0.12, distance=50)[0]
                ),
                "referenceValleys": len(
                    find_peaks(-ref, prominence=0.12, distance=50)[0]
                ),
                "minAngleRad": float(q.min()),
                "maxAngleRad": float(q.max()),
                "assistancePeakNm": float(actual[ids, 24].max()),
                "boundaryPeakNm": float(np.abs(actual[ids, 25]).max()),
            }
        )
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("subject")
    parser.add_argument("--skip-tuning", action="store_true")
    parser.add_argument("--reuse-tuning", action="store_true")
    parser.add_argument(
        "--constraint-controller",
        action="store_true",
        help="Use training-calibrated prediction guards and explicit minimum-slack torque priority.",
    )
    parser.add_argument(
        "--refine-constraint-horizon",
        action="store_true",
        help="Keep prior full-horizon training candidate and test 100/200 ms constraint horizons on training data.",
    )
    parser.add_argument(
        "--training-only",
        action="store_true",
        help="Save training options/candidates without replaying holdouts.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Reuse saved cases only when parameters, options, inputs and references match exactly.",
    )
    parser.add_argument(
        "--angle-governor",
        type=float,
        default=0.0,
        help="Predictive extension threshold, rad; zero disables the optional assistance governor.",
    )
    parser.add_argument(
        "--frozen-options",
        type=Path,
        help="Replay published MPC options without re-tuning.",
    )
    parser.add_argument(
        "--assistance-cap",
        type=float,
        default=80.0,
        help="Personalized assistance cap; hardware maximum remains 80 N m.",
    )
    args = parser.parse_args()
    out = OUT / args.subject
    start = time.monotonic()
    data = loadmat(out / "prepared.mat", simplify_cells=True)
    if args.constraint_controller and data.get("qualityGateVersion") not in (
        "clock-smoothing-gate-v2",
        "excitation-ik-three-loads-v1",
    ):
        raise RuntimeError("Prediction-guard calibration requires qualified data.")
    cache = loadmat(out / "neural_models.mat", simplify_cells=True)
    replay = loadmat(out / "hill_sequence.mat", simplify_cells=True)["result"]
    models = many(cache["models"])
    cfg = cache["controller"]
    # Keep the recovery target above the subject's unsupported static back
    # torque, so no continuing assistance is requested in a settled wait.
    steady = []
    for model in models:
        plant = Plant(model, cfg, 1.5)
        _, first = plant.init(np.zeros(2), 0.0, 1.0, np.array([1.5, 0.0]))
        steady.append(first[3])
    recovery = max(12.0, float(max(steady)) + 2.0)
    train = [c for c in many(data["cycles"]) if c["split"] == "train"]
    selected = []
    for load in (5, 10, 15):
        choices = [c for c in train if c["loadKg"] == load]
        if choices:
            selected.append(
                sorted(
                    choices,
                    key=lambda c: (not c["angleUsable"], c["angleTorqueCorrelation"]),
                )[0]
            )
    trainseq = combine(selected, args.subject, "controller_training_only")
    if args.constraint_controller:
        selected = train
        trainseq = combine(selected, args.subject, "controller_all_qualified_training")
    rec = record(trainseq, data["generator"])
    inputs = [
        rollout(
            models[s]["parameters"], models[s]["neuralParameters"], cfg, rec, s, 10
        )[2]
        for s in (0, 1)
    ]
    candidates = [
        options(64, 20, 100, recovery),
        options(128, 35, 100, recovery),
        options(256, 45, 60, recovery),
    ]
    guard_audit = None
    if args.constraint_controller:
        from calibrate_prediction_guard import estimate_guard

        guarded = options(64, 0, 100, recovery)
        guard, guard_audit = estimate_guard(data, models, cfg, guarded)
        guarded.update(torqueConstraintPriority=True, predictionErrorGuardNm=guard)
        candidates = [
            {**guarded, "torqueConstraintHorizonSeconds": h} for h in (0.1, 0.2, 0.4)
        ]
        if args.refine_constraint_horizon:
            candidates = candidates[:2]
        (out / "prediction_guard.json").write_text(
            json.dumps(guard_audit, indent=2) + "\n"
        )
        print(args.subject, "TRAINING PREDICTION GUARD", guard, flush=True)
    for candidate in candidates:
        candidate["capNm"] = args.assistance_cap
        if args.angle_governor > 0:
            candidate["angleGovernorUpperRad"] = args.angle_governor
    if args.skip_tuning:
        candidates = candidates[1:2]
    audited = []
    notified = [time.monotonic()]

    def progress(k, n):
        if time.monotonic() - notified[0] > 20:
            print(
                args.subject,
                "CONTROLLER",
                round(100 * k / n, 1),
                "%",
                "WALL",
                round(time.monotonic() - start, 1),
                flush=True,
            )
            notified[0] = time.monotonic()

    if args.reuse_tuning and (out / "mpc_candidates.json").exists():
        audited = json.loads((out / "mpc_candidates.json").read_text())
        candidates = []
    if args.refine_constraint_horizon:
        audited = json.loads((out / "mpc_candidates.json").read_text())
    if args.frozen_options is not None:
        frozen = json.loads(args.frozen_options.read_text(encoding="utf-8"))
        audited = [
            {
                "options": frozen["options"],
                "trainingScore": 0.0,
                "selection": "Frozen published controller; no tuning in this run",
            }
        ]
        candidates = []
    for opt in candidates:
        rows = []
        failures = 0
        for source in (0, 1):
            actual, nominal, flags, timings = simulate_case(
                models[source], cfg, trainseq, inputs[source], opt, progress
            )
            rows.extend(metrics(actual, nominal, trainseq, flags, source))
            failures += int(np.sum(flags < 0))
        score = float(
            np.mean(
                [
                    10 * (r["excessPeakNm"] / 20) ** 2
                    + 2 * (r["excessMeanNm"] / 5) ** 2
                    + (r["nominalAngleDeviationRad"] / 0.15) ** 2
                    + 3 * max(0, r["abdomenBackRatio"] - 0.25) ** 2
                    + 0.5 * abs(r["predictedValleys"] - r["nominalValleys"])
                    + 10 * max(0, r["maxAngleRad"] - 2.1) ** 2
                    for r in rows
                ]
            )
        )
        audited.append(
            {
                "options": opt,
                "trainingScore": score,
                "rejectedQPResults": failures,
                "cases": rows,
            }
        )
        print(
            args.subject,
            "CANDIDATE",
            opt["backTorqueWeight"],
            opt["backTorqueMarginNm"],
            "SCORE",
            score,
            flush=True,
        )
        (out / "mpc_candidates.json").write_text(json.dumps(audited, indent=2) + "\n")

    def selection(candidate):
        if not args.constraint_controller:
            return candidate["trainingScore"]
        peak = max(r["excessPeakNm"] for r in candidate["cases"])
        # First meet the training torque target (0.5 N m numerical tolerance),
        # then preserve motion. No held-out action enters horizon selection.
        valley_changes = sum(
            abs(r["predictedValleys"] - r["nominalValleys"]) for r in candidate["cases"]
        )
        return (
            peak > 0.5,
            peak if peak > 0.5 else 0.0,
            valley_changes,
            candidate["trainingScore"],
        )

    chosen = (
        audited[0]["options"]
        if args.frozen_options is not None
        else min(audited, key=selection)["options"]
    )
    if args.angle_governor > 0:
        chosen = {**chosen, "angleGovernorUpperRad": args.angle_governor}
    if args.assistance_cap != 80.0:
        chosen = {**chosen, "capNm": args.assistance_cap}
    if args.training_only:
        (out / "mpc_selected_training_options.json").write_text(
            json.dumps({"options": chosen, "trainingOnly": True}, indent=2) + "\n"
        )
        print("TRAINING SELECTION COMPLETE", args.subject, chosen, flush=True)
        return
    existing = None
    if args.resume and (out / "hill_mpc.mat").exists():
        candidate = loadmat(out / "hill_mpc.mat", simplify_cells=True)["result"]
        same_options = set(candidate["options"]) == set(chosen) and all(
            np.array_equal(
                np.asarray(candidate["options"][key]), np.asarray(chosen[key])
            )
            for key in chosen
        )
        if same_options and all(
            np.array_equal(a["parameters"], b["parameters"])
            and np.array_equal(a["neuralParameters"], b["neuralParameters"])
            for a, b in zip(many(candidate["models"]), models)
        ):
            existing = {int(g["group"]): g for g in many(candidate["groups"])}
    groups = []
    allrows = []
    checks = []
    for g in many(replay["groups"]):
        seq = g["sequence"]
        predictions = []
        nominals = []
        for source in (0, 1):
            saved = existing.get(int(g["group"])) if existing else None
            reusable = (
                saved is not None
                and np.array_equal(saved["inputs"][source], g["inputs"][source])
                and all(
                    np.array_equal(saved["sequence"][k], seq[k], equal_nan=True)
                    for k in ("q", "torque", "m", "waiting")
                )
            )
            if reusable:
                actual = np.column_stack(
                    [saved["predictions"][source][k] for k in OUTPUT_NAMES]
                )
                nominal = np.column_stack(
                    [saved["nominal"][source][k] for k in OUTPUT_NAMES]
                )
                flags = saved["predictions"][source]["solverExitFlag"]
                timings = saved["predictions"][source]["solveSeconds"]
            else:
                actual, nominal, flags, timings = simulate_case(
                    models[source], cfg, seq, g["inputs"][source], chosen, progress
                )
            baseline = g["predictions"][source]
            qdiff = float(np.max(np.abs(nominal[:, 0] - baseline["q"])))
            tdiff = float(np.max(np.abs(nominal[:, 2] - baseline["torque"])))
            assert qdiff < 1e-7 and tdiff < 1e-4, (qdiff, tdiff)
            assert actual[:, 24].min() >= -1e-8 and actual[:, 24].max() <= 80 + 1e-8
            assert np.max(np.abs(np.diff(actual[:, 24]))) <= 200 * seq["dt"] + 1e-8
            pred = prediction_dict(actual, models[source]["parameters"])
            pred["solverExitFlag"] = flags
            pred["solveSeconds"] = timings
            predictions.append(pred)
            nominals.append(prediction_dict(nominal, models[source]["parameters"]))
            rows = metrics(actual, nominal, seq, flags, source)
            allrows.extend(
                [
                    {
                        **r,
                        "group": int(g["group"]),
                        "displayOnly": seq["name"] == "continuous_available_loads",
                    }
                    for r in rows
                ]
            )
            checks.append(
                {
                    "group": int(g["group"]),
                    "source": LABELS[source],
                    "qParityRad": qdiff,
                    "torqueParityNm": tdiff,
                    "qpRejected": int(np.sum(flags < 0)),
                    "maxAngleRad": float(actual[:, 0].max()),
                    "minAngleRad": float(actual[:, 0].min()),
                    "momentumResidualNm": float(np.abs(actual[:, 16]).max()),
                    "parameterSHA256": hashlib.sha256(
                        np.r_[
                            models[source]["parameters"],
                            models[source]["neuralParameters"],
                        ].tobytes()
                    ).hexdigest(),
                }
            )
            print(args.subject, "GROUP", g["group"], LABELS[source], "DONE", flush=True)
        groups.append(
            {
                "group": g["group"],
                "sequence": seq,
                "inputs": g["inputs"],
                "predictions": predictions,
                "nominal": nominals,
            }
        )
        savemat(
            out / "hill_mpc.mat",
            {
                "result": {
                    "groups": structs(groups),
                    "models": structs(models),
                    "controller": cfg,
                    "options": chosen,
                }
            },
            do_compression=True,
            long_field_names=True,
        )
    metadata = {
        "subject": args.subject,
        "options": chosen,
        "metrics": allrows,
        "checks": checks,
        "controllerTrainingRecordIds": [r["id"] for r in selected],
        "controllerUsesTestRecords": False,
        "modelSHA256": hashlib.sha256(
            (out / "neural_models.mat").read_bytes()
        ).hexdigest(),
        "wallSeconds": time.monotonic() - start,
        "sameControllerAllCases": True,
        "independentAssistanceExperimentFitted": False,
        "torqueConstraintPriority": bool(chosen.get("torqueConstraintPriority", False)),
        "predictionGuardCalibration": guard_audit or {},
        "qualityGateVersion": data.get(
            "qualityGateVersion", "historical_unchecked_replay"
        ),
        "trainingSelectionPriority": (
            "torque feasibility, no new motion valleys, then tracking score"
            if args.constraint_controller
            else "tracking score"
        ),
    }
    (out / "mpc_evaluation.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print("COMPLETED", args.subject, "SECONDS", time.monotonic() - start, flush=True)


if __name__ == "__main__":
    main()
