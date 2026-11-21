"""Tune subject-wide NMPC with limits at 75% of each unassisted back peak.

Limits are declared offline calibration constants. The online controller sees
the current state, past activation and those constants, not future measured
EMG, reference angle or reference torque. Saved simulations are never clipped.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.io import loadmat, savemat
from scipy.signal import find_peaks

from audit_subjects import ROOT
from fit_subject import many, record, rollout
from prepare_subjects import structs
from run_hill_mpc import (
    Plant,
    simulate_case,
    prediction_dict,
    LABELS,
    configure_map_preview,
)
from fit_hill_compliant import pair_inputs, memory

BASE = ROOT / "results/nine_case_75pct"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_subject(subject):
    if subject == "240528Pilot":
        model_path = ROOT / "data/processed/hill_compliant/neural_models.mat"
        replay_path = ROOT / "results/hill_neural_sequence/hill_sequence.mat"
    else:
        model_path = BASE / subject / "neural_models.mat"
        replay_path = BASE / subject / "hill_sequence.mat"
    cache = loadmat(model_path, simplify_cells=True)
    replay = loadmat(replay_path, simplify_cells=True)["result"]
    models = many(cache["models"])
    assert len(models) == 2
    assert len({m.get("intentDecoder", "subtractive") for m in models}) == 1
    for model, used in zip(models, many(replay["models"])):
        assert np.array_equal(model["parameters"], used["parameters"])
        assert np.array_equal(model["neuralParameters"], used["neuralParameters"])
        assert model.get("intentDecoder", "subtractive") == used.get(
            "intentDecoder", "subtractive"
        )
    groups = many(replay["groups"])
    # Exp files contain nine single-action evaluations followed by optional
    # continuous demonstrations. Count each recorded condition exactly once.
    if subject != "240528Pilot":
        groups = groups[:9]
    attach_planned_inputs(groups, models)
    return cache, models, groups, model_path, replay_path


def attach_planned_inputs(groups, models):
    """Build generated activation previews with each source's fixed front end."""
    for group in groups:
        planned = []
        for source, model in enumerate(models):
            if source == 1:
                planned.append(np.asarray(group["inputs"])[1])
            else:
                seq = group["sequence"]
                base = pair_inputs(seq, 1, model["generator"])
                difference = memory(base[:, 0], 0.1) - memory(base[:, 0], 0.35)
                p = model["parameters"]
                inferred = (
                    p[15] * base[:, 0]
                    + p[16] * np.maximum(0, difference)
                    + p[17] * np.maximum(0, -difference)
                )
                front = -np.expm1(-np.maximum(0, p[14] * base[:, 1] + inferred))
                planned.append(np.c_[base[:, 0], front])
        group["plannedMAPInputs"] = np.asarray(planned)


def calibrate_guard(models, controller, groups, options):
    """Estimate unknown-input forecast errors from calibration actions only."""
    errors, normalized_errors = [], []
    by_source = []
    for source, model in enumerate(models):
        source_errors, source_normalized = [], []
        source_options = dict(options)
        if options["inputForecast"] == "hybrid":
            source_options["inputForecast"] = "hold" if source == 0 else "map"
        plant = Plant(model, controller, 1.5)
        for group in groups:
            seq = dict(group["sequence"])
            seq["plannedMAPInputs"] = group["plannedMAPInputs"][source]
            seq["activationSourceIsGenerated"] = bool(source == 1)
            inputs = np.asarray(group["inputs"])[source]
            state, first = plant.init(
                inputs[0],
                seq["m"][0],
                seq["waiting"][0],
                np.array([seq["q"][0], seq["v0"]]),
            )
            states, truth = [state.copy()], [first]
            for k in range(1, len(inputs)):
                truth.append(
                    plant.step(
                        state,
                        inputs[k],
                        seq["m"][k],
                        seq["waiting"][k],
                        0.0,
                        seq["dt"],
                        options["actualSubsteps"],
                    )
                )
                states.append(state.copy())
            truth = np.asarray(truth)
            thresholds, _ = threshold_schedule(seq, many(group["predictions"])[source])
            for k in range(0, len(inputs) - 42, 10):
                if seq["waiting"][k]:
                    continue
                h = max(0, k - round(options["inputSlopeWindow"] / seq["dt"]))
                plant.input_slope[:] = np.clip(
                    (inputs[k] - inputs[h]) / max((k - h) * seq["dt"], 1e-12),
                    -options["inputSlopeLimit"],
                    options["inputSlopeLimit"],
                )
                if source_options.get("inputForecast") in ("map", "arx"):
                    configure_map_preview(plant, seq, inputs, k, source_options)
                forecast = plant.preview(
                    states[k], np.zeros(options["blocks"]), source_options
                )
                if forecast is not None:
                    future = (
                        k
                        + np.arange(1, options["points"] + 1)
                        * options["predictionDt"]
                        / seq["dt"]
                    )
                    error = (
                        np.interp(future, np.arange(len(truth)), truth[:, 3])
                        - forecast[:, 3]
                    )
                    errors.append(error)
                    normalized_errors.append(error / (thresholds[k] / 0.75))
                    source_errors.append(error)
                    source_normalized.append(error / (thresholds[k] / 0.75))
        by_source.append(
            dict(
                inputForecast=source_options["inputForecast"],
                guardNm=(
                    2.0 + np.quantile(np.maximum(source_errors, 0.0), 0.99, axis=0)
                ).tolist(),
                relativeGuard=np.quantile(
                    np.maximum(source_normalized, 0.0), 0.99, axis=0
                ).tolist(),
            )
        )
    if not errors:
        raise RuntimeError("No calibration forecast windows")
    errors = np.asarray(errors)
    return (
        2.0 + np.quantile(np.maximum(errors, 0.0), 0.99, axis=0),
        np.quantile(np.maximum(normalized_errors, 0.0), 0.99, axis=0),
        by_source,
    )


def threshold_schedule(seq, prediction):
    limits, metadata = [], []
    schedule = np.full(len(seq["q"]), 100.0)
    for segment in many(seq["segments"]):
        ids = np.arange(int(segment["startIndex"]) - 1, int(segment["endIndex"]))
        peak = float(np.max(np.asarray(prediction["backTorque"])[ids]))
        limit = 0.75 * peak
        if not np.isfinite(limit) or limit <= 0:
            raise ValueError("Invalid unassisted back torque peak")
        schedule[ids] = limit
        limits.append(limit)
        metadata.append(
            dict(
                recordId=segment.get(
                    "recordId", f"{segment['loadKg']:g}kg_{segment['speedCode']:g}"
                ),
                loadKg=float(segment["loadKg"]),
                speedCode=int(segment["speedCode"]),
                unassistedBackPeakNm=peak,
                backThresholdNm=limit,
            )
        )
    schedule[np.asarray(seq["waiting"], dtype=bool)] = min(limits)
    return schedule, metadata


def evaluate(actual, nominal, seq, source, metadata):
    rows = []
    for segment, meta in zip(many(seq["segments"]), metadata):
        ids = np.arange(int(segment["startIndex"]) - 1, int(segment["endIndex"]))
        limit = meta["backThresholdNm"]
        q, reference = actual[ids, 0], np.asarray(seq["q"])[ids]
        excess = np.maximum(actual[ids, 3] - limit, 0)
        peak = float(actual[ids, 3].max())
        demand = nominal[ids, 3] > limit
        valleys = lambda y: len(find_peaks(-y, prominence=0.12, distance=50)[0])
        rows.append(
            dict(
                **meta,
                source=LABELS[source],
                assistedBackPeakNm=peak,
                peakMarginNm=limit - peak,
                assistedToUnassistedPeakRatio=peak / meta["unassistedBackPeakNm"],
                excessPeakNm=float(excess.max()),
                excessOver1NmSeconds=float(np.sum(excess > 1) * seq["dt"]),
                demandMeanShortfallNm=(
                    float(np.mean(np.maximum(limit - actual[ids[demand], 3], 0)))
                    if demand.any()
                    else 0.0
                ),
                thetaShapeCorrelation=(
                    float(np.corrcoef(q, reference)[0, 1]) if np.std(q) > 1e-9 else 0.0
                ),
                thetaRMSEDegrees=float(
                    np.sqrt(np.mean((q - reference) ** 2)) * 180 / np.pi
                ),
                nominalValleys=valleys(nominal[ids, 0]),
                predictedValleys=valleys(q),
                referenceValleys=valleys(reference),
                nominalAngleDeviationRad=float(
                    np.sqrt(np.mean((q - nominal[ids, 0]) ** 2))
                ),
                minimumAngleRad=float(q.min()),
                maximumAngleRad=float(q.max()),
                abdomenPeakNm=float(actual[ids, 4].max()),
                abdomenBackMeanRatio=float(
                    actual[ids, 4].sum() / max(actual[ids, 3].sum(), 1e-9)
                ),
                assistancePeakNm=float(actual[ids, 24].max()),
                assistanceImpulseNms=float(actual[ids, 24].sum() * seq["dt"]),
                boundaryPeakNm=float(np.abs(actual[ids, 25]).max()),
                peakWithinFiveNmBelow=bool(-0.05 <= limit - peak <= 5.0),
            )
        )
    return rows


def main():
    global BASE
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("subject")
    parser.add_argument("--results-root", type=Path, default=BASE)
    parser.add_argument("--candidate", default="v1")
    parser.add_argument("--guard-scale", type=float, default=0.25)
    parser.add_argument("--guard-floor", type=float, default=2.0)
    parser.add_argument("--guard-cap", type=float, default=0.0)
    parser.add_argument("--under-weight", type=float, default=12.0)
    parser.add_argument("--effort-weight", type=float, default=3.0)
    parser.add_argument("--headroom", type=float, default=2.0)
    parser.add_argument(
        "--safety-margin",
        type=float,
        default=0.0,
        help="Subject-wide numerical/prediction reserve below the displayed limit, in Nm.",
    )
    parser.add_argument("--horizon", type=float, default=0.2)
    parser.add_argument("--angle-weight", type=float, default=10.0)
    parser.add_argument("--velocity-weight", type=float, default=100.0)
    parser.add_argument(
        "--forecast", choices=("hold", "trend", "map", "hybrid", "arx"), default="hold"
    )
    parser.add_argument("--innovation-decay", type=float, default=0.25)
    parser.add_argument("--phase-observer", action="store_true")
    parser.add_argument("--forecast-time", type=float, default=0.10)
    parser.add_argument("--relative-guard", action="store_true")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--groups", type=int, nargs="*")
    args = parser.parse_args()
    BASE = args.results_root
    cache, models, groups, model_path, replay_path = load_subject(args.subject)
    controller = cache["controller"]
    out = BASE / args.subject / "nmpc" / args.candidate
    out.mkdir(parents=True, exist_ok=True)
    options = dict(
        points=16,
        blocks=4,
        predictionDt=0.025,
        predictionSubsteps=10,
        controlPeriod=0.02,
        actualSubsteps=10,
        rateNmPerSec=200.0,
        capNm=80.0,
        angleWeight=args.angle_weight,
        backTorqueWeight=64.0,
        backTorqueMarginNm=args.safety_margin,
        inputForecast=args.forecast,
        inputForecastTime=args.forecast_time,
        inputSlopeWindow=0.05,
        inputSlopeLimit=3.0,
        velocityWeight=args.velocity_weight,
        recoveryBackTargetNm=12.0,
        recoveryBackWeight=512.0,
        recoveryVelocityWeight=10.0,
        earlyRecovery=True,
        recoveryContactToleranceKg=0.5,
        recoveryIntentToleranceRad=0.10,
        recoveryAngleToleranceRad=0.025,
        recoveryVelocityToleranceRadPerSec=0.04,
        torqueConstraintPriority=True,
        torqueConstraintHorizonSeconds=args.horizon,
        predictionErrorGuardScale=args.guard_scale,
        predictionErrorGuardFloorNm=args.guard_floor,
        backTorqueUnderWeight=args.under_weight,
        trackingHeadroomNm=args.headroom,
        assistanceEffortWeight=args.effort_weight,
        angleGovernorUpperRad=1.90,
    )
    options["mapInnovationDecaySeconds"] = args.innovation_decay
    options["mapPhaseObserver"] = args.phase_observer
    options["predictionErrorGuardCapNm"] = args.guard_cap
    # Keep the settled standing equilibrium unsupported, as in the benchmark.
    steady = []
    for model in models:
        plant = Plant(model, controller, 1.5)
        _, first = plant.init(np.zeros(2), 0.0, 1.0, np.array([1.5, 0.0]))
        steady.append(first[3])
    options["recoveryBackTargetNm"] = max(12.0, max(steady) + 2.0)
    signature = dict(
        modelSHA256=digest(model_path),
        replaySHA256=digest(replay_path),
        forecastOptions=dict(
            mode=args.forecast,
            innovationDecay=args.innovation_decay,
            phaseObserver=args.phase_observer,
            forecastTime=args.forecast_time,
        ),
        predictionCalibrationVersion="source-specific-uncertainty-v2",
    )
    if args.forecast == "arx":
        from activation_forecast import fit_predictor

        options["activationResidualPredictor"] = fit_predictor(
            groups,
            options,
            signature,
            BASE / args.subject / "nmpc/activation_residual_predictor.json",
        )
    guard_path = (
        BASE
        / args.subject
        / "nmpc"
        / ("prediction_guard_" + args.forecast + str(args.forecast_time) + ".json")
    )
    saved = json.loads(guard_path.read_text()) if guard_path.exists() else None
    if saved is None or saved["signature"] != signature:
        guard, relative, by_source = calibrate_guard(
            models, controller, groups, options
        )
        guard_path.write_text(
            json.dumps(
                dict(
                    signature=signature,
                    guardNm=guard.tolist(),
                    relativeGuard=relative.tolist(),
                    bySource=by_source,
                    calibrationOnly=True,
                    quantile=0.99,
                ),
                indent=2,
            )
            + "\n"
        )
    else:
        guard = np.asarray(saved["guardNm"])
        relative = np.asarray(saved["relativeGuard"])
        by_source = saved["bySource"]
    options["predictionErrorGuardNm"] = guard.tolist()
    if args.relative_guard:
        options["predictionErrorGuardRelative"] = relative.tolist()
    if args.forecast in ("hybrid", "arx"):
        options["sourcePredictionModels"] = by_source
    selected = [g for g in groups if not args.groups or int(g["group"]) in args.groups]
    run_signature = hashlib.sha256(
        json.dumps(dict(**signature, options=options), sort_keys=True).encode()
    ).hexdigest()
    started = time.monotonic()

    def run_one(pair):
        group, source = pair
        index = int(group["group"])
        seq = dict(group["sequence"])
        seq["plannedMAPInputs"] = group["plannedMAPInputs"][source]
        seq["activationSourceIsGenerated"] = bool(source == 1)
        threshold, metadata = threshold_schedule(
            seq, many(group["predictions"])[source]
        )
        seq["backTorqueThresholdNm"] = threshold
        path = out / f"group{index}_source{source+1}.npz"
        inputs = np.asarray(group["inputs"])[source]
        notified = [time.monotonic()]
        source_options = dict(options)
        if args.forecast in ("hybrid", "arx"):
            source_options["inputForecast"] = by_source[source]["inputForecast"]
            source_options["predictionErrorGuardNm"] = by_source[source]["guardNm"]
            if args.relative_guard:
                source_options["predictionErrorGuardRelative"] = by_source[source][
                    "relativeGuard"
                ]

        def progress(k, n):
            if time.monotonic() - notified[0] > 25:
                print(
                    args.subject,
                    args.candidate,
                    index,
                    source + 1,
                    round(100 * k / n, 1),
                    "%",
                    flush=True,
                )
                notified[0] = time.monotonic()

        if path.exists() and str(np.load(path)["signature"]) == run_signature:
            z = np.load(path)
            actual, nominal, flags, timings = [
                z[k] for k in ("actual", "nominal", "flags", "timings")
            ]
        else:
            actual, nominal, flags, timings = simulate_case(
                models[source], controller, seq, inputs, source_options, progress
            )
            if not np.isfinite(actual).all():
                raise RuntimeError("Non-finite actual simulation")
            assert actual[:, 24].min() >= -1e-8 and actual[:, 24].max() <= 80 + 1e-8
            assert np.max(np.abs(np.diff(actual[:, 24]))) / seq["dt"] <= 200 + 1e-6
            assert (
                np.max(
                    abs(
                        nominal[:, 3] - many(group["predictions"])[source]["backTorque"]
                    )
                )
                < 1e-4
            )
            np.savez_compressed(
                path,
                actual=actual,
                nominal=nominal,
                flags=flags,
                timings=timings,
                threshold=threshold,
                inputs=inputs,
                signature=run_signature,
            )
        rows = evaluate(actual, nominal, seq, source, metadata)
        print(
            args.subject,
            args.candidate,
            "DONE",
            index,
            source + 1,
            "margin",
            [round(r["peakMarginNm"], 2) for r in rows],
            "valleys",
            [r["predictedValleys"] for r in rows],
            flush=True,
        )
        return index, source, seq, actual, nominal, flags, timings, rows

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        runs = list(pool.map(run_one, [(g, s) for g in selected for s in (0, 1)]))
    rows = [r for run in runs for r in run[-1]]
    result_groups = []
    for group in selected:
        predictions, nominals = [], []
        for source in (0, 1):
            run = next(
                r for r in runs if r[0] == int(group["group"]) and r[1] == source
            )
            p = prediction_dict(run[3], models[source]["parameters"])
            p.update(
                solverExitFlag=run[5],
                solveSeconds=run[6],
                backTorqueThresholdNm=run[2]["backTorqueThresholdNm"],
            )
            predictions.append(p)
            nominals.append(prediction_dict(run[4], models[source]["parameters"]))
        result_groups.append(
            dict(
                group=group["group"],
                sequence=group["sequence"],
                predictions=structs(predictions),
                nominal=structs(nominals),
                inputs=group["inputs"],
            )
        )
    savemat(
        out / "hill_mpc.mat",
        dict(
            result=dict(
                groups=structs(result_groups),
                models=structs(models),
                controller=controller,
                options=options,
            )
        ),
        long_field_names=True,
        do_compression=True,
    )
    report = dict(
        subject=args.subject,
        candidate=args.candidate,
        signature=signature,
        runSignature=run_signature,
        thresholdDefinition="0.75 * unassisted back muscle torque peak per recorded condition and activation source",
        thresholdFraction=0.75,
        options=options,
        metrics=rows,
        wallSeconds=time.monotonic() - started,
        fullNineConditions=not bool(args.groups),
        caseCountPerSource=len(rows) // 2,
        withinFiveNmBelow=sum(r["peakWithinFiveNmBelow"] for r in rows),
        maximumExcessNm=max(r["excessPeakNm"] for r in rows),
        maximumPeakShortfallNm=max(r["peakMarginNm"] for r in rows),
        addedValleyCases=sum(r["predictedValleys"] > r["nominalValleys"] for r in rows),
        changedValleyCases=sum(
            r["predictedValleys"] != r["nominalValleys"] for r in rows
        ),
        modelFrozenDuringControl=True,
        referenceSignalsUsedByOnlineController=False,
        onlineUsesFutureMeasuredActivation=False,
        independentValidation=False,
    )
    (out / "evaluation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(
        "COMPLETED",
        args.subject,
        args.candidate,
        report["withinFiveNmBelow"],
        "/",
        len(rows),
        "excess",
        report["maximumExcessNm"],
        "shortfall",
        report["maximumPeakShortfallNm"],
        flush=True,
    )


if __name__ == "__main__":
    main()
