"""Refit one fixed Hill/neural vector per subject and activation source.

Every recorded condition has equal weight in the mean Pearson criterion.
References, masks, angle scale and trained MAP profiles are immutable. This
is calibration on nine known conditions, not an independent validation claim.
The previously released models and NMPC results are never overwritten.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import json
from pathlib import Path
import threading
import time
from types import SimpleNamespace
import ctypes

import numpy as np
from scipy.io import savemat
from scipy.optimize import least_squares

import fit_subject as fit
from control_nine_cases import load_subject, digest
from prepare_subjects import structs

OUT = fit.ROOT / "results/pearson_recalibration"
TARGET_TORQUE = (0.80, 0.75)
TARGET_THETA = (0.70, 0.60)
MULTISIGNAL_LIB = None
EXTRA_LIBRARIES = {}


def extra_library(count):
    if count not in EXTRA_LIBRARIES:
        from native_library import load_library

        name = (
            "hill_neural_multisignal_core.c"
            if count == 15
            else "hill_neural_multisignal_bias_core.c"
        )
        source = fit.ROOT / "experiments/hill_compliant" / name
        library = load_library(
            source,
            [
                source,
                source.with_name("hill_sequence_core.c"),
                source.with_name("hill_pair_core.c"),
            ],
            fit.ROOT,
        )
        library.hill_neural_integrate.argtypes = fit.LIB.hill_neural_integrate.argtypes
        library.hill_neural_integrate.restype = ctypes.c_int
        EXTRA_LIBRARIES[count] = library
    return EXTRA_LIBRARIES[count]


class CriteriaMet(Exception):
    """End a local search after a full-resolution accepted-step check."""


def pearson(x, y):
    x = np.asarray(x) - np.mean(x)
    y = np.asarray(y) - np.mean(y)
    return float(x @ y / max(np.linalg.norm(x) * np.linalg.norm(y), 1e-12))


def standardized(x):
    return (x - np.mean(x)) / max(float(np.std(x)), 1e-6)


def evaluate(p, r, cfg, records, source, steps):
    rows, outputs = [], []
    for rec in records:
        if len(r) == 10:
            out, failed, inputs = fit.rollout(p, r, cfg, rec, source, steps)
        else:
            base = rec["base"][source]
            difference = rec["difference"][source]
            recruitment = (
                p[14] * base[:, 1]
                + p[15] * base[:, 0]
                + p[16] * np.maximum(difference, 0)
                + p[17] * np.maximum(-difference, 0)
            )
            inputs = fit.array(
                np.c_[base[:, 0], -np.expm1(-np.maximum(recruitment, 0))]
            )
            seq = rec["seq"]
            out = np.zeros((len(inputs), 24), order="F")
            failed = extra_library(len(r)).hill_neural_integrate(
                fit.array(p),
                inputs,
                fit.array(seq["m"]),
                fit.array(seq["waiting"]),
                len(inputs),
                fit.array([seq["q"][0], seq["v0"]]),
                seq["dt"],
                steps,
                fit.array(np.r_[seq["standingAngle"], cfg]),
                fit.array(r),
                out,
            )
        outputs.append((out, failed, inputs))
        seq = rec["seq"]
        for segment in rec["segments"]:
            ids = segment["ids"]
            reference = seq["torque"][ids]
            valid = np.isfinite(reference)
            rows.append(
                dict(
                    loadKg=float(segment["loadKg"]),
                    speedCode=int(segment["speedCode"]),
                    torqueR=pearson(out[ids, 2][valid], reference[valid]),
                    thetaR=pearson(out[ids, 0], seq["q"][ids]),
                    torqueRMSE=float(
                        np.sqrt(np.mean((out[ids, 2][valid] - reference[valid]) ** 2))
                    ),
                    referenceFraction=float(valid.mean()),
                    failed=bool(failed),
                    excursionRatio=float(
                        np.ptp(out[ids, 0]) / max(np.ptp(seq["q"][ids]), 1e-8)
                    ),
                    minimumTheta=float(out[:, 0].min()),
                    maximumTheta=float(out[:, 0].max()),
                    abdomenBackRatio=float(
                        out[ids, 4].mean() / max(out[ids, 3].mean(), 1e-8)
                    ),
                    abdomenPeakNm=float(out[ids, 4].max()),
                    backPeakNm=float(out[ids, 3].max()),
                    boundaryPeakNm=float(abs(fit.boundary_torque(p, out[:, 0])).max()),
                )
            )
    return rows, outputs


def summary(rows):
    return dict(
        meanTorqueR=float(np.mean([v["torqueR"] for v in rows])),
        meanThetaR=float(np.mean([v["thetaR"] for v in rows])),
        meanTorqueRMSE=float(np.mean([v["torqueRMSE"] for v in rows])),
        failed=sum(v["failed"] for v in rows),
        minimumTheta=min(v["minimumTheta"] for v in rows),
        maximumTheta=max(v["maximumTheta"] for v in rows),
        meanExcursionRatio=float(np.mean([v["excursionRatio"] for v in rows])),
        maximumAbdomenBackRatio=max(v["abdomenBackRatio"] for v in rows),
        boundaryPeakNm=max(v["boundaryPeakNm"] for v in rows),
    )


def fit_one(subject, source, args):
    folder = args.output / subject / f"source{source+1}"
    if args.resume and (folder / "evaluation.json").exists():
        previous = json.loads((folder / "evaluation.json").read_text())
        if previous["final"]["passed"]:
            print(
                subject,
                source + 1,
                "RETAINED passing full-resolution result",
                flush=True,
            )
            return previous
    if args.input_directory is None:
        package, models, groups, model_path, replay_path = load_subject(subject)
    else:
        from scipy.io import loadmat

        model_path = replay_path = (
            args.input_directory / subject / "calibration_input.mat"
        )
        package = loadmat(model_path, simplify_cells=True)
        models = fit.many(package["models"])
        groups = fit.many(package["groups"])[:9]
    model = models[source]
    cfg = package["controller"]
    records = [
        fit.record(g["sequence"], model["generator"], args.stride) for g in groups
    ]
    full_records = [fit.record(g["sequence"], model["generator"]) for g in groups]
    folder.mkdir(parents=True, exist_ok=True)
    original = np.r_[model["parameters"], model["neuralParameters"]]
    low = np.r_[fit.LOW, fit.NLOW].copy()
    high = np.r_[fit.HIGH, fit.NHIGH].copy()
    low[2] = 0.65
    low[32 + 6] = 0.0  # Remove the former forced 70% feedback.
    high[32 + 8] = 1.5
    high[32 + 5] = 80.0
    if args.multisignal:
        original = np.r_[original, np.zeros(5)]
        low = np.r_[low, [-12, -12, -4, -6, -8]]
        high = np.r_[high, [12, 12, 4, 6, 8]]
        if args.intent_bias:
            original = np.r_[original, 0.0]
            low = np.r_[low, -0.75]
            high = np.r_[high, 1.25]
    span = high - low
    baseline_rows, _ = evaluate(
        original[:32], original[32:], cfg, full_records, source, 10
    )
    baseline = summary(baseline_rows)
    print(subject, source + 1, "BASELINE", json.dumps(baseline), flush=True)
    lock = threading.Lock()
    best = [np.inf, original.copy()]
    counter = [0, time.monotonic()]

    def residual(z):
        vector = low + span * z
        p, r = vector[:32], vector[32:]
        rows, outputs = evaluate(p, r, cfg, records, source, args.steps)
        values = []
        for rec, (out, failed, inputs) in zip(records, outputs):
            seq = rec["seq"]
            for segment in rec["segments"]:
                ids = segment["ids"]
                n = np.sqrt(len(ids))
                reference = seq["torque"][ids]
                valid = np.isfinite(reference)
                predicted = out[ids, 2][valid]
                target = reference[valid]
                values.extend(
                    [
                        1.0
                        * (standardized(predicted) - standardized(target))
                        / np.sqrt(valid.sum()),
                        0.6
                        * (standardized(out[ids, 0]) - standardized(seq["q"][ids]))
                        / n,
                        0.45 * (predicted - target) / (80 * np.sqrt(valid.sum())),
                        0.12 * (out[ids, 0] - seq["q"][ids]) / (0.5 * n),
                        0.6
                        * np.maximum(0, out[ids, 4] - 0.28 * out[ids, 3] - 5)
                        / (25 * n),
                        0.12 * np.minimum(out[ids, 3], out[ids, 4]) / (30 * n),
                        np.array(
                            [
                                12
                                * max(
                                    0,
                                    out[ids, 4].mean() / max(out[ids, 3].mean(), 1e-8)
                                    - 0.35,
                                )
                            ]
                        ),
                        np.array(
                            [
                                1.5
                                * max(
                                    0,
                                    0.75 * np.ptp(seq["q"][ids]) - np.ptp(out[ids, 0]),
                                )
                                / 0.25
                            ]
                        ),
                        np.array(
                            [
                                max(
                                    0,
                                    float(out[ids, 0].min() - seq["q"][ids].min())
                                    - 0.20,
                                )
                                / 0.25
                            ]
                        ),
                    ]
                )
            count = np.sqrt(len(out))
            values.extend(
                [
                    4
                    * (
                        np.maximum(0, -0.10 - out[:, 0])
                        + np.maximum(0, out[:, 0] - 1.95)
                    )
                    / count,
                    0.6 * fit.boundary_torque(p, out[:, 0]) / (25 * count),
                    np.array([10 * max(0, float(out[:, 0].max()) - 1.95) / 0.2]),
                    np.array([10 * max(0, -0.10 - float(out[:, 0].min())) / 0.2]),
                    np.array(
                        [
                            10
                            * max(
                                0,
                                float(abs(fit.boundary_torque(p, out[:, 0])).max())
                                - 25,
                            )
                            / 25
                        ]
                    ),
                    np.full(1, 30.0 if failed else 0.0),
                ]
            )
            settled = np.asarray(seq["waiting"], dtype=bool).copy()
            for segment in rec["segments"]:
                end = int(segment["ids"][-1]) + 1
                settled[end : min(len(out), end + round(1.4 / seq["dt"]))] = False
            if settled.any():
                values.append(
                    0.4
                    * (out[settled, 0] - seq["standingAngle"])
                    / (0.12 * np.sqrt(settled.sum()))
                )
        stats = summary(rows)
        values.extend(
            [
                np.array(
                    [30 * max(0, TARGET_TORQUE[source] + 0.020 - stats["meanTorqueR"])]
                ),
                np.array(
                    [40 * max(0, TARGET_THETA[source] + 0.035 - stats["meanThetaR"])]
                ),
                0.025 * (vector - original) / span,
            ]
        )
        result = np.nan_to_num(np.concatenate(values), nan=1e4, posinf=1e4, neginf=-1e4)
        loss = float(result @ result)
        with lock:
            counter[0] += 1
            if loss < best[0]:
                best[:] = [loss, vector.copy()]
            if time.monotonic() - counter[1] > 25:
                print(
                    subject,
                    source + 1,
                    "calls",
                    counter[0],
                    "bestLoss",
                    round(best[0], 5),
                    "current rT/rQ",
                    round(stats["meanTorqueR"], 4),
                    round(stats["meanThetaR"], 4),
                    flush=True,
                )
                savemat(folder / "checkpoint.mat", dict(vector=best[1], loss=best[0]))
                counter[1] = time.monotonic()
        return result

    start = original.copy()
    if args.resume and (folder / "checkpoint.mat").exists():
        from scipy.io import loadmat

        start = np.asarray(
            loadmat(folder / "checkpoint.mat", simplify_cells=True)["vector"]
        )
        if len(start) < len(original):
            start = np.r_[start, np.zeros(len(original) - len(start))]
    start_vector = start.copy()
    start = np.clip((start - low) / span, 1e-5, 1 - 1e-5)
    with ThreadPoolExecutor(max_workers=args.jacobian_workers) as pool:

        def jacobian(z):
            base = residual(z)
            vector = low + span * z
            coarse_rows, _ = evaluate(
                vector[:32], vector[32:], cfg, records, source, args.steps
            )
            check = summary(coarse_rows)
            if (
                check["meanTorqueR"] > TARGET_TORQUE[source] + 0.025
                and check["meanThetaR"] > TARGET_THETA[source] + 0.025
                and quality_passed(check)
            ):
                full_rows, _ = evaluate(
                    vector[:32], vector[32:], cfg, full_records, source, 10
                )
                exact = summary(full_rows)
                if (
                    exact["meanTorqueR"] > TARGET_TORQUE[source] + 0.020
                    and exact["meanThetaR"] > TARGET_THETA[source] + 0.020
                    and quality_passed(exact)
                ):
                    best[:] = [float(base @ base), vector.copy()]
                    raise CriteriaMet

            def column(index):
                step = (
                    args.difference_step
                    if z[index] < 1 - args.difference_step
                    else -args.difference_step
                )
                changed = z.copy()
                changed[index] += step
                return (residual(changed) - base) / step

            return np.column_stack(list(pool.map(column, range(len(z)))))

        try:
            if args.evaluate_only:
                best[:] = [0.0, start_vector]
                result = SimpleNamespace(status=11)
            else:
                result = least_squares(
                    residual,
                    start,
                    jac=jacobian,
                    bounds=(np.zeros(len(low)), np.ones(len(low))),
                    max_nfev=args.max_nfev,
                    ftol=1e-6,
                    xtol=1e-6,
                    gtol=1e-5,
                )
        except CriteriaMet:
            result = SimpleNamespace(status=10)
    vector = best[1]
    savemat(folder / "checkpoint.mat", dict(vector=vector, loss=best[0]))
    updated = dict(model, parameters=vector[:32], neuralParameters=vector[32:])
    rows, _ = evaluate(vector[:32], vector[32:], cfg, full_records, source, 10)
    final = summary(rows)
    final["correlationTargetsPassed"] = bool(
        final["meanTorqueR"] > TARGET_TORQUE[source]
        and final["meanThetaR"] > TARGET_THETA[source]
    )
    final["qualityChecksPassed"] = quality_passed(final)
    final["passed"] = final["correlationTargetsPassed"] and final["qualityChecksPassed"]
    metadata = dict(
        subject=subject,
        source=source + 1,
        baseline=baseline,
        final=final,
        targetTorqueMeanR=TARGET_TORQUE[source],
        targetThetaMeanR=TARGET_THETA[source],
        conditionCount=9,
        oneFixedVectorAcrossConditions=True,
        referenceChanged=False,
        sourceModelSHA256=digest(model_path),
        sourceReplaySHA256=digest(replay_path),
        optimizationStatus=int(result.status),
        loss=best[0],
        rows=rows,
    )
    savemat(
        folder / "model.mat",
        dict(model=updated, controller=cfg),
        do_compression=True,
        long_field_names=True,
    )
    (folder / "evaluation.json").write_text(json.dumps(metadata, indent=2) + "\n")
    with (folder / "metrics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(subject, source + 1, "FINAL", json.dumps(final), flush=True)
    return metadata


def quality_passed(stats):
    return bool(
        not stats["failed"]
        and stats["minimumTheta"] > -0.15
        and stats["maximumTheta"] < 1.98
        and stats["maximumAbdomenBackRatio"] < 0.38
        and stats["boundaryPeakNm"] < 30
        and stats["meanExcursionRatio"] > 0.70
    )


def main():
    global MULTISIGNAL_LIB
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subjects", nargs="+", default=["528", "725", "730"])
    parser.add_argument("--sources", nargs="+", type=int, default=[1, 2])
    parser.add_argument("--output", type=Path, default=OUT)
    parser.add_argument("--input-directory", type=Path)
    parser.add_argument("--max-nfev", type=int, default=65)
    parser.add_argument("--stride", type=int, default=2)
    parser.add_argument("--steps", type=int, default=6)
    parser.add_argument("--difference-step", type=float, default=1e-4)
    parser.add_argument("--jacobian-workers", type=int, default=4)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--multisignal", action="store_true")
    parser.add_argument("--intent-bias", action="store_true")
    parser.add_argument("--evaluate-only", action="store_true")
    args = parser.parse_args()
    if args.multisignal:
        extra_library(16 if args.intent_bias else 15)
    catalog = {
        "528": "240528Pilot",
        "617": "240617Pilot",
        "725": "240725Exp_A01",
        "730": "240730Exp_A01",
    }
    jobs = [
        (catalog.get(s, s), source - 1)
        for s in args.subjects
        for source in args.sources
    ]
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        rows = list(pool.map(lambda pair: fit_one(*pair, args), jobs))
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "summary.json").write_text(json.dumps(rows, indent=2) + "\n")


if __name__ == "__main__":
    main()
