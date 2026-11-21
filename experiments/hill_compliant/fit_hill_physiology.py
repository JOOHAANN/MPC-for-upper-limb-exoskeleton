"""Recalibrate the unchanged Hill system for shape and restrained recruitment.

All 32 parameters remain eligible for optimization. Body/load geometry and
soft limits are shared across activation sources. Project activation anchors
are not verified %MVC; force priors below are explicit engineering priors.
"""
import argparse
import csv
import ctypes
import json
import os
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat
from scipy.optimize import least_squares
from scipy.signal import find_peaks

import fit_hill_compliant as original
from fit_hill_compliant import ROOT, LOW, HIGH, SHARED, SPECIFIC, memory, pair_inputs

# The earlier l_ref >= 1 assumption forced the back unit onto its descending
# force-length branch during flexion. Allow shorter standing fibers while
# retaining the exact same Hill/tendon equations and geometry closure.
LOW = LOW.copy()
HIGH = HIGH.copy()
LOW[2] = .65

DT = .01
PHYSIOLOGY_SCALE = 1.
BACK_ACTIVATION_WEIGHT = 0.
TRANSITION_WEIGHT = 0.
ANGLE_PROFILE = "depth"
ANGLE_PROFILES = {
    "shape": {"lower": -.05, "upper": 1.95, "sequenceLower": 0., "sequenceUpper": 1.8,
              "boundaryWeight": 1.2, "valleyTolerance": .15, "excursionRatio": .75},
    "depth": {"lower": -.3, "upper": 2.3, "sequenceLower": -.3, "sequenceUpper": 2.3,
              "boundaryWeight": .15, "valleyTolerance": .15, "excursionRatio": .75},
}
DESCRIPTION = "Same compliant Hill equations; shape-priority, activation-anchored and low-cocontraction calibration"
CONTROLLER = np.array([1.5, 200, 60, 0, 120, 2, .15], dtype=float)
SEQUENCES = loadmat(ROOT / "data/processed/hill_compliant/paired_sequences.mat", simplify_cells=True)["sequences"]
DLL = Path(__file__).with_name("hill_sequence_core.dll" if os.name == "nt" else "hill_sequence_core.so")
dependencies = [DLL.with_suffix(".c"), Path(__file__).with_name("hill_pair_core.c")]
if not DLL.exists() or DLL.stat().st_mtime < max(p.stat().st_mtime for p in dependencies):
    compiler = shutil.which("gcc") or "C:/Library/CPP/mingw64-posix/bin/gcc.exe"
    flags = ["-shared", "-O3", "-static-libgcc"] if os.name == "nt" else ["-shared", "-O3", "-fPIC"]
    subprocess.run([compiler, *flags, str(dependencies[0]), "-o", str(DLL), "-lm"], check=True)
lib = ctypes.CDLL(str(DLL))
array_type = np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS")
lib.hill_sequence_integrate.argtypes = [array_type, array_type, array_type, array_type,
                                      ctypes.c_int, array_type, ctypes.c_double,
                                      ctypes.c_int, array_type, array_type]
lib.hill_sequence_integrate.restype = ctypes.c_int


def prepared(groups, generators):
    records = []
    for group in groups:
        t = SEQUENCES[group - 1]
        record = {"trial": t, "group": group, "base": [], "difference": [], "segments": []}
        for source in (0, 1):
            base = pair_inputs(t, source, generators[source])
            record["base"].append(base)
            record["difference"].append(memory(base[:, 0], .1) - memory(base[:, 0], .35))
        for segment in t["segments"]:
            indices = np.arange(int(segment["startIndex"]) - 1, int(segment["endIndex"]))
            reference = t["q"][indices]
            valleys = find_peaks(-reference, prominence=.12, distance=50)[0]
            peaks = find_peaks(reference, prominence=.12, distance=50)[0]
            anchors = np.unique(np.r_[0, valleys, peaks, len(reference) - 1])
            direction = np.zeros(len(reference))
            for left, right in zip(anchors[:-1], anchors[1:]):
                direction[left:right] = np.sign(reference[right] - reference[left])
            record["segments"].append({"indices": indices, "reference": reference,
                                       "shape": standardized(reference), "anchors": anchors,
                                       "direction": direction, "loadKg": segment["loadKg"],
                                       "speed": segment["speedCode"], "valleys": valleys,
                                       "middlePeak": int(valleys[0] + np.argmax(reference[valleys[0]:valleys[1]+1])) if len(valleys) >= 2 else -1})
        records.append(record)
    return records


def standardized(q):
    return (q - np.mean(q)) / max(np.std(q), .08)


def valley_depths(q, segment):
    """Offline phase minima; the solver receives no reference kinematics.

    Split at the original middle rise rather than forcing a pointwise match
    to the exact reference valley time. Small phase offsets remain permitted.
    """
    middle = segment["middlePeak"]
    if middle <= 0 or middle >= len(q) - 1:
        return np.array([np.min(q), np.min(q)])
    return np.array([np.min(q[:middle + 1]), np.min(q[middle:])])


def inputs(parameters, record, source):
    base = record["base"][source]
    d = record["difference"][source]
    inferred = parameters[15] * base[:, 0] + parameters[16] * np.maximum(0, d) + parameters[17] * np.maximum(0, -d)
    front = -np.expm1(-np.maximum(0, parameters[14] * base[:, 1] + inferred))
    return np.asfortranarray(np.column_stack((base[:, 0], front))), inferred


def rollout(parameters, record, source, steps=2):
    t = record["trial"]
    a, inferred = inputs(parameters, record, source)
    out = np.zeros((len(a), 23), order="F")
    cfg = CONTROLLER.copy()
    cfg[0] = t["standingAngle"]
    failed = lib.hill_sequence_integrate(np.asfortranarray(parameters), a,
                                        np.asfortranarray(t["m"]), np.asfortranarray(t["waiting"], dtype=float),
                                        len(a), np.asfortranarray([t["q"][0], t["v0"]]), DT, steps,
                                        np.asfortranarray(cfg), out)
    return out, bool(failed), a, inferred


def residual(parameters, records, source, steps=2):
    errors = []
    angle_policy = ANGLE_PROFILES[ANGLE_PROFILE]
    for record in records:
        output, failed, a, inferred = rollout(parameters, record, source, steps)
        t = record["trial"]
        for segment in record["segments"]:
            ids = segment["indices"]
            q, v, torque, back, front = output[ids, :5].T
            n = np.sqrt(len(q))
            shape = standardized(q)
            base = record["base"][source][ids, 1]
            expected_back = record["base"][source][ids, 0] ** parameters[5]
            scale = .04 + .5 * base
            reference = segment["reference"]
            # Absolute angle has deliberately little weight; shape and phase
            # order matter, while the reference itself remains unchanged.
            frames = np.r_[
                1.5 * (shape - segment["shape"]),
                .04 * (q - reference) / .4,
                .8 * (torque - t["torque"][ids]) / 70,
                .12 * np.diff(torque) / 5,
                .8 * np.maximum(0, -segment["direction"] * v) / (3 * max(np.std(q), .08)),
                2 * (np.maximum(0, angle_policy["lower"] - q) + np.maximum(0, q - angle_policy["upper"])),
                # Recruitment is anchored to its own source, never to torque
                # or reference angle. Shared Hill mechanics are not altered.
                PHYSIOLOGY_SCALE * 1.3 * (a[ids, 1] - base) / scale,
                PHYSIOLOGY_SCALE * .8 * (output[ids, 14] - base) / scale,
                PHYSIOLOGY_SCALE * .65 * np.minimum(back, front) / 25,
                PHYSIOLOGY_SCALE * .5 * inferred[ids] / scale,
                BACK_ACTIVATION_WEIGHT * (output[ids, 13] - expected_back) / (.08 + .5 * expected_back),
                angle_policy["boundaryWeight"] * original.boundary_torque(parameters, q) / 20,
                output[ids, 9] / 15,
                .4 * output[ids, 6] / 100,
            ] / n
            anchors = segment["anchors"]
            w_errors = []
            if segment["middlePeak"] >= 0:
                middle = q[segment["middlePeak"]]
                w_errors = [max(0, .15 - middle + q[v]) / .15 for v in segment["valleys"][:2]]
            combined = np.r_[frames, .6 * (shape[anchors] - segment["shape"][anchors]) / np.sqrt(len(anchors)),
                             2 * max(0, .35 - np.ptp(q)) / .35, w_errors]
            if ANGLE_PROFILE == "depth":
                predicted_depth = valley_depths(q, segment)
                reference_depth = valley_depths(reference, segment)
                # Normalized shape alone permits an implausibly shallow W.
                # Match both flexion depths loosely, allowing deeper motion.
                shallow = np.maximum(0, predicted_depth - reference_depth - .06)
                corr = np.corrcoef(q, reference)[0, 1] if np.std(q) > 1e-6 else 0.
                left, right = segment["valleys"][:2]
                middle = np.max(q[left:right + 1])
                reference_rebound = reference[segment["middlePeak"]] - np.mean(reference_depth)
                rebound = middle - np.mean(predicted_depth)
                combined = np.r_[combined, 2.5 * shallow / .2,
                                 .25 * (predicted_depth - reference_depth) / .3,
                                 1.5 * max(0, .75 * np.ptp(reference) - np.ptp(q)) / max(np.ptp(reference), .3),
                                 2 * max(0, .75 * reference_rebound - rebound) / .25,
                                 4 * max(0, .67 - corr) / .2]
            if failed:
                combined += 15
            errors.append(np.nan_to_num(combined, nan=1e5, posinf=1e5, neginf=-1e5))
        # Evaluate every waiting sample, including the release/stance
        # transition; checking only the settled tail can hide overshoot.
        waiting = np.asarray(t["waiting"], dtype=bool)
        full_q = output[:, 0]
        transition = TRANSITION_WEIGHT * np.r_[
            (full_q[waiting] - t["standingAngle"]) / (.2 * np.sqrt(np.sum(waiting))),
            2 * max(0, np.max(full_q) - angle_policy["sequenceUpper"]),
            2 * max(0, angle_policy["sequenceLower"] - np.min(full_q)),
            .15 * np.max(np.abs(original.boundary_torque(parameters, full_q))) / 30,
            max(0, np.max(output[:, 4]) - 50) / 50,
        ]
        errors.append(np.nan_to_num(transition, nan=1e5, posinf=1e5, neginf=-1e5))
    return np.concatenate(errors) / np.sqrt(3 * len(records))


def priors(parameters, seed):
    # These broad regularizers are engineering assumptions, not an anatomical
    # calibration or a published universal abdominal/back torque ratio.
    weak = .12 * (parameters - seed) / (HIGH - LOW)
    strength = PHYSIOLOGY_SCALE * np.array([.35 * np.log((parameters[0] / parameters[20]) / 10000),
                         .35 * np.log((parameters[8] / parameters[21]) / 3000)])
    recruitment = np.array([.25 * (parameters[9] - 1), .2 * np.log(parameters[14])])
    return np.r_[weak, strength, recruitment]


def seed_pair():
    pair = np.vstack([np.array(model["parameters"]) for model in original.PACKAGE["models"]])
    pair[:, [12, 13, 18, 19]] = [ .75, .4, 27, .32 ]
    pair[:, [28, 29, 30, 31]] = [300, 80, .05, 1.9]
    for source in (0, 1):
        pair[source, [0, 1, 2, 3, 5, 7]] = [500, .4, .82, 12, 1.6, .06]
        pair[source, [8, 9, 10, 11, 14, 15, 16, 17]] = [300, .8, .45, .9, 1, .01, .02, .02]
        pair[source, [20, 21, 22, 23, 27]] = [.035, .06, .26, .42, .2]
    return np.clip(pair, LOW + 1e-6 * (HIGH - LOW), HIGH - 1e-6 * (HIGH - LOW))


def metrics(pair, records, steps=10, label="candidate"):
    rows = []
    for source in (0, 1):
        for record in records:
            o, failed, a, inferred = rollout(pair[source], record, source, steps)
            t = record["trial"]
            for seg in record["segments"]:
                ids = seg["indices"]
                q, v, torque, back, front = o[ids, :5].T
                valleys = find_peaks(-q, prominence=.12, distance=50)[0]
                reference_valleys = find_peaks(-seg["reference"], prominence=.12, distance=50)[0]
                correlation = np.corrcoef(q, seg["reference"])[0, 1] if np.std(q) > 1e-6 else 0
                base = record["base"][source][ids, 1]
                total_s = pair[source, 14] * base + inferred[ids]
                predicted_depth = valley_depths(q, seg)
                reference_depth = valley_depths(seg["reference"], seg)
                left, right = seg["valleys"][:2]
                predicted_middle = np.max(q[left:right + 1])
                reference_middle = seg["reference"][seg["middlePeak"]]
                rows.append({"model": label, "group": record["group"], "source": source + 1,
                             "loadKg": seg["loadKg"], "speedCode": seg["speed"], "failed": int(failed),
                             "shapeCorrelation": float(correlation), "predictedValleys": len(valleys),
                             "referenceValleys": len(reference_valleys),
                             "thetaRMSEDegrees": float(np.sqrt(np.mean((q - seg["reference"])**2)) * 180 / np.pi),
                             "torqueRMSENm": float(np.sqrt(np.mean((torque - t["torque"][ids])**2))),
                             "backTorqueMeanNm": float(np.mean(back)), "frontTorqueMeanNm": float(np.mean(front)),
                             "backTorquePeakNm": float(np.max(back)), "frontTorquePeakNm": float(np.max(front)),
                             "frontBackMeanTorqueRatio": float(np.sum(front) / max(np.sum(back), 1e-10)),
                             "frontActivationPeak": float(np.max(o[ids, 14])),
                             "frontActivationMean": float(np.mean(o[ids, 14])),
                             "baseFrontMean": float(np.mean(base)),
                             "frontMappingFromBackPercent": float(100 * np.sum(inferred[ids]) / max(np.sum(total_s), 1e-10)),
                             "softLimitPeakNm": float(np.max(np.abs(original.boundary_torque(pair[source], q)))),
                             "frontPassivePeakNm": float(np.max(o[ids, 9])),
                             "minimumAngleRad": float(np.min(q)), "maximumAngleRad": float(np.max(q)),
                             "referenceValley1Rad": float(reference_depth[0]), "referenceValley2Rad": float(reference_depth[1]),
                             "predictedValley1Rad": float(predicted_depth[0]), "predictedValley2Rad": float(predicted_depth[1]),
                             "maximumValleyTooShallowRad": float(np.max(predicted_depth - reference_depth)),
                             "excursionRatio": float(np.ptp(q) / max(np.ptp(seg["reference"]), .01)),
                             "predictedMiddlePeakRad": float(predicted_middle), "referenceMiddlePeakRad": float(reference_middle),
                             "middleReboundRatio": float((predicted_middle - np.mean(predicted_depth)) / max(reference_middle - np.mean(reference_depth), .01))})
    return rows


def write_metrics(out, rows, name="metrics.csv"):
    with (out / name).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def fit(records, pair, out, sweeps, max_nfev, steps):
    seed = seed_pair()
    begin = time.monotonic()

    def objective(current):
        return np.concatenate([np.r_[residual(current[s], records, s, steps), priors(current[s], seed[s])] / np.sqrt(2) for s in (0, 1)])

    best = pair.copy()
    best_loss = np.sum(objective(best)**2)
    print("START LOSS", best_loss, flush=True)
    for sweep in range(sweeps):
        initial = pair.copy()
        def source_block(source):
            count, last = 0, time.monotonic()
            best_parameters = initial[source].copy()
            best_value = np.inf
            def one(z):
                nonlocal count, last, best_parameters, best_value
                parameters = initial[source].copy()
                parameters[SPECIFIC] = LOW[SPECIFIC] + z * (HIGH[SPECIFIC] - LOW[SPECIFIC])
                r = np.r_[residual(parameters, records, source, steps), priors(parameters, seed[source])]
                count += 1
                value = r @ r
                if value < best_value:
                    best_value, best_parameters = value, parameters.copy()
                if time.monotonic() - last > 30:
                    print("SOURCE", source + 1, "SWEEP", sweep + 1, "EVALUATIONS", count, "LOSS", round(r @ r, 4), flush=True)
                    last = time.monotonic()
                    savemat(out / f"block_source{source + 1}.mat", {"parameters": best_parameters,
                                                                      "loss": best_value, "stage": "partial source block"})
                return r
            z = (initial[source, SPECIFIC] - LOW[SPECIFIC]) / (HIGH[SPECIFIC] - LOW[SPECIFIC])
            result = least_squares(one, z, bounds=(np.zeros(len(z)), np.ones(len(z))),
                                   max_nfev=max_nfev, diff_step=2e-5, ftol=2e-4)
            parameters = initial[source].copy()
            parameters[SPECIFIC] = LOW[SPECIFIC] + result.x * (HIGH[SPECIFIC] - LOW[SPECIFIC])
            savemat(out / f"block_source{source + 1}.mat", {"parameters": parameters,
                                                              "loss": result.fun @ result.fun, "stage": "completed source block"})
            print("SOURCE FINISHED", source + 1, "SWEEP", sweep + 1, "NFEV", result.nfev, flush=True)
            return parameters
        # Numerical fits run independently; common mechanics stay fixed until
        # both source blocks finish. This is not a separate agent workflow.
        with ThreadPoolExecutor(max_workers=2) as pool:
            futures = [pool.submit(source_block, s) for s in (0, 1)]
            pair = np.vstack([f.result() for f in futures])
        def common(z):
            current = pair.copy()
            current[:, SHARED] = LOW[SHARED] + z * (HIGH[SHARED] - LOW[SHARED])
            return objective(current)
        z = (pair[0, SHARED] - LOW[SHARED]) / (HIGH[SHARED] - LOW[SHARED])
        result = least_squares(common, z, bounds=(np.zeros(len(z)), np.ones(len(z))),
                               max_nfev=max_nfev, diff_step=2e-5, ftol=2e-4)
        pair[:, SHARED] = LOW[SHARED] + result.x * (HIGH[SHARED] - LOW[SHARED])
        loss = np.sum(objective(pair)**2)
        print("SWEEP", sweep + 1, "LOSS", round(loss, 5), "SECONDS", round(time.monotonic() - begin, 1), flush=True)
        if loss < best_loss:
            best, best_loss = pair.copy(), loss
            savemat(out / "progress_pair.mat", {"parameters": best, "loss": best_loss, "objective": DESCRIPTION})
            write_metrics(out, metrics(best, records, steps=10), "progress_metrics.csv")
    return best, best_loss


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="hill_physiology_v2")
    parser.add_argument("--groups", type=int, nargs="+", default=[1, 2, 3])
    parser.add_argument("--sweeps", type=int, default=2)
    parser.add_argument("--max-nfev", type=int, default=70)
    parser.add_argument("--steps", type=int, default=5)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--seed-only", action="store_true")
    args = parser.parse_args()
    out = ROOT / "results" / args.output
    out.mkdir(parents=True, exist_ok=True)
    generators = [m["generator"] for m in original.PACKAGE["models"]]
    records = prepared(args.groups, generators)
    all_records = prepared([1, 2, 3], generators)
    pair = seed_pair()
    if args.resume and (out / "progress_pair.mat").exists():
        pair = loadmat(out / "progress_pair.mat", simplify_cells=True)["parameters"]
    else:
        # A small reproducible warm-start search changes source parameters only.
        for source in (0, 1):
            candidates = []
            for strength in (350, 500, 700):
                for exponent in (1.3, 1.8, 2.2):
                    for front_strength in (150, 350, 650):
                        for front_exponent in (.7, 1):
                            p = pair[source].copy()
                            p[0], p[5] = strength, exponent
                            p[8], p[9] = front_strength, front_exponent
                            r = np.r_[residual(p, records, source, args.steps), priors(p, pair[source])]
                            candidates.append((r @ r, p))
            score, pair[source] = min(candidates, key=lambda item: item[0])
            print("SEED", source + 1, "LOSS", round(score, 4), "CAPACITY/GAMMA", pair[source, [0, 5]], flush=True)
    write_metrics(out, metrics(pair, all_records), "seed_metrics.csv")
    if args.seed_only:
        savemat(out / "seed_pair.mat", {"parameters": pair})
        return
    pair, loss = fit(records, pair, out, args.sweeps, args.max_nfev, args.steps)
    cells = np.empty((1, 2), dtype=object)
    for source in (0, 1):
        cells[0, source] = {"parameters": pair[source].reshape(1, -1), "generator": generators[source],
                            "objective": DESCRIPTION, "loss": loss, "calibrationGroups": np.array(args.groups)}
    savemat(out / "physiology_models.mat", {"models": cells, "parameters": pair,
                                             "controller": CONTROLLER[1:], "objective": DESCRIPTION,
                                             "sharedParameterIndices": SHARED + 1,
                                             "referencePolicy": "Exact original Scale; shape is prioritized without rescaling references"})
    baseline = np.vstack([m["parameters"] for m in original.PACKAGE["models"]])
    rows = metrics(baseline, all_records, label="previous") + metrics(pair, all_records)
    write_metrics(out, rows)
    (out / "policy.json").write_text(json.dumps({"objective": DESCRIPTION, "groups": args.groups,
                                                 "fitSubsteps": args.steps, "replaySubsteps": 10,
                                                 "anglePointWeight": .04, "angleShapeWeight": 1.5,
                                                 "forcePriorBackN": 10000, "forcePriorFrontN": 3000,
                                                 "forcePriorStatus": "Broad engineering regularizers, not subject measurements",
                                                 "activationStatus": "Project activation is not verified %MVC",
                                                 "evaluationStatus": "Calibrated replay; not an independent holdout"}, indent=2), encoding="utf-8")
    for source in (1, 2):
        chosen = [r for r in rows if r["model"] == "candidate" and r["source"] == source]
        print("FINAL SOURCE", source, "CORRELATION", np.mean([r["shapeCorrelation"] for r in chosen]),
              "W", sum(r["predictedValleys"] == 2 for r in chosen),
              "ABD MEAN", np.mean([r["frontTorqueMeanNm"] for r in chosen]),
              "ABD PEAK", max(r["frontTorquePeakNm"] for r in chosen), flush=True)
    print("FINISHED", out, flush=True)


if __name__ == "__main__":
    main()
