"""Personalize the existing 42-parameter neural/compliant-Hill pair.

The experimental and MAP sources use identical equations. Body/load geometry
and soft boundaries are shared; every physical and neural parameter is eligible
for bounded calibration. One fixed vector per source serves all subject cases.
Only designated training records enter the objective or model selection.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import ctypes
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.io import loadmat, savemat
from scipy.optimize import least_squares
from scipy.signal import find_peaks

from audit_subjects import ROOT, OUT
from prepare_subjects import combine, structs

sys.path.insert(0, str(ROOT / "experiments/hill_compliant"))
from fit_hill_compliant import (
    LOW,
    HIGH,
    SHARED,
    SPECIFIC,
    memory,
    pair_inputs,
    boundary_torque,
)

sys.path.insert(0, str(ROOT / "experiments/legacy_mpc"))
from run_hill_mpc import prediction_dict

LOW = LOW.copy()
LOW[2] = 0.65
NLOW = np.array([0.7, 0.3, 0.01, 30, 10, 3, 0.7, 0, 0, 0])
NHIGH = np.array([3.0, 2.0, 0.25, 250, 90, 40, 1.0, 8, 0.8, 2.5])
NATIVE = ROOT / "experiments/hill_compliant/hill_neural_core.c"
dependencies = [
    NATIVE,
    NATIVE.with_name("hill_sequence_core.c"),
    NATIVE.with_name("hill_pair_core.c"),
]
from native_library import load_library

LIB = load_library(NATIVE, dependencies, ROOT)
ARRAY = np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS")
LIB.hill_neural_integrate.argtypes = [ARRAY] * 4 + [
    ctypes.c_int,
    ARRAY,
    ctypes.c_double,
    ctypes.c_int,
    ARRAY,
    ARRAY,
    ARRAY,
]
DIVISIVE_NATIVE = NATIVE.with_name("hill_neural_divisive_core.c")
DIVISIVE_LIB = load_library(DIVISIVE_NATIVE, [DIVISIVE_NATIVE, *dependencies], ROOT)
DIVISIVE_LIB.hill_neural_integrate.argtypes = LIB.hill_neural_integrate.argtypes
SOFT_TONIC_NATIVE = NATIVE.with_name("hill_neural_soft_tonic_core.c")
SOFT_TONIC_LIB = load_library(
    SOFT_TONIC_NATIVE, [SOFT_TONIC_NATIVE, *dependencies], ROOT
)
SOFT_TONIC_LIB.hill_neural_integrate.argtypes = LIB.hill_neural_integrate.argtypes
INTEGRATORS = {
    "subtractive": LIB,
    "divisive": DIVISIVE_LIB,
    "soft_tonic": SOFT_TONIC_LIB,
}
INTENT_DECODER = "subtractive"
LABELS = ("Experimental activation", "MAP activation")
SHAPE_TAIL_WEIGHT = 0.0
ANGLE_VALUE_WEIGHT = 0.25
VALLEY_DEPTH_WEIGHT = 1.5
LANDMARK_WEIGHT = 0.0
DIRECTION_WEIGHT = 0.0
ASSISTANCE_RESPONSE_WEIGHT = 0.0
ROLLOUT_POOL = None


def absolute_jacobian(fun, step):
    """Differentiate normalized parameters without vanishing steps at bounds.

    A relative step at z=1e-5 can be smaller than the native equilibrium
    solve's resolution. An absolute normalized step also lets bound-adjacent
    muscle/feedback parameters participate in calibration reliably.
    """

    def jacobian(z):
        base = fun(z)
        columns = []
        for j in range(len(z)):
            change = step if z[j] <= 1.0 - step else -step
            trial = z.copy()
            trial[j] += change
            columns.append((fun(trial) - base) / change)
        return np.column_stack(columns)

    return jacobian


def landmark_residual(q, reference, valleys):
    """Offline, loosely match both bends and the intervening upright return."""
    if len(valleys) != 2:
        return np.zeros(7)
    left, right = map(int, valleys)
    middle = left + int(np.argmax(reference[left : right + 1]))
    predicted_depth = np.array([np.min(q[: middle + 1]), np.min(q[middle:])])
    reference_depth = reference[[left, right]]
    rebound = np.max(q[left : right + 1]) - np.mean(predicted_depth)
    reference_rebound = reference[middle] - np.mean(reference_depth)
    correlation = float(np.corrcoef(q, reference)[0, 1]) if np.std(q) > 1e-8 else 0.0
    return np.r_[
        2.0 * np.maximum(predicted_depth - reference_depth - 0.12, 0) / 0.25,
        0.3 * (predicted_depth - reference_depth) / 0.35,
        2.0 * max(0, 0.65 * reference_rebound - rebound) / 0.25,
        2.0 * max(0, 0.70 * np.ptp(reference) - np.ptp(q)) / 0.35,
        4.0 * max(0, 0.70 - correlation),
    ]


def array(value):
    return np.asfortranarray(value, dtype=float)


def many(value):
    return np.asarray(value, dtype=object).reshape(-1).tolist()


def record(seq, generator, stride=1):
    fields = (
        "q",
        "torque",
        "a",
        "front",
        "m",
        "loadKg",
        "speedCode",
        "stageId",
        "waiting",
        "referenceIsFilled",
        "t",
    )
    s = {**seq, **{k: np.asarray(seq[k])[::stride] for k in fields}}
    if "torqueReferenceValid" in seq:
        s["torqueReferenceValid"] = np.asarray(seq["torqueReferenceValid"])[::stride]
    s["dt"] = float(seq["dt"]) * stride
    segments = []
    for seg in many(seq["segments"]):
        ids = np.flatnonzero(
            (np.arange(len(s["q"])) * stride >= int(seg["startIndex"]) - 1)
            & (np.arange(len(s["q"])) * stride < int(seg["endIndex"]))
        )
        ref = s["q"][ids]
        valleys = find_peaks(
            -ref, prominence=0.12, distance=max(1, round(0.5 / s["dt"]))
        )[0]
        segments.append(
            {
                **seg,
                "ids": ids,
                "reference": ref,
                "shape": standardized(ref),
                "valleys": valleys,
            }
        )
    # pair_inputs' memories are defined at 100 Hz, so downsample after computing.
    base = [pair_inputs(seq, source, generator)[::stride] for source in (0, 1)]
    differences = [
        (
            memory(pair_inputs(seq, source, generator)[:, 0], 0.1)
            - memory(pair_inputs(seq, source, generator)[:, 0], 0.35)
        )[::stride]
        for source in (0, 1)
    ]
    return {"seq": s, "segments": segments, "base": base, "difference": differences}


def standardized(q):
    return (q - q.mean()) / max(float(q.std()), 0.08)


def rollout(p, r, cfg, record, source, steps=6, intent_decoder=None):
    s = record["seq"]
    base = record["base"][source]
    d = record["difference"][source]
    inferred = p[15] * base[:, 0] + p[16] * np.maximum(0, d) + p[17] * np.maximum(0, -d)
    a = np.c_[base[:, 0], -np.expm1(-np.maximum(0, p[14] * base[:, 1] + inferred))]
    n = len(a)
    o = np.zeros((n, 24), order="F")
    integrator = INTEGRATORS[intent_decoder or INTENT_DECODER]
    failed = integrator.hill_neural_integrate(
        array(p),
        array(a),
        array(s["m"]),
        array(s["waiting"]),
        n,
        array([s["q"][0], s["v0"]]),
        s["dt"],
        steps,
        array(np.r_[s["standingAngle"], cfg]),
        array(r),
        o,
    )
    return o, bool(failed), a


def residual(p, r, cfg, records, source, seed_p, seed_r, steps):
    values = []

    def predict(rec):
        return rollout(p, r, cfg, rec, source, steps)

    predictions = (
        list(ROLLOUT_POOL.map(predict, records))
        if ROLLOUT_POOL
        else list(map(predict, records))
    )
    for rec, (o, failed, a) in zip(records, predictions):
        s = rec["seq"]
        # The full waiting trajectory matters: a plausible isolated movement
        # can otherwise hide accumulated overshoot in a stitched sequence.
        waiting = np.asarray(s["waiting"], dtype=bool)
        settled = waiting.copy()
        for seg in rec["segments"]:
            end = int(seg["ids"][-1]) + 1
            settled[end : min(len(settled), end + round(1.6 / s["dt"]))] = False
        if np.any(settled):
            values.append(
                0.5 * (o[settled, 0] - 1.5) / (0.04 * np.sqrt(np.sum(settled)))
            )
        values.append(
            2
            * (np.maximum(0, -0.2 - o[:, 0]) + np.maximum(0, o[:, 0] - 2.1))
            / np.sqrt(len(o))
        )
        if s.get("activePerturbationStress", False):
            values.append(4 * np.maximum(0, o[:, 0] - 1.9) / np.sqrt(len(o)))
            values.append(
                0.5
                * np.maximum(0, o[:, 4] - 0.2 * o[:, 3] - 5)
                / (25 * np.sqrt(len(o)))
            )
        for seg in rec["segments"]:
            ids = seg["ids"]
            q = o[ids, 0]
            back = o[ids, 3]
            front = o[ids, 4]
            n = np.sqrt(len(ids))
            base = rec["base"][source][ids, 1]
            scale = 0.04 + 0.5 * base
            errors = [
                0.9
                * np.where(
                    np.asarray(s.get("torqueReferenceValid", np.isfinite(s["torque"])))[
                        ids
                    ],
                    o[ids, 2] - s["torque"][ids],
                    0.0,
                )
                / 70,
                0.10 * np.diff(o[ids, 2]) / 8,
                2.0 * np.maximum(0, front - 0.22 * back - 3) / 20,
                0.65 * np.minimum(front, back) / 25,
                0.7 * (a[ids, 1] - base) / scale,
                0.6 * (o[ids, 14] - base) / scale,
                float(seg.get("boundaryResidualWeight", 0.25))
                * boundary_torque(p, q)
                / 20,
                4 * (np.maximum(0, -0.2 - q) + np.maximum(0, q - 2.1)),
                0.4 * o[ids, 9] / 15,
            ]
            if bool(seg["angleUsable"]):
                ref = seg["reference"]
                shape_error = standardized(q) - seg["shape"]
                # A smooth quartic tail gives poorly fitting cases more weight
                # without introducing case-specific model parameters.
                shape_boost = np.sqrt(1 + SHAPE_TAIL_WEIGHT * np.mean(shape_error**2))
                errors.extend(
                    [
                        1.5 * shape_boost * shape_error,
                        ANGLE_VALUE_WEIGHT * (q - ref) / 0.35,
                    ]
                )
                if DIRECTION_WEIGHT and len(seg["valleys"]) == 2:
                    left, right = map(int, seg["valleys"])
                    middle = left + int(np.argmax(ref[left : right + 1]))
                    direction = np.ones(len(q))
                    direction[:left] = -1
                    direction[middle:right] = -1
                    margin = round(0.25 / s["dt"])
                    for turn in (left, middle, right):
                        direction[
                            max(0, turn - margin) : min(len(q), turn + margin + 1)
                        ] = 0
                    errors.append(
                        DIRECTION_WEIGHT
                        * np.maximum(0, -direction * o[ids, 1])
                        / (3 * max(np.std(q), 0.08))
                    )
                # Match local flexion depths loosely without moving the target.
                if len(seg["valleys"]):
                    depths = []
                    for valley in seg["valleys"]:
                        lo = max(0, valley - round(0.35 / s["dt"]))
                        hi = min(len(q), valley + round(0.35 / s["dt"]) + 1)
                        depths.append(
                            max(0, np.min(q[lo:hi]) - ref[valley] - 0.10) / 0.2
                        )
                    errors.append(
                        VALLEY_DEPTH_WEIGHT
                        * np.repeat(depths, max(1, round(len(q) / 20 / len(depths))))
                    )
            priority = np.sqrt(float(seg.get("fitWeight", 1.0)))
            v = priority * np.concatenate([np.ravel(e) for e in errors]) / n
            if failed:
                v = v + 1e2
            values.append(np.nan_to_num(v, nan=1e4, posinf=1e4, neginf=-1e4))
            if LANDMARK_WEIGHT and bool(seg["angleUsable"]):
                values.append(
                    LANDMARK_WEIGHT
                    * priority
                    * landmark_residual(q, seg["reference"], seg["valleys"])
                )
            end = int(ids[-1]) + 1
            tail = o[end : min(len(o), end + round(1.5 / s["dt"])), :]
            if len(tail):
                values.append(
                    np.r_[
                        0.2 * (tail[:, 0] - 1.5) / 0.15,
                        0.12 * np.maximum(0, tail[:, 3] - 25) / 30,
                    ]
                    * priority
                    / np.sqrt(len(tail))
                )
    # Broad priors discourage artificial antagonist strength or coactivation.
    prior = 0.08 * (p - seed_p) / (HIGH - LOW)
    neural_prior = 0.05 * (r - seed_r) / (NHIGH - NLOW)
    strength = np.array(
        [0.2 * np.log((p[0] / p[20]) / 10000), 0.2 * np.log((p[8] / p[21]) / 3000)]
    )
    return np.r_[
        np.concatenate(values) / np.sqrt(max(1, len(records))),
        prior,
        neural_prior,
        strength,
        # Engineering robustness prior, not an identified human measurement:
        # retain motion-request authority under a representative 60 Nm assist.
        ASSISTANCE_RESPONSE_WEIGHT * max(0.0, 60.0 * (1.0 - r[6]) / r[3] / 0.12 - 1.0),
    ]


def score_rows(models, cfg, records, source_label="personalized", steps=10):
    rows = []
    for source, model in enumerate(models):
        for rec in records:
            o, bad, a = rollout(
                model["parameters"],
                model["neuralParameters"],
                cfg,
                rec,
                source,
                steps,
                model.get("intentDecoder", "subtractive"),
            )
            s = rec["seq"]
            for seg in rec["segments"]:
                ids = seg["ids"]
                q = o[ids, 0]
                ref = s["q"][ids]
                torque_valid = np.asarray(
                    s.get("torqueReferenceValid", np.isfinite(s["torque"])), dtype=bool
                )[ids]
                rows.append(
                    {
                        "model": source_label,
                        "source": LABELS[source],
                        "recordId": seg["recordId"],
                        "loadKg": float(seg["loadKg"]),
                        "speedCode": int(seg["speedCode"]),
                        "split": seg["split"],
                        "angleReferenceUsable": bool(seg["angleUsable"]),
                        "failed": bool(bad),
                        "thetaRMSEDegrees": float(
                            np.sqrt(np.mean((q - ref) ** 2)) * 180 / np.pi
                        ),
                        "shapeCorrelation": (
                            float(np.corrcoef(q, ref)[0, 1])
                            if np.std(q) > 1e-8
                            else 0.0
                        ),
                        "referenceValleys": len(
                            find_peaks(-ref, prominence=0.12, distance=50)[0]
                        ),
                        "predictedValleys": len(
                            find_peaks(-q, prominence=0.12, distance=50)[0]
                        ),
                        "torqueRMSENm": (
                            float(
                                np.sqrt(
                                    np.mean(
                                        (
                                            o[ids, 2][torque_valid]
                                            - s["torque"][ids][torque_valid]
                                        )
                                        ** 2
                                    )
                                )
                            )
                            if torque_valid.any()
                            else None
                        ),
                        "torqueReferenceCoverageFraction": float(
                            np.mean(
                                np.asarray(
                                    s.get(
                                        "torqueReferenceValid", np.isfinite(s["torque"])
                                    )
                                )[ids]
                            )
                        ),
                        "backTorquePeakNm": float(o[ids, 3].max()),
                        "backTorqueMeanNm": float(o[ids, 3].mean()),
                        "abdomenTorquePeakNm": float(o[ids, 4].max()),
                        "abdomenTorqueMeanNm": float(o[ids, 4].mean()),
                        "abdomenBackRatio": float(
                            o[ids, 4].sum() / max(o[ids, 3].sum(), 1e-9)
                        ),
                        "minAngleRad": float(q.min()),
                        "maxAngleRad": float(q.max()),
                        "sequenceMinAngleRad": float(o[:, 0].min()),
                        "sequenceMaxAngleRad": float(o[:, 0].max()),
                        "wholeSequenceBoundaryPeakNm": float(
                            np.max(
                                np.abs(boundary_torque(model["parameters"], o[:, 0]))
                            )
                        ),
                        "postActionBoundaryPeakNm": (
                            float(
                                np.max(
                                    np.abs(
                                        boundary_torque(
                                            model["parameters"],
                                            o[int(ids[-1]) + 1 :, 0],
                                        )
                                    )
                                )
                            )
                            if int(ids[-1]) + 1 < len(o)
                            else 0.0
                        ),
                        "boundaryPeakNm": float(
                            np.max(np.abs(boundary_torque(model["parameters"], q)))
                        ),
                    }
                )
    return rows


def export_replay(models, cfg, sequences, generator, target):
    groups = []
    for i, seq in enumerate(sequences, 1):
        rec = record(seq, generator)
        predictions = []
        inputs = []
        for source, model in enumerate(models):
            o, bad, a = rollout(
                model["parameters"],
                model["neuralParameters"],
                cfg,
                rec,
                source,
                10,
                model.get("intentDecoder", "subtractive"),
            )
            if bad:
                raise RuntimeError(
                    f"Unstable full-resolution replay: {seq['name']}, source {source}"
                )
            p = model["parameters"]
            full = np.zeros((len(o), 27))
            full[:, :24] = o
            full[:, 25] = boundary_torque(p, o[:, 0])
            full[:, 26] = -9.81 * p[18] * p[19] * np.cos(o[:, 0]) - 9.81 * np.asarray(
                seq["m"]
            ) * p[12] * np.cos(o[:, 0] - p[13])
            predictions.append(prediction_dict(full, p))
            inputs.append(a)
        groups.append(
            {
                "group": i,
                "sequence": seq,
                "inputs": np.array(inputs),
                "predictions": predictions,
            }
        )
    result = {"groups": structs(groups), "models": structs(models), "controller": cfg}
    savemat(
        target / "hill_sequence.mat",
        {"result": result},
        do_compression=True,
        long_field_names=True,
    )


def main():
    global ROLLOUT_POOL, INTENT_DECODER
    global SHAPE_TAIL_WEIGHT, ANGLE_VALUE_WEIGHT, VALLEY_DEPTH_WEIGHT, LANDMARK_WEIGHT, DIRECTION_WEIGHT, ASSISTANCE_RESPONSE_WEIGHT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("subject")
    parser.add_argument(
        "--intent-decoder",
        choices=("subtractive", "divisive", "soft_tonic"),
        default="subtractive",
    )
    parser.add_argument("--sweeps", type=int, default=2)
    parser.add_argument("--max-nfev", type=int, default=35)
    parser.add_argument("--steps", type=int, default=6)
    parser.add_argument(
        "--rollout-workers",
        type=int,
        default=1,
        help="Parallel independent trial solves; one keeps legacy execution order.",
    )
    parser.add_argument(
        "--absolute-difference-step",
        type=float,
        default=0.0,
        help="Absolute finite-difference step in normalized parameter coordinates; zero retains legacy relative steps.",
    )
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--initial-model",
        type=Path,
        help="Explicit warm start; checkpoints still require identical prepared input hashes.",
    )
    parser.add_argument(
        "--continuous-training",
        action="store_true",
        help="Fit training cycles in their original order through standing waits.",
    )
    parser.add_argument(
        "--robust-recovery",
        action="store_true",
        help="Keep dissipative/braking reserve and include perturbed standing recovery in calibration.",
    )
    parser.add_argument("--minimum-feedback-blend", type=float, default=0.7)
    parser.add_argument(
        "--minimum-motion-stiffness",
        type=float,
        default=30.0,
        help="Engineering lower bound on the muscle-mediated motion request gain (Nm/rad).",
    )
    parser.add_argument("--angle-value-weight", type=float, default=0.25)
    parser.add_argument("--valley-depth-weight", type=float, default=1.5)
    parser.add_argument(
        "--shape-tail-weight",
        type=float,
        default=0.0,
        help="Extra smooth tail weight on poorly fitting angle shapes; zero preserves the earlier objective.",
    )
    parser.add_argument(
        "--maximum-emg-lead",
        type=float,
        default=0.8,
        help="Upper bound on causal filtered-EMG derivative gain (seconds), not future preview.",
    )
    parser.add_argument(
        "--maximum-flexion-request",
        type=float,
        default=40.0,
        help="Neural antagonist-request saturation bound in Nm; actual muscle torque remains computed by Hill.",
    )
    parser.add_argument(
        "--neural-warmup-nfev",
        type=int,
        default=0,
        help="Optional neural-only warmup before optimizing all physical/neural parameters.",
    )
    parser.add_argument("--landmark-weight", type=float, default=0.0)
    parser.add_argument("--direction-weight", type=float, default=0.0)
    parser.add_argument(
        "--assistance-response-weight",
        type=float,
        default=0.0,
        help="Soft engineering prior on muscle-mediated motion authority under assistance; zero preserves legacy calibration.",
    )
    parser.add_argument("--decoder-multistart", action="store_true")
    parser.add_argument("--active-recovery-stress", action="store_true")
    parser.add_argument(
        "--subject-body-search",
        action="store_true",
        help="Select an initial shared M/L on this subject's training data, then optimize all parameters.",
    )
    args = parser.parse_args()
    INTENT_DECODER = args.intent_decoder
    if args.rollout_workers < 1:
        parser.error("Rollout worker count must be positive.")
    if args.rollout_workers > 1:
        ROLLOUT_POOL = ThreadPoolExecutor(max_workers=args.rollout_workers)
    if not 0 <= args.absolute_difference_step < 0.01:
        parser.error("Absolute normalized finite-difference step must lie in [0,.01).")
    if not 0.0 <= args.minimum_feedback_blend < 1.0:
        parser.error(
            "The feedback blend lower bound must be nonnegative and less than 1."
        )
    if not 0 < args.minimum_motion_stiffness < NHIGH[3]:
        parser.error(
            "Motion stiffness lower bound must be positive and below its upper bound."
        )
    if args.sweeps < 0 or args.steps < 1 or args.max_nfev < 1:
        parser.error("Sweeps must be nonnegative; steps and max-nfev must be positive.")
    if (
        args.maximum_emg_lead <= 0
        or args.maximum_flexion_request <= 10
        or args.neural_warmup_nfev < 0
    ):
        parser.error(
            "EMG lead must be positive, request bound above 10 Nm and warmup nonnegative."
        )
    NHIGH[8] = args.maximum_emg_lead
    NHIGH[5] = args.maximum_flexion_request
    if args.shape_tail_weight < 0:
        parser.error("Shape-tail weight must be nonnegative.")
    SHAPE_TAIL_WEIGHT = args.shape_tail_weight
    if min(args.angle_value_weight, args.valley_depth_weight) < 0:
        parser.error("Angle-value and valley-depth weights must be nonnegative.")
    ANGLE_VALUE_WEIGHT = args.angle_value_weight
    VALLEY_DEPTH_WEIGHT = args.valley_depth_weight
    LANDMARK_WEIGHT = args.landmark_weight
    DIRECTION_WEIGHT = args.direction_weight
    ASSISTANCE_RESPONSE_WEIGHT = args.assistance_response_weight
    training_requested = bool(args.sweeps or args.neural_warmup_nfev)
    out = OUT / args.subject
    start = time.monotonic()
    data = loadmat(out / "prepared.mat", simplify_cells=True)
    if training_requested and data.get("qualityGateVersion") not in (
        "clock-smoothing-gate-v2",
        "excitation-ik-three-loads-v1",
        "excitation-ik-speed-grid-v1",
        "excitation-ik-speed-grid-torque-v1",
    ):
        raise RuntimeError(
            "Refusing calibration before timestamp/smoothness quality gate. Run quality_gate.py and prepare_subjects.py --quality-gate first."
        )
    prepared_hash = hashlib.sha256((out / "prepared.mat").read_bytes()).hexdigest()
    if args.robust_recovery:
        # Broad engineering bounds retain damping and braking authority. They
        # do not constrain the angle, prescribe an abdominal/back ratio or add
        # an external balancing torque. These parameters remain optimizable.
        LOW[3] = 2.0
        LOW[8] = 100.0
        NLOW[5] = 10.0
    NLOW[6] = args.minimum_feedback_blend
    NLOW[3] = args.minimum_motion_stiffness
    generator = data["generator"]
    if args.continuous_training:
        cycles = many(data["cycles"])
        selected = [c for c in cycles if c["split"] == "train"]
        train = [
            record(
                combine(
                    [c for c in selected if c["trialPath"] == path], args.subject, path
                ),
                generator,
                2,
            )
            for path in sorted(set(c["trialPath"] for c in selected))
        ]
    else:
        train = [record(s, generator, 2) for s in many(data["trainSequences"])]
    if args.robust_recovery:
        for initial, velocity in ((1.7, 0.6), (1.25, -0.3)):
            n = 301
            zero = np.zeros(n)
            q = np.full(n, 1.5)
            q[0] = initial
            stress = {
                "q": q,
                "torque": zero,
                "a": np.zeros((n, 2)),
                "front": zero,
                "m": zero,
                "loadKg": zero,
                "speedCode": zero,
                "stageId": zero,
                "waiting": np.ones(n),
                "referenceIsFilled": np.ones(n),
                "t": np.arange(n) * 0.02,
                "dt": 0.02,
                "standingAngle": 1.5,
                "v0": velocity,
                "segments": [],
            }
            train.append(record(stress, generator))
    if args.active_recovery_stress:
        selected = [c for c in many(data["cycles"]) if c["split"] == "train"]
        for load in sorted(set(c["loadKg"] for c in selected)):
            cycle = next(c for c in selected if c["loadKg"] == load)
            stress = combine([cycle], args.subject, "active perturbation stress")
            first = 200
            for key in (
                "q",
                "torque",
                "a",
                "front",
                "m",
                "loadKg",
                "speedCode",
                "stageId",
                "waiting",
                "referenceIsFilled",
                "t",
            ):
                stress[key] = np.asarray(stress[key])[first:].copy()
            stress["q"][0] = 1.55
            stress["v0"] = 0.6
            stress["segments"] = []
            stress["activePerturbationStress"] = True
            train.append(record(stress, generator, 2))
    full = [record(s, generator) for s in many(data["sequences"])]
    original = loadmat(
        ROOT / "data/processed/hill_compliant/neural_models.mat", simplify_cells=True
    )
    cfg = np.array(original["controller"], copy=True)
    p = np.vstack([m["parameters"] for m in original["models"]])
    r = np.vstack([m["neuralParameters"] for m in original["models"]])
    seed_p, seed_r = p.copy(), r.copy()
    if args.initial_model is not None:
        initial = loadmat(args.initial_model, simplify_cells=True)
        if initial.get("subject") != args.subject:
            raise ValueError("An explicit warm start must belong to this subject.")
        p = np.vstack([m["parameters"] for m in many(initial["models"])])
        r = np.vstack([m["neuralParameters"] for m in many(initial["models"])])
        cfg = np.array(initial["controller"], copy=True)
        previous_mode = many(initial["models"])[0].get("intentDecoder", "subtractive")
        if previous_mode != args.intent_decoder:
            # A new decoder needs a new neural seed; anatomy remains a warm
            # start and every eligible parameter is subsequently re-fitted.
            r[:] = [1.3, 0.6, 0.08, 100.0, 25.0, 20.0, 0.9, 0.2, 0.6, 0.5]
            if args.intent_decoder == "soft_tonic":
                r[:] = [1.6, 0.6, 0.06, 120.0, 25.0, 20.0, 0.9, 0.2, 0.8, 0.8]
    if args.resume and (out / "checkpoint.mat").exists():
        if training_requested:
            (out / "fit_initial_checkpoint.mat").write_bytes(
                (out / "checkpoint.mat").read_bytes()
            )
        previous = loadmat(out / "checkpoint.mat", simplify_cells=True)
        if previous.get("intentDecoder", "subtractive") != args.intent_decoder:
            raise ValueError(
                "Decoder mode changed: use an explicit initial model in a new output directory."
            )
        if training_requested and previous.get("preparedSHA256") != prepared_hash:
            raise RuntimeError(
                "A checkpoint from different/unchecked data cannot resume calibration."
            )
        p = previous["parameters"]
        r = previous["neuralParameters"]
    elif args.initial_model is None:
        if args.subject_body_search:
            # The original subject is only a numerically stable warm start for
            # muscle geometry. Body mass/lever are selected on THIS subject;
            # all source-specific and shared parameters are subsequently free.
            seed_p[:, 18] = 30.0
            seed_p[:, 19] = 0.28
            best = (np.inf, p.copy())
            for mass in (22.0, 32.0, 42.0):
                for lever in (0.23, 0.28, 0.33):
                    candidate = p.copy()
                    candidate[:, 18] = mass
                    candidate[:, 19] = lever
                    value = sum(
                        np.sum(
                            residual(
                                candidate[s],
                                r[s],
                                cfg,
                                train,
                                s,
                                seed_p[s],
                                seed_r[s],
                                args.steps,
                            )
                            ** 2
                        )
                        for s in (0, 1)
                    )
                    if value < best[0]:
                        best = (value, candidate)
            p = best[1]
            print(
                args.subject,
                "SUBJECT BODY INITIALIZATION",
                p[0, 18:20].tolist(),
                best[0],
                flush=True,
            )
        # Subjects differ in tonic activation. A small training-only decoder
        # search avoids inheriting the original subject's amplitude/exponent.
        warm = []
        seen = set()
        for rec in train:
            if not rec["segments"]:
                continue
            key = rec["segments"][0]["loadKg"]
            if key not in seen:
                seen.add(key)
                warm.append(rec)
        for source in (0, 1):
            best = (np.inf, r[source].copy())
            for amp in (0.9, 1.5, 2.3):
                for exponent in (0.6, 1.2, 1.9):
                    candidate = r[source].copy()
                    candidate[:2] = amp, exponent
                    v = residual(
                        p[source],
                        candidate,
                        cfg,
                        warm,
                        source,
                        seed_p[source],
                        seed_r[source],
                        args.steps,
                    )
                    if v @ v < best[0]:
                        best = (v @ v, candidate)
            r[source] = best[1]
    if training_requested:
        p = np.clip(p, LOW + 1e-5 * (HIGH - LOW), HIGH - 1e-5 * (HIGH - LOW))
        r = np.clip(r, NLOW + 1e-5 * (NHIGH - NLOW), NHIGH - 1e-5 * (NHIGH - NLOW))

    def loss(pp, rr):
        return sum(
            np.sum(
                residual(pp[s], rr[s], cfg, train, s, seed_p[s], seed_r[s], args.steps)
                ** 2
            )
            for s in (0, 1)
        )

    if args.decoder_multistart:
        # Clamped intent can hide the amplitude/exponent gradients entirely.
        # Search a few causal-decoder seeds jointly across ALL training cases.
        for source in (0, 1):
            best = (np.inf, r[source].copy())
            seeds = [r[source].copy()]
            for tonic in (0.0, 0.35, 0.8):
                for feedback in (0.65, 0.95):
                    for lead in (0.0, 0.25, 0.65):
                        rr = r[source].copy()
                        rr[[0, 1, 6, 8, 9]] = [2.0, 0.9, feedback, lead, tonic]
                        seeds.append(rr)
            for rr in seeds:
                rr = np.clip(
                    rr, NLOW + 1e-5 * (NHIGH - NLOW), NHIGH - 1e-5 * (NHIGH - NLOW)
                )
                v = residual(
                    p[source],
                    rr,
                    cfg,
                    train,
                    source,
                    seed_p[source],
                    seed_r[source],
                    args.steps,
                )
                value = float(v @ v)
                if value < best[0]:
                    best = value, rr.copy()
            r[source] = best[1]
            print(
                args.subject,
                "DECODER SEED",
                source + 1,
                "LOSS",
                best[0],
                "PARAMETERS",
                best[1].tolist(),
                flush=True,
            )
    current = loss(p, r)
    if training_requested and args.continuous_training and current > 100:
        candidates = []
        for blend in (1.0, 0.75, 0.5, 0.25, 0.0):
            pp = blend * p + (1 - blend) * seed_p
            rr = blend * r + (1 - blend) * seed_r
            candidates.append((loss(pp, rr), pp, rr))
        current, p, r = min(candidates, key=lambda v: v[0])
    print(
        args.subject, "INITIAL", current, "training sequences", len(train), flush=True
    )
    if args.neural_warmup_nfev:

        def warmup(source):
            low, span = NLOW, NHIGH - NLOW
            z = np.clip((r[source] - low) / span, 1e-5, 1 - 1e-5)
            best = [np.inf, z.copy()]

            def fun(z):
                v = residual(
                    p[source],
                    low + z * span,
                    cfg,
                    train,
                    source,
                    seed_p[source],
                    seed_r[source],
                    args.steps,
                )
                if v @ v < best[0]:
                    best[:] = [v @ v, z.copy()]
                return v

            least_squares(
                fun,
                z,
                jac=(
                    absolute_jacobian(fun, args.absolute_difference_step)
                    if args.absolute_difference_step
                    else "2-point"
                ),
                bounds=(np.zeros(len(z)), np.ones(len(z))),
                max_nfev=args.neural_warmup_nfev,
                diff_step=2e-5,
                ftol=5e-4,
            )
            print(
                args.subject, "NEURAL WARMUP", source + 1, "LOSS", best[0], flush=True
            )
            return low + best[1] * span

        with ThreadPoolExecutor(max_workers=2) as pool:
            r = np.vstack(list(pool.map(warmup, (0, 1))))
        current = loss(p, r)
    for sweep in range(args.sweeps):

        def fit_source(source):
            pp = p[source].copy()
            rr = r[source].copy()
            low = np.r_[LOW[SPECIFIC], NLOW]
            high = np.r_[HIGH[SPECIFIC], NHIGH]
            z = np.clip((np.r_[pp[SPECIFIC], rr] - low) / (high - low), 1e-5, 1 - 1e-5)
            best = [np.inf, z.copy()]
            last = [time.monotonic()]
            calls = [0]

            def fun(z):
                physical = pp.copy()
                neural = (
                    low[len(SPECIFIC) :]
                    + z[len(SPECIFIC) :] * (high - low)[len(SPECIFIC) :]
                )
                physical[SPECIFIC] = (
                    low[: len(SPECIFIC)]
                    + z[: len(SPECIFIC)] * (high - low)[: len(SPECIFIC)]
                )
                v = residual(
                    physical,
                    neural,
                    cfg,
                    train,
                    source,
                    seed_p[source],
                    seed_r[source],
                    args.steps,
                )
                calls[0] += 1
                if v @ v < best[0]:
                    best[:] = [v @ v, z.copy()]
                if time.monotonic() - last[0] > 30:
                    print(
                        args.subject,
                        "SWEEP",
                        sweep + 1,
                        "SOURCE",
                        source + 1,
                        "CALLS",
                        calls[0],
                        "LOSS",
                        round(best[0], 4),
                        flush=True,
                    )
                    last[0] = time.monotonic()
                return v

            result = least_squares(
                fun,
                z,
                jac=(
                    absolute_jacobian(fun, args.absolute_difference_step)
                    if args.absolute_difference_step
                    else "2-point"
                ),
                bounds=(np.zeros(len(z)), np.ones(len(z))),
                max_nfev=args.max_nfev,
                diff_step=2e-5,
                ftol=5e-4,
            )
            z = best[1]
            pp[SPECIFIC] = (
                low[: len(SPECIFIC)]
                + z[: len(SPECIFIC)] * (high - low)[: len(SPECIFIC)]
            )
            rr = (
                low[len(SPECIFIC) :]
                + z[len(SPECIFIC) :] * (high - low)[len(SPECIFIC) :]
            )
            return pp, rr

        with ThreadPoolExecutor(max_workers=2) as pool:
            fitted = list(pool.map(fit_source, (0, 1)))
        p = np.vstack([f[0] for f in fitted])
        r = np.vstack([f[1] for f in fitted])
        z = np.clip(
            (p[0, SHARED] - LOW[SHARED]) / (HIGH[SHARED] - LOW[SHARED]), 1e-5, 1 - 1e-5
        )

        def common(z):
            pp = p.copy()
            pp[:, SHARED] = LOW[SHARED] + z * (HIGH[SHARED] - LOW[SHARED])
            return np.concatenate(
                [
                    residual(
                        pp[s], r[s], cfg, train, s, seed_p[s], seed_r[s], args.steps
                    )
                    for s in (0, 1)
                ]
            )

        result = least_squares(
            common,
            z,
            jac=(
                absolute_jacobian(common, args.absolute_difference_step)
                if args.absolute_difference_step
                else "2-point"
            ),
            bounds=(np.zeros(len(z)), np.ones(len(z))),
            max_nfev=args.max_nfev,
            diff_step=2e-5,
            ftol=5e-4,
        )
        p[:, SHARED] = LOW[SHARED] + result.x * (HIGH[SHARED] - LOW[SHARED])
        current = loss(p, r)
        savemat(
            out / "checkpoint.mat",
            {
                "parameters": p,
                "neuralParameters": r,
                "intentDecoder": args.intent_decoder,
                "controller": cfg,
                "loss": current,
                "preparedSHA256": prepared_hash,
            },
            long_field_names=True,
        )
        print(
            args.subject,
            "SWEEP FINISHED",
            sweep + 1,
            "LOSS",
            current,
            "SECONDS",
            time.monotonic() - start,
            flush=True,
        )
    models = [
        {
            "parameters": p[s],
            "neuralParameters": r[s],
            "intentDecoder": args.intent_decoder,
            "generator": generator,
            "subject": args.subject,
            "objective": "Training-only personalized calibration of the same neural/Hill system",
        }
        for s in (0, 1)
    ]
    savemat(
        out / "neural_models.mat",
        {
            "models": structs(models),
            "parameters": p,
            "neuralParameters": r,
            "controller": cfg,
            "sharedParameterIndices": SHARED + 1,
            "mapModel": data["mapModel"],
            "subject": args.subject,
        },
        do_compression=True,
        long_field_names=True,
    )
    baseline = [
        {"parameters": seed_p[s], "neuralParameters": seed_r[s]} for s in (0, 1)
    ]
    metrics = score_rows(baseline, cfg, full, "original_subject_transfer") + score_rows(
        models, cfg, full
    )
    (out / "bse_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    # Display sequence uses one actually recorded held-out action per available
    # load. Absent 10/15 kg recordings are not synthesized as validation data.
    cycles = many(data["cycles"])
    display = []
    for load in (5, 10, 15):
        choices = [c for c in cycles if c["loadKg"] == load]
        if choices:
            complete = [
                c for c in choices if not c.get("rightCensored", False)
            ] or choices
            display.append(
                next(
                    (c for c in complete if str(c["split"]).startswith("test")),
                    complete[0],
                )
            )
    sequences = many(data["sequences"]) + [
        combine(display, args.subject, "continuous_available_loads")
    ]
    export_replay(models, cfg, sequences, generator, out)
    policy = {
        "subject": args.subject,
        "intentDecoder": args.intent_decoder,
        "loss": float(current),
        "wallSeconds": time.monotonic() - start,
        "runSweeps": args.sweeps,
        "resumedCheckpoint": bool(args.resume),
        "resumeInitialCheckpointSHA256": (
            hashlib.sha256(
                (out / "fit_initial_checkpoint.mat").read_bytes()
            ).hexdigest()
            if args.resume and (out / "fit_initial_checkpoint.mat").exists()
            else None
        ),
        "initialModel": (
            str(args.initial_model) if args.initial_model is not None else None
        ),
        "initialModelSHA256": (
            hashlib.sha256(args.initial_model.read_bytes()).hexdigest()
            if args.initial_model is not None
            else None
        ),
        "continuousTraining": bool(args.continuous_training),
        "robustRecoveryCalibration": bool(args.robust_recovery),
        "activePerturbationStress": bool(args.active_recovery_stress),
        "minimumFeedbackBlend": args.minimum_feedback_blend,
        "minimumMotionStiffnessNmPerRad": args.minimum_motion_stiffness,
        "absoluteNormalizedDifferenceStep": args.absolute_difference_step,
        "parallelIndependentRolloutWorkers": args.rollout_workers,
        "shapeTailWeight": args.shape_tail_weight,
        "angleValueWeight": args.angle_value_weight,
        "valleyDepthWeight": args.valley_depth_weight,
        "landmarkWeight": args.landmark_weight,
        "directionWeight": args.direction_weight,
        "assistanceResponseWeight": args.assistance_response_weight,
        "assistanceResponsePrior": dict(
            representativeAssistNm=60.0,
            nominalAngleDeflectionRad=0.12,
            independentlyMeasured=False,
        ),
        "decoderMultistart": args.decoder_multistart,
        "maximumCausalEMGLeadGainSeconds": args.maximum_emg_lead,
        "maximumFlexionRequestNm": args.maximum_flexion_request,
        "neuralWarmupMaxNfev": args.neural_warmup_nfev,
        "preparedSHA256": hashlib.sha256(
            (out / "prepared.mat").read_bytes()
        ).hexdigest(),
        "fixedVectorPerSource": True,
        "sharedBodyGeometry": True,
        "all42ParametersEligible": True,
        "qualityGateVersion": data.get(
            "qualityGateVersion", "historical_unchecked_replay"
        ),
        "qualityGateSHA256": data.get("qualityGateSHA256", ""),
        "subjectBodyInitializationSearch": bool(args.subject_body_search),
        "optimizedSourceSpecificPhysicalIndices": (SPECIFIC + 1).tolist(),
        "optimizedSharedPhysicalIndices": (SHARED + 1).tolist(),
        "optimizedNeuralParameterIndices": list(range(1, 11)),
        "parameterChangeFromOriginal": [
            np.r_[p[s] - original["models"][s]["parameters"], r[s] - seed_r[s]].tolist()
            for s in (0, 1)
        ],
        "unexcitedSoftBoundaryParametersNotIdentifiable": True,
        "calibrationUsesTestRecords": False,
        "calibrationIntegrationSeconds": 0.02 / args.steps,
        "replayIntegrationSeconds": 0.001,
        "referenceAnglePolarityChanged": False,
        "modelSHA256": hashlib.sha256(
            (out / "neural_models.mat").read_bytes()
        ).hexdigest(),
        "displayGroup": len(sequences),
        "groups": len(sequences),
        "availableLoadsKg": [c["loadKg"] for c in display],
    }
    (out / "fit_policy.json").write_text(json.dumps(policy, indent=2) + "\n")
    print("COMPLETED", json.dumps(policy), flush=True)
    if ROLLOUT_POOL is not None:
        ROLLOUT_POOL.shutdown(wait=True)


if __name__ == "__main__":
    main()
