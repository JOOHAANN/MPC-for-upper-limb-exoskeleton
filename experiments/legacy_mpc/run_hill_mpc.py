"""MPC on the selected compliant Hill pair, with a causal human-model preview.

Measured kinematics are evaluation-only. The MPC preserves the unassisted
motion predicted by the same EMG-driven model, while reducing excess back
torque. Every source retains its one calibrated parameter vector in all cases.
"""

from __future__ import annotations

import argparse
import ctypes
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat
from scipy.optimize import Bounds, LinearConstraint, minimize

ROOT = Path(__file__).resolve().parents[2]
FOLDER = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "experiments/hill_compliant"))
from native_library import load_library

LABELS = ("Experimental activation", "MAP activation")
OUTPUT_NAMES = (
    "q",
    "v",
    "torque",
    "backTorque",
    "frontTorque",
    "backActive",
    "backPassive",
    "backFiberDamping",
    "frontActive",
    "frontPassive",
    "frontFiberDamping",
    "backLengthNormalized",
    "frontLengthNormalized",
    "backActivation",
    "frontActivation",
    "equilibriumResidual",
    "momentumResidual",
    "backCommand",
    "frontCommand",
    "balanceWeight",
    "requestedHillTorque",
    "postureCorrectionRequest",
    "postureIntegral",
    "motionIntent",
    "assistanceTorque",
    "boundaryTorque",
    "gravityTorque",
)


class Plant:
    def __init__(self, model, controller, angle):
        mode = model.get("intentDecoder", "subtractive")
        if mode not in ("subtractive", "divisive", "soft_tonic"):
            raise ValueError("Unknown motion intent decoder: " + str(mode))
        source = (
            FOLDER
            / {
                "subtractive": "hill_mpc_state.c",
                "divisive": "hill_mpc_divisive_state.c",
                "soft_tonic": "hill_mpc_soft_tonic_state.c",
            }[mode]
        )
        dependencies = list(
            dict.fromkeys(
                [
                    source,
                    FOLDER / "hill_mpc_state.c",
                    *sorted((ROOT / "experiments/hill_compliant").glob("hill_*core.c")),
                ]
            )
        )
        self.library = load_library(source, dependencies, ROOT)
        array = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
        self.library.hill_mpc_init.argtypes = (
            [array] * 4 + [ctypes.c_double] * 2 + [array] * 3
        )
        self.library.hill_mpc_advance.argtypes = (
            [array] * 5 + [ctypes.c_double] * 4 + [ctypes.c_int, array]
        )
        self.library.hill_mpc_predict.argtypes = (
            [array] * 5
            + [ctypes.c_int] * 2
            + [ctypes.c_double, ctypes.c_int, ctypes.c_double, array]
        )
        self.library.hill_mpc_predict_input_trend.argtypes = (
            [array] * 6
            + [ctypes.c_int] * 2
            + [ctypes.c_double, ctypes.c_int, ctypes.c_double, ctypes.c_double, array]
        )
        self.library.hill_mpc_predict_map.argtypes = (
            [array] * 8
            + [ctypes.c_int] * 2
            + [ctypes.c_double, ctypes.c_int, ctypes.c_double, array]
        )
        self.p = np.ascontiguousarray(model["parameters"], dtype=float)
        self.r = np.ascontiguousarray(model["neuralParameters"], dtype=float)
        self.cfg = np.ascontiguousarray(np.r_[angle, controller], dtype=float)
        self.input_slope = np.zeros(2)

    def init(self, a, mass, wait, initial):
        state, out = np.zeros(34), np.zeros(27)
        failed = self.library.hill_mpc_init(
            self.p,
            self.r,
            self.cfg,
            np.ascontiguousarray(a),
            float(mass),
            float(wait),
            np.ascontiguousarray(initial),
            state,
            out,
        )
        assert not failed
        return state, out

    def step(self, state, a, mass, wait, u, dt=0.01, steps=10):
        out = np.zeros(27)
        failed = self.library.hill_mpc_advance(
            self.p,
            self.r,
            self.cfg,
            state,
            np.ascontiguousarray(a),
            float(mass),
            float(wait),
            float(u),
            float(dt),
            int(steps),
            out,
        )
        if failed:
            raise RuntimeError(
                "The coupled Hill solve failed; no angle clipping or tail replacement is allowed."
            )
        return out

    def preview(self, state, controls, options):
        out = np.zeros((options["points"], 27))
        if options.get("inputForecast", "hold") in ("map", "arx"):
            failed = self.library.hill_mpc_predict_map(
                self.p,
                self.r,
                self.cfg,
                np.ascontiguousarray(state),
                np.ascontiguousarray(controls),
                self.map_activation,
                self.map_mass,
                self.map_wait,
                len(controls),
                options["points"],
                options["predictionDt"],
                options["predictionSubsteps"],
                options["rateNmPerSec"],
                out,
            )
        elif options.get("inputForecast", "hold") == "trend":
            failed = self.library.hill_mpc_predict_input_trend(
                self.p,
                self.r,
                self.cfg,
                np.ascontiguousarray(state),
                np.ascontiguousarray(controls),
                self.input_slope,
                len(controls),
                options["points"],
                options["predictionDt"],
                options["predictionSubsteps"],
                options["rateNmPerSec"],
                options["inputForecastTime"],
                out,
            )
        else:
            failed = self.library.hill_mpc_predict(
                self.p,
                self.r,
                self.cfg,
                np.ascontiguousarray(state),
                np.ascontiguousarray(controls),
                len(controls),
                options["points"],
                options["predictionDt"],
                options["predictionSubsteps"],
                options["rateNmPerSec"],
                out,
            )
        return None if failed else out


def configure_map_preview(plant, seq, observed_inputs, current, options):
    """Only past observations correct a trained, protocol-phase MAP forecast."""
    planned = np.asarray(seq["plannedMAPInputs"], dtype=float)
    if planned.shape != np.asarray(observed_inputs).shape:
        raise ValueError("MAP activation prediction must match the protocol length")
    points = options["points"]
    elapsed = np.arange(1, points + 1) * options["predictionDt"]
    query = current + elapsed / seq["dt"]
    index = np.arange(len(planned))
    gain = np.ones(2)
    phase = 0.0
    if options.get("mapPhaseObserver", False) and not seq.get(
        "activationSourceIsGenerated", False
    ):
        start = max(0, current - round(0.8 / seq["dt"]))
        history = np.arange(start, current + 1)
        observed = np.asarray(observed_inputs)[history]
        if len(history) >= 12 and np.ptp(observed[:, 0]) > 0.025:
            best = (np.inf, 0.0, np.ones(2))
            for shift in np.arange(-0.30, 0.301, 0.03):
                candidate = np.column_stack(
                    [
                        np.interp(history + shift / seq["dt"], index, planned[:, j])
                        for j in (0, 1)
                    ]
                )
                centered = candidate - candidate.mean(0)
                amplitude = np.clip(
                    np.sum(centered * (observed - observed.mean(0)), axis=0)
                    / np.maximum(np.sum(centered**2, axis=0), 1e-6),
                    0.25,
                    3.0,
                )
                mismatch = (observed - observed.mean(0)) - amplitude * centered
                loss = float(
                    np.mean(mismatch[:, 0] ** 2)
                    + 0.20 * np.mean(mismatch[:, 1] ** 2)
                    + 2e-4 * shift**2
                )
                if loss < best[0]:
                    best = (loss, shift, amplitude)
            phase = float(best[1])
            gain = best[2]
            # A smooth observer avoids discontinuous changes of predicted phase.
            previous = getattr(plant, "map_phase_shift", 0.0)
            phase = previous + np.clip(phase - previous, -0.03, 0.03)
            plant.map_phase_shift = phase
    query = query + phase / seq["dt"]
    forecast = np.column_stack([np.interp(query, index, planned[:, j]) for j in (0, 1)])
    if options.get("inputForecast") == "arx":
        predictor = options["activationResidualPredictor"]
        lags = np.asarray(predictor["lags"], dtype=int)
        past = np.maximum(0, current - lags)
        residual = (np.asarray(observed_inputs)[past] - planned[past]).reshape(-1)
        residual_forecast = (residual / np.asarray(predictor["scale"])) @ np.asarray(
            predictor["coefficients"]
        )
        forecast += residual_forecast.reshape(points, 2)
        wait = np.interp(query, index, seq["waiting"])
        forecast *= 1 - wait[:, None]
        plant.map_activation = np.ascontiguousarray(np.clip(forecast, 0, 1))
        plant.map_mass = np.ascontiguousarray(np.interp(query, index, seq["m"]))
        plant.map_wait = np.ascontiguousarray(wait)
        return
    planned_current = np.array(
        [np.interp(current + phase / seq["dt"], index, planned[:, j]) for j in (0, 1)]
    )
    forecast = forecast * gain
    correction = observed_inputs[current] - planned_current * gain
    tau = options.get("mapInnovationDecaySeconds", 0.25)
    forecast += np.exp(-elapsed[:, None] / tau) * correction
    plant.map_activation = np.ascontiguousarray(np.clip(forecast, 0, 1))
    # These are declared protocol inputs, not measured future kinematics/forces.
    plant.map_mass = np.ascontiguousarray(np.interp(query, index, seq["m"]))
    plant.map_wait = np.ascontiguousarray(np.interp(query, index, seq["waiting"]))


def residual(
    preview,
    controls,
    current,
    nominal,
    threshold,
    q_weight,
    back_weight=2.0,
    velocity_weight=0.0,
    options=None,
):
    options = {} if options is None else options
    values = np.r_[
        np.sqrt(q_weight) * (preview[:, 0] - nominal[:, 0]) / 0.10,
        np.sqrt(back_weight) * np.maximum(preview[:, 3] - threshold, 0) / 30.0,
        np.maximum(preview[:, 4] - 12.0, 0) / 15.0,
        np.sqrt(options.get("assistanceEffortWeight", 0.01)) * controls / 80.0,
        0.20 * np.diff(np.r_[current, controls]) / 20.0,
    ]
    under_weight = options.get("backTorqueUnderWeight", 0.0)
    if under_weight > 0:
        # Only penalize over-assistance when the unassisted model predicts
        # that help is needed. Low-demand portions are never lifted to a line.
        target = np.minimum(
            nominal[:, 3], threshold - options.get("trackingHeadroomNm", 3.0)
        )
        demand = np.clip((nominal[:, 3] - threshold) / 10.0, 0.0, 1.0)
        values = np.r_[
            values,
            np.sqrt(under_weight * demand)
            * np.maximum(target - preview[:, 3], 0.0)
            / 10.0,
        ]
    if velocity_weight > 0:
        values = np.r_[
            values, np.sqrt(velocity_weight) * (preview[:, 1] - nominal[:, 1])
        ]
    return values


def mpc_command(plant, state, nominal_state, warm, threshold, options):
    size, current = len(warm), state[28]
    priority = options.get("torqueConstraintPriority", False)
    cap = options.get("capNm", 80.0)
    # Reserve torque headroom for muscle activation lag and prediction errors.
    threshold = max(0.0, threshold - options.get("backTorqueMarginNm", 0.0))
    nominal = plant.preview(nominal_state, np.zeros(size), options)
    if nominal is None:
        return (
            max(0.0, current - options["rateNmPerSec"] * options["controlPeriod"]),
            -3,
            0.0,
        )
    guard = np.asarray(
        options.get("predictionErrorGuardNm", np.zeros(options["points"])), dtype=float
    )
    relative_guard = options.get("predictionErrorGuardRelative")
    if relative_guard:
        # One normalized uncertainty policy applies to every declared limit.
        guard = 2.0 + np.asarray(relative_guard) * threshold / 0.75
    if "predictionErrorGuardScale" in options:
        floor = options.get("predictionErrorGuardFloorNm", 2.0)
        guard = floor + options["predictionErrorGuardScale"] * np.maximum(
            guard - floor, 0.0
        )
    if options.get("predictionErrorGuardCapNm", 0.0) > 0:
        guard = np.minimum(guard, options["predictionErrorGuardCapNm"])
    ceiling = threshold - guard
    constraint_points = min(
        options["points"],
        max(
            1,
            round(
                options.get(
                    "torqueConstraintHorizonSeconds",
                    options["points"] * options["predictionDt"],
                )
                / options["predictionDt"]
            ),
        ),
    )
    constraint_ceiling = ceiling[:constraint_points]
    if (
        np.max(nominal[:, 3] - ceiling) < -5
        and current < 1e-8
        and (not priority or np.max(abs(state[:2] - nominal_state[:2])) < 1e-7)
    ):
        return 0.0, 1, 0.0
    limit = options["rateNmPerSec"] * options["controlPeriod"]
    warm = np.clip(warm, 0.0, cap)
    warm[0] = np.clip(warm[0], max(0.0, current - limit), min(cap, current + limit))
    block_limit = (
        options["rateNmPerSec"] * options["predictionDt"] * options["points"] / size
    )
    for i in range(1, size):
        warm[i] = np.clip(
            warm[i],
            max(0.0, warm[i - 1] - block_limit),
            min(cap, warm[i - 1] + block_limit),
        )
    base = plant.preview(state, warm, options)
    if base is None:
        return max(0.0, current - limit), -3, 0.0
    r0 = residual(
        base,
        warm,
        current,
        nominal,
        threshold,
        options["angleWeight"],
        options["backTorqueWeight"],
        options.get("velocityWeight", 0.0),
        options,
    )
    jacobian = np.zeros((len(r0), size))
    torque_jacobian = np.zeros((options["points"], size))
    for i in range(size):
        trial = warm.copy()
        delta = 0.5 if trial[i] < cap - 0.5 else -0.5
        trial[i] += delta
        predicted = plant.preview(state, trial, options)
        if predicted is not None:
            torque_jacobian[:, i] = (predicted[:, 3] - base[:, 3]) / delta
            jacobian[:, i] = (
                residual(
                    predicted,
                    trial,
                    current,
                    nominal,
                    threshold,
                    options["angleWeight"],
                    options["backTorqueWeight"],
                    options.get("velocityWeight", 0.0),
                    options,
                )
                - r0
            ) / delta
    hessian = jacobian.T @ jacobian + 1e-8 * np.eye(size)
    gradient = jacobian.T @ r0
    difference = np.diff(np.eye(size), axis=0)
    low, high = np.zeros(size), np.full(size, cap)
    low[0], high[0] = max(0.0, current - limit), min(cap, current + limit)
    constraints = [LinearConstraint(difference, -block_limit, block_limit)]
    if priority:
        # Phase I minimizes required torque slack with amplitude/slew limits
        # unchanged. Phase II preserves motion inside that attainable slack.
        # The slack is explicit: unknown future EMG can still violate a target.
        torque_jacobian = torque_jacobian[:constraint_points]
        matrix = np.c_[torque_jacobian / 30.0, -np.ones(constraint_points)]
        rhs = (
            constraint_ceiling - base[:constraint_points, 3] + torque_jacobian @ warm
        ) / 30.0
        diff_augmented = np.c_[difference, np.zeros(size - 1)]
        phase_one = minimize(
            lambda z: z[-1] ** 2 + 1e-10 * np.sum(z[:-1] ** 2),
            np.r_[
                warm,
                max(0.0, np.max(base[:constraint_points, 3] - constraint_ceiling))
                / 30.0,
            ],
            jac=lambda z: np.r_[2e-10 * z[:-1], 2 * z[-1]],
            method="SLSQP",
            bounds=Bounds(np.r_[low, 0.0], np.r_[high, np.inf]),
            constraints=[
                LinearConstraint(diff_augmented, -block_limit, block_limit),
                LinearConstraint(matrix, -np.inf, rhs),
            ],
            options={"maxiter": 60, "ftol": 1e-9},
        )
        if phase_one.success and np.isfinite(phase_one.x).all():
            slack = max(0.0, phase_one.x[-1]) + 1e-5
            constraints.append(
                LinearConstraint(torque_jacobian / 30.0, -np.inf, rhs + slack)
            )
            initial = phase_one.x[:-1]
        else:
            return float(warm[0]), -2, float(r0 @ r0)
    else:
        initial = warm
    result = minimize(
        lambda u: 0.5 * (u - warm) @ hessian @ (u - warm) + gradient @ (u - warm),
        initial,
        jac=lambda u: hessian @ (u - warm) + gradient,
        method="SLSQP",
        bounds=Bounds(low, high),
        constraints=constraints,
        options={"maxiter": 35, "ftol": 1e-7},
    )
    # Reject an infeasible/non-finite optimizer result; apply a bounded fallback.
    feasible = (
        result.success
        and np.isfinite(result.x).all()
        and np.all(result.x >= low - 1e-7)
        and np.all(result.x <= high + 1e-7)
        and np.all(np.abs(difference @ result.x) <= block_limit + 1e-7)
    )
    if priority:
        feasible = feasible and np.all(
            torque_jacobian @ result.x <= 30.0 * (rhs + slack) + 1e-5
        )
    if not feasible:
        if priority:
            # Phase I already supplied a bounded feasible command when motion
            # optimization cannot improve it. This is a reported fallback.
            result.x = initial
        else:
            return warm[0], -2, float(r0 @ r0)
    best, value = warm.copy(), float(r0 @ r0)
    violation = float(
        np.max(np.maximum(base[:constraint_points, 3] - constraint_ceiling, 0.0))
    )
    for fraction in ((1.0, 0.5, 0.25) if priority else (1.0, 0.5)):
        candidate = warm + fraction * (result.x - warm)
        predicted = plant.preview(state, candidate, options)
        if predicted is None:
            continue
        res = residual(
            predicted,
            candidate,
            current,
            nominal,
            threshold,
            options["angleWeight"],
            options["backTorqueWeight"],
            options.get("velocityWeight", 0.0),
            options,
        )
        excess = float(
            np.max(
                np.maximum(predicted[:constraint_points, 3] - constraint_ceiling, 0.0)
            )
        )
        better = (
            (
                excess < violation - 1e-3
                or (excess <= violation + 1e-3 and res @ res < value)
            )
            if priority
            else res @ res < value
        )
        if better:
            best, value = candidate, float(res @ res)
            violation = excess
    return float(np.clip(best[0], low[0], high[0])), (1 if feasible else -2), value


def prediction_dict(out, p):
    d = {name: out[:, i] for i, name in enumerate(OUTPUT_NAMES)}
    d["effectiveActivation"] = out[:, 13:15]
    d["dampingTorque"] = -p[3] * out[:, 1]
    d["externalBalanceTorque"] = np.zeros(len(out))
    return d


def simulate_case(model, controller, seq, inputs, options, progress=None):
    """Run the same causal controller for a complete, unreset sequence."""
    plant = Plant(model, controller, seq["standingAngle"])
    inputs = np.ascontiguousarray(inputs, dtype=float)
    count, dt = len(inputs), float(seq["dt"])
    # A task may supply an offline-calibrated constant limit per action. The
    # controller receives only this declared limit, never a reference trajectory.
    threshold_schedule = seq.get("backTorqueThresholdNm")
    if threshold_schedule is not None:
        threshold_schedule = np.asarray(threshold_schedule, dtype=float)
        if threshold_schedule.ndim == 0:
            threshold_schedule = np.full(count, float(threshold_schedule))
        if (
            threshold_schedule.shape != (count,)
            or not np.isfinite(threshold_schedule).all()
            or np.any(threshold_schedule <= 0)
        ):
            raise ValueError(
                "Back torque threshold must be positive, finite and match the sequence length"
            )
    nominal_state, first = plant.init(
        inputs[0], seq["m"][0], seq["waiting"][0], np.array([seq["q"][0], seq["v0"]])
    )
    state = nominal_state.copy()
    actual, nominal = np.zeros((count, 27)), np.zeros((count, 27))
    actual[0] = nominal[0] = first
    flags, timings = np.full(count, np.nan), np.zeros(count)
    warm, command = np.zeros(options["blocks"]), 0.0
    stride = round(options["controlPeriod"] / dt)
    for k in range(1, count):
        if (k - 1) % stride == 0:
            tick = time.perf_counter()
            history_start = max(0, k - 1 - round(options["inputSlopeWindow"] / dt))
            elapsed = (k - 1 - history_start) * dt
            plant.input_slope[:] = (
                np.clip(
                    (inputs[k - 1] - inputs[history_start]) / elapsed,
                    -options["inputSlopeLimit"],
                    options["inputSlopeLimit"],
                )
                if elapsed > 0
                else 0.0
            )
            if options.get("inputForecast") in ("map", "arx"):
                configure_map_preview(plant, seq, inputs, k - 1, options)
            released = (
                options.get("earlyRecovery", False)
                and state[32] < options["recoveryContactToleranceKg"]
                and state[10]
                > seq["standingAngle"] - options["recoveryIntentToleranceRad"]
                and state[0] > seq["standingAngle"] - 0.3
                and plant.input_slope[0] <= 0.02
            )
            if seq["waiting"][k - 1] or released:
                unsettled = (
                    abs(state[0] - seq["standingAngle"])
                    > options["recoveryAngleToleranceRad"]
                    or abs(state[1]) > options["recoveryVelocityToleranceRadPerSec"]
                )
                if options["recoveryBackTargetNm"] > 0 and unsettled:
                    recovery_options = {
                        **options,
                        "backTorqueMarginNm": 0.0,
                        "backTorqueWeight": options["recoveryBackWeight"],
                        "velocityWeight": options["recoveryVelocityWeight"],
                        "torqueConstraintPriority": False,
                        "predictionErrorGuardNm": np.zeros(options["points"]).tolist(),
                        "predictionErrorGuardRelative": [],
                        "predictionErrorGuardScale": 0.0,
                        "predictionErrorGuardFloorNm": 0.0,
                        "backTorqueUnderWeight": 0.0,
                        "assistanceEffortWeight": 0.01,
                    }
                    command, flags[k], _ = mpc_command(
                        plant,
                        state,
                        nominal_state,
                        warm,
                        options["recoveryBackTargetNm"],
                        recovery_options,
                    )
                else:
                    command, flags[k] = 0.0, 1
            else:
                threshold = (
                    float(threshold_schedule[k - 1])
                    if threshold_schedule is not None
                    else float(
                        np.interp(
                            seq["loadKg"][k - 1], [0, 5, 10, 15], [90, 100, 140, 180]
                        )
                    )
                )
                command, flags[k], _ = mpc_command(
                    plant, state, nominal_state, warm, threshold, options
                )
            upper = options.get("angleGovernorUpperRad", 0.0)
            if upper > 0 and (command > 1e-6 or state[28] > 1e-6):
                # A command governor reduces positive assistance when the SAME
                # causal Hill predictor forecasts excessive backward extension.
                # It never clips/overwrites q or adds an external joint torque.
                for _ in range(8):
                    forecast = plant.preview(
                        state, np.full(options["blocks"], command), options
                    )
                    if forecast is not None and np.max(forecast[:, 0]) <= upper:
                        break
                    if command <= 1e-6:
                        break
                    command = 0.0 if command < 0.1 else 0.5 * command
            warm[:] = command
            timings[k] = time.perf_counter() - tick
        applied = float(
            np.clip(
                command,
                max(0.0, state[28] - options["rateNmPerSec"] * dt),
                min(options["capNm"], state[28] + options["rateNmPerSec"] * dt),
            )
        )
        actual[k] = plant.step(
            state,
            inputs[k],
            seq["m"][k],
            seq["waiting"][k],
            applied,
            dt,
            options["actualSubsteps"],
        )
        nominal[k] = plant.step(
            nominal_state,
            inputs[k],
            seq["m"][k],
            seq["waiting"][k],
            0.0,
            dt,
            options["actualSubsteps"],
        )
        if progress is not None:
            progress(k, count)
    return actual, nominal, flags, timings


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--groups", type=int, nargs="+", default=None)
    parser.add_argument("--angle-weight", type=float, default=10.0)
    parser.add_argument(
        "--back-weight",
        type=float,
        default=128.0,
        help="Quadratic penalty on back torque above the internal target.",
    )
    parser.add_argument(
        "--back-margin",
        type=float,
        default=35.0,
        help="Internal torque headroom below the displayed soft threshold, in N m.",
    )
    parser.add_argument(
        "--control-period", type=float, choices=[0.01, 0.02, 0.05, 0.1], default=0.02
    )
    parser.add_argument("--input-forecast", choices=["hold", "trend"], default="hold")
    parser.add_argument(
        "--velocity-weight",
        type=float,
        default=100.0,
        help="Preserve the unassisted model's angular-velocity trend.",
    )
    parser.add_argument("--output", default="hill_mpc")
    parser.add_argument(
        "--model-file",
        type=Path,
        default=ROOT / "data/processed/hill_compliant/neural_models.mat",
    )
    parser.add_argument(
        "--replay-file",
        type=Path,
        default=ROOT / "results/hill_neural_sequence/hill_sequence.mat",
    )
    parser.add_argument(
        "--recovery-back-target",
        type=float,
        default=12.0,
        help="Positive value enables state-aware post-movement support at this back-torque target (N m).",
    )
    parser.add_argument("--recovery-back-weight", type=float, default=512.0)
    parser.add_argument("--recovery-velocity-weight", type=float, default=10.0)
    parser.add_argument(
        "--recovery-early",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Permit recovery once the known load is released and the causal intent returns upright.",
    )
    args = parser.parse_args()
    if (
        not np.isfinite(
            [
                args.angle_weight,
                args.back_weight,
                args.back_margin,
                args.velocity_weight,
            ]
        ).all()
        or args.angle_weight <= 0
        or args.back_weight <= 0
        or args.back_margin < 0
        or args.velocity_weight < 0
    ):
        parser.error(
            "Angle/back weights must be finite and positive; velocity weight and torque headroom must be finite and nonnegative."
        )
    if (
        not np.isfinite(
            [
                args.recovery_back_target,
                args.recovery_back_weight,
                args.recovery_velocity_weight,
            ]
        ).all()
        or args.recovery_back_target < 0
        or args.recovery_back_weight <= 0
        or args.recovery_velocity_weight < 0
    ):
        parser.error(
            "Recovery target/velocity weight must be nonnegative and recovery back weight positive."
        )
    target = ROOT / "results" / args.output
    target.mkdir(parents=True, exist_ok=True)
    cache_path = args.model_file
    cache = loadmat(cache_path, simplify_cells=True)
    replay = loadmat(args.replay_file, simplify_cells=True)["result"]
    options = {
        "points": 16,
        "blocks": 4,
        "predictionDt": 0.025,
        "predictionSubsteps": 10,
        "controlPeriod": args.control_period,
        "actualSubsteps": 10,
        "rateNmPerSec": 200.0,
        "capNm": 80.0,
        "angleWeight": args.angle_weight,
        "backTorqueWeight": args.back_weight,
        "backTorqueMarginNm": args.back_margin,
        "inputForecast": args.input_forecast,
        "inputForecastTime": 0.10,
        "inputSlopeWindow": 0.05,
        "inputSlopeLimit": 3.0,
        "velocityWeight": args.velocity_weight,
        "recoveryBackTargetNm": args.recovery_back_target,
        "recoveryBackWeight": args.recovery_back_weight,
        "recoveryVelocityWeight": args.recovery_velocity_weight,
        "earlyRecovery": args.recovery_early,
        "recoveryContactToleranceKg": 0.5,
        "recoveryIntentToleranceRad": 0.10,
        "recoveryAngleToleranceRad": 0.025,
        "recoveryVelocityToleranceRadPerSec": 0.04,
    }
    result, summaries, checks = (
        {"groups": [], "models": cache["models"], "options": options},
        [],
        [],
    )
    started, notified = time.perf_counter(), time.perf_counter()
    for group in args.groups or range(1, len(np.atleast_1d(replay["groups"])) + 1):
        g = replay["groups"][group - 1]
        seq = g["sequence"]
        predictions, nominals = [], []
        for source in range(2):
            model = cache["models"][source]
            assert np.array_equal(
                model["parameters"], replay["models"][source]["parameters"]
            )
            assert np.array_equal(
                model["neuralParameters"], replay["models"][source]["neuralParameters"]
            )
            inputs = np.ascontiguousarray(g["inputs"][source], dtype=float)

            def progress(k, count):
                nonlocal notified
                if time.perf_counter() - notified > 15:
                    print(
                        f"Group {group}, {LABELS[source]}: {100*k/count:.1f}%; wall {time.perf_counter()-started:.1f} s",
                        flush=True,
                    )
                    notified = time.perf_counter()

            actual, nominal, flags, timings = simulate_case(
                model, cache["controller"], seq, inputs, options, progress
            )
            baseline = g["predictions"][source]
            qdiff = float(np.max(np.abs(nominal[:, 0] - baseline["q"])))
            tdiff = float(np.max(np.abs(nominal[:, 2] - baseline["torque"])))
            assert qdiff < 1e-7 and tdiff < 1e-4, (qdiff, tdiff)
            assert np.max(actual[:, 24]) <= 80.0 + 1e-8
            assert (
                np.max(np.abs(np.diff(actual[:, 24])))
                <= options["rateNmPerSec"] * float(seq["dt"]) + 1e-8
            )
            pred = prediction_dict(actual, model["parameters"])
            pred["solverExitFlag"], pred["solveSeconds"] = flags, timings
            predictions.append(pred)
            nominals.append(prediction_dict(nominal, model["parameters"]))
            movement = ~np.asarray(seq["waiting"], dtype=bool)
            row = {
                "group": group,
                "source": LABELS[source],
                "minAngleRad": float(actual[:, 0].min()),
                "maxAngleRad": float(actual[:, 0].max()),
                "assistancePeakNm": float(actual[:, 24].max()),
                "softPeakNm": float(np.abs(actual[:, 25]).max()),
                "abdomenPeakNm": float(actual[:, 4].max()),
                "thetaRMSEDegrees": float(
                    np.sqrt(np.mean((actual[movement, 0] - seq["q"][movement]) ** 2))
                    * 180
                    / np.pi
                ),
                "zeroAssistanceAngleParityRad": qdiff,
                "zeroAssistanceTorqueParityNm": tdiff,
                "rejectedQPResults": int(np.sum(flags < 0)),
                "momentumResidualNm": float(np.abs(actual[:, 16]).max()),
            }
            summaries.append(row)
            print(json.dumps(row), flush=True)
            checks.append(
                {
                    "source": LABELS[source],
                    "group": group,
                    "parameterSHA256": hashlib.sha256(
                        np.r_[model["parameters"], model["neuralParameters"]].tobytes()
                    ).hexdigest(),
                }
            )
        result["groups"].append(
            {
                "group": group,
                "sequence": seq,
                "inputs": g["inputs"],
                "predictions": predictions,
                "nominal": nominals,
            }
        )
        savemat(
            target / "hill_mpc.mat",
            {"result": result},
            do_compression=True,
            long_field_names=True,
        )
    metadata = {
        "summary": summaries,
        "options": options,
        "fixedParameterAudit": checks,
        "modelSHA256": hashlib.sha256(cache_path.read_bytes()).hexdigest(),
        "policy": "Same selected 42-parameter neural/Hill model; causal MPC; experimental reference is evaluation-only",
        "assistanceAdaptation": "Applied exoskeleton torque is subtracted from the human muscle torque request; no assistance-response experiment was fitted",
        "wallSeconds": time.perf_counter() - started,
    }
    (target / "evaluation.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
