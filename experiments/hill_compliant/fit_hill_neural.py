"""Neural-control candidate using the unchanged Hill mechanical system.

EMG drives a motion-intent estimator; proprioceptive control modifies muscle
excitation. The estimator is not the predicted angle. Both joint motion and
torque still come from forward Hill mechanics. No reference enters the solver.
"""
import argparse
import ctypes
import json
import os
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat
from scipy.optimize import least_squares

import fit_hill_physiology as fit
from fit_hill_compliant import ROOT

DLL = Path(__file__).with_name("hill_neural_core.dll" if os.name == "nt" else "hill_neural_core.so")
dependencies = [Path(__file__).with_name(name) for name in ("hill_neural_core.c", "hill_sequence_core.c", "hill_pair_core.c")]
if not DLL.exists() or DLL.stat().st_mtime < max(p.stat().st_mtime for p in dependencies):
    compiler = shutil.which("gcc") or "C:/Library/CPP/mingw64-posix/bin/gcc.exe"
    flags = ["-shared", "-O3", "-static-libgcc"] if os.name == "nt" else ["-shared", "-O3", "-fPIC"]
    subprocess.run([compiler, *flags, str(dependencies[0]), "-o", str(DLL), "-lm"], check=True)
lib = ctypes.CDLL(str(DLL))
array_type = np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS")
lib.hill_neural_integrate.argtypes = [array_type, array_type, array_type, array_type,
                                    ctypes.c_int, array_type, ctypes.c_double,
                                    ctypes.c_int, array_type, array_type, array_type]
lib.hill_neural_integrate.restype = ctypes.c_int
NEURAL_NAMES = ["intentAmplitude", "intentExponent", "intentTau", "postureKp", "postureKd",
                "maximumNegativeTorqueRequest", "neuralFeedbackBlend", "coactivationNm", "causalLeadTime",
                "loadTonicCompensation"]
NLOW = np.array([.7, .3, .01, 30, 10, 3, .7, 0., 0., 0.])
NHIGH = np.array([3., 2., .25, 250, 90, 40, 1., 8., .8, 2.5])
NEURAL = np.tile([1.6, .8, .06, 150, 35, 15, 1, 1, .35, .5], (2, 1))
ORIGINAL_ROLLOUT = fit.rollout


def parameter_rows(values):
    """Preserve nine-parameter historical models with zero tonic compensation."""
    rows = np.asarray(values, dtype=float).reshape(2, -1)
    if rows.shape[1] == 9:
        rows = np.column_stack((rows, np.zeros(2)))
    assert rows.shape == (2, len(NEURAL_NAMES))
    return rows


def rollout(parameters, record, source, steps=5):
    t = record["trial"]
    a, inferred = fit.inputs(parameters, record, source)
    o = np.zeros((len(a), 24), order="F")
    cfg = fit.CONTROLLER.copy()
    cfg[0] = t["standingAngle"]
    failed = lib.hill_neural_integrate(np.asfortranarray(parameters), a,
                                      np.asfortranarray(t["m"]), np.asfortranarray(t["waiting"], dtype=float),
                                      len(a), np.asfortranarray([t["q"][0], t["v0"]]), fit.DT, steps,
                                      np.asfortranarray(cfg), np.asfortranarray(NEURAL[source]), o)
    return o, bool(failed), a, inferred


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="hill_neural")
    parser.add_argument("--seed-only", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--sweeps", type=int, default=2)
    parser.add_argument("--max-nfev", type=int, default=55)
    parser.add_argument("--steps", type=int, default=5)
    parser.add_argument("--initial-model", type=Path)
    parser.add_argument("--physiology-scale", type=float, default=1.)
    parser.add_argument("--back-activation-weight", type=float, default=0.)
    parser.add_argument("--transition-weight", type=float, default=0.)
    parser.add_argument("--angle-profile", choices=("shape", "depth"), default="depth")
    args = parser.parse_args()
    out = ROOT / "results" / args.output
    out.mkdir(parents=True, exist_ok=True)
    generators = [m["generator"] for m in fit.original.PACKAGE["models"]]
    records = fit.prepared([1, 2, 3], generators)
    pair = fit.seed_pair()
    pair[:, [0, 1, 2, 5, 8, 9]] = [500, .35, .85, 1.4, 250, 1]
    pair[:, 15:18] = 0
    fit.PHYSIOLOGY_SCALE = args.physiology_scale
    fit.BACK_ACTIVATION_WEIGHT = args.back_activation_weight
    fit.TRANSITION_WEIGHT = args.transition_weight
    fit.ANGLE_PROFILE = args.angle_profile
    if args.initial_model:
        initial = loadmat(args.initial_model, simplify_cells=True)
        pair = np.array(initial["parameters"]) if "parameters" in initial else np.vstack([m["parameters"] for m in initial["models"]])
        values = initial["neuralParameters"] if "neuralParameters" in initial else np.vstack([m["neuralParameters"] for m in initial["models"]])
        NEURAL[:] = parameter_rows(values)
        if args.angle_profile == "depth" and np.asarray(values).shape[-1] == 9:
            NEURAL[:, 9] = .5
        if "controller" in initial:
            fit.CONTROLLER[1:] = np.asarray(initial["controller"]).reshape(6)
    if args.resume and (out / "progress_pair.mat").exists():
        pair = loadmat(out / "progress_pair.mat", simplify_cells=True)["parameters"]
        if (out / "neural_parameters.mat").exists():
            NEURAL[:] = parameter_rows(loadmat(out / "neural_parameters.mat", simplify_cells=True)["neuralParameters"])
    fit.rollout = rollout
    fit.write_metrics(out, fit.metrics(pair, records, steps=10), "seed_metrics.csv")
    print("SEED METRICS SAVED", flush=True)
    if args.seed_only:
        savemat(out / "seed_pair.mat", {"parameters": pair, "neuralParameters": NEURAL})
        return
    for sweep in range(args.sweeps):
        # Calibrate the two fixed neural policies, then all 32 Hill parameters.
        for source in (0, 1):
            count, last = 0, time.monotonic()
            best_neural, best_value = NEURAL[source].copy(), np.inf
            def one(z):
                nonlocal count, last, best_neural, best_value
                NEURAL[source] = NLOW + z * (NHIGH - NLOW)
                residual = fit.residual(pair[source], records, source, args.steps)
                count += 1
                value = float(residual @ residual)
                if value < best_value:
                    best_neural, best_value = NEURAL[source].copy(), value
                if time.monotonic() - last > 30:
                    print("NEURAL SOURCE", source + 1, "SWEEP", sweep + 1, "EVALUATIONS", count,
                          "BEST LOSS", round(best_value, 5), flush=True)
                    savemat(out / f"partial_neural_source{source + 1}.mat",
                            {"neuralParameters": best_neural, "loss": best_value})
                    last = time.monotonic()
                return residual
            start = np.clip((NEURAL[source] - NLOW) / (NHIGH - NLOW), 1e-8, 1 - 1e-8)
            result = least_squares(one, start, bounds=(np.zeros(len(NLOW)), np.ones(len(NLOW))),
                                   max_nfev=args.max_nfev, ftol=3e-4, diff_step=2e-5)
            NEURAL[source] = NLOW + result.x * (NHIGH - NLOW)
            print("NEURAL FIT", sweep + 1, source + 1, "LOSS", result.fun @ result.fun, "PARAMETERS", NEURAL[source], flush=True)
            savemat(out / "neural_parameters.mat", {"neuralParameters": NEURAL.copy()})
            savemat(out / "neural_policy_pair.mat", {"parameters": pair, "neuralParameters": NEURAL.copy()})
        fit.write_metrics(out, fit.metrics(pair, records, steps=10), "neural_policy_metrics.csv")
        pair, loss = fit.fit(records, pair, out, 1, args.max_nfev, args.steps)
        savemat(out / "neural_parameters.mat", {"neuralParameters": NEURAL.copy()})
        fit.write_metrics(out, fit.metrics(pair, records, steps=10), "progress_metrics.csv")
    cells = np.empty((1, 2), dtype=object)
    for source in (0, 1):
        cells[0, source] = {"parameters": pair[source].reshape(1, -1), "generator": generators[source],
                            "neuralParameters": NEURAL[source].reshape(1, -1),
                            "objective": "EMG-informed neural policy + same Hill mechanics, fixed across all nine conditions"}
    savemat(out / "physiology_models.mat", {"models": cells, "parameters": pair,
                                             "neuralParameters": NEURAL.copy(),
                                             "controller": fit.CONTROLLER[1:],
                                             "angleProfile": args.angle_profile,
                                             "sharedParameterIndices": fit.SHARED + 1})
    rows = fit.metrics(pair, records, steps=10)
    fit.write_metrics(out, rows)
    (out / "policy.json").write_text(json.dumps({"model": "EMG-informed intent and muscle excitation feedback",
                                                 "mechanics": "Identical compliant Hill, tendon and joint equations",
                                                 "referenceInput": False, "externalPostureTorque": 0,
                                                 "parametersAcrossConditions": "One fixed vector per source across all nine cases",
                                                 "parameterNames": NEURAL_NAMES,
                                                 "physiologyScale": args.physiology_scale,
                                                 "backActivationWeight": args.back_activation_weight,
                                                 "transitionWeight": args.transition_weight,
                                                 "angleProfile": args.angle_profile,
                                                 "anglePolicy": fit.ANGLE_PROFILES[args.angle_profile],
                                                 "physiologicalStatus": "Mechanistically constrained candidate; not independent muscle-force validation"}, indent=2), encoding="utf-8")
    print("FINISHED", out, flush=True)


if __name__ == "__main__":
    main()
