"""Calibrate the same compliant Hill equations for two activation sources.

All 32 model parameters have bounded optimization ranges. Body/load geometry
and soft-limit settings are shared between sources; muscle and activation
parameters may differ. Initial experimental states and references are unchanged.
"""
import argparse
import ctypes
import os
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat
from scipy.optimize import least_squares
from scipy.signal import find_peaks, lfilter

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/hill_compliant_geometry_fit"
DT = 0.01
PARAMETER_NAMES = (
    "backCapacity", "backArmToLength", "backReferenceLength", "jointDamping",
    "backMaximumVelocity", "backActivationExponent", "backTendonStrain",
    "backActivationTime", "frontCapacity", "frontActivationExponent",
    "frontArmToLength", "frontReferenceLength", "loadLever", "loadOffset",
    "frontInputGain", "backCoactivationGain", "risingGain", "fallingGain",
    "bodyMass", "bodyLever", "backMomentArm", "frontMomentArm",
    "backReferenceMTLength", "frontReferenceMTLength", "frontMaximumVelocity",
    "frontTendonStrain", "frontActivationTime", "referencePennation",
    "softLimitStiffness", "softLimitMaximumTorque", "softLimitLower", "softLimitUpper",
)
# The original thesis (Eq. 3.72-3.74) permits M=15..45 kg,
# L1=.20..35 m, L2=.345..85 m and phi=0..1.1 rad.
# Other bounds are engineering priors, not independently measured anatomy.
LOW = np.array([
    200, .25, 1.0, 0, 2, .4, .02, .005, 20, .4, .35, .7,
    .345, 0, .05, 0, 0, 0,
    15, .20, .020, .035, .22, .36, 2, .02, .005, .10,
    100, 40, -.10, 1.50,
])
HIGH = np.array([
    1600, 1.1, 1.3, 35, 15, 2.5, .10, .15, 1000, 1.6, .9, 1.3,
    .85, 1.1, 5, 1.5, 15, 15,
    45, .35, .040, .080, .30, .50, 15, .10, .15, .35,
    1200, 250, .30, 2.10,
])
SHARED = np.array([12, 13, 18, 19, 28, 29, 30, 31])
SPECIFIC = np.array([i for i in range(len(LOW)) if i not in SHARED])
TORQUE_WEIGHT = 2.0
TORQUE_SMOOTHNESS = .35
MONOTONIC_WEIGHT = 3.0
OBJECTIVE = "Joint geometry/force calibration; exact original reference; elastic soft limits"
DEFAULT_EXTRA = np.array([
    22.6006, .345, .0328531574526041, .075, .222, .408,
    5, .05, .04, .24086, 600, 150, .10, 1.8,
])

CORE = Path(__file__).with_name("hill_pair_core.dll" if os.name == "nt" else "hill_pair_core.so")
source_path = CORE.with_suffix(".c")
if not CORE.exists() or CORE.stat().st_mtime < source_path.stat().st_mtime:
    compiler = shutil.which("gcc")
    local_gcc = Path("C:/Library/CPP/mingw64-posix/bin/gcc.exe")
    if compiler is None and local_gcc.is_file():
        compiler = str(local_gcc)
    if compiler is None:
        raise RuntimeError("Fitting requires GCC; MATLAB replay has a compiler-free fallback.")
    flags = ["-shared", "-O3", "-static-libgcc"] if os.name == "nt" else ["-shared", "-O3", "-fPIC"]
    subprocess.run([compiler, *flags, str(source_path), "-o", str(CORE), "-lm"], check=True)
lib = ctypes.CDLL(str(CORE))
array_type = np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS")
lib.hill_pair_integrate.argtypes = [
    array_type, array_type, array_type, ctypes.c_int, ctypes.c_double,
    ctypes.c_double, ctypes.c_double, ctypes.c_int, array_type,
]
lib.hill_pair_integrate.restype = ctypes.c_int
PACKAGE = loadmat(ROOT / "data/processed/hill_compliant/calibrated_models.mat", simplify_cells=True)
trials = loadmat(ROOT / "data/processed/hill_compliant/paired_trials_original_reference.mat", simplify_cells=True)["trials"]


def upgrade(parameters):
    """Upgrade an earlier 18-parameter checkpoint without changing its inputs."""
    parameters = np.asarray(parameters).reshape(-1)
    if len(parameters) == 18:
        parameters = np.r_[parameters, DEFAULT_EXTRA]
    assert len(parameters) == 32
    return np.clip(parameters, LOW+1e-7*(HIGH-LOW), HIGH-1e-7*(HIGH-LOW))


def memory(signal, tau):
    """Causal first-order memory, initialized to the first available sample."""
    decay = np.exp(-DT/tau)
    return lfilter(1-decay, [1, -decay], signal, zi=[decay*signal[0]])[0]


def features(trial):
    """MAP-only predictors: excitation history and known load, without a clock."""
    back = trial["a"][:, 1]
    fast, slow = memory(back, .1), memory(back, .35)
    mass, nominal = trial["m"]/15, trial["loadKg"]/15
    return np.c_[
        back, back*back, fast, slow, back-slow, mass,
        np.full(len(back), nominal), back*nominal, back*mass,
        (back-fast)*nominal,
    ]


def input_generator(train):
    """Fit auxiliary front activation using training activation labels only."""
    predictors = np.vstack([features(t) for t in train])
    targets = np.concatenate([np.clip(t["front"], 1e-4, .95) for t in train])
    logits = np.log(targets/(1-targets))
    mu = predictors.mean(0)
    scale = np.maximum(predictors.std(0), .01)
    design = np.c_[np.ones(len(logits)), (predictors-mu)/scale]
    weights = np.concatenate([np.ones(len(t["q"]))/len(t["q"]) for t in train])
    penalty = np.eye(design.shape[1])*.1
    penalty[0, 0] = 0
    beta = np.linalg.solve((design.T*weights)@design+penalty, (design.T*weights)@logits)
    return {"mu": mu, "scale": scale, "beta": beta}


def pair_inputs(trial, source, generator, parameters=None):
    """Produce causal neural inputs; inferred recruitment is not a measurement."""
    back = trial["a"][:, source]
    if source == 0:
        front = trial["front"]
    else:
        predictors = features(trial)
        logits = np.c_[
            np.ones(len(back)), (predictors-generator["mu"])/generator["scale"],
        ]@np.asarray(generator["beta"]).reshape(-1)
        front = 1/(1+np.exp(-np.clip(logits, -30, 30)))
    if parameters is not None:
        difference = memory(back, .1)-memory(back, .35)
        recruitment = (
            parameters[14]*front+parameters[15]*back
            + parameters[16]*np.maximum(0, difference)
            + parameters[17]*np.maximum(0, -difference)
        )
        front = -np.expm1(-np.maximum(0, recruitment))
    result = np.asfortranarray(np.c_[back, front])
    assert result.shape == (len(back), 2), 'Activation pair must have exactly two channels.'
    return result


def boundary_torque(parameters, q):
    """Elastic soft limits: zero inside the dead zone, no measured angle input."""
    low = np.maximum(parameters[30]-q, 0)
    high = np.maximum(q-parameters[31], 0)
    low = low*low/(low+.05)
    high = high*high/(high+.05)
    return parameters[29]*(
        np.tanh(parameters[28]*low/parameters[29])
        - np.tanh(parameters[28]*high/parameters[29])
    )


def rollout(parameters, trial, source, generator, steps=2):
    """The solver sees initial state, neural input and mass, never the references."""
    inputs = pair_inputs(trial, source, generator, parameters)
    output = np.zeros((len(inputs), 17), order="F")
    failed = lib.hill_pair_integrate(
        np.asfortranarray(parameters), inputs, np.asfortranarray(trial["m"]),
        len(inputs), float(trial["q"][0]), float(trial["v0"]), DT, steps, output,
    )
    return output, failed


def training_residual(parameters, train, source, generator):
    """Angle trend and torque accuracy, with boundary force exposed and penalized."""
    residuals = []
    for trial in train:
        output, failed = rollout(parameters, trial, source, generator)
        q, velocity, torque = output[:, 0], output[:, 1], output[:, 2]
        normalizer = np.sqrt(len(q))
        reference = trial["q"]
        valleys = find_peaks(-reference, prominence=.12, distance=50)[0]
        peaks = find_peaks(reference, prominence=.12, distance=50)[0]
        anchors = np.unique(np.r_[0, valleys, peaks, len(q)-1])
        direction = np.zeros(len(q))
        for left, right in zip(anchors[:-1], anchors[1:]):
            direction[left:right] = np.sign(reference[right]-reference[left])
        reference_scale = max(np.std(reference), .12)
        prediction_scale = max(np.std(q), .12)
        reference_shape = (reference-reference.mean())/reference_scale
        prediction_shape = (q-q.mean())/prediction_scale
        correlation = np.corrcoef(q, reference)[0, 1] if np.std(q) > 1e-6 else 0
        torque_error = torque-trial["torque"]
        frames = np.r_[
            prediction_shape-reference_shape,
            .5*(q-reference)/.4,
            TORQUE_WEIGHT*torque_error/70,
            .25*(velocity/prediction_scale-np.gradient(reference, DT)/reference_scale)/3,
            8*(np.maximum(0, -q)+np.maximum(0, q-1.8)),
            TORQUE_SMOOTHNESS*np.diff(torque_error)/3,
            .12*output[:, 6]/100, .10*output[:, 9]/50,
            1.2*boundary_torque(parameters, q)/50,
        ]/normalizer
        combined = np.r_[
            frames,
            .75*(prediction_shape[anchors]-reference_shape[anchors])/np.sqrt(len(anchors)),
            1.2*(1-correlation),
            MONOTONIC_WEIGHT*np.maximum(0, -direction*velocity)/(3*prediction_scale)/normalizer,
        ]
        if failed:
            combined += 10  # Invalid states are rejected, never clipped and plotted.
        residuals.append(np.nan_to_num(combined, nan=1e5, posinf=1e5, neginf=-1e5)/np.sqrt(len(train)))
    return np.concatenate(residuals)


def priors(parameters, seed):
    """Weak regularization for parameter combinations unsupported by these data."""
    weights = np.full(32, .08)
    weights[[18, 19]] = .25
    weights[[20, 21, 22, 23, 27]] = .20
    weights[[28, 29, 30, 31]] = .20
    return weights*(parameters-seed)/(HIGH-LOW)


def fit_pair(train, generator, initial, sweeps=3, max_nfev=90):
    """Alternate source-specific and common blocks; every parameter can change."""
    pair = np.vstack([upgrade(p) for p in initial])
    pair[:, SHARED] = pair[:, SHARED].mean(0)
    prior = pair.copy()
    # A small positive onset margin protects against zero-angle crossings
    # under unseen low activation. This is a soft-spring prior, not a change
    # to the experimental reference or a hard bound on the predicted angle.
    prior[:, 30] = .10
    begin = time.time()
    status = 0

    def total(current):
        return np.r_[
            training_residual(current[0], train, 0, generator)/np.sqrt(2),
            training_residual(current[1], train, 1, generator)/np.sqrt(2),
            priors(current[0], prior[0])/np.sqrt(2),
            priors(current[1], prior[1])/np.sqrt(2),
        ]

    best_value = np.inf
    best_pair = pair.copy()
    for sweep in range(sweeps):
        for source in (0, 1):
            start = (pair[source, SPECIFIC]-LOW[SPECIFIC])/(HIGH[SPECIFIC]-LOW[SPECIFIC])
            def one(z):
                parameters = pair[source].copy()
                parameters[SPECIFIC] = LOW[SPECIFIC]+z*(HIGH[SPECIFIC]-LOW[SPECIFIC])
                return np.r_[training_residual(parameters, train, source, generator), priors(parameters, prior[source])]
            result = least_squares(one, start, bounds=(np.zeros(len(start)), np.ones(len(start))),
                                   max_nfev=max_nfev, ftol=5e-5, diff_step=1e-5)
            pair[source, SPECIFIC] = LOW[SPECIFIC]+result.x*(HIGH[SPECIFIC]-LOW[SPECIFIC])
            status = result.status
            print("SOURCE BLOCK", source+1, "SWEEP", sweep+1, "SECONDS",
                  round(time.time()-begin, 1), "NFEV", result.nfev, flush=True)
        start = (pair[0, SHARED]-LOW[SHARED])/(HIGH[SHARED]-LOW[SHARED])
        def common(z):
            current = pair.copy()
            current[:, SHARED] = LOW[SHARED]+z*(HIGH[SHARED]-LOW[SHARED])
            return total(current)
        result = least_squares(common, start, bounds=(np.zeros(len(start)), np.ones(len(start))),
                               max_nfev=max_nfev, ftol=5e-5, diff_step=1e-5)
        pair[:, SHARED] = LOW[SHARED]+result.x*(HIGH[SHARED]-LOW[SHARED])
        errors = total(pair)
        value = errors@errors
        status = result.status
        print("SWEEP", sweep+1, "TRAIN OBJECTIVE", round(value, 5),
              "SECONDS", round(time.time()-begin, 1),
              "M,L1,L2,phi", pair[0, [18, 19, 12, 13]], flush=True)
        if value < best_value:
            best_value = value
            best_pair = pair.copy()
            for source in (0, 1):
                savemat(OUT / f"progress_source{source+1}.mat", {
                    "parameters": pair[source], "generator": generator, "loss": value,
                    "objective": OBJECTIVE, "sharedParameterIndices": SHARED+1,
                })
    return best_value, best_pair, status


def engineering_pair(previous):
    """Expand an earlier same-training-set checkpoint and raise strength with load."""
    pair = np.vstack([upgrade(p) for p in previous])
    pair[:, 0] *= 1.6
    pair[:, 8] *= 1.2
    pair[:, [18, 19, 12, 13]] = [30, .345, .76, .55]
    return np.clip(pair, LOW+1e-6*(HIGH-LOW), HIGH-1e-6*(HIGH-LOW))


def main():
    global OUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="hill_compliant_geometry_fit")
    parser.add_argument("--sweeps", type=int, default=3)
    parser.add_argument("--max-nfev", type=int, default=90)
    args = parser.parse_args()
    OUT = ROOT / "results" / args.output
    OUT.mkdir(parents=True, exist_ok=True)
    generator = input_generator(trials)
    previous = [PACKAGE["models"][s]["parameters"] for s in (0, 1)]
    initial = engineering_pair(previous)
    checkpoints = [OUT / f"calibration_source{s}.mat" for s in (1, 2)]
    if all(path.exists() for path in checkpoints):
        initial = np.vstack([loadmat(path, simplify_cells=True)["parameters"] for path in checkpoints])
    value, pair, status = fit_pair(trials, generator, initial, args.sweeps, args.max_nfev)
    for source in (0, 1):
        savemat(OUT / f"calibration_source{source+1}.mat", {
            "parameters": pair[source], "loss": value, "status": status,
            "generator": generator, "objective": OBJECTIVE,
            "sharedParameterIndices": SHARED+1,
        })
    savemat(OUT / "paired_trials.mat", {"trials": trials})
    print("FINISHED", flush=True)


if __name__ == "__main__":
    main()
