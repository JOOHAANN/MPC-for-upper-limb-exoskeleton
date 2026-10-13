"""Whole-speed holdouts: jointly calibrated body geometry, no test-label fitting."""
import argparse
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.io import loadmat, savemat

import fit_hill_compliant as fitter
from fit_hill_compliant import ROOT, trials, input_generator, fit_pair, rollout, engineering_pair


def model_field(model, name):
    return model[name] if isinstance(model, dict) else getattr(model, name)


def fit_fold(speed, output, sweeps, max_nfev):
    out = ROOT / "results" / output
    out.mkdir(parents=True, exist_ok=True)
    fitter.OUT = out / f"progress_speed{speed}"
    fitter.OUT.mkdir(exist_ok=True)
    train_ids = [i for i, t in enumerate(trials) if t["speed"] != speed]
    test_ids = [i for i, t in enumerate(trials) if t["speed"] == speed]
    train = [trials[i] for i in train_ids]
    generator = input_generator(train)
    # The previous version's exact same fold is a legitimate training-only
    # warm start. Never initialize holdouts from an all-nine-trial calibration.
    baseline = fitter.PACKAGE
    fold = [35, 45, 55].index(speed)
    previous = []
    for source in (0, 1):
        old = baseline["validationModels"][source, fold]
        assert np.array_equal(model_field(old, "trainIds"), np.array(train_ids)+1)
        assert np.array_equal(model_field(old, "testIds"), np.array(test_ids)+1)
        previous.append(model_field(old, "parameters"))
    initial = engineering_pair(previous)
    checkpoints = [out / f"validation_source{s}_speed{speed}.mat" for s in (1, 2)]
    if all(path.exists() for path in checkpoints):
        old = [loadmat(path, simplify_cells=True) for path in checkpoints]
        for model in old:
            assert np.array_equal(model["trainIds"], np.array(train_ids)+1)
            assert np.array_equal(model["testIds"], np.array(test_ids)+1)
            assert model["objective"] == fitter.OBJECTIVE
        initial = np.vstack([model["parameters"] for model in old])
    print("HELD OUT SPEED", speed, "START", flush=True)
    value, pair, status = fit_pair(train, generator, initial, sweeps, max_nfev)
    for source in (0, 1):
        predictions, failures = [], []
        for i in test_ids:
            prediction, failed = rollout(pair[source], trials[i], source, generator, steps=10)
            predictions.append(prediction)
            failures.append(failed)
        savemat(checkpoints[source], {
            "parameters": pair[source], "generator": generator, "loss": value, "status": status,
            "objective": fitter.OBJECTIVE, "sharedParameterIndices": fitter.SHARED+1,
            "trainIds": np.array(train_ids)+1, "testIds": np.array(test_ids)+1,
            "predictions": predictions, "failed": failures,
        })
    print("HELD OUT SPEED", speed, "FINISHED", flush=True)
    return speed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="hill_compliant_geometry_fit")
    parser.add_argument("--speed", type=int, choices=[35, 45, 55])
    parser.add_argument("--sweeps", type=int, default=3)
    parser.add_argument("--max-nfev", type=int, default=90)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    speeds = [args.speed] if args.speed else [35, 45, 55]
    if args.workers > 1:
        with ProcessPoolExecutor(max_workers=min(args.workers, len(speeds))) as pool:
            futures = [pool.submit(fit_fold, s, args.output, args.sweeps, args.max_nfev) for s in speeds]
            for future in futures:
                future.result()
    else:
        for speed in speeds:
            fit_fold(speed, args.output, args.sweeps, args.max_nfev)
    print("VALIDATION FINISHED", flush=True)


if __name__ == "__main__":
    main()
