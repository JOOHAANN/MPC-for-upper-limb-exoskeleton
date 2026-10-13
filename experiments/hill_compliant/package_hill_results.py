"""Export MATLAB structs and fitted models for a portable replay package."""
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "results/hill_compliant_geometry_fit"
CACHE = ROOT / "data/processed/hill_compliant"


def cell_row(items):
    cells = np.empty((1, len(items)), dtype=object)
    for i, item in enumerate(items):
        cells[0, i] = item
    return cells


def struct_row(items):
    """Lists of dictionaries become MATLAB cells; explicitly export structs."""
    fields = list(items[0])
    records = np.empty((1, len(items)), dtype=[(name, object) for name in fields])
    for i, item in enumerate(items):
        for name in fields:
            value = item[name]
            if isinstance(value, (int, float, np.number)):
                value = float(value)
            if isinstance(value, np.ndarray) and value.ndim == 1:
                value = value.reshape(-1, 1)
            records[name][0, i] = value
    return records


def load_model(path):
    model = loadmat(path, simplify_cells=True)
    model = {key: value for key, value in model.items() if not key.startswith("_")}
    model.pop("predictions", None)
    model.pop("failed", None)
    model["parameters"] = np.asarray(model["parameters"]).reshape(1, -1)
    model["generator"]["mu"] = model["generator"]["mu"].reshape(1, -1)
    model["generator"]["scale"] = model["generator"]["scale"].reshape(1, -1)
    model["generator"]["beta"] = model["generator"]["beta"].reshape(-1, 1)
    return model


def main():
    from fit_hill_compliant import pair_inputs
    CACHE.mkdir(parents=True, exist_ok=True)
    trials = loadmat(OUTPUT / "paired_trials.mat", simplify_cells=True)["trials"]
    models = cell_row([load_model(OUTPUT / f"calibration_source{i}.mat") for i in (1, 2)])
    holdouts = np.empty((2, 3), dtype=object)
    for source in (1, 2):
        for fold, speed in enumerate((35, 45, 55)):
            model = load_model(OUTPUT / f"validation_source{source}_speed{speed}.mat")
            model["heldSpeed"] = speed
            holdouts[source - 1, fold] = model
    # Independent Python input traces catch omissions in the MATLAB replay
    # mapping after parameter-vector changes, before judging physics or fits.
    calibration_inputs = np.empty((2, len(trials)), dtype=object)
    validation_inputs = np.empty_like(calibration_inputs)
    for source in (0, 1):
        for k, trial in enumerate(trials):
            model = models[0, source]
            calibration_inputs[source, k] = pair_inputs(trial, source, model['generator'], model['parameters'].reshape(-1))
            fold = [35, 45, 55].index(int(trial['speed']))
            model = holdouts[source, fold]
            validation_inputs[source, k] = pair_inputs(trial, source, model['generator'], model['parameters'].reshape(-1))
    for trial in trials:
        trial["coordinateDefinition"] = "Exact original project Scale/reference: Combination_raw theta column"
        trial["initialStateDefinition"] = "Original reference start and original past-sample velocity"
    savemat(CACHE / "calibrated_models.mat", {
        "models": models, "validationModels": holdouts, "heldSpeeds": [35, 45, 55],
        "trials": struct_row(trials), "dt": 0.01, "integrationSubsteps": 10,
        "modelDescription": "Opposing compliant Hill units; 32 optimized parameters; shared body/load geometry; elastic soft limits",
        "referencePolicy": "Exact original project reference; no new rescale, filtering or seam removal",
        "angleCheckRange": np.array([0.,1.8]),
        "fiberInitialization": "Static fiber equilibrium at initial pose; original joint velocity preserved",
        "sharedParameterIndices": np.array([13, 14, 19, 20, 29, 30, 31, 32]),
        "calibrationInputAudit": calibration_inputs,
        "validationInputAudit": validation_inputs,
    }, do_compression=True)
    old = ROOT / "results/bse_torque_theta/joint_comparison.mat"
    if old.exists():
        archive = loadmat(old, simplify_cells=True)["results"]["archived"]
        savemat(CACHE / "archived_bse.mat", {
            "archived": cell_row(archive),
            "description": "Best of three archived BSE variants by all-nine-trial combined score; includes reference feedback",
        }, do_compression=True)
    print("Saved portable Hill package:", CACHE)


if __name__ == "__main__":
    main()
