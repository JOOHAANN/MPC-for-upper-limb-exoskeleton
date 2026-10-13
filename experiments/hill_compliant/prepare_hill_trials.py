"""Recover auxiliary activation while preserving the ORIGINAL reference exactly.

The portable MATLAB package already includes these arrays. This optional step
requires the local 240528Pilot recording directory. It does not fit a model.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat
from scipy.signal import lfilter

ROOT = Path(__file__).resolve().parents[2]
GROUPS = ["5_35_10_55_15_45", "5_45_10_35_15_55", "5_55_10_45_15_35"]
CONDITIONS = [(5, 35), (10, 55), (15, 45), (5, 45), (10, 35),
              (15, 55), (5, 55), (10, 45), (15, 35)]
WINDOWS = [(5321, 11740), (22037, 26059), (36732, 41109), (52761, 58381),
           (66452, 72417), (84099, 87884), (98004, 102357),
           (113394, 118437), (127817, 133717)]
CONTACT = [(6742, 10348), (22885, 25303), (37516, 40124), (53975, 57348),
           (67641, 71085), (85027, 87002), (99146, 101117),
           (114510, 117291), (128682, 132730)]


def read_sto(path):
    lines = path.read_text().splitlines()
    header = next(i for i, line in enumerate(lines) if line.strip() == "endheader")
    return np.loadtxt(lines[header+2:])


def causal_filter(signal, tau):
    decay = np.exp(-0.01/tau)
    return lfilter(1-decay, [1, -decay], signal, zi=[decay*signal[0]])[0]


def build_trials(raw_root):
    manifest = json.loads((ROOT / "config/workflows/training_map__muscle_tendon__general_model_op_scale.json").read_text())
    raw, generated = [], []
    for group in GROUPS:
        folder = ROOT / "data/processed" / f"training_map__muscle_tendon__{group}_done__build_data"
        raw.append(loadmat(folder / "Combination_raw.mat")["Combination"])
        entry = next(item for item in manifest["files"] if item["target"] == f"Combination_{group}.mat")
        generated.append(loadmat(ROOT / entry["source"])["Combination"])
        assert np.array_equal(raw[-1][:, 2], generated[-1][:, 2])
        assert np.array_equal(raw[-1][:, 0], generated[-1][:, 3])
    lengths = [len(group) for group in raw]
    combined, mapped = np.vstack(raw), np.vstack(generated)
    offsets, trials = [0, 0, 0], []
    for k, ((load, speed), (left, right), contact) in enumerate(zip(CONDITIONS, WINDOWS, CONTACT)):
        base = raw_root / f"dynamics_stoop_{load}kg_{speed}"
        lmt, activation = [read_sto(base / name) for name in ("lmt.sto", "Activations.sto")]
        start, end = (49 if k % 3 == 2 else 0), (349 if k % 3 == 1 else 299)
        count = int(np.floor((lmt[end, 0]-lmt[start, 0])/.001+.0003))+1
        group = k // 3
        indices = np.arange(left, right-1, 10)
        local = indices-sum(lengths[:group])-offsets[group]
        timestamps = lmt[start, 0]+local*.001
        # Preserve the hand-written project's final Scale/reference pipeline.
        # Select samples directly: no new rescale, smoothing or seam removal.
        q = combined[indices, 2].copy()
        previous = combined[indices[0]-10, 2]
        front = np.interp(timestamps, activation[:, 0], activation[:, 1])
        attached = load*((indices+1 > contact[0]) & (indices+1 < contact[1])).astype(float)
        # Effective load transfer is assumed, with a fixed 80 ms transition.
        mass = causal_filter(attached, .08)
        trials.append({
            "a": np.c_[combined[indices, 1], mapped[indices, 1]],
            "q": q, "qArchived": combined[indices, 2], "v0": (q[0]-previous)/.01,
            "torque": combined[indices, 0], "m": mass, "front": causal_filter(front, .12),
            "loadKg": float(load), "speed": float(speed), "t": np.arange(len(indices))*.01,
        })
        offsets[group] += count
    assert offsets == lengths
    return trials


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-root", type=Path, default=ROOT.parent / "YOUHAN LI code and project description/Data/240528Pilot")
    args = parser.parse_args()
    trials = build_trials(args.raw_root)
    out = ROOT / "results/hill_compliant_geometry_fit"
    out.mkdir(parents=True, exist_ok=True)
    from package_hill_results import struct_row
    paired = ROOT / "data/processed/hill_compliant/paired_trials_original_reference.mat"
    savemat(paired, {"trials": struct_row(trials),
        "referencePolicy": "Exact original Combination_raw theta column; no reference modification"}, do_compression=True)
    savemat(out / "paired_trials.mat", {"trials": trials}, do_compression=True)
    cache = ROOT / "data/processed/hill_compliant/calibrated_models.mat"
    if cache.exists():
        previous = loadmat(cache, simplify_cells=True)["trials"]
        for t, old in zip(trials, previous):
            for name in ("a", "m", "front", "torque"):
                assert np.max(np.abs(t[name]-old[name])) < 1e-10, name
        legacy = loadmat(ROOT / "data/processed/hill_forward/calibrated_models.mat", simplify_cells=True)["trials"]
        for t, old in zip(trials, legacy):
            assert np.array_equal(t["q"], old["q"]) and t["v0"] == old["v0"]
        print("PASS: reference and initial state equal the original archived project, sample for sample.")
    print("Saved:", out / "paired_trials.mat")


if __name__ == "__main__":
    main()
