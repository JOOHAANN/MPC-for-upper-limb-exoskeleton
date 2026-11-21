"""Publish only a screened fixed-parameter neural/Hill replay package."""
import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat

from fit_hill_compliant import ROOT, SHARED


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--evaluation", type=Path, required=True)
    args = parser.parse_args()
    check = json.loads((args.evaluation / "screening.json").read_text(encoding="utf-8"))
    required = {"entireSequenceAnglesWithinAllowedRange", "fullSequenceAbdominalPeakBelow50Nm",
                "allValleysWithinReferenceDepthTolerance", "allExcursionsAtLeast075Reference",
                "allMiddleReboundsAtLeast065Reference"}
    assert required <= check["checks"].keys(), "Full transitions and flexion-depth checks are required."
    assert check["screenPassed"] and all(check["checks"].values()), "Candidate screening failed; no promotion."
    assert check["modelSHA256"] == sha(args.model), "Model changed after evaluation."
    package = loadmat(args.model, simplify_cells=True)
    pair = package["parameters"]
    assert np.array_equal(pair[0, SHARED], pair[1, SHARED])
    cells = np.empty((1, 2), dtype=object)
    for source in (0, 1):
        model = package["models"][source].copy()
        model["parameters"] = np.asarray(model["parameters"]).reshape(1, 32)
        model["neuralParameters"] = np.asarray(model["neuralParameters"]).reshape(1, -1)
        generator = model["generator"].copy()
        generator["mu"] = generator["mu"].reshape(1, -1)
        generator["scale"] = generator["scale"].reshape(1, -1)
        generator["beta"] = generator["beta"].reshape(-1, 1)
        model["generator"] = generator
        cells[0, source] = model
    cache = ROOT / "data/processed/hill_compliant/neural_models.mat"
    savemat(cache, {"models": cells, "dt": .01, "integrationSubsteps": 10,
                   "controller": package["controller"], "sharedParameterIndices": SHARED + 1,
                   "angleProfile": check["angleProfile"], "anglePolicy": check["anglePolicy"],
                   "referencePolicy": "Exact original Scale; motion references never enter the forward solver",
                   "evaluationPolicy": "Same fixed parameter vector per source over nine calibrated conditions; not an independent holdout",
                   "modelDescription": "Load-aware EMG intent and muscle excitation feedback with unchanged compliant Hill mechanics"}, do_compression=True)
    snapshots = ROOT / "docs/hill_neural_evaluation"
    snapshots.mkdir(parents=True, exist_ok=True)
    for path in args.evaluation.iterdir():
        if path.suffix in (".png", ".pdf", ".csv", ".json"):
            shutil.copy2(path, snapshots / path.name)
    with (args.evaluation / "metrics.csv").open(newline="", encoding="utf-8") as stream:
        rows = [r for r in csv.DictReader(stream) if r["model"] == "candidate"]
    hashes = []
    controller_hash = hashlib.sha256(np.asarray(package["controller"]).astype("<f8").tobytes()).hexdigest()
    for source in (1, 2):
        selected = [r for r in rows if int(r["source"]) == source]
        assert len(selected) == 9
        vector = np.r_[pair[source - 1], package["neuralParameters"][source - 1]]
        digest = hashlib.sha256(vector.astype("<f8").tobytes()).hexdigest()
        for r in selected:
            hashes.append({"source": source, "group": r["group"], "loadKg": r["loadKg"],
                           "speedCode": r["speedCode"], "parameterCount": len(vector), "fixedParameterSHA256": digest,
                           "sharedStandingControllerSHA256": controller_hash})
    with (snapshots / "fixed_parameter_audit.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(hashes[0]))
        writer.writeheader()
        writer.writerows(hashes)
    metadata = {"fitModelSHA256": sha(args.model), "portableModelSHA256": sha(cache),
                "snapshotPolicy": "Numerically screened model predictions, not direct muscle torque measurements",
                "files": [{"file": p.name, "sha256": sha(p)} for p in sorted(snapshots.iterdir())
                          if p.suffix in (".png", ".pdf", ".csv")]}
    (snapshots / "export_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print("Saved fixed models:", cache)
    print("PASS: one fixed parameter vector per source across nine conditions; common body geometry.")
    print("Snapshots:", snapshots)


if __name__ == "__main__":
    main()
