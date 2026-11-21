"""Choose a deeper-flexion warm start across all nine cases per source.

Only fixed neural parameters are varied here. No output is rescaled and no
reference is supplied to the Hill solver. Full physical recalibration follows.
"""
import argparse
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat

import fit_hill_neural as neural
from fit_hill_compliant import ROOT


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=ROOT / "data/processed/hill_compliant/neural_models.mat")
    parser.add_argument("--output", default="hill_neural_depth")
    args = parser.parse_args()
    package = loadmat(args.model, simplify_cells=True)
    pair = np.vstack([m["parameters"] for m in package["models"]])
    start = neural.parameter_rows(np.vstack([m["neuralParameters"] for m in package["models"]]))
    records = neural.fit.prepared([1, 2, 3], [m["generator"] for m in package["models"]])
    neural.fit.rollout = neural.rollout
    neural.fit.ANGLE_PROFILE = "depth"
    neural.NEURAL[:] = start
    out = ROOT / "results" / args.output
    out.mkdir(parents=True, exist_ok=True)
    for source in (0, 1):
        best = (np.inf, start[source].copy())
        for factor in (1., 1.25, 1.5):
            for blend in (start[source, 6], .85, 1.):
                candidate = start[source].copy()
                candidate[0] = min(neural.NHIGH[0], start[source, 0] * factor)
                candidate[6] = blend
                neural.NEURAL[source] = candidate
                residual = neural.fit.residual(pair[source], records, source, steps=5)
                score = float(residual @ residual)
                if score < best[0]:
                    best = (score, candidate.copy())
                print("SOURCE", source + 1, "AMPLITUDE", round(candidate[0], 3), "BLEND", round(blend, 3), "LOSS", round(score, 4), flush=True)
        neural.NEURAL[source] = best[1]
    savemat(out / "seed_pair.mat", {"parameters": pair, "neuralParameters": neural.NEURAL.copy()})
    neural.fit.write_metrics(out, neural.fit.metrics(pair, records), "warmstart_metrics.csv")
    print("Depth warm start:", out / "seed_pair.mat", flush=True)


if __name__ == "__main__":
    main()
