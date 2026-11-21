"""Choose one common standing policy from complete continuous records.

Tune excitation feedback, not an externally applied torque or angle clip.
Transient abdominal braking and settled posture are both evaluated.
"""
import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat

import fit_hill_neural as neural
from fit_hill_compliant import ROOT, boundary_torque


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--output", default="hill_neural_depth_refined")
    args = parser.parse_args()
    package = loadmat(args.model, simplify_cells=True)
    pair = np.asarray(package["parameters"])
    neural.NEURAL[:] = neural.parameter_rows(package["neuralParameters"])
    records = neural.fit.prepared([1, 2, 3], [m["generator"] for m in package["models"]])
    out = ROOT / "results" / args.output
    out.mkdir(parents=True, exist_ok=True)
    results, best = [], (np.inf, None)
    for kp in (80., 120., 160.):
        for kd in (30., 45., 60.):
            for transition in (.05, .1):
                cfg = np.array([1.5, kp, kd, 0., 120., 2., transition])
                neural.fit.CONTROLLER[:] = cfg
                peak, error, soft, maximum, failed = 0., 0., 0., -np.inf, False
                for record in records:
                    t = record["trial"]
                    settled = np.asarray(t["waiting"], dtype=bool).copy()
                    for segment in t["segments"]:
                        stop = int(segment["endIndex"])
                        settled[stop: min(len(settled), stop + 200)] = False
                    for source in (0, 1):
                        o, bad, _, _ = neural.rollout(pair[source], record, source, steps=5)
                        failed = failed or bad
                        peak = max(peak, float(np.max(o[:, 4])))
                        error = max(error, float(np.max(np.abs(o[settled, 0] - t["standingAngle"]))))
                        soft = max(soft, float(np.max(np.abs(boundary_torque(pair[source], o[:, 0])))))
                        maximum = max(maximum, float(np.max(o[:, 0])))
                score = 5 * (max(0, peak - 45) / 20)**2 + (error / .01)**2 + .1 * (soft / 30)**2 + .001 * peak
                if failed:
                    score += 1e4
                row = {"kp": kp, "kd": kd, "transitionSeconds": transition,
                       "abdominalFullSequencePeakNm": peak, "settledAngleErrorRad": error,
                       "softPeakNm": soft, "maximumAngleRad": maximum, "failed": failed, "score": score}
                results.append(row)
                print(row, flush=True)
                if score < best[0]:
                    best = (score, cfg.copy())
    with (out / "standing_search.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)
    # Start the new decoder away from its zero bound so finite-difference
    # calibration can resolve load-dependent tonic compensation.
    neural.NEURAL[:, 9] = .5
    savemat(out / "initial_model.mat", {"parameters": pair, "neuralParameters": neural.NEURAL.copy(),
                                        "controller": best[1][1:]})
    print("Selected common standing controller:", best[1], flush=True)


if __name__ == "__main__":
    main()
