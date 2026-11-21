"""Tune one standing handoff policy across all loads and activation sources.

No muscle or mechanical parameter is changed. Candidate standing gains and
the excitation blend time are scored by the actual Hill simulation, including
recovery depth, back torque, settling, abdominal braking and movement shape.
"""

import csv
import ctypes
import json
from pathlib import Path
import sys

import numpy as np
from scipy.io import loadmat, savemat
from scipy.signal import find_peaks

from run_hill_mpc import ROOT


def main():
    out = ROOT / "results/hill_recovery_tuning"
    out.mkdir(parents=True, exist_ok=True)
    cache = loadmat(
        ROOT / "data/processed/hill_compliant/neural_models.mat", simplify_cells=True
    )
    replay = loadmat(
        ROOT / "results/hill_neural_sequence/hill_sequence.mat", simplify_cells=True
    )["result"]
    lib = ctypes.CDLL(str(ROOT / "experiments/hill_compliant/hill_neural_core.dll"))
    array = np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS")
    lib.hill_neural_integrate.argtypes = [array] * 4 + [
        ctypes.c_int,
        array,
        ctypes.c_double,
        ctypes.c_int,
        array,
        array,
        array,
    ]
    records = []
    for g in replay["groups"]:
        seq = g["sequence"]
        for source in (0, 1):
            records.append(
                (
                    seq,
                    cache["models"][source],
                    np.asfortranarray(g["inputs"][source]),
                    source,
                )
            )
    rows = []
    for kp in (15.0, 30.0, 50.0, 80.0):
        for kd in (10.0, 20.0, 30.0, 45.0):
            for transition in (0.1, 0.3, 0.6):
                cfg = np.array([1.5, kp, kd, 0.0, 120.0, 2.0, transition])
                peaks, tail_errors, settled_errors, delays = [], [], [], []
                failed = False
                abdomen, qmax, valleys = 0.0, 0.0, 0
                for seq, model, inputs, source in records:
                    n = len(inputs)
                    o = np.zeros((n, 24), order="F")
                    failed |= bool(
                        lib.hill_neural_integrate(
                            np.asfortranarray(model["parameters"]),
                            inputs,
                            np.asfortranarray(seq["m"]),
                            np.asfortranarray(seq["waiting"], dtype=float),
                            n,
                            np.asfortranarray([seq["q"][0], seq["v0"]]),
                            0.01,
                            5,
                            np.asfortranarray(cfg),
                            np.asfortranarray(model["neuralParameters"]),
                            o,
                        )
                    )
                    abdomen = max(abdomen, float(o[:, 4].max()))
                    qmax = max(qmax, float(o[:, 0].max()))
                    for seg in seq["segments"]:
                        end, start = int(seg["endIndex"]), int(seg["startIndex"]) - 1
                        stop = min(n, end + 100)
                        peaks.append(float(o[end:stop, 3].max()))
                        tail_errors.append(float(np.mean((o[end:stop, 0] - 1.5) ** 2)))
                        settled_errors.append(
                            float(
                                np.max(
                                    np.abs(
                                        o[min(n - 1, end + 150) : min(n, end + 250), 0]
                                        - 1.5
                                    )
                                )
                            )
                        )
                        valleys += (
                            len(
                                find_peaks(
                                    -o[start:end, 0], prominence=0.12, distance=50
                                )[0]
                            )
                            == 2
                        )
                score = (
                    np.mean(np.square(np.maximum(np.array(peaks) - 15, 0) / 15))
                    + 0.25 * (max(peaks) / 35) ** 2
                )
                score += (
                    0.5 * np.mean(tail_errors) / 0.1**2
                    + 2 * (max(settled_errors) / 0.025) ** 2
                )
                score += (
                    20 * (18 - valleys)
                    + max(0, abdomen - 20) ** 2 / 100
                    + max(0, qmax - 1.65) ** 2 / 0.01
                )
                if failed:
                    score += 1e5
                row = dict(
                    kp=kp,
                    kd=kd,
                    transitionSeconds=transition,
                    score=float(score),
                    tailBackPeakNm=max(peaks),
                    tailBackMeanPeakNm=float(np.mean(peaks)),
                    recoveryAngleRMSE=float(np.sqrt(np.mean(tail_errors))),
                    settledErrorRad=max(settled_errors),
                    abdomenPeakNm=abdomen,
                    maximumAngleRad=qmax,
                    twoValleyCases=valleys,
                    failed=failed,
                )
                rows.append(row)
                print(json.dumps(row), flush=True)
    rows.sort(key=lambda r: r["score"])
    with (out / "candidates.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    best = rows[0]
    cache["controller"] = np.array(
        [best["kp"], best["kd"], 0.0, 120.0, 2.0, best["transitionSeconds"]]
    )
    savemat(
        out / "recovery_models.mat",
        {k: v for k, v in cache.items() if not k.startswith("__")},
        do_compression=True,
        long_field_names=True,
    )
    print("SELECTED", json.dumps(best), flush=True)


if __name__ == "__main__":
    main()
