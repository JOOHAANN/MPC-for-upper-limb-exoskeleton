"""Escape saturated intent-decoder basins before complete Hill calibration.

Every seed is fitted jointly to all nine actions. Physical parameters remain a
warm start here; fit_subject subsequently re-optimizes every eligible parameter.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import shutil

import numpy as np
from scipy.io import loadmat, savemat
from scipy.optimize import least_squares

import fit_subject as fit
from control_nine_cases import BASE


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("subject")
    parser.add_argument(
        "--output", type=Path, default=BASE / "candidates/response_multistart"
    )
    parser.add_argument("--minimum-feedback", type=float, default=0.85)
    parser.add_argument("--minimum-stiffness", type=float, default=70.0)
    parser.add_argument("--max-nfev", type=int, default=30)
    args = parser.parse_args()
    folder = BASE / args.subject
    output = args.output / args.subject
    output.mkdir(parents=True, exist_ok=True)
    data = loadmat(folder / "prepared.mat", simplify_cells=True)
    original = loadmat(folder / "neural_models.mat", simplify_cells=True)
    models = fit.many(original["models"])
    cfg = original["controller"]
    records = [
        fit.record(s, data["generator"], 2) for s in fit.many(data["sequences"])[:9]
    ]
    fit.NLOW[6] = args.minimum_feedback
    fit.NLOW[3] = args.minimum_stiffness
    fit.NHIGH[8] = 1.5
    fit.NHIGH[5] = 80.0
    fit.SHAPE_TAIL_WEIGHT = 8.0
    fit.LANDMARK_WEIGHT = 1.0
    fit.DIRECTION_WEIGHT = 4.0
    fit.VALLEY_DEPTH_WEIGHT = 3.0
    low, span = fit.NLOW, fit.NHIGH - fit.NLOW

    def run(source):
        p = models[source]["parameters"]
        current = models[source]["neuralParameters"]
        seeds = [
            current,
            np.array([0.9, 0.5, 0.10, 100, 20, 40, 0.95, 0, 0.5, 0.2]),
            np.array([1.2, 0.6, 0.15, 150, 30, 40, 0.95, 0, 0.8, 0]),
            np.array([1.5, 1.0, 0.10, 100, 20, 40, 0.95, 0, 0.8, 0.5]),
        ]
        best = (np.inf, None)
        history = []
        for index, seed in enumerate(seeds):
            start = np.clip((seed - low) / span, 1e-5, 1 - 1e-5)
            incumbent = [np.inf, start.copy()]

            def residual(z):
                v = fit.residual(p, low + span * z, cfg, records, source, p, current, 6)
                value = float(v @ v)
                if value < incumbent[0]:
                    incumbent[:] = [value, z.copy()]
                return v

            least_squares(
                residual,
                start,
                bounds=(np.zeros(10), np.ones(10)),
                max_nfev=args.max_nfev,
                diff_step=2e-5,
                ftol=5e-4,
            )
            value, z = incumbent
            rr = low + span * z
            history.append(dict(seed=index, loss=value, neuralParameters=rr.tolist()))
            print(
                args.subject,
                "neural restart",
                source + 1,
                index,
                "loss",
                round(value, 3),
                flush=True,
            )
            if value < best[0]:
                best = (value, rr.copy())
        return best[1], history

    with ThreadPoolExecutor(max_workers=2) as pool:
        runs = list(pool.map(run, (0, 1)))
    for s in (0, 1):
        models[s]["neuralParameters"] = runs[s][0]
    savemat(
        output / "neural_seed.mat",
        dict(subject=args.subject, models=fit.structs(models), controller=cfg),
        long_field_names=True,
        do_compression=True,
    )
    for name in ("prepared.mat", "quality_gate.json", "preparation.json"):
        shutil.copy2(folder / name, output / name)
    (output / "neural_restart_audit.json").write_text(
        json.dumps(
            dict(
                subject=args.subject,
                minimumFeedback=args.minimum_feedback,
                minimumStiffness=args.minimum_stiffness,
                history=[v[1] for v in runs],
                oneNeuralVectorPerSourceAcrossNineCases=True,
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
