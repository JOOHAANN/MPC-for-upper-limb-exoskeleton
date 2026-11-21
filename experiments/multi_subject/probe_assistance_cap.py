"""Development stress probe for personalized actuator authority."""

import argparse
import json
import numpy as np
from scipy.io import loadmat
from audit_subjects import OUT
from fit_subject import many, record, rollout
from control_subject import simulate_case, metrics


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cap", type=float, required=True)
    p.add_argument("--checkpoint", action="store_true")
    a = p.parse_args()
    folder = OUT / "240723Pilot"
    cache = loadmat(folder / "neural_models.mat", simplify_cells=True)
    replay = loadmat(folder / "hill_sequence.mat", simplify_cells=True)["result"]
    meta = json.loads((folder / "mpc_evaluation.json").read_text())
    opt = {**meta["options"], "capNm": a.cap, "angleGovernorUpperRad": 1.65}
    g = many(replay["groups"])[2]
    model = dict(many(cache["models"])[0])
    inputs = g["inputs"][0]
    if a.checkpoint:
        checkpoint = loadmat(folder / "checkpoint.mat", simplify_cells=True)
        model["parameters"] = checkpoint["parameters"][0]
        model["neuralParameters"] = checkpoint["neuralParameters"][0]
        generator = model["generator"]
        rec = record(g["sequence"], generator)
        inputs = rollout(
            model["parameters"],
            model["neuralParameters"],
            cache["controller"],
            rec,
            0,
            10,
        )[2]
    actual, nominal, flags, times = simulate_case(
        model,
        cache["controller"],
        g["sequence"],
        inputs,
        opt,
    )
    rows = metrics(actual, nominal, g["sequence"], flags, 0)
    result = {
        "capNm": a.cap,
        "qmax": float(actual[:, 0].max()),
        "qmin": float(actual[:, 0].min()),
        "frontPeakNm": float(actual[:, 4].max()),
        "maximumCaseRatio": max(r["abdomenBackRatio"] for r in rows),
        "backExcessNm": max(r["excessPeakNm"] for r in rows),
        "rows": rows,
        "role": "Development stress evaluation; not an independent held-out claim",
    }
    suffix = "_checkpoint" if a.checkpoint else ""
    (folder / f"cap_probe_{a.cap:g}{suffix}.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
