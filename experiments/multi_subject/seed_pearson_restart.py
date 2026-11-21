"""Screen reproducible neural-policy restarts on unchanged nine-case data."""

from concurrent.futures import ThreadPoolExecutor
import json

import numpy as np
from scipy.io import loadmat, savemat
from scipy.stats import qmc

from recalibrate_bse_pearson import OUT, evaluate, summary, fit, load_subject


def merit(stats, source):
    torque_target = (0.82, 0.77)[source]
    angle_target = (0.735, 0.635)[source]
    return (
        (30 * max(0, torque_target - stats["meanTorqueR"])) ** 2
        + (40 * max(0, angle_target - stats["meanThetaR"])) ** 2
        + (12 * max(0, stats["maximumAbdomenBackRatio"] - 0.35)) ** 2
        + (10 * max(0, stats["boundaryPeakNm"] - 25) / 25) ** 2
        + (10 * max(0, stats["maximumTheta"] - 1.95) / 0.2) ** 2
        + (10 * max(0, -0.1 - stats["minimumTheta"]) / 0.2) ** 2
        + 1000 * stats["failed"]
    )


def main():
    subject = "240730Exp_A01"
    source = 0
    package, models, groups, _, _ = load_subject(subject)
    original = np.r_[models[source]["parameters"], models[source]["neuralParameters"]]
    folder = OUT / subject / f"source{source+1}"
    current = loadmat(folder / "checkpoint.mat", simplify_cells=True)["vector"]
    records = [
        fit.record(g["sequence"], models[source]["generator"], 2) for g in groups
    ]
    seeds = [current.copy(), original.copy()]
    for fraction in (0.25, 0.5, 0.75, 1.0):
        for rho in (0.2, 0.4, 0.6, 0.8):
            vector = current.copy()
            vector[32:] = (1 - fraction) * current[32:] + fraction * original[32:]
            vector[38] = rho
            seeds.append(vector)
    draws = qmc.LatinHypercube(d=6, seed=730).random(72)
    for draw in draws:
        vector = current.copy()
        vector[32] = 0.8 + 2.1 * draw[0]  # Intent amplitude.
        vector[33] = 0.35 + 1.5 * draw[1]  # Intent power.
        vector[38] = 0.15 + 0.7 * draw[2]  # Proprioceptive/direct-input blend.
        vector[40] = 0.2 + 1.3 * draw[3]  # Causal lead; never future EMG.
        vector[41] = 2.5 * draw[4]  # Load tonic compensation.
        vector[5] = 0.65 + 1.6 * draw[5]  # Back activation exponent.
        seeds.append(vector)

    def score(vector):
        rows, _ = evaluate(
            vector[:32], vector[32:], package["controller"], records, source, 6
        )
        stats = summary(rows)
        return merit(stats, source), vector, stats

    with ThreadPoolExecutor(max_workers=3) as pool:
        scored = list(pool.map(score, seeds))
    scored.sort(key=lambda item: item[0])
    savemat(folder / "multistart_seed.mat", dict(vector=scored[0][1]))
    (folder / "restart_screen.json").write_text(
        json.dumps(
            [
                dict(merit=cost, vector=vector.tolist(), metrics=stats)
                for cost, vector, stats in scored[:12]
            ],
            indent=2,
        )
        + "\n"
    )
    print(
        json.dumps(
            dict(seedCount=len(seeds), best=scored[0][2], merit=scored[0][0]), indent=2
        )
    )


if __name__ == "__main__":
    main()
