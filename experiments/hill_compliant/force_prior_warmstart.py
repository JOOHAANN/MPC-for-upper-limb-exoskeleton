"""Explore lower strength without changing any predicted torque curves.

This changes a parameter initialization only. The subsequent forward replay
and joint recalibration must pass the same nine-condition checks.
"""
import argparse
from pathlib import Path
import numpy as np
from scipy.io import loadmat, savemat
from fit_hill_physiology import LOW, HIGH


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    package = loadmat(args.model, simplify_cells=True)
    pair = package["parameters"].copy()
    neural = package["neuralParameters"].copy()
    for source in (0, 1):
        old_capacity, old_gamma = pair[source, 0], pair[source, 5]
        capacity = np.clip(10000 * pair[source, 20], LOW[0], HIGH[0])
        gamma = np.clip(old_gamma + np.log(old_capacity / capacity) / np.log(.3), LOW[5], HIGH[5])
        previous_suppressed_torque = old_capacity * ((1 - neural[source, 6]) * .1)**old_gamma
        new_mix = (previous_suppressed_torque / capacity)**(1 / gamma) / .1
        neural[source, 6] = np.clip(1 - new_mix, .7, .999999)
        pair[source, 0], pair[source, 5] = capacity, gamma
        print("Warm source", source + 1, "capacity/gamma/blend", capacity, gamma, neural[source, 6])
    savemat(args.out, {"parameters": pair, "neuralParameters": neural,
                      "description": "Engineering force-prior warm start, not a fitted or accepted model"})


if __name__ == "__main__":
    main()
