"""Audit fixed parameters, causal reference independence and numerical replay."""

import argparse
import copy
import hashlib
import json
import sys

import numpy as np
from scipy.io import loadmat

from audit_subjects import ROOT, OUT
from fit_subject import many

sys.path.insert(0, str(ROOT / "experiments/legacy_mpc"))
from run_hill_mpc import Plant, simulate_case


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("subject")
    args = parser.parse_args()
    path = OUT / args.subject
    cache = loadmat(path / "neural_models.mat", simplify_cells=True)
    r = loadmat(path / "hill_mpc.mat", simplify_cells=True)["result"]
    unassisted = loadmat(path / "hill_sequence.mat", simplify_cells=True)["result"]
    policy = json.loads((path / "fit_policy.json").read_text())
    meta = json.loads((path / "mpc_evaluation.json").read_text())
    assert (
        meta["modelSHA256"]
        == hashlib.sha256((path / "neural_models.mat").read_bytes()).hexdigest()
    )
    assert (
        policy["preparedSHA256"]
        == hashlib.sha256((path / "prepared.mat").read_bytes()).hexdigest()
    )
    models = many(cache["models"])
    shared = np.asarray(cache["sharedParameterIndices"], dtype=int) - 1
    assert np.array_equal(
        models[0]["parameters"][shared], models[1]["parameters"][shared]
    )
    hashes = [
        hashlib.sha256(
            np.r_[m["parameters"], m["neuralParameters"]].tobytes()
        ).hexdigest()
        for m in models
    ]
    for source in (0, 1):
        assert all(
            c["parameterSHA256"] == hashes[source]
            for c in meta["checks"]
            if c["source"] == ("Experimental activation", "MAP activation")[source]
        )
    for g, base in zip(many(r["groups"]), many(unassisted["groups"])):
        for key in ("q", "torque", "m", "waiting"):
            assert np.array_equal(
                g["sequence"][key], base["sequence"][key], equal_nan=True
            )
        assert np.array_equal(g["inputs"], base["inputs"])
    g = many(r["groups"])[-1]
    s = g["sequence"]
    n = len(s["t"])
    numerics = []
    standing = []
    for source in (0, 1):
        plant = Plant(models[source], cache["controller"], s["standingAngle"])
        a = g["inputs"][source]
        p = g["predictions"][source]
        state, first = plant.init(
            a[0], s["m"][0], s["waiting"][0], np.array([s["q"][0], s["v0"]])
        )
        fine = np.zeros((n, 27))
        fine[0] = first
        for k in range(1, n):
            fine[k] = plant.step(
                state,
                a[k],
                s["m"][k],
                s["waiting"][k],
                p["assistanceTorque"][k],
                s["dt"],
                20,
            )
        qdiff = float(np.max(np.abs(fine[:, 0] - p["q"])))
        tdiff = float(np.max(np.abs(fine[:, 2] - p["torque"])))
        numerics.append(
            {
                "source": source + 1,
                "halfStepAngleDifferenceRad": qdiff,
                "halfStepTorqueDifferenceNm": tdiff,
                "angleWithin0p02Rad": qdiff < 0.02,
                "torqueWithin3Nm": tdiff < 3,
            }
        )
        for initial, velocity in ((1.7, 0.6), (1.25, -0.3)):
            state, out = plant.init(
                np.zeros(2), 0.0, 1.0, np.array([initial, velocity])
            )
            q = [initial]
            for _ in range(3000):
                out = plant.step(state, np.zeros(2), 0.0, 1.0, 0.0, 0.01, 10)
                q.append(out[0])
            standing.append(
                {
                    "source": source + 1,
                    "initialAngleRad": initial,
                    "initialVelocity": velocity,
                    "finalErrorRad": float(abs(out[0] - 1.5)),
                    "maximumAngleRad": float(max(q)),
                    "minimumAngleRad": float(min(q)),
                    "passed": bool(
                        abs(out[0] - 1.5) < 0.01 and max(q) < 2.1 and min(q) > -0.2
                    ),
                }
            )
    # Exercise the active controller while changing ONLY evaluation references.
    seq = copy.deepcopy(s)
    count = min(len(s["t"]), 650)
    for key in (
        "t",
        "q",
        "torque",
        "m",
        "waiting",
        "loadKg",
        "speedCode",
        "stageId",
        "referenceIsFilled",
        "front",
        "a",
    ):
        seq[key] = np.asarray(seq[key])[:count]
    altered = copy.deepcopy(seq)
    altered["q"][1:] = np.linspace(-3, 3, count - 1)
    altered["torque"][:] = 999.0
    a = g["inputs"][0][:count]
    first = simulate_case(models[0], cache["controller"], seq, a, r["options"])[0]
    second = simulate_case(models[0], cache["controller"], altered, a, r["options"])[0]
    leak = float(np.max(np.abs(first - second)))
    assert leak < 1e-10
    result = {
        "subject": args.subject,
        "modelSHA256": meta["modelSHA256"],
        "fixedParameterHashes": hashes,
        "sameVectorAcrossAllCases": True,
        "sharedBodyParameters": True,
        "referenceSamplesUnchanged": True,
        "inputsMatchUnassistedReplay": True,
        "referencePerturbationOutputDifference": leak,
        "halfStepChecks": numerics,
        "perturbedStandingRecovery": standing,
        "allStandingChecksPassed": all(x["passed"] for x in standing),
        "noExternalBalanceTorque": all(
            np.max(np.abs(p["externalBalanceTorque"])) == 0
            for g in many(r["groups"])
            for p in many(g["predictions"])
        ),
    }
    (path / "validation.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
