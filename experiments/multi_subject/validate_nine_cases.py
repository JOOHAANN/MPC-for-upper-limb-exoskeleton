"""Audit frozen models, torque accounting and time-step sensitivity.

The finer-step replay applies saved actuator commands to the original Hill
plant. It never re-optimizes controls or clips an output to satisfy a limit.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat

from control_nine_cases import BASE, ROOT, digest, load_subject, many, Plant


def validate(subject, candidate):
    cache, models, groups, model_path, replay_path = load_subject(subject)
    folder = BASE / subject / "nmpc" / candidate
    report = json.loads((folder / "evaluation.json").read_text())
    saved = loadmat(folder / "hill_mpc.mat", simplify_cells=True)["result"]
    assert report["signature"]["modelSHA256"] == digest(model_path)
    assert report["signature"]["replaySHA256"] == digest(replay_path)
    assert report["fullNineConditions"] and report["caseCountPerSource"] == 9
    for model, used in zip(models, many(saved["models"])):
        np.testing.assert_array_equal(model["parameters"], used["parameters"])
        np.testing.assert_array_equal(
            model["neuralParameters"], used["neuralParameters"]
        )
    rows = []
    for group in groups:
        seq = group["sequence"]
        for source, model in enumerate(models):
            data = np.load(folder / f"group{int(group['group'])}_source{source+1}.npz")
            actual, nominal, inputs = (data[k] for k in ("actual", "nominal", "inputs"))
            assert np.isfinite(actual).all() and np.isfinite(nominal).all()
            assert actual[:, 24].min() >= -1e-8
            assert actual[:, 24].max() <= 80 + 1e-8
            assert np.abs(np.diff(actual[:, 24])).max() / seq["dt"] <= 200 + 1e-6
            accounting = float(np.abs(actual[:, 2] - actual[:, 3] + actual[:, 4]).max())
            assert accounting < 1e-8
            equilibrium = float(np.abs(actual[:, 15]).max())
            momentum = float(np.abs(actual[1:, 16]).max())
            assert equilibrium < 1e-4 and momentum < 1.0
            for segment in many(seq["segments"]):
                ids = np.arange(
                    int(segment["startIndex"]) - 1, int(segment["endIndex"])
                )
                np.testing.assert_allclose(
                    data["threshold"][ids], 0.75 * nominal[ids, 3].max(), atol=1e-5
                )
            plant = Plant(model, cache["controller"], seq["standingAngle"])
            state, first = plant.init(
                inputs[0],
                seq["m"][0],
                seq["waiting"][0],
                np.array([seq["q"][0], seq["v0"]]),
            )
            fine = [first]
            for k in range(1, len(inputs)):
                fine.append(
                    plant.step(
                        state,
                        inputs[k],
                        seq["m"][k],
                        seq["waiting"][k],
                        actual[k, 24],
                        seq["dt"],
                        20,
                    )
                )
            fine = np.asarray(fine)
            q_delta = float(np.abs(fine[:, 0] - actual[:, 0]).max())
            back_delta = float(np.abs(fine[:, 3] - actual[:, 3]).max())
            # Convergence tolerances are reported, not used to edit curves.
            assert q_delta < 0.025, (subject, group["group"], source, q_delta)
            assert back_delta < 3.0, (subject, group["group"], source, back_delta)
            fine_margins = []
            for segment in many(seq["segments"]):
                ids = np.arange(
                    int(segment["startIndex"]) - 1, int(segment["endIndex"])
                )
                fine_margins.append(
                    float(data["threshold"][ids][0] - fine[ids, 3].max())
                )
            rows.append(
                dict(
                    subject=subject,
                    group=int(group["group"]),
                    source=source + 1,
                    torqueIdentityMaximumNm=accounting,
                    fiberEquilibriumMaximum=equilibrium,
                    momentumResidualMaximumNm=momentum,
                    finerStepThetaMaximumDifferenceRad=q_delta,
                    finerStepBackMaximumDifferenceNm=back_delta,
                    appliedAssistanceMaximumNm=float(actual[:, 24].max()),
                    finerStepPeakMarginsNm=fine_margins,
                    endStandingErrorRad=float(
                        abs(actual[-1, 0] - seq["standingAngle"])
                    ),
                )
            )
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subject")
    parser.add_argument("--candidate")
    args = parser.parse_args()
    selected = (
        {args.subject: args.candidate}
        if args.subject
        else json.loads((BASE / "selected_nmpc.json").read_text())
    )
    rows = []
    for subject, candidate in selected.items():
        rows.extend(validate(subject, candidate))
    protected = json.loads((BASE / "preserved_models.json").read_text())["files"]
    for filename, expected in protected.items():
        assert digest(Path(filename)) == expected, (
            "A protected earlier model changed: " + filename
        )
    result = dict(
        passed=True,
        rows=rows,
        protectedFilesUnchanged=len(protected),
        finerStepSubsteps=20,
        nominalSubsteps=10,
        usesSavedUnmodifiedActuatorCommands=True,
    )
    path = BASE / (
        "validation_" + args.subject + ".json" if args.subject else "validation.json"
    )
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            dict(
                passed=True,
                simulations=len(rows),
                protectedFilesUnchanged=len(protected),
                maximumThetaDifferenceRad=max(
                    r["finerStepThetaMaximumDifferenceRad"] for r in rows
                ),
                maximumBackDifferenceNm=max(
                    r["finerStepBackMaximumDifferenceNm"] for r in rows
                ),
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
