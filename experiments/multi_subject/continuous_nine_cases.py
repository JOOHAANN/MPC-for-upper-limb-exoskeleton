"""Replay a true continuous 5/10/15 kg sequence with the selected controller.

State is initialized once. The subject's nine-condition controller weights,
activation forecaster and uncertainty calibration are reused without tuning.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import numpy as np
from scipy.io import loadmat, savemat

from control_nine_cases import (
    BASE,
    load_subject,
    many,
    attach_planned_inputs,
    digest,
    threshold_schedule,
    simulate_case,
    evaluate,
    prediction_dict,
    structs,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("subject")
    parser.add_argument("--candidate", required=True)
    args = parser.parse_args()
    if args.subject == "240528Pilot":
        raise ValueError(
            "The original 528 evaluations already contain three continuous sequences."
        )
    cache, models, _, model_path, replay_path = load_subject(args.subject)
    report_path = BASE / args.subject / "nmpc" / args.candidate
    report = json.loads((report_path / "evaluation.json").read_text())
    assert report["signature"]["modelSHA256"] == digest(model_path)
    groups = many(loadmat(replay_path, simplify_cells=True)["result"]["groups"])
    group = groups[-1]
    assert len(many(group["sequence"]["segments"])) == 3
    attach_planned_inputs([group], models)
    options = report["options"]
    seq = group["sequence"]

    def run(source):
        current = dict(options)
        if options["inputForecast"] in ("arx", "hybrid"):
            source_model = options["sourcePredictionModels"][source]
            current["inputForecast"] = source_model["inputForecast"]
            current["predictionErrorGuardNm"] = source_model["guardNm"]
            if "predictionErrorGuardRelative" in options:
                current["predictionErrorGuardRelative"] = source_model["relativeGuard"]
        sequence = dict(
            seq,
            plannedMAPInputs=group["plannedMAPInputs"][source],
            activationSourceIsGenerated=bool(source == 1),
        )
        threshold, meta = threshold_schedule(
            sequence, many(group["predictions"])[source]
        )
        sequence["backTorqueThresholdNm"] = threshold
        actual, nominal, flags, timings = simulate_case(
            models[source],
            cache["controller"],
            sequence,
            np.asarray(group["inputs"])[source],
            current,
        )
        assert actual[:, 24].max() <= 80 + 1e-8
        pred = prediction_dict(actual, models[source]["parameters"])
        pred.update(
            solverExitFlag=flags, solveSeconds=timings, backTorqueThresholdNm=threshold
        )
        rows = evaluate(actual, nominal, sequence, source, meta)
        print(
            args.subject,
            "continuous",
            source + 1,
            "peak margins",
            [round(r["peakMarginNm"], 2) for r in rows],
            flush=True,
        )
        return pred, prediction_dict(nominal, models[source]["parameters"]), rows

    with ThreadPoolExecutor(max_workers=2) as pool:
        runs = list(pool.map(run, (0, 1)))
    result = dict(
        groups=structs(
            [
                dict(
                    group=int(group["group"]),
                    sequence=seq,
                    predictions=structs([v[0] for v in runs]),
                    nominal=structs([v[1] for v in runs]),
                    inputs=group["inputs"],
                )
            ]
        ),
        models=structs(models),
        controller=cache["controller"],
        options=options,
    )
    savemat(
        report_path / "continuous_mpc.mat",
        dict(result=result),
        long_field_names=True,
        do_compression=True,
    )
    (report_path / "continuous_evaluation.json").write_text(
        json.dumps(
            dict(
                subject=args.subject,
                modelSHA256=digest(model_path),
                sameControllerAsNineCases=True,
                stateResetsBetweenActions=False,
                metrics=[r for v in runs for r in v[2]],
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
