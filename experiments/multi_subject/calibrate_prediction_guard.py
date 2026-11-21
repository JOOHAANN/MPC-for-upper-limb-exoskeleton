"""Measure causal hold-prediction error on qualified TRAINING signals only."""

import numpy as np
from fit_subject import record, rollout, many
from run_hill_mpc import Plant


def estimate_guard(data, models, cfg, options):
    errors = []
    details = []
    for source, model in enumerate(models):
        source_errors = []
        for seq in many(data["trainSequences"]):
            rec = record(seq, data["generator"])
            _, bad, inputs = rollout(
                model["parameters"], model["neuralParameters"], cfg, rec, source, 10
            )
            if bad:
                raise RuntimeError("Training replay unstable during guard calibration")
            plant = Plant(model, cfg, seq["standingAngle"])
            state, first = plant.init(
                inputs[0],
                seq["m"][0],
                seq["waiting"][0],
                np.array([seq["q"][0], seq["v0"]]),
            )
            states = [state.copy()]
            truth = [first]
            for k in range(1, len(inputs)):
                truth.append(
                    plant.step(
                        state, inputs[k], seq["m"][k], seq["waiting"][k], 0.0, 0.01, 10
                    )
                )
                states.append(state.copy())
            truth = np.asarray(truth)
            for k in range(0, len(inputs) - 42, 10):
                if seq["waiting"][k]:
                    continue
                history = max(0, k - round(options["inputSlopeWindow"] / 0.01))
                plant.input_slope[:] = (
                    np.clip(
                        (inputs[k] - inputs[history]) / ((k - history) * 0.01),
                        -options["inputSlopeLimit"],
                        options["inputSlopeLimit"],
                    )
                    if k > history
                    else 0.0
                )
                forecast = plant.preview(
                    states[k], np.zeros(options["blocks"]), options
                )
                if forecast is None:
                    continue
                future = (
                    k
                    + np.arange(1, options["points"] + 1)
                    * options["predictionDt"]
                    / 0.01
                )
                observed = np.interp(future, np.arange(len(truth)), truth[:, 3])
                source_errors.append(observed - forecast[:, 3])
        if not source_errors:
            raise RuntimeError("No qualified training forecast windows")
        source_errors = np.asarray(source_errors)
        errors.extend(source_errors)
        details.append(
            dict(
                source=source,
                forecastWindows=len(source_errors),
                positiveErrorP99Nm=np.quantile(
                    np.maximum(source_errors, 0), 0.99, axis=0
                ).tolist(),
                maximumPositiveErrorNm=np.max(
                    np.maximum(source_errors, 0), axis=0
                ).tolist(),
            )
        )
    errors = np.asarray(errors)
    # Engineering reserve from the positive tail of training errors, with a
    # small numerical/model floor. This is not a statistical safety guarantee.
    guard = 2.0 + np.quantile(np.maximum(errors, 0), 0.99, axis=0)
    return guard.tolist(), dict(
        trainingOnly=True,
        holdoutSignalsUsed=False,
        sourceDetails=details,
        guardNm=guard.tolist(),
        minimumReserveNm=2.0,
        quantile=0.99,
        simultaneousOrOutOfDistributionGuarantee=False,
    )
