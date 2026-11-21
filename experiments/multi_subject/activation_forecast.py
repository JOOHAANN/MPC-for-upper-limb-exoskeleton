"""Subject-specific causal residual prediction around a trained MAP generator.

One regularized, direct multi-horizon ARX matrix serves all conditions. Inputs
are only past activation residuals. No condition IDs, IK or torque labels are
features. Leave-one-condition-out activation prediction selects regularization.
"""

from __future__ import annotations

import json
from pathlib import Path
import numpy as np

LAGS = np.arange(0, 41, 2, dtype=int)


def features(observed, planned, indices):
    delta = np.asarray(observed) - np.asarray(planned)
    past = np.maximum(0, np.asarray(indices)[:, None] - LAGS[None, :])
    return delta[past].reshape(len(indices), -1)


def matrices(group, dt, points, prediction_dt):
    observed = np.asarray(group["inputs"])[0]
    planned = np.asarray(group["plannedMAPInputs"])[0]
    delta = observed - planned
    count = len(observed)
    indices = np.arange(0, count - int(np.ceil(points * prediction_dt / dt)) - 1, 2)
    x = features(observed, planned, indices)
    query = indices[:, None] + np.arange(1, points + 1)[None, :] * prediction_dt / dt
    target = np.stack(
        [
            np.interp(query.reshape(-1), np.arange(count), delta[:, j]).reshape(
                query.shape
            )
            for j in (0, 1)
        ],
        axis=2,
    )
    return x, target.reshape(len(indices), -1)


def fit(x, y, penalty):
    scale = np.maximum(np.std(x, axis=0), 0.005)
    z = x / scale
    coefficient = np.linalg.solve(
        z.T @ z / len(z) + penalty * np.eye(z.shape[1]), z.T @ y / len(z)
    )
    return scale, coefficient


def fit_predictor(groups, options, signature, path):
    if path.exists():
        saved = json.loads(path.read_text())
        if saved["signature"] == signature:
            return saved["predictor"]
    # Keep all repeated display groups out of fitting; the caller supplies the
    # nine distinct calibration conditions once (or three original sequences).
    data = []
    for group in groups:
        seq = group["sequence"]
        segments = np.asarray(seq["segments"], dtype=object).reshape(-1).tolist()
        for segment in segments:
            left = max(0, int(segment["startIndex"]) - 1 - 60)
            right = min(len(seq["q"]), int(segment["endIndex"]) + 60)
            inputs = np.stack([np.asarray(group["inputs"][s]) for s in (0, 1)])
            subset = dict(
                inputs=inputs[:, left:right],
                plannedMAPInputs=np.asarray(group["plannedMAPInputs"])[:, left:right],
            )
            data.append(
                matrices(subset, seq["dt"], options["points"], options["predictionDt"])
            )
    scores = []
    channel_weight = np.tile([1.0, 0.2], options["points"])
    for penalty in (1e-5, 1e-4, 1e-3, 0.01, 0.1, 1.0):
        fold = []
        for heldout in range(len(data)):
            x = np.vstack([v[0] for j, v in enumerate(data) if j != heldout])
            y = np.vstack([v[1] for j, v in enumerate(data) if j != heldout])
            scale, coefficient = fit(x, y, penalty)
            test_x, test_y = data[heldout]
            error = test_x / scale @ coefficient - test_y
            fold.append(float(np.mean(error**2 * channel_weight)))
        scores.append(
            dict(
                penalty=penalty,
                leaveOneConditionOutRMSE=float(np.sqrt(np.mean(fold))),
                foldRMSE=np.sqrt(fold).tolist(),
            )
        )
    penalty = min(scores, key=lambda r: r["leaveOneConditionOutRMSE"])["penalty"]
    x = np.vstack([v[0] for v in data])
    y = np.vstack([v[1] for v in data])
    scale, coefficient = fit(x, y, penalty)
    predictor = dict(
        lags=LAGS.tolist(),
        scale=scale.tolist(),
        coefficients=coefficient.tolist(),
        points=options["points"],
        predictionDt=options["predictionDt"],
        penalty=penalty,
        usesOnlyPastActivationResiduals=True,
        usesAngleOrTorqueLabels=False,
        baselineMAPResidualRMSE=float(np.sqrt(np.mean(y**2 * channel_weight))),
        calibrationRMSE=float(
            np.sqrt(np.mean((x / scale @ coefficient - y) ** 2 * channel_weight))
        ),
        crossValidation=scores,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(dict(signature=signature, predictor=predictor), indent=2) + "\n"
    )
    best = min(scores, key=lambda r: r["leaveOneConditionOutRMSE"])
    print(
        "ACTIVATION FORECAST",
        len(data),
        "conditions; MAP residual RMSE",
        predictor["baselineMAPResidualRMSE"],
        "LOCO residual RMSE",
        best["leaveOneConditionOutRMSE"],
        "penalty",
        penalty,
        flush=True,
    )
    return predictor
