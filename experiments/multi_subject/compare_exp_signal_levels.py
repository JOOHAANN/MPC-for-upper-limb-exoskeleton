"""Explain EMG/activation levels using unchanged, same-time native exports.

Saved shape factors are checked against the exported activation; no parameters
are fitted and no sensor values are rescaled to match a target amplitude.
"""

from __future__ import annotations

from collections import Counter
import json
import xml.etree.ElementTree as ET

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from audit_subjects import DATA, REPRESENTATIVES, read_sto
from quality_gate import gaussian_matlab
from review_exp_one_cycle import OUTPUT, sha

SENSOR_CHANNELS = (
    "iliocost_r",
    "iliocost_l",
    "longthor_r",
    "longthor_l",
    "longlumb_r",
    "longlumb_l",
)


def activation_curve(excitation, shape):
    """Evaluate the saved exponential mapping before channel averaging."""
    return np.expm1(np.asarray(shape) * excitation) / np.expm1(shape)


def saved_shapes(subject):
    choices = []
    for path in sorted((DATA / subject).glob("*Calibrated.xml")):
        muscles = {
            m.findtext("name"): m.findtext("shapeFactor")
            for m in ET.parse(path).iter("muscle")
        }
        if all(muscles.get(name) is not None for name in REPRESENTATIVES):
            choices.append(
                (path, np.array([float(muscles[n]) for n in REPRESENTATIVES]))
            )
    return choices


def native_pair(path):
    en, emg = read_sto(path / "emg.sto")
    an, act = read_sto(path / "Activations.sto")
    # Exact-time joins avoid fitting a delay or interpolating the verification.
    time, ei, ai = np.intersect1d(emg[:, 0], act[:, 0], return_indices=True)
    if len(time) < 30:
        raise ValueError("Insufficient exact-time EMG/activation pairs: " + str(path))
    excitation = emg[ei][:, [en.index(n) for n in REPRESENTATIVES]]
    activation = act[ai][:, [an.index(n) for n in REPRESENTATIVES]]
    return time, excitation, activation


def paired_example(row, model, shape):
    path = DATA / row["subject"] / row["trial"]
    time, excitation, activation = native_pair(path)
    predicted = activation_curve(excitation, shape)
    names, sensor = read_sto(path / "emgFilt.sto")
    lo, hi = max(time[0], sensor[0, 0]), min(time[-1], sensor[-1, 0])
    fine = np.arange(0, hi - lo, 0.001)

    # Identical physical smoothing for both stages, without peak normalization.
    def display(clock, values):
        values = np.column_stack([np.interp(fine, clock - lo, v) for v in values.T])
        return gaussian_matlab(gaussian_matlab(values, 5), 300)[::10].mean(1)

    sensor_values = sensor[:, [names.index(n) for n in SENSOR_CHANNELS]]
    displayed = [
        display(sensor[:, 0], sensor_values),
        display(time, excitation),
        display(time, activation),
    ]
    selected = row["selected"]
    relative = lo + fine[::10] - row["originEpochSeconds"]
    keep = (lo + fine[::10] >= selected["startEpochSeconds"]) & (
        lo + fine[::10] <= selected["endEpochSeconds"]
    )
    fig = plt.figure(figsize=(12, 6.6))
    layout = fig.add_gridspec(2, 2, height_ratios=[1.2, 1])
    wave = fig.add_subplot(layout[0, :])
    for values, color, label, style in zip(
        displayed,
        ("#c17a1b", "#666666", "#1469b0"),
        (
            "Filtered EMG (E)",
            "Model-input EMG, same timestamps",
            "Exported activation (A)",
        ),
        ("-", "--", "-"),
    ):
        wave.plot(
            relative[keep], values[keep], color=color, label=label, ls=style, lw=1.6
        )
    wave.set(
        xlabel="Time from filtered-EMG start (s)",
        ylabel="Mean of six back channels",
        ylim=(0, 1),
        title="Same trial and same interval; identical smoothing; no amplitude rescaling",
    )
    wave.legend(fontsize=8, ncol=3, loc="upper right")
    response, proof = fig.add_subplot(layout[1, 0]), fig.add_subplot(layout[1, 1])
    u = np.linspace(0, 1, 301)
    response.plot(
        u,
        activation_curve(u, shape[0]),
        color="#1469b0",
        lw=1.7,
        label=f"Saved shape factor = {shape[0]:.4f}",
    )
    response.plot(u, u, color=".5", ls="--", label="Equal numerical values")
    response.set(
        xlabel="EMG excitation u",
        ylabel="Activation a",
        xlim=(0, 1),
        ylim=(0, 1),
        title="a = [exp(A u) - 1] / [exp(A) - 1]",
    )
    response.legend(fontsize=8)
    proof.scatter(
        predicted[::4].ravel(), activation[::4].ravel(), s=4, alpha=0.3, color="#1469b0"
    )
    proof.plot([0, 1], [0, 1], "--", color=".4", lw=1)
    error = float(np.sqrt(np.mean((predicted - activation) ** 2)))
    proof.set(
        xlabel="Activation computed from saved parameters",
        ylabel="Exported activation",
        xlim=(0, 1),
        ylim=(0, 1),
        title=f"Native exact-time check: RMSE {error:.2e}",
    )
    for axis in (wave, response, proof):
        axis.grid(alpha=0.16)
    fig.suptitle(
        "240730Exp | 10 kg / 45: why activation is higher than filtered EMG",
        fontsize=13,
    )
    fig.tight_layout(rect=(0, 0.035, 1, 0.95))
    fig.text(
        0.5,
        0.012,
        "Both colors are experimental-data stages. Neither is MAP output. Parameters are read from the saved XML; no fit is performed.",
        ha="center",
        fontsize=8,
    )
    fig.savefig(OUTPUT / "240730Exp_EMG_activation_levels.png", dpi=155)
    plt.close(fig)
    window = (time >= selected["startEpochSeconds"]) & (
        time <= selected["endEpochSeconds"]
    )
    return dict(
        subject=row["subject"],
        trial=row["trial"],
        configuration=str(model.relative_to(DATA)),
        selectedIntervalSecondsFromEmg=[
            selected["startSecondsFromEmg"],
            selected["endSecondsFromEmg"],
        ],
        exactTimeActivationRMSE=error,
        modelInputMeanP95=float(np.percentile(excitation[window].mean(1), 95)),
        exportedActivationMeanP95=float(np.percentile(activation[window].mean(1), 95)),
    )


def main():
    review = json.loads((OUTPUT / "review.json").read_text())
    records, configs = [], {}
    for row in review["basicRecords"]:
        subject = row["subject"]
        if subject not in configs:
            configs[subject] = saved_shapes(subject)
        path = DATA / subject / row["trial"]
        time, excitation, activation = native_pair(path)
        errors = [
            (
                float(
                    np.sqrt(
                        np.mean((activation_curve(excitation, shape) - activation) ** 2)
                    )
                ),
                model,
                shape,
            )
            for model, shape in configs[subject]
        ]
        error, model, shape = min(errors, key=lambda item: item[0])
        records.append(
            dict(
                subject=subject,
                trial=row["trial"],
                exactTimeSamples=len(time),
                configuration=str(model.relative_to(DATA)),
                shapeFactors=shape.tolist(),
                activationRMSE=error,
                matchesSavedParameters=error < 1e-6,
                sourceSHA256={
                    str(p.relative_to(DATA)): sha(p)
                    for p in (path / "emg.sto", path / "Activations.sto", model)
                },
            )
        )
    example = next(
        r
        for r in review["conditions"]
        if r["subject"] == "240730Exp" and r["loadKg"] == 10 and r["speedCode"] == 45
    )
    selected_model = next(
        r
        for r in records
        if r["subject"] == example["subject"] and r["trial"] == example["trial"]
    )
    example_summary = paired_example(
        example,
        DATA / selected_model["configuration"],
        np.array(selected_model["shapeFactors"]),
    )
    payload = dict(
        equation="a = expm1(A*u) / expm1(A)",
        description="Saved exponential activation conversion verified at native exact timestamps; this is not a physiological calibration or a MAP/BSE model fit.",
        originalDataModified=False,
        parametersFitted=False,
        bothColorsAreExperimentalDataStages=True,
        MAPOutputShown=False,
        sensorAndModelInputMayDifferSlightlyInTimestampSampling=True,
        applyNonlinearityBeforeChannelAveraging=True,
        records=records,
        example=example_summary,
        summary=dict(
            records=len(records),
            matches=sum(r["matchesSavedParameters"] for r in records),
            maximumActivationRMSE=max(r["activationRMSE"] for r in records),
            configurationCounts=dict(Counter(r["configuration"] for r in records)),
        ),
    )
    (OUTPUT / "signal_level_verification.json").write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n"
    )
    print(
        json.dumps(
            {"summary": payload["summary"], "example": example_summary}, indent=2
        )
    )


if __name__ == "__main__":
    main()
