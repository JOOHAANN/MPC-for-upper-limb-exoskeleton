"""Check whether the two earlier Pilot recordings inflate EMG activation.

The +0.1 response is a numerical hypothesis checked against unchanged exports,
not a new fitted physiological parameter. Record repeats remain separate files.
"""

from __future__ import annotations

import csv
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from audit_subjects import DATA, ROOT, read_sto
from compare_exp_signal_levels import SENSOR_CHANNELS, activation_curve, native_pair
from quality_gate import gaussian_matlab
from review_exp_one_cycle import sha

OUT = ROOT / "docs/activation_level_audit"
SUBJECTS = ("240528Pilot", "240617Pilot")


def compare_trial(subject, path):
    record = dict(subject=subject, trial=str(path.relative_to(DATA / subject)))
    record["sourceSHA256"] = {
        n: sha(path / n) for n in ("emg.sto", "Activations.sto", "emgFilt.sto")
    }
    try:
        time, excitation, activation = native_pair(path)
    except ValueError as error:
        record.update(
            status="too_short_for_curve_comparison",
            error=str(error),
            samples={
                n: len(read_sto(path / n)[1]) for n in ("emg.sto", "Activations.sto")
            },
        )
        return record
    predicted = activation_curve(excitation, 0.1)
    names, sensor = read_sto(path / "emgFilt.sto")
    common = (time >= sensor[0, 0]) & (time <= sensor[-1, 0])
    sensor_at_time = np.column_stack(
        [
            np.interp(time[common], sensor[:, 0], sensor[:, names.index(n)])
            for n in SENSOR_CHANNELS
        ]
    )
    input_p95 = float(np.percentile(excitation.mean(1), 95))
    act_p95 = float(np.percentile(activation.mean(1), 95))
    record.update(
        status="checked",
        exactTimeSamples=len(time),
        checkedShapeFactor=0.1,
        activationRMSE=float(np.sqrt(np.mean((predicted - activation) ** 2))),
        maximumActivationInputDifference=float(np.max(abs(activation - excitation))),
        excitationP95=input_p95,
        activationP95=act_p95,
        activationInputP95Ratio=act_p95 / max(input_p95, 1e-12),
        sensorModelInputInterpolationRMSE=float(
            np.sqrt(np.mean((sensor_at_time - excitation[common]) ** 2))
        ),
        signalValuesClipped=False,
        inputDelaysFitted=False,
    )
    return record


def display_example(ax, subject, trial):
    path = DATA / subject / trial
    time, excitation, activation = native_pair(path)
    names, sensor = read_sto(path / "emgFilt.sto")
    origin = sensor[0, 0]
    lo, hi = max(time[0], sensor[0, 0]), min(time[-1], sensor[-1, 0])
    fine = np.arange(0, hi - lo, 0.001)
    for clock, values, color, style, label in (
        (
            sensor[:, 0],
            sensor[:, [names.index(n) for n in SENSOR_CHANNELS]],
            "#c17a1b",
            "-",
            "Filtered EMG",
        ),
        (time, excitation, ".45", "--", "Model-input EMG"),
        (time, activation, "#1469b0", "-", "Activation"),
    ):
        fine_values = np.column_stack(
            [np.interp(fine, clock - lo, v) for v in values.T]
        )
        y = gaussian_matlab(gaussian_matlab(fine_values, 5), 300)[::10].mean(1)
        ax.plot(lo + fine[::10] - origin, y, color=color, ls=style, lw=1.4, label=label)
    ax.set(
        title=subject + " | " + path.name.removeprefix("dynamics_stoop_"),
        xlabel="Time from filtered-EMG start (s)",
        ylabel="Mean of six back channels",
    )
    ax.legend(fontsize=8)


def plot(records):
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.5))
    display_example(axes[0, 0], "240528Pilot", "dynamics_stoop_5kg_45")
    display_example(axes[0, 1], "240617Pilot", "CalibrationData/dynamics_stoop_10kg_33")
    response = axes[1, 0]
    u = np.linspace(0, 1, 301)
    response.plot(u, u, "--", color=".6", label="a = u")
    for shape, label, color in (
        (0.1, "Both earlier Pilots: A = +0.1", "#1469b0"),
        (-2.13108325449385, "240725Exp: A = -2.1311", "#8b639c"),
        (-2.77886425366985, "240730Exp: A = -2.7789", "#c17a1b"),
    ):
        response.plot(u, activation_curve(u, shape), label=label, color=color, lw=1.6)
    response.set(
        xlabel="Single-channel model input u",
        ylabel="Activation a",
        title="Different saved/exported activation responses",
        xlim=(0, 1),
        ylim=(0, 1),
    )
    response.legend(fontsize=8)
    ratio_axis = axes[1, 1]
    for i, subject in enumerate(SUBJECTS):
        good = [
            r for r in records if r["subject"] == subject and r["status"] == "checked"
        ]
        ratios = np.array([r["activationInputP95Ratio"] for r in good])
        ratio_axis.scatter(
            np.linspace(i - 0.14, i + 0.14, len(good)),
            ratios,
            color="#1469b0",
            s=14,
            alpha=0.65,
        )
        ratio_axis.plot(
            [i - 0.18, i + 0.18], [np.median(ratios)] * 2, color="#993b32", lw=2
        )
    ratio_axis.axhline(1, color=".5", ls="--", label="Equal P95 values")
    ratio_axis.set(
        xticks=[0, 1],
        xticklabels=SUBJECTS,
        ylabel="P95(activation mean) / P95(EMG mean)",
        title="Every checked stoop export; each dot is one file",
        ylim=(0.94, 1.01),
    )
    ratio_axis.legend(fontsize=8)
    for ax in axes.flat:
        ax.grid(alpha=0.16)
    fig.suptitle(
        "Earlier Pilots: activation stays close to EMG; no Exp-like amplitude uplift",
        fontsize=13,
    )
    fig.tight_layout(rect=(0, 0.045, 1, 0.95))
    fig.text(
        0.5,
        0.015,
        "Exact-time verification uses native emg.sto / Activations.sto. Waveforms use identical smoothing; no peak rescaling. Four 1-2-sample files are excluded.",
        ha="center",
        fontsize=8,
    )
    fig.savefig(OUT / "pilot_EMG_activation_comparison.png", dpi=155)
    plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True, parents=True)
    records = []
    for subject in SUBJECTS:
        for file in sorted((DATA / subject).rglob("Activations.sto")):
            if file.parent.name.startswith("dynamics_stoop"):
                records.append(compare_trial(subject, file.parent))
    summary = []
    for subject in SUBJECTS:
        rows = [r for r in records if r["subject"] == subject]
        checked = [r for r in rows if r["status"] == "checked"]
        summary.append(
            dict(
                subject=subject,
                records=len(rows),
                checked=len(checked),
                tooShort=len(rows) - len(checked),
                matchesShapePointOne=sum(r["activationRMSE"] < 1e-6 for r in checked),
                maximumActivationRMSE=max(r["activationRMSE"] for r in checked),
                maximumActivationInputDifference=max(
                    r["maximumActivationInputDifference"] for r in checked
                ),
                medianActivationInputP95Ratio=float(
                    np.median([r["activationInputP95Ratio"] for r in checked])
                ),
                medianSensorInputInterpolationRMSE=float(
                    np.median([r["sensorModelInputInterpolationRMSE"] for r in checked])
                ),
            )
        )
    payload = dict(
        summary=summary,
        records=records,
        equation="a = expm1(A*u)/expm1(A)",
        checkedShapeFactor=0.1,
        parametersRetrained=False,
        originalDataModified=False,
        scope="All stoop export files, including repeats; counts are not independent experimental conditions.",
        provenance="240617Pilot/subjectUncalib.xml stores +0.1. No same-folder calibrated XML was found for 240528Pilot; its +0.1 response is verified from exports, without asserting the identity of an unavailable runtime XML.",
    )
    (OUT / "pilot_level_audit.json").write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n"
    )
    fields = [
        "subject",
        "trial",
        "status",
        "exactTimeSamples",
        "activationRMSE",
        "maximumActivationInputDifference",
        "excitationP95",
        "activationP95",
        "activationInputP95Ratio",
        "sensorModelInputInterpolationRMSE",
        "error",
    ]
    with (OUT / "pilot_level_audit.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as file:
        writer = csv.DictWriter(file, fieldnames=fields)
        writer.writeheader()
        writer.writerows({k: r.get(k) for k in fields} for r in records)
    plot(records)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
