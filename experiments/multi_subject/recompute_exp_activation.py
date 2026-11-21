"""Recompute both Exp activation inputs with shapeFactor = +0.1.

Native acquisitions and their old model outputs remain audit references. This
exports a new input profile and paired comparisons, not new muscle torques.
"""

from __future__ import annotations

from collections import Counter
import csv
import json
import xml.etree.ElementTree as ET

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import savemat

from audit_subjects import DATA, REPRESENTATIVES, ROOT, read_sto
from compare_exp_signal_levels import SENSOR_CHANNELS, activation_curve
from quality_gate import gaussian_matlab
from review_exp_one_cycle import LOADS, SPEEDS, SUBJECTS, sha

SHAPE = 0.1
PROCESSED = ROOT / "data/processed/exp_activation_A01"
FIGURES = ROOT / "docs/exp_activation_A01"
REVIEW = ROOT / "docs/exp_one_cycle_review/review.json"
MODEL_CHANNELS = (*REPRESENTATIVES, "rect_abd_r", "rect_abd_l")
FILTERED_CHANNELS = (*SENSOR_CHANNELS, "rectabd_r", "rectabd_l")


def write_sto(path, names, values, title):
    """Keep the native sample clocks and export every muscle/channel."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(
            f"{title}\nnRows={len(values)}\nnColumns={len(names)}\nendheader\n"
        )
        stream.write("\t".join(names) + "\n")
        np.savetxt(stream, values, delimiter="\t", fmt="%.17g")


def subject_configuration(subject, configuration):
    source = DATA / configuration
    tree = ET.parse(source)
    before = tree.getroot()
    shapes = {
        m.findtext("name"): float(m.findtext("shapeFactor"))
        for m in before.iter("muscle")
    }
    changed = 0
    for node in before.iter("shapeFactor"):
        node.text = str(SHAPE)
        changed += 1
    target = PROCESSED / subject / "subject_A01.xml"
    target.parent.mkdir(exist_ok=True, parents=True)
    tree.write(target, encoding="utf-8", xml_declaration=True)
    return shapes, dict(
        subject=subject,
        shapeFactor=SHAPE,
        originalConfiguration=str(source),
        originalSHA256=sha(source),
        configuration=str(target),
        changedShapeFactors=changed,
        otherMuscleParametersChanged=False,
    )


def export_native(row, shapes):
    path = DATA / row["subject"] / row["trial"]
    out = PROCESSED / row["subject"] / row["trial"]
    names, emg = read_sto(path / "emg.sto")
    activation_names, old = read_sto(path / "Activations.sto")
    if not np.isfinite(emg).all():
        raise ValueError("Non-finite native EMG: " + str(path))
    new = activation_curve(emg[:, 1:], SHAPE)
    write_sto(
        out / "Activations.sto",
        names,
        np.column_stack([emg[:, 0], new]),
        "RecomputedActivations_shapeFactor_0.1",
    )
    # Verify the original response at exact clocks; never transform old a again.
    common, ei, ai = np.intersect1d(emg[:, 0], old[:, 0], return_indices=True)
    predicted = activation_curve(emg[ei, 1:], np.array([shapes[n] for n in names[1:]]))
    reference = old[ai][:, [activation_names.index(n) for n in names[1:]]]
    error = float(np.sqrt(np.mean((predicted - reference) ** 2)))
    if len(common) < 30 or error > 1e-6:
        raise ValueError("Original response verification failed: " + str(path))
    sensor_names, sensor = read_sto(path / "emgFilt.sto")
    transformed_sensor = np.column_stack(
        [sensor[:, 0], activation_curve(sensor[:, 1:], SHAPE)]
    )
    write_sto(
        out / "sensor_activation_A01.sto",
        sensor_names,
        transformed_sensor,
        "FilteredEMG_activation_estimate_shapeFactor_0.1",
    )
    return dict(
        subject=row["subject"],
        trial=row["trial"],
        nativeMuscleSamples=len(emg),
        muscles=len(names) - 1,
        filteredSensorSamples=len(sensor),
        originalResponseRMSE=error,
        derivedActivation=str(out / "Activations.sto"),
        derivedSensorActivation=str(out / "sensor_activation_A01.sto"),
        sensorClockConflictsPreservedInNativeExport=True,
        sourceSHA256={
            n: sha(path / n) for n in ("emg.sto", "Activations.sto", "emgFilt.sto")
        },
    )


def selected_signal(row, shapes, *, additional_smoothing_ms=0):
    choice = row["selected"]
    path = DATA / row["subject"] / row["trial"]
    native = choice["source"] == "CEINMS_activation"
    name = "emg.sto" if native else "emgFilt.sto"
    names, table = read_sto(path / name)
    columns = MODEL_CHANNELS if native else FILTERED_CHANNELS
    values = table[:, [names.index(n) for n in columns]]
    clock = table[:, 0]
    if np.any(np.diff(clock) < 0):
        raise ValueError("Backward input clock: " + str(path / name))
    unique, first, counts = np.unique(clock, return_index=True, return_counts=True)
    bad = []
    for i in np.flatnonzero(counts > 1):
        repeated = values[clock == unique[i]]
        if not np.allclose(repeated, repeated[0], atol=1e-10, rtol=1e-9):
            bad.append(unique[i])
    # Preserve acquisition values in native exports; locally exclude ambiguous
    # points for this accepted window, exactly as in the preceding data review.
    clock, values = unique, values[first]
    for point in bad:
        if (
            choice["startEpochSeconds"] <= point + 0.2
            and choice["endEpochSeconds"] >= point - 0.2
        ):
            raise ValueError("Selected window overlaps ambiguous sensor samples")
    allowed = ~np.isin(clock, bad)
    clock, values = clock[allowed], values[allowed]
    ik = read_sto(path / "ik.sto")[1]
    lo, hi = max(clock[0], ik[0, 0]), min(clock[-1], ik[-1, 0])
    fine = np.arange(0, hi - lo, 0.001)
    old_shapes = np.array([shapes[n] for n in MODEL_CHANNELS])
    traces = {}
    for key, signal in (
        ("excitation", values),
        ("oldActivation", activation_curve(values, old_shapes)),
        ("activation", activation_curve(values, SHAPE)),
    ):
        sampled = np.column_stack([np.interp(fine, clock - lo, v) for v in signal.T])
        filtered = gaussian_matlab(gaussian_matlab(sampled, 5), 300)
        if additional_smoothing_ms:
            filtered = gaussian_matlab(filtered, int(additional_smoothing_ms))
        traces[key] = filtered[::10]
    times = lo + fine[::10]
    mask = (times >= choice["startEpochSeconds"]) & (times <= choice["endEpochSeconds"])
    if mask.sum() < 100:
        raise ValueError("Insufficient selected-window coverage")
    traces = {key: value[mask] for key, value in traces.items()}
    traces["timeEpochSeconds"] = times[mask]
    traces["timeSeconds"] = times[mask] - row["originEpochSeconds"]
    return traces, "model-input EMG" if native else "full filtered EMG"


def window_metrics(row, traces, source):
    old, new, u = [
        traces[k][:, :6].mean(1) for k in ("oldActivation", "activation", "excitation")
    ]
    a0, a1 = [float(np.percentile(v, 95)) for v in (old, new)]
    return dict(
        subject=row["subject"],
        loadKg=row["loadKg"],
        speedCode=row["speedCode"],
        trial=row["trial"],
        inputSource=source,
        oldShapeFactor=float(OLD_SHAPES[row["subject"]][REPRESENTATIVES[0]]),
        newShapeFactor=SHAPE,
        startSecondsFromEmg=row["selected"]["startSecondsFromEmg"],
        endSecondsFromEmg=row["selected"]["endSecondsFromEmg"],
        oldBackP95=a0,
        newBackP95=a1,
        backP95ReductionPercent=100 * (1 - a1 / a0),
        newBackInputP95Ratio=a1 / float(np.percentile(u, 95)),
        backWaveformCorrelation=float(np.corrcoef(old, new)[0, 1]),
        oldAbdomenP95=float(np.percentile(traces["oldActivation"][:, 6], 95)),
        newAbdomenP95=float(np.percentile(traces["activation"][:, 6], 95)),
        minimumNewChannelActivation=float(traces["activation"].min()),
        maximumNewChannelActivation=float(traces["activation"].max()),
        physiologicalRangePassed=bool(
            traces["activation"].min() >= -1e-6
            and traces["activation"].max() <= 1 + 1e-6
        ),
        negativeChannelSampleFraction=float(np.mean(traces["activation"] < -1e-6)),
        muscleOutputCoverage=row["selected"][
            "existingMusculoskeletalTimeCoverageFraction"
        ],
        earlierSmoothnessWarnings=row["selected"]["smoothnessWarnings"],
        musclesOrTorqueRecomputed=False,
    )


def plot_overview(rows, traces):
    fig, axes = plt.subplots(2, 2, figsize=(12, 6.8), sharex="col")
    for j, key in enumerate(((SUBJECTS[0], 5, 45), (SUBJECTS[1], 10, 45))):
        row = next(
            r for r in rows if (r["subject"], r["loadKg"], r["speedCode"]) == key
        )
        signal = traces[key]
        for i, label in enumerate(("Back: six-channel mean", "Abdomen: right channel")):
            ax = axes[i, j]
            for field, color, style, text in (
                (
                    "oldActivation",
                    "#a4473e",
                    "--",
                    f"Old A = {OLD_SHAPES[key[0]][REPRESENTATIVES[0]]:.4f}",
                ),
                ("excitation", ".5", ":", "EMG input"),
                ("activation", "#1469b0", "-", "New A = +0.1"),
            ):
                y = signal[field][:, :6].mean(1) if i == 0 else signal[field][:, 6]
                ax.plot(
                    signal["timeSeconds"], y, color=color, ls=style, lw=1.6, label=text
                )
            ax.set_ylabel(label + "\nActivation / input (0-1)")
            ax.grid(alpha=0.16)
            ax.set_ylim(bottom=0)
            if i == 0:
                ax.set_title(f"{key[0]} | {key[1]} kg / {key[2]}")
                ax.legend(fontsize=8)
            else:
                ax.set_xlabel("Time from filtered-EMG start (s)")
    fig.suptitle(
        "Both Exp groups: recomputed activation with shapeFactor = +0.1", fontsize=13
    )
    fig.tight_layout(rect=(0, 0.045, 1, 0.94))
    fig.text(
        0.5,
        0.012,
        "Same EMG, same selected window and identical smoothing. Change applied per channel before averaging. These are activation results, not new torque predictions.",
        ha="center",
        fontsize=8,
    )
    fig.savefig(FIGURES / "activation_before_after.png", dpi=155)
    plt.close(fig)


def plot_grid(subject, rows, traces):
    fig, axes = plt.subplots(5, 5, figsize=(13, 11.5), sharey=True)
    for i, load in enumerate(LOADS):
        for j, speed in enumerate(SPEEDS):
            key, ax = (subject, load, speed), axes[i, j]
            row = next(
                r for r in rows if (r["subject"], r["loadKg"], r["speedCode"]) == key
            )
            if key in traces:
                signal = traces[key]
                t = signal["timeSeconds"] - signal["timeSeconds"][0]
                ax.plot(
                    t,
                    signal["oldActivation"][:, :6].mean(1),
                    "--",
                    color="#a4473e",
                    lw=0.8,
                    alpha=0.7,
                )
                ax.plot(t, signal["activation"][:, :6].mean(1), color="#1469b0", lw=1.2)
                tag = "N" if row["selected"]["source"] == "CEINMS_activation" else "F"
                ax.set_title(f"{load:g} kg / {speed} | {tag}", fontsize=9)
                ax.set_xlabel("Time in kept window (s)", fontsize=7)
            else:
                ax.set_facecolor(".93")
                ax.text(
                    0.5,
                    0.5,
                    "No clean window\nfrom earlier review",
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                    fontsize=8,
                )
                ax.set_title(f"{load:g} kg / {speed}", fontsize=9)
            ax.grid(alpha=0.15)
    fig.suptitle(
        subject
        + " | uniform activation mapping A = +0.1\nBlue: new activation. Dashed red: old mapping, recomputed from the same input.",
        fontsize=12,
    )
    fig.text(
        0.5,
        0.015,
        "N: native model-input EMG. F: full filtered EMG estimate; paired muscle outputs may be missing. Windows and prior quality flags are unchanged.",
        ha="center",
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.94), h_pad=1)
    fig.savefig(FIGURES / (subject + "_A01_grid.png"), dpi=150)
    plt.close(fig)


OLD_SHAPES = {}


def main():
    PROCESSED.mkdir(exist_ok=True, parents=True)
    FIGURES.mkdir(exist_ok=True, parents=True)
    review = json.loads(REVIEW.read_text())
    verification = json.loads(
        (REVIEW.parent / "signal_level_verification.json").read_text()
    )
    profiles = []
    for subject in SUBJECTS:
        original = next(
            r["configuration"]
            for r in verification["records"]
            if r["subject"] == subject
        )
        OLD_SHAPES[subject], profile = subject_configuration(subject, original)
        profiles.append(profile)
    exports, skipped = [], []
    for row in review["basicRecords"] + review["supplementalNumberedRecords"]:
        if "error" in row:
            skipped.append(
                dict(subject=row["subject"], trial=row["trial"], reason=row["error"])
            )
            continue
        exports.append(export_native(row, OLD_SHAPES[row["subject"]]))
    metrics, traces = [], {}
    for row in review["conditions"]:
        if row["selected"] is None:
            continue
        signal, source = selected_signal(row, OLD_SHAPES[row["subject"]])
        key = (row["subject"], row["loadKg"], row["speedCode"])
        traces[key] = signal
        entry = window_metrics(row, signal, source)
        filename = f"{row['loadKg']:g}kg_{row['speedCode']}.mat"
        destination = PROCESSED / row["subject"] / "selected_windows" / filename
        destination.parent.mkdir(parents=True, exist_ok=True)
        savemat(
            destination,
            {
                **signal,
                "channels": np.array(MODEL_CHANNELS, dtype=object),
                "shapeFactor": SHAPE,
                "inputSource": source,
                "physiologicalRangePassed": entry["physiologicalRangePassed"],
            },
            do_compression=True,
        )
        entry["selectedWindowMAT"] = str(destination)
        metrics.append(entry)
    summary = []
    for subject in SUBJECTS:
        m = [r for r in metrics if r["subject"] == subject]
        summary.append(
            dict(
                subject=subject,
                oldShapeFactor=OLD_SHAPES[subject][REPRESENTATIVES[0]],
                newShapeFactor=SHAPE,
                exportedRecords=sum(e["subject"] == subject for e in exports),
                keptConditions=len(m),
                missingConditions=25 - len(m),
                medianBackP95ReductionPercent=float(
                    np.median([r["backP95ReductionPercent"] for r in m])
                ),
                medianBackWaveformCorrelation=float(
                    np.median([r["backWaveformCorrelation"] for r in m])
                ),
                medianNewBackInputP95Ratio=float(
                    np.median([r["newBackInputP95Ratio"] for r in m])
                ),
                medianAbdomenP95ReductionPercent=float(
                    np.median(
                        [100 * (1 - r["newAbdomenP95"] / r["oldAbdomenP95"]) for r in m]
                    )
                ),
                physiologicalRangeReviewConditions=sum(
                    not r["physiologicalRangePassed"] for r in m
                ),
                sourceCounts=dict(Counter(r["inputSource"] for r in m)),
            )
        )
        plot_grid(subject, review["conditions"], traces)
    plot_overview(review["conditions"], traces)
    payload = dict(
        shapeFactor=SHAPE,
        originalDataModified=False,
        profiles=profiles,
        mapping="a=expm1(A*u)/expm1(A), evaluated per channel before interpolation, smoothing and averaging",
        filtering="1 ms interpolation; 5 ms then 300 ms Gaussian; display/export at 100 Hz",
        windowsFixedFromReview=True,
        outputClippingApplied=False,
        MAPRetrained=False,
        BSERetrained=False,
        MPCRerun=False,
        torqueOrMuscleForcesRecomputed=False,
        nativeSensorDuplicatesPreserved=True,
        muscleOutputWarning="Old MuscleForces/Torque belong to the old activation setting. They are not new outputs for A=0.1; do not treat this activation-only directory as a complete neuromusculoskeletal replay.",
        summary=summary,
        exports=exports,
        selectedWindows=metrics,
        skipped=skipped,
        upstreamReviewSHA256=sha(REVIEW),
    )
    (PROCESSED / "profile.json").write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n"
    )
    (FIGURES / "results.json").write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n"
    )
    with (FIGURES / "condition_metrics.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as file:
        writer = csv.DictWriter(file, fieldnames=list(metrics[0]))
        writer.writeheader()
        writer.writerows(metrics)
    print(json.dumps(summary, indent=2), flush=True)
    print("New activation profile:", PROCESSED / "profile.json", flush=True)


if __name__ == "__main__":
    main()
