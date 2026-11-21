"""Plot direct native-file evidence for incomplete CEINMS export coverage.

No signal alignment, amplitude normalization, endpoint padding or model fitting
is performed. The time origin is subtracted equally from every source clock.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

from audit_subjects import DATA, REPRESENTATIVES, ROOT, read_sto

OUTPUT = ROOT / "docs/exp_truncation_evidence"
SOURCE_FILES = (
    "emgFilt.sto",
    "ik.sto",
    "emg.sto",
    "Activations.sto",
    "lmt.sto",
    "MuscleForces.sto",
    "Torque.sto",
)


def load_example(subject, load):
    path = DATA / subject / "CalibrationData" / f"dynamics_stoop_{load:g}kg_45"
    tables = {name: read_sto(path / name) for name in SOURCE_FILES}
    origin = float(tables["emgFilt.sto"][1][0, 0])
    endpoints = []
    for name, (names, data) in tables.items():
        lines = (path / name).read_text(encoding="utf-8-sig").splitlines()
        header = next(
            i for i, line in enumerate(lines) if line.strip().lower() == "endheader"
        )
        body = [
            (i + 1, line)
            for i, line in enumerate(lines)
            if i > header + 1 and line.strip()
        ]
        malformed = [i for i, line in body if len(line.split()) != len(names)]
        endpoints.append(
            dict(
                subject=subject,
                loadKg=load,
                folder=path.name,
                file=name,
                path=str(path / name),
                sha256=hashlib.sha256((path / name).read_bytes()).hexdigest(),
                sampleRows=len(data),
                malformedPhysicalRows=malformed,
                lastPhysicalDataLine=body[-1][0],
                firstNativeTimestamp=float(data[0, 0]),
                lastNativeTimestamp=float(data[-1, 0]),
                firstTimeFromEmgSeconds=float(data[0, 0] - origin),
                lastTimeFromEmgSeconds=float(data[-1, 0] - origin),
            )
        )
    return path, tables, origin, endpoints


def native_curve(ax, data, values, origin, **kwargs):
    """Break rather than interpolate across source gaps in the display."""
    values = np.asarray(values).copy()
    values[np.r_[False, np.diff(data[:, 0]) > 0.25]] = np.nan
    ax.plot(data[:, 0] - origin, values, **kwargs)


def detailed_figure(subject, load, tables, origin, output):
    emg_names, emg = tables["emgFilt.sto"]
    ac_names, activation = tables["Activations.sto"]
    ik_names, ik = tables["ik.sto"]
    _, torque = tables["Torque.sto"]
    back_names = [
        prefix + side
        for prefix in ("iliocost", "longlumb", "longthor")
        for side in ("_r", "_l")
    ]
    sensor = emg[:, [emg_names.index(n) for n in back_names]].mean(axis=1)
    drive = activation[:, [ac_names.index(n) for n in REPRESENTATIVES]].mean(axis=1)
    stop = float(activation[-1, 0] - origin)
    fig, axes = plt.subplots(
        4, 2, figsize=(12.6, 9.3), sharex="col", gridspec_kw={"width_ratios": [1.55, 1]}
    )
    for column, interval in enumerate(((0, 32), (12, 21))):
        native_curve(axes[0, column], emg, sensor, origin, color="#2166ac", lw=1)
        native_curve(
            axes[1, column], activation, drive, origin, color="#df8a16", lw=1.2
        )
        native_curve(
            axes[2, column],
            ik,
            ik[:, ik_names.index("L5_S1_Flex_Ext")],
            origin,
            color="#585858",
            lw=1.1,
        )
        native_curve(
            axes[3, column], torque, torque[:, 1], origin, color="#8c6bb1", lw=1.2
        )
        for row in range(4):
            ax = axes[row, column]
            ax.set_xlim(*interval)
            ax.axvline(stop, color="#bc302d", ls="--", lw=1.4)
            ax.grid(alpha=0.2)
            if row in (1, 3):
                last = (
                    activation[-1, 0] - origin if row == 1 else torque[-1, 0] - origin
                )
                ax.axvspan(last, interval[1], color="#c84a43", alpha=0.10)
                ax.text(
                    (max(last, interval[0]) + interval[1]) / 2,
                    0.72,
                    "NO EXPORTED DATA",
                    color="#9c312c",
                    transform=ax.get_xaxis_transform(),
                    fontsize=8,
                    ha="center",
                )
        axes[0, column].set_title(
            (
                "First 32 s of this recording"
                if column == 0
                else "Close-up of the cut activity"
            ),
            fontsize=11,
            pad=27,
        )
        axes[0, column].annotate(
            f"Activation ends at {stop:.3f} s",
            (stop, 1.02),
            xycoords=("data", "axes fraction"),
            ha="right",
            color="#a42a27",
            fontsize=9,
        )
        axes[3, column].set_xlabel("Time from first filtered-EMG sample (s)")
    # These windows are the already published EMG bout boundaries, shown as
    # screening intervals rather than asserted independent lifting repetitions.
    axes[0, 0].axvspan(13.37, 19.15, color="#e7b75b", alpha=0.18)
    axes[0, 0].axvspan(21.14, 26.54, color="#4b9b73", alpha=0.13)
    axes[0, 0].text(
        16.26,
        0.96,
        "Activity continues",
        transform=axes[0, 0].get_xaxis_transform(),
        fontsize=8,
        ha="center",
        va="top",
    )
    axes[0, 0].text(
        23.84,
        0.96,
        "Later activity",
        transform=axes[0, 0].get_xaxis_transform(),
        fontsize=8,
        ha="center",
        va="top",
    )
    labels = (
        "Filtered EMG\nback-channel mean",
        "Activation\nback-channel mean",
        "IK L5/S1\n(stored values)",
        "CEINMS Torque\n(N m)",
    )
    for row, label in enumerate(labels):
        axes[row, 0].set_ylabel(label, fontsize=10)
    fig.suptitle(
        f"{subject} | {load:g} kg | folder suffix 45\nThe activation/torque export stops while EMG and movement continue",
        fontsize=13,
    )
    fig.text(
        0.5,
        0.012,
        "Native timestamps; no fitted time shift. No values are filled after an export ends. IK is shown without rescaling.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.035, 1, 0.94), h_pad=0.65)
    fig.savefig(output / "240725Exp_5kg45_cut_detail.png", dpi=155)
    fig.savefig(output / "240725Exp_5kg45_cut_detail.svg")
    plt.close(fig)


def timeline_figure(examples, output):
    fig, axes = plt.subplots(1, len(examples), figsize=(12, 4.6), sharey=True)
    for ax, (subject, load, _, tables, origin, endpoints) in zip(axes, examples):
        for i, row in enumerate(endpoints):
            color = "#327cac" if row["file"] in ("emgFilt.sto", "ik.sto") else "#d48121"
            ax.barh(
                i,
                row["lastTimeFromEmgSeconds"] - row["firstTimeFromEmgSeconds"],
                left=row["firstTimeFromEmgSeconds"],
                height=0.58,
                color=color,
            )
            ax.text(
                row["lastTimeFromEmgSeconds"] + 0.5,
                i,
                f"{row['lastTimeFromEmgSeconds']:.3f} s",
                va="center",
                fontsize=8,
            )
        cut = float(tables["Activations.sto"][1][-1, 0] - origin)
        ax.axvline(cut, color="#bc302d", ls="--", lw=1.2)
        ax.set_title(f"{subject} | {load:g} kg | suffix 45")
        ax.set_xlabel("Time from filtered-EMG start (s)")
        ax.set_xlim(0, max(r["lastTimeFromEmgSeconds"] for r in endpoints) + 10)
        ax.grid(axis="x", alpha=0.2)
    axes[0].set_yticks(range(len(SOURCE_FILES)), SOURCE_FILES)
    axes[0].invert_yaxis()
    fig.suptitle("Native file coverage: EMG / IK continue after CEINMS exports end")
    fig.legend(
        handles=[
            Patch(color="#327cac", label="Longer measured-signal files"),
            Patch(color="#d48121", label="Shorter CEINMS exports"),
        ],
        loc="lower center",
        ncol=2,
        frameon=False,
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.08, 1, 0.94))
    fig.savefig(output / "native_file_endpoints.png", dpi=160)
    fig.savefig(output / "native_file_endpoints.svg")
    plt.close(fig)


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    examples = [
        (subject, load, *load_example(subject, load))
        for subject, load in (("240725Exp", 5), ("240730Exp", 12.5))
    ]
    subject, load, path, tables, origin, _ = examples[0]
    detailed_figure(subject, load, tables, origin, OUTPUT)
    timeline_figure(examples, OUTPUT)
    rows = [r for e in examples for r in e[-1]]
    with (OUTPUT / "native_file_endpoints.csv").open(
        "w", newline="", encoding="utf-8-sig"
    ) as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    evidence = dict(
        auditDate="2026-10-20",
        timePolicy="All original absolute timestamps minus the same first emgFilt timestamp per trial; zero fitted shifts.",
        signalPolicy="No amplitude normalization, extra smoothing, IK rescaling, extrapolation or missing activation reconstruction.",
        originNativeTimestamp=origin,
        demonstratedActivityWindowsSeconds=[[13.37, 19.15], [21.14, 26.54]],
        windowsAreEmgScreeningIntervalsNotIndependentTrials=True,
        exportStopCauseEstablished=False,
        examples=rows,
    )
    (OUTPUT / "evidence.json").write_text(json.dumps(evidence, indent=2) + "\n")
    print(
        json.dumps(
            [
                {
                    k: r[k]
                    for k in (
                        "subject",
                        "loadKg",
                        "file",
                        "sampleRows",
                        "lastPhysicalDataLine",
                        "lastTimeFromEmgSeconds",
                    )
                }
                for r in rows
            ],
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
