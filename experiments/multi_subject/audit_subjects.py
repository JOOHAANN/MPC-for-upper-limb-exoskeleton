"""Inventory real subject recordings before calibration; never trust nRows.

The acquisition software often writes nRows=0 even when samples exist.
Missing signals, duplicate exports and irregular clocks are reported explicitly.
Raw files are read only. Calibration uses the project's original angle scale.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

ROOT = Path(__file__).resolve().parents[2]
DATA = Path(
    os.environ.get(
        "EXO_SUBJECT_DATA_ROOT",
        ROOT.parent / "YOUHAN LI code and project description/Data",
    )
)
OUT = Path(os.environ.get("EXO_SUBJECT_RESULTS_ROOT", ROOT / "results/multi_subject"))
REQUIRED = (
    "ik.sto",
    "Activations.sto",
    "Torque.sto",
    "MuscleForces.sto",
    "ma_L5_S1_Flex_Ext.sto",
    "lmt.sto",
)
REPRESENTATIVES = (
    "IL_L1_r",
    "IL_L1_l",
    "LTpT_T7_r",
    "LTpT_T7_l",
    "LTpL_L1_r",
    "LTpL_L1_l",
)


def read_sto(path):
    lines = path.read_text(encoding="utf-8-sig", errors="replace").splitlines()
    end = next(i for i, s in enumerate(lines) if s.strip().lower() == "endheader")
    names = lines[end + 1].split()
    body = [s for s in lines[end + 2 :] if s.strip()]
    bad = [i for i, s in enumerate(body) if len(s.split()) != len(names)]
    # A single interrupted final write can be dropped without repairing values.
    # Interior corruption is not silently interpolated or accepted.
    if bad:
        if bad != [len(body) - 1]:
            raise ValueError(f"Malformed interior rows: {bad}")
        body = body[:-1]
    data = (
        np.loadtxt(io.StringIO("\n".join(body)), ndmin=2)
        if body
        else np.empty((0, len(names)))
    )
    return names, data


def clock_data(data):
    """Sort and combine duplicate timestamps before interpolation."""
    data = data[np.isfinite(data).all(axis=1)]
    if not len(data):
        return data
    data = data[np.argsort(data[:, 0], kind="stable")]
    times, inverse, counts = np.unique(
        data[:, 0], return_inverse=True, return_counts=True
    )
    sums = np.zeros((len(times), data.shape[1]))
    np.add.at(sums, inverse, data)
    return sums / counts[:, None]


def smooth_original_scale(raw):
    # MATLAB smoothdata(...,"gaussian",5): sigma=window/5 and truncated support.
    # Re-normalize the Gaussian at endpoints, as MATLAB shrinks the window.
    weights = np.exp(-0.5 * (np.arange(-2, 3) / 1.0) ** 2)
    y = np.convolve(raw, weights, mode="same") / np.convolve(
        np.ones(len(raw)), weights, mode="same"
    )
    return 0.5 + 0.9 * (y - y.min()) / np.ptp(y)


def extract(path):
    tables = {name: read_sto(path / name) for name in REQUIRED}
    names, lmt = tables["lmt.sto"]
    ik = clock_data(tables["ik.sto"][1])
    ac = clock_data(tables["Activations.sto"][1])
    tq = clock_data(tables["Torque.sto"][1])
    forces = clock_data(tables["MuscleForces.sto"][1])
    arms = clock_data(tables["ma_L5_S1_Flex_Ext.sto"][1])
    lmt = clock_data(lmt)
    lo = max(x[0, 0] for x in (lmt, ik, ac, tq, forces, arms))
    hi = min(x[-1, 0] for x in (lmt, ik, ac, tq, forces, arms))
    # Preserve native LMT sample locations for the original five-point smoothing.
    sample = lmt[(lmt[:, 0] >= lo) & (lmt[:, 0] <= hi), 0]
    q_raw = np.interp(
        sample, ik[:, 0], ik[:, tables["ik.sto"][0].index("L5_S1_Flex_Ext")]
    )
    q = smooth_original_scale(q_raw)
    ac_names = tables["Activations.sto"][0]
    f_names, ma_names = (
        tables["MuscleForces.sto"][0],
        tables["ma_L5_S1_Flex_Ext.sto"][0],
    )
    activations, moment_groups = [], []
    for rep in REPRESENTATIVES:
        prefix, side = rep.split("_")[0], rep[-1]
        group = [
            n
            for n in names[1:]
            if n.startswith(prefix + "_") and n.endswith("_" + side)
        ]
        moment = sum(
            np.interp(sample, forces[:, 0], forces[:, f_names.index(n)])
            * np.interp(sample, arms[:, 0], arms[:, ma_names.index(n)])
            for n in group
        )
        moment_groups.append(moment)
        activations.append(np.interp(sample, ac[:, 0], ac[:, ac_names.index(rep)]))
    moments = np.asarray(moment_groups).T
    weights = moments / np.maximum(moments.sum(axis=1, keepdims=True), 1e-10)
    six = np.asarray(activations).T
    back = np.sum(six * weights, axis=1)
    front = np.interp(sample, ac[:, 0], ac[:, ac_names.index("rect_abd_r")])
    return {
        "t": sample - sample[0],
        "q": q,
        "qRaw": q_raw,
        "back": back,
        "front": front,
        "six": six,
        "weights": weights,
        "torque": np.interp(sample, tq[:, 0], tq[:, 1]),
        "absoluteStart": float(sample[0]),
    }


def audit(data=DATA, out=OUT, subjects=None, include_derived=False):
    out.mkdir(parents=True, exist_ok=True)
    rows, signals, hashes = [], {}, {}
    for subject in sorted(data.iterdir()):
        if (
            not subject.is_dir()
            or subject.name == "240528Pilot"
            or (subjects and subject.name not in subjects)
            or (not include_derived and (subject / "derived_dataset.json").is_file())
        ):
            continue
        signals[subject.name] = []
        candidates = sorted(
            {
                p.parent
                for p in subject.rglob("Activations.sto")
                if p.parent.name.startswith("dynamics_stoop")
            }
        )
        for path in candidates:
            match = re.match(
                r"dynamics_stoop_(\d+(?:\.\d+)?)(?:kg)?(?:_(\d+))?", path.name
            )
            load = float(match[1]) if match else -1.0
            speed = int(match[2]) if match and match[2] and int(match[2]) >= 20 else -1
            row = {
                "subject": subject.name,
                "trial": str(path.relative_to(subject)),
                "loadKg": load,
                "speedCode": speed,
                "status": "usable",
                "issues": [],
                "warnings": [],
                "duplicateOf": "",
            }
            digest = hashlib.sha256()
            for name in REQUIRED:
                file = path / name
                if not file.exists():
                    row["issues"].append("missing:" + name)
                    continue
                digest.update(file.read_bytes())
                try:
                    names, a = read_sto(file)
                    file_lines = file.read_text(errors="replace").splitlines()
                    header = next(
                        i
                        for i, s in enumerate(file_lines)
                        if s.strip().lower() == "endheader"
                    )
                    body = [s for s in file_lines[header + 2 :] if s.strip()]
                    malformed = sum(len(s.split()) != len(names) for s in body)
                    if malformed:
                        row["warnings"].append(
                            f"dropped_incomplete_final_row:{name}:{malformed}"
                        )
                    row[name + "Rows"] = len(a)
                    declared = next(
                        (
                            int(s.split("=")[1])
                            for s in file_lines[:header]
                            if s.lower().startswith("nrows=")
                        ),
                        None,
                    )
                    row[name + "DeclaredRows"] = declared
                    row[name + "HeaderRowsMatch"] = declared == len(a)
                    if name == "Activations.sto" and len(a):
                        row["activationOutOfRangeFraction"] = float(
                            np.mean((a[:, 1:] < -1e-6) | (a[:, 1:] > 1 + 1e-6))
                        )
                    if len(a) < 30:
                        row["issues"].append("too_few_samples:" + name)
                    elif not np.isfinite(a).all():
                        row["issues"].append("nonfinite:" + name)
                    else:
                        dt = np.diff(a[:, 0])
                        row[name + "MedianDt"] = float(np.median(dt))
                        row[name + "MaxGap"] = float(np.max(dt))
                        row[name + "NonpositiveDt"] = int(np.sum(dt <= 0))
                        if dt.max() > 0.5:
                            row["warnings"].append("clock_gap_over_0.5s:" + name)
                except Exception as exc:
                    row["issues"].append("parse:" + name + ":" + str(exc))
            row["sha256"] = digest.hexdigest()
            for name in ("emg.sto", "emgFilt.sto"):
                file = path / name
                if not file.exists():
                    row["warnings"].append("missing_optional:" + name)
                    continue
                try:
                    names, signal = read_sto(file)
                    row[name + "Rows"] = len(signal)
                    row[name + "Channels"] = max(0, len(names) - 1)
                    row[name + "NonfiniteFraction"] = (
                        float(np.mean(~np.isfinite(signal))) if signal.size else 0.0
                    )
                    if len(signal) > 2:
                        row[name + "MedianDt"] = float(np.median(np.diff(signal[:, 0])))
                        row[name + "FlatChannels"] = int(
                            np.sum(np.ptp(signal[:, 1:], axis=0) < 1e-8)
                        )
                        row[name + "MaxGap"] = float(np.max(np.diff(signal[:, 0])))
                except Exception as exc:
                    row["warnings"].append("optional_parse:" + name + ":" + str(exc))
            if row["sha256"] in hashes:
                row["duplicateOf"] = hashes[row["sha256"]]
                row["status"] = "duplicate"
            else:
                hashes[row["sha256"]] = subject.name + "/" + row["trial"]
            if row["issues"]:
                row["status"] = "excluded"
            if row["status"] == "usable":
                try:
                    s = extract(path)
                    row.update(
                        durationSeconds=float(s["t"][-1]),
                        backMin=float(s["back"].min()),
                        backMax=float(s["back"].max()),
                        frontMax=float(s["front"].max()),
                        torqueMin=float(s["torque"].min()),
                        torqueMax=float(s["torque"].max()),
                        rawAngleExcursion=float(np.ptp(s["qRaw"])),
                    )
                    peaks = find_peaks(
                        -s["q"],
                        prominence=0.25,
                        distance=max(1, int(0.8 / np.median(np.diff(s["t"])))),
                    )[0]
                    row["flexionValleys"] = len(peaks)
                    row["angleTorqueCorrelation"] = float(
                        np.corrcoef(s["q"], s["torque"])[0, 1]
                    )
                    if row["angleTorqueCorrelation"] > -0.25:
                        row["warnings"].append(
                            "angle_polarity_or_alignment_requires_review"
                        )
                    if (
                        row["rawAngleExcursion"] < 0.015
                        or row["backMax"] < 0.005
                        or len(peaks) < 2
                    ):
                        row["status"] = "review"
                        row["issues"].append("weak_signal_or_less_than_two_flexions")
                    signals[subject.name].append((row, s))
                except Exception as exc:
                    row["status"] = "excluded"
                    row["issues"].append("extract:" + repr(exc))
            print(
                subject.name,
                row["trial"],
                row["status"],
                row.get("durationSeconds"),
                row["issues"],
                flush=True,
            )
            rows.append(row)
        if signals[subject.name]:
            n = len(signals[subject.name])
            fig, axes = plt.subplots(
                n, 3, figsize=(14, max(3, n * 1.45)), squeeze=False
            )
            for i, (row, s) in enumerate(signals[subject.name]):
                for ax, key in zip(axes[i], ("q", "back", "torque")):
                    ax.plot(s["t"], s[key], lw=1)
                    ax.grid(alpha=0.2)
                axes[i, 0].set_ylabel(
                    Path(row["trial"]).name.replace("dynamics_stoop_", ""), fontsize=8
                )
            for ax, title in zip(
                axes[0],
                (
                    "Original scale: theta (rad)",
                    "Experimental back activation",
                    "Recorded L5/S1 torque (N m)",
                ),
            ):
                ax.set_title(title)
            for ax in axes[-1]:
                ax.set_xlabel("Recording time (s)")
            fig.suptitle(subject.name + " | Raw-data audit (no fitted curves)")
            fig.tight_layout()
            fig.savefig(out / (subject.name + "_raw_audit.png"), dpi=130)
            plt.close(fig)
    serializable = [
        {**r, "issues": "; ".join(r["issues"]), "warnings": "; ".join(r["warnings"])}
        for r in rows
    ]
    fields = list(dict.fromkeys(k for r in serializable for k in r))
    with (out / "data_quality.csv").open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(serializable)
    (out / "data_quality.json").write_text(
        json.dumps(rows, indent=2) + "\n", encoding="utf-8"
    )
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DATA)
    parser.add_argument("--output", type=Path, default=OUT)
    parser.add_argument("--subjects", nargs="+")
    parser.add_argument("--include-derived", action="store_true")
    args = parser.parse_args()
    if args.include_derived and not args.subjects:
        parser.error("Select --subjects explicitly when auditing derived datasets.")
    audit(args.data, args.output, args.subjects, args.include_derived)
