"""Plot saved Hill/MPC trajectories; no fitting, simulation, or signal editing.

Each activation source keeps its selected parameters across all nine cases.
Only movement panels use local time; their original samples are unchanged.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/hill_mpc"
BUNDLE_NAME = "optimized_model_and_mpc_recovery.pdf"
CURRENT_OPTIONS = {}
LABELS = ("Experimental activation", "MAP activation")
SLUGS = ("experimental", "map")
RED, AMBER, BLUE, GREEN, PURPLE = "#cc3434", "#c8841c", "#1469b0", "#23844a", "#8054a6"
NOMINAL = {"color": AMBER, "lw": 2.6, "ls": (0, (5, 2))}
ASSISTED = {"color": BLUE, "lw": 1.35}
REFERENCE = {"color": RED, "lw": 1.25, "ls": "--"}
# Selected replay settings from run_hill_mpc.py; these penalize BACK torque.
BACK_THRESHOLDS_NM = {5: 100.0, 10: 140.0, 15: 180.0}
THRESHOLD_STYLE = {"color": "#333333", "lw": 1.35, "ls": "--"}


def back_threshold_lines(ax, seq):
    """Show load targets and the separately declared recovery target."""
    for index, segment in enumerate(seq["segments"]):
        selected = segment_slice(segment)
        start, end = seq["t"][selected.start], seq["t"][selected.stop - 1]
        threshold = BACK_THRESHOLDS_NM[int(segment["loadKg"])]
        ax.plot(
            [start, end],
            [threshold, threshold],
            **THRESHOLD_STYLE,
            label="MPC back-torque soft threshold" if index == 0 else "_nolegend_",
        )
        ax.annotate(
            f"{threshold:g} N·m",
            ((start + end) / 2, threshold),
            xytext=(0, 4),
            textcoords="offset points",
            ha="center",
            fontsize=8,
            bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none", "pad": 1},
        )
    recovery = CURRENT_OPTIONS.get("recoveryBackTargetNm", 0.0)
    if recovery > 0:
        ax.plot(
            seq["t"],
            np.where(np.asarray(seq["waiting"], dtype=bool), recovery, np.nan),
            color="#777777",
            ls="--",
            lw=0.9,
            label=f"Unloaded recovery target: {recovery:g} N·m",
        )


def segment_slice(segment):
    """MATLAB stores one-based, inclusive segment endpoints."""
    return slice(int(segment["startIndex"]) - 1, int(segment["endIndex"]))


def rmse(values, reference):
    return float(np.sqrt(np.mean((np.asarray(values) - reference) ** 2)))


def reference_line(ax, t, seq, key):
    waiting = np.asarray(seq["waiting"], dtype=bool)
    values = np.asarray(seq[key])
    ax.plot(
        t,
        np.where(~waiting, values, np.nan),
        **REFERENCE,
        label="Experimental reference (movement)",
    )
    # Filled standing targets are not independent experimental measurements.
    ax.plot(
        t,
        np.where(waiting, values, np.nan),
        color=RED,
        ls=":",
        lw=1.0,
        alpha=0.65,
        label="Filled standing reference",
    )


def shade_movements(axes, seq):
    t = seq["t"]
    for segment in seq["segments"]:
        selected = segment_slice(segment)
        start, end = t[selected.start], t[selected.stop - 1]
        for ax in axes:
            ax.axvspan(start, end, color="#f1d9ad", alpha=0.22, zorder=0)
        axes[0].text(
            (start + end) / 2,
            0.98,
            f"{int(segment['loadKg'])} kg / code {int(segment['speedCode'])}",
            transform=axes[0].get_xaxis_transform(),
            ha="center",
            va="top",
            fontsize=8,
            bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none", "pad": 2},
        )
    for ax in axes:
        ax.set_xlim(t[0], t[-1])
        ax.grid(alpha=0.2)


def sequence_figure(group, cap):
    seq, t = group["sequence"], group["sequence"]["t"]
    fig, axes = plt.subplots(
        4,
        2,
        figsize=(14.2, 11.6),
        sharex="col",
        gridspec_kw={"height_ratios": [1.15, 1.1, 0.75, 0.65]},
    )
    for col, label in enumerate(LABELS):
        nominal, assisted = group["nominal"][col], group["predictions"][col]
        ax = axes[0, col]
        reference_line(ax, t, seq, "torque")
        ax.plot(
            t, nominal["torque"], **NOMINAL, label="Optimized Hill model, no assistance"
        )
        ax.plot(
            t,
            assisted["torque"],
            **ASSISTED,
            label="Same Hill model + MPC: muscle torque",
        )
        ax.plot(
            t,
            assisted["torque"] + assisted["assistanceTorque"],
            color=GREEN,
            ls=":",
            lw=1.35,
            label="MPC: muscle + actuator torque",
        )
        back_threshold_lines(ax, seq)
        ax.set_ylabel("Torque (N·m)")
        ax.set_title(label, pad=12)
        ax.margins(y=0.21)
        ax = axes[1, col]
        reference_line(ax, t, seq, "q")
        ax.plot(t, nominal["q"], **NOMINAL)
        ax.plot(t, assisted["q"], **ASSISTED)
        ax.set_ylabel("Trunk angle (rad)")
        ax.set_title(
            f"MPC angle range: {assisted['q'].min():.3f}–{assisted['q'].max():.3f} rad",
            fontsize=10,
        )
        ax = axes[2, col]
        ax.plot(t, assisted["assistanceTorque"], color=GREEN, lw=1.4)
        ax.axhline(cap, color=RED, ls="--", lw=0.9)
        ax.set_ylabel("Assistance (N·m)")
        ax.set_ylim(-3, cap * 1.17)
        ax.text(
            0.015,
            0.97,
            f"Peak: {assisted['assistanceTorque'].max():.2f} N·m; cap: {cap:g} N·m",
            transform=ax.transAxes,
            fontsize=9,
            va="top",
            bbox={"facecolor": "white", "alpha": 0.9, "edgecolor": "none", "pad": 2},
        )
        ax = axes[3, col]
        delta = np.asarray(assisted["q"]) - nominal["q"]
        ax.axhline(0, color="#888888", lw=0.7)
        ax.plot(t, delta, color=PURPLE, lw=1.2)
        difference_limit = max(0.065, 1.25 * np.abs(delta).max())
        ax.set_ylim(-difference_limit, difference_limit)
        ax.set_ylabel("MPC − no assist\nangle (rad)")
        ax.set_xlabel("Time (s)")
        ax.set_title(
            f"Maximum |angle difference|: {np.abs(delta).max():.4f} rad",
            fontsize=9,
            loc="left",
            pad=5,
        )
        shade_movements(axes[:, col], seq)
    fig.suptitle(
        f"Optimized Hill model and MPC — continuous sequence {int(group['group'])}",
        fontsize=15,
        y=0.993,
    )
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=2,
        frameon=False,
        bbox_to_anchor=(0.5, 0.031),
        fontsize=9,
    )
    fig.text(
        0.5,
        0.009,
        "Net muscle torque = back − abdomen. Black dashed targets apply to BACK torque. Standing references are filled values.",
        ha="center",
        fontsize=8.5,
    )
    fig.tight_layout(rect=(0, 0.123, 1, 0.965), h_pad=1.45)
    return fig


def case_index(groups):
    cases = {}
    for group in groups:
        for segment in group["sequence"]["segments"]:
            key = (int(segment["loadKg"]), int(segment["speedCode"]))
            if key in cases:
                raise ValueError(f"Duplicate condition: {key}")
            cases[key] = (group, segment)
    expected = {(load, speed) for load in (5, 10, 15) for speed in (35, 45, 55)}
    if set(cases) != expected:
        raise ValueError("The full plot suite requires all nine load/speed conditions.")
    return cases


def nine_case_figure(cases, source, key):
    fig, axes = plt.subplots(3, 3, figsize=(13.8, 9.4), sharey="row")
    for row, load in enumerate((5, 10, 15)):
        for col, speed in enumerate((35, 45, 55)):
            group, segment = cases[(load, speed)]
            seq, selected = group["sequence"], segment_slice(segment)
            t = seq["t"][selected] - seq["t"][selected.start]
            nominal, assisted = group["nominal"][source], group["predictions"][source]
            ax = axes[row, col]
            ax.plot(t, seq[key][selected], **REFERENCE, label="Experimental reference")
            ax.plot(
                t,
                nominal[key][selected],
                **NOMINAL,
                label="Optimized Hill, no assistance",
            )
            ax.plot(
                t,
                assisted[key][selected],
                **ASSISTED,
                label=(
                    "Same Hill + MPC: muscle torque"
                    if key == "torque"
                    else "Same Hill + MPC"
                ),
            )
            if key == "torque":
                ax.plot(
                    t,
                    assisted["torque"][selected]
                    + assisted["assistanceTorque"][selected],
                    color=GREEN,
                    ls=":",
                    lw=1.35,
                    label="MPC: muscle + actuator torque",
                )
                threshold = BACK_THRESHOLDS_NM[load]
                ax.axhline(
                    threshold, **THRESHOLD_STYLE, label="MPC back-torque soft threshold"
                )
                ax.text(
                    0.98,
                    threshold,
                    f"{threshold:g} N·m",
                    transform=ax.get_yaxis_transform(),
                    ha="right",
                    va="bottom",
                    fontsize=8,
                    bbox={
                        "facecolor": "white",
                        "alpha": 0.85,
                        "edgecolor": "none",
                        "pad": 1,
                    },
                )
                text = (
                    f"Back peak reduction: {100*(1-assisted['backTorque'][selected].max()/nominal['backTorque'][selected].max()):.1f}%\n"
                    f"Assistance peak: {assisted['assistanceTorque'][selected].max():.1f} N·m"
                )
            else:
                text = (
                    f"RMSE no assist / MPC: {rmse(nominal['q'][selected], seq['q'][selected])*180/np.pi:.1f}° / "
                    f"{rmse(assisted['q'][selected], seq['q'][selected])*180/np.pi:.1f}°"
                )
            ax.text(
                0.025,
                0.025,
                text,
                fontsize=8,
                transform=ax.transAxes,
                va="bottom",
                bbox={
                    "facecolor": "white",
                    "alpha": 0.88,
                    "edgecolor": "none",
                    "pad": 2,
                },
            )
            ax.set_title(f"{load} kg · speed code {speed}", fontsize=11)
            ax.set_xlim(t[0], t[-1])
            if key == "q":
                ax.set_ylim(0.22, 1.8)
            else:
                ax.margins(y=0.22)
            ax.grid(alpha=0.2)
            if col == 0:
                ax.set_ylabel("Trunk angle (rad)" if key == "q" else "Torque (N·m)")
            if row == 2:
                ax.set_xlabel("Time from movement onset (s)")
    title = "Angle" if key == "q" else "Muscle and assisted total torque"
    fig.suptitle(
        f"{title}: all nine conditions — {LABELS[source]}", fontsize=15, y=0.995
    )
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=2 if key == "torque" else 3,
        frameon=False,
        bbox_to_anchor=(0.5, 0.021),
        fontsize=9,
    )
    footer = (
        "Black dashed lines are BACK-muscle soft targets, not hard limits on net or muscle + actuator torque."
        if key == "torque"
        else "One fixed parameter vector per activation source. Original reference scale; no time alignment or smoothing."
    )
    fig.text(0.5, 0.008, footer, ha="center", fontsize=8.5)
    fig.tight_layout(rect=(0, 0.135 if key == "torque" else 0.078, 1, 0.965), h_pad=2)
    return fig


def muscle_figure(group):
    seq, t = group["sequence"], group["sequence"]["t"]
    movement = ~np.asarray(seq["waiting"], dtype=bool)
    fig, axes = plt.subplots(
        3,
        2,
        figsize=(14.2, 9.4),
        sharex="col",
        gridspec_kw={"height_ratios": [1.0, 1.0, 0.65]},
    )
    for col, label in enumerate(LABELS):
        nominal, assisted = group["nominal"][col], group["predictions"][col]
        for row, key, ylabel in (
            (0, "backTorque", "Back torque (N·m)"),
            (1, "frontTorque", "Abdominal torque (N·m)"),
            (2, "boundaryTorque", "Soft-boundary torque (N·m)"),
        ):
            ax = axes[row, col]
            ax.plot(t, nominal[key], **NOMINAL, label="Optimized Hill, no assistance")
            ax.plot(t, assisted[key], **ASSISTED, label="Same Hill + MPC")
            ax.set_ylabel(ylabel)
            if row < 2:
                values = np.r_[nominal[key], assisted[key]]
                ax.set_ylim(
                    min(0.0, float(values.min())),
                    values.max() * (1.25 if row == 0 else 1.15),
                )
            else:
                limit = max(
                    1.0, 1.2 * np.max(np.abs(np.r_[nominal[key], assisted[key]]))
                )
                ax.set_ylim(-limit, limit)
                ax.set_xlabel("Time (s)")
                ax.text(
                    0.015,
                    0.86,
                    f"Peak |boundary|: {np.abs(nominal[key]).max():.3g} / "
                    f"{np.abs(assisted[key]).max():.3g} N·m (no assist / MPC)",
                    transform=ax.transAxes,
                    fontsize=9,
                )
        axes[0, col].set_title(label, pad=12)
        back_threshold_lines(axes[0, col], seq)
        ratio_n = (
            nominal["frontTorque"][movement].mean()
            / nominal["backTorque"][movement].mean()
        )
        ratio_p = (
            assisted["frontTorque"][movement].mean()
            / assisted["backTorque"][movement].mean()
        )
        axes[1, col].set_title(
            f"Action mean abdomen/back: {100*ratio_n:.1f}% / {100*ratio_p:.1f}% (no assist / MPC)",
            fontsize=10,
        )
        shade_movements(axes[:, col], seq)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.suptitle(
        f"Back, abdominal and soft-boundary torques — sequence {int(group['group'])}",
        fontsize=15,
        y=0.995,
    )
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 0.025),
    )
    fig.text(
        0.5,
        0.008,
        "Muscle torques are positive magnitudes; abdominal torque is subtracted in net joint torque. These are model predictions.",
        ha="center",
        fontsize=8.5,
    )
    fig.tight_layout(rect=(0, 0.08, 1, 0.96), h_pad=1.9)
    return fig


def heavy_load_figure(cases):
    fig, axes = plt.subplots(3, 2, figsize=(13.5, 9.8))
    for row, speed in enumerate((35, 45, 55)):
        group, segment = cases[(15, speed)]
        seq, selected = group["sequence"], segment_slice(segment)
        t = seq["t"][selected] - seq["t"][selected.start]
        for col, label in enumerate(LABELS):
            ax = axes[row, col]
            ax.plot(t, seq["q"][selected], **REFERENCE, label="Experimental reference")
            ax.plot(
                t,
                group["nominal"][col]["q"][selected],
                **NOMINAL,
                label="Optimized Hill, no assistance",
            )
            ax.plot(
                t,
                group["predictions"][col]["q"][selected],
                **ASSISTED,
                label="Same Hill + MPC",
            )
            ax.set_title(f"15 kg · speed code {speed} — {label}", fontsize=11)
            ax.set_ylabel("Trunk angle (rad)")
            ax.set_xlim(t[0], t[-1])
            ax.set_ylim(0.22, 1.8)
            ax.grid(alpha=0.2)
            if row == 2:
                ax.set_xlabel("Time from movement onset (s)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.suptitle(
        "15 kg angle detail — same optimized model before and after MPC",
        fontsize=15,
        y=0.995,
    )
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False)
    fig.tight_layout(rect=(0, 0.045, 1, 0.965), h_pad=2)
    return fig


def comparison_rows(cases):
    rows = []
    for (load, speed), (group, segment) in sorted(cases.items()):
        selected = segment_slice(segment)
        seq = group["sequence"]
        for source, label in enumerate(LABELS):
            nominal, assisted = group["nominal"][source], group["predictions"][source]
            rows.append(
                {
                    "source": label,
                    "group": int(group["group"]),
                    "loadKg": load,
                    "speedCode": speed,
                    "noAssistThetaRMSEDegrees": rmse(
                        nominal["q"][selected], seq["q"][selected]
                    )
                    * 180
                    / np.pi,
                    "mpcThetaRMSEDegrees": rmse(
                        assisted["q"][selected], seq["q"][selected]
                    )
                    * 180
                    / np.pi,
                    "maxMpcAngleChangeRad": float(
                        np.abs(assisted["q"][selected] - nominal["q"][selected]).max()
                    ),
                    "noAssistNetTorqueRMSEnm": rmse(
                        nominal["torque"][selected], seq["torque"][selected]
                    ),
                    "mpcNetTorqueRMSEnm": rmse(
                        assisted["torque"][selected], seq["torque"][selected]
                    ),
                    "mpcMusclePlusActuatorRMSEnm": rmse(
                        assisted["torque"][selected]
                        + assisted["assistanceTorque"][selected],
                        seq["torque"][selected],
                    ),
                    "noAssistBackPeakNm": float(nominal["backTorque"][selected].max()),
                    "mpcBackPeakNm": float(assisted["backTorque"][selected].max()),
                    "backPeakReductionPercent": float(
                        100
                        * (
                            1
                            - assisted["backTorque"][selected].max()
                            / nominal["backTorque"][selected].max()
                        )
                    ),
                    "noAssistAbdomenMeanNm": float(
                        nominal["frontTorque"][selected].mean()
                    ),
                    "mpcAbdomenMeanNm": float(assisted["frontTorque"][selected].mean()),
                    "assistancePeakNm": float(
                        assisted["assistanceTorque"][selected].max()
                    ),
                }
            )
    return rows


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def publish(names, metadata):
    target = ROOT / "docs/hill_mpc_evaluation"
    target.mkdir(parents=True, exist_ok=True)
    for name in names:
        shutil.copy2(OUT / name, target / name)
    # The attached historical image path is an explicitly documented current-result alias.
    legacy = ROOT / "docs/legacy_mpc_evaluation"
    if (legacy / "export_metadata.json").is_file():
        for suffix in ("png", "pdf"):
            shutil.copy2(
                OUT / f"mpc_torque_theta_group1.{suffix}",
                legacy / f"legacy_mpc_torque_theta.{suffix}",
            )
        manifest_path = legacy / "export_metadata.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        for entry in manifest["files"]:
            entry["sha256"] = sha256(legacy / entry["file"])
        manifest_path.write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
        )
    manifest = {
        "modelSHA256": metadata["modelSHA256"],
        "sourceResultSHA256": sha256(OUT / "hill_mpc.mat"),
        "activeFigureBundle": BUNDLE_NAME,
        "historicalFigureBundles": [
            name
            for name in (
                "optimized_model_and_mpc_all_results.pdf",
                "optimized_model_and_mpc_with_thresholds.pdf",
                "optimized_model_and_mpc_tuned.pdf",
                "optimized_model_and_mpc_headroom.pdf",
            )
            if (target / name).exists()
        ],
        "policy": "Saved optimized Hill / MPC trajectories only; one fixed parameter vector per activation source; original references and samples unchanged",
        "files": [
            {"file": path.name, "sha256": sha256(path)}
            for path in sorted(target.iterdir())
            if path.is_file() and path.name != "export_metadata.json"
        ],
    }
    (target / "export_metadata.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


def render_results(publish_docs=False):
    global CURRENT_OPTIONS
    result = loadmat(OUT / "hill_mpc.mat", simplify_cells=True)["result"]
    metadata = json.loads((OUT / "evaluation.json").read_text(encoding="utf-8"))
    CURRENT_OPTIONS = metadata["options"]
    groups = sorted(result["groups"], key=lambda group: int(group["group"]))
    cases = case_index(groups)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 12,
            "pdf.fonttype": 42,
        }
    )
    names = []
    bundle = BUNDLE_NAME
    with PdfPages(
        OUT / bundle,
        metadata={
            "Title": "Optimized Hill model and MPC: all nine conditions",
            "Subject": "Both activation sources; saved trajectories, no rescaling or fitting",
        },
    ) as pdf:

        def save(fig, name):
            for suffix in ("png", "pdf"):
                filename = f"{name}.{suffix}"
                fig.savefig(
                    OUT / filename, dpi=180, bbox_inches="tight", facecolor="white"
                )
                names.append(filename)
            pdf.savefig(fig, bbox_inches="tight", facecolor="white")
            plt.close(fig)

        for group in groups:
            save(
                sequence_figure(group, metadata["options"]["capNm"]),
                f"mpc_torque_theta_group{int(group['group'])}",
            )
        for source, slug in enumerate(SLUGS):
            for key, prefix in (("q", "theta"), ("torque", "torque")):
                save(nine_case_figure(cases, source, key), f"{prefix}_all9_{slug}")
        for group in groups:
            save(
                muscle_figure(group),
                f"mpc_muscle_and_boundary_torque_group{int(group['group'])}",
            )
        save(heavy_load_figure(cases), "theta_15kg_all_speeds")
    names.append(bundle)
    for suffix in ("png", "pdf"):
        alias = f"mpc_muscle_and_boundary_torque.{suffix}"
        shutil.copy2(
            OUT / f"mpc_muscle_and_boundary_torque_group1.{suffix}", OUT / alias
        )
        names.append(alias)
    rows = comparison_rows(cases)
    csv_name = "model_mpc_comparison_metrics.csv"
    with (OUT / csv_name).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    names.append(csv_name)
    if publish_docs:
        publish(names, metadata)
    print(
        f"Saved 11 figures (PNG + vector PDF), an 11-page PDF and 18 comparison rows to {OUT}"
    )
    for label in LABELS:
        selected = [row for row in rows if row["source"] == label]
        print(
            f"{label}: mean case theta RMSE, no assist / MPC = "
            f"{np.mean([row['noAssistThetaRMSEDegrees'] for row in selected]):.3f} / "
            f"{np.mean([row['mpcThetaRMSEDegrees'] for row in selected]):.3f} degrees; "
            f"max MPC angle change = {max(row['maxMpcAngleChangeRad'] for row in selected):.5f} rad"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--publish",
        action="store_true",
        help="Also refresh the tracked documentation snapshots and hashes.",
    )
    render_results(publish_docs=parser.parse_args().publish)
