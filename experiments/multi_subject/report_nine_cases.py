"""Publish nine-condition BSE and 75%-peak NMPC results without hiding failures."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat
from scipy.signal import find_peaks

from audit_subjects import ROOT
from fit_subject import many, SHARED
from fit_hill_compliant import PARAMETER_NAMES
from control_nine_cases import BASE, load_subject, LABELS, digest

DOCS = ROOT / "docs/nine_case_75pct"
SUBJECTS = ("240528Pilot", "240725Exp_A01", "240730Exp_A01")
COLORS = ("#1469b0", "#8450ad")


def save(fig, name, files):
    path = DOCS / (name + ".png")
    fig.savefig(path, dpi=155)
    plt.close(fig)
    files.append(path)


def movements(groups):
    rows = []
    for group in groups:
        for seg in many(group["sequence"]["segments"]):
            rows.append(
                (
                    group,
                    seg,
                    np.arange(int(seg["startIndex"]) - 1, int(seg["endIndex"])),
                )
            )
    return rows


def bse_scores(subject, groups):
    rows = []
    for group, seg, ids in movements(groups):
        seq = group["sequence"]
        ref = np.asarray(seq["q"])[ids]
        for source, pred in enumerate(many(group["predictions"])):
            q = np.asarray(pred["q"])[ids]
            torque = np.asarray(seq["torque"])[ids]
            valid = np.isfinite(torque)
            rows.append(
                dict(
                    subject=subject,
                    source=LABELS[source],
                    loadKg=float(seg["loadKg"]),
                    speedCode=int(seg["speedCode"]),
                    recordId=seg.get(
                        "recordId", f"{seg['loadKg']:g}kg_{seg['speedCode']:g}"
                    ),
                    shapeCorrelation=float(np.corrcoef(q, ref)[0, 1]),
                    thetaRMSEDegrees=float(
                        np.sqrt(np.mean((q - ref) ** 2)) * 180 / np.pi
                    ),
                    torqueRMSENm=(
                        float(
                            np.sqrt(
                                np.mean(
                                    (
                                        np.asarray(pred["torque"])[ids][valid]
                                        - torque[valid]
                                    )
                                    ** 2
                                )
                            )
                        )
                        if valid.any()
                        else None
                    ),
                    predictedValleys=len(
                        find_peaks(-q, prominence=0.12, distance=50)[0]
                    ),
                    referenceValleys=len(
                        find_peaks(-ref, prominence=0.12, distance=50)[0]
                    ),
                    minimumAngleRad=float(q.min()),
                    maximumAngleRad=float(q.max()),
                    backPeakNm=float(np.asarray(pred["backTorque"])[ids].max()),
                    abdomenPeakNm=float(np.asarray(pred["frontTorque"])[ids].max()),
                    abdomenBackMeanRatio=float(
                        np.asarray(pred["frontTorque"])[ids].sum()
                        / np.asarray(pred["backTorque"])[ids].sum()
                    ),
                    boundaryPeakNm=float(
                        np.abs(np.asarray(pred["boundaryTorque"])[ids]).max()
                    ),
                    torqueCoverageFraction=float(valid.mean()),
                )
            )
    return rows


def bse_figures(subject, groups, files):
    cases = movements(groups)
    assert len(cases) == 9
    for field in ("q", "torque"):
        fig, axes = plt.subplots(3, 3, figsize=(12.5, 8.7), layout="constrained")
        for ax, (group, seg, ids) in zip(axes.flat, cases):
            seq = group["sequence"]
            t = np.asarray(seq["t"])[ids]
            t = t - t[0]
            ref = np.asarray(seq[field])[ids]
            ax.plot(
                t,
                ref,
                "--",
                color="#c13b31",
                lw=1.2,
                label=(
                    "Original Scale IK" if field == "q" else "Available CEINMS estimate"
                ),
            )
            for source, pred in enumerate(many(group["predictions"])):
                ax.plot(
                    t,
                    np.asarray(pred[field])[ids],
                    color=COLORS[source],
                    lw=1.3,
                    label=LABELS[source],
                )
            ax.set_title(
                f"{seg['loadKg']:g} kg / code {seg['speedCode']:g}", fontsize=10
            )
            ax.set_xlabel("Action time (s)", fontsize=9)
            ax.set_ylabel(
                "Theta (rad)" if field == "q" else "Net muscle torque (N m)", fontsize=9
            )
            ax.grid(alpha=0.2)
            ax.tick_params(labelsize=8)
        handles, labels = axes.flat[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="outside lower center", ncol=3, fontsize=8)
        fig.suptitle(
            subject
            + " | BSE | nine recorded conditions\nOne fixed model per activation source; missing reference torque remains missing",
            fontsize=12,
        )
        save(fig, subject + "_bse_" + ("theta" if field == "q" else "torque"), files)
    fig, axes = plt.subplots(
        2, 2, figsize=(12.5, 8.7), layout="constrained", sharey=True
    )
    for source in (0, 1):
        labels = []
        reference = [[], []]
        predicted = [[], []]
        for group, seg, ids in cases:
            ref = np.asarray(group["sequence"]["q"])[ids]
            q = np.asarray(many(group["predictions"])[source]["q"])[ids]
            valley = find_peaks(-ref, prominence=0.12, distance=50)[0]
            if len(valley) >= 2:
                middle = int(valley[0] + np.argmax(ref[valley[0] : valley[-1] + 1]))
            else:
                middle = len(ref) // 2
            for j, part in enumerate((slice(None, middle + 1), slice(middle, None))):
                reference[j].append(float(ref[part].min()))
                predicted[j].append(float(q[part].min()))
            labels.append(f"{seg['loadKg']:g} kg\n{seg['speedCode']:g}")
        x = np.arange(9)
        for valley in (0, 1):
            ax = axes[valley, source]
            ax.bar(
                x - 0.18,
                reference[valley],
                0.36,
                color="#d95454",
                label="Original reference",
            )
            ax.bar(
                x + 0.18,
                predicted[valley],
                0.36,
                color=COLORS[source],
                label="Fixed-parameter BSE",
            )
            ax.set_xticks(x, labels, fontsize=8)
            ax.set_ylabel("Valley angle (rad)")
            ax.grid(axis="y", alpha=0.2)
            ax.set_title(LABELS[source] + f" | valley {valley+1}", fontsize=10)
    fig.suptitle(
        subject
        + " | both flexion depths in all nine recorded conditions\nSmaller angle = deeper flexion; reference unchanged",
        fontsize=12,
    )
    fig.legend(
        *axes[0, 0].get_legend_handles_labels(),
        loc="outside lower center",
        ncol=2,
        fontsize=9,
    )
    save(fig, subject + "_valley_depths", files)


def mpc_figures(subject, groups, files):
    cases = movements(groups)
    for source in (0, 1):
        for kind in ("torque", "theta", "abdomen_boundary"):
            fig, axes = plt.subplots(3, 3, figsize=(12.5, 8.7), layout="constrained")
            for ax, (group, seg, ids) in zip(axes.flat, cases):
                seq = group["sequence"]
                t = np.asarray(seq["t"])[ids]
                t = t - t[0]
                p = many(group["predictions"])[source]
                nom = many(group["nominal"])[source]
                limit = float(np.asarray(p["backTorqueThresholdNm"])[ids][0])
                if kind == "torque":
                    ax.plot(
                        t,
                        np.asarray(nom["backTorque"])[ids],
                        ":",
                        color="#8ba9c5",
                        lw=1.1,
                        label="Unassisted back",
                    )
                    ax.plot(
                        t,
                        np.asarray(p["backTorque"])[ids],
                        color=COLORS[source],
                        lw=1.4,
                        label="NMPC back",
                    )
                    ax.plot(
                        t,
                        np.asarray(nom["frontTorque"])[ids],
                        ":",
                        color="#bb9ed2",
                        lw=0.9,
                        label="Unassisted abdomen",
                    )
                    ax.plot(
                        t,
                        np.asarray(p["frontTorque"])[ids],
                        color="#8450ad",
                        lw=1.0,
                        label="NMPC abdomen",
                    )
                    ax.axhline(
                        limit,
                        color="#333",
                        ls="--",
                        lw=1.0,
                        label="75% back-peak threshold",
                    )
                    ax.axhspan(max(0, limit - 5), limit, color="#208552", alpha=0.10)
                    margin = limit - np.max(np.asarray(p["backTorque"])[ids])
                    ax.text(
                        0.02,
                        0.97,
                        f"Limit {limit:.1f}; peak margin {margin:.1f} N m",
                        transform=ax.transAxes,
                        va="top",
                        fontsize=7,
                        color="#c13b31" if margin < -0.05 else "#333",
                    )
                    ax.set_ylabel("Muscle torque (N m)", fontsize=9)
                elif kind == "theta":
                    ax.plot(
                        t,
                        np.asarray(seq["q"])[ids],
                        "--",
                        color="#c13b31",
                        lw=1.1,
                        label="Original Scale reference",
                    )
                    ax.plot(
                        t,
                        np.asarray(nom["q"])[ids],
                        color=".65",
                        lw=1.0,
                        label="Unassisted BSE",
                    )
                    ax.plot(
                        t,
                        np.asarray(p["q"])[ids],
                        color=COLORS[source],
                        lw=1.3,
                        label="NMPC assisted",
                    )
                    ax.set_ylabel("Theta (rad)", fontsize=9)
                else:
                    ax.plot(
                        t,
                        np.asarray(nom["frontTorque"])[ids],
                        ":",
                        color="#b395c9",
                        lw=1.2,
                        label="Unassisted abdomen",
                    )
                    ax.plot(
                        t,
                        np.asarray(p["frontTorque"])[ids],
                        color="#8450ad",
                        lw=1.3,
                        label="NMPC abdomen",
                    )
                    ax.plot(
                        t,
                        np.asarray(nom["boundaryTorque"])[ids],
                        ":",
                        color="#ddb485",
                        lw=1.0,
                        label="Unassisted soft boundary",
                    )
                    ax.plot(
                        t,
                        np.asarray(p["boundaryTorque"])[ids],
                        color="#b56c22",
                        lw=1.1,
                        label="NMPC soft boundary",
                    )
                    ax.axhline(0, color=".7", lw=0.5)
                    ax.set_ylabel("Torque (N m)", fontsize=9)
                ax.set_title(
                    f"{seg['loadKg']:g} kg / code {seg['speedCode']:g}", fontsize=10
                )
                ax.set_xlabel("Action time (s)", fontsize=9)
                ax.grid(alpha=0.2)
                ax.tick_params(labelsize=8)
            fig.suptitle(
                subject
                + " | "
                + LABELS[source]
                + " | NMPC "
                + kind
                + "\nSame subject-wide controller; 80 N m cap and 200 N m/s slew limit",
                fontsize=12,
            )
            fig.legend(
                *axes.flat[0].get_legend_handles_labels(),
                loc="outside lower center",
                ncol=3,
                fontsize=8,
            )
            save(fig, subject + "_nmpc_" + str(source + 1) + "_" + kind, files)


def continuous_figure(subject, group, files):
    seq = group["sequence"]
    t = seq["t"]
    fig, axes = plt.subplots(
        5, 2, figsize=(12.5, 10), layout="constrained", sharex=True
    )
    for source in (0, 1):
        p = many(group["predictions"])[source]
        nom = many(group["nominal"])[source]
        axes[0, source].plot(
            t, nom["backTorque"], ":", color="#8ba9c5", lw=1, label="Unassisted back"
        )
        axes[0, source].plot(
            t, p["backTorque"], color=COLORS[source], lw=1.2, label="NMPC back"
        )
        for seg in many(seq["segments"]):
            ids = np.arange(int(seg["startIndex"]) - 1, int(seg["endIndex"]))
            lim = float(np.asarray(p["backTorqueThresholdNm"])[ids][0])
            axes[0, source].plot(
                np.asarray(t)[ids], np.full(len(ids), lim), "--", color="#333", lw=1
            )
            axes[0, source].fill_between(
                np.asarray(t)[ids], lim - 5, lim, color="#208552", alpha=0.10
            )
            for ax in axes[:, source]:
                ax.axvspan(
                    float(t[ids[0]]),
                    float(t[ids[-1]]),
                    color="#f5e8d6",
                    alpha=0.25,
                    zorder=0,
                )
            axes[0, source].text(
                float(np.mean(np.asarray(t)[ids])),
                1.02,
                f"{seg['loadKg']:g} kg / {seg['speedCode']:g}",
                ha="center",
                fontsize=6,
                transform=axes[0, source].get_xaxis_transform(),
            )
        axes[1, source].plot(t, seq["q"], "--", color="#c13b31", lw=1)
        axes[1, source].plot(t, nom["q"], color=".65", lw=1)
        axes[1, source].plot(t, p["q"], color=COLORS[source], lw=1.2)
        axes[2, source].plot(t, p["assistanceTorque"], color="#25844d", lw=1)
        axes[2, source].axhline(80, color="#c13b31", ls="--", lw=0.9)
        axes[3, source].plot(t, nom["frontTorque"], ":", color="#bb9ed2", lw=1)
        axes[3, source].plot(t, p["frontTorque"], color="#8450ad", lw=1.2)
        axes[4, source].plot(t, p["boundaryTorque"], color="#b56c22", lw=1)
        axes[0, source].set_title(LABELS[source], fontsize=10)
        axes[0, source].legend(loc="upper left", fontsize=6)
        for ax, label in zip(
            axes[:, source],
            [
                "Back torque (N m)",
                "Theta (rad)",
                "Assistance (N m)",
                "Abdominal torque (N m)",
                "Soft boundary (N m)",
            ],
        ):
            ax.set_ylabel(label, fontsize=8)
            ax.grid(alpha=0.2)
            ax.tick_params(labelsize=8)
        axes[-1, source].set_xlabel("Continuous time (s)")
    fig.suptitle(
        subject
        + f" | continuous group {int(group['group'])}\nBlack dashed: 75% unassisted back-peak threshold. Green: 0-5 N m below. No state resets between actions.",
        fontsize=11,
    )
    save(fig, subject + f"_continuous_{int(group['group'])}", files)


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    from reportlab.pdfgen import canvas
    from reportlab.lib.utils import ImageReader

    DOCS.mkdir(parents=True, exist_ok=True)
    selected = json.loads((BASE / "selected_nmpc.json").read_text())
    files = []
    all_bse = []
    all_mpc = []
    audits = []
    summaries = []
    reports = {}
    parameter_values = []
    controller_values = []
    for subject in SUBJECTS:
        cache, models, groups, model_path, replay_path = load_subject(subject)
        physical = np.vstack([m["parameters"] for m in models])
        assert np.array_equal(physical[0, SHARED], physical[1, SHARED])
        conditions = movements(groups)
        assert (
            len(conditions) == 9
            and len({(s["loadKg"], s["speedCode"]) for _, s, _ in conditions}) == 9
        )
        scores = bse_scores(subject, groups)
        all_bse += scores
        for source, m in enumerate(models):
            names = list(PARAMETER_NAMES) + [
                "intentAmplitude",
                "intentExponent",
                "intentTau",
                "postureKp",
                "postureKd",
                "maximumNegativeTorqueRequest",
                "neuralFeedbackBlend",
                "coactivationNm",
                "causalLeadTime",
                "loadTonicCompensation",
            ]
            values = np.r_[m["parameters"], m["neuralParameters"]]
            parameter_values.append(
                dict(
                    subject=subject,
                    source=LABELS[source],
                    intentDecoder=m.get("intentDecoder", "subtractive"),
                    **dict(zip(names, map(float, values))),
                )
            )
            h = hashlib.sha256(
                np.r_[m["parameters"], m["neuralParameters"]].astype("<f8").tobytes()
            ).hexdigest()
            for _, seg, _ in conditions:
                audits.append(
                    dict(
                        subject=subject,
                        source=LABELS[source],
                        loadKg=seg["loadKg"],
                        speedCode=seg["speedCode"],
                        fixedParameterSHA256=h,
                        modelSHA256=digest(model_path),
                    )
                )
        folder = BASE / subject / "nmpc" / selected[subject]
        report = json.loads((folder / "evaluation.json").read_text())
        assert report["signature"]["modelSHA256"] == digest(
            model_path
        ), "MPC belongs to a different BSE model"
        assert report["fullNineConditions"] and report["caseCountPerSource"] == 9
        reports[subject] = report
        controller_keys = [
            "controlPeriod",
            "predictionDt",
            "points",
            "blocks",
            "angleWeight",
            "velocityWeight",
            "backTorqueUnderWeight",
            "trackingHeadroomNm",
            "backTorqueMarginNm",
            "predictionErrorGuardScale",
            "predictionErrorGuardFloorNm",
            "predictionErrorGuardCapNm",
            "torqueConstraintHorizonSeconds",
            "assistanceEffortWeight",
            "capNm",
            "rateNmPerSec",
            "inputForecast",
        ]
        controller_values.append(
            dict(
                subject=subject,
                candidate=selected[subject],
                **{key: report["options"][key] for key in controller_keys},
            )
        )
        all_mpc += [dict(subject=subject, **r) for r in report["metrics"]]
        actual = loadmat(folder / "hill_mpc.mat", simplify_cells=True)["result"]
        bse_figures(subject, groups, files)
        mpc_figures(subject, many(actual["groups"]), files)
        for group in many(actual["groups"]):
            if len(many(group["sequence"]["segments"])) > 1:
                continuous_figure(subject, group, files)
        extra_path = folder / "continuous_mpc.mat"
        if extra_path.exists():
            extra = loadmat(extra_path, simplify_cells=True)["result"]
            for group in many(extra["groups"]):
                continuous_figure(subject, group, files)
        for source, label in enumerate(LABELS):
            b = [r for r in scores if r["source"] == label]
            m = [r for r in report["metrics"] if r["source"] == label]
            summaries.append(
                dict(
                    subject=subject,
                    source=label,
                    conditions=9,
                    bseMeanShapeCorrelation=float(
                        np.mean([r["shapeCorrelation"] for r in b])
                    ),
                    bseMinimumShapeCorrelation=min(r["shapeCorrelation"] for r in b),
                    bseTwoValleys=sum(r["predictedValleys"] == 2 for r in b),
                    bseMeanThetaRMSEDegrees=float(
                        np.mean([r["thetaRMSEDegrees"] for r in b])
                    ),
                    bseMeanTorqueRMSENm=float(
                        np.mean(
                            [
                                r["torqueRMSENm"]
                                for r in b
                                if r["torqueRMSENm"] is not None
                            ]
                        )
                    ),
                    mpcWithinFiveNmBelow=sum(r["peakWithinFiveNmBelow"] for r in m),
                    mpcMaximumExcessNm=max(r["excessPeakNm"] for r in m),
                    mpcMaximumPeakShortfallNm=max(r["peakMarginNm"] for r in m),
                    mpcChangedValleyCases=sum(
                        r["predictedValleys"] != r["nominalValleys"] for r in m
                    ),
                    mpcAssistancePeakNm=max(r["assistancePeakNm"] for r in m),
                    mpcMeanDemandShortfallNm=float(
                        np.mean([r["demandMeanShortfallNm"] for r in m])
                    ),
                    mpcMaximumBoundaryNm=max(r["boundaryPeakNm"] for r in m),
                    bseMaximumAbdomenBackMeanRatio=max(
                        r["abdomenBackMeanRatio"] for r in b
                    ),
                    modelSHA256=digest(model_path),
                    selectedController=selected[subject],
                )
            )
    fig, ax = plt.subplots(figsize=(12.5, 6.3), layout="constrained")
    ax.axis("off")
    cells = []
    for row in summaries:
        cells.append(
            [
                row["subject"],
                row["source"].replace(" activation", ""),
                f"{row['bseMeanShapeCorrelation']:.3f}",
                f"{row['bseMinimumShapeCorrelation']:.3f}",
                f"{row['bseTwoValleys']}/9",
                f"{row['mpcWithinFiveNmBelow']}/9",
                f"{row['mpcMaximumExcessNm']:.2f}",
                f"{row['mpcMaximumPeakShortfallNm']:.2f}",
                str(row["mpcChangedValleyCases"]),
            ]
        )
    table = ax.table(
        cellText=cells,
        colLabels=[
            "Subject",
            "Source",
            "Mean r",
            "Min r",
            "BSE two\nvalleys",
            "NMPC peak\nin band",
            "Max excess\n(N m)",
            "Max gap\n(N m)",
            "NMPC valley\nchanges",
        ],
        cellLoc="center",
        bbox=[0.02, 0.28, 0.96, 0.50],
        colWidths=[0.18, 0.15, 0.075, 0.075, 0.085, 0.10, 0.085, 0.085, 0.10],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2.5)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#d7dde2")
        if row == 0:
            cell.set_facecolor("#e5edf5")
            cell.set_text_props(weight="bold")
        elif row % 2 == 0:
            cell.set_facecolor("#f6f8fa")
    ax.text(
        0,
        0.94,
        "All nine recorded conditions per subject | one BSE vector per source",
        fontsize=15,
        weight="bold",
        transform=ax.transAxes,
    )
    ax.text(
        0,
        0.87,
        "BSE r: angle trend correlation. The two-valley count also flags small additional bends.",
        fontsize=10,
        transform=ax.transAxes,
    )
    ax.text(
        0,
        0.18,
        "NMPC band: assisted back-muscle peak 0-5 N m below 75% of its unassisted peak.\n"
        "Low-demand samples are not forced to that threshold. Every condition and issue is retained.\n"
        "These are calibration/replay simulations; CEINMS is a model estimate, not an independent torque measurement.",
        fontsize=10,
        linespacing=1.6,
        va="top",
        transform=ax.transAxes,
    )
    detail = DOCS / "calibration_overview.png"
    fig.savefig(detail, dpi=170)
    plt.close(fig)
    files.insert(0, detail)
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.8), layout="constrained")
    for ax, subject in zip(axes, SUBJECTS):
        rows = [r for r in all_mpc if r["subject"] == subject]
        lo = min(r["backThresholdNm"] for r in rows) - 10
        hi = max(r["backThresholdNm"] for r in rows) + 10
        x = np.linspace(lo, hi, 100)
        ax.plot(x, x, "--", color="#333", label="Threshold")
        ax.fill_between(x, x - 5, x, color="#208552", alpha=0.12, label="0-5 N m below")
        for source, label in enumerate(LABELS):
            part = [r for r in rows if r["source"] == label]
            ax.scatter(
                [r["backThresholdNm"] for r in part],
                [r["assistedBackPeakNm"] for r in part],
                s=25,
                color=COLORS[source],
                marker="o" if source == 0 else "x",
                label=label,
            )
        count = sum(r["peakWithinFiveNmBelow"] for r in rows)
        ax.set_title(subject + f" | {count}/18 within band", fontsize=10)
        ax.set_xlabel("75% peak threshold (N m)", fontsize=9)
        ax.set_ylabel("Assisted back peak (N m)", fontsize=9)
        ax.grid(alpha=0.2)
        ax.tick_params(labelsize=8)
        ax.set_aspect("equal", adjustable="datalim")
    fig.suptitle(
        "NMPC | all 27 conditions and both activation sources\nTargets are 75% of the corresponding unassisted back-muscle torque peak",
        fontsize=12,
    )
    fig.legend(
        *axes[0].get_legend_handles_labels(),
        loc="outside lower center",
        ncol=4,
        fontsize=8,
    )
    overview = DOCS / "peak_tracking_overview.png"
    fig.savefig(overview, dpi=170)
    plt.close(fig)
    files.insert(0, overview)
    write_csv(DOCS / "bse_cases.csv", all_bse)
    write_csv(DOCS / "nmpc_cases.csv", all_mpc)
    write_csv(DOCS / "fixed_parameter_audit.csv", audits)
    write_csv(DOCS / "bse_parameter_values.csv", parameter_values)
    write_csv(DOCS / "subject_controller_parameters.csv", controller_values)
    write_csv(DOCS / "summary_metrics.csv", summaries)
    summary = dict(
        subjects=summaries,
        totalConditionSourceEvaluations=len(all_mpc),
        allPeakTrackingPassed=all(r["peakWithinFiveNmBelow"] for r in all_mpc),
        allBSEShapeIndicatorsPassed=all(
            r["shapeCorrelation"] >= 0.7 and r["predictedValleys"] == 2 for r in all_bse
        ),
        sameFixedParametersAcrossCases=True,
        independentValidation=False,
        bseIssues=[
            r
            for r in all_bse
            if r["shapeCorrelation"] < 0.7 or r["predictedValleys"] != 2
        ],
        mpcIssues=[
            r
            for r in all_mpc
            if not r["peakWithinFiveNmBelow"]
            or r["predictedValleys"] != r["nominalValleys"]
        ],
    )
    (DOCS / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    target = ROOT / "output/pdf/nine_case_bse_nmpc_75pct.pdf"
    target.parent.mkdir(parents=True, exist_ok=True)
    pdf = canvas.Canvas(str(target), pagesize=(900, 650))
    pdf.setTitle("Nine-condition BSE and 75%-peak NMPC | 528, 725, 730")
    for path in files:
        reader = ImageReader(str(path))
        w, h = reader.getSize()
        scale = min(880 / w, 620 / h)
        pdf.drawImage(
            reader,
            (900 - w * scale) / 2,
            (650 - h * scale) / 2 + 5,
            width=w * scale,
            height=h * scale,
        )
        pdf.setFont("Helvetica", 7)
        pdf.drawString(15, 10, path.stem)
        pdf.drawRightString(
            885,
            10,
            "Calibration/replay simulation; not independent physiological validation",
        )
        pdf.showPage()
    pdf.save()
    (DOCS / "report_manifest.json").write_text(
        json.dumps(dict(pdf=str(target), figures=[str(p) for p in files]), indent=2)
        + "\n"
    )
    print(
        json.dumps(
            dict(
                pdf=str(target),
                pages=len(files),
                allPeakTrackingPassed=summary["allPeakTrackingPassed"],
                allBSEShapeIndicatorsPassed=summary["allBSEShapeIndicatorsPassed"],
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
