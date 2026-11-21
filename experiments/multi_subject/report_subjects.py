"""Export reproducible scientific figures from frozen subject results."""

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

from audit_subjects import ROOT, OUT
from fit_subject import many
from fit_hill_compliant import PARAMETER_NAMES, SHARED

BLUE, AMBER, RED, GREEN, PURPLE = "#1469b0", "#c8841c", "#cc3434", "#23844a", "#8054a6"
LABELS = ("Experimental activation", "MAP activation")
PUBLISH = ROOT / "docs/multi_subject_evaluation"


def mean(rows, key):
    return float(np.mean([r[key] for r in rows])) if rows else None


def records():
    data = []
    for folder in sorted(OUT.iterdir()):
        if not folder.is_dir() or not (folder / "mpc_evaluation.json").exists():
            continue
        policy = json.loads((folder / "fit_policy.json").read_text())
        if (OUT / "quality_gate.json").exists() and policy.get(
            "qualityGateSHA256"
        ) != hashlib.sha256((OUT / "quality_gate.json").read_bytes()).hexdigest():
            print("Not publishing obsolete calibration:", folder.name, flush=True)
            continue
        assert (
            policy["preparedSHA256"]
            == hashlib.sha256((folder / "prepared.mat").read_bytes()).hexdigest()
        )
        r = loadmat(folder / "hill_mpc.mat", simplify_cells=True)["result"]
        cache = loadmat(folder / "neural_models.mat", simplify_cells=True)
        meta = json.loads((folder / "mpc_evaluation.json").read_text())
        assert (
            meta["modelSHA256"]
            == hashlib.sha256((folder / "neural_models.mat").read_bytes()).hexdigest()
        )
        data.append(
            {
                "subject": folder.name,
                "r": r,
                "cache": cache,
                "meta": meta,
                "prep": json.loads((folder / "preparation.json").read_text()),
                "bse": json.loads((folder / "bse_metrics.json").read_text()),
                "policy": policy,
                "validation": (
                    json.loads((folder / "validation.json").read_text())
                    if (folder / "validation.json").exists()
                    else None
                ),
            }
        )
    return data


def reference(ax, s, key):
    wait = np.asarray(s["waiting"], dtype=bool)
    ax.plot(
        s["t"],
        np.where(~wait, s[key], np.nan),
        color=RED,
        ls="--",
        lw=1.1,
        label="Recorded reference",
    )
    ax.plot(
        s["t"],
        np.where(wait, s[key], np.nan),
        color=RED,
        ls=":",
        lw=0.8,
        alpha=0.65,
        label="Filled waiting reference",
    )


def sequence_figure(item):
    g = many(item["r"]["groups"])[-1]
    s = g["sequence"]
    t = s["t"]
    opts = item["r"]["options"]
    fig, axes = plt.subplots(
        5,
        2,
        figsize=(13.6, 12),
        sharex=True,
        gridspec_kw={"height_ratios": [1.1, 1, 1, 0.6, 0.5]},
    )
    for source in (0, 1):
        p, p0 = many(g["predictions"])[source], many(g["nominal"])[source]
        reference(axes[0, source], s, "q")
        reference(axes[1, source], s, "torque")
        for row, key in ((0, "q"), (1, "torque")):
            axes[row, source].plot(
                t,
                p0[key],
                color=AMBER,
                lw=1.7,
                label="Personalized Hill, no assistance",
            )
            axes[row, source].plot(
                t, p[key], color=BLUE, lw=1.2, label="Same Hill + personalized MPC"
            )
        axes[1, source].plot(
            t,
            p["torque"] + p["assistanceTorque"],
            color=GREEN,
            ls=":",
            lw=1.2,
            label="MPC muscle + actuator",
        )
        axes[0, source].set_title(LABELS[source])
        if not any(bool(seg["angleUsable"]) for seg in many(s["segments"])):
            axes[0, source].text(
                0.02,
                0.04,
                "ANGLE REFERENCE FLAGGED: polarity/alignment",
                transform=axes[0, source].transAxes,
                color=RED,
                fontsize=8,
            )
        axes[2, source].plot(
            t,
            p0["backTorque"],
            color=BLUE,
            ls=":",
            alpha=0.5,
            lw=1.1,
            label="Back, no assistance",
        )
        axes[2, source].plot(t, p["backTorque"], color=BLUE, lw=1.3, label="Back, MPC")
        axes[2, source].plot(
            t,
            p0["frontTorque"],
            color=PURPLE,
            ls=":",
            alpha=0.5,
            lw=1.1,
            label="Abdomen, no assistance",
        )
        axes[2, source].plot(
            t, p["frontTorque"], color=PURPLE, lw=1.3, label="Abdomen, MPC"
        )
        for seg in many(s["segments"]):
            lo, hi = int(seg["startIndex"]) - 1, int(seg["endIndex"]) - 1
            threshold = np.interp(seg["loadKg"], [5, 10, 15], [100, 140, 180])
            axes[2, source].plot(
                t[[lo, hi]], [threshold, threshold], color="#333333", ls="--", lw=1.0
            )
            axes[0, source].annotate(
                f"{seg['loadKg']:g} kg | {seg['split']}",
                ((t[lo] + t[hi]) / 2, 0.98),
                xycoords=("data", "axes fraction"),
                ha="center",
                va="top",
                fontsize=7.5,
            )
            for ax in axes[:, source]:
                ax.axvspan(t[lo], t[hi], color=AMBER, alpha=0.045, zorder=-2)
        axes[3, source].plot(t, p["assistanceTorque"], color=GREEN, lw=1.2)
        axes[3, source].axhline(
            80, color=RED, ls="--", lw=0.9, label="80 N·m actuator cap"
        )
        if opts["capNm"] < 80:
            axes[3, source].axhline(opts["capNm"], color=BLUE, ls="--", lw=0.9)
        axes[3, source].set_ylim(-2, 88)
        axes[4, source].plot(t, p["boundaryTorque"], color=PURPLE, lw=1.2)
        peak = float(np.max(np.abs(p["boundaryTorque"])))
        boundary_label = f"Soft-boundary peak: {peak:.2f} N·m"
        boundary_scale = peak
        if item.get("showUnassistedBoundary", False):
            baseline_boundary = float(np.max(np.abs(p0["boundaryTorque"])))
            axes[4, source].plot(t, p0["boundaryTorque"], color=PURPLE, ls=":", lw=1.1)
            boundary_label = f"Boundary peak: MPC {peak:.2f}; no assistance {baseline_boundary:.2f} N·m"
            boundary_scale = max(peak, baseline_boundary)
        axes[4, source].text(
            0.015,
            0.9,
            boundary_label,
            transform=axes[4, source].transAxes,
            va="top",
            fontsize=8,
        )
        if boundary_scale < 0.1:
            axes[4, source].set_ylim(-1, 1)
        for ax, title in zip(
            axes[:, source],
            (
                "Theta (rad)",
                "Net muscle\ntorque (N·m)",
                "Back / abdomen\ntorque (N·m)",
                "Assistance\n(N·m)",
                "Soft boundary\n(N·m)",
            ),
        ):
            ax.set_ylabel(title)
            ax.grid(alpha=0.2)
        axes[-1, source].set_xlabel("Continuous time (s)")
    h, l = axes[0, 0].get_legend_handles_labels()
    h2, l2 = axes[2, 0].get_legend_handles_labels()
    h3, l3 = axes[1, 0].get_legend_handles_labels()
    fig.legend(
        h + h2 + h3[-1:],
        l + l2 + l3[-1:],
        ncol=4,
        loc="lower center",
        fontsize=8,
        frameon=False,
        bbox_to_anchor=(0.5, 0.025),
    )
    loads = " → ".join(str(int(seg["loadKg"])) for seg in many(s["segments"]))
    fig.suptitle(
        item["subject"]
        + f" | {loads} kg with standing waits | quality-screened recorded actions",
        fontsize=13,
    )
    fig.text(
        0.5,
        0.008,
        "Net muscle torque = back − abdomen. Dashed black lines: BACK soft targets (100 / 140 / 180 N·m). Missing loads are not synthesized.",
        ha="center",
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.09, 1, 0.965), h_pad=1.0)
    return fig


def activation_figure(item):
    g = many(item["r"]["groups"])[-1]
    s = g["sequence"]
    segments = many(s["segments"])
    fig, axes = plt.subplots(
        len(segments), 2, figsize=(12, 2.7 * len(segments) + 1), squeeze=False
    )
    for i, seg in enumerate(segments):
        ids = slice(int(seg["startIndex"]) - 1, int(seg["endIndex"]))
        t = s["t"][ids] - s["t"][ids.start]
        for channel in (0, 1):
            ax = axes[i, channel]
            for source, color in ((0, BLUE), (1, AMBER)):
                ax.plot(
                    t,
                    np.asarray(g["inputs"])[source, ids, channel],
                    color=color,
                    lw=1.2,
                    label=LABELS[source],
                )
            ax.set_title(
                f"{seg['loadKg']:g} kg | "
                + ("Back input" if channel == 0 else "Abdominal input")
                + f" | {seg['split']}",
                fontsize=10,
            )
            ax.set_ylabel("Activation (0–1, not %MVC)")
            ax.grid(alpha=0.2)
            ax.set_xlabel("Action time (s)")
    axes[0, 0].legend(frameon=False, fontsize=9)
    fig.suptitle(
        item["subject"]
        + " | experimental and generated neural inputs before Hill feedback"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return fig


def condition_cases(item):
    chosen = {}
    for g in many(item["r"]["groups"])[:-1]:
        for seg in many(g["sequence"]["segments"]):
            key = (float(seg["loadKg"]), int(seg["speedCode"]))
            # Prefer a held-out record when illustrating an observed condition.
            rank = 0 if str(seg["split"]).startswith("test") else 1
            if key not in chosen or rank < chosen[key][0]:
                chosen[key] = (rank, g, seg)
    if item["subject"] == "240617Pilot" and not item["prep"].get(
        "qualityCheckedBeforeTraining", False
    ):
        chosen = {k: v for k, v in chosen.items() if k[1] in (23, 33, 43)}
    return [chosen[k][1:] for k in sorted(chosen)]


def case_figure(item, source):
    cases = condition_cases(item)
    cols = min(3, len(cases))
    blocks = int(np.ceil(len(cases) / cols))
    fig, axes = plt.subplots(
        blocks * 2, cols, figsize=(13, blocks * 4.2 + 1), squeeze=False
    )
    for index, (g, seg) in enumerate(cases):
        block, col = divmod(index, cols)
        s = g["sequence"]
        ids = slice(int(seg["startIndex"]) - 1, int(seg["endIndex"]))
        t = s["t"][ids] - s["t"][ids.start]
        p = many(g["predictions"])[source]
        p0 = many(g["nominal"])[source]
        for row, key in ((2 * block, "q"), (2 * block + 1, "torque")):
            ax = axes[row, col]
            ax.plot(t, s[key][ids], color=RED, ls="--", lw=1, label="Reference")
            ax.plot(t, p0[key][ids], color=AMBER, lw=1.5, label="Personalized BSE/Hill")
            ax.plot(t, p[key][ids], color=BLUE, lw=1.1, label="BSE/Hill + MPC")
            ax.set_ylabel("Theta (rad)" if key == "q" else "Net torque (N·m)")
            ax.grid(alpha=0.2)
        speed = (
            f"code {int(seg['speedCode'])}"
            if seg["speedCode"] >= 0
            else "speed unspecified"
        )
        flag = " | theta flagged" if not bool(seg["angleUsable"]) else ""
        metadata = next(
            (c for c in item["prep"]["cycles"] if c["id"] == seg["recordId"]), {}
        )
        if metadata.get("rightCensored", False):
            flag += " | partial action"
        axes[2 * block, col].set_title(
            f"{seg['loadKg']:g} kg / {speed}\n{seg['split']}" + flag, fontsize=9
        )
        axes[2 * block + 1, col].set_xlabel("Action time (s)")
    for i in range(len(cases), blocks * cols):
        b, c = divmod(i, cols)
        axes[2 * b, c].set_visible(False)
        axes[2 * b + 1, c].set_visible(False)
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, fontsize=9, frameon=False)
    fig.suptitle(
        item["subject"]
        + " | "
        + LABELS[source]
        + (
            " | nine target load/speed conditions"
            if item["subject"] == "240617Pilot"
            and not item["prep"].get("qualityCheckedBeforeTraining", False)
            else " | all available load/speed labels"
        )
    )
    fig.tight_layout(rect=(0, 0.045, 1, 0.965))
    return fig


def summary_rows(data):
    rows = []
    for item in data:
        prep = item["prep"]
        mpc = [r for r in item["meta"]["metrics"] if not r["displayOnly"]]
        for source in LABELS:
            for split in ("train", "test"):
                selected = [
                    r
                    for r in mpc
                    if r["source"] == source
                    and (r["split"] == "train") == (split == "train")
                ]
                q = [r for r in selected if r["angleReferenceUsable"]]
                b = [
                    r
                    for r in item["bse"]
                    if r["model"] == "personalized"
                    and r["source"] == source
                    and (r["split"] == "train") == (split == "train")
                ]
                amap = [
                    r
                    for r in prep["cycles"]
                    if (r["split"] == "train") == (split == "train")
                ]
                rows.append(
                    {
                        "subject": item["subject"],
                        "source": source,
                        "split": split,
                        "cycles": len(selected),
                        "angleUsableCycles": len(q),
                        "rightCensoredCycles": sum(
                            c.get("rightCensored", False) for c in amap
                        ),
                        "mapInputRMSE": mean(amap, "mapRMSE"),
                        "bseTorqueRMSENm": mean(b, "torqueRMSENm"),
                        "mpcThetaRMSEDegrees": mean(q, "thetaRMSEDegrees"),
                        "mpcThetaCorrelation": mean(q, "shapeCorrelation"),
                        "backMeanNm": mean(selected, "backMeanNm"),
                        "abdomenMeanNm": mean(selected, "abdomenMeanNm"),
                        "abdomenBackRatio": mean(selected, "abdomenBackRatio"),
                        "worstBackExcessNm": max(r["excessPeakNm"] for r in selected),
                        "excessOver1NmSeconds": sum(
                            r["excessOver1NmSeconds"] for r in selected
                        ),
                    }
                )
    return rows


def overview(data, rows):
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.4))
    x = np.arange(len(data))
    labels = [d["subject"].replace("Pilot", "") for d in data]
    axes[0].bar(
        x,
        [
            mean([c for c in d["prep"]["cycles"] if c["split"] != "train"], "mapRMSE")
            for d in data
        ],
        color=AMBER,
    )
    axes[0].set_ylabel("MAP activation RMSE (0–1)")
    axes[0].set_title("Held-out records / cycles (includes partial actions)")
    for source, color, offset in ((LABELS[0], BLUE, -0.18), (LABELS[1], GREEN, 0.18)):
        chosen = [
            next(
                r
                for r in rows
                if r["subject"] == d["subject"]
                and r["split"] == "test"
                and r["source"] == source
            )
            for d in data
        ]
        axes[1].bar(
            x + offset,
            [
                (
                    r["mpcThetaCorrelation"]
                    if r["mpcThetaCorrelation"] is not None
                    else np.nan
                )
                for r in chosen
            ],
            0.36,
            color=color,
            label=source,
        )
        axes[2].bar(
            x + offset, [r["abdomenBackRatio"] for r in chosen], 0.36, color=color
        )
    axes[1].set_ylabel("Mean theta shape correlation")
    axes[1].set_title("Qualified angle references")
    axes[1].set_ylim(
        min(
            -0.2,
            min(
                r["mpcThetaCorrelation"]
                for r in rows
                if r["mpcThetaCorrelation"] is not None
            )
            - 0.08,
        ),
        1,
    )
    for i, item in enumerate(data):
        if not any(c["angleUsable"] for c in item["prep"]["cycles"]):
            axes[1].text(i, 0.04, "N/A", ha="center", fontsize=9, color=RED)
    axes[2].set_ylabel("Abdominal / back mean torque")
    axes[2].set_title("Same Hill mechanics; subject-specific parameters")
    for ax in axes:
        ax.set_xticks(x, labels, rotation=25)
        ax.grid(axis="y", alpha=0.2)
    axes[1].legend(fontsize=8, frameon=False)
    fig.suptitle(
        "Personalized MAP + BSE/Hill + MPC | data limitations remain visible",
        fontsize=14,
    )
    fig.text(
        0.5,
        0.02,
        "Missing conditions are not simulated as measured validation. Activation is processed CEINMS output, not independently verified %MVC.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.08, 1, 0.93))
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--partial", action="store_true")
    parser.add_argument("--publish-directory", type=Path)
    parser.add_argument("--model-directory", type=Path)
    args = parser.parse_args()
    global PUBLISH
    if args.publish_directory is not None:
        PUBLISH = args.publish_directory
    data = records()
    if not data:
        raise RuntimeError(
            "No current quality-approved, calibrated subject results available."
        )
    if not args.partial:
        assert len(data) == 5 and all(
            d["validation"] for d in data
        ), "Finish all subject runs and numerical checks first."
    PUBLISH.mkdir(parents=True, exist_ok=True)
    rows = summary_rows(data)
    figures = []
    fig = overview(data, rows)
    fig.savefig(PUBLISH / "multi_subject_overview.png", dpi=170, bbox_inches="tight")
    figures.append(fig)
    with PdfPages(PUBLISH / "multi_subject_results.pdf") as bundle:
        bundle.savefig(fig, bbox_inches="tight")
        plt.close(fig)
        for item in data:
            subject = item["subject"]
            with PdfPages(PUBLISH / (subject + "_results.pdf")) as subject_pdf:
                pages = [
                    (sequence_figure(item), "continuous"),
                    (activation_figure(item), "activation"),
                    (case_figure(item, 0), "cases_experimental"),
                    (case_figure(item, 1), "cases_map"),
                ]
                for fig, slug in pages:
                    fig.savefig(
                        PUBLISH / (subject + "_" + slug + ".png"),
                        dpi=155,
                        bbox_inches="tight",
                    )
                    subject_pdf.savefig(fig, bbox_inches="tight")
                    bundle.savefig(fig, bbox_inches="tight")
                    plt.close(fig)
    for name, values in (
        ("summary_metrics.csv", rows),
        (
            "all_mpc_cases.csv",
            [
                {"subject": d["subject"], **r}
                for d in data
                for r in d["meta"]["metrics"]
                if not r["displayOnly"]
            ],
        ),
    ):
        with (PUBLISH / name).open("w", newline="", encoding="utf-8-sig") as stream:
            w = csv.DictWriter(stream, fieldnames=list(values[0]))
            w.writeheader()
            w.writerows(values)
    summary = {
        "subjects": [d["subject"] for d in data],
        "summary": rows,
        "quality": {
            d["subject"]: {
                k: d["prep"][k]
                for k in (
                    "trainingCycles",
                    "testCycles",
                    "knownConditions",
                    "angleUsableCycles",
                )
            }
            for d in data
        },
        "parameters": {
            d["subject"]: {
                "physical": d["cache"]["parameters"].tolist(),
                "neural": d["cache"]["neuralParameters"].tolist(),
                "mpc": d["meta"]["options"],
            }
            for d in data
        },
    }
    (PUBLISH / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    for name in ("data_quality.csv", "data_quality.json"):
        if (OUT / name).exists():
            shutil.copy2(OUT / name, PUBLISH / name)
    modeldir = args.model_directory or ROOT / "data/processed/multi_subject"
    modeldir.mkdir(parents=True, exist_ok=True)
    for item in data:
        subject = item["subject"]
        shutil.copy2(
            OUT / subject / "neural_models.mat", modeldir / (subject + "_models.mat")
        )
        for name in (
            "preparation.json",
            "fit_policy.json",
            "validation.json",
            "mpc_evaluation.json",
            "mpc_candidates.json",
            "bse_metrics.json",
        ):
            file = OUT / subject / name
            if file.exists():
                shutil.copy2(file, PUBLISH / (subject + "_" + name))
    parameter_rows = []
    neural_names = (
        "intentAmplitude",
        "intentExponent",
        "intentTimeConstant",
        "postureKp",
        "postureKd",
        "negativeRequestLimit",
        "feedbackBlend",
        "coactivation",
        "causalLeadTime",
        "tonicCompensation",
    )
    for item in data:
        for source in (0, 1):
            for index, name in enumerate(PARAMETER_NAMES):
                parameter_rows.append(
                    {
                        "subject": item["subject"],
                        "source": LABELS[source],
                        "family": "BSE/Hill",
                        "parameter": name,
                        "value": float(item["cache"]["parameters"][source, index]),
                        "sharedBetweenSources": index in SHARED,
                    }
                )
            for index, name in enumerate(neural_names):
                parameter_rows.append(
                    {
                        "subject": item["subject"],
                        "source": LABELS[source],
                        "family": "Neural",
                        "parameter": name,
                        "value": float(
                            item["cache"]["neuralParameters"][source, index]
                        ),
                        "sharedBetweenSources": False,
                    }
                )
    with (PUBLISH / "all_parameters.csv").open(
        "w", newline="", encoding="utf-8-sig"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(parameter_rows[0]))
        writer.writeheader()
        writer.writerows(parameter_rows)
    manifest = {
        "bundle": "multi_subject_results.pdf",
        "subjects": {
            d["subject"]: {
                "modelSHA256": d["meta"]["modelSHA256"],
                "resultSHA256": hashlib.sha256(
                    (OUT / d["subject"] / "hill_mpc.mat").read_bytes()
                ).hexdigest(),
            }
            for d in data
        },
        "files": [
            {"file": p.name, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
            for p in sorted(PUBLISH.iterdir())
            if p.is_file() and p.name != "export_metadata.json"
        ],
    }
    (PUBLISH / "export_metadata.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
