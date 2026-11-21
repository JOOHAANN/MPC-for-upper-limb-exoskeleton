"""Export matched before/after MPC figures and complete threshold statistics."""

import csv
import hashlib
import json
from pathlib import Path
import shutil

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from scipy.io import loadmat

from audit_subjects import ROOT
from fit_subject import many
from refine_threshold import BASE, DEST, tracking_metrics
from run_hill_mpc import OUTPUT_NAMES
from run_threshold import parameter_fingerprint

PUBLISH = ROOT / "docs/mpc_threshold_tracking"
LABELS = ("Experimental activation", "MAP activation")
BLUE, OLD, RED, GREEN, PURPLE = "#1469b0", "#c8841c", "#cc3434", "#23844a", "#8054a6"


def matrix(prediction):
    return np.column_stack([prediction[key] for key in OUTPUT_NAMES])


def comparison(subject, old, new):
    before = many(old["groups"])[-1]
    after = many(new["groups"])[-1]
    seq = after["sequence"]
    t = seq["t"]
    fig, axes = plt.subplots(
        5,
        2,
        figsize=(13, 11.5),
        sharex=True,
        gridspec_kw={"height_ratios": [1.25, 1.0, 0.7, 0.6, 0.5]},
    )
    for source in (0, 1):
        p = many(before["predictions"])[source]
        q = many(after["predictions"])[source]
        nominal = many(after["nominal"])[source]
        axes[0, source].plot(
            t,
            nominal["backTorque"],
            color=BLUE,
            alpha=0.35,
            ls=":",
            lw=1.0,
            label="Without assistance",
        )
        axes[0, source].plot(
            t,
            p["backTorque"],
            color=OLD,
            ls="--",
            lw=1.3,
            label="Previous conservative MPC",
        )
        axes[0, source].plot(
            t, q["backTorque"], color=BLUE, lw=1.4, label="Revised MPC"
        )
        for seg in many(seq["segments"]):
            lo, hi = int(seg["startIndex"]) - 1, int(seg["endIndex"]) - 1
            limit = np.interp(seg["loadKg"], [5, 10, 15], [100, 140, 180])
            axes[0, source].plot(
                t[[lo, hi]], [limit, limit], color=".25", ls="--", lw=1.1
            )
            axes[0, source].fill_between(
                t[lo : hi + 1], limit - 5, limit, color=GREEN, alpha=0.14
            )
            axes[0, source].annotate(
                f"{seg['loadKg']:g} kg / {limit:g} N m",
                ((t[lo] + t[hi]) / 2, 1.01),
                xycoords=("data", "axes fraction"),
                ha="center",
                fontsize=8,
            )
        axes[0, source].set_title(LABELS[source], pad=28)
        axes[1, source].plot(
            t, seq["q"], color=RED, ls="--", lw=0.9, label="Original Scale reference"
        )
        for row, key in (
            (1, "q"),
            (2, "assistanceTorque"),
            (3, "frontTorque"),
            (4, "boundaryTorque"),
        ):
            axes[row, source].plot(t, p[key], color=OLD, ls="--", lw=1.1)
            axes[row, source].plot(
                t, q[key], color=PURPLE if row >= 3 else BLUE, lw=1.2
            )
        axes[2, source].axhline(80, color=RED, ls="--", lw=0.8)
        for row in range(5):
            axes[row, source].grid(alpha=0.2)
        axes[4, source].set_xlabel("Continuous time (s)")
        axes[4, source].set_ylim(
            -1, max(1.0, float(max(abs(q["boundaryTorque"])))) * 1.15
        )
    for row, label in enumerate(
        (
            "Back torque (N m)",
            "Theta (rad)",
            "Assistance (N m)",
            "Abdominal torque\n(N m)",
            "Soft boundary\n(N m)",
        )
    ):
        axes[row, 0].set_ylabel(label, fontsize=10)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    handles2, labels2 = axes[1, 0].get_legend_handles_labels()
    fig.legend(
        handles + handles2,
        labels + labels2,
        loc="lower center",
        ncol=4,
        fontsize=8,
        frameon=False,
        bbox_to_anchor=(0.5, 0.02),
    )
    fig.suptitle(
        subject + " | less assistance, back torque near its upper target", fontsize=13
    )
    fig.text(
        0.5,
        0.005,
        "Black dashed: back-torque limits. Green: 0-5 N m below the limit. Low-demand portions are not forced upward.",
        ha="center",
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.065, 1, 0.945), h_pad=0.7)
    return fig


def torque_panel(subject, old, new):
    before = many(old["groups"])[-1]
    after = many(new["groups"])[-1]
    seq = after["sequence"]
    t = seq["t"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 3.8), sharex=True, sharey=True)
    for source, ax in enumerate(axes):
        p = many(before["predictions"])[source]
        q = many(after["predictions"])[source]
        n = many(after["nominal"])[source]
        ax.plot(
            t,
            n["backTorque"],
            color=BLUE,
            ls=":",
            alpha=0.35,
            lw=1,
            label="Without assistance",
        )
        ax.plot(t, p["backTorque"], color=OLD, ls="--", lw=1.2, label="Previous MPC")
        ax.plot(t, q["backTorque"], color=BLUE, lw=1.7, label="Revised MPC")
        ax.plot(t, q["frontTorque"], color=PURPLE, lw=0.9, label="Abdomen, revised MPC")
        for seg in many(seq["segments"]):
            lo, hi = int(seg["startIndex"]) - 1, int(seg["endIndex"]) - 1
            limit = np.interp(seg["loadKg"], [5, 10, 15], [100, 140, 180])
            ax.plot(t[[lo, hi]], [limit, limit], color=".25", ls="--", lw=1.1)
            ax.fill_between(t[lo : hi + 1], limit - 5, limit, color=GREEN, alpha=0.16)
            ax.annotate(
                f"{seg['loadKg']:g} kg",
                ((t[lo] + t[hi]) / 2, 0.96),
                xycoords=("data", "axes fraction"),
                ha="center",
                va="top",
                fontsize=9,
            )
        ax.set_title(LABELS[source])
        ax.set_xlabel("Time (s)")
        ax.grid(alpha=0.2)
    axes[0].set_ylabel("Back / abdomen torque (N m)")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, frameon=False, fontsize=8)
    fig.suptitle(
        subject + " | black dashed = threshold; green band = 0-5 N m below", fontsize=12
    )
    fig.tight_layout(rect=(0, 0.09, 1, 0.93))
    return fig


def main():
    PUBLISH.mkdir(parents=True, exist_ok=True)
    subjects = [
        p.name
        for p in sorted(DEST.iterdir())
        if p.is_dir() and (p / "mpc_evaluation.json").exists()
    ]
    rows = []
    datasets = []
    for subject in subjects:
        folder = DEST / subject
        before = loadmat(BASE / subject / "hill_mpc.mat", simplify_cells=True)["result"]
        after = loadmat(folder / "hill_mpc.mat", simplify_cells=True)["result"]
        validation = json.loads((folder / "validation.json").read_text())
        assert validation["allStandingChecksPassed"]
        old_cache = (BASE / subject / "neural_models.mat").read_bytes()
        new_cache = (folder / "neural_models.mat").read_bytes()
        assert old_cache == new_cache
        for old_group, new_group in zip(many(before["groups"]), many(after["groups"])):
            for key in ("q", "torque", "m", "waiting"):
                assert np.array_equal(
                    old_group["sequence"][key], new_group["sequence"][key]
                )
            assert np.array_equal(old_group["inputs"], new_group["inputs"])
            if new_group["sequence"]["name"] == "continuous_available_loads":
                continue
            for source in (0, 1):
                seq = new_group["sequence"]
                nominal = matrix(many(new_group["nominal"])[source])
                old_rows = tracking_metrics(
                    matrix(many(old_group["predictions"])[source]), nominal, seq, source
                )
                new_rows = tracking_metrics(
                    matrix(many(new_group["predictions"])[source]), nominal, seq, source
                )
                for a, b in zip(old_rows, new_rows):
                    rows.append(
                        dict(
                            subject=subject,
                            source=LABELS[source],
                            recordId=b["recordId"],
                            loadKg=b["loadKg"],
                            split=b["split"],
                            thresholdNm=b["backThresholdNm"],
                            nominalPeakNm=b["nominalBackPeakNm"],
                            oldPeakNm=a["backPeakNm"],
                            newPeakNm=b["backPeakNm"],
                            oldPeakShortfallNm=a["peakShortfallNm"],
                            newPeakShortfallNm=b["peakShortfallNm"],
                            newExcessPeakNm=b["excessPeakNm"],
                            demandSeconds=b["demandSeconds"],
                            oldDemandShortfallNm=a["demandMeanShortfallNm"],
                            newDemandShortfallNm=b["demandMeanShortfallNm"],
                            oldAssistanceImpulseNms=a["assistanceImpulseNms"],
                            newAssistanceImpulseNms=b["assistanceImpulseNms"],
                            oldThetaCorrelation=a["shapeCorrelation"],
                            newThetaCorrelation=b["shapeCorrelation"],
                            referenceValleys=b["referenceValleys"],
                            nominalValleys=b["nominalValleys"],
                            newValleys=b["predictedValleys"],
                            newMinimumAngleRad=b["minAngleRad"],
                            newMaximumAngleRad=b["maxAngleRad"],
                            newBoundaryPeakNm=b["boundaryPeakNm"],
                            newAbdomenPeakNm=b["abdomenPeakNm"],
                            newAbdomenBackMeanRatio=b["abdomenBackRatio"],
                        )
                    )
        datasets.append((subject, before, after))
        tuning = json.loads((folder / "selected_threshold_options.json").read_text())
        tuning["modelParameterSHA256"] = parameter_fingerprint(folder)
        (PUBLISH / (subject + "_options.json")).write_text(
            json.dumps(tuning, indent=2) + "\n"
        )
        for name in (
            "threshold_candidates.json",
            "validation.json",
            "mpc_evaluation.json",
        ):
            shutil.copy2(folder / name, PUBLISH / (subject + "_" + name))
    with (PUBLISH / "all_cases_before_after.csv").open(
        "w", newline="", encoding="utf-8-sig"
    ) as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summaries = []
    for subject, before, after in datasets:
        chosen = [r for r in rows if r["subject"] == subject]
        demand = [r for r in chosen if r["newPeakShortfallNm"] is not None]
        old_impulse = sum(r["oldAssistanceImpulseNms"] for r in chosen)
        new_impulse = sum(r["newAssistanceImpulseNms"] for r in chosen)
        validation = json.loads((DEST / subject / "validation.json").read_text())
        summaries.append(
            dict(
                subject=subject,
                caseSourcePairs=len(chosen),
                demandCaseSourcePairs=len(demand),
                oldMeanPeakShortfallNm=float(
                    np.mean([r["oldPeakShortfallNm"] for r in demand])
                ),
                newMeanPeakShortfallNm=float(
                    np.mean([r["newPeakShortfallNm"] for r in demand])
                ),
                maximumExcessNm=max(r["newExcessPeakNm"] for r in chosen),
                assistanceImpulseReductionPercent=100
                * (1 - new_impulse / max(old_impulse, 1e-9)),
                newValleyMismatchCount=sum(
                    r["newValleys"] != r["referenceValleys"] for r in chosen
                ),
                maximumBoundaryTorqueNm=max(r["newBoundaryPeakNm"] for r in chosen),
                minimumAngleRad=min(r["newMinimumAngleRad"] for r in chosen),
                maximumAngleRad=max(r["newMaximumAngleRad"] for r in chosen),
                modelSHA256=hashlib.sha256(
                    (DEST / subject / "neural_models.mat").read_bytes()
                ).hexdigest(),
                options=after["options"],
                allStandingChecksPassed=validation["allStandingChecksPassed"],
            )
        )
    fig, axes = plt.subplots(
        1, len(subjects), figsize=(6 * len(subjects), 4.4), squeeze=False
    )
    for ax, (subject, _, _) in zip(axes[0], datasets):
        values = [
            r
            for r in rows
            if r["subject"] == subject and r["newPeakShortfallNm"] is not None
        ]
        for source, offset in ((LABELS[0], -1.7), (LABELS[1], 1.7)):
            group = [r for r in values if r["source"] == source]
            ax.scatter(
                [r["thresholdNm"] + offset for r in group],
                [r["oldPeakNm"] for r in group],
                color=OLD,
                marker="x",
                s=30,
                label="Previous MPC" if offset < 0 else None,
            )
            ax.scatter(
                [r["thresholdNm"] + offset for r in group],
                [r["newPeakNm"] for r in group],
                color=BLUE,
                s=24,
                label="Revised MPC" if offset < 0 else None,
            )
        x = np.array([90, 190])
        ax.plot(x, x, color=".25", ls="--", label="Threshold")
        ax.fill_between(x, x - 5, x, color=GREEN, alpha=0.12, label="0-5 N m below")
        ax.set(
            title=subject,
            xlabel="Back-torque threshold (N m)",
            ylabel="Action peak back torque (N m)",
            xlim=(90, 190),
            ylim=(75, 193),
        )
        ax.set_xticks([100, 140, 180])
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8, loc="upper left")
    fig.suptitle(
        "All qualified actions and both activation sources | only cases needing assistance"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(PUBLISH / "threshold_overview.png", dpi=165)
    with PdfPages(PUBLISH / "threshold_tracking_results.pdf") as pdf:
        pdf.savefig(fig)
        plt.close(fig)
        for subject, before, after in datasets:
            fig = comparison(subject, before, after)
            fig.savefig(
                PUBLISH / (subject + "_before_after.png"), dpi=150, bbox_inches="tight"
            )
            fig.savefig(PUBLISH / (subject + "_before_after.pdf"), bbox_inches="tight")
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
            fig = torque_panel(subject, before, after)
            fig.savefig(
                PUBLISH / (subject + "_torque_threshold.png"),
                dpi=170,
                bbox_inches="tight",
            )
            fig.savefig(
                PUBLISH / (subject + "_torque_threshold.pdf"), bbox_inches="tight"
            )
            plt.close(fig)

    def convert(value):
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        raise TypeError(type(value))

    summary = dict(
        subjects=subjects,
        models=summaries,
        bseAndMapUnchanged=True,
        referencesAndInputsUnchanged=True,
        controllerSelectionUsesOnlyTraining=True,
        holdoutIsDevelopmentValidation=True,
        thresholdTrackingNotAUniversalSafetyGuarantee=True,
    )
    (PUBLISH / "summary.json").write_text(
        json.dumps(summary, indent=2, default=convert) + "\n"
    )
    source_files = [
        ROOT / "experiments/legacy_mpc/run_hill_mpc.py",
        ROOT / "experiments/legacy_mpc/hill_mpc_state.c",
        ROOT / "experiments/hill_compliant/native_library.py",
        *sorted((ROOT / "experiments/hill_compliant").glob("hill_*core.c")),
        *sorted((ROOT / "experiments/multi_subject").glob("*.py")),
        ROOT / "demo_hill_subjects.m",
        ROOT / "tests/test_subject_quality.py",
    ]
    manifest = dict(
        bundle="threshold_tracking_results.pdf",
        sources=[
            dict(
                path=str(p.relative_to(ROOT)),
                sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
            )
            for p in source_files
        ],
        files=[
            dict(file=p.name, sha256=hashlib.sha256(p.read_bytes()).hexdigest())
            for p in sorted(PUBLISH.iterdir())
            if p.is_file() and p.name != "export_metadata.json"
        ],
    )
    (PUBLISH / "export_metadata.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(summary, indent=2, default=convert))


if __name__ == "__main__":
    main()
