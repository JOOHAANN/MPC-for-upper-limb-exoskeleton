"""Publish the quality-first study without overwriting the historical study."""

import hashlib
import json
import os
from pathlib import Path
import shutil

import pymupdf as fitz
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from scipy.io import loadmat

from audit_subjects import ROOT
from fit_subject import many

BASE = Path(
    os.environ.get(
        "EXO_SUBJECT_RESULTS_ROOT", ROOT / "results/multi_subject_quality_v2"
    )
)
DEST = ROOT / "docs/subject_quality_v2"
SUBJECTS = ("240617Pilot", "240627Pilot")


def save(fig, name):
    fig.tight_layout()
    fig.savefig(DEST / (name + ".png"), dpi=155, bbox_inches="tight")
    fig.savefig(DEST / (name + ".pdf"), bbox_inches="tight")
    plt.close(fig)


def main():
    global SUBJECTS
    SUBJECTS = tuple(json.loads((DEST / "summary.json").read_text())["subjects"])
    if not SUBJECTS:
        raise RuntimeError(
            "No qualified, calibrated and validated subjects to publish."
        )
    strict_dir = (
        BASE / "strict_confirmation"
        if (BASE / "strict_confirmation/quality_gate.json").exists()
        else BASE
    )
    DEST.mkdir(exist_ok=True, parents=True)
    original = json.loads((BASE / "quality_gate.json").read_text())
    strict = json.loads((strict_dir / "quality_gate.json").read_text())
    fingerprints = lambda payload: {
        (r["subject"], r["trial"]): r.get("fileSha256", {}) for r in payload["trials"]
    }
    if fingerprints(original) != fingerprints(strict):
        strict_dir = BASE
        strict = original
    subjects = set(
        r["subject"] for r in strict["summary"] if r["subject"] != "240528Pilot"
    )
    ids = lambda p: {
        (c["subject"], c["id"])
        for c in p["cycles"]
        if c["accepted"] and c["subject"] in subjects
    }
    assert ids(original) == ids(strict)
    diff = 0.0
    for row in strict["trials"]:
        if row["subject"] in subjects and row["gate"] == "pass":
            a = np.load(BASE / "signals" / row["signalCache"])
            b = np.load(strict_dir / "signals" / row["signalCache"])
            for key in a.files:
                diff = max(diff, float(np.max(abs(a[key] - b[key]))))
    assert diff == 0.0
    confirmation = dict(
        acceptedSubjectCycleIdsIdentical=True,
        maximumAuditedSignalDifference=diff,
        gateAtCalibrationSHA256=hashlib.sha256(
            (BASE / "quality_gate.json").read_bytes()
        ).hexdigest(),
        strictGateSHA256=hashlib.sha256(
            (strict_dir / "quality_gate.json").read_bytes()
        ).hexdigest(),
    )
    (DEST / "strict_gate_confirmation.json").write_text(
        json.dumps(confirmation, indent=2) + "\n"
    )
    for name in (
        "quality_gate.json",
        "data_quality_audit_v2.pdf",
        "cycle_decisions.csv",
        "trial_clocks.csv",
        "quality_overview.png",
    ):
        shutil.copy2(strict_dir / name, DEST / name)
    shutil.copy2(BASE / "quality_gate.json", DEST / "gate_at_calibration.json")
    for name in ("benchmark_confirmation.json",):
        if (BASE / name).exists():
            shutil.copy2(BASE / name, DEST / name)
    for path in (BASE / "forensics_before").glob("*"):
        if path.suffix in (".pdf", ".png", ".json"):
            shutil.copy2(path, DEST / path.name)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    kept = [c for c in strict["cycles"] if c["accepted"] and c["subject"] in subjects]
    benchmark = [c["highFrequencyFraction"] * 100 for c in strict["baseline"]]
    groups = [benchmark]
    labels = ["240528\narchived successful"]
    for subject in SUBJECTS:
        records = [c for c in kept if c["subject"] == subject]
        for key, label in (
            ("smoothnessBefore", "unfiltered"),
            ("smoothnessAfter", "matched Gaussian"),
        ):
            groups.append([100 * c[key]["highFrequencyFraction"] for c in records])
            labels.append(subject[:6] + "\n" + label)
    axes[0].boxplot(groups, tick_labels=labels, showfliers=True)
    axes[0].axhline(
        strict["smoothnessUpperLimits"]["highFrequencyFraction"] * 100,
        color="red",
        ls="--",
        label="2 x original maximum",
    )
    axes[0].set_ylabel("Activation power above 3 Hz (%)")
    axes[0].legend(fontsize=8)
    groups = [[c["normalizedSlopeP99"] for c in strict["baseline"]]]
    for subject in SUBJECTS:
        records = [c for c in kept if c["subject"] == subject]
        groups.extend(
            [
                [c[key]["normalizedSlopeP99"] for c in records]
                for key in ("smoothnessBefore", "smoothnessAfter")
            ]
        )
    axes[1].boxplot(groups, tick_labels=labels, showfliers=True)
    axes[1].axhline(
        strict["smoothnessUpperLimits"]["normalizedSlopeP99"], color="red", ls="--"
    )
    axes[1].set_ylabel("Amplitude-normalized slope, 99th percentile (1/s)")
    for ax in axes:
        ax.grid(axis="y", alpha=0.2)
        ax.tick_params(axis="x", labelsize=8)
    fig.suptitle(
        "Qualified actions: activation compared with the ORIGINAL successful subject"
    )
    save(fig, "activation_benchmark")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    trials = [
        r for r in strict["trials"] if r["subject"] in SUBJECTS and r["gate"] == "pass"
    ]
    axes[0].bar(
        np.arange(len(trials)),
        [
            1000 * r["clocks"]["Activations.sto"]["anchorNearestP95Seconds"]
            for r in trials
        ],
        color="#1469b0",
    )
    axes[0].set(
        xlabel="Qualified recording index",
        ylabel="Activation to LMT nearest-grid P95 (ms)",
        title="Asynchronous sampling grids, common exported clock",
    )
    axes[0].text(
        0.02,
        0.94,
        "Skeletal clocks: exact matching timestamps\nActivation/emg clock agreement: 100%",
        transform=axes[0].transAxes,
        va="top",
        fontsize=9,
    )
    axes[1].bar(
        np.arange(len(kept)),
        [1000 * c["maximumGapSeconds"] for c in kept],
        color="#23844a",
    )
    axes[1].axhline(250, color="red", ls="--", label="250 ms rejection threshold")
    axes[1].set(
        xlabel="Qualified action index (including the lone 240719 action)",
        ylabel="Maximum native support gap (ms)",
    )
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.grid(axis="y", alpha=0.2)
    fig.suptitle(
        "Native-clock checks BEFORE interpolation | grid distance is not hardware latency"
    )
    save(fig, "timestamp_alignment")

    models = []
    boundary_rows = []
    fig, axes = plt.subplots(
        3, len(SUBJECTS), figsize=(6 * len(SUBJECTS), 8), sharex="col", squeeze=False
    )
    for col, subject in enumerate(SUBJECTS):
        folder = BASE / subject
        result = loadmat(folder / "hill_mpc.mat", simplify_cells=True)["result"]
        cache = loadmat(folder / "neural_models.mat", simplify_cells=True)
        metadata = json.loads((folder / "mpc_evaluation.json").read_text())
        validation = json.loads((folder / "validation.json").read_text())
        assert validation["allStandingChecksPassed"]
        choices = []
        for group in many(result["groups"])[:-1]:
            for source, p in enumerate(many(group["predictions"])):
                n = many(group["nominal"])[source]
                choices.append((float(max(abs(p["boundaryTorque"]))), group, source))
                boundary_rows.append(
                    dict(
                        subject=subject,
                        group=int(group["group"]),
                        source=source,
                        assistedPeakNm=float(max(abs(p["boundaryTorque"]))),
                        nominalPeakNm=float(max(abs(n["boundaryTorque"]))),
                        assistedSecondsAbove10Nm=float(
                            sum(abs(p["boundaryTorque"]) > 10) * 0.01
                        ),
                        assistedImpulseFractionOfBack=float(
                            sum(abs(p["boundaryTorque"]))
                            / max(sum(abs(p["backTorque"])), 1e-9)
                        ),
                    )
                )
        _, group, source = max(choices, key=lambda x: x[0])
        s = group["sequence"]
        p = many(group["predictions"])[source]
        n = many(group["nominal"])[source]
        t = s["t"]
        for row, key in ((0, "q"), (1, "boundaryTorque"), (2, "backTorque")):
            axes[row, col].plot(t, n[key], color="#c8841c", label="No assistance")
            axes[row, col].plot(t, p[key], color="#1469b0", label="MPC")
            axes[row, col].grid(alpha=0.2)
        axes[0, col].plot(
            t, s["q"], color="red", ls="--", lw=0.8, label="Original Scale reference"
        )
        axes[0, col].set_title(subject + " | worst assisted soft-boundary recording")
        axes[0, col].legend(fontsize=7)
        axes[1, col].set_title(
            f"Peak: {max(abs(p['boundaryTorque'])):.2f} N m", fontsize=9
        )
        axes[2, col].set_xlabel("Continuous replay time (s)")
        real = [r for r in metadata["metrics"] if not r["displayOnly"]]
        models.append(
            dict(
                subject=subject,
                caseSourcePairs=len(real),
                maximumExcessNm=max(r["excessPeakNm"] for r in real),
                secondsAboveThresholdPlus1Nm=sum(
                    r["excessOver1NmSeconds"] for r in real
                ),
                maximumAbdomenBackMeanRatio=max(r["abdomenBackRatio"] for r in real),
                minimumAngleRad=min(r["minAngleRad"] for r in real),
                maximumAngleRad=max(r["maxAngleRad"] for r in real),
                maximumBoundaryTorqueNm=max(r["boundaryPeakNm"] for r in real),
                bodyMassKg=float(cache["parameters"][0, 18]),
                bodyLeverM=float(cache["parameters"][0, 19]),
                constraintHorizonSeconds=metadata["options"].get(
                    "torqueConstraintHorizonSeconds", 0.4
                ),
                valleyCountMismatch=sum(
                    r["predictedValleys"] != r["referenceValleys"] for r in real
                ),
                qpFallbacks=sum(r["qpRejected"] for r in metadata["checks"]),
                modelSHA256=metadata["modelSHA256"],
                qualityGateAtFit=json.loads((folder / "fit_policy.json").read_text())[
                    "qualityGateSHA256"
                ],
            )
        )
        shutil.copy2(
            folder / "prediction_guard.json",
            DEST / (subject + "_prediction_guard.json"),
        )
    for row, label in enumerate(
        ("Theta (rad)", "Soft-boundary torque (N m)", "Back torque (N m)")
    ):
        axes[row, 0].set_ylabel(label)
    fig.suptitle(
        "Soft-boundary exposure is reported separately from muscle and actuator torque"
    )
    save(fig, "boundary_audit")
    (DEST / "boundary_exposure.json").write_text(
        json.dumps(boundary_rows, indent=2) + "\n"
    )
    payload = dict(
        quality=strict["summary"],
        models=models,
        strictConfirmation=confirmation,
        completedSubjects=list(SUBJECTS),
        notCalibrated=sorted(subjects - set(SUBJECTS)),
        independentExternalValidation=False,
        hardwareSynchronizationProven=False,
        gaussianMATLABMaximumDifference=3.0531133177191805e-16,
    )
    (DEST / "quality_first_summary.json").write_text(
        json.dumps(payload, indent=2) + "\n"
    )
    with fitz.open() as document:
        with fitz.open(DEST / "data_quality_audit_v2.pdf") as part:
            document.insert_pdf(part, from_page=0, to_page=0)
        for name in (
            "timestamp_alignment.pdf",
            "activation_benchmark.pdf",
            "boundary_audit.pdf",
            "multi_subject_results.pdf",
            "240617Pilot_exceedance_forensics.pdf",
            "240627Pilot_exceedance_forensics.pdf",
            "240719Pilot_exceedance_forensics.pdf",
        ):
            if (DEST / name).exists():
                with fitz.open(DEST / name) as part:
                    document.insert_pdf(part)
        document.save(DEST / "quality_first_results.pdf", garbage=4, deflate=True)
    print(json.dumps(payload, indent=2))
    manifest = json.loads((DEST / "export_metadata.json").read_text())
    manifest["qualityBundle"] = "quality_first_results.pdf"
    manifest["sourceFiles"] = [
        {
            "path": str(p.relative_to(ROOT)),
            "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
        }
        for p in sorted((ROOT / "experiments/multi_subject").glob("*.py"))
    ]
    manifest["files"] = [
        {"file": p.name, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
        for p in sorted(DEST.iterdir())
        if p.is_file() and p.name != "export_metadata.json"
    ]
    (DEST / "export_metadata.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
