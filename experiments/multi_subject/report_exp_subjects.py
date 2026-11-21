"""Publish four-subject completion status, data exclusions and new Exp results."""

from __future__ import annotations

from collections import Counter
import json
import shutil

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

from audit_subjects import DATA, ROOT
from fit_subject import many
from quality_gate import processed_signal, sha
from run_exp_subjects import BASE, CONTROL, PUBLISH, SUBJECTS, status, completed_control


def main():
    PUBLISH.mkdir(parents=True, exist_ok=True)
    portable = ROOT / "data/processed/exp_subjects_A01"
    portable.mkdir(parents=True, exist_ok=True)
    rows = status()
    gate = json.loads((BASE / "quality_gate.json").read_text())
    maps = {}
    for subject in SUBJECTS:
        for suffix in (".png", ".json"):
            shutil.copy2(
                BASE / "map" / (subject + "_map" + suffix),
                PUBLISH / (subject + "_map" + suffix),
            )
        shutil.copy2(
            BASE / "map" / (subject + "_map.mat"), portable / (subject + "_map.mat")
        )
        report = json.loads((BASE / "map" / (subject + "_map.json")).read_text())
        maps[subject] = dict(
            trainingConditions=report["trainingConditions"],
            testConditions=report["testConditions"],
            trainingLoadsKg=report["trainingLoadsKg"],
            testLoadsKg=report["testLoadsKg"],
            testRMSE=float(
                np.mean([r["rmse"] for r in report["metrics"] if r["split"] != "train"])
            ),
            modelSHA256=sha(BASE / "map" / (subject + "_map.mat")),
        )
    for name in (
        "quality_gate.json",
        "cycle_decisions.csv",
        "trial_clocks.csv",
        "quality_overview.png",
        "data_quality_audit_v2.pdf",
    ):
        shutil.copy2(BASE / name, PUBLISH / name)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for i, subject in enumerate(SUBJECTS):
        info = next(r for r in gate["summary"] if r["subject"] == subject)
        axes[0].bar(i, info["acceptedCycles"], color="#23844a")
        axes[0].bar(
            i, info["rejectedCycles"], bottom=info["acceptedCycles"], color="#d9a17b"
        )
        axes[0].text(
            i,
            info["acceptedCycles"] + 0.3,
            str(info["acceptedCycles"]) + " qualified",
            ha="center",
            fontsize=9,
        )
        axes[1].bar(
            i - 0.15,
            maps[subject]["trainingConditions"],
            width=0.3,
            color="#1469b0",
            label="Training" if i == 0 else None,
        )
        axes[1].bar(
            i + 0.15,
            maps[subject]["testConditions"],
            width=0.3,
            color="#c8841c",
            label="Held out" if i == 0 else None,
        )
    for ax in axes:
        ax.set_xticks(range(2), [s.replace("240", "") for s in SUBJECTS])
        ax.grid(axis="y", alpha=0.2)
    axes[0].set(
        title="BSE: complete, aligned and smooth action records",
        ylabel="Actions (green qualified, tan excluded)",
    )
    axes[1].set(title="MAP: complete sensor M bouts", ylabel="Load/cadence conditions")
    axes[1].legend()
    fig.suptitle("617 / 627 already complete: skipped | New A=0.1 Exp datasets")
    fig.tight_layout()
    fig.savefig(PUBLISH / "completion_overview.png", dpi=160)
    plt.close(fig)

    # Show a complete sensor M bout whose companion force export stops early,
    # and a complete force export with weak q/torque phase agreement.
    review = json.loads((ROOT / "docs/exp_one_cycle_review/review.json").read_text())
    fig, axes = plt.subplots(2, 2, figsize=(12, 6.8))
    for i, speed in enumerate((25, 45)):
        subject = "240725Exp_A01"
        trial = f"CalibrationData/dynamics_stoop_5kg_{speed}"
        row = next(
            r
            for r in review["conditions"]
            if r["subject"] == "240725Exp"
            and r["loadKg"] == 5
            and r["speedCode"] == speed
        )
        sample = loadmat(
            DATA / subject / "selected_windows" / f"5kg_{speed}.mat",
            simplify_cells=True,
        )
        origin = row["originEpochSeconds"]
        t = sample["timeEpochSeconds"] - origin
        axes[i, 0].plot(t, sample["activation"][:, :6].mean(1), color="#1469b0")
        s = processed_signal(DATA / subject / trial)
        cutoff = s["absoluteStart"] + s["t"][-1] - origin
        if cutoff < t[-1]:
            axes[i, 0].axvspan(
                cutoff,
                t[-1],
                color="#cc3434",
                alpha=0.15,
                label="No companion CEINMS export",
            )
            axes[i, 0].axvline(cutoff, color="#cc3434", ls="--")
            axes[i, 0].legend(fontsize=8)
        axes[i, 0].set(
            title=f"5 kg / {speed}: complete observed sensor bout",
            ylabel="A=0.1 mean back activation",
            xlabel="Time since EMG start (s)",
        )
        axes[i, 1].plot(
            s["t"], s["q"], color="#1469b0", label="Original scaled IK theta"
        )
        twin = axes[i, 1].twinx()
        twin.plot(
            s["t"], s["torque"], color="#c8841c", label="A01 CEINMS torque estimate"
        )
        twin.set_ylabel("CEINMS estimate (N m)", color="#c8841c")
        axes[i, 1].set(
            title=f"Native overlap only: q/torque corr = {np.corrcoef(s['q'],s['torque'])[0,1]:.2f}",
            ylabel="Theta (rad)",
            xlabel="Time since common native start (s)",
        )
        for ax in axes[i]:
            ax.grid(alpha=0.2)
    fig.suptitle("Why 240725Exp_A01 cannot yet pass the BSE training/holdout gate")
    fig.tight_layout()
    fig.savefig(PUBLISH / "240725Exp_A01_exclusion_evidence.png", dpi=160)
    plt.close(fig)
    controls = {}
    for subject in SUBJECTS:
        if not completed_control(CONTROL / subject):
            continue
        meta = json.loads((CONTROL / subject / "mpc_evaluation.json").read_text())
        validation = json.loads((CONTROL / subject / "validation.json").read_text())
        cases = [r for r in meta["metrics"] if not r["displayOnly"]]
        prep = json.loads((BASE / subject / "preparation.json").read_text())
        controls[subject] = dict(
            trainingActions=prep["trainingCycles"],
            testActions=prep["testCycles"],
            validatedLoadsKg=sorted({r["loadKg"] for r in cases}),
            cases=len(cases),
            maximumExcessNm=max(r["excessPeakNm"] for r in cases),
            excessOver1NmSeconds=sum(r["excessOver1NmSeconds"] for r in cases),
            assistancePeakNm=max(r["assistancePeakNm"] for r in cases),
            abdomenBackMeanRatio=float(np.mean([r["abdomenBackRatio"] for r in cases])),
            validation=validation,
        )
        print(subject, json.dumps(controls[subject], indent=2))
    result = dict(
        subjects=rows,
        MAP=maps,
        quality=gate["summary"],
        NMPC=controls,
        originalDataModified=False,
        referenceTorqueIndependentlyCalibrated=False,
        referenceTorquePolicy="A01 CEINMS original-equivalent active-only estimate; not measured torque",
        referencesRescaledAgain=False,
        hardwareClosedLoopValidated=False,
    )
    (PUBLISH / "completion_status.json").write_text(json.dumps(result, indent=2) + "\n")
    manifest = dict(
        producer="experiments/multi_subject/run_exp_subjects.py",
        files={
            p.name: sha(p)
            for p in sorted(PUBLISH.iterdir())
            if p.is_file() and p.name != "export_metadata.json"
        },
    )
    (PUBLISH / "export_metadata.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
