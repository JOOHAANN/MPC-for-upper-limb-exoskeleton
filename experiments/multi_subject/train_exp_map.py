"""Fit the established MAP architecture on reviewed, complete A=0.1 bouts.

MAP can use a complete sensor bout without a complete CEINMS force export.
BSE/NMPC still require separately qualified joint, force and activation data.
Whole BSE holdout trials are excluded from MAP training to prevent leakage.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat, savemat

from audit_subjects import DATA, OUT, ROOT
from prepare_subjects import cycles, fit_map, map_predict, qualified_splits, structs
from quality_gate import smoothness, sha

SUBJECTS = ("240725Exp_A01", "240730Exp_A01")
HOLDOUT = {"240730Exp_A01": {"dynamics_stoop_12.5kg_45"}}


def native_splits(subject, gate):
    accepted = {
        r["id"] for r in gate["cycles"] if r["subject"] == subject and r["accepted"]
    }
    records = []
    for row in gate["trials"]:
        if row["subject"] != subject or row["gate"] != "pass":
            continue
        cache = OUT / "signals" / row["signalCache"]
        if sha(cache) != row["signalCacheSHA256"]:
            raise ValueError("Signal cache changed after quality review")
        records.extend(
            r for r in cycles(row, dict(np.load(cache))) if r["id"] in accepted
        )
    records = qualified_splits(subject, records)
    for r in records:
        if Path(r["trialPath"]).name in HOLDOUT.get(subject, set()):
            r["split"] = "test_file"
    return records


def fit_subject_map(subject, gate, review, destination):
    native = (
        native_splits(subject, gate)
        if gate["version"] != "excitation-only-map-v1"
        else []
    )
    native_train = [r for r in native if r["split"] == "train"]
    held_paths = {r["trialPath"] for r in native if r["split"] != "train"}
    # Freeze the back-muscle lumping on qualified training geometry only.
    weights = (
        np.mean(np.vstack([r["weights"] for r in native_train]), axis=0)
        if native_train
        else np.ones(6)
    )
    weights /= weights.sum()
    records, decisions = [], []
    for row in review["conditions"]:
        if row["subject"] + "_A01" != subject or not row.get("selected"):
            continue
        filename = f"{row['loadKg']:g}kg_{row['speedCode']}.mat"
        path = Path(
            row.get("windowPath", DATA / subject / "selected_windows" / filename)
        )
        item = loadmat(path, simplify_cells=True)
        six = item["activation"][:, :6]
        back = six @ weights
        diagnostics = smoothness(back)
        reasons = [
            "smoothness:" + k
            for k, v in gate["smoothnessUpperLimits"].items()
            if diagnostics[k] > v
        ]
        if item["activation"].min() < -1e-6 or item["activation"].max() > 1 + 1e-6:
            reasons.append("activation_outside_zero_one")
        decisions.append(
            dict(
                loadKg=row["loadKg"],
                speedCode=row["speedCode"],
                trialPath=row["trial"],
                accepted=not reasons,
                reasons=reasons,
                inputSHA256=sha(path),
                inputPath=str(path),
                smoothness=diagnostics,
            )
        )
        if reasons:
            continue
        t = item["timeEpochSeconds"] - item["timeEpochSeconds"][0]
        duration = float(t[-1])
        profile = np.column_stack(
            [np.interp(np.linspace(0, duration, 200), t, six[:, j]) for j in range(6)]
        )
        records.append(
            dict(
                id=filename[:-4],
                trialPath=row["trial"],
                loadKg=row["loadKg"],
                speedCode=row["speedCode"],
                cadenceBPM=60 / duration,
                duration=duration,
                q=np.zeros(len(t)),
                profile=profile,
                weights=np.tile(weights, (len(t), 1)),
                six=six,
                back=back,
                split="test_file" if row["trial"] in held_paths else "train",
            )
        )
    # Designate an additional whole-condition holdout per load where possible.
    # A lone usable load condition remains training-only and is reported so.
    for load in sorted({r["loadKg"] for r in records}):
        group = sorted(
            [r for r in records if r["loadKg"] == load], key=lambda r: r["speedCode"]
        )
        if len(group) >= 2 and not any(r["split"] != "train" for r in group):
            group[-1]["split"] = "test_condition"
    train = [r for r in records if r["split"] == "train"]
    if len({r["loadKg"] for r in train}) != 5:
        raise RuntimeError(subject + ": MAP training no longer covers all five loads")
    model = fit_map(train)
    metrics = []
    for r in records:
        r["prediction"] = map_predict(model, r)
        metrics.append(
            dict(
                id=r["id"],
                trialPath=r["trialPath"],
                loadKg=r["loadKg"],
                speedCode=r["speedCode"],
                split=r["split"],
                rmse=float(np.sqrt(np.mean((r["prediction"] - r["back"]) ** 2))),
                correlation=float(np.corrcoef(r["prediction"], r["back"])[0, 1]),
                outOfDomain=r["mapOutOfDomain"],
            )
        )
    destination.mkdir(parents=True, exist_ok=True)
    savemat(
        destination / (subject + "_map.mat"),
        dict(
            mapModel=model,
            trainingTrialPaths=np.array(
                sorted({r["trialPath"] for r in train}), dtype=object
            ),
            records=structs(
                [
                    {k: v for k, v in r.items() if k not in ("profile", "weights")}
                    for r in records
                ]
            ),
        ),
        long_field_names=True,
        do_compression=True,
    )
    report = dict(
        subject=subject,
        architecture=model["architecture"],
        qualityGateSHA256=(
            sha(OUT / "quality_gate.json")
            if gate["version"] != "excitation-only-map-v1"
            else ""
        ),
        trainingDataPolicy=gate["version"],
        excitationAuditPath=gate.get("excitationAuditPath", ""),
        excitationAuditSHA256=gate.get("excitationAuditSHA256", ""),
        shapeFactor=0.1,
        trainingConditions=len(train),
        testConditions=len(records) - len(train),
        trainingLoadsKg=sorted({r["loadKg"] for r in train}),
        testLoadsKg=sorted({r["loadKg"] for r in records if r["split"] != "train"}),
        bseHeldOutTrialPaths=sorted(held_paths),
        MAPTrainingTrialPaths=sorted({r["trialPath"] for r in train}),
        usesReferenceAngleOrTorque=False,
        decisions=decisions,
        metrics=metrics,
    )
    (destination / (subject + "_map.json")).write_text(
        json.dumps(report, indent=2) + "\n"
    )
    fig, axes = plt.subplots(5, 5, figsize=(15, 11), sharex=True, sharey=True)
    for i, load in enumerate((5, 7.5, 10, 12.5, 15)):
        for j, speed in enumerate((25, 35, 45, 55, 65)):
            ax = axes[i, j]
            r = next(
                (r for r in records if r["loadKg"] == load and r["speedCode"] == speed),
                None,
            )
            if r:
                phase = np.linspace(0, 1, len(r["back"]))
                ax.plot(
                    phase, r["back"], color=".55", lw=1, label="Experimental activation"
                )
                ax.plot(phase, r["prediction"], color="#1469b0", lw=1, label="MAP")
                ax.set_title(
                    f"{load:g} kg / {speed} | {'train' if r['split']=='train' else 'held out'}",
                    fontsize=8,
                )
            else:
                ax.set_facecolor(".95")
                ax.text(
                    0.5,
                    0.5,
                    "Excluded / unavailable",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                    fontsize=7,
                )
                ax.set_title(f"{load:g} kg / {speed}", fontsize=8)
            ax.grid(alpha=0.2)
    handles, labels = next(
        ax.get_legend_handles_labels() for ax in axes.ravel() if ax.lines
    )
    fig.legend(handles, labels, loc="lower center", ncol=2)
    fig.suptitle(
        subject
        + " | five-load MAP evaluation | one complete observed M bout per condition"
    )
    fig.supxlabel("Observed action phase (0-1)", y=0.035)
    fig.supylabel("Activation (dimensionless)")
    fig.tight_layout(rect=(0, 0.04, 1, 0.96))
    fig.savefig(destination / (subject + "_map.png"), dpi=150)
    plt.close(fig)
    print(
        subject,
        "MAP",
        len(train),
        "train",
        len(records) - len(train),
        "held out; RMSE",
        [
            (
                s,
                float(
                    np.mean(
                        [r["rmse"] for r in metrics if (r["split"] == "train") == s]
                    )
                ),
            )
            for s in (True, False)
        ],
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subjects", nargs="+", default=list(SUBJECTS))
    parser.add_argument("--output", type=Path, default=OUT / "map")
    args = parser.parse_args()
    for subject in args.subjects:
        if (args.output / (subject + "_map.mat")).exists():
            report = json.loads((args.output / (subject + "_map.json")).read_text())
            # BSE admission can change without changing an already-fitted MAP.
            # Validate MAP's own labels, not an unrelated IK/force quality gate.
            if (
                report.get("excitationAuditSHA256")
                and sha(Path(report["excitationAuditPath"]))
                != report["excitationAuditSHA256"]
            ):
                raise RuntimeError("MAP excitation audit changed")
            for decision in report["decisions"]:
                path = Path(
                    decision.get(
                        "inputPath",
                        (
                            DATA
                            / subject
                            / "selected_windows"
                            / f"{decision['loadKg']:g}kg_{decision['speedCode']}.mat"
                        ),
                    )
                )
                if sha(path) != decision["inputSHA256"]:
                    raise RuntimeError("MAP inputs changed; review before retraining")
            print("SKIP already fitted MAP", subject, flush=True)
            continue
        from excitation_only_map import build_review

        gate, review = build_review(subject, args.output / "excitation_review")
        fit_subject_map(subject, gate, review, args.output)


if __name__ == "__main__":
    main()
