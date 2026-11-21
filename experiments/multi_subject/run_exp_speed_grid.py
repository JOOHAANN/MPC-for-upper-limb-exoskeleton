"""Fit 725 BSE cases with torque support; retain completed MAP/MPC models."""

import argparse
import json
import shutil

from audit_subjects import ROOT
from prepare_exp_speed_grid import BASE, SUBJECT
from run_exp_three_loads import protected_files
from run_exp_subjects import run, sha, completed_model


def write_status():
    """Separate completion of fitting from meeting every trend indicator."""
    docs = ROOT / "docs/exp_speed_grid"
    folder = BASE / SUBJECT
    report = json.loads((docs / "evaluation.json").read_text())
    assert report["modelSHA256"] == sha(folder / "neural_models.mat")
    assert report["modelSHA256"] == sha(
        ROOT / "data/processed/exp_speed_grid" / (SUBJECT + "_models.mat")
    )
    issues = [
        dict(
            source=r["source"],
            loadKg=r["loadKg"],
            speedCode=r["speedCode"],
            shapeCorrelation=r["shapeCorrelation"],
            predictedValleys=r["predictedValleys"],
        )
        for r in report["metrics"]
        if r["shapeCorrelation"] < 0.7 or r["predictedValleys"] != 2
    ]
    status = dict(
        subject=SUBJECT,
        requestedConditions=15,
        retainedConditions=10,
        excludedConditions=5,
        retainedPriorityConditions=6,
        BSECalibrationCompleted=True,
        allRetainedConditionsAreCalibration=True,
        independentTestCompleted=False,
        MAPRefitted=False,
        MPCRefittedForThisModel=False,
        allConditionsMeetShapeIndicator=not issues,
        remainingIssues=issues,
        numericalValidationPassed=report["numericalValidationPassed"],
        modelSHA256=report["modelSHA256"],
    )
    (docs / "completion_status.json").write_text(json.dumps(status, indent=2) + "\n")
    for name in ("model_selection.json", "candidate_comparison.json"):
        if (BASE / name).exists():
            shutil.copy2(BASE / name, docs / name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fit", action="store_true", help="Train only if this grid model is missing."
    )
    parser.add_argument(
        "--refine",
        action="store_true",
        help="Continue the same prepared-grid checkpoint explicitly.",
    )
    parser.add_argument("--sweeps", type=int, default=2)
    parser.add_argument("--max-nfev", type=int, default=20)
    args = parser.parse_args()
    before = protected_files()
    old_model = ROOT / "results/exp_three_loads" / SUBJECT / "neural_models.mat"
    before[str(old_model)] = sha(old_model)
    folder = BASE / SUBJECT
    run("prepare_exp_speed_grid.py", result_root=BASE)
    completed = completed_model(folder)
    if (args.fit and not completed) or args.refine:
        policy_path = folder / "fit_policy.json"
        profile = json.loads(policy_path.read_text()) if policy_path.exists() else {}
        command = [
            SUBJECT,
            "--initial-model",
            old_model,
            "--continuous-training",
            "--robust-recovery",
            "--steps",
            6,
            "--sweeps",
            args.sweeps,
            "--max-nfev",
            args.max_nfev,
            "--minimum-feedback-blend",
            profile.get("minimumFeedbackBlend", 0.05),
            "--maximum-emg-lead",
            profile.get("maximumCausalEMGLeadGainSeconds", 2.5),
            "--maximum-flexion-request",
            profile.get("maximumFlexionRequestNm", 80.0),
            "--neural-warmup-nfev",
            profile.get("neuralWarmupMaxNfev", 12),
            "--shape-tail-weight",
            profile.get("shapeTailWeight", 3.0),
            "--angle-value-weight",
            profile.get("angleValueWeight", 0.25 if profile else 0.6),
            "--valley-depth-weight",
            profile.get("valleyDepthWeight", 1.5 if profile else 3.0),
        ]
        if (folder / "checkpoint.mat").exists():
            command.append("--resume")
        run("fit_subject.py", *command, result_root=BASE)
    elif not completed:
        raise RuntimeError("Grid model not fitted. Run with --fit first.")
    else:
        print("SKIP completed 725 speed-grid BSE; no MAP/MPC refit", flush=True)
    run("report_exp_speed_grid.py", result_root=BASE)
    after = protected_files()
    after[str(old_model)] = sha(old_model)
    assert before == after, "A previously completed subject/MAP/MPC model changed"
    (BASE / "preserved_models.json").write_text(
        json.dumps(dict(allUnchanged=True, files=before), indent=2) + "\n"
    )
    report = json.loads((folder / "evaluation.json").read_text())
    assert report[
        "numericalValidationPassed"
    ], "Fine-step replay/reference-independence failed"
    portable = ROOT / "data/processed/exp_speed_grid"
    portable.mkdir(parents=True, exist_ok=True)
    shutil.copy2(folder / "neural_models.mat", portable / (SUBJECT + "_models.mat"))
    for name in ("fit_policy.json", "preparation.json"):
        shutil.copy2(folder / name, ROOT / "docs/exp_speed_grid" / name)
    shutil.copy2(
        BASE / "quality_gate.json", ROOT / "docs/exp_speed_grid/quality_gate.json"
    )
    shutil.copy2(
        BASE / "preserved_models.json",
        ROOT / "docs/exp_speed_grid/preserved_models.json",
    )
    for name in ("reference_grid.png", "activation_timing_grid.png"):
        shutil.copy2(BASE / name, ROOT / "docs/exp_speed_grid" / name)
    write_status()


if __name__ == "__main__":
    main()
