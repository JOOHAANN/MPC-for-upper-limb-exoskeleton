"""Complete only missing subject-specific MAP/BSE/NMPC stages.

The successful 617/627 runs are verified and skipped. A derived Exp dataset
must pass the same strict clock/smoothing gate before its BSE can be fitted.
Existing completed stages are never silently recalibrated.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

from audit_subjects import ROOT, DATA

BASE = ROOT / "results/exp_subjects_A01"
CONTROL = ROOT / "results/exp_subjects_A01_nmpc"
PUBLISH = ROOT / "docs/exp_subjects_A01"
SUBJECTS = ("240725Exp_A01", "240730Exp_A01")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(script, *args, result_root=BASE):
    env = {
        **os.environ,
        "EXO_SUBJECT_RESULTS_ROOT": str(result_root),
        "OPENBLAS_NUM_THREADS": "1",
        "OMP_NUM_THREADS": "1",
    }
    subprocess.run(
        [sys.executable, str(Path(__file__).with_name(script)), *map(str, args)],
        cwd=ROOT,
        env=env,
        check=True,
    )


def completed_model(folder):
    needed = (
        "prepared.mat",
        "neural_models.mat",
        "hill_sequence.mat",
        "fit_policy.json",
        "bse_metrics.json",
    )
    if not all((folder / name).exists() for name in needed):
        return False
    policy = json.loads((folder / "fit_policy.json").read_text())
    if policy["preparedSHA256"] != sha(folder / "prepared.mat") or policy[
        "modelSHA256"
    ] != sha(folder / "neural_models.mat"):
        raise ValueError("Completed model fingerprint changed: " + str(folder))
    return True


def completed_control(folder):
    needed = ("hill_mpc.mat", "mpc_evaluation.json", "validation.json")
    if not all((folder / name).exists() for name in needed):
        return False
    meta = json.loads((folder / "mpc_evaluation.json").read_text())
    valid = json.loads((folder / "validation.json").read_text())
    if (
        meta["modelSHA256"] != sha(folder / "neural_models.mat")
        or valid["modelSHA256"] != meta["modelSHA256"]
    ):
        raise ValueError("Controller/model fingerprint mismatch")
    return True


def status():
    rows = []
    for subject in ("240617Pilot", "240627Pilot"):
        folder = ROOT / "results/mpc_threshold_tracking" / subject
        ready = completed_model(folder) and completed_control(folder)
        rows.append(
            dict(
                subject=subject,
                status=(
                    "already_complete_skip"
                    if ready
                    else "existing_results_need_inspection"
                ),
                MAP=ready,
                BSE=ready,
                NMPC=ready,
                modelSHA256=(
                    sha(folder / "neural_models.mat")
                    if (folder / "neural_models.mat").exists()
                    else None
                ),
            )
        )
    for subject in SUBJECTS:
        personalized_base = BASE
        personalized_control = CONTROL
        if (ROOT / "results/exp_three_loads" / subject / "prepared.mat").exists():
            personalized_base = ROOT / "results/exp_three_loads"
            personalized_control = ROOT / "results/exp_three_loads_nmpc"
        state = (
            "complete"
            if completed_model(personalized_base / subject)
            and completed_control(personalized_control / subject)
            else "pending"
        )
        if personalized_base == BASE and (BASE / "quality_gate.json").exists():
            summary = json.loads((BASE / "quality_gate.json").read_text())["summary"]
            selected = next((r for r in summary if r["subject"] == subject), None)
            if (
                selected
                and selected["acceptedCycles"] < 5
                and not completed_model(BASE / subject)
            ):
                state = "insufficient_qualified_BSE_data"
        rows.append(
            dict(
                subject=subject,
                status=state,
                MAP=(BASE / "map" / (subject + "_map.mat")).exists(),
                BSE=completed_model(personalized_base / subject),
                NMPC=completed_control(personalized_control / subject),
                modelDirectory=str(personalized_base / subject),
                controlDirectory=str(personalized_control / subject),
                requestedBSELoadsKg=[5, 10, 15],
            )
        )
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--status", action="store_true")
    parser.add_argument("--max-nfev", type=int, default=35)
    parser.add_argument("--legacy-review", action="store_true")
    args = parser.parse_args()
    if args.status:
        print(json.dumps(status(), indent=2))
        return
    if not args.legacy_review:
        run("run_exp_three_loads.py", "--max-nfev", args.max_nfev)
        return
    if not all(r["status"] == "already_complete_skip" for r in status()[:2]):
        raise RuntimeError(
            "Inspect original subjects; this runner will not overwrite them"
        )
    before = {r["subject"]: r["modelSHA256"] for r in status()[:2]}
    BASE.mkdir(parents=True, exist_ok=True)
    if not (BASE / "data_quality.json").exists():
        run("audit_subjects.py", "--subjects", *SUBJECTS, "--include-derived")
    gate_file = BASE / "quality_gate.json"
    if not gate_file.exists():
        run(
            "quality_gate.py",
            "--inventory",
            BASE / "data_quality.json",
            "--output",
            BASE,
            "--subjects",
            *SUBJECTS,
            "--loads",
            5,
            7.5,
            10,
            12.5,
            15,
        )
    gate = json.loads(gate_file.read_text())
    for row in gate["trials"]:
        for name, digest in row["fileSha256"].items():
            if sha(DATA / row["subject"] / row["trial"] / name) != digest:
                raise RuntimeError(
                    "Data changed after audit; review before new training"
                )
    run("train_exp_map.py")
    for subject in SUBJECTS:
        folder = BASE / subject
        if not (folder / "prepared.mat").exists():
            run(
                "prepare_subjects.py",
                "--quality-gate",
                gate_file,
                "--allow-training-changes",
                "--subjects",
                subject,
                "--loads",
                5,
                7.5,
                10,
                12.5,
                15,
                "--held-out-trials",
                "dynamics_stoop_12.5kg_45",
                "--map-model-directory",
                BASE / "map",
            )
        if not (folder / "prepared.mat").exists():
            print(
                "SKIP BSE/NMPC: insufficient qualified training and holdout data",
                subject,
                flush=True,
            )
            continue
        if not completed_model(folder):
            command = [
                subject,
                "--sweeps",
                2,
                "--max-nfev",
                args.max_nfev,
                "--steps",
                10,
                "--continuous-training",
                "--robust-recovery",
                "--subject-body-search",
            ]
            if (folder / "checkpoint.mat").exists():
                command.append("--resume")
            run("fit_subject.py", *command)
        if completed_control(CONTROL / subject):
            print("SKIP already complete BSE/NMPC", subject, flush=True)
            continue
        if not (folder / "mpc_selected_training_options.json").exists():
            run(
                "control_subject.py",
                subject,
                "--constraint-controller",
                "--angle-governor",
                1.75,
                "--training-only",
            )
        options = CONTROL / subject / "selected_threshold_options.json"
        if not options.exists():
            run(
                "refine_threshold.py",
                subject,
                "--base",
                BASE,
                "--destination",
                CONTROL,
                "--scales",
                0.25,
                0.5,
                "--resume",
            )
        selection = json.loads(options.read_text())
        if max(row["excessPeakNm"] for row in selection["cases"]) > 0.1:
            # The decision uses training cases only, before holdout replay.
            run(
                "refine_threshold.py",
                subject,
                "--base",
                BASE,
                "--destination",
                CONTROL,
                "--scales",
                0.625,
                0.75,
                "--resume",
            )
        if not completed_control(CONTROL / subject):
            run(
                "control_subject.py",
                subject,
                "--constraint-controller",
                "--frozen-options",
                options,
                "--resume",
                result_root=CONTROL,
            )
            run("validate_subject.py", subject, result_root=CONTROL)
    after = {r["subject"]: r["modelSHA256"] for r in status()[:2]}
    assert before == after, "Successful original subject models changed"
    if any(completed_control(CONTROL / subject) for subject in SUBJECTS):
        run(
            "report_subjects.py",
            "--partial",
            "--publish-directory",
            PUBLISH,
            "--model-directory",
            ROOT / "data/processed/exp_subjects_A01",
            result_root=CONTROL,
        )
    run("report_exp_subjects.py", result_root=CONTROL)


if __name__ == "__main__":
    main()
