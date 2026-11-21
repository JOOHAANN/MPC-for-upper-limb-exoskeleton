"""Finish only missing BSE/NMPC stages for 5/10/15 kg; keep MAP independent."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from audit_subjects import ROOT
from run_exp_subjects import completed_model, completed_control, run, sha
from prepare_exp_three_loads import BASE, SUBJECT

CONTROL = ROOT / "results/exp_three_loads_nmpc"


def protected_files():
    roots = [
        ROOT / "results/mpc_threshold_tracking" / s
        for s in ("240617Pilot", "240627Pilot")
    ]
    roots += [ROOT / "results/exp_subjects_A01_nmpc/240730Exp_A01"]
    files = [p / name for p in roots for name in ("neural_models.mat", "hill_mpc.mat")]
    files += list((ROOT / "results/exp_subjects_A01/map").glob("*_map.mat"))
    return {str(p): sha(p) for p in files}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-nfev", type=int, default=35)
    args = parser.parse_args()
    before = protected_files()
    BASE.mkdir(parents=True, exist_ok=True)
    run(
        "train_exp_map.py",
        "--output",
        ROOT / "results/exp_subjects_A01/map",
        result_root=BASE,
    )
    folder = BASE / SUBJECT
    if not (folder / "prepared.mat").exists():
        run("prepare_exp_three_loads.py", result_root=BASE)
    if not completed_model(folder):
        command = [
            SUBJECT,
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
        run("fit_subject.py", *command, result_root=BASE)
    else:
        print("SKIP fitted BSE", SUBJECT, flush=True)
    if not completed_control(CONTROL / SUBJECT):
        if not (folder / "mpc_selected_training_options.json").exists():
            run(
                "control_subject.py",
                SUBJECT,
                "--constraint-controller",
                "--angle-governor",
                1.75,
                "--training-only",
                result_root=BASE,
            )
        options = CONTROL / SUBJECT / "selected_threshold_options.json"
        if not options.exists():
            run(
                "refine_threshold.py",
                SUBJECT,
                "--base",
                BASE,
                "--destination",
                CONTROL,
                "--scales",
                0.25,
                0.5,
                "--resume",
                result_root=BASE,
            )
        chosen = json.loads(options.read_text())
        if max(r["excessPeakNm"] for r in chosen["cases"]) > 0.1:
            run(
                "refine_threshold.py",
                SUBJECT,
                "--base",
                BASE,
                "--destination",
                CONTROL,
                "--scales",
                0.625,
                0.75,
                1.0,
                "--resume",
                result_root=BASE,
            )
        run(
            "control_subject.py",
            SUBJECT,
            "--constraint-controller",
            "--frozen-options",
            options,
            "--resume",
            result_root=CONTROL,
        )
        run("validate_subject.py", SUBJECT, result_root=CONTROL)
    else:
        print("SKIP fitted NMPC", SUBJECT, flush=True)
    assert before == protected_files(), "An already-completed subject/MAP changed"
    (BASE / "preserved_models.json").write_text(
        json.dumps(dict(allUnchanged=True, files=before), indent=2) + "\n"
    )
    run("report_exp_three_loads.py", result_root=CONTROL)


if __name__ == "__main__":
    main()
