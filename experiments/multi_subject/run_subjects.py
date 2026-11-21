"""Reproduce frozen subject results or explicitly recalibrate the full pipeline."""

import argparse
import hashlib
import os
from pathlib import Path
import subprocess
import sys

from scipy.io import loadmat, savemat

from audit_subjects import ROOT, OUT

RESULTS = OUT


def run(script, *arguments):
    env = os.environ.copy()
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["OMP_NUM_THREADS"] = "1"
    env["EXO_SUBJECT_RESULTS_ROOT"] = str(RESULTS)
    subprocess.run(
        [sys.executable, str(Path(__file__).with_name(script)), *map(str, arguments)],
        check=True,
        env=env,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--retrain",
        action="store_true",
        help="Fit new MAP/BSE and tune MPC; default replays published fixed models/options.",
    )
    parser.add_argument(
        "--subjects",
        nargs="+",
        default=[
            "240617Pilot",
            "240627Pilot",
            "240719Pilot",
            "240723Pilot",
            "240725Pilot",
        ],
    )
    parser.add_argument("--max-nfev", type=int, default=35)
    parser.add_argument("--results-directory", type=Path)
    args = parser.parse_args()
    global RESULTS
    RESULTS = args.results_directory or (
        ROOT / "results/multi_subject_quality_v2" if args.retrain else OUT
    )
    RESULTS.mkdir(parents=True, exist_ok=True)
    run("audit_subjects.py", "--output", RESULTS)
    if args.retrain:
        run(
            "quality_gate.py",
            "--inventory",
            RESULTS / "data_quality.json",
            "--output",
            RESULTS,
        )
        run(
            "prepare_subjects.py",
            "--quality-gate",
            RESULTS / "quality_gate.json",
            "--allow-training-changes",
            "--subjects",
            *args.subjects,
        )
    else:
        run(
            "prepare_subjects.py",
            "--frozen-map-directory",
            ROOT / "data/processed/multi_subject",
        )
    for subject in args.subjects:
        path = RESULTS / subject
        if not (path / "prepared.mat").exists():
            print(subject, "SKIPPED: quality/coverage gate not satisfied", flush=True)
            continue
        if args.retrain:
            prepared = loadmat(path / "prepared.mat", simplify_cells=True)
            gate_hash = hashlib.sha256(
                (RESULTS / "quality_gate.json").read_bytes()
            ).hexdigest()
            if prepared.get("qualityGateSHA256") != gate_hash:
                print(
                    subject,
                    "SKIPPED: saved preparation is from an obsolete quality audit",
                    flush=True,
                )
                continue
            special = (
                ["--minimum-feedback-blend", 0.9, "--active-recovery-stress"]
                if subject == "240723Pilot"
                else []
            )
            run(
                "fit_subject.py",
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
                *special,
            )
            run(
                "control_subject.py",
                subject,
                "--angle-governor",
                1.65 if subject == "240723Pilot" else 1.75,
                "--constraint-controller",
            )
        else:
            portable = ROOT / "data/processed/multi_subject" / (subject + "_models.mat")
            model = loadmat(portable, simplify_cells=True)
            checkpoint = {
                key: model[key]
                for key in ("parameters", "neuralParameters", "controller")
            }
            path.mkdir(parents=True, exist_ok=True)
            savemat(path / "checkpoint.mat", checkpoint)
            special = (
                ["--steps", 10, "--continuous-training", "--robust-recovery"]
                if subject == "240725Pilot"
                else []
            )
            if subject == "240723Pilot":
                special = [
                    "--steps",
                    10,
                    "--continuous-training",
                    "--robust-recovery",
                    "--minimum-feedback-blend",
                    0.9,
                    "--active-recovery-stress",
                ]
            run("fit_subject.py", subject, "--sweeps", 0, "--resume", *special)
            options = (
                ROOT
                / "docs/multi_subject_evaluation"
                / (subject + "_mpc_evaluation.json")
            )
            run("control_subject.py", subject, "--frozen-options", options, "--resume")
        run("validate_subject.py", subject)
    if args.retrain:
        run(
            "report_subjects.py",
            "--partial",
            "--publish-directory",
            ROOT / "docs/subject_quality_v2",
            "--model-directory",
            ROOT / "data/processed/subject_quality_v2",
        )
        run("publish_quality.py")
    else:
        run("report_subjects.py")


if __name__ == "__main__":
    main()
