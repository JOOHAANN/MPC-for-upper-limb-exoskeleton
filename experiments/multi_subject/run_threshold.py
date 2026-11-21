"""Replay or retune near-threshold MPC on already qualified, frozen BSE models."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.io import loadmat

from audit_subjects import ROOT
from refine_threshold import BASE, DEST, setup

PUBLISH = ROOT / "docs/mpc_threshold_tracking"


def parameter_fingerprint(folder):
    model = loadmat(folder / "neural_models.mat", simplify_cells=True)
    values = np.r_[
        model["parameters"].ravel(),
        model["neuralParameters"].ravel(),
        model["controller"].ravel(),
    ]
    return hashlib.sha256(np.asarray(values, dtype="<f8").tobytes()).hexdigest()


def run(script, *arguments):
    env = os.environ.copy()
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["OMP_NUM_THREADS"] = "1"
    env["EXO_SUBJECT_RESULTS_ROOT"] = str(DEST)
    subprocess.run(
        [sys.executable, str(Path(__file__).with_name(script)), *map(str, arguments)],
        check=True,
        env=env,
        creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tune",
        action="store_true",
        help="Select controller options on training actions; never refit BSE or MAP.",
    )
    parser.add_argument("--subjects", nargs="+", default=["240617Pilot", "240627Pilot"])
    args = parser.parse_args()

    def subject_job(subject):
        if not (BASE / subject / "prepared.mat").exists():
            raise RuntimeError(
                "Qualified subject cache is missing. Run run_subjects.py --retrain first: "
                + subject
            )
        folder = setup(subject)
        if args.tune:
            run("refine_threshold.py", subject, "--scales", 0.25, 0.5)
            frozen = folder / "selected_threshold_options.json"
        else:
            frozen = PUBLISH / (subject + "_options.json")
            saved = json.loads(frozen.read_text())
            if saved["modelParameterSHA256"] != parameter_fingerprint(folder):
                raise RuntimeError(
                    "BSE parameters changed; rerun with --tune for " + subject
                )
        run(
            "control_subject.py",
            subject,
            "--constraint-controller",
            "--frozen-options",
            frozen,
            "--resume",
        )
        run("validate_subject.py", subject)

    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(subject_job, args.subjects))
    run(
        "report_subjects.py",
        "--partial",
        "--publish-directory",
        PUBLISH,
        "--model-directory",
        ROOT / "data/processed/subject_quality_v2",
    )
    run("publish_threshold.py")


if __name__ == "__main__":
    main()
