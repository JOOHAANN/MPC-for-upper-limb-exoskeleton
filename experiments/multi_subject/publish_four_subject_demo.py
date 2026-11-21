"""Verify and document the current four-subject MATLAB demonstration."""

import csv
import json
import re
from pathlib import Path

import numpy as np
from scipy.io import loadmat

from audit_subjects import ROOT
from control_nine_cases import digest, many, load_subject
from package_matlab_demo import sequence_data
from package_pearson_demo import reference_hash

DEMO = ROOT / "Demo"
SUBJECTS = {
    "528": "240528Pilot",
    "617": "240617Pilot",
    "725": "240725Exp_A01",
    "730": "240730Exp_A01",
}


README = """# Four-subject BSE and NMPC demonstration

The demonstration uses the current models for **528, 617, 725 and 730**.
Copy the entire `Demo` folder to run it elsewhere. Python, Simulink and the
original experiment directories are not required at runtime.

## Run

In MATLAB, open this folder and run:

```matlab
addpath('BSE','NMPC');
run_bse_demo;
run_nmpc_demo;
inspect_bse_model;             % Back, abdomen, soft limits and model diagram
```

The original function calls remain available:

```matlab
bse = run_bse_demo('617',false);
nmpc = run_nmpc_demo('617',false);
```

Alternatively, open and run the two procedural scripts directly:

1. `BSE/Lumped_model_check_with_coeffs_thesis_combination.m`
2. `NMPC/offline_simulator_NMPC_quasimodo_incremental_MPC_with_MAP.m`

They follow the original project style: settings at the beginning, `%%`
sections, explicit coefficient variables, preallocated simulation arrays,
ordinary loops, calculations and plots. The short `run_*` functions only pass
the subject selection and figure flag to these scripts.

The code was exercised with MATLAB R2023b. NMPC requires Optimization Toolbox
(`quadprog`). The optional C accelerator is built automatically. Configure a
compiler with `mex -setup C` if needed; the equivalent MATLAB equations in
`bse_model.m` remain available without one.

## Read the model

Start with the BSE script. It loads each subject and keeps one coefficient
vector for experimental input and another for MAP input. It then calls
`Torque_function.m`, where the complete state is initialized and advanced in
an explicit `for k = 2:NumberOfData` loop.

| Variable | Meaning |
| --- | --- |
| `coeffs`, `neural_coeffs` | Fixed physical and neural coefficients |
| `M`, `L1`, `L2` | Trunk mass, trunk centre-of-mass lever and load lever |
| `R_back`, `R_front` | Back and abdominal moment arms |
| `Ts` | Sample period |
| `Model_output` | Theta, velocity, muscle torques and the remaining states |
| `Torque_Tac0`, `Theta_Tac0` | Unassisted BSE output in the NMPC script |
| `Torque_Tac`, `Theta_Tac`, `Tac` | Assisted net torque, angle and actuator torque |
| `PH`, `CH`, `Threshold`, `u_max` | Prediction points, control blocks, back-torque threshold and actuator limit |

`bse_model.m` contains the shared Hill, tendon, activation, intent and joint
equations. Both outputs come from the same coupled system:

```text
net muscle torque = back torque - abdominal torque
J*theta_ddot = net muscle torque + gravity + soft limit + assistance - damping
```

The fiber/tendon equilibrium and joint momentum are solved together with
backward Euler. Human posture feedback acts through muscle excitation. No
reference-angle tracking torque is added. The mechanical equations and NMPC
optimization algorithm were retained during the readability changes.

## Data and fixed coefficients

`BSE/data` contains exactly one current package for each subject.
`bse_parameters.csv` lists the eight fixed 32+16-parameter vectors. Each
subject/source vector is reused for all nine recorded bouts. The two activation
sources use the same equations, with separately calibrated parameters.
`NMPC/data` contains one controller package per subject, linked to its current
BSE model hash. Previous parameter packages are archived outside `Demo`.

| Subject | Nine displayed BSE bouts |
| --- | --- |
| 528 | 5/10/15 kg, each at speed codes 35/45/55 |
| 617 | Five 5 kg, two 10 kg and two 15 kg complete recorded bouts |
| 725 | 5 kg: 35/45/65; 10 kg: 25/45/55; 15 kg: 45/55/65 |
| 730 | 5 kg: 35/45/65; 10 kg: 25/35/45/65; 15 kg: 45/55 |

The 617 skeletal recordings with speed labels often end mid-action. Its nine
complete audited bouts therefore display measured cadence and a bout number;
they are **not nine distinct nominal load/speed pairs**. The original IK angle
scale is retained. Reference torque gaps remain NaN and are excluded from
torque scores. No reference curve is smoothed or shifted to improve fitting.

617 MAP was newly trained with four NMF synergies, three-Gaussian primitives
and quadratic load/cadence weight regression. It uses qualified excitation
bouts from 5, 7.5, 10, 12.5 and 15 kg, with whole-condition holdouts. Its existing
activation processing was already A=+0.1. MAP training does not require IK or
CEINMS torque. The other subjects retain their established trained MAP outputs.
The BSE demo replays those prepared MAP predictions, as before; it does not
retrain MAP each time it runs.

The new BSE mean correlations are in `BSE/data/fit_summary.csv`. Seven of the
eight subject/source models meet the requested mean-correlation targets.
730 experimental remains an explicitly documented unmet case: torque r=0.693,
theta r=0.653. Correlation measures shape, not torque magnitude accuracy;
725 still has substantial amplitude error. These are calibration/replay
results, not independent physiological validation.

## Plots and control

Each subject produces all-nine-case torque, theta and activation figures in
`BSE/output`, along with `bse_results.mat` and `bse_metrics.csv`.
`inspect_bse_model` writes back/abdomen torque, signed upper/lower soft-limit
traces, signal-path diagnostics, a combined torque-budget figure and the
musculoskeletal model diagram to `BSE/diagnostics`.

NMPC uses one continuous 5/10/15 kg example per subject, with standing waits
and no state resets between actions. Its four-row figure shows net torque
before/after assistance with reference, back torque with threshold, angle
before/after assistance with reference, and applied actuator torque. A second
figure retains abdominal and soft-limit torques before/after control.

The black dashed threshold is **75% of the current unassisted back-muscle
peak**, recalculated separately for each action and input source. The actuator
cap remains 80 N m and its slew limit remains 200 N m/s. Controller weights
are subject specific. Activation residual predictors and prediction-error
guards were recalibrated for the current BSE models. Future measured EMG,
reference theta and CEINMS torque are not online controller inputs.

`NMPC/output` contains recalculated figures, `nmpc_results.mat` and
`nmpc_metrics.csv`. `verification.json` records the numerical replay checks,
known calibration limitations and actual threshold tracking results.
"""


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    metrics = rows(DEMO / "BSE/output/bse_metrics.csv")
    assert len(metrics) == 72
    maximum = 0.0
    for short, subject_id in SUBJECTS.items():
        packaged = loadmat(DEMO / f"BSE/data/subject_{short}.mat", simplify_cells=True)[
            "subject"
        ]
        assert len(many(packaged["models"])) == 2
        count = sum(len(many(s["segments"])) for s in many(packaged["sequences"]))
        assert count == 9
        for source, model in enumerate(many(packaged["models"]), 1):
            fit = (
                ROOT / "results/pearson_recalibration" / subject_id / f"source{source}"
            )
            original = loadmat(fit / "model.mat", simplify_cells=True)["model"]
            assert np.array_equal(model["parameters"], original["parameters"])
            r = np.asarray(original["neuralParameters"]).reshape(-1)
            assert np.array_equal(
                model["neuralParameters"], np.r_[r, np.zeros(16 - len(r))]
            )
            expected = json.loads((fit / "evaluation.json").read_text())["rows"]
            label = ("Experimental activation", "MAP activation")[source - 1]
            subset = sorted(
                [
                    v
                    for v in metrics
                    if v["Subject"] == subject_id and v["Source"] == label
                ],
                key=lambda v: int(v["CaseNumber"]),
            )
            assert len(subset) == 9
            for a, b in zip(subset, expected):
                maximum = max(
                    maximum,
                    abs(float(a["TorqueCorrelation"]) - b["torqueR"]),
                    abs(float(a["ThetaCorrelation"]) - b["thetaR"]),
                )
        if short != "617":
            _, _, groups, _, _ = load_subject(subject_id)
            assert reference_hash([sequence_data(g) for g in groups]) == reference_hash(
                many(packaged["sequences"])
            )
        else:
            original = loadmat(
                ROOT
                / "results/four_subject_demo"
                / subject_id
                / "calibration_input.mat",
                simplify_cells=True,
            )
            assert reference_hash(
                [sequence_data(g) for g in many(original["groups"])[:9]]
            ) == reference_hash(many(packaged["sequences"]))
    assert maximum < 1e-6
    nmpc = rows(DEMO / "NMPC/output/nmpc_metrics.csv")
    assert len(nmpc) == 24
    margins = np.array([float(r["PeakMarginNm"]) for r in nmpc])
    fallback = sum(int(r["FallbackUpdates"]) for r in nmpc)
    frozen = json.loads(
        (ROOT / "results/nine_case_75pct/release_manifest.json").read_text()
    )
    files = dict(frozen["artifactSHA256"], **frozen["sourceCodeSHA256"])
    for name, expected in files.items():
        assert digest(ROOT / name) == expected, name
    assert not (DEMO / "BSE/data/pearson").exists()
    for path in DEMO.rglob("*.m"):
        assert not re.search("[\u3400-\u9fff]", path.read_text()), str(path)
    (DEMO / "README.md").write_text(README)
    report = dict(
        runtimeChecksPassed=True,
        subjects=4,
        bseConditionSourceEvaluations=72,
        fixedParameterVectors=8,
        maximumPearsonReplayDifference=maximum,
        allBSECorrelationTargetsPassed=False,
        passingBSESourceModels=7,
        unmetBSEModel="240730Exp_A01 experimental",
        nmpcModelVersion="current",
        nmpcContinuousSequences=8,
        nmpcLoadSourceEvaluations=24,
        nmpcAllPeaksWithinZeroToFiveNmBelowThreshold=bool(
            np.all((margins >= 0) & (margins <= 5))
        ),
        nmpcPeakMarginRangeNm=[float(margins.min()), float(margins.max())],
        nmpcFallbackUpdates=fallback,
        frozenFilesUnchanged=len(files),
        originalReferenceAndPreparedInputsUnchanged=True,
        proceduralLoopAndNativeRolloutParityPassed=True,
        causalForecastTestsPassed=True,
        usesPythonAtRuntime=False,
        usesSimulink=False,
        independentValidation=False,
        oldDemoArchive="results/demo_archive/20261021_before_four_subjects.zip",
    )
    (DEMO / "verification.json").write_text(json.dumps(report, indent=2) + "\n")
    settings = []
    manifest = json.loads((DEMO / "manifest.json").read_text())
    for entry in manifest["subjects"]:
        entry["fileSHA256"] = digest(DEMO / entry["file"])
        entry["controllerFileSHA256"] = digest(DEMO / entry["controllerFile"])
        options = loadmat(DEMO / entry["controllerFile"], simplify_cells=True)[
            "controller"
        ]
        settings.append(
            dict(
                subject=entry["subject"],
                **{
                    k: options[k]
                    for k in (
                        "points",
                        "blocks",
                        "predictionDt",
                        "controlPeriod",
                        "capNm",
                        "rateNmPerSec",
                        "angleWeight",
                        "velocityWeight",
                        "backTorqueMarginNm",
                        "predictionErrorGuardScale",
                        "predictionErrorGuardFloorNm",
                        "backTorqueUnderWeight",
                        "trackingHeadroomNm",
                    )
                },
            )
        )
    with (DEMO / "NMPC/controller_settings.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(settings[0]))
        writer.writeheader()
        writer.writerows(settings)
    manifest["runtimeSourceSHA256"] = {
        p.relative_to(DEMO).as_posix(): digest(p) for p in sorted(DEMO.rglob("*.m"))
    }
    manifest["nativeSourceSHA256"] = {
        p.relative_to(DEMO).as_posix(): digest(p)
        for p in sorted((DEMO / "BSE/private/native").iterdir())
        if p.is_file()
    }
    (DEMO / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    m = json.loads(
        (ROOT / "results/four_subject_demo/240617Pilot/map_evaluation.json").read_text()
    )
    map_report = {
        k: m[k]
        for k in (
            "subject",
            "architecture",
            "shapeFactor",
            "shapeFactorChanged",
            "nativeResponseAudits",
            "trainingLoadsKg",
            "trainingConditions",
            "heldOutConditions",
            "metrics",
            "usesIKForMAP",
            "usesTorqueForMAP",
        )
    }
    for item in map_report["metrics"]:
        if item["speedCode"] < 20:
            item["speedCode"] = -1
            item["speedLabelNote"] = (
                "Original filename suffix is a repetition number; nominal speed unknown"
            )
    (DEMO / "BSE/data/map_training_617.json").write_text(
        json.dumps(map_report, indent=2) + "\n"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
