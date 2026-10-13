# Exoskeleton Musculoskeletal Modeling and NMPC

Current experiment: `demo_hill_compliant` couples two opposing compliant Hill
units to produce torque and angle. It compares experimental and MAP inputs,
whole-speed holdouts, the selected archived BSE, passive forces, and inferred
recruitment. All 32 parameters can be optimized, including shared body/load
geometry and elastic soft limits; the original angle reference is preserved.
[Run instructions, parameters, equations and figures](docs/HILL_COMPLIANT.zh-CN.md).
The rigid-tendon single-unit `demo_hill_forward` remains an earlier baseline.

The earlier `demo_bse_angle_priority` angle-equilibrium surrogate is historical;
its bounded angle and inverse-dynamics torque do not meet the current coupled
Hill requirement. [Historical method and limitations](docs/BSE_ANGLE_PRIORITY.zh-CN.md).

MATLAB research code for lumped muscle model calibration, muscle activation
profile (MAP) prediction, and nonlinear model predictive control (NMPC).

New experiment: run `demo_bse_generalization` for an archived BSE target audit
and nested leave-load/leave-speed-out torque estimation benchmark. See the
[method and limitations](experiments/bse_generalization/README.md) and
[Chinese optimization report](docs/BSE_OPTIMIZATION_REPORT.zh-CN.md).

The supplied project directory was named **Journal Paper of Upper-Limb
Exoskeleton**, but its project description and current equations primarily
describe a **back-support exoskeleton (BSE)** and **L5/S1** movement. This
repository preserves those models and their terminology.

## Quick start

Open MATLAB in this repository's root directory:

```matlab
setup_project;
workDir = prepare_workflow('map_predictor');
cd(workDir);
[weights, XP, MEP] = MAP_model(45, 10, 'SM', 'baseline_data.mat');
plot(MEP);
```

This uses the included fitted primitives and subject-specific baseline. The
legacy `elevation` argument is used as load in the adapted MAP model.
To reproduce the original comparison figures, run `test_MAP` in this workspace.

`prepare_workflow` copies the correct source versions and inputs to a fresh
directory under `results/`. Saved coefficients and simulation outputs stay
there. Make permanent edits in `src/`, then prepare a new workspace. Do not
recursively add `src` or `results` to the MATLAB path: experiments have
incompatible functions with the same name.

## Main workflows

| Workflow ID | Main entry script | Purpose |
| --- | --- | --- |
| `map_predictor` | `test_MAP` | Compare predicted and measured activation |
| `training_bse` | `Train_Lumped_muscle_model` | Calibrate the lumped muscle model |
| `training_bse` | `Lumped_model_check_with_coeffs_thesis_combination` | Evaluate archived model settings |
| `control_experimental` | `offline_simulator_NMPC_quasimodo_incremental_MPC` | NMPC with experimental activation |
| `control_map` | `offline_simulator_NMPC_quasimodo_incremental_MPC_with_MAP` | NMPC with predicted activation |
| `demo_calibration__geometry` | `Train_Lumped_muscle_model_Calibrate_Geo_model` | Geometric calibration example |
| `demo_calibration__muscle_tendon` | `Train_Lumped_muscle_model_EMG_and_LMT_to_MTF` | Muscle-tendon calibration example |

Use `list_workflows` for all 38 configurations, including Qualisys, damping,
scaling, and load/speed variants. [The workflow catalog](docs/WORKFLOWS.md)
maps every configuration to its original folder. Entries are preserved research
workflows, not a guarantee that every historical branch is complete.

For another experiment, return to the repository root and call
`prepare_workflow` with the desired ID, then `cd(workDir)` and run its script.
Transfer newly calibrated data between stages explicitly: a fresh preparation
always starts with archived input data.

## Repository layout

```text
config/workflows/   Experiment manifests selecting code and input files
src/shared/        Shared implementations with a single version
src/variants/      Reused implementations with experiment-specific versions
src/workflows/     Code unique to an experiment
data/processed/    Deduplicated MAT inputs and archived model coefficients
data/raw/          Local acquisition files (excluded from Git)
models/            Original Simulink models
tools/             Portable data paths and repository checks
tests/             MATLAB smoke and numerical comparison checks
docs/              Workflow catalog, migration records, known limitations
docs/reference/    Original papers, figures, and documents (local only)
results/           Generated experiment workspaces and outputs (local only)
```

Files with identical contents were consolidated; distinct implementations and
datasets were kept. [The inventory](docs/file_inventory.csv) records every
original file, destination, checksum, and disposition. The hashes describe
original files; MATLAB source was subsequently formatted and refactored.

## MATLAB requirements

Validation was performed with MATLAB R2023b. Depending on the workflow:

| Toolbox | Functions / use |
| --- | --- |
| Signal Processing Toolbox | `butter`, `filter`, `filtfilt` |
| Optimization Toolbox | `fmincon`, least-squares optimization |
| Global Optimization Toolbox | `ga` calibration |
| Parallel Computing Toolbox | Optimizers configured with `UseParallel=true` |
| Curve Fitting Toolbox | `fit`, `fittype`, primitive fitting |
| Statistics and Machine Learning Toolbox | `nnmf` synergy extraction |
| Deep Learning Toolbox | `mapminmax` in some legacy extraction variants |
| Simulink | Optional `.slx` models |

Not every workflow needs every toolbox. GA training can be expensive: some
preserved settings specify 7,000 generations and a population of 600. Select
the experiment and review its sampling, load segments, and optimizer settings
before starting a full training run.

## Data and reproducibility

Prepared MAT inputs are included. Raw recordings are retained locally but
excluded from Git. Rebuilding EMG or calibration inputs requires the original
`240528Pilot` dataset; see [data setup](data/raw/README.md). Some paths in the
historical scripts include an extra `CalibrationData` level; the path resolver
supports both layouts without renaming the recordings.

Model constants, sample ranges, and optimization bounds remain as supplied.
The geometry and muscle-tendon stages are not automatically chained.
See [known limitations](docs/KNOWN_LIMITATIONS.md) before using older branches.

## Validation

```matlab
% From the repository root:
addpath('tests');
check_syntax;
smoke_tests;
```

```shell
python tools/validate_repository.py
```

The checks cover manifest completeness, unchanged binary inputs, isolated
workspace creation, recording path resolution, and a MAP prediction.
`compare_original(originalRoot)` additionally compares representative numerical
results against the original archive. [Validation details](docs/VALIDATION.md)
state what was tested and what remains unverified.

## Publishing and attribution

The Git ignore rules exclude raw recordings, original reference documents,
generated results, editor files, and MATLAB caches. Code, experiment manifests,
processed MAT inputs, and Simulink models are included.

Original source credits, including Youhan Li and Jose Gonzalez, are preserved.
No new license or paper citation was invented during reorganization. Add the
appropriate license and final paper metadata when those have been selected.

中文使用说明：[README.zh-CN.md](README.zh-CN.md)。
