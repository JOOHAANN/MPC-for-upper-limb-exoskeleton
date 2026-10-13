# BSE target audit and grouped generalization experiment

From the repository root, run:

```matlab
demo_bse_generalization
```

This creates four PNG/PDF figures, per-trial/per-fold CSV metrics, and a MAT
result file in `results/bse_generalization`. It does not change existing BSE
equations, calibration data, or coefficients. Only base MATLAB is needed for
the new regression benchmark and the archived model audit.

## Two separate questions

1. **What does the archived BSE fit?** `audit_bse_targets` runs the archived
   `general_model_op_scale` variant with its saved coefficients on the three
   combined datasets. It compares column 1 (exactly `310*activation`) and
   column 4 (the retained recorded torque reference). The original nine cost
   windows define the diagnostic segments. This is apparent-fit analysis,
   not out-of-sample validation.
2. **Can a lightweight torque estimator generalize across conditions?** The
   benchmark starts from raw `emgFilt.sto`, `ik.sto`, and `Torque.sto` for 25
   trials, without fitted MAP activations, the combined datasets, or saved
   BSE parameters. All five loads and five speed labels are tested.

The three estimators are a fitted scalar lumped-activation gain plus intercept,
a reduced physics-inspired regression basis, and that same basis plus a ridge
residual model. **Neither baseline is the full original BSE forward model.**
No claim of a 50% improvement to the original BSE rollout or NMPC is made.

## Data handling

- Loads: 5, 7.5, 10, 12.5, 15 kg. Filename codes `75` and `125` are interpreted
  as 7.5 and 12.5; speed values 25, 35, 45, 55, 65 remain labels without
  inventing a physical unit.
- Lumbar and both arm angles are taken from `ik.sto` in their recorded units;
  no per-trial angle-range normalization is applied.
- Six back and two abdominal EMG channels remain separate in the residual.
- Timestamps are aligned on their common interval at 50 Hz. Fixed causal
  exponential smoothers and backward derivatives are used. The first 0.5 s
  is discarded consistently. Resampling does not add independent recordings.
- Linear timestamp interpolation and undocumented acquisition filtering make
  this an **offline** estimator benchmark, not a verified streaming pipeline.
- The target is the recorded `Torque.sto` reference. Its upstream generation
  is not independently established here; it must not be described as a
  verified force-sensor ground truth. A shared upstream model between EMG and
  reference torque could make this benchmark optimistic.
- The included compact cache and SHA-256 provenance permit rerunning the demo
  without shipping the large raw files. To rebuild it from raw recordings:

```matlab
setup_project;
addpath(fullfile(pwd, 'experiments', 'bse_generalization'));
run_bse_generalization(pwd, true);
```

## Evaluation protocol

Outer validation holds out all trials of one load (5 folds) or one speed
(5 folds). Each protocol therefore gives exactly one held-out prediction for
each of the 25 trials. Inner validation holds out the same type of condition
inside the outer training set. Ridge penalties are selected from
`[0.001, 0.01, 0.1, 1, 10]`; disabling the residual is also a candidate.
Feature means/scales are estimated from the current training subset only.
Every trial contributes equal total training weight, irrespective of length.
The principal score is the mean of 25 per-trial RMSEs, in Nm.

The physical basis contains mean back/abdominal activation, sine-of-angle
gravity proxies, load interactions, lumbar velocity, and acceleration terms.
It is a low-order estimator, not a uniquely identified anatomical model.
The residual adds channel-specific activation, an 80 ms activation history,
angles, velocities, and low-order interactions. Nominal load is context:
it does not imply that the object is held throughout every recording.

No torque lag, future target, trial identifier, nominal speed label, or
whole-dataset amplitude normalizer is used as a predictor. The plots display
every fold and four examples selected by condition before inspecting errors.

## Limits

This tests condition transfer within `240528Pilot`, not new subjects, new days,
different EMG electrode placements, unseen assistance torques, or closed-loop
stability. The saved fold records and predictions are not a deployable controller. The
next step is a separately validated residual inside the BSE state equations,
with uncertainty limits and multi-step rollout validation.

Run `addpath('tests'); test_bse_generalization` to check target-independent
features, equal trial weighting, and the saved split records.
