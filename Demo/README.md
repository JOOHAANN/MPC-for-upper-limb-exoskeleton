# MATLAB BSE and controller demos

## Current subject-528 workflow

The latest recalibrated two-MTU model is in `BSE/Direct_Activation_No_Limits`.
The current PID/NMPC comparison loads that exact calibrated model, with
experimental and trained-MAP inputs and fixed source-specific parameters.

From this folder, run the BSE replay:

```matlab
addpath('BSE/Direct_Activation_No_Limits')
run_dual_MTU_BSE
```

Run the controller comparison:

```matlab
addpath('NMPC')
offline_PID_NMPC_comparison_528
```

The controller entry runs nine cases per activation source and exports the
same RQ1 torque, RQ2 angle/deviation/timing, continuous-lifting and model-audit
figure types. See [NMPC instructions](NMPC/README.md) for settings and outputs.
The BSE folder is not modified by the controller-folder cleanup.

## Local inputs

A Git clone or linked worktree includes the demo source. To replay the current
calibration, copy the following directories from an authorized local research
checkout to the same relative paths within the new checkout:

- `Demo/BSE/Direct_Activation_No_Limits/data/` for the current BSE, MAP and NMPC replay
- `Demo/NMPC/data/` for the saved activation residual predictor and controller settings
- `Demo/BSE/data/` for the older multi-subject BSE demo

These directories are Git-ignored. The demos generate their outputs under their
own `output/` directories. To regenerate the continuous torque plot using a saved
run, copy `Demo/NMPC/output/current_BSE_528/current_BSE_NMPC_528.mat` to the
corresponding path and run `redraw_RQ1_528`.

## Model rules

Both joint motion and muscle torques come from the same coupled Hill-MTU
and joint-momentum equations. Unassisted movement uses the prepared activation
directly, without movement PD or angle limits. Standing/waiting feedback is
retained. During assisted movement only, human PD may correct back activation
against a causal unassisted shadow. Abdominal inputs remain measured or
separately MAP-predicted channels.

The controller uses identical actuator constraints for PID/NMPC, a threshold
at 75% of the unassisted back-muscle peak and a full upper-clipped NMPC torque
target. These are offline calibrated-model replays; the retained result
limitations are documented in `NMPC/README.md`.

## Older model revisions

The original BSE programs and other-subject packages elsewhere under `BSE`
remain available as earlier revisions. The cleaned NMPC folder supports the
latest 528 comparison; old multi-subject and no-human-PD controller entry
points were archived outside `Demo`. The previous demo guide is preserved in
`../experiments/NMPC_cleanup_20261109_214149/before_edit/Demo_README.md`.

Keep the entire `Demo` folder together: NMPC reads its physical model from the
BSE subfolder. MATLAB R2023b and Optimization Toolbox are used; the current
physical accelerator requires a configured C compiler when rebuilding.
