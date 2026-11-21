# Subject 528: legacy inputs and load-calibrated MAP

The current release refits the two-MTU BSE using the author's original
experimental selection/preprocessing and original trained back MAP. Each
activation source has one fixed 24-parameter physical vector across nine
conditions. Experimental and MAP vectors may differ. Standing balance remains
enabled outside movements; movement PD, assistance and angle limits are zero.

## Run the completed demonstration

```matlab
% From the repository root, after copying the local inputs:
run('Demo/BSE/Direct_Activation_No_Limits/run_dual_MTU_BSE.m')
```

The default evaluates all nine 528 conditions for both activation sources.
It generates torque/reference, theta/reference, back/abdominal torque,
feedback, input comparison and load-scaling figures in `output/`, together
with raw MAT outputs and CSV metrics. Figures are English PNG/PDF/FIG.
`528_overview` shows representative 5/10/15-kg bouts.

The shared `../../format_demo_figure.m` formats every figure with 16.32-pt
ordinary text, 20.4-pt titles, 24.48-pt tick numbers and line widths 1.5
times the originals.
Legends are moved outside plotted data where practical, and axis limits are
preserved when the layout is resized.

Two additional figures are produced by the same script:

- `528_continuous_torque_theta` shows net muscle torque above trunk angle,
  comparing the original CEINMS/IK references with both BSE activation
  sources. One 5-, 10- and 15-kg case is chosen with distinct speed codes.
  The selected combination maximizes the equal-weight mean of Experimental
  and MAP torque and theta Pearson correlations from the existing nine-case
  fitting replay. Selection is descriptive, not held-out validation.
- `528_mean_correlation` shows the equal-weight mean of the nine individual
  Pearson R values for torque and theta, separately for each activation
  source. Each bar reports the number of finite case correlations used.

The continuous example copies the selected recorded action inputs and
references without rescaling or time warping, joining the pieces in their
existing standing intervals. The dual-MTU model is then simulated continuously
once per activation source; muscle and joint states are not reset between
loads. Two seconds from an existing standing interval are added after the last
lift to show recovery. The dotted standing references are the original defined
waiting baselines, not independent measurements. The selected cases, their
correlations in this continuous replay, and the nine-case means are saved as
`528_continuous_example_cases.csv` and `528_correlation_summary.csv`.

The optional C accelerator is compiled on first use. A configured MATLAB C
compiler is required. Set `Use_native=false` to use the complete, equivalent
MATLAB equations in `dual_MTU_model.m`; this route is slower.

## Input preparation

`build_528_inputs.m` follows `Dynamic_kine_reference.m`, `data_combination.m`
and `NNMF_component.m` in the original `My_MAP_model - 5 10 15` project.

- Experimental back combines six representative activations using their
  muscle-group torque contributions as weights. Experimental abdomen uses
  `rect_abd_r` from `Activations.sto`. These are CEINMS activation estimates
  derived from experimental EMG, not direct measurements of internal activation.
- Original five-sample Gaussian smoothing and 1-ms interpolation are retained.
  After waiting intervals are set to zero, back uses Gaussian smoothing over
  300 samples; abdomen uses the original second-order causal Butterworth filter
  at 1/1.2 Hz with a 1000-Hz sample rate. The demo retains its 100-Hz outer grid.
- Force and moment-arm exports are aligned by timestamps before weighting.
  Exact repeated export rows are deduplicated; conflicting rows stop processing.
- Back MAP reuses the author's trained primitives, baseline and polynomial
  coefficients. Its original fusion is preserved literally, including the
  repeated IL term requested by the author:

```text
back_MAP = (MEP1 + MEP2)/12 + (MEP3 + MEP4)/4 + (MEP1 + MEP2)/6
```

Abdominal MAP is trained separately from experimental abdominal profiles using
four NMF temporal patterns, three Gaussian lobes per pattern and the original
additive quadratic load/speed increments around the 10-kg/code-45 baseline.
No back-to-abdomen recruitment, IK target, torque target or BSE output is used
to train that abdominal predictor. `predict_528_MAP.m` combines both predictors.

## MAP amplitude scaling depends on load only

For each load (5, 10, 15 kg) and MTU, the five experimental curves at codes
25/35/45/55/65 are sampled on 200 normalized action-phase points and averaged
with equal condition weights. The minimum and maximum of this mean curve
define the target range. The five raw MAP profiles are averaged in the same way.

```text
gain(load)   = (max(mean_Exp) - min(mean_Exp))
               / (max(mean_MAP) - min(mean_MAP))
offset(load) = min(mean_Exp) - gain(load) * min(mean_MAP)
scaled_MAP   = clip(gain(load) * raw_MAP + offset(load), 0, 1)
```

One affine map per load/MTU is shared by every speed. It is calculated before
BSE calibration and never optimized against theta/torque. Experimental inputs
receive no extra fitted amplitude gain. The earlier fitted abdominal gains
are removed. The old trial-specific 50:150-point valley correction is replaced
by this shared load-only range calibration; no manual per-speed shaping remains.

Two archived EMG row windows extend past the available CEINMS support: 5 kg /
codes 25 and 35. Complete earlier M pairs from the same recordings are used
for amplitude calibration, at 7.400--14.800 s and 5.569--12.309 s from EMG onset.
All five speeds remain included. Original windows, replacements and source
paths are stored in `data/legacy_MAP_528.mat`; unsupported data are not extrapolated.

## Rebuild or retrain

The three readable entry points are:

1. `build_528_inputs`: original recordings and legacy MAP coefficients ->
   experimental profiles, independent abdominal MAP and fixed load-scaling table.
2. `train_dual_MTU_BSE.m`: jointly refit all 24 physical coefficients for each
   source, with fixed input scaling and standing controller. Reference curves
   enter the objective only. Starting files are backed up before replacement.
3. `run_dual_MTU_BSE.m`: replay, audit and plot all nine conditions.

Rebuilding requires the original dataset and old project beside `GitHub_Ready`.
The completed demo uses the packaged data and does not need those raw files.
Rebuilding invalidates the previous physical fit until training is completed.
Training and evaluation here use the same nine conditions; this is a fitting
replay, not held-out BSE validation or an independently identified anatomical model.

## Mechanics and audit

Both torques and theta come from the same compliant Hill system:

```text
I(load) * theta_ddot = back_torque - abdomen_torque
                      + gravity_torque - joint_damping * theta_dot
```

Activation dynamics, force-length/velocity, passive fibers and tensile compliant
tendons are retained. Backward Euler solves fiber equilibrium and joint momentum
together at each 1-ms internal step. There is no angle clipping, boundary spring,
movement target or external assistance. Standing feedback recruits the two MTUs
only while the unchanged waiting flag is enabled. It adds no external joint torque.

`tests/test_load_scaled_528.m` verifies original MAP function parity, experimental
mean extrema, shared load scaling, unchanged reference/timeline data, fixed
physical vectors and zero movement feedback/assistance/boundary torque.
`tests/test_direct_dual_MTU.m` checks native/MATLAB/original-mechanics parity.

The 27 output columns retain the existing audit layout. Theta is column 1,
net/back/abdominal torque are 3/4/5, standing blend/request/PD request are 20:22,
and intent/assistance/boundary torque are 24:26. Columns 20:22 are exactly zero
during movement; columns 24:26 are zero everywhere.

Earlier releases are archived in `experiments/load_scaled_528/before_load_rescale`.
617/725/730 still require this revised input preparation and physical refit.
