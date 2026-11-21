# Validation record

## 2026-10-15: reduced over-assistance near the torque threshold

With BSE/MAP model bytes, references and activation inputs fixed, training-only
controller selection reduces mean peak headroom from 16.87 to 2.47 Nm in
240627 and 5.83 to 2.66 Nm in 240617, counting only actions needing assistance.
The assistance-torque time integral falls by 34.8% and 21.9%, respectively;
this quantity is not electrical energy. All 48 action/source pairs remain
below the back-torque targets with unchanged valley counts, no QP fallbacks,
and valid 80 Nm / 200 Nm/s actuator commands. Fixed-parameter, zero-assistance,
reference-perturbation, half-step and standing-recovery checks pass. Seven
Python regression checks pass, including release of unnecessary assistance and
preservation of naturally low torque. 240617's transient soft-limit torque is
65.15 Nm. These are offline model replays, not independent assisted experiments.
See [controller method](MPC_THRESHOLD_TRACKING.zh-CN.md) and
[before/after results](mpc_threshold_tracking/threshold_tracking_results.pdf).

## 2026-10-15: quality-first recalibration

Native-clock and activation checks admitted 10 actions from 240617 and 14
from 240627. MAP and all eligible BSE parameters were recalibrated separately
per folder; the other three folders did not meet quality/coverage requirements.
Across 48 action/source pairs, the final MPC replays have zero back-target
excess, preserve reference valley counts, and respect the 80 Nm / 200 Nm/s
actuator limits with zero QP fallbacks. Fixed-parameter, reference-perturbation,
half-step and perturbed standing checks pass. Five critical Python regression
checks pass, including deliberately shifted timestamps and conflicting duplicates.
The Gaussian implementation agrees with MATLAB to 3.05e-16. These checks do
not establish hardware synchronization, independent physiological validity or
real-time performance. 240617 reaches a transient soft-limit torque of 63.15 Nm.
See [method and limitations](QUALITY_FIRST.zh-CN.md) and
[full evidence](subject_quality_v2/quality_first_results.pdf).

## 2026-10-15: recovery MPC and five additional recording folders

The original subject retains its selected Hill parameters and reference Scale.
State-aware recovery support reduces the worst first-second post-action back
peak from 34.545 to 15.016 Nm and the mean from 20.685 to 13.061 Nm. All 18
main W shapes remain. The 80 Nm / 200 Nm/s actuator bounds hold; a 5 kg case
still exceeds its soft back target briefly. Recorded-command replay also
matches the new immutable native-library cache.

The added folders were audited before fitting. Their 79 recorded action
segments produce 158 activation-source case evaluations, with one fixed
parameter vector per source and one MPC setting per folder. All five pass
reference perturbation invariance, zero-assistance parity, 30-second perturbed
standing recovery and half-step numerical checks. The largest half-step
differences are 0.004820 rad and 2.105 Nm. The maximum case mean abdominal/back
torque ratio is 0.210; three rejected QP candidates use bounded fallbacks.
These checks establish implementation consistency, not universal physiological
accuracy or sealed-test generalization.

The MATLAB subject viewer was executed with all 21 images and six tab groups.
The 21-page PDF was rendered and inspected, including text-bound checks.
Missing conditions, partial recordings, 240723 angle polarity and weak 240719
MAP validation are retained explicitly. See [the method/data audit](MULTI_SUBJECT.zh-CN.md)
and [machine-readable checks](multi_subject_evaluation/validation_summary.json).

The 2026-10-14 Hill/MPC repair also passed a full run through the MATLAB
`demo_legacy_mpc` entry: both sources, all nine conditions, 18 W-shaped
movements, feasible/rate-limited actuator commands, zero external posture
torque, zero-assistance parity with the MATLAB model, half-step replays and
30-second standing extensions. See [the MPC validation](HILL_MPC.zh-CN.md).

Validated on 2026-10-12 using MATLAB R2023b on Windows.

| Check | Result |
| --- | --- |
| Workflow source/input references | All 38 manifests resolve |
| MATLAB source coverage | All 174 canonical research files referenced |
| Binary input preservation | All 132 unique MAT/Simulink assets match original SHA-256 |
| Workspace materialization | All 38 experiments prepared successfully |
| MAP smoke calculation | Finite 200-by-6 activation output |
| Raw recording resolver | Legacy and supplied folder layouts resolve to the same recording |
| Original/new MAP results | Exact equality of weights, primitives, and activation |
| Original/new low-pass filter | Exact equality for a deterministic mixed-frequency signal |
| Original/new geometry dynamics | Exact equality for 500 samples; finite outputs |
| MATLAB Code Analyzer | No syntax diagnostics; legacy quality warnings remain |
| MISS_HIT parser/formatter | All 174 research source files parsed successfully |

The regression comparison uses the subject-specific MAP mode at speed 45 and
load 10, the original second-order causal Butterworth filter, and the geometry
calibration demonstration. It does not establish equivalence for every branch
or prove closed-loop stability.

Re-run from the repository root with `addpath('tests'); check_syntax;
smoke_tests;`. Run `compare_original(originalRoot)` with the path to the backed-up
original project for the numerical comparison. Analyzer output goes to
`results/code_analyzer.json`; generated workspaces are ignored by Git.

## Selected neural/Hill model (2026-10-14)

The current `demo_hill_sequence` prefers `neural_models.mat`. Its fixed
parameters were checked across all nine conditions per activation source
(18 cases, six continuous records). These are calibrated replays, not an
independent holdout or individual physiological validation.

`test_hill_neural` passed C/readable-MATLAB parity, muscle torque balance,
future-input isolation and absence of external joint PD torque. The previous
`test_hill_compliant` regression also passed after accelerator rebuilding;
cached solver comparisons allow floating-point compiler rounding, while
reference-input isolation compares two fresh executions exactly.

The selected depth-priority model has 42 parameters per activation source,
plus a common standing controller. Both flexion valleys span 0.488–0.596 rad
for experimental inputs and 0.470–0.654 rad for MAP inputs. All nine cases
retain two major valleys and the middle rise. The relaxed −0.3–2.3 rad range
is a screening criterion, not an angle clip.

Complete MATLAB replays matched the C calibration solver within 5.1e-10 rad
and 1.2e-7 Nm. Halving the internal time step changed angle by at most 0.00309
rad and torque by 1.900 Nm. Adding 30 seconds of standing to each record,
without state resets, gave at most 8.4e-10 rad error in the final ten seconds.
The original angle and torque references were preserved exactly; enabled
elastic soft boundaries produced zero torque throughout these replays.

Reports, equations, fixed-parameter hashes, numerical checks and figures are
in [HILL_NEURAL.zh-CN.md](HILL_NEURAL.zh-CN.md) and
`docs/hill_neural_evaluation/`. The old speed-holdout calibration is a
different model and does not establish the new model's generalization.

## Optimized-model / MPC figure suite (2026-10-15)

The plot-only entry `demo_hill_mpc_plots` passed in MATLAB R2023b: Python
rendered all 11 figures, and the MATLAB viewer created 11 tabs and 11 image
objects. The full simulation entry now uses the same renderer and viewer.
The Python plotting and evaluation modules also passed syntax compilation.

The combined PDF has 11 nonempty pages, with 18 unique source/load/speed rows
in `model_mpc_comparison_metrics.csv`. PNGs and publication file hashes were
checked. The saved MPC trajectories and selected parameter cache retained
their original SHA-256 hashes; plotting did not rerun fitting or simulation.
Both inputs cover all nine cases. Filled standing references are marked,
and muscle torque is distinguished from muscle-plus-actuator torque.

## Stronger MPC back-torque penalty (2026-10-15)

Both weights 8 and 16 were replayed on all six continuous records with the
same calibrated Hill parameters. Weight 16 is now the default (previously 2);
angle weight 10, thresholds 100/140/180 Nm, actuator cap 80 Nm and slew limit
200 Nm/s are unchanged. All 18 conditions reduce both mean positive threshold
excess and peak excess, and retain two principal angle valleys. Exceedances
are reduced rather than eliminated; no output clipping is used.

The selected replay passed zero-assistance parity, original-reference equality,
actuator amplitude/rate bounds, muscle dominance and all-six numerical checks.
Halving the integration step changes angle by at most 0.003245 rad and net
muscle torque by 2.065 Nm. All extended 30-second standing errors are below
9e-10 rad. There are zero rejected QP results and zero elastic-boundary or
external posture-balance torques. Before/after comparison checks identical
physical/neural parameters, activations, time and reference arrays.

## MPC headroom and angular-velocity preservation (2026-10-15)

Eleven settings including the preceding controller were compared on all nine
conditions per source. The selected shared controller uses a 20-ms update,
35-Nm headroom, back weight 128, angle weight 10 and velocity weight 100.
Actuator cap/rate remain 80 Nm and 200 Nm/s. The selected input forecast is
causal hold; a past-input trend option was tested but not selected.

Compared with weight 16, mean positive threshold excess changes from
1.38498 to 0.01310 Nm (experimental) and 1.62566 to zero (MAP). Seven/nine
experimental cases and all MAP cases have zero exceedance; eight/nine
experimental cases stay within 1 Nm of the threshold. Remaining time above
threshold + 1 Nm totals 0.12 s / 0 s. The worst 7.84336-Nm peak coincides with
80-Nm actuator saturation. All 18 angle curves retain two principal valleys.

The full numerical evaluation passed. Half-step replay differs by at most
0.003641 rad and 2.633 Nm; added 30-second standing errors are below 9e-10
rad. All zero-assistance parity, reference equality, actuator bounds/rate,
muscle dominance and momentum checks pass, with zero rejected QP results.
The new optional C forecast has exact zero-slope/hold parity, leaves the
input state unchanged, and responds to nonzero past-input slopes. Human
muscle dynamics and calibrated parameter hashes are unchanged.

## Post-movement transient diagnosis (2026-10-15)

All 18 endings were analyzed with matched end-state ablations of standing P,
standing D, standing PD, movement feedback and both feedback layers. Every
counterfactual uses the same subsequent activations, mass and recorded
actuator command for up to one second; it does not rerun MPC optimization.
The baseline replay splits the last internal step to recover the exact q/v
used for feedback. It passes angle/back-torque parity and reconstruction of
the recorded blended correction, without altering the selected trajectories
or parameter cache. Diagnostics and the four-page report are linked from
[HILL_MPC.zh-CN.md](HILL_MPC.zh-CN.md).
