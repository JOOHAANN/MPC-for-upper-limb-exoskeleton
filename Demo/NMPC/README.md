# Subject 528: PID versus NMPC

## Run

Open `offline_PID_NMPC_comparison_528.m` and click **Run**, or enter:

```matlab
cd('D:\matlab\Journal Paper of Upper-Limb Exoskeleton\GitHub_Ready\Demo\NMPC')
offline_PID_NMPC_comparison_528
```

The script loads the latest calibrated dual-MTU package from
`../BSE/Direct_Activation_No_Limits/data/subject_528.mat`. Experimental and MAP
each use one fixed physical vector over all nine load/speed conditions.
The trained BSE, human feedback and NMPC settings are unchanged. PID gains
are tuned separately for the Experimental and MAP activation sources.
MATLAB R2023b, Optimization Toolbox and a configured C compiler are required.
The supplied Windows MEX accelerator rebuilds from source if necessary.

The three settings at the beginning of the script are optional:

```matlab
Source_list = 1:2;       % 1: Experimental; 2: MAP
Show_figures = true;    % Display figures as well as saving them
Reuse_cache = true;     % false recomputes the complete closed-loop simulations
offline_PID_NMPC_comparison_528
```

## Files

| File or folder | Purpose |
| --- | --- |
| `offline_PID_NMPC_comparison_528.m` | Single main script: load, simulate, check and plot |
| `current_528_settings.m` | PID, source-specific NMPC and human-PD settings |
| `current_bse_plant.m` | Latest BSE with actuator and assistance-gated human PD |
| `simulate_threshold_pid.m` | Reactive PID on current back-muscle torque |
| `simulate_threshold_nmpc.m` | NMPC following the complete clipped torque reference |
| `nmpc_forecast.m` | Trained MAP preview plus past activation residuals |
| `nmpc_file_hash.m` | Cache/model identity check |
| `../format_demo_figure.m` | Shared font sizes, line widths and readable legend layout |
| `private/` | Plotting, metrics and the physical MEX implementation |
| `data/` | Trained residual predictor and current controller-setting exports |
| `output/current_BSE_528/` | One current result set, figures, metrics and replay caches |

## Retained figures

Both sources retain the same figure types as before:

- `RQ1_torque_all_nine_*`: unassisted/PID/NMPC muscle torque, PID/NMPC assistance,
  upper-clipped target, muscle threshold and common actuator cap.
- `RQ1_continuous_example_*`: continuous 5/10/15-kg example with standing waits.
- `RQ2_theta_all_nine_*` and `RQ2_continuous_example_*`: absolute angle trajectories.
- `RQ2_theta_deviation_all_nine_*`: angle changes relative to unassisted BSE.
- `RQ2_timing_all_nine_*`: demand/assistance timing without temporal shifts.
- `528_nmpc`: torque and angle versus experimental references, thresholds and assistance.
- `528_muscle_components`: abdominal and angle-boundary torques.
- `528_assistance_PD_audit`: assistance, back-input corrections, abdominal torque
  and exact standing/movement feedback gates.
- `PID_NMPC_band_entry_528`: the continuous example with threshold bands and
  the first qualifying post-peak entries marked. The event and summary CSVs
  include every eligible interval, including those without a sustained entry.

Figures use 16.32 pt for ordinary text, 20.4 pt for most titles,
and 24.48 pt for axis tick numbers. The MAP continuous RQ1 and RQ2
titles use 30.6 pt and 40.8 pt, respectively. Existing line widths are multiplied
by 1.5 once. Legends are placed outside the data region where practical.
The formatter preserves plotted values and pre-format axis limits.

There are 16 figure types/files in each PNG/FIG set, with PDF exports for the
RQ figures, band-entry report and PD audit. MAT files preserve the raw trajectories; CSV files
contain all 54 controller/source/case metric rows. Cached runs are reused only
when model, physical/solver code, inputs and active settings match.

## Unchanged model and comparison rules

The reactive PID is implemented in incremental form, with the next reachable
actuator command retained as its memory. This prevents a rate-limit anti-windup
rule from discarding the whole integral increment. Its only control inputs are
current back torque, the threshold and past controller state; it has no future
activation preview. `integralNm` is the equivalent integral contribution
after saturation tracking and can be negative.

| Activation input | Kp | Ki | Kd | Derivative filter | Withdrawal filter |
| --- | ---: | ---: | ---: | ---: | ---: |
| Experimental | 0.83 | 3.4 | 0.12 | 0.027 s | 0.065 s |
| MAP | 0.3901 | 7.8257 | 0.1786 | 0.0197 s | 0.0334 s |

Each source uses its single fixed tuning in all nine cases, and RQ1/RQ2 use the
same resulting trajectory. No BSE, human-PD, NMPC or actuator settings changed.
The following comparison uses the saved MATLAB trajectories from before and
after this PID tuning, with exactly the same unassisted and NMPC trajectories:

| Source | Sustained entries | Settled by demand end | Total excess area | Mean angle RMSE |
| --- | ---: | ---: | ---: | ---: |
| Experimental | 5/14 -> 7/14 | 4/14 -> 6/14 | 83.26 -> 82.20 N m s | 2.85 -> 3.02 deg |
| MAP | 6/18 -> 7/18 | 0/18 -> 0/18 | 175.72 -> 178.72 N m s | 4.34 -> 2.83 deg |

These are action intervals, not counts of the nine load/speed conditions.
The MAP tuning reduces the average fast variation of muscle torque, but a few
intervals still oscillate and its total excess area increases slightly. It is
not a claim that every action settles at the threshold. The gains were checked
with three illustrative synthetic noise sequences, not measured sensor noise.

For the band-entry report, an eligible interval is one where unassisted back
torque exceeds the threshold for at least 0.30 s. After the largest overshoot,
the report finds the first period of at least 0.15 s within +/-max(5 N m, 5% of
the threshold). `StayedUntilDemandEnd` means it stays there after that first
entry. `SettledByDemandEnd` also allows later re-entry and requires at least
0.15 s inside the band at the end of the recorded demand interval. Entry and
final-settling times are measured from demand onset; confirmation includes the
0.15-s observation. Unsuccessful intervals remain in the CSV. These are
finite-action measurements, not standard fixed-setpoint step-response settling.
The report counts separate excursions above or below the band using a 0.05-s
moving mean only for counting; plotted torques and all band-entry times remain
unfiltered simulation outputs. `AtMostTwoOvershootsAndSettled` counts events
that meet the finite-interval settling criterion after at most two overshoots.

The current calibrated BSE's back-activation time constant is 0.20 s; the
archived model used about 0.005 s and included a different assistance-dependent
neural recruitment rule. Consequently the earlier controller and current PID
do not operate on the same plant. Changing the BSE time constant solely to
improve the controller plot would require a new calibration of the BSE itself.

Movement PD uses the original P=5, D=3 structure to correct back activation only
while actual assistance is nonzero during a movement. It tracks the current
causal unassisted shadow, not recorded or future IK. At zero assistance,
movement PD is off and the plant matches the standalone calibrated BSE.
Standing/waiting feedback remains active. Abdominal movement input remains
the experimental or separately trained MAP channel. Angle-boundary torque is zero.

PID and NMPC share assistance in [0,80] N m, positive rise <=200 N m/s and
instantaneous withdrawal capability. PID includes withdrawal-only command
shaping. The threshold is 75% of the unassisted back-muscle peak per case.
NMPC tracks `min(unassisted back torque, threshold)` with a 0.5-N-m feasibility
reserve. Control effort is reported, not matched. Actual limits and separate
Experimental/MAP horizons are listed in `data/current_controller_settings_528.csv`.

The unchanged NMPC results retain their limits: Experimental maximum excess
is 13.7069 N m (15 kg / code 35); MAP maximum
is 0.4074 N m. Mean angle deviations are 1.3911 and 1.0335 deg. The recorded
solver does not consistently meet the 20-ms update budget. These are offline
calibration-set replays, not real-time hardware or unseen-condition validation.
Full details remain in `output/current_BSE_528/release_audit.json`.

## Historical files

Old entry points, unused data, no-human-PD variants, unused native solvers and
30 historical result folders were moved outside `Demo` to:

`../../experiments/NMPC_cleanup_20261109_214149/retired/`

`moves.json` in that archive records every move. The former main-script runner
was merged into the retained entry script; plotting and metrics are private
helpers. The current settings no longer load a historical controller package.
