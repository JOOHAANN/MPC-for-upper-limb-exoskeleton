# Torque and trajectory control without human posture PD

Run the separate subject-528 experiment:

```matlab
run('D:/matlab/Journal Paper of Upper-Limb Exoskeleton/GitHub_Ready/Demo/NMPC/offline_PID_NMPC_no_human_PD_528.m')
```

The previous four-subject and PID/NMPC demonstrations are not modified.
This script uses experimental activation and trained MAP activation, each
with its existing subject-specific Hill parameters across nine cases:
5/10/15 kg at speed codes 35/45/55. Three continuous sequences include standing
waits and retain their complete muscle and joint states between actions.

## What is removed

The actual human plant and every NMPC prediction have movement posture gains
`neuralParameters(4:5) = 0` and standing posture gains
`balanceController(1:3) = 0`. The inactive standing integral memory is zeroed.
Output columns 22 and 23 are checked to be exactly zero. Neither trajectory
error nor desired angle enters a muscle command in this plant.

The same compliant Hill muscle, tendon and joint equations are retained.
Gravity compensation, excitation filtering, assistance-dependent excitation
response, passive mechanical damping and the original soft boundaries remain.
They are separate from posture-error feedback. Muscle torque is output column
4, abdominal torque column 5 and signed soft-boundary torque column 26.
The ablated plant is not retrained and cannot inherit the original model's
fit or validation claims. Its uncontrolled motion can drift substantially.

## Shared external actuator and references

PID and NMPC both command a signed exoskeleton torque in **[-80, +80] N m**.
Positive torque extends the trunk; negative torque provides active flexion
or extension braking. Increasing torque magnitude is limited to 200 N m/s.
Magnitude withdrawal is immediate. Reversal first withdraws to zero, then
builds torque in the opposite direction at the same limited rate. Both
controllers update every 0.02 s. Plant output is every 0.01 s, with 0.001-s
internal integration. This bidirectional device assumption differs from the
previous one-way demonstration; it is not measured hardware behavior.

`data/no_human_PD_task_528.mat` contains a fixed desired movement for each
sequence, generated **offline** from the previously calibrated unassisted MAP
model. Its original posture gains are used only in that reference-generation
step. There is no posture PD inside either new actual or prediction plant.
Both controllers receive the same task; PID uses its current sample, NMPC
uses its planned horizon. No future experimental activation or IK is used.
This is a declared known-task reference, not online unknown-intention inference.

Muscle thresholds remain 75% of the previously calibrated, successful-task
unassisted peak within each condition and source. They match the preceding
demonstration. The drifting no-human-PD ablation is not used to redefine the
required-task threshold. Both controllers use the same thresholds. The
frozen MAP task also supplies the desired torque, upper-clipped at that limit.

## Controllers and the two research questions

The new external PID includes the previously damped threshold-reactive torque
branch plus a current-angle PID branch. Only the external exoskeleton applies
this angle correction. Giving both methods the same torque/trajectory goals
avoids comparing a trajectory-aware NMPC with a torque-only PID. The torque
branch reacts to threshold exceedance; the angle branch can act below that
threshold and remains active during standing waits. One PID tuning is shared
across all loads, speeds, sources and both research questions.

NMPC predicts the no-human-PD Hill model using the trained MAP-plus-past-ARX
activation forecast. Its objective includes task-angle and velocity errors,
tracking the upper-clipped frozen MAP task torque, threshold exceedance,
and small effort/move penalties. Exceedance has a soft penalty because finite
actuation may not satisfy both goals simultaneously. The 0.5-N-m reserve is
part of this penalty, not a hard guarantee on actual torque. Angle and muscle
torque are never clipped or overwritten after simulation. PID and NMPC have
matched equipment constraints; squared-torque effort is reported, not matched.

RQ1 reports peak, area and duration of threshold exceedance. RQ2 reports
tracking error relative to the fixed task trajectory and also records change
from the uncontrolled baseline. The uncontrolled no-PD trajectory is not the
tracking target. Positive-assistance/demand lag is reported separately and
can be ambiguous when active braking is important. Smaller angle error alone
does not establish reduced hardware latency, synchronization or comfort.

Outputs are separate in `output/no_human_PD_528/`: raw MAT trajectories,
per-case CSV metrics, all-nine torque and theta comparisons, continuous
examples, and abdominal/soft-boundary torque audits. Figures are exported as
PNG, vector PDF and editable FIG. Parameter choices are calibrated and checked
on these recorded cases, not an independent subject or hardware test.

`Source_list`, `Group_list`, `Show_figures` and `Reuse_cache` optionally select
the run. Leave `Max_samples` unset for the full experiment. Development probes
use separate output labels; do not mix them with the complete nine-case result.
