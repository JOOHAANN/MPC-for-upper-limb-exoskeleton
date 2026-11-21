# Direct-activation dual-MTU BSE

This separate demonstration retains the existing back and abdominal MTUs and
their coupled tendon/fiber/joint mechanics. It removes the movement-intent
decoder, movement posture feedback and both angle boundaries. The previous
BSE and NMPC demonstrations are not replaced.

**The physical/recruitment coefficients are transferred from the current
calibration. They have not been refitted after removing movement feedback.**
A completed numerical replay is not a claim of accurate fitting or validation.
The actual accuracy and any failed trajectories are reported without clipping.

## Run

```matlab
run('D:/matlab/Journal Paper of Upper-Limb Exoskeleton/GitHub_Ready/Demo/BSE/Direct_Activation_No_Limits/run_dual_MTU_BSE.m')
```

The default evaluates 528, 617, 725 and 730, with nine recorded conditions and
both experimental and frozen trained-MAP activation per subject. To select one:

```matlab
Subject_list = {'528'};
Show_figures = true;
run('D:/matlab/Journal Paper of Upper-Limb Exoskeleton/GitHub_Ready/Demo/BSE/Direct_Activation_No_Limits/run_dual_MTU_BSE.m')
```

The default compiles a small C accelerator on first use. A configured MATLAB C
compiler is required for that route. Set `Use_native=false` to execute the full
readable MATLAB equations instead; this is slower. The two implementations are
checked against one another and against the existing dual-MTU physics.

## Files to read

- `run_dual_MTU_BSE.m`: subject/case loops, simulation, checks and metrics.
- `dual_MTU_model.m`: complete MATLAB equations and optional accelerator calls.
- `prepare_activation.m`: unchanged input preparation and abdominal MAP generator.
- `plot_direct_BSE.m`: torque, theta, muscle torque and standing-feedback figures.
- `data/`: four subject packages and source provenance.
- `private/native/`: equivalent accelerator; not a second model or controller.

## Feedback and mechanics

External assistance is zero throughout this BSE-only demonstration. During
movement, each MTU receives its prepared experimental/MAP excitation directly:

```text
command = prepared_input
activation_target = command ^ activation_exponent
d(activation)/dt = (activation_target - activation) / activation_time
```

There is no target movement angle and no pose-error-driven change to those
commands. Excitation/activation remain physiologically bounded to [0,1]; this
does not bound the joint angle. Input preparation retains the existing fitted
abdominal recruitment terms and MAP generator. It uses input history/load,
not the current simulated angle or an angle reference.

Standing feedback uses the packaged subject balance coefficients and fixed
standing angle. The recorded waiting flag is its only mode gate. Waiting
feedback blends in with the original time constant and exits immediately at
the first movement sample, so it cannot track reference motion during actions.
This explicit flag replaces the old low-activation proxy for standing.

The standing PD request and gravity compensation request neural recruitment
of the two MTUs. They do not add a hidden mechanical joint torque. The code
does not reuse the old experiment's direct external `T_balance` joint torque;
it retains the existing dual-MTU standing-recruitment implementation.

Both muscle torques and theta are outputs of the same compliant-Hill system:

```text
I(m) * theta_ddot = back_torque - abdomen_torque
                    + gravity_torque - joint_damping * theta_dot
```

The tendon is tensile; normalized fiber force includes active force-length,
force-velocity, passive stretch and fiber damping. Backward Euler solves the
fiber equilibrium and joint momentum together at each 1-ms internal step.
No angle spring, boundary damping, angle saturation or angle wrapping exists.
Newton-step bounds and fiber root brackets protect numerical solution; they
do not constrain the physical joint angle. Failed equilibrium solves stop and
leave subsequent samples missing rather than repairing the plotted output.

One fixed 28-coefficient vector per subject/source is used across all nine
conditions. All 16 prior intent/neural parameters and four angle-boundary
parameters are absent from these new model vectors. No parameter is switched
by load or speed. Muscle and joint states persist between loads within each
recorded sequence. Separately recorded sequences initialize independently.

## Outputs

`output/` contains raw MAT results, all 72 case/source metric rows, and four
figure types per subject in PNG/PDF/FIG: net muscle torque versus reference,
theta versus reference, back/abdominal torques, and standing-feedback requests.
Original reference scales and unsupported CEINMS intervals are retained.

The 27 output columns retain the existing audit layout:

```text
1 theta             2 angular velocity        3 net muscle torque
4 back torque       5 abdomen torque          6:11 force contributions
12:13 fiber length  14:15 activation           16 fiber residual
17 momentum residual 18:19 excitation command  20 standing blend
21 standing request 22 standing PD correction 23 integral memory
24 intent: zero     25 assistance: zero       26 angle boundary: zero
27 gravity torque
```

Column 22 is a neural-recruitment correction request, not an added external
torque. Its value and the standing blend must be exactly zero during movement.
Columns 24:26 must be zero everywhere. The demo checks these requirements and
that movement commands equal the prepared inputs sample by sample.
