# Preserved research limitations

This reorganization preserves experiment-specific scientific choices. It does
not certify every archived optimization or controller branch.

The multi-subject pipeline audits five additional recording folders and fits
MAP/Hill/MPC parameters per folder. Only 240617 contains the nine named target
conditions, with some truncated records; other folders have fewer speeds or
missing loads. The 240723 angle polarity and 240719 10 kg angle alignment are
flagged, not silently corrected. Holdouts were inspected during development
and are not a sealed final test. Predicted assisted muscle torque still relies
on a modeled human response rather than an independently fitted assistance
experiment. [Full data audit and method](MULTI_SUBJECT.zh-CN.md).

The current `demo_hill_compliant` preserves the exact archived angle reference
and evaluates a coupled compliant Hill pair. Its auxiliary flexor recruitment,
effective load geometry and muscle priors are estimated and not independently
validated. Its 32 parameters include shared body/load geometry and elastic
soft limits; these values do not establish anatomical identifiability.
Whole-speed holdouts condition on the archived MAP inputs; they do
not establish end-to-end or cross-subject generalization. See
[the current method and limitations](HILL_COMPLIANT.zh-CN.md).

- The main project description identifies a back-support exoskeleton and
  L5/S1 model despite the outer folder's upper-limb name.
- `training_map__muscle_tendon__general_model_op_damp` and
  `training_map__muscle_tendon__general_model_op_scale` contain a top-level
  `MAP_output_model_fusion` script without local `baseline_data.mat` and
  `Combination_raw.mat`. Use the matching `__build_data` workflow to prepare
  those inputs. The static missing-load inventory is in
  `missing_archived_inputs.json`. The manifests do not guess which case data
  to use for an incomplete historical directory.
- Legacy `getWeights` implementations use speed/load polynomial coefficients
  outside the subject-specific branch even though those coefficients are
  loaded inside that branch. The demonstrated and tested mode is `SM`.
  The historical `GM` branch requires further review.
- The Qualisys `from_lmt_to_mtf/Cost_function_from5to15` implementation has a
  possible unset output reported by MATLAB Code Analyzer. It was preserved
  rather than assigning an invented objective value.
- Some geometry calibration scripts have unequal lower/upper bound lengths
  and seeds outside the specified bounds. These are original settings;
  review the intended parameter constraints before training.
- Sample boundaries, channel indices, load transitions, and model constants
  remain specific to the supplied experiments. Normalization helpers retain
  their original behavior for degenerate inputs.
- Original scripts retain some workspace-style MAT loads, global variables,
  unused exploratory variables, and dynamically grown arrays. Code Analyzer
  warnings remain; there were no syntax diagnostics in the validation run.
- The original project description notes slow controller computation, limited
  fidelity of activation response to assistance torque, and EMG noise. An
  offline simulation is not evidence of real-time hardware performance.

The archived experimental and MAP MPC simulators now have complete 45.4-second
checked replays. They show angle excursions, failed optimizer exits, and
assistance exceeding the configured cap. They use the archived BSE model,
not the selected neural/Hill model. See [the MPC test report](LEGACY_MPC_TEST.zh-CN.md).
The default MPC entry now uses the selected neural/Hill model and preserves
the W-shaped angle pattern in all nine cases per source. It assumes human
excitation adaptation to known assistance torque; this response is not fitted
to independent assisted EMG measurements. Remaining phase/torque mismatch is
documented in [the current MPC report](HILL_MPC.zh-CN.md).
Full GA calibration, every NMPC branch, Simulink execution, and every raw
recording preprocessing route have not all been rerun.
