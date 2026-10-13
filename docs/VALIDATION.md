# Validation record

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
