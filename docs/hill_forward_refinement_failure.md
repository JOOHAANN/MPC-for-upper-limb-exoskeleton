# Rejected exploratory refinement

After the 180-seed screening and six-start Hill calibration, an additional
experiment initialized each source's optimizer with parameters fitted to the
other source on the same training conditions. No held-out labels entered the
optimizer, but an improved training objective did not ensure a valid forward
prediction on held-out conditions. At least one prediction crossed the numerical
failure threshold, so the experiment aborted before replacing saved parameters.

The shipped default retains the preceding calibration. This choice was informed
by development-condition results; the reported leave-speed-out evaluation is
therefore exploratory, not a locked independent benchmark. The optional third
argument to `run_hill_forward` reproduces this refinement experiment and may
abort. It is disabled by default and is not required to reproduce the baseline.

This failure is evidence against claiming that the single lumped rigid-tendon
Hill model is already a stable, validated predictor for arbitrary activations.
The previous angle-equilibrium surrogate has not been restored as a fallback.
