# MATLAB source organization

- `shared`: identical implementations shared by several experiments.
- `variants`: reused implementations with names that also have other versions.
- `workflows`: implementations unique to one experiment.

Do not use `addpath(genpath(...))` on this tree: different experiments contain
incompatible functions with the same name. `config/workflows/*.json` maps the
correct canonical files to each experiment. `prepare_workflow` creates an
isolated working directory with the expected original filenames.

Identical files were consolidated by SHA-256, not by filename. Distinct model
equations, calibration settings, and case-specific data were retained.
Original author credits remain in the source files.
