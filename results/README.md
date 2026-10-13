# Experiment outputs

`prepare_workflow` creates a fresh subdirectory for each run. It materializes
the exact source variants and input assets selected by a workflow manifest.
Run scripts there so saved coefficients, simulations, and figures do not
overwrite the repository inputs. These generated folders are ignored by Git.

Edit the canonical files listed in the manifest to make lasting code changes.
Re-prepare a workflow after editing source. Transfer calibrated parameters to
the next stage deliberately; preparing another workflow uses its archived inputs.
