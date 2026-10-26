# Workflow entry points

Install the package in the unified PyAutoLens 2025.1.18.7 environment first.

- `prepare_run.py` delegates to `disk-lens prepare-main`.
- `run_workflow.py` delegates to `disk-lens run-main`.
- `recalculate_selection.py` delegates to `disk-lens selection-radii` for the
  separate particle-mass screening calculation.
- `recalculate_stellar_ml.py` delegates to `disk-lens stellar-ml` for the
  real-stellar aperture measurement and mock normalization.

All scientific implementation is in `src/disk_lens/`. See
[the complete workflow](../docs/WORKFLOW.md) for inputs, environment selection
and output paths.
