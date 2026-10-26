# Updated production mock-model Einstein radii

All 37 radii were recalculated on 2026-10-08 from the actual SPH stellar
template and analytic Ludlow NFW model, using central deflection radius
zero. The generated source images retain radius 0.2 arcsec. Each radius
passed fixed-model contour-sampling convergence and boundary checks.

The primary radius is the largest closed critical-curve area-equivalent
radius. Additional loop areas and the historical sum of individual radii
remain diagnostics. This is distinct from the TNG particle-mass selection
radius and from inherited `_legacy` mock radii. The differences table is
a historical comparison, not a controlled single-change experiment.

`mock_truth_37.csv` and `measurements_37.json` provide radii, actual model
parameters and all evaluation-grid results; input and method hashes are
recorded in the manifests. Full contour coordinates and unmasked native
deflections are saved with each external run's `truth/` bundle.
See [method](../../../docs/MOCK_TRUTH.md).
