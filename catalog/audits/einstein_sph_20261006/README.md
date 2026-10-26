# Component-specific SPH selection-radius audit

Computed 2026-10-06 for snapshot 84 of TNG50-1. Mass/kinematic preselection
uses the adopted catalogue-frame audit, with 216 passing objects from 260.
The strict three-view maximum radius cut is >0.5 arcsec.

- `selection_216.csv`: all preselected objects, XY/XZ/YZ radii, numerical
  checks, BH correction, multiple-curve diagnostics and result/code hashes.
- `production_37.csv`: exactly the user-approved 37 members.
- `historical_crossers_10.csv`: the available earlier NGP versus new SPH
  numerical comparison; membership is not inferred from those ten candidates.
- `summary.json`: counts, refinements and classification diagnostics.
- `method_checks.json`: recorded analytic,
  normalization, representative-object and execution-equivalence findings. Validation
  scripts and bulk per-object run outputs remain outside the public release.
- `threshold_support_sensitivity.json`: explicit support variations for the
  near-threshold objects 508423 and 509383, distinct from grid/FOV convergence.
- `software_versions.json` and `manifest.json`: environment and checksums.

173 of 216 pass the numerical radius cut;
37 of 37 confirmed members pass. All reported
adopted classifications pass selection-maximum grid convergence, FOV and
threshold-stability checks. The stricter all-view flag remains false for
507525: its weak XY curve is about 0.041 arcsec; its selection maximum is
0.334495 arcsec and is stably below the cut.
Flags control numerical stability within the stated reconstruction; they
do not imply a unique physical density below simulation resolution.

`combined_observed_delta_arcsec` is the sum of observed changes between
grids and fields, not a statistical uncertainty or rigorous error bound.
Initial finest pixels are 0.025 arcsec, with flagged objects refined further.
`method_hash` identifies the batch producer or the production API used for
refinement; the manifest lists both. Source formatting and refactoring
change these hashes even when the mass/lensing numerics agree. The bundled
API additionally reports convergence of the actual screening maximum, while
retaining the stricter all-view flag. Particle loading, view rendering and
mass/lensing mathematical statements are unchanged from the frozen audit
producer, as verified by an AST comparison. The real
production API reproduced all nine batch/field radii for 83999 exactly.

The main result includes bound stars plus wind, DM and gas. Compact BH
dynamical particle mass is a separately labelled diagnostic. Other
structures along the LOS and FoF fuzz are outside this mass domain.

The BH-inclusive column is diagnostic, with its own numerical stability
flag in the tables. Objects flagged there do not define a confirmed
four-component classification. Nominal BH correction changes the cut for
[508423]; BH diagnostic stability is unresolved
for [456102, 502255, 509383]. These are distinct from
the adopted, converged three-component screening.

The approved 37-member list, inherited mock radii and masks, and fit-status
records are independent of this numerical screening. See
[the scientific method](../../../docs/EINSTEIN_RADIUS.md).
