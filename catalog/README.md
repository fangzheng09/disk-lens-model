# Sample membership and per-system inputs

The active sample is `production_37` in `sample_sets.json`: the original
33 objects plus exactly **4, 83993, 83999 and 297779**, confirmed by the user
on 2026-10-06. `additional_4` records these additions. The older
`additional_6`, ten velocity-threshold crossers and historical membership
review remain provenance records, not active sample definitions.
Validate the structure with `disk-lens validate-catalog`.

`production_systems.json` contains inputs for exactly these 37 objects and
is the bundled fallback for the fitting pipeline. The source position,
image positions, fitting mask and historical normalization values are preserved
from `system_parameters.provisional.json`. **Membership is confirmed; stellar
M/L and adopted normalization have been recalculated. Historical mock grids
have been audited, and all members adopt 0.2 arcsec central source suppression.
The remaining inherited values are adopted fiducial inputs; experiments
record explicit per-run overrides.** The original
41-system file is retained as the historical input record.

`fit_status.production.json` records the user-confirmed historical paper
20/13 success/failure split and successful additions 4, 83993 and 83999.
The [priority audit](audits/fitting_priority_20261009/manifest.json) matches
the successful 20 to `finish/` and the paper analysis entry points.
`main_fit_priority_23` is the adopted first queue. `main_fit_deferred_14`
contains the original unsuccessful 13 plus unclassified addition 297779.
All 37 still require a uniform review of the new SPH fits; these historical
labels determine scheduling rather than certify the updated fit outcomes.

Small selection-audit CSVs and their input/code checksums remain under
`audits/selection_20260910/`. Historical Einstein-radius estimates used
nearest-grid-point mass deposition. The new component-specific SPH results
are in [`audits/einstein_sph_20261006/`](audits/einstein_sph_20261006/README.md),
with all 216 preselected objects and a separate 37-member table. These are
particle-mass selection radii. The [stellar M/L update](audits/stellar_ml_20261008/README.md)
synchronizes each production record with its YZ aperture M/L and mock coefficient.
Other adopted mock inputs retain their historical provenance.

The adopted kinematic selection keeps `SubhaloPos/SubhaloVel`, corrected
particle velocities (`Velocities * sqrt(a)`), a physical 30 kpc aperture and
`kappa_rot > 0.45`. Alternative stellar-centroid and iterative-frame columns
remain sensitivity evidence. See [scientific conventions](../docs/CONVENTIONS.md).

The [2026-10-08 mock-grid audit](audits/mock_grid_20261008/README.md)
separates the historical ray-tracing crop from the inherited fitting/plotting
grid. `grid_shape` now stores the adopted native mock grid; each record keeps
its previous inherited grid and evidence in `mock_grid_provenance`.
`central_deflection_radius_arcsec` is the adopted central source-suppression
radius, independent of `mask_radius_arcsec`. Preparation reads it by default;
historical external records lacking it require an explicit radius.

`mock_einstein_radius_arcsec` records the physical radius of the updated
SPH-template plus analytic halo mock, measured without central source
suppression. `mock_einstein_radius_arcsec_legacy` and particle-mass selection
radii remain distinct. Each updated truth radius links to its saved
convergence/provenance record. See [mock truth](../docs/MOCK_TRUTH.md).

Image positions are expressed in the native PyAutoLens fitting frame. The
production `flip_for_ds9: true` setting flips raw FITS arrays on input;
raw-FITS `origin='upper'` overlays must not be used to diagnose their y signs.
All 37 original position lists passed the corrected brightness check. See
[the coordinate audit](audits/positions_manual_20261008/README.md).
