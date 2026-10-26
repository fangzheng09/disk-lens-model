# Human-reviewed image positions

Each per-system record stores an explicit user selection, previous and adopted
(y,x) coordinates in arcseconds, and arc image / detector provenance hashes.
Candidates use the noiseless PSF-convolved lensed source component without
foreground light. Existing frozen fit inputs and mock images are unchanged.
Only confirmed systems are included; detection alone is not confirmation.

## FITS-coordinate correction

The production AutoArray setting `flip_for_ds9: true` flips raw FITS pixels
on loading. Schema-1 raw-array position checks omitted that conversion and
misclassified 15 systems. The detector now loads arrays exactly as fitting
does. All 37 original position lists pass the corrected brightness checks.
`positions_fit_frame_37.csv` and `fits_frame_correction.json` record the
corrected check. The user's C1/C2 choice for 83999 keeps the same two image
pixels; their y coordinates were converted into the native fitting frame.
No mock files, truth radii or active fitting inputs were changed.
