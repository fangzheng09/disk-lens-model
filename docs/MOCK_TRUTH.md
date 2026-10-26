# Physical Einstein radii of the production mock model

The mock truth is computed from the same unblurred SPH stellar light template,
adopted stellar coefficient, native deflection grid, bilinear interpolation
and analytic spherical Ludlow NFW halo used to generate each mock.
New mock generation and truth calculations require PyAutoLens 2025.1.18.7
and lenstronomy 1.13.2 in the same environment as fitting. The historical
Ludlow implementation retains its Planck15 convention; photometry remains
Planck18. No particle-mass selection radius is substituted.

Stellar deflections are rebuilt from the original light image with the
central-deflection radius set to **zero**. Source-suppression radius 0.2
arcsec remains in the generated image. Source brightness, PSF, Poisson
noise and the observational fitting mask do not enter the truth.

The primary radius is `sqrt(A/pi)` for the largest closed tangential
critical curve of the exact mock tracer. The area/radius of every loop and
the historical PyAutoLens sum of individual radii are recorded separately.
Use the same definition for truth and posterior-derived comparisons.

The tangential eigenvalue is computed through the tracer's deflection
Jacobian, not the placeholder zero convergence in its input-deflection
profile. Evaluation stays inside the native interpolation domain. Sampling
at native pixel width, half and quarter width is checked; an eighth width
is added if needed. The last radius difference must be no greater than
`max(0.005 arcsec, 0.01 * radius)`. This checks contour sampling of the fixed
discrete mock model; it does not improve the physical SPH resolution. Open
curves or a critical region on the evaluation boundary are rejected.

`run-main --through mock` now saves `truth/` after the imaging bundle;
fit/summary stages require that truth bundle to validate against their
mock configuration and input hashes. If an older mock lacks truth, compute
it in a new unified-version plan/output directory. Existing immutable plans
without the truth-path field use a sibling `truth/` directory. Historical
bundles retain their recorded software versions and are not relabelled as
new 2025 products.

```bash
/path/to/python -m disk_lens mock-truth \
  --config /work/run_3/config/mock.json \
  --mock /work/run_3/mock --output /work/run_3/truth
```

The truth bundle contains `truth.json`, unmodified unit-stellar deflections
and a manifest with model/input/file hashes. These files remain separate
from the 0.2 arcsec-modified deflections used for the source image. The
production catalogue's `mock_einstein_radius_arcsec` is a synchronized
value from that truth record. Historical `_legacy` radii are preserved.
Prepared catalogue snapshots are not mutated during a run.

The unified version is validated against all 37 existing mocks and saved
unmasked radii; see [version compatibility](VERSION_COMPATIBILITY.md).
