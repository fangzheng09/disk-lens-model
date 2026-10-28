# Production scientific conventions

The adopted stellar renderer is SPH with projected spherical M4 kernels,
TNG StellarHsml support and 8 subpixels. Centring, axis iteration, YZ rotation,
nominal cube and integer-pixel foreground aperture retain the tested geometry.
Historical mock crops are audited separately from fitting/plotting grids
and must fit within that unchanged foreground; preparation rejects oversized
crops. No production foreground padding or aperture expansion is used.
Particle coordinates and masses use physical kpc and Msun; extraction handles
periodic offsets and excludes wind particles. The snapshot scale factor is
independent of the adopted observational redshifts.

Stellar photometry uses intrinsic SDSS i magnitudes, Planck18 distance modulus
and the numerical F814W zero point 25.948. No dust, K-correction or filter
response conversion is introduced. The stellar template coefficient from the
supplied systems file is distinct from physical M/L and fitted MGE normalization;
see [photometry](PHOTOMETRY.md).

The mock retains lens/source redshifts 0.15/0.65, source intensity 0.4,
pixel-centre Gaussian PSF and noise derived from realised Poisson counts.
There is no target-S/N rescaling. New mock generation and truth evaluation use
PyAutoLens 2025.1.18.7; mock ray tracing explicitly keeps one sample per
native image pixel. This is separate from SPH subpixel integration and
from the fitting stages' own oversampling settings. The spherical Ludlow halo uses the inspected
historical PyAutoLens convention. Foreground and source ray-grid dimensions
are separate. The central pseudo-core radius and fitting mask are distinct,
per-system inputs. All34 adopted members use the user-confirmed 0.2 arcsec
central source suppression, with an explicit override available.

The five-stage SLaM model, priors, meshes and image-position constraints are
preserved. New main-run live-point defaults are documented in [fitting](FITTING.md);
historical saved configurations determine the settings of existing runs.
The Pix2 diagnostic guard skips invalid intermediate
plot/runtime diagnostics, while likelihood and final-result exceptions retain
their original behavior. Final parameter summaries use the recorded posterior
weights. Enclosed mass and other derived paper quantities are outside this
production workflow.

## Adopted kinematic selection

For production sample membership, retain the historical catalogue frame:
positions are relative to `SubhaloPos`, and particle velocities are measured
relative to `SubhaloVel`. Convert coordinates to physical kpc with `a/h` and
snapshot particle velocities to peculiar km/s with `sqrt(a)`; catalogue
`SubhaloVel` is already in km/s. Use the physical 30 kpc spherical aperture,
the existing mass-weighted statistic `sum(m_i * v_phi_i^2 / v_i^2) / sum(m_i)`
and the strict threshold `kappa_rot > 0.45`.

The adopted catalogue contains 34 systems, with20 included in the main
recovery statistics. The runner consumes these fixed inputs and does not
automatically reselect galaxies. The iterative centre used for the SPH
foreground is a separate rendering operation.

## Particle-mass selection Einstein radii

The adopted reconstruction uses the LOS-projected spherical M4 kernel,
with separate stored stellar and dark-matter supports and cell-volume-derived
gas supports. All components share physical units, periodic offsets, the
iterative rendering centre, principal axes and lens plane. The historical
implementation labelled CIC was nearest-grid-point deposition.

Use the effective radius `sqrt(A/pi)` of the largest closed tangential
critical curve. Record additional curves, the historical sum of curve radii
and a compact BH dynamical-mass correction separately. Selection retains
its historical three-component domain (stars plus winds, DM and gas),
the three-view maximum and strict `R_max > 0.5 arcsec` cut. Each object
requires convergence of its three-view maximum, enlarged-field convergence
and a stable threshold decision. The stricter individual-view grid flags
remain separate, including unresolved tiny weak-view curves.
The mass map is a bound-subhalo reconstruction, without unrelated LOS
structures or parent-halo fuzz; gas smoothing is an approximation to AREPO
cells. These numerical checks do not remove physical resolution limitations.

This selection radius is distinct from the production mock's SPH stellar
light template with its adopted coefficient and analytic Ludlow NFW halo.
The mock does not project TNG gas or dark-matter particles. Inherited
`_legacy` mock radii, central-deflection radii and fitting masks remain
separate per-system inputs. No artificial central deflection is used for
the particle-mass selection radius.

See [method, component scales and invocation](EINSTEIN_RADIUS.md). Numerical
screening does not automatically change the adopted 34 members or their
fit classifications.

## Adopted stellar aperture M/L

The 2026-10-08 update measures real-stellar mass and intrinsic SDSS i light
with identical native StellarHsml projected M4 aperture weights, at each
updated YZ selection radius. Solar magnitude 4.53 is used consistently for
physical M/L and conversion to the counts/s/pixel mock coefficient. New plans
prefer the recalculated coefficient, with historical-file fallback. Old
values and existing outputs remain preserved. Other adopted inputs retain
their historical provenance; run-specific overrides are recorded separately. See [stellar M/L](STELLAR_ML.md).

## Physical mock-model truth

The mock stage rebuilds the same native stellar deflections at central
radius zero and evaluates the actual stellar-template plus Ludlow-halo
tracer. It saves the largest closed critical-curve effective radius and
records all loop areas and the historical sum-of-radii diagnostic. Contour
sampling convergence is checked within the native interpolation domain.
The source-suppression image retains the user-confirmed 0.2 arcsec radius.
See [mock truth](MOCK_TRUTH.md).
