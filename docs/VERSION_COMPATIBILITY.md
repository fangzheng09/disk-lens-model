# Unified AutoLens version and mock compatibility

New mock generation, physical mock truth and fitting use **AutoLens
2025.1.18.7** in the same Python environment, with lenstronomy 1.13.2.
The original 37 SPH mock images were generated with AutoLens 2024.11.13.2;
their saved manifests retain that actual provenance.

## Preserving the original mock calculation

The old 2024 grid default used native-pixel ray tracing. The 2025 grid
default uses 4×4 oversampling. New mock and truth code explicitly sets
`over_sample_size=1`, preserving the original numerical prescription.
This is independent of SPH subpixel integration and of the fitting stages'
own oversampling settings. The mass model, source, PSF prescription,
Poisson seed and noise definition are retained.

A direct switch to the 2025 default grid would alter source images. The
explicit native sampling is therefore part of the version migration.

## Numerical comparison

On 26 October 2026, all 37 production systems were regenerated from their
existing immutable SPH lightmaps and mock configurations using the unified
2025 code. The comparison included the current NumPy 2.1.3/SciPy 1.16.0
numerical stack and lenstronomy 1.13.2; the reference mock stack used NumPy
1.26.4/SciPy 1.14.0. Both used the same physical inputs and random seeds.
The verification required no new fitting or particle projection.

| Quantity | Comparison over all 37 systems |
|---|---|
| Observed data | Every pixel exactly equal in all 37 systems |
| Realised noise map | Every pixel exactly equal in all 37 systems |
| Stellar deflections, central source suppression included | Every element exactly equal in all 37 systems |
| Unmodified stellar deflections for truth | Every element exactly equal in all 37 systems |
| Gaussian PSF | Maximum absolute difference 6.94e-18 |
| Model and convolved source images | Maximum absolute difference 2.49e-13 counts/s/pixel |
| Unconvolved source image | Maximum absolute difference 6.83e-13 counts/s/pixel |
| Convolved foreground | Maximum absolute difference 4.98e-14 counts/s/pixel |
| Unmasked mock Einstein radius | Maximum absolute difference 5.29e-13 arcsec |

The PSF is not bitwise identical: NumPy's `exp` changes last-place rounding
in 46 of its 225 pixels. Coordinates, exponents and the unnormalized sum
match in the isolated comparison. The PSF's largest difference relative to
its peak is 4.36e-17, below one float64 machine epsilon; normalization differs
by one machine epsilon. The strict bitwise-PSF criterion remains recorded
as unmet. This comparison establishes numerical equivalence of the
preserved mock prescription, rather than claiming every output file has
identical bytes.

Detailed per-system differences, source/input checksums, actual import
paths, PSF rounding diagnostics and scheduler records remain in the
research validation archive. A subsequent source-formatting cleanup was
verified to leave the Python AST unchanged. Historical input files and
results were not overwritten.

## Plans and existing products

New plans record their required package versions and native mock sampling;
reusing a bundle from another version or sampling prescription is rejected.
Historical schema-2 plans without this version field may reuse their complete,
verified 2024 mock/truth products. Missing historical products require a new
plan and output directory; the runner does not silently regenerate them
in place. A fit that reads an existing 2024 observation retains that
observation's recorded provenance.
