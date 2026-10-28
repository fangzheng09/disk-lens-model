# Adopted stellar SPH rendering

The foreground uses the line-of-sight projection of a spherical 3D M4
cubic-spline kernel, integrated over image subpixels. Its support radius is
TNG StellarHsml in physical kpc, representing approximately 32 same-type
neighbours. A missing input support array is computed from the distance to
the 32nd other stellar particle in the full input before aperture selection.

The production settings use 32 neighbours and 8 subpixels. The historical
centre iteration, reduced-inertia axes, YZ rotation and nominal cube selection
match the adopted SPH rendering. The grid spans the integer-pixel portion of
the 10 R_half aperture. SPH is the particle renderer; it is distinct from
PSF convolution and from the mock-observation noise prescription.

Small unresolved kernels use the existing numerical sampling floor. Wide
kernels retain the existing analytic normalization. Field-edge losses are
recorded rather than renormalized into the image:

```text
selected flux = deposited flux + edge kernel flux + outside-field flux
```

`templates/lightmap.template.json` is the SPH foreground configuration.
Manifests record support source, kernel settings, geometry and flux accounting.

The production ray grid is the audited historical mock crop and must fit
inside the nominal foreground. `prepare-main` rejects oversized crops rather
than padding the image. The full SPH image contributes to stellar deflections
before the ray grid is cropped; no foreground aperture change is introduced.
The adopted grids are included in `catalog/production_systems.json`.
