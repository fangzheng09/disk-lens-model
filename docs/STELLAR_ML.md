# Stellar aperture M/L and the mock coefficient

The35 adopted records use each object's YZ particle-mass selection Einstein
radius for the input stellar aperture measurement. The three-view maximum
is the selection statistic; it does not replace that YZ aperture.

## Measurement

Only real stars (`GFM_StellarFormationTime > 0`) enter either stellar mass
or luminosity. Current particle masses are in Msun and intrinsic SDSS i
luminosities use `10**[0.4*(4.53 - M_i)]`. Wind, gas and dark-matter masses
contribute to the radius calculation but not to stellar M/L.

The centre and principal axes use the same real-star geometry as the
foreground and mass audit. The aperture is a centred projected YZ circle,
with physical radius `theta_E_YZ * kpc_per_arcsec(z_l=0.15)`, retaining the
full bound-subhalo LOS as in the historical aperture measurement.

`disk_lens.stellar_ml` integrates each native StellarHsml projected M4 kernel
over that circle. Identical kernel fractions weight mass and luminosity:

```text
Upsilon_i = sum(w_j * m_j) / sum(w_j * L_i,j)
```

The radial kernel table is shared with the SPH foreground renderer. Angular
overlap of a kernel-centred ring with the aperture is integrated by Gaussian
quadrature, including kernels whose centres lie outside the circle. This
continuum aperture integration has no image pixels or sampling floor; it
does not alter the foreground's 8-by-8 image sampling. The published values
use quadrature order 128, checked against 64 for every production object.
Independent integration of the underlying 3D M4 kernel and a Cartesian
renderer check are recorded in the audit summary.

The full bound-subhalo aperture is distinct from the foreground's nominal
5 R_half cube selection. The measurements also record the cube-only M/L
as a diagnostic; its difference from the adopted aperture value is below
1% for these 35 systems. Neither an aperture average nor a constant mock
coefficient represents the galaxy's full spatial M/L variation.

## Calibration and parameter use

Solar i magnitude **4.53** is now used consistently in both luminosities
and the counts/s/pixel-to-convergence conversion. At Planck18, zl=0.15,
zs=0.65, ZP=25.948 and pixel scale 0.05 arcsec:

```text
n = Upsilon_i * 0.19887930016415595
```

Each `production_systems.json` record now includes:

| Field | Meaning |
|---|---|
| `selection_einstein_radius_yz_arcsec` | updated YZ aperture radius |
| `selection_einstein_radius_max_arcsec` | updated screening maximum |
| `stellar_mass_to_light_msun_per_lsun_i` | adopted real-stellar aperture M/L |
| `stellar_mass_normalization` | adopted mock coefficient n |
| `stellar_aperture` | aperture masses, luminosities, geometry and calibration |
| `stellar_mass_normalization_legacy` | unchanged historical n |
| `mock_einstein_radius_arcsec_legacy` | unchanged historical mock-model radius |

`prepare-main` prefers `stellar_mass_normalization`; a historical systems
file lacking this field falls back to `stellar_mass_normalization_legacy`.
An explicit `--stellar-normalization` overrides either. The chosen source
is recorded in `plan.json`. Existing prepared plans and completed mocks or
fits are preserved; new parameters require preparing a new run directory.
Masks, source centres, image-position constraints and halo masses remain
inherited per-system inputs. The selection radius does not replace the
mock-model radius, which must be measured again after generating a new mock.

## Recompute one aperture

Use an extracted wind-free stellar NPZ and the chosen YZ selection radius:

```bash
disk-lens stellar-ml --particles /work/inputs/stars_3.npz \
  --radius-arcsec 0.9694408726417415 --output /work/ml_3.json
```

The CLI records explicit calibration inputs, geometry and the particle-file
checksum. It refuses to replace an existing output. Run particle calculations
on an allocated compute node. Other band, pixel or redshift conventions
must be supplied consistently through the command options.
