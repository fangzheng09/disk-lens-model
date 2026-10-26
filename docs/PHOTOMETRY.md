# Photometry and stellar convergence

Implemented in `disk_lens.photometry`. All scientific calibration inputs are
explicit. No solar magnitude, zero point, redshift or pixel scale is hard-coded.

## Definitions

For an absolute magnitude M and a supplied cosmology:

```text
m = M + 5 log10(D_L / 10 pc)
F = 10^[0.4 (ZP - m)]                          [counts/s]
L / L_sun = 10^[0.4 (M_sun - M)]
A_pix = [D_A(z_l) * pixel_angle_in_radians]^2  [physical kpc²]
Sigma_crit = c²/(4 pi G) * D_s/(D_l D_ls)      [Msun/kpc²]
```

Let F_sun be the count rate of one solar luminosity at the lens distance, in
the same band/convention as the physical M/L. If I_pix contains counts/s/pixel:

```text
L_pix / L_sun = I_pix / F_sun
M_pix = (M/L) * I_pix / F_sun
kappa_pix = M_pix / (A_pix * Sigma_crit)
          = I_pix * stellar_kappa_per_count_rate(...)
```

The inverse pixel-area factor in the original conversion notebook is therefore
appropriate **for this template convention**. It is not a universal conversion
for every light profile or a counts/s/arcsec² image. Applying this factor to a
PyAutoLens analytic MGE requires confirming the units of that fitted intensity;
that interface has not yet been validated.

## Historical coefficient check

For Planck18, zl=0.15, zs=0.65, pixel scale 0.05 arcsec, ZP=25.948,
M_sun=4.52 and M/L=1:

| Quantity | Value |
|---|---:|
| Solar count rate | 6.860562596652988e-8 counts/s |
| Physical Sigma_crit | 4.054264102341438e9 Msun/kpc² |
| New factor using exact angle units | 0.19705596378056342 |
| Saved notebook factor | 0.1970563339860656 |

The approximately 1.9 ppm relative difference is accounted for by the notebook's
rounded 206265 arcsec/radian, rather than a reversed pixel-area factor.

Using M_sun=4.53 instead changes the factor at fixed numerical M/L to
0.198879300164156. But changing the solar convention **consistently** also changes
the numerical M/L inferred from a fixed physical mass and absolute magnitudes;
these effects cancel. A mismatch between the conventions in separate steps does
not cancel. The production stellar-aperture update now adopts 4.53 consistently for
both steps; the generic conversion functions still require an explicit
solar magnitude. See [stellar M/L](STELLAR_ML.md).

## Bandpass scope

The input notebook takes column 6 of `GFM_StellarPhotometrics` (SDSS i), adds a
distance modulus, and uses the ACS/F814W numerical zero point 25.948. This is a
specified HST-like mock prescription, not a synthetic transformation of the
stellar spectrum into redshifted F814W. No K-correction, extinction or spectral
response integration is implemented here. Keep this approximation explicit.

References: [TNG photometric fields](https://www.tng-project.org/data/docs/specifications/),
[Astropy cosmology distances](https://docs.astropy.org/en/stable/cosmology/index.html).

## Example and verification

After installing this project into an environment with its dependencies:

```python
from astropy.cosmology import Planck18
from disk_lens.photometry import stellar_kappa_per_count_rate

factor = stellar_kappa_per_count_rate(
    mass_to_light_msun_per_lsun=1.0,
    solar_absolute_magnitude=4.52,
    zeropoint=25.948,
    pixel_scale_arcsec=0.05,
    lens_redshift=0.15,
    source_redshift=0.65,
    cosmology=Planck18,
)
# kappa = light_map_counts_per_second_per_pixel * factor
```

The fitted MGE normalization is distinct from the SPH light-template coefficient.
Parameter summaries retain the original model units and do not automatically
convert fitted MGE amplitudes to physical mass-to-light ratios.
