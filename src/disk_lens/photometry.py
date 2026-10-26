"""Magnitude, count-rate and stellar-convergence conversions.

Cosmology, band conventions and pixel scale are explicit inputs. No dust,
K-correction or bandpass transformation is inferred. Counts must use the same
instrument convention as the supplied magnitude zero point. These functions
do not define a universal PyAutoLens mass-to-light parameter: the returned
factor applies to a numerical counts/s/pixel template at the specified scale.
See docs/PHOTOMETRY.md for the derivation and limitations.
"""

import numpy as np
from astropy import constants as const
from astropy import units as u
from astropy.cosmology import FLRW

from ._validation import finite_array, scalar


def absolute_magnitude_to_count_rate(
    absolute_magnitude,
    *,
    redshift: float,
    zeropoint: float,
    cosmology: FLRW,
) -> np.ndarray:
    """Return counts/s from absolute magnitudes and luminosity distance.

    Uses m = M + distance modulus and rate = 10**(0.4 * (ZP - m)).
    The distance modulus is the only correction applied. For a particle array,
    each output is the integrated count rate of that particle, not a surface
    brightness; deposition into pixels is a separate operation.
    """
    z = scalar(redshift, "redshift", positive=True)
    zp = scalar(zeropoint, "zeropoint")
    magnitudes = finite_array(absolute_magnitude, "absolute_magnitude")
    distance_modulus = cosmology.distmod(z).to_value(u.mag)
    with np.errstate(over="raise", invalid="raise"):
        return np.power(10.0, 0.4 * (zp - magnitudes - distance_modulus))


def absolute_magnitude_to_solar_luminosity(
    absolute_magnitude, *, solar_absolute_magnitude: float
) -> np.ndarray:
    """Return luminosity in solar units in the SAME band/magnitude convention."""
    solar_mag = scalar(solar_absolute_magnitude, "solar_absolute_magnitude")
    magnitudes = finite_array(absolute_magnitude, "absolute_magnitude")
    with np.errstate(over="raise", invalid="raise"):
        return np.power(10.0, 0.4 * (solar_mag - magnitudes))


def critical_surface_density_msun_kpc2(
    *, lens_redshift: float, source_redshift: float, cosmology: FLRW
) -> float:
    """Return physical Sigma_crit for a single lens/source plane in Msun/kpc²."""
    zl = scalar(lens_redshift, "lens_redshift", positive=True)
    zs = scalar(source_redshift, "source_redshift", positive=True)
    if zs <= zl:
        raise ValueError("source_redshift must exceed lens_redshift")
    dl = cosmology.angular_diameter_distance(zl)
    ds = cosmology.angular_diameter_distance(zs)
    dls = cosmology.angular_diameter_distance_z1z2(zl, zs)
    sigma = const.c**2 * ds / (4 * np.pi * const.G * dl * dls)
    return scalar(sigma.to_value(u.Msun / u.kpc**2), "Sigma_crit", positive=True)


def physical_pixel_area_kpc2(
    *, pixel_scale_arcsec: float, redshift: float, cosmology: FLRW
) -> float:
    """Return physical area of a square angular pixel at the given redshift."""
    pixel_scale = scalar(pixel_scale_arcsec, "pixel_scale_arcsec", positive=True)
    z = scalar(redshift, "redshift", positive=True)
    theta = (pixel_scale * u.arcsec).to_value(u.rad)
    side_kpc = cosmology.angular_diameter_distance(z).to_value(u.kpc) * theta
    return scalar(side_kpc**2, "pixel area", positive=True)


def stellar_kappa_per_count_rate(
    *,
    mass_to_light_msun_per_lsun: float,
    solar_absolute_magnitude: float,
    zeropoint: float,
    pixel_scale_arcsec: float,
    lens_redshift: float,
    source_redshift: float,
    cosmology: FLRW,
) -> float:
    """Return C such that kappa = C * light_map_counts_per_second_per_pixel.

    C = (M/L) / (solar_count_rate * pixel_area * Sigma_crit).
    M/L is a global, non-negative physical value in the specified solar band.
    Its solar-magnitude convention must match that used to derive M/L itself.
    This factor cannot be applied unchanged to counts/s/arcsec² or total counts.
    """
    ml = scalar(mass_to_light_msun_per_lsun, "mass_to_light_msun_per_lsun")
    if ml < 0:
        raise ValueError("mass_to_light_msun_per_lsun must be non-negative")
    solar_mag = scalar(solar_absolute_magnitude, "solar_absolute_magnitude")
    sigma = critical_surface_density_msun_kpc2(
        lens_redshift=lens_redshift,
        source_redshift=source_redshift,
        cosmology=cosmology,
    )
    area = physical_pixel_area_kpc2(
        pixel_scale_arcsec=pixel_scale_arcsec,
        redshift=lens_redshift,
        cosmology=cosmology,
    )
    solar_rate = scalar(
        absolute_magnitude_to_count_rate(
            solar_mag, redshift=lens_redshift, zeropoint=zeropoint, cosmology=cosmology
        ),
        "solar count rate",
        positive=True,
    )
    return ml / (solar_rate * area * sigma)


def physical_mass_to_light_from_kappa_per_count_rate(
    kappa_per_count_rate: float,
    *,
    solar_absolute_magnitude: float,
    zeropoint: float,
    pixel_scale_arcsec: float,
    lens_redshift: float,
    source_redshift: float,
    cosmology: FLRW,
) -> float:
    """Convert the legacy/instrumental light normalization to physical M/L.

    ``kappa_per_count_rate`` must multiply a counts/s/pixel light template to
    produce convergence. This is the exact inverse of
    :func:`stellar_kappa_per_count_rate` for the same band, grid and cosmology.
    """
    normalization = scalar(kappa_per_count_rate, "kappa_per_count_rate")
    if normalization < 0:
        raise ValueError("kappa_per_count_rate must be non-negative")
    unit_normalization = stellar_kappa_per_count_rate(
        mass_to_light_msun_per_lsun=1.0,
        solar_absolute_magnitude=solar_absolute_magnitude,
        zeropoint=zeropoint,
        pixel_scale_arcsec=pixel_scale_arcsec,
        lens_redshift=lens_redshift,
        source_redshift=source_redshift,
        cosmology=cosmology,
    )
    return normalization / unit_normalization
