"""Intrinsic stellar M/L in a circular aperture of the principal YZ projection."""

import numpy as np
from scipy.integrate import cumulative_simpson

from ._validation import finite_array, scalar
from .geometry import principal_axes, stellar_centre
from .photometry import absolute_magnitude_to_solar_luminosity
from .sph import _PROJECTED_KERNEL, _PROJECTED_RADIUS


_RADIAL_DENSITY = 16 * _PROJECTED_RADIUS * _PROJECTED_KERNEL
_CDF = cumulative_simpson(_RADIAL_DENSITY, x=_PROJECTED_RADIUS, initial=0)
_NORMALIZATION = float(_CDF[-1])
_CDF /= _NORMALIZATION


def projected_m4_aperture_weights(
    projected_distance, support_radius, aperture_radius, *, quadrature_order=64
):
    """Fraction of each LOS-projected spherical M4 kernel inside a circle.

    All distances share one physical or angular unit. Native supports are used
    without a pixel-dependent sampling floor. Kernels wholly inside or outside
    have exact weights 1 or 0; crossing kernels are integrated over radial rings
    using their analytic angular overlap. The same kernel table as sph_project
    is used, with its continuum normalization. This is an aperture integral,
    not a change to the foreground image renderer.
    """
    distance = finite_array(projected_distance, "projected_distance")
    support = finite_array(support_radius, "support_radius")
    radius = scalar(aperture_radius, "aperture_radius", positive=True)
    if distance.ndim != 1 or distance.shape != support.shape:
        raise ValueError("distances and supports must be matching 1D arrays")
    if np.any(distance < 0) or np.any(support <= 0):
        raise ValueError("distances must be non-negative and supports positive")
    if (
        isinstance(quadrature_order, bool)
        or not isinstance(quadrature_order, int)
        or quadrature_order < 16
    ):
        raise ValueError("quadrature_order must be an integer >= 16")
    weights = np.zeros_like(distance)
    inside = distance + support <= radius
    weights[inside] = 1
    crossing = np.flatnonzero(~inside & (distance - support < radius))
    nodes, gauss_weights = np.polynomial.legendre.leggauss(quadrature_order)
    t = (nodes + 1) / 2
    # This substitution regularizes the square-root angular-overlap endpoints.
    fraction = np.sin(np.pi * t / 2) ** 2
    derivative = np.pi / 2 * np.sin(np.pi * t)
    for start in range(0, len(crossing), 4096):
        indices = crossing[start : start + 4096]
        d, h = distance[indices], support[indices]
        r, offset = radius / h, d / h
        full_upper = np.clip(r - offset, 0, 1)
        values = np.interp(full_upper, _PROJECTED_RADIUS, _CDF)
        lower = np.minimum(np.abs(r - offset), 1)
        upper = np.minimum(r + offset, 1)
        span = np.maximum(upper - lower, 0)
        q = lower[:, None] + span[:, None] * fraction
        denominator = 2 * offset[:, None] * q
        cosine = np.divide(
            r[:, None] ** 2 - offset[:, None] ** 2 - q**2,
            denominator,
            out=np.ones_like(q),
            where=denominator > 0,
        )
        angular_fraction = 1 - np.arccos(np.clip(cosine, -1, 1)) / np.pi
        kernel = np.interp(q, _PROJECTED_RADIUS, _PROJECTED_KERNEL)
        partial = (
            np.sum(
                16 * q * kernel * angular_fraction * derivative * gauss_weights / 2,
                axis=1,
            )
            * span
            / _NORMALIZATION
        )
        weights[indices] = np.clip(values + partial, 0, 1)
    return weights


def stellar_aperture_mass_to_light(
    particles,
    radius_kpc,
    *,
    solar_absolute_magnitude=4.53,
    quadrature_order=64,
):
    """Return real-stellar mass/luminosity using identical SPH aperture weights.

    The caller supplies wind-free StellarParticles and a physical circular
    YZ radius. The full bound-subhalo LOS is retained, matching the historical
    aperture measurement and particle-mass selection. Dark matter, gas and
    wind masses never enter the stellar numerator. No PSF or dust is applied.
    """
    if particles.smoothing_lengths_kpc is None:
        raise ValueError("aperture M/L requires stored physical StellarHsml")
    radius = scalar(radius_kpc, "radius_kpc", positive=True)
    centre, iterations = stellar_centre(particles.positions_kpc, particles.masses_msun)
    centred = particles.positions_kpc - centre
    rotation, shape = principal_axes(
        centred, particles.masses_msun, 5 * particles.halfmass_radius_kpc
    )
    aligned = centred @ rotation.T
    distance = np.linalg.norm(aligned[:, 1:3], axis=1)
    weights = projected_m4_aperture_weights(
        distance,
        particles.smoothing_lengths_kpc,
        radius,
        quadrature_order=quadrature_order,
    )
    luminosity = absolute_magnitude_to_solar_luminosity(
        particles.absolute_magnitudes_i,
        solar_absolute_magnitude=solar_absolute_magnitude,
    )
    mass = float(np.dot(weights, particles.masses_msun))
    light = float(np.dot(weights, luminosity))
    if mass <= 0 or light <= 0:
        raise ValueError("stellar aperture mass and luminosity must be positive")
    point = distance < radius
    cube = np.all(np.abs(aligned) < 5 * particles.halfmass_radius_kpc, axis=1)
    cube_mass = float(np.dot(weights[cube], particles.masses_msun[cube]))
    cube_light = float(np.dot(weights[cube], luminosity[cube]))
    return {
        "radius_kpc": radius,
        "stellar_mass_msun": mass,
        "stellar_luminosity_lsun_i": light,
        "mass_to_light_msun_per_lsun": mass / light,
        "point_particle_mass_to_light": float(
            particles.masses_msun[point].sum() / luminosity[point].sum()
        ),
        "foreground_cube_mass_to_light": cube_mass / cube_light,
        "foreground_cube_mass_fraction": cube_mass / mass,
        "foreground_cube_light_fraction": cube_light / light,
        "solar_absolute_magnitude": solar_absolute_magnitude,
        "quadrature_order": quadrature_order,
        "centre_offset_kpc": centre.tolist(),
        "centre_iterations": iterations,
        "rotation_rows": rotation.tolist(),
        "shape": shape,
        "scope": "real stars; full bound-subhalo LOS; intrinsic SDSS i",
        "kernel": "continuum circular-aperture integral of LOS-projected 3D M4",
    }
