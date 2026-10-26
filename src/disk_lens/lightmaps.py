"""Render the adopted stellar SPH light map in the validated YZ geometry."""

import numpy as np
from astropy import units as u

from ._validation import scalar
from .geometry import principal_axes, stellar_centre
from .grid import ImageGrid
from .particles import StellarParticles
from .photometry import absolute_magnitude_to_count_rate

SPH_SETTINGS = {"method": "sph", "neighbors": 32, "workers": 1, "subpixels": 8}


def aperture_grid(
    halfmass_radius_kpc,
    *,
    pixel_scale_arcsec,
    lens_redshift,
    cosmology,
    aperture_halfmass_multiple=5.0,
):
    """Use the tested floor(10 R_half / angular pixel width) foreground size."""
    radius = scalar(halfmass_radius_kpc, "halfmass_radius_kpc", positive=True)
    multiplier = scalar(
        aperture_halfmass_multiple, "aperture_halfmass_multiple", positive=True
    )
    pixel = scalar(pixel_scale_arcsec, "pixel_scale_arcsec", positive=True)
    scale = cosmology.kpc_proper_per_arcmin(lens_redshift).to_value(u.kpc / u.arcsec)
    size = int(np.floor(2 * multiplier * radius / scale / pixel))
    return ImageGrid((size, size), pixel)


def build_lightmap(
    particles: StellarParticles,
    grid: ImageGrid,
    *,
    lens_redshift,
    zeropoint,
    cosmology,
    view="yz",
    aperture_halfmass_multiple=5.0,
    deposition=None,
):
    """Preserve the centring, principal-frame rotation and cube from the SPH run."""
    from .sph import neighbour_support_radii, sph_project

    if view != "yz":
        raise ValueError("the production SPH workflow uses the validated yz view")
    if grid.shape[0] != grid.shape[1]:
        raise ValueError("the production SPH foreground grid must be square")
    deposition = dict(SPH_SETTINGS) if deposition is None else dict(deposition)
    if set(deposition) != set(SPH_SETTINGS) or deposition["method"] != "sph":
        raise ValueError(
            "production deposition requires sph, neighbors, workers and subpixels"
        )
    if deposition["neighbors"] != 32:
        raise ValueError("the adopted SPH support uses 32 neighbours")
    if (
        isinstance(deposition["workers"], bool)
        or not isinstance(deposition["workers"], int)
        or deposition["workers"] < 1
    ):
        raise ValueError("SPH workers must be a positive integer")
    multiplier = scalar(
        aperture_halfmass_multiple, "aperture_halfmass_multiple", positive=True
    )
    centre, iterations = stellar_centre(particles.positions_kpc, particles.masses_msun)
    centred = particles.positions_kpc - centre
    half_width = multiplier * particles.halfmass_radius_kpc
    rotation, shape_info = principal_axes(centred, particles.masses_msun, half_width)
    aligned = centred @ rotation.T
    rotated = aligned @ np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]]).T
    scale = cosmology.kpc_proper_per_arcmin(lens_redshift).to_value(u.kpc / u.arcsec)
    angular = rotated / scale
    half_fov = grid.shape[0] * grid.pixel_scale_arcsec / 2
    selected = np.all(np.abs(rotated) < half_width, axis=1)
    selected &= np.all((angular >= -half_fov) & (angular <= half_fov), axis=1)
    rates = absolute_magnitude_to_count_rate(
        particles.absolute_magnitudes_i,
        redshift=lens_redshift,
        zeropoint=zeropoint,
        cosmology=cosmology,
    )
    if particles.smoothing_lengths_kpc is None:
        support = neighbour_support_radii(
            particles.positions_kpc, neighbours=32, workers=deposition["workers"]
        )
        support_source = "recomputed exact 32nd other stellar neighbour in 3D"
    else:
        support = particles.smoothing_lengths_kpc
        support_source = particles.metadata.get(
            "smoothing", "stored physical StellarHsml"
        )
    image, flux = sph_project(
        angular[selected, :2],
        rates[selected],
        support[selected] / scale,
        grid,
        subpixels=deposition["subpixels"],
    )
    flux.update(
        {
            "all_particle_flux": float(rates.sum()),
            "outside_aperture_flux": float(rates[~selected].sum()),
            "support_radius_source": support_source,
        }
    )
    metadata = {
        "particle_input": particles.metadata,
        "view": "yz",
        "geometry": "validated paper-era centre, reduced inertia axes and YZ rotation",
        "principal_axes": "X=major,Y=intermediate,Z=minor",
        "centre_offset_kpc": centre.tolist(),
        "rotation_rows": rotation.tolist(),
        "centre": {"iterations": iterations},
        "shape": shape_info,
        "flux": flux,
        "aperture_halfmass_multiple": multiplier,
        "lens_redshift": lens_redshift,
        "zeropoint": zeropoint,
        "cosmology": cosmology.name,
        "deposition": "line-of-sight projected spherical M4 kernel; no PSF",
        "deposition_parameters": deposition,
        "intrinsic_light": True,
        "psf_convolved": False,
        "calibration": "SDSS i plus distance modulus and supplied zero point; no K-correction",
    }
    return image, metadata
