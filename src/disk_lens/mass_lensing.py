"""Mass projection and isolated single-plane tangential critical curves.

The renderer uses the production LOS-projected spherical M4 kernel, not
GLAMER's 2014 truncated Gaussian or its newer two-dimensional B-spline.
GLAMER motivates three-dimensional, separately defined component supports.
"""

import numpy as np
from scipy.signal import fftconvolve
from skimage.measure import find_contours

from disk_lens.grid import ImageGrid
from disk_lens.sph import sph_project


def projected_mass(positions_yx, masses, supports, grid: ImageGrid, *, subpixels=2):
    image, info = sph_project(positions_yx, masses, supports, grid, subpixels=subpixels)
    accounting = {
        key.replace("flux", "mass_msun"): value for key, value in info.items()
    }
    total = accounting["input_mass_msun"]
    residual = total - (
        accounting["deposited_mass_msun"]
        + accounting["outside_fov_mass_msun"]
        + accounting["edge_kernel_mass_msun"]
    )
    accounting["accounting_relative_error"] = abs(residual) / total if total else 0.0
    if accounting["accounting_relative_error"] > 1e-9:
        raise ValueError("Mass accounting failed")
    return image, accounting


def point_mass_pixels(positions_yx, masses, grid):
    """Pixel-integrated point masses via four-centre linear assignment.

    Used only for BH dynamical particle masses with no resolved spatial
    density support. This is a numerical representation of compact mass,
    not a BH_Hsml accretion-kernel model.
    """
    image = np.zeros(grid.shape)
    row = (grid.shape[0] - 1) / 2 - positions_yx[:, 0] / grid.pixel_scale_arcsec
    col = (grid.shape[1] - 1) / 2 + positions_yx[:, 1] / grid.pixel_scale_arcsec
    r0, c0 = np.floor(row).astype(int), np.floor(col).astype(int)
    for dr in [0, 1]:
        for dc in [0, 1]:
            r, c = r0 + dr, c0 + dc
            weights = (1 - abs(row - r)) * (1 - abs(col - c))
            valid = (r >= 0) & (r < grid.shape[0]) & (c >= 0) & (c < grid.shape[1])
            np.add.at(image, (r[valid], c[valid]), masses[valid] * weights[valid])
    return image


def shear_from_kappa(kappa, grid: ImageGrid):
    """Isolated linear convolution with the point-mass shear Green functions.

    Local convergence is the input mass density. Green functions are
    integrated analytically over each square source pixel; its self shear
    is zero by symmetry. Unlike a periodic FFT Poisson
    solve, this linear convolution introduces no periodic images or
    compensating mean-density sheet. Spatial convergence remains required.
    """
    kappa = np.asarray(kappa, dtype=float)
    if kappa.shape != grid.shape or not np.all(np.isfinite(kappa)) or np.any(kappa < 0):
        raise ValueError("Invalid convergence map")
    ny, nx = grid.shape
    p = grid.pixel_scale_arcsec
    dy, dx = np.meshgrid(
        -np.arange(-ny + 1, ny) * p,
        np.arange(-nx + 1, nx) * p,
        indexing="ij",
    )
    x0, x1, y0, y1 = dx - p / 2, dx + p / 2, dy - p / 2, dy + p / 2
    kernel1 = (
        np.arctan(y1 / x1)
        - np.arctan(y0 / x1)
        - np.arctan(y1 / x0)
        + np.arctan(y0 / x0)
    ) / np.pi
    kernel2 = (
        0.5
        / np.pi
        * np.log(
            ((x1 * x1 + y1 * y1) * (x0 * x0 + y0 * y0))
            / ((x0 * x0 + y1 * y1) * (x1 * x1 + y0 * y0))
        )
    )
    kernel1[ny - 1, nx - 1] = 0
    kernel2[ny - 1, nx - 1] = 0
    gamma1 = fftconvolve(kappa, kernel1, mode="same")
    gamma2 = fftconvolve(kappa, kernel2, mode="same")
    return gamma1, gamma2


def critical_radius(kappa, grid):
    """Largest closed tangential critical curve, radius sqrt(area/pi).

    Additional curves and the historical sum-of-radii convention are
    reported separately. A component reaching the map edge is invalid.
    Absence of a critical curve is a valid zero radius, not an exception.
    No source-position mask or artificial central deflection is used.
    """
    gamma1, gamma2 = shear_from_kappa(kappa, grid)
    eigenvalue = 1 - kappa - np.hypot(gamma1, gamma2)
    if np.any(eigenvalue[[0, -1], :] <= 0) or np.any(eigenvalue[:, [0, -1]] <= 0):
        raise ValueError("Tangential critical region reaches map boundary")
    areas = []
    for curve in find_contours(eigenvalue, 0, fully_connected="low"):
        if not np.allclose(curve[0], curve[-1], atol=1e-8):
            raise ValueError("Open tangential critical curve")
        y, x = curve[:, 0], curve[:, 1]
        area = abs(np.dot(x[:-1], y[1:]) - np.dot(y[:-1], x[1:])) / 2
        # Report subpixel loops too; do not silently add them to the main radius.
        if area > 0:
            areas.append(float(area * grid.pixel_scale_arcsec**2))
    areas.sort(reverse=True)
    radii = np.sqrt(np.array(areas) / np.pi)
    return {
        "einstein_radius_arcsec": float(radii[0]) if len(radii) else 0.0,
        "critical_curve_areas_arcsec2": areas,
        "critical_curve_count": len(areas),
        "historical_sum_of_radii_arcsec": float(radii.sum()),
        "second_to_first_area": areas[1] / areas[0] if len(areas) > 1 else 0.0,
        "boundary_clear": True,
    }
