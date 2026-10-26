"""Adaptive stellar-light rendering with a projected 3D cubic-spline kernel.

The neighbour distance is measured in 3D and becomes the support radius of
a spherical M4 kernel. The kernel is integrated along the line of sight
before it is sampled on the image plane. This is an SPH-like rendering
choice, not a claim that collisionless stars were evolved with SPH.
"""

import math

import numpy as np
from scipy.spatial import cKDTree

try:
    from numba import njit
except ImportError:
    _NUMBA = False

    def njit(function=None, **options):
        """Accept both bare and configured Numba decorator syntax."""
        return (lambda function: function) if function is None else function
else:
    _NUMBA = True

from ._validation import finite_array
from .grid import ImageGrid


def neighbour_support_radii(
    positions_kpc, *, neighbours=32, workers=1, chunk_size=100_000
):
    """3D distance to the requested number of *other* stellar particles.

    ``cKDTree.query`` includes each query particle itself at distance zero,
    so the desired column is obtained with ``k=neighbours + 1``. The query
    is chunked to avoid an N-by-neighbours distance allocation.
    """
    positions = finite_array(positions_kpc, "positions_kpc")
    if positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError("positions_kpc must have shape (N, 3)")
    if (
        isinstance(neighbours, bool)
        or not isinstance(neighbours, int)
        or not 1 <= neighbours < len(positions)
    ):
        raise ValueError("neighbours must be an integer between 1 and N-1")
    if isinstance(workers, bool) or not isinstance(workers, int) or workers < 1:
        raise ValueError("workers must be a positive integer")
    if (
        isinstance(chunk_size, bool)
        or not isinstance(chunk_size, int)
        or chunk_size < 1
    ):
        raise ValueError("chunk_size must be a positive integer")
    tree = cKDTree(positions)
    radii = np.empty(len(positions), dtype=float)
    for start in range(0, len(positions), chunk_size):
        end = min(start + chunk_size, len(positions))
        distances = tree.query(positions[start:end], k=neighbours + 1, workers=workers)[
            0
        ]
        radii[start:end] = distances[:, -1]
    return radii


def _projected_cubic_spline_table(size=4097, quadrature_order=64):
    """Tabulate the LOS integral of a spherical M4 kernel.

    The 3D kernel uses support radius H and the Gadget/Arepo convention
    ``f(u)=1-6u^2+6u^3`` below u=1/2 and ``2(1-u)^3`` up to u=1. Its
    dimensional normalization cancels when each discrete particle stencil
    is flux-normalized, but the line-of-sight integral sets the radial shape.
    """
    projected_radius = np.linspace(0.0, 1.0, size)
    nodes, weights = np.polynomial.legendre.leggauss(quadrature_order)
    height = np.sqrt(np.maximum(0.0, 1.0 - projected_radius**2))
    radius_3d = np.sqrt(
        projected_radius[:, None] ** 2 + (height[:, None] * nodes[None, :]) ** 2
    )
    shape = np.where(
        radius_3d < 0.5,
        1.0 - 6.0 * radius_3d**2 + 6.0 * radius_3d**3,
        2.0 * (1.0 - radius_3d) ** 3,
    )
    return projected_radius, height * np.sum(weights[None, :] * shape, axis=1)


_PROJECTED_RADIUS, _PROJECTED_KERNEL = _projected_cubic_spline_table()


@njit(cache=True)
def _projected_kernel_value(q, table):
    """Linearly interpolate the dimensionless projected-kernel table."""
    if q >= 1.0:
        return 0.0
    coordinate = q * (len(table) - 1)
    lower = int(coordinate)
    fraction = coordinate - lower
    return table[lower] * (1.0 - fraction) + table[lower + 1] * fraction


@njit(cache=True)
def _accumulate(positions, values, effective, ny, nx, pixel, kernel_table):
    """Normalize each kernel, then deposit only its image-intersecting stencil."""
    image = np.zeros((ny, nx), dtype=np.float64)
    outside = 0.0
    clipped = 0.0
    contributing = 0
    for i in range(len(values)):
        y, x = positions[i, 0], positions[i, 1]
        value = values[i]
        if value == 0:
            continue
        row = (ny - 1) / 2 - y / pixel
        col = x / pixel + (nx - 1) / 2
        extent = effective[i] / pixel
        r0, r1 = math.ceil(row - extent), math.floor(row + extent)
        c0, c1 = math.ceil(col - extent), math.floor(col + extent)
        if r1 < 0 or r0 >= ny or c1 < 0 or c0 >= nx:
            outside += value
            continue
        # The line-of-sight integral of the unnormalised 3D M4 kernel has
        # plane integral pi/8. For a wide, well-sampled kernel, its full
        # discrete stencil sum is (pi/8)*extent**2 to ~1e-7 relative accuracy.
        # This avoids scanning millions of off-image fine pixels for the few
        # broad stellar kernels; only the image-intersecting stencil is drawn.
        if extent >= 32.0:
            total = (math.pi / 8.0) * extent * extent
        else:
            total = 0.0
            for r in range(r0, r1 + 1):
                for c in range(c0, c1 + 1):
                    q = math.hypot(r - row, c - col) / extent
                    total += _projected_kernel_value(q, kernel_table)
        if total <= 0:
            raise RuntimeError("SPH kernel has no sampled pixel centres")
        delivered = 0.0
        for r in range(max(r0, 0), min(r1, ny - 1) + 1):
            for c in range(max(c0, 0), min(c1, nx - 1) + 1):
                q = math.hypot(r - row, c - col) / extent
                contribution = value * _projected_kernel_value(q, kernel_table) / total
                image[r, c] += contribution
                delivered += contribution
        clipped += value - delivered
        contributing += 1
    return image, outside, clipped, contributing


def sph_project(
    positions_yx_arcsec, flux, support_radii_arcsec, grid: ImageGrid, *, subpixels=8
):
    """Deposit particle flux with a projected spherical cubic-spline kernel.

    The 3D neighbour distance is the spherical kernel support radius H. Its
    line-of-sight integral is tabulated by Gaussian quadrature. Projected
    values are evaluated on an oversampled grid, normalized over the complete,
    unbounded local stencil, and block-integrated onto the requested pixels.
    A minimal numerical floor below one fine pixel only guarantees that an
    unresolved kernel is sampled; no gravitational-softening floor is imposed.

    Particles outside the image may still contribute when their kernels
    overlap it. Flux missed by the finite field is recorded, never silently
    renormalized into the image. The result is counts/s/pixel, not surface
    brightness per square arcsecond and not a PSF-convolved observation.
    """
    positions = finite_array(positions_yx_arcsec, "positions_yx_arcsec")
    values = finite_array(flux, "flux")
    radii = finite_array(support_radii_arcsec, "support_radii_arcsec")
    if positions.ndim != 2 or positions.shape[1] != 2:
        raise ValueError("positions_yx_arcsec must have shape (N, 2)")
    if values.shape != (len(positions),) or radii.shape != values.shape:
        raise ValueError("flux and support_radii_arcsec must have shape (N,)")
    if np.any(values < 0) or np.any(radii < 0):
        raise ValueError("flux and support radii must be non-negative")
    if isinstance(subpixels, bool) or not isinstance(subpixels, int) or subpixels < 1:
        raise ValueError("subpixels must be a positive integer")
    fine_pixel = grid.pixel_scale_arcsec / subpixels
    sampling_floor = np.sqrt(0.5) * fine_pixel * (1.0 + 1e-12)
    effective = np.maximum(radii, sampling_floor)
    if len(effective) and np.any(effective == 0):
        raise ValueError("effective SPH support radii must be positive")
    ny, nx = grid.shape
    fine_image, outside, clipped, contributing = _accumulate(
        positions,
        values,
        effective,
        ny * subpixels,
        nx * subpixels,
        fine_pixel,
        _PROJECTED_KERNEL,
    )
    image = fine_image.reshape(ny, subpixels, nx, subpixels).sum(axis=(1, 3))
    accounting = {
        "input_flux": float(values.sum()),
        "outside_fov_flux": outside,
        "edge_kernel_flux": clipped,
        "deposited_flux": float(image.sum()),
        "input_particles": len(values),
        "contributing_particles": contributing,
        "renderer": "numba" if _NUMBA else "python",
        "neighbour_support_radius_arcsec_min": float(radii.min())
        if len(radii)
        else None,
        "neighbour_support_radius_arcsec_median": float(np.median(radii))
        if len(radii)
        else None,
        "neighbour_support_radius_arcsec_max": float(radii.max())
        if len(radii)
        else None,
        "kernel": "line-of-sight projected 3D cubic spline (M4)",
        "subpixels": subpixels,
        "sampling_floor_arcsec": sampling_floor,
        "sampling_floor_particles": int(np.count_nonzero(radii < sampling_floor)),
    }
    return image, accounting
