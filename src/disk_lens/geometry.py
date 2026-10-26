"""Particle centring and mass-weighted reduced-inertia principal axes."""

import numpy as np

from ._validation import finite_array, scalar


def _particles(positions, masses):
    x = finite_array(positions, "positions")
    m = finite_array(masses, "masses")
    if x.ndim != 2 or x.shape[1] != 3 or m.shape != (len(x),) or len(x) < 4:
        raise ValueError("need at least four positions (N,3) and masses (N,)")
    if np.any(m <= 0):
        raise ValueError("particle masses must be positive")
    return x, m


def periodic_offsets(coordinates, centre, *, box_size: float):
    """Minimum-image offsets; inputs all use the same original length units."""
    x = finite_array(coordinates, "coordinates")
    c = finite_array(centre, "centre")
    size = scalar(box_size, "box_size", positive=True)
    if x.ndim != 2 or x.shape[1] != 3 or c.shape != (3,):
        raise ValueError("coordinates must be (N,3), centre (3,)")
    return (x - c + size / 2) % size - size / 2


def stellar_centre(positions, masses, fraction=0.2, tolerance=0.01, max_iterations=100):
    """Literal centre iteration of the supplied image_light.ipynb."""
    centre = np.average(positions, axis=0, weights=masses)
    for iteration in range(1, max_iterations + 1):
        radius = np.linalg.norm(positions - centre, axis=1)
        selected = radius < np.percentile(radius, 100 * fraction)
        updated = np.average(positions[selected], axis=0, weights=masses[selected])
        if np.linalg.norm(updated - centre) <= tolerance:
            return updated, iteration
        centre = updated
    return centre, max_iterations


def principal_axes(positions, masses, major_radius, tolerance=0.01, max_iterations=100):
    """Literal axis iteration of the supplied image_light.ipynb."""
    s, q = 1.0, 1.0
    transform = np.eye(3)
    order = [2, 1, 0]
    for iteration in range(1, max_iterations + 1):
        y = positions @ transform.T
        radius = np.sqrt(
            y[:, order[2]] ** 2
            + y[:, order[1]] ** 2 / q**2
            + y[:, order[0]] ** 2 / s**2
        )
        selected = (radius < major_radius) & (radius > 0)
        points, weights, r = y[selected], masses[selected], radius[selected]
        tensor = (points.T * (weights / r**2)) @ points
        eigenvalues, eigenvectors = np.linalg.eig(tensor)
        order = np.argsort(eigenvalues)
        new_q = float(np.sqrt(eigenvalues[order[1]] / eigenvalues[order[2]]))
        new_s = float(np.sqrt(eigenvalues[order[0]] / eigenvalues[order[2]]))
        transform = np.linalg.inv(eigenvectors) @ transform
        if max(abs(q - new_q), abs(s - new_s)) <= tolerance:
            q, s = new_q, new_s
            break
        q, s = new_q, new_s
    transform = transform[order[::-1], :]
    if np.dot(transform[2], np.cross(transform[0], transform[1])) < 0:
        transform[1] *= -1
    return transform, {"b_over_a": q, "c_over_a": s, "iterations": iteration}
