"""Convert raw TNG snapshot values to physical units.

Inputs are numerical values in the TNG code units named in each function.
No centring, periodic wrapping, particle selection or bulk subtraction is done.
In particular, SubhaloVel is already in peculiar km/s and must not be passed
to the particle-velocity conversion.

Reference: https://www.tng-project.org/data/docs/specifications/
"""

import numpy as np

from ._validation import finite_array, scalar


def comoving_ckpc_h_to_physical_kpc(
    values, *, scale_factor: float, hubble_h: float
) -> np.ndarray:
    """Convert ckpc/h coordinates or lengths to physical kpc using a/h."""
    a = scalar(scale_factor, "scale_factor", positive=True)
    h = scalar(hubble_h, "hubble_h", positive=True)
    return finite_array(values, "values") * (a / h)


def snapshot_mass_to_msun(values, *, hubble_h: float) -> np.ndarray:
    """Convert masses from 1e10 Msun/h to Msun; scale factor is not involved."""
    h = scalar(hubble_h, "hubble_h", positive=True)
    masses = finite_array(values, "values")
    if np.any(masses < 0):
        raise ValueError("masses must be non-negative")
    return masses * (1.0e10 / h)


def snapshot_velocity_to_peculiar_kms(values, *, scale_factor: float) -> np.ndarray:
    """Convert snapshot particle Velocities to peculiar km/s using sqrt(a)."""
    a = scalar(scale_factor, "scale_factor", positive=True)
    return finite_array(values, "values") * np.sqrt(a)
