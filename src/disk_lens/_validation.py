"""Small numerical input checks shared by the physical conversion functions."""

import numpy as np


def finite_array(value, name: str) -> np.ndarray:
    """Convert real, finite numeric input to a floating array (including scalars)."""
    if not np.isrealobj(value):
        raise ValueError(f"{name} must be real")
    array = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite")
    return array


def scalar(value, name: str, *, positive: bool = False) -> float:
    """Require a finite scalar, optionally strictly positive."""
    array = finite_array(value, name)
    if array.ndim != 0:
        raise ValueError(f"{name} must be a scalar")
    result = float(array)
    if positive and result <= 0:
        raise ValueError(f"{name} must be positive")
    return result
