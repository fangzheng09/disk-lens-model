"""Angular image grids and FITS coordinate headers."""

from dataclasses import dataclass

import numpy as np
from astropy.io import fits

from ._validation import scalar


@dataclass(frozen=True)
class ImageGrid:
    """Pixel centres: x increases with column, y decreases with row.

    The origin is the geometrical centre, including for even-sized images.
    Positions and deflections use (y, x), in arcseconds. Plot with origin='upper'.
    """

    shape: tuple[int, int]
    pixel_scale_arcsec: float

    def __post_init__(self):
        if len(self.shape) != 2 or any(
            isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 2
            for n in self.shape
        ):
            raise ValueError("shape must contain two integers >= 2")
        object.__setattr__(self, "shape", tuple(int(n) for n in self.shape))
        object.__setattr__(
            self,
            "pixel_scale_arcsec",
            scalar(self.pixel_scale_arcsec, "pixel_scale_arcsec", positive=True),
        )

    def coordinates(self):
        ny, nx = self.shape
        y = ((ny - 1) / 2 - np.arange(ny)) * self.pixel_scale_arcsec
        x = (np.arange(nx) - (nx - 1) / 2) * self.pixel_scale_arcsec
        return np.stack(np.meshgrid(y, x, indexing="ij"), axis=-1)

    def header(self, *, bunit: str) -> fits.Header:
        ny, nx = self.shape
        return fits.Header(
            {
                "CTYPE1": "LINEAR",
                "CTYPE2": "LINEAR",
                "CUNIT1": "arcsec",
                "CUNIT2": "arcsec",
                "CRPIX1": (nx + 1) / 2,
                "CRPIX2": (ny + 1) / 2,
                "CRVAL1": 0.0,
                "CRVAL2": 0.0,
                "CDELT1": self.pixel_scale_arcsec,
                "CDELT2": -self.pixel_scale_arcsec,
                "PIXSCALE": self.pixel_scale_arcsec,
                "BUNIT": bunit,
                "YXORDER": "row-down; y-up",
            }
        )
