"""Per-system inputs kept separate from sample-set membership."""

from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class SystemInputs:
    subhalo_id: int
    grid_shape: tuple[int, int]
    pixel_scale_arcsec: float
    mask_radius_arcsec: float
    source_centre_yx_arcsec: tuple[float, float] | None
    positions_yx_arcsec: tuple[tuple[float, float], ...]
    stellar_mass_normalization_legacy: float
    halo_mass_at_200_msun_legacy: float
    mock_einstein_radius_arcsec_legacy: float | None
    status: str
    stellar_mass_normalization: float | None = None
    central_deflection_radius_arcsec: float | None = None

    @property
    def adopted_stellar_mass_normalization(self):
        """Prefer the recalculated coefficient, with historical-file fallback."""
        return (
            self.stellar_mass_normalization_legacy
            if self.stellar_mass_normalization is None
            else self.stellar_mass_normalization
        )


def load_system_inputs(path, subhalo_id):
    """Load one provisional or validated system record from schema version 1."""
    raw = json.loads(Path(path).read_text())
    if raw.get("schema_version") != 1 or not isinstance(raw.get("systems"), dict):
        raise ValueError(
            "system catalogue must use schema_version 1 and contain systems"
        )
    key = str(int(subhalo_id))
    if key not in raw["systems"]:
        raise KeyError(f"subhalo {subhalo_id} is absent from the system catalogue")
    item = raw["systems"][key]
    required = {
        "grid_shape",
        "pixel_scale_arcsec",
        "mask_radius_arcsec",
        "source_centre_yx_arcsec",
        "positions_yx_arcsec",
        "stellar_mass_normalization_legacy",
        "halo_mass_at_200_msun_legacy",
        "mock_einstein_radius_arcsec_legacy",
        "status",
    }
    missing = required - item.keys()
    if missing:
        raise ValueError(f"system {subhalo_id} is missing {sorted(missing)}")
    shape = tuple(item["grid_shape"])
    positions = np.asarray(item["positions_yx_arcsec"], dtype=float)
    centre = item["source_centre_yx_arcsec"]
    if (
        len(shape) != 2
        or any(not isinstance(n, int) or n < 2 for n in shape)
        or positions.ndim != 2
        or positions.shape[1:] != (2,)
        or not np.all(np.isfinite(positions))
    ):
        raise ValueError(f"system {subhalo_id} has invalid grid or positions")
    if centre is not None and (len(centre) != 2 or not np.all(np.isfinite(centre))):
        raise ValueError(f"system {subhalo_id} has an invalid source centre")
    for name in [
        "pixel_scale_arcsec",
        "mask_radius_arcsec",
        "stellar_mass_normalization_legacy",
        "halo_mass_at_200_msun_legacy",
    ]:
        if not np.isfinite(item[name]) or item[name] <= 0:
            raise ValueError(f"system {subhalo_id}.{name} must be positive and finite")
    normalization = item.get("stellar_mass_normalization")
    if normalization is not None and (
        isinstance(normalization, bool)
        or not isinstance(normalization, (int, float))
        or not np.isfinite(normalization)
        or normalization <= 0
    ):
        raise ValueError(
            f"system {subhalo_id}.stellar_mass_normalization must be positive and finite"
        )
    einstein = item["mock_einstein_radius_arcsec_legacy"]
    central_radius = item.get("central_deflection_radius_arcsec")
    if central_radius is not None and (
        isinstance(central_radius, bool)
        or not isinstance(central_radius, (int, float))
        or not np.isfinite(central_radius)
        or central_radius < 0
    ):
        raise ValueError(
            f"system {subhalo_id}.central_deflection_radius_arcsec must be finite and non-negative"
        )
    if einstein is not None and (not np.isfinite(einstein) or einstein <= 0):
        raise ValueError(f"system {subhalo_id} has an invalid Einstein radius")
    return SystemInputs(
        subhalo_id=int(key),
        grid_shape=(int(shape[0]), int(shape[1])),
        pixel_scale_arcsec=float(item["pixel_scale_arcsec"]),
        mask_radius_arcsec=float(item["mask_radius_arcsec"]),
        source_centre_yx_arcsec=None if centre is None else tuple(map(float, centre)),
        positions_yx_arcsec=tuple(tuple(map(float, value)) for value in positions),
        stellar_mass_normalization_legacy=float(
            item["stellar_mass_normalization_legacy"]
        ),
        halo_mass_at_200_msun_legacy=float(item["halo_mass_at_200_msun_legacy"]),
        mock_einstein_radius_arcsec_legacy=None
        if einstein is None
        else float(einstein),
        status=str(item["status"]),
        stellar_mass_normalization=None
        if normalization is None
        else float(normalization),
        central_deflection_radius_arcsec=None
        if central_radius is None
        else float(central_radius),
    )
