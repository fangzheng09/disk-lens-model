"""Portable stellar-particle input and an optional IllustrisTNG reader adapter."""

from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np

from .geometry import _particles, periodic_offsets
from .units import comoving_ckpc_h_to_physical_kpc, snapshot_mass_to_msun


@dataclass
class StellarParticles:
    """Physical offsets, current stellar masses, and intrinsic SDSS i magnitudes.

    Wind particles must already be removed. Coordinates are not yet aligned.
    The metadata records how this input was prepared.
    """

    positions_kpc: np.ndarray
    masses_msun: np.ndarray
    absolute_magnitudes_i: np.ndarray
    halfmass_radius_kpc: float
    metadata: dict
    smoothing_lengths_kpc: np.ndarray | None = None

    def __post_init__(self):
        from ._validation import finite_array, scalar

        self.positions_kpc, self.masses_msun = _particles(
            self.positions_kpc, self.masses_msun
        )
        self.absolute_magnitudes_i = finite_array(
            self.absolute_magnitudes_i, "magnitudes"
        )
        if self.absolute_magnitudes_i.shape != self.masses_msun.shape:
            raise ValueError("magnitudes must have one value per particle")
        if self.smoothing_lengths_kpc is not None:
            self.smoothing_lengths_kpc = finite_array(
                self.smoothing_lengths_kpc, "smoothing_lengths_kpc"
            )
            if self.smoothing_lengths_kpc.shape != self.masses_msun.shape:
                raise ValueError("smoothing lengths must have one value per particle")
            if np.any(self.smoothing_lengths_kpc <= 0):
                raise ValueError("smoothing lengths must be positive")
        self.halfmass_radius_kpc = scalar(
            self.halfmass_radius_kpc, "halfmass_radius_kpc", positive=True
        )

    def save(self, filename):
        # Exclusive creation prevents accidental replacement of a particle input.
        arrays = {
            "positions_kpc": self.positions_kpc,
            "masses_msun": self.masses_msun,
            "absolute_magnitudes_i": self.absolute_magnitudes_i,
            "halfmass_radius_kpc": self.halfmass_radius_kpc,
            "metadata_json": json.dumps(self.metadata, allow_nan=False),
        }
        if self.smoothing_lengths_kpc is not None:
            arrays["smoothing_lengths_kpc"] = self.smoothing_lengths_kpc
        with Path(filename).open("xb") as handle:
            np.savez_compressed(handle, **arrays)

    @classmethod
    def load(cls, filename):
        with np.load(filename, allow_pickle=False) as data:
            smoothing = (
                data["smoothing_lengths_kpc"]
                if "smoothing_lengths_kpc" in data
                else None
            )
            return cls(
                data["positions_kpc"],
                data["masses_msun"],
                data["absolute_magnitudes_i"],
                float(data["halfmass_radius_kpc"]),
                json.loads(str(data["metadata_json"])),
                smoothing,
            )


def load_tng_stars(base_path, *, snapshot: int, subhalo_id: int) -> StellarParticles:
    """Read one subhalo through the installed illustris_python package.

    The reader's path conventions must match the data layout. This adapter
    does not alter sys.path, monkeypatch the reader or copy upstream code.
    """
    try:
        import illustris_python as il
    except ImportError as exc:
        raise ImportError(
            "TNG extraction requires an installed illustris_python reader; see docs/INSTALLATION.md"
        ) from exc
    if snapshot < 0 or subhalo_id < 0:
        raise ValueError("snapshot and subhalo_id must be non-negative")
    base_path = str(Path(base_path).resolve())
    header = il.groupcat.loadHeader(base_path, snapshot)
    subhalo = il.groupcat.loadSingle(base_path, snapshot, subhaloID=subhalo_id)
    data = il.snapshot.loadSubhalo(
        base_path,
        snapshot,
        subhalo_id,
        partType=4,
        fields=[
            "Coordinates",
            "Masses",
            "GFM_StellarPhotometrics",
            "GFM_StellarFormationTime",
            "StellarHsml",
        ],
    )
    if data.get("count", 0) < 4:
        raise ValueError("subhalo has insufficient PartType4 particles")
    stars = data["GFM_StellarFormationTime"] > 0
    h = float(header["HubbleParam"])
    a = float(header.get("Time", 1 / (1 + header["Redshift"])))
    offsets = periodic_offsets(
        data["Coordinates"][stars],
        subhalo["SubhaloPos"],
        box_size=float(header["BoxSize"]),
    )
    return StellarParticles(
        comoving_ckpc_h_to_physical_kpc(offsets, scale_factor=a, hubble_h=h),
        snapshot_mass_to_msun(data["Masses"][stars], hubble_h=h),
        data["GFM_StellarPhotometrics"][stars, 6],
        float(
            comoving_ckpc_h_to_physical_kpc(
                subhalo["SubhaloHalfmassRadType"][4], scale_factor=a, hubble_h=h
            )
        ),
        {
            "simulation": "TNG50-1",
            "snapshot": snapshot,
            "subhalo_id": subhalo_id,
            "snapshot_redshift": float(header["Redshift"]),
            "scale_factor": a,
            "hubble_h": h,
            "wind_filter": "GFM_StellarFormationTime > 0",
            "excluded_wind_particles": int((~stars).sum()),
            "origin": "SubhaloPos; periodic minimum-image offsets",
            "subhalo_position_ckpc_h": np.asarray(subhalo["SubhaloPos"]).tolist(),
            "photometry": "GFM_StellarPhotometrics[:,6], intrinsic SDSS i AB",
            "smoothing": "PartType4/StellarHsml; physical kpc; 32+/-1 same-type neighbours",
        },
        comoving_ckpc_h_to_physical_kpc(
            data["StellarHsml"][stars], scale_factor=a, hubble_h=h
        ),
    )
