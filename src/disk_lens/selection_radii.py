"""Production SPH particle-mass selection radii, independent of mock fitting."""

from datetime import datetime, timezone
from importlib.metadata import version
import hashlib
import json
from pathlib import Path
import time

import h5py
import numpy as np
from astropy.cosmology import Planck18

from .geometry import periodic_offsets, principal_axes, stellar_centre
from .grid import ImageGrid
from .photometry import critical_surface_density_msun_kpc2
from .mass_lensing import critical_radius, point_mass_pixels, projected_mass

ZL, ZS = 0.15, 0.65


def _method_hash():
    digest = hashlib.sha256()
    for name in [
        "selection_radii.py",
        "mass_lensing.py",
        "sph.py",
        "geometry.py",
        "photometry.py",
    ]:
        digest.update((Path(__file__).parent / name).read_bytes())
    return digest.hexdigest()


def _load_components(tng_root, snapshot, sid):
    import illustris_python as il

    header = il.groupcat.loadHeader(tng_root, snapshot)
    h, a, box = (
        float(header["HubbleParam"]),
        float(header["Time"]),
        float(header["BoxSize"]),
    )
    catalog = il.groupcat.loadSingle(tng_root, snapshot, subhaloID=sid)
    origin = catalog["SubhaloPos"]
    rhalf = float(catalog["SubhaloHalfmassRadType"][4]) * a / h

    def offsets(raw):
        return periodic_offsets(raw, origin, box_size=box) * a / h

    star = il.snapshot.loadSubhalo(
        tng_root,
        snapshot,
        sid,
        4,
        ["Coordinates", "Masses", "StellarHsml", "GFM_StellarFormationTime"],
    )
    star_xyz = offsets(star["Coordinates"])
    star_m = star["Masses"].astype(float) * 1e10 / h
    real = star["GFM_StellarFormationTime"] > 0
    centre, iterations = stellar_centre(star_xyz[real], star_m[real])
    rotation, shape = principal_axes(star_xyz[real] - centre, star_m[real], 5 * rhalf)
    components = {
        "stars_and_winds": (star_xyz, star_m, star["StellarHsml"].astype(float) * a / h)
    }
    del star

    dm = il.snapshot.loadSubhalo(
        tng_root, snapshot, sid, 1, ["Coordinates", "SubfindHsml"]
    )
    with h5py.File(il.snapshot.snapPath(tng_root, snapshot), "r") as handle:
        dm_mass = float(handle["Header"].attrs["MassTable"][1]) * 1e10 / h
    components["dark_matter"] = (
        offsets(dm["Coordinates"]),
        np.full(dm["count"], dm_mass),
        dm["SubfindHsml"].astype(float) * a / h,
    )
    del dm

    gas = il.snapshot.loadSubhalo(
        tng_root, snapshot, sid, 0, ["Coordinates", "Masses", "Density"]
    )
    if gas["count"]:
        volume = gas["Masses"].astype(float) / gas["Density"].astype(float)
        support = 2 * np.cbrt(3 * volume / (4 * np.pi)) * a / h
        components["gas"] = (
            offsets(gas["Coordinates"]),
            gas["Masses"].astype(float) * 1e10 / h,
            support,
        )
    del gas
    bh = il.snapshot.loadSubhalo(
        tng_root, snapshot, sid, 5, ["Coordinates", "Masses", "BH_Mass"]
    )
    bh_component = (
        (offsets(bh["Coordinates"]), bh["Masses"].astype(float) * 1e10 / h)
        if bh["count"]
        else (np.empty((0, 3)), np.empty(0))
    )
    bh_physical_mass = float(np.sum(bh.get("BH_Mass", []))) * 1e10 / h
    # Each stored support is defined in the full snapshot, before a FOV cut.
    # Keep all bound particles along the LOS, as in the historical selection.
    rotated = {}
    for name, (xyz, masses, support) in components.items():
        if not np.all(np.isfinite(support)) or np.any(support <= 0):
            raise ValueError(f"Invalid stored/cell-derived support in {name}")
        rotated[name] = ((xyz - centre) @ rotation.T, masses, support)
    bh_component = ((bh_component[0] - centre) @ rotation.T, bh_component[1])
    return (
        rotated,
        bh_component,
        {
            "halfmass_radius_kpc": rhalf,
            "scale_factor": a,
            "hubble_h": h,
            "dm_particle_mass_msun": dm_mass,
            "centre_offset_from_subhalo_kpc": centre.tolist(),
            "centre_iterations": iterations,
            "shape": shape,
            "stellar_particles": int(real.sum()),
            "wind_particles": int((~real).sum()),
            "bh_dynamical_mass_msun": float(bh_component[1].sum()),
            "bh_physical_mass_msun": bh_physical_mass,
        },
    )


def _render_views(
    sid, components, bh, meta, pixels, subpixels, fov_multiplier, supports_scale
):
    begin = time.monotonic()
    kpc_arcsec = float(Planck18.kpc_proper_per_arcmin(ZL).value) / 60
    sigma_crit = critical_surface_density_msun_kpc2(
        lens_redshift=ZL, source_redshift=ZS, cosmology=Planck18
    )
    result = {"maps": {}}
    for pixel in pixels:
        n = max(
            32,
            int(
                np.ceil(
                    fov_multiplier * meta["halfmass_radius_kpc"] / kpc_arcsec / pixel
                )
            ),
        )
        n += n % 2
        grid = ImageGrid((n, n), pixel)
        for view, axes in [("xy", [1, 0]), ("xz", [2, 0]), ("yz", [2, 1])]:
            mass_map = np.zeros(grid.shape)
            accounting = {}
            for name, (xyz, masses, support) in components.items():
                component, info = projected_mass(
                    xyz[:, axes] / kpc_arcsec,
                    masses,
                    support * supports_scale / kpc_arcsec,
                    grid,
                    subpixels=subpixels,
                )
                mass_map += component
                accounting[name] = info
            kappa = mass_map / ((pixel * kpc_arcsec) ** 2 * sigma_crit)
            radius = critical_radius(kappa, grid)
            # Preserve the historical three-component scope for selection.
            # Show the compact BH dynamical mass as a distinct physical correction.
            bh_map = point_mass_pixels(bh[0][:, axes] / kpc_arcsec, bh[1], grid)
            with_bh = critical_radius(
                kappa + bh_map / ((pixel * kpc_arcsec) ** 2 * sigma_crit), grid
            )
            result["maps"][f"{pixel:g}_{view}"] = {
                **radius,
                "pixel_scale_arcsec": pixel,
                "grid_pixels": n,
                "full_width_kpc": n * pixel * kpc_arcsec,
                "mass_accounting": accounting,
                "with_bh_radius_arcsec": with_bh["einstein_radius_arcsec"],
                "bh_deposited_mass_msun": float(bh_map.sum()),
            }
            print(
                f"{sid} pixel={pixel:g} view={view} R={radius['einstein_radius_arcsec']:.6f} elapsed={time.monotonic() - begin:.1f}s",
                flush=True,
            )
    return result["maps"]


def calculate_selection_radius(
    tng_root,
    *,
    snapshot=84,
    subhalo_id,
    output_dir,
    pixel_scales=(0.05, 0.025),
    subpixels=2,
    fov_multiplier=10.0,
    support_scale=1.0,
):
    """Calculate one particle-mass selection radius in three principal views.

    All parameter and output checks precede particle loading. An existing
    result is never overwritten. This calculation does not change sample
    membership, inherited mock radii or mock/fitting mask inputs.
    """
    sid = subhalo_id
    pixels = tuple(pixel_scales)
    supports_scale = support_scale
    output = Path(output_dir)
    if isinstance(sid, bool) or not isinstance(sid, int) or sid < 0:
        raise ValueError("subhalo_id must be a non-negative integer")
    if isinstance(snapshot, bool) or not isinstance(snapshot, int) or snapshot < 0:
        raise ValueError("snapshot must be a non-negative integer")
    if not pixels or any(not np.isfinite(p) or p <= 0 for p in pixels):
        raise ValueError("pixel_scales must be positive and finite")
    if isinstance(subpixels, bool) or not isinstance(subpixels, int) or subpixels < 1:
        raise ValueError("subpixels must be a positive integer")
    if (
        not np.isfinite(fov_multiplier)
        or fov_multiplier <= 0
        or not np.isfinite(support_scale)
        or support_scale <= 0
    ):
        raise ValueError("FOV multiple and support scale must be positive and finite")
    if (output / f"{sid}.json").exists():
        raise FileExistsError(f"selection result already exists: {output}/{sid}.json")
    begin = time.monotonic()
    components, bh, meta = _load_components(tng_root, snapshot, sid)
    result = {
        "id": sid,
        "method_hash": _method_hash(),
        "software_versions": {
            name: version(name)
            for name in ["numpy", "scipy", "astropy", "numba", "scikit-image", "h5py"]
        },
        "snapshot": snapshot,
        "tng_root": str(Path(tng_root).resolve()),
        "metadata": meta,
        "maps": {},
        "subpixels": subpixels,
        "fov_full_width_halfmass_multiple": fov_multiplier,
        "support_scale": supports_scale,
        "lens_redshift": ZL,
        "source_redshift": ZS,
        "scope": "bound subhalo stars+winds, DM and gas; BH correction reported separately",
        "radius_definition": "largest closed lambda_t=0 curve sqrt(A/pi)",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }
    result["maps"] = _render_views(
        sid, components, bh, meta, pixels, subpixels, fov_multiplier, supports_scale
    )
    fine, coarse = min(pixels), max(pixels)
    result["fine_max_radius_arcsec"] = max(
        result["maps"][f"{fine:g}_{v}"]["einstein_radius_arcsec"]
        for v in ["xy", "xz", "yz"]
    )
    result["coarse_max_radius_arcsec"] = max(
        result["maps"][f"{coarse:g}_{v}"]["einstein_radius_arcsec"]
        for v in ["xy", "xz", "yz"]
    )
    result["passes_einstein_strict_0p5"] = result["fine_max_radius_arcsec"] > 0.5
    # The screening statistic is the maximum of the three physical views.
    # Preserve the stricter all-view flag, including weak/subpixel curves.
    result["selection_max_resolution_delta_arcsec"] = abs(
        result["fine_max_radius_arcsec"] - result["coarse_max_radius_arcsec"]
    )
    result["max_resolution_delta_arcsec"] = max(
        abs(
            result["maps"][f"{fine:g}_{v}"]["einstein_radius_arcsec"]
            - result["maps"][f"{coarse:g}_{v}"]["einstein_radius_arcsec"]
        )
        for v in ["xy", "xz", "yz"]
    )
    # Classify only when resolution error is small and cannot move it across .5.
    tolerance = max(0.005, 0.01 * result["fine_max_radius_arcsec"])
    result["resolution_converged"] = result["max_resolution_delta_arcsec"] <= tolerance
    result["selection_resolution_converged"] = (
        result["selection_max_resolution_delta_arcsec"] <= tolerance
    )
    result["threshold_stable"] = (
        abs(result["fine_max_radius_arcsec"] - 0.5)
        > result["max_resolution_delta_arcsec"] + 0.0025
    )
    if len(set(pixels)) < 2:
        result["max_resolution_delta_arcsec"] = None
        result["resolution_converged"] = False
        result["selection_max_resolution_delta_arcsec"] = None
        result["selection_resolution_converged"] = False
        result["threshold_stable"] = False
        result["resolution_note"] = "Single-grid calculation; convergence not checked"
    field_maps = _render_views(
        sid,
        components,
        bh,
        meta,
        (coarse,),
        subpixels,
        2 * fov_multiplier,
        supports_scale,
    )
    field_delta = max(
        abs(
            field_maps[f"{coarse:g}_{view}"]["einstein_radius_arcsec"]
            - result["maps"][f"{coarse:g}_{view}"]["einstein_radius_arcsec"]
        )
        for view in ["xy", "xz", "yz"]
    )
    result["field_converged"] = field_delta <= tolerance
    result["field_check"] = {
        "full_width_halfmass_multiple": 2 * fov_multiplier,
        "pixel_scale_arcsec": coarse,
        "max_radius_delta_arcsec": field_delta,
        "radii_arcsec": {
            view: field_maps[f"{coarse:g}_{view}"]["einstein_radius_arcsec"]
            for view in ["xy", "xz", "yz"]
        },
    }
    if result["max_resolution_delta_arcsec"] is not None:
        result["numerical_error_arcsec"] = (
            result["max_resolution_delta_arcsec"] + field_delta
        )
        result["threshold_stable"] = (
            result["selection_resolution_converged"]
            and result["field_converged"]
            and abs(result["fine_max_radius_arcsec"] - 0.5)
            > result["numerical_error_arcsec"] + 0.1 * fine
        )
    else:
        result["numerical_error_arcsec"] = None
    result["elapsed_seconds"] = time.monotonic() - begin
    output.mkdir(parents=True, exist_ok=True)
    path = output / f"{sid}.json"
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)
    return result
