"""Bright-arc candidates and coordinate checks for generated mock images."""

import csv
import json
from pathlib import Path

import numpy as np
from importlib.metadata import version
from scipy.ndimage import map_coordinates, maximum_filter

from ._io import checksum, validate_bundle, write_json
from ._validation import finite_array, scalar
from .grid import ImageGrid
from ._paths import resource_path
from .systems import load_system_inputs


def _sample(image, position, pixel):
    row = (image.shape[0] - 1) / 2 - position[0] / pixel
    col = (image.shape[1] - 1) / 2 + position[1] / pixel
    inside = 0 <= row <= image.shape[0] - 1 and 0 <= col <= image.shape[1] - 1
    value = float(map_coordinates(image, [[row], [col]], order=1, mode="constant")[0])
    return value, bool(inside)


def _peak_offset(image, row, col):
    """Estimate a local maximum with a negative-definite quadratic patch."""
    patch = image[row - 1 : row + 2, col - 1 : col + 2]
    gradient = np.array(
        [(patch[2, 1] - patch[0, 1]) / 2, (patch[1, 2] - patch[1, 0]) / 2]
    )
    mixed = (patch[2, 2] - patch[2, 0] - patch[0, 2] + patch[0, 0]) / 4
    hessian = np.array(
        [
            [patch[2, 1] - 2 * patch[1, 1] + patch[0, 1], mixed],
            [mixed, patch[1, 2] - 2 * patch[1, 1] + patch[1, 0]],
        ]
    )
    if np.max(np.linalg.eigvalsh(hessian)) >= 0 or np.linalg.cond(hessian) > 1e6:
        return np.zeros(2), False
    offset = -np.linalg.solve(hessian, gradient)
    if not np.all(np.isfinite(offset)) or np.max(np.abs(offset)) > 0.75:
        return np.zeros(2), False
    return offset, True


def find_arc_peaks(
    image,
    *,
    pixel_scale_arcsec,
    min_peak_fraction=0.3,
    min_separation_arcsec=0.2,
    max_candidates=12,
    min_radius_arcsec=0.0,
):
    """Find separated brightness maxima; do not infer distinct lensed images."""
    values = finite_array(image, "source image")
    if values.ndim != 2 or min(values.shape) < 3 or values.max() <= 0:
        raise ValueError("source image must be 2D, at least 3x3, with positive flux")
    pixel = scalar(pixel_scale_arcsec, "pixel_scale_arcsec", positive=True)
    fraction = scalar(min_peak_fraction, "min_peak_fraction", positive=True)
    separation = scalar(min_separation_arcsec, "min_separation_arcsec", positive=True)
    radius = scalar(min_radius_arcsec, "min_radius_arcsec")
    if fraction > 1 or radius < 0:
        raise ValueError("peak fraction must be <= 1 and minimum radius >= 0")
    if (
        isinstance(max_candidates, bool)
        or not isinstance(max_candidates, int)
        or max_candidates < 1
    ):
        raise ValueError("max_candidates must be a positive integer")
    grid = ImageGrid(values.shape, pixel)
    coordinates = grid.coordinates()
    maxima = values == maximum_filter(values, size=3, mode="constant", cval=-np.inf)
    maxima &= values >= fraction * values.max()
    maxima &= np.linalg.norm(coordinates, axis=-1) >= radius
    maxima[[0, -1], :] = False
    maxima[:, [0, -1]] = False
    pixels = list(zip(*np.nonzero(maxima)))
    pixels.sort(key=lambda p: (-values[p], p[0], p[1]))
    candidates = []
    for row, col in pixels:
        offset, refined = _peak_offset(values, row, col)
        position = coordinates[row, col] + np.array([-offset[0], offset[1]]) * pixel
        if np.linalg.norm(position) < radius:
            continue
        if any(
            np.linalg.norm(position - np.array([p["y_arcsec"], p["x_arcsec"]]))
            < separation
            for p in candidates
        ):
            continue
        value, _ = _sample(values, position, pixel)
        candidates.append(
            {
                "rank": len(candidates) + 1,
                "y_arcsec": float(position[0]),
                "x_arcsec": float(position[1]),
                "source_cps_pixel": value,
                "peak_fraction": value / float(values.max()),
                "pixel_peak_cps": float(values[row, col]),
                "row": int(row),
                "column": int(col),
                "subpixel_refined": refined,
            }
        )
        if len(candidates) == max_candidates:
            break
    return candidates


def review_mock_positions(
    mock_dir,
    systems_file,
    *,
    subhalo_id,
    output,
    min_peak_fraction=0.3,
    min_separation_arcsec=0.2,
    max_candidates=12,
    max_arc_distance_arcsec=0.1,
):
    """Report peaks and old-point brightness without modifying fit inputs."""
    mock_dir = Path(mock_dir).resolve()
    systems_file = Path(systems_file).resolve()
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    validate_bundle(mock_dir, "mock")
    config = json.loads((mock_dir / "config.json").read_text())
    if config["subhalo_id"] != subhalo_id:
        raise ValueError("mock subhalo differs from requested system")
    system = load_system_inputs(systems_file, subhalo_id)
    # Match the actual fitting loader, including flip_for_ds9. Raw FITS row
    # coordinates are not necessarily the native PyAutoLens (y, x) frame.
    from autoconf import conf

    config_dir = resource_path("config")
    conf.instance.push(new_path=config_dir)
    import autolens as al

    pixel = system.pixel_scale_arcsec
    source = np.asarray(
        al.Array2D.from_fits(
            file_path=mock_dir / "source_convolved.fits", pixel_scales=pixel
        ).native,
        dtype=float,
    )
    data = np.asarray(
        al.Array2D.from_fits(
            file_path=mock_dir / "data.fits", pixel_scales=pixel
        ).native,
        dtype=float,
    )
    loader = {
        "method": "PyAutoLens Array2D.from_fits, matching production Imaging.from_fits",
        "autolens_version": version("autolens"),
        "flip_for_ds9": bool(conf.instance["general"]["fits"]["flip_for_ds9"]),
        "fits_configuration_sha256": checksum(config_dir / "general.yaml"),
        "coordinates": "native fitting grid: y decreases with row; display origin='upper'",
    }
    if pixel != config["pixel_scale_arcsec"] or data.shape != source.shape:
        raise ValueError("mock source/data grid differs from system configuration")
    distance_limit = scalar(
        max_arc_distance_arcsec, "max_arc_distance_arcsec", positive=True
    )
    candidates = find_arc_peaks(
        source,
        pixel_scale_arcsec=pixel,
        min_peak_fraction=min_peak_fraction,
        min_separation_arcsec=min_separation_arcsec,
        max_candidates=max_candidates,
        min_radius_arcsec=config["central_deflection_radius_arcsec"],
    )
    for candidate in candidates:
        candidate["inside_fitting_mask"] = bool(
            np.hypot(candidate["y_arcsec"], candidate["x_arcsec"])
            <= system.mask_radius_arcsec
        )
    coordinates = ImageGrid(source.shape, pixel).coordinates()
    bright = coordinates[source >= min_peak_fraction * source.max()]
    peak = float(source.max())
    rows = []
    for index, point in enumerate(system.positions_yx_arcsec):
        value, inside = _sample(source, point, pixel)
        bright_distance = float(np.min(np.linalg.norm(bright - point, axis=-1)))
        within_mask = bool(np.linalg.norm(point) <= system.mask_radius_arcsec)
        passes = (
            inside
            and within_mask
            and value / peak >= min_peak_fraction
            and bright_distance <= distance_limit
        )
        closest = (
            min(
                candidates,
                key=lambda p: np.linalg.norm(
                    np.array([p["y_arcsec"], p["x_arcsec"]]) - point
                ),
            )
            if candidates
            else None
        )
        flipped, _ = _sample(source, (-point[0], point[1]), pixel)
        rows.append(
            {
                "point": index + 1,
                "y_arcsec": point[0],
                "x_arcsec": point[1],
                "source_peak_fraction": value / peak,
                "distance_to_bright_arc_arcsec": bright_distance,
                "inside_image": inside,
                "inside_fitting_mask": within_mask,
                "passes_brightness_check": bool(passes),
                "flipped_y_source_peak_fraction": flipped / peak,
                "nearest_candidate_rank": closest["rank"] if closest else None,
                "nearest_candidate_distance_arcsec": float(
                    np.linalg.norm(
                        np.array([closest["y_arcsec"], closest["x_arcsec"]]) - point
                    )
                )
                if closest
                else None,
            }
        )
    possible_y_sign = bool(rows) and all(
        not row["passes_brightness_check"]
        and row["flipped_y_source_peak_fraction"] >= min_peak_fraction
        for row in rows
    )
    report = {
        "schema_version": 2,
        "subhalo_id": subhalo_id,
        "image_loader": loader,
        "coordinate_order": "y,x arcsec; y increases towards upper image rows",
        "candidate_source": "noiseless PSF-convolved lensed source component",
        "settings": {
            "min_peak_fraction": min_peak_fraction,
            "min_separation_arcsec": min_separation_arcsec,
            "max_candidates": max_candidates,
            "max_arc_distance_arcsec": distance_limit,
            "excluded_central_radius_arcsec": config[
                "central_deflection_radius_arcsec"
            ],
        },
        "candidates": candidates,
        "suggested_positions_yx_arcsec": [
            [p["y_arcsec"], p["x_arcsec"]]
            for p in candidates
            if p["inside_fitting_mask"]
        ],
        "catalogue_points": rows,
        "all_catalogue_points_on_bright_arc": bool(rows)
        and all(row["passes_brightness_check"] for row in rows),
        "possible_y_sign_mismatch": possible_y_sign,
        "catalogue_modified": False,
        "fit_submitted": False,
        "limitations": "Bright maxima need review: several may belong to one continuous arc; brightness alone does not establish a common source-plane position.",
        "input_sha256": {
            "mock_manifest": checksum(mock_dir / "manifest.json"),
            "source_convolved.fits": checksum(mock_dir / "source_convolved.fits"),
            "systems": checksum(systems_file),
        },
        "method_source_sha256": checksum(__file__),
    }
    output.mkdir(parents=True)
    write_json(output / "report.json", report)
    for name, records in [
        ("candidates.csv", candidates),
        ("catalogue_points.csv", rows),
    ]:
        with (output / name).open("w", newline="") as stream:
            fields = list(records[0]) if records else ["y_arcsec", "x_arcsec"]
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(records)
    _plot_review(data, source, coordinates, report, output / "positions.png", pixel)
    return output


def _plot_review(data, source, coordinates, report, path, pixel):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from astropy.visualization import ImageNormalize, AsinhStretch

    points = report["catalogue_points"]
    candidates = report["candidates"]
    radii = [np.hypot(p["y_arcsec"], p["x_arcsec"]) for p in points + candidates]
    half = min(min(source.shape) * pixel / 2, max(radii, default=0.8) + 0.35)
    window = np.all(np.abs(coordinates) <= half, axis=-1)
    bounds = [
        -source.shape[1] * pixel / 2,
        source.shape[1] * pixel / 2,
        -source.shape[0] * pixel / 2,
        source.shape[0] * pixel / 2,
    ]
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), layout="constrained")
    for ax, image, title in zip(
        axes, [source, data], ["PSF lensed source", "Noisy mock"]
    ):
        vmax = source.max() if image is source else np.percentile(image[window], 99.8)
        norm = ImageNormalize(
            vmin=0,
            vmax=max(float(vmax), 1e-15),
            stretch=AsinhStretch(0.07 if image is source else 0.025),
        )
        ax.imshow(
            image,
            origin="upper",
            extent=bounds,
            cmap="magma",
            norm=norm,
            interpolation="nearest",
        )
        for p in candidates:
            ax.plot(p["x_arcsec"], p["y_arcsec"], "o", mfc="none", mec="cyan", ms=8)
            ax.annotate(
                f"C{p['rank']}",
                (p["x_arcsec"], p["y_arcsec"]),
                xytext=(5, -12),
                textcoords="offset points",
                color="cyan",
                fontsize=8,
            )
        for p in points:
            color = "white" if p["passes_brightness_check"] else "red"
            ax.plot(p["x_arcsec"], p["y_arcsec"], "+", color=color, ms=12, mew=2)
            ax.annotate(
                f"P{p['point']}",
                (p["x_arcsec"], p["y_arcsec"]),
                xytext=(5, 5),
                textcoords="offset points",
                color=color,
                fontsize=8,
            )
        ax.set(
            xlim=(-half, half),
            ylim=(-half, half),
            xlabel="x [arcsec]",
            ylabel="y [arcsec]",
            title=title,
        )
    fig.suptitle(
        f"ID {report['subhalo_id']}: cyan = peak candidates; white/red = stored points pass/fail\nCandidates require review; input catalogue unchanged"
    )
    fig.savefig(path, dpi=160)
    plt.close(fig)
