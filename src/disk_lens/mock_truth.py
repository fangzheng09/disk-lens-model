"""Unmasked effective Einstein radii of the exact production mock lens model."""

import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from skimage.measure import find_contours

from ._io import checksum, validate_bundle, write_json
from .simulation import deflections_from_light, mock_tracer
from ._environment import (
    MOCK_OVER_SAMPLE_SIZE,
    REQUIRED_PACKAGES,
    package_versions,
    require_environment,
    validate_recorded_packages,
)


def _critical_curves(eigenvalue, grid):
    values = np.asarray(eigenvalue.native, dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError("non-finite mock tangential eigenvalues")
    if np.any(values[[0, -1], :] <= 0) or np.any(values[:, [0, -1]] <= 0):
        raise ValueError("mock tangential critical region reaches evaluation boundary")
    pixel = float(grid.pixel_scales[0])
    origin = np.asarray(grid.native)[0, 0]
    curves = []
    for contour in find_contours(values, 0, fully_connected="low"):
        if not np.allclose(contour[0], contour[-1], atol=1e-8):
            raise ValueError("open mock tangential critical curve")
        row, col = contour.T
        area = abs(np.dot(col[:-1], row[1:]) - np.dot(row[:-1], col[1:])) / 2
        if area > 0:
            curves.append(
                (
                    float(area * pixel**2),
                    np.column_stack((origin[0] - row * pixel, origin[1] + col * pixel)),
                )
            )
    curves.sort(key=lambda item: item[0], reverse=True)
    areas = [item[0] for item in curves]
    radii = np.sqrt(np.asarray(areas) / np.pi)
    return {
        "pixel_scale_arcsec": pixel,
        "shape": list(values.shape),
        "einstein_radius_arcsec": float(radii[0]) if len(radii) else 0.0,
        "critical_curve_areas_arcsec2": areas,
        "critical_curve_count": len(areas),
        "historical_sum_of_radii_arcsec": float(radii.sum()),
        "boundary_clear": True,
        "largest_curve_yx_arcsec": curves[0][1].tolist() if curves else [],
    }


def run_mock_truth(config_file, mock_bundle, output):
    """Rebuild radius-zero stellar deflections and measure the actual mock tracer.

    Finer grids sample the same bilinearly interpolated native deflections;
    they do not change the mock's SPH resolution, mass model or halo.
    """
    require_environment()

    import autolens as al

    config_path = Path(config_file).resolve()
    config = json.loads(config_path.read_text())
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    mock_bundle = Path(mock_bundle).resolve()
    mock = validate_bundle(mock_bundle, "mock", config_file=config_path)
    validate_recorded_packages(mock, REQUIRED_PACKAGES, stage="mock")
    if (
        mock.get("required_packages") != REQUIRED_PACKAGES
        or mock.get("raytrace_grid", {}).get("over_sample_size")
        != MOCK_OVER_SAMPLE_SIZE
    ):
        raise ValueError("mock truth requires a verified native-sampled 2025 mock")
    light_bundle = (config_path.parent / config["input_lightmap"]).resolve()
    validate_bundle(light_bundle, "lightmap")
    light_path = light_bundle / "light.fits"
    light_hash = checksum(light_path)
    if mock["input_manifest_sha256"] != checksum(light_bundle / "manifest.json"):
        raise ValueError(
            "mock input light-map manifest differs from the truth template"
        )
    if mock["input_light_sha256"] != light_hash:
        raise ValueError("mock input light differs from the truth template")
    light = fits.getdata(light_path).astype(float)
    shape = tuple(config["raytrace_grid_shape"])
    native_pixel = float(config["pixel_scale_arcsec"])
    unmasked = deflections_from_light(
        light,
        raytrace_shape=shape,
        pixel_scale_arcsec=native_pixel,
        central_radius_arcsec=0.0,
    )
    tracer = mock_tracer(
        unmasked,
        pixel_scale_arcsec=native_pixel,
        stellar_mass_normalization=config["stellar_mass_normalization"],
        halo=config["halo"],
        source=config["source"],
    )
    # Stay within the interpolation domain, leaving a native-pixel edge margin.
    width = (shape[0] - 2) * native_pixel
    checks = []
    for pixel in (native_pixel, native_pixel / 2, native_pixel / 4):
        count = int(round(width / pixel))
        grid = al.Grid2D.uniform(
            shape_native=(count, count),
            pixel_scales=pixel,
            over_sample_size=MOCK_OVER_SAMPLE_SIZE,
        )
        checks.append(
            _critical_curves(tracer.tangential_eigen_value_from(grid=grid), grid)
        )
    radius = checks[-1]["einstein_radius_arcsec"]
    tolerance = max(0.005, 0.01 * radius)
    delta = abs(radius - checks[-2]["einstein_radius_arcsec"])
    converged = delta <= tolerance
    if not converged:
        pixel = native_pixel / 8
        count = int(round(width / pixel))
        grid = al.Grid2D.uniform(
            shape_native=(count, count),
            pixel_scales=pixel,
            over_sample_size=MOCK_OVER_SAMPLE_SIZE,
        )
        checks.append(
            _critical_curves(tracer.tangential_eigen_value_from(grid=grid), grid)
        )
        radius = checks[-1]["einstein_radius_arcsec"]
        tolerance = max(0.005, 0.01 * radius)
        delta = abs(radius - checks[-2]["einstein_radius_arcsec"])
        converged = delta <= tolerance
    if not converged:
        raise ValueError(f"mock critical-curve sampling not converged: delta={delta}")
    packages = package_versions(
        ["numpy", "scipy", "astropy", "autolens", "lenstronomy"]
    )
    truth = {
        "schema_version": 1,
        "subhalo_id": config["subhalo_id"],
        "einstein_radius_arcsec": radius,
        "definition": "sqrt(A/pi) for largest closed tangential critical curve",
        "historical_sum_of_radii_arcsec": checks[-1]["historical_sum_of_radii_arcsec"],
        "central_deflection_radius_arcsec": 0.0,
        "image_central_deflection_radius_arcsec": config[
            "central_deflection_radius_arcsec"
        ],
        "mass_model": "actual mock SPH stellar template coefficient + analytic Ludlow NFW",
        "stellar_mass_normalization": config["stellar_mass_normalization"],
        "halo": config["halo"],
        "halo_cosmology": (
            "NFWMCRLudlowSph default in PyAutoLens " + packages["autolens"]
        ),
        "native_pixel_scale_arcsec": native_pixel,
        "native_raytrace_shape": list(shape),
        "over_sample_size": MOCK_OVER_SAMPLE_SIZE,
        "convergence": {
            "passed": converged,
            "last_grid_change_arcsec": delta,
            "tolerance_arcsec": tolerance,
            "scope": "critical-curve sampling of the fixed native interpolated mock model",
        },
        "evaluations": checks,
        "note": "No PSF, noise, source-suppression or fitting mask enters the truth. This is not the particle-mass selection radius.",
    }
    output.mkdir(parents=True)
    write_json(output / "config.json", config)
    write_json(output / "truth.json", truth)
    np.save(output / "unit_stellar_deflections_unmasked_yx.npy", unmasked)
    names = ["truth.json", "unit_stellar_deflections_unmasked_yx.npy"]
    write_json(
        output / "manifest.json",
        {
            "schema_version": 1,
            "stage": "mock_truth",
            "status": "complete",
            "subhalo_id": config["subhalo_id"],
            "required_packages": dict(REQUIRED_PACKAGES),
            "packages": packages,
            "over_sample_size": MOCK_OVER_SAMPLE_SIZE,
            "config_sha256": checksum(config_path),
            "input_mock_manifest_sha256": checksum(mock_bundle / "manifest.json"),
            "input_manifest_sha256": checksum(light_bundle / "manifest.json"),
            "input_light_sha256": light_hash,
            "files_sha256": {name: checksum(output / name) for name in names},
            "method_source_sha256": {
                name: checksum(Path(__file__).parent / name)
                for name in ["mock_truth.py", "simulation.py", "_environment.py"]
            },
            "autolens_version": packages["autolens"],
            "lenstronomy_version": packages["lenstronomy"],
        },
    )
    return output
