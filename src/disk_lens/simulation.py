"""SPH-light mock observations using the validated historical numerical settings."""

import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import zoom
from scipy.signal import convolve2d

from ._validation import finite_array, scalar
from .grid import ImageGrid
from ._io import checksum, write_json, validate_bundle
from ._environment import (
    LEGACY_MOCK_PACKAGES,
    MOCK_OVER_SAMPLE_SIZE,
    REQUIRED_PACKAGES,
    package_versions,
    require_environment,
)


def _require_keys(mapping, expected, name):
    if not isinstance(mapping, dict) or set(mapping) != set(expected):
        raise ValueError(f"{name} has missing or unexpected fields")


def centre_crop(array, shape):
    """Return the historical integer-centred crop."""
    values = finite_array(array, "array")
    target = tuple(int(value) for value in shape)
    if len(target) != 2 or any(value < 2 for value in target):
        raise ValueError("shape must contain two integers >= 2")
    if target[0] > values.shape[0] or target[1] > values.shape[1]:
        raise ValueError("ray-tracing grid cannot exceed the deflection grid")
    row = (values.shape[0] - target[0]) // 2
    col = (values.shape[1] - target[1]) // 2
    return values[row : row + target[0], col : col + target[1], ...]


def centre_embed(background, addition):
    """Add a smaller image to the historical integer-centred slice."""
    result = finite_array(background, "background").copy()
    small = finite_array(addition, "addition")
    if result.ndim != 2 or small.ndim != 2:
        raise ValueError("background and addition must be two-dimensional")
    if small.shape[0] > result.shape[0] or small.shape[1] > result.shape[1]:
        raise ValueError("ray-traced source cannot exceed the foreground image")
    row = (result.shape[0] - small.shape[0]) // 2
    col = (result.shape[1] - small.shape[1]) // 2
    result[row : row + small.shape[0], col : col + small.shape[1]] += small
    return result


def gaussian_psf(*, sigma_arcsec, pixel_scale_arcsec, size=15):
    """Gaussian sampled at pixel centres, as in ``sim_util.gauss_2d``."""
    sigma = scalar(sigma_arcsec, "sigma_arcsec", positive=True)
    pixel = scalar(pixel_scale_arcsec, "pixel_scale_arcsec", positive=True)
    if isinstance(size, bool) or not isinstance(size, int) or size < 3 or size % 2 == 0:
        raise ValueError("size must be an odd integer >= 3")
    coordinate = np.arange(size) * pixel
    coordinate -= coordinate.mean()
    x, y = np.meshgrid(coordinate, coordinate)
    psf = np.exp(-0.5 * (x * x + y * y) / sigma**2)
    return psf / psf.sum()


def poisson_observation(expected_image, *, sky_cps_pixel, exposure_seconds, seed=None):
    """Draw Poisson counts and derive noise from that realised count map."""
    expected = finite_array(expected_image, "expected_image")
    sky = scalar(sky_cps_pixel, "sky_cps_pixel", positive=True)
    exposure = scalar(exposure_seconds, "exposure_seconds", positive=True)
    if np.any(expected < 0):
        raise ValueError("expected image must be non-negative")
    if seed is not None and (
        isinstance(seed, bool) or not isinstance(seed, int) or seed < 0
    ):
        raise ValueError("seed must be null or a non-negative integer")
    expected_counts = (expected + sky) * exposure
    if seed is None:
        counts = np.random.poisson(expected_counts)
        rng = "numpy global RandomState; unseeded, matching the notebooks"
    else:
        counts = np.random.RandomState(seed).poisson(expected_counts)
        rng = "numpy.random.RandomState/MT19937"
    return (
        counts / exposure - sky,
        np.sqrt(counts) / exposure,
        {"seed": seed, "rng": rng},
    )


def deflections_from_light(
    light, *, raytrace_shape, pixel_scale_arcsec, central_radius_arcsec
):
    """Reproduce zoom + lenstronomy adaptive kernel + crop + central mask."""
    try:
        from lenstronomy.LensModel.convergence_integrals import (
            deflection_from_kappa_grid_adaptive,
        )
    except ImportError as error:
        raise RuntimeError(
            "SPH mock generation requires the 'mock' dependencies"
        ) from error
    foreground = finite_array(light, "light")
    if foreground.ndim != 2 or foreground.shape[0] != foreground.shape[1]:
        raise ValueError("SPH light map must be a square two-dimensional array")
    pixel = scalar(pixel_scale_arcsec, "pixel_scale_arcsec", positive=True)
    radius = scalar(central_radius_arcsec, "central_radius_arcsec")
    if radius < 0:
        raise ValueError("central radius must be non-negative")
    high_resolution = zoom(foreground, zoom=2, order=1)
    kernel_size = int(len(high_resolution) + 1)
    alpha_x, alpha_y = deflection_from_kappa_grid_adaptive(
        high_resolution,
        grid_spacing=pixel / 2,
        low_res_factor=2,
        high_res_kernel_size=kernel_size,
    )
    alpha_x = centre_crop(alpha_x, raytrace_shape)
    alpha_y = centre_crop(alpha_y, raytrace_shape)
    if radius > 0:
        height, width = alpha_x.shape
        yy, xx = np.meshgrid(np.arange(height), np.arange(width), indexing="ij")
        central = ((yy - (height - 1) / 2) * pixel) ** 2 + (
            (xx - (width - 1) / 2) * pixel
        ) ** 2 <= radius**2
        alpha_x = alpha_x.copy()
        alpha_y = alpha_y.copy()
        alpha_x[central] = 1e6
        alpha_y[central] = 1e6
    return np.stack((np.flipud(alpha_y), np.flipud(alpha_x)), axis=-1)


def mock_tracer(
    stellar_deflections_yx,
    *,
    pixel_scale_arcsec,
    stellar_mass_normalization,
    halo,
    source,
):
    """Build the exact historical stellar-template plus Ludlow-halo tracer."""
    try:
        import autoarray as aa
        import autolens as al
        from autogalaxy.profiles.mass.abstract.abstract import MassProfile
    except ImportError as error:
        raise RuntimeError(
            "SPH mock generation requires the 'mock' dependencies"
        ) from error
    deflections = finite_array(stellar_deflections_yx, "stellar_deflections_yx")
    if deflections.ndim != 3 or deflections.shape[-1] != 2:
        raise ValueError("stellar deflections must have shape (ny,nx,2)")
    grid_size = deflections.shape[0]
    if deflections.shape[1] != grid_size:
        raise ValueError("historical ray-tracing grid must be square")
    grid = al.Grid2D.uniform(
        shape_native=(grid_size, grid_size),
        pixel_scales=pixel_scale_arcsec,
        over_sample_size=MOCK_OVER_SAMPLE_SIZE,
    )
    native = grid.reshape((grid_size, grid_size, 2))
    y_axis = native[:, 0, 0]
    x_axis = native[0, :, 1]

    class InputDeflections(MassProfile):
        def __init__(self, mass_to_light_ratio):
            super().__init__()
            self.mass_to_light_ratio = mass_to_light_ratio
            zeros = np.zeros(deflections.shape[:2])
            axes = (y_axis, x_axis)
            settings = {"fill_value": 0.0, "bounds_error": False, "method": "linear"}
            self.interp_defl_y = RegularGridInterpolator(
                axes, deflections[..., 0], **settings
            )
            self.interp_defl_x = RegularGridInterpolator(
                axes, deflections[..., 1], **settings
            )
            self.interp_kappa = RegularGridInterpolator(axes, zeros, **settings)
            self.interp_psi = RegularGridInterpolator(axes, zeros, **settings)

        @aa.grid_dec.to_vector_yx
        def convergence_2d_from(self, grid):
            return self.interp_kappa(grid) * self.mass_to_light_ratio

        @aa.grid_dec.to_vector_yx
        def potential_2d_from(self, grid):
            return self.interp_psi(grid) * self.mass_to_light_ratio

        @aa.grid_dec.to_vector_yx
        def deflections_yx_2d_from(self, grid):
            values = np.stack(
                (self.interp_defl_y(grid), self.interp_defl_x(grid)), axis=-1
            )
            return values * self.mass_to_light_ratio

    mass = InputDeflections(
        scalar(stellar_mass_normalization, "stellar_mass_normalization", positive=True)
    )
    dark = al.mp.NFWMCRLudlowSph(
        centre=tuple(halo["centre_yx_arcsec"]),
        mass_at_200=scalar(halo["mass_at_200_msun"], "mass_at_200_msun", positive=True),
        redshift_object=scalar(halo["lens_redshift"], "lens_redshift", positive=True),
        redshift_source=scalar(
            halo["source_redshift"], "source_redshift", positive=True
        ),
    )
    lens = al.Galaxy(redshift=halo["lens_redshift"], mass=mass, dark=dark)
    bulge = al.lp.Sersic(
        centre=tuple(source["centre_yx_arcsec"]),
        ell_comps=al.convert.ell_comps_from(
            axis_ratio=source["axis_ratio"], angle=source["angle_deg"]
        ),
        intensity=source["intensity"],
        effective_radius=source["effective_radius_arcsec"],
        sersic_index=source["sersic_index"],
    )
    source_galaxy = al.Galaxy(redshift=halo["source_redshift"], bulge=bulge)
    return al.Tracer(galaxies=[lens, source_galaxy])


def lensed_source(
    stellar_deflections_yx,
    *,
    pixel_scale_arcsec,
    stellar_mass_normalization,
    halo,
    source,
):
    """Evaluate the historical PyAutoLens tracer on its native pixel grid."""
    import autolens as al

    tracer = mock_tracer(
        stellar_deflections_yx,
        pixel_scale_arcsec=pixel_scale_arcsec,
        stellar_mass_normalization=stellar_mass_normalization,
        halo=halo,
        source=source,
    )
    grid_size = np.asarray(stellar_deflections_yx).shape[0]
    grid = al.Grid2D.uniform(
        shape_native=(grid_size, grid_size),
        pixel_scales=pixel_scale_arcsec,
        # 2024's default used native pixels; 2025 defaults to 4x4 sub-pixels.
        over_sample_size=MOCK_OVER_SAMPLE_SIZE,
    )
    return np.flipud(
        np.asarray(tracer.image_2d_from(grid=grid)).reshape(grid_size, grid_size)
    )


def run_simulation(config_file, output):
    """Generate a production mock from a complete SPH light-map bundle."""
    config_path = Path(config_file).resolve()
    config = json.loads(config_path.read_text())
    required = {
        "schema_version",
        "mode",
        "subhalo_id",
        "input_lightmap",
        "pixel_scale_arcsec",
        "raytrace_grid_shape",
        "stellar_mass_normalization",
        "halo",
        "source",
        "central_deflection_radius_arcsec",
        "psf",
        "observation",
    }
    if not isinstance(config, dict) or set(config) != required:
        raise ValueError("invalid schema-1 sph configuration")
    if config["schema_version"] != 1 or config["mode"] != "sph":
        raise ValueError("invalid schema-1 sph configuration")
    _require_keys(
        config["halo"],
        {"mass_at_200_msun", "centre_yx_arcsec", "lens_redshift", "source_redshift"},
        "halo",
    )
    _require_keys(
        config["source"],
        {
            "centre_yx_arcsec",
            "axis_ratio",
            "angle_deg",
            "intensity",
            "effective_radius_arcsec",
            "sersic_index",
        },
        "source",
    )
    _require_keys(config["psf"], {"sigma_arcsec", "size"}, "psf")
    _require_keys(
        config["observation"],
        {"sky_cps_pixel", "exposure_seconds", "seed"},
        "observation",
    )
    if (
        isinstance(config["subhalo_id"], bool)
        or not isinstance(config["subhalo_id"], int)
        or config["subhalo_id"] < 0
    ):
        raise ValueError("subhalo_id must be a non-negative integer")
    historical_values = {
        "pixel_scale_arcsec": (config["pixel_scale_arcsec"], 0.05),
        "halo.lens_redshift": (config["halo"]["lens_redshift"], 0.15),
        "halo.source_redshift": (config["halo"]["source_redshift"], 0.65),
        "source.axis_ratio": (config["source"]["axis_ratio"], 0.7),
        "source.angle_deg": (config["source"]["angle_deg"], 0.0),
        "source.intensity": (config["source"]["intensity"], 0.4),
        "source.effective_radius_arcsec": (
            config["source"]["effective_radius_arcsec"],
            0.15,
        ),
        "source.sersic_index": (config["source"]["sersic_index"], 1.0),
        "psf.sigma_arcsec": (config["psf"]["sigma_arcsec"], 0.05),
        "psf.size": (config["psf"]["size"], 15),
        "observation.sky_cps_pixel": (config["observation"]["sky_cps_pixel"], 0.1),
        "observation.exposure_seconds": (
            config["observation"]["exposure_seconds"],
            2400,
        ),
    }
    changed = [
        name
        for name, (actual, expected) in historical_values.items()
        if not np.isclose(actual, expected)
    ]
    if changed:
        raise ValueError(
            "sph fixes the paper-wide settings; changed: " + ", ".join(changed)
        )
    # The central pseudo-core exclusion was adjusted by hand for some systems.
    # It is an explicit per-system input, not a paper-wide constant.
    central_radius = config["central_deflection_radius_arcsec"]
    if (
        isinstance(central_radius, bool)
        or not isinstance(central_radius, (int, float))
        or not np.isfinite(central_radius)
        or central_radius < 0
    ):
        raise ValueError(
            "central_deflection_radius_arcsec must be finite and non-negative"
        )
    if tuple(config["halo"]["centre_yx_arcsec"]) != (0.0, 0.0):
        raise ValueError("sph fixes halo.centre_yx_arcsec to [0.0, 0.0]")
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    bundle = (config_path.parent / config["input_lightmap"]).resolve()
    light_manifest = validate_bundle(bundle, "lightmap")
    if light_manifest.get("deposition_parameters", {}).get("method") != "sph":
        raise ValueError("production mock requires an SPH light-map bundle")
    if light_manifest["grid"]["pixel_scale_arcsec"] != config["pixel_scale_arcsec"]:
        raise ValueError("SPH light and mock pixel scales differ")
    if light_manifest["lens_redshift"] != config["halo"]["lens_redshift"]:
        raise ValueError("SPH light and mock lens redshifts differ")
    if (
        light_manifest["zeropoint"] != 25.948
        or light_manifest["cosmology"] != "Planck18"
    ):
        raise ValueError("production mock requires the validated stellar calibration")
    require_environment()
    light_path = bundle / "light.fits"
    light = fits.getdata(light_path).astype(float)
    image_grid = ImageGrid(light.shape, config["pixel_scale_arcsec"])
    raytrace_grid = ImageGrid(
        tuple(config["raytrace_grid_shape"]), config["pixel_scale_arcsec"]
    )
    if raytrace_grid.shape[0] != raytrace_grid.shape[1]:
        raise ValueError("sph requires a square ray-tracing grid")
    if raytrace_grid.shape[0] % 2:
        raise ValueError("sph requires an even ray-tracing grid size")
    if (
        raytrace_grid.shape[0] > image_grid.shape[0]
        or raytrace_grid.shape[1] > image_grid.shape[1]
    ):
        raise ValueError("ray-tracing grid cannot exceed the SPH foreground image")
    halo = dict(config["halo"])
    stellar_deflections = deflections_from_light(
        light,
        raytrace_shape=raytrace_grid.shape,
        pixel_scale_arcsec=image_grid.pixel_scale_arcsec,
        central_radius_arcsec=config["central_deflection_radius_arcsec"],
    )
    source_native = lensed_source(
        stellar_deflections,
        pixel_scale_arcsec=image_grid.pixel_scale_arcsec,
        stellar_mass_normalization=config["stellar_mass_normalization"],
        halo=halo,
        source=config["source"],
    )
    source_full = centre_embed(np.zeros(image_grid.shape), source_native)
    psf = gaussian_psf(
        sigma_arcsec=config["psf"]["sigma_arcsec"],
        pixel_scale_arcsec=image_grid.pixel_scale_arcsec,
        size=config["psf"]["size"],
    )
    foreground_convolved = convolve2d(light, psf, mode="same")
    source_convolved = convolve2d(source_full, psf, mode="same")
    model_image = foreground_convolved + source_convolved
    data, noise_map, rng = poisson_observation(
        model_image,
        sky_cps_pixel=config["observation"]["sky_cps_pixel"],
        exposure_seconds=config["observation"]["exposure_seconds"],
        seed=config["observation"]["seed"],
    )

    output.mkdir(parents=True)
    write_json(output / "config.json", config)
    arrays = {
        "data.fits": data,
        "noise_map.fits": noise_map,
        "psf.fits": psf,
        "model_image.fits": model_image,
        "foreground_convolved.fits": foreground_convolved,
        "source_convolved.fits": source_convolved,
        "source_native.fits": source_native,
    }
    files = {}
    for filename, array in arrays.items():
        if filename == "psf.fits":
            grid_for_header = ImageGrid(array.shape, image_grid.pixel_scale_arcsec)
            unit = "1"
        elif filename == "source_native.fits":
            grid_for_header = raytrace_grid
            unit = "electron / s"
        else:
            grid_for_header = image_grid
            unit = "electron / s"
        fits.writeto(
            output / filename, array, header=grid_for_header.header(bunit=unit)
        )
        files[filename] = checksum(output / filename)
    np.save(output / "stellar_deflections_yx.npy", stellar_deflections)
    files["stellar_deflections_yx.npy"] = checksum(
        output / "stellar_deflections_yx.npy"
    )
    packages = package_versions(
        ["numpy", "scipy", "astropy", "autolens", "lenstronomy"]
    )
    reference_packages = LEGACY_MOCK_PACKAGES
    expected_noise = np.sqrt(
        (model_image + config["observation"]["sky_cps_pixel"])
        / config["observation"]["exposure_seconds"]
    )
    manifest = {
        "schema_version": 1,
        "stage": "mock",
        "mode": "sph",
        "status": "complete",
        "subhalo_id": config["subhalo_id"],
        "required_packages": dict(REQUIRED_PACKAGES),
        "historical_reference_packages": reference_packages,
        "historical_package_match": all(
            packages[name] == required_version
            for name, required_version in reference_packages.items()
        ),
        "input_light_fits": str(light_path),
        "input_manifest_sha256": checksum(bundle / "manifest.json"),
        "deposition": "sph",
        "input_light_sha256": checksum(light_path),
        "image_grid": {
            "shape": list(image_grid.shape),
            "pixel_scale_arcsec": image_grid.pixel_scale_arcsec,
        },
        "raytrace_grid": {
            "shape": list(raytrace_grid.shape),
            "pixel_scale_arcsec": raytrace_grid.pixel_scale_arcsec,
            "over_sample_size": MOCK_OVER_SAMPLE_SIZE,
        },
        "source_peak_snr_expected": float(np.max(source_convolved / expected_noise)),
        "noise": rng
        | {
            "definition": "sqrt(realised Poisson counts) / exposure",
            "sky_cps_pixel": config["observation"]["sky_cps_pixel"],
            "exposure_seconds": config["observation"]["exposure_seconds"],
        },
        "central_deflection": {
            "radius_arcsec": config["central_deflection_radius_arcsec"],
            "replacement_value_arcsec": 1e6,
        },
        "packages": packages,
        "method_source_sha256": {
            name: checksum(Path(__file__).parent / name)
            for name in ["simulation.py", "_environment.py"]
        },
        "config_sha256": checksum(config_path),
        "files_sha256": files,
    }
    write_json(output / "manifest.json", manifest)
    return output
