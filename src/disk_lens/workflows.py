"""SPH light-map and validated mock-observation file workflows."""

from importlib.metadata import PackageNotFoundError, version
import json
from pathlib import Path

import numpy as np
from astropy.cosmology import Planck18
from astropy.io import fits

from ._io import checksum, write_json
from .grid import ImageGrid
from .lightmaps import aperture_grid, build_lightmap
from .particles import StellarParticles


def run_lightmap(config_file, output):
    """Render a new SPH foreground bundle and record its scientific inputs."""
    output = Path(output)
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    path = Path(config_file).resolve()
    config = json.loads(path.read_text())
    required = {
        "schema_version",
        "input_particles",
        "grid",
        "lens_redshift",
        "zeropoint",
        "cosmology",
        "view",
        "aperture_halfmass_multiple",
        "deposition",
    }
    if (
        not isinstance(config, dict)
        or set(config) != required
        or config["schema_version"] != 1
        or config["cosmology"] != "Planck18"
    ):
        raise ValueError("invalid production SPH lightmap configuration")
    if not isinstance(config["grid"], dict) or set(config["grid"]) not in (
        {"shape", "pixel_scale_arcsec"},
        {"shape", "pixel_scale_arcsec", "render_shape"},
    ):
        raise ValueError("grid requires shape and pixel_scale_arcsec")
    particle_path = (path.parent / config["input_particles"]).resolve()
    particles = StellarParticles.load(particle_path)
    if config["grid"]["shape"] is None:
        grid = aperture_grid(
            particles.halfmass_radius_kpc,
            pixel_scale_arcsec=config["grid"]["pixel_scale_arcsec"],
            lens_redshift=config["lens_redshift"],
            cosmology=Planck18,
            aperture_halfmass_multiple=config["aperture_halfmass_multiple"],
        )
    else:
        grid = ImageGrid(
            tuple(config["grid"]["shape"]), config["grid"]["pixel_scale_arcsec"]
        )
    render_grid = grid
    if "render_shape" in config["grid"]:
        render_grid = ImageGrid(
            tuple(config["grid"]["render_shape"]), grid.pixel_scale_arcsec
        )
        if any(
            output_size < render_size or (output_size - render_size) % 2
            for output_size, render_size in zip(grid.shape, render_grid.shape)
        ):
            raise ValueError("foreground padding must be non-negative and symmetric")
    image, metadata = build_lightmap(
        particles,
        render_grid,
        lens_redshift=config["lens_redshift"],
        zeropoint=config["zeropoint"],
        cosmology=Planck18,
        view=config["view"],
        aperture_halfmass_multiple=config["aperture_halfmass_multiple"],
        deposition=config["deposition"],
    )
    if render_grid.shape != grid.shape:
        padding = tuple(
            ((output_size - render_size) // 2,) * 2
            for output_size, render_size in zip(grid.shape, render_grid.shape)
        )
        image = np.pad(image, padding, mode="constant")
        metadata["image_padding"] = {
            "method": "symmetric zero padding before deflection calculation",
            "render_shape": list(render_grid.shape),
            "output_shape": list(grid.shape),
            "padding_yx": [list(pair) for pair in padding],
            "pixel_centres_preserved": True,
        }
    packages = {}
    for name in ["numpy", "scipy", "astropy", "numba"]:
        try:
            packages[name] = version(name)
        except PackageNotFoundError:
            packages[name] = None
    output.mkdir(parents=True)
    write_json(output / "config.json", config)
    fits.writeto(output / "light.fits", image, header=grid.header(bunit="electron / s"))
    metadata.update(
        {
            "schema_version": 1,
            "stage": "lightmap",
            "status": "complete",
            "input_particles_sha256": checksum(particle_path),
            "files_sha256": {"light.fits": checksum(output / "light.fits")},
            "grid": {
                "shape": list(grid.shape),
                "pixel_scale_arcsec": grid.pixel_scale_arcsec,
            },
            "provenance": {
                "config_sha256": checksum(path),
                "packages": packages,
                "source_sha256": {
                    p.name: checksum(p)
                    for p in sorted(Path(__file__).parent.glob("*.py"))
                },
            },
        }
    )
    write_json(output / "manifest.json", metadata)
    return output
