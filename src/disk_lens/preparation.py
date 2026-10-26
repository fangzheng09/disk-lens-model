"""Prepare and run the adopted SPH → mock → five-stage fit in one environment."""

import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
from astropy.cosmology import Planck18

from ._fitting_config import load_fitting_config
from ._environment import (
    LEGACY_MOCK_PACKAGES,
    MOCK_OVER_SAMPLE_SIZE,
    REQUIRED_PACKAGES,
    validate_recorded_packages,
)
from ._io import checksum, validate_bundle, write_json
from ._validation import scalar
from .lightmaps import SPH_SETTINGS, aperture_grid
from .systems import load_system_inputs

MOCK_PACKAGES = dict(REQUIRED_PACKAGES)
FIT_PACKAGES = dict(REQUIRED_PACKAGES)


def _relative(target, parent):
    return os.path.relpath(Path(target).resolve(), Path(parent).resolve())


def prepare_main_run(
    *,
    systems_file,
    subhalo_id,
    particles_file,
    fitting_template,
    output_dir,
    seed,
    central_deflection_radius_arcsec=None,
    stellar_mass_normalization=None,
    mass_at_200_msun=None,
    sph_workers=1,
    sph_subpixels=8,
):
    """Resolve the checked system record into one immutable SPH production plan."""
    systems = Path(systems_file).resolve()
    particles = Path(particles_file).resolve()
    fitting = Path(fitting_template).resolve()
    output = Path(output_dir).resolve()
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    if not particles.is_file():
        raise FileNotFoundError(particles)
    _, fitting_config = load_fitting_config(fitting)
    system = load_system_inputs(systems, subhalo_id)
    if system.source_centre_yx_arcsec is None or not system.positions_yx_arcsec:
        raise ValueError(
            "the system record must provide a source centre and image positions"
        )
    if system.pixel_scale_arcsec != 0.05:
        raise ValueError("the validated historical mock uses 0.05 arcsec pixels")
    if system.grid_shape[0] != system.grid_shape[1] or system.grid_shape[0] % 2:
        raise ValueError("the ray-tracing grid must be square and even-sized")
    for name, value in [("sph_workers", sph_workers), ("sph_subpixels", sph_subpixels)]:
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"{name} must be a positive integer")
    if isinstance(seed, bool) or not isinstance(seed, int) or not 0 <= seed < 2**32:
        raise ValueError("seed must be an integer between 0 and 2**32-1")
    adopted_radius = (
        system.central_deflection_radius_arcsec
        if central_deflection_radius_arcsec is None
        else central_deflection_radius_arcsec
    )
    if adopted_radius is None:
        raise ValueError(
            "provide central_deflection_radius_arcsec in the system record or explicitly"
        )
    central_radius = scalar(adopted_radius, "central_deflection_radius_arcsec")
    if central_radius < 0:
        raise ValueError("central_deflection_radius_arcsec must be non-negative")
    normalization = scalar(
        system.adopted_stellar_mass_normalization
        if stellar_mass_normalization is None
        else stellar_mass_normalization,
        "stellar_mass_normalization",
        positive=True,
    )
    halo_mass = scalar(
        system.halo_mass_at_200_msun_legacy
        if mass_at_200_msun is None
        else mass_at_200_msun,
        "mass_at_200_msun",
        positive=True,
    )
    # Read only the small radius entry, not millions of particles during preparation.
    with np.load(particles, allow_pickle=False) as archive:
        halfmass = float(archive["halfmass_radius_kpc"])
    foreground = aperture_grid(
        halfmass,
        pixel_scale_arcsec=0.05,
        lens_redshift=fitting_config["lens_redshift"],
        cosmology=Planck18,
    )
    render_grid = foreground
    if foreground.shape[0] < system.grid_shape[0]:
        raise ValueError(
            f"system {subhalo_id}: mock grid {system.grid_shape} exceeds the "
            f"nominal SPH foreground {foreground.shape}; check the historical "
            "mock crop rather than the fitting/plotting grid"
        )
    config_dir = output / "config"
    lightmap = {
        "schema_version": 1,
        "input_particles": _relative(particles, config_dir),
        "grid": {"shape": list(foreground.shape), "pixel_scale_arcsec": 0.05},
        "lens_redshift": 0.15,
        "zeropoint": 25.948,
        "cosmology": "Planck18",
        "view": "yz",
        "aperture_halfmass_multiple": 5.0,
        "deposition": dict(SPH_SETTINGS, workers=sph_workers, subpixels=sph_subpixels),
    }
    mock = {
        "schema_version": 1,
        "mode": "sph",
        "subhalo_id": int(subhalo_id),
        "input_lightmap": "../lightmap",
        "pixel_scale_arcsec": 0.05,
        "raytrace_grid_shape": list(system.grid_shape),
        "stellar_mass_normalization": normalization,
        "halo": {
            "mass_at_200_msun": halo_mass,
            "centre_yx_arcsec": [0.0, 0.0],
            "lens_redshift": 0.15,
            "source_redshift": 0.65,
        },
        "source": {
            "centre_yx_arcsec": list(system.source_centre_yx_arcsec),
            "axis_ratio": 0.7,
            "angle_deg": 0.0,
            "intensity": 0.4,
            "effective_radius_arcsec": 0.15,
            "sersic_index": 1.0,
        },
        "central_deflection_radius_arcsec": central_radius,
        "psf": {"sigma_arcsec": 0.05, "size": 15},
        "observation": {"sky_cps_pixel": 0.1, "exposure_seconds": 2400, "seed": seed},
    }
    config_dir.mkdir(parents=True)
    for name, content in [
        ("lightmap", lightmap),
        ("mock", mock),
        ("fitting", fitting_config),
    ]:
        write_json(config_dir / f"{name}.json", content)
    plan = {
        "schema_version": 2,
        "workflow": "sph",
        "status": "prepared",
        "subhalo_id": int(subhalo_id),
        "input_system_status": system.status,
        "required_packages": dict(REQUIRED_PACKAGES),
        "paths": {
            "systems": _relative(systems, output),
            "particles": _relative(particles, output),
            "lightmap_config": "config/lightmap.json",
            "mock_config": "config/mock.json",
            "fitting_config": "config/fitting.json",
            "lightmap_output": "lightmap",
            "mock_output": "mock",
            "mock_truth_output": "truth",
            "fit_output": "fit",
            "summary_output": "summary/parameters.json",
        },
        "input_sha256": {
            "systems": checksum(systems),
            "particles": checksum(particles),
            **{
                name + "_config": checksum(config_dir / f"{name}.json")
                for name in ["lightmap", "mock", "fitting"]
            },
        },
        "resolved_parameters": {
            "stellar_mass_normalization": normalization,
            "stellar_normalization_source": (
                "explicit override"
                if stellar_mass_normalization is not None
                else "systems.stellar_mass_normalization"
                if system.stellar_mass_normalization is not None
                else "systems.stellar_mass_normalization_legacy"
            ),
            "mass_at_200_msun": halo_mass,
            "central_deflection_radius_arcsec": central_radius,
            "central_deflection_radius_source": (
                "systems.central_deflection_radius_arcsec"
                if central_deflection_radius_arcsec is None
                else "explicit override"
            ),
            "seed": seed,
            "deposition": lightmap["deposition"],
            "foreground_shape": list(foreground.shape),
            "foreground_render_shape": list(render_grid.shape),
            "raytrace_shape": list(system.grid_shape),
            "raytrace_over_sample_size": MOCK_OVER_SAMPLE_SIZE,
        },
    }
    path = output / "plan.json"
    write_json(path, plan)
    return path


def _load_plan(plan_file):
    path = Path(plan_file).resolve()
    plan = json.loads(path.read_text())
    if (
        not isinstance(plan, dict)
        or plan.get("schema_version") != 2
        or plan.get("workflow") != "sph"
    ):
        raise ValueError(
            "expected a production SPH schema-2 plan; prepare a new output directory"
        )
    if "required_packages" in plan and plan["required_packages"] != REQUIRED_PACKAGES:
        raise ValueError(
            "prepared plan required_packages differs from the supported workflow "
            f"{REQUIRED_PACKAGES}; prepare a new output directory"
        )
    if "required_packages" in plan and (
        plan.get("resolved_parameters", {}).get("raytrace_over_sample_size")
        != MOCK_OVER_SAMPLE_SIZE
    ):
        raise ValueError("prepared plan requires the preserved native mock sampling")
    paths = {
        key: (path.parent / value).resolve() for key, value in plan["paths"].items()
    }
    for name in [
        "systems",
        "particles",
        "lightmap_config",
        "mock_config",
        "fitting_config",
    ]:
        if checksum(paths[name]) != plan["input_sha256"][name]:
            raise ValueError(f"prepared input or configuration changed: {name}")
    load_fitting_config(paths["fitting_config"])
    if (
        json.loads(paths["lightmap_config"].read_text())["deposition"]["method"]
        != "sph"
    ):
        raise ValueError("production plans require SPH light maps")
    return path, plan, paths


def _check_interpreter(executable, packages):
    """Check environment pins without importing the modeling stack or starting work."""
    code = (
        "from importlib.metadata import version\n"
        "import json, sys\n"
        "required=json.loads(sys.argv[1])\n"
        "actual={name:version(name) for name in required}\n"
        "if actual != required:\n"
        "    raise RuntimeError(f'wrong workflow environment: {actual}; expected {required}')\n"
    )
    subprocess.run([str(executable), "-c", code, json.dumps(packages)], check=True)


def _complete(directory, stage, config_file, input_hashes, *, packages=None):
    if not directory.exists():
        return False
    manifest = validate_bundle(
        directory, stage, config_file=config_file, input_sha256=input_hashes
    )
    if packages is not None:
        validate_recorded_packages(manifest, packages, stage=stage)
        if packages == REQUIRED_PACKAGES:
            recorded_sampling = (
                manifest.get("raytrace_grid", {}).get("over_sample_size")
                if stage == "mock"
                else manifest.get("over_sample_size")
            )
            if (
                manifest.get("required_packages") != REQUIRED_PACKAGES
                or recorded_sampling != MOCK_OVER_SAMPLE_SIZE
            ):
                raise ValueError(
                    f"{stage} bundle requirements or native sampling differ from "
                    "the prepared workflow; prepare a new output directory"
                )
    return True


def _prepared_products(plan_path, plan, paths):
    """Check every existing bundle before any computation or interpreter dispatch."""
    legacy = "required_packages" not in plan
    packages = LEGACY_MOCK_PACKAGES if legacy else plan["required_packages"]
    lightmap_complete = _complete(
        paths["lightmap_output"],
        "lightmap",
        paths["lightmap_config"],
        {"input_particles_sha256": plan["input_sha256"]["particles"]},
    )
    truth_path = paths.get("mock_truth_output", plan_path.parent / "truth")
    if not lightmap_complete:
        if legacy or paths["mock_output"].exists() or truth_path.exists():
            raise ValueError(
                "plan requires a verified lightmap, mock and truth chain; "
                "prepare a new output directory"
            )
        return False, False, False, truth_path
    hashes = {
        "input_manifest_sha256": checksum(paths["lightmap_output"] / "manifest.json"),
        "input_light_sha256": checksum(paths["lightmap_output"] / "light.fits"),
    }
    mock_complete = _complete(
        paths["mock_output"],
        "mock",
        paths["mock_config"],
        hashes,
        packages=packages,
    )
    if not mock_complete:
        if legacy or truth_path.exists():
            raise ValueError(
                "historical or incomplete mock chain cannot be regenerated in "
                "place; prepare a new output directory"
            )
        return True, False, False, truth_path
    truth_complete = _complete(
        truth_path,
        "mock_truth",
        paths["mock_config"],
        {
            **hashes,
            "input_mock_manifest_sha256": checksum(
                paths["mock_output"] / "manifest.json"
            ),
        },
        packages=packages,
    )
    if legacy and not truth_complete:
        raise ValueError(
            "historical plan has no required_packages and may only reuse a "
            "complete verified 2024 mock and truth; prepare a new output directory"
        )
    return lightmap_complete, mock_complete, truth_complete, truth_path


def run_prepared_main(
    plan_file,
    *,
    through="summary",
    number_of_cores=1,
    iterations_per_update=5000,
    n_live_override=None,
    mock_python=None,
    fit_python=None,
):
    """Execute with the current Python; optional interpreter overrides are supported.

    Historical schema-2 plans without requirements may reuse their verified
    2024 mock and truth, but cannot generate missing products in place.
    """
    if through not in {"lightmap", "mock", "fit", "summary"}:
        raise ValueError("through must be lightmap, mock, fit or summary")
    plan_path, plan, paths = _load_plan(plan_file)
    lightmap_complete, mock_complete, truth_complete, truth_path = _prepared_products(
        plan_path, plan, paths
    )
    if mock_python is not None or fit_python is not None:
        mock_exe = str(mock_python or sys.executable)
        fit_exe = str(fit_python or sys.executable)
        if through != "lightmap" and not (mock_complete and truth_complete):
            _check_interpreter(mock_exe, MOCK_PACKAGES)
        if through in {"fit", "summary"}:
            _check_interpreter(fit_exe, FIT_PACKAGES)
        if not lightmap_complete or (
            through != "lightmap" and not (mock_complete and truth_complete)
        ):
            subprocess.run(
                [
                    mock_exe,
                    "-m",
                    "disk_lens",
                    "run-main",
                    "--plan",
                    str(plan_path),
                    "--through",
                    "lightmap" if through == "lightmap" else "mock",
                ],
                check=True,
            )
        if through in {"fit", "summary"}:
            command = [
                fit_exe,
                "-m",
                "disk_lens",
                "run-main",
                "--plan",
                str(plan_path),
                "--through",
                through,
                "--cores",
                str(number_of_cores),
                "--iterations-per-update",
                str(iterations_per_update),
            ]
            if n_live_override is not None:
                command += ["--n-live", str(n_live_override)]
            subprocess.run(command, check=True)
        return paths[
            {
                "lightmap": "lightmap_output",
                "mock": "mock_output",
                "fit": "fit_output",
                "summary": "summary_output",
            }[through]
        ]
    # Fail before expensive rendering if any requested computation needs the pins.
    if through in {"fit", "summary"} or (
        through != "lightmap" and not (mock_complete and truth_complete)
    ):
        _check_interpreter(sys.executable, REQUIRED_PACKAGES)
    from .workflows import run_lightmap

    if not lightmap_complete:
        run_lightmap(paths["lightmap_config"], paths["lightmap_output"])
    if through == "lightmap":
        return paths["lightmap_output"]
    if not mock_complete:
        from .simulation import run_simulation

        run_simulation(paths["mock_config"], paths["mock_output"])
    if not truth_complete:
        from .mock_truth import run_mock_truth

        run_mock_truth(paths["mock_config"], paths["mock_output"], truth_path)
    if not (mock_complete and truth_complete):
        _prepared_products(plan_path, plan, paths)
    if through == "mock":
        return paths["mock_output"]
    from .pipeline import run_main_slam

    run_main_slam(
        paths["fitting_config"],
        paths["systems"],
        subhalo_id=plan["subhalo_id"],
        dataset_dir=paths["mock_output"],
        output_dir=paths["fit_output"],
        number_of_cores=number_of_cores,
        iterations_per_update=iterations_per_update,
        n_live_override=n_live_override,
        reuse_completed=True,
    )
    if through == "fit":
        return paths["fit_output"]
    return summarize_main_run(plan_path)


def summarize_main_run(plan_file, parameters=None):
    """Summarize the verified final posterior; infer exact model columns by default."""
    from .analysis import summarize_samples_csv
    from .pipeline import _validate_completed_result
    from ._result_files import result_checksum

    plan_path, plan, paths = _load_plan(plan_file)
    _prepared_products(plan_path, plan, paths)
    fit_manifest = json.loads(
        (paths["fit_output"] / "pipeline_manifest.json").read_text()
    )
    if fit_manifest.get("status") != "complete":
        raise ValueError("the five-stage fit has not completed")
    required_stages = {
        "source_lp",
        "source_pix_1",
        "source_pix_2",
        "light",
        "mass_light_dark",
    }
    if not required_stages.issubset(fit_manifest.get("completed_stages", [])):
        raise ValueError("completed fit is missing stage completion records")
    expected = {
        "subhalo_id": plan["subhalo_id"],
        "config_sha256": plan["input_sha256"]["fitting_config"],
        "systems_sha256": plan["input_sha256"]["systems"],
        "dataset_sha256": {
            name: checksum(paths["mock_output"] / f"{name}.fits")
            for name in ["data", "psf", "noise_map"]
        },
    }
    if any(fit_manifest.get(key) != value for key, value in expected.items()):
        raise ValueError("completed fit inputs differ from the SPH plan")
    _validate_completed_result(paths["fit_output"], fit_manifest)
    samples = (
        paths["fit_output"] / fit_manifest["result_directory"] / "files/samples.csv"
    )
    summary = summarize_samples_csv(samples, paths["summary_output"], parameters)
    write_json(
        paths["summary_output"].parent / "manifest.json",
        {
            "schema_version": 1,
            "stage": "summary",
            "status": "complete",
            "plan_sha256": checksum(plan_path),
            "samples_sha256": result_checksum(samples),
            "parameters": list(summary["parameters"]),
            "files_sha256": {
                paths["summary_output"].name: checksum(paths["summary_output"])
            },
        },
    )
    return paths["summary_output"]
