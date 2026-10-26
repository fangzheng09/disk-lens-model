"""Original-style, five-stage SLaM fitting main program.

The fit body follows the paper-era pipeline.  Paths and per-system inputs are
explicit for fitting production SPH mock FITS.
"""

from importlib.metadata import version
import json
import os
from pathlib import Path

# Workflows imports NumPy/SciPy, so apply the original single-thread setting
# before importing it, including when using this module as the main program.
for _thread_variable in (
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OMP_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread_variable] = "1"

from disk_lens.systems import load_system_inputs
from disk_lens._paths import resource_path
from disk_lens._io import checksum, write_json as _write_manifest
from disk_lens._result_files import result_checksum
from disk_lens._fitting_config import (
    STAGE_LIVE_POINTS as STAGE_LIVE_POINTS,
    live_points as _live_points,
    load_fitting_config as _load_config,
)

SUPPORTED_AUTOLENS_VERSION = "2025.1.18.7"


def fit(
    dataset_name: str,
    mask_radius: float = None,
    number_of_cores: int = 1,
    iterations_per_update: int = 5000,
    *,
    input_dir,
    output_dir,
    systems_file=None,
    n_live=None,
    on_stage_complete=None,
    on_dataset_loaded=None,
):
    import numpy as np
    import os

    # Load the historical PyAutoLens configuration before any model defaults
    # are instantiated; this must work when launched outside the project root.
    from autoconf import conf

    conf.instance.push(new_path=resource_path("config"))

    import autofit as af
    import autolens as al

    from disk_lens import slam
    from disk_lens.systems import load_system_inputs

    if systems_file is None:
        systems_file = resource_path("catalog", "production_systems.json")
    system = load_system_inputs(systems_file, int(dataset_name))
    dataset_name = str(dataset_name)
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    stage_n_live = _live_points(n_live)

    """
    __Dataset__ 

    """

    dataset_path = Path(input_dir).resolve()

    dataset = al.Imaging.from_fits(
        data_path=dataset_path / "data.fits",
        psf_path=dataset_path / "psf.fits",
        noise_map_path=dataset_path / "noise_map.fits",
        pixel_scales=system.pixel_scale_arcsec,
    )
    if on_dataset_loaded is not None:
        on_dataset_loaded(dataset.shape_native)

    dataset_centre = (0.0, 0.0)

    if mask_radius is None:
        mask_radius = system.mask_radius_arcsec
    if not np.isfinite(mask_radius) or mask_radius <= 0:
        raise ValueError("mask_radius must be positive and finite")

    mask = al.Mask2D.circular(
        shape_native=dataset.shape_native,
        pixel_scales=dataset.pixel_scales,
        radius=mask_radius,
    )

    dataset = dataset.apply_mask(mask=mask)

    over_sample_size = al.util.over_sample.over_sample_size_via_radial_bins_from(
        grid=dataset.grid,
        sub_size_list=[4, 2, 1],
        radial_list=[0.1, 0.3],
        centre_list=[dataset_centre],
    )

    dataset = dataset.apply_over_sampling(over_sample_size_lp=over_sample_size)

    """
    __Settings AutoFit__
    
    The settings of autofit, which controls the output paths, parallelization, database use, etc.
    """

    settings_search = af.SettingsSearch(
        path_prefix=str(output_dir / "searches"),
        unique_tag=dataset_name,
        info=None,
        number_of_cores=number_of_cores,
        session=None,
    )

    positions = al.Grid2DIrregular(values=list(system.positions_yx_arcsec))

    """
    __Redshifts__
    
    The redshifts of the lens and source galaxies.
    
    These are placeholders for now given we probably don't know the redshifts of the lens and source galaxies,
    but amnually input via command line can be added.
    """
    redshift_lens = 0.15
    redshift_source = 0.65

    """
    __HPC Mode__
    
    When running in parallel via Python `multiprocessing`, display issues with the `matplotlib` backend can arise
    and cause the code to crash.
    
    HPC mode sets the backend to mitigate this issue and is set to run throughout the entire pipeline below.
    
    The `iterations_per_update` below specifies the number of iterations performed by the non-linear search between
    output, where visuals of the maximum log likelihood model, lens model parameter estimates and other information
    are output to hard-disk.
    
    There are a number of environment variables which must be set to ensure parallelization is efficient, which
    are set below in this script to ensure the pipeline always runs efficiently even if you have not manually set them.
    """
    from autoconf import conf

    conf.instance["general"]["hpc"]["hpc_mode"] = True
    conf.instance["general"]["hpc"]["iterations_per_update"] = iterations_per_update

    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"

    """
    __SOURCE LP PIPELINE__
    
    The SOURCE LP PIPELINE uses one search to initialize a robust model for the source galaxy's light, which in 
    this example:
    
     - Uses a multi Gaussian expansion with 2 sets of 15 Gaussians each for the bulge and disk (4 sets total)
       for the lens galaxy's light.
    
     - Uses an `Isothermal` model for the lens's total mass distribution with an `ExternalShear`.
    
     - Uses a multi Gaussian expansion with 1 set of 10 Gaussians for the source galaxy's light.
    
     __Settings__:
    
     - Mass Centre: Fix the mass profile centre to (0.0, 0.0) (this assumption will be relaxed in the MASS TOTAL PIPELINE).
    """
    analysis = al.AnalysisImaging(
        dataset=dataset,
        positions_likelihood=al.PositionsLHPenalty(threshold=0.1, positions=positions),
    )

    # Lens Light

    total_gaussians = 15
    gaussian_per_basis = 2

    log10_sigma_list = np.linspace(-3, np.log10(mask_radius), total_gaussians)

    centre_0 = af.UniformPrior(
        lower_limit=dataset_centre[0] - 0.05, upper_limit=dataset_centre[0] + 0.05
    )
    centre_1 = af.UniformPrior(
        lower_limit=dataset_centre[1] - 0.05, upper_limit=dataset_centre[1] + 0.05
    )

    bulge_gaussian_list = []

    for j in range(gaussian_per_basis):
        ell_comps_0 = af.UniformPrior(lower_limit=-0.7, upper_limit=0.7)
        ell_comps_1 = af.UniformPrior(lower_limit=-0.7, upper_limit=0.7)

        gaussian_list = af.Collection(
            af.Model(al.lp_linear.Gaussian) for _ in range(total_gaussians)
        )

        for i, gaussian in enumerate(gaussian_list):
            gaussian.centre.centre_0 = centre_0
            gaussian.centre.centre_1 = centre_1
            gaussian.ell_comps.ell_comps_0 = ell_comps_0
            gaussian.ell_comps.ell_comps_1 = ell_comps_1
            gaussian.sigma = 10 ** log10_sigma_list[i]

        bulge_gaussian_list += gaussian_list

    lens_bulge = af.Model(
        al.lp_basis.Basis,
        profile_list=bulge_gaussian_list,
    )

    log10_sigma_list = np.linspace(-3, np.log10(mask_radius), total_gaussians)

    disk_gaussian_list = []

    for j in range(gaussian_per_basis):
        ell_comps_0 = af.UniformPrior(lower_limit=-0.7, upper_limit=0.7)
        ell_comps_1 = af.UniformPrior(lower_limit=-0.7, upper_limit=0.7)

        gaussian_list = af.Collection(
            af.Model(al.lp_linear.Gaussian) for _ in range(total_gaussians)
        )

        for i, gaussian in enumerate(gaussian_list):
            gaussian.centre.centre_0 = centre_0
            gaussian.centre.centre_1 = centre_1
            gaussian.ell_comps.ell_comps_0 = ell_comps_0
            gaussian.ell_comps.ell_comps_1 = ell_comps_1
            gaussian.sigma = 10 ** log10_sigma_list[i]

        disk_gaussian_list += gaussian_list

    lens_disk = af.Model(
        al.lp_basis.Basis,
        profile_list=disk_gaussian_list,
    )

    # Source Light

    centre_0 = af.GaussianPrior(mean=0.0, sigma=0.3)
    centre_1 = af.GaussianPrior(mean=0.0, sigma=0.3)

    total_gaussians = 10
    gaussian_per_basis = 1

    log10_sigma_list = np.linspace(-3, np.log10(1.0), total_gaussians)

    bulge_gaussian_list = []

    for j in range(gaussian_per_basis):
        gaussian_list = af.Collection(
            af.Model(al.lp_linear.Gaussian) for _ in range(total_gaussians)
        )

        for i, gaussian in enumerate(gaussian_list):
            gaussian.centre.centre_0 = centre_0
            gaussian.centre.centre_1 = centre_1
            gaussian.ell_comps = gaussian_list[0].ell_comps
            gaussian.sigma = 10 ** log10_sigma_list[i]

        bulge_gaussian_list += gaussian_list

    source_bulge = af.Model(
        al.lp_basis.Basis,
        profile_list=bulge_gaussian_list,
    )

    source_lp_result = slam.source_lp.run(
        settings_search=settings_search,
        analysis=analysis,
        lens_bulge=lens_bulge,
        lens_disk=lens_disk,
        mass=af.Model(al.mp.Isothermal),
        shear=af.Model(al.mp.ExternalShear),
        source_bulge=source_bulge,
        mass_centre=dataset_centre,
        redshift_lens=redshift_lens,
        redshift_source=redshift_source,
        n_live=stage_n_live["source_lp_n_live"],
    )
    if on_stage_complete is not None:
        on_stage_complete("source_lp")

    """
    __SOURCE PIX PIPELINE__
    
    The SOURCE PIX PIPELINE uses two searches to initialize a robust model for the `Pixelization` that
    reconstructs the source galaxy's light. 
    
    This pixelization adapts its source pixels to the morphology of the source, placing more pixels in its 
    brightest regions. To do this, an "adapt image" is required, which is the lens light subtracted image meaning
    only the lensed source emission is present.
    
    The SOURCE LP Pipeline result is not good enough quality to set up this adapt image (e.g. the source
    may be more complex than a simple light profile). The first step of the SOURCE PIX PIPELINE therefore fits a new
    model using a pixelization to create this adapt image.
    
    The first search, which is an initialization search, fits an `Overlay` image-mesh, `Delaunay` mesh 
    and `AdaptiveBrightnessSplit` regularization.
    
    __Adapt Images / Image Mesh Settings__
    
    If you are unclear what the `adapt_images` and `SettingsInversion` inputs are doing below, refer to the 
    `autolens_workspace/*/imaging/advanced/chaining/pix_adapt/start_here.py` example script.
    
    __Settings__:
    
     - Positions: We update the positions and positions threshold using the previous model-fitting result (as described 
     in `chaining/examples/parametric_to_pixelization.py`) to remove unphysical solutions from the `Inversion` model-fitting.
    """
    analysis = al.AnalysisImaging(
        dataset=dataset,
        adapt_image_maker=al.AdaptImageMaker(result=source_lp_result),
        positions_likelihood=source_lp_result.positions_likelihood_from(
            factor=3.0, minimum_threshold=0.1, positions=positions
        ),
    )

    source_pix_result_1 = slam.source_pix.run_1(
        settings_search=settings_search,
        analysis=analysis,
        source_lp_result=source_lp_result,
        mesh_init=al.mesh.Delaunay,
        image_mesh_init_shape=(30, 30),
        n_live=stage_n_live["source_pix_1_n_live"],
    )
    if on_stage_complete is not None:
        on_stage_complete("source_pix_1")

    """
    __SOURCE PIX PIPELINE 2 (with lens light)__
    
    The second search, which uses the mesh and regularization used throughout the remainder of the SLaM pipelines,
    fits the following model:
    
    - Uses a `Hilbert` image-mesh. 
    
    - Uses a `Delaunay` mesh.
    
     - Uses an `AdaptiveBrightnessSplit` regularization.
    
     - Carries the lens redshift, source redshift and `ExternalShear` of the SOURCE LP PIPELINE through to the
     SOURCE PIX PIPELINE.
    
    The `Hilbert` image-mesh and `AdaptiveBrightness` regularization adapt the source pixels and regularization weights
    to the source's morphology.
    
    Below, we therefore set up the adapt image using this result.
    """
    from disk_lens.fit_visualization import DiagnosticSafeAnalysisImaging

    analysis = DiagnosticSafeAnalysisImaging(
        dataset=dataset,
        adapt_image_maker=al.AdaptImageMaker(result=source_pix_result_1),
        settings_inversion=al.SettingsInversion(
            image_mesh_min_mesh_pixels_per_pixel=3,
            image_mesh_min_mesh_number=5,
            image_mesh_adapt_background_percent_threshold=0.1,
            image_mesh_adapt_background_percent_check=0.8,
        ),
    )

    source_pix_result_2 = slam.source_pix.run_2(
        settings_search=settings_search,
        analysis=analysis,
        source_lp_result=source_lp_result,
        source_pix_result_1=source_pix_result_1,
        image_mesh=al.image_mesh.Hilbert,
        mesh=al.mesh.Delaunay,
        regularization=al.reg.AdaptiveBrightnessSplit,
        image_mesh_pixels_fixed=1000,
        n_live=stage_n_live["source_pix_2_n_live"],
    )
    if on_stage_complete is not None:
        on_stage_complete("source_pix_2")

    """
    __LIGHT LP PIPELINE__
    
    The LIGHT LP PIPELINE uses one search to fit a complex lens light model to a high level of accuracy, using the
    lens mass model and source light model fixed to the maximum log likelihood result of the SOURCE PIX PIPELINE.
    In this example it:
    
     - Uses a multi Gaussian expansion with 2 sets of 15 Gaussians each for the bulge and disk (4 sets total)
       for the lens galaxy's light.
    
     - Uses an `Isothermal` mass model with `ExternalShear` for the lens's total mass distribution [fixed from SOURCE PIX PIPELINE].
    
     - Uses a `Pixelization` for the source's light [fixed from SOURCE PIX PIPELINE].
    
     - Carries the lens redshift and source redshift of the SOURCE PIPELINE through to the MASS PIPELINE [fixed values].   
    """
    analysis = al.AnalysisImaging(
        dataset=dataset,
        adapt_image_maker=al.AdaptImageMaker(result=source_pix_result_1),
    )

    centre_0 = af.GaussianPrior(mean=dataset_centre[0], sigma=0.1)
    centre_1 = af.GaussianPrior(mean=dataset_centre[1], sigma=0.1)

    total_gaussians = 15
    gaussian_per_basis = 2

    log10_sigma_list = np.linspace(-2, np.log10(mask_radius), total_gaussians)

    bulge_gaussian_list = []

    for j in range(gaussian_per_basis):
        ell_comps_0 = af.UniformPrior(lower_limit=-0.7, upper_limit=0.7)
        ell_comps_1 = af.UniformPrior(lower_limit=-0.7, upper_limit=0.7)

        gaussian_list = af.Collection(
            af.Model(al.lp_linear.Gaussian) for _ in range(total_gaussians)
        )

        for i, gaussian in enumerate(gaussian_list):
            gaussian.centre.centre_0 = centre_0
            gaussian.centre.centre_1 = centre_1
            gaussian.ell_comps.ell_comps_0 = ell_comps_0
            gaussian.ell_comps.ell_comps_1 = ell_comps_1
            gaussian.sigma = 10 ** log10_sigma_list[i]

        bulge_gaussian_list += gaussian_list

    lens_bulge = af.Model(
        al.lp_basis.Basis,
        profile_list=bulge_gaussian_list,
    )

    log10_sigma_list = np.linspace(-2, np.log10(mask_radius), total_gaussians)

    disk_gaussian_list = []

    for j in range(gaussian_per_basis):
        ell_comps_0 = af.UniformPrior(lower_limit=-0.7, upper_limit=0.7)
        ell_comps_1 = af.UniformPrior(lower_limit=-0.7, upper_limit=0.7)

        gaussian_list = af.Collection(
            af.Model(al.lp_linear.Gaussian) for _ in range(total_gaussians)
        )

        for i, gaussian in enumerate(gaussian_list):
            gaussian.centre.centre_0 = centre_0
            gaussian.centre.centre_1 = centre_1
            gaussian.ell_comps.ell_comps_0 = ell_comps_0
            gaussian.ell_comps.ell_comps_1 = ell_comps_1
            gaussian.sigma = 10 ** log10_sigma_list[i]

        disk_gaussian_list += gaussian_list

    lens_disk = af.Model(
        al.lp_basis.Basis,
        profile_list=disk_gaussian_list,
    )

    light_result = slam.light_lp.run(
        settings_search=settings_search,
        analysis=analysis,
        source_result_for_lens=source_pix_result_1,
        source_result_for_source=source_pix_result_2,
        lens_bulge=lens_bulge,
        lens_disk=lens_disk,
        n_live=stage_n_live["light_n_live"],
    )
    if on_stage_complete is not None:
        on_stage_complete("light")

    """
    __MASS LIGHT DARK PIPELINE (with lens light)__

    The MASS LIGHT DARK PIPELINE (with lens light) uses one search to fit a complex lens mass model to a high level of 
    accuracy, using the source model of the SOURCE PIX PIPELINE and the lens light model of the LIGHT LP PIPELINE to 
    initialize the model priors. In this example it:

     - Uses the MGE light model (via lp_chain_tracer) from the LIGHT LP PIPELINE as the stellar mass distribution,
       scaled by a free mass-to-light ratio.

     - The lens galaxy's dark matter mass distribution is a `NFWMCRLudlowSph` whose centre is aligned with the bulge of 
       the light model [1 free parameter: M200].

     - Uses an `Inversion` for the source's light [priors fixed from SOURCE PIX PIPELINE].

     - Carries the lens redshift, source redshift and `ExternalShear` of the SOURCE LP PIPELINE through to the MASS 
       LIGHT DARK PIPELINE.
    """

    analysis = al.AnalysisImaging(
        dataset=dataset,
        adapt_image_maker=al.AdaptImageMaker(result=source_pix_result_1),
        positions_likelihood=source_pix_result_2.positions_likelihood_from(
            factor=3.0, minimum_threshold=0.1, positions=positions
        ),
    )

    lp_chain_tracer = al.util.chaining.lp_chain_tracer_from(
        light_result=light_result, settings_search=settings_search
    )

    dark = af.Model(al.mp.NFWMCRLudlowSph)
    dark.redshift_object = 0.15
    dark.redshift_source = 0.65

    mass_result = slam.mass_light_dark.run(
        settings_search=settings_search,
        analysis=analysis,
        lp_chain_tracer=lp_chain_tracer,
        source_result_for_lens=source_pix_result_1,
        source_result_for_source=source_pix_result_2,
        light_result=light_result,
        dark=dark,
        n_live=stage_n_live["mass_n_live"],
    )
    if on_stage_complete is not None:
        on_stage_complete("mass_light_dark")

    analysis = al.AnalysisImaging(
        dataset=dataset,
        adapt_image_maker=al.AdaptImageMaker(result=source_pix_result_1),
        positions_likelihood=source_pix_result_2.positions_likelihood_from(
            factor=3.0, minimum_threshold=0.1, positions=positions
        ),
    )

    return mass_result


def _fit_result_checksums(output, result_directory):
    """Record the final sampler's completion marker and scientific outputs."""
    result_directory = Path(result_directory).resolve()
    result_directory.relative_to(output)
    files = {}
    for name in (".completed", "model.results", "files/samples.csv"):
        path = result_directory / name
        try:
            files[str(path.relative_to(output))] = result_checksum(path)
        except FileNotFoundError as exc:
            raise ValueError(f"completed fit is missing result file: {path}") from exc
    return files


def _validate_completed_result(output, manifest):
    relative = manifest.get("result_directory")
    if not isinstance(relative, str) or not manifest.get("result_files_sha256"):
        raise ValueError(
            "completed fit has no output checksums; keep the existing results "
            "and use a new output directory for a verified run"
        )
    actual = _fit_result_checksums(output, output / relative)
    if actual != manifest["result_files_sha256"]:
        raise ValueError("completed fit result checksum mismatch")


def run_main_slam(
    config_file,
    systems_file,
    *,
    subhalo_id,
    dataset_dir,
    output_dir,
    number_of_cores=1,
    iterations_per_update=5000,
    n_live_override=None,
    mask_radius_override=None,
    reuse_completed=False,
):
    """Record inputs and call the original-style fit used by both CLIs."""
    config_path, config = _load_config(config_file)
    systems_file = Path(systems_file).resolve()
    system = load_system_inputs(systems_file, subhalo_id)
    if not system.positions_yx_arcsec:
        raise ValueError(f"subhalo {subhalo_id} has no image-plane positions")
    if number_of_cores < 1 or iterations_per_update < 1:
        raise ValueError("number_of_cores and iterations_per_update must be positive")
    if mask_radius_override is not None:
        import math

        if not math.isfinite(mask_radius_override) or mask_radius_override <= 0:
            raise ValueError("mask_radius_override must be positive and finite")
    effective_mask_radius = (
        system.mask_radius_arcsec
        if mask_radius_override is None
        else float(mask_radius_override)
    )
    live = _live_points(
        n_live_override if n_live_override is not None else config["searches"]
    )
    dataset_dir = Path(dataset_dir).resolve()
    files = {
        name: dataset_dir / f"{name}.fits" for name in ("data", "psf", "noise_map")
    }
    missing = [str(path) for path in files.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing imaging files: {missing}")
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    cache_dir = output / ".cache"
    cache_dir.mkdir(exist_ok=True)
    os.environ.setdefault("NUMBA_CACHE_DIR", str(cache_dir / "numba"))
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir / "matplotlib"))
    for name in (
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OMP_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[name] = "1"
    if version("autolens") != SUPPORTED_AUTOLENS_VERSION:
        raise RuntimeError(
            f"validated PyAutoLens version is {SUPPORTED_AUTOLENS_VERSION}; "
            f"found {version('autolens')}"
        )
    # Keep the complete original AutoConf settings with each resumed/new run.
    # Older runs did not record them; the 2026-09-27 restoration is documented.
    autolens_config_dir = resource_path("config")
    config_files = {
        str(path.relative_to(autolens_config_dir)): checksum(path)
        for path in sorted(autolens_config_dir.rglob("*"))
        if path.is_file()
    }
    if not config_files:
        raise FileNotFoundError(
            f"missing historical AutoConf configuration: {autolens_config_dir}"
        )
    config_snapshot_path = output / "autolens_config_snapshot.json"
    if config_snapshot_path.is_file():
        snapshot = json.loads(config_snapshot_path.read_text())
        if snapshot["files_sha256"] != config_files:
            raise ValueError("cannot resume fit after AutoConf configuration changed")
    else:
        _write_manifest(
            config_snapshot_path,
            {
                "schema_version": 1,
                "files_sha256": config_files,
                "restoration_note": "Original configuration restored on 2026-09-27; see docs/FITTING.md",
            },
        )
    manifest_path = output / "pipeline_manifest.json"
    current_inputs = {
        "schema_version": 1,
        "status": "running",
        "completed_stages": [],
        "subhalo_id": int(subhalo_id),
        "experiment": "main",
        "input_system_status": system.status,
        "config_sha256": checksum(config_path),
        "systems_sha256": checksum(systems_file),
        "dataset_sha256": {name: checksum(path) for name, path in files.items()},
        "dataset_paths": {name: str(path.resolve()) for name, path in files.items()},
        "live_points": live,
        "autolens_version": version("autolens"),
        "autofit_version": version("autofit"),
    }
    current_inputs["mask_radius_arcsec"] = effective_mask_radius
    if manifest_path.is_file():
        previous = json.loads(manifest_path.read_text())
        # Older records stored only explicit overrides. The system hash below
        # must still match before its default radius can be trusted.
        previous.setdefault(
            "mask_radius_arcsec",
            previous.get("mask_radius_override_arcsec", system.mask_radius_arcsec),
        )
        comparable = set(current_inputs) - {"status", "completed_stages"}
        changed = [
            key for key in comparable if previous.get(key) != current_inputs[key]
        ]
        if changed:
            raise ValueError(
                f"cannot resume fit after inputs changed: {sorted(changed)}"
            )
        if previous.get("status") == "complete":
            if not reuse_completed:
                raise FileExistsError(f"fit is already complete: {output}")
            required_stages = {
                "source_lp",
                "source_pix_1",
                "source_pix_2",
                "light",
                "mass_light_dark",
            }
            if not required_stages.issubset(previous.get("completed_stages", [])):
                raise ValueError("completed fit is missing stage completion records")
            _validate_completed_result(output, previous)
            return output
        manifest = current_inputs | {
            "completed_stages": previous.get("completed_stages", []),
            "status": "running",
        }
    else:
        manifest = current_inputs
    manifest.pop("error", None)
    _write_manifest(manifest_path, manifest)

    def stage_done(stage):
        if stage not in manifest["completed_stages"]:
            manifest["completed_stages"].append(stage)
        _write_manifest(manifest_path, manifest)

    def dataset_loaded(shape):
        manifest["dataset_shape"] = list(shape)
        _write_manifest(manifest_path, manifest)

    try:
        result = fit(
            dataset_name=str(subhalo_id),
            mask_radius=effective_mask_radius,
            number_of_cores=number_of_cores,
            iterations_per_update=iterations_per_update,
            input_dir=dataset_dir,
            output_dir=output,
            systems_file=systems_file,
            n_live=live,
            on_stage_complete=stage_done,
            on_dataset_loaded=dataset_loaded,
        )
        result_directory = Path(result.paths.output_path).resolve()
        manifest["result_files_sha256"] = _fit_result_checksums(
            output, result_directory
        )
        manifest["result_directory"] = str(result_directory.relative_to(output))
        manifest["status"] = "complete"
        _write_manifest(manifest_path, manifest)
        return result
    except Exception as exc:
        manifest["status"] = "failed"
        manifest["error"] = f"{type(exc).__name__}: {exc}"
        _write_manifest(manifest_path, manifest)
        raise


def main(argv=None):
    import argparse

    parser = argparse.ArgumentParser(
        description="Fit one mock with the original five-stage SLaM pipeline"
    )
    parser.add_argument("--dataset", required=True, help="TNG subhalo ID")
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--systems", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--mask-radius", "--mask_radius", type=float)
    parser.add_argument("--number-of-cores", "--number_of_cores", type=int, default=1)
    parser.add_argument(
        "--iterations-per-update", "--iterations_per_update", type=int, default=5000
    )
    parser.add_argument("--n-live", "--n_live", type=int)
    args = parser.parse_args(argv)
    return run_main_slam(
        args.config or resource_path("configs", "fitting.template.json"),
        args.systems or resource_path("catalog", "production_systems.json"),
        subhalo_id=int(args.dataset),
        dataset_dir=args.input_dir,
        output_dir=args.output,
        number_of_cores=args.number_of_cores,
        iterations_per_update=args.iterations_per_update,
        n_live_override=args.n_live,
        mask_radius_override=args.mask_radius,
    )


if __name__ == "__main__":
    main()
