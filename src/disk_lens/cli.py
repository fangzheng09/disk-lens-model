"""Production SPH rendering, particle-mass screening and five-stage lens fitting."""

import argparse
import os

from ._paths import resource_path


def main(argv=None):
    # Set numerical thread limits before any scientific imports. Samplers use cores.
    for name in [
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OMP_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ]:
        os.environ[name] = "1"
    parser = argparse.ArgumentParser(prog="disk-lens", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    extract = sub.add_parser(
        "extract-stars", help="extract one TNG subhalo to stellar NPZ"
    )
    extract.add_argument("--tng-root", required=True)
    extract.add_argument("--snapshot", type=int, default=84)
    extract.add_argument("--subhalo", type=int, required=True)
    extract.add_argument("--output", required=True)
    for command, help_text in [
        ("lightmap", "render the production SPH foreground"),
        ("simulate", "generate the validated SPH-light mock"),
    ]:
        stage = sub.add_parser(command, help=help_text)
        stage.add_argument("--config", required=True)
        stage.add_argument("--output", required=True)
    catalog = sub.add_parser(
        "validate-catalog", help="validate stored sample membership"
    )
    catalog.add_argument(
        "--catalog", default=str(resource_path("catalog", "sample_sets.json"))
    )
    prepare = sub.add_parser(
        "prepare-main", help="prepare a resolved SPH production plan"
    )
    prepare.add_argument("--systems", required=True)
    prepare.add_argument("--subhalo", type=int, required=True)
    prepare.add_argument("--particles", required=True)
    prepare.add_argument(
        "--fitting-template",
        default=str(resource_path("configs", "fitting.template.json")),
    )
    prepare.add_argument("--output", required=True)
    prepare.add_argument("--seed", type=int, required=True)
    prepare.add_argument(
        "--central-deflection-radius",
        type=float,
        help="override the adopted central deflection radius in the system record",
    )
    prepare.add_argument("--stellar-normalization", type=float)
    prepare.add_argument("--mass-at-200", type=float)
    prepare.add_argument("--sph-workers", type=int, default=1)
    prepare.add_argument("--sph-subpixels", type=int, default=8)
    run = sub.add_parser(
        "run-main", help="execute SPH → mock → fit → parameter summary"
    )
    run.add_argument("--plan", required=True)
    run.add_argument(
        "--through", choices=["lightmap", "mock", "fit", "summary"], default="summary"
    )
    run.add_argument(
        "--mock-python",
        help="optional mock interpreter (default: current Python; AutoLens 2025.1.18.7)",
    )
    run.add_argument(
        "--fit-python",
        help="optional fitting interpreter (default: current Python; AutoLens 2025.1.18.7)",
    )
    fit = sub.add_parser(
        "fit-main", help="fit an existing SPH mock with five-stage SLaM"
    )
    fit.add_argument(
        "--config", default=str(resource_path("configs", "fitting.template.json"))
    )
    fit.add_argument("--systems", required=True)
    fit.add_argument("--subhalo", type=int, required=True)
    fit.add_argument("--dataset", required=True)
    fit.add_argument("--output", required=True)
    for stage in [run, fit]:
        stage.add_argument("--cores", type=int, default=1)
        stage.add_argument("--iterations-per-update", type=int, default=5000)
        stage.add_argument(
            "--n-live", type=int, help="override live points in all fitting stages"
        )
    summarize = sub.add_parser(
        "summarize-main", help="summarize the plan's completed final posterior"
    )
    summarize.add_argument("--plan", required=True)
    summarize.add_argument("--parameter", action="append")
    samples = sub.add_parser(
        "summarize-samples", help="summarize an explicit posterior CSV"
    )
    samples.add_argument("--samples", required=True)
    samples.add_argument("--output", required=True)
    samples.add_argument("--parameter", action="append")
    selection = sub.add_parser(
        "selection-radii", help="calculate SPH particle-mass selection radii"
    )
    selection.add_argument("--tng-root", required=True)
    selection.add_argument("--snapshot", type=int, default=84)
    selection.add_argument("--subhalo", type=int, required=True)
    selection.add_argument("--output", required=True)
    selection.add_argument(
        "--pixel-scales", type=float, nargs="+", default=[0.05, 0.025]
    )
    selection.add_argument("--subpixels", type=int, default=2)
    selection.add_argument("--fov-multiple", type=float, default=10.0)
    selection.add_argument("--support-scale", type=float, default=1.0)
    stellar_ml = sub.add_parser(
        "stellar-ml", help="measure real-stellar SPH M/L inside a YZ circular aperture"
    )
    stellar_ml.add_argument("--particles", required=True)
    stellar_ml.add_argument("--radius-arcsec", type=float, required=True)
    stellar_ml.add_argument("--output", required=True)
    stellar_ml.add_argument("--solar-magnitude", type=float, default=4.53)
    stellar_ml.add_argument("--lens-redshift", type=float, default=0.15)
    stellar_ml.add_argument("--source-redshift", type=float, default=0.65)
    stellar_ml.add_argument("--zeropoint", type=float, default=25.948)
    stellar_ml.add_argument("--pixel-scale", type=float, default=0.05)
    stellar_ml.add_argument("--quadrature-order", type=int, default=128)
    truth = sub.add_parser(
        "mock-truth", help="save the unmasked Einstein radius of an existing SPH mock"
    )
    truth.add_argument("--config", required=True)
    truth.add_argument("--mock", required=True)
    truth.add_argument("--output", required=True)
    positions = sub.add_parser(
        "check-positions",
        help="report bright-arc candidates and check stored positions",
    )
    positions.add_argument("--mock", required=True)
    positions.add_argument("--systems", required=True)
    positions.add_argument("--subhalo", type=int, required=True)
    positions.add_argument("--output", required=True)
    positions.add_argument("--min-peak-fraction", type=float, default=0.3)
    positions.add_argument("--min-separation", type=float, default=0.2)
    positions.add_argument("--max-candidates", type=int, default=12)
    positions.add_argument("--max-arc-distance", type=float, default=0.1)
    args = parser.parse_args(argv)
    if args.command == "check-positions":
        from .positions import review_mock_positions

        result = review_mock_positions(
            args.mock,
            args.systems,
            subhalo_id=args.subhalo,
            output=args.output,
            min_peak_fraction=args.min_peak_fraction,
            min_separation_arcsec=args.min_separation,
            max_candidates=args.max_candidates,
            max_arc_distance_arcsec=args.max_arc_distance,
        )
        print(f"Position review: {result}")
        return
    if args.command == "mock-truth":
        from .mock_truth import run_mock_truth

        result = run_mock_truth(args.config, args.mock, args.output)
        print(f"Mock truth: {result}")
        return
    if args.command == "extract-stars":
        from .particles import load_tng_stars
        from pathlib import Path

        destination = Path(args.output)
        if destination.exists():
            raise FileExistsError(f"output already exists: {destination}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        load_tng_stars(
            args.tng_root, snapshot=args.snapshot, subhalo_id=args.subhalo
        ).save(destination)
    elif args.command == "validate-catalog":
        from .catalog import load_sample_catalog

        result = load_sample_catalog(args.catalog)
        sizes = ", ".join(
            f"{name}={len(item.ids)}" for name, item in result.sets.items()
        )
        print(f"Valid: {result.simulation} snapshot {result.snapshot}; {sizes}")
        return
    elif args.command == "prepare-main":
        from .preparation import prepare_main_run

        result = prepare_main_run(
            systems_file=args.systems,
            subhalo_id=args.subhalo,
            particles_file=args.particles,
            fitting_template=args.fitting_template,
            output_dir=args.output,
            seed=args.seed,
            central_deflection_radius_arcsec=args.central_deflection_radius,
            stellar_mass_normalization=args.stellar_normalization,
            mass_at_200_msun=args.mass_at_200,
            sph_workers=args.sph_workers,
            sph_subpixels=args.sph_subpixels,
        )
        print(f"Prepared: {result}")
        return
    elif args.command == "run-main":
        from .preparation import run_prepared_main

        result = run_prepared_main(
            args.plan,
            through=args.through,
            number_of_cores=args.cores,
            iterations_per_update=args.iterations_per_update,
            n_live_override=args.n_live,
            mock_python=args.mock_python,
            fit_python=args.fit_python,
        )
        print(f"Completed through {args.through}: {result}")
        return
    elif args.command == "fit-main":
        from .pipeline import run_main_slam

        run_main_slam(
            args.config,
            args.systems,
            subhalo_id=args.subhalo,
            dataset_dir=args.dataset,
            output_dir=args.output,
            number_of_cores=args.cores,
            iterations_per_update=args.iterations_per_update,
            n_live_override=args.n_live,
        )
    elif args.command == "summarize-main":
        from .preparation import summarize_main_run

        print(f"Saved: {summarize_main_run(args.plan, args.parameter)}")
        return
    elif args.command == "summarize-samples":
        from .analysis import summarize_samples_csv

        summarize_samples_csv(args.samples, args.output, args.parameter)
    elif args.command == "selection-radii":
        from .selection_radii import calculate_selection_radius

        result = calculate_selection_radius(
            args.tng_root,
            snapshot=args.snapshot,
            subhalo_id=args.subhalo,
            output_dir=args.output,
            pixel_scales=args.pixel_scales,
            subpixels=args.subpixels,
            fov_multiplier=args.fov_multiple,
            support_scale=args.support_scale,
        )
        print(f"Selection radius: {result['fine_max_radius_arcsec']:.6f} arcsec")
        print(
            f"All views converged: {result['resolution_converged']}; "
            f"selection maximum converged: {result['selection_resolution_converged']}; "
            f"threshold stable: {result['threshold_stable']}"
        )
    elif args.command == "stellar-ml":
        from pathlib import Path
        import hashlib
        import json
        from astropy.cosmology import Planck18

        from ._validation import scalar
        from .particles import StellarParticles
        from .photometry import stellar_kappa_per_count_rate
        from .stellar_ml import stellar_aperture_mass_to_light

        radius = scalar(args.radius_arcsec, "radius_arcsec", positive=True)
        destination = Path(args.output)
        if destination.exists():
            raise FileExistsError(destination)
        factor = stellar_kappa_per_count_rate(
            mass_to_light_msun_per_lsun=1.0,
            solar_absolute_magnitude=args.solar_magnitude,
            zeropoint=args.zeropoint,
            pixel_scale_arcsec=args.pixel_scale,
            lens_redshift=args.lens_redshift,
            source_redshift=args.source_redshift,
            cosmology=Planck18,
        )
        result = stellar_aperture_mass_to_light(
            StellarParticles.load(args.particles),
            radius * Planck18.kpc_proper_per_arcmin(args.lens_redshift).value / 60,
            solar_absolute_magnitude=args.solar_magnitude,
            quadrature_order=args.quadrature_order,
        )
        result.update(
            {
                "radius_arcsec": radius,
                "stellar_mass_normalization": result["mass_to_light_msun_per_lsun"]
                * factor,
                "unit_ml_template_factor": factor,
                "lens_redshift": args.lens_redshift,
                "source_redshift": args.source_redshift,
                "zeropoint": args.zeropoint,
                "pixel_scale_arcsec": args.pixel_scale,
                "cosmology": "Planck18",
                "particle_input_sha256": hashlib.sha256(
                    Path(args.particles).read_bytes()
                ).hexdigest(),
            }
        )
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("x") as handle:
            json.dump(result, handle, indent=2, allow_nan=False)
            handle.write("\n")
        print(f"Stellar M/L: {result['mass_to_light_msun_per_lsun']:.6f} Msun/Lsun")
    elif args.command == "lightmap":
        from .workflows import run_lightmap

        run_lightmap(args.config, args.output)
    else:
        from .simulation import run_simulation

        run_simulation(args.config, args.output)
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
