# Repository scope

This release provides the main modelling workflow: TNG50 stellar-particle extraction,
foreground rendering, mock imaging, five-stage SLaM fitting and weighted
model-parameter summaries. `src/disk_lens/` contains the implementation;
`scripts/` contains thin command entry points. The root `pipeline.py` links
to `src/disk_lens/pipeline.py` for compatibility with direct script
invocations. The implementation is maintained inside the package;
`python -m disk_lens.pipeline` provides its module entry point.

`config/` preserves the AutoConf settings and priors. `templates/` contains
portable workflow settings, and `catalog/` contains the adopted 37-object
sample, fiducial system inputs and dated scientific provenance. These
resources are included in installed distributions.

Independent true-light experiments and Figure 1–8 postprocessing are maintained
in the research archive. They are outside this package's supported
particle-to-parameter-summary interface. Runtime position overlays, fitting
diagnostics, input validation and resume guards are part of this interface.

## Fitting configuration

New main runs default to **800/800/1000/1000/1500** live points, with
Nautilus / Nautilus / DynestyStatic / Nautilus / Nautilus in stage order.
Historical runs also used 2000 live points in the final stage. Each run's
saved manifest determines its actual configuration. Updating release
defaults does not alter prepared plans or frozen production runs.

## Results

The final manuscript sample and adopted runs will be entered under
[`results/`](../results/README.md) after selection and physical postprocessing.
The existing `catalog/fit_status.production.json` is explicitly a dated
historical scheduling record, not a table of new SPH outcomes. Computation
completion, numerical coefficient diagnostics and scientific fit acceptance
are separate fields in the result contract.

Large TNG data, particle files, FITS observations and sampler outputs belong
in external run storage. Compact published tables must identify the complete
adopted run and its configuration; independent posterior samples are never
concatenated across reruns.

## Installation and redistribution

See [installation](INSTALLATION.md) for the unified PyAutoLens 2025.1.18.7
environment and its recorded dependency limitations. Package installation and resource
checks are distinct from a new full scientific sampling run. Third-party
attribution and redistribution terms are recorded in
[`THIRD_PARTY.md`](../THIRD_PARTY.md).
