# Main five-stage fitting

`pipeline.py` links to `src/disk_lens/pipeline.py`. Its complete `fit()`
contains dataset setup, explicit MGE construction, analysis objects and
five SLaM calls. The stage modules live under `src/disk_lens/slam/`.
The package pipeline is the sole maintained implementation of the main fit.

| Stage | Model / search | Default live points |
|---|---|---:|
| Source LP | lens/source MGE, Isothermal, shear; Nautilus | 800 |
| Source Pix 1 | fixed lens light, Delaunay; Nautilus | 800 |
| Source Pix 2 | Hilbert image mesh, 1000 source pixels; DynestyStatic | 1000 |
| Light | refit foreground MGE with previous mass/source; Nautilus | 1000 |
| Mass Light Dark | MGE stellar mass plus spherical Ludlow NFW; Nautilus | 1500 |

## Inputs and invocation

A dataset directory contains `data.fits`, `noise_map.fits`, `psf.fits`.
A schema-1 systems JSON supplies the mask, pixel scale and image positions;
see `systems.py` and `catalog/production_systems.json` for the
schema. The observation FITS size is independent of the smaller ray grid.
The bundled catalogue provides fiducial inputs. Mask experiments use explicit
per-run overrides and independent outputs. Historical final-stage 2000-live
results retain their original saved settings.

```bash
python pipeline.py --dataset 3 --input-dir /work/mock_3       --systems /work/systems.json --output /work/fit_3 --number-of-cores 30
```

Equivalently, from anywhere after wheel installation:

```bash
disk-lens fit-main --subhalo 3 --dataset /work/mock_3       --systems /work/systems.json --config /work/fitting.json       --output /work/fit_3 --cores 30
```

Use the unified scientific environment with PyAutoLens 2025.1.18.7. The standard
`templates/fitting.template.json` is bundled; custom files must preserve the
fixed model settings, while sampler live points can be explicitly supplied.
Production runs use the full live-point settings in the fitting template.

## Original configuration and resume records

All 94 original AutoConf files are preserved, loaded before creating model
defaults and shipped in the wheel. In particular,
`general.yaml:test.disable_positions_lh_inversion_check: true` is a
historical setting. It is retained together with the original priors and
stage-specific image-position constraints. Pix2 uses the diagnostic guard
from the tested research workspace: an `InversionException` in an intermediate
plot or optional runtime profile is logged and skips that diagnostic only.
The likelihood, mesh thresholds and priors are unchanged; final-plot errors
and unrelated exceptions still propagate.

`pipeline_manifest.json` records FITS/config/system hashes, library
versions, effective mask radius, live points and completed stages. `autolens_config_snapshot.json`
records the configuration hashes. Resume uses native sampler checkpoints
and rejects changed inputs/settings, including removing an override that
changes the effective mask radius. Older records without an effective-radius
field are interpreted using their explicit override or the unchanged systems
catalogue. `run-main` also validates inputs, live points and stage completion
records before reusing a completed fit. New completed runs additionally record
and verify the final sampler's `.completed`, `model.results` and
`files/samples.csv` checksums. Older completed runs without output checksums
remain intact but cannot be automatically reused as verified results; use a
new output directory for a verified run. This does not prevent resuming an
unfinished older run whose inputs still match. Direct `fit-main` still refuses an
already completed output. Different scientific experiments need separate
output directories.

AutoFit may replace each completed search directory with a sibling `.zip`
archive. Completion checks and parameter summaries read the same logical
result paths directly from either storage format; checksums refer to the
uncompressed file contents. Archives need not be extracted. Missing result
members, duplicate members, corrupt contents and changed checksums remain
errors. Automatic compression does not indicate a fitting failure.

## Cluster execution

The repository CLI and `pipeline.py` accept explicit interpreter, dataset,
systems, configuration and output paths. The workflow runner uses the same interpreter for mock and fitting by default;
optional interpreter-path flags must satisfy the saved unified version requirements. Execute
long rendering and sampling work on allocated compute resources. Limit numerical-library threads before importing NumPy/SciPy;
the sampler can still use its allocated worker processes.
