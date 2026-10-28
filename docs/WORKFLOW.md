# SPH production workflow

## Inputs

Use TNG50-1 snapshot 84 and an installed `illustris_python` reader. One
schema-1 systems JSON supplies the adopted ray-tracing grid, 0.05 arcsec
pixels, fitting mask, source centre, image positions, stellar normalization
and halo mass. The bundled `catalog/production_systems.json` contains
exactly the adopted 34 systems. The central source-suppression radius is
0.2 arcsec; actual masks and settings of the 20 adopted results are recorded
in their metadata. See [data interfaces](DATA.md) and [sample records](../catalog/README.md).

## Extract and prepare

The production noise seed is `2026100800 + subhalo_id`; the example below
uses system 3. A new seed creates a different noise realization.

```bash
/path/to/python -m disk_lens extract-stars \
  --tng-root /data/TNG50-1 --snapshot 84 --subhalo 3 \
  --output /work/inputs/stars_3.npz

/path/to/python -m disk_lens prepare-main \
  --systems catalog/production_systems.json --subhalo 3 \
  --particles /work/inputs/stars_3.npz \
  --seed 2026100803 --central-deflection-radius 0.2 \
  --output /work/run_3
```

Preparation validates the fixed model, ray grid and scalar settings before
creating output. It records SPH settings, resolved masses, source geometry
and input/configuration hashes in `plan.json`. Stellar normalization and
halo mass come from the supplied system record; `--stellar-normalization`
and `--mass-at-200` are explicit overrides. The normalization is the adopted
light-template coefficient, distinct from a fitted MGE physical M/L. The
recalculated `stellar_mass_normalization` is preferred when present; historical
files fall back to `_legacy`. See [stellar M/L](STELLAR_ML.md).

The foreground size is `floor(10 R_half / angular pixel width)`, preserving
the tested particle centring, principal-axis rotation and nominal 3D cube.
The foreground may be odd-sized; the separate ray-tracing grid must be
square and even-sized, and must fit within the nominal foreground. The
adopted grid is the historical mock crop, not a fitting/plotting grid.
Preparation rejects an oversized crop before creating output. No foreground
padding or aperture expansion is used by the production workflow. Deflections
use the full SPH foreground before being cropped to the ray grid; the source
is then embedded in that foreground before PSF convolution and noise. The
final mock therefore retains the original foreground size. The central
pseudo-core radius defaults to the adopted system record (0.2 arcsec for
all production members), independently of the fitting mask; an explicit
command-line value overrides it.
The default SPH support uses stored StellarHsml; a missing support array is
computed from the 32nd other stellar neighbour before aperture selection.

## Run all stages

```bash
/path/to/python -m disk_lens run-main \
  --plan /work/run_3/plan.json \
  --cores 30
```

The same Python environment executes lightmap, mock, truth, fit and summary.
New plans record PyAutoLens 2025.1.18.7 for both mock and fitting, with
lenstronomy 1.13.2 for stellar deflections; the runner checks these requirements
before expensive work. The optional `--mock-python` and `--fit-python` path
arguments remain available and must satisfy the same saved version requirements.
Stage live points come from the fitting template. Run expensive TNG rendering and sampling on the allocated
compute node; `--cores` controls sampler processes. `--sph-workers` is set
when preparing the plan and controls neighbour queries.

Stages can also be run separately:

```bash
/path/to/python -m disk_lens run-main --plan /work/run_3/plan.json --through mock
/path/to/python -m disk_lens run-main --plan /work/run_3/plan.json --through summary --cores 30
```

A repeated invocation verifies completed products and resumes unfinished
sampler stages. Changed inputs or damaged files cause an error and remain
intact. Schema-1 plans from the earlier generic/CIC workflow are not converted;
prepare a new SPH run directory. Older completed fits lacking output hashes
remain available as archived results but are not automatically reused. Historical schema-2
plans without saved package requirements may reuse their verified existing
mock/truth bundles, preserving the actual recorded versions. They do not
silently regenerate missing historical products under the unified version;
prepare a new plan/output directory for a new complete run.

## Outputs

| Path under the run directory | Contents |
|---|---|
| `config/`, `plan.json` | resolved inputs and scientific settings |
| `lightmap/light.fits` | intrinsic SPH foreground counts/s/pixel |
| `mock/` | data, noise map, PSF, source and foreground components, unit stellar deflections |
| `truth/` | physical mock-model Einstein radius, unmasked stellar deflections and provenance |
| `fit/searches/` | five-stage SLaM posterior and diagnostic outputs |
| `fit/pipeline_manifest.json` | fitting inputs, effective mask, stage status and final-file hashes |
| `summary/parameters.json` | q16, median and q84 for every final model parameter |
| `summary/manifest.json` | plan/posterior hashes and summary columns |

The summary reads exact CSV column names, including chained MGE parameter
paths, and strips AutoFit's whitespace around statistical headers. It uses
posterior weights and excludes zero-weight samples. These are model-parameter
summaries, not posterior-derived enclosed masses, Einstein radii or Figure
1–8 quantities. For a subset of parameters:

```bash
/path/to/python -m disk_lens summarize-main --plan /work/run_3/plan.json \
  --parameter galaxies.lens.dark.mass_at_200
```

Standalone `lightmap`, `simulate` and `fit-main` also accept explicit paths.
The workflow scripts in `scripts/` delegate to the same package entry points.

## Particle-mass selection radii

Use the separate [selection-radius entry point](EINSTEIN_RADIUS.md) to calculate
the TNG particle-mass screening radius. It is independent of the production
mock mass model and its inherited Einstein-radius inputs.

## Review a central suppression radius

The adopted production default is 0.2 arcsec central source suppression.
`--central-deflection-radius` explicitly overrides the per-system value. It replaces stellar deflections inside the circle with 1e6 arcsec
before ray tracing; it is independent of the observational fitting mask.
The PSF is applied afterwards, so light can spread back into the circle.
For the physical mock-model Einstein radius, rebuild stellar deflections
with radius zero and retain the same SPH template, adopted coefficient and
analytic Ludlow halo. Do not derive this truth from the modified deflections
or replace it with the TNG particle-mass screening radius.

The mock stage includes [physical mock truth](MOCK_TRUTH.md). Its Einstein
radius is saved before fitting, with central source suppression disabled
only for this truth calculation. The 0.2 arcsec source-suppression image is
retained as the fitting dataset.

The adopted fitting template uses 1500 live points for the final mass
stage, with 800, 800, 1000 and 1000 for the preceding stages. This is also
the default stage mapping for direct fitting. `--n-live 1500` would override
every stage and is therefore not equivalent. Use a separate output directory
when changing settings; existing prepared or running fits retain their
frozen configuration. Historical main runs used 2000 live points in the final
stage; their saved manifests remain the authority for those runs.

After generating a mock, use [`check-positions`](POSITIONS.md) to inspect
bright-arc candidates and validate stored manual coordinates. The command
saves a report and overlays; it does not modify coordinates or submit fits.
