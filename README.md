# TNG50 disc lens modeling

Research workflow for **What can strong lensing recover in edge-on
disc galaxies? A mock study based on TNG50 simulations**, by Zheng Fang,
Xiaoyue Cao, Ran Li and Weipeng Lin.

The package extracts stellar particles, renders an SPH foreground, generates
a mock observation, runs the five-stage SLaM fit and summarizes its weighted
posterior. The adopted renderer integrates the three-dimensional M4 kernel,
using TNG StellarHsml, 32 neighbours and 8 subpixels. Mock images retain
source intensity 0.4, the sampled Gaussian PSF and the realised-Poisson
noise convention.

## Install

Use one Python 3.11 environment with **PyAutoLens 2025.1.18.7** for both
mock generation and five-stage fitting. From a checkout:

```bash
/path/to/python -m pip install -c requirements/bootstrap.txt '.[mock,lensing]'
/path/to/python -m pip install --no-deps -r requirements/production.txt
```

Both commands operate on the same environment. The second restores the
recorded production NumPy 2.1.3 and SciPy 1.16.0 versions after a
resolver-compatible bootstrap; their known upstream metadata conflicts and
installation details are documented in [installation](docs/INSTALLATION.md).
Lenstronomy 1.13.2 retains the mock deflection implementation. The
[37-system compatibility check](docs/VERSION_COMPATIBILITY.md) verifies the
preserved native-pixel mock prescription under the unified version. The requirement
files are direct constraints rather than a complete transitive lockfile.
Particle extraction also requires the external `illustris_python` reader
and local TNG50 data.

## Run one system

The bundled `catalog/production_systems.json` contains the adopted 37-member
sample: the original 33 plus 4, 83993, 83999 and 297779. It supplies source
and image positions, ray-grid size, fiducial fitting mask, stellar
normalization and halo mass. Stellar aperture M/L, native mock grids and
image-coordinate conventions have dated provenance in `catalog/audits/`.
Run-specific mask experiments retain separate inputs and manifests.

```bash
/path/to/python -m disk_lens extract-stars \
  --tng-root /data/TNG50-1 --snapshot 84 --subhalo 3 \
  --output /work/inputs/stars_3.npz

/path/to/python -m disk_lens prepare-main \
  --systems catalog/production_systems.json --subhalo 3 \
  --particles /work/inputs/stars_3.npz \
  --seed 2026100803 --central-deflection-radius 0.2 \
  --output /work/run_3

/path/to/python -m disk_lens run-main \
  --plan /work/run_3/plan.json \
  --cores 30
```

The production seed is `2026100800 + subhalo_id`. The 0.2 arcsec central
source-suppression radius is distinct from the observational fitting mask.
Physical mock truth uses unsuppressed stellar deflections. New main fits
use **800/800/1000/1000/1500** live points and
**Nautilus / Nautilus / DynestyStatic / Nautilus / Nautilus**. Historical
final-stage 2000-live runs retain their saved configurations.

The complete invocation runs through the parameter summary. Execute long
rendering and sampling on allocated compute resources; `--cores` must match
the allocation. See [workflow](docs/WORKFLOW.md), [fitting](docs/FITTING.md),
[mock truth](docs/MOCK_TRUTH.md) and [coordinate review](docs/POSITIONS.md).

## Layout and scope

```text
src/disk_lens/   rendering, selection radii, mock, fitting and summary
pipeline.py     compatibility link to src/disk_lens/pipeline.py
scripts/        thin command entry points
config/         preserved AutoConf settings and priors
templates/      portable lightmap, mock and fitting configurations
catalog/        sample membership, fiducial inputs and dated provenance
docs/           running instructions and scientific conventions
requirements/   unified production constraints and historical provenance
results/        interface for the final selected compact manuscript tables
```

The main workflow ends at weighted model-parameter summaries. Independent
true-light experiments and Figure 1–8 physical postprocessing are maintained
in the research archive. The separate `selection-radii` command measures
particle-mass selection radii; see [its definitions](docs/EINSTEIN_RADIUS.md).

Final adopted runs and figure samples await author selection under
[`results/`](results/README.md). The old `catalog/fit_status.production.json`
is a dated scheduling record. A completed run or an MGE coefficient close
to the mock normalization does not alone certify physical mass accuracy.
Large particles, FITS observations and posterior searches use external run
storage. See [repository scope](docs/REPOSITORY.md) and [data](docs/DATA.md).

Project-authored code is MIT licensed. Permission for the copied workspace
SLaM/config material remains unconfirmed and must be resolved before public
distribution of this candidate.

Citation metadata is in [CITATION.cff](CITATION.cff); third-party source and
redistribution status are in [THIRD_PARTY.md](THIRD_PARTY.md). This is a
0.1.0 release candidate hosted at
[github.com/fangzheng09/disk-lens-model](https://github.com/fangzheng09/disk-lens-model).
Final manuscript results are pending selection; an archival DOI has not yet
been assigned.
