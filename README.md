# TNG50 disc lens modeling

Code and results for **What can strong lensing recover in edge-on disc
galaxies? A mock study based on TNG50 simulations**, by Zheng Fang,
Xiaoyue Cao, Ran Li and Weipeng Lin.

The adopted sample contains **34 edge-on TNG50 galaxies**. The main
mass-recovery analysis uses **20 adopted decomposed fits**.

## Results

The [20-system results directory](results/README.md) contains each system's
complete five-stage fitting archives: SP, SPix1, SPix2, L and M. Models,
posterior samples, imaging data and diagnostics are preserved inside the
original ZIPs. The adopted catalogue contains 34 galaxies.

## Install

Use one Python 3.11 environment with **PyAutoLens 2025.1.18.7** for mock
generation, physical mock truth and five-stage fitting:

```bash
python -m pip install -c requirements/bootstrap.txt '.[mock,lensing]'
python -m pip install --no-deps -r requirements/production.txt
```

The second command restores the recorded production numerical stack after
a resolver-compatible bootstrap. See [installation](docs/INSTALLATION.md)
for dependency constraints and the external IllustrisTNG reader.

## Workflow

The package extracts stellar particles, renders the foreground, generates
mock observations, runs the five-stage SLaM fit and summarizes its posterior.
The foreground uses the projected three-dimensional M4 kernel with
StellarHsml support, 32 neighbours and 8 subpixels. See [SPH rendering](docs/SPH.md).

```bash
python -m disk_lens extract-stars \
  --tng-root /data/TNG50-1 --snapshot 84 --subhalo 3 \
  --output /work/inputs/stars_3.npz

python -m disk_lens prepare-main \
  --systems catalog/production_systems.json --subhalo 3 \
  --particles /work/inputs/stars_3.npz --seed 2026100803 \
  --central-deflection-radius 0.2 --output /work/run_3

python -m disk_lens run-main --plan /work/run_3/plan.json --cores 30
```

Long rendering and sampling run on allocated compute resources. The fitting
pipeline is `src/disk_lens/pipeline.py`; the root `pipeline.py` is its
compatibility link. [Workflow](docs/WORKFLOW.md), [fitting](docs/FITTING.md),
[data](docs/DATA.md) and [conventions](docs/CONVENTIONS.md) describe the inputs.

```text
results/        complete five-stage results for 20 adopted fits
catalog/        the 34-system sample and fiducial inputs
src/disk_lens/  scientific workflow implementation
scripts/        command entry points
config/         AutoConf settings and priors
templates/      portable workflow configurations
requirements/   unified scientific environment constraints
docs/           methods, units and running instructions
```

Project-authored code is MIT licensed; see [LICENSE](LICENSE),
[NOTICE](NOTICE) and [THIRD_PARTY.md](THIRD_PARTY.md).
Citation metadata is in [CITATION.cff](CITATION.cff).
