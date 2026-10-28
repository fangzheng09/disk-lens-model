# TNG50 disc lens modeling

Code and results for **What can strong lensing recover in edge-on disc
galaxies? A mock study based on TNG50 simulations**, by Zheng Fang,
Xiaoyue Cao, Ran Li and Weipeng Lin.

The adopted sample contains **34 edge-on TNG50 galaxies**. The main
mass-recovery analysis uses **20 adopted decomposed fits**; a Sérsic
true-light diagnostic covers the 34-system sample.

## Results

- [Browse the 20 adopted fits](results/GALLERY.md)
- [Parameter summaries](results/main_fits.csv) · [All fitted parameters](results/posterior_parameters.csv)
- [Enclosed masses](results/mass_recovery.csv) · [Radial structure](results/radial_structure.csv)
- [True-light masses](results/true_light_mass_recovery.csv) · [Summary statistics](results/summary_statistics.json)
- [Sample membership](results/sample_membership.json) · [Figure index](results/figures/README.md)

Each system has its final model report, model definition, weighted parameter
intervals, fitting diagnostics and settings under `results/systems/ID/`.
The [result definitions](results/README.md) specify apertures, units and
posterior summaries. Full sampler archives and large simulation inputs are
kept separately and are available on reasonable request.

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
results/        adopted fits, tables and manuscript figures
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
