# Installation

Use one Python 3.11 environment with **PyAutoLens 2025.1.18.7** for
particle preparation, SPH light rendering, mock imaging, physical mock
truth and five-stage lens fitting. The following commands run from a source
checkout. Store simulation inputs and run products separately.

## Unified scientific environment

```bash
python3.11 -m venv /work/envs/disk-lens
/work/envs/disk-lens/bin/python -m pip install -c requirements/bootstrap.txt '.[mock,lensing]'
/work/envs/disk-lens/bin/python -m pip install --no-deps -r requirements/production.txt
/work/envs/disk-lens/bin/disk-lens validate-catalog
```

Both installation steps use this same environment. All five PyAuto packages
use release 2025.1.18.7. Mock stellar deflections retain lenstronomy 1.13.2;
Nautilus 1.0.4, Dynesty 2.1.4 and Colossus 1.3.1 are explicitly constrained.
The primary `production.txt` file records the numerical production versions;
`bootstrap.txt` provides the initial resolver-compatible installation.

The observed fitting environment used NumPy 2.1.3 and SciPy 1.16.0, above
limits declared by some upstream AutoConf/anesthetic and AutoFit/AutoArray
packages. A single normal resolver invocation with `production.txt` cannot
install this combination. The second command deliberately restores these
observed direct versions with `--no-deps`. `pip check` then reports the known
upstream NumPy/SciPy metadata conflicts. This is an explicit reconstruction
of the recorded numerical stack, rather than a claim that upstream metadata
supports that combination. Keeping only the bootstrap gives a different
numerical environment whose scientific equivalence has not been established.

The requirement files are direct constraints, not a full transitive lockfile.
Use a new environment for this installation; running production jobs retain
their original environment. Install the external
[IllustrisTNG reader](https://www.tng-project.org/data/docs/scripts/) following
its own dependency and data-layout instructions. The portable particle input
boundary is a normalized NPZ; see [data interfaces](DATA.md).

## Running the complete workflow

The same interpreter prepares and executes the plan:

```bash
/work/envs/disk-lens/bin/python -m disk_lens prepare-main \
  --systems catalog/production_systems.json --subhalo 3 \
  --particles /work/inputs/stars_3.npz --seed 2026100803 \
  --output /work/run_3
/work/envs/disk-lens/bin/python -m disk_lens run-main \
  --plan /work/run_3/plan.json --cores 30
```

The two optional interpreter-path flags remain compatible with existing
scripts, but new plans require the unified release for both stages. Long
SPH rendering and fitting run on allocated compute nodes. See
[workflow](WORKFLOW.md) for model, resume and output definitions.

## Base package and distributions

The base package installs with `python -m pip install .` for numerical
utilities and parameter summaries. It does not include the full modeling
extras. Build distributions with `python -m build`. For a wheel installation,
replace `'.[mock,lensing]'` in the first installation command with the wheel
path and both extras; use the source checkout's constraint files. From a
working directory outside the checkout, verify:

```bash
disk-lens --help
disk-lens validate-catalog
python -m disk_lens.pipeline --help
```

Installed resources include the catalogue, workflow templates and AutoConf
snapshot. Standalone scripts require the package installed.

## Particle-mass selection radii

Install `.[selection]` for h5py and scikit-image, together with the external
`illustris_python` reader. Selection radii use the base numerical
implementation and do not require a separate PyAutoLens environment. See
[selection-radius instructions](EINSTEIN_RADIUS.md).
