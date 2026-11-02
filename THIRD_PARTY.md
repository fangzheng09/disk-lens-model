# Third-party provenance

## Project-authored code

The authors have selected the [MIT License](LICENSE) for this project's
original code and documentation, copyright (c) 2026 Zheng Fang and
contributors. The standard MIT text is reproduced unchanged apart from the
copyright holder and year. [NOTICE](NOTICE) records the scope of this grant.

## Copied workspace material

`src/disk_lens/slam/` and `config/` contain material adapted from
[PyAutoLens Workspace](https://github.com/PyAutoLabs/autolens_workspace).
The content was compared with the official workspace tag
[`2025.1.18.7`](https://github.com/PyAutoLabs/autolens_workspace/tree/2025.1.18.7),
which resolves to commit
[`cb5a94f41273c0df26ad9bc269d25014594ebd3d`](https://github.com/PyAutoLabs/autolens_workspace/tree/cb5a94f41273c0df26ad9bc269d25014594ebd3d).
This is a verified comparison baseline, not a recovered record of the original
checkout used to make the historical local copies.

Of the 94 files in this repository's `config/` snapshot, 88 are byte-for-byte
identical to that tag. The following six files contain local changes:

- `config/general.yaml`
- `config/grids.yaml`
- `config/priors/mass/dark/nfw.yaml`
- `config/priors/mass/stellar/gaussian.yaml`
- `config/visualize/mat_wrap_2d.yaml`
- `config/visualize/plots.yaml`

The SLaM stages also contain local changes, including configurable live-point
counts and preserved research-model choices. See
[the stage provenance record](src/disk_lens/slam/ORIGIN.md) for the differences
and file hashes. These are adapted modules, not an unmodified upstream release.

## Separately installed software and simulation inputs

Particle extraction uses the separately installed
[IllustrisTNG `illustris_python` reader](https://github.com/illustristng/illustris_python);
it is not vendored here. The simulation uses lenstronomy, and fitting uses the
PyAutoLens ecosystem, Nautilus and Dynesty. These dependencies retain their
own licenses and scientific citation requests; installing a dependency does
not transfer ownership of its source to this project.

IllustrisTNG simulation snapshots and extracted particle inputs are not bundled.
The Git repository and source distribution include the complete five-stage
fitting archives for the 21 adopted main-analysis systems under `results/`,
together with their SHA256 checksums. These result archives are not included
in the Python wheel. Data access and repository scope are described in
[DATA.md](docs/DATA.md) and [REPOSITORY.md](docs/REPOSITORY.md).
