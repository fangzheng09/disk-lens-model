# Third-party provenance and release status

## Project-authored code

The authors have selected the [MIT License](LICENSE) for this project's
original code and documentation, copyright (c) 2026 Zheng Fang and
contributors. The standard MIT text is reproduced unchanged apart from the
copyright holder and year. [NOTICE](NOTICE) records the scope of this grant.

The grant does not relicense the upstream-derived SLaM Python modules under
`src/disk_lens/slam/` or the copied `config/` snapshot. Their applicable
upstream terms remain a separate release requirement, described below.

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

## Workspace permission still to confirm

No `LICENSE`, `LICENSE.txt` or `COPYING` was found in the root of the pinned
workspace tree or in the available local history for those names. The
[current workspace repository](https://github.com/PyAutoLabs/autolens_workspace)
also did not identify a repository license when checked on 26 October 2026.
The [PyAutoLabs organization page](https://github.com/PyAutoLabs) describes its
packages as MIT licensed, but the applicable terms for this separately copied
historical workspace material have not been established by this audit.

The separately installed
[PyAutoLens library has an MIT license](https://github.com/PyAutoLabs/PyAutoLens/blob/main/LICENSE).
That library license is not being substituted for missing workspace-specific
evidence. Obtain an applicable upstream license notice or permission for the
copied SLaM/config material before distributing it publicly, and preserve the
required notice with that material.

The project-authored code's MIT license does not replace the third-party
permission step above. This candidate therefore makes no claim that all its
contents are covered by a single verified MIT grant.

## Separately installed software and simulation inputs

Particle extraction uses the separately installed
[IllustrisTNG `illustris_python` reader](https://github.com/illustristng/illustris_python);
it is not vendored here. The simulation uses lenstronomy, and fitting uses the
PyAutoLens ecosystem, Nautilus and Dynesty. These dependencies retain their
own licenses and scientific citation requests; installing a dependency does
not transfer ownership of its source to this project.

IllustrisTNG particle data and generated fit products are not bundled in the
code distribution. Data access and the intended code-release scope are
described in [DATA.md](docs/DATA.md) and [REPOSITORY.md](docs/REPOSITORY.md).
