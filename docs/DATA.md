# Data interfaces and availability

TNG50 is available through the IllustrisTNG data portal under its access
and citation terms. This source release contains the main workflow, the 34-system catalogue and the adopted result
tables and images. It does not bundle simulation snapshots, particles, observations or
posterior archives. No public archive URL or DOI for those generated data
has been assigned. The current manuscript directs requests for generated
mock images and modelling outputs to the corresponding author.

Obtain [IllustrisTNG data](https://www.tng-project.org/data/) under the
project's access and citation terms. This workflow uses TNG50-1 snapshot 84;
the snapshot scale factor and observed lens redshift are distinct.
Consult the [field specifications](https://www.tng-project.org/data/docs/specifications/)
and [reader instructions](https://www.tng-project.org/data/docs/scripts/).

## Normalized stellar NPZ

| Array | Unit / meaning |
|---|---|
| `positions_kpc` | `(N,3)`, physical kpc relative to reference point; periodic wrapping handled |
| `masses_msun` | stellar masses in Msun, wind removed |
| `absolute_magnitudes_i` | intrinsic SDSS i absolute AB magnitudes |
| `halfmass_radius_kpc` | physical stellar half-mass radius |
| `metadata_json` | source, selection and unit metadata, read without pickle |
| `smoothing_lengths_kpc` | optional physical 3D support radii |

## Imaging bundle

`data.fits` and `noise_map.fits` use count-rate units per pixel; the PSF is
normalized. Source, foreground and model-component FITS retain their units
and pixel coordinate headers. Manifests record the model route, seed/RNG,
both foreground and ray-grid dimensions, input hashes and actual calibration.

## Adopted results

The repository includes compact results and fitting diagnostics for20 main
systems, the 34-system Sérsic diagnostic and the final figures. See
[results](../results/README.md) for units, apertures, intervals and checksums.
Complete sampler archives, FITS observations and simulation data remain in
separate storage and are available on reasonable request.
