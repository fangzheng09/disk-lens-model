# Particle-mass selection Einstein radii

This calculation measures the selection radius of the TNG particle mass
distribution. It is separate from the mock's SPH stellar light-template
coefficient plus analytic Ludlow NFW halo, and from fitted lens parameters.
Do not substitute a particle-mass selection radius for a `_legacy` mock
radius or a central-deflection/fitting mask input.

## Mass reconstruction

[GLAMER Part I, section 3.1](https://arxiv.org/pdf/1312.1128) motivates
three-dimensional adaptive smoothing and independent particle species.
The production mass renderer shares the foreground's line-of-sight-integrated
spherical M4 kernel. GLAMER's original common truncated Gaussian and its
newer 2D B-spline have different shapes and size conventions; this code
does not claim to reproduce a GLAMER tree-code calculation.

| Component | Mass and outer kernel support H |
|---|---|
| Stars and wind, PartType4 | `Masses`; stored `StellarHsml`, approximately 32 same-type neighbours |
| Dark matter, PartType1 | exact `Header/MassTable[1]`; stored `SubfindHsml`, approximately 64 DM neighbours |
| Gas, PartType0 | `Masses`; `H = 2 (3 V / 4 pi)^(1/3)`, with `V = Masses / Density` |

Lengths become physical kpc using `a/h`; masses become Msun using `1e10/h`.
The gas support follows a documented TNG cell-volume rendering convention,
an approximation to AREPO cells. Gas `SubfindHsml` follows DM neighbours
and is not a gas-cell length. Stored stellar/DM supports come from the full
snapshot, before any image cut. See the
[TNG specifications](https://www.tng-project.org/data/docs/specifications/).

Wind particles contribute to mass even though stellar foreground rendering
excludes them. The main result retains the historical stars+wind, DM and gas
scope. A separately labelled `with_bh_radius_arcsec` also includes compact BH
**dynamical particle Masses**. `BH_Mass` is recorded for interpretation and is
not added again; the accretion length `BH_Hsml` is not a lensing support.
The API numerical flags apply to the primary three-component radius. A BH
diagnostic very close to the cut needs its own grid/FOV assessment before
it can define a different, four-component selected sample.

The particle domain contains the bound subhalo, including its full LOS depth.
It does not add FoF fuzz or unrelated foreground/background structures.
All components share the production iterative stellar centre, principal-axis
rotation and periodic offsets. This rendering centre is independent of the
adopted `SubhaloPos/SubhaloVel` kinematic screening frame.

## Lensing and screening definition

Use Planck18 and lens/source redshifts 0.15/0.65. The physical snapshot scale
factor comes from its Header. Deposit mass per pixel and calculate
`kappa = M_pixel / (A_pixel Sigma_crit)`.

The shear is an isolated linear FFT convolution with Green functions integrated
analytically over square source pixels. Self shear is zero by symmetry and
local convergence is the input mass density. No periodic images, compensating
mean-density sheet, central artificial deflections or source/fitting mask enter
this calculation. Finite-field losses are recorded, without renormalizing
clipped particle mass into the field. Unresolved supports have a numerical
sampling floor of `pixel_scale / (sqrt(2) * subpixels)` (plus a negligible
round-off margin), ensuring at least one fine-pixel sample. This is a grid
regularization, not a gravitational softening scale or a replacement for
the component support definition. Per-component floor counts are recorded
in each run, and grid/subpixel refinement checks its effect on the radius.

Trace `lambda_t = 1 - kappa - |gamma| = 0` by marching squares. The main
effective radius is `sqrt(A/pi)` for the **largest closed tangential critical
curve**. Save other curve areas and the historical sum of individual radii as
diagnostics. The latter can exceed the main galaxy radius when several loops
exist. A critical region touching the boundary is invalid; no curve is a valid
zero radius, distinct from a computation error.

Render XY, XZ and YZ principal-axis views and apply the historical strict cut
`max(R_XY, R_XZ, R_YZ) > 0.5 arcsec`. Screening results do not automatically
change the user-approved `production_37` membership or classify fit success.

## Invocation and numerical checks

Install the selection dependencies in addition to an `illustris_python`
reader compatible with the local TNG directory layout:

```bash
python -m pip install '.[selection]'

python -m disk_lens selection-radii \
  --tng-root /data/TNG50-1 --snapshot 84 --subhalo 4 \
  --output /work/selection_radii
```

The same entry point is available through `scripts/recalculate_selection.py`.
Each run loads particles once, evaluates 0.05 and 0.025 arcsec grids with two
integration subpixels, and automatically checks a doubled FOV at the coarser
pixel scale. The default main full width is 10 stellar half-mass radii; this is
different from the 5-half-mass-radius bound used for the shape calculation.
Foreground light rendering retains its separate eight-subpixel setting.

The output directory contains `{subhalo_id}.json` with three views, all grid
radii, component mass accounting, BH corrections, FOV check, numerical flags,
software versions and a method-source hash. Existing results are not overwritten.

The grid change of the actual screening statistic `R_max` and the largest
per-view FOV change must be no larger than
`max(0.005 arcsec, 0.01 * R_max)`. The stricter all-view grid flag is also
retained: a tiny weak-view curve can remain under-resolved while the maximum
and its screening classification are stable. Threshold stability also requires the
distance from 0.5 arcsec to exceed the combined observed grid/FOV change plus
`0.1 * finest_pixel_scale`. A single requested pixel scale cannot establish
either grid-convergence flag and is labelled accordingly. These are numerical checks;
they do not claim a unique mass reconstruction below simulation resolution.

If flags remain false, prepare a new output directory with finer grids or a
larger field, retaining the same component supports:

```bash
python -m disk_lens selection-radii \
  --tng-root /data/TNG50-1 --snapshot 84 --subhalo 4 \
  --pixel-scales 0.025 0.0125 --fov-multiple 20 \
  --output /work/selection_radii_refined
```

`--support-scale` is an explicit sensitivity setting and changes the reconstructed
mass density. Record it separately from numerical grid refinement. The scientific
default is 1.0. Subpixel refinement is controlled by `--subpixels`.

## Recorded 2026-10-06 audit

All 216 mass/kinematically preselected objects were recalculated in three
views with independent grid and doubled-field checks. After refinement,
173 satisfy the strict 0.5 arcsec radius cut;
37 of the 37 confirmed production members pass.
The confirmed membership is unchanged. See the
[results and provenance](../catalog/audits/einstein_sph_20261006/README.md).
Inherited mock parameters and fit classifications are separate records.
