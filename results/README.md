# Adopted manuscript results

This release contains **20 adopted main fits from a 34-system parent sample**,
the 34-system Sérsic true-light diagnostic, and the final manuscript figures.

| File | Contents |
|---|---|
| `sample_membership.json` |34 parent IDs,20 main IDs and figure subsets |
| `main_fits.csv` | Input and recovered stellar normalization, intervals and actual fit masks |
| `posterior_parameters.csv` | Complete weighted summaries for every sampled final parameter |
| `mass_recovery.csv` | Mock and fitted total, stellar and dark-matter aperture masses, in solar masses |
| `radial_structure.csv` | Einstein radius, curvature, convergence and slope summaries |
| `true_light_mass_recovery.csv` |34-system Sérsic control mass summaries |
| `summary_statistics.json` | Cross-system median signed biases and NMAD scatters, in per cent |
| `systems/ID/` | Adopted model, original report, diagnostic images and metadata |
| `figures/` | Manuscript figures in PDF and PNG |

## Parameter summaries

Parameter intervals use the complete final posterior with positive saved
weights, at quantiles 0.15865,0.5 and 0.84135. The numerical stellar template
coefficient `n` is distinct from an enclosed stellar mass and a physical M/L.
Run metadata gives the actual mask, live points, archive identifier and
posterior/data checksums. System 414684 reuses its completed earlier stages
and adopts the separately sampled final mass stage; its posterior is not
combined with another run. [Browse the fitting diagnostics](GALLERY.md).

## Derived masses and radial quantities

Main-analysis derived quantities propagate 200 joint equal-weight posterior
draws per system. Mock masses use the mock Einstein aperture; fitted masses
use each draw's own critical-curve Einstein aperture. Units are solar masses.
The Sérsic true-light diagnostic uses the fixed mock Einstein aperture.
Quoted mass intervals are the 15.865th and 84.135th percentiles.

For radial structure, `R_E` is in arcsec and `alpha2_E` in inverse arcsec;
`xi`, `kappa_E`, `mw_slope` and `annular_slope` are dimensionless. The last
two are the enclosed and narrow-annulus projected mass-weighted slopes.
The annulus is 0.95–1.05 times the effective Einstein radius. Each fit draw
uses its own critical-curve radius, 300 radial samples and 64 azimuths.

For any quantity, the signed residual is `(fit-mock)/mock`. The cross-system
bias is its median; NMAD is 1.48 times the median absolute deviation about
that median. NMAD describes scatter, and is separate from the signed bias.

Large TNG snapshots, particles, FITS observations and complete sampler
archives are available separately on reasonable request. The repository
contains the selected outcomes rather than intermediate runs.
