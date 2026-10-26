# Adopted mock grids and central source suppression

On 2026-10-08 the user confirmed the historical cropped mock grid and a
0.2 arcsec central deflection radius for all 37 production members.

`grid_comparison_37.csv` records original light sizes, previously inherited
fit/plot grids, adopted native ray grids, and read-only evidence hashes.
The original 33 use the archived native lensed-source FITS shapes. The four
additions retain their tested mock grids. In particular, the archived
468189 arc is 100 pixels; a historical override table says 160 and is not
used for this object.

All adopted grids are even and fit within the original foreground. No
foreground padding or aperture expansion is needed. Stellar deflections
are calculated from the full original SPH foreground before cropping; the
source is then embedded in that foreground before PSF and Poisson noise,
as in the historical mock. The full mock dimensions and native ray grid
remain separate. Fit masks, sample membership and masses are unchanged.
Physical mock-model Einstein radii require deflections rebuilt at central
radius zero; this audit does not recalculate them.
