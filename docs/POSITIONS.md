# Review bright-arc positions

`check-positions` reports peak candidates and assesses stored manual positions
on the actual generated mock. It uses `source_convolved.fits`: the noiseless,
PSF-convolved lensed source. Detecting peaks on the total image would often
select the foreground galaxy's nucleus. Source surface brightness before PSF
is not used to pretend that arbitrary pixels are observed image peaks.

```bash
/path/to/mock/python -m disk_lens check-positions \
  --mock /work/run_83999/mock \
  --systems /work/production_systems.json --subhalo 83999 \
  --output /work/position_review_83999
```

The same command is available through `scripts/check_positions.py`.
No image, catalogue, source centre or fitting job is changed by this command.
Use a fresh output directory for a new review.

## Detection and coordinates

The tool finds 3x3 local maxima, keeps pixels at least 30% of the global
source peak, excludes the image's central source-suppression circle, and
selects up to 12 peaks separated by at least 0.2 arcsec. Strongest peaks are
considered first. Where the local curvature is well-defined, a quadratic
3x3 patch estimates a subpixel peak; otherwise the pixel centre is retained.
Those parameters can be changed using `--min-peak-fraction`,
`--min-separation` and `--max-candidates`.

Reported positions are `(y,x)` in arcseconds with positive y towards upper
image rows, matching `Grid2D` / `ImageGrid`. Odd and even foreground sizes
are handled explicitly. The checker uses `Array2D.from_fits` with the same
production AutoConf settings as fitting. The adopted `flip_for_ds9: true`
flips raw FITS arrays during loading. Display the **loaded native** array
with `origin='upper'`; for these production files, displaying the raw
`fits.getdata` array with `origin='lower'` gives the same orientation.
The `ImageGrid` FITS headers alone do not control AutoArray's row conversion.
Do not confuse raw-storage coordinates with fitting coordinates or globally
flip stored positions based only on a raw-array overlay.

An earlier checker used raw FITS arrays and misclassified y signs. Those
schema-1 reports are superseded by schema-2 reports recording the actual
image loader, `flip_for_ds9`, package version and configuration checksum.

## Reports and limitations

- `report.json`: candidate coordinates, old-point brightness/distance checks,
  coordinate-direction diagnostics and input/code hashes.
- `candidates.csv`: all peak locations and flags indicating whether they lie
  within the fitting mask.
- `catalogue_points.csv`: current manual coordinates, their source brightness
  relative to peak, distance to the bright arc and nearest detected candidate.
- `positions.png`: source-only and total-mock overlays; cyan circles are
  candidates, white/red crosses are stored points passing/failing the check.

A stored point passes the brightness check when it is inside the image and
fitting mask, reaches the chosen source-brightness fraction and lies within
0.1 arcsec of the bright region. `--max-arc-distance` changes that distance.
These are review thresholds, not a declaration that a fit will succeed.
When all stored points are faint but their y-reflections are bright, the
report flags a possible y-sign mismatch. It never applies the reflection.

Multiple local peaks can belong to the same continuous arc. Brightness alone
cannot establish distinct lensed images or a shared source-plane location.
Review and select a suitable subset before adopting candidate coordinates;
existing good manual points need not be moved to every local maximum.
For real observations without a known source component, foreground removal
and noise-aware peak selection require a separate workflow; this command
is specifically for generated mock bundles.
