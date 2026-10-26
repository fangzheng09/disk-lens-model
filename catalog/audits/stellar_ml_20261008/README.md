# Production stellar M/L update, 2026-10-08

All 37 approved systems were measured at their updated YZ selection radius.
`comparison_37.csv` compares saved historical coefficients and implied
physical M/L with the new SPH aperture measurement and consistently calibrated
coefficient. Old implied M/L is `legacy_n / 0.1970563339860656`, retaining
historical rounding; it is not a new measurement. Direct old M/L and aperture
records are available for 15 objects and blank for the other 22.

`measurements_37.json` retains per-object aperture measurements, checks,
particle-extraction provenance and source hashes. `summary.json` records
the comparison ranges and independent numerical checks. `manifest.json`
records checksums of these audit products. Historical values are preserved
in `production_systems.json` and the original 41-system provisional file.

The physical M/L changes range from -0.05623% to +0.01537%. The mock n
changes range from +0.86835% to +0.94061%, dominated by correcting the old
mixed solar-magnitude calibration. The maximum relative M/L difference
between quadrature orders 64 and 128 is below 4.2e-11; the Cartesian image
check agrees to about 4e-6 and independent 3D kernel checks to below 1e-7
in aperture fraction. These are integration checks, not physical errors.

All per-object centres and axis ratios were checked against the completed
mass-radius audit. No sample membership, fit classifications, image data
or posteriors were changed. See [measurement and use](../../../docs/STELLAR_ML.md).
