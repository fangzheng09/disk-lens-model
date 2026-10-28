# Numerical environment

New mock generation, mock truth and five-stage fitting use one Python 3.11
environment with PyAutoLens 2025.1.18.7 and lenstronomy 1.13.2.
See [installation](INSTALLATION.md) for the constrained numerical stack.

Mock ray tracing explicitly sets `over_sample_size=1` to retain the
adopted native-pixel prescription. This is independent of the SPH subpixel
integration and the fitting stages' own oversampling settings. Source,
PSF and realised-Poisson noise prescriptions are preserved.

Prepared plans record their version requirements and validate existing
products before reuse. Saved results retain their original software
provenance; installing the current package does not relabel earlier data.
The published sample contains 35 systems. Each of the 21 adopted fits records
its available software and input checksums in its result metadata.
