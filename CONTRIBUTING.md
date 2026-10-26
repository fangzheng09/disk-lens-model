# Maintaining the production workflow

Preserve the adopted SPH geometry/kernel, historical mock prescription,
five-stage model, priors and the documented sampler defaults. Changes to scientific methods
must state their units and numerical impact. Keep per-system values and
catalogue provenance and per-run overrides explicit.

Use Ruff 0.16.9 to format and check `src scripts`. Copied SLaM modules retain
their source formatting and provenance. Build the sdist/wheel and inspect
installed entry points and resources. Development verification material is
maintained outside the production repository.

Keep generated particles, observations, posteriors and private job metadata
in external or ignored run storage. Add small data only with its provenance.
