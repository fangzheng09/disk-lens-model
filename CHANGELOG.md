# Changelog

## 0.1.0 — release candidate

- Package the main lens-modelling workflow from particle extraction through five-stage
  lens fitting and weighted posterior-parameter summaries.
- Include the adopted 37-system catalogue, component-specific SPH provenance,
  mock-truth definitions, coordinate conventions and portable resources.
- Default new main fits to 800/800/1000/1000/1500 live points; retain saved
  configuration records for historical final-stage 2000 runs.
- Bound NumPy below 2.2 so a base installation remains compatible with
  the Astropy 6.1 interfaces used during run preparation.
- Unify new mock generation, truth calculation and fitting on PyAutoLens
  2025.1.18.7 in one Python environment, retaining the observed production
  NumPy/SciPy versions and their documented upstream metadata conflicts.
- Record new plans' package requirements and preserve historical bundles'
  original provenance when reusing existing products.

Final manuscript result tables are awaiting run selection. This candidate
is hosted at [github.com/fangzheng09/disk-lens-model](https://github.com/fangzheng09/disk-lens-model);
an archive DOI has not yet been assigned.
