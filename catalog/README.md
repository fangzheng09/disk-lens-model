# Manuscript sample

The adopted sample contains **35 TNG50 galaxies**. `sample_sets.json` defines
`paper_parent_35`, the 21-object `paper_main_21`, and the 14 other parent systems.
These last 14 are outside the main component-recovery statistics; membership
does not assign a uniform fitting outcome to them.

`production_systems.json` provides the 35 systems' fiducial mock and fitting
inputs. Values include the native ray grid, source position, image-position
constraints, stellar normalization, halo mass and central source suppression.
Fields ending in `_legacy` retain the runtime input schema and their original
units. The actual adopted settings are retained in each system's five-stage archives.

Validate membership with `disk-lens validate-catalog`. See the
[fitting results](../results/README.md) and [scientific conventions](../docs/CONVENTIONS.md).
