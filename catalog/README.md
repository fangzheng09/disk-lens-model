# Manuscript sample

The adopted sample contains **34 TNG50 galaxies**. `sample_sets.json` defines
`paper_parent_34`, the 20-object `paper_main_20`, and the 14 other parent systems.
These last14 are outside the main component-recovery statistics; membership
does not assign a uniform fitting outcome to them.

`production_systems.json` provides the 34 systems' fiducial mock and fitting
inputs. Values include the native ray grid, source position, image-position
constraints, stellar normalization, halo mass and central source suppression.
Fields ending in `_legacy` retain the runtime input schema and their original
units. The actual settings used for each adopted result, including its fitting
mask and live points, are recorded in `results/systems/ID/metadata.json`.

Validate membership with `disk-lens validate-catalog`. See the
[result tables](../results/README.md) and [scientific conventions](../docs/CONVENTIONS.md).
