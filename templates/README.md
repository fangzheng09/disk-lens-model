# Production configurations

| Template | Role |
|---|---|
| `lightmap.template.json` | adopted SPH foreground, 32 neighbours and 8 subpixels |
| `main.template.json` | SPH-light mock with source intensity 0.4 and historical PSF/noise |
| `fitting.template.json` | preserved five-stage model and full sampler settings |

Use `prepare-main` to resolve per-system inputs and write the three configs
plus `plan.json`. For standalone commands, fill every null field. A null
lightmap grid shape requests the integer-pixel 10 R_half aperture. Relative
input paths resolve against each JSON directory; output CLI paths resolve
against the launch directory. Filled private copies may use `local.*.json`.

The separate `config/` tree is the preserved AutoConf settings and priors;
it is packaged along with these workflow templates and the catalogue.

New five-stage runs default to 800/800/1000/1000/1500 live points. Historical
completed runs may have used 2000 in the mass stage; their stored configurations
remain the authority for those results. The sampler sequence is Nautilus,
Nautilus, DynestyStatic, Nautilus, Nautilus.
