# Repository scope

The adopted sample contains 35 TNG50 galaxies. The main analysis uses 21
systems, each supplied with its complete five-stage fitting archives under
[`results/`](../results/README.md). The archives preserve model files,
posterior samples, imaging data and fitting diagnostics.

`src/disk_lens/` implements the scientific workflow; `pipeline.py` links to
its implementation. `config/`, `templates/` and `catalog/` provide runtime
settings and inputs. See [installation](INSTALLATION.md) and
[workflow](WORKFLOW.md) for running instructions.
