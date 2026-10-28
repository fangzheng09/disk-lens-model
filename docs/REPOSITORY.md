# Repository scope

This repository presents the 34-system manuscript sample and20 adopted
main-analysis fits. [Results](../results/README.md) provides the fitted
models, weighted parameter summaries, diagnostics, derived mass and radial
tables, the Sérsic control, and the manuscript figures.

`src/disk_lens/` implements particle extraction, foreground rendering,
mock imaging, five-stage SLaM fitting and parameter summarization. `scripts/`
provides entry points; the root `pipeline.py` links to the implementation.
`config/`, `templates/` and `catalog/` are portable runtime resources.

New runs use PyAutoLens 2025.1.18.7. Their default stage live points are
800/800/1000/1000/1500; each adopted result records its actual settings.
Raw simulation inputs and complete sampler outputs remain in separate run
storage. The published tables and per-system settings describe the adopted results.

See [installation](INSTALLATION.md), [workflow](WORKFLOW.md),
[result definitions](../results/README.md) and [third-party terms](../THIRD_PARTY.md).
