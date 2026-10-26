# Final manuscript results

Final adopted fits and figure samples are awaiting author selection. This
directory defines the publication interface; it currently contains no claimed
final fit outcomes. Large posterior searches and observations remain outside
the repository.

When the manuscript runs are selected, add small CSV/JSON tables and a
provenance manifest containing:

- System ID, analysis inclusion and review status.
- Run ID, experiment kind, complete five-stage rerun versus final-stage
  resampling, and any reused stages.
- Mock/input identity and checksums, mask radius, MGE width-grid limit,
  stage live points and sampler classes, and software/configuration identity.
- Completed-stage and posterior checksums, adopted-run selection rule,
  parameter paths, posterior-weight and quantile definitions.
- Physical aperture and unit definitions for derived masses and structure
  quantities, with joint-posterior propagation and the matching mock truth.
- Separate flags for computational completion, coefficient diagnostics and
  scientific acceptance. A coefficient within 20% of the mock normalization
  alone does not establish physical mass accuracy or convergence.

For independent repeats, retain the run identities and select a whole run's
posterior according to the documented rule. A result based on shared earlier
stages must identify that dependence. Missing or unaccepted systems remain
explicit; no temporary running-stage summary enters the final table.

The dated historical scheduling categories in
[`catalog/fit_status.production.json`](../catalog/fit_status.production.json)
remain separate from these new SPH results.
