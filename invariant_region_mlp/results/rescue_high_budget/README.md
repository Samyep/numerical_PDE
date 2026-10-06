# Rescue high-budget result bundle

This bundle contains per-repetition metrics, aggregate summaries, immutable raw
NPZ artifacts, fixed test sets, manifests, provenance, validation, and five
work-normalized PDF figures for the Funding and Rosenbrock HJB rescue study.

Key files:

- `funding_repetitions.csv`: one row per Funding stochastic root.
- `funding_summary.csv`: aggregate and paired Raw-vs-IR Funding metrics.
- `funding_shrink_validation.csv`: separate-seed diagnostic shrinkage tuning.
- `hjb_repetitions.csv`: one row per HJB method/repetition.
- `hjb_summary.csv`: H1/H2 aggregate value, gradient, generator, constraint,
  correction-variance, work, and timing metrics.
- `work_summary.csv`: common compact work/error table.
- `full_summary.json`: machine-readable aggregate results and verdicts.
- `provenance.json`: exact grids, source hashes, commits, and compute totals.
- `validation_audit.json`: bitwise and completeness audit.
- `raw/`: resumable immutable scientific result artifacts.
- `test_points/`: fixed HJB geometries and quadrature references.
- `figures/`: visually verified PDF plots.

Verdicts:

- Funding: **F-B, partially rescued**.
- HJB: **H-B, trending but not reached**.
