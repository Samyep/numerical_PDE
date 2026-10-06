# Allen-Cahn recovery result bundle

This directory contains source-faithful per-repetition results, pathwise
equivalence metrics, aggregate summaries, resumable task artifacts,
provenance, validation, and four work/mechanism figures.

- `repetition_metrics.csv`: every production method/repetition.
- `pathwise_equivalence.csv`: direct Beck versus generic interval IR pairs.
- `summary.csv`: aggregate accuracy, work, timing, and activation metrics.
- `provenance.json`: sources, exact grids, hashes, environment, and runtime.
- `validation_audit.json`: scientific and artifact integrity checks.
- `raw/`: resumable task-level JSON artifacts.
- `figures/`: rendered and visually checked PDF figures.

Verdict: **B**. Mathematical and code containment are exact, but no
source-faithful tested state leaves even the tighter `[0,1]` interval, so no
active finite-budget truncation rescue is observed.
