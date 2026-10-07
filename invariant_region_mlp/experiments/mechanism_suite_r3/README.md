# Mechanism suite round 3

This directory contains the isolated implementation of the frozen round-3
protocol in `ROUND3_PREREGISTRATION.md`.  It subclasses or wraps the earlier
mechanism-suite code and never modifies `FullHistoryMLP`.

Run from the repository root:

```text
python -m invariant_region_mlp.experiments.mechanism_suite_r3.run --stage all --workers 6
python -m invariant_region_mlp.experiments.mechanism_suite_r3.analyze
```

The runner is resumable and refuses to overwrite an artifact whose protocol
metadata does not match.  The registered order is S0, S2, S1, S3, S4.  S1
runs d=20/100 before the d=400 cutoff pilots and d=400 grid.

