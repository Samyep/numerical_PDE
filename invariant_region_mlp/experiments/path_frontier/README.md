# Pathwise-only equal-cost frontier

The frozen specification is `PREREGISTRATION.md`. This study imports the
unchanged double-estimator implementation, reuses every applicable frozen
result, and computes only the missing `path` rows.

From the repository root:

```powershell
python -m invariant_region_mlp.experiments.path_frontier.run run --workers 8
python -m invariant_region_mlp.experiments.path_frontier.run audit --workers 8
python -m invariant_region_mlp.experiments.path_frontier.analyze
```

All task artifacts are atomic and the run is resumable.
