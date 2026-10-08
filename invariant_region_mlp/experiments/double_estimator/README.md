# Double-estimator study

The frozen specification is `PREREGISTRATION.md`.  The implementation adds a
double generator estimator without changing `FullHistoryMLP` or any prior
experiment.

Run the registered priority blocks from the repository root:

```powershell
python -m invariant_region_mlp.experiments.double_estimator.run core --workers 6
python -m invariant_region_mlp.experiments.double_estimator.run c2 --workers 6
python -m invariant_region_mlp.experiments.double_estimator.run ec --workers 6
python -m invariant_region_mlp.experiments.double_estimator.run p4 --workers 6
python -m invariant_region_mlp.experiments.double_estimator.analyze
```

Every task writes an atomic JSON artifact, so every stage is resumable.  The
collector builds `results/double_estimator/rows.csv`; the analyzer writes all
pre-registered tables, figures, verdicts, and `docs/DOUBLE_ESTIMATOR_REPORT.md`.

