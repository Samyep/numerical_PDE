# Correction-level baseline at n=3

This directory implements Study F without modifying the frozen
`experiments/double_estimator/` implementation.  `CorrectionBaselineMLP`
subclasses `DoubleEstimatorMLP`; the only estimator change is subtraction of
the node-local plain-generator baseline from the level-1 gradient integrand.

Typical commands from the repository root:

```bash
python -m invariant_region_mlp.experiments.baseline_n3.run priority1 --workers 8
python -m invariant_region_mlp.experiments.baseline_n3.run priority2 --workers 8
python -m invariant_region_mlp.experiments.baseline_n3.run floor256 --workers 8
python -m invariant_region_mlp.experiments.baseline_n3.run floor512 --workers 8
python -m invariant_region_mlp.experiments.baseline_n3.run audit --workers 8
python -m invariant_region_mlp.experiments.baseline_n3.analyze
```
