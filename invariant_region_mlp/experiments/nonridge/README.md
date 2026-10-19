# Non-ridge spectral HJB benchmark

This isolated experiment extends the frozen CPU and GPU MLP solvers to the
pre-registered DCT-separable spectral terminal family. The GPU DCT/IDCT uses
the Makhoul FFT construction in float64.

Run from the repository root:

```bash
python -m invariant_region_mlp.experiments.nonridge.run validate
python -m invariant_region_mlp.experiments.nonridge.run calibrate
python -m invariant_region_mlp.experiments.nonridge.run stage_a
python -m invariant_region_mlp.experiments.nonridge.analyze
```

Commit and push the complete Stage-A results and interim report before any
Stage-B command. The runner verifies that the clean local HEAD equals its
upstream commit before it permits Stage B:

```bash
python -m invariant_region_mlp.experiments.nonridge.run stage_b
python -m invariant_region_mlp.experiments.nonridge.analyze
```

All stages are resumable. Use `run progress` and `run collect` for audits.
