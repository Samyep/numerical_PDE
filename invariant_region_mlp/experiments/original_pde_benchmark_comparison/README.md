# Original-PDE benchmark comparison

This experiment compares Terminal-Aware Tangent-Actor with direct martingale
control variates against unenhanced methods from the official DeepBSDE,
DeepPicardIteration, and SDGD-PINN repositories.  It separates native
single-point tasks from amortized function approximation and freezes Actor
selection before reading final-test references.

The generated test banks, checkpoints, official working directories, and
Slurm logs remain local.  Their hashes and scalar outcomes are recorded in the
tracked result manifests and tables.

Run correctness and source-independent reference preparation with:

```bash
python -m pytest \
  invariant_region_mlp/experiments/original_pde_benchmark_comparison/test_original_benchmarks.py -q
python -m invariant_region_mlp.experiments.original_pde_benchmark_comparison.prepare
```

GPU launchers are in `slurm/`.  `actor_broad.slurm` requires a `BENCHMARK`
environment variable and runs the four Stage-B candidates on four GH200 GPUs.

