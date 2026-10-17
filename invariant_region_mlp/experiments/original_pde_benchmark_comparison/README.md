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

GPU launchers are in `slurm/`.  The cluster permits one active job for this
allocation, so the recorded study uses the `*_serial.slurm` launchers on one
GH200.  Actor jobs take `BENCHMARK=hjb_log_d100`, `dpi_ou_gmm_d100`, or
`allen_cahn_d100` through `--export`:

1. `actor_broad_serial.slurm` runs the four broad candidates; `select_trials.py`
   writes the frozen successive-halving decision.
2. `actor_stage_c_serial.slurm` continues the two reference-free winners and
   runs the architecture, data, iteration, and plain-label ablations.
3. `actor_max_quality_serial.slurm` is a conditional one-update upper-budget
   check at 16,384 points, 256 paths/label, and 64 label steps.
4. `actor_final_seeds_serial.slurm` and
   `actor_final_evaluation_serial.slurm` reproduce the frozen selection with
   three seeds and 32/64/128-step evaluation.

Official baseline launchers retain their source algorithms.  The Deep BSDE
launcher is an audited PyTorch compatibility port; DPI, PINN-HTE, and D-DBSDE
run through the authors' package; SDGD-PINN executes the official definitions
and a numerically validated checkpointable transcription of its monolithic
training loop.  No baseline receives the Tangent-Actor terminal construction
or MCV.

After all jobs finish, regenerate the auditable tables and figures with:

```bash
python -m invariant_region_mlp.experiments.original_pde_benchmark_comparison.collect_slurm_accounting \
  --start 2026-10-12
python -m invariant_region_mlp.experiments.original_pde_benchmark_comparison.aggregate_results
```

Native single-point and amortized function results are intentionally written
to different CSV files.  Raw checkpoints, common-bank predictions, official
working directories, and Slurm logs remain local; compact scalar JSON, CSV,
source audits, and figures are tracked.

