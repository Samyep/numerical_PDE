# Original-PDE benchmark comparison

This experiment compares Terminal-Aware Tangent-Actor with direct martingale
control variates against unenhanced methods from the official DeepBSDE,
DeepPicardIteration, and SDGD-PINN repositories.  It separates native
single-point tasks from amortized function approximation and freezes Actor
selection before reading final-test references.

The generated test banks, checkpoints, official working directories, and
Slurm logs remain local.  Their hashes and scalar outcomes are recorded in the
tracked result manifests and tables.  `reproducibility_metadata.json` records
the software/hardware environment, source commits, bank hashes, Slurm resource
contract, and seed policy; `protocol.json` records the frozen scientific
protocol.

Run correctness and source-independent reference preparation with:

```bash
source /sw/user/python/miniforge3-pytorch-2.10.0/etc/profile.d/conda.sh
conda activate base
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
   three seeds and 32/64/128-step evaluation.  The seed launcher accepts a
   space-separated `SEEDS` export so an upper-budget two-round Actor can be
   reproduced in separate allocations without changing its configuration.
   The evaluation launcher likewise accepts `SEEDS`, `DIAGNOSTIC_STEPS`, and
   `RUN_SCREEN`; their defaults reproduce the full frozen protocol and only
   permit lossless partitioning across allocations.

Official baseline launchers retain their source algorithms.  The Deep BSDE
launcher is an audited PyTorch compatibility port; DPI, PINN-HTE, and D-DBSDE
run through the authors' package; SDGD-PINN executes the official definitions
and a numerically validated checkpointable transcription of its monolithic
training loop.  No baseline receives the Tangent-Actor terminal construction
or MCV.

The Full-History launcher keeps the frozen six-configuration, three-seed grid
as its default.  Its `BENCHMARKS`, `CONFIGURATIONS` (for example `3:24`), and
`SEEDS` exports only partition that same grid across allocations when the
two-hour queue limit requires it.  Completed atomic result files are skipped
by default, so a timed-out partition can be resumed without rerunning finished
configurations; set `SKIP_EXISTING=0` only for an intentional reproduction.

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

