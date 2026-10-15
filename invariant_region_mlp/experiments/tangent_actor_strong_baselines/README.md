# Tangent-Actor strong-baseline pilot

This experiment compares Terminal-Aware Tangent-Actor Iteration with amortized Deep BSDE, SDGD-PINN, Deep Picard Iteration, SCaSML-L1, and inference-only Full-History MLP on two semilinear PDEs with independent two-dimensional references.

The frozen configuration is in [`protocol.json`](protocol.json), and the inspected official sources are in [`OFFICIAL_SOURCE_REVIEW.md`](OFFICIAL_SOURCE_REVIEW.md). Heavy references, checkpoints, raw traces, predictions, and Slurm logs stay in the workspace; compact CSVs, figures, and the report are versioned.

From the repository `source/` directory:

```bash
python -m invariant_region_mlp.experiments.tangent_actor_strong_baselines.prepare
pytest invariant_region_mlp/experiments/tangent_actor_strong_baselines/test_strong_baselines.py -q
sbatch invariant_region_mlp/experiments/tangent_actor_strong_baselines/slurm/calibrate.slurm
sbatch invariant_region_mlp/experiments/tangent_actor_strong_baselines/slurm/main_array.slurm
sbatch invariant_region_mlp/experiments/tangent_actor_strong_baselines/slurm/evaluate.slurm
python -m invariant_region_mlp.experiments.tangent_actor_strong_baselines.aggregate
python -m invariant_region_mlp.experiments.tangent_actor_strong_baselines.make_report
```

`main_array.slurm` and `evaluate.slurm` each request one four-GH200 node and launch four isolated one-GPU steps for the two benchmarks and two seeds. Every timing inside a method uses explicit CUDA synchronization.
