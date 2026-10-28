# Final component ablation

This package implements the preregistered ICML 2026 Tangent-Actor component
ablation. Frozen declarations live in `protocol.py`; measured results must not
be used to change them.

Reproduce the CPU gates:

```bash
python -m invariant_region_mlp.experiments.final_component_ablation.prepare_protocol
python -m invariant_region_mlp.experiments.final_component_ablation.preflight
python -m pytest -q invariant_region_mlp/experiments/final_component_ablation/test_correctness.py
```

The immutable four-GPU queue is `formal_queue.csv`. The validated Slurm entry
point is `slurm/formal_four_gpu.slurm`; each slot runs independent cells and
resumes only at completed Actor-update boundaries. A formal learned cell is:

```bash
python -m invariant_region_mlp.experiments.final_component_ablation.run_ablation \
  --pde joint_concave --variant a1_no_mcv --seed 0 --replicate 0
python -m invariant_region_mlp.experiments.final_component_ablation.run_ablation \
  --pde joint_concave --variant a1_no_mcv --seed 0 --replicate 1 \
  --evaluate-only
```

After all 80 evaluation records exist, build the audited tables, figures, and
report with the exact Slurm allocation IDs:

```bash
python -m invariant_region_mlp.experiments.final_component_ablation.build_artifacts \
  --job-ids INTERACTIVE_JOB_ID,REGULAR_JOB_ID
```

The historical 512-state oscillatory metrics are never reused in the primary
matrix. `preflight.py` freezes the mandated 960-state bank from the audited
ICML-final-suite generation rule before evaluation.
