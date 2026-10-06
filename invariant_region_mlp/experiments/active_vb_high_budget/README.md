# Active-gradient VB high-budget experiment

This directory implements the published `sigma=sqrt(2)` SCaSML
gradient-dependent nonlinear benchmark with `z=sigma*grad(u)`. It intentionally
does not use the public repository's `sigma=0.25`, float16 arithmetic, terminal
EBL normalization bug, or heuristic output clipping.

Key commands (run from the `numerical_PDE` repository root):

```powershell
# Equation/projection/RNG/base-case tests
python invariant_region_mlp/experiments/active_vb_high_budget/test_vb_equation.py

# Signal gate and full d=20 target-grid pilot
python invariant_region_mlp/experiments/active_vb_high_budget/run_vb_pilot.py --workers 8 --resume

# Resumable plan-driven main run
python invariant_region_mlp/experiments/active_vb_high_budget/run_vb_high_budget.py main `
  --plan invariant_region_mlp/results/active_vb_high_budget/main_plan.json `
  --skip-signal-gate --workers 8 --resume

# Aggregate main plus expensive consistency extension
python invariant_region_mlp/experiments/active_vb_high_budget/analyze_vb_high_budget.py `
  --stage main,high_budget_extension `
  --selection-stage tuning,tuning_refine --verdict A

# Pairing/integrity audit and compact raw archives
python invariant_region_mlp/experiments/active_vb_high_budget/validate_vb_results.py
python invariant_region_mlp/experiments/active_vb_high_budget/compact_vb_results.py --stage main
```

Batch-IR is box-compatible: it first maps each negative z coordinate to zero,
then uses one common alpha per parent sibling group so the largest remaining
coordinate is at most `sigma/4`. The u coordinates are clipped independently.
This rule is identity on a feasible batch and is not a radial Euclidean box
projection.

See `docs/VB_HIGH_BUDGET_REPORT.md` for the full protocol and verdict.
