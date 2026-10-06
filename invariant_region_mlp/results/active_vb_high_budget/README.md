# Active VB result artifacts

The experiment was evaluated in float64. The committed result set contains:

- `repetition_metrics.csv`: one scalar-diagnostic row for every final main or
  high-budget-extension repetition;
- `repetition_metrics_{pilot,tuning,tuning_refine}.csv`: raw per-repetition
  pilot and validation-only tuning diagnostics;
- `compact_repetitions/*.npz`: float64 per-point value predictions, nonlinear
  value corrections, and gradient-error norms for every final repetition;
- `work_summary.csv` and `full_summary.json`: aggregate results;
- `signal_strength.json`, `tuning_choices.json`, manifests, fixed main/pilot
  point sets, figures, and `validation_audit.json`.

The detailed local `repetitions/` tree is approximately 2.6 GB and additionally
contains every predicted gradient coordinate, exact state, and root terminal
state. It is intentionally retained locally and ignored by Git; it was not
deleted. It can be reproduced from the committed scripts and manifests. The
compact archives preserve the final raw data needed to recompute value metrics,
nonlinear-correction statistics, and pointwise gradient-error norms without
exceeding ordinary Git hosting limits.

Reproduce or resume the target-grid main run with:

```powershell
python invariant_region_mlp/experiments/active_vb_high_budget/run_vb_high_budget.py main `
  --plan invariant_region_mlp/results/active_vb_high_budget/main_plan.json `
  --skip-signal-gate --workers 8 --resume
```

Run the integrity audit with:

```powershell
python invariant_region_mlp/experiments/active_vb_high_budget/validate_vb_results.py
```
