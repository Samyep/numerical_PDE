# Legacy benchmark high-budget rescue

This directory contains the float64, corrected-EBL rescue study for the 100D
nonlinear Funding benchmark and the 100--160D Rosenbrock HJB benchmark.

The primary comparison is always Raw full-history MLP versus Samplewise
certified IR-MLP. Suppression and tighter/shrinkage controls are restricted to
representative diagnostics. No method applies heuristic output clipping, and
returned root states are never projected.

Main entry points:

- `funding_rescue.py`: resumable block runner for paired Funding roots.
- `hjb_rescue.py`: resumable HJB stage runner and exact legacy reproduction.
- `analyze_rescue.py`: CSV aggregation, work-normalized figures, and reports.
- `validate_rescue.py`: bitwise pairing and artifact integrity audit.
- `test_rescue.py`: equation, certificate, projection, EBL, and RNG tests.

The HJB coefficient construction deliberately uses JAX `PRNGKey(0/1)` because
that is part of the validated historical benchmark. SciPy supplies the stable
Gauss--Laguerre reference. Scientific recursion and results use NumPy float64.

From the repository root:

```powershell
python invariant_region_mlp/experiments/rescue_high_budget/test_rescue.py
python invariant_region_mlp/experiments/rescue_high_budget/analyze_rescue.py
python invariant_region_mlp/experiments/rescue_high_budget/validate_rescue.py
```

See `docs/FUNDING_HIGH_BUDGET_RESCUE.md`,
`docs/HJB_HIGH_BUDGET_RESCUE.md`, and
`docs/OLD_BENCHMARK_RESCUE_SUMMARY.md` for conclusions.
