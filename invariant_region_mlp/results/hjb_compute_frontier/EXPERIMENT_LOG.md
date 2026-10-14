# HJB Compute Frontier Experiment Log

Created: 2026-10-14 America/New_York

Branch: `ir-mlp-hjb-strong-baselines-compute-frontier`

Starting commit: `56f6448e1557118f8cc6e51ed038ff03512c84d2`

## Preflight

- CPU correctness audit: passed.
- Focused tests: 6 passed.
- B1 analytic residual maximum: `2.22e-16`.
- HJB-Log quadrature versus independent adaptive integration: Value `2.30e-13`, z `1.53e-15` maximum absolute error.
- New 2,048-point test-bank hashes and 256-point policy-validation-bank hashes recorded in `validation.json`.
- External baseline settings frozen before GPU results; one configuration per method, no test-based baseline tuning.
- Actor selection override recorded as `test_rmse_oracle_allowed_by_user`.

## Slurm jobs

No job submitted yet. The next entry will record the GH200 smoke/profile job ID, script hash, state, and accounting.
