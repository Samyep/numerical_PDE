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

### 3377537 — `hjbfront_smoke`

- Script: `experiments/hjb_compute_frontier/slurm/smoke_profile.slurm`
- Allocation: one NVIDIA GH200, 8 CPUs, 64 GB host memory, 30-minute limit.
- State: `FAILED`, exit `1:0`, elapsed 79 GPU-seconds (`0.02194` GPU-hours), MaxRSS 661,312 KB.
- Completed before failure: B1 PINN+SDGD, Deep BSDE, and DPI smoke runs, including checkpoints and independently recomputed test metrics.
- Failure: Actor smoke declared 8 paths but inherited a 64-path chunk. The explicit divisibility guard stopped before allocating the invalid tensor. HJB-Log was not started.
- Diagnosis: configuration plumbing only; no OOM, NaN, or training instability.
- Fix: every Actor label/evaluation call now passes `path_chunk=min(64, paths)`, with a regression test for the 8-path setting. Focused suite now has 7 passing tests.
- Recovery script prepared: `experiments/hjb_compute_frontier/slurm/smoke_continue.slurm`; it skips the three completed B1 baseline cells.
