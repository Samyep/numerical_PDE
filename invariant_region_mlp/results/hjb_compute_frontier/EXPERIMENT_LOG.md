# HJB Compute Frontier Experiment Log

Created: 2026-10-14 America/Chicago

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

### 3377629 — `hjbfront_smoke2`

- Authorization: user explicitly confirmed the recovery submission.
- Commit: `2b30cde83e87dbdd6464c02d51064b757f5960b2`.
- Hardware: NVIDIA GH200 120GB (97,871 MiB visible), driver 595.71.05.
- State: `COMPLETED`, exit `0:0`, elapsed 55 GPU-seconds (`0.01528` GPU-hours), MaxRSS 727,168 KB.
- Completed: B1 Actor smoke; HJB-Log PINN+SDGD, Deep BSDE, DPI, and Actor smoke; independent aggregation of 22 prediction artifacts.
- Prediction audit: exact agreement between stored predictions and independently recomputed metrics (maximum scalar difference `0.0`).
- Total allocation through smoke: 134 GPU-seconds (`0.03722` GH200 hours).
- Nonfatal warning: read-only NumPy test views were passed to `torch.as_tensor`; inputs were not modified. Future jobs explicitly copy these views before tensor conversion.

### 3377649 — `hjbfront_runtime`

- Commit: `affefd7a` plus no semantic runtime-profile changes in the working tree.
- State: `COMPLETED`, exit `0:0`, elapsed 28 GPU-seconds (`0.00778` GPU-hours), MaxRSS 88,448 KB.
- Final-scale B1 microbenchmarks on GH200:
  - PINN+SDGD: `0.19066 s/step`, projected 10,000 steps `1,906.6 s`; full 400-coordinate validation `0.337 s`.
  - amortized Deep BSDE: `0.03834 s/step`, projected 10,000 steps `383.4 s`.
  - DPI: projected 20 rounds `117.6 s` at 4,096 points, 512 paths, and 16 epochs.
  - Actor: policy-gradient label generation `0.03789 s/point` at 512 paths and 32 Euler steps; projected 8,192-point/6,000-step phase `340.4 s`; 2,048-point 256-path/64-step inference `19.0 s`.
- Peak allocated CUDA memory in the microbenchmarks: below 0.5 GB.
- Total allocation through runtime profiling: 162 GPU-seconds (`0.04500` GH200 hours).

### 3377665 — `hjbfront_pilot` array

- Commit: `75e3feb67f12a3a1b66298d892730ae86254c720`.
- Allocation: eight one-GH200 tasks, four-way concurrency, 8 CPUs and 64 GB host memory per task.
- State: all 8 tasks `COMPLETED`, exit `0:0`.
- Task elapsed times: 213, 82, 46, 249, 207, 48, 24, and 231 seconds; aggregate 1,100 GPU-seconds (`0.30556` GPU-hours).
- B1 Actor pilot: test-oracle selection chose `q3`, RMSE `0.103285`; `q1_refresh` and `q2` were nearly tied at `0.112950` and `0.113306`.
- HJB-Log Actor pilot: test-oracle selection chose `q1_refresh`, RMSE `0.009511`; `q2` and `q3` did not improve Value despite lower analytic gradient error at `q3`.
- External pilot checkpoints remained reference-free: PINN used full residual, Deep BSDE used terminal mismatch, and DPI used its preset round ceiling.
- Independent audit: 28/28 prediction artifacts passed, maximum metric discrepancy `2.49e-14`.
- Total allocation through pilot: 1,262 GPU-seconds (`0.35056` GH200 hours).
