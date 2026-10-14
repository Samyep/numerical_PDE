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

### 3377699 — `hjbfront_final` array

- Commits: tasks 0–1 recorded `3d1bb40c`; tasks 2–23 recorded `c3b992ac`. The only intervening commit added the high-budget Slurm launcher and did not change training code.
- Allocation: 24 one-GH200 tasks for 3 seeds × 2 PDEs × 4 methods, at most four concurrent, 8 CPUs and 64 GB host memory per task.
- State: all 24 tasks `COMPLETED`, exit `0:0`; no retrials.
- Aggregate allocation: 21,197 GPU-seconds (`5.88806` GPU-hours). Maximum Python-step MaxRSS was 2,268,288 KB.
- The first 18 tasks requested 75 minutes. After measured maxima were 32.0 minutes for PINN and 25.3 minutes for Actor, pending tasks 18–23 were safely reduced to 45 minutes to satisfy the group billing-minute scheduler limit. Algorithms and numerical budgets were unchanged.
- Independent aggregation audited 96 prediction artifacts; maximum recomputed-metric discrepancy was `1.42e-14`.
- B1 mean RMSE across seeds: Actor `0.06301`, DPI `0.06841`, PINN+SDGD `0.08254`, amortized Deep BSDE `0.31655`.
- HJB-Log mean RMSE across seeds: Actor `0.004974`, amortized Deep BSDE `0.02035`, DPI `0.06378`, PINN+SDGD `0.84205`.
- Total allocation through the frozen final study: 22,459 GPU-seconds (`6.23861` GH200 hours), including smoke, recovery, profiling, pilot, and final runs.

### 3377705 — `hjbfront_actorhi` array

- Submitted as two high-quality Actor tasks, one per PDE.
- HJB-Log task 1 was cancelled while still pending and consumed zero GPU-seconds: the final three-seed study showed RMSE near `0.0050` for every Actor round with per-point MC standard error near `0.0044`, so more Actor training would not address the observed bottleneck.
- B1 task 0 was also cancelled while pending, with zero GPU allocation.  The
  array reservation continued to trigger `QOSGrpBillingMinutes` even after its
  wall limit was reduced.  It was replaced by phase-resumable single jobs; no
  completed computation was discarded.
- The replacement keeps `q2/q3` at the validated final budget and concentrates
  the much larger 16,384-point/1,024-path/64-step/12,000-step budget on `q1`
  and matched `q1_refresh`.

### 3378026–3378028 — staged `hjbfront_hi_stage` launch

- Jobs 3378026, 3378027, and 3378028 were initially submitted as a three-stage
  dependency chain.  Slurm counted the two unfulfilled 100-minute dependency
  jobs against the group billing-minute limit and blocked stage 1.
- Jobs 3378027 and 3378028 were cancelled before allocation, with zero elapsed
  GPU-seconds.  Job 3378026 was retained; later stages are submitted only after
  their predecessor completes.  This scheduling change does not alter any
  numerical setting or recompute a completed phase.

### 3378042 — `hjbfront_evalmc` array

- Commit: `68f1e297`; three serial one-GH200 tasks evaluated frozen HJB-Log
  Actor candidates.  No training was performed.
- State: all 3 tasks `COMPLETED`, exit `0:0`; elapsed times 507, 488, and 515
  seconds, totaling 1,510 GPU-seconds (`0.41944` GPU-hours).  Maximum Python
  MaxRSS was 980,224 KB.
- Each seed evaluated all five candidates at 1,024 paths/64 Euler steps and at
  256 paths/128 steps, using candidate-level common random numbers.
- Independent audit: 45/45 prediction files reproduced every stored scalar
  exactly.  At 1,024 paths, `q3` was selected in all seeds with RMSE
  `0.003094 +/- 0.000092`; fixed `q0` had `0.003233 +/- 0.000037`.
- Increasing paths reduced RMSE materially, while doubling time steps at 256
  paths did not.  This confirms final Monte Carlo variance as the main visible
  HJB-Log evaluation bottleneck.
- Total allocated compute through this diagnostic: 23,969 GPU-seconds
  (`6.65806` GH200 hours).
