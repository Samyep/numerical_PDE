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

### 3378026 — high-quality B1 Actor `q0`

- Script: `experiments/hjb_compute_frontier/slurm/high_actor_continue.slurm`.
- State: `COMPLETED`, exit `0:0`, elapsed 88 GPU-seconds; Python step 41
  seconds, MaxRSS 1,645,952 KB.
- Produced the width-512 terminal-pathwise `q0` checkpoint.  RMSE was
  `0.170516`; path-training time was 13.37 seconds.

### 3378260 — high-quality B1 Actor `q1`

- State: `COMPLETED`, exit `0:0`, elapsed 1,338 GPU-seconds (`0.37167`
  GPU-hours); Python step 1,294 seconds, MaxRSS 1,611 MiB.
- Configuration: 16,384 training points, 1,024 validation points, 1,024 paths
  per label, 64 label Euler steps, 12,000 optimizer steps, width 512.
- The selected optimizer checkpoint was step 12,000/12,000.  `q1` RMSE was
  `0.053064`, 40.8% below the same-seed frozen-budget `q1`.

### 3378329 — matched fresh-label `q1_refresh`

- State: `COMPLETED`, exit `0:0`, elapsed 1,378 GPU-seconds (`0.38278`
  GPU-hours); Python step 1,347 seconds, MaxRSS 1,776,448 KB.
- Used the same high-quality budget and training-state distribution as the
  planned true `q2`, with fresh independent `q0` policy-evaluation labels.
- The selected optimizer checkpoint was step 12,000/12,000.  RMSE was
  `0.033069`, a further 37.7% reduction from high-quality `q1`.

### 3379826 — matched true `q2`

- State: `COMPLETED`, exit `0:0`, elapsed 1,342 GPU-seconds (`0.37278`
  GPU-hours); Python step 1,327 seconds, MaxRSS 2,220,480 KB.
- Used the same 16,384/1,024/1,024/64/12,000 budget and matched states as
  `q1_refresh`, but generated labels from the frozen `q1` policy evaluation.
- RMSE was `0.040800`, 23.4% above `q1_refresh`, despite improving direct
  gradient relative error from `0.2450` to `0.2305`.

### 3379966 — high-quality true `q3` and final candidate evaluation

- The QOS had only 5.691 one-GH200 minutes available at submission.  The job
  used the declared smaller `q3` budget (8,192 points, 512 paths, 32 steps,
  6,000 optimizer steps) and a five-minute wall limit.
- The numerical `srun` completed in 283 seconds with exit `0:0`.  It atomically
  wrote the `q3` checkpoint, all five prediction arrays, and the final result
  JSON.  Slurm marked the outer job `TIMEOUT` at 306 allocated seconds because
  batch cleanup ended six seconds after the wall limit; no numerical artifact
  was partial and no retry was performed.
- `q3` RMSE was `0.036694`: 10.1% below `q2` but 11.0% above `q1_refresh`.
  The test oracle therefore retained `q1_refresh`; policy validation and the
  direct-gradient diagnostic favored `q3`.  This disagreement did not justify
  a `q4` run.
- Independent recomputation of the five high-profile prediction artifacts
  agreed with the run JSON to `1.67e-16` maximum scalar discrepancy.

### Final accounting and seed-extension status

- The five high-quality B1 stages used 4,452 allocated GH200-seconds
  (`1.23667` GPU-hours).
- Total allocation across smoke, recovery, profiling, pilot, frozen final,
  HJB-Log evaluation sensitivity, and high-quality B1 diagnostics was 28,421
  GH200-seconds (`7.89472` GPU-hours).
- After job 3379966, QOS accounting showed 209,999,409 of 210,000,000 billing
  minutes consumed, leaving 591 billing-minutes, or 0.591 one-GH200 wall
  minutes.  The validated seeds 3--4 extension was therefore not submitted.
  The main study remains three complete seeds for every method/PDE cell, as
  permitted by the original cost-contingent protocol.
