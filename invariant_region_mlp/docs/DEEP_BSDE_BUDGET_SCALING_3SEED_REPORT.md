# Deep BSDE hyperparameter-budget scaling (three seeds)

## Outcome

- Formal terminal records: **1920 / 1920**.
- Numerical failures: **0**.
- Allocation preflight: **ALLOCATION_GATE_PASS**. Its 92.18 H100-hour projection covered the original larger grid; the active scope was subsequently narrowed before truth opening.
- Actual allocation (including recorded failures/retries): **30.34 H100 GPU-hours**.
- Frozen DPI/PINN predictions on the paired 64-point subset: **COMPLETE**, **108 / 108** seed-level artifacts (no retraining).
- Scope amendment: **C0/C1/C2 only**, 30 tasks / 1920 terminal records. Joint J1/J2/J3 scaling is excluded from analysis and tables at the user's request before the truth gate opened.
- Excluded work is retained for audit and charged to total allocation: one completed J1 bundle and four terminal records from an interrupted J3 bundle; no excluded accuracy is computed.
- Infrastructure incidents: **3**. Two CUDA-unavailable incidents affected task launches; a later ICE quota incident interrupted eight partially complete tasks. Exact frozen-task retries preserve completed query records, and every failed/cancelled allocation is retained in accounting.
- Resume recovery: **PASS**. The 17 superseded infrastructure markers are retained in an audit archive; the final numerical-failure count is **0**.

Deep BSDE is a native scalar-initial-value method here: every one of the 64 query states receives a fresh model. Reported Deep BSDE compute is therefore the cumulative active training time for all 64 models, averaged across training seeds—not a one-query runtime.

## Individual sweeps

| PDE | Sweep | Level | Hyperparameter | RMSE mean ± seed SD | RelMAE | Total 64-query min | Historical TA-reference status |
|---|---|---:|---|---:|---:|---:|---|
| b6_fullrank | C0 | 1 | updates=1000 | 0.164782 ± 0.000876 | 0.13149 | 21.5 | TARGET_NOT_REACHED_WITHIN_EVALUATED_BUDGET |
| b6_fullrank | C0 | 2 | updates=2000 | 0.160081 ± 0.00101 | 0.12721 | 43.0 | TARGET_NOT_REACHED_WITHIN_EVALUATED_BUDGET |
| b6_fullrank | C0 | 3 | updates=4000 | 0.144701 ± 0.00061 | 0.11441 | 85.8 | TARGET_NOT_REACHED_WITHIN_EVALUATED_BUDGET |
| b6_fullrank | C1 | 1 | batch_size=32 | 0.161231 ± 0.00151 | 0.12897 | 42.2 | TARGET_NOT_REACHED_WITHIN_EVALUATED_BUDGET |
| b6_fullrank | C1 | 2 | batch_size=64 | 0.160081 ± 0.00101 | 0.12721 | 43.0 | TARGET_NOT_REACHED_WITHIN_EVALUATED_BUDGET |
| b6_fullrank | C1 | 3 | batch_size=128 | 0.158102 ± 0.000256 | 0.12602 | 45.9 | TARGET_NOT_REACHED_WITHIN_EVALUATED_BUDGET |
| b6_fullrank | C2 | 1 | time_steps=10 | 0.147119 ± 0.000689 | 0.11717 | 22.1 | TARGET_NOT_REACHED_WITHIN_EVALUATED_BUDGET |
| b6_fullrank | C2 | 2 | time_steps=20 | 0.160081 ± 0.00101 | 0.12721 | 43.0 | TARGET_NOT_REACHED_WITHIN_EVALUATED_BUDGET |
| b6_fullrank | C2 | 3 | time_steps=40 | 0.167064 ± 0.000825 | 0.13363 | 81.9 | TARGET_NOT_REACHED_WITHIN_EVALUATED_BUDGET |
| dense_fourier | C0 | 1 | updates=1000 | 0.0766677 ± 0.00328 | 0.046615 | 22.6 | TARGET_NOT_REACHED_WITHIN_EVALUATED_BUDGET |
| dense_fourier | C0 | 2 | updates=2000 | 0.0741138 ± 0.000809 | 0.044534 | 45.1 | TARGET_NOT_REACHED_WITHIN_EVALUATED_BUDGET |
| dense_fourier | C0 | 3 | updates=4000 | 0.0594877 ± 0.00139 | 0.034757 | 89.9 | TARGET_NOT_REACHED_WITHIN_EVALUATED_BUDGET |
| dense_fourier | C1 | 1 | batch_size=32 | 0.0746 ± 0.00191 | 0.04456 | 44.5 | TARGET_NOT_REACHED_WITHIN_EVALUATED_BUDGET |
| dense_fourier | C1 | 2 | batch_size=64 | 0.0741138 ± 0.000809 | 0.044534 | 45.1 | TARGET_NOT_REACHED_WITHIN_EVALUATED_BUDGET |
| dense_fourier | C1 | 3 | batch_size=128 | 0.0669422 ± 0.000758 | 0.04083 | 44.3 | TARGET_NOT_REACHED_WITHIN_EVALUATED_BUDGET |
| dense_fourier | C2 | 1 | time_steps=10 | 0.0650344 ± 0.00158 | 0.038846 | 23.2 | TARGET_NOT_REACHED_WITHIN_EVALUATED_BUDGET |
| dense_fourier | C2 | 2 | time_steps=20 | 0.0741138 ± 0.000809 | 0.044534 | 45.1 | TARGET_NOT_REACHED_WITHIN_EVALUATED_BUDGET |
| dense_fourier | C2 | 3 | time_steps=40 | 0.0767117 ± 0.00287 | 0.046433 | 88.2 | TARGET_NOT_REACHED_WITHIN_EVALUATED_BUDGET |

C2 changes Euler time-discretization steps; it is not an optimizer-step sweep. C0 checkpoints 1000/2000/4000 come from a single continuous trajectory. The common 2000/64/20 checkpoint is reused exactly for C0-level2, C1-level2, and C2-level2.

## Observed budget trends

- **B6 Full-Rank HJB:** increasing updates from 1000 to 4000 reduced RMSE by 12.2%; batch 32 to 128 reduced it by 1.9%. Increasing Euler steps from 10 to 40 changed RMSE by +13.6% (nonmonotonic/worse). The best observed setting was C0-level 3 (updates=4000), RMSE 0.144701 at 85.8 cumulative minutes for 64 queries; it remained 201.0x the descriptive historical TA-reference RMSE.
- **Dense Fourier:** increasing updates from 1000 to 4000 reduced RMSE by 22.4%; batch 32 to 128 reduced it by 10.3%. Increasing Euler steps from 10 to 40 changed RMSE by +18.0% (nonmonotonic/worse). The best observed setting was C0-level 3 (updates=4000), RMSE 0.0594877 at 89.9 cumulative minutes for 64 queries; it remained 33.8x the descriptive historical TA-reference RMSE.

## Interpretation boundaries

- A lower observed RMSE at a larger setting is a genuine observed improvement only at the evaluated configurations; nonmonotonic sequences are reported without smoothing.
- `TARGET_NOT_REACHED_WITHIN_EVALUATED_BUDGET` refers only to the pre-existing audited TA reference RMSE used to frame the reviewer question. That TA value was measured on 960 points, whereas Deep BSDE uses the frozen 64-point subset, so it is descriptive rather than a paired hypothesis test.
- A configuration with any final `numerical_failure` record is labelled numerical instability and is not recomputed on a finite subset. Infrastructure failures are separately retained and retried under the unchanged frozen configuration.
- The experiment does not establish that Deep BSDE cannot converge under arbitrarily larger budgets.
- The main ICML table keeps historical DPI/PINN values unchanged and states the 960-vs-64 protocol limitation. The supplementary table uses their saved predictions on the identical 64 points when available, without retraining.

## Provenance and outputs

- Report-generation commit: `325792538bdf41d6857068718d781f57064b8c3b`
- Main table: `results/deep_bsde_budget_scaling_3seed/tables/deep_bsde_budget_scaling_main.tex`
- Seed metrics: `results/deep_bsde_budget_scaling_3seed/deep_bsde_seed_metrics.csv`
- Resource accounting: `results/deep_bsde_budget_scaling_3seed/resource_accounting.json`
- Artifact audit: `results/deep_bsde_budget_scaling_3seed/artifact_audit.json`
