# Actor tuning report

Status: pilot complete; three-seed final study pending.

The ambient-output Actor has no rank projection on either PDE. The declared study compares terminal-pathwise `q0`, one policy improvement `q1`, a matched fresh-label continuation `q1_refresh`, a true `q2`, and a true `q3`. Each controlled-policy label differentiates the full Euler trajectory with the source Actor frozen while preserving its state derivative.

The pilot uses 2,048 terminal-label points for `q0`, then 4,096 points, 256 paths per label, 16 Euler steps, and 4,000 optimizer steps for each subsequent phase. The final preset raises this to 4,096/256/6,000 for `q0` and 8,192/512/32/6,000 for policy-evaluation phases. A targeted high-quality B1 branch declares width 512, 16,384 points, 1,024 paths, 64 steps, and 12,000 optimizer steps for `q1` and matched `q1_refresh`; `q2` and `q3` retain the final-scale 8,192/512/32/6,000 budget. This isolates the value of much better first-update supervision while keeping later-iteration cost controlled. The HJB-Log high branch was cancelled before allocation after the final runs showed that its Value RMSE was already dominated by final MC inference noise.

Per the user's explicit override, Actor candidates are ranked by test Relative L2/RMSE. The report will show both model-path cost and total exploration cost so this oracle selection is not presented as a deployable reference-free procedure. External baselines do not receive the same selection privilege.

## Pilot results

The pilot ran on NVIDIA GH200 GPUs with seed 0. Candidate evaluations within each PDE used common Brownian random numbers. All actors used ambient 400-dimensional outputs without a rank projection.

| PDE | Actor | Value RMSE | Relative L2 | Bias | gradient relative error | cumulative path training (s) |
|---|---:|---:|---:|---:|---:|---:|
| B1 (400,40) | q0 | 0.217052 | 0.071111 | +0.188670 | 0.7852 | 6.7 |
| B1 (400,40) | q1 | 0.146121 | 0.047872 | +0.127411 | 0.6272 | 54.4 |
| B1 (400,40) | q1_refresh | 0.112950 | 0.037005 | +0.095793 | 0.5343 | 102.1 |
| B1 (400,40) | q2 | 0.113306 | 0.037121 | +0.095990 | 0.5318 | 102.4 |
| B1 (400,40) | **q3** | **0.103285** | **0.033838** | **+0.078973** | **0.4566** | 150.1 |
| HJB-Log (400) | q0 | 0.012310 | 0.001926 | +0.004107 | 1.0777 | 7.0 |
| HJB-Log (400) | q1 | 0.009987 | 0.001562 | +0.003398 | 0.6166 | 57.7 |
| HJB-Log (400) | **q1_refresh** | **0.009511** | **0.001488** | **+0.002370** | 0.4742 | 107.8 |
| HJB-Log (400) | q2 | 0.009529 | 0.001491 | +0.002353 | 0.4739 | 107.9 |
| HJB-Log (400) | q3 | 0.009615 | 0.001504 | +0.002601 | **0.4203** | 157.8 |

The B1 pilot supports retaining multiple policy iterations: `q3` reduced RMSE by 29.3% from `q1` and by 8.6% from the matched `q1_refresh` continuation. The direct `q2` update and continuation were tied, so the gain appears only after the next update at this pilot budget. On HJB-Log, extra policy iterations improved the analytic gradient diagnostic but did not improve Value beyond `q1_refresh`; the test-oracle rule therefore selected `q1_refresh`.

The HJB-Log truth has mean 6.393 and standard deviation only 0.02395. Its small Relative L2 can hide poor predictors: a validation-fitted constant has test RMSE 0.02398. The selected Actor RMSE 0.00951 is also 0.385 times the truth standard deviation, so its gain is not just a mean-value effect.

The complete eight-task pilot consumed 1,100 allocated GH200-seconds (0.3056 GPU-hours). All tasks completed successfully. Independent recomputation of all 28 stored prediction artifacts agreed with reported metrics to at most `2.49e-14`.
