# B1 actor bottleneck diagnostics

## Numerical diagnosis

The oracle actor produces a large and broad accuracy improvement, so learned-actor error is a major current bottleneck. Oracle Actor beats the learned actor in 15/15 cells, with geometric-mean Oracle/Learned RMSE ratio 0.2629.

The learned actor improves geometric-mean RMSE over Pathwise + Baseline without an actor, although the cell count shows how consistent that benefit is. Learned Actor beats the no-actor variant in 13/15 cells, with geometric-mean Learned/No-Actor RMSE ratio 0.5593. Oracle Actor beats the no-actor variant in 15/15 cells, with geometric-mean Oracle/No-Actor ratio 0.1471.

## Scope and pairing

These are supplemental diagnostics and do not replace any primary six-method result. Every cell uses the frozen Stage B PDE, 960 held-out points, `n=3`, `M=16`, and 10 repetitions. Learned Actor, Oracle Actor, and Pathwise + Baseline without Actor use the same random seed for each `(d,k,repetition)`, enabling paired comparisons.

Oracle Actor evaluates `q(t,x)=z*(t,x)` at every recursive child and in the matching level-1 baseline. It uses forbidden oracle information and is only a diagnostic upper bound. The no-actor variant uses the original quadratic driver with the pathwise terminal gradient and original-driver level-matched baseline.

## Value RMSE

| d | k | Learned Actor | Oracle Actor | Pathwise + Baseline, no Actor | Oracle/Learned | Learned/No Actor | Actor relative error |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 | 1 | 0.1807 | **0.1231** | 0.3885 | 0.681 | 0.465 | 0.205 |
| 100 | 5 | 0.1719 | **0.06439** | 0.272 | 0.375 | 0.632 | 0.369 |
| 100 | 10 | 0.1747 | **0.06599** | 0.2628 | 0.378 | 0.665 | 0.410 |
| 100 | 20 | 0.1317 | **0.04231** | 0.133 | 0.321 | 0.990 | 0.469 |
| 100 | 40 | 0.1268 | **0.03234** | 0.0914 | 0.255 | 1.387 | 0.589 |
| 200 | 1 | 0.2216 | **0.1219** | 0.7253 | 0.550 | 0.305 | 0.251 |
| 200 | 5 | 0.2728 | **0.06653** | 0.5565 | 0.244 | 0.490 | 0.474 |
| 200 | 10 | 0.3237 | **0.06544** | 0.5167 | 0.202 | 0.627 | 0.532 |
| 200 | 20 | 0.1836 | **0.0418** | 0.2721 | 0.228 | 0.675 | 0.599 |
| 200 | 40 | 0.1864 | **0.03068** | 0.1716 | 0.165 | 1.086 | 0.741 |
| 400 | 1 | 0.3302 | **0.1215** | 1.42 | 0.368 | 0.233 | 0.314 |
| 400 | 5 | 0.3426 | **0.06807** | 1.118 | 0.199 | 0.306 | 0.611 |
| 400 | 10 | 0.3544 | **0.06593** | 1.053 | 0.186 | 0.337 | 0.697 |
| 400 | 20 | 0.2778 | **0.0407** | 0.5213 | 0.146 | 0.533 | 0.778 |
| 400 | 40 | 0.2472 | **0.03067** | 0.3556 | 0.124 | 0.695 | 0.883 |

## Actor quality relation

Across the 15 cells, the Spearman correlation between held-out learned-actor relative gradient error and the Learned/Oracle value-RMSE ratio is 0.968 (two-sided p=3.6e-09). For each cell, `actor_diagnostic_comparison.csv` also reports the implied mean tangent-driver gap `0.5 * actor_gradient_rmse^2` at the exact gradient.

## Training-budget diagnosis

Across all 15 actors, noisy-label generation took only 0.572 seconds and network training took 10.016 seconds. Per actor, network training had median 0.561 seconds (range 0.382--1.144) and median 425 executed optimizer steps (range 275--925); the median selected checkpoint was step 225.

The resulting held-out actor gradient RMSE had median 0.958; relative gradient error had median 0.532 and maximum 0.883. The noisy training labels themselves had median RMSE 2.342 against the post-training exact diagnostic. Together with the oracle gap, these measurements show that the current actor data/training allocation is too small for the accuracy demanded by the recursive solver.

## Measured cost

The existing 15 learned actors used 10.588 seconds for label generation plus training and 35.356 seconds for the common MLP workload. Oracle diagnostic inference used 35.163 seconds; Pathwise + Baseline without Actor used 34.015 seconds. Oracle cost is descriptive only because exact-gradient access is unavailable in a practical solver.

The complete diagnostic allocation, including 18 GPU-side correctness tests, report generation, and scheduler overhead, used 108 GH200 seconds (0.0300 GH200 hours).

## Artifacts

- `actor_diagnostic_results.csv`: all 45 learned/oracle/no-actor method-cells.
- `actor_diagnostic_comparison.csv`: paired ratios, differences, and actor-quality metrics.
- `summary.json`: aggregate counts, geometric means, correlation, and runtime.
- `figures/actor_diagnostic_value_rmse.*` and `figures/actor_quality_vs_oracle_gap.*`.
