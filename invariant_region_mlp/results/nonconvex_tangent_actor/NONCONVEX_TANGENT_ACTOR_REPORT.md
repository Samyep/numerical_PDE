# Nonconvex–Nonconcave Tangent-Actor Report

## Decision: **CONDITIONAL GO**

The experiment evaluates the unchanged undamped Terminal-Aware Tangent-Actor update with alpha=1 martingale control variates. Checkpoints are selected by an independently replicated, noise-corrected fixed-point consistency score before the numerical reference is loaded.

## Frozen selected results

| family | d | lambda | seed | Actor | Value RMSE | Relative L2 | Actor gradient-vector RMSE | MCV variance reduction | ESS |
|---|---:|---:|---:|---|---:|---:|---:|---:|---:|
| primary | 100 | 0.2 | 0 | a2 | 0.00150872 | 0.00303075 | 0.0607912 | 9.76× | 0.9736 |

## Interpretation

Completed 1 training runs and 1 post-freeze evaluations. The rank correlation between the reference-free consistency score and post-freeze true Value error is 1.0.

No monotone policy-improvement theorem or nonnegative Bregman-regret identity is used for these drivers. Linear-policy values are diagnostics only. The experiment supports claims only for smooth bounded-curvature drivers in the tested locally stable regime.

The Actor always receives the full ambient state. Active coordinates and reference gradients are used only in the post-freeze audit. These are ambient 100D/400D tests with a two-dimensional reference structure, not evidence for arbitrary effective dimension.

The independent comparison uses raw Full-History MLP and a plain soft-terminal PINN. Native scalar Deep BSDE does not estimate the required random-time solution function, the available Deep Picard adaptation embeds our terminal-aware Actor, and the available DPI implementation is tied to a control-form HJB; they are excluded rather than modified into favorable but method-mismatched baselines.

## Artifacts

- [Actor iteration behavior](figures/actor_iteration_behavior.png)
- [Fixed-point consistency versus true error](figures/fixed_point_vs_true_error.png)
- [Curvature-region coverage](figures/curvature_region_coverage.png)
- [MCV variance reduction](figures/mcv_variance_reduction.png)
- `ACTOR_CONSISTENCY_DIAGNOSTICS.csv`
- `ACTOR_ITERATION_RESULTS.csv`
- `MCV_VARIANCE_RESULTS.csv`
- `BASELINE_COMPARISON.csv`
- `COMPUTE_ACCOUNTING.csv`
- `REPRODUCIBILITY_MANIFEST.json`
