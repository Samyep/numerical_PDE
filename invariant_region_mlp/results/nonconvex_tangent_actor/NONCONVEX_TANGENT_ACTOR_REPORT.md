# Nonconvex–Nonconcave Tangent-Actor Report

## Decision: **CONDITIONAL GO**

The experiment evaluates the unchanged undamped Terminal-Aware Tangent-Actor update with alpha=1 martingale control variates. Checkpoints are selected by an independently replicated, noise-corrected fixed-point consistency score before the numerical reference is loaded.

## Frozen selected results

| family | d | lambda | seed | Actor | Value RMSE | Relative L2 | Actor gradient-vector RMSE | MCV variance reduction | ESS |
|---|---:|---:|---:|---|---:|---:|---:|---:|---:|
| coupled | 100 | 0.2 | 0 | a2 | 0.00198365 | 0.00422741 | 0.0630903 | 7.44× | 0.9504 |
| primary | 100 | 0.0 | 0 | a2 | 0.00278866 | 0.00556405 | 0.0525457 | 11.99× | 1.0000 |
| primary | 100 | 0.05 | 0 | a2 | 0.00141132 | 0.00283456 | 0.0519344 | 14.65× | 0.9984 |
| primary | 100 | 0.2 | 0 | a4 | 0.000684943 | 0.00137593 | 0.007139 | 51.75× | 0.9733 |
| primary | 100 | 0.2 | 0 | a2 | 0.00150872 | 0.00303075 | 0.0607912 | 9.76× | 0.9736 |
| primary | 100 | 0.5 | 0 | a2 | 0.00954707 | 0.018995 | 0.0970053 | 2.22× | 0.8513 |
| primary | 400 | 0.2 | 0 | a2 | 0.00629103 | 0.0124982 | 0.153572 | 2.82× | 0.9781 |

## Frozen endpoint comparison

This table uses the reference-free selected main Actor, the predeclared highest Full-History endpoint `n=3,M=64`, and the PINN final step. No baseline row is selected by test RMSE.

| method | frozen configuration | seeds | Value RMSE mean ± sd | mean total seconds |
|---|---|---:|---:|---:|
| Tangent-Actor + MCV | reference-free FP-selected Actor | 1 | 0.000684943 ± 0 | 2805.0 |

## Ambient-dimension scaling diagnosis

| d | params / train point | a1 selected step | a1 label-replication Gradient RMSE | selected Value RMSE | selected active / inactive Gradient RMSE |
|---:|---:|---:|---:|---:|---:|
| 100 | 254.7 | 250 | 0.015087 | 0.00150872 | 0.060729 / 0.0027515 |
| 400 | 367.3 | 1 | 0.017067 | 0.00629103 | 0.15318 / 0.010965 |

The d=400 Gradient-label replication noise is only 1.13× the d=100 value, while the selected error remains almost entirely in the two active directions and the inactive-output energy is negligible. Training after the first step nearly interpolates the training bank but worsens validation. The observed failure is therefore most consistent with full-ambient sample complexity/generalization under many irrelevant coordinates, rather than exponential-weight instability, MC-label noise alone, or spurious inactive Gradient output.

## Interpretation

Completed 7 training runs and 7 post-freeze evaluations. The rank correlation between the reference-free consistency score and post-freeze true Value error is 0.933794466403162.

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
