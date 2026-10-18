# Tangent-Actor High-Accuracy Frontier

## Decision: **GO**

This follow-up tests the maximum numerical accuracy of frozen Tangent-Actor models, a stronger raw Full-History frontier, a three-seed nonseparable problem, and matched-budget 400D training. Ground Truth is absent from Actor training and checkpoint selection throughout.

## Decisive findings

- The independent 512-step primary Tangent-Actor mean RMSE is `8.14161e-05` across three frozen training seeds.
- The strongest completed Full-History endpoint is `n=4, M=96` at `0.000219475` mean RMSE (95% t interval `[0.00020048, 0.00023847]`); the Actor value is lower.
- Tangent-Actor first reaches the frozen `5.0e-04` RMSE threshold at 128 Euler steps, with `3039.8` mean production GPU-seconds including training and selection.
- Full-History first reaches that threshold at `n=4, M=64` with `173.8` mean GPU-seconds.
- A primary accuracy advantage that clears MC and reference uncertainty is `True`.
- The three-seed high-budget nonseparable support criterion is `True`; every selected Actor improves over its own `a0`: `True`.
- Matched-budget full-ambient training fixes the 400D degradation across the required replications: `True` (3 completed seed(s)).

## Answers to the research questions

1. **128-step replication.** The independent three-seed mean is `0.000258083`; it beats the strongest completed Full-History mean: **False**, and clears the conservative uncertainty test: **False**.
2. **256 and 512 steps.** Mean RMSEs are `0.000138283` and `8.14161e-05`; their robust strongest-baseline comparisons are **True** and **True**.
3. **Full-History ceiling.** `STRONG_FULL_HISTORY_COMPARISON.md` retains every completed or failed depth/budget row.
4. **Strongest primary comparison.** Robust Tangent-Actor advantage: **True**.
5. **Nonseparable improvement.** Three-seed supporting result: **True**; lower mean error than the strongest coupled Full-History endpoint: **True**.
6. **Seed stability.** The primary 512-step seed standard deviation is `2.25e-06`; nonseparable trajectories are reported without selecting the best seed.
7. **400D matched budget.** Seed 0 crossed the predeclared 10% gate: **True**. Replicated success criterion met: **True**.
8. **Learned low rank.** The conditional branch was not triggered because matched-budget full-ambient 400D training crossed the 10% gate and replicated across three seeds. No oracle subspace is used, and no low-rank outcome is claimed.
9. **Remaining error.** The largest 512-step diagnostic component is **time_discretization_proxy**. The reported component magnitudes are `{'actor_approximation_proxy': 3.45957104580319e-05, 'finite_path_mc': 4.3253977118129776e-05, 'reference_uncertainty': 8.516656394984202e-06, 'time_discretization_proxy': 5.905697805414628e-05}`; the Actor term is a residual proxy after removing measured MC, reference, and CRN time-step contributions, rather than an exact decomposition.
10. **Defensible claim.** Tangent-Actor is an accurate learned linearization method for the tested smooth bounded-curvature, locally stable semilinear PDEs. Claims of universal nonconvex convergence, a global reference-free certificate, or general speed superiority are not supported.

## Artifacts

- `HIGH_PRECISION_EULER_AUDIT.md`
- `STRONG_FULL_HISTORY_COMPARISON.md`
- `NONSEPARABLE_HIGH_BUDGET_REPORT.md`
- `AMBIENT_400D_GENERALIZATION_AUDIT.md`
- `REFERENCE_FREE_ERROR_THEORY_NOTE.md`
- `accuracy_frontier.csv`, `euler_convergence.csv`, `full_history_extended_grid.csv`
- `nonseparable_three_seed_results.csv`, `ambient_scaling_comparison.csv`, `compute_accounting.csv`
- `reproducibility_manifest.json`
- `figures/primary_accuracy_frontier.{png,pdf}`
- `figures/nonseparable_high_budget.{png,pdf}`
- `figures/ambient_400d_generalization.{png,pdf}`
