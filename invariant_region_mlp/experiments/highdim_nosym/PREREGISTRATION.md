# High-effective-dimension, symmetry-free benchmark: Tangent-Actor vs Full-History MLP on S05 and B6 (pre-registered)

Repository `Samyep/numerical_PDE`, subproject `invariant_region_mlp`. Create branch `ir-mlp-highdim-nosym` from the
current head of `ir-mlp-tangent-actor-icml-final-suite` (`72d7c09a` or later). Bring in the validated S-family code from
branch `ir-mlp-nonridge` (commit `d384f4e8`, directory `experiments/nonridge/`) with `git checkout d384f4e8 -- <paths>`
and verify byte identity. Do not modify the final suite's frozen files, its `FROZEN_PROTOCOL.json`, or
`REVIEW_GATE_STATUS.json`; do not touch any other running branch.
You are authorized to run long experiments locally and use all available CPU/GPU resources reasonably,
but do not delete existing results or rewrite the manuscript.

## Why this study

Every PDE in the final suite has low intrinsic complexity: k=1 or k=2 ridges (joint, nonconvex, MLP failure), radial
symmetry (HJB-Log, Allen–Cahn), or an isotropic quadratic plus a 5-direction ridge (OU/GMM). Three open questions are
answered only by a problem without these structures:

1. **Non-triviality of learning.** In the final suite, the untrained terminal tangent `a0` with MCV inference is within
   2–4x of DPI/PINN-HTE, and almost all of the Actor's gain comes from the first update. Is the learned Actor still
   decisive, and do later updates matter, when the solution has no low-dimensional structure?
2. **Tangent-Actor vs Full-History MLP in the regime where Full-History MLP is supposed to fail.** On the final-suite
   PDEs, raw Full-History MLP is cheap and accurate. On S05, raw MLP at `n=3` is 10^2–10^4 times worse than its
   pathwise variant (`ir-mlp-nonridge`, criterion N-1).
3. **Does the method depend on symmetry or a fixed active subspace?**

## Hard rules

1. New code goes in `invariant_region_mlp/experiments/highdim_nosym/`, results in `results/highdim_nosym/`, and the
   report in `docs/HIGHDIM_NOSYM_REPORT.md`. Import the mature Tangent-Actor, MCV, labels, and selector code and the
   GPU Full-History solver; do not edit them. Add new equations through adapters in the new directory.
2. Freeze this document first: commit it verbatim as `experiments/highdim_nosym/PREREGISTRATION.md` before any
   implementation or run, and record the hash in the report.
3. Analytic or quadrature truth never enters Actor labels, fitting, early stopping, iteration choice, checkpoint
   selection, or the B6 parameter choice. The B6 calibration below is reference-only and must be frozen before any
   solver touches B6.
4. Tangent-Actor uses the **frozen mature configuration** of the final suite: terminal-aware Gradient Actor
   `q = sqrt(2) grad g(x) + (T - t) Q_theta(t, x)` (no Value head is needed because `f` does not depend on `u`), 8,192
   training states, 128 paths per label, 32 label steps, 4x384 SiLU, 5,000 optimizer steps per update, four updates,
   MCV training and inference, and the final suite's reference-free selector. No new hyperparameter search. If a run
   diverges or is non-finite, apply one predeclared fallback (halve the learning rate and restart that update from the
   previous Actor), record it, and report it as a deviation.
5. Full-History MLP grids are frozen below. Report every cell, including failures. Never select a Full-History cell
   with truth for a headline; report the full grid and the descriptive best.
6. Verdicts are mechanical; report every failure; no post-hoc thresholds; exploratory analyses are labelled.

## PDEs

Common equation and convention (as in P1/S-family): `u_t + Laplacian(u) + f(z) = 0`, `z = sqrt(2) grad u`,
`f(z) = -|z|^2 / 2`, i.e. `u_t + Laplacian(u) - |grad u|^2 = 0`, `sigma = sqrt(2)`, `mu = 0`, terminal `u(T, x) = g(x)`.
Hopf–Cole gives `u = -log w`, `w(t, x) = E[exp(-g(x + sqrt(2 (T - t)) G))]`.
Dimensions `d in {100, 400}`.

### S05 (already implemented and validated on `ir-mlp-nonridge`)

- `phi(y) = log 3 - log(exp(0.5 y) + exp(1.5 y) + exp(3 y))`;
  `Q x = DCT(s * x)` (orthonormal DCT-II, `s = 2 * default_rng(SeedSequence([20261201, d, 31])).integers(0, 2, d) - 1`);
  `c_j = j^(-1/2) / sqrt(sum_i i^(-1))`; `g(x) = sum_j c_j phi((Q x)_j)`; `T = 0.25`.
- Exact solution: `u(t, x) = sum_j v(c_j, T - t, (Q x)_j)`, `v(c, tau, y) = -log E[exp(-c phi(y + sqrt(2 tau) xi))]`,
  `z = sqrt(2) Q^T (dv/dy)_j`; Gauss–Hermite with 96 nodes (validated to 1e-15; PDE residual 1e-7).
- Test bank: the frozen `ir-mlp-nonridge` bank (`fixed_points` with `pde_id = "S05"`: `t ~ U[0, T)`,
  `x ~ U[-1, 1]^d`, 1,200 states, the same 960-state test split), so that existing Full-History rows pair exactly.
  Verify the bank hash against `ir-mlp-nonridge`.
- Known calibration (`ir-mlp-nonridge`): Hessian effective rank 16 (d=100) / 27 (d=400), centred `z*` variance outside
  the top-10 directions 0.41 / 0.52, heat-only skill 0.41 / 0.32.

### B6 (new): full-rank, non-separable log-sum-exp HJB

- `g(x) = -log[(1/K) sum_{m=1..K} exp(b_m . x)]`, with `K = 2d` and
  `b_m = gamma * v + beta_m * e_m`. Here `e_m` are i.i.d. uniform on the unit sphere,
  `beta_m = beta * U_m` with `U_m ~ U[0.5, 1.5]`, and `v` is a fixed unit vector. All are drawn from
  `default_rng(SeedSequence([20261201, d, 66]))` in that order (`e_m`, then `U_m`, then `v`). `T = 0.25`.
- Exact solution: `u(t, x) = -log[(1/K) sum_m exp(b_m . x + (T - t) |b_m|^2)]`,
  `z*(t, x) = -sqrt(2) sum_m pi_m b_m`, `pi_m` the softmax of those logits,
  Hessian `grad^2 u = -Cov_pi(b_m)` (full rank whenever `{b_m}` affinely span `R^d`).
- Test bank: same law as S05 (`t ~ U[0, T)`, `x ~ U[-1, 1]^d`, 1,200 states, 20% validation split), built with the
  same `fixed_points` construction and `pde_id = "B6"`. Training and validation states for the Actor use that law,
  as for S05.
- Unequal `|b_m|` are required: with equal norms the time term is a common shift and `z*` would not depend on `t`.

**B6 calibration (reference-only, before any solver run; separately for each d).**
Grid: `beta in {1.0, 1.5, 2.0, 3.0, 4.0}` x `gamma in {0, 1}` (10 configurations). On the 960 test states compute:
- `ESS_med`: median over states of `exp(entropy(pi))`, the effective number of active components;
- `z_rms`: RMS `|z*|`;
- `f_cv`: spatial standard deviation of `f(z*)` divided by `|mean f(z*)|`;
- heat-only skill: RMSE between `u*` and the `f = 0` solution `E[g(x + sqrt(2 tau) G)]`, divided by `std(u*)`, using
  4,096 Monte Carlo samples per state;
- `r_glob`: effective rank `(tr C)^2 / tr(C^2)` of the centred covariance `C` of `z*` over the test states, and the
  fraction of `tr C` outside its top-10 eigen-directions;
- `r_loc`: median local Hessian effective rank of `Cov_pi(b_m)`.

**Selection rule (frozen now):** among configurations with `ESS_med >= 5`, `0.5 <= z_rms <= 4`, heat-only
skill `>= 0.15`, and `f_cv >= 0.05`, choose the one with the largest `r_glob`; break ties by the smaller `beta`. If none
qualifies, choose the configuration that minimizes the sum of normalized constraint violations, label B6
"calibration-limited", and continue.

Validate B6 before any solver run:
- terminal identity;
- autodiff PDE residual at most 1e-10 relative on 200 states;
- `z*` against finite differences of `u*`;
- Hopf–Cole Monte Carlo with 1e6 samples at 20 states, within 4 standard errors;
- GPU equation adapters equal numpy to 1e-12;
- a draw-tape replay of the GPU Full-History solver for `path` and `double_path_base` at d=100, cells (2,8) and (3,4),
  agreeing to 1e-10.

Report the calibration table in full.

## Methods

| ID | Method | Notes |
|---|---|---|
| `TA` | Tangent-Actor, mature frozen config | selected checkpoint by the final suite's reference-free selector |
| `TA_a0` | untrained terminal tangent `q = sqrt(2) grad g` | same MCV inference as `TA`; zero training |
| `TA_ai` | post-freeze per-iteration diagnostics `a0 ... a4` | common Brownian numbers; descriptive only, never used to select |
| `FH_raw` | Full-History MLP, Bismut terminal and Bismut levels | GPU solver, `raw` |
| `FH_path` | Full-History MLP with pathwise terminal gradient | GPU solver, `path` |
| `FH_rep` | repaired Full-History MLP: pathwise + level-1 baseline + product double | GPU solver, `double_path_base` (exact for quadratic `f`) |

Full-History frozen grid for every PDE and d:
- `raw` and `path`: `(n, M) in {(2,16), (2,64), (2,256), (3,8), (3,16), (3,24), (4,6), (4,10)}`;
- `double_path_base`: `(n, M) in {(3,16), (3,24), (3,32), (4,6), (4,10), (4,16)}`;
- 3 estimator seeds per cell (the S05 rows from `ir-mlp-nonridge` with 10 repetitions may be reused only where the
  method, cell, d, bank and code revision match exactly; record the reuse).

Tangent-Actor inference:
- primary protocol: 64 Euler steps, 256 paths, 4 independent replicates;
- post-freeze audit: 128 and 256 steps, labelled numerical audits.

Training seeds:
- seed 0 for all four (PDE, d) cells;
- seeds 1 and 2 at d=100 for both PDEs if the budget allows.

## Metrics and accounting

For every row, on the 960 test states:
- value RMSE, relative L2, `RMSE / std(u*)`, bias, gradient-vector RMSE where available;
- MC standard-error contribution, separated from the non-MC proxy (as in the final suite);
- nonfinite counts.

Time: synchronized GH200 wall time.
- Tangent-Actor: label generation, fitting, reference-free selection and final inference are reported separately and
  as a total.
- Full-History: inference only.

Also report the descriptive accuracy-cost envelope `E*(C) = min_{i: C_i <= C} E_i` for the `TA` family (iterations,
inference budgets) and for each Full-History family, and the GPU seconds each family needs to first reach
`RMSE <= 2 x RMSE(TA, seed 0)` (or "not reached" within its grid).

## Pre-registered criteria (evaluated on seed 0 at the primary 64-step protocol)

- **H1 (learning is decisive):** `RMSE(TA) <= 0.5 x RMSE(TA_a0)` in all 4 (PDE, d) cells.
- **H2 (more than one update matters):** `RMSE(TA) <= 0.8 x RMSE(a1)`, both under the same post-freeze common-random-number
  diagnostic, in at least 2 of the 4 cells. A failure means "one policy-improvement step suffices here" and must be
  stated plainly.
- **H3 (beats raw Full-History MLP):** `RMSE(TA) <= 0.5 x` the minimum RMSE over the entire frozen `FH_raw` grid, in
  all 4 cells.
- **H4 (beats every Full-History variant):** `RMSE(TA) <=` the minimum RMSE over all frozen `FH_raw`, `FH_path` and
  `FH_rep` cells (any cost), in at least 3 of the 4 cells.
- **H5 (efficiency vs the best Full-History variant):** in at least 2 of the 4 cells, the best Full-History family
  does not reach `2 x RMSE(TA)` within its frozen grid, or needs more GPU seconds than `TA`'s total to do so.
- **H6 (no dependence on a fixed subspace):** for B6, `r_glob >= d/4` at the frozen calibration, and H1 holds on B6 at
  both d. Report `r_glob` for S05 alongside.
- **Reported, no verdict:** seeds 1–2 at d=100; 128/256-step audits; MCV variance-reduction factors (never as RMSE
  factors); the reference-free selector choice against the post-freeze per-iteration truth; peak memory.

Predictions (no verdict): H3 passes easily (`FH_raw` collapses at `n >= 3`); H4 and H5 are the informative criteria;
H2 may fail on S05 (as on the joint PDEs) and pass on B6 if the nonlinearity is stronger there.

## Optional Stage C (only after all of the above, labelled exploratory)

1. **Noise-replicate test of the residual Jensen term:**
   - setting: S05 at d=100;
   - from the frozen `a3`, run one update with 4 independent label-noise seeds;
   - keep states, initialization, optimizer seed and minibatch order identical;
   - report `J = mean_k v(A_k) - v(mean_k A_k)` on common Brownian paths, and `J_pred = 0.5 E int mean_k |q_k - q_bar|^2 ds`.
2. **One neural baseline:** DPI with its official configuration transferred to S05 and B6 at d=100 (its full published
   budget), on the same banks.

## Outputs

```
results/highdim_nosym/  validation.csv  b6_calibration.csv  main_rows.csv  fh_grid.csv  iteration_diagnostics.csv
                        frontier.csv  time_to_threshold.csv  resource_accounting.csv  analysis_summary.json  figures/
docs/HIGHDIM_NOSYM_REPORT.md
```

Figures:
- RMSE against GPU seconds per (PDE, d), with `TA`, `TA_a0`, the per-iteration points and all three Full-History
  families;
- RMSE by Actor iteration;
- the B6 calibration table as a heat map.

The report starts with an outcome table H1–H6, then each criterion verbatim with its verdict, the calibration, failures
and deviations, and provenance (preregistration commit, code commits, GPU and software versions, bank hashes,
validation results, reuse audit).

## Priority if compute is short

1. Validation and B6 calibration.
2. S05 d=100: `TA`, `TA_a0`, then the Full-History grid.
3. B6 d=100: same order.
4. S05 d=400.
5. B6 d=400.
6. Per-iteration diagnostics.
7. Seeds 1–2.
8. Stage C.

If a single Full-History cell exceeds 30 minutes per seed, run one seed and label it. Never drop a criterion silently;
mark it "not evaluated (compute)".
