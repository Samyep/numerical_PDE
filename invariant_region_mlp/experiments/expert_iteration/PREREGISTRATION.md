# Noise-aware defect correction and expert iteration: pre-registered study

You are working in the repository `Samyep/numerical_PDE`, subproject `invariant_region_mlp`.
Branch from `ir-mlp-candidate-screen` (latest commit includes `experiments/effective_dim_screen/defect_probe.py`,
which is the numpy prototype of the defect-corrected MLP used below). Create branch `ir-mlp-expert-iteration`.
You are authorized to run long experiments locally and use all available CPU/GPU resources reasonably,
but do not delete existing results or rewrite the manuscript.

## Hard rules

1. Do **not** modify the manuscript, `FullHistoryMLP`, any existing module, or any existing result.
   New code goes in `invariant_region_mlp/experiments/expert_iteration/`, results in `results/expert_iteration/`,
   the report in `docs/EXPERT_ITERATION_REPORT.md`. Subclass or wrap existing classes (as `defect_probe.py` does).
2. **Freeze this document first.** Commit it verbatim as `experiments/expert_iteration/PREREGISTRATION.md`
   before writing implementation code or running anything. Record that commit hash in the report.
3. Seeds: base seed `20261101`. Every stochastic component (network init, collocation points, MC samples, test points)
   gets its own `SeedSequence([20261101, study_id, problem_id, d, seed, component_id, ...])`. Record the scheme.
4. Verdicts are mechanical. Report every failure. No post-hoc thresholds. Exploratory analyses are labelled as such.
   Hyperparameter tuning is allowed only where this document says so, only on the validation split, and with an
   equal tuning budget for every method that is tuned. Record every tuning choice in `tuning_choices.json`.
5. Monte Carlo recursion and references in float64. Networks may train in float32, but every network evaluation inside
   the MLP recursion and every reported error uses float64 (cast the trained weights).
6. Record for every row: wall time, device, network forward/backward counts, generator calls, terminal samples,
   non-finite counts.

## Background (why these experiments)

MLP estimates z = sigma grad u with Bismut-Elworthy-Li (Malliavin) weights. Each estimate carries isotropic noise in
all d directions; plugging it into a nonlinear generator f gives a bias of about 0.5 tr(Hess_z f * Cov z_hat), growing
like (d-1) c / M. In SCaSML-style defect correction (MLP on e = u - u_theta), the constant c scales with the surrogate's
gradient error, so the correction helps only when the surrogate is already good. Replacing the terminal Bismut weight
by the pathwise gradient sigma grad(g - u_theta(T, .))(X) (possible because g and the network are differentiable)
removes the d-dependence in controlled tests (`results/effective_dim_screen/defect_probe_summary.md`).
Study A tests this on SCaSML's own LQG benchmark with real networks. Study B tests whether the corrected solver
can be used as a teacher (expert iteration) and whether that beats Deep Picard Iteration (DPI) at equal cost.

## Common definitions

- Generator convention: `u_t + Lap u + f(u, z) = 0`, `z = sigma grad u`, `sigma = sqrt(2)`, `mu = 0`.
- Defect equation for a surrogate u_th (z_th = sigma grad u_th by autodiff):
  `e_t + Lap e + f_e = 0`, `e(T) = g - u_th(T)`,
  `f_e(t, x, e, z_e) = f(u_th + e, z_th + z_e) - f(u_th, z_th) + r_th`,
  `r_th = d_t u_th + Lap u_th + f(u_th, z_th)` (surrogate residual, autodiff).
  Laplacian: exact for d <= 50; for d > 50 use Hutchinson with Rademacher probes (number of probes = ceil(d/4),
  as in SCaSML). r_th enters f_e linearly, so Hutchinson noise does not enter the nonlinearity; additionally run the
  exact Laplacian at d = 100 on one seed to confirm (report the difference).
- Prediction after correction: `u_th + e_hat`.
- The MLP recursion is the existing full-history MLP (Beta(0.5) time sampling, as in `mechanism_mlp.py`).
  Gradient estimators:
  - `bismut`: unchanged base code (terminal and level terms with Bismut weights).
  - `path`: terminal gradient replaced by `mean_k sigma grad(g - u_th(T, .))(X_k)`; level terms unchanged (Bismut).
- Error metrics on test points: relative L2 of u (primary, as in SCaSML), skill = RMSE / std(u*), relative L2 of z (reported only).

---

## Study A: SCaSML's LQG benchmark: reproduce the failure, then fix it

**PDE.** `u_t + Lap u - lam |grad u|^2 = 0` (so `f(z) = -(lam/2)|z|^2`), terminal
`g(x) = log( (1 + sum_{i=1}^{d-1} [c1_i (x_i - x_{i+1})^2 + c2_i x_{i+1}^2]) / 2 )`, `c1_i, c2_i ~ U[0.5, 1.5]` (fixed seed).
SCaSML leaves some details unspecified; our reading (record it as such in the report): `lam = 1`, no drift (the HJB as
written has none), `T = 0.5`, test points `t ~ U[0, 0.5)`, `x` uniform in the unit ball of R^d, 1,000 test points + 200
validation points.

**Reference.** Hopf-Cole: `u = -(1/lam) log E exp(-lam g(x + sqrt2 W_{T-t}))`,
`grad u = E[exp(-lam g(X)) grad g(X)] / E[exp(-lam g(X))]`. Monte Carlo with antithetic pairs; increase samples until the
estimated standard error of u is <= 1e-4 at every test point (report max SE). This is gate A-G1.

**Gate A-G0 (how nonlinear is this benchmark?).** Report the relative L2 error of the f = 0 solution `E g(x + sqrt2 W)`
(same MC). Prediction: small. Whatever it is, include `f_zero` as a baseline in every table below.

**Surrogates.** PINN as in SCaSML: 5 hidden layers x 50, tanh, Adam (lr 1e-3, betas (0.9, 0.99)), 2,500 iterations,
100 interior + 1,000 terminal points per iteration, Hutchinson d/4 Laplacian in training, unit loss weights.
Save checkpoints at 500, 1,000 and 2,500 iterations (three surrogate quality levels). d in {20, 50, 100, 120, 140, 160};
3 network seeds.

**Methods** (all on the same surrogate and test points):

| name | description |
|---|---|
| `surrogate` | u_th alone |
| `f_zero` | f = 0 solution (no surrogate) |
| `mlp` / `mlp_clip` | plain MLP (bismut), without / with SCaSML's clipping (threshold 10 on u and each z coordinate, per level) |
| `scasml` | defect MLP, bismut, clipping 0.1 (faithful SCaSML) |
| `scasml_noclip` | defect MLP, bismut, no clipping |
| `path` | defect MLP, pathwise terminal gradient, no clipping (**primary**) |
| `path_clip` | `path` + clipping 0.1 |
| `oracle_state` | defect MLP with exact (e, z_e) fed to f_e (lower bound) |

MLP cells `(n, M)` in {(2, 10), (2, 32), (3, 6)}; primary cell (2, 10) (SCaSML's VB setting). 5 MC repetitions per cell.

**Pre-registered criteria** (final 2,500-iteration surrogate unless stated; median over network seeds):
- **A-1 (reproduction):** `mlp` relative L2 >= 1 at d >= 100, and `scasml` / `surrogate` error ratio >= 0.5 at d = 160 in the
  primary cell. (If not reproduced, report and continue.)
- **A-2 (fix):** in the primary cell, `path` <= 0.2 x `surrogate` and `path` <= 0.5 x `scasml` at every d in {100, ..., 160}.
- **A-3 (dimension):** the ratio `path / surrogate` at d = 160 is <= 1.3 x its value at d = 20, while the ratio
  `scasml_noclip / surrogate` at d = 160 is >= 1.5 x its value at d = 20.
- **A-4 (mechanism):** mean generator bias (mean of f_e(estimate) - f_e(truth) over all generator calls) for
  `scasml_noclip` is negative at every d and its magnitude at d = 160 is >= 3 x that at d = 20; for `path` the magnitude at
  d = 160 is <= 1.5 x that at d = 20.
- **A-5 (surrogate quality, exploratory):** plot `path` and `scasml_noclip` error against surrogate error across the three
  checkpoints; report whether `scasml_noclip` is ever worse than the surrogate alone.

---

## Study B: expert iteration (EI) versus Deep Picard Iteration (DPI)

**Problems.**
- B1: Study A's LQG, d = 100.
- B2: multi-ridge log-sum-exp HJB, `k = 10, d = 100, scale = 4, T = 0.1` (`MultiRidgeHJB` in `experiments/effective_dim_screen/effdim_mlp.py`; closed form).
- B3: multi-ridge, `k = 10, d = 100, scale = 2, T = 0.5` (longer horizon; DPI's advantage was reported to grow with T).
- B4 (prediction test, lower priority): multi-ridge, `k = 50, d = 400, scale = 5, T = 0.1`.
Test points as in the respective source (multi-ridge: `t ~ U[0,T)`, `x ~ U[-1,1]^d`); 1,000 test + 200 validation.

**Network (identical for all methods).** Scalar u_theta(t, x): 4 hidden layers x 128, tanh; z = sigma grad u_theta by autodiff.
Loss for every regression step: `mean |u - y|^2 + w_z mean |sigma grad u - z_label|^2`, `w_z` tuned (see below).
Every method warm-starts from the previous round's weights.

**Round 0 (common to all methods).** Fit the network to f = 0 Feynman-Kac labels (value and pathwise gradient of
`E g(x + sqrt2 W)`) at N label points. All methods start from this same checkpoint.

**Methods.**
- `EI-path` (**primary**): each round, draw N fresh label points from the training distribution; labels
  `(u_th + e_hat, z_th + z_e_hat)` from defect MLP with the `path` estimator, cell (2, 10) unless tuned; refit for E epochs.
- `EI-bismut`: identical but `bismut` estimator (prediction: fails to contract on B4).
- `DPI`: Han, Hu, Long, Zhao (arXiv 2409.08526). Each round, value labels
  `y = mean_m [ g(X_T) + (T - t) f(u_k, sigma grad u_k)(s_m, X_{s_m}) ]` with `s_m ~ U[t, T]` and f evaluated on the current
  network; gradient labels with their Bismut-Elworthy-Li estimator including their control variates (implement from their
  Eqs. 19-20; cite equation numbers in code). Refit for E epochs.
- For every method, the final network is also evaluated with inference-time correction (`path`, cell (2, 32)).

**Budget and tuning.** Primary cost axis: wall-clock on the same device, including label generation and training.
Secondary: network evaluations + generator calls. Tuning (validation split, B2 only, equal budget of 12 configurations
per method): learning rate, epochs per round E, labels per round N, w_z, and MC samples per label (DPI: M_DPI; EI: MLP cell
among (2, 10), (2, 32), (3, 6)). Freeze the chosen configuration and use it unchanged on B1, B3, B4.
Run EI for R = 6 rounds. Run DPI for as many rounds as fit in 4x the EI total wall-clock.
3 seeds (network init + sampling).

**Pilot first.** One seed on B2 with the tuned configuration. Go/no-go (pre-registered): proceed to the full study only if
`EI-path` test error after round 3 is <= 0.5 x after round 0. If not, stop Study B, write the report, and include the
pilot curves. Do not redesign the method after seeing pilot results.

**Pre-registered criteria** (median over seeds; "error" = test relative L2 of the network alone unless stated):
- **B-1 (contraction):** `EI-path` error after round 3 <= 0.2 x its error after round 0, on B1, B2, B3.
- **B-2 (divergence prediction):** on B4, `EI-bismut` error after round 3 >= its error after round 1, while `EI-path`
  satisfies B-1 on B4.
- **B-3 (headline, equal cost):** at the wall-clock of `EI-path` round 6, `EI-path` error <= 0.5 x the best DPI error
  reached within the same wall-clock, on at least 2 of B1-B3.
- **B-4 (DPI given more compute):** report whether DPI at 4x the wall-clock reaches `EI-path`'s round-6 error (no verdict).
- **B-5 (inference-time correction):** for the final `EI-path` network, corrected error <= 0.3 x network error on B1-B3.

---

## Outputs

```
results/expert_iteration/
  A_reference.json  A_gates.json  A_rows.csv  A_summary.csv
  B_tuning.csv  tuning_choices.json  B_pilot.csv  B_rounds.csv  B_summary.csv
  analysis_summary.json  figures/
docs/EXPERT_ITERATION_REPORT.md
```
Figures: (A) error vs d for all methods, primary cell; generator bias vs d; error vs surrogate error (A-5).
(B) error vs wall-clock for EI-path, EI-bismut, DPI per problem (log-log), with round markers.
The report starts with an outcome table of all criteria and gates, then each criterion verbatim with its verdict,
then all failures and deviations (including every ambiguity in SCaSML/DPI that you had to resolve), then provenance
(prereg commit, code commit, environment, device, totals).

Priority if compute is short: A-G0/G1 -> Study A at d in {20, 100, 160} -> B pilot -> rest of Study A -> Study B full
(B2, B1, B3) -> B4. Never drop a criterion silently; mark it "not evaluated (compute)".
