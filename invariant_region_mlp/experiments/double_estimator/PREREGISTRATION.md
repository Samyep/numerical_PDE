# Double estimator for the generator: pre-registered study

You are working in the repository `Samyep/numerical_PDE`, subproject `invariant_region_mlp`.
Branch from `ir-mlp-expert-iteration` (commit `2db332ea` or later) and create branch `ir-mlp-double-estimator`.
You are authorized to run long experiments locally and use all available CPU/GPU resources reasonably,
but do not delete existing results or rewrite the manuscript.

## Hard rules

1. Do **not** modify the manuscript, `FullHistoryMLP`, any existing module, or any existing result.
   New code goes in `invariant_region_mlp/experiments/double_estimator/`, results in `results/double_estimator/`,
   the report in `docs/DOUBLE_ESTIMATOR_REPORT.md`. Subclass or wrap existing classes (see
   `experiments/effective_dim_screen/cv_probe.py`, which already overrides `solve` and `_terminal_estimate`).
2. **Freeze this document first**: commit it verbatim as `experiments/double_estimator/PREREGISTRATION.md` before any
   implementation or run, and record the commit hash in the report.
3. Base seed `20261201`; paired seeds `SeedSequence([20261201, d, n, M, rep, chunk_index])`, identical chunk size for all
   methods in a study, so all methods see the same random numbers wherever they draw the same variables.
4. float64. Value accuracy (skill = RMSE/std(u*)) decides every verdict. Verdicts are mechanical; report every failure;
   no post-hoc thresholds; exploratory analyses labelled as such.
5. Record per row: skill, mean generator bias, generator RMSE, generator calls, recursive calls, wall time, non-finite counts.

## Background

Inside MLP, the generator f receives a Monte Carlo gradient estimate. For a curved f this produces a bias
(approximately 0.5 tr(Hess f Cov z_hat), growing like (d-1)c/M), the PDE analogue of the maximization bias of Q-learning.
The pathwise terminal gradient removes the noise only in the terminal block; correction-level integrands are Monte Carlo
outputs and cannot be differentiated, so deep cells still diverge. Double Q-learning removes maximization bias by using
two independent estimates. The analogue here: evaluate f with two **independent** recursive estimates of the same state.

## Methods

At every generator call in the correction levels (Algorithm: `D_j = f(U_l) - f(U_{l-1})` at the child point
`(S_j, X_j)`), compute **two independent recursive calls** `U_l^(1)`, `U_l^(2)` (fresh randomness, same child point), and
likewise `U_{l-1}^(1)`, `U_{l-1}^(2)`, and replace `f(U)` by `f_double(U^(1), U^(2))`:

- quadratic generators `f(z) = -(lam/2)|z|^2` (P1, multi-ridge): `f_double = -(lam/2) z^(1) . z^(2)`
  (exactly unbiased for -(lam/2)|E z_hat|^2 given independence);
- block game `f(z) = -(a/2)|z_A|^2 + (b/2)|z_B|^2`: `f_double = -(a/2) z_A^(1).z_A^(2) + (b/2) z_B^(1).z_B^(2)`;
- norm generator `f(z) = -lam |z|` (P4; Double-Q form, since |z| = max_{|a|<=1} a.z):
  `f_double = -lam (z^(1)/|z^(1)|) . z^(2)` (one-sided bias, bounded), with `z^(1)=0` mapped to 0.

The two estimates are used **only** inside f. The value and gradient contributions of the level otherwise use the
average `(U^(1)+U^(2))/2` wherever the base code would use U (state this choice in the report).
`double` therefore costs about twice the recursive work of `raw` at the same (n, M); all comparisons at equal cost use
measured generator calls and wall time.

Methods per PDE: `raw`; `path` (pathwise terminal gradient, as in `cv_probe.py`); `double`; `double_path`;
`oracle_state`; `centre` (data-free certificate centre); `f_zero`; and the certified `box` (P1, C2, multi-ridge
`sub_box`) as the certificate-based reference.

## PDEs and grid

| id | PDE | reference | d |
|---|---|---|---|
| P1 | ridge LSE HJB (round-1 parameters) | closed form | 20, 100, 400 |
| MR | multi-ridge LSE HJB, k=10, scale=4, T=0.1 (`effdim_mlp.py`) | closed form | 100, 400 |
| C2 | game, convex / cancel / flip (round-3 parameters) | P1 closed form | 20, 100, 400 |
| P4 | norm HJB (round-1 parameters, cached 1-D reference) | 1-D reference | 20, 100 |

Cells `(n, M)` in `{(2,32), (3,6), (3,10), (4,6)}`; 10 paired repetitions; 1,200 test points with a fixed 20% validation
split (as in round 3). Equal-cost study (S-EC) for P1 and MR at d=100 and d=400: the round-3 S1 grid
(`n=2: M in {8,16,32,64,128}`, `n=3: M in {4,6,10,16}`, `n=4: M in {3,4,6}`), frontier = lower envelope of mean skill versus
mean generator calls, evaluated at 10 log-spaced cost levels where both frontiers exist; 1,000-draw bootstrap over reps.

## Pre-registered predictions and criteria

- **D-1 (bias removed):** for P1 and MR, in cells (3,6) and (4,6), |mean generator bias| of `double` <= 0.05 x that of
  `raw`, at every d.
- **D-2 (deep cells repaired):** for P1 and MR, `double` skill <= 0.1 x `raw` skill in cells (3,6) and (4,6) at every d,
  and `double` skill at the largest d <= 1.5 x its value at the smallest d in each of those cells.
- **D-3 (complements pathwise):** `double_path` skill <= `path` skill in cells (3,6) and (4,6) for at least 75% of the
  (PDE, d, cell) combinations over P1 and MR.
- **D-4 (equal cost):** for P1 and MR at d=100 and d=400, the frontier of the better of `double`/`double_path` is <= 0.8 x
  the `raw` frontier at >= 8 of 10 cost levels.
- **D-5 (sign-independent):** on C2, |bias(`double`)| <= 0.1 x |bias(`raw`)| for convex and flip, at cells (2,32), (3,6),
  every d.
- **D-6 (exploratory, P4):** report bias and skill of the Double-Q norm form; no verdict.
- **Stated risk (prediction, not a criterion):** `double` removes the bias of f but not the variance of the Bismut-weighted
  level terms; in the deepest cells the variance (roughly d c^2 / M^2 per call) may still dominate. Report generator RMSE and
  across-repetition spread so that bias and variance can be separated.

## Outputs

```
results/double_estimator/
  rows.csv  summary.csv  frontier.csv  frontier_levels.csv  analysis_summary.json  figures/
docs/DOUBLE_ESTIMATOR_REPORT.md
```
Figures: skill vs d per cell (raw, path, double, double_path, box, oracle); bias vs d; equal-cost frontiers.
The report starts with an outcome table of all criteria, then each criterion verbatim with its verdict, then failures and
deviations, then provenance (prereg commit, code commit, environment, totals).
Priority if compute is short: P1 and MR core grid -> C2 -> S-EC -> P4. Never drop a criterion silently; mark it
"not evaluated (compute)".
