# Round 3: equal-cost (Pareto) study, new PDE families, and tests of the one-step theory

You are working in the repository `Samyep/numerical_PDE`. Branch from `ir-mlp-candidate-screen`
(it contains rounds 1-2, the one-step theory note `docs/theory/`, and the candidate screen).
Create branch `ir-mlp-mechanism-suite-r3`. You may run long experiments and use the machine's CPUs.

## Hard rules

1. Do **not** modify the manuscript, `FullHistoryMLP`, any round-1/round-2 module, or any existing result.
   Put new code in `invariant_region_mlp/experiments/mechanism_suite_r3/` (subclass or wrap existing classes;
   e.g. generalise `_project_segment` in your own module, as `experiments/candidate_screen/screen.py` does).
   Results go to `results/mechanism_suite_r3/`, the report to `docs/MECHANISM_SUITE_R3_REPORT.md`.
2. **Freeze this document first.** Commit it verbatim as `experiments/mechanism_suite_r3/ROUND3_PREREGISTRATION.md`
   before writing implementation code or running anything. Record that commit hash in the report.
3. Fresh randomness: base seed `20261207`; new 1,200-point test sets with a fixed 20% validation split.
   Tree seeds `SeedSequence([20261207, d, n, M, rep, chunk_index])`; identical chunk size for all methods in a study.
4. Value accuracy decides every verdict (skill = RMSE/std(u*) on the test split). Gradients are reported only.
5. Verdicts are mechanical. Report every failure. No post-hoc thresholds. Exploratory analyses must be labelled as such.
6. float64 everywhere. Record work counters (generator calls, recursively evaluated states, stochastic samples, wall time)
   and non-finite counts for every row.

## 1. PDEs

All use `z = sigma * grad u`, `sigma = sqrt(2)`, `mu = 0`, test points `t ~ U[0,T)`, `x ~ U[-1,1]^d` unless stated.
Implementations of C1-C3 exist in `experiments/candidate_screen/candidates.py`; reuse them, but C3 needs a production reference (below).

| id | PDE | parameters | reference | certificate (all PDE-derived) | primary method |
|---|---|---|---|---|---|
| P1 | `u_t + Lap u - |grad u|^2 = 0`, `g = -log mean_k exp(lam_k w.x)` | lam=(0.5,1.5,3), T=0.25 | closed form | segment (Hopf-Cole hull), box hull, ball, span | box |
| P4 | `u_t + Lap u - |grad u| = 0`, `g = (log cosh(2 w.x) - log 2)/2` | T=0.5 | cached 1-D reference | tight comparison interval (round 2), round-1 segment, box, ball, span | tight_segment |
| C1 | l1-control HJB `u_t + Lap u - lam ||grad u||_1 = 0`, `lam = 1/||w||_1` | as P4 | P4 reference (identical ridge solution) | segment `|psi_s|<=1`, box, ball, span | segment |
| C2-convex | LQ game `u_t + Lap u - a|grad_A u|^2 + b|grad_B u|^2 = 0`, w mass 1/2 per block | a=2, b=0, d_A=d/2 | P1 closed form | P1 segment, box, ball, span | box |
| C2-cancel | same | a=3, b=1, d_A=round(d/4) | P1 closed form | same | box |
| C2-flip | same | a=3, b=1, d_A=round(d/8) | P1 closed form | same | box |
| C3 | cubic viscous HJ `u_t + Lap u - 0.5 |grad u|^3 = 0`, `g = 2 (log cosh(2 w.x) - log 2)/2` | T=0.25 | production 1-D reference (below) | segment `|psi_s|<=2`, box, ball, span | segment |
| LQG (exploratory) | `u_t + Lap u - |grad u|^2 = 0`, `g = gamma |x|^2` | gamma=1, T=0.5, x = r*theta with r ~ U[0,1], theta uniform on the sphere | `u = (d/2) log(1+4 gamma h) + gamma |x|^2/(1+4 gamma h)` | `z = 2 sqrt2 p x`, `p in [gamma/(1+4 gamma T), gamma]` (Riccati monotonicity; label "Riccati-derived") | segment |

Notes.
- C2 one-step prediction (docs/theory, Prop. 2 applied blockwise): orthogonal Jensen bias = `(c/M) K` with
  `K = 0.5*(-a(d_A-1/2) + b(d_B-1/2))`; negative for convex, approximately zero for cancel, positive for flip.
- C3 reference (gate G1 for C3): monotone Godunov scheme on `s in [-L, L]` with `L in {6, 8}` and
  `ds in {0.02, 0.01, 0.005, 0.0025}`; Richardson extrapolation of the first-order scheme. Pass if successive
  differences decrease, the finest Richardson difference is <= 1e-6 on `|s| <= 3`, the L=6 vs L=8 difference is <= 1e-7 there,
  and an independent finite-difference residual of the interpolated reference is <= 1e-5. Otherwise stop C3 and report.
- LQG prediction stated in advance: it may fail G2 because the d-scaled time-only term dominates std(u), and its raw
  gradient noise contains a `sigma^2 h |G|^2` term of size O(d h). Run E0 only; continue to S1 only if all gates pass.
  If you can implement antithetic terminal sampling as a subclass without touching `FullHistoryMLP`, add `raw_antithetic`
  as a control; otherwise state that it was not run.

## 2. Methods

raw; the primary certified method; the other valid certified geometries (segment/box/ball/span as listed);
`oracle_state`; `centre` (data-free certificate centre); `f_zero`; best `shrink_c`, c in {0.1, 0.25, 0.5, 0.75},
tuned on the validation split per cell; `oracle_z`, `oracle_u` only in S2.

## 3. Studies (run in this order)

**S0 Gates.** G1-G6 exactly as in round 1 for C1, C2 (all three), C3, LQG; carry P1/P4 gates from round 1.
Configurations `(3,6)` and `(4,6)`, 3 reps, d in {20, 100}. Write `gates.json`; stop any PDE that fails.

**S1 Equal-cost frontier (headline).** PDEs: P1, P4, C1, C2-convex, C2-flip, C3 (and LQG only if admitted).
d in {20, 100, 400}. Grid:
`n=1: M in {16, 64, 256, 1024}`, `n=2: M in {8, 16, 32, 64, 128, 256}`, `n=3: M in {4, 6, 10, 16, 24}`, `n=4: M in {3, 4, 6, 8}`.
10 reps per cell. Drop cells whose single-repetition wall time exceeds 30 min at d=400 and list them.
For each method, the frontier is the lower envelope of mean test skill against mean generator calls over all cells
(secondary axes: wall time, stochastic samples). Evaluate frontiers at 10 log-spaced cost levels spanning the range where
both frontiers exist (piecewise-constant envelope: best skill at or below the cost level). Bootstrap repetitions (1,000 draws) for 95% intervals.

**S2 Theory predictions on C2.** All three C2 configurations, d in {20, 100, 400}, cells `(2,8), (2,32), (3,6), (4,6)`, 10 reps.
Also extend `experiments/theory_onestep/verify_onestep_theory.py` to the C2 generator (terminal block, d in {10, 50, 200, 1000}, M in {4, 16}).

**S3 Dimension sweep.** P1, C1, C2-convex, C3 at d in {20, 50, 100, 200, 400}, cells `(2,32)` and `(3,10)`, 10 reps.

**S4 Certificate geometry.** C1 and C3: segment vs box vs ball vs span at cells `(2,32), (3,6), (4,6)`, d in {20, 100}, 10 reps (reuse S1 rows where identical).

## 4. Pre-registered criteria

- **PF-1 (frontier):** for each PDE in S1, at d=100 and at d=400, the primary method's frontier has skill <= 0.8 x the raw frontier
  at >= 8 of the 10 cost levels.
- **PF-2 (dimension robustness at equal cost):** for P1, C1, C2-convex, C3: at the median cost level, raw-frontier skill at d=400
  is >= 2 x its value at d=20, while the primary-method frontier at d=400 is <= 1.3 x its value at d=20.
- **PF-3 (data matters):** in the S1 cells, the primary method wins >= 7/10 paired reps against `centre` in >= 75% of cells,
  for P1, C1, C2-convex, C3. Prediction for P4 (tight certificate): fails PF-3 (round-2 replication).
- **PF-4 (beyond generic suppression):** same as PF-3 against the best tuned `shrink_c`.
- **T-1 (sign):** in S2 at cells (2,8) and (2,32), mean raw generator bias is < 0 for C2-convex and > 0 for C2-flip at every d,
  and |bias(C2-cancel)| <= 0.25 x min(|bias(convex)|, |bias(flip)|) at the same d and cell.
- **T-2 (where the curse appears):** best raw skill over the S2 cells grows by >= 2x from d=20 to d=400 for C2-convex and C2-flip,
  and by <= 1.3x for C2-cancel.
- **T-3 (one-step quantitative):** in the extended terminal-block check, measured Jensen gap / predicted `(c/M)K` lies in [0.9, 1.1]
  for every C2 configuration with |K| >= 5.
- **DIM (S3):** for P1, C1, C2-convex, C3, primary-method skill at d=400 / d=20 <= 1.3 at both cells.
- **LQG-G (exploratory, no claim):** report the gate outcome against the stated prediction.

## 5. Outputs

```
results/mechanism_suite_r3/
  gates.json  s1_frontier.csv  s1_frontier_levels.csv  s2_theory.csv  s2_onestep_check.json
  s3_dimension.csv  s4_geometry.csv  tuning_choices.json  analysis_summary.json  figures/
docs/MECHANISM_SUITE_R3_REPORT.md
```
Figures: one frontier panel per PDE and d (raw, primary, centre, best shrink, f_zero; log cost axis);
C2 bias-sign panel; dimension panel. The report starts with an outcome table of all criteria, then lists each criterion verbatim
with its verdict, then all failures, then provenance (prereg commit, code commit, environment, work totals).

Priority if compute is short: S0 -> S2 -> S1 (d=20, 100) -> S3 -> S1 (d=400) -> S4. Never drop a pre-registered criterion
silently; mark it "not evaluated (compute)" instead.
