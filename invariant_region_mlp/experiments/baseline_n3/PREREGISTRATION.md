# Correction-level baseline at n = 3, and a GPU port: pre-registered study

Repository `Samyep/numerical_PDE`, subproject `invariant_region_mlp`. Branch from `ir-mlp-path-frontier`
(commit `a42e2c5c` or later) and create branch `ir-mlp-baseline-gpu`.
You are authorized to run long experiments locally and use all available CPU/GPU resources reasonably,
but do not delete existing results or rewrite the manuscript.

The work has two studies. **Study F** (CPU, scientific) comes first. **Study G** (GPU port, engineering) may be
implemented in parallel, but its exploratory runs (G-5) start only after its validation criteria pass.

## Hard rules

1. Do not modify any existing module or result. New code goes in `invariant_region_mlp/experiments/baseline_n3/`
   (Study F) and `invariant_region_mlp/experiments/gpu_mlp/` (Study G). Study F must **import or subclass** the frozen
   `experiments/double_estimator/solver.py` (`DoubleEstimatorMLP`) unchanged. Results go in `results/baseline_n3/` and
   `results/gpu_mlp/`; reports in `docs/BASELINE_N3_REPORT.md` and `docs/GPU_MLP_REPORT.md`.
2. Freeze this document first: commit it verbatim as `experiments/baseline_n3/PREREGISTRATION.md` before any
   implementation or run, and record the commit hash in both reports.
3. Study F uses exactly the double-estimator protocol: base seed `20261201`, primary seed
   `SeedSequence([20261201, d, n, M, rep, chunk_index])`, the same chunk size, the same fixed 1,200 test points with the
   20% validation split (verdicts on the 960 test points), float64, 10 repetitions. Reuse every existing row (`path` and
   `double_path` at (3,16); `path` at n=2, M <= 128); do not rerun them.
4. Verdicts are mechanical; report every failure; no post-hoc thresholds; exploratory analyses labelled as such.

## Background

In MLP the terminal block is centred by g(x): it uses (g(X_T) - g(x)) G/sqrt(h), which has the same mean because
E[g(x) G] = 0 (a control variate, exactly the baseline of REINFORCE). The correction levels are **not** centred: with
Beta(1/2,1) time sampling the level-1 gradient increment is mean_j 2 sqrt(h) D_j G_j with D_j = f(U^1(R_j, X_j)), so the
whole magnitude of f multiplies the d-dimensional Gaussian weight. This puts variance of order h E[D^2]/M into every
coordinate of z_hat, and at the next level the curvature of f turns that variance into a bias proportional to d. The
pathwise terminal gradient cannot remove it (the integrand is a Monte Carlo output, not a differentiable function).
Evidence (path-frontier study): `path` at n=3 on P1 has skill ~ 0.42 d/M (M x skill = 42-46 for M = 4..16; ratio
d=400/d=100 = 3.9 at every M). `double_path` removes the bias, but at d=400 the residual variance keeps n=3 above the
n=2 floor (0.165 at (3,16) versus 0.127).

## Method: correction-level baseline

At every recursive call `solve(n, t, x)` with n >= 2, after the terminal estimate T(t,x) has been computed (pathwise in
the path variants), compute `b(t,x) = f(T(t,x))` with the plain generator (one generator call per parent; count it).
In **level 1 only**, replace the gradient increment by

    mean_j  w_j * (D_j - b) * G_j / sqrt(dt_j)        (w_j = the existing importance weight)

Leave the value increment and all levels l >= 2 unchanged. b depends only on this node's terminal draws, which are
independent of the level-1 transition draws (R_j, G_j), so the expectation is unchanged. The baseline consumes no random
draws, so the primary draw stream of each baseline method must be identical to that of its non-baseline counterpart
(audit the fingerprints).

| method | definition |
|---|---|
| `path_base` | `path` + baseline |
| `double_path_base` | `double_path` + baseline (D_j is the double-generator value; b uses the plain f) |
| `path_base_oracle` | diagnostic only: b = f(exact state at (t,x)) |
| `path`, `double_path` | existing implementations (reuse rows at M=16) |

The root value at n=2 is unaffected by the baseline (only z changes), so no n=2 baseline runs are needed.

Required unit tests (in addition to pairing audits): (i) replacing b by any constant leaves the mean of the level-1
gradient increment unchanged (statistical test, 5 standard errors); (ii) on a linear-terminal HJB
(g(x) = x_1, f(z) = -|z|^2/2, sigma = sqrt(2), d = 100) the level-1 gradient increment of `path_base` is exactly zero
at every node.

## Study F: grid

- PDEs: P1 and MR (k=10, scale=4, T=0.1); d in {100, 400}; 10 repetitions.
- n=3, M in {16, 24, 32}: `path`, `double_path` (new rows at M = 24, 32), `path_base`, `double_path_base`, `path_base_oracle`.
- n=2 floor reference: `path` at M in {256, 512}.

Record per row, in addition to the existing fields:
- root gradient skill RMSE(z_hat - z*)/RMS(z*) on the test points;
- Gaussian vectors drawn (terminal + transition, including auxiliary trees) and the same multiplied by d;
- for the root level-1 gradient increment: the per-coordinate sample variance of the integrand over siblings, averaged
  over coordinates and test points.

Cost axes: (a) mean generator calls (primary, for continuity with D-4/E-1); (b) Gaussian vectors x d, the dominant
arithmetic for large d (generator calls undercount the terminal block); (c) wall time.

## Study F: pre-registered criteria

- **F-1 (baseline repairs `path` at n=3):** `path_base` mean skill <= 0.5 x `path` mean skill in at least 10 of the 12
  (PDE, d, M) n=3 combinations.
- **F-2 (baseline complements double):** `double_path_base` mean skill <= `double_path` mean skill in at least 10 of the
  12 combinations.
- **F-3 (n=3 below the n=2 floor at d=400):** for each PDE at d=400, the minimum over n=3 cells and methods
  {`double_path`, `double_path_base`} of mean skill is <= 0.9 x the minimum over all n=2 `path` cells
  (M = 8 ... 512, including reused rows). PASS if it holds for both PDEs. Report the same quantity at d=100 (no verdict).
- **F-4 (reported, no verdict):** frontiers (lower envelope of mean skill vs cost) of the n=2-only `path` family and of
  the combined family {`path`, `path_base`, `double_path`, `double_path_base`} over n in {2,3}, on each of the three cost
  axes, at 10 log-spaced cost levels where both exist; the attaining cell at each level; 1,000-draw bootstrap over reps.
- **Predictions (no verdict):** `path_base_oracle` <= `path_base`; the root level-1 gradient-increment variance of
  `path_base` is far below that of `path`; for `path`, M x skill / d is roughly constant.

## Study G: GPU port

- **G-0 (environment):** record GPU model, driver, CUDA, framework versions, FP64 throughput class. If no GPU is
  available, implement and validate on the framework's CPU backend and mark G-4 and G-5 "not evaluated (no GPU)".
- **G-1 (implementation):** PyTorch (or JAX), **float64 only**, breadth-first batched like the numpy code, supporting
  `raw`, `path`, `double`, `double_path`, `path_base`, `double_path_base` for P1 and MR; points processed in chunks to fit
  memory. Re-implement the equations; verify terminal, exact state and generator against the numpy equations to relative
  1e-12 on random inputs. Seeding: a device generator seeded from
  `SeedSequence([20261201, d, n, M, rep, chunk_index]).generate_state(2)`, identical for all methods of a cell, so GPU
  methods are paired with each other (they are not paired with CPU rows; say so).
- **G-2 (exactness, verdict):** with a draw-tape interface, record every primary and auxiliary draw of the numpy solver
  for one chunk and replay it in the GPU solver. For P1 and MR at d=100, cells (2,8), (3,4), (3,6), all six methods:
  every output (u and z) agrees with numpy to relative 1e-10. PASS only if all agree.
- **G-3 (statistical agreement, verdict):** with native GPU randomness, P1 and MR, d in {100, 400}, cells (2,32) and
  (3,10), methods `path`, `double_path`, `path_base`, `double_path_base`: 30 GPU repetitions versus 30 new CPU repetitions
  (stored in `results/gpu_mlp/`). PASS if |difference of mean skill| <= 3 x Welch standard error in every cell.
- **G-4 (speed, reported):** wall time per repetition (1,200 points) for `double_path_base`, P1, d=400, cells (3,16),
  (3,32), (4,6), (4,10), (4,16), on GPU and on CPU (numpy, same machine, report core count); peak GPU memory. CPU times
  that would exceed 2 hours may be extrapolated from one chunk; label them as extrapolated.
- **G-5 (exploratory, only if G-2 and G-3 pass):** on GPU, n=4 with M in {6, 10, 16}, methods `path`, `double_path`,
  `double_path_base`; P1 and MR; d in {100, 400}; 10 repetitions. Report skill and all three costs next to the best n=2
  and n=3 results from Study F. No verdict.

## Outputs

```
results/baseline_n3/  rows.csv  summary.csv  frontier.csv  frontier_levels.csv  analysis_summary.json  figures/
results/gpu_mlp/      validation.csv  agreement.csv  speed.csv  n4_rows.csv  n4_summary.csv  figures/
docs/BASELINE_N3_REPORT.md   docs/GPU_MLP_REPORT.md
```

Figures: (F) skill versus M at n=3 for the five methods, one panel per (PDE, d), with the n=2 floor as a horizontal
line; frontiers on the three cost axes. (G) GPU/CPU speed-up versus cell; n=4 skill versus M.
Each report starts with an outcome table of its criteria, then each criterion verbatim with its verdict, then failures
and deviations, then provenance (prereg commit, code commit, environment, totals, pairing audit).

## Priority if compute is short

F at M=16 and 24 -> F at M=32 -> G-1/G-2 -> F n=2 M=256 -> G-3 -> G-4 -> G-5 -> F n=2 M=512.
Never drop a criterion silently; mark it "not evaluated (compute)".
