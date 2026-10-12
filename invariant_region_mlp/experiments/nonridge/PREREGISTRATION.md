# Non-ridge benchmarks: do pathwise + baseline + double survive full-rank terminals? (pre-registered, GPU)

Repository `Samyep/numerical_PDE`, subproject `invariant_region_mlp`. Branch from `ir-mlp-baseline-gpu`
(commit `7cb74835` or later; it contains the GPU solver and the baseline methods) and create branch `ir-mlp-nonridge`.
Another machine may be running `ir-mlp-actor-critic`; do not push to, rebase onto or modify it.
You are authorized to run long experiments locally and use all available CPU/GPU resources reasonably,
but do not delete existing results or rewrite the manuscript.

**Speed matters.** Run Stage A first, then commit and push Stage A results with an interim report before Stage B starts.

## Hard rules

1. Do not modify any existing module or result. New code in `invariant_region_mlp/experiments/nonridge/`, results in
   `results/nonridge/`, report in `docs/NONRIDGE_REPORT.md`. Import or subclass the existing GPU solver
   (`experiments/gpu_mlp/solver.py`) and the frozen CPU solvers (`experiments/baseline_n3/solver.py`,
   `experiments/double_estimator/solver.py`) without editing them. If the GPU equation wrapper rejects new families,
   subclass it in the new directory.
2. Freeze this document first: commit it verbatim as `experiments/nonridge/PREREGISTRATION.md` before any
   implementation or run; record the hash in the report.
3. Equation parameters below are fixed now. Do not tune them after seeing any MLP result. If calibration shows a
   degenerate problem, report it and continue with the frozen parameters.
4. Verdict rows come from the GPU solver (float64 only, native GPU randomness, device seed from
   `SeedSequence([20261201, d, n, M, rep, chunk_index]).generate_state(2)`, identical for all methods of a cell so methods
   are paired). 10 repetitions. Test points and the 20% validation split are built exactly as
   `double_estimator.equations.fixed_points` does, with the new `pde_id` strings; verdicts use the 960 test points.
5. Verdicts are mechanical; report every failure; no post-hoc thresholds; exploratory analyses labelled as such.

## Why

Every benchmark so far (P1, MR) is exactly low-rank: z* lies in a fixed k-dimensional subspace, so the d-k remaining
directions carry pure noise and zero signal. A reviewer can say the repairs only work because of this. Here the
terminal function is curved in **all** d directions (full-rank Hessian, z* not confined to any subspace), and an
approximately low-rank member interpolates between the ridge case and the full-rank case.

## Equations: spectral separable HJB family

Same PDE and conventions as P1: sigma = sqrt(2), mu = 0, T = 0.25, z = sigma grad u,
f(z) = -|z|^2 / 2, i.e. u_t + Laplacian(u) - |grad u|^2 = 0, u(T, x) = g(x).

- Profile (the P1 profile): phi(y) = log 3 - log(exp(0.5 y) + exp(1.5 y) + exp(3 y)).
- Rotation: Q x = DCT(s * x), with DCT = `scipy.fft.dct(type=2, norm="ortho")` along the last axis and
  s in {-1, +1}^d drawn as `2 * np.random.default_rng(np.random.SeedSequence([20261201, d, 31])).integers(0, 2, d) - 1`.
  Q is orthogonal; Q^T v = s * IDCT(v) with `scipy.fft.idct(type=2, norm="ortho")`.
- Spectrum: c_j = j^(-alpha) / sqrt(sum_{i=1..d} i^(-2 alpha)), j = 1..d (so sum c_j^2 = 1, matching P1's |grad g| scale).
- Terminal: g(x) = sum_j c_j phi(y_j), y = Q x. Terminal gradient: sigma grad g = sigma Q^T (c_j phi'(y_j))_j.

| pde_id | alpha | character |
|---|---|---|
| `S0` | 0 | full rank, flat spectrum (no low-dimensional structure at all) |
| `S05` | 0.5 | approximately low-rank, effective rank ~ 16 (d=100) / ~ 26 (d=400) |
| `S1` | 1 | near-ridge with a full-rank tail, effective rank ~ 2.5 |

**Exact solution (Hopf-Cole, separable).** u = -log w with w the backward heat flow of exp(-g). Because Q is orthogonal,
the coordinates of Q(x + sqrt(2 tau) G) are independent, so with tau = T - t:

    u(t, x) = sum_j v(c_j, tau, y_j),        v(c, tau, y) = -log E[ exp(-c phi(y + sqrt(2 tau) xi)) ],  xi ~ N(0,1)
    z(t, x) = sigma Q^T ( dv/dy (c_j, tau, y_j) )_j,
    dv/dy  = E[ c phi'(Y) exp(-c phi(Y)) ] / E[ exp(-c phi(Y)) ],   Y = y + sqrt(2 tau) xi.

Evaluate the 1-D expectations with Gauss-Hermite quadrature (K = 96 nodes) in the log domain (logsumexp). The exact
solution is needed only at test points (and, on the CPU reference solver, at generator calls for diagnostics). The heat
(Picard-0) iterate is u0(t, x) = sum_j c_j E[phi(y_j + sqrt(2 tau) xi)] by the same quadrature.

**Equation validation (must pass before any MLP run; record in `results/nonridge/validation.csv`):**
- V-1 quadrature convergence: K = 96 versus K = 192, max relative difference in u and z <= 1e-11 on 200 test points per
  (family, d).
- V-2 special case: with c = e_1 (all weight on j = 1) the quadrature u, z equal the closed-form P1 formula with
  w = Q^T e_1 to relative 1e-11.
- V-3 PDE residual: for the 1-D function v, d v/d tau - d^2 v/dy^2 + (dv/dy)^2 = 0 by centred finite differences
  (step 1e-4) to <= 1e-6 at 100 random (c, tau, y); and z equals sigma times the finite-difference gradient of u along 5
  random directions per point to <= 1e-6.
- V-4 Monte Carlo: Hopf-Cole u by 10^6 samples at 20 points per (family, d) agrees with quadrature within 4 standard
  errors.
- V-5 GPU equations: terminal, terminal gradient, generator and double generator on GPU equal numpy to relative 1e-12
  on random inputs. Implement the DCT on GPU (e.g. via FFT, Makhoul's algorithm) or with a dense orthogonal matrix;
  report which and its time per call.
- V-6 draw-tape exactness: replay recorded numpy draws in the GPU solver for S0 and S05, d = 100, cells (2,8) and (3,4),
  methods `path` and `double_path_base`: all outputs agree to relative 1e-10. For the frozen CPU solver supply
  diagnostic-only attributes that cannot affect estimates (box_low = -inf, box_high = +inf, z_ball_radius = inf,
  no u certificate).

**Calibration table (reported, no verdict), on the 960 test points for each family and d, and for P1 from its closed
form:** std of u*; mean and std of f(z*); RMS |z*| and |mean z*|; skill of the heat iterate u0 as a predictor of u*
(nonlinearity strength); fraction of the **centred** variance of z* (covariance over test points) outside its top-1 and
top-10 principal directions; for the terminal Hessian of g at the test points (diagonal in y-coordinates:
c_j phi''(y_j)), the mean Frobenius norm and the mean effective rank (sum lambda^2)^2 / sum lambda^4.

The design keeps the total curvature fixed and varies only its rank. A pre-check on 600 other random points gave:
Hessian Frobenius norm ~0.89 for every family (P1: 0.87); Hessian effective rank ~86 / 344 (S0, d = 100 / 400),
~16 / 26 (S05), ~2.7 (S1), 1 (P1); centred z* variance outside the top-10 directions ~0.82 / 0.92 (S0), ~0.41 / 0.52
(S05), ~0.05 (S1); heat-iterate skill ~0.38 / 0.27 (S0), ~0.41 / 0.32 (S05), ~0.49 / 0.43 (S1). Every family keeps a
common mean gradient (|mean z*| ~ 2.5, because phi is monotone, as in P1), and for S0 f(z*) is nearly constant over
the test points (std ~0.1-0.2); this is expected for a flat spectrum and is reported, not corrected. Large deviations
from these values indicate an implementation error.

## Methods

The existing GPU methods, unchanged: `raw` (Bismut terminal), `path` (pathwise terminal), `path_base` (path + level-1
baseline b = f(T(t,x))), `double_path` (path + product double generator, f quadratic so -z1.z2/2), `double_path_base`.
Add the root level-1 gradient-integrand variance diagnostic (per-coordinate sample variance over siblings, averaged over
coordinates and test points, as in Study F) to every n >= 2 row if it is not already recorded.

Record per row: skill (value RMSE / std u* on the test points; the verdict metric), root gradient relative L2, generator
calls, Gaussian vectors drawn and the same times d, GPU wall time, the level-1 variance diagnostic, non-finite counts.

## Stage A (verdicts; push results and an interim report immediately when done)

- Families `S0`, `S05`; d in {100, 400}; 10 repetitions.
- n = 2: `path` at M in {16, 64, 256}; `raw` at M = 64.
- n = 3, M = 16: `raw`, `path`, `path_base`, `double_path`, `double_path_base`.
- n = 3, M = 24: `double_path`, `double_path_base`.

## Stage A pre-registered criteria

- **N-1 (pathwise matters off-ridge):** `path` skill <= 0.5 x `raw` skill at (2,64) and at (3,16), in all 8
  (family, d, cell) combinations. A non-finite `raw` mean counts as satisfied.
- **N-2 (baseline off-ridge):** `path_base` skill <= 0.5 x `path` skill at (3,16), in all 4 (family, d) combinations.
- **N-3 (double complements baseline off-ridge):** `double_path_base` skill <= `double_path` skill at (3,16) and (3,24),
  in at least 7 of the 8 (family, d, M) combinations.
- **N-4 (n = 3 below the n = 2 floor at d = 400):** for each family at d = 400, the minimum over the n = 3 cells and
  methods {`double_path`, `double_path_base`} of mean skill is <= 0.9 x the minimum over the n = 2 `path` cells
  (M = 16, 64, 256). PASS if it holds for both families. Report the same quantity at d = 100 (no verdict).
- **N-5 (dimension robustness):** `double_path_base` skill ratio d = 400 / d = 100 <= 1.3 at (3,16) and (3,24), for both
  families (4 comparisons; PASS only if all hold).
- **Reported, no verdict:**
  - the ratios path_base/path, double_path_base/double_path and best-n3/n2-floor next to the same ratios for P1 and MR
    from `results/baseline_n3/summary.csv` (CPU rows, read-only; ratio-level comparison only, rows are not paired);
  - for `path` and `path_base` at (3,16) and (3,24), skill versus V1 * d / M (V1 = level-1 variance diagnostic) on
    log-log axes with a fitted slope per family;
  - equal-cost frontiers on the generator-call and Gaussian-coordinate axes, with the attaining cells.

## Stage B (after the interim push)

1. `S1` on the full Stage-A grid, with N-1 ... N-5 evaluated for `S1` as a separate, also pre-registered, block (same
   thresholds; N-4 and N-5 for `S1` alone).
2. All three families: n = 3, M = 32 for `double_path` and `double_path_base`; n = 2 `path` at M = 512 (extends the
   floor; re-evaluate N-4 including it and report both versions).
3. Exploratory (no verdict): n = 4, M in {6, 10}, `double_path_base`, d = 400, all three families.

## Outputs

```
results/nonridge/  validation.csv  calibration.csv  rows.csv  summary.csv  frontier.csv  analysis_summary.json  figures/
docs/NONRIDGE_REPORT.md
```

Figures: skill versus M per method with the n = 2 floor as a horizontal line, one panel per (family, d); the V1-collapse
plot; skill versus d for `double_path_base`; frontiers. The report starts with an outcome table (N-1 ... N-5 for Stage A,
then the S1 block), then each criterion verbatim with its verdict, the calibration table, failures and deviations, and
provenance (prereg commit, code commit, GPU/driver/CUDA/PyTorch versions, totals, validation results, timings).

## Priority if compute is short

V-1 ... V-6 -> Stage A at d = 400 (both families) -> Stage A at d = 100 -> interim push -> S1 -> M = 32 and M = 512 ->
n = 4. If one Stage-A repetition of a cell takes more than 30 minutes on the GPU, run 5 repetitions for that cell and
label it. Never drop a criterion silently; mark it "not evaluated (compute)".
