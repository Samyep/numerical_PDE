# ICML final benchmark suite report

**Execution status:** staged execution in progress on 2026-10-18. This document
records only completed, audited evidence. It does not assign the final
Scientific GO, Benchmark GO, or compute-efficiency verdict while supported
benchmark cells remain pending.

## Frozen protocol

The suite starts from commit
`d36ab9e79d0e6ffe49c43d60cdb28a064dc9ac09` and uses seven PDE families, a
six-dimension grid, five solver entries, and three requested training seeds.
The machine-readable contract is `FROZEN_PROTOCOL.json`. All common test banks
were generated before final training, and every completed result records its
test-bank hash.

The amortized whole-function task contains Tangent-Actor, DPI, PINN-HTE, and
SDGD-PINN. Deep BSDE remains exclusively in the native single-initial-state
task, where its full solve cost is charged. The external baselines retain the
authors' source algorithms and budgets; they receive no Tangent-Actor terminal
residual, learned tangent update, martingale control variate, rank projection,
or oracle coordinate input.

The mature Tangent-Actor configuration uses 8,192 training states, 128 paths
per label, 32 label steps, a four-layer width-384 SiLU Actor, 5,000 optimizer
updates per Actor iteration, and four learned updates. Actor fitting and final
reference-free model-selection time are accounted for separately. Analytic or
numerical nonterminal truth is opened only after every candidate checkpoint
and selection decision are frozen.

## Current coverage

The current canonical aggregation contains 12 completed PDE–dimension–method
cells out of 210 requested cells, 18 whole-function rows, 30 native-point rows,
36 accuracy-frontier rows, and 111 ablation rows. These counts include exact,
hash-audited historical reuse. Smoke runs are excluded. The execution manifest
records the actual seed count for every cell rather than manufacturing
replication intervals.

Completed 100D evidence currently includes:

- Tangent-Actor on HJB-Log, OU/GMM, joint concave, joint convex, nonconvex coupled, the
  nonsmooth MLP-failure PDE, and Allen–Cahn.
- Native SDGD-PINN on HJB-Log.
- Native DPI and PINN-HTE on OU/GMM.
- Native Deep BSDE on HJB-Log and Allen–Cahn, reported only in the separate
  single-point table.

The mature 100D joint-convex Tangent-Actor cell and its common-bank
all-iteration diagnostic are complete. The remaining external 100D transfers,
dimension extremes, intermediate
dimensions, and missing independent seeds remain pending.

## Audited numerical findings to date

The following values are measured results, not extrapolated endpoints. Main
rows use the frozen 64-step protocol; finer Euler results remain in
`ACCURACY_COMPUTE_FRONTIER.csv`.

| PDE and method | Seeds | Value RMSE | Relative L2 | Evidence status |
|---|---:|---:|---:|---|
| HJB-Log, Tangent-Actor + MCV | 3 | mean `8.535e-4` | mean `1.862e-4` | exact historical reuse |
| HJB-Log, SDGD-PINN | 1 | `2.992e-2` | `6.527e-3` | native official reuse |
| OU/GMM, Tangent-Actor + MCV | 3 | mean `2.210e-1` | mean `7.126e-4` | exact historical reuse |
| OU/GMM, DPI | 1 | `1.766e1` | `5.697e-2` | native official reuse |
| OU/GMM, PINN-HTE | 1 | `3.867e1` | `1.247e-1` | native official reuse |
| Joint concave, Tangent-Actor + MCV | 1 | `3.231e-4` | `5.441e-4` | new final-suite run |
| Joint convex, Tangent-Actor + MCV | 1 | `2.116e-4` | `4.830e-4` | new final-suite run |
| Nonconvex coupled, Tangent-Actor + MCV | 3 | mean `8.093e-4` | mean `1.725e-3` | exact historical reuse |
| MLP failure, Tangent-Actor + MCV | 1 | `6.662e-3` | `6.144e-3` | adapted minimum-norm subgradient |
| Allen–Cahn, value-only Tangent-Actor | 3 | mean `1.666e-4` | mean `3.636e-3` | exact historical reuse |

The table is incomplete and cannot yet support a full cross-method ranking.
In particular, the current external comparisons cover only their native 100D
anchors. A low error against an unfinished competitor grid is not recorded as
a benchmark win.

### Joint-concave mature run

The new mature joint-concave run completed in Slurm job `3406194`. Four Actor
updates reduced held-out noisy tangent-label loss as follows:

`0.002238 -> 0.0008226 -> 0.0005869 -> 0.0004275`.

The independent noise-corrected fixed-point selector chose `a3`. Its total
inconsistency was `3.32377e-4`, versus `3.33057e-4` for `a4`; the selection did
not use numerical truth. At 128 Euler steps and 1,024 paths, the selected Actor
gave:

- plain policy RMSE `1.85162e-3`;
- MCV policy RMSE `2.13470e-4` and relative L2 `3.59451e-4`;
- MCV finite-path contribution `1.72335e-4`;
- variance reduction `117.58x`;
- exponential-weight ESS fraction `0.98804`;
- zero nonfinite values.

The independent 2D reference uncertainty is `8.454e-5`. Thus the 128-step
result is already close to the combined MC/reference resolution; the
predeclared 256/512-step audit is still needed before assigning a remaining
time-discretization floor.

The common-random-number post-freeze iteration diagnostic (job `3406317`)
found MCV policy RMSEs

`a0 0.003361, a1 0.0004013, a2 0.0003897, a3 0.0003855, a4 0.0003858`.

Most policy-value benefit arrived in the first learned update. Later updates
continued to reduce Actor value and gradient errors, while the final policy
estimate became MC limited. The reference-free choice of `a3` agrees with this
post-hoc truth diagnostic to the available resolution.

### Joint-convex mature run

The new mature joint-convex run completed in Slurm job `3406335`. Four Actor
updates reduced held-out noisy tangent-label loss as follows:

`0.004118 -> 0.001517 -> 0.0008143 -> 0.0005582`.

The independent noise-corrected fixed-point selector chose `a4`, with total
inconsistency decreasing from `0.8132` at `a0` to `3.78926e-4` at `a4`.
Numerical truth was not used for this decision. At 128 Euler steps and 1,024
paths, the frozen Actor gave:

- plain policy RMSE `1.50891e-3`;
- MCV policy RMSE `1.43167e-4` and relative L2 `3.26811e-4`;
- MCV finite-path contribution `1.45160e-4`;
- variance reduction `117.51x`;
- exponential-weight ESS fraction `0.98982`;
- zero nonfinite values.

The independent 2D reference uncertainty is `7.784e-5`. The MCV result is
therefore already at the combined MC/reference resolution. The primary frozen
64-step row is RMSE `2.11576e-4`; finer-step values are numerical audits and
did not participate in Actor selection.

The common-random-number post-freeze diagnostic (job `3406472`, seed
`683824272`) found 64-step MCV policy RMSEs

`a0 0.005115, a1 0.0003368, a2 0.0003062, a3 0.0002991, a4 0.0002984`.

The reference-free selector therefore chose the smallest observed error, but
the `a3`--`a4` difference is far below the finite-path contribution of about
`2.9e-4` and is not statistically resolved. Actor value RMSE decreased from
`0.1007` at `a0` to `0.001969` at `a4`, and Actor gradient-vector RMSE
decreased from `0.06865` to `0.002562`; all ten estimator rows were finite.
The evidence supports a decisive first update and progressively better Actor
fits, while the policy-value metric becomes MC limited after `a1`.

### MLP dimensional-failure mechanism

The reproduced raw Full-History MLP curve uses the audited `n=m=3` recursion
and 512 independent repetitions per dimension. Its RMSEs at
`d=20,50,100,200,400,1000` are respectively

`7.5654, 19.4756, 39.2967, 79.0731, 158.0477, 396.9779`.

All samples were finite. At the newly completed 400D point, the 95% bootstrap
interval is `[154.1189, 162.0082]`. The separately trained 100D Tangent-Actor
has 64-step MCV RMSE
`0.0066623` and 128-step MCV RMSE `0.0063182`. This Tangent-Actor result uses
the documented minimum-norm subgradient at the nonsmooth origin and is labeled
`adapted`; it is not a smooth-PDE convergence theorem. The remaining
Tangent-Actor dimensions are still pending. No claim of fundamental
Full-History MLP nonconvergence is made; this is the frozen low-level mechanism
experiment rather than an asymptotic convergence study.

### MCV and high-accuracy evidence

MCV reduces variance without changing the corresponding discrete policy
value. The observed reduction factors include about `18--20x` on HJB-Log,
`44x` on OU/GMM, `59.95x` at 64 steps and `117.58x` at 128 steps on the new
joint-concave run, and about `50x` on the mature nonconvex coupled runs. It
must not be described as an equal-factor reduction in RMSE or as a correction
for Actor bias.

The mature nonconvex 100D frontier reaches mean MCV RMSE around `1.36e-4` at
512 steps with an independently audited reference uncertainty of `8.08e-6`.
The method remains empirical for the indefinite driver; these experiments do
not establish global policy monotonicity or universal convergence.

## Numerical-reference audit

All seven references are independent of learned solvers. The current
uncertainty proxies are `5.33e-15` for HJB quadrature, `1.19e-12` for the
OU/GMM residual audit, `8.45e-5` for joint concave, `7.78e-5` for joint convex,
`8.08e-6` for nonconvex coupled, floating-point analytic evaluation for the
MLP counterexample, and `9.64e-6` for the 100D Allen–Cahn native anchor.
Dimension-specific Allen references are generated by two-level radial
refinement and remain pending until their own audit succeeds.

## Compute accounting

Through the completed joint-convex iteration diagnostic, the final-suite
incremental allocation is `2.0978` GH200 GPU-hours. The separately labeled
historical source experiments account for `21.9911` GPU-hours and are not
double counted. Failed and cancelled scheduler attempts remain visible in
`RESOURCE_ACCOUNTING.csv`; the corrected pre-allocation cancellation consumed
zero GPU time.

The Full-History implementation is CPU-only, but the DeltaAI GH partition
requires every submitted job to reserve a GPU. Its 11-second scheduler
allocation is therefore charged as `0.0031` GH200 GPU-hours, bringing the
incremental total to `2.1008` hours even though the algorithm did not use GPU
arithmetic.

No speed claim is currently supported. The final comparison will separate
Actor fitting, reference-free model selection, frozen-model inference, and
total allocated compute. Published baselines that are still improving at
their fixed budget will be labeled compute limited.

## Independent artifact checks

`FINAL_ARTIFACT_AUDIT.json` currently passes twelve executable consistency
checks: Deep BSDE is absent from the whole-function table, Deep BSDE costs are
charged in the native table, smoke artifacts are absent, main scientific keys
are unique, metrics and nonfinite counts are valid, reference uncertainties
and test-bank hashes are complete, source provenance and source artifacts are
present, and Allen–Cahn has no artificial MCV row. The suite status remains
`in_progress` because supported method cells are still pending.

## Remaining work before a final verdict

1. Run the source-faithful 100D Deep BSDE native solves and DPI, PINN-HTE, and
   SDGD-PINN whole-function transfers.
2. Execute dimension extremes, then the intermediate dimension grid.
3. Complete missing independent training seeds and the frozen
   64/128/256/512 evaluation frontier.
4. Regenerate all tables and figures, quantify incomplete or resource-limited
   cells, and only then assign Scientific GO, Benchmark GO, and compute
   efficiency separately.

The canonical numerical sources are `MAIN_RESULTS.csv`,
`NATIVE_POINT_RESULTS.csv`, `ACCURACY_COMPUTE_FRONTIER.csv`,
`DIMENSION_SCALING.csv`, `ABLATION_RESULTS.csv`, and
`RESOURCE_ACCOUNTING.csv`. Raw JSON, checkpoints, prediction arrays, cached
labels, and Slurm logs remain in the shared workspace and are hash referenced
by the reproducibility manifest.
