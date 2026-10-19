# ICML final benchmark suite report

**Execution status:** the corrected single-seed 100D comparison and all native
paper-fidelity anchors are complete as of 2026-10-19. The post-100D dimension
sweep is intentionally blocked pending independent review. This document does
not assign a full-grid Scientific GO, Benchmark GO, or compute-efficiency
verdict while that gate remains closed.

## Frozen protocol

The suite starts from commit
`d36ab9e79d0e6ffe49c43d60cdb28a064dc9ac09` and uses seven PDE families, a
six-dimension grid, five solver entries, and one fixed primary training seed
(`0`) per method cell. The single-seed amendment preserves earlier multi-seed
artifacts as historical evidence but excludes them from rankings and
completion counts. The machine-readable contract is `FROZEN_PROTOCOL.json`.
All common test banks were generated before final training, and every
completed result records its test-bank hash.

The 2026-10-19 paired-bank remediation subsequently opened every frozen
prediction artifact and separated numerical array hashes from NPZ file hashes.
It established that the reported nonconvex 960/512 mismatch was a metadata
error: all four methods had used the same 512 states and reference. Allen–Cahn
was also recomputed against one shared high-accuracy radial reference. The
strict artifact-level audit now passes on all seven 100D PDEs; full details and
pre-repair values are preserved in `PAIRWISE_BANK_REMEDIATION_20261019.md`.

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

The canonical aggregation contains 35 completed PDE-dimension-method cells out
of 210 requested cells, 29 whole-function rows, 44 native-point rows, 41
accuracy-frontier rows, 28 dimension-scaling rows, and 111 ablation rows.
These counts include exact, hash-audited historical reuse. Smoke runs and
superseded baseline profiles are excluded. The extra whole-function row is the
source-native SDGD 250D fidelity anchor; it is deliberately excluded from the
predeclared `10/20/50/100/200/400` dimension-scaling table.

Completed 100D evidence includes:

- Tangent-Actor on HJB-Log, OU/GMM, joint concave, joint convex, nonconvex
  coupled, the nonsmooth MLP-failure PDE, and Allen–Cahn.
- Corrected DPI, PINN-HTE, and SDGD-PINN whole-function results on all seven
  families at their complete published budgets.
- Native SDGD-PINN on HJB-Log at 250D and native DPI/PINN-HTE on OU/GMM at
  100D for the paper-fidelity review.
- Native Deep BSDE on all seven 100D families, reported only in the separate
  single-point table.

The remaining supported cells are the post-review dimension sweep. No such job
has been launched because `REVIEW_GATE_STATUS.json` remains
`pending_independent_review`.

## Audited numerical findings to date

The following seed-0 Value RMSEs are measured endpoints. Tangent-Actor rows use
the frozen 64-step primary protocol; finer Euler results remain numerical
audits in `ACCURACY_COMPUTE_FRONTIER.csv`.

| PDE, 100D | Tangent-Actor | PINN-HTE | SDGD-PINN | DPI |
|---|---:|---:|---:|---:|
| HJB-Log | `0.000814` | `0.053621` | `0.039884` | `0.133935` |
| OU/GMM | `0.219763` | `38.6698` | `62.4838` | `17.6638` |
| Joint concave | `0.000323` | `0.003226` | `0.023332` | `0.002142` |
| Joint convex | `0.000212` | `0.003267` | `0.012780` | `0.001857` |
| Nonconvex coupled | `0.000804` | `0.003068` | `0.008877` | `0.001699` |
| MLP failure | `0.006662` | `0.039939` | `0.066871` | `0.010021` |
| Allen-Cahn | `0.000173` | `0.000948` | `0.002692` | `0.001067` |

Tangent-Actor has the lowest measured seed-0 whole-function RMSE on all seven
100D PDEs under the frozen budgets. This is evidence for a 100D accuracy
advantage, with important limits: most external rows are faithful transfers
rather than native paper tasks; MLP-failure rows are adapted to a nonsmooth
driver; Tangent-Actor usually spends more compute than DPI/PINN-HTE; and one
seed provides no training-seed uncertainty. The result does not establish a
universal or full-dimension ranking.

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

### Deep BSDE native 100D task

The five missing source-faithful transfers completed in Slurm job `3406569`;
the HJB-Log and Allen–Cahn native anchors were already available as exact
audited reuse. Every new run used the pinned DeepBSDE commit
`e76fed80b1995daf4a0dd6d3e2cf64a0931dabda`, 20 BSDE time intervals, two
width-110 hidden layers, 2,000 optimizer iterations, float64 arithmetic, and
the full training cost for its single initial state. All runs had zero
nonfinite values.

| 100D native PDE | Deep BSDE absolute error | Final clipped terminal MSE |
|---|---:|---:|
| OU/GMM | `71.0813` | `287.936` |
| Joint concave | `3.812e-4` | `4.391e-3` |
| Joint convex | `1.252e-3` | `4.486e-3` |
| Nonconvex coupled | `2.180e-3` | `1.337e-2` |
| MLP failure | `6.646e-2` | `1.976e-1` |

The joint-concave Deep BSDE transfer is a meaningful competitor win on the
native point: its `3.812e-4` error is below the selected Tangent-Actor's
64-step MCV native error `5.596e-4`, while using about 61 seconds of training
instead of about 1,895 seconds of fitting plus selection. On joint convex, the
ordering reverses: Tangent-Actor has 64-step MCV native error `3.278e-4` and
128-step audit error `4.187e-5`, at much higher training cost. OU/GMM remained
far from its native reference at the fixed Deep BSDE budget and is labeled
compute limited or poorly matched to this transfer rather than converged.
These scalar comparisons remain separate from the amortized whole-function
table.

### DPI whole-function 100D task

The six missing DPI transfers completed in Slurm job `3406624`; the native
OU/GMM result remains exact audited reuse. Every new run used the pinned
upstream commit `e8616c2f531c9f95299b246a13c1180b8f02822a`, all 40 Picard
rounds, 16 optimizer epochs per round, 4,096 training states, the full source
path budget, and the source PISGradNet. PyTorch 2.6 initially rejected the
source's trusted Lightning checkpoint objects after round one. Job `3406600`
preserves that 56-second failed attempt. The corrected runner restores the
pre-2.6 `torch.load` semantics without changing the pinned checkout or any
numerical update; a real failed-run checkpoint and all later round reloads
passed the compatibility audit.

| 100D whole-function PDE | DPI RMSE | Relative L2 | Training seconds |
|---|---:|---:|---:|
| HJB-Log | `0.133935` | `0.029215` | `447.5` |
| Joint concave | `0.002142` | `0.003606` | `1244.3` |
| Joint convex | `0.001857` | `0.004238` | `1234.0` |
| Nonconvex coupled | `0.001699` | `0.003621` | `4011.6` |
| MLP failure | `0.010021` | `0.009241` | `261.4` |
| Allen-Cahn | `0.001067` | `0.023284` | `414.5` |

Together with the reused OU/GMM DPI RMSE `17.66`, Tangent-Actor has lower
measured 100D whole-function error on all seven families. This is a single-seed
DPI comparison outside its native OU/GMM problem, and Tangent-Actor usually
spends substantially more fitting time. The nonconvex DPI transfer is a cost
exception: it used `4011.6` seconds, more than Tangent-Actor's `3098.4` seconds
of training plus `25.1` seconds of primary inference.

### Corrected PINN-HTE and SDGD-PINN transfers

The corrected PINN-HTE runs use the source HTE-16 configuration only for its
native OU/GMM problem and the paper's generic semilinear PINN profile on the
other PDEs. Their 100D RMSEs are `0.053621`, `38.6698`, `0.003226`, `0.003267`,
`0.003068`, `0.039939`, and `0.000948` in the PDE order of the main table.
Training took between `281` and `1241` seconds. The native OU/GMM result remains
CONDITIONAL in the paper-fidelity review because its rMAE and gradient rMAE are
worse than the approximate paper bars.

SDGD-PINN completed every 10,000-epoch fixed budget with exact RNG-preserving
checkpoint continuation. Its corresponding RMSEs are `0.039884`, `62.4838`,
`0.023332`, `0.012780`, `0.008877`, `0.066871`, and `0.002692`; training took
`5396--14701` seconds. The source-native 250D HJB-Log anchor reached relative
L2 `0.009072` in `204.39` minutes, versus the paper's five-seed value
`0.006147` in `210` minutes. Its source fidelity and runtime match, but the
accuracy gap gives a CONDITIONAL verdict. The broad-test-law OU/GMM transfer
is reported as a distribution/profile limitation and is not used to discredit
SDGD on its native published HJB benchmark.

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

The final-suite incremental Slurm ledger records `35.2136` allocated GH200
GPU-hours since the suite's accounting start. Of these, `11.6939` hours belong
to jobs whose final Slurm state is COMPLETED and `23.5197` to jobs marked
TIMEOUT/CANCELLED. The latter is not discarded compute: the two-hour SDGD jobs
checkpointed every 100 epochs and resumed the exact model, optimizer,
scheduler, collocation, and RNG states in the next allocation. The accounting
classification follows Slurm's final state while the scientific work is
retained.

Separately labeled historical source experiments account for `21.9911`
GPU-hours and are not double counted. The interactive queue's allocation-credit
multiplier is not folded into `RESOURCE_ACCOUNTING.csv`; that file reports
physical allocated GPU-hours. Failed pre-allocation attempts remain visible
with zero usage.

No speed claim is currently supported. The final comparison will separate
Actor fitting, reference-free model selection, frozen-model inference, and
total allocated compute. Published baselines that are still improving at
their fixed budget will be labeled compute limited.

Using the completed seed-0 100D wall times, the blocked 175-cell post-gate
dimension sweep is estimated at `139.53` GPU-hours. The deliberately wide
extrapolation range is `81.67--284.82` GPU-hours. SDGD dominates the point
estimate at `94.14` hours, followed by Tangent-Actor `21.35`, DPI `13.87`,
PINN-HTE `6.29`, and Deep BSDE `3.89`. Queue wait is excluded.

## Independent artifact checks

`FINAL_ARTIFACT_AUDIT.json` passes the executable consistency
checks: Deep BSDE is absent from the whole-function table, Deep BSDE costs are
charged in the native table, smoke artifacts are absent, main scientific keys
are unique, metrics and nonfinite counts are valid, reference uncertainties
and test-bank hashes are complete, source provenance and source artifacts are
present, and Allen–Cahn has no artificial MCV row. The suite status remains
`in_progress` because the independently gated dimension sweep is pending.

## Review-gated next step

The 100D comparison, native fidelity anchors, convergence plot, corrected
tables, and single-seed compute estimate are complete. The next action is an
independent review of `BASELINE_PAPER_FIDELITY_REVIEW.md`. The remaining
dimension extremes and intermediate grid must not be launched until
`REVIEW_GATE_STATUS.json` is explicitly changed to `approved`. If approved,
the sweep will retain one seed per method cell and the original fixed budgets;
if not approved, the frozen 100D results remain available without being used
to justify a broader ranking.

The canonical numerical sources are `MAIN_RESULTS.csv`,
`NATIVE_POINT_RESULTS.csv`, `ACCURACY_COMPUTE_FRONTIER.csv`,
`DIMENSION_SCALING.csv`, `ABLATION_RESULTS.csv`, and
`RESOURCE_ACCOUNTING.csv`. Raw JSON, checkpoints, prediction arrays, cached
labels, and Slurm logs remain in the shared workspace and are hash referenced
by the reproducibility manifest.
