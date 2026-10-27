# Addendum: matched-accuracy amortization on a joint nonseparable PDE

This addendum extends the completed Allen--Cahn 100D study on branch
`ir-mlp-single-point-break-even`.  It must be committed unchanged before any
new joint-PDE bank is generated or any joint-PDE solver is run.  Allen--Cahn
code, banks, raw records, decisions, and verdicts remain frozen.

## 1. Why this PDE is chosen

The second PDE is the already validated 100D `joint_nonseparable` benchmark
from branch `ir-mlp-icml-joint-pde`, not a problem selected after inspecting
the Allen--Cahn timing curves.  Before this addendum, its independent 32-query
study reported RMSE `0.00103796` for Tangent-Actor, `0.0515340` for the
preregistered FH-MLP cell `(n=3,M=16)`, `0.0347129` for FH-MLP `(n=3,M=24)`,
and `0.0348293` for Deep BSDE.  These values may motivate targets and candidate
budgets, but none of those queries or stochastic draws may enter this study.

This PDE is scientifically useful because it has genuine joint value-gradient
dependence and a nonzero mixed derivative, whereas Allen--Cahn has a
value-only driver.  It tests whether the Allen--Cahn amortization conclusion
persists under a materially different semilinear mechanism.  A failure to
observe an FH-MLP crossing is retained as a negative result.

## 2. PDE and reference

At `d=100` and `T=0.5`, solve

\[
  \partial_t u+\Delta u+f(u,z)=0,\qquad z=\sqrt2\nabla u,
\]

with

\[
  f(y,z)=y-y^3-\tfrac12\lvert z+yv\rvert^2,
\]

and

\[
  g(x)=0.5+0.15\sin(\beta_1^\top x)
            +0.15\cos(\beta_2^\top x).
\]

The vectors `beta_1`, `beta_2`, and `v=beta_1` are exactly those defined by
the validated joint benchmark: sign-canonicalized thin QR of the matrix drawn
from `default_rng(SeedSequence([20261232, 100, 17]))`.  The driver satisfies
`nabla_z f_y=-v != 0`.

The ambient reference is obtained only from the independent two-dimensional
periodic Fourier solver used by the validated joint study.  A fresh primary
`256 x 256` reference and independent `128 x 128` audit must satisfy:

- primary/audit value RMSE at most `5e-5` and maximum difference at most
  `2e-4`;
- PDE-residual RMSE at most `2e-6`;
- terminal interpolation maximum error at most `1e-4`;
- solution values in `[0,1]` up to `1e-8`.

Reference failure stops the extension as `BLOCKED_REFERENCE_VALIDATION`.
Truth is forbidden from labels, fitting, stochastic stopping, Actor selection,
candidate selection after calibration, timing, and construction of the
truth-free timing bank.

## 3. Fresh banks and randomness

Use revision `single_point_break_even_joint_v1` and base seed `20261282`.
Every stochastic seed is derived by SHA-256 from revision, stage, method,
configuration, query, replicate, and training seed.

Native queries are `(t,x)=(0,x)` with `x in R^100`:

| bank | size | use | truth |
|---|---:|---|---|
| calibration | 32 | accuracy configuration selection | independent reference |
| final | 64 | frozen accuracy verification | independent reference |
| timing | 1,048,576 | actual cumulative timing prefixes | absent |

Calibration and final each retain the historically mandated exact origin as
row zero; all remaining Gaussian rows are mutually disjoint iid `N(0,I_100)`.
The shared origin is disclosed as the sole overlap.  Timing is entirely fresh
iid Gaussian and contains no reference array.

The timing prefixes are
`1,2,4,...,1048576`.  Deep BSDE uses
`1,2,4,8,16,24,32,48,64` and stops at its active-compute cap if necessary.

## 4. Accuracy targets

- primary: `epsilon_1 = 6.0e-2`;
- secondary: `epsilon_2 = 2.0e-2`.

The primary threshold was chosen because all three methods met it in the
pre-existing study.  The secondary threshold was chosen because Tangent-Actor
and at least one frozen FH-MLP cell met it, while Deep BSDE did not.  A method
is eligible only after passing on calibration and again on the disjoint final
bank.  A more accurate method is never deliberately degraded.

## 5. Candidate configurations and selection

All methods import the audited implementations from
`single_point_joint_extension` and `single_point_32_ice`; adapters for this
addendum live under `experiments/single_point_break_even_joint/`.  Existing
solver modules are not edited.

### 5.1 Tangent-Actor

Use the mature joint terminal-aware Actor unchanged: 8,192 training states,
512 validation states, 128 paths per label, 32 label steps, 4x384 SiLU, 5,000
AdamW steps per update, four updates, joint value-gradient MCV labels, and the
reference-free selector.  Calibration training seeds are `{0,1,2}`; fresh
final/timing seeds are `{10,11,12}`.

Evaluate frozen inference candidates `direct`, `mcv_32`, `mcv_128`, and
`mcv_512`, all with 64 Euler steps and two independent inference replicates.
For each target, select the accuracy-eligible candidate with the smallest
mean measured query time; ties within 2% choose the smaller path budget, with
`direct` smallest.  Common training cost is not removed from break-even
timing.

### 5.2 Full-History MLP

Use the validated float64 recursion with pathwise terminal gradient and
untruncated Bismut level gradients.  Calibration candidates are

`(n,M) in {(2,8),(2,16),(2,32),(3,8),(3,12),(3,16),(3,24)}`,

with five estimator replicates.  For each target, choose the eligible cell
with the smallest mean measured single-query time; ties within 2% choose the
smaller generator-call count and then lexicographic `(n,M)`.  Final accuracy
uses five fresh replicates.

### 5.3 Deep BSDE

Use the audited native scalar-`Y_0` solver: 20 time intervals, batch size 64,
validation size 256, float64, and a fresh network for every query.  Candidate
checkpoints are `{500,1000,2000,4000}` updates from one frozen 4,000-update
trajectory per query and seed.  Calibration and final use three training
seeds per query.  Select the least-update eligible checkpoint.  If no
candidate reaches a target, report `TARGET_NOT_REACHED`; do not substitute an
amortized decoder or another training schedule.

After calibration, freeze all target-specific choices and their canonical
hash in a separate commit before final errors are computed.

## 6. Final accuracy and timing

Final RMSE is computed independently for every stochastic run and pooled over
the frozen runs; nonfinite predictions make the corresponding configuration
ineligible.  No final error can change a configuration.

All runtime comparisons use NVIDIA H100 80GB GPUs only, synchronized at timing
boundaries, with three independent timing replicates.  The measured total at
query count `Q` includes:

- Tangent-Actor: fresh Actor construction, labels, fitting, reference-free
  selection, checkpoint loading, and actual evaluation of the first `Q`
  timing queries;
- FH-MLP: actual cumulative recursive solution of the first `Q` queries;
- Deep BSDE: actual cumulative training of `Q` fresh scalar-`Y_0` models.

Queue time, reference construction, truth evaluation, plotting, hashing, and
upload are excluded from plotted method-active time but remain documented.
No time at `Q` may be obtained by multiplying a one-query measurement.

Each timing replicate has a one-hour active-compute cap.  A time-capped run
reports the last complete frozen prefix.  Primary-target timing has priority;
secondary timing may be marked `NOT_EVALUATED_COMPUTE_CAP` rather than reducing
replicates or silently raising the total cap.

## 7. Break-even and failure rules

For each eligible baseline and target, the observed break-even is the smallest
common executed prefix `Q` for which mean baseline total time is no smaller
than mean Tangent-Actor total time.  Report the interval from the preceding
prefix to `Q`, every per-replicate crossing, and mean plus standard deviation
at every point.  Never interpolate or extrapolate a headline crossing.

If no crossing occurs, report `NOT_REACHED_WITHIN_TESTED_BUDGET`, the largest
common executed `Q`, both measured mean times, and their ratio.  OOM, timeout,
nonfinite values, infrastructure failures, and exact reruns are all retained.

## 8. Compute ceiling and outputs

The joint addendum has an independent ceiling of **18 allocated H100
GPU-hours**, including failed GPU allocations and reruns.  It may not be raised
after results are observed.  If the ceiling prevents a criterion, mark it
`NOT_EVALUATED_COMPUTE_CAP`.

New outputs are isolated under

```text
invariant_region_mlp/experiments/single_point_break_even_joint/
invariant_region_mlp/results/single_point_break_even_joint/
invariant_region_mlp/docs/SINGLE_POINT_BREAK_EVEN_JOINT_REPORT.md
```

Required compact outputs mirror the Allen--Cahn study: bank/reference audits,
calibration rows, frozen choices, final accuracy rows, executed timing rows,
break-even summary, resource accounting, artifact audit, report, PNG, and PDF.
Large checkpoints and the truth-free timing bank remain in the ICE archive,
bound to compact artifacts by hashes.  After this addendum is audited and
pushed, stop; do not start a third PDE.
