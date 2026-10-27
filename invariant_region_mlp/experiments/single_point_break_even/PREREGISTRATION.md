# Matched-Accuracy Amortization Break-Even

Status: **FROZEN BEFORE NEW BANK GENERATION, CALIBRATION, OR SOLVER EXECUTION**  
Revision: `single_point_break_even_v1`  
Freeze date: 2026-10-27  
Branch: `ir-mlp-single-point-break-even`  
Base branch: `ir-mlp-icml-single-point-32-ice`  
Base commit: `9ce44b4393d2e54ed461b69a214af0d2019bb304`

## 1. Scientific question

At a common, prespecified numerical accuracy, how many distinct initial-state
queries must be solved before the measured total computation of a native
single-point method exceeds that of one globally trained Tangent-Actor?

The primary benchmark is Allen-Cahn in dimension 100.  Query time is always
`t=0`, and states follow the previous audited single-point law: the first row
is the exact origin and all remaining rows have independent `N(0,I_100)`
entries.  The equation, horizon, terminal condition, diffusion convention,
and reference implementation are unchanged from `single_point_32_ice`.

The methods are:

1. `tangent_actor`: one global Actor is trained once; every later state is a
   query to that frozen Actor, optionally followed by the selected MCV
   refinement.
2. `fh_mlp_standard`: the audited Full-History MLP recursion is rerun for each
   state; it has no training phase.
3. `deep_bsde`: the audited native scalar-`Y_0` Deep BSDE model is independently
   initialized and trained for every state.  It is not converted into a
   whole-function learner.

Previous 32-query errors (`5.2279e-5`, `2.1472e-4`, and `1.7769e-4` for the
three methods, respectively) are feasibility evidence only.  Those queries,
predictions, and errors are forbidden from new configuration selection and
formal accuracy evaluation.

## 2. Accuracy targets and interpretation

The common targets are frozen as

- primary: `epsilon_1 = 2.5e-4`;
- secondary: `epsilon_2 = 1.0e-4`.

A method is eligible for a target only when the configuration selected on the
calibration bank also achieves pooled RMSE at or below that target on the
disjoint final accuracy bank.  A more accurate method is not degraded to
force equal errors.  If no calibration candidate reaches a target, report
`TARGET_NOT_REACHED` together with the best measured calibration RMSE and its
actual cost.  If a frozen configuration passes calibration but fails final
accuracy, report `FINAL_ACCURACY_FAILED`; do not select a replacement.

Only accuracy-eligible methods receive a break-even verdict for that target.

## 3. Frozen randomness and query banks

The namespace seed is `20271027`.  Seeds are derived with NumPy
`SeedSequence`, recorded in manifests, and never reused across stages.

Three immutable banks are generated before any solver touches them:

| Bank | Size | Purpose | Truth access |
|---|---:|---|---|
| calibration | 32 | select the cheapest eligible numerical configuration | evaluation after each fixed candidate completes |
| final accuracy | 64 | one-shot qualification of the frozen configuration | evaluation only |
| timing | 524,288 | actual prefix timing | no truth is generated or read |

All banks contain distinct states and are disjoint from the previous 32-query
bank.  The calibration and final banks each contain the origin as row zero and
independent Gaussian rows thereafter.  The timing bank contains independent
Gaussian rows; it is used only for computation measurement.

The frozen timing query counts are

`1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384,
32768, 65536, 131072, 262144, 524288`.

Deep BSDE is measured only through
`1, 2, 4, 8, 16, 24, 32, 48, 64`, because it performs a fresh training at
each state.  These are actual cumulative runs, not a line fitted to one point.

Calibration training seeds are `(0,1,2)`.  Final training seeds are
`(10,11,12)` and are disjoint.  FH-MLP uses five independent estimator
replicates in calibration and final accuracy, and three independent timing
replicates.  Tangent-Actor uses three training seeds and, for stochastic MCV
inference, two independent inference replicates per Actor during calibration
and final accuracy.  Direct Actor evaluation is deterministic conditional on
the trained Actor.

## 4. Reference gate

The convergence-audited radial Allen-Cahn solver is copied through an adapter;
the source implementation is not edited.  The calibration and final
references are recomputed with at least three refinement levels and a domain
expansion check.  Before calibration:

- all reference values and uncertainties must be finite;
- reference uncertainty RMS must be at most `5e-6`;
- maximum pointwise uncertainty must be at most `1e-5`;
- terminal identity and a 200-state PDE residual audit must pass;
- hashes of `t`, `x`, truth, and uncertainty arrays are frozen.

Failure stops the study as `BLOCKED_REFERENCE_GATE`.  Truth is prohibited from
training, Actor selection, stochastic stopping, timing-bank construction, and
runtime stopping.

## 5. Candidate configurations and mechanical selection

### 5.1 Tangent-Actor

Training is the frozen mature Allen-Cahn value-Actor configuration from the
audited source: four Actor updates, 8,192 training states, 128 paths per label,
32 label steps, a 4x384 SiLU network, 5,000 optimizer steps per update, MCV
labels, and the reference-free fixed-point selector.  No training
hyperparameter is searched.

The only calibration choice is query evaluation:

- `direct`: the selected frozen Actor value head, with no Feynman-Kac paths;
- `mcv_P`, with `P in {32,64,128,256,512,1024}`, 64 Euler steps and the
  audited value-Actor martingale control variate.

Direct evaluation is a legitimate Tangent-Actor query, not an oracle: the
Actor is the learned global approximation and receives only `(t,x)`.  For
each target, select the eligible candidate with the smallest synchronized
mean query time.  Ties within 2% choose the smaller path budget, with `direct`
ordered first.

### 5.2 Full-History MLP

Use the audited standard estimator unchanged: float64, pathwise terminal
gradient, untruncated Bismut level gradients, beta-time sampling with
`alpha=0.5`, and independent high/low recursive subtrees.  Candidate `(n,M)`
configurations, in deterministic tie order, are

`(3,12), (3,16), (3,24), (3,32), (3,48), (3,64), (3,96),
 (4,6), (4,8), (4,10), (4,12), (4,16)`.

For each target, select the eligible candidate with the smallest synchronized
mean inference time over the calibration bank and five estimator replicates.
No truth-selected envelope is allowed after final evaluation.

### 5.3 Deep BSDE

Use the audited fixed-initial-state official-compatible PyTorch port: 20 time
intervals, float64, two hidden layers of width `d+10`, batch 64, validation
size 256, learning rate `5e-4`, and the native trainable scalar `Y_0`.
Every query and seed starts a fresh model.

One trajectory is trained through 8,000 updates.  Only the prespecified
checkpoints at `2,000`, `4,000`, and `8,000` updates are accuracy candidates;
their `Y_0` values are recorded during training.  For each target, select the
earliest checkpoint whose pooled calibration RMSE meets the target.  If none
does, report `TARGET_NOT_REACHED`.  Validation loss may be recorded but does
not select a checkpoint.

## 6. Freeze barrier

Calibration produces `calibration_rows.csv` and
`FROZEN_SELECTIONS.json`.  The latter records every candidate's RMSE, time,
failure count, selected configuration per target, array hashes, code commit,
GPU model, and a canonical SHA-256.  It is committed before final jobs are
submitted.  After this barrier, configuration changes are forbidden.

## 7. Final accuracy and actual runtime measurement

All formal computation uses one NVIDIA H100 80GB per run.  GH200, RTX 5070 Ti,
CPU, and mixed-hardware times are excluded from the primary curves.

Primary time is synchronized method-active wall time on the H100:

- Tangent-Actor: Actor initialization, all label generation, fitting, the
  reference-free selector, checkpoint load, and selected query evaluation;
- FH-MLP: solver construction and the complete recursive computation for
  every query;
- Deep BSDE: fresh model construction, all frozen training updates, and
  extraction of the scalar `Y_0` for every query.

Queue time, Slurm startup, reference construction, truth evaluation,
aggregation, plotting, hashing, and artifact upload are excluded from primary
time but task wall time and allocated GPU time are also reported.  CUDA is
synchronized at timing boundaries.  A fixed untimed CUDA warm-up is allowed
and recorded.  Failed attempts and reruns remain in resource accounting.

Final accuracy uses the 64-state final bank.  No final error can change a
configuration.  Runtime uses timing-bank prefixes and records every requested
prefix independently or as an exact cumulative checkpoint in a single
sequential run.  It is forbidden to multiply one single-point time by `Q` and
present that extrapolation as measured.

For Tangent-Actor, total cost at `Q` is the measured fresh training and
selection cost plus measured frozen-Actor query cost.  For FH-MLP it is the
measured cumulative recursive inference cost.  For Deep BSDE it is the
measured cumulative cost of `Q` independent trainings.

Each timing replicate has a one-hour active-compute cap.  A method that reaches
the cap before its next prefix records `TIME_CAP` and no partial-prefix RMSE or
time.  FH-MLP and Tangent-Actor run until the first measured crossing and one
additional frozen grid point, or the largest grid point.  Deep BSDE runs
through 64 queries unless numerical failure or its time cap intervenes.

## 8. Break-even rule

For each eligible baseline and target, the observed break-even is the smallest
common measured query count `Q` for which

`mean total baseline time(Q) >= mean total Tangent-Actor time(Q)`.

Report the previous tested count and the first crossing count as an interval.
Do not interpolate a headline crossing.  If no crossing is observed, report
`NOT_REACHED_WITHIN_TESTED_BUDGET` and the largest common measured `Q`, both
actual times, and their ratio.  Also report per-replicate crossings and
mean/standard deviation of time at every measured point.

If Deep BSDE is accuracy-ineligible, numerically nonfinite, or unable to
complete the frozen schedule within its cap, state this briefly and make the
FH-MLP comparison the primary baseline result.  Do not replace Deep BSDE with
an amortized decoder.

## 9. Compute limits and failure policy

Allen-Cahn has an 18 H100 GPU-hour total allocation cap, subdivided as a
planning limit of 3 hours for Tangent-Actor, 4 hours for FH-MLP, 10 hours for
Deep BSDE, and 1 hour for audits/contingency.  A subdivision may be rebalanced
without changing numerical configurations, but the total cap may not be
raised post hoc.

OOM, timeout, nonfinite values, missing artifacts, and scheduler failures are
recorded.  An infrastructure-only rerun may use the identical numerical seed
and configuration; both allocations count.  Numerical failures are not
silently retried with new hyperparameters.

## 10. Outputs

New code and results are isolated under

```
invariant_region_mlp/experiments/single_point_break_even/
invariant_region_mlp/results/single_point_break_even/
invariant_region_mlp/docs/SINGLE_POINT_BREAK_EVEN_REPORT.md
```

Required artifacts are

```
bank_manifest.json
reference_audit.json
calibration_rows.csv
FROZEN_SELECTIONS.json
final_accuracy_rows.csv
timing_rows.csv
break_even_summary.csv
resource_accounting.csv
artifact_audit.json
analysis_summary.json
figures/matched_accuracy_break_even.png
figures/matched_accuracy_break_even.pdf
```

The figure has one panel per accuracy target, log query count on the x-axis,
and measured total H100 seconds on the y-axis.  It shows actual points with
mean +/- standard deviation, labels ineligible methods directly, and never
draws an extrapolated crossing.

## 11. Optional second PDE

No second PDE is allowed to start before Allen-Cahn's final artifact audit is
complete.  If compute remains and a second PDE is scientifically useful, a
separate addendum must freeze its equation, targets, candidate budgets, banks,
and compute cap before generating its new banks or running any solver.  The
choice must use only evidence predating this study and cannot be selected to
make the Allen-Cahn outcome look stronger.

