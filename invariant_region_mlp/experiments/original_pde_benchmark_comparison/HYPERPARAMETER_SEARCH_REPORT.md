# Tangent-Actor hyperparameter search log

Actor tuning uses no analytic value or gradient labels and never opens the
final-test reference.  Each candidate is ranked by held-out tangent-label
error and a common reference-free linear-policy value bank.  Positive
`source_minus_new` means the new Actor lowers the tangent-policy value.

## Stage B candidates

| Trial | Width x depth | Activation | Points | Paths/label | Label steps | Optimizer steps | Initial status |
|---|---:|---|---:|---:|---:|---:|---|
| `pilot192` | 192 x 3 | SiLU | 2,048 | 64 | 16 | 2,000 | Completed |
| `wide384` | 384 x 4 | SiLU | 4,096 | 64 | 24 | 3,000 | Completed |
| `deep_elu` | 256 x 5 | ELU | 4,096 | 64 | 24 | 3,000 | Completed |
| `quality128` | 256 x 4 | SiLU | 4,096 | 128 | 32 | 4,000 | Completed |

All four receive one update at seed 0.  The best two reference-free candidates
continue before Stage C.  This is a successive-halving study, not a Cartesian
grid.

The frozen Stage B decisions were:

| Benchmark | Continued trials | Best one-update policy value | Held-out / replication RMSE for best |
|---|---|---:|---:|
| HJB-Log | `deep_elu`, `quality128` | `4.5755605` (`deep_elu`) | `0.011321 / 0.000741` |
| OU/GMM HJB | `wide384`, `quality128` | `309.81450` (`wide384`) | `0.284557 / 0.019439` |
| Allen–Cahn | `deep_elu`, `quality128` | `0.046354659` (`deep_elu`) | `0.003222 / 0.001038` |

These policy values are reference-free linear-policy estimates on a common
bank.  They are comparable within a benchmark and have no cross-benchmark
meaning.  In HJB-Log and OU/GMM, approximation error against held-out tangent
labels is much larger than independent label-replication noise.  In
Allen–Cahn, all four one-update policy values differ by less than `4e-7`, so
architecture ranking is practically tied at this stage.

## Stage C candidates

`data8192` isolates more spatial points at 64 paths/label.  `focused` combines
8,192 points, 128 paths, 32 label steps, a 384 x 4 SiLU Actor, and 5,000 fitting
steps.  `focused_plain` is the identical plain-FK-label ablation.  Later Actor
rounds warm-start only when architecture is unchanged.

The additional architecture controls are `compact128`, `wide512`, and
`tanh256`.  A conditional `max_quality` trial uses 16,384 points, 256
paths/label, 64 label steps, and 6,000 optimizer steps.  It receives one update
only when the staged run is stable and either still improves the replicated
policy value or leaves held-out Actor error above five times label-replication
RMSE.  Further updates require a new improvement on the replicated common
policy bank.  This gate prevents a nominal upper-budget option from becoming
an automatic compute sweep.

Plain-FK labels are retained as a reported training ablation.  They are not
eligible to replace the declared MCV-training method in the primary
reference-free selection.  All candidates may appear later in a descriptive
postselection test-error plot, clearly marked as not used for selection.

## HJB-Log frozen search outcome

The three-replicate, 512-path, 64-step common-bank ranking selected
`focused a4` after `993.72` seconds of cumulative training.  Its policy value
was `4.575605754`, MC standard error `2.00e-5`, held-out gradient-label RMSE
`0.006573`, and independent label-replication RMSE `0.000847`.  Relative to
the selected checkpoint:

| Candidate | Policy minus selected | Paired state SE | Interpretation |
|---|---:|---:|---|
| `focused a3` | `+6.70e-5` | `2.91e-5` | the fourth update gives a small replicated improvement |
| `focused_plain a3` | `+1.03e-4` | `6.99e-5` | plain-FK training is competitive but does not beat MCV training |
| `max_quality a1` | `+3.07e-4` | `1.01e-4` | twice the states and paths plus 64 steps is detectably worse |

`max_quality a1` took `1977.89` seconds, used its 5,500-step checkpoint, and
reduced label-replication RMSE to `0.000426`; nevertheless, its held-out RMSE
remained `0.009486` and its policy value worsened.  The conditional continuation
gate therefore stopped it before `a2`.  This result separates cleaner labels
from better policies: supervision noise fell, but the chosen policy did not
improve.

## OU/GMM HJB staged search outcome

The replicated common-bank ranking selected `focused a4` after `1539.59`
seconds of cumulative training.  Its policy value was `308.7257452`, MC
standard error `0.002914`, held-out gradient-label RMSE `0.059290`, and
independent label-replication RMSE `0.009644`.  All exponential weights were
finite, their mean ESS fraction was one, and MCV reduced inference-estimator
variance by a factor of `43.57` on this ranking bank.

| Candidate | Policy minus selected | Paired state SE | Interpretation |
|---|---:|---:|---|
| `focused a3` | `+0.027388` | `0.003114` | the fourth update gives a clearly resolved improvement |
| `focused_plain a3` | `+0.028678` | `0.004092` | matched three-round plain training is close to MCV training at `a3` |
| `data8192 a3` | `+0.097540` | `0.010293` | more states alone does not match the combined focused budget |
| `wide384 a3` | `+0.106273` | `0.011713` | the smaller Stage-B budget plateaus earlier |

At the matched third update, MCV and plain-FK training produced nearly equal
policy values (`308.75313` and `308.75442`).  MCV did make the supervision
substantially cleaner: label-replication RMSE was `0.010465`, versus `0.029624`
for plain FK, and held-out label RMSE was `0.071546`, versus `0.077291`.
Thus the DPI evidence supports a training-variance benefit but does not by
itself show a material matched-round policy gain from MCV.  The fourth MCV
update is the material improvement.

Because `a4` remained better by about `8.8` paired standard errors and its
held-out error was still `6.15` times the label-replication RMSE, the
predeclared gate admits one `max_quality` update.  Further upper-budget rounds
remain conditional on that checkpoint improving the same replicated bank.

