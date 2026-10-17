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

