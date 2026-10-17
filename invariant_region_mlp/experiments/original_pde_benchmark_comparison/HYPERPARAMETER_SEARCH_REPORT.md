# Tangent-Actor hyperparameter search log

Actor tuning uses no analytic value or gradient labels and never opens the
final-test reference.  Each candidate is ranked by held-out tangent-label
error and a common reference-free linear-policy value bank.  Positive
`source_minus_new` means the new Actor lowers the tangent-policy value.

## Stage B candidates

| Trial | Width x depth | Activation | Points | Paths/label | Label steps | Optimizer steps | Initial status |
|---|---:|---|---:|---:|---:|---:|---|
| `pilot192` | 192 x 3 | SiLU | 2,048 | 64 | 16 | 2,000 | Scheduled |
| `wide384` | 384 x 4 | SiLU | 4,096 | 64 | 24 | 3,000 | Scheduled |
| `deep_elu` | 256 x 5 | ELU | 4,096 | 64 | 24 | 3,000 | Scheduled |
| `quality128` | 256 x 4 | SiLU | 4,096 | 128 | 32 | 4,000 | Scheduled |

All four receive one update at seed 0.  The best two reference-free candidates
continue before Stage C.  This is a successive-halving study, not a Cartesian
grid.

## Stage C candidates

`data8192` isolates more spatial points at 64 paths/label.  `focused` combines
8,192 points, 128 paths, 32 label steps, a 384 x 4 SiLU Actor, and 5,000 fitting
steps.  `focused_plain` is the identical plain-FK-label ablation.  Later Actor
rounds warm-start only when architecture is unchanged.

The tables in this file will be replaced with measured selection decisions,
failed trials, and cumulative costs after the staged jobs finish.

