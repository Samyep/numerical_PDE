# Current-round experiments

## Conservative HJB bound

| d | heuristic | exact A | uniform sqrt(15) | box |
|---:|---:|---:|---:|---:|
| 100 | 1.527 | 0.792 | **0.830** | 1.052 |
| 120 | 1.564 | 0.812 | **0.851** | 1.050 |
| 140 | 1.568 | 0.828 | **0.864** | 1.048 |
| 160 | 1.549 | 0.849 | **0.879** | 1.044 |

The uniform radius uses only the public coefficient ranges and remains close to the exact-matrix projector.

## Higher-level 100D funding

| n | M | baseline MAE | hard MAE | box MAE |
|---:|---:|---:|---:|---:|
| 4 | 2 | 6.483 | **1.341** | 4.791 |
| 4 | 3 | 3.025 | **0.877** | 2.708 |
| 5 | 2 | 5.196 | **1.459** | 4.592 |

## Counterparty-credit-risk negative control

| n | M | baseline MAE | hard MAE | hard violation rate |
|---:|---:|---:|---:|---:|
| 2 | 6 | 0.721 | 0.721 | 0.00000 |
| 2 | 10 | 0.477 | 0.477 | 0.00000 |
| 3 | 4 | 0.595 | 0.595 | 0.00024 |
| 3 | 6 | 0.288 | 0.288 | 0.00020 |
| 4 | 3 | 0.428 | 0.428 | 0.00000 |

Hard projection is effectively inactive here, giving a clean negative control.
