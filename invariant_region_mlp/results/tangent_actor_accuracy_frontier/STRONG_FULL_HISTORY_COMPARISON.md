# Strong Full-History Comparison

The baseline is the raw recursive Full-History MLP estimator. It has no Tangent-Actor, MCV, terminal-aware learned correction, neural surrogate, value projection, or oracle active subspace. Large terminal Monte Carlo tensors are accumulated in blocks to fit device memory; this changes storage and summation order only, not the estimator or sample count.

| family | n | M | role | finite seeds | Value RMSE mean ± sd | wall seconds mean | peak GPU GB |
|---|---:|---:|---|---:|---:|---:|---:|
| coupled | 3 | 256 | frontier | 1/1 | 0.00305039 ± 0 | 82.0 | 1.04 |
| coupled | 3 | 384 | frontier | 1/1 | 0.0038804 ± 0 | 274.1 | 1.19 |
| coupled | 3 | 512 | endpoint | 3/3 | 0.00431753 ± 0.00011 | 671.1 | 1.31 |
| coupled | 4 | 16 | frontier | 1/1 | 0.02219 ± 0 | 8.0 | 0.15 |
| coupled | 4 | 32 | frontier | 1/1 | 0.00665695 ± 0 | 15.9 | 1.03 |
| coupled | 4 | 64 | frontier | 1/1 | 0.00340733 ± 0 | 175.2 | 1.32 |
| coupled | 4 | 96 | endpoint | 3/3 | 0.00163057 ± 7.1e-05 | 886.5 | 3.05 |
| primary | 3 | 384 | endpoint | 3/3 | 0.000534724 ± 1.8e-05 | 271.3 | 1.19 |
| primary | 3 | 512 | endpoint | 3/3 | 0.000579189 ± 1e-05 | 668.8 | 1.31 |
| primary | 4 | 2 | profile | 1/1 | 1.5157 ± 0 | 0.7 | 0.03 |
| primary | 4 | 16 | frontier | 1/1 | 0.00242686 ± 0 | 7.4 | 0.15 |
| primary | 4 | 32 | frontier | 1/1 | 0.000885548 ± 0 | 15.0 | 1.03 |
| primary | 4 | 64 | frontier | 1/1 | 0.000366765 ± 0 | 173.8 | 1.32 |
| primary | 4 | 96 | endpoint | 3/3 | 0.000219475 ± 7.6e-06 | 876.4 | 3.05 |

## Fair conclusion

The strongest completed primary endpoint was `n=4, M=96` with mean RMSE `0.000219475` across 3 seeds (small-sample 95% t interval `[0.00020048, 0.00023847]`). Configuration grids were frozen from resource feasibility and algorithmic budget, and failed or worsening rows remain in the table.
Runtime is kept visible, while the primary comparison concerns maximum attainable numerical accuracy. Historical rows were recomputed against the refined reference whenever stored predictions were available.
