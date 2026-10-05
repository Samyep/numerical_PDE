# Extended Batch-IR experiment suite

Priority order requested: **1 -> 2 -> 3 -> 6 -> 4**.

All comparisons use paired stochastic budgets within each configuration. Batch-IR uses one common feasible scale for a Monte Carlo sibling batch before the nonlinear generator is evaluated. Final root outputs are not clipped.

## 1. Headline 100--160D Rosenbrock HJB

Protocol: corrected SCaSML-scale standalone MLP reproduction, (n=2), (M=10), 1,000 interior unit-ball points + 200 unit-sphere points, 10 paired repetitions, fixed Rosenbrock matrix per dimension, and the non-oracle certified radius (|z|_2lesqrt{15}). The child-gradient estimator uses the corrected centered Elworthy--Bismut--Li form. The reference uses scaled Gauss--Laguerre quadrature with log-sum-exp stabilization.

| d | Raw rel-L2 | Samplewise IR | Batch-IR | Batch vs sample reduction | Batch wins |
|---:|---:|---:|---:|---:|---:|
| 100 | 2.27995 | 0.62077 | **0.26530** | 57.3% | 10/10 |
| 120 | 2.66715 | 0.62851 | **0.25586** | 59.3% | 10/10 |
| 140 | 3.09462 | 0.62791 | **0.24556** | 60.9% | 10/10 |
| 160 | 3.47939 | 0.64044 | **0.24287** | 62.1% | 10/10 |

Paired t-test p-values for Batch-IR versus Samplewise IR are (4.9	imes10^{-20}), (8.8	imes10^{-19}), (3.8	imes10^{-20}), and (1.9	imes10^{-19}), respectively.

The Samplewise values reproduce the magnitude of the archived corrected (sqrt{15})-ball results (about 0.63--0.64). The Raw column here means **no correction** and should not be conflated with the older heuristic-clipping column.

Batch-IR is aggressive in this regime: its sibling-batch activation rate is 99.87--99.98%, and its mean common scale decreases from about 0.330 at (d=100) to 0.260 at (d=160). This strengthens the finite-budget evidence while also emphasizing the finite-sample bias/coupling caveat.

## 6. Batchwise mechanism diagnostics on the same HJB protocol

### Generator MSE

| d | Raw | Samplewise IR | Batch-IR |
|---:|---:|---:|---:|
| 100 | 4427.03 | 38.42 | **10.50** |
| 120 | 9992.13 | 43.17 | **10.98** |
| 140 | 12730.85 | 42.62 | **10.81** |
| 160 | 11312.72 | 42.99 | **10.30** |

Batch-IR reduces generator MSE versus raw by about **99.76--99.92%**.

### Variance across paired runs of the actual nonlinear level correction

| d | Raw | Samplewise IR | Batch-IR |
|---:|---:|---:|---:|
| 100 | 20.1419 | 0.2221 | **0.1403** |
| 120 | 29.9815 | 0.2156 | **0.1393** |
| 140 | 45.2710 | 0.2058 | **0.1409** |
| 160 | 60.7858 | 0.2109 | **0.1526** |

Batch-IR reduces the actual level-correction variance versus raw by about **99.30--99.75%** and is lower than Samplewise IR in all four dimensions.

The variance of *individual* nonlinear increments is not uniformly lower than Samplewise IR: Batch-IR is comparable and is slightly higher in some dimensions. The stronger empirical statement is therefore about generator MSE and the variance of the actual averaged level correction, not every possible variance statistic.

## 2. Deeper 100D nonlinear-funding recursion

Reference value: 21.299. The state-dependent certified ellipsoid is enforced in normalized delta coordinates.

| n | M | Raw MAE | Samplewise IR | Batch-IR | Batch vs sample reduction | Batch wins | paired p |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 10 | 1.52683 | 0.52147 | **0.45658** | 12.4% | 73% | 4.7e-6 |
| 3 | 8 | 1.10714 | 0.38394 | **0.30012** | 21.8% | 75% | 3.2e-10 |
| 4 | 3 | 3.20564 | 0.77079 | **0.68970** | 10.5% | 60% | 1.4e-3 |

Thus the Batch-IR advantage survives deeper recursion; it is not confined to the earlier shallow (n=2) sweep.

## 3. Neufeld--Wu at M=2

Thirty repetitions per dimension, (n=M=2).

| d | Raw MAE | Samplewise IR | Batch-IR |
|---:|---:|---:|---:|
| 100 | 0.06066 | **0.02131** | **0.02131** |
| 200 | 0.03878 | **0.01682** | **0.01682** |
| 300 | 0.04039 | **0.01737** | **0.01737** |

Samplewise IR and Batch-IR are exactly tied in value error, just as in the previous (n=M=3) experiment. Either correction moves the offending gradient below the threshold where the driver activates, so both remove the same spurious nonlinear term. A strict “win fraction” is therefore not meaningful here; these are ties.

## 4. Batchwise structural sanity check on the dimensionality counterexample

The batchwise structural correction is
[
z_i mapsto alpha_B Q z_i,
]
where (Q) is the certified projection onto (operatorname{span}(e_1)). Pure radial scaling without (Q) is **not** an exact rescue.

At (n=m=3), 512 repetitions:

| d | Raw value RMSE | Samplewise subspace | Batch structural | max |Batch-Sample| |
|---:|---:|---:|---:|---:|
| 10 | 3.59502 | **0.110083** | **0.110083** | 0 |
| 100 | 39.12143 | **0.110083** | **0.110083** | 0 |
| 1000 | 396.59147 | **0.110083** | **0.110083** | 0 |

The two corrected value estimators are identical realization-by-realization. The common scale is active in about 59% of the tested sibling groups, so equality is not caused by an inactive Batch-IR step; it follows from returning every gradient to the generator-compatible certified subspace, where the offending generator vanishes.

## Overall verdict

The extended suite strengthens the empirical picture:

- on the headline HJB protocol, Batch-IR beats Samplewise IR in all paired repetitions at all four dimensions;
- the gain tracks a much smaller generator MSE and a smaller variance of the actual nonlinear level correction;
- the benefit survives deeper funding recursion;
- on a thresholded driver, Batch-IR and Samplewise IR tie when both remove the same spurious activation;
- on the structural counterexample, adding a common batch scale after the certified subspace map preserves exact rescue.

These experiments do **not** establish that Batch-IR is universally superior. The common scale couples samples and moves the empirical mean, and in the HJB/funding regimes it activates almost always. Sharp consistency and complexity theory for batch-coupled retractions remains open.

## Reproducibility files

- `experiments/batch_contraction/hjb_batchwise_headline.py`
- `experiments/batch_contraction/funding_batchwise_deep.py`
- `experiments/batch_contraction/neufeld_batch_m2.py`
- `experiments/batch_contraction/counterexample_batchwise_sanity.py`
- `results/batch_contraction/hjb_batchwise_headline_results.json`
- `results/batch_contraction/batchwise_extended_suite_summary.json`
