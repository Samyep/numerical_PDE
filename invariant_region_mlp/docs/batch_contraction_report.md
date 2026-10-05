# Batchwise contraction versus samplewise invariant-region projection

We tested the user's proposal that an entire Monte Carlo batch can be shrunk by a common factor determined by the largest constraint violation. We also tested a mean-preserving version.

## Methods

For a batch z_1,...,z_N and convex ball constraints:
- Raw: no correction.
- Samplewise IR: Euclidean projection of each z_i independently.
- Mean-preserving batch contraction (MPBC): z_i' = zbar + alpha(z_i-zbar), with the largest alpha in [0,1] making the entire batch feasible. If zbar is outside the feasible set, an all-feasible mean-preserving correction is mathematically impossible; the implementation performs no correction.
- Uniform batch shrinkage: z_i' = alpha z_i, where alpha is the largest common factor that makes every sample feasible. This changes the empirical mean but sharply reduces batch variance.

Corrections occur only before the nonlinear generator F; final root outputs are never clipped. All comparisons use paired stochastic paths.

## HJB

For the 20D/40D quadratic HJB with n=2, M=6 and the certified ball ||z||<=sqrt(2):

| d | Raw rel-L2 | Samplewise IR | MPBC | Uniform shrink |
|---:|---:|---:|---:|---:|
| 20 | 314.136 | 0.2412 | 314.136 | **0.0966** |
| 40 | 1034.262 | 0.1993 | 1034.262 | **0.0833** |

MPBC is effectively unavailable: the empirical batch mean itself is almost always outside the certified ball. Uniform shrinkage beats samplewise IR in every paired run.

## 100D nonlinear funding

| M | Raw MAE | Samplewise IR | MPBC | Uniform shrink | MPBC feasible |
|---:|---:|---:|---:|---:|---:|
| 10 | 1.4042 | 0.5102 | 1.4042 | **0.4807** | 0% |
| 20 | 0.8848 | 0.2965 | 0.3350 | **0.2539** | 100% |
| 40 | 0.5221 | 0.1993 | 0.1773 | **0.1486** | 100% |

When the batch mean becomes feasible, MPBC is genuinely useful and at M=40 beats samplewise IR. Uniform shrinkage remains best. Paired uniform-versus-samplewise improvements are significant at M=10,20,40 (p approximately 1.3e-13, 2.9e-35, 8.4e-54).

## Neufeld--Wu 100--300D

For M=n=3:

| d | Raw MAE | Samplewise IR | MPBC | Uniform shrink |
|---:|---:|---:|---:|---:|
| 100 | 0.03029 | **0.00665** | 0.03029 | **0.00665** |
| 200 | 0.01441 | **0.00711** | 0.01441 | **0.00711** |
| 300 | 0.04889 | **0.00679** | 0.04889 | **0.00679** |

Samplewise IR and uniform shrink are identical in value error because either one moves the offending gradients below the generator threshold 25, eliminating the same spurious nonlinear activation. MPBC never contracts an actually harmful violating batch.

## Interpretation

The experiment changes the methodological picture:
1. Variance reduction alone can be highly effective.
2. Mean preservation is not sufficient, and may be impossible in the hardest regimes because the empirical batch mean itself is infeasible.
3. Samplewise Euclidean projection is not empirically optimal among all feasible shrinkage rules.
4. Uniform batch shrinkage is biased at finite sample size because it changes the empirical mean and couples samples through the batch maximum; its strong performance motivates a bias--variance theory rather than an immediate replacement of IR.
5. Samplewise IR retains distinct theoretical advantages: locality, nearest-point geometry, deterministic Fejer improvement, and no cross-sample coupling.

A natural next question is: among constraint-preserving transformations of a stochastic Picard batch, what bias--variance tradeoff is optimal before a nonlinear generator?


## Extended suite (2026-10-05)

We subsequently ran the requested extension in priority order 1 -> 2 -> 3 -> 6 -> 4:
headline 100--160D HJB, deeper funding, Neufeld--Wu at M=2, matched Batch-IR mechanism diagnostics, and a batchwise structural sanity check on the dimensionality counterexample.

The main new result is that, on the corrected headline HJB protocol (n=2, M=10), Batch-IR reduces mean relative-L2 error from Samplewise IR values 0.621--0.640 to 0.243--0.265 across d=100--160, winning all 10/10 paired repetitions at every dimension. On the same states, Batch-IR also has lower generator MSE and lower variance of the actual nonlinear level correction. Deeper funding tests retain a 10.5--21.8% MAE reduction relative to Samplewise IR. Neufeld--Wu remains an exact value-error tie between the two feasible corrections, and the structural counterexample confirms realization-by-realization equality after applying the certified subspace map before the common batch scale.

See `docs/batchwise_extended_suite_report.md` and `results/batch_contraction/batchwise_extended_suite_summary.json` for the complete tables and caveats.
