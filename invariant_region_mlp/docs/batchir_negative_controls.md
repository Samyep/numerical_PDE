# Batch-IR negative-control check

All methods use paired random seeds. Corrections are applied to recursive child states before the nonlinear generator when such a generator exists. The purpose is to test whether Batch-IR creates an artificial benefit or degradation when the proposed mechanism is absent.

## Allen--Cahn 100D

Certified interval: `0 <= u <= 1`. 100 repetitions per budget, n=3.

| M | Raw rel-MAE | Samplewise rel-MAE | Batch-IR rel-MAE | Violation | Batch activation |
|---:|---:|---:|---:|---:|---:|
| 2 | 3.3035% | 3.3035% | 3.3035% | 0.0000% | 0.0000% |
| 3 | 1.9590% | 1.9590% | 1.9590% | 0.0000% | 0.0000% |
| 4 | 1.4215% | 1.4215% | 1.4215% | 0.0000% | 0.0000% |

The three methods are identical realization-by-realization. The constraint is never violated, hence Batch-IR has alpha=1 for every tested recursive batch.

## Counterparty credit risk 100D

Standard benchmark parameters: T=2, sigma=.2, beta=.03, K1=30, K2=60, L=15; reference value 2.626. The certified value interval used here is [-15,15]. 100 paired repetitions per budget.

| n | M | Raw MAE | Samplewise MAE | Batch-IR MAE | Violation | Batch activation |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 10 | 0.450018 | 0.450018 | 0.450018 | 0.00000% | 0.00000% |
| 3 | 6 | 0.297810 | 0.297810 | 0.297810 | 0.00000% | 0.00000% |
| 4 | 3 | 0.516940 | 0.516940 | 0.516940 | 0.00000% | 0.00000% |

For these paired seeds no recursive child state left the certified interval, so Samplewise IR and Batch-IR are exactly inert.

## Linear convection--diffusion 100D

Here f=0, so the constrained gradient z is not reused by a nonlinear generator. We deliberately estimate z with finite Monte Carlo noise and apply the gradient constraint, but the value estimator is structurally independent of that correction.

| M | Raw MAE | Samplewise MAE | Batch-IR MAE | z violation | Batch activation | Max value difference |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | 0.114882 | 0.114882 | 0.114882 | 100.00% | 100.00% | 0.0e+00 |
| 8 | 0.050744 | 0.050744 | 0.050744 | 97.00% | 97.00% | 0.0e+00 |
| 16 | 0.026660 | 0.026660 | 0.026660 | 8.50% | 8.50% | 0.0e+00 |

This is the stronger negative control: the gradient constraint can be active, but because z is not fed into f, neither Samplewise IR nor Batch-IR changes the value estimate.

## Verdict

- inactive constraint -> alpha=1 and exact overlap (Allen--Cahn, credit risk);
- active but generator-irrelevant gradient correction -> exact overlap in value (linear convection--diffusion).

Thus the positive Batch-IR gains in HJB/funding are not explained by a generic effect of common-factor shrinkage on every MLP computation.
