# 1D SWE Phase 2: Fully Discrete Safety and Multi-Step Training

This checkpoint extends the first 1D shallow-water HCFL pilot with fully discrete safety layers, multi-step training, a MUSCL baseline, and a near-dry stress test.

## Safety construction

Starting from an entropy-feasible learned candidate F_H and a low-order Rusanov flux F_L, use

F_{i+1/2} = F^L_{i+1/2} + alpha_{i+1/2}(F^H_{i+1/2}-F^L_{i+1/2}),  0 <= alpha <= 1.

For each cell, the low-order update h_i^L is assumed positive under an explicitly checked CFL condition. A cellwise budget limits the sum of the two adverse height corrections; the shared interface coefficient is the minimum of the two neighboring cell factors. This guarantees h_i^{n+1} >= h_floor.

Because Tadmor's interface condition is affine in the flux, and both Rusanov and the hard-projected learned flux are entropy feasible, the positivity blend remains inside the same entropy half-space.

A second per-trajectory scalar blend

F(beta) = F_L + beta(F_pos-F_L),  0 <= beta <= 1,

is used when necessary to enforce non-increase of the periodic-domain total physical entropy after the explicit step. The update is affine in beta and the shallow-water energy is convex for h>0, so a bisection finds the largest feasible beta. Entropy is accumulated in float64 and a small endpoint safety margin is used.

Adaptive CFL substepping based on |u|+sqrt(g h) is required in near-dry OOD rollouts so that the low-order Rusanov endpoint remains feasible.

## Verified stress-test behavior

On 4,000 random 64-cell states with h in [0.03,2.8] and u in [-2,2], entropy-only HCFL produced a minimum next-step h of -0.0875. The positivity blend raised the minimum to 1e-4. Only about 0.021% of interfaces were limited. The global entropy limiter was active on about 0.25% of trajectories, with minimum beta about 0.993. After the refined float64 bisection, no total-entropy increases above 1e-8 were observed in this test.

## Four-step trajectory fine-tuning

Three-seed rollout NRMSE:

| split | one-step HCFL-safe | four-step HCFL-safe |
|---|---:|---:|
| ID | 0.04741 +/- 0.00132 | **0.04542 +/- 0.00095** |
| OOD | 0.11456 +/- 0.00992 | **0.11038 +/- 0.00909** |

Thus multi-step fine-tuning gives a modest but consistent improvement.

## MUSCL baseline

On the original non-near-dry tests:

| split | HCFL-safe | MUSCL-Rusanov |
|---|---:|---:|
| ID | **0.04738 +/- 0.00184** | 0.06663 +/- 0.00396 |
| OOD | **0.11443 +/- 0.00990** | 0.14061 +/- 0.01114 |

## Near-dry stress test

Reference trajectories include depths as low as 0.05. With adaptive CFL substepping, the fully safe HCFL rollout maintains h >= 1e-4.

Three-seed rollout NRMSE:

| method | NRMSE |
|---|---:|
| HCFL-safe one-step | 0.18591 +/- 0.03041 |
| **HCFL-safe four-step** | **0.17712 +/- 0.02587** |
| MUSCL-Rusanov | 0.20718 +/- 0.01813 |
| HLL | 0.26825 +/- 0.02765 |
| Rusanov | 0.29097 +/- 0.03004 |

For the fine-tuned HCFL model, the positivity limiter is active on about 0.17% of interfaces and the global entropy line search on about 2.5% of substeps. The stress test uses roughly 2.6 adaptive substeps per saved-data interval on average.

## Runtime

Approximate CPU cost for 100 simultaneous 64-cell trajectories:

| method | ms / saved-data step |
|---|---:|
| Rusanov | 0.23 |
| HLL | 0.38 |
| MUSCL-Rusanov | 0.86 |
| HCFL entropy-only | 2.97 |
| HCFL fully-safe | 3.87 |

These timings are not optimized.

## Precise current guarantees / limitations

Supported in the current 1D flat-bottom SWE implementation:
1. exact finite-volume conservation from shared interface fluxes;
2. hard Tadmor interface half-space feasibility to numerical precision;
3. fully discrete h >= h_floor under the checked low-order CFL premise;
4. global periodic-domain total entropy non-increase to bisection/numerical precision;
5. multi-step trajectory training improves rollout accuracy in all three tested seeds.

Limitations:
- the fully discrete entropy condition is global, not a local cellwise entropy-flux inequality;
- the near-dry solver uses a positive floor (1e-4), not a true wet/dry front treatment;
- bathymetry / well-balancing are not included;
- runtime is several times higher than the classical baselines in the current Python code;
- Euler is not yet implemented.
