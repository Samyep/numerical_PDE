# Seed-0 result: direct 512-cell SWE training

## Question

The 64-cell-trained closure was accurate after deployment on 512 cells, but
showed large shock-local ringing.  This audit tests the hypothesis that the
ringing is primarily a resolution-transfer error: the closure learned the
amount of antidiffusion needed by a 64-cell operator and applied that same
correction to an already less-diffusive 512-cell operator.

## Controlled experiment

Only one network was retrained:

- central flux + nonnegative Roe dissipation + proposal-feasibility loss;
- symmetric 4-cell stencil;
- 6,050 parameters;
- 580 training trajectories and 136 independently seeded validation
  trajectories;
- HLL-2048 labels conservatively restricted to 512 cells;
- eight differentiable substeps per `5e-4` supervision interval, preserving
  the original stable ratio `dt/dx = 0.032`;
- hard Tadmor projection in every substep;
- a 50,000-update cap, with selection and stopping controlled only by
  independent validation rollout error.

The best checkpoint occurred at update 30,600.  Training stopped at update
33,000 after a validation plateau at the minimum learning rate `3e-6`.
Best validation rollout NRMSE was `0.00461806`; all validation depths remained
positive.

GPU-float32 HLL-2048 label generation was checked against an independent
CPU-float64 implementation on two trajectories.  The normalized RMS
difference after restriction was `2.58e-6`, far below the learned-model error
scale.

## Held-out 512-cell deployment

All methods use the same hard projection, positivity safety, fully discrete
entropy safety, test initial conditions, and HLL-2048 scoring reference.
`Roe-512` is a zero-network-output control (`m=1`), not a learned checkpoint.

| boundary | method | mean rollout NRMSE | positive final TV excess | excess extrema |
|---|---:|---:|---:|---:|
| periodic | HLL-512 | 0.023165 | 0.09% | 0 |
| periodic | classical Roe-512 | 0.015394 | 0.05% | 0 |
| periodic | HCFL trained at 64 | 0.011638 | 40.03% | 129 |
| periodic | **HCFL trained at 512** | **0.005618** | **0.88%** | **7** |
| transmissive | HLL-512 | 0.010140 | 0.06% | 0 |
| transmissive | classical Roe-512 | 0.006503 | 0.03% | 0 |
| transmissive | HCFL trained at 64 | 0.005775 | 39.75% | 48 |
| transmissive | **HCFL trained at 512** | **0.002775** | **1.55%** | **3** |

Relative to the identical 64-trained architecture, direct 512 training:

- lowers mean rollout NRMSE by 51.7% on periodic cases and 51.9% on
  transmissive cases;
- lowers positive final TV excess by 97.8% and 96.1%, respectively;
- lowers excess-extrema counts from 129 to 7 and from 48 to 3.

The improvement is not merely a switch from HLL to an unlearned Roe flux.
Relative to classical Roe-512, the 512-trained network lowers NRMSE by 63.5%
on periodic cases and 57.3% on transmissive cases.  Its mean absolute learned
Roe-multiplier change from one is about 0.48, so the network has not collapsed
to the zero-output state.

## Safety and remaining limitation

- Minimum depth remained positive in every test.
- The local positivity limiter never intervened.
- Maximum post-projection entropy residual was `1.89e-13`.
- Fully discrete entropy intervention was zero on periodic tests and 1.94% on
  transmissive tests; mean transmissive beta was 0.9845.
- The raw proposal violation rate at tolerance `1e-8` was zero.

The resolution-transfer hypothesis is therefore strongly supported for this
single seed and this retained architecture.  It is not proven to be the only
source of ringing: the 512-trained model still has seven periodic and three
transmissive excess extrema, concentrated in collision, supercritical, and
transcritical cases.  A shock-local monotonicity control remains useful, but
the dominant 64-to-512 ringing should not be attributed to the HCFL form
itself.

Primary figures:

- `periodic_and_nonperiodic_train64_vs_train512_seed0.png` (combined)
- `periodic_train64_vs_train512_seed0.png`
- `nonperiodic_train64_vs_train512_seed0.png`
- `stability_accuracy_train64_vs_train512_seed0.png`
- `validation_convergence_train512_seed0.png`
