# Entropy/OOD benchmark design and implementation record

This is an **exploratory design record**, not an independently timestamped
preregistration.  Four of the five Euler initial conditions had appeared in an
earlier local benchmark; this run extends their horizon and adds a transonic
case.  The detailed SWE split was recorded before its formal model training,
but the repository had no immutable commit at that point.  Consequently these
results must not be described as a blind confirmatory test.

## Scientific question

Does the 64-cell HCFL solver retain accuracy while enforcing its stated
conservation, positivity, interface-entropy, and fully-discrete entropy
safeguards on wave patterns and rollout horizons outside its training task?

The purpose of the suite is to stress the claimed mechanism.  It is not a
claim that every competing method is universally inferior.

## Euler suite

All methods start from the same 64 finite-volume cell averages.  Fine
references start from the piecewise-constant prolongation of those same
averages, are evolved by SSP-RK2/HLLC on 2048 cells, and are conservatively
restricted to 64 cells.  The domain is periodic for this matched first phase.
The final time is 63 saved intervals, four times the 15-interval training
horizon.

The named cases are fixed as follows (primitive variables are `(rho, u, p)`):

| case | left state | right state | role |
|---|---|---|---|
| Sod | `(1, 0, 1)` | `(0.125, 0, 0.1)` | ordinary control |
| transonic rarefaction | `(0.1, -2, 0.1)` | `(1, -1, 1)` | entropy-solution stress |
| near-vacuum expansion | `(1, -2, 0.4)` | `(1, 2, 0.4)` | positivity/entropy stress |
| collision | `(1, 2, 1)` | `(1, -2, 1)` | compressive-wave stress |
| strong pressure | `(1, 0, 5)` | `(1, 0, 0.05)` | strong-shock stress |

The main comparison uses:

1. HLLC-64;
2. TVD MUSCL-HLLC-64;
3. the previously trained vanilla residual FNO-64;
4. the three previously trained HCFL-64 seeds;
5. the same-data RoeNet architecture adaptation;
6. the authors' official data-free PINN checkpoints on their native Sod and
   Lax tasks, reported separately rather than pooled into this five-case mean.

RoeNet and the PINN are not to be called exact paper reproductions unless
their original task, data, optimizer, and architecture are unchanged.  Any
adaptation is labelled in every table.

## SWE suite

The secondary comparison uses the official clawFNO implementation on the 2D
radial dam-break task, alongside vanilla FNO and a 2D HCFL solver.  Tests are
split into:

- in-distribution radius and height ratio;
- radius outside the training interval;
- stronger, strictly positive height ratios;
- a doubled rollout horizon.

Dry-bed and bathymetric cases are excluded because the current HCFL solver has
neither a wet/dry treatment nor a well-balanced bathymetry discretization.

The implemented suite uses the following details.  The square
domain is `[-2.5, 2.5]^2`, gravity is one, the deployed grid is `32 x 32`, and
the independent reference grid is `128 x 128`.  The reference is a
second-order MC-limited HLL finite-volume solver with SSP-RK2 time stepping and
the same constant-extrapolation boundary condition as the PDEBench generator.
It is not represented as a byte-for-byte PyClaw reproduction.  Fine states are
conservatively averaged by `4 x 4` blocks before any method sees them.

The data split is 100 training trajectories, 20 validation trajectories,
and disjoint test sets of 20 in-distribution radii, 20 radii outside the
training interval, 20 stronger positive height ratios, and 10 doubled-horizon
trajectories.  Training and ordinary testing use 25 saved states from
`t = 0` through `t = 0.96`; the doubled-horizon set uses 49 saved states
through `t = 1.92`.  Training radii are uniform on `[0.3, 0.7]`; radius OOD is
the equal mixture `[0.15, 0.25] U [0.75, 0.85]`; strong-height OOD uses an
inner height uniform on `[2.5, 3.0]` instead of the training value two.  The
outer height is one and the initial velocity is zero throughout.

The vanilla FNO and clawFNO use the authors' source architecture and published
radial-dam hyperparameters (`modes=8`, `width=20`, one-shot 24-frame output).
Because the locally generated data and split differ from the downloadable
author dataset, both are labelled **official-architecture adaptations**, not
exact reproductions of published error numbers.  Both operator models and
HCFL receive the same 100 trajectories.  Hyperparameter selection uses only
the 20 validation trajectories.  A second 24-frame operator call, initialized
from the first block's final prediction, is used for the doubled-horizon test.

### Actual SWE training configuration

The HCFL comparison uses four- and six-cell stencils, width 72, seed zero,
Adam with learning rate `5e-4`, batch size four, validation every 500 updates,
an 8,000-update arm cap, proposal-feasibility weight `1e-3`, and a numerical
fully-discrete entropy tolerance of `5e-7`.  Checkpoint selection is
lexicographic: physical validation completion first, completed-only NRMSE
second.  The six-cell arm is selected using validation only.

FNO and clawFNO use seed zero, batch size ten, Adam learning rate `1e-2`,
weight decay `1e-4`, at most 500 epochs, and patience 100.  Their checkpoint
selection likewise ranks validation completion before completed-only error; a
checkpoint with no physical validation trajectory is retained only as an
explicitly inadmissible failure diagnostic.

Dataset RNG seeds are `61001`, `62001`, `63001`, `64001`, `65001`, and `66001`
for train, validation, ID test, radius-OOD test, height-OOD test, and long test,
respectively.  These values document the run; they are not evidence of an
externally preregistered protocol.

## Metrics and failure accounting

The analysis uses the following evidence order:

1. physical completion rate (finite state, positive density/pressure/depth);
2. global discrete entropy-balance violations;
3. conservation drift;
4. rollout and final-time errors on completed trajectories;
5. total variation, range overshoot, and wave-location diagnostics;
6. offline training and online rollout cost, reported separately.

Interface Tadmor residuals are reported only for flux-form methods for which
that quantity is defined.  A failed trajectory is never removed from the
completion denominator.  Conditional error and completion rate are reported
separately; no arbitrary finite error is assigned to NaNs.

After the predictions and checkpoints were frozen, normalized mean absolute
error (NMAE) was added as a second accuracy view at the user's request.  This
post-hoc metric did not select or retrain any model.  It uses the same
training-channel standardization and the same completed trajectories as
NRMSE:

`NMAE = mean(abs((prediction - reference) / training_channel_std))`.

Euler and the native PINN table use the conserved channels `(rho, rho*u, E)`,
matching their reported NRMSE; 2-D SWE uses primitive channels `(h, u, v)`,
again matching its NRMSE.  The machine-readable tables additionally retain
unnormalized per-primitive-channel MAE so that the aggregate cannot hide which
physical variable dominates.  Model selection remains based on the original
completion-first, NRMSE-second rule.

## Tuning history rule

The first run evaluated existing 64-cell checkpoints.  Subsequent tuning used
validation rather than the named result tables.  The maximum training cap was
50,000 optimizer updates.  A higher grid was to be attempted only if 64-cell
HCFL was materially deficient while satisfying the safety audit; that trigger
was not met.  This history is retained for transparency, but it does not turn
the current suite into a confirmatory benchmark.
