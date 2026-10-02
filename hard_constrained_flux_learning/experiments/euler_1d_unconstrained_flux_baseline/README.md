# Unconstrained DNN flux baseline

## Pre-run hypothesis

> With the direct HLLC + five-point MLP flux proposal held fixed, removing all
> hard output post-processing may slightly improve ordinary validation error,
> but it will increase Tadmor violations and can lose physical admissibility on
> severe shock-tube rollouts.

This experiment revisits the historical plain-flux control under the current
validation-convergence protocol. The earlier 1,100-update plain/hard results do
not establish the converged comparison.

## Single changed factor

The unconstrained arm emits the existing `DirectFlux` proposal directly as the
shared interface flux. It does **not** apply:

- the hard Tadmor projection;
- the conservative local density/pressure limiter;
- the fully-discrete total-entropy line search.

The proposal itself is deliberately unchanged for a one-factor ablation: it is
the same HLLC-centered five-point MLP correction used by `direct_broad`.
Conservation remains exact because neighboring cells use the same interface
flux. “Unconstrained” here means no physics-feasibility post-processing of the
proposal; it does not mean changing the proposal architecture at the same time.

## Held fixed

- seed-0 broad 580-trajectory training tensor;
- independent 136-trajectory validation tensor;
- 512-cell strict Rusanov + SSP-RK2 teacher restricted to 64 cells;
- width 72, initialization seed, minibatch sequence, Adam, gradient clipping,
  and state normalization;
- validation every 100 updates, learning-rate schedule, plateau rule, and
  20,000-update cap;
- ordinary ID, broad-random in-support, moderate high-frequency OOD, and five
  canonical final tests;
- periodic boundary treatment.

Raw rollouts are never repaired. A trajectory is marked failed at the first
nonfinite state, density below `1e-5`, or pressure below `1e-5`; it is not
silently replaced by a classical step. Checkpoint selection gives every fully
admissible validation rollout priority over a checkpoint with any failed
validation trajectory, then minimizes rollout NRMSE.

## Run

```powershell
python run_unconstrained_baseline.py --seed 0
python plot_results.py --seed 0
python plot_unconstrained512_vs_fvm.py --seed 0
```

The 512-cell comparison deploys the validation-selected 64-cell checkpoint
with eight raw learned updates per saved interval, matching the training
`dt/dx`. It plots native FVM-512 on the same cell centers and keeps native
FVM-2048 as a light background reference. No curve is averaged or sampled,
and a learned curve is omitted if its raw rollout becomes inadmissible before
the final time.
