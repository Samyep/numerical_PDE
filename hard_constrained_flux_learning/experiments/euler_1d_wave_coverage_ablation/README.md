# 1D Euler wave-coverage ablation

## Pre-run hypothesis

> The remaining severe held-out failures of direct-vector HLLC-HCFL are
> primarily caused by insufficient coverage of distinct Euler wave regimes
> (contacts, compressions, expansions, pressure jumps, and collisions), rather
> than insufficient network capacity. Replacing only the 100-trajectory
> random-extreme portion of the current 580-trajectory broad training set with
> a balanced, parameter-randomized wave-regime set should improve several
> canonical held-out tests without materially harming ordinary-ID accuracy.

This hypothesis was recorded before running the experiment.

## Audit motivating the experiment

- The current local baseline is the five-point direct-vector correction on top
  of HLLC, followed by hard Tadmor projection.
- Historical three-seed results identify canonical shock-tube accuracy, not
  entropy feasibility, as the main remaining 1D bottleneck.
- HLLC proposal, local admissibility limiting, global fully-discrete entropy
  limiting, and four-step fine-tuning have already been tested.
- The local ablation runner's split named `OOD` is sampled from the same
  distribution included in broad training. It is therefore reported here as
  `broad_random_in_support`, not as OOD.
- The exploratory reference generator contains an emergency state-floor
  repair. Instrumentation showed zero repair events for the seed-0 baseline
  training set. This experiment uses a strict generator that fails instead of
  repairing a state.
- PyClaw is not installed in the local environment. Both arms therefore use
  the same existing 512-cell Rusanov + SSP-RK2 reference infrastructure.

## The single changed factor

Only the distribution of the final 100 training trajectories changes:

- `broad_random` (baseline): 220 ordinary + 260 broad-random + 100
  random-extreme trajectories;
- `wave_coverage`: the same first 480 trajectories + 100 balanced structured
  wave-regime trajectories.

The structured trajectories are parameter-randomized and do not contain any
exact canonical evaluation initial condition.

## Held fixed

- direct-vector HLLC correction architecture and width;
- parameter initialization and parameter count;
- baseline-derived input and loss normalization for both arms;
- Adam optimizer, learning rate, batch size, gradient clipping, and 1,100
  update steps;
- hard Tadmor projection;
- local conservative density/pressure limiter;
- trajectory-wise fully-discrete entropy safeguard;
- evaluation data, time horizon, coarse/fine grids, and error normalization.

For extreme states, the evaluator recomputes the projection residual in
float64 and applies a machine-precision backoff before returning float32 fluxes.
It likewise validates the fully-discrete entropy line search after float32
flux reconstruction. This fixes two roundoff leaks found during pre-commit QA;
the same safety implementation is used for both arms.

## Evaluation labels

- `ordinary_id`: ordinary held-out random trajectories;
- `broad_random_in_support`: the historical broad-random split, explicitly not
  called OOD because both arms train on this distribution;
- `moderate_ood_high_frequency`: smooth frequencies 4--6, while training uses
  frequencies 1--3;
- five fixed held-out canonical tests: Sod, Lax, collision, strong-pressure
  jump, and near-vacuum expansion.

Ordinary and moderate suites use the existing 16 snapshots (`t_final=0.006`).
Canonical tests use 64 snapshots (`t_final=0.0252`), long enough for the known
collision/near-vacuum admissibility and fully-discrete entropy safeguards to be
exercised. The historical repository records stress-test CSVs but not the
original evaluator, so this horizon is declared explicitly for reproducibility
and both arms are rerun under it.

The reported Tadmor violation rate uses the repository's established `1e-5`
residual tolerance, with the residual itself evaluated in float64.

## Run

```powershell
python run_coverage_ablation.py --seed 0
```

The screening decision is made from seed 0. Additional seeds are run only if
seed 0 shows a meaningful, non-benchmark-specific improvement without material
ordinary-ID degradation.

The pre-specified positive-screen criterion is:

1. ordinary-ID and moderate-OOD NRMSE each degrade by no more than 5%;
2. at least three of the five canonical tests improve;
3. mean canonical NRMSE improves by at least 5%;
4. no canonical test degrades by more than 10%; and
5. all rollouts remain admissible with zero measured Tadmor violations.
