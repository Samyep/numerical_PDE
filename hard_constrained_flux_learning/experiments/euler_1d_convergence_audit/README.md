# 1D Euler validation-convergence audit

## Motivation

Earlier 1D Euler architecture and wave-coverage comparisons used a fixed budget
of 1,100 Adam updates. Their logged minibatch losses did not establish a
validation plateau. Consequently, conclusions such as "Conv1d is worse" or
"wave coverage does not help" are provisional until all arms receive the same
validation-controlled opportunity to converge.

## Pre-run hypothesis

> Some rankings obtained at 1,100 updates are optimization-budget artifacts.
> Training every candidate with a shared validation-based learning-rate and
> early-stopping rule may materially change at least one architecture or
> training-distribution conclusion.

This hypothesis and the convergence rule below were recorded before the run.

## Arms

Architecture comparison on the same broad-random training tensor:

- `direct_broad`;
- `invariant_broad`;
- `characteristic_broad`;
- `dissipation_broad`;
- `conv_broad`.

Training-distribution comparison with the direct architecture fixed:

- `direct_broad`;
- `direct_wave`.

The wave arm is the previous balanced contact/compression/expansion/
pressure-jump/collision replacement of the final 100 trajectories. The common
direct-broad arm serves as the comparator in both audits.

## Shared validation set

Validation trajectories are independent of training and final testing. The
common validation set contains:

- 44 ordinary trajectories;
- 52 broad-random trajectories;
- 20 random-extreme trajectories;
- 20 structured wave-regime trajectories.

This union mirrors the shared 480-trajectory core and gives equal validation
weight to the two alternative 100-trajectory tails. No canonical benchmark is
used for optimization, scheduling, checkpoint selection, or stopping.

## Convergence protocol

- reproduce the original first 1,100 updates exactly;
- validate every 100 updates;
- choose checkpoints by safe validation rollout NRMSE;
- after update 1,100, reduce learning rate by `0.3` after ten validation checks
  without at least `0.1%` relative improvement;
- continue down to a minimum learning rate of `3e-6`;
- declare convergence after another ten non-improving checks at the minimum
  learning rate;
- cap at 20,000 updates and report explicitly if the cap is reached;
- retain both the exact 1,100-update checkpoint and the best validation
  checkpoint for every arm.

The deterministic full-training one-step loss, validation one-step loss, and
safe validation rollout NRMSE are recorded at every validation check.

## Held fixed

- seed-0 512-cell strict Rusanov + SSP-RK2 references restricted to 64 cells;
- 580 training trajectories per arm and the previous distributions;
- width, initialization seed, minibatch-index sequence, Adam, gradient clipping,
  and initial learning rate;
- baseline-derived normalization for every arm;
- hard Tadmor projection and the complete shared rollout safety stack;
- final test suite and metrics.

Boundary treatment remains periodic during this optimization audit. Physical
non-periodic boundary closures are a separate change and are not mixed into
this run.

## Run

```powershell
python run_convergence_audit.py --seed 0
python plot_results.py --seed 0
python plot_best_vs_fvm.py --seed 0
```

## Result

The seed-0 run completed with all six arms satisfying the validation-plateau
criterion before the update cap. See [`RESULTS.md`](RESULTS.md) for the revised
architecture and wave-coverage conclusions. The FVM comparison plot shows the
HCFL-64 and an unaveraged 64-point subsample of FVM-2048 only. The fine-grid
series takes the near-center cell (`16 + 32*i`) from each consecutive block of
32 and connects those values with an ordinary line. Like-for-like numerical
error metrics still use the conservatively averaged reference.

The same command also writes
`results/fvm512_vs_hcfl512_with_fvm2048_seed0.png`: native FVM-512 and the
64-grid-trained HCFL flux deployed on 512 cells are emphasized, while native
FVM-2048 remains as a light background curve. HCFL uses eight CFL-matched
substeps per saved interval. This is a resolution-transfer visualization, not
a separately trained 512-cell checkpoint.
