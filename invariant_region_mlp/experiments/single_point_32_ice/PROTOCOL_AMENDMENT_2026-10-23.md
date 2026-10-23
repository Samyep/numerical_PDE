# Stage-C Protocol Amendment — 2026-10-23

The original frozen protocol remains immutable.  This amendment changes only
the formal repetition policy for Tangent-Actor and the Slurm dispatch
granularity.  Frozen PDEs, parameters, query banks, references, mathematical
sample counts, solver implementations, and method-specific hyperparameters do
not change.

## Primary repetition policy

- Tangent-Actor uses exactly one fresh training run (`seed=0`) for each
  `(PDE, dimension)` cell.  That one global Actor is evaluated at all 32 frozen
  query points with the four preregistered independent MCV inference
  replicates.  The primary analysis estimates inference-MC variability, not
  between-training-seed variability.
- Deep BSDE retains three independent training seeds at every frozen point.
  Its pointwise training randomness is part of the baseline comparison and is
  not assumed negligible.
- FH-MLP retains all six frozen `(n, M)` configurations and ten independent MC
  replicates per configuration.

The already completed Tangent-Actor seeds 1 and 2 at `d=20` are retained as
nonprimary supplementary evidence.  They are not selected using query truth
and are excluded from the amended primary tables so that every dimension uses
the same one-seed Tangent-Actor protocol.

## Allocation bundling

One H100 allocation now executes a complete natural cell before releasing the
GPU:

- Tangent-Actor: one `(PDE, dimension, seed=0)` training followed by all 32
  query evaluations;
- Deep BSDE: one `(PDE, dimension, seed)` allocation sequentially trains the 32
  independent single-initial-state models;
- FH-MLP: one `(PDE, dimension, n, M)` allocation sequentially executes all ten
  MC replicates, each of which evaluates all 32 points.

This produces 40 Slurm allocations per dimension: 4 Tangent-Actor, 12 Deep
BSDE, and 24 FH-MLP.  The underlying scientific task counts are 4, 384, and
240, respectively.  Task-level atomic output, writer locks, exact checkpoint
resume, deterministic seeds, and task-level accounting remain unchanged.  A
retried bundle skips already complete tasks and resumes the first incomplete
task, so bundling changes scheduling overhead but not numerical results.

The amended primary experiment therefore contains 20 Tangent-Actor trainings,
1,920 Deep-BSDE trainings, and 1,200 FH-MLP configuration-replicate estimates
over the same 640 frozen queries.
