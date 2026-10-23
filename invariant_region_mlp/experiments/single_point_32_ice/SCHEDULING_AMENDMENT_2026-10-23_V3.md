# Stage-C Scheduling Amendment V3 — 2026-10-23

This is a scheduling-only refinement of
`single_point_32_ice_v2_ta_seed0_bundled`.  It does not change the primary
seed policy, mathematical tasks, frozen inputs, solver configurations, random
seeds, output records, or analysis.

Starting with `d=100`, all sixty FH-MLP tasks for one `(PDE, dimension)` cell
(six frozen `(n, M)` configurations times ten MC replicates) run sequentially
inside one H100 allocation.  Deep BSDE remains one allocation per
`(PDE, dimension, seed)` containing all 32 pointwise trainings.  Tangent-Actor
remains one allocation per `(PDE, dimension, seed=0)` containing training and
all 32 query evaluations.

The resulting allocation count is 20 per dimension:

- 4 Tangent-Actor allocations;
- 12 Deep-BSDE allocations;
- 4 FH-MLP allocations.

Task-level writer locks, immutable output records, exact checkpoint resume,
and skip-on-complete behavior remain active.  A retry therefore resumes only
the incomplete scientific task inside the larger allocation.

The already running `d=50` queue retains the V2 40-allocation layout so that
no live array index changes meaning.  V3 begins only after the validated
`d=50` audit.
