# Scope amendment: C0/C1/C2 only

Effective 2026-10-27 at 21:12 UTC, before the final-reference truth gate was
opened, the user instructed the agent to execute and report only the new Deep
BSDE C0/C1/C2 sweeps and not the joint J1/J2/J3 scaling study.

The frozen preregistration, original 42-task manifest, and already-created
artifacts remain unchanged for provenance.  The active analysis scope contains
30 tasks (1,920 terminal query records): five unique C0/C1/C2 trajectories for
each of two PDEs and three training seeds.  Original Slurm array indices are
preserved so pending tasks cannot be reinterpreted.

At the instruction time, task 5 (`b6_fullrank`, seed 0, `joint_j1`) had already
completed 64 queries.  Task 6 (`b6_fullrank`, seed 0, `joint_j3`) was running
and was cancelled with four completed query records.  The remaining ten J1/J3
array tasks had not started and were cancelled.  These excluded artifacts are
retained, never enter accuracy aggregation or tables, and their actual GPU
allocation remains included in resource accounting.  J2 required no separate
training because it coincided with the common C0/C1/C2 reference checkpoint.

No numerical hyperparameter, seed, query bank, stopping rule, threshold, or
evaluation rule changed.  The final-reference truth remained unopened when
this amendment was applied.
