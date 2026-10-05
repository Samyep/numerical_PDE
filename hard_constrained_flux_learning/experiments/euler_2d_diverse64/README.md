# Diverse fixed-64 2-D Euler HCFL training

This experiment tests whether the frozen `flat18` HCFL architecture benefits
from physically diverse high-resolution supervision while deployment remains
strictly fixed at `64 x 64`.

## Locked design

- Reference solver: classic PyClaw 5.10 four-wave 2-D Euler Roe solver with
  `transverse_waves=2` and constant-extrapolation boundaries.
- Reference grid: `512 x 512`, selected after a six-family
  `256/512/1024 -> 64` convergence audit.
- Training/deployment grid: `64 x 64` only.  No cross-grid claim is made.
- Restriction: block average of conserved variables
  `(rho, rho*u, rho*v, E)`.
- Model: the already-selected direct `3 x 6` (`flat18`) face patch, HLLC base,
  signed Roe correction, hard Tadmor projection, and proposal-feasibility
  loss.  The architecture is not changed in this experiment.
- Data: 192 train, 36 validation, and 36 held-out test trajectories, balanced
  over oblique Riemann, contact/shear, oblique quadrant, radial interface,
  colliding-wave, and smooth-packet initial conditions.
- Saved horizon: `t=0.1` at intervals of `0.005`.
- Checkpoint selection: physical completion first and validation rollout NMAE
  second.  Test trajectories are not read during training.

The six deterministic reference-audit cases use seed `84001`.  Train,
validation, and test use disjoint seeds `86001`, `87001`, and `88001`.

## Commands

Dataset generation runs inside the pinned `hcfl-clawpack:5.10.0` container.
After generation:

```text
python run_diverse64.py reference-audit
python run_diverse64.py audit-data
python run_diverse64.py train --seed 0
python run_diverse64.py train --seed 1
python run_diverse64.py train --seed 2
python run_diverse64.py evaluate --seeds 0 1 2
```

Large generated `.npz` files are reproducible and ignored by Git.  Selected
checkpoints, scalar curves, audits, and the final report are retained.
