# Prospective protocol: diverse fixed-64 2-D Euler confirmation

This protocol was written after the reference/data audits and after observing
validation curves, but before reading any metric from `test.npz`.  The test
archive was generated at the same time as the other splits but is not loaded
by the training command.

## Question

With architecture and deployment grid fixed, does physically diverse
high-resolution supervision allow the 18-cell HCFL numerical flux to improve a
`64 x 64` Euler rollout relative to HLLC-64, while retaining conservation,
positivity, and the hard interface-entropy condition?  PyClaw Roe-64 is the
strong conventional comparator, not a result HCFL is assumed to beat.

## Frozen model

The model is `flat18` from the preceding transverse-stencil ablation: a shared
MLP observes an oriented `3 x 6` primitive-variable patch and predicts four
signed Roe-wave multipliers.  The proposal is HLLC plus the learned Roe-basis
correction.  The update uses the hard Tadmor half-space projection.  Proposal
feasibility is penalized during training.  HLL is available only in the
deployment safety wrapper and its activation rate must be reported.

No architecture, loss weight, stencil, or test-dependent hyperparameter is
changed in this confirmation.  Training uses three seeds and a 50,000-update
cap, with validation-driven checkpoint selection, learning-rate reductions,
and stopping at a minimum-learning-rate plateau.

## Reference and splits

- Fine reference: classic PyClaw 5.10 four-wave 2-D Euler Roe with two
  transverse corrections on `512 x 512`.
- Coarse representation and deployment: `64 x 64`.
- Restriction: exact block averaging of fine-cell conserved variables.
- Train: 192 trajectories, seed `86001`.
- Validation: 36 trajectories, seed `87001`.
- Test: 36 trajectories, seed `88001`.
- End time: `0.1`; saved interval: `0.005`.

Every split is balanced over six families: oblique Riemann, contact/shear,
oblique quadrant, radial interface, colliding waves, and smooth wave packets.
Serialized parameter specifications are exactly disjoint between splits.

A one-case-per-family convergence audit compares `256/512/1024 -> 64`.
The observed `512` versus `1024` rollout NMAE is `0.0012335`, small relative
to the untrained HLLC validation NMAE `0.0283969`; this supports `512` as the
reference used for this fixed-grid experiment without calling it exact.

## Model selection and test use

Selection is lexicographic: full validation completion first, then aggregate
validation rollout NMAE.  Per-family validation NMAE is diagnostic and cannot
replace the locked aggregate selection criterion.  The held-out test archive
is evaluated once after all three seed checkpoints have stopped.

Every test method starts from the same `512 -> 64` restricted initial state.
The native Roe-64 baseline is rerun from that state; it does not receive a
cell-centre reconstruction of the analytic initial condition.

## Required outputs

For HLLC-64, native PyClaw Roe-64, and HCFL over all three seeds, report:

1. completion and minimum density/pressure;
2. aggregate and per-family NMAE/NRMSE plus primitive-variable MAE;
3. signed density-TV error and positive TV excess;
4. maximum interface-projection residual and fully discrete entropy balance;
5. conservation and any positivity/entropy fallback activation;
6. nonzero learned-correction diagnostics, to rule out convergence to HLLC.

The result is a fixed-64 learned-closure result.  It does not support a
cross-resolution generalization claim.
