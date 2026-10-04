# Prospective 2-D Euler / PyClaw evaluation protocol

This protocol was frozen before training or evaluating a 2-D HCFL model.  The
public gallery initial condition and a small PyClaw installation/reference-grid
smoke test were inspected first; no learned 2-D Euler prediction existed when
this file was committed.  This is therefore a prospective model-evaluation
record, not a blind or externally registered study.

## Question

Can the same HCFL construction used in one dimension -- an HLLC base flux, a
learned Roe-wave correction, a hard Tadmor projection, and a proposal
feasibility loss -- improve a 64 by 64 finite-volume rollout on a genuinely
two-dimensional shock interaction without losing conservation, positivity, or
entropy admissibility?

## Locked public test

The main out-of-distribution test is the Clawpack 5.10 PyClaw gallery example
`2-dimensional Euler equations`, also called the Liska--Wendroff quadrant
problem:

<https://www.clawpack.org/gallery/pyclaw/gallery/quadrants.html>

The domain is `[0,1] x [0,1]`, `gamma=1.4`, the two discontinuities are at
`x=0.8` and `y=0.8`, and all four boundaries use constant extrapolation.  The
four primitive states `(rho,u,v,p)` are, in `(upper-right, upper-left,
lower-left, lower-right)` order,

1. `(1.5, 0, 0, 1.5)`;
2. `(0.532258064516129, 1.206045378311055, 0, 0.3)`;
3. `(0.137992831541219, 1.206045378311055, 1.206045378311055,
   0.029032258064516)`;
4. `(0.532258064516129, 0, 1.206045378311055, 0.3)`.

The locked end time is `t=0.8`.  This exact case, its rotations/reflections,
and its states are excluded from training and validation.

## Reference and common initial state

The reference generator is Clawpack 5.10.0, classic PyClaw, the four-wave 2-D
Euler Roe solver, two transverse-wave corrections, and extrapolation boundary
conditions.  It is built by `Dockerfile.clawpack` so the executable environment
is reproducible.

The formal gallery reference uses `1024 x 1024` cells and is conservatively
block-averaged to `64 x 64`.  A `256/512/1024` reference-grid audit is reported
separately.  Every coarse method starts from the *same* restricted
`1024 x 1024` initial finite-volume state.  Thus a method cannot benefit from a
different cell-centre discretization of the discontinuity at `0.8`.

The native PyClaw Roe-64 comparator is re-run from that same common coarse
initial state.  No fine-grid values after `t=0` are supplied to any learned or
coarse-grid method.

## Training and validation data

Training data are also produced by the same isolated PyClaw configuration, on
`256 x 256` grids and conservatively restricted to `64 x 64`.

- training: 48 trajectories, RNG seed `81001`;
- validation: 12 disjoint trajectories, RNG seed `82001`;
- held-out in-distribution test: 12 trajectories, RNG seed `83001`;
- saved interval: `0.005`;
- training/validation/test end time: `0.10`.

For every trajectory, four quadrant primitive states are sampled independently:
`log(rho)` is uniform between `log(0.5)` and `log(2.0)`, `log(p)` is uniform
between `log(0.4)` and `log(2.0)`, and both velocity components are uniform on
`[-0.6,0.6]`.  The x and y split locations are uniform on `[0.30,0.70]`.
The seeds fix both the samples and their order.  A generated trajectory that is
non-finite or has non-positive density or pressure is a recorded generation
failure; it is not silently replaced by a new random draw.

## HCFL model and training

The state ordering is `(rho,rho*u,rho*v,E)`.  For each coordinate direction it
is re-oriented to `(rho,m_normal,m_tangent,E)`.  One network is shared by both
directions.  It predicts four signed multipliers in the local four-wave Roe
basis from an interface-centred primitive-variable stencil.  Its proposal is

`F_tilde = F_HLLC - 0.5 R (d * abs(lambda) * alpha)`.

The final learned flux used by the update is the Euclidean projection of this
proposal onto the Tadmor entropy half-space.  A squared feasibility loss on the
*unprojected* proposal has weight `1e-3`.  The HLL safety flux is not used in the
training forward pass.

Four- and six-cell stencils are compared using only the validation set and seed
zero.  Selection is lexicographic: full physical completion first, then lower
validation NMAE.  The selected stencil is trained with seeds `0,1,2`.  Adam,
width 72, batch size 4, initial learning rate `5e-4`, validation every 500
updates, and a hard cap of 50,000 updates are fixed.  Learning-rate reductions
and stopping are driven only by validation plateaus.

At deployment, the hard interface projection is always active.  A global
convex blend toward a positivity/entropy-stable HLL step is allowed only as a
safety wrapper when the projected proposal would violate positive density,
positive pressure, or the fully discrete entropy inequality.  Its activation
rate and minimum blend coefficient are mandatory outputs; a result that hides
fallback use is invalid.

## Comparators

1. native PyClaw Roe-64 from the common coarse initial state;
2. first-order HLLC-64;
3. HLLC plus the selected HCFL Roe correction, reported over three seeds.

This first experiment tests the HCFL mechanism against auditable finite-volume
baselines.  FNO/PINN/learned-flux literature adaptations are not added until
this physical audit passes; otherwise architecture differences would obscure
the numerical question.

## Metrics and failure accounting

Evidence is ordered as follows:

1. completion and finite/positive density and pressure;
2. maximum hard-projection interface entropy residual;
3. fully discrete entropy-balance violations;
4. boundary-aware conservation closure;
5. NMAE (primary accuracy metric), NRMSE, and raw MAE in `(rho,u,v,p)`;
6. density/pressure range overshoot, total-variation excess, and wall time.

NMAE is

`mean(abs((primitive_prediction - primitive_reference) / training_std))`,

where `training_std` is fixed from the training trajectories.  NRMSE uses the
same channels and scale.  Failed trajectories remain in the completion
denominator; conditional errors are never presented without completion.

No hyperparameter, checkpoint, or method is selected using the official
quadrant result.  Any later change to this protocol must be recorded as a
dated amendment before the affected run.

## Amendments

### 2026-10-04, before formal data generation or model training

The `0.005` saved interval is restricted to the random training, validation,
and ID-test trajectories.  The official `t=0.8` gallery case is saved at the
gallery's 10 equal output intervals (`0.08`, giving 11 states including the
initial state).  All coarse solvers still integrate with adaptive internal
substeps and audit every internal step.  This change avoids retaining 161
full `1024 x 1024 x 4` snapshots solely for plotting/error quadrature; it does
not alter the reference evolution, endpoint, model selection, or safety audit.
