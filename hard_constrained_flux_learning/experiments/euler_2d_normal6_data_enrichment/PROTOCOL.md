# Prospective protocol: training-data-only enrichment

## Question

Does adding more genuinely two-dimensional reference trajectories reduce the
accuracy gap and excess total variation of the retained six-cell normal HCFL
plus fixed Roe transverse transport method?

## Frozen factors

The following are identical to the original experiment:

- model architecture and 7,348 trainable parameters;
- shared x/y six-cell normal stencil;
- fixed parameter-free Roe transverse increment transport;
- feasibility, trajectory, and positivity losses and their weights;
- Adam optimizer, batch size, learning-rate schedule, update cap, validation
  cadence, checkpoint rule, and plateau stopping rule;
- original-training input normalization and conserved-variable loss scales;
- original-training validation and test NMAE scales;
- the 36-trajectory validation and 36-trajectory test archives.

The test archive is not loaded by training and is evaluated only after all
three validation-selected checkpoints have stopped.

## Changed factor

The training archive grows from 192 to 384 trajectories.  It preserves the
original 192 trajectories exactly and appends 32 trajectories from each of
six new families: crossed Riemann, corrugated contact, shock/contact triple
point, shock-vortex, double blast, and four-way collision.

Every added target is generated with PyClaw 5.10 Roe,
`transverse_waves=2`, at 512 by 512 and conservatively restricted to 64 by 64.

## Required decision metrics

The primary comparison is paired held-out rollout NMAE on the unchanged 36
test cases, after averaging each model over three seeds.  Signed density-TV
error, primitive-variable MAEs, positivity, fully discrete entropy balance,
interface projection residual, and conservation are required diagnostics.

An improvement may be attributed only to the combined effect of greater data
quantity and the six added families.  This experiment does not separate those
two aspects of enrichment.
