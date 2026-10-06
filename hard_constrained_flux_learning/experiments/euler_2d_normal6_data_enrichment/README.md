# Data-only 2-D Euler enrichment ablation

This experiment changes only the training archive used by the converged
six-cell normal HCFL plus fixed Roe transverse transport experiment.  The
network, flux, projection, losses, optimizer, update cap, stopping rule,
validation split, test split, and deployment grid remain unchanged.
The original training archive also supplies the frozen input normalization,
conserved-variable loss scales, validation NMAE scales, and test NMAE scales;
the enlarged archive changes only which trajectories are sampled.

The original 192 trajectories are retained.  The enriched archive appends 32
fine-grid trajectories from each of six genuinely two-dimensional families:

- crossed Riemann problems with non-orthogonal discontinuities;
- corrugated contact/shear interfaces;
- shock/contact triple points;
- shock-vortex interactions;
- interacting double blasts;
- four-way collisions.

All extra targets use the same PyClaw 5.10 Roe solver with
`transverse_waves=2` on 512 by 512, followed by conservative restriction to
64 by 64.  The final training archive therefore has 384 trajectories.  The
locked 36-trajectory validation and 36-trajectory test archives are reused
byte-for-byte.

This design tests the data-coverage hypothesis.  It does not test a new loss
or a new numerical method.
