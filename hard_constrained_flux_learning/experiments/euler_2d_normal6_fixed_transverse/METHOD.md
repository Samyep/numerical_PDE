# Dimension-consistent 2-D construction

For an oriented face, the learned normal flux is exactly the retained 1-D
construction with the extra tangential momentum component:

\[
F^N_{i+1/2,j}
=\frac{F(U_{i,j})+F(U_{i+1,j})}{2}
-\frac12 R_x
\left(d_\theta\odot |\lambda_x|_{\rm ef}\odot\alpha_x\right),
\qquad 0\le d_\theta\le2.
\]

The single network sees the six cells along the face normal and outputs four
nonnegative Roe multipliers.  The same weights are used after orienting every
y face; there is no x-specific or y-specific network.

## Fixed corner transport

Every conservative numerical face flux defines consistent normal
fluctuations

\[
D^- = F^N-F(U_L),\qquad D^+=F(U_R)-F^N,
\]

so that \(D^-+D^+=F(U_R)-F(U_L)\).  A parameter-free transverse Roe solver
splits each fluctuation,

\[
D\mapsto B^-D+ B^+D,
\]

using the two transverse acoustic waves and the repeated contact/shear speed.
For an x-normal sweep, the down- and up-going pieces modify the two neighboring
y-face fluxes by the standard corner coefficient

\[
-\frac12\frac{\Delta t}{\Delta x} B^\pm D.
\]

The y-normal sweep supplies the symmetric correction to the x-face fluxes.
The final update remains a single conservative flux divergence.

## Constraint ordering

The fixed transverse correction is part of the raw proposal
\(\widetilde F\).  The feasibility loss is evaluated on this complete proposal,
and the flux used by the PDE update is

\[
F^H=P(\widetilde F),
\]

where \(P\) is the same hard interface entropy projection used by the 1-D
method.  Thus transverse transport cannot bypass the hard face constraint.
The low-order HLL flux remains confined to the deployment safety wrapper; it is
not used to generate training proposals or training targets.

## Exact 1-D reduction

If the state is constant in the transverse direction and the transverse
velocity is zero, every transverse face receives the same fixed correction.
Its discrete transverse divergence is therefore exactly zero, while the
orthogonal normal fluctuation is zero.  The complete 2-D update reduces to the
same six-cell 1-D HCFL update (verified by an automated numerical invariant
test).

## Why the implementation is level 1, not level 2

This implements Clawpack's standard transverse **increment-wave** transport,
corresponding to `transverse_waves=1`.  In Clawpack,
`transverse_waves=2` additionally transports the explicit second-order normal
correction wave `cqxx`.  That correction is a component of Clawpack's separate
second-order reconstruction/limiting algorithm and is absent from the retained
1-D HCFL update.  Adding it only in 2-D would break the method-consistency goal.

Reference implementation: Clawpack 5.10
[`flux2.f90`](https://github.com/clawpack/pyclaw/blob/v5.10.0/src/pyclaw/classic/flux2.f90)
and
[`rpt2_euler.f90`](https://github.com/clawpack/riemann/blob/v5.10.0/src/rpt2_euler.f90).
