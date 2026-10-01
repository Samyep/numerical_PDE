# 1D Shallow-Water HCFL Pilot

## Setup

Flat-bottom 1D shallow-water equations

\[
U=(h,m)^T,\qquad
f(U)=\left(m,\;m^2/h+\tfrac12 g h^2\right)^T,\qquad g=9.81.
\]

Physical entropy and entropy variables:

\[
\eta(U)=\frac{m^2}{2h}+\frac12gh^2,
\qquad
v(U)=\left(gh-\frac12u^2,\;u\right)^T,
\qquad
\psi(U)=\frac12ghm.
\]

For fixed left/right states, Tadmor's interface condition is the affine half-space

\[
(v_R-v_L)^T F \le \psi_R-\psi_L.
\]

HCFL-P takes a learned flux proposal and applies the exact Euclidean projection onto this half-space.

## Data

- Reference solver: fine-grid finite volume with Rusanov flux + SSP-RK2.
- Fine grid: 1024 cells.
- Learned/coarse grid: 64 cells (exact cell averaging from the fine grid).
- Snapshot spacing: \`5e-4\`.
- 18 snapshots per trajectory.
- Training: 500 random trajectories per seed.
- ID validation: 100 trajectories per seed.
- OOD validation: 100 trajectories per seed with wider depth/velocity ranges.
- Seeds: 0, 1, 2.
- Neural proposal: five-cell stencil, primitive variables \`(h,u)\`, two hidden layers of width 64.
- Training supervision: one-step coarse trajectory targets only; no numerical-flux labels.

The learned proposal is initialized from coarse Rusanov and learns a bounded local correction. All learned solvers retain the shared-flux conservative finite-volume update.

## Methods

- \`plain\`: learned flux, no entropy constraint.
- \`soft\`: same learned flux plus a soft squared penalty on positive Tadmor residuals.
- \`hard\`: exact half-space projection at every interface.
- \`rusanov\`, \`hll\`: classical coarse-grid baselines.

## Three-seed results

Mean +/- standard deviation rollout normalized RMSE:

| Split | Plain | Soft | Hard (HCFL-P) | Rusanov | HLL |
|---|---:|---:|---:|---:|---:|
| ID | 0.0383 +/- 0.0019 | 0.0422 +/- 0.0020 | 0.0474 +/- 0.0013 | 0.0966 +/- 0.0030 | 0.0940 +/- 0.0029 |
| OOD | 0.1027 +/- 0.0091 | 0.1099 +/- 0.0098 | 0.1146 +/- 0.0099 | 0.2042 +/- 0.0175 | 0.1957 +/- 0.0169 |

Mean entropy-violation rate during rollout:

| Split | Plain | Soft | Hard |
|---|---:|---:|---:|
| ID | 2.61% | 1.93% | 0% |
| OOD | 2.64% | 1.99% | 0% |

For \`hard\`, the maximum positive Tadmor residual over the three runs is about \`2e-6\`, consistent with float32 numerical precision. The classical baselines are also entropy stable by construction.

The hard projection modifies about 5--6% of interface fluxes during rollout:

- ID active fraction: roughly 5.5--6.2% across seeds.
- OOD active fraction: roughly 4.7--6.2% across seeds.

The mean correction norm is small (about \`0.0021--0.0024\` ID and \`0.0033--0.0056\` OOD), so the hard layer acts as a sparse safety correction rather than replacing the learned solver.

All tested rollouts remained positive in water depth. However, entropy stability alone is not a proof of positivity for arbitrary states, so positivity remains an open theoretical/methodological requirement.

## Interpretation

This is the first system-level signal for the project:

1. A five-point learned flux trained only from trajectories can substantially outperform coarse Rusanov/HLL on fine-grid-reference rollouts.
2. Soft entropy regularization reduces but does not eliminate violations.
3. HCFL-P removes all measured violations while retaining most of the learned accuracy advantage over classical coarse solvers.
4. There is a real accuracy/guarantee tradeoff: Plain is more accurate than Hard on this first setup.
5. The hard layer intervenes only on a small fraction of interfaces, supporting the intended \`learn freely, correct only when necessary\` interpretation.

## Implementation note

An early multi-step-training NaN was traced to an autograd implementation issue, not to the mathematical projection: using \`torch.where\` around \`residual / ||a||^2\` still eagerly formed near-\`0/0\` expressions on inactive states. The committed implementation performs the division only on the active mask. With this correction, the projection backward pass is stable.

## Open work

- trajectory multi-step training and longer horizons;
- direct HCFL-F parameterization for vector fluxes;
- formal positivity/invariant-domain enforcement;
- stronger classical baselines (MUSCL, Roe/HLLC where appropriate);
- deterministic reference-grid convergence checks;
- dam-break-specific benchmark suite;
- 1D Euler after SWE stabilizes.
