# 2D shallow-water HCFL proof of concept

This experiment tests whether hard-constrained conservative flux learning
extends beyond 1D without changing the basic mechanism.

## Setup

Periodic 2D shallow water on a Cartesian grid:

[
U=(h,hu,hv)^T.
]

Reference trajectories are generated on a 48x48 grid with SSP-RK2/Rusanov and
cell-averaged to a 12x12 training grid. The learned solver uses a shared
orientation-aware 3x3 local network to propose both x- and y-face fluxes.

For each face, the candidate flux is projected onto the Tadmor half-space

[
(v_R-v_L)^T F_n le psi_{n,R}-psi_{n,L}.
]

For shallow water,

[
v=(gh-	frac12(u^2+v^2),u,v)^T,
]

and for a face with normal direction n,

[
psi_n=	frac12 g h^2 (ucdot n).
]

The network is trained only from coarse solution snapshot pairs; no numerical
flux labels are used.

## Three-seed results

| split | Plain | Soft | Hard | Rusanov |
|---|---:|---:|---:|---:|
| ID rollout NRMSE | 0.06190 ± 0.00673 | 0.06191 ± 0.00673 | 0.06200 ± 0.00689 | 0.08747 ± 0.01467 |
| OOD rollout NRMSE | 0.08145 ± 0.01652 | 0.08122 ± 0.01678 | 0.08127 ± 0.01652 | 0.10319 ± 0.02089 |

Entropy violation rates:

- Plain ID: 5.98% mean, up to 8.50%.
- Soft ID: 5.46% mean, up to 7.79%.
- Hard ID: 0%.
- Plain OOD: 8.48% mean, up to 10.93%.
- Soft OOD: 7.77% mean, up to 9.99%.
- Hard OOD: 0%.

The hard layer therefore removes all measured face entropy violations while
leaving rollout accuracy essentially unchanged. The learned models remain
substantially more accurate than the coarse Rusanov baseline.

## Interpretation

This is the first direct evidence in the project that the HCFL mechanism
extends naturally across spatial dimension:

- conservation remains exact because each face flux is shared;
- the same affine Tadmor projection is applied independently on every x/y face;
- the hard guarantee survives 2D rollout;
- no special 2D penalty or retraining objective is required.

The current 2D test is intentionally modest. It does not yet include dry states,
bathymetry, well-balancedness, or high-order classical baselines.

## Next step

The next strongest experiment is 2D Euler, reusing the HLLC-based proposal and
local admissibility machinery developed in 1D.
