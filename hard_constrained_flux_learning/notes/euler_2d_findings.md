# 2D Euler HCFL proof of concept

This experiment extends the HLLC-based HCFL construction to the two-dimensional
ideal-gas Euler equations.

## Model

The coarse solver is periodic on a 12x12 Cartesian grid. Reference trajectories
are generated on a 48x48 grid and cell-averaged to the coarse grid.

A single orientation-aware 3x3 network is shared between x- and y-faces. For
each face it proposes a correction to an HLLC base flux:

[
widetilde F_{	heta,n}
=
F_{m HLLC,n}
+
Delta F_	heta(	ext{3x3 local patch}).
]

For the standard mathematical entropy

[
eta(U)=-ho s/(gamma-1),
]

the entropy potential in the coordinate directions is

[
psi_x=m_x,qquad psi_y=m_y.
]

Each proposed face flux is therefore projected onto the affine Tadmor
half-space

[
(v_R-v_L)^T F_n le psi_{n,R}-psi_{n,L}.
]

No flux labels are used in training; supervision is only from coarse
solution-snapshot pairs.

## Three-seed results

| split | Plain | Soft | Hard | coarse HLLC |
|---|---:|---:|---:|---:|
| ID rollout NRMSE | 0.006662 ± 0.000848 | 0.006661 ± 0.000847 | 0.006670 ± 0.000846 | 0.012755 ± 0.001831 |
| OOD rollout NRMSE | 0.018555 ± 0.003391 | 0.018556 ± 0.003391 | 0.018555 ± 0.003401 | 0.029553 ± 0.007609 |

Entropy violation rates:

- Plain ID: 0.67% mean, up to 0.93%.
- Hard ID: 0%.
- Plain OOD: 1.50% mean, up to 2.09%.
- Hard OOD: 0%.

All three learned variants completed the tested rollouts with positive density
and pressure. The hard projection changes the guarantee but essentially not
the rollout error on this moderate 2D distribution.

## Interpretation

Together with the 2D shallow-water experiment, this is direct evidence that the
core HCFL mechanism is not tied to 1D:

- shared face fluxes retain exact conservation in 2D;
- the Tadmor constraint remains affine in each normal face flux;
- one orientation-aware network can be shared between coordinate directions;
- the hard projection removes all measured entropy violations without an
  observable accuracy penalty.

This is still a proof-of-concept distribution. The next necessary experiment is
a strong-shock / implosion-style 2D Euler test, where admissibility limiting
should matter.
