# Results: zero-shot nonperiodic Euler audit

Seed: 0.  The checkpoints were trained only on periodic 64-cell trajectories
and were not retrained for this experiment.

## Boundary treatment

- Boundary condition: transmissive / constant extrapolation.
- Neural fluxes: only the 511 interior interfaces of the 512-cell deployment.
- Domain boundary fluxes: the physical Euler flux of the endpoint state.
- Near-boundary stencil: replicated ghost cells, with no circular wrapping.
- Interior learned proposals retain the Tadmor hard projection.
- The fully-discrete entropy target includes the physical left/right boundary
  entropy flux.

The operator self-tests verify equal-state consistency, zero left/right wrap
leakage, exact classical boundary fluxes, constant-state preservation, and the
finite-volume boundary conservation identity.

## Test suite

The common reference is strict nonperiodic HLLC-2048 + SSP-RK2,
conservatively restricted to 512 cells only for metrics.  All runs end at
`t=0.0252`.

- Five centered cases: Sod, Lax, collision, strong pressure, near-vacuum
  expansion.
- Three boundary-interaction cases: a contact exiting each side and a strong
  pressure wave exiting the right side.

## Aggregate result

| Method | 8-case rollout NRMSE | 8-case final NRMSE | Centered rollout | Boundary-interaction rollout | Normalized TV excess | Excess extrema |
|---|---:|---:|---:|---:|---:|---:|
| Native HLLC-512 | 0.0512534 | 0.0534821 | 0.0515424 | 0.0507717 | **0.01514** | **0** |
| HLLC + Roe correction | **0.0347194** | **0.0308995** | **0.0334210** | **0.0368833** | 0.04533 | 51 |
| Central + nonnegative Roe, no feasibility loss | 0.0356543 | 0.0320357 | 0.0336011 | 0.0390761 | 0.04912 | 58 |
| Central + nonnegative Roe + feasibility (`1e-3`) | 0.0368331 | 0.0360070 | 0.0358374 | 0.0384924 | 0.06632 | **31** |

Relative to native HLLC-512, HLLC + Roe correction lowers the eight-case
rollout/final error by 32.26% / 42.22%.  The feasibility-trained nonnegative
Roe model lowers them by 28.14% / 32.67%.  Both improvements remain present
on the three cases that actually interact with a boundary.

## Safety and interpretation

All learned runs completed 8/8 cases.  No local positivity limiter activation
or low-order time-step halving was needed.  The maximum recorded interior
Tadmor residual is below `6.2e-13`, and every boundary-aware fully-discrete
entropy balance has nonpositive residual.

HLLC + Roe correction is the strongest nonperiodic method in this seed: it has
the lowest overall, centered, and boundary-interaction rollout error.  The
feasibility-trained nonnegative Roe method remains useful and has fewer excess
extrema than the other learned arms, but it is less accurate, has larger
normalized TV excess, and reaches a smaller pressure margin (`8.75e-4` versus
`9.02e-3` for HLLC + Roe correction).

This experiment supports the proposed design principle: for transmissive
boundaries, the network need not learn a boundary flux.  A periodic-trained
local interior flux can be transferred zero-shot when the boundary operator is
handled analytically.

Limitations: one seed, short final time, and only transmissive boundaries.
Reflective-wall and prescribed inflow boundaries require separate
boundary-entropy operators and are not established by this result.

## Artifacts

- `results/nonperiodic_transmissive512_seed0.png`
- `results/nonperiodic_transmissive512_seed0.json`
- `results/nonperiodic_transmissive512_seed0.csv`
- `results/scientific_integrity_audit_seed0.json`
