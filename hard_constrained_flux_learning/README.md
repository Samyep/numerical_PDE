# Hard-Constrained Conservative Flux Learning

Research prototype for learning numerical fluxes from data while enforcing conservation and entropy admissibility by construction.

## Core idea

For a 1D hyperbolic conservation law

\[
U_t + f(U)_x = 0,
\]

we retain the conservative finite-volume update

\[
U_i^{n+1}=U_i^n-\frac{\Delta t}{\Delta x}\left(F_{i+1/2}-F_{i-1/2}\right),
\]

but replace a hand-designed numerical flux by a learned proposal

\[
\widetilde F_{i+1/2}=\mathrm{DNN}_\theta(\text{local stencil}).
\]

A differentiable hard layer then projects the proposal into an entropy-admissible set before the conservative update is applied.

For systems with a convex entropy pair \((\eta,q)\), entropy variables \(v=\nabla_U\eta\), and entropy potential \(\psi=v^T f-q\), Tadmor's interface condition is

\[
(v_R-v_L)^T F \le \psi_R-\psi_L.
\]

For fixed left/right states, this is an affine half-space constraint on the flux vector. The Euclidean hard projection is

\[
F^H=\widetilde F-\frac{[(v_R-v_L)^T\widetilde F-(\psi_R-\psi_L)]_+}{\|v_R-v_L\|_2^2}(v_R-v_L).
\]

Thus an arbitrary learned proposal can be minimally corrected to satisfy the interface entropy inequality to machine precision.

## What is already validated

1. **Burgers flux sanity check.** A poor central-flux proposal becomes much more robust after hard entropy correction. Finite Kruzkov constraint sets progressively suppress shock overshoot. In the continuum-constraint limit for convex Burgers, the construction recovers the Godunov flux, so Burgers is primarily a sanity/theory case rather than the main novelty.
2. **Learned Burgers flux prototype.** A small DNN learned a local flux rule and was rolled out autoregressively through the conservative update. Hard Kruzkov layers improved several shock/rarefaction tests and strongly reduced nonphysical overshoot. More constraints were not always lower-error because stronger enforcement can increase dissipation.
3. **Shallow-water hard projection.** For 200,000 random physically admissible interfaces, a deliberately corrupted flux proposal violated Tadmor's inequality on about half the samples. Closed-form half-space projection reduced all violations below `1e-10`; the maximum residual after projection was about `8.5e-14`.

See `STATUS.md` for exact boundaries between verified facts and open work.

## Repository layout

- `plan.md` — research plan and go/no-go criteria.
- `STATUS.md` — verified results, theoretical observations, and unresolved claims.
- `experiments/burgers_hardnet_sanity.py` — non-learned Burgers flux projection sanity check.
- `experiments/burgers_dnn_flux.py` — compact learned-flux prototype.
- `experiments/shallow_water_projection.py` — closed-form shallow-water Tadmor projection test.
- `results/` — snapshots from the exploratory runs.
- `notes/literature_gap.md` — current positioning and nearest-neighbor checklist.

## Current intended positioning

**Hard-constrained conservative flux learning**: learn numerical fluxes from trajectory data, preserve conservation through flux form, and enforce entropy stability with an architecture-level hard projection rather than a soft penalty.

The main target is **systems** (starting with 1D shallow water, then 1D Euler), not scalar Burgers.


## 1D SWE Phase 2

The SWE prototype now includes a fully discrete safety wrapper:
- adaptive CFL substepping;
- conservative positivity limiting for water height;
- global periodic-domain entropy non-increase via a cheap convex line search;
- four-step trajectory fine-tuning;
- MUSCL-Rusanov and near-dry baselines.

See:
- `experiments/swe_1d_phase2_safe.py`
- `experiments/swe_1d_phase2_multistep.py`
- `notes/swe_1d_phase2_findings.md`
- `results/swe_1d_phase2/`

On the near-dry stress test, the 3-seed multi-step HCFL-safe model reaches NRMSE `0.177 +/- 0.026`, versus `0.207 +/- 0.018` for MUSCL-Rusanov, while maintaining `h >= 1e-4` under adaptive CFL substepping.


## 1D Euler checkpoint

The project now includes a 1D ideal-gas Euler HCFL checkpoint.

Validated in 3-seed trajectory experiments:
- five-point learned flux trained only from downsampled high-resolution trajectories;
- exact Tadmor half-space projection using the standard mathematical entropy;
- hard model has zero measured interface entropy violations to numerical precision;
- conservative convex admissibility limiter enforces `rho>0` and `p>0` relative to an admissible Rusanov low-order endpoint;
- a second line search enforces nonincrease of total mathematical entropy over the fully discrete periodic step;
- on the original random ID distribution, hard HCFL NRMSE is about `0.0229`, versus `0.0323` for MUSCL-HLLC;
- broad-distribution fine-tuning reaches about `0.0730` on the widened random test distribution, versus `0.0795` for MUSCL-HLLC.

Important limitation: on canonical severe shock-tube tests (Sod, Lax, opposing-flow collision, strong pressure jump, near-vacuum expansion), MUSCL-HLLC remains substantially more accurate. The full HCFL safety wrapper prevents negative-density/pressure failures that occur for unconstrained or entropy-only learned rollouts.

Files:
- `notes/euler_1d_findings.md`
- `euler_1d_phase/euler_1d_phase.zip` (all Euler source code + summary CSVs)
