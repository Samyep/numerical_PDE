# Research Plan: Hard-Constrained Conservative Flux Learning

## 1. Research question

Can a learned numerical flux remain expressive and data-adaptive while satisfying conservation and entropy-stability constraints by construction?

We want a reusable learned finite-volume time stepper, not a PINN that directly represents one solution. The model is rolled out in time:

\[
U^n \to \widetilde F_\theta \to F_\theta^H \to U^{n+1} \to \cdots.
\]

## 2. Proposed model

For each interface, a DNN receives a local stencil and proposes a flux vector

\[
\widetilde F_{i+1/2}=G_\theta(U_{i-r},\ldots,U_{i+r}).
\]

The physical update keeps the standard conservative form

\[
U_i^{n+1}=U_i^n-\lambda(F_{i+1/2}^H-F_{i-1/2}^H),\qquad \lambda=\Delta t/\Delta x.
\]

### Hard entropy layer

For a convex entropy pair, define

\[
a(U_L,U_R)=v_R-v_L,\qquad b(U_L,U_R)=\psi_R-\psi_L.
\]

Project the proposal onto

\[
\mathcal C(U_L,U_R)=\{F:a^T F\le b\}
\]

with

\[
F^H=\widetilde F-\frac{[a^T\widetilde F-b]_+}{\|a\|^2}a.
\]

This is cheap, differentiable almost everywhere, and exact to floating-point precision.

## 3. Training regime

Primary goal: **trajectory-only supervision**.

Train on solution snapshots, not flux labels:

\[
\mathcal L_{\rm traj}=\sum_{k=1}^K w_k\|\widehat U^{n+k}-U_{\rm data}^{n+k}\|^2.
\]

Backpropagate through the unrolled finite-volume steps and the hard projection.

A consistency-preserving parameterization should be used so that

\[
F(U,U)=f(U)
\]

holds by construction or to a very tight tolerance. One candidate is

\[
\widetilde F=F_{\rm base}(U_L,U_R)+s(U_L,U_R)G_\theta(\text{stencil}),\qquad s(U,U)=0.
\]

## 4. Phase A: 1D shallow water

### PDE

\[
U=(h,m)^T,\qquad f(U)=\left(m,\frac{m^2}{h}+\frac12gh^2\right)^T.
\]

Use physical energy

\[
\eta(U)=\frac{m^2}{2h}+\frac12gh^2
\]

with entropy variables

\[
v=(gh-u^2/2,u)^T
\]

and entropy potential

\[
\psi=\frac12ghm.
\]

### Tasks

- Dam break / Riemann problems.
- Smooth periodic waves.
- OOD depths and velocities.
- Long-rollout stability.

### Baselines

- Rusanov / local Lax-Friedrichs.
- HLL / HLLC where appropriate.
- Roe with entropy fix.
- Plain learned flux with no hard layer.
- Soft entropy-penalty learned flux.
- Closest learned conservative/entropy-stable flux methods identified in literature review.

### Metrics

- trajectory `L1/L2` error;
- shock location/speed;
- conservation drift;
- interface Tadmor residual (max/mean/violation rate);
- positivity failures (`h<=0`);
- rollout failure rate;
- runtime and fraction/magnitude of hard corrections.

## 5. Positivity/admissibility

Entropy stability alone does not guarantee admissible states. For shallow water we need `h>0`; for Euler we need `rho>0` and `p>0`.

Investigate whether positivity can be expressed as a cheap additional constraint on flux/update and jointly enforced with the entropy half-space. Candidate strategies:

1. project flux into intersection of entropy half-space and update-level positivity constraints;
2. use an analytic scaling limiter after entropy projection;
3. parameterize a guaranteed-safe baseline plus learned correction and constrain the correction coefficient.

A main theoretical target is a projection/limiter with an explicit guarantee and low inference cost.

## 6. Semi-discrete versus fully discrete entropy

The interface Tadmor condition gives a clean semi-discrete entropy argument. The paper must not overclaim a fully discrete theorem without handling time integration.

Investigate SSP-RK / forward-Euler conditions and CFL restrictions under which the hard flux layer yields a fully discrete entropy or invariant-domain guarantee.

## 7. Phase B: 1D Euler

After shallow water works, move to Euler:

- Sod shock tube;
- Lax problem;
- strong blast wave;
- Shu-Osher interaction;
- OOD Riemann states.

The hard layer should enforce entropy stability and, if possible, density/pressure admissibility.

## 8. Theory targets

Minimal theorem package:

1. **Projection theorem.** Hard layer is the Euclidean projection onto the Tadmor half-space and enforces the entropy interface condition exactly.
2. **Conservation proposition.** Shared projected interface flux preserves the finite-volume conservation identity.
3. **Semi-discrete entropy theorem.** Under the standard convex entropy assumptions, projected fluxes yield a semi-discrete entropy inequality.
4. **Expressivity/noninterference statement.** If the proposal already lies in the admissible half-space, the layer is identity; the correction is minimum norm otherwise.
5. If possible, a fully discrete entropy/admissibility result under an explicit time-stepping/CFL condition.

## 9. What Burgers is for

Burgers is a diagnostic/theory example, not the main contribution.

For a finite Kruzkov set, hard projection improves entropy behavior. As the constraint set approaches the full Kruzkov family for convex Burgers, the admissible boundary recovers the Godunov flux. This is useful for interpretation but also shows why the main method should be tested on systems.

## 10. Go/no-go criteria

Continue aggressively if, on shallow water:

- hard projection produces essentially zero entropy violations;
- trajectory accuracy is competitive with or better than a plain learned flux and standard low-order solvers;
- the method is less dissipative than a strongly safe baseline in smooth regions;
- positivity can be guaranteed without an expensive per-interface nonlinear solve;
- training from trajectory data remains stable over multi-step rollout.

Reframe or stop if the hard layer collapses almost all learned proposals to a standard classical flux, or if positivity/fully-discrete stability requires an optimization problem expensive enough to remove the learned-solver advantage.


## Phase A update: 1D shallow water pilot completed

The first system-level pilot is positive enough to continue.

Current empirical tradeoff:
- Plain has the best raw accuracy but 2--3% entropy-violating interfaces.
- A soft penalty reduces but does not remove violations.
- HCFL-P removes measured violations to float32 precision while retaining a
  large accuracy advantage over coarse HLL/Rusanov.
- The projection changes only about 5--6% of interface fluxes, supporting the
  intended sparse-safety-filter interpretation.

Next priority order:
1. add/derive a positivity or invariant-domain mechanism and a fully discrete
   guarantee compatible with the entropy half-space;
2. implement a direct vector HCFL-F parameterization and compare it against
   projection;
3. move training to stable multi-step rollout and add longer-horizon tests;
4. add MUSCL/high-order entropy-stable classical baselines and reference-grid
   convergence checks;
5. only after the SWE method is stable, move to 1D Euler.


## Phase A.2 update: fully discrete SWE safety completed

The 1D SWE signal remains positive after strengthening the solver.

Implemented:
1. adaptive CFL substepping;
2. an interface-wise conservative flux-correction limiter that guarantees h >= h_floor when the Rusanov low-order endpoint is positive;
3. preservation of the Tadmor half-space under the positivity blend because the constraint is affine in flux;
4. a per-trajectory global convex line search that enforces total periodic-domain entropy non-increase;
5. four-step trajectory fine-tuning;
6. MUSCL-Rusanov and severe near-dry stress tests.

The next main technical target is now **1D Euler**.  The key question is whether the same safe-baseline + learned-correction + convex-limiting construction can jointly preserve density and pressure/internal-energy admissibility while keeping the entropy guarantee cheap.  Do not move to 2D until this is resolved.


## Phase B update: 1D Euler checkpoint

The basic 1D Euler mechanism works and exposes the next bottleneck.

Positive findings:
1. the Euler entropy condition is again an affine half-space in the 3-vector numerical flux;
2. trajectory-only hard HCFL is markedly better than coarse HLL/Rusanov on random training-like trajectories;
3. exact entropy projection does not materially damage accuracy there;
4. conservative admissibility and fully discrete entropy limiters prevent catastrophic density/pressure failure.

Negative/limiting finding:
- MUSCL-HLLC is substantially stronger on canonical severe shock tubes. Broad random fine-tuning improves some stress cases but does not close this gap.

Next priority:
1. replace the Rusanov-centered proposal with a stronger entropy-stable / HLLC-like learned proposal;
2. use multi-step training with broader wave-pattern coverage while keeping benchmark instances held out;
3. derive a local cell/interface admissibility limiter rather than a single per-trajectory convex coefficient;
4. only then consider 2D Euler/SWE.

Do not claim universal solver superiority from the current 1D results.


## Phase C update: 2D systems and OOD trust

The 1D limiter design is now considered structurally complete for the current paper:
- use hard face-level Tadmor projection;
- use local interface-wise admissibility limiting;
- retain the rare global fully-discrete total-entropy line search;
- do not force a cell-local fully-discrete entropy condition, because tested local constructions were significantly more restrictive.

### 2D shallow water

Completed a 3-seed 2D proof of concept with a shared orientation-aware 3x3 face network.
Hard projection removes all measured x/y-face entropy violations while preserving essentially the same rollout accuracy as the unconstrained learned solver.

### 2D Euler

Completed a 3-seed 2D HLLC-based HCFL proof of concept.
Hard entropy projection again gives zero measured violations with negligible accuracy cost on moderate ID/OOD data.

### Severe OOD and trust fallback

Strong 2D Euler blast tests show that hard physical feasibility is not sufficient for predictive reliability: a learned correction can remain entropy-feasible and admissible while being badly wrong far outside training support.

Current trust design:
1. compute a training-calibrated state score from standardized primitive variables;
2. use the 99th-percentile training score as the support threshold;
3. set
   [
   \tau(U)=\operatorname{clip}[(q_{0.99}/z(U))^2,0.05,1];
   ]
4. blend between entropy-projected HLLC and entropy-projected learned HCFL.

The trust layer is nearly identity in-distribution, partially active under moderate OOD, and strongly falls back under extreme blast states.

### Current next priorities

1. **Literature re-audit now.**
   Search specifically for:
   - learned numerical flux + direct Tadmor projection;
   - learned Riemann solver + entropy half-space projection;
   - certified learned hyperbolic solver + invariant-domain/positivity limiting;
   - neural solver + training-calibrated fallback/trust region to classical flux;
   - 2025--2026 work combining learned corrections with convex entropy/admissibility limiting.

2. **Paper-level theory package.**
   Formalize:
   - face-projection theorem;
   - conservation theorem in arbitrary dimension;
   - semi-discrete entropy inequality;
   - convex preservation of entropy feasibility under local admissibility limiting and trust blending;
   - admissibility result conditional on safe low-order update/CFL;
   - global fully-discrete total-entropy guarantee for the final scalar line search;
   - Burgers two-point identifiability theorem.

3. **Publication-quality baselines.**
   Add or verify:
   - high-order classical entropy-stable schemes;
   - HLLC/MUSCL-HLLC in 1D and appropriate 2D baselines;
   - CFN/ESCFN/NESCFN where reproducible;
   - recent hard-constrained neural Riemann solvers.

4. **Reference-resolution checks.**
   Repeat key results with at least one finer reference grid to verify that gains are not artifacts of the reference solver/resolution.

5. **Runtime/memory profile.**
   Report both learned inference and safety-layer overhead.

### Paper framing at this stage

Do not claim a universally better PDE solver.

Preferred framing:

> **Hard-Constrained Conservative Flux Learning: learn accurate coarse-grid flux corrections from trajectories, enforce conservation and entropy by construction, locally recover admissibility, and fall back toward a trusted classical flux only when the state leaves the training support.**

The current empirical evidence spans:
- Burgers theory/identification;
- 1D shallow water;
- 1D Euler;
- 2D shallow water;
- 2D Euler;
- severe OOD stress tests with a training-calibrated trust fallback.

### 1D training-coverage decision update (revised after convergence audit)

A validation-converged seed-0 audit reproduces the earlier 1,100-update
checkpoints exactly and shows that neither broad nor wave training was then
converged. At their validation-selected optima, wave and broad training are
effectively tied: wave is 0.12% worse on validation, 0.52% better on ordinary
ID, 0.42% better on broad-random ID, 0.58% worse on moderate OOD, and 0.98%
worse on the five-case canonical mean.

Do not describe wave learning as a failed direction. The narrower conclusion
is that replacing random extremes with this coarse balanced-wave slice offers
no material overall gain. If coverage is revisited, use additive coverage or
interface-level sampling while retaining random extremes.

### 1D architecture decision update (validation-converged seed 0)

The shared 1,100-update architecture ranking is retired. Best checkpoints now
occur at 6,100--12,500 updates. Learned dissipation leads seed 0, improving
validation by 36.16% and the canonical mean by 28.94% relative to converged
direct correction. Characteristic correction improves them by 23.49% and
14.18%. CNN is effectively tied on validation (+0.34%) and improves the
canonical mean by 9.98%, so it remains a valid control rather than a rejected
architecture.

Next, run multi-seed confirmation for direct, characteristic, dissipation, and
CNN under the same convergence rule, and report every seed against the common
FVM-2048 reference. The seed-0 FVM-2048 re-evaluation preserves the ranking
(dissipation, characteristic, CNN, direct) but raises the best model's absolute
canonical mean by 45.76% relative to FVM-512. Keep the named canonical problems
as final tests; add an independently parameterized stress-validation family so
ordinary validation convergence is not confused with severe-OOD optimality.

### 1D reference-precision decision update

A seed-0 paired ablation increased only the 1D Euler Rusanov + SSP-RK2 teacher
grid from 512 to 2048 cells. Both arms retained the same 580 initial conditions,
64-cell learning grid, model, training schedule, safety stack, and common
2048-cell evaluation references.

The finer teacher improved moderate OOD by 1.46%, collision by 3.39%, and
near-vacuum expansion by 4.51%, but worsened Sod by 3.53% and improved the
five-case canonical mean by only 1.56% at the historical 1,100-update budget.
The new audit shows that this budget is insufficient, so these predictive
differences are provisional rather than a reason to stop the direction.

Reference quality remains a publication requirement. For final results,
replace or verify the exploratory Rusanov teacher with converged high-order
references. Any future 512-versus-2048 or sharper-teacher comparison must train
both arms to the same validation-convergence rule and include an explicit
target-fit/learning-capacity diagnostic.
