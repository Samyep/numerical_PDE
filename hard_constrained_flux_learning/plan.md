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
