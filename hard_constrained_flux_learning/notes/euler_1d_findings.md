# 1D Euler findings

## Setup

We study the ideal-gas Euler equations

\[
U=(\rho,m,E)^T,\qquad
f(U)=\left(m,\frac{m^2}{\rho}+p,\frac{m}{\rho}(E+p)\right)^T,
\]

with

\[
p=(\gamma-1)\left(E-\frac{m^2}{2\rho}\right),\qquad \gamma=1.4.
\]

The mathematical entropy is

\[
\eta(U)=-\frac{\rho s}{\gamma-1},\qquad
s=\log p-\gamma\log\rho,
\]

with entropy variables

\[
v=\left(\frac{\gamma-s}{\gamma-1}-\frac{\rho u^2}{2p},\frac{\rho u}{p},-\frac{\rho}{p}\right)^T
\]

and entropy potential \(\psi=m\). Hence, for fixed left/right states, the Tadmor condition

\[
(v_R-v_L)^T F\le \psi_R-\psi_L
\]

is an affine half-space constraint on the three-component numerical flux.

## Learned model

A five-point stencil in primitive variables is mapped to a three-component learned correction on top of a Rusanov proposal. Training uses only coarse solution snapshots obtained by cell-averaging high-resolution Rusanov/SSP-RK2 reference trajectories. No numerical-flux labels are used.

The hard model applies the closed-form projection onto the Tadmor half-space. A full safe rollout additionally uses:

1. an adaptive-CFL Rusanov low-order update;
2. a conservative per-trajectory convex limiter between the low-order and learned entropy-feasible flux to enforce \(\rho>0\) and \(p>0\);
3. a second convex line search to enforce nonincrease of total mathematical entropy over the fully discrete step.

Because the ideal-gas admissible set \(\{\rho>0,p>0\}\) is convex, the conservative line search is well-defined whenever the low-order update is admissible.

## Three-seed random-trajectory results

| split | method | rollout NRMSE |
|---|---|---:|
| ID | Plain | 0.02233 +/- 0.00159 |
| ID | Soft entropy penalty | 0.02244 +/- 0.00159 |
| ID | HCFL hard entropy | 0.02293 +/- 0.00157 |
| ID | HLL | 0.04761 +/- 0.00240 |
| ID | Rusanov | 0.05025 +/- 0.00263 |
| ID | MUSCL-HLLC | 0.03228 +/- 0.00179 |
| OOD (original model) | Plain | 0.08382 +/- 0.00453 |
| OOD (original model) | Soft | 0.08402 +/- 0.00454 |
| OOD (original model) | HCFL hard | 0.08389 +/- 0.00479 |
| OOD (original model) | HLL | 0.11855 +/- 0.00250 |
| OOD (original model) | Rusanov | 0.12843 +/- 0.00321 |
| OOD (original model) | MUSCL-HLLC | 0.07946 +/- 0.00228 |

The hard model has zero measured Tadmor violations to numerical precision. Plain and soft models violate the interface inequality on roughly 2--3% of evaluated interfaces.

## Broad-distribution fine-tuning

Fine-tuning the hard model on a wider random distribution reduces the error on that wider test distribution from 0.08389 to 0.07297 while retaining the hard entropy projection. This comparison should not be described as OOD after fine-tuning, because the widened training distribution overlaps the original OOD range.

## Canonical/strong stress tests

The following tests are deliberately far from the initial random-training distribution: Sod, Lax, opposing-flow collision, strong pressure jump, and near-vacuum expansion.

The full safe limiter is essential. Plain and entropy-only learned fluxes can generate negative pressure on Sod, collision, and near-vacuum tests. The full safe method completes all tests while preserving admissibility.

However, MUSCL-HLLC is substantially more accurate on these classical stress tests. Representative NRMSEs:

| case | broad HCFL-safe | MUSCL-HLLC |
|---|---:|---:|
| Sod | 0.1050 | 0.0447 |
| Lax | 0.2289 | 0.1462 |
| collision | 0.2519 | 0.1522 |
| strong pressure | 0.3149 | 0.1644 |
| near-vacuum expansion | 0.3182 | 0.1626 |

This is an important limitation: the current learned solver is competitive on the trajectory distribution it was trained to approximate, but is not yet a universal replacement for mature high-order shock-capturing methods.

## Current interpretation

The strongest supported claim is not that HCFL beats classical solvers everywhere. It is:

> A learned coarse-grid conservative flux can be trained from trajectories, projected exactly into the entropy-stable half-space, and combined with conservative admissibility/fully-discrete limiters so that catastrophic learned-solver failures are prevented. On the training-like distribution the learned coarse solver can outperform standard coarse classical baselines, while extreme shock-tube generalization remains a limitation.

## Next research questions

1. Replace Rusanov-centered learned corrections with a stronger entropy-stable or HLLC-like proposal while retaining hard guarantees.
2. Train with multi-step objectives and broader wave-pattern coverage without explicitly training on the final benchmark instances.
3. Derive a less global admissibility limiter (cell/interface local rather than one scalar per trajectory) to reduce unnecessary fallback toward the low-order flux.
4. Move to 2D only after the 1D proposal/limiter tradeoff is satisfactory.
