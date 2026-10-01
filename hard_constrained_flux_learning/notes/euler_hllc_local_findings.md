# 1D Euler: HLLC proposal and local admissibility limiting

## What changed

The original Euler HCFL used a Rusanov base plus a learned five-point correction.
This is robust but unnecessarily diffusive on contact/shock structure.

This phase makes two changes:

1. **HLLC-based proposal**
   \[
   \widetilde F_\theta = F_{\rm HLLC} + \Delta F_\theta(\text{five-point stencil}),
   \]
   followed by the same hard Tadmor projection.

2. **Local conservative admissibility limiter**
   A global trajectory-wise blending coefficient was replaced by shared
   interface coefficients \(\alpha_{i+1/2}\in[0,1]\). The final flux is
   \[
   F_{i+1/2}
   =
   F^{\rm lo}_{i+1/2}
   +
   \alpha_{i+1/2}
   \left(F^{\rm hi}_{i+1/2}-F^{\rm lo}_{i+1/2}\right).
   \]
   Because the same interface coefficient is used by both neighboring cells,
   conservation is preserved. Because the Tadmor feasible set is affine/convex
   in the interface flux, convex interpolation between two entropy-feasible
   fluxes remains entropy feasible.

The low-order reference is entropy-projected Rusanov. The local limiter reduces
only interfaces adjacent to inadmissible cells and iterates until
\(\rho>0,\ p>0\). A trajectory-wise fully-discrete total-entropy line search is
retained as a final safeguard.

## Main stress-test results

Three-seed mean NRMSE, using broad training:

| case | Rusanov-base HCFL, global safe | HLLC-base HCFL, global safe | HLLC-base HCFL, local safe | MUSCL-HLLC |
|---|---:|---:|---:|---:|
| Sod | 0.1050 | 0.0765 | **0.0765** | 0.0447 |
| Lax | 0.2289 | 0.1958 | **0.1958** | 0.1462 |
| collision | 0.2519 | 0.2347 | **0.0952** | 0.1522 |
| strong pressure | 0.3149 | 0.2498 | **0.2498** | 0.1644 |
| near vacuum | 0.3182 | 0.2816 | **0.1791** | 0.1626 |

The local limiter is the major gain on the difficult admissibility-dominated
cases:

- collision: 0.2347 -> **0.0952**;
- near vacuum: 0.2816 -> **0.1791**.

It modifies only about 2--3% of interfaces on these tests, instead of scaling
the whole trajectory back toward Rusanov.

## Multi-step fine-tuning

A four-step fine-tuning pass was also tested on the HLLC-based model:

| case | local-safe HLLC HCFL | + 4-step fine-tune |
|---|---:|---:|
| Sod | 0.0765 | **0.0720** |
| Lax | 0.1958 | **0.1940** |
| collision | 0.0952 | **0.0903** |
| strong pressure | **0.2498** | 0.2514 |
| near vacuum | **0.1791** | 0.1905 |

Thus multi-step fine-tuning is a modest, non-uniform optimization improvement.
The structural gains come primarily from the HLLC proposal and local limiter.

## Current interpretation

The strongest current claim is not that HCFL uniformly beats a high-order
classical solver. MUSCL-HLLC remains substantially stronger on Sod, Lax, and
strong-pressure canonical tests. However:

- the learned HLLC-based HCFL beats first-order HLLC on Sod and strong-pressure;
- it strongly beats MUSCL-HLLC on the collision test;
- it comes close to MUSCL-HLLC on near-vacuum expansion;
- it preserves exact conservative flux form;
- interface Tadmor feasibility is enforced by construction;
- the local admissibility limiter prevents catastrophic negative density or
  pressure while intervening only where needed.

This supports the narrower framing of HCFL as a **certified learned coarse
solver** rather than a universal replacement for classical high-order methods.

## Next decision

Do not move to 2D yet. The remaining technical question is whether the
fully-discrete entropy safeguard can also be localized, avoiding the remaining
trajectory-wise scalar beta. If that can be done without losing the guarantee,
then the 1D method is structurally mature enough for 2D.
