# Benchmark compatibility audit

This audit was completed before GPU training.  All signs use
`z = sigma^T grad(u)`.

| Benchmark | Exact equation used | Tangent eligibility | Reference | Decision |
|---|---|---|---|---|
| DeepBSDE / SDGD HJB-Log | `u_t + Delta u - ||grad u||^2 = 0`, `T=1`, `g=log((1+||x||^2)/2)` | Globally concave `f(z)=-||z||^2/2`; Gradient Actor only | Independent Hopf-Cole 1D quadrature | Include |
| DPI OU/GMM HJB | `u_t + alpha/2 Delta u - F.grad u - alpha/2||grad u||^2 - d theta=0`, `F=theta(mu-x)`, `alpha=theta=1`, `T=1` | Known linear drift plus globally concave quadratic gradient driver | Exact evolving five-component GMM | Include |
| DeepBSDE Allen-Cahn | `u_t + Delta u + u-u^3=0`, `T=0.3`, `g=1/(2+0.4||x||^2)` | `f''(u)=-6u<=0` on the comparison-principle invariant interval `[0,1]` | Independent squared-radius PDE with mesh refinement and Richardson audit | Include |
| SDGD HJB-Rosenbrock | Same concave HJB driver with a Rosenbrock terminal | Driver is compatible | Official repository omits generated reference files; its generator uses only 1,000 effective MC paths per point and gives no reference-error audit | Exclude from compute study |

For the DPI equation, the solver places `-F.grad(u)` in the known forward
generator and tangent-linearizes only `-||z||^2/2-d theta`.  Expanding the
controlled generator recovers the exact source equation; no coefficient is
changed.

The DPI source has a documentation inconsistency.  Its display string writes a
constant proportional to `alpha/2 div(F)`, while `ff()` uses `-d theta`.  The
latter is required by the repository's analytic evolving-GMM solution.  Direct
automatic differentiation gives a maximum PDE residual of about `1.2e-12`
for the implemented equation, recorded in `reference_accuracy_audit.json`.

The HJB-Log reference is independent of every learned method.  At the native
point it gives `u(0,0)=4.590161724604883`; quadrature orders 64, 128, and 256
and adaptive integration agree to floating-point precision.  The Allen-Cahn
reference is numerical rather than exact and is labelled accordingly.

