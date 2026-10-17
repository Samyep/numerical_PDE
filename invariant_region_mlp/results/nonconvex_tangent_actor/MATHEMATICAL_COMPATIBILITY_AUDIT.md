# Mathematical Compatibility Audit

The implementation uses the existing generic tangent-coefficient and weighted Feynman–Kac code. No PDE-specific linear solver is introduced.

| Driver | max tangent error | max Hessian error | f_uu>0 fraction | f_uu<0 fraction | invariant [0,1] | gate |
|---|---:|---:|---:|---:|---|---|
| coupled | 4.441e-16 | 0.000e+00 | 0.136 | 0.864 | True | PASS |
| primary | 4.441e-16 | 0.000e+00 | 0.457 | 0.543 | True | PASS |

The curvature audit samples `t ~ Uniform(0,0.8T)` and `x ~ N(0,I_d)` from the independent reference, so it excludes the terminal boundary where sign changes could be inherited trivially from the terminal condition.

| Driver | Value mean ± std | Value q05 / q50 / q95 | strong f_uu>0 | strong f_uu<0 | joint Hessian indefinite |
|---|---:|---:|---:|---:|---:|
| coupled | 0.4596 ± 0.0607 | 0.3207 / 0.4577 / 0.6325 | 0.124 | 0.849 | 0.207 |
| primary | 0.4955 ± 0.0673 | 0.3297 / 0.4959 / 0.6679 | 0.418 | 0.489 | 0.457 |

## Controlled curvature-strength references

All fractions use the same independent interior distribution and are frozen before Actor training.

| lambda | Value mean ± std | Value q05 / q50 / q95 | f_uu>0 | f_uu<0 | strong positive | strong negative |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.4927 ± 0.0864 | 0.3026 / 0.4930 / 0.6933 | 0.000 | 0.000 | 0.000 | 0.000 |
| 0.05 | 0.4935 ± 0.0811 | 0.3097 / 0.4938 / 0.6867 | 0.453 | 0.547 | 0.348 | 0.414 |
| 0.2 | 0.4955 ± 0.0673 | 0.3297 / 0.4959 / 0.6679 | 0.457 | 0.543 | 0.418 | 0.489 |
| 0.5 | 0.4979 ± 0.0477 | 0.3637 / 0.4982 / 0.6353 | 0.464 | 0.536 | 0.433 | 0.495 |

The interval check uses the forward-time form `w_tau = Delta w + f(w,sqrt(2) grad w)`: 0 is a subsolution and 1 is a supersolution for both drivers. The numerical reference independently remains inside this interval.

The joint Hessian is bounded on the invariant Value interval: the Value curvature has magnitude at most `4*pi^2*lambda + eta^2`, the Gradient block is `-I`, and the mixed block has norm `eta`. Therefore the tangent remainder is quadratically bounded whenever the Actor Gradient remains in a bounded neighborhood.

A local fixed-point argument would additionally require a stable linear solution operator (uniformly controlled exponential weights), a neighborhood in which the tangent-solution map is Lipschitz with product constant below one, bounded trajectory moments, and Actor/MC/discretization errors smaller than that neighborhood's radius. The experiment audits weight ESS, finite trajectories, fixed-point contraction, and step refinement as numerical evidence for these conditions. It does not prove that the map is globally contractive.

Accordingly, the local Taylor bound is used as a consistency mechanism only. It does not imply monotone policy improvement, a one-sided tangent bound, or global convergence for a nonconvex driver.
