# Frozen benchmark definitions

All equations use the stated diffusion convention. Tangent-Actor records `z = sigma^T grad(u)` explicitly; baseline equations are supplied in each method's native variables without adding Tangent-Actor features.

| Key | Equation | Terminal | T | Reference | Dimensions |
|---|---|---|---:|---|---|
| `hjb_log` | `u_t + Delta u - ||grad u||² = 0` | `log((1+||x||²)/2)` | 1 | Hopf–Cole 1D quadrature | 10,20,50,100,200,400 |
| `ou_gmm` | `u_t + 0.5 Delta u -(mu-x)·grad u -0.5||grad u||²-d=0` | frozen five-component GMM negative log density | 1 | exact evolving GMM | 10,20,50,100,200,400 |
| `joint_concave` | `u_t+Delta u+u-u³-0.5||sqrt(2)grad u||²=0` | `0.5+0.15(sin(b1·x)+sin(b2·x))` | 0.5 | independent 2D spectral solver | 10,20,50,100,200,400 |
| `joint_convex` | `u_t+Delta u+u³-u+0.5||sqrt(2)grad u||²=0` | same two-direction terminal | 0.5 | independent 2D spectral solver | 10,20,50,100,200,400 |
| `nonconvex_coupled` | `u_t+Delta u+0.2 sin(2pi u)-0.5||sqrt(2)grad u+u e||²=0` | `0.5+0.3 sin(b1·x)cos(b2·x)` | 0.5 | independent 2D Fourier solver | 10,20,50,100,200,400 |
| `mlp_failure` | `u_t+0.5 Delta u+||grad(u)_{2:d}||=0` | `|x1|` | 1 | exact Gaussian absolute moment | 20,50,100,200,400,1000 |
| `allen_cahn` | `u_t+Delta u+u-u³=0` | `1/(2+0.4||x||²)` | 0.3 | dimension-specific radial PDE | 10,20,50,100,200,400 |

The dense basis is generated deterministically and checked for orthonormality. Every learned model receives `(t,x)` in the full ambient dimension. Reduced coordinates are used only by independent references.

The default whole-function bank contains 960 held-out states. HJB-Log uses `t~U[0,T), x~N(0,I)`. OU/GMM reproduces `x|t~N(0,(4+t)I)`. The joint and nonconvex PDEs use `t~U[0,T), x~N(0,I)`. Allen–Cahn uses the natural Brownian marginal `x|t~N(0,2tI)`. Exact historical 100D banks are retained when the audit verifies identity; each bank's SHA-256 is recorded.

For the nonsmooth MLP example, Tangent-Actor uses the explicitly declared minimum-norm subgradient of the inactive norm: `b(q)=q_inactive/||q_inactive||` away from zero and `b(0)=0`. This is an adapted extension and supplies no global smooth-tangent or convergence theorem. Any external solver using an automatic subgradient is also labeled adapted.
