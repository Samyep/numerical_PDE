# Allen-Cahn truncated-MLP source notes

## Primary theory source

Christian Beck, Fabian Hornung, Martin Hutzenthaler, Arnulf Jentzen, and
Thomas Kruse, *Overcoming the curse of dimensionality in the numerical
approximation of Allen-Cahn partial differential equations via truncated
full-history recursive multilevel Picard approximations*, Journal of
Numerical Mathematics 28(4), 197-222 (2020), DOI
[`10.1515/jnma-2019-0074`](https://doi.org/10.1515/jnma-2019-0074),
[arXiv:1907.06729](https://arxiv.org/abs/1907.06729).

The inspected arXiv PDF has SHA-256
`521ed4d065ed24349f629f80794ca094ffb1e87664afd5341882a60bfd10aec5`.

This paper is a theory and complexity paper. It defines and analyzes
truncated full-history recursive MLP for locally Lipschitz coercive reaction
nonlinearities, but it does **not** contain a concrete finite-budget numerical
Allen-Cahn benchmark or a table of simulation results.

### Equation and time convention

The Allen-Cahn specialization in Corollary 5.3 is written in forward form:

```text
partial_t u_d(t,x) = Delta_x u_d(t,x) + u_d(t,x) - u_d(t,x)^3,
```

with data prescribed at `t=0`. Its diffusion samples are
`x + sqrt(2)(W_t-W_s)`. Earlier Setting 3.1 also states the time-reversed
terminal-value convention with generator `f_r` and standard Brownian
diffusion; Proposition 3.5 explicitly relates the two orientations.

### Truncated driver and radius schedule

Setting 3.1, equation (76), defines

```text
f_r(t,x,u) = f(t,x,min(r,max(-r,u))).
```

Every nonlinear evaluation in the recursive formula (77) uses this truncated
driver. Corollaries 5.1-5.3 allow a radius sequence `rho_M` satisfying

```text
rho_M -> infinity,
limsup rho_M / log(log M) < infinity
```

asymptotically; the finitely many small indices can be assigned separately.
The explicit introductory specialization uses

```text
rho_M = log(1 + log M).
```

The algorithmic subscript is important: the recursion uses `f_M`, so this is
a sampling-parameter radius, not a dimension radius. Corollary 5.3 then takes
`M=n` for the headline complexity statement.

For `f(u)=u-u^3`, Beck et al. verify local Lipschitz continuity and the
coercivity inequality `u f(u)=u^2-u^4 <= 1+u^2`. Their conclusion is that for
every positive slack exponent there are levels achieving RMS error at most
`epsilon` with work bounded by `c_delta d epsilon^(-(2+delta))` after
renaming the slack parameter. This is an asymptotic existence/complexity
statement, not a prescription that the finite-budget numerical radius must be
one.

## Published numerical companion

Sebastian Becker, Ramon Braunwarth, Martin Hutzenthaler, Arnulf Jentzen, and
Philippe von Wurstemberger, *Numerical simulations for full history recursive
multilevel Picard approximations for systems of high-dimensional partial
differential equations*, Communications in Computational Physics 28(5),
2109-2138 (2020), DOI
[`10.4208/cicp.OA-2020-0130`](https://doi.org/10.4208/cicp.OA-2020-0130),
[arXiv:2005.10206](https://arxiv.org/abs/2005.10206).

The inspected arXiv PDF has SHA-256
`bff46fd1c1b99ceac7a0204bcf97eb9ea4f8d93b18739307cd4bd1fef925060f`.

Section 3.1 and the accompanying C++ code provide the concrete benchmark used
in this recovery study:

```text
T = 1,
partial_t u + Delta u + u - u^3 = 0,
u(T,x) = g(x) = 1 / (2 + (2/5)||x||^2),
X_{t,s}^x = x + sqrt(2)(W_s-W_t),
x_root = 0,
d in {10,100,1000},
n in {1,...,8}, M=n,
fixed truncation radius r=4.
```

The random time is uniform on the remaining interval. The published table
reports one realization at each cell and estimates relative L2 error from five
independent runs. The unknown exact root values are themselves approximated:

| d | Deep-splitting reference | Five-run `V_(8,8,4)` reference |
| ---: | ---: | ---: |
| 10 | 0.29614 | 0.29555 |
| 100 | 0.03376 | 0.03373 |
| 1000 | 0.00339 | 0.00340 |

These are published numerical references, not analytic truth. The companion
paper's `r=4` is justified by the a priori estimate in its equation (6), which
gives a solution magnitude no larger than `sqrt(5e)/2 <= 4`. It is distinct
from the asymptotically growing theorem schedule above.

## Solution-side invariant interval

For this specific companion instance, the terminal function satisfies
`0 < g <= 1/2`. With `v(s,x)=u(1-s,x)`, the PDE becomes

```text
partial_s v = Delta v + v - v^3.
```

The constants zero and one are stationary sub- and supersolutions because
`f(0)=f(1)=0`. The parabolic comparison principle therefore gives
`0 <= u(t,x) <= 1`. Consequently `[0,1]` is a rigorous, tighter
solution-side invariant interval for this particular PDE instance.

This interval is not Beck et al.'s general theorem schedule and is not the
fixed `[-4,4]` truncation used in the published numerical table. The recovery
experiment records all three concepts separately.

## Novelty boundary

Beck et al. established scalar truncated MLP and its local-to-global
Lipschitz modification for reaction nonlinearities. The containment result in
this repository does not claim those ideas as new. It shows that the generic
Samplewise IR map reduces exactly to their scalar operation when the feasible
set is `[-r,r] x R^d`. Any broader novelty claim must concern genuinely
structured value-gradient geometry in gradient-dependent PDEs, not scalar
Allen-Cahn truncation.
