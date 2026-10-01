# Fully-discrete entropy localization ablation

The 1D Euler method already has:

- shared conservative interface fluxes;
- hard Tadmor projection at every interface;
- local interface-wise admissibility limiting for \\(\\rho>0\\) and \\(p>0\\);
- a final trajectory-wise scalar \\(\\beta\\in[0,1]\\) that guarantees
  \\(\\sum_i \\eta(U_i^{n+1})\\le \\sum_i\\eta(U_i^n)\\).

We tested whether the last global scalar could be replaced by a purely local
fully-discrete entropy limiter.

## Attempt 1: cell-wise entropy budgets

For the entropy-projected low-order flux \\(F^{lo}\\), define

\\[
Q^{lo}_{i+1/2}
=
\\{v\\}_{i+1/2}^T F^{lo}_{i+1/2}
-
\\{\\psi\\}_{i+1/2},
\\]

and the forward-Euler cell budget

\\[
B_i
=
\\eta(U_i^n)
-
\\lambda
\\left(
Q^{lo}_{i+1/2}-Q^{lo}_{i-1/2}
\\right).
\\]

A local correction can then be constrained to keep
\\(\\eta(U_i^{n+1})\\le B_i\\).

This is rigorous when the low-order forward-Euler step itself satisfies the
budget, but it is too restrictive in practice. On the canonical seed-0 Euler
tests it limited roughly 44--94% of interfaces and degraded accuracy:

| case | strict local-budget NRMSE |
|---|---:|
| Sod | 0.0715 |
| Lax | 0.2377 |
| collision | 0.2006 |
| strong pressure | 0.2716 |
| near vacuum | 0.2245 |

The low-order local entropy condition can also require a substantially smaller
time step on extreme random states.

## Attempt 2: sequential interface spending of global entropy slack

Start from the low-order update and its total entropy slack

\\[
S
=
\\sum_i\\eta(U_i^n)
-
\\sum_i\\eta(U_i^{lo})
\\ge0.
\\]

Process interfaces sequentially. A correction at interface \\(i+1/2\\)
changes only cells \\(i\\) and \\(i+1\\). Choose the largest
\\(\\beta_{i+1/2}\\in[0,1]\\) that preserves both Euler admissibility
and the remaining global entropy slack.

This is conservative and fully-discrete entropy safe, but order dependent. On
seed 0 it gave:

| case | sequential-slack NRMSE |
|---|---:|
| Sod | 0.0745 |
| Lax | 0.1960 |
| collision | 0.1872 |
| strong pressure | 0.2418 |
| near vacuum | 0.2460 |

It is less restrictive than cell-wise budgets but still materially worse than
the current method on collision and near-vacuum expansion.

## Why we keep the global scalar entropy safeguard

We directly ablated the final global entropy scalar after the local
admissibility limiter.

Three-seed HLLC-HCFL results:

| case | local admissibility only | + global entropy scalar | mean global beta |
|---|---:|---:|---:|
| Sod | 0.07653 | 0.07653 | 1.0000 |
| Lax | 0.19578 | 0.19578 | 1.0000 |
| collision | 0.09474 | 0.09518 | 0.9984 |
| strong pressure | 0.24975 | 0.24975 | 1.0000 |
| near vacuum | 0.17587 | 0.17909 | 0.9828 |

Thus the global fully-discrete entropy safeguard:

- is exactly inactive on Sod, Lax, and strong-pressure cases;
- changes collision accuracy by about 0.5%;
- changes near-vacuum accuracy by about 1.8%;
- converts large positive one-step total-entropy excursions into guaranteed
  nonincrease.

This is a better accuracy/guarantee tradeoff than either local entropy
construction tested so far.

## Current design decision

Keep:

1. hard interface Tadmor projection;
2. local interface-wise admissibility limiter;
3. rare global fully-discrete entropy scalar.

Do **not** force the fully-discrete entropy condition to be cell-local in the
current version. The local alternatives are mathematically clean but
unnecessarily restrictive.

This closes the main 1D limiter-design question and supports moving to 2D.
