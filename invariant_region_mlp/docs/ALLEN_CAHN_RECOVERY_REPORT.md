# Allen-Cahn truncated-MLP recovery report

## Outcome

**VERDICT B: Mathematical/code equivalence is verified, but the tested source-faithful Allen-Cahn configurations do not activate truncation enough to show a numerical rescue.**

The experiment proves containment and verifies it pathwise. Direct Beck-style `f(P_r(u))` and generic Samplewise projection onto `[-r,r] x R^d` agree exactly in two independently written recursions. The published Allen-Cahn convergence regime is also reproduced numerically. However, Raw MLP never leaves even the tighter certified interval `[0,1]` anywhere in the completed source, pilot, or 30-repetition final grids. Consequently Raw, Beck-Truncated, same-radius Interval IR, and certified `[0,1]` IR are identical in every paired production run. Claiming an active rescue would be unsupported.

The detailed source audit is in `docs/ALLEN_CAHN_SOURCE_NOTES.md`. Beck et al. (2020) provide the truncation theory and complexity result but no finite-budget numerical table. The concrete benchmark comes from Becker et al. (2020), their published companion simulation paper.

## Exact source setting

The reproduced terminal-value problem is

```text
partial_t u(t,x) + Delta u(t,x) + u(t,x) - u(t,x)^3 = 0,
T = 1,
u(T,x) = 1 / (2 + (2/5)||x||^2),
X_(t,s)^x = x + sqrt(2)(W_s-W_t),
x_root = 0.
```

The companion benchmark uses `d in {10,100,1000}`, `n=M in {1,...,8}`, a uniform random time, and fixed radius `r=4`. This study reproduced `n=M=1,...,5` with five independent repetitions per cell, matching the paper's five-run error protocol. The published `V_(8,8,4)` reference values (`0.29555`, `0.03373`, `0.00340`) are numerical rather than exact.

Beck's theory separately permits a growing `rho_M` with `rho_M -> infinity` and `rho_M=O(log log M)`; its explicit example is `log(1+log M)`. That theorem schedule, the companion's fixed `r=4`, and this instance's solution-side invariant interval `[0,1]` are not conflated.

## Mathematical containment proposition

Let `H=R x R^d`, `Y=(u,z)`, and `C_r=[-r,r] x R^d`. Because `C_r` is a Cartesian product of a closed interval and the full gradient space, Euclidean distance separates and

```text
Pi_(C_r)(u,z) = (P_r(u),z),
P_r(u) = min(r,max(-r,u)).
```

For a reaction driver independent of `z`, Samplewise IR therefore evaluates

```text
bar f(u,z) = f(P_r(u)),
```

which is exactly Beck's truncated driver.

For pathwise equality, fix the random tree, random-time variables, diffusion increments, sample allocation, and radius schedule. At Picard depth zero both recursions return zero. Assume equality at every depth below `n`. At every level in the depth-`n` estimator, the fine and coarse child values are equal by the induction hypothesis; applying the identities above gives identical fine and coarse generator values. Terminal samples, level weights, additions, and subtractions are then identical. Hence the complete depth-`n` estimates are equal on that random tree. Induction proves `U_IR(n,M,r)=U_Beck(n,M,r)` pathwise, not only in law.

## Independent code-path verification

The direct implementation performs its own scalar clamp and cubic evaluation. The IR implementation runs a separately written recursion, calls the generic interval projector, and then calls the original reaction. They share only equation primitives and an immutable keyed random tree.

- dedicated equivalence grid: 240 paired configurations;
- dimensions: `1,10,100,1000`;
- depths/sample sizes: `(0,2),(1,2),(2,2),(3,2),(3,3)`;
- three seeds, zero and nonzero ramp points;
- fixed `r=4` and `rho_M=log(1+log M)` tracks;
- maximum root discrepancy: `0`;
- exact root matches: `240/240`;
- maximum saved-correction discrepancy: `0`;
- exact saved-correction matches: `9408/9408`.

Across all source/pilot/final tasks as well, the maximum Beck-vs-IR discrepancy is `0` and all draw fingerprints match.

## Published-regime reproduction

| d | n | M | reps | mean value | MAE | random draws | max reused u |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 10 | 1 | 1 | 5 | 0.124723 | 0.170827 | 21 | 0 |
| 10 | 2 | 2 | 5 | 0.210621 | 0.084929 | 190 | 0.139602 |
| 10 | 3 | 3 | 5 | 0.273542 | 0.030285 | 2688 | 0.23676 |
| 10 | 4 | 4 | 5 | 0.286311 | 0.01771 | 5.1772e+04 | 0.360778 |
| 10 | 5 | 5 | 5 | 0.295944 | 0.004549 | 1.2770e+06 | 0.354488 |
| 100 | 1 | 1 | 5 | 0.011083 | 0.022647 | 201 | 0 |
| 100 | 2 | 2 | 5 | 0.025026 | 0.008704 | 1810 | 0.015156 |
| 100 | 3 | 3 | 5 | 0.03092 | 0.003746 | 2.5638e+04 | 0.024543 |
| 100 | 4 | 4 | 5 | 0.032766 | 0.001592 | 4.9421e+05 | 0.028819 |
| 100 | 5 | 5 | 5 | 0.033438 | 0.000474 | 1.2197e+07 | 0.033357 |
| 1000 | 1 | 1 | 5 | 0.001254 | 0.002146 | 2001 | 0 |
| 1000 | 2 | 2 | 5 | 0.002487 | 0.000913 | 1.8010e+04 | 0.001295 |
| 1000 | 3 | 3 | 5 | 0.003014 | 0.000386 | 2.5514e+05 | 0.002417 |
| 1000 | 4 | 4 | 5 | 0.00331 | 0.000147 | 4.9186e+06 | 0.003052 |
| 1000 | 5 | 5 | 5 | 0.00335 | 6.6119e-05 | 1.2139e+08 | 0.00328 |

At `n=M=5`, the means are `0.295944` (d=10), `0.033438` (d=100), and `0.003350` (d=1000), close to the companion's numerical references `0.29555`, `0.03373`, and `0.00340`. The reference is not treated as analytic truth.

## Controlled depth/budget sweep and final cells

Stage B kept the PDE, terminal function, horizon, diffusion, and evaluation point unchanged. It explored 26 pilot cells covering high sampling at depth two and increasingly deep/under-sampled allocations through `(n,M)=(8,2)`. Nine representative cells were then frozen at 30 repetitions:

| d | n | M | reps | Raw MAE | Beck MAE | IR r=4 MAE | IR [0,1] MAE | max reused u |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 10 | 2 | 16 | 30 | 0.069767 | 0.069767 | 0.069767 | 0.069767 | 0.255201 |
| 10 | 5 | 3 | 30 | 0.007794 | 0.007794 | 0.007794 | 0.007794 | 0.365371 |
| 10 | 8 | 2 | 30 | 0.009748 | 0.009748 | 0.009748 | 0.009748 | 0.407293 |
| 100 | 3 | 3 | 30 | 0.003489 | 0.003489 | 0.003489 | 0.003489 | 0.025816 |
| 100 | 5 | 3 | 30 | 0.000895 | 0.000895 | 0.000895 | 0.000895 | 0.036736 |
| 100 | 6 | 2 | 30 | 0.001619 | 0.001619 | 0.001619 | 0.001619 | 0.037945 |
| 1000 | 3 | 3 | 30 | 0.000292 | 0.000292 | 0.000292 | 0.000292 | 0.002494 |
| 1000 | 5 | 2 | 30 | 0.000198 | 0.000198 | 0.000198 | 0.000198 | 0.003715 |
| 1000 | 6 | 2 | 30 | 9.4877e-05 | 9.4877e-05 | 9.4877e-05 | 9.4877e-05 | 0.003561 |

The largest pre-projection child value in all production artifacts is `0.407293`. The minimum is zero. There are no `[-4,4]` violations, no `[0,1]` violations, no activations, and no nonfinite values or generators. Thus this is a strong inactive sanity check, not an active stabilization result.

## Driver regularity

For `f(u)=u-u^3`, `f'(u)=1-3u^2`. On `[-r,r]`, the clipped driver is globally Lipschitz with sharp constant

```text
max(1, |1-3r^2|),
```

bounded by `1+3r^2` as in a standard local estimate. At `r=4` the sharp constant is `47` (the loose bound is 49); on `[0,1]` it is 2. A dense numerical slope check is included in the tests.

## Work and reproducibility

The production grid contains 1692 method-repetitions plus 240 dedicated equivalence pairs. Recorded parallel orchestration wall time is `156.625` seconds; summed task wall time is `1128.200` seconds. Work columns contain terminal evaluations, generator evaluations, recursive states, normal scalar draws, uniform draws, total stochastic samples, and per-root wall time.

No manuscript source was modified. The result establishes backward compatibility with scalar truncated MLP; it does not claim truncation or the local-to-global argument as new. The broader IR contribution, if used later, must concern structured value-gradient geometry for gradient-dependent nonlinear reuse.

**VERDICT B: Mathematical/code equivalence is verified, but the tested source-faithful Allen-Cahn configurations do not activate truncation enough to show a numerical rescue.**
