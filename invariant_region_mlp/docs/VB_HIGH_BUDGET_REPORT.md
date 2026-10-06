# Active-gradient viscous-Burgers high-budget MLP study

Generated from the `main,high_budget_extension` repetition artifacts on 2026-10-06T05:51:25.500195+00:00. No manuscript source was modified.

## Benchmark and exact solution

The main experiment uses the published SCaSML equation with `T=0.5`, `sigma=sqrt(2)`, coordinatewise drift `mu=-(1/d+sigma^2/2)`, generator `f(u,z)=sigma*u*sum_i z_i`, and state convention `z=sigma*grad(u)`. The exact solution is `u*(t,x)=sigmoid(t+sum_i x_i)` and `z_i*=sigma*u*(1-u)`.

Writing `q=u(1-u)`, the analytic residual is

`q + d*mu*q + (sigma^2/2)d*q(1-2u) + d*sigma^2*u*q = 0`.

The unit test evaluates this identity in dimensions 2, 20, and 80 and requires maximum absolute residual below `2e-13`. A separate centered finite-difference check of `u_t`, every gradient coordinate, and the Laplacian requires residual below `2e-7`; both checks pass. The 100,000-point signal gate attained a maximum analytic residual of `1.42e-14`.

## Exact certificates and corrections

Because `0<u<1` and `0<=u(1-u)<=1/4`, every exact state lies in `[0,1] x [0,sigma/4]^d` and in the z-ball of radius `sigma*sqrt(d)/4`. Samplewise box, z-only box, and ball corrections are the corresponding Euclidean projections. Batch-IR first replaces negative z coordinates by zero, then applies one common sibling factor `alpha=min(1,(sigma/4)/max_{i,j}(z_ij)_+)`; u is clipped samplewise. Thus Batch-IR is box-feasible and identity on feasible batches, but it is explicitly not called a coordinatewise box projection.

## Published/public-repository discrepancy and solver provenance

The branch was created from `origin/ir-mlp-life-or-death` at commit `fd2eece949b3a72ab24986cb8cba00e4adcc863c`; its corrected full-history Funding/HJB implementations supplied the local recursion pattern. The [published report](https://2prime.github.io/scasml_techreport.pdf) specifies `sigma=sqrt(2)`. [Public SCaSML](https://github.com/Francis-Fan-create/SCaSML) commit `c2412b942b1720084accf9ac9651c6085e073919` instead returns `sigma=0.25` for this equation. The public full-history terminal estimator also divides a standard normal by `T-t`; this study uses the corrected EBL normalization `xi/sqrt(T-t)`. All scientific runs use float64, no SCaSML surrogate clipping, and Beta(1/2,1) time importance sampling. The algebraically zero level-0 generator summand is elided and recorded as such. The optional `sigma=0.25` diagnostic was not run, so no repository-sigma result is mixed into the main tables.

## Nonlinear-signal gate

| d | mean u* | median ||z*|| | mean |f*| | max |f*| | f=0 relative L2 |
|---:|---:|---:|---:|---:|---:|
| 20 | 0.5465 | 1.302 | 3.957 | 5.926 | 0.729 |
| 40 | 0.5392 | 1.554 | 6.61 | 11.85 | 0.8494 |
| 60 | 0.5341 | 1.611 | 8.667 | 17.78 | 0.8953 |
| 80 | 0.5327 | 1.598 | 10.48 | 23.7 | 0.9171 |

## Grid actually completed

The study was progressive rather than a single rectangular run:

- Signal gate: 100,000 fixed samples per dimension.
- Stage-1 pilot, `d=20`: the complete requested 16-point grid, namely `n=2, M={2,4,8,16,24,32}`, `n=3, M={2,3,4,6,8,10}`, `n=4, M={2,3,4}`, and `n=5, M=2`; 240 points, 2 paired repetitions, and 17 non-duplicate methods/controls. (`c=0/1` and `a=0/1` were evaluated through their exactly equivalent named methods.)
- Validation tuning, all four dimensions: 10 representative grid points `[(2,2),(2,8),(2,24),(2,32),(3,3),(3,6),(3,10),(4,2),(4,3),(5,2)]`, 240 validation points, and 3 paired repetitions. A second independent refinement added `c={0.375,0.625,0.875}` and `a={0.875,0.9}`.
- Held-out main run, all four dimensions: the same 10 representative grid points on 1,000 interior plus 200 boundary points. Low/high/deep headline cells have 10 paired repetitions; the other cells have 5. This produced 2,400 method-repetitions.
- High-budget consistency extension, all four dimensions: `(n,M)=(3,16)` and `(4,6)`, with Raw, certified box, `z=0`, and `f=0`, 5 paired repetitions on all 1,200 points. This added 160 method-repetitions.

Thus the final held-out aggregate covers dimensions `[20,40,60,80]`, 12 `(n,M)` configurations, and 2,560 method-repetitions. The complete pilot, tuning, and refinement manifests and raw scalar repetition files are committed alongside the final aggregate.

Fixed points are sampled on `t in [0,0.5]`, `x in [-0.5,0.5]^d`, with boundary points forcing one random coordinate to `+/-0.5`. Hyperparameters are selected only on the stratified 20% validation split; reported headline errors use the 80% test split.

## Work accounting

Every repetition records terminal g evaluations, nonlinear f evaluations, recursively evaluated states, transition and terminal stochastic samples, scalar normal draws, and wall time. `z=0` traverses the complete recursive tree and evaluates the zeroed generator; `f=0` deletes recursion and therefore has zero f calls. Total stochastic samples—not f calls alone—are used for comparisons involving `f=0`. Peak memory was not measured reliably under the shared eight-process runner and is therefore not reported rather than approximated.

## Validation-tuned controls

The full machine-readable choices, candidate losses, and aliases (`c=0` is `z=0`, `c=1` is Raw, tighter-box `a=0` is value-equivalent to `z=0`, and `a=1` is the certified box) are in `results/active_vb_high_budget/tuning_choices.json`. Selection used only the fixed 20% validation split. The held-out 80% test split was opened only after `c`, illegal `a<1`, and relaxed valid `a>1` had been frozen for each dimension/budget cell.

## Main empirical findings

### The nonlinear term is genuinely active

The mean exact generator magnitude grows from `3.96` at `d=20` to `10.48` at `d=80`, and deleting it incurs relative L2 error from `0.729` to `0.917`. This is qualitatively different from the earlier weak-signal HJB/funding life-or-death cases. In the stochastic runs, `z=0` and `f=0` are numerically identical under paired roots, but `z=0` pays the full recursive work. Their high-sample held-out floors are about `0.738`, `0.848`, `0.898`, and `0.920` for `d=20,40,60,80`.

### High-sampling crossover and consistency

The following held-out value relative L2 errors show the relevant high-sampling sequence. Each `M=10` cell has 10 repetitions; each expensive `M=16` cell has 5.

| d | Raw, n=3 M=10 | box, n=3 M=10 | z=0 floor | Raw, n=3 M=16 | box, n=3 M=16 |
|---:|---:|---:|---:|---:|---:|
| 20 | 0.466 | 0.469 | 0.737 | 0.439 | 0.461 |
| 40 | 0.696 | 0.669 | 0.848 | 0.677 | 0.670 |
| 60 | 0.808 | 0.780 | 0.898 | 0.793 | 0.781 |
| 80 | 0.859 | 0.834 | 0.920 | 0.847 | 0.833 |

Raw and certified box therefore approach one another as M rises at fixed `n=3`; both remain separated from the wrong-PDE floor. At `M=16`, Raw is better at `d=20`, while box is slightly better at `d=40,60,80`. The paired box-vs-Raw win counts are `0/5`, `5/5`, `5/5`, and `5/5`, respectively. This is the requested high-budget consistency behavior, although the remaining value error shows that this is not an asymptotic-accuracy study.

Picard depth without enough per-level sampling is not equivalent to high effective accuracy. At `(n,M)=(4,6)`, which has about 1.27 million f calls and 14.65 million stochastic samples per repetition, Raw still amplifies recursive noise while box remains stable:

| d | Raw | certified box | z=0/f=0 |
|---:|---:|---:|---:|
| 20 | 2.68 | 0.393 | 0.737 |
| 40 | 4.22 | 0.608 | 0.848 |
| 60 | 6.07 | 0.737 | 0.898 |
| 80 | 9.63 | 0.797 | 0.920 |

Box beats both Raw and the wrong-PDE floor in every one of the 20 paired repetitions in this deep/high-work diagnostic. Conversely, at the deliberately under-sampled `n=5,M=2` cells, validation-tuned illegal tighter boxes beat the exact box; those cells are a noise-dominated depth regime, not evidence of convergence.

### Certificate versus generic suppression

On the original 40 held-out dimension/budget cells, the primary untuned certified box beats `z=0`, `f=0`, validation-tuned constant shrinkage, and the validation-tuned illegal tighter box in 29 cells: `7/10`, `8/10`, `8/10`, and `6/10` from `d=20` through `80`. Including the 8 expensive consistency cells, the count is 37/48. At `n=3,M=10`, box beats both wrong-PDE baselines in all 40 paired repetitions; it beats Raw in `2/10`, `10/10`, `10/10`, and `10/10` repetitions across the four dimensions. It also beats the tuned illegal tighter control in all 40 repetitions. This is an empirical regime of advantage that constant or aggressive gradient suppression does not explain.

The exact tight bound is nevertheless not accuracy-optimal everywhere. A validation-selected *looser but still valid* outer box beats the tight box in 22/40 target-grid cells, commonly choosing `a=2` for high-M `n=2` cells. This result argues against interpreting every gain as “more shrinkage is better”: in a large part of the grid, less shrinkage improves value error. It also means the paper must not claim that the tightest exact certificate is universally optimal.

### Gradient, generator, and correction diagnostics

The positive result is specific to the value solution. On gradient relative L2, a certified method beats the suppression controls in only 3/40 target-grid cells. For example, at `n=3,M=16`, certified-box gradient errors are `0.919,1.274,1.498,1.676`, while `z=0/f=0` gives `0.749,0.880,0.922,0.946`. Projection substantially improves the very noisy Raw gradients (`1.49,2.12,2.64,2.97`) but does not beat the biased suppression baseline. Any paper claim must state this limitation explicitly.

Generator diagnostics do support the value-error mechanism. At `n=3,M=16`, Raw-to-box generator MSE changes from `5.50→5.29`, `18.49→17.82`, `34.34→33.30`, and `54.33→53.05`; the corresponding `z=0` MSEs are `9.71,26.83,44.81,68.06`. At `n=3,M=10`, samplewise box reduces the across-repetition nonlinear-correction variance by roughly a factor of five to six relative to Raw. Batch-IR reduces that variance still further, but its shared alpha introduces enough bias that it is never the overall value-error winner in the main cells.

Intermediate-state constraints remain active even at the largest budgets: the certified box pre-violation and activation rates are approximately `0.995–1.000`. For Batch-IR at `n=3,M=10`, mean alpha rises from `0.508` at `d=20` to `0.643` at `d=80`; the batch method is consequently a strong common contraction and is more biased than samplewise box. The u constraint matters mainly in deep/noisy settings: at `d=80,n=5,M=2`, full box has value error `0.950` versus `1.481` for z-only box, though the tuned tighter control is better than both.

## Complete aggregate result table

The exact full table is `results/active_vb_high_budget/work_summary.csv`. The compact table below gives held-out value error, gradient error, work, wall time, generator bias, and pre-correction box violations for every completed cell.

| d | n | M | method | reps | value rel L2 | grad rel L2 | f calls | samples | sec | gen bias | box viol. |
|---:|---:|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| 20 | 2 | 2 | f_zero | 10 | 0.7568 | 2.185 | 0 | 4800 | 0.009358 | NA | 0 |
| 20 | 2 | 8 | f_zero | 5 | 0.7388 | 0.9024 | 0 | 76800 | 0.02965 | NA | 0 |
| 20 | 2 | 24 | f_zero | 5 | 0.7375 | 0.766 | 0 | 691200 | 0.297 | NA | 0 |
| 20 | 2 | 32 | f_zero | 5 | 0.7376 | 0.7571 | 0 | 1.2288e+06 | 0.586 | NA | 0 |
| 20 | 3 | 3 | f_zero | 5 | 0.7398 | 1.083 | 0 | 32400 | 0.01591 | NA | 0 |
| 20 | 3 | 6 | f_zero | 5 | 0.7376 | 0.7963 | 0 | 259200 | 0.1082 | NA | 0 |
| 20 | 3 | 10 | f_zero | 10 | 0.7374 | 0.7571 | 0 | 1.2e+06 | 0.6548 | NA | 0 |
| 20 | 3 | 16 | f_zero | 5 | 0.7375 | 0.7489 | 0 | 4.9152e+06 | 3.029 | NA | 0 |
| 20 | 4 | 2 | f_zero | 5 | 0.7425 | 1.262 | 0 | 19200 | 0.01183 | NA | 0 |
| 20 | 4 | 3 | f_zero | 5 | 0.7388 | 0.8724 | 0 | 97200 | 0.0318 | NA | 0 |
| 20 | 4 | 6 | f_zero | 5 | 0.7372 | 0.7543 | 0 | 1.5552e+06 | 0.9501 | NA | 0 |
| 20 | 5 | 2 | f_zero | 10 | 0.7403 | 1.037 | 0 | 38400 | 0.01726 | NA | 0 |
| 20 | 2 | 2 | batch_box | 10 | 0.6975 | 2.384 | 2400 | 12000 | 0.03572 | -2.186 | 1 |
| 20 | 2 | 2 | box_factor_a0.875 | 10 | 0.6608 | 2.905 | 2400 | 12000 | 0.03101 | -1.865 | 1 |
| 20 | 2 | 2 | raw | 10 | 0.7812 | 3.997 | 2400 | 12000 | 0.03177 | -1.744 | 1 |
| 20 | 2 | 2 | sample_ball | 10 | 0.7254 | 2.507 | 2400 | 12000 | 0.03209 | -2.211 | 1 |
| 20 | 2 | 2 | sample_box | 10 | 0.6576 | 3.052 | 2400 | 12000 | 0.02897 | -1.798 | 1 |
| 20 | 2 | 2 | shrink_c0.25 | 10 | 0.7326 | 2.333 | 2400 | 12000 | 0.03188 | -2.299 | 1 |
| 20 | 2 | 2 | z_only_box | 10 | 0.6576 | 3.052 | 2400 | 12000 | 0.03162 | -1.798 | 1 |
| 20 | 2 | 2 | z_zero | 10 | 0.7568 | 2.185 | 2400 | 12000 | 0.03063 | -2.485 | 1 |
| 20 | 2 | 8 | batch_box | 5 | 0.6722 | 0.9486 | 9600 | 163200 | 0.09194 | -2.207 | 1 |
| 20 | 2 | 8 | box_factor_a0.9 | 5 | 0.601 | 1.262 | 9600 | 163200 | 0.08353 | -1.797 | 1 |
| 20 | 2 | 8 | box_factor_a2 | 5 | 0.5562 | 1.71 | 9600 | 163200 | 0.087 | -1.427 | 1 |
| 20 | 2 | 8 | raw | 5 | 0.6157 | 1.503 | 9600 | 163200 | 0.08505 | -1.674 | 1 |
| 20 | 2 | 8 | sample_ball | 5 | 0.6419 | 1.166 | 9600 | 163200 | 0.08611 | -1.942 | 1 |
| 20 | 2 | 8 | sample_box | 5 | 0.5932 | 1.317 | 9600 | 163200 | 0.08444 | -1.746 | 1 |
| 20 | 2 | 8 | z_only_box | 5 | 0.5932 | 1.317 | 9600 | 163200 | 0.08676 | -1.746 | 1 |
| 20 | 2 | 8 | z_zero | 5 | 0.7388 | 0.9024 | 9600 | 163200 | 0.08432 | -2.499 | 1 |
| 20 | 2 | 24 | batch_box | 5 | 0.6646 | 0.7353 | 28800 | 1.4112e+06 | 0.6471 | -2.18 | 0.9987 |
| 20 | 2 | 24 | box_factor_a0.9 | 5 | 0.5947 | 0.8592 | 28800 | 1.4112e+06 | 0.6377 | -1.768 | 0.9987 |
| 20 | 2 | 24 | box_factor_a2 | 5 | 0.5646 | 0.9876 | 28800 | 1.4112e+06 | 0.6336 | -1.536 | 0.9987 |
| 20 | 2 | 24 | raw | 5 | 0.5983 | 0.9396 | 28800 | 1.4112e+06 | 0.6217 | -1.663 | 0.9987 |
| 20 | 2 | 24 | sample_ball | 5 | 0.6055 | 0.8886 | 28800 | 1.4112e+06 | 0.6345 | -1.75 | 0.9987 |
| 20 | 2 | 24 | sample_box | 5 | 0.5885 | 0.8792 | 28800 | 1.4112e+06 | 0.6286 | -1.726 | 0.9987 |
| 20 | 2 | 24 | z_only_box | 5 | 0.5885 | 0.8792 | 28800 | 1.4112e+06 | 0.6306 | -1.726 | 0.9987 |
| 20 | 2 | 24 | z_zero | 5 | 0.7375 | 0.766 | 28800 | 1.4112e+06 | 0.6343 | -2.505 | 0.9987 |
| 20 | 3 | 3 | batch_box | 5 | 0.6581 | 1.297 | 28800 | 165600 | 0.1391 | -1.964 | 1 |
| 20 | 3 | 3 | box_factor_a0.9 | 5 | 0.5405 | 2.226 | 28800 | 165600 | 0.1227 | -1.537 | 1 |
| 20 | 3 | 3 | box_factor_a1.25 | 5 | 0.5183 | 3.127 | 28800 | 165600 | 0.1282 | -1.336 | 1 |
| 20 | 3 | 3 | raw | 5 | 1.275 | 10.82 | 28800 | 165600 | 0.1153 | -1.315 | 1 |
| 20 | 3 | 3 | sample_ball | 5 | 0.634 | 1.845 | 28800 | 165600 | 0.1258 | -1.842 | 1 |
| 20 | 3 | 3 | sample_box | 5 | 0.5292 | 2.471 | 28800 | 165600 | 0.1222 | -1.476 | 1 |
| 20 | 3 | 3 | shrink_c0.5 | 5 | 0.662 | 2.303 | 28800 | 165600 | 0.1174 | -1.817 | 1 |
| 20 | 3 | 3 | z_only_box | 5 | 0.5301 | 2.52 | 28800 | 165600 | 0.1193 | -1.476 | 1 |
| 20 | 3 | 3 | z_zero | 5 | 0.7398 | 1.083 | 28800 | 165600 | 0.1137 | -2.246 | 1 |
| 20 | 2 | 32 | batch_box | 5 | 0.6604 | 0.7163 | 38400 | 2.496e+06 | 1.268 | -2.165 | 0.9958 |
| 20 | 2 | 32 | box_factor_a0.9 | 5 | 0.5913 | 0.8099 | 38400 | 2.496e+06 | 1.203 | -1.759 | 0.9958 |
| 20 | 2 | 32 | box_factor_a2 | 5 | 0.5659 | 0.8959 | 38400 | 2.496e+06 | 1.22 | -1.564 | 0.9958 |
| 20 | 2 | 32 | raw | 5 | 0.5931 | 0.8691 | 38400 | 2.496e+06 | 1.262 | -1.662 | 0.9958 |
| 20 | 2 | 32 | sample_ball | 5 | 0.5975 | 0.8423 | 38400 | 2.496e+06 | 1.256 | -1.719 | 0.9958 |
| 20 | 2 | 32 | sample_box | 5 | 0.5855 | 0.8246 | 38400 | 2.496e+06 | 1.266 | -1.72 | 0.9958 |
| 20 | 2 | 32 | z_only_box | 5 | 0.5855 | 0.8246 | 38400 | 2.496e+06 | 1.253 | -1.72 | 0.9958 |
| 20 | 2 | 32 | z_zero | 5 | 0.7376 | 0.7571 | 38400 | 2.496e+06 | 1.23 | -2.504 | 0.9958 |
| 20 | 4 | 2 | batch_box | 5 | 0.6629 | 1.747 | 67200 | 247200 | 0.3216 | -1.727 | 1 |
| 20 | 4 | 2 | box_factor_a0.875 | 5 | 0.5662 | 3.552 | 67200 | 247200 | 0.2855 | -1.323 | 1 |
| 20 | 4 | 2 | raw | 5 | 148.8 | 1136 | 67200 | 247200 | 0.254 | -0.3867 | 1 |
| 20 | 4 | 2 | sample_ball | 5 | 0.6714 | 2.717 | 67200 | 247200 | 0.2622 | -1.673 | 1 |
| 20 | 4 | 2 | sample_box | 5 | 0.5772 | 4.228 | 67200 | 247200 | 0.2688 | -1.24 | 1 |
| 20 | 4 | 2 | shrink_c0.375 | 5 | 0.8291 | 4.777 | 67200 | 247200 | 0.2494 | -1.689 | 1 |
| 20 | 4 | 2 | z_only_box | 5 | 0.6054 | 4.752 | 67200 | 247200 | 0.2647 | -1.237 | 1 |
| 20 | 4 | 2 | z_zero | 5 | 0.7425 | 1.262 | 67200 | 247200 | 0.2488 | -2.026 | 1 |
| 20 | 3 | 6 | batch_box | 5 | 0.644 | 0.8493 | 100800 | 1.1736e+06 | 0.632 | -1.952 | 0.9995 |
| 20 | 3 | 6 | box_factor_a0.9 | 5 | 0.4993 | 1.337 | 100800 | 1.1736e+06 | 0.6233 | -1.5 | 0.9996 |
| 20 | 3 | 6 | box_factor_a1.5 | 5 | 0.4255 | 2.316 | 100800 | 1.1736e+06 | 0.6133 | -1.216 | 0.9996 |
| 20 | 3 | 6 | raw | 5 | 0.558 | 3.8 | 100800 | 1.1736e+06 | 0.5868 | -1.325 | 0.9996 |
| 20 | 3 | 6 | sample_ball | 5 | 0.557 | 1.311 | 100800 | 1.1736e+06 | 0.6078 | -1.679 | 0.9996 |
| 20 | 3 | 6 | sample_box | 5 | 0.4804 | 1.484 | 100800 | 1.1736e+06 | 0.6214 | -1.444 | 0.9996 |
| 20 | 3 | 6 | shrink_c0.875 | 5 | 0.5372 | 2.793 | 100800 | 1.1736e+06 | 0.6094 | -1.446 | 0.9996 |
| 20 | 3 | 6 | z_only_box | 5 | 0.4801 | 1.503 | 100800 | 1.1736e+06 | 0.6006 | -1.444 | 0.9996 |
| 20 | 3 | 6 | z_zero | 5 | 0.7376 | 0.7963 | 100800 | 1.1736e+06 | 0.6052 | -2.223 | 0.9995 |
| 20 | 4 | 3 | batch_box | 5 | 0.6483 | 1.096 | 190800 | 1.0728e+06 | 0.6037 | -1.714 | 0.9999 |
| 20 | 4 | 3 | box_factor_a0.9 | 5 | 0.4731 | 2.355 | 190800 | 1.0728e+06 | 0.5603 | -1.271 | 0.9999 |
| 20 | 4 | 3 | box_factor_a1.25 | 5 | 0.4554 | 3.797 | 190800 | 1.0728e+06 | 0.5561 | -1.065 | 1 |
| 20 | 4 | 3 | raw | 5 | 33.45 | 342.8 | 190800 | 1.0728e+06 | 0.5385 | -0.9077 | 1 |
| 20 | 4 | 3 | sample_ball | 5 | 0.5835 | 1.968 | 190800 | 1.0728e+06 | 0.5535 | -1.559 | 0.9999 |
| 20 | 4 | 3 | sample_box | 5 | 0.4567 | 2.736 | 190800 | 1.0728e+06 | 0.5691 | -1.209 | 1 |
| 20 | 4 | 3 | shrink_c0.375 | 5 | 0.662 | 1.905 | 190800 | 1.0728e+06 | 0.56 | -1.658 | 0.9999 |
| 20 | 4 | 3 | z_only_box | 5 | 0.4622 | 2.935 | 190800 | 1.0728e+06 | 0.5437 | -1.208 | 1 |
| 20 | 4 | 3 | z_zero | 5 | 0.7388 | 0.8724 | 190800 | 1.0728e+06 | 0.5223 | -1.988 | 0.9999 |
| 20 | 3 | 10 | batch_box | 10 | 0.6321 | 0.7511 | 264000 | 5.172e+06 | 3.092 | -1.925 | 0.9966 |
| 20 | 3 | 10 | box_factor_a0.9 | 10 | 0.4884 | 1.028 | 264000 | 5.172e+06 | 3.031 | -1.466 | 0.9972 |
| 20 | 3 | 10 | box_factor_a2 | 10 | 0.4002 | 2.226 | 264000 | 5.172e+06 | 2.993 | -1.098 | 0.9974 |
| 20 | 3 | 10 | raw | 10 | 0.4662 | 2.206 | 264000 | 5.172e+06 | 3.014 | -1.307 | 0.9972 |
| 20 | 3 | 10 | sample_ball | 10 | 0.5175 | 1.114 | 264000 | 5.172e+06 | 3.091 | -1.553 | 0.9972 |
| 20 | 3 | 10 | sample_box | 10 | 0.4694 | 1.126 | 264000 | 5.172e+06 | 3.053 | -1.414 | 0.9973 |
| 20 | 3 | 10 | z_only_box | 10 | 0.4688 | 1.137 | 264000 | 5.172e+06 | 3.112 | -1.414 | 0.9973 |
| 20 | 3 | 10 | z_zero | 10 | 0.7374 | 0.7571 | 264000 | 5.172e+06 | 3.079 | -2.198 | 0.9953 |
| 20 | 5 | 2 | batch_box | 10 | 0.6595 | 1.56 | 302400 | 1.1016e+06 | 0.9128 | -1.508 | 1 |
| 20 | 5 | 2 | box_factor_a0.75 | 10 | 0.5525 | 3.037 | 302400 | 1.1016e+06 | 0.7968 | -1.193 | 1 |
| 20 | 5 | 2 | raw | 10 | 1.804e+07 | 1.89e+08 | 302400 | 1.1016e+06 | 0.7492 | -1553 | 1 |
| 20 | 5 | 2 | sample_ball | 10 | 0.6564 | 2.969 | 302400 | 1.1016e+06 | 0.7649 | -1.434 | 1 |
| 20 | 5 | 2 | sample_box | 10 | 0.5784 | 4.721 | 302400 | 1.1016e+06 | 0.7534 | -1.027 | 1 |
| 20 | 5 | 2 | shrink_c0.25 | 10 | 0.7053 | 2.19 | 302400 | 1.1016e+06 | 0.7621 | -1.578 | 1 |
| 20 | 5 | 2 | z_only_box | 10 | 0.6456 | 5.683 | 302400 | 1.1016e+06 | 0.769 | -1.024 | 1 |
| 20 | 5 | 2 | z_zero | 10 | 0.7403 | 1.037 | 302400 | 1.1016e+06 | 0.7545 | -1.799 | 1 |
| 20 | 3 | 16 | raw | 5 | 0.4391 | 1.489 | 652800 | 2.06016e+07 | 13.64 | -1.297 | 0.9949 |
| 20 | 3 | 16 | sample_box | 5 | 0.4612 | 0.9187 | 652800 | 2.06016e+07 | 13.65 | -1.393 | 0.9949 |
| 20 | 3 | 16 | z_zero | 5 | 0.7375 | 0.7489 | 652800 | 2.06016e+07 | 13.64 | -2.18 | 0.9896 |
| 20 | 4 | 6 | raw | 5 | 2.679 | 32.66 | 1.2672e+06 | 1.4652e+07 | 10.34 | -1.017 | 0.9988 |
| 20 | 4 | 6 | sample_box | 5 | 0.3927 | 1.54 | 1.2672e+06 | 1.4652e+07 | 10.42 | -1.17 | 0.9987 |
| 20 | 4 | 6 | z_zero | 5 | 0.7372 | 0.7543 | 1.2672e+06 | 1.4652e+07 | 10.4 | -1.949 | 0.9976 |
| 40 | 2 | 2 | f_zero | 10 | 0.8573 | 2.743 | 0 | 4800 | 0.0118 | NA | 0 |
| 40 | 2 | 8 | f_zero | 5 | 0.8477 | 1.087 | 0 | 76800 | 0.04245 | NA | 0 |
| 40 | 2 | 24 | f_zero | 5 | 0.8477 | 0.903 | 0 | 691200 | 0.784 | NA | 0 |
| 40 | 2 | 32 | f_zero | 5 | 0.8477 | 0.8909 | 0 | 1.2288e+06 | 1.293 | NA | 0 |
| 40 | 3 | 3 | f_zero | 5 | 0.8489 | 1.329 | 0 | 32400 | 0.02606 | NA | 0 |
| 40 | 3 | 6 | f_zero | 5 | 0.8479 | 0.9437 | 0 | 259200 | 0.2458 | NA | 0 |
| 40 | 3 | 10 | f_zero | 10 | 0.8478 | 0.8913 | 0 | 1.2e+06 | 1.339 | NA | 0 |
| 40 | 3 | 16 | f_zero | 5 | 0.8476 | 0.8799 | 0 | 4.9152e+06 | 5.639 | NA | 0 |
| 40 | 4 | 2 | f_zero | 5 | 0.8498 | 1.566 | 0 | 19200 | 0.01612 | NA | 0 |
| 40 | 4 | 3 | f_zero | 5 | 0.8486 | 1.051 | 0 | 97200 | 0.05104 | NA | 0 |
| 40 | 4 | 6 | f_zero | 5 | 0.848 | 0.888 | 0 | 1.5552e+06 | 1.845 | NA | 0 |
| 40 | 5 | 2 | f_zero | 10 | 0.8489 | 1.27 | 0 | 38400 | 0.02333 | NA | 0 |
| 40 | 2 | 2 | batch_box | 10 | 0.8223 | 2.964 | 2400 | 12000 | 0.03791 | -3.396 | 1 |
| 40 | 2 | 2 | box_factor_a0.9 | 10 | 0.7971 | 3.842 | 2400 | 12000 | 0.03443 | -2.938 | 1 |
| 40 | 2 | 2 | box_factor_a1.5 | 10 | 0.8035 | 4.857 | 2400 | 12000 | 0.03294 | -2.577 | 1 |
| 40 | 2 | 2 | raw | 10 | 0.8865 | 5.02 | 2400 | 12000 | 0.0328 | -2.866 | 1 |
| 40 | 2 | 2 | sample_ball | 10 | 0.8436 | 3.047 | 2400 | 12000 | 0.03134 | -3.461 | 1 |
| 40 | 2 | 2 | sample_box | 10 | 0.7965 | 4.015 | 2400 | 12000 | 0.03413 | -2.872 | 1 |
| 40 | 2 | 2 | shrink_c0.5 | 10 | 0.844 | 3.455 | 2400 | 12000 | 0.0303 | -3.28 | 1 |
| 40 | 2 | 2 | z_only_box | 10 | 0.7965 | 4.015 | 2400 | 12000 | 0.03201 | -2.872 | 1 |
| 40 | 2 | 2 | z_zero | 10 | 0.8573 | 2.743 | 2400 | 12000 | 0.03204 | -3.695 | 1 |
| 40 | 2 | 8 | batch_box | 5 | 0.8096 | 1.159 | 9600 | 163200 | 0.1122 | -3.415 | 1 |
| 40 | 2 | 8 | box_factor_a0.9 | 5 | 0.7576 | 1.635 | 9600 | 163200 | 0.1227 | -2.881 | 1 |
| 40 | 2 | 8 | box_factor_a2 | 5 | 0.7269 | 2.356 | 9600 | 163200 | 0.1167 | -2.352 | 1 |
| 40 | 2 | 8 | raw | 5 | 0.7791 | 1.898 | 9600 | 163200 | 0.1074 | -2.783 | 1 |
| 40 | 2 | 8 | sample_ball | 5 | 0.8014 | 1.322 | 9600 | 163200 | 0.1182 | -3.209 | 1 |
| 40 | 2 | 8 | sample_box | 5 | 0.7525 | 1.715 | 9600 | 163200 | 0.1189 | -2.815 | 1 |
| 40 | 2 | 8 | z_only_box | 5 | 0.7525 | 1.715 | 9600 | 163200 | 0.1217 | -2.815 | 1 |
| 40 | 2 | 8 | z_zero | 5 | 0.8477 | 1.087 | 9600 | 163200 | 0.1052 | -3.703 | 1 |
| 40 | 2 | 24 | batch_box | 5 | 0.8075 | 0.9117 | 28800 | 1.4112e+06 | 1.678 | -3.387 | 1 |
| 40 | 2 | 24 | box_factor_a0.9 | 5 | 0.7569 | 1.111 | 28800 | 1.4112e+06 | 1.691 | -2.86 | 1 |
| 40 | 2 | 24 | box_factor_a2 | 5 | 0.7323 | 1.346 | 28800 | 1.4112e+06 | 1.697 | -2.503 | 1 |
| 40 | 2 | 24 | raw | 5 | 0.7679 | 1.202 | 28800 | 1.4112e+06 | 1.655 | -2.763 | 1 |
| 40 | 2 | 24 | sample_ball | 5 | 0.7782 | 1.082 | 28800 | 1.4112e+06 | 1.654 | -2.952 | 1 |
| 40 | 2 | 24 | sample_box | 5 | 0.7524 | 1.143 | 28800 | 1.4112e+06 | 1.662 | -2.805 | 1 |
| 40 | 2 | 24 | z_only_box | 5 | 0.7524 | 1.143 | 28800 | 1.4112e+06 | 1.662 | -2.805 | 1 |
| 40 | 2 | 24 | z_zero | 5 | 0.8477 | 0.903 | 28800 | 1.4112e+06 | 1.666 | -3.692 | 1 |
| 40 | 3 | 3 | batch_box | 5 | 0.799 | 1.572 | 28800 | 165600 | 0.2036 | -2.919 | 1 |
| 40 | 3 | 3 | box_factor_a0.9 | 5 | 0.7118 | 3.119 | 28800 | 165600 | 0.1789 | -2.361 | 1 |
| 40 | 3 | 3 | box_factor_a1.25 | 5 | 0.7077 | 4.57 | 28800 | 165600 | 0.1807 | -2.102 | 1 |
| 40 | 3 | 3 | raw | 5 | 1.366 | 22.58 | 28800 | 165600 | 0.1693 | -2.155 | 1 |
| 40 | 3 | 3 | sample_ball | 5 | 0.7959 | 1.945 | 28800 | 165600 | 0.1797 | -2.855 | 1 |
| 40 | 3 | 3 | sample_box | 5 | 0.7062 | 3.513 | 28800 | 165600 | 0.1888 | -2.284 | 1 |
| 40 | 3 | 3 | shrink_c0.5 | 5 | 0.8019 | 3.8 | 28800 | 165600 | 0.1763 | -2.724 | 1 |
| 40 | 3 | 3 | z_only_box | 5 | 0.7103 | 3.698 | 28800 | 165600 | 0.1817 | -2.281 | 1 |
| 40 | 3 | 3 | z_zero | 5 | 0.8489 | 1.329 | 28800 | 165600 | 0.167 | -3.198 | 1 |
| 40 | 2 | 32 | batch_box | 5 | 0.8063 | 0.8919 | 38400 | 2.496e+06 | 2.836 | -3.373 | 1 |
| 40 | 2 | 32 | box_factor_a0.9 | 5 | 0.7574 | 1.047 | 38400 | 2.496e+06 | 2.728 | -2.855 | 1 |
| 40 | 2 | 32 | box_factor_a2 | 5 | 0.7367 | 1.212 | 38400 | 2.496e+06 | 2.692 | -2.543 | 1 |
| 40 | 2 | 32 | raw | 5 | 0.7668 | 1.116 | 38400 | 2.496e+06 | 2.822 | -2.761 | 1 |
| 40 | 2 | 32 | sample_ball | 5 | 0.7731 | 1.05 | 38400 | 2.496e+06 | 2.847 | -2.894 | 1 |
| 40 | 2 | 32 | sample_box | 5 | 0.7534 | 1.071 | 38400 | 2.496e+06 | 2.854 | -2.803 | 1 |
| 40 | 2 | 32 | z_only_box | 5 | 0.7534 | 1.071 | 38400 | 2.496e+06 | 2.881 | -2.803 | 1 |
| 40 | 2 | 32 | z_zero | 5 | 0.8477 | 0.8909 | 38400 | 2.496e+06 | 2.804 | -3.692 | 1 |
| 40 | 4 | 2 | batch_box | 5 | 0.799 | 2.099 | 67200 | 247200 | 0.3631 | -2.477 | 1 |
| 40 | 4 | 2 | box_factor_a0.75 | 5 | 0.7309 | 4.31 | 67200 | 247200 | 0.3202 | -2.047 | 1 |
| 40 | 4 | 2 | raw | 5 | 248.4 | 3061 | 67200 | 247200 | 0.301 | 0.01226 | 1 |
| 40 | 4 | 2 | sample_ball | 5 | 0.8145 | 2.752 | 67200 | 247200 | 0.3103 | -2.478 | 1 |
| 40 | 4 | 2 | sample_box | 5 | 0.7594 | 6.343 | 67200 | 247200 | 0.3121 | -1.834 | 1 |
| 40 | 4 | 2 | shrink_c0.1 | 5 | 0.8403 | 1.637 | 67200 | 247200 | 0.3036 | -2.68 | 1 |
| 40 | 4 | 2 | z_only_box | 5 | 0.8583 | 8.735 | 67200 | 247200 | 0.3171 | -1.826 | 1 |
| 40 | 4 | 2 | z_zero | 5 | 0.8498 | 1.566 | 67200 | 247200 | 0.3049 | -2.771 | 1 |
| 40 | 3 | 6 | batch_box | 5 | 0.7939 | 1.029 | 100800 | 1.1736e+06 | 1.29 | -2.864 | 1 |
| 40 | 3 | 6 | box_factor_a0.9 | 5 | 0.6866 | 1.851 | 100800 | 1.1736e+06 | 1.299 | -2.291 | 1 |
| 40 | 3 | 6 | box_factor_a1.5 | 5 | 0.6504 | 3.414 | 100800 | 1.1736e+06 | 1.31 | -1.912 | 1 |
| 40 | 3 | 6 | raw | 5 | 0.7614 | 5.554 | 100800 | 1.1736e+06 | 1.233 | -2.136 | 1 |
| 40 | 3 | 6 | sample_ball | 5 | 0.7544 | 1.453 | 100800 | 1.1736e+06 | 1.239 | -2.652 | 1 |
| 40 | 3 | 6 | sample_box | 5 | 0.6743 | 2.082 | 100800 | 1.1736e+06 | 1.244 | -2.219 | 1 |
| 40 | 3 | 6 | shrink_c0.875 | 5 | 0.7343 | 4.01 | 100800 | 1.1736e+06 | 1.252 | -2.27 | 1 |
| 40 | 3 | 6 | z_only_box | 5 | 0.6746 | 2.128 | 100800 | 1.1736e+06 | 1.253 | -2.218 | 1 |
| 40 | 3 | 6 | z_zero | 5 | 0.8479 | 0.9437 | 100800 | 1.1736e+06 | 1.234 | -3.129 | 1 |
| 40 | 4 | 3 | batch_box | 5 | 0.792 | 1.33 | 190800 | 1.0728e+06 | 0.9046 | -2.449 | 1 |
| 40 | 4 | 3 | box_factor_a0.9 | 5 | 0.6613 | 3.671 | 190800 | 1.0728e+06 | 0.8631 | -1.875 | 1 |
| 40 | 4 | 3 | raw | 5 | 116.7 | 1238 | 190800 | 1.0728e+06 | 0.8124 | -1.371 | 1 |
| 40 | 4 | 3 | sample_ball | 5 | 0.7661 | 2.055 | 190800 | 1.0728e+06 | 0.8403 | -2.352 | 1 |
| 40 | 4 | 3 | sample_box | 5 | 0.657 | 4.34 | 190800 | 1.0728e+06 | 0.8565 | -1.796 | 1 |
| 40 | 4 | 3 | shrink_c0.25 | 5 | 0.817 | 1.417 | 190800 | 1.0728e+06 | 0.8518 | -2.481 | 1 |
| 40 | 4 | 3 | z_only_box | 5 | 0.6966 | 4.987 | 190800 | 1.0728e+06 | 0.864 | -1.793 | 1 |
| 40 | 4 | 3 | z_zero | 5 | 0.8486 | 1.051 | 190800 | 1.0728e+06 | 0.9005 | -2.719 | 1 |
| 40 | 3 | 10 | batch_box | 10 | 0.7879 | 0.926 | 264000 | 5.172e+06 | 6.491 | -2.842 | 0.9999 |
| 40 | 3 | 10 | box_factor_a0.9 | 10 | 0.682 | 1.405 | 264000 | 5.172e+06 | 6.247 | -2.265 | 0.9999 |
| 40 | 3 | 10 | box_factor_a2 | 10 | 0.6489 | 3.381 | 264000 | 5.172e+06 | 6.241 | -1.737 | 0.9999 |
| 40 | 3 | 10 | raw | 10 | 0.6962 | 3.188 | 264000 | 5.172e+06 | 6.27 | -2.118 | 0.9999 |
| 40 | 3 | 10 | sample_ball | 10 | 0.7308 | 1.299 | 264000 | 5.172e+06 | 6.339 | -2.52 | 0.9999 |
| 40 | 3 | 10 | sample_box | 10 | 0.6692 | 1.554 | 264000 | 5.172e+06 | 6.381 | -2.198 | 0.9999 |
| 40 | 3 | 10 | z_only_box | 10 | 0.6691 | 1.587 | 264000 | 5.172e+06 | 6.438 | -2.198 | 0.9999 |
| 40 | 3 | 10 | z_zero | 10 | 0.8478 | 0.8913 | 264000 | 5.172e+06 | 6.345 | -3.107 | 0.9999 |
| 40 | 5 | 2 | batch_box | 10 | 0.7949 | 1.894 | 302400 | 1.1016e+06 | 1.187 | -2.1 | 1 |
| 40 | 5 | 2 | box_factor_a0.75 | 10 | 0.7192 | 4.911 | 302400 | 1.1016e+06 | 1.101 | -1.671 | 1 |
| 40 | 5 | 2 | raw | 10 | 2.807e+07 | 2.452e+08 | 302400 | 1.1016e+06 | 1.024 | 7592 | 1 |
| 40 | 5 | 2 | sample_ball | 10 | 0.7933 | 3.029 | 302400 | 1.1016e+06 | 1.075 | -2.082 | 1 |
| 40 | 5 | 2 | sample_box | 10 | 0.7911 | 7.875 | 302400 | 1.1016e+06 | 1.045 | -1.461 | 1 |
| 40 | 5 | 2 | shrink_c0.1 | 10 | 0.838 | 1.333 | 302400 | 1.1016e+06 | 1.054 | -2.295 | 1 |
| 40 | 5 | 2 | z_only_box | 10 | 1.08 | 12.68 | 302400 | 1.1016e+06 | 1.096 | -1.453 | 1 |
| 40 | 5 | 2 | z_zero | 10 | 0.8489 | 1.27 | 302400 | 1.1016e+06 | 1.052 | -2.39 | 1 |
| 40 | 3 | 16 | raw | 5 | 0.6772 | 2.116 | 652800 | 2.06016e+07 | 25.88 | -2.102 | 0.9992 |
| 40 | 3 | 16 | sample_box | 5 | 0.6705 | 1.273 | 652800 | 2.06016e+07 | 26.08 | -2.18 | 0.9992 |
| 40 | 3 | 16 | z_zero | 5 | 0.8476 | 0.8799 | 652800 | 2.06016e+07 | 26.26 | -3.086 | 0.9985 |
| 40 | 4 | 6 | raw | 5 | 4.225 | 56.45 | 1.2672e+06 | 1.4652e+07 | 21.33 | -1.6 | 0.9999 |
| 40 | 4 | 6 | sample_box | 5 | 0.6082 | 2.361 | 1.2672e+06 | 1.4652e+07 | 21.45 | -1.735 | 0.9999 |
| 40 | 4 | 6 | z_zero | 5 | 0.848 | 0.888 | 1.2672e+06 | 1.4652e+07 | 21.28 | -2.646 | 0.9998 |
| 60 | 2 | 2 | f_zero | 10 | 0.9059 | 3.204 | 0 | 4800 | 0.01152 | NA | 0 |
| 60 | 2 | 8 | f_zero | 5 | 0.8978 | 1.195 | 0 | 76800 | 0.06225 | NA | 0 |
| 60 | 2 | 24 | f_zero | 5 | 0.8982 | 0.9515 | 0 | 691200 | 1.116 | NA | 0 |
| 60 | 2 | 32 | f_zero | 5 | 0.8982 | 0.9365 | 0 | 1.2288e+06 | 1.932 | NA | 0 |
| 60 | 3 | 3 | f_zero | 5 | 0.8995 | 1.492 | 0 | 32400 | 0.02815 | NA | 0 |
| 60 | 3 | 6 | f_zero | 5 | 0.8983 | 1.008 | 0 | 259200 | 0.3548 | NA | 0 |
| 60 | 3 | 10 | f_zero | 10 | 0.8982 | 0.937 | 0 | 1.2e+06 | 2.003 | NA | 0 |
| 60 | 3 | 16 | f_zero | 5 | 0.8982 | 0.9219 | 0 | 4.9152e+06 | 8.684 | NA | 0 |
| 60 | 4 | 2 | f_zero | 5 | 0.8982 | 1.781 | 0 | 19200 | 0.02062 | NA | 0 |
| 60 | 4 | 3 | f_zero | 5 | 0.8984 | 1.142 | 0 | 97200 | 0.08211 | NA | 0 |
| 60 | 4 | 6 | f_zero | 5 | 0.8981 | 0.932 | 0 | 1.5552e+06 | 2.757 | NA | 0 |
| 60 | 5 | 2 | f_zero | 10 | 0.8994 | 1.418 | 0 | 38400 | 0.03136 | NA | 0 |
| 60 | 2 | 2 | batch_box | 10 | 0.8846 | 3.417 | 2400 | 12000 | 0.04085 | -3.99 | 1 |
| 60 | 2 | 2 | box_factor_a0.75 | 10 | 0.8674 | 4.115 | 2400 | 12000 | 0.03457 | -3.619 | 1 |
| 60 | 2 | 2 | raw | 10 | 0.9376 | 5.386 | 2400 | 12000 | 0.03319 | -3.562 | 1 |
| 60 | 2 | 2 | sample_ball | 10 | 0.8986 | 3.388 | 2400 | 12000 | 0.03599 | -4.089 | 1 |
| 60 | 2 | 2 | sample_box | 10 | 0.8649 | 4.592 | 2400 | 12000 | 0.03704 | -3.447 | 1 |
| 60 | 2 | 2 | shrink_c0.5 | 10 | 0.9032 | 3.869 | 2400 | 12000 | 0.03565 | -3.913 | 1 |
| 60 | 2 | 2 | z_only_box | 10 | 0.8649 | 4.592 | 2400 | 12000 | 0.03508 | -3.447 | 1 |
| 60 | 2 | 2 | z_zero | 10 | 0.9059 | 3.204 | 2400 | 12000 | 0.03412 | -4.265 | 1 |
| 60 | 2 | 8 | batch_box | 5 | 0.8739 | 1.276 | 9600 | 163200 | 0.1717 | -4.077 | 1 |
| 60 | 2 | 8 | box_factor_a0.9 | 5 | 0.8373 | 1.856 | 9600 | 163200 | 0.1494 | -3.527 | 1 |
| 60 | 2 | 8 | box_factor_a2 | 5 | 0.8173 | 2.757 | 9600 | 163200 | 0.1618 | -2.969 | 1 |
| 60 | 2 | 8 | raw | 5 | 0.8542 | 2.081 | 9600 | 163200 | 0.1324 | -3.481 | 1 |
| 60 | 2 | 8 | sample_ball | 5 | 0.8715 | 1.401 | 9600 | 163200 | 0.1546 | -3.932 | 1 |
| 60 | 2 | 8 | sample_box | 5 | 0.8339 | 1.953 | 9600 | 163200 | 0.1542 | -3.461 | 1 |
| 60 | 2 | 8 | z_only_box | 5 | 0.8339 | 1.953 | 9600 | 163200 | 0.1536 | -3.461 | 1 |
| 60 | 2 | 8 | z_zero | 5 | 0.8978 | 1.195 | 9600 | 163200 | 0.1714 | -4.335 | 1 |
| 60 | 2 | 24 | batch_box | 5 | 0.8733 | 0.9758 | 28800 | 1.4112e+06 | 2.376 | -4.03 | 1 |
| 60 | 2 | 24 | box_factor_a0.9 | 5 | 0.8374 | 1.235 | 28800 | 1.4112e+06 | 2.415 | -3.489 | 1 |
| 60 | 2 | 24 | box_factor_a2 | 5 | 0.819 | 1.556 | 28800 | 1.4112e+06 | 2.419 | -3.086 | 1 |
| 60 | 2 | 24 | raw | 5 | 0.8465 | 1.315 | 28800 | 1.4112e+06 | 2.358 | -3.428 | 1 |
| 60 | 2 | 24 | sample_ball | 5 | 0.8572 | 1.138 | 28800 | 1.4112e+06 | 2.382 | -3.673 | 1 |
| 60 | 2 | 24 | sample_box | 5 | 0.8342 | 1.275 | 28800 | 1.4112e+06 | 2.373 | -3.431 | 1 |
| 60 | 2 | 24 | z_only_box | 5 | 0.8342 | 1.275 | 28800 | 1.4112e+06 | 2.38 | -3.431 | 1 |
| 60 | 2 | 24 | z_zero | 5 | 0.8982 | 0.9515 | 28800 | 1.4112e+06 | 2.397 | -4.305 | 1 |
| 60 | 3 | 3 | batch_box | 5 | 0.8675 | 1.737 | 28800 | 165600 | 0.229 | -3.417 | 1 |
| 60 | 3 | 3 | box_factor_a0.9 | 5 | 0.8088 | 3.666 | 28800 | 165600 | 0.2118 | -2.834 | 1 |
| 60 | 3 | 3 | raw | 5 | 1.939 | 31.1 | 28800 | 165600 | 0.2136 | -2.683 | 1 |
| 60 | 3 | 3 | sample_ball | 5 | 0.8661 | 2.019 | 28800 | 165600 | 0.2084 | -3.4 | 1 |
| 60 | 3 | 3 | sample_box | 5 | 0.8061 | 4.13 | 28800 | 165600 | 0.2088 | -2.754 | 1 |
| 60 | 3 | 3 | shrink_c0.5 | 5 | 0.8849 | 4.899 | 28800 | 165600 | 0.1977 | -3.231 | 1 |
| 60 | 3 | 3 | z_only_box | 5 | 0.8157 | 4.553 | 28800 | 165600 | 0.2028 | -2.75 | 1 |
| 60 | 3 | 3 | z_zero | 5 | 0.8995 | 1.492 | 28800 | 165600 | 0.1903 | -3.673 | 1 |
| 60 | 2 | 32 | batch_box | 5 | 0.8723 | 0.9514 | 38400 | 2.496e+06 | 4.324 | -4.036 | 1 |
| 60 | 2 | 32 | box_factor_a0.9 | 5 | 0.8369 | 1.153 | 38400 | 2.496e+06 | 4.135 | -3.501 | 1 |
| 60 | 2 | 32 | box_factor_a2 | 5 | 0.8205 | 1.383 | 38400 | 2.496e+06 | 4.12 | -3.142 | 1 |
| 60 | 2 | 32 | raw | 5 | 0.8452 | 1.216 | 38400 | 2.496e+06 | 4.357 | -3.434 | 1 |
| 60 | 2 | 32 | sample_ball | 5 | 0.8529 | 1.106 | 38400 | 2.496e+06 | 4.347 | -3.625 | 1 |
| 60 | 2 | 32 | sample_box | 5 | 0.8339 | 1.184 | 38400 | 2.496e+06 | 4.381 | -3.447 | 1 |
| 60 | 2 | 32 | z_only_box | 5 | 0.8339 | 1.184 | 38400 | 2.496e+06 | 4.389 | -3.447 | 1 |
| 60 | 2 | 32 | z_zero | 5 | 0.8982 | 0.9365 | 38400 | 2.496e+06 | 4.259 | -4.32 | 1 |
| 60 | 4 | 2 | batch_box | 5 | 0.8655 | 2.319 | 67200 | 247200 | 0.4488 | -2.849 | 1 |
| 60 | 4 | 2 | box_factor_a0.75 | 5 | 0.8226 | 5.209 | 67200 | 247200 | 0.395 | -2.399 | 1 |
| 60 | 4 | 2 | raw | 5 | 381.5 | 9260 | 67200 | 247200 | 0.3657 | 2.413 | 1 |
| 60 | 4 | 2 | sample_ball | 5 | 0.8711 | 2.742 | 67200 | 247200 | 0.3892 | -2.886 | 1 |
| 60 | 4 | 2 | sample_box | 5 | 0.8532 | 7.667 | 67200 | 247200 | 0.3723 | -2.182 | 1 |
| 60 | 4 | 2 | shrink_c0.25 | 5 | 0.8807 | 3.164 | 67200 | 247200 | 0.3822 | -2.9 | 1 |
| 60 | 4 | 2 | z_only_box | 5 | 0.9668 | 10.26 | 67200 | 247200 | 0.3845 | -2.172 | 1 |
| 60 | 4 | 2 | z_zero | 5 | 0.8982 | 1.781 | 67200 | 247200 | 0.367 | -3.118 | 1 |
| 60 | 3 | 6 | batch_box | 5 | 0.8631 | 1.109 | 100800 | 1.1736e+06 | 1.777 | -3.33 | 1 |
| 60 | 3 | 6 | box_factor_a0.9 | 5 | 0.7876 | 2.163 | 100800 | 1.1736e+06 | 1.822 | -2.737 | 1 |
| 60 | 3 | 6 | box_factor_a1.5 | 5 | 0.7732 | 4.112 | 100800 | 1.1736e+06 | 1.814 | -2.339 | 1 |
| 60 | 3 | 6 | raw | 5 | 0.889 | 6.469 | 100800 | 1.1736e+06 | 1.715 | -2.626 | 1 |
| 60 | 3 | 6 | sample_ball | 5 | 0.8422 | 1.501 | 100800 | 1.1736e+06 | 1.727 | -3.18 | 1 |
| 60 | 3 | 6 | sample_box | 5 | 0.78 | 2.449 | 100800 | 1.1736e+06 | 1.707 | -2.663 | 1 |
| 60 | 3 | 6 | shrink_c0.875 | 5 | 0.8443 | 4.698 | 100800 | 1.1736e+06 | 1.788 | -2.755 | 1 |
| 60 | 3 | 6 | z_only_box | 5 | 0.7816 | 2.534 | 100800 | 1.1736e+06 | 1.766 | -2.662 | 1 |
| 60 | 3 | 6 | z_zero | 5 | 0.8983 | 1.008 | 100800 | 1.1736e+06 | 1.748 | -3.571 | 1 |
| 60 | 4 | 3 | batch_box | 5 | 0.861 | 1.43 | 190800 | 1.0728e+06 | 1.505 | -2.766 | 1 |
| 60 | 4 | 3 | box_factor_a0.9 | 5 | 0.7749 | 4.541 | 190800 | 1.0728e+06 | 1.271 | -2.169 | 1 |
| 60 | 4 | 3 | raw | 5 | 37.18 | 489.5 | 190800 | 1.0728e+06 | 1.329 | -1.772 | 1 |
| 60 | 4 | 3 | sample_ball | 5 | 0.8518 | 2.015 | 190800 | 1.0728e+06 | 1.432 | -2.717 | 1 |
| 60 | 4 | 3 | sample_box | 5 | 0.7742 | 5.348 | 190800 | 1.0728e+06 | 1.489 | -2.088 | 1 |
| 60 | 4 | 3 | shrink_c0.375 | 5 | 0.8674 | 3.34 | 190800 | 1.0728e+06 | 1.29 | -2.67 | 1 |
| 60 | 4 | 3 | z_only_box | 5 | 0.8174 | 6.699 | 190800 | 1.0728e+06 | 1.454 | -2.082 | 1 |
| 60 | 4 | 3 | z_zero | 5 | 0.8984 | 1.142 | 190800 | 1.0728e+06 | 1.426 | -3.014 | 1 |
| 60 | 3 | 10 | batch_box | 10 | 0.8614 | 0.9899 | 264000 | 5.172e+06 | 9.811 | -3.291 | 1 |
| 60 | 3 | 10 | box_factor_a0.9 | 10 | 0.7877 | 1.661 | 264000 | 5.172e+06 | 9.589 | -2.699 | 1 |
| 60 | 3 | 10 | box_factor_a1.5 | 10 | 0.7638 | 3.015 | 264000 | 5.172e+06 | 9.559 | -2.342 | 1 |
| 60 | 3 | 10 | raw | 10 | 0.8081 | 4.023 | 264000 | 5.172e+06 | 9.615 | -2.594 | 1 |
| 60 | 3 | 10 | sample_ball | 10 | 0.8295 | 1.368 | 264000 | 5.172e+06 | 9.461 | -3.044 | 1 |
| 60 | 3 | 10 | sample_box | 10 | 0.7797 | 1.858 | 264000 | 5.172e+06 | 9.588 | -2.631 | 1 |
| 60 | 3 | 10 | z_only_box | 10 | 0.781 | 1.927 | 264000 | 5.172e+06 | 9.566 | -2.63 | 1 |
| 60 | 3 | 10 | z_zero | 10 | 0.8982 | 0.937 | 264000 | 5.172e+06 | 9.577 | -3.53 | 1 |
| 60 | 5 | 2 | batch_box | 10 | 0.8638 | 2.064 | 302400 | 1.1016e+06 | 1.575 | -2.33 | 1 |
| 60 | 5 | 2 | box_factor_a0.75 | 10 | 0.82 | 6.019 | 302400 | 1.1016e+06 | 1.432 | -1.877 | 1 |
| 60 | 5 | 2 | raw | 10 | 5.798e+06 | 1.297e+08 | 302400 | 1.1016e+06 | 1.357 | 8411 | 1 |
| 60 | 5 | 2 | sample_ball | 10 | 0.8701 | 2.885 | 302400 | 1.1016e+06 | 1.444 | -2.348 | 1 |
| 60 | 5 | 2 | sample_box | 10 | 0.9001 | 9.65 | 302400 | 1.1016e+06 | 1.403 | -1.662 | 1 |
| 60 | 5 | 2 | shrink_c0.25 | 10 | 1.085 | 9.613 | 302400 | 1.1016e+06 | 1.397 | -2.366 | 1 |
| 60 | 5 | 2 | z_only_box | 10 | 1.314 | 19.72 | 302400 | 1.1016e+06 | 1.415 | -1.651 | 1 |
| 60 | 5 | 2 | z_zero | 10 | 0.8994 | 1.418 | 302400 | 1.1016e+06 | 1.367 | -2.597 | 1 |
| 60 | 3 | 16 | raw | 5 | 0.7927 | 2.64 | 652800 | 2.06016e+07 | 39.28 | -2.573 | 0.9999 |
| 60 | 3 | 16 | sample_box | 5 | 0.7807 | 1.498 | 652800 | 2.06016e+07 | 39.47 | -2.611 | 0.9999 |
| 60 | 3 | 16 | z_zero | 5 | 0.8982 | 0.9219 | 652800 | 2.06016e+07 | 39.55 | -3.505 | 0.9998 |
| 60 | 4 | 6 | raw | 5 | 6.071 | 159.5 | 1.2672e+06 | 1.4652e+07 | 32.28 | -1.942 | 1 |
| 60 | 4 | 6 | sample_box | 5 | 0.7375 | 2.87 | 1.2672e+06 | 1.4652e+07 | 32.84 | -2.035 | 1 |
| 60 | 4 | 6 | z_zero | 5 | 0.8981 | 0.932 | 1.2672e+06 | 1.4652e+07 | 33.05 | -2.94 | 1 |
| 80 | 2 | 2 | f_zero | 10 | 0.9244 | 3.608 | 0 | 4800 | 0.01291 | NA | 0 |
| 80 | 2 | 8 | f_zero | 5 | 0.9208 | 1.283 | 0 | 76800 | 0.08584 | NA | 0 |
| 80 | 2 | 24 | f_zero | 5 | 0.9205 | 0.9842 | 0 | 691200 | 1.551 | NA | 0 |
| 80 | 2 | 32 | f_zero | 5 | 0.9204 | 0.9647 | 0 | 1.2288e+06 | 2.568 | NA | 0 |
| 80 | 3 | 3 | f_zero | 5 | 0.9197 | 1.632 | 0 | 32400 | 0.03876 | NA | 0 |
| 80 | 3 | 6 | f_zero | 5 | 0.9206 | 1.053 | 0 | 259200 | 0.5612 | NA | 0 |
| 80 | 3 | 10 | f_zero | 10 | 0.9205 | 0.9654 | 0 | 1.2e+06 | 2.73 | NA | 0 |
| 80 | 3 | 16 | f_zero | 5 | 0.9204 | 0.9462 | 0 | 4.9152e+06 | 11.64 | NA | 0 |
| 80 | 4 | 2 | f_zero | 5 | 0.9199 | 1.977 | 0 | 19200 | 0.02508 | NA | 0 |
| 80 | 4 | 3 | f_zero | 5 | 0.9206 | 1.218 | 0 | 97200 | 0.1317 | NA | 0 |
| 80 | 4 | 6 | f_zero | 5 | 0.9205 | 0.9598 | 0 | 1.5552e+06 | 3.682 | NA | 0 |
| 80 | 5 | 2 | f_zero | 10 | 0.9205 | 1.548 | 0 | 38400 | 0.03996 | NA | 0 |
| 80 | 2 | 2 | batch_box | 10 | 0.9084 | 3.818 | 2400 | 12000 | 0.04584 | -4.719 | 1 |
| 80 | 2 | 2 | box_factor_a0.9 | 10 | 0.8941 | 4.978 | 2400 | 12000 | 0.03844 | -4.194 | 1 |
| 80 | 2 | 2 | raw | 10 | 0.9398 | 5.943 | 2400 | 12000 | 0.03625 | -4.255 | 1 |
| 80 | 2 | 2 | sample_ball | 10 | 0.9196 | 3.762 | 2400 | 12000 | 0.0392 | -4.842 | 1 |
| 80 | 2 | 2 | sample_box | 10 | 0.894 | 5.203 | 2400 | 12000 | 0.03858 | -4.121 | 1 |
| 80 | 2 | 2 | shrink_c0.375 | 10 | 0.9188 | 4.017 | 2400 | 12000 | 0.03613 | -4.72 | 1 |
| 80 | 2 | 2 | z_only_box | 10 | 0.894 | 5.203 | 2400 | 12000 | 0.03768 | -4.121 | 1 |
| 80 | 2 | 2 | z_zero | 10 | 0.9244 | 3.608 | 2400 | 12000 | 0.03718 | -4.999 | 1 |
| 80 | 2 | 8 | batch_box | 5 | 0.9035 | 1.378 | 9600 | 163200 | 0.2558 | -4.765 | 1 |
| 80 | 2 | 8 | box_factor_a0.9 | 5 | 0.8751 | 2.117 | 9600 | 163200 | 0.1935 | -4.151 | 1 |
| 80 | 2 | 8 | box_factor_a2 | 5 | 0.8615 | 3.253 | 9600 | 163200 | 0.2226 | -3.501 | 1 |
| 80 | 2 | 8 | raw | 5 | 0.89 | 2.453 | 9600 | 163200 | 0.2 | -4.108 | 1 |
| 80 | 2 | 8 | sample_ball | 5 | 0.9025 | 1.505 | 9600 | 163200 | 0.2324 | -4.647 | 1 |
| 80 | 2 | 8 | sample_box | 5 | 0.8725 | 2.237 | 9600 | 163200 | 0.2 | -4.076 | 1 |
| 80 | 2 | 8 | z_only_box | 5 | 0.8725 | 2.237 | 9600 | 163200 | 0.2232 | -4.076 | 1 |
| 80 | 2 | 8 | z_zero | 5 | 0.9208 | 1.283 | 9600 | 163200 | 0.2463 | -5.027 | 1 |
| 80 | 2 | 24 | batch_box | 5 | 0.9023 | 1.019 | 28800 | 1.4112e+06 | 3.283 | -4.722 | 1 |
| 80 | 2 | 24 | box_factor_a0.9 | 5 | 0.8744 | 1.36 | 28800 | 1.4112e+06 | 3.316 | -4.119 | 1 |
| 80 | 2 | 24 | box_factor_a2 | 5 | 0.8599 | 1.786 | 28800 | 1.4112e+06 | 3.244 | -3.633 | 1 |
| 80 | 2 | 24 | raw | 5 | 0.8832 | 1.437 | 28800 | 1.4112e+06 | 3.249 | -4.062 | 1 |
| 80 | 2 | 24 | sample_ball | 5 | 0.8927 | 1.193 | 28800 | 1.4112e+06 | 3.29 | -4.386 | 1 |
| 80 | 2 | 24 | sample_box | 5 | 0.8719 | 1.412 | 28800 | 1.4112e+06 | 3.258 | -4.054 | 1 |
| 80 | 2 | 24 | z_only_box | 5 | 0.8719 | 1.412 | 28800 | 1.4112e+06 | 3.269 | -4.054 | 1 |
| 80 | 2 | 24 | z_zero | 5 | 0.9205 | 0.9842 | 28800 | 1.4112e+06 | 3.277 | -4.997 | 1 |
| 80 | 3 | 3 | batch_box | 5 | 0.8957 | 1.897 | 28800 | 165600 | 0.3064 | -3.988 | 1 |
| 80 | 3 | 3 | box_factor_a0.9 | 5 | 0.8549 | 4.218 | 28800 | 165600 | 0.2739 | -3.341 | 1 |
| 80 | 3 | 3 | box_factor_a1.25 | 5 | 0.8701 | 6.295 | 28800 | 165600 | 0.2706 | -3.041 | 1 |
| 80 | 3 | 3 | raw | 5 | 2.738 | 37.01 | 28800 | 165600 | 0.2655 | -3.182 | 1 |
| 80 | 3 | 3 | sample_ball | 5 | 0.8986 | 2.094 | 28800 | 165600 | 0.2772 | -3.996 | 1 |
| 80 | 3 | 3 | sample_box | 5 | 0.8558 | 4.776 | 28800 | 165600 | 0.2749 | -3.253 | 1 |
| 80 | 3 | 3 | shrink_c0.75 | 5 | 1.448 | 16.48 | 28800 | 165600 | 0.273 | -3.498 | 1 |
| 80 | 3 | 3 | z_only_box | 5 | 0.868 | 5.208 | 28800 | 165600 | 0.2839 | -3.247 | 1 |
| 80 | 3 | 3 | z_zero | 5 | 0.9197 | 1.632 | 28800 | 165600 | 0.2655 | -4.251 | 1 |
| 80 | 2 | 32 | batch_box | 5 | 0.9016 | 0.9896 | 38400 | 2.496e+06 | 5.686 | -4.744 | 1 |
| 80 | 2 | 32 | box_factor_a0.9 | 5 | 0.8734 | 1.254 | 38400 | 2.496e+06 | 5.437 | -4.149 | 1 |
| 80 | 2 | 32 | box_factor_a2 | 5 | 0.8598 | 1.565 | 38400 | 2.496e+06 | 5.456 | -3.716 | 1 |
| 80 | 2 | 32 | raw | 5 | 0.8811 | 1.323 | 38400 | 2.496e+06 | 5.735 | -4.088 | 1 |
| 80 | 2 | 32 | sample_ball | 5 | 0.889 | 1.156 | 38400 | 2.496e+06 | 5.735 | -4.348 | 1 |
| 80 | 2 | 32 | sample_box | 5 | 0.8711 | 1.293 | 38400 | 2.496e+06 | 5.729 | -4.088 | 1 |
| 80 | 2 | 32 | z_only_box | 5 | 0.8711 | 1.293 | 38400 | 2.496e+06 | 5.768 | -4.088 | 1 |
| 80 | 2 | 32 | z_zero | 5 | 0.9204 | 0.9647 | 38400 | 2.496e+06 | 5.575 | -5.026 | 1 |
| 80 | 4 | 2 | batch_box | 5 | 0.8997 | 2.514 | 67200 | 247200 | 0.5722 | -3.229 | 1 |
| 80 | 4 | 2 | box_factor_a0.75 | 5 | 0.8761 | 5.722 | 67200 | 247200 | 0.5143 | -2.729 | 1 |
| 80 | 4 | 2 | raw | 5 | 4.78e+04 | 1.235e+06 | 67200 | 247200 | 0.4878 | 286.8 | 1 |
| 80 | 4 | 2 | sample_ball | 5 | 0.9079 | 2.787 | 67200 | 247200 | 0.5142 | -3.292 | 1 |
| 80 | 4 | 2 | sample_box | 5 | 0.91 | 8.574 | 67200 | 247200 | 0.5118 | -2.492 | 1 |
| 80 | 4 | 2 | shrink_c0.1 | 5 | 0.9162 | 2.089 | 67200 | 247200 | 0.4969 | -3.414 | 1 |
| 80 | 4 | 2 | z_only_box | 5 | 1.007 | 12.06 | 67200 | 247200 | 0.5207 | -2.474 | 1 |
| 80 | 4 | 2 | z_zero | 5 | 0.9199 | 1.977 | 67200 | 247200 | 0.4778 | -3.503 | 1 |
| 80 | 3 | 6 | batch_box | 5 | 0.8953 | 1.181 | 100800 | 1.1736e+06 | 2.864 | -3.862 | 1 |
| 80 | 3 | 6 | box_factor_a0.9 | 5 | 0.8393 | 2.554 | 100800 | 1.1736e+06 | 2.883 | -3.212 | 1 |
| 80 | 3 | 6 | box_factor_a1.25 | 5 | 0.8303 | 3.874 | 100800 | 1.1736e+06 | 2.866 | -2.94 | 1 |
| 80 | 3 | 6 | raw | 5 | 0.8966 | 7.65 | 100800 | 1.1736e+06 | 2.805 | -3.114 | 1 |
| 80 | 3 | 6 | sample_ball | 5 | 0.8836 | 1.551 | 100800 | 1.1736e+06 | 2.794 | -3.743 | 1 |
| 80 | 3 | 6 | sample_box | 5 | 0.8347 | 2.906 | 100800 | 1.1736e+06 | 2.824 | -3.13 | 1 |
| 80 | 3 | 6 | shrink_c0.75 | 5 | 0.8691 | 3.897 | 100800 | 1.1736e+06 | 2.8 | -3.381 | 1 |
| 80 | 3 | 6 | z_only_box | 5 | 0.8407 | 3.183 | 100800 | 1.1736e+06 | 2.864 | -3.128 | 1 |
| 80 | 3 | 6 | z_zero | 5 | 0.9206 | 1.053 | 100800 | 1.1736e+06 | 2.769 | -4.107 | 1 |
| 80 | 4 | 3 | batch_box | 5 | 0.8932 | 1.51 | 190800 | 1.0728e+06 | 2.135 | -3.171 | 1 |
| 80 | 4 | 3 | box_factor_a0.875 | 5 | 0.8236 | 5.008 | 190800 | 1.0728e+06 | 1.986 | -2.538 | 1 |
| 80 | 4 | 3 | raw | 5 | 54.17 | 2860 | 190800 | 1.0728e+06 | 2.004 | -1.994 | 1 |
| 80 | 4 | 3 | sample_ball | 5 | 0.8871 | 2.071 | 190800 | 1.0728e+06 | 2.041 | -3.153 | 1 |
| 80 | 4 | 3 | sample_box | 5 | 0.8241 | 6.199 | 190800 | 1.0728e+06 | 2.026 | -2.427 | 1 |
| 80 | 4 | 3 | shrink_c0.375 | 5 | 0.8955 | 5.525 | 190800 | 1.0728e+06 | 1.89 | -3.066 | 1 |
| 80 | 4 | 3 | z_only_box | 5 | 0.8773 | 8.46 | 190800 | 1.0728e+06 | 1.993 | -2.418 | 1 |
| 80 | 4 | 3 | z_zero | 5 | 0.9206 | 1.218 | 190800 | 1.0728e+06 | 1.973 | -3.424 | 1 |
| 80 | 3 | 10 | batch_box | 10 | 0.8942 | 1.033 | 264000 | 5.172e+06 | 13.25 | -3.793 | 1 |
| 80 | 3 | 10 | box_factor_a0.9 | 10 | 0.839 | 1.873 | 264000 | 5.172e+06 | 13.28 | -3.145 | 1 |
| 80 | 3 | 10 | box_factor_a1.25 | 10 | 0.8264 | 2.752 | 264000 | 5.172e+06 | 13.2 | -2.895 | 1 |
| 80 | 3 | 10 | raw | 10 | 0.8593 | 4.528 | 264000 | 5.172e+06 | 13.18 | -3.048 | 1 |
| 80 | 3 | 10 | sample_ball | 10 | 0.874 | 1.413 | 264000 | 5.172e+06 | 13.13 | -3.578 | 1 |
| 80 | 3 | 10 | sample_box | 10 | 0.8338 | 2.106 | 264000 | 5.172e+06 | 13.19 | -3.069 | 1 |
| 80 | 3 | 10 | shrink_c0.875 | 10 | 0.8562 | 3.336 | 264000 | 5.172e+06 | 13.02 | -3.177 | 1 |
| 80 | 3 | 10 | z_only_box | 10 | 0.8368 | 2.241 | 264000 | 5.172e+06 | 13.36 | -3.068 | 1 |
| 80 | 3 | 10 | z_zero | 10 | 0.9205 | 0.9654 | 264000 | 5.172e+06 | 13.06 | -4.035 | 1 |
| 80 | 5 | 2 | batch_box | 10 | 0.8963 | 2.282 | 302400 | 1.1016e+06 | 1.901 | -2.627 | 1 |
| 80 | 5 | 2 | box_factor_a0.5 | 10 | 0.8674 | 3.757 | 302400 | 1.1016e+06 | 1.745 | -2.367 | 1 |
| 80 | 5 | 2 | raw | 10 | 4.662e+06 | 2.294e+08 | 302400 | 1.1016e+06 | 1.676 | 2.047e+05 | 1 |
| 80 | 5 | 2 | sample_ball | 10 | 0.8993 | 2.909 | 302400 | 1.1016e+06 | 1.766 | -2.675 | 1 |
| 80 | 5 | 2 | sample_box | 10 | 0.9503 | 11.4 | 302400 | 1.1016e+06 | 1.722 | -1.883 | 1 |
| 80 | 5 | 2 | shrink_c0.1 | 10 | 0.9155 | 1.668 | 302400 | 1.1016e+06 | 1.704 | -2.809 | 1 |
| 80 | 5 | 2 | z_only_box | 10 | 1.481 | 23.97 | 302400 | 1.1016e+06 | 1.763 | -1.863 | 1 |
| 80 | 5 | 2 | z_zero | 10 | 0.9205 | 1.548 | 302400 | 1.1016e+06 | 1.696 | -2.903 | 1 |
| 80 | 3 | 16 | raw | 5 | 0.8472 | 2.968 | 652800 | 2.06016e+07 | 54.72 | -3.04 | 1 |
| 80 | 3 | 16 | sample_box | 5 | 0.8333 | 1.676 | 652800 | 2.06016e+07 | 55.39 | -3.063 | 1 |
| 80 | 3 | 16 | z_zero | 5 | 0.9204 | 0.9462 | 652800 | 2.06016e+07 | 54.85 | -4.019 | 1 |
| 80 | 4 | 6 | raw | 5 | 9.633 | 207.4 | 1.2672e+06 | 1.4652e+07 | 41.04 | -2.255 | 1 |
| 80 | 4 | 6 | sample_box | 5 | 0.7967 | 3.354 | 1.2672e+06 | 1.4652e+07 | 41.07 | -2.331 | 1 |
| 80 | 4 | 6 | z_zero | 5 | 0.9205 | 0.9598 | 1.2672e+06 | 1.4652e+07 | 39.29 | -3.309 | 1 |

## Generator, projection, and crossover diagnostics

Generator MSE/bias/MAE, nonlinear-correction variance across paired repetitions, box and ball violation rates, overshoot energies, activation rates, Batch-alpha quantiles, per-method win fractions, and paired win fractions against Raw/box/`z=0`/`f=0` are columns in the CSV/JSON summaries. The six PDFs in the result figure directory show error against f calls, wall time, and total stochastic work, plus generator bias, violation rate, and the winner sequence by dimension.

There is no single crossover ordered only by scalar work: how work is allocated between Picard depth and sibling sampling is decisive. The observed regimes are:

1. Low-depth/low-sample (`n=2,M=2`): the certified box already wins or ties the best control; pure `z=0` is not best.
2. Balanced medium/high sampling (`n=3,M=6–16`): box beats suppression and, except at `d=20`, beats or approaches Raw; Raw and box draw together as M rises.
3. Deep/under-sampled (`n=4/5` with small M): Raw can explode and illegal tighter bounds can beat the exact box by suppressing recursion noise.
4. Wrong-PDE baselines: `z=0/f=0` converge to a clear nonzero bias floor and never recover the nonlinear solution.

The integrity audit covers all 2,560 final repetition files. It confirms bitwise-identical root terminal estimates across paired methods, bitwise equality of `z=0` and `f=0` outputs, zero nonlinear correction/f calls for `f=0`, equal recursive f work for all other methods, finite predictions/truth, and complete manifests. Full coordinatewise arrays (about 2.6 GB) remain locally preserved; committed compact float64 archives retain every value prediction, nonlinear value correction, and pointwise gradient-error norm.

## Life-or-death conclusion and paper recommendation

Across 48 comparable dimension/budget cells, the primary untuned certified box beats all available suppression controls in 37; a suppression control wins in 11, concentrated in deep/under-sampled cells. The main-paper recommendation is therefore to use this benchmark as qualified positive evidence that certified geometry can regularize value estimation beyond generic gradient suppression. The claim must be limited to value accuracy and to the observed balanced-sampling regime. The paper should also show the relaxed-valid-box result, disclose that tighter illegal bounds still win in some depth-noise cells, state that gradient error does not improve over `z=0`, and avoid claiming that Batch-IR or the tightest certificate is universally optimal. No manuscript source was modified in this branch.

This satisfies the strong-signal decision rule for value estimation: certified box beats `z=0`, `f=0`, validation-tuned constant shrinkage, and validation-tuned invalid tighter bounds over meaningful medium/high-budget regimes, while the high-sampling extension shows Raw moving toward box and both remaining below the wrong-PDE floor.

VERDICT A:
Certified geometry has an empirical regime of advantage beyond generic shrinkage.
