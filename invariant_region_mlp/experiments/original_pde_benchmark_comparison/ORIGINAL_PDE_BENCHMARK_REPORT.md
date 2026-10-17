# Original-PDE benchmark report

## Scope and experimental contract

This study tests Terminal-Aware Tangent-Actor with direct martingale control
variates (MCV) on PDEs taken from the official
[Deep BSDE](https://github.com/frankhan91/DeepBSDE),
[Deep Picard Iteration](https://github.com/DeepOptimalControl/DeepPicardIteration),
and [SDGD-PINN](https://github.com/zheyuanhu01/SDGD_PINN) repositories.  It keeps the competitors' published
architectures and algorithms intact.  In particular, no baseline receives our
terminal-aware Actor construction, tangent labels, invariant-region layer, or
MCV.  An original terminal construction or variance-reduction device remains
in place when it is part of that baseline.

Actor hyperparameters and iteration count are selected without analytic
nonterminal values or gradients.  Selection uses held-out noisy tangent labels
and a replicated common-state policy-value bank.  Final test references are
opened only after the checkpoint is frozen.  A later screen of all checkpoints
is descriptive and is never used to change the selected method.

All GPU runs use one NVIDIA GH200 120 GB GPU.  Reported training and inference
times are synchronized wall times.  Failed and restarted allocations remain in
the resource audit.

## Benchmarks

| Benchmark | PDE and terminal condition | Native task | Independent reference | Eligibility |
|---|---|---|---|---|
| HJB-Log, `d=100`, `T=1` | `u_t + Delta u - ||grad u||^2 = 0`; `g(x)=log((1+||x||^2)/2)` | `u(0,0)` for Deep BSDE; distributed function task for SDGD | Hopf–Cole plus 1D Gauss–Laguerre quadrature | globally concave in `z=sqrt(2) grad u` |
| OU/GMM HJB, `d=100`, `T=1` | `u_t + (alpha/2) Delta u - F.grad(u) - (alpha/2)||grad u||^2 - d theta = 0`, with `F=theta(mu-x)` and `alpha=theta=1`; terminal value is the negative log density of the official five-component GMM | official `t~U[0,1]`, `x|t~N(0,(4+t)I)` function task | exact evolving Gaussian mixture | known linear drift plus globally concave gradient driver |
| Allen–Cahn, `d=100`, `T=0.3` | `u_t + Delta u + u-u^3 = 0`; `g(x)=1/(2+0.4||x||^2)` | `u(0,0)` | independent squared-radius PDE with mesh refinement | concave on the comparison-principle interval `[0,1]` |

HJB-Rosenbrock was not used.  Its driver is compatible, but the official
repository omits its generated test data, and the supplied 1,000-path reference
generator has no numerical-error audit.  Changing that PDE or treating this
reference as exact would violate the eligibility rule.

## Reference validation

The HJB-Log quadrature gives `u(0,0)=4.590161724604883`; orders 64, 128,
256 and adaptive integration agree within `5.4e-15`.  The independent NumPy
and Torch implementations of the DPI analytic solution agree within
`1.9e-13`, and its maximum automatic-differentiation PDE residual is
`1.2e-12`.  For Allen–Cahn, Richardson extrapolation gives
`u(0,0)=0.05278135060749693` with a `9.64e-6` refinement-difference proxy.
Allen–Cahn results near that scale are therefore labelled reference-limited.

## Baseline reproduction before comparison

Deep BSDE was ported from its official TensorFlow implementation to PyTorch
because TensorFlow is absent on the GH200 image.  The port retains the scalar
free initial value, free initial gradient, 19 nonshared time-slice networks,
normalization, 20-step BSDE rollout, loss, float64 arithmetic, and published
iteration schedules.  It remains a native single-point solver.

| Fixed official method | Evaluation task | Seeds | Error | Mean training time |
|---|---|---:|---:|---:|
| Deep BSDE | HJB-Log `u(0,0)` | 3 | `0.0096785 ± 0.0027611` absolute | 65.15 s |
| Deep BSDE | Allen–Cahn `u(0,0)` | 3 | `0.00016911 ± 0.00017982` absolute | 127.69 s |
| DPI | OU/GMM common 960-state bank | 1 | `17.6638` RMSE | 627.42 s |
| PINN-HTE | OU/GMM common 960-state bank | 1 | `38.6698` RMSE | 767.21 s |
| D-DBSDE | OU/GMM common 960-state bank | 1 | `30.7247` RMSE | 453.05 s |

The `±` entries are sample standard deviations across independent seeds.  The
DPI-family entries have one fixed official run and provide no seed-stability
estimate.  DPI's value at `x=0` is reported only as a diagnostic because its
published benchmark is the distributed function task.

## Fair task and cost separation

Native single-point and amortized function approximation are kept in separate
tables.  Deep BSDE is compared at its published initial point and would require
an independent training run for every additional state.  Tangent-Actor, DPI,
PINN-HTE, D-DBSDE, and SDGD are evaluated on common held-out state banks where
their native models are amortized.  Query-count cost curves charge all Actor
training and measured inference; they do not amortize training away at one
query.

For Tangent-Actor, the main solve cost includes label generation, Actor fitting,
and final Monte Carlo policy evaluation.  Analytic-gradient Actor diagnostics,
Bregman checks, and artifact I/O are recorded separately and also included in
an auxiliary all-diagnostics total.  MCV factors are reported as estimator
variance-reduction factors, never as RMSE-reduction factors.

## Reference-free Actor tuning before final-test evaluation

The staged search used no analytic nonterminal values or gradients.  The
selection statistic was a replicated common-state linear-policy value with
common Brownian paths across candidates; held-out tangent-label error broke
ties.  Plain-FK training was retained as an ablation but was not eligible to
replace the declared MCV-training method.

The ultimately selected `focused` configuration uses 8,192 spatial training
points, 128 paths per label, 32 label-time steps, a four-hidden-layer
width-384 SiLU Actor, and 5,000 optimizer steps in each of four tangent updates.
Successive Actors warm-start from the preceding round.  Direct MCV with
coefficient one is used for eligible gradient-policy labels and final
inference; Allen–Cahn retains its Value-only Actor and does not receive an
extra Gradient head merely to enable MCV.

| Benchmark | Frozen staged checkpoint | Training time | Common-bank policy value | Held-out label RMSE | Label-replication RMSE |
|---|---|---:|---:|---:|---:|
| HJB-Log | `focused a4` | 993.72 s | 4.575605754 | 0.006573 | 0.000847 |
| OU/GMM HJB | `focused a4` | 1539.59 s | 308.7257452 | 0.059290 | 0.009644 |
| Allen–Cahn | `focused a4` | 635.78 s | 0.04632316949 | 0.001436 | 0.000603 |

For HJB-Log, `focused a4` beat `focused a3` by `6.70e-5` with paired
state SE `2.91e-5`.  A nominally larger `max_quality a1` used 16,384 points,
256 paths per label, and 64 label steps, yet was worse by `3.07e-4` with
paired SE `1.01e-4`; it was therefore stopped before a second round.

For OU/GMM, the fourth focused update beat the third by `0.027388` with
paired state SE `0.003114`.  At three matched rounds, MCV and plain-FK
training had essentially the same policy value (`308.75313` versus
`308.75442`).  MCV nevertheless reduced independent label-replication RMSE
from `0.029624` to `0.010465` and modestly reduced held-out error from
`0.077291` to `0.071546`.  This is evidence for cleaner supervision, while
the matched-round policy gain is negligible.  One upper-budget OU/GMM update
was admitted by the predeclared gate because iteration four still improved
and Actor error remained well above label noise.  It was worse than
`focused a4` by `1.124919 ± 0.124412` on the replicated bank and worsened
held-out RMSE to `0.250733`, despite reducing label-replication RMSE to
`0.008390`.  The max-budget branch therefore stopped before `a2`; the final
OU/GMM checkpoint remains `focused a4`.

For Allen–Cahn, `focused a4` beat `focused a3` by only `1.47e-8` with
paired state SE `2.31e-9`.  The leading `focused`, `quality128`, `data8192`,
and `deep_elu` third/fourth-round checkpoints all lie within `4.5e-8` in
common-bank policy value.  The selected held-out error was only `2.38` times
the label-replication RMSE, below the predeclared factor-five upper-budget
trigger.  The `max_quality` branch was therefore not run, and `focused a4`
was frozen before the final-test reference was opened.

## Frozen Tangent-Actor evaluation

The table reports the mean and sample standard deviation across three fully
independent `focused a4` training seeds.  Function inference uses 64 Euler
steps, 256 paths per state, and four independent replicates.  Actor error is
vector RMSE against `z=sigma^T grad(u)` for the two gradient problems and
direct Value-Actor RMSE for Allen–Cahn.

| Benchmark | Value RMSE | Relative L2 | Bias | Actor error | MCV variance reduction | Training | Function inference |
|---|---:|---:|---:|---:|---:|---:|---:|
| HJB-Log | `0.0008535 ± 0.0000365` | `0.0001862 ± 0.0000080` | `+0.0003478 ± 0.0000218` | `0.04716 ± 0.00206` | `18.83 ± 0.93 x` | `1030.3 ± 75.9 s` | `101.6 ± 6.1 s` |
| OU/GMM HJB | `0.22095 ± 0.00149` | `0.0007127 ± 0.0000048` | `-0.18458 ± 0.00128` | `0.59796 ± 0.00654` | `44.37 ± 0.35 x` | `1382.3 ± 136.8 s` | `134.2 ± 1.9 s` |
| Allen–Cahn | `0.0001666 ± 0.0000083` | `0.003636 ± 0.000182` | `-0.0000005 ± 0.0000200` | `0.001347 ± 0.000041` | not used | `580.2 ± 63.8 s` | `44.3 ± 1.7 s` |

There were no nonfinite values.  Mean exponential-weight ESS fractions were
`1.000`, `1.000`, and `0.9924`, respectively.  On the official OU/GMM bank,
the frozen 64-step Actor RMSE is far below DPI (`17.6638`), PINN-HTE
(`38.6698`), and D-DBSDE (`30.7247`).  Those competitors have one fixed
official seed, so the magnitude is compelling but is not presented as a
multi-seed significance test.

The native published-point comparison remains a different task:

| Benchmark and method | Native absolute error | Training | Native inference | Interpretation |
|---|---:|---:|---:|---|
| HJB-Log Tangent-Actor | `0.003150 ± 0.000207` | `1030.3 ± 75.9 s` | `25.84 ± 1.14 s` | more accurate, substantially slower at one query |
| HJB-Log Deep BSDE | `0.009678 ± 0.002761` | `65.15 s` | included in training output | native single-point baseline |
| Allen–Cahn Tangent-Actor | `0.0000295 ± 0.0000241` | `580.2 ± 63.8 s` | `26.20 ± 1.04 s` | mean MC SE is `0.0000467`; reference-limited |
| Allen–Cahn Deep BSDE | `0.0001691 ± 0.0001798` | `127.69 s` | included in training output | seed spread exceeds its mean error |

Thus Tangent-Actor is more accurate on the HJB-Log native point but does not
win on single-query runtime.  Its training cost is amortized only when the
learned function is queried repeatedly.  On Allen–Cahn the Actor mean error is
lower, but its MC standard error and the independent-reference uncertainty
prevent a fine-grained ranking of very small point errors.

### MC, time discretization, and Actor error

At 64 steps the mean sampling-variance MSE components are `5.37e-7` for
HJB-Log, `4.26e-4` for OU/GMM, and `2.76e-8` for Allen–Cahn.  Subtracting
the sampling component from each seed's total MSE, truncating at zero, and
then averaging gives non-MC RMSE proxies `0.000435`, `0.21999`, and
`0.0000284`.  The Allen–Cahn decomposition is itself noisy: one seed has an
estimated sampling component just above its observed MSE, and pooling the
three seed MSEs before subtraction instead gives `0.0000127`.  HJB-Log and
Allen–Cahn are substantially MC-limited at the declared path budget, whereas
OU/GMM is dominated by systematic error.

The frozen seed-0 discretization diagnostic is:

| Benchmark | 32-step RMSE | 64-step RMSE | 128-step RMSE | 32 / 64 / 128 function inference |
|---|---:|---:|---:|---:|
| HJB-Log | `0.000911` | `0.000814` | `0.000706` | `50.7 / 99.3 / 183.0 s` |
| OU/GMM HJB | `0.47401` | `0.21976` | `0.08768` | `73.4 / 136.4 / 253.9 s` |
| Allen–Cahn | `0.0001717` | `0.0001734` | `0.0001727` | `22.2 / 43.9 / 81.9 s` |

OU/GMM bias changes from `-0.39888` to `-0.18354` to `-0.06888`; its current
accuracy ceiling is therefore controlled mainly by Euler discretization.
The 128-step result is reported as a post-freeze numerical diagnostic, not as
a replacement selected after seeing test error.  HJB-Log benefits modestly
from finer steps.  Allen–Cahn is unchanged within MC noise.

For HJB-Log, the 64-step Bregman diagnostic has observed mean policy gap
`0.000279`, integrated quadratic gap `0.000241`, and identity residual
`3.8e-5`.  For OU/GMM these are `-0.1860`, `0.04987`, and `-0.2359` at 64
steps; the residual contracts from `-0.4443` at 32 steps to `-0.1187` at 128
steps.  The sign violation is a finite-step bias diagnostic rather than a
failure of the continuous-time identity.  The quadratic-gradient identity
does not apply to Allen–Cahn's cubic value driver.

### Postselection audit

The descriptive, lower-path test screen agrees with the frozen HJB-Log choice:
`focused a4` is best.  On OU/GMM it ranks `quality128 a3` first with RMSE
`0.1367`, while frozen `focused a4` has `0.2266` under that noisy screen.  On
Allen–Cahn it ranks `data8192 a1` first, but all leading candidates agree to
about `1e-8` in screen RMSE.  These outcomes do not revise the declared
reference-free selection.  The OU/GMM mismatch is an honest limitation of the
policy-value selection statistic and is retained in the descriptive compute
frontier rather than hidden by selecting on test RMSE.

## Numerical integrity and leakage controls

Every Actor consumes time and all 100 ambient state coordinates.  The model
receives no norm, active-subspace projection, GMM responsibility, exact
nonterminal value, or exact nonterminal gradient.  The terminal value or
terminal gradient appears only in our declared boundary-aware residual
parameterization.  Competing methods retain their own original terminal
construction and receive none of ours.

The HJB implementations consistently use `z=sigma^T grad(u)`.  For HJB-Log,
`sigma=sqrt(2) I`, so the tangent driver is `-|z|^2/2`; for OU/GMM the known
linear drift stays in the forward SDE and only the concave quadratic driver is
linearized.  Allen–Cahn uses a Value Actor constrained by its comparison-
principle interval `[0,1]`; no Gradient head is added just to manufacture an
MCV.  Nine automated tests cover terminal data, diffusion normalization,
sign conventions, exact-reference residuals, and tangent identities.

Ground Truth is absent from label generation, fitting, early stopping, Actor
iteration, and frozen checkpoint ranking.  Analytic values and gradients are
opened only for the final evaluation and Bregman diagnostics.  The later
all-checkpoint test screen is marked descriptive and cannot revise the frozen
selection.

## Uncertainty and comparison rules

Headline Tangent-Actor results use three independently trained Actors when
the frozen configuration is affordable.  Each frozen Actor evaluation uses
four independent MC replicates on the common held-out states.  The report
separates the replicate sampling-variance component from the residual
non-MC error proxy, and it reports seed variation separately from conditional
MC standard error.

Deep BSDE has three independent seeds on its native single-point tasks.
DPI, PINN-HTE, D-DBSDE, and the full official SDGD run each have one fixed
official seed, so their differences cannot be presented as multi-seed
statistical dominance.  Allen–Cahn comparisons within roughly `1e-5` of the
reference are labelled reference-limited.  A variance-reduction factor only
compares estimator variances on shared paths; it is never interpreted as a
factor reduction in Value RMSE or systematic bias.
