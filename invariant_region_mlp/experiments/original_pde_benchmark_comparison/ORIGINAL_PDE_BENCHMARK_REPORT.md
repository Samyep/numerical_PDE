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

All completed numerical runs use one NVIDIA GH200 120 GB GPU.  Reported
training and inference times are synchronized wall times.  Failed and
restarted allocations remain in the resource audit, including one accidental
four-GPU request cancelled after four seconds before numerical work.

## Benchmarks

### Table 1 — PDE benchmarks

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
| SDGD-PINN | HJB-Log common 960-state bank | 1 | `0.0299212` RMSE | 11,402.80 s |

The `±` entries are sample standard deviations across independent seeds.  The
DPI-family and SDGD-PINN entries have one fixed official run and provide no
seed-stability estimate.  DPI's value at `x=0` is reported only as a diagnostic
because its published benchmark is the distributed function task.

The official SDGD-PINN run completed all 10,000 published epochs.  Its
otherwise monolithic loop was split at epoch 5,000 while preserving model,
Adam, scheduler, collocation points, and every random-number state; a 100-epoch
audit matched the uninterrupted official loop bit-for-bit.  The 10,000-epoch
checkpoint contains all 10,000 loss and L2 entries.  On the full official
20,000-state test it obtains RMSE `0.0276622`, relative L2 `0.00606570`, and
bias `+0.0157126`.  Direct evaluation on the common 960-state bank gives RMSE
`0.0299212`, relative L2 `0.00652670`, bias `+0.0177013`, and Gradient-Actor
vector RMSE `0.0247201`.  Its native `u(0,0)` prediction is `2.5220323` versus
the independent reference `4.5901617`, an absolute error of `2.06813`.
Distributed and single-point results are therefore reported separately.

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

### Table 2 — Native-task accuracy

The native published-point comparison remains a different task.  The
Full-History row uses the member with the lowest function-bank RMSE in its
frozen grid; all six configurations are reported below, so this descriptive
label is not presented as a preregistered selection.

| Benchmark and method | Native absolute error | Training | Native inference | Interpretation |
|---|---:|---:|---:|---|
| HJB-Log Tangent-Actor | `0.003150 ± 0.000207` | `1030.3 ± 75.9 s` | `25.84 ± 1.14 s` | more accurate, substantially slower at one query |
| HJB-Log Deep BSDE | `0.009678 ± 0.002761` | `65.15 s` | included in training output | native single-point baseline |
| HJB-Log SDGD-PINN | `2.06813` | `11402.80 s` | `0.00068 s` | official amortized model evaluated at the native point |
| HJB-Log Full-History `n3/M24` | `0.001647 ± 0.00175` | none | `0.00383 ± 0.00007 s` | inference-only; large three-seed point-estimator spread |
| Allen–Cahn Tangent-Actor | `0.0000295 ± 0.0000241` | `580.2 ± 63.8 s` | `26.20 ± 1.04 s` | mean MC SE is `0.0000467`; reference-limited |
| Allen–Cahn Deep BSDE | `0.0001691 ± 0.0001798` | `127.69 s` | included in training output | seed spread exceeds its mean error |
| Allen–Cahn Full-History `n3/M24` | `0.0001969 ± 0.000220` | none | `0.00309 ± 0.00007 s` | inference-only; point estimate is reference- and sampling-limited |

Thus Tangent-Actor is more accurate than Deep BSDE and SDGD-PINN on the
HJB-Log native point but does not win on single-query runtime.  The displayed
Full-History mean is lower still, with a standard deviation slightly larger
than its mean.  Actor training cost is amortized only when the learned
function is queried repeatedly.  On Allen–Cahn the Actor mean error is lower,
but its MC standard error and the independent-reference uncertainty prevent a
fine-grained ranking of very small point errors.

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

### What tuning changed

Iteration count was the most consequential Actor hyperparameter on OU/GMM:
the replicated paired improvement fell from `18.77` on the first update to
`1.077`, `0.0879`, and `0.0277` on updates two through four, but remained
positive for all three seeds.  On HJB-Log, later-round gains were much smaller
(`1.42e-4 ± 8.07e-5` for update four), and Allen–Cahn was already on a
practical plateau after update one.

Width/depth mattered mainly during the low-budget HJB searches.  Increasing
only the number of spatial points or only the label-path count did not
reliably improve the reference-free ranking.  The upper-budget experiments
make this especially clear: doubling states and paths and using 64 label steps
reduced label-replication noise on HJB-Log and OU/GMM, yet worsened their
policy rankings.  Cleaner labels alone therefore do not solve Actor
approximation or optimization error.  MCV training likewise reduced label
noise without a resolved matched-round policy advantage over plain-FK
training.

Extensive tuning did materially lower error relative to the initial pilots,
especially on OU/GMM, but it did not remove every systematic limitation.  The
frozen OU/GMM checkpoint is sensitive to the reference-free selection
statistic, and its final 64-step bias is primarily an Euler error.  Those two
effects are more important than further increasing label paths at the tested
budgets.

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
MCV.  Eleven automated tests cover terminal data, diffusion normalization,
sign conventions, exact-reference residuals, tangent identities, the
value-only Full-History path, and exact SDGD random-state restoration.

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

## Frozen Full-History MLP control

The Full-History control has no learned parameters or training phase.  The
grid below was fixed before results were opened.  Entries are mean ± sample
standard deviation across three independent estimator seeds.  Function timing
is synchronized GPU time for the full common bank: 960 HJB-Log states and 410
Allen–Cahn states.  It excludes Python and `srun` startup, as do the timed
training and inference kernels for the neural methods.  The complete
36-process allocation, including every interpreter and scheduler launch, took
132 GPU seconds.

| Benchmark | `(n,M)` | Value RMSE | Relative L2 | Bias | Function inference | Native absolute error |
|---|---:|---:|---:|---:|---:|---:|
| HJB-Log | `(2,16)` | `0.007426 ± 0.000125` | `0.001620 ± 0.000027` | `-0.000123 ± 0.000385` | `0.1109 ± 0.0002 s` | `0.00235 ± 0.00231` |
| HJB-Log | `(2,32)` | `0.003819 ± 0.000052` | `0.000833 ± 0.000011` | `-0.000136 ± 0.000089` | `0.1086 ± 0.0031 s` | `0.00614 ± 0.00286` |
| HJB-Log | `(3,8)` | `0.005371 ± 0.000147` | `0.001172 ± 0.000032` | `-0.000825 ± 0.000111` | `0.1844 ± 0.0039 s` | `0.00936 ± 0.00729` |
| HJB-Log | `(3,12)` | `0.002917 ± 0.000113` | `0.000636 ± 0.000025` | `-0.000439 ± 0.000026` | `0.1825 ± 0.0045 s` | `0.00127 ± 0.00092` |
| HJB-Log | `(3,16)` | `0.001888 ± 0.000031` | `0.000412 ± 0.000007` | `-0.000404 ± 0.000060` | `0.2585 ± 0.0213 s` | `0.00249 ± 0.00212` |
| HJB-Log | `(3,24)` | `0.001073 ± 0.000016` | `0.000234 ± 0.000003` | `-0.000273 ± 0.000023` | `0.7689 ± 0.0161 s` | `0.00165 ± 0.00175` |
| Allen–Cahn | `(2,16)` | `0.001323 ± 0.000072` | `0.02887 ± 0.00158` | `-0.000607 ± 0.000043` | `0.0773 ± 0.0014 s` | `0.00216 ± 0.00099` |
| Allen–Cahn | `(2,32)` | `0.001124 ± 0.000020` | `0.02454 ± 0.00044` | `-0.000646 ± 0.000037` | `0.0780 ± 0.0004 s` | `0.00201 ± 0.00078` |
| Allen–Cahn | `(3,8)` | `0.000541 ± 0.000018` | `0.01182 ± 0.00039` | `-0.000036 ± 0.000022` | `0.1036 ± 0.0005 s` | `0.000612 ± 0.000540` |
| Allen–Cahn | `(3,12)` | `0.000389 ± 0.000019` | `0.00850 ± 0.00041` | `-0.000059 ± 0.000013` | `0.1072 ± 0.0052 s` | `0.000350 ± 0.000280` |
| Allen–Cahn | `(3,16)` | `0.000291 ± 0.000002` | `0.00635 ± 0.00004` | `-0.000042 ± 0.000013` | `0.1316 ± 0.0082 s` | `0.000343 ± 0.000502` |
| Allen–Cahn | `(3,24)` | `0.000205 ± 0.000005` | `0.00448 ± 0.00012` | `-0.000046 ± 0.000013` | `0.3171 ± 0.0081 s` | `0.000197 ± 0.000220` |

All 36 estimators were finite.  The best observed grid member is `(3,24)` on
both PDEs, but it is labelled descriptively because the protocol commits to
reporting the complete grid rather than selecting a configuration with Ground
Truth.  Tangent-Actor lowers its RMSE by only 20.5% on HJB-Log and 18.8% on
Allen–Cahn.  That accuracy increment costs orders of magnitude more: for
HJB-Log, `0.7689` seconds of Full-History inference becomes `1030.3 + 101.6`
seconds of Actor training and inference.  Even the entire two-PDE,
three-seed Full-History grid consumed only 132 allocated GPU seconds.  These
data contradict an efficiency advantage over Full-History MLP on the two
supported original PDEs.

## Consolidated fixed comparison

### Table 3 — Amortized function accuracy

The table uses the frozen Actor and fixed official baselines.  The
Full-History entry is the highest-accuracy observed member of the fully
reported frozen grid.  Total time is training plus inference on the indicated
bank and excludes diagnostic Ground-Truth gradients.

| Benchmark | Method | Seeds | States | Value RMSE | Relative L2 | Training | Inference | Total |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| HJB-Log | Full-History MLP `(3,24)` | 3 | 960 | `0.001073 ± 0.000016` | `0.000234 ± 0.000003` | none | `0.769 s` | `0.769 s` |
| HJB-Log | Tangent-Actor + MCV | 3 | 960 | `0.000854 ± 0.000037` | `0.000186 ± 0.000008` | `1030.3 s` | `101.6 s` | `1132.0 s` |
| HJB-Log | SDGD-PINN | 1 | 960 | `0.029921` | `0.006527` | `11402.8 s` | `0.187 s` | `11403.0 s` |
| OU/GMM HJB | Tangent-Actor + MCV | 3 | 960 | `0.22095 ± 0.00149` | `0.000713 ± 0.000005` | `1382.3 s` | `134.2 s` | `1516.4 s` |
| OU/GMM HJB | DPI | 1 | 960 | `17.6638` | `0.05697` | `627.4 s` | `0.119 s` | `627.5 s` |
| OU/GMM HJB | PINN-HTE | 1 | 960 | `38.6698` | `0.12472` | `767.2 s` | `0.120 s` | `767.3 s` |
| OU/GMM HJB | D-DBSDE | 1 | 960 | `30.7247` | `0.09910` | `453.0 s` | `0.117 s` | `453.2 s` |
| Allen–Cahn | Full-History MLP `(3,24)` | 3 | 410 | `0.000205 ± 0.000005` | `0.00448 ± 0.00012` | none | `0.317 s` | `0.317 s` |
| Allen–Cahn | Tangent-Actor | 3 | 410 | `0.000167 ± 0.000008` | `0.00364 ± 0.00018` | `580.2 s` | `44.3 s` | `624.6 s` |

Tangent-Actor is 35.1 times more accurate than SDGD-PINN on HJB-Log while
using about one tenth of its total time.  On OU/GMM it is 79.9 times more
accurate than DPI, 175.0 times more accurate than PINN-HTE, and 139.1 times
more accurate than D-DBSDE, but takes 2.4, 2.0, and 3.3 times their total time,
respectively.  These are one-run official neural baselines, so the ratios
describe the observed fixed runs rather than seed-level significance.

## Compute, memory, and stability

| Benchmark and method | Parameters | Peak GPU bytes | Nonfinite count |
|---|---:|---:|---:|
| HJB-Log Tangent-Actor | 521,188 | 138,283,520 | 0 |
| HJB-Log Full-History `(3,24)` | 0 | 138,526,720 | 0 |
| HJB-Log SDGD-PINN | 2,204,673 | 893,657,088 | 0 |
| HJB-Log Deep BSDE | 665,861 | 103,439,872 | 0 |
| OU/GMM Tangent-Actor | 521,188 | 155,096,064 | 0 |
| OU/GMM DPI | 967,624 | 98,693,219,840 | 0 |
| OU/GMM PINN-HTE | 967,624 | 2,464,704,512 | 0 |
| OU/GMM D-DBSDE | 967,624 | 3,565,502,464 | 0 |
| Allen–Cahn Tangent-Actor | 483,073 | 109,887,488 | 0 |
| Allen–Cahn Full-History `(3,24)` | 0 | 85,722,112 | 0 |
| Allen–Cahn Deep BSDE | 665,861 | 102,528,000 | 0 |

The complete study used 32 Slurm allocation records: 10.752 completed GH200
GPU-hours and 0.078 failed or cancelled GPU-hours.  The latter includes every
calibration failure, an accidental four-GPU launch cancelled after four
seconds, a failed evaluation launch, and the 20-second SDGD RNG-device resume
failure.  None is removed from `runtime_resource_accounting.csv`.  The
headline method times in Tables 2 and 3 are per solve; the allocation total
also includes broad Actor tuning, conditional trials, reference audits,
independent seeds, diagnostics, and failed starts.

For amortized workloads, Figure 2 evaluates `K=1,16,128,960` and uses
`E*(C)=min_{i:C_i<=C} E_i`.  Full-History owns the low-budget envelope on both
supported PDEs.  The frozen Actor enters the HJB-Log envelope only around
1,132 seconds, lowering RMSE from `0.001073` to `0.000854`; on Allen–Cahn it
enters around 684 seconds for `K=960`, lowering `0.000205` to `0.000167`.
The dashed Actor curve is the explicitly descriptive postselection screen and
does not alter the primary envelope.  Figure 1 places zero-training
Full-History points at a `0.01`-second display floor because its horizontal axis
is logarithmic; Figure 2 uses their actual measured total inference time.

Deep BSDE remains much faster for one native point: approximately 67 total
seconds for HJB-Log and 130 seconds for Allen–Cahn, versus roughly 1,056 and
606 seconds for Actor training plus native inference.  If one only extrapolates
cost and assumes independent Deep BSDE retraining at every point, batched Actor
cost crosses below Deep BSDE at about 16 HJB queries and five Allen–Cahn
queries.  This is a cost crossover, not a multi-query accuracy claim, because
the native Deep BSDE architecture supplies no common-bank function prediction.

## Answers to the research questions

1. **Eligible original PDEs.** HJB-Log and OU/GMM have globally concave
   quadratic gradient drivers; Allen–Cahn is concave on its rigorously
   invariant `[0,1]` interval.  All three are semilinear with a linear
   highest-order operator.  HJB-Rosenbrock has a compatible driver but was
   excluded because its official test data are missing and the available
   reference is not accurate enough for a headline comparison.

2. **Deep BSDE.** At the published point, Tangent-Actor has lower observed
   error on HJB-Log (`0.003150` versus `0.009678`) but costs much more for one
   query.  Its Allen–Cahn mean is also lower (`0.0000295` versus `0.0001691`),
   but Actor MC error, Deep BSDE seed spread, and the `9.64e-6` reference proxy
   make this a qualified result rather than a precise ranking.

3. **Deep Picard Iteration.** Yes on the official 960-state OU/GMM function
   task: `0.22095` versus `17.6638` RMSE.  The Actor run is about 2.4 times
   slower and its 64-step error is still dominated by Euler bias; 128 steps
   reduce its seed-0 RMSE further to `0.08768`.

4. **PINNs.** The frozen Actor beats SDGD-PINN on HJB-Log (`0.000854` versus
   `0.029921` RMSE) with about one tenth of the total time.  It also beats
   PINN-HTE on OU/GMM (`0.22095` versus `38.6698`) but costs about twice as
   much.  There is no original amortized PINN result for the selected
   Allen–Cahn protocol.

5. **GPU-time budgets.** Below one second of timed inference, Full-History
   reaches `0.001073` on HJB-Log and `0.000205` on Allen–Cahn.  The frozen
   Actor improves those errors modestly only after hundreds to more than one
   thousand seconds.  On OU/GMM the official neural baselines occupy the
   453–767 second range with RMSE 17.7–38.7; the frozen Actor uses about 1,516
   seconds and reaches `0.221`.  Lower-budget Actor checkpoints look strong in
   the descriptive test screen, but they were not selected with an admissible
   final-test rule and cannot replace the frozen headline.

6. **Single point versus multiple queries.** Deep BSDE wins single-query
   runtime; Tangent-Actor gives lower observed point error but does not hide
   its training overhead.  Its batched cost becomes preferable to retraining
   Deep BSDE near 16 HJB or five Allen–Cahn queries.  Full-History remains much
   cheaper across all tested query counts on the two PDEs it supports.

7. **Important Actor hyperparameters.** Iteration count matters most on
   OU/GMM, where all four updates give independently replicated positive
   policy improvements.  Width and depth matter mainly in the low-budget HJB
   search.  More points, paths, or label steps alone do not reliably improve
   the policy.  MCV sharply reduces label and inference variance, while its
   matched-round policy benefit over plain-FK training is unresolved.

8. **Systematic error after tuning.** Tuning produces stable, strong neural
   baseline results but does not remove systematic error.  HJB-Log retains a
   non-MC RMSE proxy of `0.000435`.  OU/GMM retains `0.21999`, driven mainly by
   time discretization and reference-free checkpoint ranking; increasing from
   64 to 128 Euler steps cuts seed-0 RMSE from `0.21976` to `0.08768`.  Cleaner
   upper-budget labels by themselves do not improve the chosen Actor.

9. **Efficiency relative to Full-History MLP.** The evidence is negative on
   these benchmarks.  Actor is only 20.5% and 18.8% more accurate than the
   strongest frozen Full-History result on HJB-Log and Allen–Cahn, while one
   Actor solve costs far more than the complete 36-run Full-History grid.  The
   current zero-drift recursive wrapper does not cover OU/GMM, so that PDE
   supplies no countervailing Full-History evidence.

10. **Defensible ICML-level claim.** The data support a conditional claim:
    terminal-aware Tangent-Actor with direct MCV is stable across three
    original 100D semilinear PDEs and substantially improves on fixed official
    amortized neural baselines in Value accuracy, with verified iterative
    policy improvement and useful variance reduction.  They do not support a
    general efficiency claim over recursive nonlinear Monte Carlo, nor a
    universal single-query advantage.  A paper centered on replacing
    Full-History MLP would need benchmarks where effective dimension,
    nonseparability, or recursion depth makes the matched Full-History control
    genuinely costly.  With the present evidence, that stronger positioning is
    a no-go; the neural-baseline and Actor-iteration result is a conditional go.

## Research conclusion

The study changes the proposed story.  Tangent-Actor + MCV is a strong
high-dimensional neural function approximator on the three eligible original
benchmarks.  It clearly outperforms the tested official DPI, PINN-HTE,
D-DBSDE, and SDGD-PINN runs in Value RMSE, and it gives lower mean native-point
error than Deep BSDE.  Its advantages have explicit costs and limitations:
single-point runtime is worse, Allen–Cahn is reference-limited, OU/GMM is
Euler-limited, and the raw Full-History estimator is already extremely strong
on HJB-Log and Allen–Cahn.  The appropriate outcome is therefore
**CONDITIONAL GO** for the Tangent-Actor method and **NO-GO** for the current
claim that it is a more efficient replacement for Full-History MLP.

## Generated artifacts

The publication figures are
[accuracy versus training time](figures/accuracy_vs_training_gpu_time.png),
[the true total-compute envelope](figures/total_compute_budget_envelope.png),
[Actor hyperparameter scaling](figures/tangent_actor_hyperparameter_scaling.png),
and [error/compute decomposition](figures/error_and_compute_decomposition.png).
Their vector PDF versions are stored beside the PNG files.  The underlying
rows are in `native_task_comparison.csv`,
`function_approximation_comparison.csv`, `hyperparameter_trials.csv`,
`training_compute_frontier.csv`, `total_compute_frontier.csv`, and
`runtime_resource_accounting.csv`.
