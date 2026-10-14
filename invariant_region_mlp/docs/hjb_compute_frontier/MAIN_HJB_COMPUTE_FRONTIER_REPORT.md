# Teacher-free Actor Iteration versus strong HJB baselines

## Scope and status

This experiment compares ambient-output Teacher-free Actor Iteration with a
martingale control variate against PINN+SDGD, an amortized Deep BSDE solver,
and Deep Picard Iteration (DPI).  The two frozen benchmarks are B1
`(d=400,k=40,T=0.25)` and HJB-Log `(d=400,T=1)`.  Every reported main result
uses a new shared 2,048-state test bank; the historical 960-state B1 bank is
not reused.

The main comparison consists of three independent training seeds for every
method and PDE.  All 24 final jobs completed on NVIDIA GH200 GPUs without a
retry.  The external baselines use one configuration each and never use an
analytic solution for tuning or checkpoint selection.  Following the user's
explicit override, the Actor may use final-test RMSE to rank the five declared
Actor candidates.  Those rows are therefore an oracle development result,
and search cost is reported separately.  B1 happens to select `q3` under both
the test oracle and the independent policy-value validation rule in all three
seeds.  HJB-Log selection is not stable because its policy differences are
smaller than final Monte Carlo noise.

The main numerical result is mixed and useful:

- On B1, the final `q3` Actor has the lowest mean Value RMSE, `0.06301`, but
  its small `0.00540` advantage over DPI is uncertain with three seeds and
  costs much more computation.  DPI is the strongest compute-efficient B1
  baseline.
- On HJB-Log, even the fixed, inexpensive `q0` Actor reaches mean RMSE
  `0.005017`, about 4.06 times lower than Deep BSDE and 12.7 times lower than
  DPI.  Later Actor iterations improve the direct gradient but barely change
  Value accuracy.  The final policy evaluation, rather than Actor training,
  is the visible accuracy bottleneck.
- A separate high-quality B1 seed-0 diagnostic reaches RMSE `0.03307` with a
  refreshed one-step Actor.  A matched true `q2` and a later `q3` do not beat
  it, so the best tuned result comes from better supervision and continued
  fitting rather than additional policy iterations.  This test-oracle,
  single-seed result is evidence of attainable accuracy, not a replacement
  for the frozen comparison.

## Frozen methods and selection

The exact settings were committed before the final array in
[`protocol_frozen.json`](../../experiments/hjb_compute_frontier/protocol_frozen.json).
The external adaptations and their official source links are documented in
[`baseline_protocols.md`](baseline_protocols.md).

| Method | Frozen final configuration | Selection without analytic test reference |
|---|---|---|
| PINN+SDGD | 4x512 Tanh, coordinate budget 32, 10,000 Adam steps | minimum full 400-coordinate PDE-residual MSE on 8 frozen states |
| Amortized Deep BSDE | shared time-dependent 3x256 SiLU `z`, 32 time steps, 10,000 Adam steps | minimum terminal-mismatch MSE on frozen paths |
| DPI | terminal-aware 4x512 ELU joint Value/gradient net, 20 rounds, 4,096 points and 512 paths per round | preset round-20 ceiling; held-out round-local target loss is diagnostic |
| Actor + CV | ambient 400-output 4x384 SiLU, `q0`, `q1`, matched `q1_refresh`, true `q2`, true `q3` | user-authorized minimum test RMSE; independent policy validation also recorded |

There is no rank projection in either Actor benchmark.  Actor labels use only
the terminal condition and pathwise derivatives of a frozen controlled
trajectory.  They do not use a nonlinear gradient teacher or an analytic
nonterminal gradient.  The final estimator is

\[
\widehat u^q = \operatorname{MC}\!\left[
g(X_T^q)+\frac12\int_0^T\|q(s,X_s^q)\|^2ds
-\int_0^Tq(s,X_s^q)^\top dW_s\right],
\]

with 256 antithetic paths and 64 Euler steps in the frozen comparison.

## Final three-seed accuracy

Values are mean +/- sample standard deviation over three training seeds.  A
single 2,048-point reference bank is intentionally shared across methods, so
the standard deviation measures training and evaluation-seed variability,
not test-state resampling.

### B1 `(400,40)`

| Method | selected checkpoint | Relative L2 | Value RMSE | bias | RMSE / std(`u*`) | gradient relative error |
|---|---:|---:|---:|---:|---:|---:|
| **Actor + CV** | `q3` | **0.02059 +/- 0.00170** | **0.06301 +/- 0.00521** | +0.04155 +/- 0.00127 | **0.1032 +/- 0.0085** | 0.3043 +/- 0.0039 |
| PINN+SDGD | 10,000 | 0.02698 +/- 0.00380 | 0.08254 +/- 0.01163 | +0.02208 +/- 0.01065 | 0.1352 +/- 0.0190 | **0.1342 +/- 0.0147** |
| Amortized Deep BSDE | 10,000 | 0.10345 +/- 0.00237 | 0.31655 +/- 0.00724 | +0.12827 +/- 0.00139 | 0.5184 +/- 0.0119 | 0.5848 +/- 0.0060 |
| DPI | round 20 | 0.02236 +/- 0.00061 | 0.06841 +/- 0.00186 | -0.00063 +/- 0.00099 | 0.1120 +/- 0.0030 | 0.5394 +/- 0.0143 |

The paired Actor-minus-DPI RMSE difference is `-0.00540`.  Its exact
small-sample 95% Student-t interval is `[-0.02023,+0.00944]`, so three seeds do
not establish a reliable advantage over DPI.  The Actor is lower than
PINN+SDGD in every seed, but that paired interval is also marginal at this
sample size: `-0.01952 +/- 0.01984`.  The comparison supports competitive
accuracy, not a broad statistical superiority claim.

### HJB-Log `(400)`

| Method | selected checkpoint | Relative L2 | Value RMSE | bias | RMSE / std(`u*`) | gradient relative error |
|---|---:|---:|---:|---:|---:|---:|
| **Actor + CV, test oracle** | `q3,q0,q3` | **0.0007780 +/- 0.0000070** | **0.004974 +/- 0.000045** | +0.002210 +/- 0.000133 | **0.2077 +/- 0.0019** | 0.3901 +/- 0.0509 |
| PINN+SDGD | `5k,5k,10k` | 0.13171 +/- 0.01129 | 0.84205 +/- 0.07219 | +0.83831 +/- 0.06962 | 35.16 +/- 3.01 | 5.033 +/- 3.597 |
| Amortized Deep BSDE | 10,000 | 0.003183 +/- 0.000159 | 0.02035 +/- 0.00102 | +0.00128 +/- 0.00304 | 0.8497 +/- 0.0425 | 1.0080 +/- 0.0025 |
| DPI | round 20 | 0.009977 +/- 0.000807 | 0.06378 +/- 0.00516 | -0.03243 +/- 0.01026 | 2.663 +/- 0.216 | 11.948 +/- 0.247 |

The fixed `q0` Actor, before any policy-improvement iteration, has RMSE
`0.005017 +/- 0.000051`.  This differs little from the test-oracle result and
avoids candidate-selection ambiguity.  Relative to fixed `q0`, Deep BSDE has
about 4.06 times higher RMSE and DPI about 12.7 times higher RMSE.

HJB-Log has truth mean about `6.393` and truth standard deviation only
`0.02395`, so Relative L2 alone is misleading.  A constant fitted on the
independent validation reference has test RMSE `0.02398`.  The Actor's
`RMSE/std(u*)` near `0.21` and its absolute RMSE show that it learns state
variation rather than merely predicting the mean.

### Frozen-Actor evaluation sensitivity

The five final Actor candidates were subsequently frozen and evaluated with
two independent budgets.  This diagnostic performs no training.  Candidate
comparisons within each seed and budget use common Brownian random numbers.

| evaluation | fixed `q0` RMSE | fixed `q3` RMSE | test-oracle actors | oracle RMSE | `q3` inference, 2,048 states |
|---|---:|---:|---|---:|---:|
| 256 paths, 64 steps (frozen main) | 0.005017 +/- 0.000051 | 0.004996 +/- 0.000051 | `q3,q0,q3` | 0.004974 +/- 0.000045 | 16.86 s |
| **1,024 paths, 64 steps** | 0.003233 +/- 0.000037 | **0.003094 +/- 0.000092** | `q3,q3,q3` | **0.003094 +/- 0.000092** | 65.60 s |
| 256 paths, 128 steps | 0.005067 +/- 0.000155 | 0.004963 +/- 0.000082 | `q2,q1_refresh,q3` | 0.004937 +/- 0.000087 | 32.36 s |

Quadrupling paths lowers fixed-`q3` RMSE by 38.1%, whereas doubling Euler
steps at fixed path count changes it by only 0.66%.  At 1,024 paths, the
estimated MC variance still accounts for 51.9% of observed `q3` MSE; the
variance-subtracted RMSE diagnostic is `0.002147`.  The corresponding fixed
`q0` values are 46.6% and `0.002361`.  Thus `q3` has a small detectable policy
benefit under more accurate evaluation, but most of the practical gain comes
from evaluating a cheap Actor with more paths.  The high-path `q3` RMSE is
only 4.3% below `q0` despite roughly 75 times more Actor path-training time.
These more expensive inference points are a sensitivity analysis and are not
substituted into the frozen four-method compute curves.

## Computation

These are measured method-internal GH200 wall times, mean +/- sample standard
deviation.  Total cost adds inference for 960 states to the complete path cost
for the selected model.  Actor label generation, supervised optimization,
policy validation, and Monte Carlo inference are included.

### Selected model path

| PDE | Method | training / model-path seconds | inference, 960 states | total seconds |
|---|---|---:|---:|---:|
| B1 | Actor `q3` | 978.7 +/- 49.0 | 8.83 +/- 1.74 | 987.6 +/- 50.7 |
| B1 | PINN+SDGD | 1554.6 +/- 40.9 | 0.00145 | 1554.6 +/- 40.9 |
| B1 | Deep BSDE | 345.5 +/- 23.1 | 0.00057 | 345.5 +/- 23.1 |
| B1 | DPI | **36.74 +/- 1.12** | 0.00097 | **36.74 +/- 1.12** |
| HJB-Log | fixed Actor `q0` | **13.21** | 7.86 | **21.07** |
| HJB-Log | Actor test oracle | 637.3 +/- 542.1 | 7.95 +/- 0.56 | 645.2 +/- 542.4 |
| HJB-Log | PINN+SDGD | 1718.4 +/- 130.2 | 0.00131 | 1718.4 +/- 130.2 |
| HJB-Log | Deep BSDE | 339.6 +/- 3.9 | 0.00052 | 339.6 +/- 3.9 |
| HJB-Log | DPI | 25.42 +/- 0.66 | 0.00063 | 25.43 +/- 0.66 |

The HJB-Log fixed `q0` result is therefore both more accurate and cheaper than
the three baselines for a 960-query workload.  On B1, DPI reaches nearly the
same final accuracy roughly 27 times faster than the Actor model path.  The
Actor's B1 advantage is final accuracy, not compute efficiency.

The Actor oracle requires evaluating all candidates.  Mean total exploration
plus candidate-evaluation cost is `1395.4` seconds per B1 seed and `1351.8`
seconds per HJB-Log seed.  This search overhead is separate from the cost of a
single selected model and must be counted in any end-to-end tuning comparison.

The complete frozen final array consumed 21,197 allocated GH200-seconds
(`5.8881` GPU-hours).  Smoke tests, one diagnosed smoke failure, profiling,
the pilot, and the final array together consumed 22,459 seconds (`6.2386`
GPU-hours).  Frozen HJB-Log evaluation sensitivity added 1,510 seconds, and
the staged high-quality B1 diagnostic added 4,452 seconds.  Total recorded
allocation for the complete study is 28,421 GH200-seconds (`7.8947`
GPU-hours), including failed-job and scheduler overhead.

## Compute-accuracy curves

The figures use validation-eligible checkpoints at each budget for the three
external methods.  They never take a test-error lower envelope.  Actor points
show the declared discrete candidates and therefore retain the explicit test
oracle label.

| PDE | method/checkpoint | mean training seconds | mean Relative L2 | mean RMSE |
|---|---|---:|---:|---:|
| B1 | DPI round 2 | 4.71 | 0.02938 | 0.08990 |
| B1 | DPI round 20 | 36.74 | 0.02236 | 0.06841 |
| B1 | Actor `q1` | 335.3 | 0.02923 | 0.08945 |
| B1 | Actor `q1_refresh` | 655.9 | 0.02329 | 0.07126 |
| B1 | Actor `q3` | 978.7 | 0.02059 | 0.06301 |
| B1 | PINN 2k / 10k | 311.3 / 1554.6 | 0.02850 / 0.02698 | 0.08720 / 0.08254 |
| HJB-Log | Actor `q0` | 13.21 | 0.0007848 | 0.005017 |
| HJB-Log | DPI round 20 | 25.42 | 0.009977 | 0.06378 |
| HJB-Log | Deep BSDE 10k | 339.6 | 0.003183 | 0.02035 |
| HJB-Log | Actor `q3` | 952.2 | 0.0007814 | 0.004996 |

The rendered plots are:

- [combined 2x2 frontier](../../results/hjb_compute_frontier/figures/hjb_compute_frontier_2x2.png)
- [B1 training-time curve](../../results/hjb_compute_frontier/figures/b1_d400_k40_training_time.png)
- [B1 total-cost curve](../../results/hjb_compute_frontier/figures/b1_d400_k40_total_cost.png)
- [HJB-Log training-time curve](../../results/hjb_compute_frontier/figures/hjb_log_d400_training_time.png)
- [HJB-Log total-cost curve](../../results/hjb_compute_frontier/figures/hjb_log_d400_total_cost.png)

## What the Actor iterations show

| PDE | Actor | Value RMSE, mean +/- SD | gradient relative error | path training seconds |
|---|---|---:|---:|---:|
| B1 | `q0` | 0.16866 +/- 0.00211 | 0.6929 | 12.6 |
| B1 | `q1` | 0.08945 +/- 0.00033 | 0.4621 | 334.7 |
| B1 | `q1_refresh` | 0.07126 +/- 0.00186 | 0.3593 | 655.3 |
| B1 | true `q2` | 0.08610 +/- 0.00345 | 0.3587 | 655.8 |
| B1 | true `q3` | **0.06301 +/- 0.00521** | **0.3043** | 978.1 |
| HJB-Log | `q0` | 0.005017 +/- 0.000051 | 0.4500 | 12.7 |
| HJB-Log | `q1` | 0.005161 +/- 0.000016 | 0.3872 | 325.7 |
| HJB-Log | `q1_refresh` | 0.005058 +/- 0.000041 | 0.3700 | 639.5 |
| HJB-Log | true `q2` | 0.005060 +/- 0.000046 | 0.3701 | 638.1 |
| HJB-Log | true `q3` | 0.004996 +/- 0.000051 | **0.3604** | 951.7 |

On B1, `q2` is worse than the matched fresh-label continuation
`q1_refresh`, while `q3` then improves substantially in all seeds.  Approximate
policy iteration is therefore nonmonotone, and one low-quality `q2` result
would have led to the wrong stopping conclusion.  On HJB-Log, gradient error
falls steadily but Value RMSE remains near `0.005`; extra policy iterations
do not provide a worthwhile accuracy-cost tradeoff.

### Targeted high-quality B1 diagnostic

One seed was then used to raise supervision quality without changing the
ambient-output algorithm.  The width-512 `q1`, `q1_refresh`, and true `q2`
phases each use 16,384 training points, 1,024 paths per label, 64 label Euler
steps, and 12,000 optimizer steps.  `q1_refresh` and `q2` use the same states
and compute budget; only the policy-evaluation source differs.  `q3` uses the
smaller 8,192/512/32/6,000 follow-up budget.

| Actor | Value RMSE | Relative L2 | gradient relative error | model-path seconds |
|---|---:|---:|---:|---:|
| `q1` | 0.053064 | 0.017343 | 0.3339 | 1,260.1 |
| **`q1_refresh`** | **0.033069** | **0.010808** | 0.2450 | 2,536.8 |
| true `q2` | 0.040800 | 0.013334 | 0.2305 | 2,493.6 |
| true `q3` | 0.036694 | 0.011992 | **0.2151** | 2,667.8 |

The matched true `q2` has 23.4% higher Value RMSE than continued fitting with
fresh `q0` labels.  `q3` recovers 10.1% from `q2` but remains 11.0% above
`q1_refresh`.  MC-variance-subtracted diagnostics preserve this ordering.  The
test oracle therefore selects `q1_refresh`; reference-free policy validation
instead ranks `q3` first, so the small single-seed gap is selection-rule
sensitive.  A `q4` was not justified.  The high-quality `q1_refresh` is 42.1%
below the same-seed frozen `q3` RMSE and demonstrates that label quality and
continued fitting are the larger B1 lever.

## Numerical and artifact audit

The HJB-Log reference uses a float64 one-dimensional Hopf-Cole quadrature after
the rescaling `y=d*s`.  Independent adaptive integration agrees within
`2.30e-13` for Value and `1.53e-15` for `z`; B1's analytic PDE residual is at
most `2.22e-16`.  Seven focused tests cover terminal constraints, HJB signs,
SDGD second derivatives, Deep BSDE drift, controlled pathwise differentiation,
and the martingale correction.

Independent aggregation reloaded all 96 final prediction artifacts and
recomputed the scalar metrics.  The maximum discrepancy from the stored run
records is `1.42e-14`.  The five high-quality Actor prediction arrays were
audited separately with maximum discrepancy `1.67e-16`.  Reproducibility
artifacts include:

- [`summary.csv`](../../results/hjb_compute_frontier/summary.csv)
- [`compute_frontier.csv`](../../results/hjb_compute_frontier/compute_frontier.csv)
- [`prediction_audit.json`](../../results/hjb_compute_frontier/prediction_audit.json)
- [`evaluation_sensitivity_summary.csv`](../../results/hjb_compute_frontier/evaluation_sensitivity_summary.csv)
- [`evaluation_sensitivity_aggregate.json`](../../results/hjb_compute_frontier/evaluation_sensitivity_aggregate.json)
- [`high_actor_summary.csv`](../../results/hjb_compute_frontier/high_actor_summary.csv)
- [`high_actor_prediction_audit.json`](../../results/hjb_compute_frontier/high_actor_prediction_audit.json)
- [`validation.json`](../../results/hjb_compute_frontier/validation.json)
- [`scheduler_accounting.json`](../../results/hjb_compute_frontier/scheduler_accounting.json)
- [`EXPERIMENT_LOG.md`](../../results/hjb_compute_frontier/EXPERIMENT_LOG.md)

## Research conclusion and limits

The strongest paper evidence is benchmark-dependent.  HJB-Log shows a clear
accuracy-and-cost win for terminal-pathwise Actor initialization plus
martingale-CV evaluation: the fixed `q0` result is already substantially
better than all three baselines.  B1 shows that high-quality pathwise
supervision can attain much lower single-seed Value error and that true policy
iterations can recover after a poor intermediate update.  The matched high
budget comparison also shows that a true `q2` does not beat refreshed
one-step training.  DPI remains far more compute-efficient and statistically
tied to the frozen Actor result at three seeds.

The main limits are the explicit test-oracle Actor selection, only three final
training seeds, one frozen configuration per external baseline, a scaled
rather than author-scale PINN/DPI budget, and the amortized adaptation required
for Deep BSDE.  The poor HJB-Log PINN result should be read as the result of
this documented single preset, not as a universal statement about PINNs.
Likewise, one targeted high-budget B1 seed is a tuning diagnostic rather than
a replacement for the three-seed frozen comparison.  A prepared seeds 3--4
extension was not submitted because the project QOS had only 0.591 one-GH200
minutes left after the final diagnostic; no five-seed claim is made.
