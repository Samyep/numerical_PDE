# High-Budget MCV Training Validation

**Decision: GO.** Direct \(\alpha=1\) martingale control variates should be the
default Tangent-Actor label estimator for this solver. Across both benchmarks
and both training seeds, MCV substantially reduced Value and Gradient label
noise, improved held-out Actor accuracy, increased reference-free policy
improvement, and lowered final PDE error at matched GH200 time. Complete-update
runtime was unchanged within run-to-run timing variation; peak label-generation
memory rose by 1.9% in \(d=100\) and 6.4% in \(d=400\).

The final PDE gain is smaller than the label and Actor-fit gain. At the B5
operating point, mean systematic RMSE fell by 1.36% on
`nonseparable_d100_k2` and 10.85% on `joint_concave_d400_k2` under one common
plain inference protocol. This distinction matters: MCV removes much more
supervision noise than the present Actor can convert into PDE accuracy, so
Actor approximation and time discretization become the next bottlenecks.

## Method and controls

For a frozen source Actor \(A_j=(r_j,q_j)\), training uses

\[
\widetilde P_j = P_j-\sum_{m=0}^{N-1} C_m q_j(t_m,X_m)^\top\Delta W_m,
\qquad
\widehat z_j=\sqrt2\nabla_x\mathbb E_M[\widetilde P_j].
\]

The implementation differentiates the complete corrected payoff. Source-Actor
parameters are frozen, while derivatives through \(C_m\), \(q_j(t_m,X_m)\),
and \(X_m\) are preserved. The correction is evaluated at predictable left
endpoints before either \(C\) or \(X\) is advanced. There is one source-Actor
forward pass per time step, just as in plain tangent labels.

All paired training branches start from the same saved `a3`. Architecture,
optimizer, loss, validation bank, checkpoint selection, 128 paths/label, and
16 label time steps are unchanged. Stage 1 uses 2,048 points; Stage 2 uses
8,192 points. The primary Actor-quality comparison applies the same plain
weighted tangent Feynman--Kac estimator with 512 paths and 32 steps to every
Actor. Ground Truth is used only after training and never for labels, stopping,
or checkpoint selection.

## Stage 0: correctness and label noise

The mandatory gate used eight independent replications on 128 states, common
state and Brownian banks, and two frozen checkpoints.

| Source Actor | paths | Value variance MCV/plain | Value reduction | Gradient variance MCV/plain | Gradient reduction | time MCV/plain | memory MCV/plain | max CRN-FD relative error |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `a3` | 128 | 0.07926 | 12.6x | 0.07455 | 13.4x | 1.030 | 1.019 | 0.00610 |
| late | 128 | 0.06226 | 16.1x | 0.06072 | 16.5x | 1.046 | 1.019 | 0.01003 |

The paired MCV-minus-plain Value and directional-Gradient confidence intervals
both contain zero for both checkpoints. Raw Gradient labels are correlated at
0.9993, while their replication noise correlation is only 0.17--0.20. This is
the expected pattern for two unbiased estimators of the same signal with very
different noise. ESS is 0.9804--0.9806 and all nonfinite counts are zero.

The four correctness tests cover predictable left endpoints, differentiation
of the complete correction, exact reduction to plain labels when \(q=0\),
parameter freezing with state-graph preservation, and common-random-number
finite differences. All four pass.

## Stage 1: does lower noise generalize?

The short 2,048-point comparison uses a separate 512-path target bank for the
common Actor diagnostic.

| Benchmark | Value RMSE plain | Value RMSE MCV | Gradient-vector RMSE plain | Gradient-vector RMSE MCV | plain-minus-MCV policy value | post-hoc PDE RMSE MCV/plain |
|---|---:|---:|---:|---:|---:|---:|
| nonseparable, \(d=100\) | 0.008830 | 0.008250 | 0.012120 | 0.010594 | +2.665e-5 | 0.9923 |
| joint concave, \(d=400\) | 0.013712 | 0.013974 | 0.028576 | 0.026845 | +5.154e-5 | 0.9871 |

Both reference-free policy differences are positive with MC confidence
intervals excluding zero. The \(d=400\) Value head is 1.9% worse in this short
run, but its Gradient head, policy value, and post-hoc PDE accuracy improve.
The predeclared Stage-1 gate therefore passed before Stage 2 was launched.

## Stage 2: Actor quality at high budget

The table reports means over two training seeds. B4 corresponds to the first
complete 8,192-point update (`a4`). B5 uses the validation-selected complete
checkpoint available below 960 seconds: `a5` for the nonseparable problem and
`a6` for the \(d=400\) problem.

| Benchmark | budget point | time plain / MCV (s) | Actor Value RMSE plain / MCV | reduction | Actor Gradient-vector RMSE plain / MCV | reduction |
|---|---|---:|---:|---:|---:|---:|
| nonseparable, \(d=100\) | B4 | 443.9 / 448.6 | 0.007266 / 0.005054 | 30.4% | 0.009630 / 0.006856 | 28.8% |
| nonseparable, \(d=100\) | B5 | 718.9 / 715.9 | 0.006244 / 0.003441 | 44.9% | 0.008680 / 0.004841 | 44.2% |
| joint concave, \(d=400\) | B4 | 348.4 / 350.9 | 0.014223 / 0.012217 | 14.1% | 0.019060 / 0.015650 | 17.9% |
| joint concave, \(d=400\) | B5 | 776.0 / 771.6 | 0.009185 / 0.006888 | 25.0% | 0.013690 / 0.008931 | 34.8% |

These validation errors use each branch's independently generated labels, so
Stage 1 supplies the stricter common-target generalization check. The Stage-2
trend is nevertheless consistent across every seed and iteration.

Reference-free cumulative policy improvement from `a3` to the B5 checkpoint
is 0.0001935 with MCV versus 0.0001681 with plain labels on the nonseparable
problem (+15.1%), and 0.0011532 versus 0.0010416 on the \(d=400\) problem
(+10.7%). Thus the gain is visible without Ground Truth.

## Matched-compute PDE accuracy

These are validation-selected operating points evaluated with the identical
plain 512-path, 32-step estimator. Plain candidates include the full
compute-reallocation envelope: extra iterations, 4,096/8,192 points, and
256/512 paths per label. At B4 and B5, the selected plain competitor is the
8,192-point branch.

| Benchmark | budget | Relative L2 plain | Relative L2 MCV | reduction | systematic RMSE plain | systematic RMSE MCV | reduction |
|---|---:|---:|---:|---:|---:|---:|---:|
| nonseparable, \(d=100\) | B4 | 0.006202 | 0.006181 | 0.34% | 0.002219 | 0.002196 | 1.06% |
| nonseparable, \(d=100\) | B5 | 0.006179 | 0.006151 | 0.46% | 0.002195 | 0.002165 | 1.36% |
| joint concave, \(d=400\) | B4 | 0.006521 | 0.006451 | 1.07% | 0.002128 | 0.002009 | 5.57% |
| joint concave, \(d=400\) | B5 | 0.006172 | 0.006124 | 0.77% | 0.001442 | 0.001285 | 10.85% |

At B5, both seeds improve separately. The \(d=400\) systematic RMSE changes
from 0.002178 to 0.002084 in seed 0 and from 0.000706 to 0.000486 in seed 1.
The corresponding nonseparable reductions are 0.002279 to 0.002247 and
0.002112 to 0.002084.

B3 is not a symmetric two-seed comparison: the staged protocol ran the
2,048-point MCV update only for seed 0, while the pre-existing plain envelope
contains cheap 2,048-point updates for both seeds. The B3 CSV is retained, but
the research conclusion is based on the predeclared 480--960 second region.

## Separating training and inference gains

The B5 Actors were evaluated in all four combinations. Values below are
two-seed means at 512 paths and 32 steps.

| Benchmark | training | inference | Relative L2 | systematic RMSE | MC fraction of observed MSE | inference seconds |
|---|---|---|---:|---:|---:|---:|
| nonseparable, \(d=100\) | plain | plain | 0.006179 | 0.002195 | 0.677 | 45.6 |
| nonseparable, \(d=100\) | MCV | plain | 0.006151 | 0.002165 | 0.683 | 46.9 |
| nonseparable, \(d=100\) | plain | MCV | 0.003565 | 0.002140 | 0.077 | 48.0 |
| nonseparable, \(d=100\) | MCV | MCV | 0.003512 | 0.002110 | 0.076 | 48.1 |
| joint concave, \(d=400\) | plain | plain | 0.006172 | 0.001442 | 0.863 | 33.5 |
| joint concave, \(d=400\) | MCV | plain | 0.006124 | 0.001285 | 0.880 | 37.1 |
| joint concave, \(d=400\) | plain | MCV | 0.002399 | 0.001406 | 0.253 | 38.7 |
| joint concave, \(d=400\) | MCV | MCV | 0.002273 | 0.001318 | 0.268 | 38.7 |

MCV inference lowers MC variance by 22.4--27.7x. More importantly for this
study, MCV training remains better under both inference estimators. Under MCV
inference, its Relative-L2 gain is 1.48% on the nonseparable problem and 5.27%
on \(d=400\); systematic-RMSE gains are 1.42% and 6.22%, respectively.

## Euler sensitivity

Both validation-selected training variants were evaluated with MCV inference
at 16, 32, and 64 Euler steps. Entries are two-seed means.

| Benchmark | steps | Relative L2 plain / MCV training | MCV-training reduction | systematic RMSE plain / MCV training | MCV-training reduction |
|---|---:|---:|---:|---:|---:|
| nonseparable, \(d=100\) | 16 | 0.006958 / 0.006925 | 0.47% | 0.004264 / 0.004245 | 0.44% |
| nonseparable, \(d=100\) | 32 | 0.003565 / 0.003512 | 1.48% | 0.002140 / 0.002110 | 1.42% |
| nonseparable, \(d=100\) | 64 | 0.001908 / 0.001853 | 2.89% | 0.001100 / 0.001071 | 2.62% |
| joint concave, \(d=400\) | 16 | 0.004174 / 0.004010 | 3.93% | 0.002614 / 0.002499 | 4.39% |
| joint concave, \(d=400\) | 32 | 0.002399 / 0.002273 | 5.27% | 0.001406 / 0.001318 | 6.22% |
| joint concave, \(d=400\) | 64 | 0.001532 / 0.001413 | 7.78% | 0.000817 / 0.000735 | 9.98% |

The absolute error roughly halves as the grid doubles, so 32 steps is not the
numerical accuracy limit. The training advantage persists at every grid and
grows at 64 steps. It therefore is not an artifact of the original time
discretization; the coarser grid partly hides the higher-quality Actor.

## Cost, stability, and remaining bottleneck

Across matched 8,192-point updates, the geometric-mean MCV/plain ratios are
0.997 for complete-update time and 0.994 for label time; these differences are
timing noise. Stage-0 isolated timing gives the conservative overhead estimate
of 3.0--4.6%. Peak label memory changes from 86.2 to 87.9 MB in \(d=100\) and
from 102.8 to 109.4 MB in \(d=400\). Training-label ESS ranges from 0.9799 to
0.9891, there are no nonfinites, and no Value labels require invariant-region
projection.

Label replication RMSE falls by about 3.4--4.3x, while Actor validation RMSE
falls by 14--45%. Consequently, the Actor error divided by label noise is often
higher under MCV. This is evidence that supervision noise has ceased to be the
only bottleneck; Actor capacity, optimization, and the Euler grid now limit how
much of the variance reduction becomes final PDE accuracy.

## Answers to the eight research questions

1. **Does MCV reduce Gradient-label variance?** Yes. The independently
   replicated component/vector variance falls by 13.4x at `a3` and 16.5x at
   the late Actor, with confidence intervals far below one.
2. **Does lower noise improve Actor generalization?** Yes. The independent
   Stage-1 target bank improves Gradient-vector error on both benchmarks, and
   the high-budget held-out Actor errors improve for every benchmark, seed,
   and selected iteration.
3. **Does MCV improve training accuracy at matched GPU time?** Yes. B4/B5
   complete-update times are matched, while Actor, reference-free policy, and
   post-hoc PDE metrics improve.
4. **Does it help in the 480--960 second region?** Yes. Both benchmarks and
   both seeds improve at B4 and B5.
5. **Is it useful with 8,192 points?** Yes; this is where the strongest Actor
   gains are observed.
6. **Does it beat reallocating plain compute?** Yes in the target B4/B5
   region. Validation selects the 8,192-point plain branch over extra paths,
   yet the matched MCV branch is still better.
7. **Does it lower the systematic error floor?** Yes, consistently but by
   different amounts: 1.36% on the nonseparable problem and 10.85% on
   \(d=400\) at B5 under common plain inference. Euler refinement shows that
   the reported 32-step floor still contains material discretization error.
8. **Is the gain large enough for the default algorithm?** Yes. The gain is
   consistent, reference-free policy improvement is 10.7--15.1% larger, and
   the runtime cost is negligible. The correction should be enabled by
   default, with plain labels retained as an ablation.

## Limits

The study covers two PDEs and two training seeds; the Stage-0 replication gate
uses one PDE and one seed. Per-seed systematic-RMSE confidence intervals can be
wide under plain inference because 68--88% of its observed MSE is MC variance;
the MCV-inference rows give much tighter diagnostics. The 32-step evaluator is
not Euler-converged. These facts limit claims about universal effect size, but
not the observed matched-protocol direction of improvement.

## Reproducibility and artifacts

The experiment ran on NVIDIA GH200 120 GB GPUs with PyTorch 2.10.0+cu129.
Jobs were `3392994` (Stage 0), `3393041` (Stage 1), `3393082` (Stage 2),
`3393224` (common evaluation), and `3393481` (symmetric Euler extension).
Their elapsed times were 3:06, 7:16, 58:28, 33:34, and 7:35,
respectively: 1.15 GH200-hours for the staged training study and 1.83
GH200-hours including all final diagnostics.
Configuration and selection rules are fixed in
[`protocol.json`](../../experiments/tangent_actor_mcv_training/protocol.json).

Primary artifacts:

- [`label_variance_comparison.csv`](label_variance_comparison.csv)
- [`actor_training_frontier.csv`](actor_training_frontier.csv)
- [`actor_validation_diagnostics.csv`](actor_validation_diagnostics.csv)
- [`budget_operating_points.csv`](budget_operating_points.csv)
- [`systematic_error_comparison.csv`](systematic_error_comparison.csv)
- [`euler_sensitivity.csv`](euler_sensitivity.csv)
- [`runtime_and_resource_accounting.csv`](runtime_and_resource_accounting.csv)
- [`reproducibility_manifest.json`](reproducibility_manifest.json)

Figures:

- [`relative_l2_vs_training_seconds.png`](figures/relative_l2_vs_training_seconds.png)
- [`systematic_rmse_vs_training_seconds.png`](figures/systematic_rmse_vs_training_seconds.png)
- [`label_variance_reduction.png`](figures/label_variance_reduction.png)
- [`actor_validation_by_iteration.png`](figures/actor_validation_by_iteration.png)
