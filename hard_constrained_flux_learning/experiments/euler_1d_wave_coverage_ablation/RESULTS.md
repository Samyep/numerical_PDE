# Seed-0 result: balanced wave coverage

## Conclusion

The pre-registered screen is **negative**. Replacing the 100 random-extreme
training trajectories with 100 balanced contact/compression/expansion/
pressure-jump/collision trajectories does not produce a sufficiently large or
consistent gain. Do not spend compute on multi-seed confirmation of this exact
formulation.

## What was inspected

Before the run, the repository status, plan, recent history, Euler notes,
committed CSVs, local ablation runner, PyClaw generator, HLLC proposal, local
admissibility limiter, and fully-discrete entropy ablations were reviewed.

The current 1D baseline is the five-point **direct-vector HLLC correction**
followed by hard Tadmor projection, local conservative admissibility limiting,
and the rare trajectory-wise fully-discrete entropy line search. Prior work had
already tested invariant inputs, characteristic corrections, learned
dissipation, Conv1d, HLLC versus Rusanov proposals, local versus global safety,
and four-step fine-tuning. The remaining empirical weakness was accuracy on
held-out canonical wave regimes.

Two evaluation issues were corrected for both arms:

1. the historical local runner's `OOD` generator is also present in broad
   training, so it is named `broad_random_in_support` here;
2. extreme float32 states exposed roundoff leaks after a one-pass Tadmor
   projection and after the double-precision entropy line search was converted
   back to float32. A float64 residual recheck/backoff and a reconstructed-
   update recheck close these numerical leaks symmetrically for both arms.

The committed historical stress CSVs do not include their evaluator. This run
therefore declares a reproducible 64-snapshot canonical horizon
(`t_final=0.0252`) and reruns the baseline rather than mixing the new values
with old numbers.

## Hypothesis and controlled change

**Hypothesis:** severe held-out errors are primarily a wave-regime coverage
problem rather than a capacity problem.

Only the last 100 of 580 training trajectories changed. The first 220 ordinary
and 260 broad-random trajectories were bit-for-bit shared. Both arms used the
same baseline-derived normalization, 6,627-parameter network initialization,
minibatch-index sequence, Adam settings, 1,100 updates, hard/safety layers,
reference solver, and evaluation data.

PyClaw was unavailable. Both arms used the existing 512-cell Rusanov + SSP-RK2
reference, with silent state repair disabled.

## NRMSE comparison

Negative change means the wave-coverage arm is better.

| evaluation | broad-random baseline | wave coverage | relative change |
|---|---:|---:|---:|
| ordinary ID | 0.012356 | **0.012196** | -1.29% |
| broad random (in support) | **0.034333** | 0.034355 | +0.06% |
| moderate OOD, frequencies 4--6 | **0.006759** | 0.006798 | +0.58% |
| Sod | **0.017851** | 0.021613 | +21.08% |
| Lax | **0.062785** | 0.062863 | +0.12% |
| collision | 0.101875 | **0.096482** | -5.29% |
| strong pressure | 0.119305 | **0.118283** | -0.86% |
| near-vacuum expansion | 0.089056 | **0.086370** | -3.02% |
| canonical mean | 0.078174 | **0.077122** | -1.35% |

The arm improves three of five canonical cases, but only collision clears 5%,
the mean gain is much smaller than the pre-specified 5%, and Sod regresses by
21.08%. Ordinary ID is not harmed; moderate OOD is effectively unchanged.

## Constraint and intervention metrics

Tadmor violation rate is zero for every row at the repository-standard `1e-5`
tolerance. Minimum pressure remains above the `1e-5` floor. The cases where
the safety layers matter are:

| case / arm | min density | min pressure | local limiter rate | FD entropy rate |
|---|---:|---:|---:|---:|
| collision / baseline | 0.044719 | 1.0054e-5 | 0.198% | 4.762% |
| collision / coverage | 0.057314 | 1.0064e-5 | 0.198% | 6.349% |
| near vacuum / baseline | 0.007578 | 1.0014e-5 | 1.587% | 23.810% |
| near vacuum / coverage | 0.010462 | 1.0014e-5 | 1.438% | 22.222% |

The largest measured Tadmor residual is `1.39e-8`, below the reporting
tolerance. The maximum fully-discrete entropy change is nonpositive within the
declared `1e-8` numerical tolerance for every rollout.

## Compute

Both models have 6,627 trainable parameters. CPU training took 5.32 seconds for
the baseline and 5.19 seconds for wave coverage; full evaluation took 1.51 and
1.55 seconds, respectively. The experiment does not change model compute.

## Interpretation and decision

This result rules out the simple claim that swapping a small random-extreme
slice for uniformly balanced wave templates is enough to solve the remaining
canonical gap. The collision and near-vacuum improvements are real within this
seed, but they are accompanied by a much larger Sod regression and no moderate
OOD gain. Hard feasibility and admissibility remain intact, but they do not
imply predictive accuracy; zero measured Tadmor violations also do not establish
entropy-solution uniqueness.

**Decision: stop this exact replacement strategy; do not run more seeds.** If
training coverage is revisited later, it should be modified rather than merely
confirmed—for example, test additive coverage or interface-level regime
balancing without removing random extremes. No such second experiment was run
here.

Machine-readable results are in [`results/`](results/).
