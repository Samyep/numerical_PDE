# Results: direct complete-flux baseline

## Verdict

The direct full-flux MLP is easy to optimize on the random validation
distribution, and it is **not** a zero-output solution.  However, it loses the
canonical-shock robustness and 64-to-512 grid transfer supplied by the
HLLC/Roe-dissipation parameterization.  It is therefore not the better HCFL
architecture for the current objective.

This is a seed-0 controlled result, not a multi-seed superiority claim.

## Convergence

The run used a 50,000-update cap but stopped automatically at update 23,000
after reaching the minimum learning rate and satisfying the validation
plateau rule.  The validation-selected checkpoint is update 20,000.

| model | parameters | selected update | stop update | validation rollout NRMSE |
|---|---:|---:|---:|---:|
| HLLC + direct vector correction | 6,627 | 8,300 | 10,000 | 0.0494303 |
| HLLC + Roe-dissipation correction | 6,627 | 12,500 | 16,500 | 0.0315568 |
| direct complete flux | 6,627 | 20,000 | 23,000 | **0.0280469** |

On the validation selection metric, full flux is 43.26% below direct
correction and 11.12% below Roe-dissipation.  This makes the held-out tests
essential: validation alone would select the wrong architecture for robust
shock transfer.

## 64-cell held-out tests

| split | direct correction | Roe dissipation | complete flux |
|---|---:|---:|---:|
| ordinary ID | 0.011437 | 0.009132 | **0.009032** |
| broad random in-support | 0.033078 | 0.023035 | **0.019202** |
| moderate OOD high-frequency | 0.006677 | 0.006798 | **0.005754** |
| Sod | 0.020117 | **0.014678** | 0.036873 |
| Lax | 0.065987 | **0.055243** | 0.070964 |
| collision | 0.116844 | **0.075152** | 0.129348 |
| strong pressure | 0.124339 | **0.088452** | 0.111991 |
| near-vacuum expansion | 0.112997 | **0.079350** | 0.138356 |

Across the three random-distribution tests, complete flux has mean NRMSE
0.011329, 12.78% below Roe-dissipation.  Across the five canonical Riemann
tests, it has mean NRMSE 0.097506, 55.82% above Roe-dissipation and loses all
five cases.

## The network did not collapse to zero

Diagnostics on independent validation ground-truth states give:

- raw learned flux RMS: 4.1786;
- HLLC flux RMS: 4.5641;
- projected learned flux-divergence RMS: 1.0043;
- HLLC flux-divergence RMS: 0.9281;
- learned/HLLC divergence RMS ratio: 1.0821.

A zero or spatially constant proposal would have near-zero flux divergence,
so the measured ratio rules out that failure mode.

The hard Tadmor projection changes 9.87% of validation interfaces, but its RMS
change is only 1.61% of the raw learned-flux RMS.  The local positivity limiter
is inactive in every reported 64-cell test.  Only the broad-random split uses
the fully-discrete entropy limiter, at a 0.0741% intervention rate.  Thus the
good random-distribution result is not primarily a low-order fallback result.

The model is not exactly classically consistent: on sampled constant states,
the normalized error in `F_NN(U,U) = F_physical(U)` is 0.1120.  This is an
expected identifiability weakness because periodic trajectory loss observes
flux differences, not the additive flux gauge.  Constant states are still
preserved to numerical precision because their spatial flux divergence is
approximately zero (maximum measured value `3.81e-6`).

## Zero-shot deployment on 512 cells

Errors below use native HLLC-2048 restricted conservatively to 512 cells only
for the metric.  Every plotted curve itself remains on its native grid.

| case | native HLLC-512 | Roe HCFL-512 | complete-flux HCFL-512 |
|---|---:|---:|---:|
| Sod | 0.010811 | **0.006820** | 0.026413 |
| Lax | 0.051048 | **0.030757** | 0.067519 |
| collision | 0.068022 | **0.051155** | 0.162742 |
| strong pressure | 0.171412 | **0.129729** | 0.242031 |
| near-vacuum expansion | 0.069891 | **0.062885** | **failed** |

For the four completed cases, complete flux has mean NRMSE 0.124676, which is
128.28% above Roe-dissipation and 65.52% above native HLLC-512.  Its hard
projection intervention rate rises to 43.9--51.5%, and the final profiles show
visible overshoots/oscillations near discontinuities.

The near-vacuum rollout reaches the explicit computational guard at saved
snapshot 26.  A nominal update was allowed up to 64 internal safe substeps
(the matched expected count is one).  Before failure:

- minimum pressure reached `1.0000001e-5`, essentially the safety floor;
- minimum density was `4.369e-4`;
- the internal step fell to `7.147e-10`;
- 1,188 low-order step halvings had occurred.

The maximum characteristic speed was only 2.94, so this is not a wave-speed
explosion.  It is a time-step collapse caused by approaching the admissibility
boundary while establishing the low-order positivity/entropy premise.

## Scientific interpretation

Directly learning all three flux components gives the MLP more freedom to fit
the random trajectory distribution, but removes three useful inductive biases:

1. exact physical consistency inherited from an analytic base flux;
2. an established shock-capturing Riemann flux outside the learned correction;
3. Roe-wave coordinates that localize learned changes to characteristic
   dissipation.

The controlled experiment therefore refines the conclusion rather than merely
declaring complete-flux learning “bad”: it wins the random validation objective
but fails the more demanding canonical and resolution-transfer criteria.  For
the current scientific claim, HLLC plus learned Roe-dissipation remains the
defensible method.
