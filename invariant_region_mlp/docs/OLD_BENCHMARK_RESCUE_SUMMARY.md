# Old benchmark high-budget rescue summary

## Executive result

The two older benchmarks do not join active VB as clean positive benchmarks, although both show that Samplewise certified IR can stabilize Raw full-history MLP.

| PDE | Verdict | Raw vs IR | Suppression floor crossed? | Paper role |
| --- | --- | --- | --- | --- |
| 100D Funding | F-B, partially rescued | IR consistently improves Raw and Raw approaches IR at high sampling | No; `z=0` remains much better | Secondary stabilization/mechanism result |
| 100--160D Rosenbrock HJB | H-B, trending but not reached | IR strongly stabilizes Raw; both improve along n=2 high sampling | No; `f=0` remains 60--160x better at M=96 | Weak-nonlinearity limitation/boundary case |

## Scientific recommendation

Active VB should remain the principal positive nonlinear benchmark. Funding may be retained as realistic evidence that certified projection regularizes Raw recursion, but not as proof that the correct gradient channel is useful at attainable budget. Rosenbrock HJB should be explicitly reframed as a mechanism/limitation example: correctness-preserving geometry can be statistically unhelpful when the true nonlinear correction is much smaller than gradient Monte Carlo noise.

The primary comparison was Raw versus Samplewise certified IR throughout. Suppression, fixed shrinkage, and an invalid tighter envelope were restricted to representative diagnostics. Batch-IR was omitted from the new study.

## Compute and reproducibility

The recorded parallel orchestration wall time is 25.031 minutes; summed per-task wall time is 1.629 CPU-task hours. The experiment contains 4410 Funding root repetitions and 213 HJB method-repetitions. All scientific arrays are float64. Legacy headline reproduction, paired terminal blocks, random-tree fingerprints, finite-state checks, work accounting, and manifest completeness are covered by the integrity audit.

Detailed reports:

- `docs/FUNDING_HIGH_BUDGET_RESCUE.md`
- `docs/HJB_HIGH_BUDGET_RESCUE.md`

Results and figures are under `results/rescue_high_budget/`. No manuscript source was modified.
