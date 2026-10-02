# Seed-0 result

## Outcome

The unconstrained DNN flux reached its best fully admissible validation
checkpoint at update 6,500 and satisfied the common plateau rule at update
10,300. Its validation rollout NRMSE was `0.04933`, only `0.19%` below the
matched constrained direct model's `0.04943`. The parameter count was identical
(`6,627`) and training took 89.1 seconds on CPU.

Thus, removing the hard output safeguards produced no material accuracy gain
on the ordinary validation distribution. It did remove reliability: raw
rollouts developed negative pressure on two of the five canonical tests.

## Controlled comparison

| split | constrained direct NRMSE | unconstrained NRMSE | change | raw completion |
|---|---:|---:|---:|---:|
| ordinary ID | 0.01144 | 0.01140 | -0.31% | 90/90 |
| broad in-support | 0.03308 | 0.03304 | -0.11% | 90/90 |
| moderate OOD | 0.00668 | 0.00668 | -0.02% | 90/90 |
| Sod | 0.02012 | 0.01987 | -1.22% | 1/1 |
| Lax | 0.06599 | 0.06620 | +0.32% | 1/1 |
| collision | 0.11684 | **failed** | n/a | 0/1 |
| strong pressure | 0.12434 | 0.12573 | +1.12% | 1/1 |
| near-vacuum expansion | 0.11300 | **failed** | n/a | 0/1 |

Negative values in the change column favor the unconstrained model. All six
completed comparisons differ by less than 1.3%, so there is no meaningful
accuracy advantage in either direction.

## Failure and entropy evidence

- Collision first became inadmissible at saved snapshot 32, with minimum
  pressure `-0.01032`.
- Near-vacuum expansion first became inadmissible at snapshot 14, with minimum
  pressure `-0.00922`.
- Raw Tadmor violation rates were `0.32%` on ordinary ID, `0.68%` on broad
  in-support data, `4.09%` on Sod, `2.10%` before collision failed, `3.40%` on
  strong pressure, and `1.56%` before near-vacuum failed.
- The constrained direct comparator completed every test with zero measured
  post-projection Tadmor violations.

The 1,100-update checkpoint already failed collision and near-vacuum, so
training longer did not repair the severe-rollout reliability problem. During
convergence, many later checkpoints also lost one of 136 validation
trajectories even while their one-step and partial-rollout losses decreased.

## Decision

Stop this direction as a candidate production solver and retain it as the
plain-DNN negative control. The seed-0 evidence is already decisive for the
structural question: the hard layers have negligible cost on completed
rollouts but prevent catastrophic physical-state failure. Multi-seed expansion
is not justified unless the question changes to quantifying failure
probability rather than selecting a solver.

See `results/unconstrained_flux_baseline_seed0.png`,
`results/metrics_seed0.csv`, and
`results/comparison_vs_constrained_direct_seed0.csv`.
