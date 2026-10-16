# Tangent-Actor martingale control-variate quick validation

This inference-only experiment freezes existing Terminal-Aware Tangent-Actor
checkpoints and augments their weighted linear Feynman--Kac rollout with

\[
H=\sum_m C_m q(t_m,X_m)^\top\Delta W_m.
\]

The accumulator uses the same Actor call, the left endpoint before the state
increment, and the unscaled Brownian increment. Since `q` already represents
`sqrt(2) grad u`, it is not scaled again. The original evaluator remains
unchanged and is executed on identical random streams as an implementation
cross-check and runtime baseline.

Calibration estimates one scalar coefficient from a separate random stream.
It centers payoff and martingale pair samples within each initial state and
uses antithetic pair averages as the independent units. Numerical Ground Truth
is opened only after all coefficients and predictions are frozen.

The staged launchers are:

- `slurm/pilot.slurm`: `d=100, seed=0`, checkpoints `a1` and `a6`.
- `slurm/expand.slurm`: the missing `a3`, the second seed, and both `d=400`
  seeds, conditional on a useful pilot.

Rebuild lightweight outputs after all runs with:

```bash
python -m invariant_region_mlp.experiments.tangent_actor_mcv_quick.aggregate
pytest invariant_region_mlp/experiments/tangent_actor_mcv_quick/test_mcv.py -q
```
