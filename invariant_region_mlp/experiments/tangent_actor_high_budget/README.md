# Terminal-Aware Tangent-Actor high-budget scaling

This experiment preserves the original three-update TA_A records and branches
from each saved `a3` checkpoint. It isolates extra iterations, nested training
point counts, and MC paths per tangent label. Held-out tangent-label error and
paired tangent-policy value are the only selection signals. The independent 2D
reference is read only in posthoc evaluation.

The comparison also repairs a protocol leak in the preceding strong-baseline
pilot. `deep_bsde_original` uses raw amortized `Y0(x0)` and independent raw
time-slice `Z_n(x)` networks. `sdgd_pinn_original` uses an ordinary tanh PINN
with the authors' soft terminal value/gradient penalty. Neither competitor uses
the Terminal-Aware Tangent-Actor ansatz. The earlier terminal-aware competitor
records remain untouched and are labelled historical in the final report.

Run one workload locally with:

```bash
python -m invariant_region_mlp.experiments.tangent_actor_high_budget.run \
  --method actor_pilot --benchmark nonseparable_d100_k2 --seed 0
```

