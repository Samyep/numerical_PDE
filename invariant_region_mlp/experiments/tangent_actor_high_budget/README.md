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

The staged GPU workflow is recorded in `slurm/actor_pilot.slurm`,
`slurm/stage2_and_baselines.slurm`, and `slurm/evaluate.slurm`.  The second
script continues only the 8192-point branch selected by the frozen stage-one
rule, then trains the corrected original-form Deep BSDE and SDGD-PINN
baselines.  The evaluation script uses a common reference-free validation bank
to select an Actor before opening the independent numerical reference.

After all four evaluation tasks finish, rebuild every lightweight artifact with:

```bash
python -m invariant_region_mlp.experiments.tangent_actor_high_budget.aggregate
python -m invariant_region_mlp.experiments.tangent_actor_high_budget.make_report
python -m pytest \
  invariant_region_mlp/experiments/tangent_actor_high_budget/test_high_budget.py -q
```

The results directory contains the execution manifest, protocol and
reproducibility metadata, scalar CSV files, independent metric replay and
protocol audits, figures in PNG/PDF form, and the final Markdown report.  Raw
predictions and model checkpoints are intentionally excluded from Git because
they are large; their paths and scheduler accounting remain in the committed
metadata.

`slurm/sdgd_recovery.slurm` is a narrowly scoped recovery launcher.  It runs
only the four SDGD tasks if a completed Deep-BSDE run prevents the chained
command from reaching SDGD.  It never retrains or tunes Deep BSDE; the
`recover_deep_bsde` utility records already-written timed checkpoints and the
observed instability with unavailable diagnostics represented as null.

