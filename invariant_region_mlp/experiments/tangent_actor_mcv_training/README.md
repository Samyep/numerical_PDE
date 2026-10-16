# MCV training validation

This experiment tests whether the predictable martingale correction improves
the **training labels** and the resulting Tangent-Actor, separately from its
already established inference benefit.

The mandatory Stage 0 command is:

```bash
python -m invariant_region_mlp.experiments.tangent_actor_mcv_training.stage0
```

The source Actor's parameters are frozen, while the complete corrected payoff
`P-H` retains its computational graph through `C`, `q(t,X)`, and the controlled
state trajectory. Ground Truth is never used to create labels or select a
checkpoint.
