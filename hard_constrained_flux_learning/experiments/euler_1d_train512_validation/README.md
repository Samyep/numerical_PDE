# Euler direct-512 HCFL validation

This experiment trains the retained symmetric four-cell
`Central + nonnegative Roe + proposal-feasibility` model directly on 512-cell
trajectories. Training labels are native periodic HLLC-2048 trajectories
conservatively restricted to 512 finite-volume cell averages.

The seed-0 data split contains 580 training trajectories and 136 independent
validation trajectories. The training split contains 220 ordinary, 260 broad,
and 100 extreme trajectories. The validation split contains 44 ordinary, 52
broad, 20 extreme, and 20 wave-coverage trajectories. Checkpoint selection
uses autoregressive validation rollout NRMSE through the complete deployment
safety stack; the named deployment cases are not used for training or
selection. The update cap is 50,000, and a checkpoint is retained only after
the validation-plateau rule is satisfied at the minimum learning rate.

Run the differentiable-forward smoke test:

```powershell
python run_train512.py --self-test --device cuda
```

Generate data, train to convergence, evaluate, and make the combined figure:

```powershell
python run_train512.py --phase all --seed 0 --device cuda --resume
```

The comparison uses only three plotted curves: native HLLC-2048, a zero-NN
entropy-fixed Roe-512 control, and HCFL trained directly at 512. Both 512-cell
methods use the same hard entropy projection, local positivity limiter, and
fully-discrete entropy limiter during deployment.

The evaluation also runs the trained checkpoint with both `F_low` blends
removed. This projection-only ablation records admissibility failures,
fully-discrete entropy increases, and its state difference from the complete
hard-safe solver. It is not treated as hard-positive merely because a finite
test set happens to remain admissible.

Seed-0 results and the interpretation of the `F_low` ablation are in
[`RESULTS.md`](RESULTS.md). The deterministic 70 MB trajectory cache is local
and git-ignored; the converged checkpoint, metrics, audit metadata, and figures
are tracked.
