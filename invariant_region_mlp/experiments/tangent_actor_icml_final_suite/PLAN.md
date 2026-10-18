# ICML final benchmark execution plan

The suite is frozen from base commit `d36ab9e79d0e6ffe49c43d60cdb28a064dc9ac09` on branch `ir-mlp-tangent-actor-icml-final-suite`. It contains 42 PDE–dimension configurations and 210 method cells. The machine-readable contract is `results/tangent_actor_icml_final_suite/FROZEN_PROTOCOL.json`.

Execution is staged to prevent a bad transfer from consuming the full allocation:

1. **Stage A — audit and preflight.** Verify pinned sources, equations, diffusion and `z` conventions, references, invariant intervals, deterministic bases, test-bank hashes, and all historical reuse candidates. Freeze method-transfer rules before opening new test errors.
2. **Stage B — 100D anchors.** Reuse only byte-audited exact matches. Run every missing supported 100D cell with the published external budget and the mature four-update Tangent-Actor budget. Keep Deep BSDE in the separate native single-point table.
3. **Stage C — dimension extremes.** Run 10D and 400D for smooth PDEs and 20D/1000D for the MLP counterexample. Diagnose nonfinite values, memory, gradient scaling, and generalization before filling intermediate dimensions.
4. **Stage D — remaining dimensions.** Use resumable one-GPU jobs. A failed official job remains in the manifest. Published baseline budgets are never shortened to fit a queue window.
5. **Stage E — replications and precision.** Complete missing seeds and the frozen 64/128/256/512 Euler audit. Refinement is a post-freeze numerical audit and cannot select an Actor.

The cluster permits one submitted or running job for this user, so jobs execute serially. Checkpoints are idempotent, and scheduler retries do not alter numerical updates. If quota or allocation ends first, the final report will quantify every unfinished cell rather than substitute a cheaper run.

Primary selection rules are fixed in advance. External methods use their published stopping budget. Tangent-Actor uses only held-out tangent-label error and paired reference-free policy values. Ground Truth is opened only after the Actor is frozen. Deep BSDE results are never inserted into the amortized whole-function table.
