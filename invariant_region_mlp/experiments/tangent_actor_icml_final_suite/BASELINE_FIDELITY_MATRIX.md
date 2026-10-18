# Baseline fidelity matrix

| Method | Pinned revision | Native anchor | Frozen algorithm/budget | Transfer rule | Prohibited additions |
|---|---|---|---|---|---|
| Deep BSDE | `e76fed80b1995daf4a0dd6d3e2cf64a0931dabda` | HJB-Log 100D; Allen–Cahn 100D | scalar free `Y0`, free `Z0`, 19 nonshared time networks, 20 BSDE steps, float64, official schedules and terminal loss | tensor shapes and PDE coefficients only; each initial state is a new solve | amortized `Y0(x)`, terminal-aware residuals, MCV |
| DPI | `e8616c2f531c9f95299b246a13c1180b8f02822a` | OU/GMM 100D | OU/GMM: PISGradNet, 40 Picard rounds; other PDEs: paper §5.2 4×128 ELU, 80 Picard rounds, 4096 states and paths, 16 epochs/round | select the published profile from PDE structure before evaluation; no PISGradNet outside OU/GMM | our tangent update, MCV, oracle coordinates |
| PINN-HTE | same DPI revision | OU/GMM 100D | OU/GMM: source PISGradNet, 20k epochs, HTE-16; other PDEs: paper §5.2 4×128 ELU, soft terminal weight 10, 30k epochs, HTE-16 | no hard-terminal/PIS specialization outside OU/GMM | Tangent-Actor terminal residual, MCV, extra tuning |
| PINN + SDGD | `c597cc4fe9caac0b3ca118ea84453405f25ff833` | HJB-Log 250D (paper anchor; suite comparison is a 100D transfer) | HJB: 3×1024 tanh, hard terminal, PGD-20, 100 training/10 adversarial coordinates, 10k epochs; generic: 3×1024 tanh, 100×20 SDE points, soft terminal losses, 10 coordinates, 10k epochs | HJB profile only for HJB-structured PDEs; paper §5.2 profile elsewhere | Tangent-Actor terminal features, MCV, reduced epoch count |
| Tangent-Actor | base `d36ab9e7…` | method under study | 8192 states, 128 paths/label, 32 label steps, 4x384 SiLU, 5000 optimizer steps/update, four updates; frozen reference-free selection | only mathematically required driver/value-head changes | oracle reduced coordinates or analytic nonterminal labels |

Status labels in `DIMENSION_COMPATIBILITY_MATRIX.csv` have the following meaning:

- `native_official`: original PDE, dimension and source configuration.
- `faithful_transfer`: unchanged algorithm and per-run budget with only PDE/dimension tensors changed.
- `adapted`: a disclosed mathematical or technical change beyond those tensors.
- `unsupported`: no scientifically defensible execution path; the reason must be recorded.

Deep BSDE is a native single-point solver. Its rows belong only in `NATIVE_POINT_RESULTS.csv`; its full training cost is charged for each independently solved state. Earlier amortized Deep BSDE pilots remain useful diagnostics but are not relabeled as the official baseline.

The audited PyTorch Deep BSDE compatibility port and exact-resume SDGD transcription may be reused. Framework and checkpointing fixes must preserve every numerical update. Earlier compute-matched strong-baseline pilots used adapted architectures/budgets and therefore cannot populate the final official baseline cells.

## Implemented final-suite runners

The transfer runners are deliberately thin wrappers around the pinned source
algorithms. Smoke mode exists only to exercise equations, shapes, checkpoint
I/O, and evaluation; every smoke artifact is marked `smoke_nonresult` and is
excluded from final aggregation.

| Runner | Actual final configuration | Allowed transfer changes | Resume/status |
|---|---|---|---|
| `official_dpi_transfer.py --method dpi` | `hjb_pis` on OU/GMM; `semilinear_paper` elsewhere. The latter fixes paper Table-1 `K=80`, `M=N=4096`, `E=16`, 4×128 ELU, and source gradient supervision | registered equation class, horizon, diffusion, tensor dimension, coefficients and terminal data | source checkpoints plus per-round RNG snapshots |
| `official_dpi_transfer.py --method pinn_hte` | `hjb_pis` on OU/GMM; `semilinear_paper` elsewhere, with 4×128 ELU, soft terminal weight 10, HTE-16 and 30,000 epochs | the same equation/tensor substitutions as DPI | source final checkpoint and RNG metadata |
| `deep_bsde_transfer.py` | scalar trainable `Y0`, free `Z0`, 19 nonshared subnetworks, 20 time intervals, float64, clipped terminal loss; 2,000 HJB-style or 4,000 Allen–Cahn-style updates and the published learning-rate schedule | state dimension, drift/diffusion, driver and terminal tensors; hidden width follows the source `d+10` rule | exact model, Adam, data-generator and RNG state at logging boundaries |
| `sdgd_transfer.py` | `hjb_adversarial`: hard terminal, 100 collocation points, PGD-20, 100 training and 10 PGD coordinates, linear schedule; `semilinear_sde`: 100 trajectories ×20 steps, soft terminal value weight 20 plus gradient weight 1, 10 coordinates, exponential 0.9995 schedule. Both use 3×1024 tanh and 10,000 epochs | state dimension and PDE coefficients/terminal function only; two fresh coordinate subsets are sampled independently at every trajectory time slice as in the source | exact model, optimizer, scheduler, collocation state where applicable, and RNG state |

For high-dimensional DPI transfers, the pinned upstream implementation's own
`GPUMemoryTracker` may split one data-generation call into smaller resident
chunks. This does not change `DATA_SIZE=4096`, the 4,096 terminal/integral
paths per target, the configured 40 OU/GMM or 80 semilinear Picard rounds, or
the 16 optimizer epochs per round. It is source-native memory scheduling
rather than a reduced-budget adaptation.

The pinned DPI source calls `torch.load` without a `weights_only` argument.
PyTorch 2.6 changed that default and rejects the source's own Lightning
checkpoint objects. The final runner sets
`TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1`, restoring the pre-2.6 trusted-checkpoint
behavior without editing the pinned checkout, network state, optimizer, data,
or any numerical update. This is an environment compatibility fix, not an
algorithmic adaptation.

The DPI/PINN OU/GMM reproduction and SDGD HJB profile preserve their authors'
own source-native hard terminal embeddings. General semilinear transfers use
the authors' soft-terminal profiles. None of the external
runners uses the Tangent-Actor residual, learned tangent coefficients, MCV,
oracle coordinates, rank projection, or analytic nonterminal labels.

## Preflight evidence

CPU smoke executions completed for all distinct mathematical paths:

- DPI: joint concave, Allen–Cahn, the nonsmooth MLP-failure PDE, and OU/GMM.
- PINN-HTE: the nonconvex coupled PDE.
- Deep BSDE: HJB-Log, OU/GMM, joint convex, MLP failure, and Allen–Cahn.
- SDGD-PINN: HJB-Log, joint concave, MLP failure, and Allen–Cahn.
- Mature Tangent-Actor runners: standard gradient-Actor, value-only
  Allen–Cahn, joint-driver, convex sign-transform, nonconvex, and nonsmooth
  minimum-norm-subgradient paths.

The smoke checks found and corrected two transfer-only defects before final
execution: PINN-HTE requires a finite evaluation frequency, and the
MLP-failure terminal must accept the unbatched samples used by `functorch`.
They also verified that Allen–Cahn final evaluation uses the independently
refined radial reference. The tracked `RUNNER_PREFLIGHT_AUDIT.json` records
the coverage and source hashes. Smoke errors are never numerical results.

## 2026-10-18 fidelity correction

An audit against the paper text found that the earlier suite transcription
used 10 trace coordinates for both HJB training and PGD. The paper and pinned
HJB source prescribe 100 for training and 10 for PGD. It also applied the
specialized DPI/PINN HJB profile and SDGD adversarial profile to every PDE.
Those completed artifacts are preserved as nonprimary diagnostics. Corrected
profile-tagged runs are the only external-baseline rows eligible for the main
table. New training uses seed 0 only; historical extra seeds remain in
`HISTORICAL_AND_NONPRIMARY_RESULTS.csv` and are never selected by accuracy.
