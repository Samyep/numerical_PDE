# Baseline fidelity matrix

| Method | Pinned revision | Native anchor | Frozen algorithm/budget | Transfer rule | Prohibited additions |
|---|---|---|---|---|---|
| Deep BSDE | `e76fed80b1995daf4a0dd6d3e2cf64a0931dabda` | HJB-Log 100D; Allen–Cahn 100D | scalar free `Y0`, free `Z0`, 19 nonshared time networks, 20 BSDE steps, float64, official schedules and terminal loss | tensor shapes and PDE coefficients only; each initial state is a new solve | amortized `Y0(x)`, terminal-aware residuals, MCV |
| DPI | `e8616c2f531c9f95299b246a13c1180b8f02822a` | OU/GMM 100D | official PISGradNet, recursive targets, MC supervision, 40 Picard rounds, native data and optimizer budget | dimension/PDE tensors only; official terminal embedding remains because it is DPI's own feature | our tangent update, MCV, oracle coordinates |
| PINN-HTE | same DPI revision | OU/GMM 100D | official hard terminal embedding, architecture, stochastic derivative samples, residual, optimizer and epoch count | dimension/PDE tensors only | Tangent-Actor terminal residual, MCV, extra tuning |
| PINN + SDGD | `c597cc4fe9caac0b3ca118ea84453405f25ff833` | HJB-Log 100D | width 1024, three tanh layers, official coordinate subsampling, adversarial loop, collocation bank, Adam schedule, 10,000 epochs | dimension/PDE tensors and source-defined dimension scaling only | terminal-aware Actor features, MCV, reduced epoch count |
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
| `official_dpi_transfer.py --method dpi` | pinned `PISGradNet`, source data generator and objectives; 40 Picard rounds, 20,000 optimizer epochs per round, official MC label budgets and optimizer | registered equation class, tensor dimension, PDE coefficients/terminal data, and frozen dimension-specific OU/GMM tensors | source checkpoints plus per-round RNG snapshots; full results require the complete published loop |
| `official_dpi_transfer.py --method pinn_hte` | pinned PINN/PISGradNet hard-terminal model, residual objective, stochastic derivative estimator, architecture and optimizer budget | the same equation/tensor substitutions as DPI | source final checkpoint and RNG metadata |
| `deep_bsde_transfer.py` | scalar trainable `Y0`, free `Z0`, 19 nonshared subnetworks, 20 time intervals, float64, clipped terminal loss; 2,000 HJB-style or 4,000 Allen–Cahn-style updates and the published learning-rate schedule | state dimension, drift/diffusion, driver and terminal tensors; hidden width follows the source `d+10` rule | exact model, Adam, data-generator and RNG state at logging boundaries |
| `sdgd_transfer.py` | source hard terminal embedding, width 1024, three tanh hidden layers, 100 collocation points, 20 PGD steps of size 0.05, `min(10,d)` sampled trace coordinates, Adam 1e-3 with the source linear schedule, 10,000 epochs | state dimension and PDE coefficients/terminal function only | exact model, optimizer, scheduler, adversarial states and RNG state |

The DPI/PINN transfer preserves the official method's own terminal embedding;
the SDGD transfer preserves its source-native hard terminal embedding. These
are baseline features, not Tangent-Actor additions. None of the external
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
