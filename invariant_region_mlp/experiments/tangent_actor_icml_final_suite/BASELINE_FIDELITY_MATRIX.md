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
