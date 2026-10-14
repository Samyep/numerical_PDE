# HJB compute-frontier preflight report

## Status

CPU audit passed before GPU submission. The authoritative machine-readable record is `results/hjb_compute_frontier/validation.json`.

The experiment runs on branch `ir-mlp-hjb-strong-baselines-compute-frontier`, starting from commit `56f6448e`. It creates independent experiment, results, documentation, cache, checkpoint, prediction, and Slurm-log locations; existing B1 artifacts are read-only.

## Equation checks

For B1 `(400,40)`, the exact analytic PDE residual was at most `2.22e-16`. Terminal Value and (z=\sqrt2\nabla u) agree exactly between the terminal and nonterminal-reference implementations, and the PyTorch terminal derivative agrees to `1.39e-16`.

For HJB-Log `(400)`, the Hopf-Cole integral is evaluated after the change of variable (y=400s), using 256-node float64 Gauss-Legendre quadrature on `[0,80]` and log-domain summation. On independent points at `t=0`, `t=0.53`, and `t=T`, comparison with adaptive infinite-interval `scipy.quad` gave maximum absolute errors:

- Value: `2.30e-13`
- (z): `1.53e-15`
- terminal derivative/autodiff consistency: `3.47e-18`

The (y=ds) rescaling prevents high-dimensional mass near (s=0) from being missed by generic adaptive integration.

## Data separation

Each PDE has a new 256-point policy-validation bank and a 2,048-point sealed test bank drawn from independent deterministic namespaces. The historical 960-point B1 bank is not read. The frozen SHA-256 identifiers are recorded in `validation.json`, including points and exact references.

## Sign and implementation checks

Automated tests cover:

- B1 PDE residual and HJB-Log quadrature against an independent integrator.
- terminal Value and terminal (z) in NumPy and PyTorch.
- controlled drift `-sqrt(2) q`, cost `+0.5 ||q||^2`, and martingale correction `-q^T dW` at the terminal boundary.
- freezing source-Actor parameters while retaining derivatives through the controlled state path.
- hard terminal constraints for PINN and DPI.
- Deep BSDE forward sign `dY=+0.5||Z||^2 dt+Z dW`.
- finite SDGD second-derivative objective and backward propagation.

Seven focused tests pass. Baseline source code contains no `exact_u` or `exact_z` call. Actor label generation also receives no nonterminal analytic reference.

## Frozen selection asymmetry

External baselines use only full PDE residual, terminal mismatch, or preset Picard round rules. Per the user's latest instruction, our declared Actor candidates may be ranked by final test RMSE. This supersedes the reference-free Actor-selection clause in the original task and is recorded as `test_rmse_oracle_allowed_by_user`. It must be described as oracle development selection in every final comparison.

## GPU preflight gate

The smoke suite completed all four methods on both PDEs. Twenty-two stored prediction artifacts reproduced every reported scalar exactly. A final-scale GH200 microbenchmark then measured projected per-run costs: 31.8 minutes for 10,000-step PINN+SDGD, 6.4 minutes for 10,000-step Deep BSDE, 2.0 minutes for 20-round DPI, and 5.7 minutes for one 8,192-point high-quality Actor policy-improvement phase. A 2,048-point 256-path/64-step Actor evaluation projects to 19 seconds. Peak allocated CUDA memory was below 0.5 GB.

These measurements place the planned two-PDE, three-seed comparison and targeted Actor exploration below the 18-hour allocation ceiling. Complete pilot accounting, rather than the linear projection, will determine the final launch grouping. The external configurations remain unchanged.
