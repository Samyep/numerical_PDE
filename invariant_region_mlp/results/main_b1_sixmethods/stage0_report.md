# B1 Stage 0 correctness report

**Verdict:** reference construction PASS; MLP draw replay deferred (reference-only invocation).

The pre-method analytic calibration selected coefficient scale **2**. The rule and all rejected candidate metrics are saved before solver evaluation.

## Checks

- Stable float64 terminal and exact-solution log-sum-exp.
- Analytic PDE residual and terminal identity.
- Closed-form gradient checked independently by PyTorch autodiff in the unit suite.
- Dense coefficient affine rank and active Hessian rank.
- `k=1` exactly reproduces the established P1 direction and coefficients.
- Diffusion and gradient-state convention: `sigma=sqrt(2)`, `z=sqrt(2) grad u`.
- CPU/Torch recursion on identical draws: deferred (reference-only invocation).

## Calibration

Candidate rows: 20. Frozen-grid rows: 15.
Configurations flagged as weak by the preregistered diagnostic thresholds: none.

Artifacts: `scale_calibration_candidates.csv`, `pde_calibration.csv`, `frozen_pde_config.json`, and `stage0_mlp_tape_checks.csv` when run with a GPU.
