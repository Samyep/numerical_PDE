# Method Fidelity Audit

This audit belongs to protocol `single_point_32_ice_v1`. It describes the
algorithmic identity of every formal method before the pilot and formal runs.

## Tangent-Actor + MCV

- Every PDE-dimension-seed cell is trained from a fresh deterministic `a0`;
  previous whole-function checkpoints are not formal inputs.
- P1 and P3 use the mature terminal-aware gradient Actor and its direct-alpha
  martingale control variate.
- P2 remains a terminal-aware **value Actor**. Its MCV does not add a learned
  gradient head: the predictable integrand is
  `sqrt(2) * grad_x(r_theta)` obtained by automatic differentiation of the
  frozen value Actor at each left endpoint. The stochastic integral is
  zero-mean and is subtracted on the same paths. This is an audited extension
  of the mature value-Actor runner, not an oracle feature.
- P4 uses the complete mature joint value-gradient Actor, including the
  tangent reaction coefficient, drift, and source integral.
- All Actors receive the full ambient `(t,x)` input. Exact nonterminal values,
  query truth, low-dimensional oracle projections, and test-point selection are
  absent from training and selection.

Classification: `native_project_method`, with the P2 autograd martingale
explicitly tagged `value_actor_mcv_extension`.

## Full-History Multilevel Picard

- The method is an untrained recursive Monte Carlo point estimator.
- Its value component uses the standard multilevel telescoping sum. When a
  gradient is consumed by the driver, the state also carries `z=sqrt(2) grad u`:
  the terminal contribution uses `sqrt(2) grad g`, and driver contributions use
  the Bismut weight `N/sqrt(delta)`.
- The time randomization is Beta(`alpha=1/2`, 1), with exact reciprocal-density
  weight `(T-t) * fraction^(1-alpha) / alpha`.
- High- and low-level recursive subtrees share the sampled evaluation
  time/state but use independent subsequent random streams.
- No value clipping or `delta=1e-7` Bismut truncation is used in the formal
  standard estimator. Nonfinite recursion is retained as a numerical failure.
- P4 uses this same coupled `(u,z)` recursion with its actual
  `0.5*tanh(u)*||z||^2` driver. The B6 quadratic double-estimator is not reused
  for P4.
- The standard nonlinear plug-in estimator is not described as unbiased.

Classification: `standard_fh_mlp`. A future B6 repaired estimator, if emitted,
must remain a separately named secondary variant.

## Deep BSDE

- Scalar trainable `Y0`, free `Z0`, 19 nonshared time-slice Z networks, batch
  normalization, 20 Euler intervals, float64, clipped terminal loss, and the
  source-family optimization schedules are unchanged.
- The sole structural generalization is that forward state row zero is the
  frozen requested `x_i` instead of a hard-coded zero. Each point and seed gets
  a newly initialized and independently optimized model.
- Checkpoints preserve the model, optimizer, completed update count, Python,
  NumPy, CPU/CUDA PyTorch RNG states, the training generator, and fixed
  validation Brownian increments.
- P1 and P2 are `native_source_family` after the 100D zero-state rollout anchor.
- P3 and P4 are `faithful_transfer`: only equation data, horizon, dimension, and
  fixed initial state change. Their main runs retain the non-Allen source
  schedule and `[0,1]` Y0 initialization. Any later stabilization change must
  be emitted as `adapted_transfer` and cannot overwrite the faithful result.
- Exact solution values are loaded only after the fixed training schedule ends.

## Comparison boundary

Tangent-Actor amortizes one global training run over 32 queries; FH-MLP and Deep
BSDE preserve their point-solver identities. Accuracy is compared on identical
frozen points. Compute accounting therefore reports both per-query inference
and full training/allocation cost rather than pretending the methods have the
same amortization model.
