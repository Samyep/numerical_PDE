# Baseline Fidelity Audit

DPI, PINN-HTE, and SDGD-PINN are **faithful method transfers** to the two new PDEs. They are not claimed as native experiments from the source papers.

- No terminal-aware Tangent-Actor feature or MCV estimator is added to a baseline.
- DPI/PINN-HTE use the pinned generic semilinear paper profile and its source budget.
- SDGD-PINN uses the source-native adversarial HJB profile for B6 and the source-native general semilinear SDE profile for Dense Fourier.
- Only the equation, horizon, ambient dimension, and evaluation bank are adapted.
- Exact solutions are opened for frozen post-training evaluation only.
- Smoke runs are non-result artifacts and are excluded from every table and figure.

Completed transferred cells: 0 / 30.

| Method | Source commit | Classification |
|---|---|---|
| DPI | `e8616c2f531c9f95299b246a13c1180b8f02822a` | faithful_transfer |
| PINN-HTE | `e8616c2f531c9f95299b246a13c1180b8f02822a` | faithful_transfer |
| SDGD-PINN | `c597cc4fe9caac0b3ca118ea84453405f25ff833` | faithful_transfer |
