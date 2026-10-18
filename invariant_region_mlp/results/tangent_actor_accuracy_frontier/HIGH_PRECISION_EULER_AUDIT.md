# High-Precision Euler Audit

This is a follow-up precision experiment on the three already selected frozen `a4` Actors. No Actor was retrained or reselected. Fine Brownian increments were regenerated on an independent stream and aggregated to construct common random numbers across 64, 128, 256, and 512 Euler steps.

The numerical reference was refined to `n=256`, `dt=3.125e-05`. The complete previous-to-current discrepancy, `8.517e-06`, is propagated as reference uncertainty.

| Euler steps | seeds | Value RMSE mean ± sd | per-seed 95% CI envelope | Relative L2 | residual non-MC proxy | MC contribution | mean bias | MCV VR | ESS | inference s | peak GB | weak order |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 64 | 3 | 0.000500862 ± 3.7e-06 | [0.000449, 0.000553] | 1.006e-03 | 0.000491007 | 9.8499e-05 | +3.645e-04 | 51.3× | 0.9733 | 90.6 | 1.60 | 1.01 |
| 128 | 3 | 0.000258083 ± 2.5e-06 | [0.000224, 0.000294] | 5.184e-04 | 0.000247609 | 7.22766e-05 | +1.851e-04 | 95.2× | 0.9733 | 178.6 | 1.60 | 1.03 |
| 256 | 3 | 0.000138283 ± 2.1e-06 | [0.000112, 0.000164] | 2.778e-04 | 0.000126737 | 5.46547e-05 | +9.524e-05 | 166.6× | 0.9733 | 354.6 | 1.60 | — |
| 512 | 3 | 8.14161e-05 ± 2.3e-06 | [5.83e-05, 0.000104] | 1.636e-04 | 6.84441e-05 | 4.3254e-05 | +5.066e-05 | 266.1× | 0.9734 | 705.0 | 1.60 | — |

The reported raw RMSE includes finite-path Monte Carlo error. The residual proxy removes the estimated antithetic-pair sampling contribution and the reference discrepancy in quadrature; it is a diagnostic rather than an alternative headline metric. Confidence intervals use a state bootstrap with Gaussian perturbations calibrated to per-state MC variance and the reference proxy.

The observed order is reported only when two successive common-random-number refinement differences remain positive after noise correction. It is not forced to first order.
