# Tangent-Actor whole-function benchmark

**Integrity audit:** PASS. **Execution:** in_progress. Generated 2026-10-19T21:04:15.284785+00:00.

This benchmark contains the predeclared 7 PDE x 5 dimension x 4 method grid. Every training seed is 0. Deep BSDE is excluded because it is a separate single-point task.

## Execution coverage

- COMPLETE: 28
- OUT_OF_TIME: 0
- OUT_OF_MEMORY: 0
- FAILED: 0
- PENDING: 112
- New measured GPU time: 0.000 GH200 hours (120-hour cap).

Incomplete cells have blank accuracy fields and are absent from all numerical curves.

## Strictly complete rankings

All four methods are complete in 7 PDE/dimension settings. Tangent-Actor has the smallest Value RMSE in 7/7 of them.

- Allen-Cahn, d=100: tangent_actor (0.0001734) < pinn_hte (0.000948) < dpi (0.001067) < sdgd_pinn (0.002692)
- HJB-Log, d=100: tangent_actor (0.000814) < sdgd_pinn (0.03988) < pinn_hte (0.05362) < dpi (0.1339)
- Jointly concave value-gradient PDE, d=100: tangent_actor (0.0003231) < dpi (0.002142) < pinn_hte (0.003226) < sdgd_pinn (0.02333)
- Jointly convex value-gradient PDE, d=100: tangent_actor (0.0002116) < dpi (0.001857) < pinn_hte (0.003267) < sdgd_pinn (0.01278)
- Hutzenthaler-Nguyen MLP failure example, d=100: tangent_actor (0.006662) < dpi (0.01002) < pinn_hte (0.03994) < sdgd_pinn (0.06687)
- Nonconvex-nonconcave coupled PDE, d=100: tangent_actor (0.0008041) < dpi (0.001699) < pinn_hte (0.003068) < sdgd_pinn (0.008877)
- OU/GMM HJB, d=100: tangent_actor (0.2198) < dpi (17.66) < pinn_hte (38.67) < sdgd_pinn (62.48)

## Interpretation limits

The baseline rows are source-faithful method transfers unless explicitly labeled native. The table therefore supports comparison on this frozen whole-function suite, rather than a claim of universal superiority over each method's original paper task. Results are single-seed training endpoints; numerical-reference and finite-MC uncertainty remain in the CSV.
