# Tangent-Actor whole-function benchmark

**Integrity audit:** PASS. **Execution:** in_progress. Generated 2026-10-21T22:53:28.898784+00:00.

This benchmark contains the predeclared 7 PDE x 5 dimension x 4 method grid. Every training seed is 0. Deep BSDE is excluded because it is a separate single-point task.

## Execution coverage

- COMPLETE: 113
- OUT_OF_TIME: 0
- OUT_OF_MEMORY: 0
- FAILED: 0
- PENDING: 27
- New measured GPU time: 67.371 GH200 hours (120-hour cap).
- Slurm account charge: 134.790 credit-hours.

Incomplete cells have blank accuracy fields and are absent from all numerical curves.
Cell time is summed from completed attempts, including checkpoint resumes and failed attempts. Suite allocation time also includes reserved idle capacity, so the two totals need not agree.

## 100D repair and strict pairing

The apparent 960-versus-512 Nonconvex coupled mismatch was traced to stale metadata and confusion between file hashes and hashes of the actual arrays. Direct inspection showed that all four frozen predictions use the same 512-state `t` and `x` arrays and the same truth array. Their metrics were recomputed into canonical prediction artifacts without retraining or checkpoint reselection. Allen--Cahn was likewise recomputed against one common highest-accuracy frozen reference; the metric changes were below the recorded reference uncertainty and did not change the ranking.

The independent repair review and the content-hash audit both pass all seven 100D PDEs. The audit opens every prediction artifact, hashes the actual `t`, `x`, and truth arrays, and recomputes the reported errors; matching CSV metadata alone is insufficient.

## Strictly complete rankings

All four methods are complete in 18 PDE/dimension settings. Tangent-Actor has the smallest Value RMSE in 18/18 of them.

- Allen-Cahn, d=20: tangent_actor (0.001026) < dpi (0.00145) < pinn_hte (0.001718) < sdgd_pinn (0.01187)
- Allen-Cahn, d=100: tangent_actor (0.0001734) < pinn_hte (0.000948) < dpi (0.001067) < sdgd_pinn (0.002692)
- Allen-Cahn, d=400: tangent_actor (2.153e-05) < pinn_hte (0.0004262) < dpi (0.0004546) < sdgd_pinn (0.000616)
- HJB-Log, d=20: tangent_actor (0.0009483) < sdgd_pinn (0.006094) < pinn_hte (0.0131) < dpi (0.01486)
- HJB-Log, d=50: tangent_actor (0.000568) < dpi (0.03614) < pinn_hte (0.04226) < sdgd_pinn (0.04659)
- HJB-Log, d=100: tangent_actor (0.000814) < sdgd_pinn (0.03988) < pinn_hte (0.05362) < dpi (0.1339)
- HJB-Log, d=400: tangent_actor (0.0009245) < sdgd_pinn (0.02925) < dpi (0.3496) < pinn_hte (0.4749)
- Jointly concave value-gradient PDE, d=20: tangent_actor (0.0003143) < dpi (0.0008811) < pinn_hte (0.003834) < sdgd_pinn (0.02592)
- Jointly concave value-gradient PDE, d=50: tangent_actor (0.0003246) < dpi (0.0009667) < pinn_hte (0.002799) < sdgd_pinn (0.01795)
- Jointly concave value-gradient PDE, d=100: tangent_actor (0.0003231) < dpi (0.002142) < pinn_hte (0.003226) < sdgd_pinn (0.02333)
- Jointly convex value-gradient PDE, d=20: tangent_actor (0.0002138) < dpi (0.0007283) < pinn_hte (0.001622) < sdgd_pinn (0.01326)
- Jointly convex value-gradient PDE, d=100: tangent_actor (0.0002116) < dpi (0.001857) < pinn_hte (0.003267) < sdgd_pinn (0.01278)
- Hutzenthaler-Nguyen MLP failure example, d=20: tangent_actor (0.004033) < dpi (0.005611) < pinn_hte (0.006659) < sdgd_pinn (0.03903)
- Hutzenthaler-Nguyen MLP failure example, d=100: tangent_actor (0.006662) < dpi (0.01002) < pinn_hte (0.03994) < sdgd_pinn (0.06687)
- Hutzenthaler-Nguyen MLP failure example, d=400: tangent_actor (0.009286) < pinn_hte (0.0101) < dpi (0.01823) < sdgd_pinn (0.244)
- Nonconvex-nonconcave coupled PDE, d=100: tangent_actor (0.0008041) < dpi (0.001699) < pinn_hte (0.003068) < sdgd_pinn (0.008877)
- OU/GMM HJB, d=20: tangent_actor (0.06178) < pinn_hte (0.6489) < dpi (1.572) < sdgd_pinn (9.932)
- OU/GMM HJB, d=100: tangent_actor (0.2198) < dpi (17.66) < pinn_hte (38.67) < sdgd_pinn (62.48)

## Measured compute by method

The total below is measured cell GPU time. Training, selection, and inference components come from completed canonical rows and may omit setup or wrapper overhead.

| Method | Complete cells | Cell GPU h | Training h | Selection h | Inference h |
|---|---:|---:|---:|---:|---:|
| Tangent-Actor + MCV | 23 | 14.006 | 9.947 | 0.389 | 0.870 |
| DPI | 28 | 6.360 | 6.172 | 0.000 | 0.000 |
| PINN-HTE | 34 | 5.501 | 5.310 | 0.000 | 0.000 |
| SDGD-PINN | 28 | 59.830 | 58.040 | 0.000 | 0.000 |

## Stability and capped runs

All 113 completed rows have finite primary metrics, and the audit found no reported nonfinite outputs. Checkpoint resumes are charged to the same six-hour cell budget. A cell is marked OUT_OF_TIME only after its cumulative attempts reach that cap; capped checkpoints never enter the complete ranking.

No cell has yet ended as OUT_OF_TIME, OUT_OF_MEMORY, or FAILED.

## Recommended paper placement

The ICML main text can use the strict complete-setting ranking, the dimension-scaling figure, and the accuracy--compute figure, with the single-seed and transfer wording stated in their captions. The full 140-cell table, per-cell cost components, reference uncertainty, source commits, checkpoint-resume ledger, capped statuses, and 100D repair provenance belong in the appendix and reproducibility package.

## Reviewer-facing limitations

- Every new training result uses one predeclared seed (seed 0), so the suite does not estimate training-seed variability.
- DPI, PINN-HTE, and SDGD-PINN are source-faithful transfers to this common whole-function task. These rows do not establish universal superiority over results on the baselines' original paper-specific tasks.
- Numerical references and finite-path Monte Carlo evaluations have nonzero uncertainty; the corresponding fields remain in the canonical CSV.
- Known low-dimensional structure is used only by independent reference solvers. The trained methods receive full ambient coordinates; this suite therefore tests ambient scaling for these PDE families, while broader high-effective-dimension evidence remains a separate claim.
- Incomplete or capped cells cannot be counted as accuracy wins. Any remaining PENDING, OUT_OF_TIME, OUT_OF_MEMORY, or FAILED entries must remain visible in tables and plots.

## Interpretation limits

The baseline rows are source-faithful method transfers unless explicitly labeled native. The table therefore supports comparison on this frozen whole-function suite, rather than a claim of universal superiority over each method's original paper task. Results are single-seed training endpoints; numerical-reference and finite-MC uncertainty remain in the CSV.
