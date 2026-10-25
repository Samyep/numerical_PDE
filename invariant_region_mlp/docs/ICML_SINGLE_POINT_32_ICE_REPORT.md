# ICML Single-Point 32 ICE Report

Revision: single_point_32_ice_v2_ta_seed0_bundled (source protocol: single_point_32_ice_v1)

This report compares fresh Tangent-Actor plus MCV training, standard full-history MLP, and independent single-initial-state Deep BSDE on the same 640 frozen time-zero queries.
Tangent-Actor uses one fresh training seed per PDE-dimension cell and four independent MCV inference replicates, as documented in the 2026-10-23 protocol amendment.

## Main RMSE over points and formal repeats

| PDE | d | Tangent-Actor | FH-MLP | Deep BSDE |
|---|---:|---:|---:|---:|
| hjb_log | 20 | 2.000e-03 | 6.669e-03 | 2.553e-02 |
| hjb_log | 50 | 6.601e-03 | 3.461e-03 | 1.535e-02 |
| hjb_log | 100 | 5.592e-03 | 2.319e-03 | 7.952e-03 |
| hjb_log | 200 | 5.892e-03 | 1.502e-03 | 3.739e-03 |
| hjb_log | 400 | 3.576e-03 | 1.062e-03 | 1.450e-03 |
| allen_cahn | 20 | 1.475e-04 | 1.036e-03 | 3.123e-04 |
| allen_cahn | 50 | 7.148e-05 | 4.541e-04 | 2.216e-04 |
| allen_cahn | 100 | 5.228e-05 | 2.147e-04 | 1.777e-04 |
| allen_cahn | 200 | 4.967e-05 | 1.164e-04 | 1.439e-04 |
| allen_cahn | 400 | 8.938e-06 | 5.287e-05 | 1.089e-04 |
| b6_fullrank | 20 | 7.126e-04 | 8.560e-02 | 1.128e-01 |
| b6_fullrank | 50 | 1.419e-03 | 2.714e-01 | 2.526e-01 |
| b6_fullrank | 100 | 1.591e-03 | 5.893e-01 | 2.816e-01 |
| b6_fullrank | 200 | 4.168e-03 | 1.925e+00 | 3.722e-01 |
| b6_fullrank | 400 | 3.630e-03 | 3.190e+00 | 3.242e-01 |
| dense_fourier | 20 | 4.253e-03 | 3.598e-02 | 4.376e-02 |
| dense_fourier | 50 | 4.320e-03 | 4.226e-02 | 7.790e-02 |
| dense_fourier | 100 | 3.717e-03 | 1.331e-01 | 1.130e-01 |
| dense_fourier | 200 | 8.688e-03 | 1.427e-01 | 7.519e-02 |
| dense_fourier | 400 | 1.399e-02 | 2.933e-01 | 7.908e-02 |

## Relative MAE for the two newly added PDEs

For the 40 frozen Tangent-Actor inference runs on B6 Full-Rank and Dense
Fourier (2 PDEs x 5 dimensions x 4 independent inference replicates), we
additionally report
`rMAE = sum |prediction - reference| / sum |reference|`. Each run uses the
same 32 frozen queries as its RMSE. The table below pools the four formal
replicates in each PDE-dimension cell; all 40 per-run values are retained in
`results/single_point_32_ice/b6_dense_ta_40_run_metrics.csv`.

| PDE | d | RMSE | rMAE |
|---|---:|---:|---:|
| b6_fullrank | 20 | 7.126e-04 | 7.844e-04 |
| b6_fullrank | 50 | 1.419e-03 | 8.759e-04 |
| b6_fullrank | 100 | 1.591e-03 | 9.041e-04 |
| b6_fullrank | 200 | 4.168e-03 | 1.331e-03 |
| b6_fullrank | 400 | 3.630e-03 | 1.254e-03 |
| dense_fourier | 20 | 4.253e-03 | 6.063e-03 |
| dense_fourier | 50 | 4.320e-03 | 7.251e-03 |
| dense_fourier | 100 | 3.717e-03 | 2.333e-03 |
| dense_fourier | 200 | 8.688e-03 | 1.979e-02 |
| dense_fourier | 400 | 1.399e-02 | 3.147e-02 |

The complete pooled table is
`results/single_point_32_ice/b6_dense_ta_aggregate_metrics.csv`; the
derivation audit is
`results/single_point_32_ice/b6_dense_ta_rmae_audit.json`.

## Required conclusions

1. The complete dimension trends are given in the table above and in figure1_accuracy_vs_dimension.pdf. A failed numerical outcome is never ranked using a finite-only subset.

2. Tangent-Actor has the lowest finite RMSE in 16 of 20 complete PDE-dimension cells. Exact ties, if any, are counted for every tied method.

3. No formal FH-MLP configuration produced a nonfinite prediction. The full six-configuration data remain in fh_mlp_results.csv; the preregistered main comparison uses n=3, M=16.

4. Deep BSDE stability is measured from three independent trainings at every initial point. Its median and maximum cell-level repeat-RMSE standard deviations are 4.780e-04 and 7.419e-03.

5. B6 and Dense Fourier retain their full-rank or dense frozen parameters. Their rows therefore measure effective-dimensional difficulty rather than a changed low-dimensional surrogate.

6. For Tangent-Actor, randomness_decomposition.csv reports within-Actor FK Monte Carlo noise for the single fresh training run. Between-training-seed variance is intentionally not estimated by the amended primary protocol. The median inferred MC fraction is 1.000e+00.

7. The accuracy-spread rule identifies b6_fullrank, dense_fourier as the first PDEs for a separate future accuracy-compute frontier.

8. HJB-Log and Allen-Cahn Deep BSDE runs are native source-family reproductions; B6 and Dense Fourier are faithful transfers. Allen-Cahn Tangent-Actor uses the explicitly audited value-Actor autograd martingale-CV extension. Standard FH-MLP is retained without truth-based configuration selection.

## Audit and reproducibility

Final artifact audit: VALIDATED.

All query hashes, references, source and amended dispatch manifests, algorithm hashes, raw predictions, formal repeats, and Slurm/task timing records are checked by artifact_audit.json. Numerical failures and incomplete artifacts are reported rather than replaced or omitted.
