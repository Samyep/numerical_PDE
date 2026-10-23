# ICML Main High-Dimensional Expansion: B6 and Dense Fourier

**Execution state:** COMPLETE
**Cells:** 40 complete, 0 pending, 0 out of time, 0 OOM, 0 failed, 0 budget-truncated.
**Expansion allocation:** 26.731 GPU hours under the separate 120-hour cap.
**Artifact audit:** PASS.

The primary Tangent-Actor row is the reference-free selected Actor evaluated once with 1,024 paths and 64 Euler steps (replicate 0). Replicate 1 is reported separately and is not averaged into a fictitious higher-path estimator. Exact solutions are used only after selection for evaluation and for the pre-frozen B6 calibration/benchmark characterization.

## Current formal results

| PDE | d | Method | Status | Value RMSE | Relative L2 | GPU h |
|---|---:|---|---|---:|---:|---:|
| B6 Full-Rank HJB | 20 | Tangent-Actor + MCV | COMPLETE | 0.000356648 | 0.00050307 | 0.646 |
| B6 Full-Rank HJB | 20 | DPI | COMPLETE | 0.0084997 | 0.0119893 | 0.100 |
| B6 Full-Rank HJB | 20 | PINN-HTE | COMPLETE | 0.00732977 | 0.010339 | 0.065 |
| B6 Full-Rank HJB | 20 | SDGD-PINN | COMPLETE | 0.000243414 | 0.000343347 | 1.524 |
| B6 Full-Rank HJB | 50 | Tangent-Actor + MCV | COMPLETE | 0.000796705 | 0.000745553 | 0.629 |
| B6 Full-Rank HJB | 50 | DPI | COMPLETE | 0.047715 | 0.0446515 | 0.117 |
| B6 Full-Rank HJB | 50 | PINN-HTE | COMPLETE | 0.0335646 | 0.0314096 | 0.067 |
| B6 Full-Rank HJB | 50 | SDGD-PINN | COMPLETE | 0.00222581 | 0.0020829 | 1.610 |
| B6 Full-Rank HJB | 100 | Tangent-Actor + MCV | COMPLETE | 0.000703699 | 0.000622747 | 0.641 |
| B6 Full-Rank HJB | 100 | DPI | COMPLETE | 0.142568 | 0.126168 | 0.210 |
| B6 Full-Rank HJB | 100 | PINN-HTE | COMPLETE | 0.192089 | 0.169992 | 0.084 |
| B6 Full-Rank HJB | 100 | SDGD-PINN | COMPLETE | 0.0113061 | 0.0100055 | 1.989 |
| B6 Full-Rank HJB | 200 | Tangent-Actor + MCV | COMPLETE | 0.00202258 | 0.00115139 | 0.635 |
| B6 Full-Rank HJB | 200 | DPI | COMPLETE | 0.478418 | 0.272349 | 0.203 |
| B6 Full-Rank HJB | 200 | PINN-HTE | COMPLETE | 0.752889 | 0.428597 | 0.099 |
| B6 Full-Rank HJB | 200 | SDGD-PINN | COMPLETE | 0.0749329 | 0.0426571 | 1.976 |
| B6 Full-Rank HJB | 400 | Tangent-Actor + MCV | COMPLETE | 0.00177859 | 0.000997321 | 0.661 |
| B6 Full-Rank HJB | 400 | DPI | COMPLETE | 0.846719 | 0.474787 | 0.350 |
| B6 Full-Rank HJB | 400 | PINN-HTE | COMPLETE | 1.37024 | 0.768345 | 0.103 |
| B6 Full-Rank HJB | 400 | SDGD-PINN | COMPLETE | 0.0941086 | 0.0527703 | 1.897 |
| Dense Fourier Joint | 20 | Tangent-Actor + MCV | COMPLETE | 0.00235343 | 0.00346857 | 1.116 |
| Dense Fourier Joint | 20 | DPI | COMPLETE | 0.0206245 | 0.0303971 | 0.104 |
| Dense Fourier Joint | 20 | PINN-HTE | COMPLETE | 0.048961 | 0.0721605 | 0.096 |
| Dense Fourier Joint | 20 | SDGD-PINN | COMPLETE | 0.0134925 | 0.0198857 | 0.902 |
| Dense Fourier Joint | 50 | Tangent-Actor + MCV | COMPLETE | 0.00221841 | 0.00333184 | 1.086 |
| Dense Fourier Joint | 50 | DPI | COMPLETE | 0.0624759 | 0.093833 | 0.144 |
| Dense Fourier Joint | 50 | PINN-HTE | COMPLETE | 0.0655093 | 0.0983888 | 0.071 |
| Dense Fourier Joint | 50 | SDGD-PINN | COMPLETE | 0.051308 | 0.0770599 | 0.887 |
| Dense Fourier Joint | 100 | Tangent-Actor + MCV | COMPLETE | 0.00164706 | 0.00119031 | 1.114 |
| Dense Fourier Joint | 100 | DPI | COMPLETE | 0.0646384 | 0.0467136 | 0.194 |
| Dense Fourier Joint | 100 | PINN-HTE | COMPLETE | 0.0771032 | 0.0557218 | 0.129 |
| Dense Fourier Joint | 100 | SDGD-PINN | COMPLETE | 0.0741177 | 0.0535642 | 0.916 |
| Dense Fourier Joint | 200 | Tangent-Actor + MCV | COMPLETE | 0.00426578 | 0.00778318 | 1.078 |
| Dense Fourier Joint | 200 | DPI | COMPLETE | 0.24618 | 0.44917 | 0.345 |
| Dense Fourier Joint | 200 | PINN-HTE | COMPLETE | 0.262188 | 0.478378 | 0.102 |
| Dense Fourier Joint | 200 | SDGD-PINN | COMPLETE | 0.208589 | 0.380583 | 0.905 |
| Dense Fourier Joint | 400 | Tangent-Actor + MCV | COMPLETE | 0.005091 | 0.00941985 | 1.184 |
| Dense Fourier Joint | 400 | DPI | COMPLETE | 0.390364 | 0.722289 | 0.562 |
| Dense Fourier Joint | 400 | PINN-HTE | COMPLETE | 0.400522 | 0.741084 | 0.153 |
| Dense Fourier Joint | 400 | SDGD-PINN | COMPLETE | 0.23285 | 0.43084 | 0.920 |

## Interpretation guardrails

DPI, PINN-HTE, and SDGD-PINN rows are faithful transfers of pinned source protocols to new equations. They support a controlled comparison in this suite; they do not establish that Tangent-Actor universally outperforms every native experiment reported in those papers.

Four-method settings currently complete: 10 / 10. Tangent-Actor currently has the lowest RMSE in 9 of those settings (ties split fractionally). This statement will update automatically as cells finish.

The effective-rank, outside-top-10 variance, nonlinear-skill, iteration, MC-noise, and accuracy-compute evidence is stored in the accompanying CSVs and Figures A–C. Missing or capped cells stay explicitly non-numerical.

## Reproducibility

- Frozen training seed: 0; frozen whole-function bank size: 960 per PDE/dimension.
- Dimensions: 20, 50, 100, 200, 400.
- Per-cell cap: 6 GPU hours; expansion cap: 120 GPU hours.
- Heavy banks, checkpoints, predictions, and logs remain in the workspace; committed manifests bind them by SHA-256.
- Baseline classification and source commits are recorded in `baseline_fidelity_audit.md`.
