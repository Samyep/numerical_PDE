# ICML Main High-Dimensional Expansion: B6 and Dense Fourier

**Execution state:** IN_PROGRESS_OR_TERMINALLY_INCOMPLETE
**Cells:** 0 complete, 40 pending, 0 out of time, 0 OOM, 0 failed.
**Expansion allocation:** 0.020 GH200 hours under the separate 120-hour cap.
**Artifact audit:** PASS.

The primary Tangent-Actor row is the reference-free selected Actor evaluated once with 1,024 paths and 64 Euler steps (replicate 0). Replicate 1 is reported separately and is not averaged into a fictitious higher-path estimator. Exact solutions are used only after selection for evaluation and for the pre-frozen B6 calibration/benchmark characterization.

## Current formal results

| PDE | d | Method | Status | Value RMSE | Relative L2 | GH200 h |
|---|---:|---|---|---:|---:|---:|
| B6 Full-Rank HJB | 20 | Tangent-Actor + MCV | PENDING | — | — | 0.000 |
| B6 Full-Rank HJB | 20 | DPI | PENDING | — | — | 0.000 |
| B6 Full-Rank HJB | 20 | PINN-HTE | PENDING | — | — | 0.000 |
| B6 Full-Rank HJB | 20 | SDGD-PINN | PENDING | — | — | 0.000 |
| B6 Full-Rank HJB | 50 | Tangent-Actor + MCV | PENDING | — | — | 0.000 |
| B6 Full-Rank HJB | 50 | DPI | PENDING | — | — | 0.000 |
| B6 Full-Rank HJB | 50 | PINN-HTE | PENDING | — | — | 0.000 |
| B6 Full-Rank HJB | 50 | SDGD-PINN | PENDING | — | — | 0.000 |
| B6 Full-Rank HJB | 100 | Tangent-Actor + MCV | PENDING | — | — | 0.000 |
| B6 Full-Rank HJB | 100 | DPI | PENDING | — | — | 0.000 |
| B6 Full-Rank HJB | 100 | PINN-HTE | PENDING | — | — | 0.000 |
| B6 Full-Rank HJB | 100 | SDGD-PINN | PENDING | — | — | 0.000 |
| B6 Full-Rank HJB | 200 | Tangent-Actor + MCV | PENDING | — | — | 0.000 |
| B6 Full-Rank HJB | 200 | DPI | PENDING | — | — | 0.000 |
| B6 Full-Rank HJB | 200 | PINN-HTE | PENDING | — | — | 0.000 |
| B6 Full-Rank HJB | 200 | SDGD-PINN | PENDING | — | — | 0.000 |
| B6 Full-Rank HJB | 400 | Tangent-Actor + MCV | PENDING | — | — | 0.000 |
| B6 Full-Rank HJB | 400 | DPI | PENDING | — | — | 0.000 |
| B6 Full-Rank HJB | 400 | PINN-HTE | PENDING | — | — | 0.000 |
| B6 Full-Rank HJB | 400 | SDGD-PINN | PENDING | — | — | 0.000 |
| Dense Fourier Joint | 20 | Tangent-Actor + MCV | PENDING | — | — | 0.000 |
| Dense Fourier Joint | 20 | DPI | PENDING | — | — | 0.000 |
| Dense Fourier Joint | 20 | PINN-HTE | PENDING | — | — | 0.000 |
| Dense Fourier Joint | 20 | SDGD-PINN | PENDING | — | — | 0.000 |
| Dense Fourier Joint | 50 | Tangent-Actor + MCV | PENDING | — | — | 0.000 |
| Dense Fourier Joint | 50 | DPI | PENDING | — | — | 0.000 |
| Dense Fourier Joint | 50 | PINN-HTE | PENDING | — | — | 0.000 |
| Dense Fourier Joint | 50 | SDGD-PINN | PENDING | — | — | 0.000 |
| Dense Fourier Joint | 100 | Tangent-Actor + MCV | PENDING | — | — | 0.000 |
| Dense Fourier Joint | 100 | DPI | PENDING | — | — | 0.000 |
| Dense Fourier Joint | 100 | PINN-HTE | PENDING | — | — | 0.000 |
| Dense Fourier Joint | 100 | SDGD-PINN | PENDING | — | — | 0.000 |
| Dense Fourier Joint | 200 | Tangent-Actor + MCV | PENDING | — | — | 0.000 |
| Dense Fourier Joint | 200 | DPI | PENDING | — | — | 0.000 |
| Dense Fourier Joint | 200 | PINN-HTE | PENDING | — | — | 0.000 |
| Dense Fourier Joint | 200 | SDGD-PINN | PENDING | — | — | 0.000 |
| Dense Fourier Joint | 400 | Tangent-Actor + MCV | PENDING | — | — | 0.000 |
| Dense Fourier Joint | 400 | DPI | PENDING | — | — | 0.000 |
| Dense Fourier Joint | 400 | PINN-HTE | PENDING | — | — | 0.000 |
| Dense Fourier Joint | 400 | SDGD-PINN | PENDING | — | — | 0.000 |

## Interpretation guardrails

DPI, PINN-HTE, and SDGD-PINN rows are faithful transfers of pinned source protocols to new equations. They support a controlled comparison in this suite; they do not establish that Tangent-Actor universally outperforms every native experiment reported in those papers.

Four-method settings currently complete: 0 / 10. Tangent-Actor currently has the lowest RMSE in 0 of those settings (ties split fractionally). This statement will update automatically as cells finish.

The effective-rank, outside-top-10 variance, nonlinear-skill, iteration, MC-noise, and accuracy-compute evidence is stored in the accompanying CSVs and Figures A–C. Missing or capped cells stay explicitly non-numerical.

## Reproducibility

- Frozen training seed: 0; frozen whole-function bank size: 960 per PDE/dimension.
- Dimensions: 20, 50, 100, 200, 400.
- Per-cell cap: 6 GH200 hours; expansion cap: 120 GH200 hours.
- Heavy banks, checkpoints, predictions, and logs remain in the workspace; committed manifests bind them by SHA-256.
- Baseline classification and source commits are recorded in `baseline_fidelity_audit.md`.
