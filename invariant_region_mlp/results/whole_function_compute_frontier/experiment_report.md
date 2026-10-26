# Two-PDE 100D Accuracy--Compute Frontier

## Outcome

All entries below are measured checkpoints. No RMSE values were interpolated, smoothed, or repeated to create overlap.

| PDE | Method | Original-budget output | Training + selection (min) | RMSE | RelMAE |
|---|---|---:|---:|---:|---:|
| B6 Full-Rank HJB | Tangent-Actor | a4 | 11.569 | 0.000731862 | 0.000492141 |
| B6 Full-Rank HJB | DPI | picard_0080 | 4.769 | 0.124786 | 0.120237 |
| B6 Full-Rank HJB | PINN-HTE | epoch_0030000 | 3.362 | 0.156696 | 0.153139 |
| Dense Fourier | Tangent-Actor | a3 | 23.441 | 0.00174699 | 0.000822212 |
| Dense Fourier | DPI | picard_0080 | 5.071 | 0.115137 | 0.0682635 |
| Dense Fourier | PINN-HTE | epoch_0030000 | 3.515 | 0.0809444 | 0.0482655 |

## Headline measured comparisons

- **B6 Full-Rank HJB original-budget outputs:** Tangent-Actor RMSE `0.000731862` is **170.5x** lower than DPI and **214.1x** lower than PINN-HTE. These checkpoints have their own measured method-specific costs and are not an equal-time comparison.
- **B6 Full-Rank HJB observed C0 overlap (2.04--18.51 min):** the worst measured Tangent-Actor RMSE is `0.00075947`, while the best measured DPI/PINN-HTE RMSE is `0.10075`--a conservative **132.7x** separation using measured checkpoints only (no time interpolation).
- **Dense Fourier original-budget outputs:** Tangent-Actor RMSE `0.00174699` is **65.9x** lower than DPI and **46.3x** lower than PINN-HTE. These checkpoints have their own measured method-specific costs and are not an equal-time comparison.
- **Dense Fourier observed C0 overlap (2.02--15.13 min):** the worst measured Tangent-Actor RMSE is `0.00177371`, while the best measured DPI/PINN-HTE RMSE is `0.0669744`--a conservative **37.8x** separation using measured checkpoints only (no time interpolation).

## C0 anytime trajectories

The primary curves are C0. C1/C2 are sensitivity analyses and are not collapsed into a test-selected envelope.

| PDE | Method | First measured RMSE (min) | Final measured RMSE (min) | Stop reason |
|---|---|---:|---:|---|
| B6 Full-Rank HJB | Tangent-Actor | 0.0010319 (0.11) | 0.00072793 (23.19) | maximum_predeclared_actor_updates |
| B6 Full-Rank HJB | DPI | 0.15809 (2.04) | 0.23642 (18.51) | maximum_predeclared_picard_iterations |
| B6 Full-Rank HJB | PINN-HTE | 1.0926 (0.06) | 0.17125 (30.22) | predeclared_validation_plateau |
| Dense Fourier | Tangent-Actor | 0.0027963 (0.24) | 0.0017146 (34.54) | active_training_cap |
| Dense Fourier | DPI | 0.11493 (2.02) | 0.053313 (18.91) | maximum_predeclared_picard_iterations |
| Dense Fourier | PINN-HTE | 1.4289 (0.12) | 0.10433 (15.13) | predeclared_validation_plateau |

## Configuration sensitivity

Each row is the last genuinely trained checkpoint of that trajectory; this table is descriptive, not a selection rule.

| PDE | Method | Config | Final time (min) | Final RMSE |
|---|---|---|---:|---:|
| B6 Full-Rank HJB | Tangent-Actor | C0 | 23.186 | 0.000727926 |
| B6 Full-Rank HJB | Tangent-Actor | C1 | 22.234 | 0.000733735 |
| B6 Full-Rank HJB | Tangent-Actor | C2 | 11.869 | 0.000730319 |
| B6 Full-Rank HJB | DPI | C0 | 18.509 | 0.236418 |
| B6 Full-Rank HJB | DPI | C1 | 6.527 | 0.160876 |
| B6 Full-Rank HJB | DPI | C2 | 5.091 | 0.123611 |
| B6 Full-Rank HJB | PINN-HTE | C0 | 30.220 | 0.171247 |
| B6 Full-Rank HJB | PINN-HTE | C1 | 30.034 | 0.11941 |
| B6 Full-Rank HJB | PINN-HTE | C2 | 30.122 | 0.187293 |
| Dense Fourier | Tangent-Actor | C0 | 34.538 | 0.00171465 |
| Dense Fourier | Tangent-Actor | C1 | 31.464 | 0.00174204 |
| Dense Fourier | Tangent-Actor | C2 | 22.118 | 0.00175797 |
| Dense Fourier | DPI | C0 | 18.910 | 0.0533129 |
| Dense Fourier | DPI | C1 | 6.494 | 0.0838375 |
| Dense Fourier | DPI | C2 | 5.252 | 0.0834959 |
| Dense Fourier | PINN-HTE | C0 | 15.131 | 0.104331 |
| Dense Fourier | PINN-HTE | C1 | 20.116 | 0.0731126 |
| Dense Fourier | PINN-HTE | C2 | 20.167 | 0.0685889 |

## Observed cost overlap

- **B6 Full-Rank HJB:** C0 observed-range intersection = 2.04--18.51 min.
- **Dense Fourier:** C0 observed-range intersection = 2.02--15.13 min.

## Fidelity and accounting

- Non-result smoke preflight: **PASS**.
- Artifact/fidelity audit: **FAIL**.
- Allocated GPU time: **7.669 H100 GPU-hours**.
- Active measured training time: **5.786 GPU-hours**.
- 18-GPU-hour cap respected: **True**.
- Protocol SHA-256: `54a76a89c475f1c5ef11b1ba8279eb82ad78e2a53d0def320a2dcbabe73e1422`.
- Complete checkpoint, prediction, test-bank, source-revision, and scheduler details are in the companion CSV/JSON artifacts.

### Failed audit checks

- **Check 6 -- original-budget C0 consistency: FAIL.**
  - Dense Fourier / DPI: new original-budget RMSE `0.115137165` versus historical `0.0646384146`; absolute difference `0.0504987502` exceeds the frozen tolerance `0.05`. The failure is retained without changing the threshold.

## Limitations

This is a single-training-seed study (seed 0). Tangent-Actor uses two fixed Monte Carlo inference replicates, but those replicates do not measure training-seed variability. C1/C2 are sensitivity runs, not test-selected hyperparameters. Timing is H100-specific and separates active numerical work from Slurm allocation time. Existing SDGD-PINN fixed-budget results were not rerun and do not appear on these scaling curves.

## Figures

- `figures/b6_100d_accuracy_compute.pdf`
- `figures/dense_fourier_100d_accuracy_compute.pdf`
- `figures/two_pde_accuracy_compute.pdf` and `.png`
- `figures/configuration_ablation.pdf`
- `figures/accuracy_vs_total_cost.pdf`
- `figures/checkpoint_time_coverage.pdf`
- `figures/cost_breakdown.pdf`
