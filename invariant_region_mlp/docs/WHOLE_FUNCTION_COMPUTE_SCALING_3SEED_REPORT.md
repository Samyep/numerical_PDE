# Three-seed hyperparameter-specific accuracy--compute scaling

## Completeness and audit

| Check | Result |
|---|---:|
| Artifact audit | **PASS** |
| Formal runs | 168 / 168 |
| Aggregated curve points | 76 / 76 |
| Training seeds per point | 3 (0, 1, 2) |
| Frozen calibration grid | FROZEN_BEFORE_FORMAL_TRAINING |
| Resource cap | PASS |

Each point is the mean of three separately computed per-seed RMSE values; error bars are their sample standard deviation. Predictions are not ensembled.

## Descriptive best observed points

These post-evaluation minima summarize the measured grid only; they were not used for training, stopping, checkpoint selection, or a test-selected envelope.

| PDE | Method | Sweep | Parameter | Mean time (min) | Mean RMSE | Ratio to TA best |
|---|---|---|---:|---:|---:|---:|
| b6_fullrank | tangent_actor | C0 | 12 | 32.410 | 0.00071568086 | 1.00x |
| b6_fullrank | dpi | C2 | 128 | 18.054 | 0.047256172 | 66.03x |
| b6_fullrank | pinn_hte | C0 | 240000 | 26.396 | 0.092138954 | 128.74x |
| dense_fourier | tangent_actor | C2 | 40000 | 26.126 | 0.0016853092 | 1.00x |
| dense_fourier | dpi | C0 | 450 | 25.999 | 0.061796779 | 36.67x |
| dense_fourier | pinn_hte | C2 | 4096 | 25.822 | 0.068057781 | 40.38x |

## Main results

| PDE | Method | Sweep | Parameter | Mean time (min) | Mean RMSE | Seed SD |
|---|---|---|---:|---:|---:|---:|
| b6_fullrank | dpi | C0 | 80 | 6.152 | 0.13743992 | 0.0129 |
| b6_fullrank | dpi | C0 | 160 | 13.129 | 0.14890907 | 0.0658 |
| b6_fullrank | dpi | C0 | 320 | 24.411 | 0.1843084 | 0.0679 |
| b6_fullrank | dpi | C0 | 450 | 32.351 | 0.14979788 | 0.0578 |
| b6_fullrank | dpi | C1 | 512 | 5.105 | 0.16159003 | 0.0396 |
| b6_fullrank | dpi | C1 | 16384 | 7.615 | 0.14388911 | 0.0166 |
| b6_fullrank | dpi | C1 | 32768 | 10.870 | 0.13275198 | 0.0259 |
| b6_fullrank | dpi | C1 | 65536 | 30.334 | 0.13977332 | 0.0383 |
| b6_fullrank | dpi | C1 | 131072 | 22.452 | 0.12416244 | 0.0338 |
| b6_fullrank | dpi | C2 | 4 | 3.311 | 0.14269588 | 0.018 |
| b6_fullrank | dpi | C2 | 32 | 6.222 | 0.15527597 | 0.00486 |
| b6_fullrank | dpi | C2 | 64 | 9.637 | 0.10816922 | 0.0209 |
| b6_fullrank | dpi | C2 | 128 | 18.054 | 0.047256172 | 0.00474 |
| b6_fullrank | dpi | C2 | 256 | 28.803 | 0.1129611 | 0.0305 |
| b6_fullrank | pinn_hte | C0 | 30000 | 3.307 | 0.21792051 | 0.0543 |
| b6_fullrank | pinn_hte | C0 | 60000 | 6.606 | 0.11995371 | 0.0218 |
| b6_fullrank | pinn_hte | C0 | 120000 | 13.209 | 0.18623126 | 0.0524 |
| b6_fullrank | pinn_hte | C0 | 240000 | 26.396 | 0.092138954 | 0.04 |
| b6_fullrank | pinn_hte | C1 | 64 | 7.096 | 0.15585441 | 0.0882 |
| b6_fullrank | pinn_hte | C1 | 128 | 8.578 | 0.20591067 | 0.148 |
| b6_fullrank | pinn_hte | C1 | 256 | 26.567 | 0.19963294 | 0.0181 |
| b6_fullrank | pinn_hte | C1 | 512 | 27.384 | 0.1849748 | 0.0944 |
| b6_fullrank | pinn_hte | C2 | 512 | 3.404 | 0.21792051 | 0.0543 |
| b6_fullrank | pinn_hte | C2 | 2048 | 12.799 | 0.12333936 | 0.0467 |
| b6_fullrank | pinn_hte | C2 | 4096 | 25.577 | 0.13612064 | 0.0284 |
| b6_fullrank | pinn_hte | C2 | 8192 | 50.806 | 0.16952913 | 0.0302 |
| b6_fullrank | pinn_hte | C2 | 16384 | 100.716 | 0.10268974 | 0.0122 |
| b6_fullrank | tangent_actor | C0 | 2 | 5.400 | 0.00072765051 | 8.79e-06 |
| b6_fullrank | tangent_actor | C0 | 4 | 10.757 | 0.0007200539 | 1.3e-05 |
| b6_fullrank | tangent_actor | C0 | 6 | 16.184 | 0.00071805609 | 1.19e-05 |
| b6_fullrank | tangent_actor | C0 | 8 | 21.655 | 0.00071666207 | 1.38e-05 |
| b6_fullrank | tangent_actor | C0 | 12 | 32.410 | 0.00071568086 | 1.36e-05 |
| b6_fullrank | tangent_actor | C1 | 64 | 5.640 | 0.00072228849 | 1.35e-05 |
| b6_fullrank | tangent_actor | C1 | 128 | 10.765 | 0.0007200539 | 1.3e-05 |
| b6_fullrank | tangent_actor | C1 | 192 | 15.852 | 0.00071927927 | 1.31e-05 |
| b6_fullrank | tangent_actor | C1 | 256 | 20.951 | 0.00072061531 | 1.34e-05 |
| b6_fullrank | tangent_actor | C1 | 384 | 31.292 | 0.00071960252 | 1.32e-05 |
| b6_fullrank | tangent_actor | C2 | 1250 | 10.466 | 0.00073507493 | 1.36e-05 |
| b6_fullrank | tangent_actor | C2 | 40000 | 14.263 | 0.00071787336 | 1.25e-05 |
| b6_fullrank | tangent_actor | C2 | 80000 | 18.145 | 0.00071655214 | 1.31e-05 |
| b6_fullrank | tangent_actor | C2 | 160000 | 25.895 | 0.00071576974 | 1.4e-05 |
| dense_fourier | dpi | C0 | 80 | 4.749 | 0.11078381 | 0.0333 |
| dense_fourier | dpi | C0 | 160 | 9.291 | 0.099395396 | 0.0611 |
| dense_fourier | dpi | C0 | 320 | 18.493 | 0.063192384 | 0.0146 |
| dense_fourier | dpi | C0 | 450 | 25.999 | 0.061796779 | 0.0147 |
| dense_fourier | dpi | C1 | 512 | 4.071 | 0.076695677 | 0.00376 |
| dense_fourier | dpi | C1 | 16384 | 8.193 | 0.093650119 | 0.0302 |
| dense_fourier | dpi | C1 | 32768 | 8.570 | 0.11515494 | 0.0611 |
| dense_fourier | dpi | C1 | 65536 | 13.221 | 0.07397991 | 0.00383 |
| dense_fourier | dpi | C1 | 131072 | 22.379 | 0.12234705 | 0.0401 |
| dense_fourier | dpi | C2 | 16 | 5.008 | 0.11078381 | 0.0333 |
| dense_fourier | dpi | C2 | 64 | 9.441 | 0.097404466 | 0.0237 |
| dense_fourier | dpi | C2 | 128 | 16.188 | 0.063945332 | 0.00879 |
| dense_fourier | dpi | C2 | 256 | 29.599 | 0.07765453 | 0.0246 |
| dense_fourier | pinn_hte | C0 | 30000 | 3.374 | 0.070337809 | 0.00948 |
| dense_fourier | pinn_hte | C0 | 60000 | 6.743 | 0.077262657 | 0.0134 |
| dense_fourier | pinn_hte | C0 | 120000 | 13.475 | 0.087662238 | 0.00856 |
| dense_fourier | pinn_hte | C0 | 240000 | 26.946 | 0.10519939 | 0.0484 |
| dense_fourier | pinn_hte | C1 | 4 | 3.354 | 0.10353622 | 0.038 |
| dense_fourier | pinn_hte | C1 | 64 | 7.322 | 0.087511977 | 0.0113 |
| dense_fourier | pinn_hte | C1 | 128 | 8.621 | 0.076908924 | 0.00544 |
| dense_fourier | pinn_hte | C2 | 256 | 3.331 | 0.093718502 | 0.0201 |
| dense_fourier | pinn_hte | C2 | 2048 | 13.402 | 0.07787646 | 0.0197 |
| dense_fourier | pinn_hte | C2 | 4096 | 25.822 | 0.068057781 | 0.0103 |
| dense_fourier | tangent_actor | C0 | 1 | 5.235 | 0.00181545 | 9.02e-05 |
| dense_fourier | tangent_actor | C0 | 2 | 10.626 | 0.0017591236 | 9.37e-05 |
| dense_fourier | tangent_actor | C0 | 3 | 15.976 | 0.0017693691 | 8.8e-05 |
| dense_fourier | tangent_actor | C0 | 4 | 21.314 | 0.0017584167 | 9.63e-05 |
| dense_fourier | tangent_actor | C0 | 6 | 32.053 | 0.0017270626 | 7.43e-05 |
| dense_fourier | tangent_actor | C1 | 32 | 5.862 | 0.0018050703 | 8.69e-05 |
| dense_fourier | tangent_actor | C1 | 64 | 10.889 | 0.0017680528 | 9.99e-05 |
| dense_fourier | tangent_actor | C1 | 128 | 21.112 | 0.0017584167 | 9.63e-05 |
| dense_fourier | tangent_actor | C1 | 192 | 31.388 | 0.0017626905 | 9.19e-05 |
| dense_fourier | tangent_actor | C2 | 1250 | 20.546 | 0.0017838204 | 0.000131 |
| dense_fourier | tangent_actor | C2 | 40000 | 26.126 | 0.0016853092 | 4.43e-05 |
| dense_fourier | tangent_actor | C2 | 80000 | 31.092 | 0.0016951761 | 7.26e-05 |

## Interpretation boundary

These are one-parameter scaling curves, not globally optimized Pareto frontiers. C0 points are checkpoints from one trajectory per seed; C1/C2 points are independent fresh runs. No running minimum, smoothing, interpolation, or test-selected envelope is used.
The plotted x-axis is synchronized active training time, as frozen. Reference-free TA selector time, inference time, scheduler wall time, and allocated GPU time are retained separately rather than silently folded into that axis.

## Timing calibration and frozen grids

Calibration used no accuracy and projected 51.23 formal H100 GPU-hours before formal training.

| PDE | Method | Sweep | Frozen parameter values |
|---|---|---|---|
| b6_fullrank | dpi | C0 | 80, 160, 320, 450 |
| b6_fullrank | dpi | C1 | 512, 16384, 32768, 65536, 131072 |
| b6_fullrank | dpi | C2 | 4, 32, 64, 128, 256 |
| b6_fullrank | pinn_hte | C0 | 30000, 60000, 120000, 240000 |
| b6_fullrank | pinn_hte | C1 | 64, 128, 256, 512 |
| b6_fullrank | pinn_hte | C2 | 512, 2048, 4096, 8192, 16384 |
| b6_fullrank | tangent_actor | C0 | 2, 4, 6, 8, 12 |
| b6_fullrank | tangent_actor | C1 | 64, 128, 192, 256, 384 |
| b6_fullrank | tangent_actor | C2 | 1250, 40000, 80000, 160000 |
| dense_fourier | dpi | C0 | 80, 160, 320, 450 |
| dense_fourier | dpi | C1 | 512, 16384, 32768, 65536, 131072 |
| dense_fourier | dpi | C2 | 16, 64, 128, 256 |
| dense_fourier | pinn_hte | C0 | 30000, 60000, 120000, 240000 |
| dense_fourier | pinn_hte | C1 | 4, 64, 128 |
| dense_fourier | pinn_hte | C2 | 256, 2048, 4096 |
| dense_fourier | tangent_actor | C0 | 1, 2, 3, 4, 6 |
| dense_fourier | tangent_actor | C1 | 32, 64, 128, 192 |
| dense_fourier | tangent_actor | C2 | 1250, 40000, 80000 |

## Resource accounting and failures

Total allocated H100 time (including calibration failures and retries): **69.424 GPU-hours**.

Scheduler rows not marked COMPLETED: 66.

Audit failures and recorded deviations:

- Final artifact audit: no failures.
- Scheduler history retained and counted: CANCELLED=5, FAILED=58, TIMEOUT=3.
- The FAILED and CANCELLED scheduler rows are calibration/preflight attempts (mostly zero- or few-second startup failures on subsequently excluded nodes, plus curtailed calibration arrays); they do not represent missing formal runs.
- Formal deviation: the three B6 PINN-HTE C2 runs at batch size 16384 timed out at 70:12, then completed with the identical numerical configuration after only the Slurm wall-time limit was raised; both attempts are included in resource accounting.
- Timing-calibration deviation: points beyond the intended 25--35 minute endpoint are retained (b6_fullrank/pinn_hte/C2=50.8 min, b6_fullrank/pinn_hte/C2=100.7 min).
- A pre-fix artifact audit failed because the evaluator skipped pre-registered DPI C0 checkpoints at 80/160/320. The FAIL record is preserved in `artifact_audit_pre_dpi_c0_evaluation_fix.json`; the evaluator-only fix did not retrain models, and the six C0 trajectories were re-evaluated in Slurm job 6619464.
- 3 pre-final DPI chunk-timing directories were moved to `calibration/archived_pre_final_chunk/` and preserved; they were not used to freeze the final grid.

## DPI large-sample chunking validation

Status: **PASS**; all nominal samples used: `True`; one-chunk maximum absolute difference `0.0`; true multi-chunk maximum absolute difference `2.7755575615628914e-17`.

## Dense Fourier DPI provenance

The preceding study's historical C0 consistency check remains FAIL. Its RNG-only diagnostic explained 96.72% of the discrepancy but was not bitwise. This study disables the historical exact-solution callback (`EVAL.FREQ=null`) and averages new truth-gated runs over seeds 0, 1, and 2.

## Provenance

Calibration manifest SHA-256: `6a6d2bef704f35c3976400bf59b50182d2038429f8d95cc353f2a72a6665f1d5`.
Run manifest SHA-256: `69a8c185e5917383c7b110935ff876b37e97d69989d757686b4fec7790a1c16f`.
Checked-out training revision: `425f8245a4b8a4c66696934a18917a615052006f`.
Frozen input hashes: `{"b6_fullrank": {"parameter_file_sha256": "21b508c30ad63de912979cc9c85ecda0d81aedab1357c88a40e04692495972a4", "test_bank_file_sha256": "6ee6c56ca3d4f413764f28865cb8c00823ea2881ebc90f2643ff4a36fe878730", "test_input_sha256": "0ee05fff7ce16f30d0d31563f60df7db3eb23a9c46bbe7bc9be8eb3dbee1311b", "truth_value_sha256": "6262db7f6b7ee2c41cb6011dadb7288226157533488f708be7d6ea6773f2e13a"}, "dense_fourier": {"parameter_file_sha256": "042588d7c3cc8175cc7b9344284725443650837f72008806c668d35b9a005307", "test_bank_file_sha256": "998fde8b2c0556de4b8fa9a55dec5e489224828cf74e5f3b36f77cb1a8e2c011", "test_input_sha256": "d409296128f48c5ae00c5de2c6ddd998768c01dff18a89c385245a50ff581b70", "truth_value_sha256": "3b33ea36e21595b5c20e5d91097d1fadfcbf3f7884c0328f20579e69995fc816"}}`.
Seed metrics SHA-256: `d93d5e5472d6103b5ded2228c60b35ebec15170d92b57100f3b65d0ca3a151e9`.
Scaling summary SHA-256: `fab0b32170ffa92fc0e0842f29f35bea864b62bba3fc078e6d56ca045377b853`.
Figure PDF SHA-256: `4be9d6c0ad9c3719b60ec933939843de35c258d56d9a54291c4dfca8201ac7f7`.
DPI chunking validation SHA-256: `bba17d298ab8a271925ef446442dc5ed76c5404d46c09a903227ebc05a3a9b00`.
