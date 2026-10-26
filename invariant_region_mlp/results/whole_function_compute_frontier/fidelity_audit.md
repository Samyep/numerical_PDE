# Fidelity audit

Overall: **FAIL**

The checks below are mechanical; a failed check is not silently waived.

| # | Check | Result |
|---:|---|---|
| 1 | monotone checkpoint times | PASS |
| 2 | checkpoint weights match recorded hashes | PASS |
| 3 | predictions bind to recorded checkpoints/artifacts | PASS |
| 4 | TA timing includes label generation and fitting | PASS |
| 5 | truth-free training/selection and global truth gate | PASS |
| 6 | original-budget C0 consistency | FAIL |
| 7 | completion and numerical-failure records agree with predictions | PASS |
| 8 | checkpoint timing components and scheduler accounting reconcile | PASS |
| 9 | common frozen protocols and input hashes on both PDEs | PASS |
| 10 | no new SDGD-PINN job submitted | PASS |

## Post-freeze diagnostics

These diagnostics explain failures but never change a frozen verdict or enter a frontier.

- Dense Fourier DPI historical-callback RNG replay: historical RMSE `0.0646384146`, frontier RMSE `0.115137165`, RNG-replay RMSE `0.0662925394`; frozen Check 6 remains unchanged.

```json
{
  "checks": [
    {
      "failures": [],
      "id": 1,
      "name": "monotone checkpoint times",
      "pass": true
    },
    {
      "failures": [],
      "id": 2,
      "name": "checkpoint weights match recorded hashes",
      "pass": true
    },
    {
      "failures": [],
      "id": 3,
      "name": "predictions bind to recorded checkpoints/artifacts",
      "pass": true
    },
    {
      "failures": [],
      "id": 4,
      "name": "TA timing includes label generation and fitting",
      "pass": true
    },
    {
      "failures": [],
      "id": 5,
      "name": "truth-free training/selection and global truth gate",
      "pass": true
    },
    {
      "details": [
        {
          "absolute_difference": 2.816314953419394e-05,
          "historical_rmse": 0.0007036985506070463,
          "historical_training_seconds": 669.0896890976001,
          "method": "tangent_actor",
          "new_inference_replicates": 2,
          "new_rmse": 0.0007318617001412402,
          "new_rmse_sd_across_inference_replicates": 3.9828708030398504e-05,
          "new_training_plus_selection_seconds": 694.1233984371647,
          "pass": true,
          "pde": "b6_fullrank",
          "relative_difference_to_historical": 0.04002161083023259,
          "tolerance": 0.05,
          "uncertainty_scope": "inference-only; one training seed in both studies"
        },
        {
          "absolute_difference": 0.017782676567985906,
          "historical_rmse": 0.1425682861836468,
          "historical_training_seconds": 537.9383489890024,
          "method": "dpi",
          "new_inference_replicates": 1,
          "new_rmse": 0.1247856096156609,
          "new_rmse_sd_across_inference_replicates": 0.0,
          "new_training_plus_selection_seconds": 286.14510092069395,
          "pass": true,
          "pde": "b6_fullrank",
          "relative_difference_to_historical": 0.12473094152986777,
          "tolerance": 0.0712841430918234,
          "uncertainty_scope": "no replicate uncertainty estimate; one training seed in both studies"
        },
        {
          "absolute_difference": 0.03539253099748854,
          "historical_rmse": 0.19208899694342052,
          "historical_training_seconds": 235.67849746206775,
          "method": "pinn_hte",
          "new_inference_replicates": 1,
          "new_rmse": 0.15669646594593198,
          "new_rmse_sd_across_inference_replicates": 0.0,
          "new_training_plus_selection_seconds": 201.73992082942277,
          "pass": true,
          "pde": "b6_fullrank",
          "relative_difference_to_historical": 0.18425069400468236,
          "tolerance": 0.09604449847171026,
          "uncertainty_scope": "no replicate uncertainty estimate; one training seed in both studies"
        },
        {
          "absolute_difference": 9.99293211790288e-05,
          "historical_rmse": 0.0016470558953415148,
          "historical_training_seconds": 1357.8923614919186,
          "method": "tangent_actor",
          "new_inference_replicates": 2,
          "new_rmse": 0.0017469852165205436,
          "new_rmse_sd_across_inference_replicates": 0.00014132140129011935,
          "new_training_plus_selection_seconds": 1083.3854689649306,
          "pass": true,
          "pde": "dense_fourier",
          "relative_difference_to_historical": 0.06067148143646248,
          "tolerance": 0.05,
          "uncertainty_scope": "inference-only; one training seed in both studies"
        },
        {
          "absolute_difference": 0.05049875017250341,
          "historical_rmse": 0.0646384146456871,
          "historical_training_seconds": 633.6548566620331,
          "method": "dpi",
          "new_inference_replicates": 1,
          "new_rmse": 0.11513716481819052,
          "new_rmse_sd_across_inference_replicates": 0.0,
          "new_training_plus_selection_seconds": 304.25195059599355,
          "pass": false,
          "pde": "dense_fourier",
          "relative_difference_to_historical": 0.7812498256541454,
          "tolerance": 0.05,
          "uncertainty_scope": "no replicate uncertainty estimate; one training seed in both studies"
        },
        {
          "absolute_difference": 0.003841257402348086,
          "historical_rmse": 0.07710316643174087,
          "historical_training_seconds": 290.03436847403646,
          "method": "pinn_hte",
          "new_inference_replicates": 1,
          "new_rmse": 0.08094442383408895,
          "new_rmse_sd_across_inference_replicates": 0.0,
          "new_training_plus_selection_seconds": 210.8713668207638,
          "pass": true,
          "pde": "dense_fourier",
          "relative_difference_to_historical": 0.049819710137958294,
          "tolerance": 0.05,
          "uncertainty_scope": "no replicate uncertainty estimate; one training seed in both studies"
        }
      ],
      "id": 6,
      "name": "original-budget C0 consistency",
      "pass": false
    },
    {
      "details": {
        "incomplete_trajectories": [],
        "nonfinite_prediction_or_metric_rows": [],
        "recorded_failure_payloads": [],
        "trajectories_without_reported_metrics": []
      },
      "id": 7,
      "name": "completion and numerical-failure records agree with predictions",
      "pass": true
    },
    {
      "details": {
        "active_seconds_from_accounting": 20828.27771383268,
        "active_seconds_from_trajectories": 20828.27771383268,
        "active_seconds_reconciled": true,
        "allocated_seconds_cover_active_training": true,
        "checkpoint_timing_failures": [],
        "failed_or_cancelled_retry_gpu_seconds_included": 678.0,
        "live_gpu_jobs_at_accounting_snapshot": [],
        "metric_timing_failures": [],
        "totals": [
          {
            "active_training_hours": 5.785632698286855,
            "active_training_seconds": 20828.27771383268,
            "allocated_gpu_hours": 8.156944444444445,
            "allocated_gpu_seconds": 29365.0,
            "elapsed": "",
            "gpu_count": "",
            "job_id": "TOTAL",
            "job_name": "all_wcf_jobs",
            "job_wall_seconds": "",
            "state": "SUMMARY",
            "within_18_gpu_hour_cap": true
          }
        ],
        "within_18_gpu_hour_cap": true
      },
      "id": 8,
      "name": "checkpoint timing components and scheduler accounting reconcile",
      "pass": true
    },
    {
      "details": {
        "configuration_snapshot_failures": [],
        "direct_input_hash_rows": 79,
        "embedded_prediction_array_hash_rows": 35,
        "input_hash_failures": [],
        "missing_cells": [],
        "unexpected_cells": []
      },
      "id": 9,
      "name": "common frozen protocols and input hashes on both PDEs",
      "pass": true,
      "protocol_sha256": "54a76a89c475f1c5ef11b1ba8279eb82ad78e2a53d0def320a2dcbabe73e1422"
    },
    {
      "details": {
        "accounted_job_mentions_sdgd": false,
        "manifest_no_sdgd_jobs": true,
        "submitted_job_mentions_sdgd": false
      },
      "id": 10,
      "name": "no new SDGD-PINN job submitted",
      "pass": true
    }
  ],
  "completed_utc": "2026-10-26T17:04:16.749257+00:00",
  "overall_pass": false,
  "post_freeze_diagnostics": [
    {
      "changes_frozen_check_6_verdict": false,
      "classification": "POST_FREEZE_DIAGNOSTIC_EXCLUDED_FROM_FRONTIER",
      "excluded_from_frontier": true,
      "frontier_rmse": 0.11513716481819052,
      "historical_rmse": 0.0646384146456871,
      "name": "Dense Fourier DPI historical-callback RNG replay",
      "prediction_artifact_hash_valid": true,
      "relative_rmse_difference_from_historical": 0.025590428389247428,
      "rng_replay_callback_calls": 1280,
      "rng_replay_rmse": 0.06629253936687204,
      "status": "COMPLETE"
    }
  ],
  "revision": "whole_function_compute_frontier_v1"
}
```
