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
| 7 | every numerical failure recorded | PASS |
| 8 | scheduler accounting reconciled and within cap | PASS |
| 9 | common frozen protocols on both PDEs | PASS |
| 10 | no new SDGD-PINN job submitted | PASS |

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
          "historical_rmse": 0.0007036985506070463,
          "method": "tangent_actor",
          "new_rmse": 0.0007318617001412402,
          "pass": true,
          "pde": "b6_fullrank",
          "tolerance": 0.05
        },
        {
          "historical_rmse": 0.1425682861836468,
          "method": "dpi",
          "new_rmse": 0.1247856096156609,
          "pass": true,
          "pde": "b6_fullrank",
          "tolerance": 0.0712841430918234
        },
        {
          "historical_rmse": 0.19208899694342052,
          "method": "pinn_hte",
          "new_rmse": 0.15669646594593198,
          "pass": true,
          "pde": "b6_fullrank",
          "tolerance": 0.09604449847171026
        },
        {
          "historical_rmse": 0.0016470558953415148,
          "method": "tangent_actor",
          "new_rmse": 0.0017469852165205436,
          "pass": true,
          "pde": "dense_fourier",
          "tolerance": 0.05
        },
        {
          "historical_rmse": 0.0646384146456871,
          "method": "dpi",
          "new_rmse": 0.11513716481819052,
          "pass": false,
          "pde": "dense_fourier",
          "tolerance": 0.05
        },
        {
          "historical_rmse": 0.07710316643174087,
          "method": "pinn_hte",
          "new_rmse": 0.08094442383408895,
          "pass": true,
          "pde": "dense_fourier",
          "tolerance": 0.05
        }
      ],
      "id": 6,
      "name": "original-budget C0 consistency",
      "pass": false
    },
    {
      "id": 7,
      "name": "every numerical failure recorded",
      "pass": true,
      "recorded_failures": []
    },
    {
      "details": [
        {
          "active_training_hours": 5.785632698286855,
          "active_training_seconds": 20828.27771383268,
          "allocated_gpu_hours": 7.668611111111111,
          "allocated_gpu_seconds": 27607.0,
          "elapsed": "",
          "gpu_count": "",
          "job_id": "TOTAL",
          "job_name": "all_wcf_jobs",
          "job_wall_seconds": "",
          "state": "SUMMARY",
          "within_18_gpu_hour_cap": true
        }
      ],
      "id": 8,
      "name": "scheduler accounting reconciled and within cap",
      "pass": true
    },
    {
      "id": 9,
      "name": "common frozen protocols on both PDEs",
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
  "completed_utc": "2026-10-26T13:36:58.010210+00:00",
  "overall_pass": false,
  "revision": "whole_function_compute_frontier_v1"
}
```
