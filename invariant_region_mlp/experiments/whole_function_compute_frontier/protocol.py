"""Frozen declarations for the whole-function accuracy--compute frontier."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Final


REVISION: Final = "whole_function_compute_frontier_v1"
BRANCH: Final = "ir-mlp-whole-function-compute-frontier"
BASE_REVISION: Final = "d43b45bdcf068f8fb1e09cf2d9a7e1d5bfd5c038"
UPSTREAM_DPI_REVISION: Final = "e8616c2f531c9f95299b246a13c1180b8f02822a"

PDES: Final = ("b6_fullrank", "dense_fourier")
METHODS: Final = ("tangent_actor", "dpi", "pinn_hte")
CONFIGURATIONS: Final = ("C0", "C1", "C2")
DIMENSION: Final = 100
TRAINING_SEED: Final = 0
TEST_STATES: Final = 960
TARGET_MINUTES: Final = (2, 5, 8, 10, 15, 20, 25, 30)
PRIMARY_MINUTES: Final = (5, 20)
ACTIVE_TRAINING_CAP_SECONDS: Final = 30 * 60
GLOBAL_GPU_SECONDS_CAP: Final = 18 * 3600
C0_CONSISTENCY_ABSOLUTE_RMSE: Final = 0.05
C0_CONSISTENCY_RELATIVE_RMSE: Final = 0.50

TA_EVALUATION_PATHS: Final = 1024
TA_EVALUATION_STEPS: Final = 64
TA_INFERENCE_REPLICATES: Final = 2

# PINN stopping uses a fixed, truth-free validation bank.  These values are
# frozen before formal training and are shared across C0/C1/C2 and both PDEs.
PINN_VALIDATION_STATES: Final = 256
PINN_VALIDATION_HTE_VECTORS: Final = 16
PINN_VALIDATION_SEED: Final = 20270117
PINN_CHECK_INTERVAL_EPOCHS: Final = 100

PLATEAU_MIN_SECONDS: Final = 15 * 60
PLATEAU_WINDOW: Final = 3
PLATEAU_RELATIVE_IMPROVEMENT: Final = 0.01

EXPECTED_INPUTS: Final = {
    "b6_fullrank": {
        "test_input_sha256": "0ee05fff7ce16f30d0d31563f60df7db3eb23a9c46bbe7bc9be8eb3dbee1311b",
        "test_bank_file_sha256": "6ee6c56ca3d4f413764f28865cb8c00823ea2881ebc90f2643ff4a36fe878730",
        "truth_value_sha256": "6262db7f6b7ee2c41cb6011dadb7288226157533488f708be7d6ea6773f2e13a",
        "parameter_file_sha256": "21b508c30ad63de912979cc9c85ecda0d81aedab1357c88a40e04692495972a4",
    },
    "dense_fourier": {
        "test_input_sha256": "d409296128f48c5ae00c5de2c6ddd998768c01dff18a89c385245a50ff581b70",
        "test_bank_file_sha256": "998fde8b2c0556de4b8fa9a55dec5e489224828cf74e5f3b36f77cb1a8e2c011",
        "truth_value_sha256": "3b33ea36e21595b5c20e5d91097d1fadfcbf3f7884c0328f20579e69995fc816",
        "parameter_file_sha256": "042588d7c3cc8175cc7b9344284725443650837f72008806c668d35b9a005307",
    },
}

CONFIGURATION_GRID: Final = {
    "tangent_actor": {
        "C0": {"paths_per_label": 128, "optimizer_steps": 5000, "actor_updates": 4,
               "continuation_actor_updates_max": 8},
        "C1": {"paths_per_label": 256, "optimizer_steps": 5000, "actor_updates": 4,
               "continuation_actor_updates_max": 4},
        "C2": {"paths_per_label": 128, "optimizer_steps": 10000, "actor_updates": 4,
               "continuation_actor_updates_max": 4},
    },
    "dpi": {
        "C0": {"epochs_per_picard": 16, "integral_mc_samples": 4096,
               "original_picard_iterations": 80, "picard_iterations_max": 320},
        "C1": {"epochs_per_picard": 32, "integral_mc_samples": 4096,
               "original_picard_iterations": 80, "picard_iterations_max": 80},
        "C2": {"epochs_per_picard": 16, "integral_mc_samples": 8192,
               "original_picard_iterations": 80, "picard_iterations_max": 80},
    },
    "pinn_hte": {
        "C0": {"hutchinson_vectors": 16, "batch_size": 512, "original_epochs": 30000},
        "C1": {"hutchinson_vectors": 32, "batch_size": 512, "original_epochs": 30000},
        "C2": {"hutchinson_vectors": 16, "batch_size": 1024, "original_epochs": 30000},
    },
}

PACKAGE_DIR = Path(__file__).resolve().parent
SOURCE_ROOT = PACKAGE_DIR.parents[2]
RESULTS_DIR = SOURCE_ROOT / "invariant_region_mlp" / "results" / "whole_function_compute_frontier"
RAW_DIR = RESULTS_DIR / "raw"
CHECKPOINT_DIR = RESULTS_DIR / "checkpoints"
PREDICTION_DIR = RESULTS_DIR / "predictions"
FIGURE_DIR = RESULTS_DIR / "figures"
LOG_DIR = RESULTS_DIR / "slurm"
STATUS_DIR = RESULTS_DIR / "status"
WORK_DIR = RESULTS_DIR / "work"
REPORT_PATH = RESULTS_DIR / "experiment_report.md"


def protocol_payload() -> dict[str, Any]:
    payload: dict[str, Any] = {
        "revision": REVISION,
        "branch": BRANCH,
        "base_revision": BASE_REVISION,
        "upstream_dpi_revision": UPSTREAM_DPI_REVISION,
        "pdes": list(PDES),
        "methods": list(METHODS),
        "configurations": list(CONFIGURATIONS),
        "dimension": DIMENSION,
        "training_seed": TRAINING_SEED,
        "test_states": TEST_STATES,
        "target_minutes": list(TARGET_MINUTES),
        "primary_minutes": list(PRIMARY_MINUTES),
        "active_training_cap_seconds": ACTIVE_TRAINING_CAP_SECONDS,
        "global_gpu_seconds_cap": GLOBAL_GPU_SECONDS_CAP,
        "historical_C0_consistency": {
            "rule": "abs(new_rmse-old_rmse) <= max(0.05, 0.50*old_rmse)",
            "absolute_rmse_tolerance": C0_CONSISTENCY_ABSOLUTE_RMSE,
            "relative_rmse_tolerance": C0_CONSISTENCY_RELATIVE_RMSE,
            "role": "fidelity audit only; never used for selection or stopping",
        },
        "configuration_grid": CONFIGURATION_GRID,
        "expected_inputs": EXPECTED_INPUTS,
        "ta_evaluation": {
            "paths": TA_EVALUATION_PATHS,
            "euler_steps": TA_EVALUATION_STEPS,
            "independent_replicates": TA_INFERENCE_REPLICATES,
        },
        "pinn_validation": {
            "states": PINN_VALIDATION_STATES,
            "hutchinson_vectors": PINN_VALIDATION_HTE_VECTORS,
            "seed": PINN_VALIDATION_SEED,
            "check_interval_epochs": PINN_CHECK_INTERVAL_EPOCHS,
            "metric": "sqrt(mean(PDE_residual^2) + 10*mean(terminal_error^2))",
            "uses_nonterminal_truth": False,
        },
        "plateau_rule": {
            "minimum_training_seconds": PLATEAU_MIN_SECONDS,
            "window": PLATEAU_WINDOW,
            "maximum_relative_improvement": PLATEAU_RELATIVE_IMPROVEMENT,
            "tangent_actor_metric": "independent noise-corrected fixed-point inconsistency",
            "pinn_hte_metric": "independent fixed-bank PDE residual plus terminal objective",
            "dpi": "not applied; no source-compatible independent validation metric",
        },
        "timing": {
            "main_axis": "cumulative synchronized active training plus causal selection",
            "training_includes": {
                "tangent_actor": "label generation and actor fitting",
                "dpi": "Picard data generation and network fitting",
                "pinn_hte": "HTE residual, terminal loss, and optimizer update",
            },
            "separate_fields": [
                "training_seconds", "selection_seconds", "inference_seconds",
                "allocated_gpu_seconds", "job_wall_seconds",
            ],
        },
        "truth_policy": (
            "frozen nonterminal truth is opened only after all 18 trajectories have stopped; "
            "it is forbidden from training, configuration choice, selection, stopping, and extension"
        ),
        "sdgd_policy": "no SDGD-PINN job may be submitted",
        "historical_reuse": {
            "tangent_actor_C0": "reuse audited actor checkpoints and phase timings; continuation is actor-boundary exact",
            "dpi_C0": "rerun because the historical artifact lacks complete continuation/timing state",
            "pinn_hte_C0": "rerun because the historical artifact lacks optimizer/RNG continuation state",
        },
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    payload["protocol_sha256"] = sha256(canonical.encode()).hexdigest()
    return payload


def ensure_directories() -> None:
    for path in (
        RESULTS_DIR, RAW_DIR, CHECKPOINT_DIR, PREDICTION_DIR, FIGURE_DIR,
        LOG_DIR, STATUS_DIR, WORK_DIR,
    ):
        path.mkdir(parents=True, exist_ok=True)


__all__ = [name for name in globals() if name.isupper()] + ["ensure_directories", "protocol_payload"]

