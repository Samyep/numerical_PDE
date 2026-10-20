"""Frozen constants and artifact routing for the high-dimensional pilot."""

from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path
from typing import Any, Final

from invariant_region_mlp.experiments.double_estimator.equations import (
    BASE_SEED,
    PACKAGE_ROOT,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.tangent_solver import (
    ActorConfig,
)
from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis.training import (
    FitConfig,
)


REVISION: Final = "main_highdim_pilot_5070_v1"
DIMENSION: Final = 100
TERMINAL_TIME: Final = 0.5
TRAINING_SEED: Final = 0
UPDATES: Final = 2

# The preregistration fixes the distributions but intentionally leaves the
# concrete random seeds to the implementation.  These tuples are frozen here
# before any parameter draw, calibration, or solver run.
DENSE_PARAMETER_SEED: Final = (BASE_SEED, DIMENSION, 41)
LQG_COEFFICIENT_SEED: Final = (BASE_SEED, DIMENSION, 73)
STATE_SEED_NAMESPACE: Final = "main_highdim_pilot_5070_v1"

TEST_STATES: Final = 512
VALIDATION_STATES: Final = 256
SELECTOR_STATES: Final = 64
HEAT_MC_PATHS: Final = 16384
HEAT_MC_REPLICATES: Final = 2
INFERENCE_PATHS: Final = 1024
INFERENCE_STEPS: Final = 64
INFERENCE_REPLICATES: Final = 2
REFERENCE_QUADRATURE_ORDER: Final = 192
REFERENCE_AUDIT_ORDERS: Final = (96, 144, 192, 256)

JOINT_CONFIG = FitConfig(
    training_points=2048,
    validation_points=256,
    paths_per_label=64,
    label_time_steps=16,
    optimizer_steps=1500,
    width=256,
    depth=3,
    learning_rate=1.5e-4,
    batch_size=256,
    parent_batch=4,
    path_chunk=32,
    checkpoint_interval=150,
    gradient_clip=10.0,
)

GRADIENT_CONFIG = ActorConfig(
    training_points=2048,
    validation_points=256,
    paths_per_label=64,
    label_time_steps=16,
    optimizer_steps=1500,
    width=256,
    depth=3,
    activation="silu",
    learning_rate=1.5e-4,
    batch_size=256,
    parent_batch=4,
    path_chunk=32,
    checkpoint_interval=150,
    gradient_clip=10.0,
    training_distribution="frozen_pde_specific",
    label_estimator="mcv",
)

RESULTS_ROOT = PACKAGE_ROOT / "results" / "main_highdim_pilot"
INPUT_ROOT = RESULTS_ROOT / "inputs"
RAW_ROOT = RESULTS_ROOT / "raw"
CACHE_ROOT = RESULTS_ROOT / "cache"
CHECKPOINT_ROOT = RESULTS_ROOT / "checkpoints"
PREDICTION_ROOT = RESULTS_ROOT / "predictions"
FIGURES_ROOT = RESULTS_ROOT / "figures"
REPORT_PATH = PACKAGE_ROOT / "docs" / "MAIN_HIGHDIM_PILOT_5070_REPORT.md"
VALIDATION_PATH = RESULTS_ROOT / "validation.json"
MANIFEST_PATH = RESULTS_ROOT / "test_bank_manifest.json"
CALIBRATION_PATH = RESULTS_ROOT / "benchmark_calibration.csv"


def ensure_directories() -> None:
    for path in (
        RESULTS_ROOT,
        INPUT_ROOT,
        RAW_ROOT,
        CACHE_ROOT,
        CHECKPOINT_ROOT,
        PREDICTION_ROOT,
        FIGURES_ROOT,
    ):
        path.mkdir(parents=True, exist_ok=True)


def write_json_atomic(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def frozen_configuration() -> dict[str, Any]:
    return {
        "revision": REVISION,
        "base_seed": BASE_SEED,
        "dimension": DIMENSION,
        "terminal_time": TERMINAL_TIME,
        "training_seed": TRAINING_SEED,
        "updates": UPDATES,
        "dense_parameter_seed": list(DENSE_PARAMETER_SEED),
        "lqg_coefficient_seed": list(LQG_COEFFICIENT_SEED),
        "state_seed_namespace": STATE_SEED_NAMESPACE,
        "joint_actor": asdict(JOINT_CONFIG),
        "gradient_actor": asdict(GRADIENT_CONFIG),
        "test_states": TEST_STATES,
        "validation_states": VALIDATION_STATES,
        "selector_states": SELECTOR_STATES,
        "heat_mc_paths": HEAT_MC_PATHS,
        "heat_mc_replicates": HEAT_MC_REPLICATES,
        "inference_paths": INFERENCE_PATHS,
        "inference_steps": INFERENCE_STEPS,
        "inference_replicates": INFERENCE_REPLICATES,
        "reference_quadrature_order": REFERENCE_QUADRATURE_ORDER,
        "reference_audit_orders": list(REFERENCE_AUDIT_ORDERS),
        "per_new_pde_gpu_seconds_cap": 3600,
        "total_gpu_seconds_cap": 10800,
    }


__all__ = [name for name in globals() if name.isupper()] + [
    "ensure_directories",
    "frozen_configuration",
    "write_json_atomic",
]

