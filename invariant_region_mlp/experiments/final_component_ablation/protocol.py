"""Immutable declarations for the final component ablation.

This module contains declarations only.  Formal runners must refuse settings
that differ from this protocol instead of accepting ad-hoc command-line
overrides.
"""

from __future__ import annotations

from dataclasses import asdict
import hashlib
import json
from pathlib import Path
from typing import Any

from invariant_region_mlp.experiments.original_pde_benchmark_comparison.tangent_solver import (
    ActorConfig,
)
from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis.training import (
    FitConfig,
)


PACKAGE_DIR = Path(__file__).resolve().parent
PROJECT_DIR = PACKAGE_DIR.parents[1]
RESULTS_DIR = PROJECT_DIR / "results" / "final_component_ablation"
RAW_DIR = RESULTS_DIR / "raw"
CACHE_DIR = RESULTS_DIR / "cache"
CHECKPOINT_DIR = RESULTS_DIR / "checkpoints"
PREDICTION_DIR = RESULTS_DIR / "predictions"
INPUT_DIR = RESULTS_DIR / "inputs"
FIGURE_DIR = RESULTS_DIR / "figures"
TABLE_DIR = RESULTS_DIR / "tables"
SLURM_DIR = RESULTS_DIR / "slurm"

BRANCH = "ir-mlp-final-component-ablation"
DIMENSION = 100
PDES = ("hjb_log", "joint_concave", "joint_convex", "nonconvex_coupled")
VARIANTS = ("a0_full", "a1_no_mcv", "a2_no_terminal_anchor", "a3_no_updates")
TRAINING_SEEDS = (0, 1, 2)
UPDATES = 4
EVALUATION_STATES = 960
EVALUATION_STEPS = 64
EVALUATION_PATHS = 1024
EVALUATION_REPLICATES = (0, 1)
ALLOCATED_GPU_HOUR_CEILING = 36.0

# HJB preserves the published primary configuration's mixed state distribution.
HJB_CONFIG = ActorConfig(
    training_points=8192,
    validation_points=512,
    paths_per_label=128,
    label_time_steps=32,
    optimizer_steps=5000,
    width=384,
    depth=4,
    activation="silu",
    learning_rate=1.5e-4,
    batch_size=512,
    parent_batch=8,
    path_chunk=32,
    checkpoint_interval=250,
    gradient_clip=10.0,
    training_distribution="mixed",
    label_estimator="mcv",
)

JOINT_CONFIG = FitConfig(
    training_points=8192,
    validation_points=512,
    paths_per_label=128,
    label_time_steps=32,
    optimizer_steps=5000,
    width=384,
    depth=4,
    learning_rate=1.5e-4,
    batch_size=512,
    parent_batch=8,
    path_chunk=32,
    checkpoint_interval=250,
    gradient_clip=10.0,
)


def stable_seed(*parts: object) -> int:
    payload = "::".join(str(part) for part in parts).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little") % (2**31 - 1)


def ensure_dirs() -> None:
    for path in (
        RESULTS_DIR,
        RAW_DIR,
        CACHE_DIR,
        CHECKPOINT_DIR,
        PREDICTION_DIR,
        INPUT_DIR,
        FIGURE_DIR,
        TABLE_DIR,
        SLURM_DIR,
    ):
        path.mkdir(parents=True, exist_ok=True)


def frozen_protocol() -> dict[str, Any]:
    body: dict[str, Any] = {
        "schema_version": 1,
        "branch": BRANCH,
        "dimension": DIMENSION,
        "pdes": list(PDES),
        "variants": list(VARIANTS),
        "training_seeds": list(TRAINING_SEEDS),
        "actor_updates": UPDATES,
        "hjb_configuration": asdict(HJB_CONFIG),
        "joint_configuration": asdict(JOINT_CONFIG),
        "evaluation": {
            "states": EVALUATION_STATES,
            "steps": EVALUATION_STEPS,
            "paths": EVALUATION_PATHS,
            "replicates": list(EVALUATION_REPLICATES),
            "primary_replicate": 0,
            "metric": "whole-function value RMSE",
            "rmae": "sum(abs(prediction-reference))/sum(abs(reference))",
        },
        "estimators": {
            "a0_full": {"training_labels": "mcv", "selector": "mcv", "inference": "mcv"},
            "a1_no_mcv": {"training_labels": "plain", "selector": "plain", "inference": "plain"},
            "a2_no_terminal_anchor": {"training_labels": "mcv", "selector": "mcv", "inference": "mcv"},
            "a3_no_updates": {"training_labels": None, "selector": None, "inference": "mcv"},
        },
        "selection": {
            "joint": "minimum independent noise-corrected normalized fixed-point inconsistency",
            "hjb": "published four-update protocol with held-out noisy-label optimizer checkpoints",
            "ground_truth_used": False,
        },
        "nonconvex_960_bank_amendment": {
            "reason": "historical nonconvex manuscript bank has 512 states but this task mandates 960",
            "rule": "icml_final_test_bank/nonconvex_coupled/d100, 960 states, T=0.5",
            "status": "frozen_before_evaluation",
            "historical_512_results_reused_as_primary": False,
        },
        "allocated_gpu_hour_ceiling": ALLOCATED_GPU_HOUR_CEILING,
        "no_test_reference_for_training_selection_or_tuning": True,
        "formal_timing_gpu": "NVIDIA GH200",
    }
    encoded = json.dumps(body, sort_keys=True, separators=(",", ":")).encode("utf-8")
    body["protocol_sha256"] = hashlib.sha256(encoded).hexdigest()
    return body


__all__ = [
    "ALLOCATED_GPU_HOUR_CEILING",
    "BRANCH",
    "CACHE_DIR",
    "CHECKPOINT_DIR",
    "DIMENSION",
    "EVALUATION_PATHS",
    "EVALUATION_REPLICATES",
    "EVALUATION_STATES",
    "EVALUATION_STEPS",
    "FIGURE_DIR",
    "HJB_CONFIG",
    "INPUT_DIR",
    "JOINT_CONFIG",
    "PACKAGE_DIR",
    "PDES",
    "PREDICTION_DIR",
    "PROJECT_DIR",
    "RAW_DIR",
    "RESULTS_DIR",
    "SLURM_DIR",
    "TABLE_DIR",
    "TRAINING_SEEDS",
    "UPDATES",
    "VARIANTS",
    "ensure_dirs",
    "frozen_protocol",
    "stable_seed",
]
