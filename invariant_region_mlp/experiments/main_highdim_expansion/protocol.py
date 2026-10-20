"""Frozen protocol for the two-PDE high-dimensional main-suite expansion.

Definitions in this module are fixed before any production training.  Runtime
artifacts may record progress, but they must never mutate these declarations.
"""

from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Final

from invariant_region_mlp.experiments.double_estimator.equations import BASE_SEED
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.tangent_solver import (
    ActorConfig,
)
from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis.training import (
    FitConfig,
)


REVISION: Final = "main_highdim_expansion_v1"
BRANCH: Final = "ir-mlp-icml-two-highdim-expansion"
DIMENSIONS: Final = (20, 50, 100, 200, 400)
PDES: Final = ("b6_fullrank", "dense_fourier")
METHODS: Final = ("tangent_actor", "dpi", "pinn_hte", "sdgd_pinn")
TRAINING_SEED: Final = 0
TEST_STATES: Final = 960
VALIDATION_STATES: Final = 512
ACTOR_UPDATES: Final = 4
EVALUATION_PATHS: Final = 1024
PRIMARY_EULER_STEPS: Final = 64
AUDIT_EULER_STEPS: Final = (128, 256)
INFERENCE_REPLICATES: Final = 2

CELL_GPU_SECONDS_CAP: Final = 6 * 3600
GLOBAL_GPU_SECONDS_CAP: Final = 120 * 3600

B6_BETAS: Final = (1.0, 1.5, 2.0, 3.0, 4.0)
B6_GAMMAS: Final = (0.0, 1.0)
B6_HEAT_SAMPLES: Final = 4096
B6_FIXED_SELECTIONS: Final = {
    100: {"beta": 1.5, "gamma": 1.0},
    400: {"beta": 2.0, "gamma": 1.0},
}

GRADIENT_ACTOR_CONFIG: Final = ActorConfig(
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
    training_distribution="function",
    label_estimator="mcv",
)

JOINT_ACTOR_CONFIG: Final = FitConfig(
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

SOURCE_COMMITS: Final = {
    "icml_final_suite": "f762109160a5c0bb6bba928370491bc9bbd04232",
    "b6_validated": "281977883e0a5ce4f05895878680b582620461be",
    "dense_fourier_pilot": "3d3d4c80f0b29af39e7d509a143812a1f5a702f0",
    "dpi_pinn_hte": "e8616c2f531c9f95299b246a13c1180b8f02822a",
    "dpi": "e8616c2f531c9f95299b246a13c1180b8f02822a",
    "pinn_hte": "e8616c2f531c9f95299b246a13c1180b8f02822a",
    "sdgd_pinn": "c597cc4fe9caac0b3ca118ea84453405f25ff833",
}


PACKAGE_DIR = Path(__file__).resolve().parent
SOURCE_ROOT = PACKAGE_DIR.parents[2]
RESULTS_DIR = SOURCE_ROOT / "invariant_region_mlp" / "results" / "main_highdim_expansion"
INPUT_DIR = RESULTS_DIR / "inputs"
RAW_DIR = RESULTS_DIR / "raw"
CACHE_DIR = RESULTS_DIR / "cache"
CHECKPOINT_DIR = RESULTS_DIR / "checkpoints"
PREDICTION_DIR = RESULTS_DIR / "predictions"
FIGURE_DIR = RESULTS_DIR / "figures"
LOG_DIR = RESULTS_DIR / "slurm"
STATUS_DIR = RESULTS_DIR / "cell_status"
MANIFEST_PATH = RESULTS_DIR / "parameter_and_test_bank_manifest.json"
CALIBRATION_PATH = RESULTS_DIR / "b6_calibration_freeze.json"
VALIDATION_PATH = RESULTS_DIR / "mathematical_validation.json"
CHARACTERIZATION_PATH = RESULTS_DIR / "benchmark_characterization.csv"
EXECUTION_STATUS_PATH = RESULTS_DIR / "execution_status.json"
REPORT_PATH = SOURCE_ROOT / "invariant_region_mlp" / "docs" / "ICML_MAIN_HIGHDIM_EXPANSION_REPORT.md"


def stable_seed(*parts: object) -> int:
    digest = sha256("::".join(map(str, parts)).encode("utf-8")).digest()
    return int.from_bytes(digest[:4], "little")


def ensure_directories() -> None:
    for path in (
        RESULTS_DIR,
        INPUT_DIR,
        RAW_DIR,
        CACHE_DIR,
        CHECKPOINT_DIR,
        PREDICTION_DIR,
        FIGURE_DIR,
        LOG_DIR,
        STATUS_DIR,
    ):
        path.mkdir(parents=True, exist_ok=True)


def frozen_protocol() -> dict[str, Any]:
    payload: dict[str, Any] = {
        "revision": REVISION,
        "branch": BRANCH,
        "base_seed": BASE_SEED,
        "dimensions": list(DIMENSIONS),
        "pdes": list(PDES),
        "methods": list(METHODS),
        "training_seed": TRAINING_SEED,
        "test_states": TEST_STATES,
        "validation_states": VALIDATION_STATES,
        "actor_updates": ACTOR_UPDATES,
        "evaluation_paths": EVALUATION_PATHS,
        "primary_euler_steps": PRIMARY_EULER_STEPS,
        "audit_euler_steps": list(AUDIT_EULER_STEPS),
        "inference_replicates": INFERENCE_REPLICATES,
        "cell_gpu_seconds_cap": CELL_GPU_SECONDS_CAP,
        "global_gpu_seconds_cap": GLOBAL_GPU_SECONDS_CAP,
        "gradient_actor": asdict(GRADIENT_ACTOR_CONFIG),
        "joint_actor": asdict(JOINT_ACTOR_CONFIG),
        "b6_calibration": {
            "betas": list(B6_BETAS),
            "gammas": list(B6_GAMMAS),
            "heat_samples": B6_HEAT_SAMPLES,
            "fixed_existing_selections": B6_FIXED_SELECTIONS,
            "selection_rule": (
                "qualify iff ESS_med>=5, 0.5<=z_rms<=4, "
                "heat_only_skill>=0.15, f_cv>=0.05; maximize r_glob; "
                "tie-break smaller beta then gamma"
            ),
        },
        "test_distributions": {
            "b6_fullrank": "t~Uniform[0,0.25), x~Uniform([-1,1]^d)",
            "dense_fourier": "t~Uniform[0,0.5), x~Uniform([-1,1]^d)",
        },
        "selection": {
            "tangent_actor": "reference-free held-out noisy tangent labels and fixed-point consistency",
            "external_baselines": "source-faithful fixed protocol; no final-test model selection",
            "formal_seed": 0,
        },
        "source_commits": SOURCE_COMMITS,
        "ground_truth_policy": (
            "exact solutions are evaluation/calibration only and are forbidden "
            "from labels, optimization, checkpoint selection, and stopping"
        ),
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    payload["protocol_sha256"] = sha256(canonical).hexdigest()
    return payload


def write_json_atomic(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


__all__ = [name for name in globals() if name.isupper()] + [
    "ensure_directories",
    "frozen_protocol",
    "stable_seed",
    "write_json_atomic",
]
