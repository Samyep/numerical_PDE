"""Frozen constants and artifact helpers for the joint-PDE extension."""

from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Final

import numpy as np

from invariant_region_mlp.experiments.single_point_32_ice.protocol import (
    JOINT_ACTOR_CONFIG,
)


REVISION: Final = "single_point_joint_nonseparable_v1"
BRANCH: Final = "ir-mlp-icml-joint-pde"
BASE_COMMIT: Final = "32ef8b074"
PREREGISTRATION_COMMIT: Final = "0f1485703"
BASE_SEED: Final = 20261232
PDE: Final = "joint_nonseparable"
DIMENSIONS: Final = (20, 50, 100, 200, 400)
METHODS: Final = ("tangent_actor", "fh_mlp_standard", "deep_bsde")
QUERY_COUNT: Final = 32
TRAINING_SEEDS: Final = (0, 1, 2)
PRIMARY_TA_TRAINING_SEEDS: Final = (0,)
ACTOR_UPDATES: Final = 4
TA_EVALUATION_PATHS: Final = 1024
TA_PRIMARY_EULER_STEPS: Final = 64
TA_INFERENCE_REPLICATES: Final = 4
FH_CONFIGURATIONS: Final = (
    (2, 16), (2, 32), (3, 8), (3, 12), (3, 16), (3, 24),
)
FH_MAIN_CONFIGURATION: Final = (3, 16)
FH_REPLICATES: Final = 10
DEEP_BSDE_INTERVALS: Final = 20
DEEP_BSDE_BATCH_SIZE: Final = 64
DEEP_BSDE_VALIDATION_SIZE: Final = 256
DEEP_BSDE_CHECKPOINT_INTERVAL: Final = 100

PACKAGE_DIR = Path(__file__).resolve().parent
SOURCE_ROOT = PACKAGE_DIR.parents[2]
RESULTS_DIR = SOURCE_ROOT / "invariant_region_mlp" / "results" / "single_point_joint_extension"
PARAMETER_DIR = RESULTS_DIR / "parameters"
QUERY_DIR = RESULTS_DIR / "frozen_query_banks"
PREDICTION_DIR = RESULTS_DIR / "predictions"
CHECKPOINT_DIR = RESULTS_DIR / "checkpoints"
CACHE_DIR = RESULTS_DIR / "cache"
RAW_DIR = RESULTS_DIR / "raw"
STATUS_DIR = RESULTS_DIR / "task_status"
LOG_DIR = RESULTS_DIR / "slurm"
FIGURE_DIR = RESULTS_DIR / "figures"
TABLE_DIR = RESULTS_DIR / "tables"
FROZEN_PROTOCOL_JSON = RESULTS_DIR / "FROZEN_PROTOCOL.json"
PARAMETER_MANIFEST = RESULTS_DIR / "parameter_manifest.json"
REFERENCE_VALIDATION = RESULTS_DIR / "reference_validation.json"
PREFLIGHT_AUDIT = RESULTS_DIR / "preflight_audit.json"
SMOKE_AUDIT = RESULTS_DIR / "smoke_audit.json"
EXECUTION_STATUS = RESULTS_DIR / "execution_status.json"


def stable_seed(*parts: object) -> int:
    payload = "::".join((REVISION, *map(str, parts))).encode("utf-8")
    return int.from_bytes(sha256(payload).digest()[:4], "little")


def array_sha256(*arrays: np.ndarray) -> str:
    digest = sha256()
    for array in arrays:
        value = np.ascontiguousarray(array)
        digest.update(str(value.dtype).encode("ascii"))
        digest.update(str(value.shape).encode("ascii"))
        digest.update(value.tobytes())
    return digest.hexdigest()


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json_atomic(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def write_npz_atomic(path: Path, **arrays: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp.npz")
    np.savez_compressed(temporary, **arrays)
    temporary.replace(path)


def ensure_directories() -> None:
    for path in (
        RESULTS_DIR,
        PARAMETER_DIR,
        QUERY_DIR,
        PREDICTION_DIR,
        CHECKPOINT_DIR,
        CACHE_DIR,
        RAW_DIR,
        STATUS_DIR,
        LOG_DIR,
        FIGURE_DIR,
        TABLE_DIR,
    ):
        path.mkdir(parents=True, exist_ok=True)


def query_path(dimension: int) -> Path:
    return QUERY_DIR / f"{PDE}_d{int(dimension)}"


def load_query_bank(
    dimension: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    root = query_path(dimension)
    return (
        np.load(root / "t.npy", allow_pickle=False),
        np.load(root / "x.npy", allow_pickle=False),
        np.load(root / "reference.npy", allow_pickle=False),
        np.load(root / "reference_z.npy", allow_pickle=False),
    )


def _require(path: Path, label: str) -> dict[str, Any]:
    if not path.is_file():
        raise RuntimeError(f"missing {label}: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("status") != "VALIDATED" or payload.get("overall_pass") is not True:
        raise RuntimeError(f"{label} has not passed")
    return payload


def require_preflight() -> dict[str, Any]:
    return _require(PREFLIGHT_AUDIT, "joint-extension preflight")


def require_smoke() -> dict[str, Any]:
    require_preflight()
    return _require(SMOKE_AUDIT, "joint-extension smoke gate")


def deep_bsde_configuration(dimension: int, *, pilot: bool = False) -> dict[str, Any]:
    d = int(dimension)
    return {
        "benchmark": f"{PDE}_d{d}",
        "total_time": 0.5,
        "intervals": DEEP_BSDE_INTERVALS,
        "hidden_widths": (d + 10, d + 10),
        "learning_rate": 1e-2,
        "iterations": 4 if pilot else 2000,
        "batch_size": 8 if pilot else DEEP_BSDE_BATCH_SIZE,
        "validation_size": 16 if pilot else DEEP_BSDE_VALIDATION_SIZE,
        "logging_frequency": 1 if pilot else DEEP_BSDE_CHECKPOINT_INTERVAL,
        "y_init_range": (0.0, 1.0),
        "dtype": "float64",
    }


def frozen_protocol() -> dict[str, Any]:
    payload: dict[str, Any] = {
        "revision": REVISION,
        "branch": BRANCH,
        "base_commit": BASE_COMMIT,
        "preregistration_commit": PREREGISTRATION_COMMIT,
        "base_seed": BASE_SEED,
        "pde": PDE,
        "dimensions": list(DIMENSIONS),
        "methods": list(METHODS),
        "query_count": QUERY_COUNT,
        "query_construction": {
            "origin_row": 0,
            "remaining_law": "N(0,I_d)",
            "seed_sequence": "[20261232, d, 32]",
        },
        "parameter_construction": {
            "law": "sign-canonicalized thin QR of d x 2 N(0,1)",
            "seed_sequence": "[20261232, d, 17]",
            "v": "beta_1",
        },
        "training_seeds": list(TRAINING_SEEDS),
        "primary_ta_training_seeds": list(PRIMARY_TA_TRAINING_SEEDS),
        "actor_updates": ACTOR_UPDATES,
        "joint_actor": asdict(JOINT_ACTOR_CONFIG),
        "ta_evaluation": {
            "paths": TA_EVALUATION_PATHS,
            "primary_euler_steps": TA_PRIMARY_EULER_STEPS,
            "inference_replicates": TA_INFERENCE_REPLICATES,
        },
        "fh_mlp": {
            "configurations": [list(item) for item in FH_CONFIGURATIONS],
            "main_configuration": list(FH_MAIN_CONFIGURATION),
            "replicates": FH_REPLICATES,
        },
        "deep_bsde": {
            "intervals": DEEP_BSDE_INTERVALS,
            "batch_size": DEEP_BSDE_BATCH_SIZE,
            "validation_size": DEEP_BSDE_VALIDATION_SIZE,
            "checkpoint_interval": DEEP_BSDE_CHECKPOINT_INTERVAL,
        },
        "reference": {
            "primary_grid": 256,
            "primary_dt_max": 1.25e-4,
            "audit_grid": 128,
            "audit_dt_max": 2.5e-4,
            "value_refinement_rmse_max": 5e-5,
            "value_refinement_max_abs_max": 2e-4,
            "pde_residual_rmse_max": 2e-6,
            "terminal_interpolation_max_abs_max": 1e-4,
        },
        "judgments": {
            "J-1": "TA lowest RMSE in at least 4 of 5 dimensions",
            "J-2": "TA <= 0.5 * FH(n=3,M=16) in at least 4 of 5 dimensions",
            "J-3": "TA <= 0.5 * Deep BSDE in at least 3 of 5 dimensions",
            "J-4": "TA RMSE(d=400)/RMSE(d=20) <= 3.0",
        },
        "ground_truth_policy": (
            "evaluation and validation only; forbidden from labels, fitting, "
            "stopping, checkpoint selection, and method-budget selection"
        ),
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    normalized = json.loads(canonical)
    normalized["protocol_sha256"] = sha256(canonical.encode("utf-8")).hexdigest()
    return normalized


__all__ = [name for name in globals() if name.isupper()] + [
    "array_sha256",
    "deep_bsde_configuration",
    "ensure_directories",
    "file_sha256",
    "frozen_protocol",
    "load_query_bank",
    "query_path",
    "require_preflight",
    "require_smoke",
    "stable_seed",
    "write_json_atomic",
    "write_npz_atomic",
]

