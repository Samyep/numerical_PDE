"""Machine-readable constants frozen by :mod:`FROZEN_PROTOCOL.md`."""

from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Final

import numpy as np

from invariant_region_mlp.experiments.original_pde_benchmark_comparison.tangent_solver import (
    ActorConfig,
)
from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis.training import (
    FitConfig,
)


REVISION: Final = "single_point_32_ice_v1"
BRANCH: Final = "ir-mlp-icml-single-point-32-ice"
BASE_COMMIT: Final = "027240a72fe3028b4308275e5ec8f06722b79c33"
DIMENSIONS: Final = (20, 50, 100, 200, 400)
PDES: Final = ("hjb_log", "allen_cahn", "b6_fullrank", "dense_fourier")
METHODS: Final = ("tangent_actor", "fh_mlp_standard", "deep_bsde")
QUERY_COUNT: Final = 32
TRAINING_SEEDS: Final = (0, 1, 2)
ACTOR_UPDATES: Final = 4
TA_EVALUATION_PATHS: Final = 1024
TA_PRIMARY_EULER_STEPS: Final = 64
TA_AUDIT_EULER_STEPS: Final = (128, 256)
TA_INFERENCE_REPLICATES: Final = 4
FH_CONFIGURATIONS: Final = ((2, 16), (2, 32), (3, 8), (3, 12), (3, 16), (3, 24))
FH_MAIN_CONFIGURATION: Final = (3, 16)
FH_REPLICATES: Final = 10
DEEP_BSDE_INTERVALS: Final = 20
DEEP_BSDE_BATCH_SIZE: Final = 64
DEEP_BSDE_VALIDATION_SIZE: Final = 256
DEEP_BSDE_CHECKPOINT_INTERVAL: Final = 100

QUERY_SEEDS: Final = {
    "hjb_log": {20: 1014295083, 50: 1870289811, 100: 1898610954,
                200: 2768728099, 400: 1832156210},
    "allen_cahn": {20: 1716816384, 50: 259137838, 100: 4074736349,
                   200: 1881011736, 400: 255990953},
    "b6_fullrank": {20: 214171249, 50: 697137073, 100: 1737294624,
                    200: 2575111948, 400: 3771049365},
    "dense_fourier": {20: 3392201883, 50: 1193042412, 100: 2884838781,
                      200: 2659431886, 400: 2229501467},
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

B6_SELECTIONS: Final = {
    20: {"beta": 1.0, "gamma": 1.0},
    50: {"beta": 1.5, "gamma": 1.0},
    100: {"beta": 1.5, "gamma": 1.0},
    200: {"beta": 2.0, "gamma": 1.0},
    400: {"beta": 2.0, "gamma": 1.0},
}

SOURCE_COMMITS: Final = {
    "base_and_highdim_results": BASE_COMMIT,
    "mature_whole_function_suite": "136e892ba6ffd5fcf33d5d3dc618828584675f44",
    "validated_b6": "281977883e0a5ce4f05895878680b582620461be",
    "dense_fourier_pilot": "3d3d4c80f0b29af39e7d509a143812a1f5a702f0",
    "strong_baselines": "66495efb3510ab99514942ed7f9ee62259b2d21d",
}

PACKAGE_DIR = Path(__file__).resolve().parent
SOURCE_ROOT = PACKAGE_DIR.parents[2]
RESULTS_DIR = SOURCE_ROOT / "invariant_region_mlp" / "results" / "single_point_32_ice"
SOURCE_EXPANSION_RESULTS = (
    SOURCE_ROOT / "invariant_region_mlp" / "results" / "main_highdim_expansion"
)
QUERY_DIR = RESULTS_DIR / "frozen_query_banks"
PARAMETER_DIR = RESULTS_DIR / "parameters"
PREDICTION_DIR = RESULTS_DIR / "predictions"
CHECKPOINT_DIR = RESULTS_DIR / "checkpoints"
CACHE_DIR = RESULTS_DIR / "cache"
RAW_DIR = RESULTS_DIR / "raw"
STATUS_DIR = RESULTS_DIR / "task_status"
LOG_DIR = RESULTS_DIR / "slurm"
FIGURE_DIR = RESULTS_DIR / "figures"
TABLE_DIR = RESULTS_DIR / "tables"
QUEUE_DIR = RESULTS_DIR / "dispatch_queues"
FROZEN_PROTOCOL_JSON = RESULTS_DIR / "FROZEN_PROTOCOL.json"
PARAMETER_MANIFEST = RESULTS_DIR / "parameter_manifest.json"
REFERENCE_VALIDATION = RESULTS_DIR / "reference_validation.json"
EXECUTION_STATUS = RESULTS_DIR / "execution_status.json"
STAGE_A_AUDIT = RESULTS_DIR / "stage_a_audit.json"
STAGE_B_AUDIT = RESULTS_DIR / "stage_b_audit.json"
METHOD_FIDELITY_AUDIT = RESULTS_DIR / "method_fidelity_audit.md"


def stable_seed(*parts: object) -> int:
    """Return a revision-namespaced deterministic uint32 seed."""
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
        RESULTS_DIR, QUERY_DIR, PARAMETER_DIR, PREDICTION_DIR, CHECKPOINT_DIR,
        CACHE_DIR, RAW_DIR, STATUS_DIR, LOG_DIR, FIGURE_DIR, TABLE_DIR,
        QUEUE_DIR,
    ):
        path.mkdir(parents=True, exist_ok=True)


def query_path(pde: str, dimension: int) -> Path:
    return QUERY_DIR / f"{pde}_d{int(dimension)}"


def load_query_bank(pde: str, dimension: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    root = query_path(pde, dimension)
    return (
        np.load(root / "t.npy", allow_pickle=False),
        np.load(root / "x.npy", allow_pickle=False),
        np.load(root / "reference.npy", allow_pickle=False),
    )


def require_stage_a() -> dict[str, Any]:
    if not STAGE_A_AUDIT.is_file():
        raise RuntimeError(f"missing Stage-A audit: {STAGE_A_AUDIT}")
    payload = json.loads(STAGE_A_AUDIT.read_text(encoding="utf-8"))
    if payload.get("status") != "VALIDATED" or payload.get("overall_pass") is not True:
        raise RuntimeError("Stage-A CUDA mathematical preflight has not passed")
    return payload


def require_stage_b() -> dict[str, Any]:
    """Require the nonformal 100D end-to-end pilot before formal dispatch."""
    require_stage_a()
    if not STAGE_B_AUDIT.is_file():
        raise RuntimeError(f"missing Stage-B audit: {STAGE_B_AUDIT}")
    payload = json.loads(STAGE_B_AUDIT.read_text(encoding="utf-8"))
    if payload.get("status") != "VALIDATED" or payload.get("overall_pass") is not True:
        raise RuntimeError("Stage-B representative pilot has not passed")
    return payload


def deep_bsde_configuration(pde: str, dimension: int, *, pilot: bool = False) -> dict[str, Any]:
    allen = pde == "allen_cahn"
    horizon = {"hjb_log": 1.0, "allen_cahn": 0.3,
               "b6_fullrank": 0.25, "dense_fourier": 0.5}[pde]
    return {
        "benchmark": f"{pde}_d{int(dimension)}",
        "total_time": horizon,
        "intervals": DEEP_BSDE_INTERVALS,
        "hidden_widths": (int(dimension) + 10, int(dimension) + 10),
        "learning_rate": 5e-4 if allen else 1e-2,
        "iterations": 4 if pilot else (4000 if allen else 2000),
        "batch_size": 8 if pilot else DEEP_BSDE_BATCH_SIZE,
        "validation_size": 16 if pilot else DEEP_BSDE_VALIDATION_SIZE,
        "logging_frequency": 1 if pilot else DEEP_BSDE_CHECKPOINT_INTERVAL,
        "y_init_range": (0.3, 0.5) if allen else (0.0, 1.0),
        "dtype": "float64",
    }


def frozen_protocol() -> dict[str, Any]:
    payload: dict[str, Any] = {
        "revision": REVISION,
        "branch": BRANCH,
        "base_commit": BASE_COMMIT,
        "dimensions": list(DIMENSIONS),
        "pdes": list(PDES),
        "methods": list(METHODS),
        "query_count": QUERY_COUNT,
        "query_seeds": QUERY_SEEDS,
        "training_seeds": list(TRAINING_SEEDS),
        "actor_updates": ACTOR_UPDATES,
        "gradient_actor": asdict(GRADIENT_ACTOR_CONFIG),
        "joint_actor": asdict(JOINT_ACTOR_CONFIG),
        "ta_evaluation": {
            "paths": TA_EVALUATION_PATHS,
            "primary_euler_steps": TA_PRIMARY_EULER_STEPS,
            "audit_euler_steps": list(TA_AUDIT_EULER_STEPS),
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
        "b6_selections": B6_SELECTIONS,
        "source_commits": SOURCE_COMMITS,
        "ground_truth_policy": (
            "evaluation and reference validation only; forbidden from training, "
            "checkpoint selection, stopping, and solver-budget selection"
        ),
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    normalized = json.loads(canonical)
    normalized["protocol_sha256"] = sha256(canonical.encode("utf-8")).hexdigest()
    return normalized


__all__ = [name for name in globals() if name.isupper()] + [
    "array_sha256", "deep_bsde_configuration", "ensure_directories",
    "file_sha256", "frozen_protocol", "load_query_bank", "query_path", "require_stage_a",
    "require_stage_b",
    "stable_seed", "write_json_atomic", "write_npz_atomic",
]

