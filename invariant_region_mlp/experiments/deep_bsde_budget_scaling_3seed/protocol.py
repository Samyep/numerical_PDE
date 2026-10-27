"""Frozen declarations and artifact helpers for the Deep-BSDE study."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from hashlib import sha256
import json
import os
from pathlib import Path
import tempfile
from typing import Any, Final

import numpy as np


REVISION: Final = "deep_bsde_budget_scaling_3seed_v1"
BRANCH: Final = "ir-mlp-deep-bsde-budget-scaling-3seed"
PREREGISTRATION_COMMIT: Final = "a1ff02e15f70c6c04a6450ca15c839ee1ea0ca3f"
BASE_REVISION: Final = "9cb0d6f1d69c47882c497807283700fc60c770de"
OFFICIAL_PORT_SOURCE_REVISION: Final = "29518ced89b3a8c0e89cdfa50b2a1cec999904aa"
AUDITED_QUERY_SOLVER_BLOB: Final = "bb2179ffdcd81282829f30569e545e0252c4b5d3"

PDES: Final = ("b6_fullrank", "dense_fourier")
DIMENSION: Final = 100
TRAINING_SEEDS: Final = (0, 1, 2)
CALIBRATION_COUNT: Final = 32
FINAL_COUNT: Final = 64
SOURCE_BANK_COUNT: Final = 960
SUBSET_SEED: Final = 20261309
SUBSET_STRATA: Final = 32
MAX_H100_GPU_HOURS: Final = 120.0
PREFLIGHT_UPDATES: Final = 50
PREFLIGHT_CONTINGENCY: Final = 1.10

DTYPE: Final = "float64"
HIDDEN_WIDTHS: Final = (DIMENSION + 10, DIMENSION + 10)
LEARNING_RATE: Final = 1.0e-2
VALIDATION_SIZE: Final = 256
Y_INIT_RANGE: Final = (0.0, 1.0)
TERMINAL_LOSS_CLIP: Final = 50.0

EXPECTED_INPUTS: Final = {
    "b6_fullrank": {
        "test_input_sha256": "0ee05fff7ce16f30d0d31563f60df7db3eb23a9c46bbe7bc9be8eb3dbee1311b",
        "test_bank_file_sha256": "6ee6c56ca3d4f413764f28865cb8c00823ea2881ebc90f2643ff4a36fe878730",
        "truth_value_sha256": "6262db7f6b7ee2c41cb6011dadb7288226157533488f708be7d6ea6773f2e13a",
        "parameter_file_sha256": "21b508c30ad63de912979cc9c85ecda0d81aedab1357c88a40e04692495972a4",
        "T": 0.25,
        "driver": "-0.5*sum(z**2)",
    },
    "dense_fourier": {
        "test_input_sha256": "d409296128f48c5ae00c5de2c6ddd998768c01dff18a89c385245a50ff581b70",
        "test_bank_file_sha256": "998fde8b2c0556de4b8fa9a55dec5e489224828cf74e5f3b36f77cb1a8e2c011",
        "truth_value_sha256": "3b33ea36e21595b5c20e5d91097d1fadfcbf3f7884c0328f20579e69995fc816",
        "parameter_file_sha256": "042588d7c3cc8175cc7b9344284725443650837f72008806c668d35b9a005307",
        "T": 0.5,
        "driver": "0.5*tanh(y)*sum(z**2)",
    },
}

# Pre-existing audited 960-bank TA reference-configuration means (four actor
# updates, 128 paths/label, 5000 fit steps/update).  These are frozen context
# for the reviewer question, not new TA rows and not a paired 64-point target.
HISTORICAL_TA_REFERENCE_RMSE: Final = {
    "b6_fullrank": 0.0007200539,
    "dense_fourier": 0.0017584167,
}


@dataclass(frozen=True)
class TrajectorySpec:
    """One independently trained trajectory before query and seed expansion."""

    trajectory_id: str
    updates: int
    batch_size: int
    time_steps: int
    capture_updates: tuple[int, ...]

    def payload(self) -> dict[str, Any]:
        return {**asdict(self), "capture_updates": list(self.capture_updates)}

    @property
    def signature(self) -> str:
        return canonical_sha256(self.payload())


TRAJECTORIES: Final = (
    TrajectorySpec("reference", 4000, 64, 20, (1000, 2000, 4000)),
    TrajectorySpec("c1_batch32", 2000, 32, 20, (2000,)),
    TrajectorySpec("c1_batch128", 2000, 128, 20, (2000,)),
    TrajectorySpec("c2_steps10", 2000, 64, 10, (2000,)),
    TrajectorySpec("c2_steps40", 2000, 64, 40, (2000,)),
    TrajectorySpec("joint_j1", 1000, 32, 10, (1000,)),
    TrajectorySpec("joint_j3", 4000, 128, 40, (4000,)),
)
TRAJECTORY_BY_ID: Final = {spec.trajectory_id: spec for spec in TRAJECTORIES}

# Each displayed point maps to a genuine checkpoint.  Repeated tuples are
# intentional exact reuse, not independently restarted approximations.
DISPLAY_POINTS: Final = (
    {"family": "sweep", "sweep": "C0", "level": 1, "parameter_name": "updates", "parameter_value": 1000, "trajectory_id": "reference", "checkpoint_updates": 1000},
    {"family": "sweep", "sweep": "C0", "level": 2, "parameter_name": "updates", "parameter_value": 2000, "trajectory_id": "reference", "checkpoint_updates": 2000},
    {"family": "sweep", "sweep": "C0", "level": 3, "parameter_name": "updates", "parameter_value": 4000, "trajectory_id": "reference", "checkpoint_updates": 4000},
    {"family": "sweep", "sweep": "C1", "level": 1, "parameter_name": "batch_size", "parameter_value": 32, "trajectory_id": "c1_batch32", "checkpoint_updates": 2000},
    {"family": "sweep", "sweep": "C1", "level": 2, "parameter_name": "batch_size", "parameter_value": 64, "trajectory_id": "reference", "checkpoint_updates": 2000},
    {"family": "sweep", "sweep": "C1", "level": 3, "parameter_name": "batch_size", "parameter_value": 128, "trajectory_id": "c1_batch128", "checkpoint_updates": 2000},
    {"family": "sweep", "sweep": "C2", "level": 1, "parameter_name": "time_steps", "parameter_value": 10, "trajectory_id": "c2_steps10", "checkpoint_updates": 2000},
    {"family": "sweep", "sweep": "C2", "level": 2, "parameter_name": "time_steps", "parameter_value": 20, "trajectory_id": "reference", "checkpoint_updates": 2000},
    {"family": "sweep", "sweep": "C2", "level": 3, "parameter_name": "time_steps", "parameter_value": 40, "trajectory_id": "c2_steps40", "checkpoint_updates": 2000},
    {"family": "joint", "sweep": "J", "level": 1, "parameter_name": "joint", "parameter_value": 1, "trajectory_id": "joint_j1", "checkpoint_updates": 1000},
    {"family": "joint", "sweep": "J", "level": 2, "parameter_name": "joint", "parameter_value": 2, "trajectory_id": "reference", "checkpoint_updates": 2000},
    {"family": "joint", "sweep": "J", "level": 3, "parameter_name": "joint", "parameter_value": 3, "trajectory_id": "joint_j3", "checkpoint_updates": 4000},
)

# The original frozen preregistration and task-index mapping remain intact.
# The user narrowed the executable/reporting scope on 2026-10-27, before the
# truth gate opened, to the individual C0/C1/C2 sweeps only.  Keeping the
# original declarations here prevents pending Slurm array indices from being
# reinterpreted while these filtered views drive completion and analysis.
SCALING_DISPLAY_POINTS: Final = tuple(
    point for point in DISPLAY_POINTS if point["family"] == "sweep"
)
EXCLUDED_JOINT_TRAJECTORIES: Final = ("joint_j1", "joint_j3")


PACKAGE_DIR = Path(__file__).resolve().parent
SOURCE_ROOT = PACKAGE_DIR.parents[2]
RESULTS_DIR = SOURCE_ROOT / "invariant_region_mlp" / "results" / "deep_bsde_budget_scaling_3seed"
INPUT_DIR = RESULTS_DIR / "inputs"
RAW_DIR = RESULTS_DIR / "raw"
CHECKPOINT_DIR = RESULTS_DIR / "checkpoints"
STATUS_DIR = RESULTS_DIR / "status"
LOG_DIR = RESULTS_DIR / "slurm"
TABLE_DIR = RESULTS_DIR / "tables"
REPORT_PATH = SOURCE_ROOT / "invariant_region_mlp" / "docs" / "DEEP_BSDE_BUDGET_SCALING_3SEED_REPORT.md"


def ensure_directories() -> None:
    for path in (RESULTS_DIR, INPUT_DIR, RAW_DIR, CHECKPOINT_DIR, STATUS_DIR, LOG_DIR, TABLE_DIR):
        path.mkdir(parents=True, exist_ok=True)


def canonical_sha256(payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return sha256(encoded).hexdigest()


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
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_seed(*parts: object) -> int:
    digest = sha256("::".join(map(str, parts)).encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "little") % (2**31 - 1)


def write_json_atomic(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as stream:
        temporary = Path(stream.name)
        stream.write(text)
    temporary.replace(path)


def write_npz_atomic(path: Path, **arrays: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as stream:
        np.savez_compressed(stream, **arrays)
    temporary.replace(path)


def historical_results_root() -> Path:
    configured = os.environ.get("DBBS_HISTORICAL_RESULTS") or os.environ.get("WCF_HISTORICAL_RESULTS")
    candidates = [
        Path(configured) if configured else None,
        SOURCE_ROOT / "invariant_region_mlp" / "results" / "main_highdim_expansion",
    ]
    for candidate in candidates:
        if candidate is not None and (candidate / "inputs" / "b6_fullrank_d100_test_bank.npz").is_file():
            return candidate.resolve()
    raise FileNotFoundError("set DBBS_HISTORICAL_RESULTS to the audited main_highdim_expansion result root")


def protocol_payload() -> dict[str, Any]:
    payload = {
        "revision": REVISION,
        "branch": BRANCH,
        "preregistration_commit": PREREGISTRATION_COMMIT,
        "base_revision": BASE_REVISION,
        "official_port_source_revision": OFFICIAL_PORT_SOURCE_REVISION,
        "audited_query_solver_blob": AUDITED_QUERY_SOLVER_BLOB,
        "pdes": list(PDES),
        "dimension": DIMENSION,
        "training_seeds": list(TRAINING_SEEDS),
        "calibration_count": CALIBRATION_COUNT,
        "final_count": FINAL_COUNT,
        "subset_seed": SUBSET_SEED,
        "subset_strata": SUBSET_STRATA,
        "max_h100_gpu_hours": MAX_H100_GPU_HOURS,
        "profile": {
            "dtype": DTYPE,
            "hidden_widths": list(HIDDEN_WIDTHS),
            "learning_rate": LEARNING_RATE,
            "validation_size": VALIDATION_SIZE,
            "y_init_range": list(Y_INIT_RANGE),
            "terminal_loss_clip": TERMINAL_LOSS_CLIP,
        },
        "expected_inputs": EXPECTED_INPUTS,
        "historical_ta_reference_rmse_960_bank": HISTORICAL_TA_REFERENCE_RMSE,
        "trajectories": [spec.payload() | {"signature": spec.signature} for spec in TRAJECTORIES],
        "display_points": list(DISPLAY_POINTS),
        "deep_bsde_fresh_model_per_query": True,
        "truth_used_for_selection_training_or_stopping": False,
    }
    payload["protocol_sha256"] = canonical_sha256(payload)
    return payload


def task_specs() -> list[dict[str, Any]]:
    manifest = json.loads((RESULTS_DIR / "query_manifest.json").read_text(encoding="utf-8"))
    tasks: list[dict[str, Any]] = []
    for pde in PDES:
        final_rows = manifest["pdes"][pde]["final_indices"]
        for seed in TRAINING_SEEDS:
            for trajectory in TRAJECTORIES:
                tasks.append({
                    "task_index": len(tasks),
                    "pde": pde,
                    "seed": seed,
                    "trajectory_id": trajectory.trajectory_id,
                    "trajectory_signature": trajectory.signature,
                    "query_count": len(final_rows),
                })
    return tasks


def scaling_task_specs() -> list[dict[str, Any]]:
    """Return C0/C1/C2 tasks with their original frozen array indices."""

    return [
        task for task in task_specs()
        if task["trajectory_id"] not in EXCLUDED_JOINT_TRAJECTORIES
    ]


__all__ = [
    "CHECKPOINT_DIR", "DISPLAY_POINTS", "SCALING_DISPLAY_POINTS",
    "EXCLUDED_JOINT_TRAJECTORIES", "EXPECTED_INPUTS", "INPUT_DIR",
    "MAX_H100_GPU_HOURS", "PDES", "PREFLIGHT_CONTINGENCY",
    "PREFLIGHT_UPDATES", "RAW_DIR", "REPORT_PATH", "RESULTS_DIR",
    "REVISION", "STATUS_DIR", "TABLE_DIR", "TRAINING_SEEDS",
    "TRAJECTORIES", "TRAJECTORY_BY_ID", "TrajectorySpec", "array_sha256",
    "canonical_sha256", "ensure_directories", "file_sha256",
    "historical_results_root", "protocol_payload", "stable_seed",
    "scaling_task_specs", "task_specs", "write_json_atomic", "write_npz_atomic",
]
