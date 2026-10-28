"""Frozen constants and artifact helpers for the high-accuracy study."""

from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
import json
import math
import os
from pathlib import Path
from typing import Any, Final

import numpy as np

from invariant_region_mlp.experiments.single_point_32_ice.protocol import (
    JOINT_ACTOR_CONFIG,
)


REVISION: Final = "high_accuracy_fhmlp_break_even_v1"
PREREGISTRATION_COMMIT: Final = "60c04015682a5f77ff03fd829f575d059620a1bd"
SOURCE_BASE_COMMIT: Final = "49b48e646a8ffa246e4b780c2c54f2a56f24f3a4"
SOURCE_VALIDATED_COMMIT: Final = "b7dd884225936a0d67f770efce1207bc478d4f7c"
BASE_SEED: Final = 20261327
PDE: Final = "joint_nonseparable"
DIMENSION: Final = 100
HORIZON: Final = 0.5

TARGETS: Final = {"primary": 1.0e-2, "secondary": 5.0e-3}
BANK_SIZES: Final = {"calibration": 64, "final": 128, "timing": 1_048_576}
TA_TRAINING_SEEDS: Final = (0, 1, 2)
TA_PATHS: Final = (32, 128, 512)
TA_STEPS: Final = 64
TA_MCV_ACCURACY_REPLICATES: Final = 2
TA_TIMING_REPLICATES: Final = 3
FH_CONFIGS: Final = (
    (2, 32), (2, 48), (2, 64), (2, 96), (2, 128),
    (2, 192), (2, 256), (2, 384), (2, 512),
    (3, 24), (3, 32), (3, 48), (3, 64), (3, 96),
)
FH_ACCURACY_REPLICATES: Final = 5
FH_TIMING_REPLICATES: Final = 3
GEOMETRIC_TIMING_COUNTS: Final = tuple(2**power for power in range(21))
PROFILE_QUERIES: Final = 4
PROFILE_REPLICATES: Final = 2
PROFILE_TASK_CAP_SECONDS: Final = 900.0
CALIBRATION_TASK_CAP_SECONDS: Final = 7200.0
FORMAL_TIMING_TASK_CAP_SECONDS: Final = 10_800.0
MCV_TIMING_TASK_CAP_SECONDS: Final = 3600.0
ALLOCATED_GPU_HOUR_CAP: Final = 30.0
TIMING_GPU_HOUR_RESERVE: Final = 12.0
TA_AND_FINAL_ACCURACY_GPU_HOUR_RESERVE: Final = 5.0
REFERENCE_THRESHOLDS: Final = {
    "value_refinement_rmse": 2.0e-5,
    "value_refinement_max_abs": 1.0e-4,
}
REFINEMENT_MULTIPLIERS: Final = tuple(1.0 + index / 8.0 for index in range(1, 8))

EXPECTED_TA_ARTIFACTS: Final = {
    0: {
        "training_sha256": "ee0876f16cac7df73ab762737f0ea3582d15c3b4e693d4a2ed815390cbe02f76",
        "checkpoint_sha256": "092af2949b7cf4ae7803f62266622f7eb8a749e23ed8816500712d310485080d",
        "selected_actor": "a4",
        "charged_training_seconds": 1360.5530377796385,
    },
    1: {
        "training_sha256": "d395199f578cf8896262ebe5fb24ba3326db4eba518052361107e4f8a122052c",
        "checkpoint_sha256": "65f57b3d47314a975a9918f617cd330515e79ac589b11af69864b241435f5be9",
        "selected_actor": "a4",
        "charged_training_seconds": 1360.9297719779424,
    },
    2: {
        "training_sha256": "e7cd46dcd1503c6a93165bc66d52cfa7623d785e816ab703ae999fe65f0ce149",
        "checkpoint_sha256": "a12b7c4c00fb7d2ddb4223cae15beb14b412ec3c4ee042d53ae051d5023b1428",
        "selected_actor": "a4",
        "charged_training_seconds": 1365.9744462522212,
    },
}

PACKAGE_DIR = Path(__file__).resolve().parent
SOURCE_ROOT = PACKAGE_DIR.parents[2]
RESULTS_DIR = SOURCE_ROOT / "invariant_region_mlp" / "results" / "high_accuracy_fhmlp_break_even"
BANK_DIR = RESULTS_DIR / "banks"
RAW_DIR = RESULTS_DIR / "raw"
PREDICTION_DIR = RESULTS_DIR / "predictions"
FIGURE_DIR = RESULTS_DIR / "figures"
SLURM_DIR = RESULTS_DIR / "slurm"
MANIFEST_DIR = RESULTS_DIR / "manifests"
TASK_STATUS_DIR = RAW_DIR / "task_status"
FROZEN_SELECTIONS = RESULTS_DIR / "FROZEN_SELECTIONS.json"
FROZEN_TIMING_GRID = RESULTS_DIR / "FROZEN_TIMING_GRID.json"
REFERENCE_AUDIT = RESULTS_DIR / "reference_audit.json"
PREFLIGHT_AUDIT = RESULTS_DIR / "preflight_audit.json"
PROFILE_RESULTS = RESULTS_DIR / "performance_profile.json"
COST_GATE = RESULTS_DIR / "cost_gate.json"


def source_results_dir() -> Path:
    override = os.environ.get("IRMLP_JOINT_SOURCE_RESULTS")
    if override:
        return Path(override).expanduser().resolve()
    return SOURCE_ROOT / "invariant_region_mlp" / "results" / "single_point_break_even_joint"


def ensure_directories() -> None:
    for path in (
        RESULTS_DIR, BANK_DIR, RAW_DIR, PREDICTION_DIR, FIGURE_DIR,
        SLURM_DIR, MANIFEST_DIR, TASK_STATUS_DIR,
    ):
        path.mkdir(parents=True, exist_ok=True)


def stable_seed(*parts: object) -> int:
    payload = "::".join((REVISION, str(BASE_SEED), *map(str, parts))).encode("utf-8")
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


def canonical_sha256(payload: Any) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return sha256(encoded).hexdigest()


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


def bank_path(name: str) -> Path:
    if name not in BANK_SIZES:
        raise KeyError(name)
    return BANK_DIR / f"{name}.npz"


def load_bank(name: str, *, require_truth: bool = False) -> dict[str, np.ndarray]:
    with np.load(bank_path(name), allow_pickle=False) as data:
        result = {key: np.asarray(data[key]) for key in data.files}
    if require_truth and "reference" not in result:
        raise RuntimeError(f"bank {name!r} has no reference")
    if name == "timing" and "reference" in result:
        raise RuntimeError("timing bank must never contain truth")
    return result


def _require_validated(path: Path, name: str) -> dict[str, Any]:
    if not path.is_file():
        raise RuntimeError(f"{name} is missing")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("status") != "VALIDATED" or payload.get("overall_pass") is not True:
        raise RuntimeError(f"{name} has not passed")
    return payload


def require_reference() -> dict[str, Any]:
    return _require_validated(REFERENCE_AUDIT, "reference audit")


def require_preflight() -> dict[str, Any]:
    require_reference()
    return _require_validated(PREFLIGHT_AUDIT, "H100 preflight")


def require_cost_gate() -> dict[str, Any]:
    require_preflight()
    if not COST_GATE.is_file():
        raise RuntimeError("cost gate is missing")
    payload = json.loads(COST_GATE.read_text(encoding="utf-8"))
    if payload.get("status") != "PASS":
        raise RuntimeError("cost gate has not passed")
    return payload


def _require_canonical(path: Path, revision: str) -> dict[str, Any]:
    if not path.is_file():
        raise RuntimeError(f"frozen artifact missing: {path.name}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    expected = payload.get("canonical_sha256")
    bare = dict(payload)
    bare.pop("canonical_sha256", None)
    if expected != canonical_sha256(bare):
        raise RuntimeError(f"{path.name} canonical hash mismatch")
    if payload.get("revision") != revision:
        raise RuntimeError(f"{path.name} revision mismatch")
    return payload


def require_frozen_selections() -> dict[str, Any]:
    return _require_canonical(FROZEN_SELECTIONS, REVISION)


def require_frozen_timing_grid() -> dict[str, Any]:
    return _require_canonical(FROZEN_TIMING_GRID, REVISION)


def refined_timing_counts(estimated_crossing: float) -> tuple[int, ...]:
    """Octile refinement of the estimated bracket and its two neighbours."""
    estimate = min(max(float(estimated_crossing), 1.0), float(BANK_SIZES["timing"]))
    values = set(GEOMETRIC_TIMING_COUNTS)
    center_power = int(math.floor(math.log2(estimate)))
    for lower_power in (center_power - 1, center_power, center_power + 1):
        if lower_power < 0 or lower_power >= 20:
            continue
        lower = 2**lower_power
        upper = min(2 * lower, BANK_SIZES["timing"])
        quantum = max(1, lower // 256)
        for multiplier in REFINEMENT_MULTIPLIERS:
            value = int(round((lower * multiplier) / quantum) * quantum)
            values.add(min(max(value, lower + 1), upper - 1))
    return tuple(sorted(value for value in values if 1 <= value <= BANK_SIZES["timing"]))


def protocol_payload() -> dict[str, Any]:
    return {
        "revision": REVISION,
        "preregistration_commit": PREREGISTRATION_COMMIT,
        "source_base_commit": SOURCE_BASE_COMMIT,
        "source_validated_commit": SOURCE_VALIDATED_COMMIT,
        "base_seed": BASE_SEED,
        "pde": PDE,
        "dimension": DIMENSION,
        "horizon": HORIZON,
        "targets": TARGETS,
        "bank_sizes": BANK_SIZES,
        "ta_training_seeds": list(TA_TRAINING_SEEDS),
        "joint_actor_config": asdict(JOINT_ACTOR_CONFIG),
        "ta_paths": list(TA_PATHS),
        "ta_steps": TA_STEPS,
        "ta_mcv_accuracy_replicates": TA_MCV_ACCURACY_REPLICATES,
        "fh_configs": [list(value) for value in FH_CONFIGS],
        "fh_accuracy_replicates": FH_ACCURACY_REPLICATES,
        "geometric_timing_counts": list(GEOMETRIC_TIMING_COUNTS),
        "refinement_multipliers": list(REFINEMENT_MULTIPLIERS),
        "profile_queries": PROFILE_QUERIES,
        "profile_replicates": PROFILE_REPLICATES,
        "caps": {
            "profile_task_seconds": PROFILE_TASK_CAP_SECONDS,
            "calibration_task_seconds": CALIBRATION_TASK_CAP_SECONDS,
            "formal_timing_task_seconds": FORMAL_TIMING_TASK_CAP_SECONDS,
            "mcv_timing_task_seconds": MCV_TIMING_TASK_CAP_SECONDS,
            "allocated_h100_gpu_hours": ALLOCATED_GPU_HOUR_CAP,
            "timing_gpu_hour_reserve": TIMING_GPU_HOUR_RESERVE,
            "ta_and_final_accuracy_gpu_hour_reserve": TA_AND_FINAL_ACCURACY_GPU_HOUR_RESERVE,
        },
        "expected_ta_artifacts": EXPECTED_TA_ARTIFACTS,
        "reference_thresholds": REFERENCE_THRESHOLDS,
    }


__all__ = [name for name in globals() if name.isupper()] + [
    "array_sha256", "bank_path", "canonical_sha256", "ensure_directories",
    "file_sha256", "load_bank", "protocol_payload", "refined_timing_counts",
    "require_cost_gate", "require_frozen_selections", "require_frozen_timing_grid",
    "require_preflight", "require_reference", "source_results_dir", "stable_seed",
    "write_json_atomic", "write_npz_atomic",
]
