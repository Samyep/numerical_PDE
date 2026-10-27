"""Frozen constants and artifact helpers for the joint break-even addendum."""

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


REVISION: Final = "single_point_break_even_joint_v1"
ADDENDUM_COMMIT: Final = "8a18360be907da36b8b875c4fdf6ba0b0cc1007d"
SOURCE_IMPORT_COMMIT: Final = "5d4bc017c"
SOURCE_VALIDATED_COMMIT: Final = "b7dd884225936a0d67f770efce1207bc478d4f7c"
BASE_SEED: Final = 20261282
PDE: Final = "joint_nonseparable"
DIMENSION: Final = 100
HORIZON: Final = 0.5

TARGETS: Final = {"epsilon_1": 6.0e-2, "epsilon_2": 2.0e-2}
BANK_SIZES: Final = {"calibration": 32, "final": 64, "timing": 1_048_576}
CALIBRATION_SEEDS: Final = (0, 1, 2)
FINAL_SEEDS: Final = (10, 11, 12)
TA_PATHS: Final = (32, 128, 512)
TA_STEPS: Final = 64
TA_INFERENCE_REPLICATES: Final = 2
FH_CONFIGS: Final = ((2, 8), (2, 16), (2, 32), (3, 8), (3, 12), (3, 16), (3, 24))
FH_REPLICATES: Final = 5
FH_TIMING_REPLICATES: Final = 3
DB_CANDIDATE_UPDATES: Final = (500, 1000, 2000, 4000)
TIMING_COUNTS: Final = tuple(2**power for power in range(21))
DB_TIMING_COUNTS: Final = (1, 2, 4, 8, 16, 24, 32, 48, 64)
ACTIVE_CAP_SECONDS: Final = 3600.0
ALLOCATED_GPU_HOUR_CAP: Final = 18.0

PACKAGE_DIR = Path(__file__).resolve().parent
SOURCE_ROOT = PACKAGE_DIR.parents[2]
RESULTS_DIR = SOURCE_ROOT / "invariant_region_mlp" / "results" / "single_point_break_even_joint"
BANK_DIR = RESULTS_DIR / "banks"
RAW_DIR = RESULTS_DIR / "raw"
CACHE_DIR = RESULTS_DIR / "cache"
CHECKPOINT_DIR = RESULTS_DIR / "checkpoints"
PREDICTION_DIR = RESULTS_DIR / "predictions"
LOG_DIR = RESULTS_DIR / "logs"
FIGURE_DIR = RESULTS_DIR / "figures"
SLURM_DIR = RESULTS_DIR / "slurm"
FROZEN_SELECTIONS = RESULTS_DIR / "FROZEN_SELECTIONS.json"
REFERENCE_AUDIT = RESULTS_DIR / "reference_audit.json"
PREFLIGHT_AUDIT = RESULTS_DIR / "preflight_audit.json"


def ensure_directories() -> None:
    for path in (
        RESULTS_DIR, BANK_DIR, RAW_DIR, CACHE_DIR, CHECKPOINT_DIR,
        PREDICTION_DIR, LOG_DIR, FIGURE_DIR, SLURM_DIR,
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
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    return sha256(encoded).hexdigest()


def write_json_atomic(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
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


def require_reference() -> dict[str, Any]:
    if not REFERENCE_AUDIT.is_file():
        raise RuntimeError("reference audit is missing")
    payload = json.loads(REFERENCE_AUDIT.read_text(encoding="utf-8"))
    if payload.get("status") != "VALIDATED" or payload.get("overall_pass") is not True:
        raise RuntimeError("reference audit has not passed")
    return payload


def require_preflight() -> dict[str, Any]:
    require_reference()
    if not PREFLIGHT_AUDIT.is_file():
        raise RuntimeError("H100 preflight is missing")
    payload = json.loads(PREFLIGHT_AUDIT.read_text(encoding="utf-8"))
    if payload.get("status") != "VALIDATED" or payload.get("overall_pass") is not True:
        raise RuntimeError("H100 preflight has not passed")
    return payload


def require_frozen_selections() -> dict[str, Any]:
    if not FROZEN_SELECTIONS.is_file():
        raise RuntimeError("calibration selections are not frozen")
    payload = json.loads(FROZEN_SELECTIONS.read_text(encoding="utf-8"))
    expected = payload.get("canonical_sha256")
    bare = dict(payload)
    bare.pop("canonical_sha256", None)
    if expected != canonical_sha256(bare):
        raise RuntimeError("FROZEN_SELECTIONS canonical hash mismatch")
    if payload.get("revision") != REVISION:
        raise RuntimeError("FROZEN_SELECTIONS revision mismatch")
    return payload


def protocol_payload() -> dict[str, Any]:
    return {
        "revision": REVISION,
        "addendum_commit": ADDENDUM_COMMIT,
        "source_import_commit": SOURCE_IMPORT_COMMIT,
        "source_validated_commit": SOURCE_VALIDATED_COMMIT,
        "base_seed": BASE_SEED,
        "pde": PDE,
        "dimension": DIMENSION,
        "horizon": HORIZON,
        "targets": TARGETS,
        "bank_sizes": BANK_SIZES,
        "calibration_seeds": list(CALIBRATION_SEEDS),
        "final_seeds": list(FINAL_SEEDS),
        "joint_actor_config": asdict(JOINT_ACTOR_CONFIG),
        "ta_paths": list(TA_PATHS),
        "ta_steps": TA_STEPS,
        "ta_inference_replicates": TA_INFERENCE_REPLICATES,
        "fh_configs": [list(value) for value in FH_CONFIGS],
        "fh_replicates": FH_REPLICATES,
        "deep_bsde_candidate_updates": list(DB_CANDIDATE_UPDATES),
        "timing_counts": list(TIMING_COUNTS),
        "deep_bsde_timing_counts": list(DB_TIMING_COUNTS),
        "active_cap_seconds": ACTIVE_CAP_SECONDS,
        "allocated_gpu_hour_cap": ALLOCATED_GPU_HOUR_CAP,
    }


__all__ = [name for name in globals() if name.isupper()] + [
    "array_sha256", "bank_path", "canonical_sha256", "ensure_directories",
    "file_sha256", "load_bank", "protocol_payload", "require_frozen_selections",
    "require_preflight", "require_reference", "stable_seed", "write_json_atomic",
    "write_npz_atomic",
]

