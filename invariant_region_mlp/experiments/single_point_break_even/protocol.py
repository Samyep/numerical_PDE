"""Frozen constants and artifact helpers for the break-even study."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Final

import numpy as np


REVISION: Final = "single_point_break_even_v1"
PREREGISTRATION_COMMIT: Final = "498589cc2"
SOURCE_COMMIT: Final = "9ce44b4393d2e54ed461b69a214af0d2019bb304"
NAMESPACE_SEED: Final = 20271027
PDE: Final = "allen_cahn"
DIMENSION: Final = 100
HORIZON: Final = 0.3

TARGETS: Final = {"epsilon_1": 2.5e-4, "epsilon_2": 1.0e-4}
BANK_SIZES: Final = {"calibration": 32, "final": 64, "timing": 524_288}
CALIBRATION_SEEDS: Final = (0, 1, 2)
FINAL_SEEDS: Final = (10, 11, 12)
FH_REPLICATES: Final = 5
FH_TIMING_REPLICATES: Final = 3
TA_INFERENCE_REPLICATES: Final = 2
TA_PATHS: Final = (32, 64, 128, 256, 512, 1024)
TA_STEPS: Final = 64
FH_CONFIGS: Final = (
    (3, 12), (3, 16), (3, 24), (3, 32), (3, 48), (3, 64), (3, 96),
    (4, 6), (4, 8), (4, 10), (4, 12), (4, 16),
)
DB_CANDIDATE_UPDATES: Final = (2000, 4000, 8000)
TIMING_COUNTS: Final = (
    1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096,
    8192, 16384, 32768, 65536, 131072, 262144, 524288,
)
DB_TIMING_COUNTS: Final = (1, 2, 4, 8, 16, 24, 32, 48, 64)

PACKAGE_DIR = Path(__file__).resolve().parent
SOURCE_ROOT = PACKAGE_DIR.parents[2]
RESULTS_DIR = SOURCE_ROOT / "invariant_region_mlp" / "results" / "single_point_break_even"
BANK_DIR = RESULTS_DIR / "banks"
RAW_DIR = RESULTS_DIR / "raw"
CACHE_DIR = RESULTS_DIR / "cache"
CHECKPOINT_DIR = RESULTS_DIR / "checkpoints"
PREDICTION_DIR = RESULTS_DIR / "predictions"
LOG_DIR = RESULTS_DIR / "logs"
FIGURE_DIR = RESULTS_DIR / "figures"
SLURM_DIR = RESULTS_DIR / "slurm"
FROZEN_SELECTIONS = RESULTS_DIR / "FROZEN_SELECTIONS.json"


def ensure_directories() -> None:
    for path in (
        RESULTS_DIR, BANK_DIR, RAW_DIR, CACHE_DIR, CHECKPOINT_DIR,
        PREDICTION_DIR, LOG_DIR, FIGURE_DIR, SLURM_DIR,
    ):
        path.mkdir(parents=True, exist_ok=True)


def _spawn_word(part: object) -> int:
    return int.from_bytes(sha256(str(part).encode("utf-8")).digest()[:4], "little")


def stable_seed(*parts: object) -> int:
    """Deterministic uint32 seed derived through NumPy SeedSequence."""
    spawn_key = tuple(_spawn_word(part) for part in (REVISION, *parts))
    sequence = np.random.SeedSequence(NAMESPACE_SEED, spawn_key=spawn_key)
    return int(sequence.generate_state(1, dtype=np.uint32)[0])


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
        "preregistration_commit": PREREGISTRATION_COMMIT,
        "source_commit": SOURCE_COMMIT,
        "namespace_seed": NAMESPACE_SEED,
        "pde": PDE,
        "dimension": DIMENSION,
        "horizon": HORIZON,
        "targets": TARGETS,
        "bank_sizes": BANK_SIZES,
        "calibration_seeds": list(CALIBRATION_SEEDS),
        "final_seeds": list(FINAL_SEEDS),
        "fh_replicates": FH_REPLICATES,
        "fh_timing_replicates": FH_TIMING_REPLICATES,
        "ta_inference_replicates": TA_INFERENCE_REPLICATES,
        "ta_paths": list(TA_PATHS),
        "ta_steps": TA_STEPS,
        "fh_configs": [list(value) for value in FH_CONFIGS],
        "deep_bsde_candidate_updates": list(DB_CANDIDATE_UPDATES),
        "timing_counts": list(TIMING_COUNTS),
        "deep_bsde_timing_counts": list(DB_TIMING_COUNTS),
    }


__all__ = [name for name in globals() if name.isupper()] + [
    "array_sha256", "bank_path", "canonical_sha256", "ensure_directories",
    "file_sha256", "load_bank", "protocol_payload", "require_frozen_selections",
    "stable_seed", "write_json_atomic", "write_npz_atomic",
]
