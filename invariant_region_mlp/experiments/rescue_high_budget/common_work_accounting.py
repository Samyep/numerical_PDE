"""Shared work accounting and immutable-result helpers for rescue experiments."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
from typing import Any

import numpy as np


@dataclass
class WorkCounters:
    terminal_g_evals: int = 0
    terminal_samples: int = 0
    f_evals: int = 0
    recursively_evaluated_states: int = 0
    transition_samples: int = 0
    standard_normal_variates: int = 0
    time_uniform_variates: int = 0
    correction_opportunities: int = 0
    constraint_violations: int = 0
    activated_states: int = 0
    overshoot_energy_sum: float = 0.0
    nonfinite_states: int = 0
    nonfinite_generators: int = 0

    def merge(self, other: "WorkCounters") -> None:
        for name in (
            "terminal_g_evals",
            "terminal_samples",
            "f_evals",
            "recursively_evaluated_states",
            "transition_samples",
            "standard_normal_variates",
            "time_uniform_variates",
            "correction_opportunities",
            "constraint_violations",
            "activated_states",
            "nonfinite_states",
            "nonfinite_generators",
        ):
            setattr(self, name, int(getattr(self, name) + getattr(other, name)))
        self.overshoot_energy_sum = float(
            self.overshoot_energy_sum + other.overshoot_energy_sum
        )

    def summary(self) -> dict[str, Any]:
        opportunities = max(self.correction_opportunities, 1)
        payload = asdict(self)
        payload.update(
            {
                "total_stochastic_samples": self.terminal_samples
                + self.transition_samples,
                "constraint_violation_rate": self.constraint_violations
                / opportunities,
                "projection_activation_rate": self.activated_states / opportunities,
                "mean_overshoot_energy": self.overshoot_energy_sum / opportunities,
            }
        )
        return payload


class DrawRecorder:
    """RNG wrapper that makes paired-tree equality directly auditable."""

    def __init__(self, rng: np.random.Generator, *, enabled: bool = True) -> None:
        self.rng = rng
        self.enabled = bool(enabled)
        self._hash = hashlib.sha256() if enabled else None

    def _record(self, values: np.ndarray) -> np.ndarray:
        values = np.asarray(values, dtype=np.float64)
        if self._hash is not None:
            contiguous = np.ascontiguousarray(values)
            self._hash.update(np.asarray(contiguous.shape, dtype=np.int64).tobytes())
            self._hash.update(contiguous.tobytes())
        return values

    def normal(self, shape: tuple[int, ...]) -> np.ndarray:
        return self._record(self.rng.standard_normal(shape, dtype=np.float64))

    def uniform(self, shape: tuple[int, ...]) -> np.ndarray:
        return self._record(self.rng.random(shape, dtype=np.float64))

    def power(self, alpha: float, shape: tuple[int, ...]) -> np.ndarray:
        return self._record(self.rng.power(alpha, size=shape))

    @property
    def fingerprint(self) -> str | None:
        return None if self._hash is None else self._hash.hexdigest()


def combine_fingerprints(fingerprints: list[str]) -> str:
    digest = hashlib.sha256()
    for value in fingerprints:
        digest.update(value.encode("ascii"))
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def git_output(repository_root: Path, *args: str) -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(repository_root), *args],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def environment_summary() -> dict[str, str]:
    return {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "platform": platform.platform(),
    }


def atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp-{os.getpid()}")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False),
        encoding="utf-8",
    )
    os.replace(temporary, path)


def save_npz_atomic(path: Path, **arrays: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.stem + f".tmp-{os.getpid()}.npz")
    np.savez_compressed(temporary, **arrays)
    os.replace(temporary, path)


def metadata_array(payload: dict[str, Any]) -> np.ndarray:
    return np.asarray(json.dumps(payload, sort_keys=True, allow_nan=False))


def load_metadata(data: Any) -> dict[str, Any]:
    return json.loads(str(data["metadata_json"]))
