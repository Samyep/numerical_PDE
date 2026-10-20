"""Freeze equation parameters and all state banks before training."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Callable

import numpy as np

from invariant_region_mlp.experiments.tangent_actor_icml_final_suite.protocol import (
    stable_seed,
)

from .equations import array_sha256, sample_unit_ball
from .protocol import (
    DENSE_PARAMETER_SEED,
    DIMENSION,
    INPUT_ROOT,
    LQG_COEFFICIENT_SEED,
    MANIFEST_PATH,
    SELECTOR_STATES,
    STATE_SEED_NAMESPACE,
    TERMINAL_TIME,
    TEST_STATES,
    VALIDATION_STATES,
    ensure_directories,
    frozen_configuration,
    write_json_atomic,
)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_npz(path: Path, **arrays: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp.npz")
    np.savez_compressed(temporary, **arrays)
    temporary.replace(path)


def _time(rng: np.random.Generator, count: int) -> np.ndarray:
    return rng.uniform(
        0.0, np.nextafter(TERMINAL_TIME, 0.0), int(count)
    ).astype(np.float64)


def _bank(
    pde: str,
    sampler: Callable[[np.random.Generator, int, int], np.ndarray],
) -> tuple[dict[str, np.ndarray], dict[str, object]]:
    arrays: dict[str, np.ndarray] = {}
    splits = {
        "test": TEST_STATES,
        "validation": VALIDATION_STATES,
        "selector": SELECTOR_STATES,
        "train_a1": 2048,
        "train_a2": 2048,
    }
    seed_manifest: dict[str, int] = {}
    for split, count in splits.items():
        seed = stable_seed(STATE_SEED_NAMESPACE, pde, split, count, DIMENSION)
        rng = np.random.default_rng(seed)
        arrays[f"{split}_t"] = _time(rng, count)
        arrays[f"{split}_x"] = sampler(rng, count, DIMENSION)
        seed_manifest[split] = seed
    return arrays, {
        "seeds": seed_manifest,
        "counts": splits,
        "array_sha256": array_sha256(*arrays.values()),
    }


def _cube(rng: np.random.Generator, count: int, d: int) -> np.ndarray:
    return rng.uniform(-1.0, 1.0, (int(count), int(d))).astype(np.float64)


def prepare() -> dict[str, object]:
    ensure_directories()
    expected_files = (
        INPUT_ROOT / "dense_fourier_parameters.npz",
        INPUT_ROOT / "coupled_lqg_parameters.npz",
        INPUT_ROOT / "dense_fourier_banks.npz",
        INPUT_ROOT / "coupled_lqg_banks.npz",
    )
    if MANIFEST_PATH.is_file() and all(path.is_file() for path in expected_files):
        existing = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
        if existing.get("configuration") != frozen_configuration():
            raise RuntimeError("frozen input configuration no longer matches code")
        if existing.get("status") != "FROZEN_BEFORE_TRAINING":
            raise RuntimeError("unexpected input manifest status")
        return existing
    dense_rng = np.random.default_rng(np.random.SeedSequence(DENSE_PARAMETER_SEED))
    omega = dense_rng.standard_normal((2 * DIMENSION, DIMENSION))
    omega /= np.linalg.norm(omega, axis=1, keepdims=True)
    dense_signs = (2 * dense_rng.integers(0, 2, 2 * DIMENSION) - 1).astype(
        np.float64
    )
    phases = dense_rng.uniform(0.0, 2.0 * np.pi, 2 * DIMENSION)
    dense_path = INPUT_ROOT / "dense_fourier_parameters.npz"
    _atomic_npz(
        dense_path,
        omega=omega.astype(np.float64),
        signs=dense_signs,
        phases=phases.astype(np.float64),
        seed=np.asarray(DENSE_PARAMETER_SEED, dtype=np.int64),
    )

    lqg_rng = np.random.default_rng(np.random.SeedSequence(LQG_COEFFICIENT_SEED))
    c1 = lqg_rng.uniform(0.5, 1.5, DIMENSION - 1).astype(np.float64)
    c2 = lqg_rng.uniform(0.5, 1.5, DIMENSION - 1).astype(np.float64)
    lqg_path = INPUT_ROOT / "coupled_lqg_parameters.npz"
    _atomic_npz(
        lqg_path,
        c1=c1,
        c2=c2,
        seed=np.asarray(LQG_COEFFICIENT_SEED, dtype=np.int64),
    )

    dense_bank, dense_meta = _bank("dense_fourier", _cube)
    dense_bank_path = INPUT_ROOT / "dense_fourier_banks.npz"
    _atomic_npz(dense_bank_path, **dense_bank)
    lqg_bank, lqg_meta = _bank("coupled_lqg", sample_unit_ball)
    lqg_bank_path = INPUT_ROOT / "coupled_lqg_banks.npz"
    _atomic_npz(lqg_bank_path, **lqg_bank)

    manifest: dict[str, object] = {
        "status": "FROZEN_BEFORE_TRAINING",
        "configuration": frozen_configuration(),
        "dense_fourier": {
            "parameter_file": str(dense_path),
            "parameter_file_sha256": file_sha256(dense_path),
            "parameter_array_sha256": array_sha256(omega, dense_signs, phases),
            "omega_sha256": array_sha256(omega),
            "signs_sha256": array_sha256(dense_signs),
            "phases_sha256": array_sha256(phases),
            "bank_file": str(dense_bank_path),
            "bank_file_sha256": file_sha256(dense_bank_path),
            "state_law": "t~Uniform[0,0.5), x~Uniform[-1,1]^100",
            **dense_meta,
        },
        "coupled_lqg": {
            "benchmark_label": "RECREATED_PUBLISHED_BENCHMARK",
            "parameter_file": str(lqg_path),
            "parameter_file_sha256": file_sha256(lqg_path),
            "parameter_array_sha256": array_sha256(c1, c2),
            "bank_file": str(lqg_bank_path),
            "bank_file_sha256": file_sha256(lqg_bank_path),
            "state_law": (
                "t~Uniform[0,0.5), x uniform by volume on the Euclidean unit ball"
            ),
            "paper_url": (
                "https://proceedings.iclr.cc/paper_files/paper/2026/file/"
                "cf8353a0efc3644454a5591a7cb55c17-Paper-Conference.pdf"
            ),
            "official_code_url": "https://github.com/Francis-Fan-create/SCaSML",
            "official_code_commit_audited": "c2412b942b1720084accf9ac9651c6085e073919",
            "source_resolution": (
                "ICLR Section 3.3 fixes [0,0.5] x B^d. The repository defaults "
                "to T=1 and creates coefficient arrays as a function of batch shape; "
                "we follow the paper interval and freeze one global coefficient draw."
            ),
            **lqg_meta,
        },
    }
    write_json_atomic(MANIFEST_PATH, manifest)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.parse_args()
    manifest = prepare()
    print(f"frozen inputs: {MANIFEST_PATH}")
    print(manifest["dense_fourier"]["parameter_array_sha256"])
    print(manifest["coupled_lqg"]["parameter_array_sha256"])


if __name__ == "__main__":
    main()


__all__ = ["file_sha256", "prepare"]

