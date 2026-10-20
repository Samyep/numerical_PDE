"""Freeze parameters, common 960-state banks, and exact references."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .benchmarks import (
    array_sha256,
    dense_parameters,
    exact_state,
    fixed_b6_points,
    make_benchmark,
    metadata,
)
from .protocol import (
    BASE_SEED,
    CALIBRATION_PATH,
    DIMENSIONS,
    INPUT_DIR,
    MANIFEST_PATH,
    PDES,
    REVISION,
    TEST_STATES,
    ensure_directories,
    frozen_protocol,
    stable_seed,
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


def _bank(pde: str, d: int) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    if pde == "b6_fullrank":
        all_t, all_x, is_validation = fixed_b6_points(int(d))
        t, x = all_t[~is_validation].copy(), all_x[~is_validation].copy()
        provenance = {
            "seed_namespace": [BASE_SEED, d, len(all_t), "B6 stable point seed"],
            "full_bank_sha256": array_sha256(all_t, all_x, is_validation),
            "compatibility": "exact highdim_nosym test subset",
        }
    elif pde == "dense_fourier":
        seed = stable_seed(REVISION, "test_bank", pde, d, TEST_STATES)
        rng = np.random.default_rng(seed)
        t = rng.uniform(0.0, np.nextafter(0.5, 0.0), TEST_STATES).astype(np.float64)
        x = rng.uniform(-1.0, 1.0, (TEST_STATES, d)).astype(np.float64)
        provenance = {"seed": seed, "seed_namespace": [REVISION, "test_bank", pde, d, TEST_STATES]}
    else:
        raise KeyError(pde)
    if len(t) != TEST_STATES or x.shape != (TEST_STATES, d):
        raise AssertionError((pde, d, t.shape, x.shape))
    return t, x, provenance


def _expected_paths() -> list[Path]:
    return [
        INPUT_DIR / f"{pde}_d{d}_parameters.npz"
        for pde in PDES for d in DIMENSIONS
    ] + [
        INPUT_DIR / f"{pde}_d{d}_test_bank.npz"
        for pde in PDES for d in DIMENSIONS
    ]


def _verify_existing(manifest: dict[str, Any]) -> None:
    if manifest.get("status") != "FROZEN_BEFORE_TRAINING":
        raise RuntimeError("input manifest has an unexpected status")
    if manifest.get("protocol") != frozen_protocol():
        raise RuntimeError("input manifest protocol drift")
    for item in manifest["cells"]:
        for key in ("parameter_file", "test_bank_file"):
            path = Path(item[key])
            if not path.is_file():
                raise RuntimeError(f"missing frozen artifact: {path}")
            expected = item[f"{key}_sha256"]
            if file_sha256(path) != expected:
                raise RuntimeError(f"frozen artifact hash drift: {path}")


def prepare() -> dict[str, Any]:
    ensure_directories()
    if not CALIBRATION_PATH.is_file():
        raise RuntimeError("B6 calibration must be frozen before input preparation")
    calibration = json.loads(CALIBRATION_PATH.read_text(encoding="utf-8"))
    if calibration.get("status") != "FROZEN_BEFORE_TRAINING":
        raise RuntimeError("B6 calibration was not frozen before training")
    if MANIFEST_PATH.is_file():
        existing = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
        _verify_existing(existing)
        return existing
    if any(path.exists() for path in _expected_paths()):
        raise RuntimeError("partial frozen inputs exist without a manifest; refusing overwrite")

    cells: list[dict[str, Any]] = []
    for pde in PDES:
        for d in DIMENSIONS:
            benchmark = make_benchmark(pde, d)
            if pde == "dense_fourier":
                omega, signs, phases = dense_parameters(d)
                parameter_arrays = {
                    "omega": omega,
                    "signs": signs,
                    "phases": phases,
                    "seed": np.asarray([BASE_SEED, d, 41], dtype=np.int64),
                }
            else:
                equation = benchmark.equation
                parameter_arrays = {
                    "directions": equation.e,
                    "multipliers": equation.U,
                    "common_vector": equation.v,
                    "vectors": equation.b,
                    "squared_norms": equation.n2,
                    "beta": np.asarray(equation.beta),
                    "gamma": np.asarray(equation.gamma),
                    "seed": np.asarray([BASE_SEED, d, 66], dtype=np.int64),
                }
            parameter_path = INPUT_DIR / f"{pde}_d{d}_parameters.npz"
            _atomic_npz(parameter_path, **parameter_arrays)

            t, x, bank_provenance = _bank(pde, d)
            truth = exact_state(benchmark, t, x)
            if truth.shape != (TEST_STATES, d + 1) or not np.all(np.isfinite(truth)):
                raise FloatingPointError((pde, d, truth.shape))
            bank_path = INPUT_DIR / f"{pde}_d{d}_test_bank.npz"
            _atomic_npz(
                bank_path,
                t=t,
                x=x,
                truth_u=truth[:, 0],
                truth_z=truth[:, 1:],
            )
            cells.append(
                {
                    "pde": pde,
                    "dimension": d,
                    "parameter_file": str(parameter_path),
                    "parameter_file_sha256": file_sha256(parameter_path),
                    "parameter_array_sha256": array_sha256(*parameter_arrays.values()),
                    "test_bank_file": str(bank_path),
                    "test_bank_file_sha256": file_sha256(bank_path),
                    "test_input_sha256": array_sha256(t, x),
                    "truth_sha256": array_sha256(truth[:, 0], truth[:, 1:]),
                    "truth_value_sha256": array_sha256(truth[:, 0]),
                    "truth_z_sha256": array_sha256(truth[:, 1:]),
                    "test_states": len(t),
                    "state_law": f"t~Uniform[0,{benchmark.T}), x~Uniform([-1,1]^{d})",
                    "benchmark": metadata(benchmark),
                    "bank_provenance": bank_provenance,
                    "nonfinite_parameters": int(
                        sum(np.size(a) - np.count_nonzero(np.isfinite(a)) for a in parameter_arrays.values())
                    ),
                    "nonfinite_bank": int(
                        np.size(t) + np.size(x) + np.size(truth)
                        - np.count_nonzero(np.isfinite(t))
                        - np.count_nonzero(np.isfinite(x))
                        - np.count_nonzero(np.isfinite(truth))
                    ),
                }
            )
    manifest = {
        "status": "FROZEN_BEFORE_TRAINING",
        "revision": REVISION,
        "protocol": frozen_protocol(),
        "calibration_file": str(CALIBRATION_PATH),
        "calibration_file_sha256": file_sha256(CALIBRATION_PATH),
        "cells": cells,
    }
    write_json_atomic(MANIFEST_PATH, manifest)
    return manifest


def load_test_bank(pde: str, d: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Path]:
    manifest = prepare()
    rows = [item for item in manifest["cells"] if item["pde"] == pde and item["dimension"] == int(d)]
    if len(rows) != 1:
        raise KeyError((pde, d))
    path = Path(rows[0]["test_bank_file"])
    with np.load(path) as data:
        t = data["t"].copy()
        x = data["x"].copy()
        truth_u = data["truth_u"].copy()
        truth_z = data["truth_z"].copy()
    return t, x, truth_u, truth_z, path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.parse_args()
    result = prepare()
    print(f"frozen {len(result['cells'])} parameter/test-bank cells at {MANIFEST_PATH}")


if __name__ == "__main__":
    main()
