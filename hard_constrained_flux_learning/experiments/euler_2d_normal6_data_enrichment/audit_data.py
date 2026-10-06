"""Audit that the enrichment changes only the training trajectories."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import numpy as np


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as archive:
        return {key: archive[key] for key in archive.files}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--baseline-data-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = args.data_dir.resolve()
    baseline = args.baseline_data_dir.resolve()

    enriched = load(data / "train.npz")
    original = load(baseline / "train.npz")
    q = enriched["q"]
    rho = q[..., 0]
    pressure = 0.4 * (
        q[..., 3]
        - 0.5 * (q[..., 1] ** 2 + q[..., 2] ** 2) / rho
    )
    family_names = [str(value) for value in enriched["family_names"].tolist()]
    metadata = json.loads(str(enriched["metadata_json"].item()))
    checks = {
        "original_training_prefix_is_byte_exact": bool(
            np.array_equal(q[: original["q"].shape[0]], original["q"])
        ),
        "saved_times_unchanged": bool(
            np.array_equal(enriched["times"], original["times"])
        ),
        "validation_file_unchanged": sha256(data / "validation.npz")
        == sha256(baseline / "validation.npz"),
        "test_file_unchanged": sha256(data / "test.npz")
        == sha256(baseline / "test.npz"),
        "all_values_finite": bool(np.isfinite(q).all()),
        "all_densities_positive": float(rho.min()) > 0.0,
        "all_pressures_positive": float(pressure.min()) > 0.0,
        "metadata_declares_no_method_change": not bool(
            metadata["model_loss_optimizer_changed"]
        ),
    }
    result = {
        "passed": all(checks.values()),
        "checks": checks,
        "training_sha256": sha256(data / "train.npz"),
        "baseline_training_sha256": sha256(baseline / "train.npz"),
        "validation_sha256": sha256(data / "validation.npz"),
        "test_sha256": sha256(data / "test.npz"),
        "shape": list(q.shape),
        "family_counts": dict(sorted(Counter(family_names).items())),
        "minimum_density": float(rho.min()),
        "minimum_pressure": float(pressure.min()),
        "maximum_density": float(rho.max()),
        "maximum_pressure": float(pressure.max()),
        "metadata": metadata,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True), encoding="utf-8"
    )
    print(json.dumps(result, indent=2), flush=True)
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
