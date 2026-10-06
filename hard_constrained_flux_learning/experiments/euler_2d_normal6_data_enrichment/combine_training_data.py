"""Append targeted 2-D trajectories to the locked original training split."""

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


def scalar_json(value: np.ndarray) -> object:
    return json.loads(str(value.item()))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--extra", type=Path, nargs="+", required=True)
    parser.add_argument("--validation", type=Path, required=True)
    parser.add_argument("--test", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    base = load(args.base)
    extras = [load(path) for path in args.extra]
    reference_times = base["times"]
    expected_shape = base["q"].shape[1:]
    for path, extra in zip(args.extra, extras):
        if not np.array_equal(extra["times"], reference_times):
            raise RuntimeError(f"Time mismatch in {path}")
        if extra["q"].shape[1:] != expected_shape:
            raise RuntimeError(f"Trajectory-shape mismatch in {path}")

    base_specs = scalar_json(base["specs_json"])
    extra_specs = [
        spec for extra in extras for spec in scalar_json(extra["specs_json"])
    ]
    family_names = np.concatenate(
        [base["family_names"], *(extra["family_names"] for extra in extras)]
    )
    q = np.concatenate([base["q"], *(extra["q"] for extra in extras)], axis=0)
    counts = Counter(str(value) for value in family_names.tolist())
    base_metadata = scalar_json(base["metadata_json"])
    extra_metadata = [scalar_json(extra["metadata_json"]) for extra in extras]
    metadata = {
        "purpose": "data-only enrichment ablation",
        "model_loss_optimizer_changed": False,
        "base_training_archive": str(args.base.resolve()),
        "base_training_sha256": sha256(args.base),
        "extra_archives": [
            {
                "path": str(path.resolve()),
                "sha256": sha256(path),
                "metadata": item,
            }
            for path, item in zip(args.extra, extra_metadata)
        ],
        "validation_sha256": sha256(args.validation),
        "test_sha256": sha256(args.test),
        "validation_and_test_unchanged": True,
        "base_count": int(base["q"].shape[0]),
        "extra_count": int(sum(extra["q"].shape[0] for extra in extras)),
        "count": int(q.shape[0]),
        "family_counts": dict(sorted(counts.items())),
        "fine_cells": int(base_metadata["fine_cells"]),
        "coarse_cells": int(base_metadata["coarse_cells"]),
        "tfinal": float(base_metadata["tfinal"]),
        "output_intervals": int(base_metadata["output_intervals"]),
        "restriction": base_metadata["restriction"],
    }
    result = {
        "q": q,
        "times": reference_times,
        "family_names": family_names,
        "specs_json": np.asarray(
            json.dumps([*base_specs, *extra_specs], sort_keys=True)
        ),
        "metadata_json": np.asarray(json.dumps(metadata, sort_keys=True)),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, **result)
    print(
        json.dumps(
            {
                "output": str(args.output),
                "sha256": sha256(args.output),
                "shape": list(q.shape),
                "metadata": metadata,
            },
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
