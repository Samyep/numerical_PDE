"""Independently recompute all Value metrics from immutable saved arrays."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from .common import PREDICTION_DIR, RESULTS_DIR, atomic_json, model_metrics, write_csv


PROBLEMS = ("allen_cahn", "joint_concave", "joint_convex")


def _read_csv(path: Path) -> dict[tuple[str, ...], dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if path.name == "results.csv":
        return {(row["pde"], row["iteration"]): row for row in rows}
    return {(row["pde"], row["actor"], row["variant"]): row for row in rows}


def _hash_arrays(arrays: np.lib.npyio.NpzFile) -> str:
    digest = hashlib.sha256()
    for key in sorted(arrays.files):
        value = np.ascontiguousarray(arrays[key])
        digest.update(key.encode())
        digest.update(str(value.shape).encode())
        digest.update(value.tobytes())
    return digest.hexdigest()


def main() -> None:
    reported = _read_csv(RESULTS_DIR / "results.csv")
    reported_floor = _read_csv(RESULTS_DIR / "mc_floor_discretization.csv")
    rows = []
    differences = []
    hashes = {}
    for problem in PROBLEMS:
        with np.load(PREDICTION_DIR / f"{problem}.npz") as arrays:
            hashes[problem] = _hash_arrays(arrays)
            truth = arrays["reference_original_value"]
            actors = tuple(
                actor for actor in ("a0", "a1", "a2", "a3", "exact_tangent")
                if (problem, actor) in reported
            )
            for actor in actors:
                metrics = model_metrics(arrays[f"{actor}_original_value"], truth)
                source = reported[(problem, actor)]
                delta = max(
                    abs(metrics[key] - float(source[key]))
                    for key in ("rmse", "relative_l2", "bias", "rmse_over_truth_std")
                )
                differences.append(delta)
                rows.append({"pde": problem, "actor": actor, "variant": "main", **metrics,
                             "max_abs_difference_from_report": delta})
            for variant in ("main", "high_paths", "high_steps"):
                for actor in ("a2", "a3", "exact_tangent"):
                    if (problem, actor, variant) not in reported_floor:
                        continue
                    metrics = model_metrics(
                        arrays[f"floor_{actor}_{variant}_original_value"], truth
                    )
                    source = reported_floor[(problem, actor, variant)]
                    delta = max(
                        abs(metrics[key] - float(source[key]))
                        for key in ("rmse", "relative_l2", "bias", "rmse_over_truth_std")
                    )
                    differences.append(delta)
    write_csv(RESULTS_DIR / "recomputed_metrics.csv", rows)
    maximum = max(differences, default=0.0)
    atomic_json(
        RESULTS_DIR / "independent_recomputation.json",
        {
            "method": "metrics recomputed only from committed NPZ prediction arrays",
            "prediction_array_sha256": hashes,
            "maximum_absolute_scalar_difference": maximum,
            "tolerance": 1e-12,
            "passed": maximum <= 1e-12,
        },
    )
    if maximum > 1e-12:
        raise RuntimeError(f"independent recomputation mismatch: {maximum}")


if __name__ == "__main__":
    main()
