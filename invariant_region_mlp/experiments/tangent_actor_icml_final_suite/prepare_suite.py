"""Write the frozen protocol, compatibility grid, and initial execution manifest."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess

from .protocol import (
    BASE_COMMIT, BRANCH, METHODS, PDES, RESULTS_DIR, SOURCE_COMMITS,
    TRAINING_SEEDS, compatibility, dimensions_for, frozen_protocol, stable_seed,
)


def sha256_file(path: Path) -> str:
    block = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            block.update(chunk)
    return block.hexdigest()


def git_head(path: Path) -> str | None:
    try:
        return subprocess.check_output(
            ["git", "-C", str(path), "rev-parse", "HEAD"], text=True
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    tmp.replace(path)


def historical_candidates(source_root: Path) -> list[dict[str, object]]:
    paths = [
        ("original_100d", source_root / "invariant_region_mlp/results/original_pde_benchmark_comparison/final_summary.json"),
        ("original_100d_table", source_root / "invariant_region_mlp/experiments/original_pde_benchmark_comparison/function_approximation_comparison.csv"),
        ("nonconvex_100d", source_root / "invariant_region_mlp/results/tangent_actor_accuracy_frontier/FINAL_SUMMARY.json"),
        ("nonconvex_rows", source_root / "invariant_region_mlp/results/tangent_actor_accuracy_frontier/nonseparable_three_seed_results.csv"),
        ("mlp_failure", source_root / "invariant_region_mlp/results/counterexample_rescue/summary_all.csv"),
        ("joint_pilot", source_root / "invariant_region_mlp/results/convex_concave_semilinear_extended/ground_truth_results.csv"),
    ]
    records = []
    for role, path in paths:
        records.append({
            "role": role,
            "path": str(path.relative_to(source_root)),
            "exists": path.exists(),
            "sha256": sha256_file(path) if path.exists() else None,
            "reuse_status": "candidate_requires_exact_protocol_match",
        })
    return records


def main() -> None:
    source_root = Path(__file__).resolve().parents[3]
    results = RESULTS_DIR
    results.mkdir(parents=True, exist_ok=True)
    protocol = frozen_protocol()
    write_json(results / "FROZEN_PROTOCOL.json", protocol)

    rows: list[dict[str, object]] = []
    cells: list[dict[str, object]] = []
    for pde in PDES:
        for dimension in dimensions_for(pde.key):
            for method in METHODS:
                classification, state, reason = compatibility(pde.key, method, dimension)
                row = {
                    "pde": pde.key,
                    "dimension": dimension,
                    "method": method,
                    "method_status": classification,
                    "execution_status": state,
                    "reason": reason,
                    "native_task": "single_point" if method == "deep_bsde" else "whole_function",
                    "source_commit": SOURCE_COMMITS[method],
                    "target_training_seeds": ";".join(map(str, TRAINING_SEEDS)),
                    "test_bank_seed": stable_seed("icml_final_test_bank", pde.key, dimension),
                    "reference_seed": stable_seed("icml_final_reference", pde.key, dimension),
                }
                rows.append(row)
                cells.append({**row, "attempts": [], "completed_seeds": [], "artifacts": []})

    with (results / "DIMENSION_COMPATIBILITY_MATRIX.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    externals = source_root.parent / "external"
    source_audit = {
        name: {
            "expected": expected,
            "actual": git_head(externals / directory),
            "path": str(externals / directory),
        }
        for name, expected, directory in (
            ("deep_bsde", SOURCE_COMMITS["deep_bsde"], "DeepBSDE"),
            ("dpi", SOURCE_COMMITS["dpi"], "DeepPicardIteration"),
            ("sdgd_pinn", SOURCE_COMMITS["sdgd_pinn"], "SDGD_PINN"),
        )
    }
    for item in source_audit.values():
        item["matches"] = item["actual"] == item["expected"]

    manifest = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "base_commit": BASE_COMMIT,
        "branch": BRANCH,
        "protocol_sha256": protocol["protocol_sha256"],
        "source_audit": source_audit,
        "historical_candidates": historical_candidates(source_root),
        "configurations": [
            {"pde": p.key, "dimension": d}
            for p in PDES for d in dimensions_for(p.key)
        ],
        "method_cells": cells,
        "counts": {
            "pde_dimension_configurations": sum(len(dimensions_for(p.key)) for p in PDES),
            "method_cells": len(cells),
        },
        "status": "stage_a_frozen",
    }
    manifest_path = results / "EXECUTION_MANIFEST.json"
    if manifest_path.exists():
        existing = json.loads(manifest_path.read_text(encoding="utf-8"))
        if existing.get("protocol_sha256") != protocol["protocol_sha256"]:
            raise RuntimeError("existing execution manifest belongs to a different frozen protocol")
    else:
        write_json(manifest_path, manifest)
    write_json(results / "SOURCE_REVISION_AUDIT.json", source_audit)
    print(json.dumps(manifest["counts"], sort_keys=True))


if __name__ == "__main__":
    main()
