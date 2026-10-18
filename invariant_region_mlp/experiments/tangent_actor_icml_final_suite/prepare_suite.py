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
    PROTOCOL_AMENDMENT, TRAINING_SEEDS, compatibility, dimensions_for,
    frozen_protocol, stable_seed,
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
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
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
        "schema_version": 2,
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
        "seed_policy": {
            "primary_training_seeds": list(TRAINING_SEEDS),
            "historical_extra_seeds": "preserved outside primary aggregation",
            "best_seed_selection": False,
        },
        "review_gate": {
            "status": "pending_independent_review",
            "required_before": "any post-100D dimension sweep",
            "required_artifacts": [
                "BASELINE_PAPER_FIDELITY_REVIEW.md",
                "REMAINING_GPU_HOUR_ESTIMATE.json",
                "100D_RESULTS.csv",
            ],
        },
        "protocol_amendments": [PROTOCOL_AMENDMENT],
        "status": "stage_b_fidelity_gate_pending",
    }
    manifest_path = results / "EXECUTION_MANIFEST.json"
    if manifest_path.exists():
        existing = json.loads(manifest_path.read_text(encoding="utf-8"))
        old_hash = existing.get("protocol_sha256")
        if old_hash != protocol["protocol_sha256"]:
            # This is the one explicitly authorized post-freeze amendment:
            # reduce future training to seed 0 and add the review gate. Merge
            # by scientific cell so no attempts, hashes, or completed runs are
            # lost or relabeled.
            old_cells = {
                (cell["pde"], int(cell["dimension"]), cell["method"]): cell
                for cell in existing.get("method_cells", [])
            }
            for cell in manifest["method_cells"]:
                old = old_cells.get((cell["pde"], int(cell["dimension"]), cell["method"]))
                if old is None:
                    continue
                for key in ("attempts", "completed_seeds", "artifacts"):
                    cell[key] = old.get(key, cell[key])
                if old.get("execution_status", "").startswith("completed"):
                    cell["execution_status"] = old["execution_status"]
            manifest["created_utc"] = existing.get("created_utc", manifest["created_utc"])
            manifest["protocol_history"] = list(existing.get("protocol_history", [])) + [{
                "protocol_sha256": old_hash,
                "superseded_by": protocol["protocol_sha256"],
                "reason": PROTOCOL_AMENDMENT["name"],
                "preserved_completed_runs": True,
            }]
            write_json(manifest_path, manifest)
        else:
            # Refresh declarations while retaining mutable execution evidence.
            manifest = existing
    else:
        write_json(manifest_path, manifest)
    gate_path = results / "REVIEW_GATE_STATUS.json"
    if not gate_path.exists():
        write_json(gate_path, {
            "status": "pending_independent_review",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "scope": "post-100D dimension sweep",
            "approval_must_be_explicit": True,
        })
    write_json(results / "SOURCE_REVISION_AUDIT.json", source_audit)
    print(json.dumps(manifest["counts"], sort_keys=True))


if __name__ == "__main__":
    main()
