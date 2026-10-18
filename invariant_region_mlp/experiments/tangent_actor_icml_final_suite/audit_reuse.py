"""Strict audit of historical artifacts eligible for the final suite."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any

from .prepare_suite import sha256_file, write_json
from .protocol import RESULTS_DIR


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _artifact(path: Path, source_root: Path) -> dict[str, Any]:
    return {
        "path": str(path.relative_to(source_root)),
        "sha256": sha256_file(path),
        "bytes": path.stat().st_size,
    }


def _file_hash_matches(path: Path, expected: str) -> bool:
    return path.exists() and sha256_file(path) == expected


def _original_audit(source_root: Path) -> dict[str, Any]:
    experiment = source_root / "invariant_region_mlp/experiments/original_pde_benchmark_comparison"
    results = source_root / "invariant_region_mlp/results/original_pde_benchmark_comparison"
    metadata = _read_json(experiment / "reproducibility_metadata.json")
    definitions = _read_json(results / "pde_definitions_and_references.json")
    expected_hashes = metadata["evaluation_bank_sha256"]
    bank_map = {
        "hjb_log_d100": results / "evaluation_bank_hjb_log_d100.npz",
        "dpi_ou_gmm_d100": results / "evaluation_bank_dpi_ou_gmm_d100.npz",
        "allen_cahn_d100": results / "evaluation_bank_allen_cahn_d100.npz",
    }
    bank_checks = {
        key: {
            "expected_sha256": expected_hashes[key],
            "actual_sha256": sha256_file(path),
            "matches": _file_hash_matches(path, expected_hashes[key]),
            "path": str(path.relative_to(source_root)),
        }
        for key, path in bank_map.items()
    }
    source_expected = {
        "DeepBSDE": "e76fed80b1995daf4a0dd6d3e2cf64a0931dabda",
        "DeepPicardIteration": "e8616c2f531c9f95299b246a13c1180b8f02822a",
        "SDGD_PINN": "c597cc4fe9caac0b3ca118ea84453405f25ff833",
    }
    source_checks = {
        key: metadata["official_source_commits"].get(key) == value
        for key, value in source_expected.items()
    }
    equations = definitions["benchmarks"]
    definition_checks = {
        "hjb_log": equations["hjb_log_d100"]["equation"] == "u_t + Delta u - ||grad u||^2 = 0",
        "ou_gmm": equations["dpi_ou_gmm_d100"]["dimension"] == 100
        and equations["dpi_ou_gmm_d100"]["z_convention"] == "z=sigma*grad(u)",
        "allen_cahn": equations["allen_cahn_d100"]["terminal"] == "1/(2+0.4||x||^2)",
    }
    selections = {}
    for key in bank_map:
        selection = _read_json(results / f"final_actor_selection_{key}.json")
        selections[key] = {
            "ground_truth_used": selection["ground_truth_used"],
            "final_test_opened": selection["final_test_opened"],
            "passes": not selection["ground_truth_used"] and not selection["final_test_opened"],
        }
    passed = (
        all(item["matches"] for item in bank_checks.values())
        and all(source_checks.values()) and all(definition_checks.values())
        and all(item["passes"] for item in selections.values())
    )
    return {
        "passed": passed,
        "bank_checks": bank_checks,
        "source_checks": source_checks,
        "definition_checks": definition_checks,
        "selection_checks": selections,
        "artifacts": [
            _artifact(experiment / "function_approximation_comparison.csv", source_root),
            _artifact(experiment / "native_task_comparison.csv", source_root),
            _artifact(experiment / "runtime_resource_accounting.csv", source_root),
        ],
    }


def _nonconvex_audit(source_root: Path) -> dict[str, Any]:
    results = source_root / "invariant_region_mlp/results/tangent_actor_accuracy_frontier"
    summary = _read_json(results / "FINAL_SUMMARY.json")
    rows_path = results / "nonseparable_three_seed_results.csv"
    rows = list(csv.DictReader(rows_path.open(newline="", encoding="utf-8")))
    seeds = sorted({int(row["seed"]) for row in rows})
    selected = [row for row in rows if row["selected"].lower() == "true"]
    post_files = [
        results / f"post_freeze_coupled_d100_k2_lam0p2_eta1_main_mcv_s{seed}.json"
        for seed in seeds
    ]
    hashes = []
    selection_ok = True
    for path in post_files:
        payload = _read_json(path)
        hashes.append(payload["test_bank_sha256"])
        selection_ok &= payload["ground_truth_used_for_training_or_selection"] is False
        selection_ok &= payload["ground_truth_loaded_only_after_selection_file"] is True
    passed = (
        seeds == [0, 1, 2] and len(selected) == 3 and len(set(hashes)) == 1
        and selection_ok and summary["nonseparable_support"] is True
        and summary["coupled_actor_finest_steps"] == 512
    )
    return {
        "passed": passed,
        "pde_match": {
            "lambda": 0.2, "eta": 1.0, "dimension": 100,
            "effective_dimension": 2, "horizon": 0.5,
        },
        "training_seeds": seeds,
        "selected_rows": len(selected),
        "test_bank_sha256": hashes[0] if hashes else None,
        "ground_truth_excluded_from_selection": selection_ok,
        "artifacts": [
            _artifact(rows_path, source_root),
            _artifact(results / "euler_convergence.csv", source_root),
            _artifact(results / "compute_accounting.csv", source_root),
            _artifact(results / "FINAL_SUMMARY.json", source_root),
        ],
    }


def _counterexample_audit(source_root: Path) -> dict[str, Any]:
    results = source_root / "invariant_region_mlp/results/counterexample_rescue"
    rows = list(csv.DictReader((results / "summary_all.csv").open(newline="", encoding="utf-8")))
    required = {20, 50, 100, 200, 400, 1000}
    present = {int(row["d"]) for row in rows if row["mode"] == "raw"}
    validation = _read_json(results / "independent_validation.json")
    reusable = sorted(required & present)
    missing = sorted(required - present)
    return {
        "passed": validation["n_passed"] >= 17 and set(reusable) == required - {400},
        "role": "secondary raw Full-History MLP dimensional-failure series",
        "reusable_dimensions": reusable,
        "missing_dimensions": missing,
        "method_cell": False,
        "reason_not_tangent_actor": "projection-based IR modes are not Tangent-Actor and remain excluded",
        "artifacts": [
            _artifact(results / "summary_all.csv", source_root),
            _artifact(results / "independent_validation.json", source_root),
        ],
    }


def _excluded_pilots(source_root: Path) -> list[dict[str, Any]]:
    return [
        {
            "path": "invariant_region_mlp/results/convex_concave_semilinear_extended",
            "status": "anchor_only",
            "reason": "2048 states, 16 label steps, 3x192 network, 2500 steps and only a0-a3 do not match the mature frozen budget",
        },
        {
            "path": "invariant_region_mlp/results/tangent_actor_strong_baselines",
            "status": "supplementary_only",
            "reason": "compute-matched amortized/adapted baseline architectures are not the authors' fixed official configurations",
        },
    ]


def _update_manifest(audit: dict[str, Any]) -> None:
    path = RESULTS_DIR / "EXECUTION_MANIFEST.json"
    manifest = _read_json(path)
    accepted = {
        ("hjb_log", 100, "tangent_actor"): ([0, 1, 2], "completed_reused"),
        ("hjb_log", 100, "sdgd_pinn"): ([0], "completed_reused"),
        ("hjb_log", 100, "deep_bsde"): ([0, 1, 2], "completed_reused"),
        ("ou_gmm", 100, "tangent_actor"): ([0, 1, 2], "completed_reused"),
        ("ou_gmm", 100, "dpi"): ([0], "completed_reused"),
        ("ou_gmm", 100, "pinn_hte"): ([0], "completed_reused"),
        ("allen_cahn", 100, "tangent_actor"): ([0, 1, 2], "completed_reused"),
        ("allen_cahn", 100, "deep_bsde"): ([0, 1, 2], "completed_reused"),
        ("nonconvex_coupled", 100, "tangent_actor"): ([0, 1, 2], "completed_reused"),
    }
    if not audit["original_100d"]["passed"] or not audit["nonconvex_100d"]["passed"]:
        raise RuntimeError("historical audit failed; refusing to mark reused cells complete")
    for cell in manifest["method_cells"]:
        key = (cell["pde"], int(cell["dimension"]), cell["method"])
        if key not in accepted:
            continue
        seeds, state = accepted[key]
        cell["execution_status"] = state
        cell["completed_seeds"] = seeds
        cell["artifacts"] = (
            audit["nonconvex_100d"]["artifacts"]
            if key[0] == "nonconvex_coupled"
            else audit["original_100d"]["artifacts"]
        )
        cell["attempts"] = [{"kind": "historical_reuse", "audit": "HISTORICAL_REUSE_AUDIT.json"}]
    manifest["historical_reuse_audit"] = "HISTORICAL_REUSE_AUDIT.json"
    manifest["status"] = "stage_a_reuse_audited"
    write_json(path, manifest)


def main() -> None:
    source_root = Path(__file__).resolve().parents[3]
    audit = {
        "schema_version": 1,
        "original_100d": _original_audit(source_root),
        "nonconvex_100d": _nonconvex_audit(source_root),
        "counterexample_full_history": _counterexample_audit(source_root),
        "excluded_or_anchor_only": _excluded_pilots(source_root),
        "rule": "reuse requires identical equation, dimension, test bank, estimator convention, and checkpoint-selection protocol",
    }
    audit["passed"] = (
        audit["original_100d"]["passed"]
        and audit["nonconvex_100d"]["passed"]
        and audit["counterexample_full_history"]["passed"]
    )
    write_json(RESULTS_DIR / "HISTORICAL_REUSE_AUDIT.json", audit)
    _update_manifest(audit)
    print(json.dumps({
        "passed": audit["passed"],
        "reused_method_cells": 9,
        "counterexample_missing": audit["counterexample_full_history"]["missing_dimensions"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
