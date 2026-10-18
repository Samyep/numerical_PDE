"""Fail-closed audit of final frontier artifacts and scientific boundaries."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from .common import PACKAGE_DIR, RESULTS_DIR


REQUIRED_REPORTS = (
    "ACCURACY_FRONTIER_REPORT.md",
    "HIGH_PRECISION_EULER_AUDIT.md",
    "STRONG_FULL_HISTORY_COMPARISON.md",
    "NONSEPARABLE_HIGH_BUDGET_REPORT.md",
    "AMBIENT_400D_GENERALIZATION_AUDIT.md",
    "REFERENCE_FREE_ERROR_THEORY_NOTE.md",
)
REQUIRED_TABLES = (
    "accuracy_frontier.csv",
    "euler_convergence.csv",
    "full_history_extended_grid.csv",
    "nonseparable_three_seed_results.csv",
    "ambient_scaling_comparison.csv",
    "compute_accounting.csv",
)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def run() -> dict:
    missing = [name for name in (*REQUIRED_REPORTS, *REQUIRED_TABLES, "reproducibility_manifest.json", "FINAL_SUMMARY.json") if not (RESULTS_DIR / name).exists()]
    if missing:
        raise RuntimeError(f"missing final artifacts: {missing}")
    precision_all = sorted((RESULTS_DIR / "raw").glob("precision_primary_a4_seed*.json"))
    precision_by_seed = {}
    for path in precision_all:
        payload = json.loads(path.read_text(encoding="utf-8"))
        seed = int(payload["training_seed"])
        paths = int(payload["configuration"]["paths"])
        if seed not in precision_by_seed or paths > precision_by_seed[seed][0]:
            precision_by_seed[seed] = (paths, path)
    precision = [precision_by_seed[seed][1] for seed in sorted(precision_by_seed)]
    if len(precision) != 3:
        raise RuntimeError(f"expected three primary precision seeds, found {len(precision)}")
    hashes, primary_rows = set(), 0
    for path in precision:
        payload = json.loads(path.read_text(encoding="utf-8"))
        hashes.add(payload["test_bank_sha256"])
        if not payload["checkpoint_selection_reused_unchanged"]:
            raise RuntimeError(path)
        if not payload["ground_truth_used_for_training_or_selection"] is False:
            raise RuntimeError(path)
        if [row["steps"] for row in payload["rows"]] != [64, 128, 256, 512]:
            raise RuntimeError(path)
        if int(payload["configuration"]["paths"]) < 8192:
            raise RuntimeError(f"insufficient primary precision paths: {path}")
        primary_rows += len(payload["rows"])
    if len(hashes) != 1:
        raise RuntimeError("primary precision test banks differ")
    coupled = sorted(RESULTS_DIR.glob("post_freeze_coupled_d100*_main_mcv_s*.json"))
    if len(coupled) != 3:
        raise RuntimeError(f"expected three coupled main runs, found {len(coupled)}")
    coupled_hashes = set()
    for path in coupled:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload["ground_truth_used_for_training_or_selection"]:
            raise RuntimeError(path)
        if len(payload["iteration_rows"]) != 5:
            raise RuntimeError(path)
        coupled_hashes.add(payload["test_bank_sha256"])
    if len(coupled_hashes) != 1:
        raise RuntimeError("coupled Actor test banks differ")
    endpoint_seeds = {}
    for path in sorted((RESULTS_DIR / "raw").glob("full_history_*.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload["configuration_role"] != "endpoint":
            continue
        key = (payload["spec"]["family"], int(payload["n"]), int(payload["M"]))
        endpoint_seeds.setdefault(key, set()).add(int(payload["baseline_seed"]))
    required_endpoints = {
        ("primary", 3, 384), ("primary", 3, 512), ("primary", 4, 96),
        ("coupled", 3, 512), ("coupled", 4, 96),
    }
    for endpoint in required_endpoints:
        if endpoint_seeds.get(endpoint) != {0, 1, 2}:
            raise RuntimeError(
                f"Full-History endpoint {endpoint} has seeds {sorted(endpoint_seeds.get(endpoint, set()))}"
            )
    ambient = sorted(RESULTS_DIR.glob("post_freeze_primary_d400*_main_mcv_s0.json"))
    if len(ambient) != 1:
        raise RuntimeError(f"expected one matched-budget d400 seed-0 run, found {len(ambient)}")
    ambient_payload = json.loads(ambient[0].read_text(encoding="utf-8"))
    a0_rmse = float(ambient_payload["iteration_rows"][0]["mcv_value_rmse"])
    selected_rmse = float(ambient_payload["selected_row"]["mcv_value_rmse"])
    ambient_improvement = (a0_rmse - selected_rmse) / a0_rmse
    ambient_training = json.loads(
        (RESULTS_DIR / f"training_primary_d400_k2_lam0p2_eta0_main_mcv_s0.json").read_text(encoding="utf-8")
    )
    protocol = json.loads((PACKAGE_DIR / "protocol.json").read_text(encoding="utf-8"))
    gap_rule = protocol["ambient_400"]["persistent_generalization_gap_definition"]
    gap_phases = 0
    for phase in ambient_training["phases"]:
        fit = phase["training"]
        train = float(fit["selected_train_gradient_vector_rmse"])
        valid = float(fit["selected_validation_gradient_vector_rmse"])
        replication = float(phase["labels"]["gradient_label_replication_vector_rmse"])
        if (
            valid >= float(gap_rule["validation_over_train_gradient_rmse_ratio"]) * max(train, 1e-300)
            and valid - train >= float(gap_rule["minimum_excess_as_fraction_of_label_replication_rmse"]) * replication
        ):
            gap_phases += 1
    persistent_gap = gap_phases >= int(gap_rule["minimum_update_phases"])
    if ambient_improvement >= float(protocol["ambient_400"]["meaningful_value_improvement_fraction"]):
        ambient_replicates = sorted(RESULTS_DIR.glob("post_freeze_primary_d400*_main_mcv_s*.json"))
        if len(ambient_replicates) != 3:
            raise RuntimeError("meaningful d400 seed-0 result was not replicated with seeds 1 and 2")
    elif persistent_gap:
        low_rank = sorted(RESULTS_DIR.glob("post_freeze_primary_d400*_lowrank_r*.json"))
        ranks = {int(json.loads(path.read_text(encoding="utf-8"))["rank"]) for path in low_rank}
        if ranks != set(protocol["ambient_400"]["learned_low_rank_ranks"]):
            raise RuntimeError(f"conditional low-rank grid incomplete: {sorted(ranks)}")
    for table in REQUIRED_TABLES:
        rows = _read_csv(RESULTS_DIR / table)
        if not rows:
            raise RuntimeError(f"empty table {table}")
        for row in rows:
            for key, value in row.items():
                if value and value.lower() in {"nan", "inf", "-inf"}:
                    raise RuntimeError(f"nonfinite marker in {table}: {key}={value}")
    reference = json.loads((RESULTS_DIR / "references" / "primary_refinement_audit.json").read_text())
    if reference["reference_uncertainty_value_rmse"] >= 1e-5:
        raise RuntimeError("reference uncertainty did not reach the precision gate")
    if not reference["inside_invariant_interval"]:
        raise RuntimeError("refined reference left the invariant interval")
    manifest = json.loads((RESULTS_DIR / "reproducibility_manifest.json").read_text())
    if manifest["oracle_subspace_used"] or manifest["historical_results_overwritten"]:
        raise RuntimeError("manifest boundary violation")
    summary = json.loads((RESULTS_DIR / "FINAL_SUMMARY.json").read_text())
    result = {
        "passed": True,
        "primary_precision_runs": len(precision),
        "primary_precision_rows": primary_rows,
        "coupled_training_runs": len(coupled),
        "matched_budget_d400_runs": len(ambient),
        "d400_seed0_value_improvement_fraction": ambient_improvement,
        "d400_persistent_generalization_gap": persistent_gap,
        "shared_primary_test_bank_sha256": next(iter(hashes)),
        "reference_uncertainty_value_rmse": reference["reference_uncertainty_value_rmse"],
        "decision": summary["decision"],
    }
    (RESULTS_DIR / "FINAL_ARTIFACT_AUDIT.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
