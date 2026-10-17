"""Fail-closed audit of final frontier artifacts and scientific boundaries."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from .common import RESULTS_DIR


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
        primary_rows += len(payload["rows"])
    if len(hashes) != 1:
        raise RuntimeError("primary precision test banks differ")
    coupled = sorted(RESULTS_DIR.glob("post_freeze_coupled_d100*_main_mcv_s*.json"))
    if len(coupled) != 3:
        raise RuntimeError(f"expected three coupled main runs, found {len(coupled)}")
    for path in coupled:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload["ground_truth_used_for_training_or_selection"]:
            raise RuntimeError(path)
        if len(payload["iteration_rows"]) != 5:
            raise RuntimeError(path)
    ambient = sorted(RESULTS_DIR.glob("post_freeze_primary_d400*_main_mcv_s0.json"))
    if len(ambient) != 1:
        raise RuntimeError(f"expected one matched-budget d400 run, found {len(ambient)}")
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
