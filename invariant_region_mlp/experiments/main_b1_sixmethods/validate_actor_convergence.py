"""Validate the completed protocol-v3 actor convergence study."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable

from .common import PROJECT_DIR, RESULTS_DIR, atomic_json, load_config


PROTOCOL_DIR = RESULTS_DIR / "actor_convergence" / "protocol_v3"
RAW_DIR = PROTOCOL_DIR / "raw"
CHECKPOINT_DIR = PROTOCOL_DIR / "checkpoints"
PREDICTION_DIR = PROTOCOL_DIR / "predictions"
VALIDATION_PATH = PROTOCOL_DIR / "validation.json"


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _finite(values: Iterable[Any]) -> bool:
    try:
        return all(math.isfinite(float(value)) for value in values)
    except (TypeError, ValueError):
        return False


def _geometric_mean(values: Iterable[float]) -> float:
    data = list(values)
    return math.exp(sum(math.log(value) for value in data) / len(data))


def _close(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=1.0e-10, abs_tol=1.0e-12)


def validate() -> dict[str, Any]:
    config = load_config()
    expected_cells = {
        (int(d), int(k))
        for d in config["pde"]["dimensions"]
        for k in config["pde"]["effective_dimensions"]
    }
    expected_repetitions = int(config["stage_b"]["mlp_seeds"])
    failures: list[str] = []
    selected_levels: dict[str, int] = {}
    test_hashes: dict[str, str] = {}
    rows_checked = 0

    for d, k in sorted(expected_cells):
        cell = f"d{d}_k{k}"
        path = RAW_DIR / f"actor_converged_{cell}.json"
        if not path.is_file():
            failures.append(f"{cell}: missing converged raw record")
            continue
        payload = _load(path)
        if payload.get("status") != "complete":
            failures.append(f"{cell}: status={payload.get('status')}")
        if payload.get("d") != d or payload.get("k") != k:
            failures.append(f"{cell}: dimension metadata mismatch")
        if payload.get("paired_stage_b_seed_namespace") != "B1|actor_mlp":
            failures.append(f"{cell}: wrong paired Stage B RNG namespace")
        if payload.get("protocol", {}).get("protocol_version") != 3:
            failures.append(f"{cell}: wrong protocol version")
        for flag in (
            "nonterminal_oracle_reference_used_for_training",
            "nonterminal_oracle_reference_used_for_selection",
        ):
            if payload.get(flag) is not False:
                failures.append(f"{cell}: {flag} is not false")
        if payload.get("exact_reference_first_used_after_checkpoint_frozen") is not True:
            failures.append(f"{cell}: exact reference timing is invalid")

        diagnostic = payload.get("postselection_exact_diagnostic", {})
        for flag in (
            "used_for_training",
            "used_for_early_stopping",
            "used_for_budget_convergence",
            "used_for_checkpoint_selection",
        ):
            if diagnostic.get(flag) is not False:
                failures.append(f"{cell}: post-selection diagnostic {flag} is not false")

        subspace = payload.get("label_subspace", {})
        if (
            subspace.get("status") != "complete"
            or subspace.get("rank") != k
            or subspace.get("rank_check_passed") is not True
            or subspace.get("nonterminal_oracle_reference_used") is not False
        ):
            failures.append(f"{cell}: invalid data-derived label subspace certificate")
        residual = subspace.get("relative_projection_residual")
        if not _finite([residual]) or float(residual) > 1.0e-10:
            failures.append(f"{cell}: label-subspace residual={residual}")
        if payload.get("audit", {}).get("nonterminal_oracle_reference_used") is not False:
            failures.append(f"{cell}: stochastic audit used a nonterminal oracle")

        convergence = payload.get("budget_convergence", {})
        levels = payload.get("levels", [])
        if convergence.get("converged") is not True:
            failures.append(f"{cell}: convergence certificate is false")
        if int(convergence.get("stable_transitions", -1)) < 2:
            failures.append(f"{cell}: fewer than two stable transitions")
        if len(levels) < 3:
            failures.append(f"{cell}: insufficient budget levels")
        else:
            selected = levels[-1].get("budget", {}).get("name")
            if payload.get("selected_level") != selected:
                failures.append(f"{cell}: selected level is not the final fitted level")
            selected_levels[str(selected)] = selected_levels.get(str(selected), 0) + 1
            tolerance = payload["protocol"]["convergence"]
            audit_tolerance = float(
                tolerance["audit_excess_rmse_relative_change_tolerance"]
            )
            prediction_tolerance = float(
                tolerance["prediction_relative_change_tolerance"]
            )
            final_stable: list[bool] = []
            for previous, current in zip(levels[-3:-1], levels[-2:]):
                previous_audit = float(previous["audit_noise_corrected_rmse"])
                current_audit = float(current["audit_noise_corrected_rmse"])
                audit_change = abs(current_audit - previous_audit) / max(
                    previous_audit, float.fromhex("0x1.0p-1022")
                )
                final_stable.append(
                    current.get("optimizer_converged") is True
                    and audit_change <= audit_tolerance
                    and float(current["prediction_relative_change"])
                    <= prediction_tolerance
                )
            if final_stable != [True, True]:
                failures.append(f"{cell}: final two transitions are not both stable")

        for index, level in enumerate(levels, start=1):
            if level.get("status") != "complete":
                failures.append(f"{cell}: level {index} is incomplete")
            for flag in (
                "nonterminal_oracle_reference_used_for_training",
                "nonterminal_oracle_reference_used_for_selection",
            ):
                if level.get(flag) is not False:
                    failures.append(f"{cell}: level {index} {flag} is not false")

        rows = payload.get("rows", [])
        if len(rows) != expected_repetitions:
            failures.append(
                f"{cell}: Stage B rows={len(rows)}, expected={expected_repetitions}"
            )
        repetitions = {int(row.get("repetition", -1)) for row in rows}
        if repetitions != set(range(expected_repetitions)):
            failures.append(f"{cell}: Stage B repetition set is invalid")
        cell_hashes = {str(payload.get("test_point_sha256"))}
        for row in rows:
            rows_checked += 1
            cell_hashes.add(str(row.get("test_point_sha256")))
            if row.get("method") != "actor_mlp_converged":
                failures.append(f"{cell}: wrong Stage B method label")
            if int(row.get("test_points", -1)) != 960:
                failures.append(f"{cell}: Stage B test point count is not 960")
            if row.get("test_distribution") != "x~N(0,I_d)":
                failures.append(f"{cell}: wrong Stage B test distribution")
            if row.get("value_nonfinite") != 0 or row.get("gradient_nonfinite") != 0:
                failures.append(f"{cell}: non-finite Stage B prediction")
            if not _finite(
                row.get(key)
                for key in ("value_rmse", "value_bias", "gradient_rmse")
            ):
                failures.append(f"{cell}: non-finite Stage B scalar metric")

        comparison_paths = (
            RESULTS_DIR / "raw" / "stage_b" / f"actor_mlp_{cell}.json",
            RESULTS_DIR
            / "actor_bottleneck_diagnostics"
            / "raw"
            / f"oracle_actor_{cell}.json",
            RESULTS_DIR
            / "actor_bottleneck_diagnostics"
            / "raw"
            / f"pathwise_baseline_no_actor_{cell}.json",
        )
        for comparison_path in comparison_paths:
            if not comparison_path.is_file():
                failures.append(f"{cell}: missing {comparison_path.name}")
                continue
            comparison = _load(comparison_path)
            cell_hashes.update(
                str(row.get("test_point_sha256"))
                for row in comparison.get("rows", [])
            )
        if len(cell_hashes) != 1:
            failures.append(f"{cell}: paired held-out hashes do not match")
        else:
            test_hashes[cell] = next(iter(cell_hashes))

        for artifact in (
            CHECKPOINT_DIR / f"actor_converged_{cell}.pt",
            PREDICTION_DIR / f"actor_converged_{cell}.npz",
        ):
            if not artifact.is_file() or artifact.stat().st_size == 0:
                failures.append(f"{cell}: missing local artifact {artifact.name}")

    comparison_path = PROTOCOL_DIR / "actor_convergence_comparison.csv"
    with comparison_path.open(newline="") as handle:
        comparison_rows = list(csv.DictReader(handle))
    comparison_cells = {(int(row["d"]), int(row["k"])) for row in comparison_rows}
    if comparison_cells != expected_cells or len(comparison_rows) != len(expected_cells):
        failures.append("comparison CSV does not contain exactly the 15 expected cells")
    ratio_old = [float(row["converged_over_original"]) for row in comparison_rows]
    ratio_oracle = [float(row["converged_over_oracle"]) for row in comparison_rows]
    ratio_no_actor = [float(row["converged_over_no_actor"]) for row in comparison_rows]
    actor_ratios = [float(row["actor_error_new_over_old"]) for row in comparison_rows]
    old_errors = [
        float(row["original_actor_relative_error_same_audit"])
        for row in comparison_rows
    ]
    new_errors = [
        float(row["converged_actor_relative_error_same_audit"])
        for row in comparison_rows
    ]
    if not _finite(ratio_old + ratio_oracle + ratio_no_actor + actor_ratios):
        failures.append("comparison CSV contains non-finite ratios")

    summary = _load(PROTOCOL_DIR / "summary.json")
    expected_summary = {
        "cells": len(expected_cells),
        "converged_cells": len(expected_cells),
        "converged_better_than_original_cells": sum(value < 1.0 for value in ratio_old),
        "converged_better_than_no_actor_cells": sum(
            value < 1.0 for value in ratio_no_actor
        ),
        "converged_within_10pct_oracle_cells": sum(
            value <= 1.10 for value in ratio_oracle
        ),
        "actor_error_improved_cells": sum(value < 1.0 for value in actor_ratios),
    }
    for key, expected in expected_summary.items():
        if int(summary.get(key, -1)) != expected:
            failures.append(f"summary {key}={summary.get(key)}, expected={expected}")
    summary_floats = {
        "geometric_mean_converged_over_original": _geometric_mean(ratio_old),
        "geometric_mean_converged_over_oracle": _geometric_mean(ratio_oracle),
        "geometric_mean_converged_over_no_actor": _geometric_mean(ratio_no_actor),
        "geometric_mean_actor_error_new_over_old": _geometric_mean(actor_ratios),
        "original_actor_relative_error_median": sorted(old_errors)[len(old_errors) // 2],
        "converged_actor_relative_error_median": sorted(new_errors)[len(new_errors) // 2],
        "converged_actor_relative_error_max": max(new_errors),
    }
    for key, expected in summary_floats.items():
        actual = float(summary.get(key, math.nan))
        if not _close(actual, expected):
            failures.append(f"summary {key}={actual}, expected={expected}")
    if summary.get("selected_level_counts") != selected_levels:
        failures.append("summary selected-level counts do not match raw cells")
    for flag in (
        "analytic_reference_used_for_training",
        "analytic_reference_used_for_selection",
    ):
        if summary.get(flag) is not False:
            failures.append(f"summary {flag} is not false")

    accounting = _load(PROTOCOL_DIR / "execution_accounting.json")
    records = accounting.get("records", [])
    accounted_seconds = sum(float(record.get("gpu_seconds", math.nan)) for record in records)
    if (
        len(records) != 19
        or any(record.get("state") != "COMPLETED" for record in records)
        or any(record.get("exit_code") != "0:0" for record in records)
        or not _close(accounted_seconds, float(accounting.get("actual_gpu_seconds", math.nan)))
        or not _close(accounted_seconds, float(summary.get("scheduler_actual_gpu_seconds", math.nan)))
    ):
        failures.append("scheduler accounting is inconsistent")

    manifest = _load(PROTOCOL_DIR / "artifact_manifest.json")
    manifest_checked = 0
    for record in manifest.get("artifacts", []):
        artifact = PROJECT_DIR / record["path"]
        if not artifact.is_file():
            failures.append(f"manifest artifact missing: {record['path']}")
            continue
        manifest_checked += 1
        digest = hashlib.sha256(artifact.read_bytes()).hexdigest()
        if digest != record.get("sha256") or artifact.stat().st_size != record.get("bytes"):
            failures.append(f"manifest mismatch: {record['path']}")

    result = {
        "schema_version": 1,
        "status": "PASS" if not failures else "FAIL",
        "cells_checked": len(expected_cells),
        "stage_b_rows_checked": rows_checked,
        "manifest_artifacts_checked": manifest_checked,
        "scheduler_allocations_checked": len(records),
        "selected_level_counts": selected_levels,
        "test_point_hashes": test_hashes,
        "failures": failures,
    }
    atomic_json(VALIDATION_PATH, result)
    if failures:
        raise RuntimeError("Actor convergence validation failed:\n" + "\n".join(failures))
    return result


def main() -> None:
    print(json.dumps(validate(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
