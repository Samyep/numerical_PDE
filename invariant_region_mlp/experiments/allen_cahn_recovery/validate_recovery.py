"""Validate the Allen--Cahn recovery result bundle from raw artifacts.

The checks deliberately start from the task-level JSON files instead of
trusting the aggregate tables.  A machine-readable audit is written next to
the result bundle and the process exits nonzero if any acceptance condition
fails.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import subprocess
from typing import Any


HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[1]
REPOSITORY_ROOT = PROJECT_ROOT.parent
RESULT_ROOT = PROJECT_ROOT / "results" / "allen_cahn_recovery"
RAW_ROOT = RESULT_ROOT / "raw"
DOC_ROOT = PROJECT_ROOT / "docs"

EXPECTED_TASKS = {
    "equivalence": 240,
    "source": 75,
    "pilot": 78,
    "final": 270,
}
PRODUCTION_STAGES = ("source", "pilot", "final")
PRODUCTION_METHODS = {
    "raw",
    "beck_truncated",
    "interval_ir",
    "certified_ir_0_1",
}
FIGURES = (
    "raw_vs_truncated_error.pdf",
    "violation_vs_budget.pdf",
    "pathwise_equivalence.pdf",
    "error_vs_work.pdf",
)


def _finite_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and math.isfinite(float(value))


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _git_status_paths() -> list[str]:
    try:
        output = subprocess.check_output(
            ["git", "-C", str(REPOSITORY_ROOT), "status", "--porcelain=v1"],
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return []
    return [line[3:].replace("\\", "/") for line in output.splitlines() if line]


def main() -> None:
    checks: list[dict[str, Any]] = []

    def check(name: str, condition: bool, detail: Any) -> None:
        checks.append({"name": name, "passed": bool(condition), "detail": detail})

    payloads: dict[str, list[tuple[Path, dict[str, Any]]]] = {}
    manifests: dict[str, dict[str, Any]] = {}
    for stage, expected in EXPECTED_TASKS.items():
        files = sorted((RAW_ROOT / stage).glob("*.json"))
        rows = [(path, json.loads(path.read_text(encoding="utf-8"))) for path in files]
        payloads[stage] = rows
        check(
            f"raw task count: {stage}",
            len(rows) == expected,
            {"expected": expected, "observed": len(rows)},
        )

        manifest_path = RESULT_ROOT / f"{stage}_manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifests[stage] = manifest
        manifest_paths = {
            str(Path(item["path"])).replace("\\", "/")
            for item in manifest.get("tasks", [])
        }
        raw_paths = {
            str(path.relative_to(RESULT_ROOT)).replace("\\", "/")
            for path, _ in rows
        }
        check(
            f"complete manifest: {stage}",
            manifest.get("status") == "complete"
            and manifest.get("tasks_total") == expected
            and manifest.get("tasks_completed") == expected
            and manifest_paths == raw_paths,
            {
                "status": manifest.get("status"),
                "tasks_total": manifest.get("tasks_total"),
                "tasks_completed": manifest.get("tasks_completed"),
                "listed_paths_match": manifest_paths == raw_paths,
            },
        )

    all_rows = [item for rows in payloads.values() for item in rows]
    production = [
        item
        for stage in PRODUCTION_STAGES
        for item in payloads[stage]
    ]
    check("total raw tasks", len(all_rows) == 663, len(all_rows))
    check("production raw tasks", len(production) == 423, len(production))

    malformed_methods: list[str] = []
    invalid_task_metadata: list[str] = []
    nonfinite_fields: list[str] = []
    negative_counters: list[str] = []
    equivalence_failures: list[str] = []
    production_pairing_failures: list[str] = []
    activation_failures: list[str] = []
    correction_count = 0
    correction_exact = 0
    maximum_root_discrepancy = 0.0
    maximum_correction_discrepancy = 0.0
    maximum_production_state = -math.inf
    minimum_production_state = math.inf

    finite_fields = (
        "value",
        "wall_clock_seconds",
        "pre_truncation_violation_rate",
        "truncation_activation_rate",
        "mean_squared_overshoot",
        "pre_projection_min",
        "pre_projection_max",
        "generator_abs_before_mean",
        "generator_abs_after_mean",
        "nonlinear_correction_variance",
    )
    counter_fields = (
        "terminal_g_evals",
        "f_evals",
        "recursive_states",
        "normal_scalar_draws",
        "uniform_draws",
        "total_stochastic_samples",
        "state_transitions",
        "checked_states",
        "violating_states",
        "projection_activations",
        "nonfinite_states",
        "nonfinite_corrections",
        "generator_before_nonfinite",
        "generator_after_nonfinite",
    )

    for path, payload in all_rows:
        relative = str(path.relative_to(RESULT_ROOT)).replace("\\", "/")
        stage = path.parent.name
        task = payload.get("task", {})
        methods = payload.get("methods", [])
        expected_methods = (
            {"beck_truncated", "interval_ir"}
            if stage == "equivalence"
            else PRODUCTION_METHODS
        )
        observed_methods = {row.get("method") for row in methods}
        if observed_methods != expected_methods or len(methods) != len(expected_methods):
            malformed_methods.append(relative)
        if task.get("stage") != stage or any(
            row.get(key) != task.get(key)
            for row in methods
            for key in (
                "stage",
                "dimension",
                "depth",
                "sample_size",
                "repetition",
                "seed",
                "test_point",
                "radius",
                "radius_label",
            )
        ):
            invalid_task_metadata.append(relative)

        for row in methods:
            for field in finite_fields:
                if not _finite_number(row.get(field)):
                    nonfinite_fields.append(f"{relative}:{row.get('method')}:{field}")
            if not row.get("finite_output"):
                nonfinite_fields.append(f"{relative}:{row.get('method')}:finite_output")
            for field in counter_fields:
                value = row.get(field)
                if not isinstance(value, int) or value < 0:
                    negative_counters.append(f"{relative}:{row.get('method')}:{field}")

        equivalence = payload.get("equivalence", {})
        root_discrepancy = float(equivalence.get("absolute_value_discrepancy", math.inf))
        correction_discrepancy = float(
            equivalence.get("max_correction_discrepancy", math.inf)
        )
        maximum_root_discrepancy = max(maximum_root_discrepancy, root_discrepancy)
        maximum_correction_discrepancy = max(
            maximum_correction_discrepancy, correction_discrepancy
        )
        if not (
            root_discrepancy <= 1e-12
            and correction_discrepancy <= 1e-12
            and equivalence.get("exact_value_match") is True
            and equivalence.get("draw_fingerprints_match") is True
            and equivalence.get("work_counters_match") is True
        ):
            equivalence_failures.append(relative)

        if stage == "equivalence":
            correction_count += int(equivalence["correction_count"])
            correction_exact += int(equivalence["exact_correction_matches"])
        else:
            values = [float(row["value"]) for row in methods]
            fingerprints = {row["draw_fingerprint"] for row in methods}
            work_signatures = {
                tuple(row[field] for field in counter_fields[:7]) for row in methods
            }
            if not (
                len(set(values)) == 1
                and len(fingerprints) == 1
                and len(work_signatures) == 1
            ):
                production_pairing_failures.append(relative)
            for row in methods:
                minimum_production_state = min(
                    minimum_production_state, float(row["pre_projection_min"])
                )
                maximum_production_state = max(
                    maximum_production_state, float(row["pre_projection_max"])
                )
                if (
                    row["violating_states"] != 0
                    or row["projection_activations"] != 0
                    or row["pre_truncation_violation_rate"] != 0.0
                    or row["truncation_activation_rate"] != 0.0
                    or row["mean_squared_overshoot"] != 0.0
                    or row["nonfinite_states"] != 0
                    or row["nonfinite_corrections"] != 0
                    or row["generator_before_nonfinite"] != 0
                    or row["generator_after_nonfinite"] != 0
                ):
                    activation_failures.append(f"{relative}:{row['method']}")

    check("method sets and cardinalities", not malformed_methods, malformed_methods[:10])
    check("task metadata consistency", not invalid_task_metadata, invalid_task_metadata[:10])
    check("finite outputs and diagnostics", not nonfinite_fields, nonfinite_fields[:10])
    check("nonnegative integer work counters", not negative_counters, negative_counters[:10])
    check(
        "Beck versus IR pathwise tolerance",
        not equivalence_failures and maximum_root_discrepancy <= 1e-12,
        {
            "failures": equivalence_failures[:10],
            "max_root_discrepancy": maximum_root_discrepancy,
            "max_correction_discrepancy": maximum_correction_discrepancy,
        },
    )
    check(
        "dedicated saved corrections are bitwise exact",
        correction_count == 9408 and correction_exact == correction_count,
        {"exact": correction_exact, "total": correction_count},
    )
    check(
        "inactive production methods coincide pathwise",
        not production_pairing_failures,
        production_pairing_failures[:10],
    )
    check(
        "zero activation and no nonfinite production states",
        not activation_failures,
        {
            "failures": activation_failures[:10],
            "minimum_reused_state": minimum_production_state,
            "maximum_reused_state": maximum_production_state,
        },
    )
    check(
        "production states remain inside certified interval",
        minimum_production_state >= 0.0 and maximum_production_state <= 1.0,
        {"minimum": minimum_production_state, "maximum": maximum_production_state},
    )

    repetition_rows = _read_csv(RESULT_ROOT / "repetition_metrics.csv")
    equivalence_rows = _read_csv(RESULT_ROOT / "pathwise_equivalence.csv")
    summary_rows = _read_csv(RESULT_ROOT / "summary.csv")
    check("repetition CSV row count", len(repetition_rows) == 1692, len(repetition_rows))
    check("equivalence CSV row count", len(equivalence_rows) == 663, len(equivalence_rows))
    check("summary CSV row count", len(summary_rows) == 200, len(summary_rows))

    n5_means: dict[str, float] = {}
    published = {10: 0.29555, 100: 0.03373, 1000: 0.00340}
    for dimension, reference in published.items():
        values = [
            float(row["value"])
            for row in repetition_rows
            if row["stage"] == "source"
            and int(row["dimension"]) == dimension
            and int(row["depth"]) == 5
            and int(row["sample_size"]) == 5
            and row["method"] == "raw"
        ]
        n5_means[str(dimension)] = sum(values) / len(values)
    check(
        "source n=M=5 means reproduce published references",
        all(abs(n5_means[str(d)] - ref) < 1e-3 for d, ref in published.items()),
        {"observed": n5_means, "published_numerical_references": published},
    )

    figure_details: dict[str, Any] = {}
    figures_valid = True
    for name in FIGURES:
        path = RESULT_ROOT / "figures" / name
        valid = path.is_file() and path.stat().st_size > 1000
        if valid:
            valid = path.read_bytes()[:5] == b"%PDF-"
        figures_valid &= valid
        figure_details[name] = {
            "exists": path.is_file(),
            "bytes": path.stat().st_size if path.is_file() else 0,
            "pdf_header": valid,
        }
    check("four nonempty PDF figures", figures_valid, figure_details)

    report_path = DOC_ROOT / "ALLEN_CAHN_RECOVERY_REPORT.md"
    source_notes_path = DOC_ROOT / "ALLEN_CAHN_SOURCE_NOTES.md"
    report = report_path.read_text(encoding="utf-8")
    source_notes = source_notes_path.read_text(encoding="utf-8")
    check(
        "report records honest Verdict B",
        "VERDICT B" in report and "do not activate truncation" in report,
        str(report_path.relative_to(REPOSITORY_ROOT)).replace("\\", "/"),
    )
    check(
        "source notes distinguish theory and companion numerics",
        "1907.06729" in source_notes
        and "2005.10206" in source_notes
        and "does **not** contain a concrete finite-budget numerical" in source_notes
        and "fixed truncation radius r=4" in source_notes,
        str(source_notes_path.relative_to(REPOSITORY_ROOT)).replace("\\", "/"),
    )

    status_paths = _git_status_paths()
    manuscript_like = [
        path
        for path in status_paths
        if path.lower().endswith((".tex", ".bib"))
        or "/manuscript" in path.lower()
        or "/paper/" in path.lower()
    ]
    check("no manuscript files modified", not manuscript_like, manuscript_like)

    elapsed = sum(float(manifests[stage]["elapsed_seconds"]) for stage in EXPECTED_TASKS)
    task_seconds = sum(
        float(manifests[stage]["task_seconds_sum"]) for stage in EXPECTED_TASKS
    )
    passed = all(item["passed"] for item in checks)
    audit = {
        "schema_version": 1,
        "status": "passed" if passed else "failed",
        "verdict": "B",
        "checks_passed": sum(item["passed"] for item in checks),
        "checks_total": len(checks),
        "statistics": {
            "raw_tasks": len(all_rows),
            "production_method_repetitions": len(repetition_rows),
            "dedicated_equivalence_pairs": len(payloads["equivalence"]),
            "exact_saved_corrections": correction_exact,
            "saved_corrections_total": correction_count,
            "maximum_root_discrepancy": maximum_root_discrepancy,
            "maximum_correction_discrepancy": maximum_correction_discrepancy,
            "minimum_reused_state": minimum_production_state,
            "maximum_reused_state": maximum_production_state,
            "orchestration_wall_seconds": elapsed,
            "task_wall_seconds": task_seconds,
        },
        "checks": checks,
    }
    output = RESULT_ROOT / "validation_audit.json"
    output.write_text(json.dumps(audit, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({key: audit[key] for key in ("status", "checks_passed", "checks_total", "statistics")}, indent=2))
    if not passed:
        failed = [item["name"] for item in checks if not item["passed"]]
        raise SystemExit("validation failed: " + "; ".join(failed))


if __name__ == "__main__":
    main()
