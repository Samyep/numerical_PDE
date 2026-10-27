"""Render the final publication-oriented Markdown report from audited artifacts."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from .protocol import RESULTS_DIR, SOURCE_ROOT, TARGETS


REPORT_PATH = (
    SOURCE_ROOT / "invariant_region_mlp" / "docs" / "SINGLE_POINT_BREAK_EVEN_REPORT.md"
)


def _read(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _fmt(value: Any, *, digits: int = 4) -> str:
    if value is None:
        return "--"
    if isinstance(value, float):
        return f"{value:.{digits}g}"
    return str(value)


def build() -> Path:
    analysis = _read(RESULTS_DIR / "analysis_summary.json")
    frozen = _read(RESULTS_DIR / "FROZEN_SELECTIONS.json")
    reference = _read(RESULTS_DIR / "reference_audit.json")
    preflight = _read(RESULTS_DIR / "preflight_audit.json")
    bank = _read(RESULTS_DIR / "bank_manifest.json")
    artifact = _read(RESULTS_DIR / "artifact_audit.json") if (
        RESULTS_DIR / "artifact_audit.json"
    ).is_file() else {"status": "PENDING"}
    timing_coverage = _read(RESULTS_DIR / "timing_coverage.json")
    resource_path = RESULTS_DIR / "resource_accounting.csv"
    resources = []
    if resource_path.is_file():
        with resource_path.open(newline="", encoding="utf-8") as handle:
            resources = list(csv.DictReader(handle))
    gpu_hours = sum(
        float(row["task_wall_seconds"]) * float(row.get("gpu_count") or 0)
        for row in resources
    ) / 3600.0
    failed_tasks = sum(row["status"] != "complete" for row in resources)

    outcome_rows = []
    for target_name, threshold in TARGETS.items():
        for method, label in (
            ("tangent_actor", "Tangent-Actor"),
            ("fh_mlp_standard", "Full-History MLP"),
            ("deep_bsde", "Deep BSDE"),
        ):
            item = analysis["targets"][target_name][method]
            if method == "tangent_actor":
                crossing = (
                    "reference curve"
                    if timing_coverage["targets"][target_name][method] == "RUN"
                    else timing_coverage["targets"][target_name][method]
                )
            else:
                cross = analysis["break_even"][target_name].get(method, {})
                if cross.get("status") == "OBSERVED_BREAK_EVEN":
                    crossing = f"({cross.get('previous_queries')}, {cross['crossing_queries']}] queries"
                else:
                    crossing = cross.get("status", "not evaluated")
            outcome_rows.append(
                f"| {target_name} ({threshold:.1e}) | {label} | "
                f"{item.get('selected_candidate') or '--'} | {item['status']} | "
                f"{_fmt(item.get('final_rmse'))} | {crossing} |"
            )

    selection_lines = []
    for target_name in TARGETS:
        values = []
        for method in ("tangent_actor", "fh_mlp_standard", "deep_bsde"):
            item = frozen["targets"][target_name][method]
            values.append(
                f"`{method}` = `{item.get('selected_candidate')}` "
                f"(calibration RMSE {_fmt(item.get('calibration_rmse', item.get('best_calibration_rmse')))})"
            )
        selection_lines.append(f"- **{target_name}:** " + "; ".join(values) + ".")

    report = f"""# Matched-Accuracy Single-Point Amortization Break-Even

## Outcome

| Target | Method | Frozen configuration | Final status | Final pooled RMSE | Measured break-even |
|---|---|---|---|---:|---|
{chr(10).join(outcome_rows)}

The headline crossing is the first **actually executed** common query-count
prefix at which mean baseline time is no smaller than mean Tangent-Actor total
time.  No line fitted to one query and no interpolated crossing is used.
The strict-target accuracy result is retained even when its runtime comparison
is marked `NOT_EVALUATED_COMPUTE_CAP`; no reduced-replicate curve is substituted.

## Frozen calibration choices

{chr(10).join(selection_lines)}

Calibration used 32 new queries.  Final accuracy used a disjoint 64-query
bank and could not change any configuration.  Tangent-Actor is trained once;
Full-History MLP reruns its recursion at every point; native Deep BSDE trains
a new scalar-$Y_0$ model at every point.

## Reference and correctness gates

- Reference gate: **{reference['status']}**.
- Reference uncertainty RMS: `{reference['uncertainty_rms']:.6e}` (limit `5e-6`).
- Maximum reference uncertainty: `{reference['uncertainty_max']:.6e}` (limit `1e-5`).
- 200-state PDE/terminal audit: **{reference['pde_and_terminal_audit']['status']}**.
- H100 numerical/source preflight: **{preflight['status']}** on `{preflight['gpu']}`.
- Final artifact audit: **{artifact['status']}**.

The frozen prose contained one unavoidable set-theoretic conflict: it required
both calibration and final banks (and the previous bank) to contain the exact
origin while also calling all banks disjoint.  The mandated origin is retained
and reported as the sole overlap; every stochastic row is disjoint.  Bank
audit status: `{bank['separation_audit']['status']}`.

## Timing and accounting

Primary times are synchronized method-active times on NVIDIA H100 80GB GPUs.
Queue time, reference construction, truth evaluation, plotting, and artifact
upload are excluded from the curves.  The task records currently account for
`{gpu_hours:.3f}` allocated H100 task-hours across `{len(resources)}` tasks,
including `{failed_tasks}` recorded failed tasks.  Detailed rows are in
`results/single_point_break_even/resource_accounting.csv`.

For Tangent-Actor, every plotted total contains its measured Actor training,
reference-free selection, checkpoint loading, and the executed query prefix.
For Full-History MLP it contains all recursive inference.  For Deep BSDE it
contains all fresh per-query training.  The timing bank has 524,288 Gaussian
states and deliberately contains no truth array.

## Provenance

- Preregistration commit: `498589cc2`.
- Frozen implementation commit: `{frozen['code_commit']}`.
- Frozen selection hash: `{frozen['canonical_sha256']}`.
- Study revision: `{analysis['revision']}`.
- Figure: `results/single_point_break_even/figures/matched_accuracy_break_even.pdf`.

The optional second PDE is not part of this report.  It may begin only after
this Allen--Cahn artifact audit is validated and under a separately frozen
addendum.
"""
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    temporary = REPORT_PATH.with_suffix(REPORT_PATH.suffix + ".tmp")
    temporary.write_text(report, encoding="utf-8")
    temporary.replace(REPORT_PATH)
    return REPORT_PATH


if __name__ == "__main__":
    print(build())

