"""Write the concise final Markdown report from audited tables."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import json
import subprocess

from .protocol import REPORT_PATH, RESULTS_DIR, SOURCE_ROOT, write_json_atomic


def _rows(name: str) -> list[dict[str, str]]:
    with (RESULTS_DIR / name).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_report() -> str:
    scaling = _rows("deep_bsde_scaling_summary.csv")
    accounting = json.loads((RESULTS_DIR / "resource_accounting.json").read_text(encoding="utf-8"))
    completeness = json.loads((RESULTS_DIR / "training_completeness.json").read_text(encoding="utf-8"))
    preflight = json.loads((RESULTS_DIR / "preflight_allocation_gate.json").read_text(encoding="utf-8"))
    matched = json.loads((RESULTS_DIR / "matched_existing_baselines_status.json").read_text(encoding="utf-8"))
    scope = json.loads((RESULTS_DIR / "scope_amendment.json").read_text(encoding="utf-8"))
    incidents = json.loads(
        (RESULTS_DIR / "infrastructure_incidents.json").read_text(encoding="utf-8")
    )
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=SOURCE_ROOT, text=True).strip()
    lines = [
        "# Deep BSDE hyperparameter-budget scaling (three seeds)",
        "",
        "## Outcome",
        "",
        f"- Formal terminal records: **{completeness['terminal_records']} / {completeness['expected_terminal_records']}**.",
        f"- Numerical failures: **{completeness['numerical_failure_count']}**.",
        f"- Allocation preflight: **{preflight['status']}**. Its {preflight['projected_total_gpu_hours_with_contingency']:.2f} H100-hour projection covered the original larger grid; the active scope was subsequently narrowed before truth opening.",
        f"- Actual allocation (including recorded failures/retries): **{accounting['total_allocated_h100_gpu_hours']:.2f} H100 GPU-hours**.",
        f"- Frozen DPI/PINN predictions on the paired 64-point subset: **{matched['status']}**, **{matched.get('row_count', 0)} / 108** seed-level artifacts (no retraining).",
        f"- Scope amendment: **C0/C1/C2 only**, {scope['expected_active_tasks']} tasks / {scope['expected_terminal_records']} terminal records. Joint J1/J2/J3 scaling is excluded from analysis and tables at the user's request before the truth gate opened.",
        "- Excluded work is retained for audit and charged to total allocation: one completed J1 bundle and four terminal records from an interrupted J3 bundle; no excluded accuracy is computed.",
        f"- Infrastructure incidents: **{len(incidents.get('incidents', []))}**. A CUDA-unavailable H100 node affected nine task launches; exact frozen-task retries exclude that node, and every failed/cancelled allocation is retained in accounting.",
        "",
        "Deep BSDE is a native scalar-initial-value method here: every one of the 64 query states receives a fresh model. Reported Deep BSDE compute is therefore the cumulative active training time for all 64 models, averaged across training seeds—not a one-query runtime.",
        "",
        "## Individual sweeps",
        "",
        "| PDE | Sweep | Level | Hyperparameter | RMSE mean ± seed SD | RelMAE | Total 64-query min | Historical TA-reference status |",
        "|---|---|---:|---|---:|---:|---:|---|",
    ]
    for row in scaling:
        if row["status"] == "complete":
            metric = f"{float(row['mean_rmse']):.6g} ± {float(row['std_rmse']):.3g}"
            relmae = f"{float(row['mean_relmae']):.5g}"
            minutes = f"{float(row['mean_total_active_training_minutes_64_queries']):.1f}"
        else:
            metric, relmae, minutes = "NUMERICAL INSTABILITY", "NA", "NA"
        lines.append(
            f"| {row['pde']} | {row['sweep']} | {row['level']} | "
            f"{row['parameter_name']}={row['parameter_value']} | {metric} | {relmae} | "
            f"{minutes} | {row['historical_ta_reference_target_status']} |"
        )
    lines.extend([
        "",
        "C2 changes Euler time-discretization steps; it is not an optimizer-step sweep. C0 checkpoints 1000/2000/4000 come from a single continuous trajectory. The common 2000/64/20 checkpoint is reused exactly for C0-level2, C1-level2, and C2-level2.",
    ])
    lines.extend([
        "",
        "## Interpretation boundaries",
        "",
        "- A lower observed RMSE at a larger setting is a genuine observed improvement only at the evaluated configurations; nonmonotonic sequences are reported without smoothing.",
        "- `TARGET_NOT_REACHED_WITHIN_EVALUATED_BUDGET` refers only to the pre-existing audited TA reference RMSE used to frame the reviewer question. That TA value was measured on 960 points, whereas Deep BSDE uses the frozen 64-point subset, so it is descriptive rather than a paired hypothesis test.",
        "- A failed 64-query configuration is numerical instability and is not recomputed on a finite subset.",
        "- The experiment does not establish that Deep BSDE cannot converge under arbitrarily larger budgets.",
        "- The main ICML table keeps historical DPI/PINN values unchanged and states the 960-vs-64 protocol limitation. The supplementary table uses their saved predictions on the identical 64 points when available, without retraining.",
        "",
        "## Provenance and outputs",
        "",
        f"- Report-generation commit: `{commit}`",
        "- Main table: `results/deep_bsde_budget_scaling_3seed/tables/deep_bsde_budget_scaling_main.tex`",
        "- Seed metrics: `results/deep_bsde_budget_scaling_3seed/deep_bsde_seed_metrics.csv`",
        "- Resource accounting: `results/deep_bsde_budget_scaling_3seed/resource_accounting.json`",
        "- Artifact audit: `results/deep_bsde_budget_scaling_3seed/artifact_audit.json`",
        "",
    ])
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text("\n".join(lines), encoding="utf-8")
    return str(REPORT_PATH)


if __name__ == "__main__":
    print(write_report())

