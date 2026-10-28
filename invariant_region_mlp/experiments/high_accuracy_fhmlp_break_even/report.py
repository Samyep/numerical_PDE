"""Generate the concise scientific report without changing experiment state."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess

from .protocol import RESULTS_DIR, SOURCE_ROOT, TARGETS


REPORT = SOURCE_ROOT / "invariant_region_mlp" / "docs" / "HIGH_ACCURACY_FHMLP_BREAK_EVEN_REPORT.md"


def _read(name: str):
    return json.loads((RESULTS_DIR / name).read_text(encoding="utf-8"))


def _candidate(item: dict) -> str:
    return str(item.get("selected_candidate") or "--")


def generate() -> Path:
    frozen = _read("FROZEN_SELECTIONS.json")
    final = _read("final_accuracy_summary.json")
    analysis = _read("analysis_summary.json")
    resources = _read("resource_accounting.json")
    audit = _read("artifact_audit.json")
    preflight = _read("preflight_audit.json")
    commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=SOURCE_ROOT, text=True
    ).strip()

    lines = [
        "# High-Accuracy FH-MLP Amortization Break-Even",
        "",
        "## Outcome",
        "",
        "| Target | TA-direct final RMSE | TA-FK-MCV final RMSE | FH-MLP config | FH final RMSE | Measured break-even |",
        "|---|---:|---:|---|---:|---|",
    ]
    for target, threshold in TARGETS.items():
        local = final["targets"][target]
        crossing = analysis["break_even"][target]
        direct_rmse = local["ta_direct"].get("final_rmse")
        mcv_rmse = local["ta_fk_mcv"].get("final_rmse")
        fh_rmse = local["fh_mlp_standard"].get("final_rmse")
        if crossing["status"] == "OBSERVED_BREAK_EVEN":
            outcome = f"({crossing['previous_queries']:,}, {crossing['crossing_queries']:,}]"
        else:
            outcome = crossing["status"]
        lines.append(
            f"| {target} ({threshold:.3g}) | "
            f"{('--' if direct_rmse is None else format(direct_rmse, '.6g'))} | "
            f"{('--' if mcv_rmse is None else format(mcv_rmse, '.6g'))} | "
            f"{_candidate(local['fh_mlp_standard'])} | "
            f"{('--' if fh_rmse is None else format(fh_rmse, '.6g'))} | {outcome} |"
        )

    lines += [
        "",
        "The previous audited RMSE <= 0.06 result remains unchanged: its observed bracket was (524,288, 1,048,576] queries. The new brackets above use only actually executed common prefixes; no interpolation or single-query extrapolation enters the verdict.",
        "",
        "## Frozen calibration",
        "",
        "| Target | Method | Calibration status | Selected configuration | Calibration RMSE |",
        "|---|---|---|---|---:|",
    ]
    for target in TARGETS:
        for method in ("ta_direct", "ta_fk_mcv", "fh_mlp_standard"):
            item = frozen["targets"][target][method]
            rmse = item.get("calibration_rmse", item.get("best_calibration_rmse"))
            lines.append(
                f"| {target} | {method} | {item['status']} | {_candidate(item)} | "
                f"{('--' if rmse is None else format(rmse, '.6g'))} |"
            )

    lines += [
        "",
        "## Measured cumulative timing",
        "",
    ]
    for target in TARGETS:
        crossing = analysis["break_even"][target]
        if crossing["status"] == "OBSERVED_BREAK_EVEN":
            lines += [
                f"- **{target}:** observed bracket `({crossing['previous_queries']:,}, {crossing['crossing_queries']:,}]`; at the first crossing TA used `{crossing['ta_mean_seconds_at_crossing'] / 60.0:.3f} +/- {crossing['ta_std_seconds_at_crossing'] / 60.0:.3f}` minutes and FH-MLP used `{crossing['fh_mean_seconds_at_crossing'] / 60.0:.3f} +/- {crossing['fh_std_seconds_at_crossing'] / 60.0:.3f}` minutes.",
                f"  At `{crossing['representative_beyond_queries']:,}` queries, TA used `{100.0 * crossing['ta_relative_computation_saved_beyond']:.1f}%` less measured cumulative computation.",
            ]
        else:
            lines.append(f"- **{target}:** `{crossing['status']}`.")

    lines += [
        "",
        "## Correctness, accounting, and interpretation",
        "",
        f"- Artifact audit: **{audit['status']}**.",
        f"- Actual allocation: **{resources['total_allocated_h100_gpu_hours']:.3f} / 30 H100 GPU-hours**, including failed allocations and retries recorded by Slurm.",
        f"- Reference refinement RMSE: `{_read('reference_audit.json')['fresh_bank_value_refinement_rmse']:.3e}`.",
        f"- FH batching audit: **{preflight['fh_batched_execution_audit']['status']}**. {preflight['fh_batched_execution_audit']['reason']}",
        "- TA training and reference-free checkpoint selection are charged in every cumulative curve even though byte-identical audited checkpoints were reused; new experimental allocation is reported separately.",
        "- The main timing curve uses `TA-direct`. The current manuscript describes FK-MCV as the final evaluation estimator, so FK-MCV accuracy and timing are explicitly separated rather than relabelled as direct inference.",
        "- A missing crossing means `NOT_REACHED_WITHIN_TESTED_BUDGET`; it is not converted into an interpolated claim.",
        "",
        "## Provenance",
        "",
        "- Preregistration commit: `60c04015682a5f77ff03fd829f575d059620a1bd`.",
        f"- Report-generation commit: `{commit}`.",
        f"- Frozen selection hash: `{frozen['canonical_sha256']}`.",
        "- Raw timing: `results/high_accuracy_fhmlp_break_even/timing_replicates.csv`.",
        "- Main figure: `results/high_accuracy_fhmlp_break_even/figures/high_accuracy_break_even_primary.pdf`.",
    ]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return REPORT


if __name__ == "__main__":
    print(generate())
