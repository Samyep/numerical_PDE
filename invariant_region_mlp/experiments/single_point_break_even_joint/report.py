"""Build the compact joint-PDE matched-accuracy report."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from .protocol import RESULTS_DIR, SOURCE_ROOT, TARGETS


REPORT_PATH = SOURCE_ROOT / "invariant_region_mlp" / "docs" / "SINGLE_POINT_BREAK_EVEN_JOINT_REPORT.md"


def _read(path: Path) -> dict[str, Any]: return json.loads(path.read_text(encoding="utf-8"))
def _fmt(value: Any) -> str: return "--" if value is None else (f"{value:.5g}" if isinstance(value, float) else str(value))


def build() -> Path:
    analysis = _read(RESULTS_DIR / "analysis_summary.json"); frozen = _read(RESULTS_DIR / "FROZEN_SELECTIONS.json")
    reference = _read(RESULTS_DIR / "reference_audit.json"); preflight = _read(RESULTS_DIR / "preflight_audit.json")
    artifact = _read(RESULTS_DIR / "artifact_audit.json") if (RESULTS_DIR / "artifact_audit.json").is_file() else {"status": "PENDING", "allocated_h100_gpu_hours": None}
    coverage = _read(RESULTS_DIR / "timing_coverage.json")
    rows = []
    for target, threshold in TARGETS.items():
        for method, label in (("tangent_actor", "Tangent-Actor"), ("fh_mlp_standard", "Full-History MLP"), ("deep_bsde", "Deep BSDE")):
            item = analysis["targets"][target][method]
            if method == "tangent_actor": crossing = "reference curve" if coverage["targets"][target][method] == "RUN" else coverage["targets"][target][method]
            else: crossing = analysis["break_even"][target].get(method, {}).get("status", "NOT_EVALUATED")
            rows.append(f"| {target} ({threshold:.1e}) | {label} | {item.get('selected_candidate') or '--'} | {item['status']} | {_fmt(item.get('final_rmse'))} | {crossing} |")
    statements = []
    for target in TARGETS:
        for method, label in (("fh_mlp_standard", "Full-History MLP"), ("deep_bsde", "Deep BSDE")):
            item = analysis["break_even"][target].get(method, {}); status = item.get("status", "NOT_EVALUATED")
            if status == "OBSERVED_BREAK_EVEN": statements.append(f"- **{target}, {label}:** observed in `({item.get('previous_queries')}, {item['crossing_queries']}]` queries; at the first measured crossing TA used `{item['tangent_actor_mean_seconds_at_crossing']:.3f}` s and the baseline `{item['baseline_mean_seconds_at_crossing']:.3f}` s.")
            elif status == "NOT_REACHED_WITHIN_TESTED_BUDGET": statements.append(f"- **{target}, {label}:** `NOT_REACHED_WITHIN_TESTED_BUDGET`; at `Q={item['largest_common_queries']}`, TA used `{item['tangent_actor_mean_seconds']:.3f}` s and the baseline `{item['baseline_mean_seconds']:.3f}` s (baseline/TA `{item['baseline_over_tangent_actor_time']:.4f}`).")
            else: statements.append(f"- **{target}, {label}:** `{status}`.")
    choices = []
    for target in TARGETS:
        values = []
        for method in ("tangent_actor", "fh_mlp_standard", "deep_bsde"):
            item = frozen["targets"][target][method]; values.append(f"`{method}` = `{item.get('selected_candidate')}` (calibration RMSE {_fmt(item.get('calibration_rmse', item.get('best_calibration_rmse')))})")
        choices.append(f"- **{target}:** " + "; ".join(values) + ".")
    report = f"""# Joint Nonseparable Matched-Accuracy Amortization Break-Even

## Outcome

| Target | Method | Frozen configuration | Final status | Final pooled RMSE | Measured break-even |
|---|---|---|---|---:|---|
{chr(10).join(rows)}

The crossing is the first actually executed common query prefix at which the
mean baseline total is no smaller than the mean Tangent-Actor total.  No
single-query extrapolation, curve fit, interpolation, ensemble prediction, or
post-final configuration change is used.

## Measured amortization result

{chr(10).join(statements)}

## Frozen calibration choices

{chr(10).join(choices)}

Calibration used 32 new queries and final accuracy used a disjoint 64-query
bank, except for the mandated shared origin.  Tangent-Actor is trained once;
FH-MLP re-executes its recursion per point; native Deep BSDE trains a fresh
scalar-$Y_0$ model per point.  Secondary timing is deliberately omitted when
it would require a different selected configuration, preserving the frozen
18 H100-hour ceiling.

## Correctness and accounting

- Independent reference: **{reference['status']}**; fresh-bank primary/audit RMSE `{reference['fresh_bank_value_refinement_rmse']:.6e}`.
- H100/source preflight: **{preflight['status']}** on `{preflight['gpu']}`.
- Artifact audit: **{artifact['status']}**.
- Allocated H100 task-hours, including recorded failures: `{_fmt(artifact.get('allocated_h100_gpu_hours'))}` / `18`.

Method-active time includes Actor label construction, fitting, reference-free
selection and inference for TA; recursive inference for FH-MLP; and complete
fresh per-query training for Deep BSDE.  The 1,048,576-state timing bank has
no truth array.  Queueing, reference construction, truth evaluation, hashing,
plotting and upload are excluded from method curves.

## Provenance

- Frozen addendum commit: `8a18360be907da36b8b875c4fdf6ba0b0cc1007d`.
- Imported validated-source commit: `b7dd884225936a0d67f770efce1207bc478d4f7c`.
- Frozen implementation commit: `{frozen['code_commit']}`.
- Frozen selection hash: `{frozen['canonical_sha256']}`.
- Analysis commit: `{analysis['analysis_code_commit']}`.
- Figure: `results/single_point_break_even_joint/figures/joint_matched_accuracy_break_even.pdf`.

Large checkpoints and the truth-free timing bank remain in the ICE archive;
their hashes and every reported prediction hash are retained in compact
records.  This is the final preregistered extension; no third PDE is started.
"""
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True); temporary = REPORT_PATH.with_suffix(REPORT_PATH.suffix + ".tmp"); temporary.write_text(report, encoding="utf-8"); temporary.replace(REPORT_PATH); return REPORT_PATH


if __name__ == "__main__": print(build())
