"""Render the audited experiment report from final CSV/JSON artifacts."""

from __future__ import annotations

import csv
import json

from .protocol import PDES, REPORT_PATH, RESULTS_DIR


LABEL = {
    "b6_fullrank": "B6 Full-Rank HJB",
    "dense_fourier": "Dense Fourier",
    "tangent_actor": "Tangent-Actor",
    "dpi": "DPI",
    "pinn_hte": "PINN-HTE",
}


def _csv(name: str) -> list[dict[str, str]]:
    with (RESULTS_DIR / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _bool(value: str) -> bool:
    return value.lower() == "true"


def _original(rows: list[dict[str, str]], pde: str, method: str) -> tuple[dict, float, str]:
    boundary = [
        row for row in rows
        if row["pde"] == pde and row["method"] == method
        and row["configuration_id"] == "C0" and _bool(row["is_original_budget"])
    ]
    if len(boundary) != 1:
        raise RuntimeError((pde, method, len(boundary)))
    item = boundary[0]
    selected_id = item["checkpoint_id"]
    rmse = float(item["rmse"])
    if method == "tangent_actor" and item.get("causal_selected_checkpoint"):
        selected_id = item["causal_selected_checkpoint"]
        selected = [
            row for row in rows
            if row["pde"] == pde and row["method"] == method
            and row["configuration_id"] == "C0"
            and row["checkpoint_id"] == selected_id
        ]
        if len(selected) != 1:
            raise RuntimeError((pde, method, selected_id))
        rmse = float(selected[0]["rmse"])
    return item, rmse, selected_id


def _overlap_separation(
    rows: list[dict[str, str]], pde: str, lower: float, upper: float,
) -> tuple[float, float, float]:
    def in_interval(row: dict[str, str]) -> bool:
        minutes = float(row["main_axis_seconds"]) / 60.0
        return (
            row["pde"] == pde
            and row["configuration_id"] == "C0"
            and lower <= minutes <= upper
        )

    ta = [float(row["rmse"]) for row in rows if in_interval(row) and row["method"] == "tangent_actor"]
    baselines = [float(row["rmse"]) for row in rows if in_interval(row) and row["method"] != "tangent_actor"]
    if not ta or not baselines:
        raise RuntimeError((pde, lower, upper, len(ta), len(baselines)))
    worst_ta = max(ta)
    best_baseline = min(baselines)
    return worst_ta, best_baseline, best_baseline / worst_ta


def render() -> None:
    rows = _csv("checkpoint_metrics.csv")
    trajectories = _csv("trajectory_summary.csv")
    accounting = _csv("compute_accounting.csv")
    overlap = json.loads((RESULTS_DIR / "overlap_audit.json").read_text(encoding="utf-8"))
    audit = json.loads((RESULTS_DIR / "artifact_audit.json").read_text(encoding="utf-8"))
    protocol = json.loads((RESULTS_DIR / "EXPERIMENT_PROTOCOL.json").read_text(encoding="utf-8"))
    smoke = json.loads((RESULTS_DIR / "raw" / "smoke_nonresults.json").read_text(encoding="utf-8"))
    diagnostic_path = RESULTS_DIR / "diagnostics" / "dense_dpi_callback_replay.json"
    callback_diagnostic = (
        json.loads(diagnostic_path.read_text(encoding="utf-8"))
        if diagnostic_path.is_file() else None
    )
    totals = [row for row in accounting if row["job_id"] == "TOTAL"]
    total = totals[0] if totals else {}
    diagnostic_jobs = [row for row in accounting if row.get("job_name") == "wcf_dpi_audit"]
    diagnostic_gpu_seconds = sum(
        float(row.get("allocated_gpu_seconds") or 0.0) for row in diagnostic_jobs
    )
    lines = [
        "# Two-PDE 100D Accuracy--Compute Frontier", "",
        "## Outcome", "",
        "All entries below are measured checkpoints. No RMSE values were interpolated, smoothed, or repeated to create overlap.", "",
        "| PDE | Method | Prespecified C0 checkpoint | Training + selection (min) | RMSE | RelMAE |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for pde in PDES:
        for method in ("tangent_actor", "dpi", "pinn_hte"):
            boundary, rmse, selected_id = _original(rows, pde, method)
            selected_row = next(
                row for row in rows
                if row["pde"] == pde and row["method"] == method
                and row["configuration_id"] == "C0" and row["checkpoint_id"] == selected_id
            )
            lines.append(
                f"| {LABEL[pde]} | {LABEL[method]} | {selected_id} | "
                f"{float(boundary['main_axis_seconds']) / 60:.3f} | {rmse:.6g} | "
                f"{float(selected_row['relmae']):.6g} |"
            )
    lines += ["", "## Headline measured comparisons", ""]
    for pde in PDES:
        _ta_boundary, ta_rmse, _ta_id = _original(rows, pde, "tangent_actor")
        _dpi_boundary, dpi_rmse, _dpi_id = _original(rows, pde, "dpi")
        _pinn_boundary, pinn_rmse, _pinn_id = _original(rows, pde, "pinn_hte")
        interval = overlap[pde]["intersection_minutes"]
        worst_ta, best_baseline, conservative_gap = _overlap_separation(
            rows, pde, float(interval[0]), float(interval[1])
        )
        lines += [
            f"- **{LABEL[pde]} prespecified C0 checkpoints:** Tangent-Actor RMSE "
            f"`{ta_rmse:.6g}` is **{dpi_rmse / ta_rmse:.1f}x** lower than DPI and "
            f"**{pinn_rmse / ta_rmse:.1f}x** lower than PINN-HTE. These checkpoints "
            "have their own measured method-specific costs and are not an equal-time comparison.",
            f"- **{LABEL[pde]} observed C0 overlap ({float(interval[0]):.2f}--{float(interval[1]):.2f} min):** "
            f"the worst measured Tangent-Actor RMSE is `{worst_ta:.6g}`, while the best measured "
            f"DPI/PINN-HTE RMSE is `{best_baseline:.6g}`--a conservative **{conservative_gap:.1f}x** "
            "separation using measured checkpoints only (no time interpolation).",
        ]
    lines += ["", "## C0 anytime trajectories", "",
              "The primary curves are C0. C1/C2 are sensitivity analyses and are not collapsed into a test-selected envelope.", "",
              "| PDE | Method | First measured RMSE (min) | Final measured RMSE (min) | Stop reason |",
              "|---|---|---:|---:|---|" ]
    for pde in PDES:
        for method in ("tangent_actor", "dpi", "pinn_hte"):
            candidates = sorted(
                [row for row in rows if row["pde"] == pde and row["method"] == method and row["configuration_id"] == "C0"],
                key=lambda row: float(row["main_axis_seconds"]),
            )
            trajectory = next(
                row for row in trajectories
                if row["pde"] == pde and row["method"] == method and row["configuration_id"] == "C0"
            )
            first, final = candidates[0], candidates[-1]
            lines.append(
                f"| {LABEL[pde]} | {LABEL[method]} | "
                f"{float(first['rmse']):.5g} ({float(first['main_axis_seconds'])/60:.2f}) | "
                f"{float(final['rmse']):.5g} ({float(final['main_axis_seconds'])/60:.2f}) | "
                f"{trajectory['stop_reason']} |"
            )
    lines += [
        "", "### Why some curves stop before 30 minutes", "",
        "Thirty minutes was a maximum active-training window, not a requirement to manufacture a point at 30 minutes. "
        "The frozen protocol also imposed algorithmic iteration caps and truth-free validation stopping. "
        "Consequently, an endpoint before 30 minutes is expected when one of those rules fires first.", "",
        "| PDE | Method | C0 endpoint (min) | Mechanical reason |",
        "|---|---|---:|---|",
    ]
    reason_text = {
        "maximum_predeclared_actor_updates": "Reached the frozen maximum of 8 actor updates",
        "maximum_predeclared_picard_iterations": "Reached the frozen maximum of 320 Picard rounds",
        "predeclared_validation_plateau": "Triggered the frozen truth-free validation plateau rule",
        "active_training_cap": "Crossed the 30-minute cap at an algorithmic boundary",
        "active_training_cap_at_actor_boundary": "Crossed the 30-minute cap at a complete actor-update boundary",
        "configured_actor_updates": "Completed the configured actor updates",
    }
    for pde in PDES:
        for method in ("tangent_actor", "dpi", "pinn_hte"):
            trajectory = next(
                row for row in trajectories
                if row["pde"] == pde and row["method"] == method
                and row["configuration_id"] == "C0"
            )
            endpoint = (
                float(trajectory["training_seconds"])
                + float(trajectory["selection_seconds"])
            ) / 60.0
            reason = reason_text.get(trajectory["stop_reason"], trajectory["stop_reason"])
            lines.append(
                f"| {LABEL[pde]} | {LABEL[method]} | {endpoint:.2f} | {reason} |"
            )
    lines += [
        "",
        "The main paper figure is restricted to the preregistered 0--30 minute window and contains no "
        "original-budget star markers. Checkpoints beyond 30 minutes are not deleted: the Dense Fourier "
        "Tangent-Actor boundary points at 31.46 and 34.54 minutes remain in `checkpoint_metrics.csv` and "
        "the configuration/total-cost appendix figures. Lines are not extended by carrying their last value forward.",
    ]
    lines += ["", "## Configuration sensitivity", "",
              "Each row is the last genuinely trained checkpoint of that trajectory; this table is descriptive, not a selection rule.", "",
              "| PDE | Method | Config | Final time (min) | Final RMSE |",
              "|---|---|---|---:|---:|" ]
    for pde in PDES:
        for method in ("tangent_actor", "dpi", "pinn_hte"):
            for config in ("C0", "C1", "C2"):
                candidates = sorted(
                    [row for row in rows if row["pde"] == pde and row["method"] == method and row["configuration_id"] == config],
                    key=lambda row: float(row["main_axis_seconds"]),
                )
                final = candidates[-1]
                lines.append(
                    f"| {LABEL[pde]} | {LABEL[method]} | {config} | "
                    f"{float(final['main_axis_seconds'])/60:.3f} | {float(final['rmse']):.6g} |"
                )
    lines += ["", "## Observed cost overlap", ""]
    for pde in PDES:
        item = overlap[pde]
        intersection = item["intersection_minutes"]
        wording = (
            f"{intersection[0]:.2f}--{intersection[1]:.2f} min"
            if intersection is not None else "none"
        )
        lines.append(f"- **{LABEL[pde]}:** C0 observed-range intersection = {wording}.")
    lines += [
        "", "## Fidelity and accounting", "",
        f"- Non-result smoke preflight: **{smoke['status']}**.",
        f"- Artifact/fidelity audit: **{'PASS' if audit['overall_pass'] else 'FAIL'}**.",
        f"- Allocated GPU time including post-freeze diagnostics: **{float(total.get('allocated_gpu_hours') or 0):.3f} H100 GPU-hours**.",
        f"- Post-freeze diagnostic allocation: **{diagnostic_gpu_seconds / 3600:.3f} H100 GPU-hours** across {len(diagnostic_jobs)} attempts.",
        f"- Active measured training time: **{float(total.get('active_training_hours') or 0):.3f} GPU-hours**.",
        f"- 18-GPU-hour cap respected: **{total.get('within_18_gpu_hour_cap')}**.",
        f"- Protocol SHA-256: `{protocol['protocol_sha256']}`.",
        "- Complete checkpoint, prediction, test-bank, source-revision, and scheduler details are in the companion CSV/JSON artifacts.",
    ]
    failed_checks = [item for item in audit["checks"] if not item["pass"]]
    if failed_checks:
        lines += ["", "### Failed audit checks", ""]
        for item in failed_checks:
            lines.append(f"- **Check {item['id']} -- {item['name']}: FAIL.**")
            if item["id"] == 6:
                for detail in item.get("details", []):
                    if detail.get("pass"):
                        continue
                    difference = abs(float(detail["new_rmse"]) - float(detail["historical_rmse"]))
                    relative_difference = float(
                        detail.get(
                            "relative_difference_to_historical",
                            difference / abs(float(detail["historical_rmse"])),
                        )
                    )
                    lines.append(
                        f"  - {LABEL[detail['pde']]} / {LABEL[detail['method']]}: new original-budget "
                        f"RMSE `{float(detail['new_rmse']):.9g}` versus historical "
                        f"`{float(detail['historical_rmse']):.9g}`; absolute difference "
                        f"`{difference:.9g}` exceeds the frozen tolerance "
                        f"`{float(detail['tolerance']):.9g}`. Relative to the historical RMSE, the discrepancy is "
                        f"`{100 * relative_difference:.2f}%`. The failure is retained "
                        "without changing the threshold; with one training seed, inference uncertainty cannot adjudicate it."
                    )
    lines += ["", "### Dense Fourier DPI targeted consistency diagnosis", ""]
    lines += [
        "The expanded C0 YAML was compared byte-for-byte at the semantic level. Apart from the output directory, "
        "the only numerical difference is `EVAL.FREQ: 8` historically versus `EVAL.FREQ: null` in this study. "
        "The historical upstream callback samples 10,000 nonterminal states with the global Torch RNG and evaluates "
        "the exact solution during training. It therefore both adds substantial diagnostic work to the historical elapsed time and advances "
        "the same RNG stream later used for stochastic Picard training. The new run disabled that callback to enforce "
        "the frozen truth gate. Thus the 80-round models are not expected to be bitwise or numerically identical, and "
        "the old 10.56-minute and new 5.07-minute figures use different work definitions.",
    ]
    if callback_diagnostic is None:
        lines.append(
            "A post-freeze RNG-only callback replay has not been completed; the original consistency check remains FAIL."
        )
    else:
        failed_detail = next(
            detail for item in audit["checks"] if item["id"] == 6
            for detail in item.get("details", [])
            if detail["pde"] == "dense_fourier" and detail["method"] == "dpi"
        )
        frontier_difference = abs(
            float(failed_detail["new_rmse"]) - float(failed_detail["historical_rmse"])
        )
        replay_difference = abs(
            float(callback_diagnostic["replay_rmse"])
            - float(callback_diagnostic["historical_rmse"])
        )
        discrepancy_reduction = 1.0 - replay_difference / frontier_difference
        lines += [
            "The post-freeze replay reproduces only the historical callback's global-RNG draws, without evaluating "
            "nonterminal truth during training. It is diagnostic only and is excluded from every frontier curve.", "",
            "| Quantity | Value |", "|---|---:|",
            f"| Historical 80-round RMSE | {float(callback_diagnostic['historical_rmse']):.9g} |",
            f"| Truth-gated frontier 80-round RMSE | {float(failed_detail['new_rmse']):.9g} |",
            f"| RNG-only replay RMSE | {float(callback_diagnostic['replay_rmse']):.9g} |",
            f"| Relative RMSE difference | {100 * float(callback_diagnostic['relative_rmse_difference']):.3f}% |",
            f"| Reduction of the historical-RMSE discrepancy | {100 * discrepancy_reduction:.2f}% |",
            f"| Historical training time | {float(callback_diagnostic['historical_training_seconds']) / 60:.3f} min |",
            f"| RNG-only replay time (diagnostic instrumentation; not comparable) | {float(callback_diagnostic['replay_training_seconds']) / 60:.3f} min |",
            f"| RNG replay callback calls | {int(callback_diagnostic['rng_only_callback_calls'])} |",
            f"| Prediction-array RMSE (replay vs historical) | {float(callback_diagnostic['prediction_rmse_vs_historical']):.9g} |",
            "",
            "The endpoint RMSE moves from a 78.12% discrepancy to 2.56%, reducing the discrepancy by "
            f"{100 * discrepancy_reduction:.2f}%. This is strong evidence that callback-induced RNG advancement "
            "explains most of the failed consistency check. It is not a bitwise reproduction: the prediction arrays "
            "still differ, and there is only one training seed. Check 6 therefore remains FAIL rather than being waived.",
        ]
    lines += [
        "", "## Limitations", "",
        "This is a single-training-seed study (seed 0). Tangent-Actor uses two fixed Monte Carlo inference replicates, but those replicates do not measure training-seed variability. C1/C2 are sensitivity runs, not test-selected hyperparameters. Timing is H100-specific and separates active numerical work from Slurm allocation time. Existing SDGD-PINN fixed-budget results were not rerun and do not appear on these scaling curves.",
        "", "## Figures", "",
        "- `figures/b6_100d_accuracy_compute.pdf`",
        "- `figures/dense_fourier_100d_accuracy_compute.pdf`",
        "- `figures/two_pde_accuracy_compute.pdf` and `.png`",
        "- `figures/configuration_ablation.pdf`",
        "- `figures/accuracy_vs_total_cost.pdf`",
        "- `figures/checkpoint_time_coverage.pdf`",
        "- `figures/cost_breakdown.pdf`",
        "",
    ]
    REPORT_PATH.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    render()
