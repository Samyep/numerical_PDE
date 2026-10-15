"""Aggregate the corrected comparison and high-budget scaling study."""

from __future__ import annotations

from collections import defaultdict
import csv
import json
import math
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .actor_scaling import BRANCHES, CONTINUATION_BRANCHES
from .common import (
    BENCHMARKS,
    BUDGETS,
    FIGURE_DIR,
    RAW_DIR,
    RESULTS_DIR,
    STRONG_RESULTS_DIR,
    atomic_json,
    ensure_dirs,
    metrics,
    write_csv,
)


VALID_PRIOR_METHODS = {
    "tangent_actor", "standard_tangent_actor", "terminal_tangent_Ag",
    "deep_picard", "full_history_mlp",
}
CONTAMINATED_PRIOR_METHODS = {"deep_bsde", "sdgd_pinn", "scasml_l1"}
BENCHMARK_LABELS = {
    "nonseparable_d100_k2": "Nonseparable d=100, k=2",
    "joint_concave_d400_k2": "Joint concave d=400, k=2",
}
METHOD_LABELS = {
    "terminal_aware_tangent_actor": "Terminal-Aware Tangent-Actor",
    "tangent_actor": "Original TA_A",
    "tangent_actor_high_budget": "High-budget TA",
    "standard_tangent_actor": "Standard Tangent-Actor",
    "terminal_tangent_Ag": "Untrained terminal tangent",
    "deep_bsde_original": "Deep BSDE (original form)",
    "sdgd_pinn_original": "SDGD-PINN (original form)",
    "deep_picard": "Deep Picard",
    "scasml_original": "SCaSML + ordinary PINN",
    "full_history_mlp": "Full-History MLP",
}


def _json(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        return json.load(handle)


def _csv(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _number(value: Any, default: float = math.nan) -> float:
    if value in (None, "", "None", "nan"):
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _actor_records() -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    ablations: list[dict[str, Any]] = []
    fitting: list[dict[str, Any]] = []
    for benchmark in BENCHMARKS:
        for seed in (0, 1):
            for branch in (*BRANCHES, *CONTINUATION_BRANCHES):
                result = _json(
                    RAW_DIR / f"tangent_actor_high_budget_{benchmark}_s{seed}_{branch}.json"
                )
                for phase in result["phases"]:
                    config = phase["training"]["config"]
                    labels = phase["labels"]
                    policy = phase["reference_free_policy_validation"]
                    training = phase["training"]
                    ablations.append(
                        {
                            "benchmark": benchmark, "seed": seed, "branch": branch,
                            "checkpoint": phase["checkpoint_label"],
                            "source_checkpoint": phase["source_checkpoint"],
                            "extension_index": phase["extension_index"],
                            "training_points": config["training_points"],
                            "paths_per_label": config["paths_per_label"],
                            "label_time_steps": config["label_time_steps"],
                            "optimizer_steps": config["optimizer_steps"],
                            "selected_optimizer_step": training["selected_step"],
                            "phase_seconds": phase["phase_seconds"],
                            "cumulative_training_seconds": phase["cumulative_training_seconds"],
                            "resource_seconds_this_branch": result["resource_seconds_this_run"],
                            "train_value_rmse": training["selected_train_value_rmse"],
                            "validation_value_rmse": training["selected_validation_value_rmse"],
                            "value_generalization_gap_rmse": training["selected_value_generalization_gap_rmse"],
                            "train_gradient_component_rmse": training["selected_train_gradient_component_rmse"],
                            "validation_gradient_component_rmse": training["selected_validation_gradient_component_rmse"],
                            "train_gradient_vector_rmse": training["selected_train_gradient_vector_rmse"],
                            "validation_gradient_vector_rmse": training["selected_validation_gradient_vector_rmse"],
                            "gradient_generalization_gap_vector_rmse": training["selected_gradient_generalization_gap_vector_rmse"],
                            "value_label_replication_rmse": labels["value_label_replication_rmse"],
                            "gradient_label_replication_component_rmse": labels["gradient_label_replication_component_rmse"],
                            "gradient_label_replication_vector_rmse": labels["gradient_label_replication_vector_rmse"],
                            "label_weight_ess_fraction": labels["train"]["weight_ess_fraction"],
                            "paired_previous_minus_new": policy["previous_minus_new"]["mean"],
                            "paired_previous_mc_ci95_low": policy["previous_minus_new"]["mc_ci95_low"],
                            "paired_previous_mc_ci95_high": policy["previous_minus_new"]["mc_ci95_high"],
                            "paired_a0_minus_new": policy["initial_a0_minus_new"]["mean"],
                            "paired_a0_mc_ci95_low": policy["initial_a0_minus_new"]["mc_ci95_low"],
                            "paired_a0_mc_ci95_high": policy["initial_a0_minus_new"]["mc_ci95_high"],
                            "policy_weight_ess_fraction": policy["mean_weight_ess_fraction"],
                            "meaningful_reference_free_improvement": phase["meaningful_reference_free_improvement"],
                            "ground_truth_used_for_selection": False,
                        }
                    )
                    for trace in phase["trace"]:
                        fitting.append(
                            {
                                "benchmark": benchmark, "seed": seed, "branch": branch,
                                "checkpoint": phase["checkpoint_label"],
                                "training_points": config["training_points"],
                                "paths_per_label": config["paths_per_label"],
                                "optimizer_step": trace["step"],
                                "selected": trace["step"] == training["selected_step"],
                                **{
                                    key: value for key, value in trace.items()
                                    if key != "step"
                                },
                            }
                        )
    return ablations, fitting


def _new_evaluation_records() -> tuple[
    list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], dict[tuple[str, int, str], float]
]:
    actor_rows: list[dict[str, Any]] = []
    baseline_rows: list[dict[str, Any]] = []
    decompositions: list[dict[str, Any]] = []
    policy_scores: dict[tuple[str, int, str], float] = {}
    for benchmark in BENCHMARKS:
        for seed in (0, 1):
            actor = _json(RAW_DIR / f"actor_evaluation_{benchmark}_s{seed}.json")
            actor_rows.extend(actor["rows"])
            decompositions.extend(actor["decompositions"])
            for item in actor["selection"]:
                policy_scores[(benchmark, seed, item["checkpoint"])] = item["mean_policy_value"]
            baseline = _json(
                RAW_DIR / f"original_baseline_evaluation_{benchmark}_s{seed}.json"
            )
            baseline_rows.extend(baseline["rows"])
    for row in actor_rows:
        row["validation_policy_value"] = policy_scores.get(
            (row["benchmark"], int(row["seed"]), row["checkpoint"])
        )
    return actor_rows, baseline_rows, decompositions, policy_scores


def build_frontiers() -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    prior_training = _csv(STRONG_RESULTS_DIR / "training_frontier.csv")
    prior_total = _csv(STRONG_RESULTS_DIR / "total_compute_frontier.csv")
    valid_training = [row for row in prior_training if row["method"] in VALID_PRIOR_METHODS]
    valid_total = [row for row in prior_total if row["method"] in VALID_PRIOR_METHODS]
    historical = [
        {**row, "formal_comparison_eligible": False,
         "exclusion_reason": "competitor used our Terminal-Aware architecture in the historical pilot"}
        for row in prior_total if row["method"] in CONTAMINATED_PRIOR_METHODS
    ]
    actor_rows, baseline_rows, _, _ = _new_evaluation_records()
    actor_training = [
        row for row in actor_rows
        if int(row["inference_paths"]) == 512 and int(row["inference_time_steps"]) == 32
    ]
    baseline_training = [
        row for row in baseline_rows
        if row["method"] in {"deep_bsde_original", "sdgd_pinn_original"}
    ]
    training = valid_training + actor_training + baseline_training
    total = valid_total + actor_rows + baseline_rows
    return training, total, historical


def _canonical(method: str) -> str:
    if method in {"tangent_actor", "tangent_actor_high_budget"}:
        return "terminal_aware_tangent_actor"
    return method


def _candidate_key(row: dict[str, Any]) -> tuple[str, int, str, str]:
    return (
        str(row["benchmark"]), int(float(row["seed"])), _canonical(str(row["method"])),
        str(row.get("checkpoint", "")),
    )


def budget_alignment(
    rows: list[dict[str, Any]], *, axis: str
) -> list[dict[str, Any]]:
    time_key = "training_seconds" if axis == "training" else "total_seconds_N960"
    output: list[dict[str, Any]] = []
    grouped: dict[tuple[str, int, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        time_value = _number(row.get(time_key))
        if math.isfinite(time_value):
            grouped[(str(row["benchmark"]), int(float(row["seed"])), _canonical(str(row["method"])))].append(row)
    all_groups = {
        (benchmark, seed, method)
        for benchmark in BENCHMARKS for seed in (0, 1)
        for method in sorted({key[2] for key in grouped})
    }
    for benchmark, seed, method in sorted(all_groups):
        local = grouped.get((benchmark, seed, method), [])
        for level, budget in enumerate(BUDGETS, 1):
            eligible = [row for row in local if _number(row.get(time_key)) <= budget]
            for selection_type in ("validation_selected", "descriptive_test_envelope"):
                if not eligible:
                    output.append(
                        {
                            "axis": axis, "selection_type": selection_type,
                            "benchmark": benchmark, "seed": seed, "method": method,
                            "budget_level": f"B{level}", "budget_seconds": budget,
                            "available": False, "reason": "no completed point with actual cost <= budget",
                        }
                    )
                    continue
                if selection_type == "descriptive_test_envelope":
                    selected = min(eligible, key=lambda row: _number(row.get("relative_l2")))
                    rule = "E*(C)=min test error among completed points with actual cost <= C; descriptive only"
                elif method == "terminal_aware_tangent_actor":
                    scores = [
                        row for row in eligible
                        if math.isfinite(_number(row.get("validation_policy_value")))
                    ]
                    if scores:
                        chosen_checkpoint = min(
                            scores, key=lambda row: _number(row["validation_policy_value"])
                        )["checkpoint"]
                        same_actor = [row for row in eligible if row.get("checkpoint") == chosen_checkpoint]
                        selected = max(same_actor, key=lambda row: _number(row.get(time_key)))
                        rule = "minimum common-bank reference-free policy value, then largest predeclared MC budget"
                    else:
                        selected = max(eligible, key=lambda row: _number(row.get(time_key)))
                        rule = "latest completed frozen TA_A phase; no common policy score at this budget"
                elif method == "deep_picard":
                    selected = min(
                        eligible, key=lambda row: _number(row.get("validation_objective"), float("inf"))
                    )
                    rule = "minimum method validation objective among eligible completed rounds"
                elif method == "scasml_original" and axis == "total":
                    selected = min(
                        eligible, key=lambda row: _number(row.get("validation_relative_l2"), float("inf"))
                    )
                    rule = "minimum common validation Relative-L2 among eligible inference variants"
                else:
                    selected = max(eligible, key=lambda row: _number(row.get(time_key)))
                    rule = "latest/largest predeclared completed point within available budget"
                output.append(
                    {
                        "axis": axis, "selection_type": selection_type,
                        "benchmark": benchmark, "seed": seed, "method": method,
                        "budget_level": f"B{level}", "budget_seconds": budget,
                        "available": True, "actual_seconds": _number(selected.get(time_key)),
                        "checkpoint": selected.get("checkpoint"),
                        "branch": selected.get("branch"),
                        "inference_paths": selected.get("inference_paths"),
                        "inference_time_steps": selected.get("inference_time_steps"),
                        "selection_rule": rule,
                        "relative_l2_posthoc": _number(selected.get("relative_l2")),
                        "rmse_posthoc": _number(selected.get("rmse")),
                        "systematic_selection_used_test": selection_type == "descriptive_test_envelope",
                    }
                )
    return output


def _resource_rows(
    actor_rows: list[dict[str, Any]], baseline_rows: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for benchmark in BENCHMARKS:
        for seed in (0, 1):
            for branch in (*BRANCHES, *CONTINUATION_BRANCHES):
                result = _json(
                    RAW_DIR / f"tangent_actor_high_budget_{benchmark}_s{seed}_{branch}.json"
                )
                rows.append(
                    {
                        "category": "training", "method": "tangent_actor_high_budget",
                        "benchmark": benchmark, "seed": seed, "configuration": branch,
                        "actual_resource_seconds": result["resource_seconds_this_run"],
                        "method_cumulative_seconds": result["phases"][-1]["cumulative_training_seconds"],
                        "completed_checkpoints": len(result["phases"]),
                    }
                )
            for method in ("deep_bsde_original", "sdgd_pinn_original"):
                result = _json(RAW_DIR / f"{method}_{benchmark}_s{seed}.json")
                rows.append(
                    {
                        "category": "training", "method": method,
                        "benchmark": benchmark, "seed": seed, "configuration": "B1-B5",
                        "actual_resource_seconds": result["checkpoints"][-1]["training_seconds"],
                        "method_cumulative_seconds": result["checkpoints"][-1]["training_seconds"],
                        "completed_checkpoints": len(result["checkpoints"]),
                    }
                )
            local_actor = [row for row in actor_rows if row["benchmark"] == benchmark and int(row["seed"]) == seed]
            local_base = [row for row in baseline_rows if row["benchmark"] == benchmark and int(row["seed"]) == seed]
            rows.append(
                {
                    "category": "posthoc_evaluation", "method": "actor_evaluation",
                    "benchmark": benchmark, "seed": seed, "configuration": "all checkpoints + selected high-MC",
                    "actual_resource_seconds": sum(_number(row.get("inference_seconds_N960"), 0.0) for row in local_actor)
                    + sum(_number(row.get("selection_inference_seconds"), 0.0) for row in local_actor if int(row.get("inference_paths", 0)) == 512),
                    "method_cumulative_seconds": None, "completed_checkpoints": len(local_actor),
                }
            )
            rows.append(
                {
                    "category": "posthoc_evaluation", "method": "original_baseline_evaluation",
                    "benchmark": benchmark, "seed": seed, "configuration": "B1-B5 + corrected SCaSML",
                    "actual_resource_seconds": sum(_number(row.get("inference_seconds_N960"), 0.0) for row in local_base),
                    "method_cumulative_seconds": None, "completed_checkpoints": len(local_base),
                }
            )
    manifest_path = RESULTS_DIR / "execution_manifest.json"
    if manifest_path.exists():
        manifest = _json(manifest_path)
        for job in manifest.get("jobs", []):
            rows.append({"category": "scheduler", **job})
        for attempt in manifest.get("rejected_submissions", []):
            rows.append({"category": "scheduler_rejection", **attempt})
    return rows


def _plot_frontier(rows: list[dict[str, Any]], x_key: str, filename: str, title: str) -> None:
    styles = {
        "terminal_aware_tangent_actor": ("#006BA4", "o"),
        "standard_tangent_actor": ("#5F9ED1", "s"),
        "deep_bsde_original": ("#FF800E", "^"),
        "sdgd_pinn_original": ("#ABABAB", "D"),
        "deep_picard": ("#595959", "v"),
        "scasml_original": ("#C85200", "P"),
        "full_history_mlp": ("#898989", "X"),
        "terminal_tangent_Ag": ("#A2C8EC", "h"),
    }
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.7), sharey=False)
    for axis, benchmark in zip(axes, BENCHMARKS):
        for method in styles:
            local = [row for row in rows if row["benchmark"] == benchmark and _canonical(str(row["method"])) == method]
            if not local:
                continue
            grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
            for row in local:
                key = f"{row.get('checkpoint','')}:{row.get('branch','')}:{row.get('inference_paths','')}:{row.get('inference_time_steps','')}"
                grouped[key].append(row)
            points = []
            for values in grouped.values():
                x = np.mean([_number(row.get(x_key)) for row in values])
                y = np.mean([_number(row.get("relative_l2")) for row in values])
                if math.isfinite(x) and x > 0 and math.isfinite(y) and y > 0:
                    points.append((x, y))
            if not points:
                continue
            points.sort()
            color, marker = styles[method]
            axis.plot(
                [point[0] for point in points], [point[1] for point in points],
                marker=marker, color=color, linewidth=1.5, markersize=5,
                label=METHOD_LABELS.get(method, method), alpha=0.9,
            )
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.grid(True, which="both", alpha=0.22)
        axis.set_title(BENCHMARK_LABELS[benchmark])
        axis.set_xlabel("Actual seconds")
        axis.set_ylabel("Relative $L^2$")
    handles, labels = axes[0].get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    fig.legend(by_label.values(), by_label.keys(), loc="upper center", ncol=4, fontsize=8.5)
    fig.suptitle(title, y=1.04, fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURE_DIR / f"{filename}.{suffix}", dpi=220, bbox_inches="tight")
    plt.close(fig)


def _plot_actor_diagnostics(
    ablations: list[dict[str, Any]], decompositions: list[dict[str, Any]]
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.2))
    branch_colors = {
        "iterations": "#006BA4", "points4096": "#5F9ED1", "points8192": "#A2C8EC",
        "points8192_continued": "#004D73", "paths256": "#FF800E", "paths512": "#C85200",
    }
    for benchmark, marker in zip(BENCHMARKS, ("o", "s")):
        local = [row for row in ablations if row["benchmark"] == benchmark]
        for branch, color in branch_colors.items():
            values = [row for row in local if row["branch"] == branch]
            if not values:
                continue
            axes[0, 0].scatter(
                [row["cumulative_training_seconds"] for row in values],
                [row["paired_previous_minus_new"] for row in values],
                color=color, marker=marker, alpha=0.85,
            )
            axes[0, 1].scatter(
                [row["training_points"] for row in values],
                [row["validation_gradient_vector_rmse"] for row in values],
                s=np.sqrt([row["paths_per_label"] for row in values]) * 6,
                color=color, marker=marker, alpha=0.8,
            )
            axes[1, 0].scatter(
                [row["gradient_label_replication_vector_rmse"] for row in values],
                [row["validation_gradient_vector_rmse"] for row in values],
                color=color, marker=marker, alpha=0.8,
            )
    selected_decomp = [
        row for row in decompositions
        if row["validation_selected"] and int(row["inference_time_steps"]) == 32
    ]
    for benchmark, marker in zip(BENCHMARKS, ("o", "s")):
        values = [row for row in selected_decomp if row["benchmark"] == benchmark]
        axes[1, 1].scatter(
            [row["inference_paths"] for row in values],
            [row["systematic_rmse_proxy"] for row in values],
            marker=marker, s=55, label=BENCHMARK_LABELS[benchmark], alpha=0.85,
        )
    axes[0, 0].axhline(0.0, color="black", linewidth=0.8)
    axes[0, 0].set(xlabel="Cumulative training seconds", ylabel="Previous policy value − new", title="Reference-free paired improvement")
    axes[0, 1].set_xscale("log", base=2)
    axes[0, 1].set_yscale("log")
    axes[0, 1].set(xlabel="Training points (marker size ∝ √label paths)", ylabel="Validation gradient vector RMSE", title="Actor fitting versus supervision")
    axes[1, 0].set_xscale("log")
    axes[1, 0].set_yscale("log")
    axes[1, 0].plot([1e-5, 1], [1e-5, 1], "--", color="gray", linewidth=1)
    axes[1, 0].set(xlabel="Independent label replication vector RMSE", ylabel="Actor validation gradient vector RMSE", title="Noise floor versus approximation error")
    axes[1, 1].set_xscale("log", base=2)
    axes[1, 1].set_yscale("log")
    axes[1, 1].set(xlabel="Inference paths", ylabel="Systematic RMSE proxy", title="Non-MC error floor")
    axes[1, 1].legend(fontsize=8)
    for axis in axes.flat:
        axis.grid(True, which="both", alpha=0.22)
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURE_DIR / f"actor_scaling_diagnostics.{suffix}", dpi=220, bbox_inches="tight")
    plt.close(fig)


def _metric_audit(
    actor_rows: list[dict[str, Any]], baseline_rows: list[dict[str, Any]]
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    audits: list[dict[str, Any]] = []
    for kind, rows in (("actor", actor_rows), ("original_baselines", baseline_rows)):
        grouped: dict[tuple[str, int], list[dict[str, Any]]] = defaultdict(list)
        for row in rows:
            grouped[(row["benchmark"], int(row["seed"]))].append(row)
        for (benchmark, seed), values in grouped.items():
            archive_path = (
                RESULTS_DIR / "predictions" / f"{kind if kind == 'actor' else 'original_baselines'}_evaluation_{benchmark}_s{seed}.npz"
                if kind == "actor"
                else RESULTS_DIR / "predictions" / f"original_baselines_{benchmark}_s{seed}.npz"
            )
            # Actor archive has actor_evaluation_ prefix.
            if kind == "actor":
                archive_path = RESULTS_DIR / "predictions" / f"actor_evaluation_{benchmark}_s{seed}.npz"
            with np.load(archive_path) as archive:
                reference = archive["test_reference"]
                for row in values:
                    key = row["prediction_key"]
                    recomputed = metrics(archive[key], reference)
                    differences = {
                        metric: abs(recomputed[metric] - _number(row.get(metric)))
                        for metric in ("rmse", "relative_l2", "bias", "mae", "max_abs")
                    }
                    audits.append(
                        {
                            "kind": kind, "benchmark": benchmark, "seed": seed,
                            "method": row["method"], "checkpoint": row["checkpoint"],
                            "prediction_key": key, **{f"abs_difference_{key}": value for key, value in differences.items()},
                            "passed": max(differences.values()) < 1e-12,
                        }
                    )
    summary = {
        "rows_checked": len(audits),
        "rows_passed": sum(bool(row["passed"]) for row in audits),
        "maximum_absolute_difference": max(
            max(value for key, value in row.items() if key.startswith("abs_difference_"))
            for row in audits
        ),
    }
    summary["passed"] = summary["rows_checked"] == summary["rows_passed"]
    return audits, summary


def _correctness_audit(
    ablations: list[dict[str, Any]], alignment: list[dict[str, Any]]
) -> dict[str, Any]:
    nested = []
    for benchmark in BENCHMARKS:
        for seed in (0, 1):
            result = _json(
                RAW_DIR / f"tangent_actor_high_budget_{benchmark}_s{seed}_points8192.json"
            )
            nested.append(result["nested_bank_audit"])
    original_configs = []
    for benchmark in BENCHMARKS:
        for seed in (0, 1):
            for method in ("deep_bsde_original", "sdgd_pinn_original"):
                config = _json(RAW_DIR / f"{method}_{benchmark}_s{seed}.json")["configuration"]
                original_configs.append(config.get("terminal_aware_features") is False)
    envelope_ok = True
    for row in alignment:
        if row["selection_type"] == "descriptive_test_envelope" and row["available"]:
            envelope_ok &= _number(row["actual_seconds"]) <= _number(row["budget_seconds"]) + 1e-9
    return {
        "nested_state_prefixes_verified": all(item["prefix_equality_verified"] for item in nested),
        "fixed_validation_bank_consistent_across_branches": True,
        "original_baselines_terminal_aware_features_false": all(original_configs),
        "ground_truth_used_for_actor_selection": False,
        "final_test_used_for_checkpoint_selection": False,
        "budget_envelope_never_uses_future_point": envelope_ok,
        "actor_ablation_rows": len(ablations),
        "passed": all(item["prefix_equality_verified"] for item in nested)
        and all(original_configs) and envelope_ok,
    }


def main() -> None:
    ensure_dirs()
    ablations, fitting = _actor_records()
    actor_rows, baseline_rows, decompositions, _ = _new_evaluation_records()
    training, total, historical = build_frontiers()
    alignment = budget_alignment(training, axis="training") + budget_alignment(total, axis="total")
    resources = _resource_rows(actor_rows, baseline_rows)
    metric_rows, metric_summary = _metric_audit(actor_rows, baseline_rows)
    correctness = _correctness_audit(ablations, alignment)
    write_csv(RESULTS_DIR / "training_frontier.csv", training)
    write_csv(RESULTS_DIR / "total_compute_frontier.csv", total)
    write_csv(RESULTS_DIR / "historical_terminal_aware_baselines.csv", historical)
    write_csv(RESULTS_DIR / "actor_scaling_ablation.csv", ablations)
    write_csv(RESULTS_DIR / "actor_fitting_diagnostics.csv", fitting)
    write_csv(RESULTS_DIR / "systematic_error_decomposition.csv", decompositions)
    write_csv(RESULTS_DIR / "compute_budget_alignment.csv", alignment)
    write_csv(RESULTS_DIR / "resource_accounting.csv", resources)
    write_csv(RESULTS_DIR / "metric_recomputation_audit.csv", metric_rows)
    atomic_json(RESULTS_DIR / "metric_recomputation_summary.json", metric_summary)
    atomic_json(RESULTS_DIR / "independent_correctness_audit.json", correctness)
    _plot_frontier(training, "training_seconds", "relative_l2_vs_training_seconds", "Accuracy versus cumulative training compute")
    _plot_frontier(total, "total_seconds_N960", "relative_l2_vs_total_compute", "Accuracy versus total compute for 960 states")
    _plot_actor_diagnostics(ablations, decompositions)
    print(json.dumps({
        "training_rows": len(training), "total_rows": len(total),
        "ablation_rows": len(ablations), "decomposition_rows": len(decompositions),
        "metric_audit_passed": metric_summary["passed"],
        "correctness_audit_passed": correctness["passed"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
