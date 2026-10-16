"""Build compact CSVs, figures, and reproducibility metadata."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import platform
import subprocess
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .common import (
    FIGURE_DIR,
    HIGH_BUDGET_RESULTS_DIR,
    RAW_DIR,
    RESULTS_DIR,
    STRONG_RESULTS_DIR,
    atomic_json,
    ensure_dirs,
    write_csv,
)


BENCHMARKS = ("nonseparable_d100_k2", "joint_concave_d400_k2")
LABELS = {
    "nonseparable_d100_k2": "Nonseparable (d=100)",
    "joint_concave_d400_k2": "Joint concave (d=400)",
}
COLORS = {"plain": "#6B6B6B", "mcv": "#006BA4"}


def _json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _csv(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _float(value: Any, default: float = math.nan) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _plain_training_runs() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for benchmark in BENCHMARKS:
        for seed in (0, 1):
            for path in sorted((HIGH_BUDGET_RESULTS_DIR / "raw").glob(
                f"tangent_actor_high_budget_{benchmark}_s{seed}_*.json"
            )):
                result = _json(path)
                branch = result["branch"]
                config = result["configuration"]
                for phase in result["phases"]:
                    labels = phase["labels"]
                    training = phase["training"]
                    value_noise = labels["value_label_replication_rmse"] / np.sqrt(2.0)
                    gradient_noise = (
                        labels["gradient_label_replication_vector_rmse"] / np.sqrt(2.0)
                    )
                    rows.append({
                        "diagnostic_kind": "selected_training_trace",
                        "training_variant": "plain",
                        "benchmark": benchmark,
                        "seed": seed,
                        "branch": branch,
                        "iteration": phase["extension_index"] + 3,
                        "training_points": config["training_points"],
                        "paths_per_label": config["paths_per_label"],
                        "cumulative_training_seconds": phase["cumulative_training_seconds"],
                        "selected_optimizer_step": training["selected_step"],
                        "actor_validation_value_rmse": training[
                            "selected_validation_value_rmse"
                        ],
                        "actor_validation_gradient_component_rmse": training[
                            "selected_validation_gradient_component_rmse"
                        ],
                        "actor_validation_gradient_vector_rmse": training[
                            "selected_validation_gradient_vector_rmse"
                        ],
                        "value_label_noise_rmse": value_noise,
                        "gradient_label_vector_noise_rmse": gradient_noise,
                        "value_fit_to_label_noise_ratio": training[
                            "selected_validation_value_rmse"
                        ] / max(value_noise, np.finfo(float).tiny),
                        "gradient_fit_to_label_noise_ratio": training[
                            "selected_validation_gradient_vector_rmse"
                        ] / max(gradient_noise, np.finfo(float).tiny),
                        "validation_target_paths": config["paths_per_label"],
                        "common_high_quality_target": False,
                    })
    return rows


def _mcv_training_runs() -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    diagnostics: list[dict[str, Any]] = []
    runtime: list[dict[str, Any]] = []
    for path in sorted(RAW_DIR.glob("mcv_points*_*.json")):
        result = _json(path)
        config = result["configuration"]
        for phase in result["phases"]:
            labels = phase["labels"]
            training = phase["training"]
            value_noise = labels["value_label_replication_rmse"] / np.sqrt(2.0)
            gradient_noise = labels["gradient_label_replication_vector_rmse"] / np.sqrt(2.0)
            diagnostics.append({
                "diagnostic_kind": "selected_training_trace",
                "training_variant": "mcv",
                "benchmark": result["benchmark"],
                "seed": result["seed"],
                "branch": result["branch"],
                "iteration": phase["iteration"],
                "training_points": config["training_points"],
                "paths_per_label": config["paths_per_label"],
                "cumulative_training_seconds": phase["cumulative_training_seconds"],
                "selected_optimizer_step": training["selected_step"],
                "actor_validation_value_rmse": training[
                    "selected_validation_value_rmse"
                ],
                "actor_validation_gradient_component_rmse": training[
                    "selected_validation_gradient_component_rmse"
                ],
                "actor_validation_gradient_vector_rmse": training[
                    "selected_validation_gradient_vector_rmse"
                ],
                "value_label_noise_rmse": value_noise,
                "gradient_label_vector_noise_rmse": gradient_noise,
                "value_fit_to_label_noise_ratio": training[
                    "selected_validation_value_rmse"
                ] / max(value_noise, np.finfo(float).tiny),
                "gradient_fit_to_label_noise_ratio": training[
                    "selected_validation_gradient_vector_rmse"
                ] / max(gradient_noise, np.finfo(float).tiny),
                "validation_target_paths": config["paths_per_label"],
                "common_high_quality_target": False,
            })
            label_peak = max(
                item.get("peak_gpu_bytes") or 0
                for item in (labels["train"], labels["validation"], labels["replication"])
            )
            runtime.append({
                "record_kind": "complete_actor_update",
                "training_variant": "mcv",
                "benchmark": result["benchmark"],
                "seed": result["seed"],
                "branch": result["branch"],
                "iteration": phase["iteration"],
                "training_points": config["training_points"],
                "paths_per_label": config["paths_per_label"],
                "base_a3_training_seconds": phase["base_a3_training_seconds"],
                "label_seconds": labels["wall_seconds"],
                "fit_seconds": training["wall_seconds"],
                "auxiliary_evaluation_checkpoint_seconds": (
                    phase["phase_seconds"] - labels["wall_seconds"] - training["wall_seconds"]
                ),
                "phase_seconds": phase["phase_seconds"],
                "additional_continuation_seconds": phase[
                    "additional_continuation_seconds"
                ],
                "cumulative_training_seconds": phase["cumulative_training_seconds"],
                "label_peak_gpu_bytes": label_peak,
                "fit_peak_gpu_bytes": training.get("peak_gpu_bytes"),
                "weight_ess_fraction_train": labels["train"]["weight_ess_fraction"],
                "nonfinite_count": sum(
                    item["nonfinite_count"]
                    for item in (labels["train"], labels["validation"], labels["replication"])
                ),
            })
    return diagnostics, runtime


def _plain_runtime() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for benchmark in BENCHMARKS:
        for seed in (0, 1):
            for path in sorted((HIGH_BUDGET_RESULTS_DIR / "raw").glob(
                f"tangent_actor_high_budget_{benchmark}_s{seed}_*.json"
            )):
                result = _json(path)
                config = result["configuration"]
                for phase in result["phases"]:
                    labels, training = phase["labels"], phase["training"]
                    rows.append({
                        "record_kind": "complete_actor_update",
                        "training_variant": "plain",
                        "benchmark": benchmark,
                        "seed": seed,
                        "branch": result["branch"],
                        "iteration": phase["extension_index"] + 3,
                        "training_points": config["training_points"],
                        "paths_per_label": config["paths_per_label"],
                        "base_a3_training_seconds": phase["base_a3_training_seconds"],
                        "label_seconds": labels["wall_seconds"],
                        "fit_seconds": training["wall_seconds"],
                        "auxiliary_evaluation_checkpoint_seconds": (
                            phase["phase_seconds"] - labels["wall_seconds"]
                            - training["wall_seconds"]
                        ),
                        "phase_seconds": phase["phase_seconds"],
                        "additional_continuation_seconds": (
                            phase["cumulative_training_seconds"]
                            - phase["base_a3_training_seconds"]
                        ),
                        "cumulative_training_seconds": phase[
                            "cumulative_training_seconds"
                        ],
                        "label_peak_gpu_bytes": max(
                            item.get("peak_gpu_bytes") or 0
                            for item in (
                                labels["train"], labels["validation"], labels["replication"]
                            )
                        ),
                        "fit_peak_gpu_bytes": training.get("peak_gpu_bytes"),
                        "weight_ess_fraction_train": labels["train"][
                            "weight_ess_fraction"
                        ],
                        "nonfinite_count": sum(
                            item["nonfinite_count"]
                            for item in (
                                labels["train"], labels["validation"], labels["replication"]
                            )
                        ),
                    })
    return rows


def _prefix_frontier() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for benchmark in BENCHMARKS:
        for seed in (0, 1):
            source = _json(
                STRONG_RESULTS_DIR / "raw" / f"evaluation_{benchmark}_s{seed}.json"
            )
            for row in source["training_frontier"]:
                if row["method"] != "tangent_actor" or row["checkpoint"] not in ("a1", "a2"):
                    continue
                rows.append({
                    "training_variant": "shared_prefix",
                    "benchmark": benchmark,
                    "seed": seed,
                    "evaluation_id": row["checkpoint"],
                    "checkpoint": row["checkpoint"],
                    "branch": "shared_prefix",
                    "training_points": 2048,
                    "paths_per_label": 128,
                    "training_seconds": row["training_seconds"],
                    "inference_variant": "plain",
                    "inference_paths": 512,
                    "inference_time_steps": 32,
                    "rmse": row["rmse"],
                    "relative_l2": row["relative_l2"],
                    "bias": row["bias"],
                    "systematic_rmse_proxy": "",
                    "final_test_used_for_selection": False,
                })
    return rows


def _envelope(rows: list[dict[str, Any]], key: str) -> tuple[np.ndarray, np.ndarray]:
    ordered = sorted(rows, key=lambda row: _float(row["training_seconds"]))
    x = np.asarray([_float(row["training_seconds"]) for row in ordered])
    y = np.asarray([_float(row[key]) for row in ordered])
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    return x, np.minimum.accumulate(y)


def _plot_frontier(
    frontier: list[dict[str, Any]], budgets: list[dict[str, Any]], key: str, filename: str
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.8), constrained_layout=True)
    for axis, benchmark in zip(axes, BENCHMARKS, strict=True):
        for variant in ("plain", "mcv"):
            for seed in (0, 1):
                prefix = [
                    row for row in frontier
                    if row["benchmark"] == benchmark and int(row["seed"]) == seed
                    and row["training_variant"] == "shared_prefix"
                ]
                body = [
                    row for row in frontier
                    if row["benchmark"] == benchmark and int(row["seed"]) == seed
                    and row["training_variant"] == variant
                ]
                x, y = _envelope(prefix + body, key)
                axis.plot(
                    x, y, color=COLORS[variant], alpha=0.35, linewidth=1.2,
                    linestyle="--", marker=".", markersize=3,
                )
            operating = [
                row for row in budgets
                if row["benchmark"] == benchmark
                and row["training_variant"] == variant
            ]
            for budget in (240.0, 480.0, 960.0):
                local = [row for row in operating if _float(row["budget_seconds"]) == budget]
                if not local:
                    continue
                values = [_float(row[key]) for row in local]
                times = [_float(row["actual_training_seconds"]) for row in local]
                axis.scatter(
                    [np.mean(times)], [np.mean(values)], s=42,
                    color=COLORS[variant], edgecolor="white", linewidth=0.7,
                    label=("Plain training" if variant == "plain" else "MCV training")
                    if budget == 240.0 else None,
                    zorder=3,
                )
        axis.set_title(LABELS[benchmark])
        axis.set_xlabel("Cumulative GH200 training seconds")
        axis.grid(alpha=0.2)
        axis.set_yscale("log")
    axes[0].set_ylabel("Relative L2" if key == "relative_l2" else "Systematic RMSE proxy")
    axes[0].legend(frameon=False, fontsize=9)
    fig.savefig(FIGURE_DIR / f"{filename}.png", dpi=220)
    fig.savefig(FIGURE_DIR / f"{filename}.pdf")
    plt.close(fig)


def _plot_label_variance() -> None:
    rows = [
        row for row in _csv(RESULTS_DIR / "label_variance_comparison.csv")
        if int(row["paths"]) == 128
    ]
    labels = [row["checkpoint"] for row in rows]
    value = [_float(row["value_variance_ratio_mcv_over_plain"]) for row in rows]
    gradient = [
        _float(row["gradient_component_variance_ratio_mcv_over_plain"]) for row in rows
    ]
    x = np.arange(len(rows))
    width = 0.34
    fig, axis = plt.subplots(figsize=(5.8, 3.8), constrained_layout=True)
    axis.bar(x - width / 2, value, width, label="Value label", color="#5DA5DA")
    axis.bar(x + width / 2, gradient, width, label="Gradient label", color="#006BA4")
    axis.axhline(1.0, color="black", linewidth=1.0, linestyle="--")
    for position, ratio in zip(x - width / 2, value, strict=True):
        axis.text(position, ratio + 0.012, f"{1 / ratio:.1f}×", ha="center", fontsize=8)
    for position, ratio in zip(x + width / 2, gradient, strict=True):
        axis.text(position, ratio + 0.012, f"{1 / ratio:.1f}×", ha="center", fontsize=8)
    axis.set_xticks(x, labels)
    axis.set_ylabel("MCV / Plain estimator variance")
    axis.set_ylim(0.0, 1.05)
    axis.legend(frameon=False)
    axis.grid(axis="y", alpha=0.2)
    fig.savefig(FIGURE_DIR / "label_variance_reduction.png", dpi=220)
    fig.savefig(FIGURE_DIR / "label_variance_reduction.pdf")
    plt.close(fig)


def _plot_actor_validation(rows: list[dict[str, Any]]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.0), constrained_layout=True)
    for column, benchmark in enumerate(BENCHMARKS):
        for variant in ("plain", "mcv"):
            local = [
                row for row in rows
                if row.get("diagnostic_kind") == "selected_training_trace"
                and row["benchmark"] == benchmark
                and row["training_variant"] == variant
                and int(row["training_points"]) == 8192
            ]
            for seed in (0, 1):
                seed_rows = sorted(
                    [row for row in local if int(row["seed"]) == seed],
                    key=lambda row: int(row["iteration"]),
                )
                if not seed_rows:
                    continue
                x = [int(row["iteration"]) for row in seed_rows]
                axes[0, column].plot(
                    x, [_float(row["actor_validation_value_rmse"]) for row in seed_rows],
                    color=COLORS[variant], alpha=0.55, marker="o", markersize=3,
                    label=variant.upper() if seed == 0 else None,
                )
                axes[1, column].plot(
                    x, [_float(row["actor_validation_gradient_vector_rmse"]) for row in seed_rows],
                    color=COLORS[variant], alpha=0.55, marker="o", markersize=3,
                    label=variant.upper() if seed == 0 else None,
                )
        axes[0, column].set_title(LABELS[benchmark])
        axes[1, column].set_xlabel("Actor iteration")
        for row in range(2):
            axes[row, column].grid(alpha=0.2)
            axes[row, column].set_yscale("log")
    axes[0, 0].set_ylabel("Validation Value RMSE")
    axes[1, 0].set_ylabel("Validation Gradient vector RMSE")
    axes[0, 0].legend(frameon=False)
    fig.savefig(FIGURE_DIR / "actor_validation_by_iteration.png", dpi=220)
    fig.savefig(FIGURE_DIR / "actor_validation_by_iteration.pdf")
    plt.close(fig)


def run() -> None:
    ensure_dirs()
    frontier_path = RESULTS_DIR / "actor_training_frontier.csv"
    budget_path = RESULTS_DIR / "budget_operating_points.csv"
    frontier = _csv(frontier_path)
    prefix = [row for row in frontier if row["training_variant"] == "shared_prefix"]
    if not prefix:
        prefix = _prefix_frontier()
    frontier = prefix + [row for row in frontier if row["training_variant"] != "shared_prefix"]
    write_csv(frontier_path, frontier)

    stage1_common = [
        row for row in _csv(RESULTS_DIR / "actor_validation_diagnostics.csv")
        if row.get("stage") == "1" and row.get("variant") in ("plain", "mcv")
    ]
    for row in stage1_common:
        row["diagnostic_kind"] = "common_high_quality_stage1"
        row["common_high_quality_target"] = True
    plain_diagnostics = _plain_training_runs()
    mcv_diagnostics, mcv_runtime = _mcv_training_runs()
    diagnostics = stage1_common + plain_diagnostics + mcv_diagnostics
    write_csv(RESULTS_DIR / "actor_validation_diagnostics.csv", diagnostics)
    runtime = _plain_runtime() + mcv_runtime
    write_csv(RESULTS_DIR / "runtime_and_resource_accounting.csv", runtime)

    budgets = _csv(budget_path)
    _plot_frontier(frontier, budgets, "relative_l2", "relative_l2_vs_training_seconds")
    _plot_frontier(
        frontier, budgets, "systematic_rmse_proxy", "systematic_rmse_vs_training_seconds"
    )
    _plot_label_variance()
    _plot_actor_validation(diagnostics)

    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True
    ).stdout.strip()
    manifest = {
        "experiment": "tangent_actor_mcv_training_high_budget_v1",
        "git_commit_at_aggregation": commit,
        "branch": "ir-mlp-tangent-actor-mcv-training",
        "python": platform.python_version(),
        "platform": platform.platform(),
        "benchmarks": list(BENCHMARKS),
        "training_seeds": [0, 1],
        "label_estimator": "direct alpha=1 P-H with full pathwise derivative",
        "common_actor_quality_inference": {"kind": "plain", "paths": 512, "steps": 32},
        "secondary_inference": {"kind": "MCV", "paths": 512, "steps": [16, 32, 64]},
        "budget_boundaries_seconds": [240, 480, 960],
        "ground_truth_used_for_training_or_selection": False,
        "source_results": {
            "strong_baselines": str(STRONG_RESULTS_DIR),
            "high_budget_plain": str(HIGH_BUDGET_RESULTS_DIR),
        },
        "deliverables": sorted(
            str(path.relative_to(RESULTS_DIR))
            for path in RESULTS_DIR.rglob("*")
            if path.is_file() and "slurm" not in path.parts
            and "checkpoints" not in path.parts and "labels" not in path.parts
            and "predictions" not in path.parts
        ),
    }
    atomic_json(RESULTS_DIR / "reproducibility_manifest.json", manifest)


if __name__ == "__main__":
    run()
