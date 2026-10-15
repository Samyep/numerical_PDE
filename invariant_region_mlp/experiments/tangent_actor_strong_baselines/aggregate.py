"""Aggregate measured frontiers, independently recompute metrics, and plot."""

from __future__ import annotations

from collections import defaultdict
import json
import math
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .common import (
    BUDGETS,
    FIGURE_DIR,
    RAW_DIR,
    RESULTS_DIR,
    TRAINING_SEEDS,
    atomic_json,
    ensure_dirs,
    evaluation_data,
    metrics,
    write_csv,
)


BENCHMARKS = ("nonseparable_d100_k2", "joint_concave_d400_k2")
TRAINING_METHODS = (
    "tangent_actor",
    "standard_tangent_actor",
    "terminal_tangent_Ag",
    "deep_bsde",
    "sdgd_pinn",
    "deep_picard",
    "scasml_l1",
    "full_history_mlp",
)
LABELS = {
    "tangent_actor": "Terminal-aware Tangent Actor",
    "standard_tangent_actor": "Standard Tangent Actor",
    "terminal_tangent_Ag": "Untrained terminal tangent",
    "deep_bsde": "Amortized Deep BSDE",
    "sdgd_pinn": "PINN + SDGD",
    "deep_picard": "Deep Picard",
    "scasml_l1": "SCaSML-L1",
    "full_history_mlp": "Full-History MLP",
}
COLORS = {
    "tangent_actor": "#0072B2",
    "standard_tangent_actor": "#56B4E9",
    "terminal_tangent_Ag": "#999999",
    "deep_bsde": "#D55E00",
    "sdgd_pinn": "#009E73",
    "deep_picard": "#CC79A7",
    "scasml_l1": "#E69F00",
    "full_history_mlp": "#000000",
}


def _read(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        return json.load(handle)


def load_frontiers() -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    training: list[dict[str, Any]] = []
    total: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    for benchmark in BENCHMARKS:
        for seed in TRAINING_SEEDS:
            payload = _read(RAW_DIR / f"evaluation_{benchmark}_s{seed}.json")
            training.extend(payload["training_frontier"])
            total.extend(payload["total_compute_frontier"])
            failures.extend(
                {"benchmark": benchmark, "seed": seed, **row}
                for row in payload["failures"]
            )
    return training, total, failures


def metric_recomputation_audit(
    training: list[dict[str, Any]], total: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    seen: set[tuple[str, int, str]] = set()
    for row in training + total:
        key = (row["benchmark"], int(row["seed"]), row.get("prediction_key", ""))
        if not key[2] or key in seen:
            continue
        seen.add(key)
        prediction_path = RESULTS_DIR / "predictions" / f"frontier_{key[0]}_s{key[1]}.npz"
        with np.load(prediction_path) as archive:
            prediction = archive[key[2]]
            reference = archive["test_reference"]
        recomputed = metrics(prediction, reference)
        audit = {
            "benchmark": key[0],
            "seed": key[1],
            "prediction_key": key[2],
            "method": row["method"],
            "stored_rmse": row["rmse"],
            "recomputed_rmse": recomputed["rmse"],
            "rmse_abs_difference": abs(row["rmse"] - recomputed["rmse"]),
            "stored_relative_l2": row["relative_l2"],
            "recomputed_relative_l2": recomputed["relative_l2"],
            "relative_l2_abs_difference": abs(
                row["relative_l2"] - recomputed["relative_l2"]
            ),
            "point_count": len(prediction),
        }
        audit["passed_1e-12"] = (
            audit["rmse_abs_difference"] <= 1e-12
            and audit["relative_l2_abs_difference"] <= 1e-12
        )
        rows.append(audit)
    return rows


def budget_alignment(training: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str, int], list[dict[str, Any]]] = defaultdict(list)
    for row in training:
        grouped[(row["benchmark"], row["method"], int(row["seed"]))].append(row)
    output: list[dict[str, Any]] = []
    for benchmark in BENCHMARKS:
        for method in TRAINING_METHODS:
            for seed in TRAINING_SEEDS:
                candidates_all = grouped.get((benchmark, method, seed), [])
                for level, target in enumerate(BUDGETS, 1):
                    if method == "full_history_mlp":
                        output.append(
                            {
                                "benchmark": benchmark,
                                "method": method,
                                "seed": seed,
                                "budget_level": f"B{level}",
                                "target_seconds": target,
                                "available": False,
                                "reason": "inference-only method; training budget is not applicable",
                            }
                        )
                        continue
                    eligible = [
                        row
                        for row in candidates_all
                        if 0.75 * target <= float(row["training_seconds"]) <= 1.25 * target
                    ]
                    if not eligible:
                        output.append(
                            {
                                "benchmark": benchmark,
                                "method": method,
                                "seed": seed,
                                "budget_level": f"B{level}",
                                "target_seconds": target,
                                "available": False,
                                "reason": "no complete real checkpoint within +/-25%",
                            }
                        )
                        continue
                    # Selection uses elapsed-time proximity only; the final-test metric
                    # is copied after the checkpoint identity has been fixed.
                    selected = min(
                        eligible,
                        key=lambda row: (
                            abs(float(row["training_seconds"]) - target),
                            float(row["validation_objective"])
                            if row.get("validation_objective") is not None
                            else math.inf,
                        ),
                    )
                    output.append(
                        {
                            "benchmark": benchmark,
                            "method": method,
                            "seed": seed,
                            "budget_level": f"B{level}",
                            "target_seconds": target,
                            "available": True,
                            "actual_seconds": selected["training_seconds"],
                            "deviation_percent": 100.0
                            * (float(selected["training_seconds"]) - target)
                            / target,
                            "checkpoint": selected["checkpoint"],
                            "selection_rule": "closest real checkpoint within +/-25%; validation objective tie-break",
                            "validation_objective": selected.get("validation_objective"),
                            "relative_l2_posthoc": selected["relative_l2"],
                            "rmse_posthoc": selected["rmse"],
                            "rmse_over_reference_std_posthoc": selected[
                                "rmse_over_reference_std"
                            ],
                        }
                    )
    return output


def total_budget_alignment(total: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str, int], list[dict[str, Any]]] = defaultdict(list)
    for row in total:
        grouped[(row["benchmark"], row["method"], int(row["seed"]))].append(row)
    output: list[dict[str, Any]] = []
    for benchmark in BENCHMARKS:
        for method in TRAINING_METHODS:
            for seed in TRAINING_SEEDS:
                method_rows = grouped.get((benchmark, method, seed), [])
                for level, target in enumerate(BUDGETS, 1):
                    eligible = [
                        row
                        for row in method_rows
                        if 0.75 * target <= float(row["total_seconds_N960"]) <= 1.25 * target
                    ]
                    if not eligible:
                        output.append(
                            {
                                "benchmark": benchmark,
                                "method": method,
                                "seed": seed,
                                "budget_level": f"B{level}",
                                "target_total_seconds": target,
                                "available": False,
                                "reason": "no real training-plus-inference point within +/-25%",
                            }
                        )
                        continue
                    # Tune inference variants only *within the same trained
                    # checkpoint* using the common value-validation bank.  A
                    # tangent-label fit loss and a value Relative-L2 have
                    # different units and must never rank different Actor
                    # checkpoints.  Once each checkpoint has one variant, use
                    # total-time proximity across checkpoints.
                    by_checkpoint: dict[str, list[dict[str, Any]]] = defaultdict(list)
                    for row in eligible:
                        by_checkpoint[str(row["checkpoint"])].append(row)
                    representatives: list[dict[str, Any]] = []
                    for checkpoint_rows in by_checkpoint.values():
                        validated = [
                            row
                            for row in checkpoint_rows
                            if row.get("validation_relative_l2") is not None
                            and math.isfinite(float(row["validation_relative_l2"]))
                        ]
                        if validated:
                            representatives.append(
                                min(
                                    validated,
                                    key=lambda row: (
                                        float(row["validation_relative_l2"]),
                                        abs(float(row["total_seconds_N960"]) - target),
                                    ),
                                )
                            )
                        else:
                            representatives.append(
                                min(
                                    checkpoint_rows,
                                    key=lambda row: abs(
                                        float(row["total_seconds_N960"]) - target
                                    ),
                                )
                            )
                    selected = min(
                        representatives,
                        key=lambda row: abs(float(row["total_seconds_N960"]) - target),
                    )
                    output.append(
                        {
                            "benchmark": benchmark,
                            "method": method,
                            "seed": seed,
                            "budget_level": f"B{level}",
                            "target_total_seconds": target,
                            "available": True,
                            "actual_training_seconds": selected["training_seconds"],
                            "actual_inference_seconds_N960": selected[
                                "inference_seconds_N960"
                            ],
                            "actual_total_seconds_N960": selected["total_seconds_N960"],
                            "deviation_percent": 100.0
                            * (float(selected["total_seconds_N960"]) - target)
                            / target,
                            "checkpoint": selected["checkpoint"],
                            "total_variant": selected.get("total_variant"),
                            "selection_rule": "within-checkpoint inference budget selected by common validation Relative-L2; closest real total time across checkpoints",
                            "relative_l2_posthoc": selected["relative_l2"],
                            "rmse_posthoc": selected["rmse"],
                            "rmse_over_reference_std_posthoc": selected[
                                "rmse_over_reference_std"
                            ],
                        }
                    )
    return output


def seed_summary(
    rows: list[dict[str, Any]], *, time_key: str, include_variant: bool
) -> list[dict[str, Any]]:
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        key: tuple[Any, ...] = (
            row["benchmark"],
            row["method"],
            row["checkpoint"],
        )
        if include_variant:
            key = key + (row.get("total_variant"),)
        grouped[key].append(row)
    output: list[dict[str, Any]] = []
    for key, values in sorted(grouped.items(), key=lambda item: tuple(map(str, item[0]))):
        relative = np.asarray([float(row["relative_l2"]) for row in values])
        rmse = np.asarray([float(row["rmse"]) for row in values])
        times = np.asarray([float(row[time_key]) for row in values])
        result = {
            "benchmark": key[0],
            "method": key[1],
            "checkpoint": key[2],
            "seeds": len(values),
            f"{time_key}_mean": float(np.mean(times)),
            f"{time_key}_std": float(np.std(times, ddof=1)) if len(times) > 1 else 0.0,
            "relative_l2_mean": float(np.mean(relative)),
            "relative_l2_std": float(np.std(relative, ddof=1)) if len(relative) > 1 else 0.0,
            "rmse_mean": float(np.mean(rmse)),
            "rmse_std": float(np.std(rmse, ddof=1)) if len(rmse) > 1 else 0.0,
        }
        if include_variant:
            result["total_variant"] = key[3]
        output.append(result)
    return output


def convergence_diagnostics(training: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    grouped: dict[tuple[str, str, int], list[dict[str, Any]]] = defaultdict(list)
    for row in training:
        grouped[(row["benchmark"], row["method"], int(row["seed"]))].append(row)
    for (benchmark, method, seed), values in sorted(grouped.items()):
        values = sorted(values, key=lambda row: float(row["training_seconds"]))
        first, last = values[0], values[-1]
        previous = values[-2] if len(values) > 1 else None
        if previous is None:
            rmse_change = relative_change = None
            improving = None
        else:
            rmse_change = float(last["rmse"] - previous["rmse"])
            relative_change = rmse_change / max(float(previous["rmse"]), 1e-12)
            improving = bool(relative_change < -0.05)
        nonfinite = max(
            int(value.get("nonfinite_count_training", 0) or 0)
            + int(value.get("nonfinite_count_inference", 0) or 0)
            for value in values
        )
        if method in {"terminal_tangent_Ag"}:
            status = "untrained ablation"
        elif previous is None:
            status = "single measured point"
        elif nonfinite:
            status = "numerically unstable"
        elif abs(relative_change or 0.0) <= 0.05:
            status = "reasonably stable over last measured checkpoint"
        elif improving:
            status = "still improving; compute-limited"
        else:
            status = "not monotone / unsettled"
        rows.append(
            {
                "benchmark": benchmark,
                "method": method,
                "seed": seed,
                "checkpoints": len(values),
                "first_training_seconds": first["training_seconds"],
                "last_training_seconds": last["training_seconds"],
                "first_rmse_posthoc": first["rmse"],
                "last_rmse_posthoc": last["rmse"],
                "last_checkpoint_rmse_change_posthoc": rmse_change,
                "last_checkpoint_relative_rmse_change_posthoc": relative_change,
                "validation_objective_first": first.get("validation_objective"),
                "validation_objective_last": last.get("validation_objective"),
                "nonfinite_count_max": nonfinite,
                "still_improving_at_largest_budget_posthoc": improving,
                "stability_assessment": status,
                "posthoc_test_used_for_selection": False,
            }
        )
    return rows


def runtime_accounting(
    training: list[dict[str, Any]], total: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for benchmark in BENCHMARKS:
        for method in TRAINING_METHODS:
            for seed in TRAINING_SEEDS:
                train_rows = [
                    row
                    for row in training
                    if row["benchmark"] == benchmark
                    and row["method"] == method
                    and int(row["seed"]) == seed
                ]
                total_rows = [
                    row
                    for row in total
                    if row["benchmark"] == benchmark
                    and row["method"] == method
                    and int(row["seed"]) == seed
                ]
                if not train_rows and not total_rows:
                    continue
                rows.append(
                    {
                        "benchmark": benchmark,
                        "method": method,
                        "seed": seed,
                        "maximum_training_seconds": max(
                            [float(row["training_seconds"]) for row in train_rows]
                            or [0.0]
                        ),
                        "number_training_checkpoints_evaluated": len(train_rows),
                        "number_total_compute_variants_evaluated": len(total_rows),
                        "sum_measured_test_inference_seconds": sum(
                            float(row["inference_seconds_N960"]) for row in total_rows
                        ),
                        "minimum_total_seconds_N960": min(
                            [float(row["total_seconds_N960"]) for row in total_rows]
                            or [math.nan]
                        ),
                        "maximum_total_seconds_N960": max(
                            [float(row["total_seconds_N960"]) for row in total_rows]
                            or [math.nan]
                        ),
                    }
                )
    return rows


def _plot_frontier(
    rows: list[dict[str, Any]], *, x_key: str, filename: str, training_guides: bool
) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(13.2, 5.0), sharey=False)
    for axis, benchmark in zip(axes, BENCHMARKS, strict=True):
        local = [row for row in rows if row["benchmark"] == benchmark]
        for method in TRAINING_METHODS:
            method_rows = [row for row in local if row["method"] == method]
            if not method_rows:
                continue
            for seed in TRAINING_SEEDS:
                seed_rows = sorted(
                    [row for row in method_rows if int(row["seed"]) == seed],
                    key=lambda row: float(row[x_key]),
                )
                if not seed_rows:
                    continue
                x = np.asarray([max(float(row[x_key]), 0.1) for row in seed_rows])
                y = np.asarray([float(row["relative_l2"]) for row in seed_rows])
                if method == "full_history_mlp":
                    # These points are distinct predeclared (level, sample-base)
                    # configurations, rather than checkpoints on one training
                    # trajectory.  Connecting them would suggest an ordering
                    # that the algorithm does not have.
                    axis.scatter(
                        x,
                        y,
                        marker="D",
                        s=18,
                        alpha=0.55,
                        color=COLORS[method],
                        label=LABELS[method] if seed == 0 else None,
                    )
                else:
                    axis.plot(
                        x,
                        y,
                        marker="o",
                        markersize=3.3,
                        linewidth=1.25,
                        alpha=0.55,
                        color=COLORS[method],
                        label=LABELS[method] if seed == 0 else None,
                    )
            method_groups: dict[tuple[str, str | None], list[dict[str, Any]]] = defaultdict(list)
            for row in method_rows:
                method_groups[(row["checkpoint"], row.get("total_variant"))].append(row)
            for group in method_groups.values():
                if len(group) < 2:
                    continue
                mean_x = float(np.mean([max(float(row[x_key]), 0.1) for row in group]))
                mean_y = float(np.mean([float(row["relative_l2"]) for row in group]))
                std_y = float(np.std([float(row["relative_l2"]) for row in group], ddof=1))
                axis.errorbar(
                    [mean_x],
                    [mean_y],
                    yerr=[std_y],
                    fmt="D",
                    markersize=3.5,
                    capsize=2,
                    linewidth=0.9,
                    color=COLORS[method],
                    alpha=0.9,
                )
        if training_guides:
            for budget in BUDGETS:
                axis.axvline(budget, color="#bbbbbb", linewidth=0.7, linestyle=":", zorder=0)
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.grid(True, which="both", alpha=0.22)
        axis.set_title(benchmark.replace("_", " "))
        axis.set_xlabel("Actual GH200 seconds")
        axis.set_ylabel("Relative $L^2$ on 960 held-out states")
    handles, labels = axes[0].get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    figure.legend(
        by_label.values(), by_label.keys(), loc="upper center", ncol=4, frameon=False
    )
    figure.tight_layout(rect=(0, 0, 1, 0.88))
    figure.savefig(FIGURE_DIR / filename, dpi=220)
    figure.savefig(FIGURE_DIR / filename.replace(".png", ".pdf"))
    plt.close(figure)


def _plot_method_convergence(training: list[dict[str, Any]]) -> None:
    for benchmark in BENCHMARKS:
        methods = sorted({row["method"] for row in training if row["benchmark"] == benchmark})
        for method in methods:
            rows = [
                row
                for row in training
                if row["benchmark"] == benchmark and row["method"] == method
            ]
            figure, axes = plt.subplots(1, 2, figsize=(9.2, 3.7))
            for seed in TRAINING_SEEDS:
                local = sorted(
                    [row for row in rows if int(row["seed"]) == seed],
                    key=lambda row: float(row["training_seconds"]),
                )
                if not local:
                    continue
                time_values = [max(float(row["training_seconds"]), 0.1) for row in local]
                axes[0].plot(
                    time_values,
                    [row["rmse"] for row in local],
                    marker="o",
                    label=f"seed {seed}",
                )
                validation = [row.get("validation_objective") for row in local]
                if all(value is not None and float(value) > 0 for value in validation):
                    axes[1].plot(time_values, validation, marker="o", label=f"seed {seed}")
            axes[0].set_ylabel("Test Value RMSE (post hoc)")
            axes[1].set_ylabel("Reference-free validation objective")
            for axis in axes:
                axis.set_xscale("log")
                axis.set_yscale("log")
                axis.set_xlabel("Actual training seconds")
                axis.grid(True, which="both", alpha=0.25)
                handles, labels = axis.get_legend_handles_labels()
                if handles:
                    axis.legend(frameon=False)
            figure.suptitle(f"{LABELS[method]} — {benchmark}")
            figure.tight_layout()
            stem = f"convergence_{benchmark}_{method}"
            figure.savefig(FIGURE_DIR / f"{stem}.png", dpi=200)
            plt.close(figure)


def main() -> None:
    ensure_dirs()
    training, total, failures = load_frontiers()
    audit = metric_recomputation_audit(training, total)
    if not all(row["passed_1e-12"] for row in audit):
        raise AssertionError("independent metric recomputation failed")
    alignment = budget_alignment(training)
    total_alignment = total_budget_alignment(total)
    convergence = convergence_diagnostics(training)
    runtime = runtime_accounting(training, total)
    training_seed_summary = seed_summary(
        training, time_key="training_seconds", include_variant=False
    )
    total_seed_summary = seed_summary(
        total, time_key="total_seconds_N960", include_variant=True
    )
    write_csv(RESULTS_DIR / "training_frontier.csv", training)
    write_csv(RESULTS_DIR / "total_compute_frontier.csv", total)
    write_csv(RESULTS_DIR / "compute_budget_alignment.csv", alignment)
    write_csv(RESULTS_DIR / "total_compute_budget_alignment.csv", total_alignment)
    write_csv(RESULTS_DIR / "convergence_diagnostics.csv", convergence)
    write_csv(RESULTS_DIR / "runtime_accounting.csv", runtime)
    write_csv(RESULTS_DIR / "metric_recomputation_audit.csv", audit)
    write_csv(RESULTS_DIR / "evaluation_failures.csv", failures)
    write_csv(RESULTS_DIR / "training_seed_summary.csv", training_seed_summary)
    write_csv(RESULTS_DIR / "total_compute_seed_summary.csv", total_seed_summary)
    atomic_json(
        RESULTS_DIR / "metric_recomputation_summary.json",
        {
            "rows": len(audit),
            "all_passed_1e-12": all(row["passed_1e-12"] for row in audit),
            "maximum_rmse_abs_difference": max(row["rmse_abs_difference"] for row in audit),
            "maximum_relative_l2_abs_difference": max(
                row["relative_l2_abs_difference"] for row in audit
            ),
        },
    )
    _plot_frontier(
        training,
        x_key="training_seconds",
        filename="accuracy_vs_training_time.png",
        training_guides=True,
    )
    _plot_frontier(
        total,
        x_key="total_seconds_N960",
        filename="accuracy_vs_total_compute_N960.png",
        training_guides=False,
    )
    _plot_method_convergence(training)
    print(
        json.dumps(
            {
                "training_rows": len(training),
                "total_rows": len(total),
                "alignment_rows": len(alignment),
                "total_alignment_rows": len(total_alignment),
                "audit_rows": len(audit),
                "failures": len(failures),
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
