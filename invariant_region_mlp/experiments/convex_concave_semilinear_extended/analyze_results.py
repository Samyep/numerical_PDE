"""Aggregate complete CSVs and figures without introducing new model selection."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .common import FIGURE_DIR, RESULTS_DIR, atomic_json, write_csv


FAMILIES = ("allen_cahn", "joint_concave", "joint_convex")
LABELS = {
    "allen_cahn": "Allen–Cahn",
    "joint_concave": "Joint concave",
    "joint_convex": "Convex via sign",
    "nonseparable": "Nonseparable",
}


def _load_jsons(prefix: str) -> list[dict[str, Any]]:
    result = []
    for path in sorted(RESULTS_DIR.glob(f"{prefix}_*.json")):
        with path.open() as handle:
            result.append(json.load(handle))
    return result


def _configuration_summary(payload: dict[str, Any]) -> dict[str, Any]:
    spec, rows, signal = payload["spec"], payload["iteration_rows"], payload["terminal_signal"]
    by_actor = {row["iteration"]: row for row in rows}
    prediction_path = RESULTS_DIR.parent.parent / payload["prediction_path"]
    with np.load(prediction_path) as arrays:
        total_delta = arrays["a3_working_value"] - arrays["a0_working_value"]
    a3 = by_actor["a3"]
    phase_seconds = sum(row["label_seconds"] + row["training_seconds"] for row in rows)
    return {
        **spec,
        "run_id": a3["run_id"],
        "terminal_variance": signal["terminal_variance"],
        "terminal_std": signal["terminal_std"],
        "terminal_gradient_energy": signal["gradient_energy"],
        "terminal_gradient_rank": signal["ambient_gradient_matrix_rank"],
        "a0_working_value_mean": by_actor["a0"]["reference_free_working_value_mean"],
        "a3_working_value_mean": a3["reference_free_working_value_mean"],
        "a3_minus_a0_paired_mean": float(np.mean(total_delta)),
        "a3_minus_a0_paired_se_across_states": float(
            np.std(total_delta, ddof=1) / np.sqrt(len(total_delta))
        ),
        "improvement_over_terminal_std": float(
            -np.mean(total_delta) / max(signal["terminal_std"], np.finfo(float).tiny)
        ),
        "a3_actor_validation_value_rmse": a3["actor_validation_value_rmse"],
        "a3_actor_validation_gradient_vector_rmse": a3[
            "actor_validation_gradient_vector_rmse"
        ],
        "a3_value_label_replication_rmse": a3["value_label_replication_rmse"],
        "a3_gradient_label_replication_vector_rmse": a3[
            "gradient_label_replication_vector_rmse"
        ],
        "a3_policy_mc_mean_se": a3["reference_free_value_mc_se"],
        "a3_label_weight_ess_fraction": a3["label_weight_ess_fraction"],
        "a3_evaluation_weight_ess_fraction": a3["evaluation_weight_ess_fraction"],
        "a3_evaluation_weight_max": a3["evaluation_weight_max"],
        "total_label_and_training_seconds": phase_seconds,
        "total_reference_free_inference_seconds": sum(row["inference_seconds"] for row in rows),
        "all_rollouts_finite": all(row["evaluation_nonfinite_count"] == 0 for row in rows),
    }


def _plot_effective_dimension(summaries: list[dict[str, Any]]) -> None:
    rows = [row for row in summaries if row["family"] in FAMILIES]
    figure, axes = plt.subplots(2, 2, figsize=(10.8, 7.8))
    styles = {100: "o-", 400: "s--"}
    for family in FAMILIES:
        for d in sorted({row["d"] for row in rows}):
            subset = sorted(
                (row for row in rows if row["family"] == family and row["d"] == d),
                key=lambda row: row["k"],
            )
            if not subset:
                continue
            x = [row["k"] for row in subset]
            label = f"{LABELS[family]}, d={d}"
            style = styles.get(d, "o-")
            axes[0, 0].plot(
                x, [row["improvement_over_terminal_std"] for row in subset], style, label=label
            )
            axes[0, 1].plot(
                x,
                [row["a3_actor_validation_value_rmse"] / row["terminal_std"] for row in subset],
                style,
                label=label,
            )
            gradient = [row["a3_actor_validation_gradient_vector_rmse"] for row in subset]
            if any(value is not None for value in gradient):
                axes[1, 0].plot(x, gradient, style, label=label)
            axes[1, 1].plot(
                x, [row["a3_evaluation_weight_ess_fraction"] for row in subset], style, label=label
            )
    axes[0, 0].set_title("Paired a0→a3 improvement / terminal std")
    axes[0, 1].set_title("a3 Value-Actor validation RMSE / terminal std")
    axes[1, 0].set_title("a3 Gradient-Actor validation vector RMSE")
    axes[1, 1].set_title("a3 evaluation exponential-weight ESS")
    for axis in axes.flat:
        axis.set_xlabel("effective dimension k")
        axis.set_xscale("log", base=2)
        axis.set_xticks((2, 5, 10, 20))
        axis.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        axis.grid(alpha=0.25)
    axes[0, 0].axhline(0.0, color="0.4", lw=1)
    axes[1, 1].set_ylim(0.0, 1.03)
    axes[0, 1].legend(frameon=False, fontsize=7, ncol=2)
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "accuracy_stability_vs_effective_dimension.png", dpi=180)
    plt.close(figure)


def _plot_signal(summaries: list[dict[str, Any]]) -> None:
    rows = [
        row
        for row in summaries
        if row["family"] == "allen_cahn" and row["d"] == 100 and row["training_seed"] == 0
    ]
    rows.sort(key=lambda row: row["k"])
    figure, axes = plt.subplots(1, 2, figsize=(8.8, 3.6))
    x = [row["k"] for row in rows]
    axes[0].plot(x, [row["terminal_variance"] for row in rows], "o-")
    axes[1].plot(x, [row["terminal_gradient_energy"] for row in rows], "o-")
    axes[0].set_title("Terminal signal variance")
    axes[1].set_title(r"Terminal gradient energy $E\|\nabla g\|^2$")
    for axis in axes:
        axis.set_xlabel("effective dimension k")
        axis.set_xscale("log", base=2)
        axis.set_yscale("log")
        axis.set_xticks((2, 5, 10, 20))
        axis.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        axis.grid(alpha=0.25)
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "terminal_signal_scaling.png", dpi=180)
    plt.close(figure)


def _plot_nonseparable(
    iteration_rows: list[dict[str, Any]], ground_truth_rows: list[dict[str, Any]]
) -> None:
    rf = [row for row in iteration_rows if row["family"] == "nonseparable"]
    gt = [
        row
        for row in ground_truth_rows
        if row["family"] == "nonseparable" and row["iteration"] != "exact_tangent"
    ]
    if not rf or not gt:
        return
    order = ("a0", "a1", "a2", "a3")
    figure, axes = plt.subplots(1, 2, figsize=(9.4, 3.7))
    for seed in sorted({row["training_seed"] for row in rf}):
        subset = {row["iteration"]: row for row in rf if row["training_seed"] == seed}
        base = subset["a0"]["reference_free_working_value_mean"]
        axes[0].plot(
            order,
            [subset[name]["reference_free_working_value_mean"] - base for name in order],
            "o-",
            label=f"seed {seed}",
        )
        exact = {row["iteration"]: row for row in gt if row["training_seed"] == seed}
        axes[1].plot(order, [exact[name]["rmse"] for name in order], "o-", label=f"seed {seed}")
    axes[0].set_title("Reference-free paired tangent value")
    axes[0].set_ylabel("change from a0 (working variable)")
    axes[1].set_title("Independent 2D Ground Truth")
    axes[1].set_ylabel("Value RMSE")
    axes[1].set_yscale("log")
    for axis in axes:
        axis.set_xlabel("Actor iteration")
        axis.grid(alpha=0.25)
        axis.legend(frameon=False, fontsize=8)
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "nonseparable_driver_iterations.png", dpi=180)
    plt.close(figure)


def _plot_k2_ground_truth(rows: list[dict[str, Any]]) -> None:
    subset = [row for row in rows if row["iteration"] != "exact_tangent"]
    if not subset:
        return
    figure, axis = plt.subplots(figsize=(7.4, 4.2))
    order = ("a0", "a1", "a2", "a3")
    for family in ("nonseparable", *FAMILIES):
        family_rows = [row for row in subset if row["family"] == family]
        for d in sorted({row["d"] for row in family_rows}):
            selected = [row for row in family_rows if row["d"] == d]
            values = []
            spreads = []
            for actor in order:
                samples = [row["rmse"] for row in selected if row["iteration"] == actor]
                values.append(float(np.mean(samples)))
                spreads.append(float(np.std(samples, ddof=1)) if len(samples) > 1 else 0.0)
            axis.errorbar(
                order,
                values,
                yerr=spreads,
                marker="o",
                capsize=3,
                label=f"{LABELS[family]}, d={d}",
            )
    axis.set_yscale("log")
    axis.set_xlabel("Actor iteration")
    axis.set_ylabel("Value RMSE against independent 2D reference")
    axis.grid(alpha=0.25)
    axis.legend(frameon=False, fontsize=8, ncol=2)
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "k2_ground_truth_accuracy.png", dpi=180)
    plt.close(figure)


def main() -> None:
    training_payloads = _load_jsons("training")
    ground_truth_payloads = _load_jsons("ground_truth")
    if not training_payloads:
        raise RuntimeError("no completed training results")
    iteration_rows = [row for payload in training_payloads for row in payload["iteration_rows"]]
    summaries = [_configuration_summary(payload) for payload in training_payloads]
    ground_truth_rows = [row for payload in ground_truth_payloads for row in payload["rows"]]
    bregman_rows = [
        {**payload["spec"], **payload["bregman"]} for payload in ground_truth_payloads
    ]
    sort_key = lambda row: (
        ("nonseparable", *FAMILIES).index(row["family"]),
        row["d"],
        row["k"],
        row["training_seed"],
        ("a0", "a1", "a2", "a3", "exact_tangent").index(row.get("iteration", "a0")),
    )
    iteration_rows.sort(key=sort_key)
    summaries.sort(key=sort_key)
    ground_truth_rows.sort(key=sort_key)
    bregman_rows.sort(key=sort_key)
    write_csv(RESULTS_DIR / "iteration_results.csv", iteration_rows)
    write_csv(RESULTS_DIR / "configuration_summary.csv", summaries)
    if ground_truth_rows:
        write_csv(RESULTS_DIR / "ground_truth_results.csv", ground_truth_rows)
        write_csv(RESULTS_DIR / "bregman_identity.csv", bregman_rows)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    _plot_effective_dimension(summaries)
    _plot_signal(summaries)
    _plot_nonseparable(iteration_rows, ground_truth_rows)
    _plot_k2_ground_truth(ground_truth_rows)

    gate = {}
    for family in FAMILIES:
        matches = [
            row
            for row in summaries
            if row["family"] == family and row["d"] == 100 and row["k"] == 10
        ]
        if matches:
            row = matches[0]
            upper = row["a3_minus_a0_paired_mean"] + 1.96 * row[
                "a3_minus_a0_paired_se_across_states"
            ]
            gate[family] = {
                "paired_delta_95pct_upper": upper,
                "significant_improvement": upper < 0.0,
                "finite": row["all_rollouts_finite"],
                "ess_above_half": row["a3_evaluation_weight_ess_fraction"] > 0.5,
                "passes": upper < 0.0
                and row["all_rollouts_finite"]
                and row["a3_evaluation_weight_ess_fraction"] > 0.5,
            }
    atomic_json(
        RESULTS_DIR / "analysis_summary.json",
        {
            "completed_configurations": len(summaries),
            "completed_ground_truth_runs": len(ground_truth_payloads),
            "initial_stage_complete": len(training_payloads) >= 12,
            "conditional_expansion_complete": len(training_payloads) >= 27,
            "reference_free_expansion_gate": gate,
            "all_terminal_gradient_ranks_verified": all(
                row["terminal_gradient_rank"] == row["k"] for row in summaries
            ),
            "all_rollouts_finite": all(row["all_rollouts_finite"] for row in summaries),
            "plots": sorted(path.name for path in FIGURE_DIR.glob("*.png")),
        },
    )
    print(
        f"aggregated {len(summaries)} configurations, "
        f"{len(ground_truth_payloads)} Ground Truth runs"
    )


if __name__ == "__main__":
    main()
