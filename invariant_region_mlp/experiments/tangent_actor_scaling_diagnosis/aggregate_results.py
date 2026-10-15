"""Create complete CSVs and diagnostic plots from frozen run artifacts."""

from __future__ import annotations

from collections import defaultdict
import json
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from .common import FIGURE_DIR, RESULTS_DIR, atomic_json, ensure_dirs, write_csv


def _payloads() -> list[dict[str, Any]]:
    result = []
    for path in sorted(RESULTS_DIR.glob("run_*.json")):
        with path.open() as handle:
            result.append(json.load(handle))
    if not result:
        raise RuntimeError("no completed scaling-diagnosis runs")
    return result


def _noise_corrected(fit: float | None, replication: float | None) -> float | None:
    if fit is None or replication is None:
        return None
    return float(np.sqrt(max(fit**2 - 0.5 * replication**2, 0.0)))


def _ablation_row(payload: dict[str, Any]) -> dict[str, Any]:
    rows = payload["iteration_rows"]
    final = rows[-1]
    phase = payload["phases"][-1]
    total_label = sum(row["label_seconds"] for row in rows)
    total_training = sum(row["training_seconds"] for row in rows)
    total_inference = sum(row["inference_seconds"] for row in rows)
    improvement = final["improvement_from_a0_mean"]
    se = final["delta_from_a0_se_across_states"]
    value_fit = final["validation_value_rmse"]
    value_noise = final["value_label_replication_rmse"]
    gradient_fit = final["validation_gradient_vector_rmse"]
    gradient_noise = final["gradient_label_replication_vector_rmse"]
    best = max(rows, key=lambda row: row["improvement_from_a0_mean"])
    return {
        **payload["spec"],
        **payload["variant"],
        "run_id": final["run_id"],
        "terminal_variance": payload["terminal_signal"]["terminal_variance"],
        "terminal_gradient_energy": payload["terminal_signal"]["gradient_energy"],
        "final_iteration": final["iteration"],
        "final_improvement_from_a0": improvement,
        "final_improvement_se_across_states": se,
        "final_improvement_ci95_low": improvement - 1.96 * se,
        "final_improvement_ci95_high": improvement + 1.96 * se,
        "final_improvement_state_paired_significant": improvement - 1.96 * se > 0.0,
        "best_reported_iteration": best["iteration"],
        "best_reported_improvement_from_a0": best["improvement_from_a0_mean"],
        "final_working_value_mean": final["working_value_mean"],
        "final_selected_step": final["selected_step"],
        "final_train_value_rmse": final["train_value_rmse"],
        "final_validation_value_rmse": value_fit,
        "final_value_generalization_gap_rmse": final["value_generalization_gap_rmse"],
        "final_value_label_replication_rmse": value_noise,
        "final_value_fit_over_replication_noise": value_fit
        / max(value_noise, np.finfo(float).tiny),
        "final_value_noise_corrected_fit_proxy": _noise_corrected(value_fit, value_noise),
        "final_train_gradient_component_rmse": final["train_gradient_component_rmse"],
        "final_validation_gradient_component_rmse": final[
            "validation_gradient_component_rmse"
        ],
        "final_train_gradient_vector_rmse": final["train_gradient_vector_rmse"],
        "final_validation_gradient_vector_rmse": gradient_fit,
        "final_gradient_generalization_gap_vector_rmse": final[
            "gradient_generalization_gap_vector_rmse"
        ],
        "final_gradient_label_replication_component_rmse": final[
            "gradient_label_replication_component_rmse"
        ],
        "final_gradient_label_replication_vector_rmse": gradient_noise,
        "final_gradient_fit_over_replication_noise": None
        if gradient_fit is None
        else gradient_fit / max(gradient_noise, np.finfo(float).tiny),
        "final_gradient_noise_corrected_fit_proxy": _noise_corrected(
            gradient_fit, gradient_noise
        ),
        "final_label_weight_ess_fraction": final["label_weight_ess_fraction"],
        "final_evaluation_weight_ess_fraction": final[
            "evaluation_weight_ess_fraction"
        ],
        "final_mc_mean_se": final["mc_mean_se"],
        "selected_lower_saturation_fraction": phase["training"][
            "selected_lower_saturation_fraction"
        ],
        "selected_upper_saturation_fraction": phase["training"][
            "selected_upper_saturation_fraction"
        ],
        "actor_parameters": final["actor_parameters"],
        "attributed_label_seconds": total_label,
        "training_seconds": total_training,
        "inference_seconds": total_inference,
        "attributed_total_gpu_seconds": total_label + total_training + total_inference,
        "observed_restart_aware_run_wall_seconds": payload["run_wall_seconds"],
        "improvement_per_attributed_gpu_minute": improvement
        / max((total_label + total_training + total_inference) / 60.0, np.finfo(float).tiny),
        "all_rollouts_finite": all(row["evaluation_nonfinite_count"] == 0 for row in rows),
        "exact_test_used_for_selection": payload["selection_used_exact_test"],
    }


def _flatten_curves(payloads: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for payload in payloads:
        spec, variant = payload["spec"], payload["variant"]
        for phase in payload["phases"]:
            labels = phase["labels"]
            for point in phase["trace"]:
                rows.append(
                    {
                        "run_id": phase["run_id"],
                        **spec,
                        **variant,
                        "iteration": f"a{phase['iteration']}",
                        "value_label_replication_rmse": labels[
                            "value_replication_rmse"
                        ],
                        "gradient_label_replication_component_rmse": labels[
                            "gradient_component_replication_rmse"
                        ],
                        "gradient_label_replication_vector_rmse": labels[
                            "gradient_vector_replication_rmse"
                        ],
                        "train_value_label_std": phase["training"][
                            "train_value_label_std"
                        ],
                        "validation_value_label_std": phase["training"][
                            "validation_value_label_std"
                        ],
                        "train_gradient_label_component_rms": phase["training"][
                            "train_gradient_label_component_rms"
                        ],
                        "validation_gradient_label_component_rms": phase["training"][
                            "validation_gradient_label_component_rms"
                        ],
                        "selected_step": phase["training"]["selected_step"],
                        **point,
                    }
                )
    return rows


def _compute_rows(payloads: list[dict[str, Any]]) -> list[dict[str, Any]]:
    result = []
    for payload in payloads:
        base = {**payload["spec"], **payload["variant"]}
        for row in payload["iteration_rows"]:
            result.append(
                {
                    **base,
                    "run_id": row["run_id"],
                    "iteration": row["iteration"],
                    "actor_parameters": row["actor_parameters"],
                    "attributed_label_seconds": row["label_seconds"],
                    "training_seconds": row["training_seconds"],
                    "inference_seconds": row["inference_seconds"],
                    "peak_training_gpu_bytes": None
                    if row["iteration"] == "a0"
                    else payload["phases"][int(row["iteration"][1:]) - 1]["training"][
                        "peak_gpu_bytes"
                    ],
                }
            )
    return result


def _variant_summary(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[tuple[str, int, int, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[(row["family"], row["d"], row["k"], row["variant"])].append(row)
    result = []
    for (family, d, k, variant), samples in sorted(groups.items()):
        improvements = np.asarray(
            [sample["final_improvement_from_a0"] for sample in samples]
        )
        result.append(
            {
                "family": family,
                "d": d,
                "k": k,
                "variant": variant,
                "seeds": [sample["training_seed"] for sample in samples],
                "seed_count": len(samples),
                "mean_final_improvement": float(np.mean(improvements)),
                "sd_final_improvement_across_seeds": float(np.std(improvements, ddof=1))
                if len(samples) > 1
                else 0.0,
                "reliable_seed_count": sum(
                    bool(sample["final_improvement_state_paired_significant"])
                    for sample in samples
                ),
                "mean_validation_value_rmse": float(
                    np.mean([sample["final_validation_value_rmse"] for sample in samples])
                ),
                "mean_validation_gradient_vector_rmse": None
                if samples[0]["final_validation_gradient_vector_rmse"] is None
                else float(
                    np.mean(
                        [sample["final_validation_gradient_vector_rmse"] for sample in samples]
                    )
                ),
                "mean_attributed_gpu_seconds": float(
                    np.mean([sample["attributed_total_gpu_seconds"] for sample in samples])
                ),
                "mean_improvement_per_gpu_minute": float(
                    np.mean(
                        [sample["improvement_per_attributed_gpu_minute"] for sample in samples]
                    )
                ),
            }
        )
    return result


def _plot_training(curves: list[dict[str, Any]]) -> None:
    stage_a = [
        row
        for row in curves
        if row["family"] == "joint_concave" and row["d"] == 400 and row["k"] == 5
    ]
    if not stage_a:
        return
    figure, axes = plt.subplots(2, 2, figsize=(11.2, 7.8), sharex=True)
    variant_names = sorted({row["variant"] for row in stage_a})
    colors = dict(zip(variant_names, plt.get_cmap("tab10").colors, strict=False))
    for iteration, column in (("a1", 0), ("a3", 1)):
        for variant in variant_names:
            subset = [
                row for row in stage_a
                if row["variant"] == variant and row["iteration"] == iteration
            ]
            if not subset:
                continue
            steps = sorted({row["step"] for row in subset})
            for prefix, linestyle in (("train", "--"), ("validation", "-")):
                value = [
                    np.mean(
                        [row[f"{prefix}_value_rmse"] for row in subset if row["step"] == step]
                    )
                    for step in steps
                ]
                axes[0, column].plot(
                    steps, value, linestyle, color=colors[variant]
                )
                gradient = [
                    np.mean(
                        [
                            row[f"{prefix}_gradient_vector_rmse"]
                            for row in subset
                            if row["step"] == step
                        ]
                    )
                    for step in steps
                ]
                axes[1, column].plot(
                    steps, gradient, linestyle, color=colors[variant]
                )
        axes[0, column].set_title(f"{iteration}: Value Actor")
        axes[1, column].set_title(f"{iteration}: Gradient Actor")
    for axis in axes.flat:
        axis.set_yscale("log")
        axis.grid(alpha=0.25)
        axis.set_xlabel("optimizer step")
    axes[0, 0].set_ylabel("RMSE")
    axes[1, 0].set_ylabel("vector RMSE")
    handles = [
        Line2D([0], [0], color=colors[name], lw=2, label=name)
        for name in variant_names
    ]
    handles.extend(
        (
            Line2D([0], [0], color="0.25", lw=2, ls="--", label="train"),
            Line2D([0], [0], color="0.25", lw=2, ls="-", label="validation"),
        )
    )
    axes[0, 1].legend(handles=handles, frameon=False, fontsize=7, ncol=2)
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "training_validation_loss.png", dpi=180)
    plt.close(figure)


def _plot_accuracy_compute(rows: list[dict[str, Any]]) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(10.5, 4.1))
    stage_a = [
        row for row in rows
        if row["family"] == "joint_concave" and row["d"] == 400 and row["k"] == 5
    ]
    stage_groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in stage_a:
        stage_groups[row["variant"]].append(row)
    colors = dict(zip(sorted(stage_groups), plt.get_cmap("tab10").colors, strict=False))
    for variant, samples in sorted(stage_groups.items()):
        x = np.asarray([sample["attributed_total_gpu_seconds"] / 60.0 for sample in samples])
        y = np.asarray([sample["final_improvement_from_a0"] for sample in samples])
        axes[0].errorbar(
            np.mean(x),
            np.mean(y),
            xerr=np.std(x, ddof=1) if len(x) > 1 else 0.0,
            yerr=np.std(y, ddof=1) if len(y) > 1 else 0.0,
            marker="s" if samples[0]["actor_kind"] == "terminal_aware" else "o",
            color=colors[variant],
            capsize=3,
            label=variant,
        )
        axes[0].annotate(
            variant,
            (np.mean(x), np.mean(y)),
            fontsize=8,
            xytext=(5, 4),
            textcoords="offset points",
        )
    axes[0].axhline(0.0, color="0.4", lw=1)
    axes[0].set_xlabel("attributed GPU minutes")
    axes[0].set_ylabel("paired improvement, a0 − final")
    axes[0].set_title("Accuracy–compute tradeoff")

    transfer = [row for row in rows if row["variant"] == "TA_A"]
    transfer_groups: dict[tuple[str, int], list[dict[str, Any]]] = defaultdict(list)
    for row in transfer:
        transfer_groups[(row["family"], row["d"])].append(row)
    transfer_colors = {
        ("joint_concave", 400): "#d62728",
        ("allen_cahn", 400): "#1f77b4",
        ("joint_concave", 100): "#2ca02c",
    }
    for (family, d), samples in sorted(transfer_groups.items()):
        by_k: dict[int, list[dict[str, Any]]] = defaultdict(list)
        for sample in samples:
            by_k[sample["k"]].append(sample)
        k_values = sorted(by_k)
        means = [
            np.mean([sample["final_improvement_from_a0"] for sample in by_k[k]])
            for k in k_values
        ]
        deviations = [
            np.std(
                [sample["final_improvement_from_a0"] for sample in by_k[k]], ddof=1
            )
            if len(by_k[k]) > 1
            else 0.0
            for k in k_values
        ]
        axes[1].errorbar(
            k_values,
            means,
            yerr=deviations,
            marker="o" if d == 400 else "D",
            color=transfer_colors[(family, d)],
            capsize=3,
            label=f"{family}, d={d}",
        )
    axes[1].axhline(0.0, color="0.4", lw=1)
    axes[1].set_xlabel("effective dimension k")
    axes[1].set_ylabel("paired improvement, a0 − final")
    axes[1].set_title("Transfer across effective dimension")
    axes[1].legend(frameon=False, fontsize=6)
    for axis in axes:
        axis.grid(alpha=0.25)
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "scaling_accuracy_compute.png", dpi=180)
    plt.close(figure)


def _plot_iterations(iteration_rows: list[dict[str, Any]]) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(11.2, 4.4), sharey=True)
    stage_a = [
        row for row in iteration_rows
        if row["family"] == "joint_concave" and row["d"] == 400 and row["k"] == 5
    ]
    transfer = [row for row in iteration_rows if row["variant"] == "TA_A"]

    def plot_groups(axis: Any, source: list[dict[str, Any]], group_keys: tuple[str, ...]) -> None:
        grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
        for row in source:
            grouped[tuple(row[key] for key in group_keys)].append(row)
        for key, samples in sorted(grouped.items(), key=lambda item: str(item[0])):
            by_iteration: dict[str, list[float]] = defaultdict(list)
            for sample in samples:
                by_iteration[sample["iteration"]].append(
                    sample["improvement_from_a0_mean"]
                )
            order = ("a0", "a1", "a2", "a3")
            means = [np.mean(by_iteration[name]) for name in order]
            deviations = [
                np.std(by_iteration[name], ddof=1) if len(by_iteration[name]) > 1 else 0.0
                for name in order
            ]
            label = key[0] if len(key) == 1 else ", ".join(map(str, key))
            axis.errorbar(order, means, yerr=deviations, marker="o", capsize=3, label=label)

    plot_groups(axes[0], stage_a, ("variant",))
    plot_groups(axes[1], transfer, ("family", "d", "k"))
    axes[0].set_title("Stage A: controlled ablation")
    axes[1].set_title("Terminal-aware transfer")
    for axis in axes:
        axis.axhline(0.0, color="0.4", lw=1)
        axis.set_xlabel("Actor iteration")
        axis.grid(alpha=0.25)
        axis.legend(frameon=False, fontsize=7)
    axes[0].set_ylabel("paired improvement from a0")
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "iteration_improvement.png", dpi=180)
    plt.close(figure)


def main() -> None:
    ensure_dirs()
    payloads = _payloads()
    iteration_rows = [row for payload in payloads for row in payload["iteration_rows"]]
    curves = _flatten_curves(payloads)
    ablation = [_ablation_row(payload) for payload in payloads]
    compute = _compute_rows(payloads)
    summary = _variant_summary(ablation)
    write_csv(RESULTS_DIR / "reference_free_improvement.csv", iteration_rows)
    write_csv(RESULTS_DIR / "training_curves.csv", curves)
    write_csv(RESULTS_DIR / "ablation_results.csv", ablation)
    write_csv(RESULTS_DIR / "compute_accounting.csv", compute)
    _plot_training(curves)
    _plot_accuracy_compute(ablation)
    _plot_iterations(iteration_rows)
    atomic_json(
        RESULTS_DIR / "analysis_summary.json",
        {
            "completed_runs": len(payloads),
            "variant_summary": summary,
            "selection_constraints": {
                "exact_test_used": False,
                "ground_truth_available_for_k_gt_2": False,
                "checkpoint_metric": "held-out noisy tangent-label loss",
            },
        },
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

