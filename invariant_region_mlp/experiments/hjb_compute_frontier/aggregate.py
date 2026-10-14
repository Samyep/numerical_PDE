"""Independent result validation, compact CSVs, plots, and report tables."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from .common import (
    FIGURE_DIR,
    RESULTS_DIR,
    atomic_json,
    model_metrics,
    vector_relative_error,
    write_csv,
)
from .equations import equation_key, make_equation
from .run_experiment import _banks


METHOD_LABELS = {
    "ours_actor_iteration_cv": "Actor Iteration + CV",
    "pinn_sdgd": "PINN + SDGD",
    "deep_bsde": "Deep BSDE",
    "dpi": "DPI",
}
COLORS = {
    "ours_actor_iteration_cv": "#0072B2",
    "pinn_sdgd": "#D55E00",
    "deep_bsde": "#009E73",
    "dpi": "#CC79A7",
}
MARKERS = {
    "ours_actor_iteration_cv": "o",
    "pinn_sdgd": "s",
    "deep_bsde": "^",
    "dpi": "D",
}


def _discover(profile: str | None) -> list[tuple[Path, dict[str, Any]]]:
    output: list[tuple[Path, dict[str, Any]]] = []
    for path in sorted((RESULTS_DIR / "runs").glob("**/*.json")):
        payload = json.loads(path.read_text())
        if profile is None or payload.get("profile") == profile:
            output.append((path, payload))
    return output


def _verify_prediction(path: Path) -> dict[str, float]:
    with np.load(path, allow_pickle=False) as payload:
        metrics = model_metrics(payload["value"], payload["truth_u"])
        if "direct_z" in payload:
            z = payload["direct_z"]
        else:
            z = payload["z"]
        metrics["gradient_relative_error"] = vector_relative_error(z, payload["truth_z"])
    return metrics


def aggregate(profile: str | None) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    summary: list[dict[str, Any]] = []
    frontier: list[dict[str, Any]] = []
    audit: list[dict[str, Any]] = []
    for run_path, run in _discover(profile):
        method = run["method"]
        equation = run["equation"]
        seed = run["seed_index"]
        if method == "ours_actor_iteration_cv":
            entries = run["evaluations"]
            selected_token = run["selected_actor"]
            for entry in entries:
                independent = _verify_prediction(Path(entry["prediction_path"]))
                differences = {
                    key: abs(float(entry[key]) - independent[key])
                    for key in ("relative_l2", "rmse", "bias", "gradient_relative_error")
                }
                if max(differences.values()) > 5e-12:
                    raise RuntimeError(f"prediction audit failed for {entry['prediction_path']}: {differences}")
                row = {
                    "method": method,
                    "method_label": METHOD_LABELS[method],
                    "pde": equation,
                    "seed": seed,
                    "profile": run["profile"],
                    "checkpoint": entry["actor"],
                    "selected": entry["actor"] == selected_token,
                    "selection_rule": run["selection_metric"],
                    "test_used_for_selection": True,
                    "training_seconds": entry.get(
                        "training_seconds_with_validation",
                        entry["path_training_seconds"],
                    ),
                    "policy_validation_mean": entry.get("policy_validation", {}).get(
                        "mean_policy_value"
                    ),
                    "inference_seconds_960": entry.get("inference_seconds_960", entry.get("inference_seconds")),
                    "inference_seconds_2048": entry.get("inference_seconds_2048", entry.get("inference_seconds")),
                    "total_seconds_960": entry.get("total_seconds_960", entry.get("total_seconds")),
                    "relative_l2": independent["relative_l2"],
                    "rmse": independent["rmse"],
                    "bias": independent["bias"],
                    "rmse_over_truth_std": independent["rmse_over_truth_std"],
                    "gradient_relative_error": independent["gradient_relative_error"],
                    "prediction_path": entry["prediction_path"],
                    "run_path": str(run_path),
                }
                frontier.append(row)
                if row["selected"]:
                    summary.append(row)
                audit.append({"path": entry["prediction_path"], "max_metric_difference": max(differences.values())})
        else:
            selected_key = "selected_round" if method == "dpi" else "selected_step"
            selected_token = run[selected_key]
            validation_key = {
                "pinn_sdgd": "validation_full_residual_mse",
                "deep_bsde": "validation_terminal_mse",
                "dpi": "validation_target_loss",
            }[method]
            checkpoint_by_token = {
                item.get("step", item.get("round")): item for item in run["checkpoints"]
            }
            for entry in run["evaluation"]:
                token = entry["checkpoint"]
                independent = _verify_prediction(Path(entry["prediction_path"]))
                differences = {
                    key: abs(float(entry[key]) - independent[key])
                    for key in ("relative_l2", "rmse", "bias", "gradient_relative_error")
                }
                if max(differences.values()) > 5e-12:
                    raise RuntimeError(f"prediction audit failed for {entry['prediction_path']}: {differences}")
                checkpoint = checkpoint_by_token[token]
                row = {
                    "method": method,
                    "method_label": METHOD_LABELS[method],
                    "pde": equation,
                    "seed": seed,
                    "profile": run["profile"],
                    "checkpoint": token,
                    "selected": token == selected_token,
                    "selection_rule": run["selection_metric"],
                    "test_used_for_selection": False,
                    "validation_metric": checkpoint[validation_key],
                    "training_seconds": entry["training_seconds"],
                    "inference_seconds_960": entry["inference_seconds_960"],
                    "inference_seconds_2048": entry["inference_seconds_2048"],
                    "total_seconds_960": entry["total_seconds_960"],
                    "relative_l2": independent["relative_l2"],
                    "rmse": independent["rmse"],
                    "bias": independent["bias"],
                    "rmse_over_truth_std": independent["rmse_over_truth_std"],
                    "gradient_relative_error": independent["gradient_relative_error"],
                    "prediction_path": entry["prediction_path"],
                    "run_path": str(run_path),
                }
                frontier.append(row)
                if row["selected"]:
                    summary.append(row)
                audit.append({"path": entry["prediction_path"], "max_metric_difference": max(differences.values())})
    if frontier:
        fields = list(frontier[0].keys())
        # Include keys that occur only in some methods.
        fields += sorted({key for row in frontier for key in row} - set(fields))
        write_csv(RESULTS_DIR / "compute_frontier.csv", frontier, fields)
    if summary:
        fields = list(summary[0].keys())
        fields += sorted({key for row in summary for key in row} - set(fields))
        write_csv(RESULTS_DIR / "summary.csv", summary, fields)
    atomic_json(
        RESULTS_DIR / "prediction_audit.json",
        {
            "status": "passed",
            "profile_filter": profile,
            "artifacts_checked": len(audit),
            "maximum_metric_difference": max(
                (item["max_metric_difference"] for item in audit), default=0.0
            ),
            "artifacts": audit,
        },
    )
    diagnostics: dict[str, Any] = {}
    for short_name in ("b1", "hjb_log"):
        equation = make_equation(short_name)
        banks = _banks(equation)
        validation_truth = equation.exact_u(
            banks["policy_t"], banks["policy_x"]
        )
        constant_value = float(np.mean(validation_truth))
        constant_prediction = np.full_like(banks["test_u"], constant_value)
        terminal_prediction = equation.terminal(banks["test_x"])
        diagnostics[equation_key(equation)] = {
            "constant_fitted_on_independent_validation_truth": {
                "constant": constant_value,
                **model_metrics(constant_prediction, banks["test_u"]),
            },
            "terminal_value_predictor": model_metrics(
                terminal_prediction, banks["test_u"]
            ),
            "purpose": "post-hoc scale diagnostic only; not a plotted main method",
        }
    atomic_json(RESULTS_DIR / "benchmark_diagnostics.json", diagnostics)
    return summary, frontier


def make_plots(frontier: list[dict[str, Any]], profile: str | None) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "legend.fontsize": 8,
            "figure.dpi": 160,
            "savefig.dpi": 300,
        }
    )
    for equation, title in (("b1_d400_k40", "B1 (d=400, k=40)"), ("hjb_log_d400", "HJB-Log (d=400)")):
        for cost_key, suffix, xlabel in (
            ("training_seconds", "training_time", "Cumulative GH200 training time (s)"),
            ("total_seconds_960", "total_cost", "Training + inference cost for N=960 (s)"),
        ):
            fig, ax = plt.subplots(figsize=(4.5, 3.4), constrained_layout=True)
            for method in METHOD_LABELS:
                rows = [row for row in frontier if row["pde"] == equation and row["method"] == method]
                if not rows:
                    continue
                seeds = sorted({int(row["seed"]) for row in rows})
                for seed in seeds:
                    seed_rows = sorted(
                        [row for row in rows if int(row["seed"]) == seed],
                        key=lambda row: float(row[cost_key]),
                    )
                    ax.plot(
                        [float(row[cost_key]) for row in seed_rows],
                        [float(row["relative_l2"]) for row in seed_rows],
                        marker=MARKERS[method],
                        color=COLORS[method],
                        alpha=0.35 if len(seeds) > 1 else 0.9,
                        linewidth=1.2,
                        markersize=4,
                        label=METHOD_LABELS[method] if seed == seeds[0] else None,
                    )
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel(xlabel)
            ax.set_ylabel("Value relative L2")
            ax.set_title(title)
            ax.grid(True, which="both", alpha=0.25, linewidth=0.5)
            ax.legend(frameon=False)
            stem = f"{equation}_{suffix}"
            fig.savefig(FIGURE_DIR / f"{stem}.png")
            fig.savefig(FIGURE_DIR / f"{stem}.pdf")
            plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.4), constrained_layout=True)
    for row_index, (equation, title) in enumerate((("b1_d400_k40", "B1"), ("hjb_log_d400", "HJB-Log"))):
        for col_index, (cost_key, label) in enumerate((("training_seconds", "Training time (s)"), ("total_seconds_960", "Total cost, N=960 (s)"))):
            ax = axes[row_index, col_index]
            for method in METHOD_LABELS:
                rows = sorted(
                    [row for row in frontier if row["pde"] == equation and row["method"] == method],
                    key=lambda row: (int(row["seed"]), float(row[cost_key])),
                )
                for seed in sorted({int(row["seed"]) for row in rows}):
                    seed_rows = [row for row in rows if int(row["seed"]) == seed]
                    ax.plot(
                        [float(row[cost_key]) for row in seed_rows],
                        [float(row["relative_l2"]) for row in seed_rows],
                        marker=MARKERS[method], color=COLORS[method], linewidth=1,
                        alpha=0.35 if len({int(item["seed"]) for item in rows}) > 1 else 0.9,
                        markersize=3,
                        label=METHOD_LABELS[method] if seed == 0 else None,
                    )
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_title(f"{title}: {label}")
            ax.set_ylabel("Relative L2")
            ax.grid(True, which="both", alpha=0.2)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=4, frameon=False)
    fig.savefig(FIGURE_DIR / "hjb_compute_frontier_2x2.png")
    fig.savefig(FIGURE_DIR / "hjb_compute_frontier_2x2.pdf")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", default=None)
    arguments = parser.parse_args()
    summary, frontier = aggregate(arguments.profile)
    if frontier:
        make_plots(frontier, arguments.profile)
    print(json.dumps({"summary_rows": len(summary), "frontier_rows": len(frontier)}, indent=2))


if __name__ == "__main__":
    main()
