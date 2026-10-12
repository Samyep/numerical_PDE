"""Create final tables, figures, diagnostics, and the paper-facing report."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import shutil
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .common import FIGURES_DIR, PROJECT_DIR, RESULTS_DIR, atomic_json, load_config
from .runner import METHODS


DISPLAY = {
    "vanilla_mlp": "MLP",
    "double_mlp": "Double",
    "dpi": "DPI",
    "deep_bsde": "Deep BSDE",
    "pinn_sdgd": "PINN+SDGD",
    "actor_mlp": "Ours",
}
COLORS = {
    "vanilla_mlp": "#7f7f7f",
    "double_mlp": "#9467bd",
    "dpi": "#1f77b4",
    "deep_bsde": "#ff7f0e",
    "pinn_sdgd": "#2ca02c",
    "actor_mlp": "#d62728",
}
MARKERS = {
    "vanilla_mlp": "o",
    "double_mlp": "s",
    "dpi": "^",
    "deep_bsde": "D",
    "pinn_sdgd": "v",
    "actor_mlp": "*",
}


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str] | None = None) -> None:
    if fields is None:
        fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _number(value: str | None) -> float | None:
    if value in {None, "", "None"}:
        return None
    return float(value)


def _format_rmse(row: dict[str, str] | None) -> str:
    if row is None:
        return "N/A"
    value = float(row["value_rmse"])
    ci = _number(row.get("per_seed_value_rmse_ci95"))
    return f"{value:.4g}" if ci is None else f"{value:.4g} ± {ci:.2g}"


def _load_actor_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    raw = RESULTS_DIR / "raw" / "stage_b"
    for path in sorted(raw.glob("actor_mlp_d*_k*.json")):
        with path.open() as handle:
            payload = json.load(handle)
        actor = payload.get("shared", {}).get("actor", {})
        if not actor:
            continue
        rows.append(
            {
                "d": payload["d"],
                "k": payload["k"],
                "best_step": actor["best_step"],
                "steps_executed": actor["steps_executed"],
                "validation_noisy_label_loss": actor["best_noisy_label_validation_loss"],
                "training_label_rmse_against_exact": actor["training_label_rmse_against_exact"],
                "actor_gradient_rmse": actor["actor_validation_gradient_rmse"],
                "actor_relative_gradient_error": actor["actor_validation_relative_error"],
                "actor_bias_norm": actor["actor_validation_bias_norm"],
                "actor_norm_mean": actor["actor_validation_norm_mean"],
                "label_generation_seconds": actor["label_generation_seconds"],
                "training_seconds": actor["training_seconds"],
                "total_actor_seconds": actor["total_actor_seconds"],
                "parameter_count": actor["parameter_count"],
                "analytic_reference_used_for_training": actor["analytic_reference_used_for_training"],
            }
        )
    return rows


def _accuracy_figure(rows: list[dict[str, str]]) -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.2), sharey=True)
    for axis, d in zip(axes, (100, 200, 400)):
        for method in METHODS:
            selected = sorted(
                (row for row in rows if int(row["d"]) == d and row["method"] == method),
                key=lambda row: int(row["k"]),
            )
            if not selected:
                continue
            axis.plot(
                [int(row["k"]) for row in selected],
                [float(row["value_rmse"]) for row in selected],
                label=DISPLAY[method],
                color=COLORS[method],
                marker=MARKERS[method],
                linewidth=1.7,
                markersize=6 if method != "actor_mlp" else 9,
            )
        axis.set_title(f"d = {d}")
        axis.set_xlabel("effective dimension k")
        axis.set_xticks([1, 5, 10, 20, 40])
        axis.set_yscale("log")
        axis.grid(True, which="both", alpha=0.25)
    axes[0].set_ylabel("value RMSE")
    handles, labels = axes[-1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=6, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURES_DIR / f"accuracy_vs_effective_dimension.{suffix}", dpi=220)
    plt.close(fig)


def _compute_figure(
    accuracy: list[dict[str, str]], costs: list[dict[str, str]]
) -> None:
    accuracy_map = {
        (int(row["d"]), int(row["k"]), row["method"]): float(row["value_rmse"])
        for row in accuracy
    }
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.2), sharey=True)
    for axis, d in zip(axes, (100, 200, 400)):
        for method in METHODS:
            selected = [
                row for row in costs if int(row["d"]) == d and row["method"] == method
            ]
            x, y = [], []
            for row in selected:
                key = (d, int(row["k"]), method)
                if key in accuracy_map:
                    x.append(max(float(row["total_gpu_seconds"]), 1e-5))
                    y.append(accuracy_map[key])
            if x:
                axis.scatter(
                    x,
                    y,
                    label=DISPLAY[method],
                    color=COLORS[method],
                    marker=MARKERS[method],
                    s=55 if method != "actor_mlp" else 90,
                    alpha=0.85,
                )
        axis.set_title(f"d = {d}")
        axis.set_xlabel("measured total GPU seconds")
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.grid(True, which="both", alpha=0.25)
    axes[0].set_ylabel("value RMSE")
    handles, labels = axes[-1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=6, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURES_DIR / f"accuracy_vs_total_compute.{suffix}", dpi=220)
    plt.close(fig)


def _gradient_figure(rows: list[dict[str, str]]) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.2), sharey=True)
    for axis, d in zip(axes, (100, 200, 400)):
        for method in METHODS:
            selected = sorted(
                (
                    row
                    for row in rows
                    if int(row["d"]) == d
                    and row["method"] == method
                    and row.get("gradient_rmse") not in {None, "", "None"}
                ),
                key=lambda row: int(row["k"]),
            )
            if selected:
                axis.plot(
                    [int(row["k"]) for row in selected],
                    [float(row["gradient_rmse"]) for row in selected],
                    label=DISPLAY[method],
                    color=COLORS[method],
                    marker=MARKERS[method],
                    linewidth=1.7,
                )
        axis.set_title(f"d = {d}")
        axis.set_xlabel("effective dimension k")
        axis.set_xticks([1, 5, 10, 20, 40])
        axis.set_yscale("log")
        axis.grid(True, which="both", alpha=0.25)
    axes[0].set_ylabel("gradient-state RMSE")
    handles, labels = axes[-1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=6, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURES_DIR / f"gradient_accuracy.{suffix}", dpi=220)
    plt.close(fig)


def make_report(accounting_path: Path | None = None) -> dict[str, Any]:
    accuracy = _read_csv(RESULTS_DIR / "stage_b_results.csv")
    costs = _read_csv(RESULTS_DIR / "stage_b_costs.csv")
    config = load_config()
    dpi_profile_path = (
        PROJECT_DIR
        / "experiments"
        / "main_b1_sixmethods"
        / "dpi_final_config.json"
    )
    with dpi_profile_path.open() as handle:
        dpi_final_profile = json.load(handle)
    expected = 15 * len(METHODS)
    if len(accuracy) != expected:
        raise RuntimeError(f"need {expected} complete method-cells, found {len(accuracy)}")

    indexed = {
        (int(row["d"]), int(row["k"]), row["method"]): row for row in accuracy
    }
    wide_rows: list[dict[str, Any]] = []
    markdown_rows: list[str] = []
    ours_best = 0
    pairwise_wins = {method: 0 for method in METHODS if method != "actor_mlp"}
    ratios = {method: [] for method in METHODS if method != "actor_mlp"}
    for d in config["pde"]["dimensions"]:
        for k in config["pde"]["effective_dimensions"]:
            cell = {method: indexed.get((d, k, method)) for method in METHODS}
            numeric = {
                method: float(row["value_rmse"])
                for method, row in cell.items()
                if row is not None
            }
            best = min(numeric, key=numeric.get)
            ours_best += int(best == "actor_mlp")
            ours_value = numeric["actor_mlp"]
            for method in ratios:
                pairwise_wins[method] += int(ours_value < numeric[method])
                ratios[method].append(ours_value / numeric[method])
            wide = {"d": d, "k": k}
            formatted: list[str] = []
            for method in (
                "vanilla_mlp",
                "double_mlp",
                "dpi",
                "deep_bsde",
                "pinn_sdgd",
                "actor_mlp",
            ):
                text = _format_rmse(cell[method])
                wide[DISPLAY[method]] = text
                formatted.append(f"**{text}**" if method == best else text)
            wide_rows.append(wide)
            markdown_rows.append(
                "| " + " | ".join([str(d), str(k), *formatted]) + " |"
            )
    _write_csv(
        RESULTS_DIR / "main_accuracy_table.csv",
        wide_rows,
        ["d", "k", "MLP", "Double", "DPI", "Deep BSDE", "PINN+SDGD", "Ours"],
    )

    cost_rows: list[dict[str, Any]] = []
    for row in costs:
        cost_rows.append(
            {
                "d": row["d"],
                "k": row["k"],
                "method": DISPLAY[row["method"]],
                "training_seconds": row["training_seconds"],
                "full_evaluation_seconds": row["full_evaluation_seconds"],
                "total_gpu_seconds": row["total_gpu_seconds"],
                "parameter_count": row["parameter_count"],
                "peak_gpu_bytes": row["peak_gpu_bytes"],
                "optimizer_steps_total": row.get("optimizer_steps_total", ""),
                "training_parent_samples": row.get("training_parent_samples", ""),
                "validation_parent_samples": row.get("validation_parent_samples", ""),
                "actor_label_parent_points": row.get("actor_label_parent_points", ""),
                "terminal_samples": row.get("terminal_samples", ""),
                "transition_samples": row.get("transition_samples", ""),
                "generator_calls": row.get("generator_calls", ""),
                "gaussian_coordinate_work": row.get("gaussian_coordinate_work", ""),
                "sdgd_coordinate_second_derivatives": row.get(
                    "sdgd_coordinate_second_derivatives", ""
                ),
                "bsde_time_steps_simulated": row.get("bsde_time_steps_simulated", ""),
            }
        )
    _write_csv(RESULTS_DIR / "main_cost_table.csv", cost_rows)
    cost_summary_rows: list[dict[str, Any]] = []
    cost_markdown_rows: list[str] = []
    for method in METHODS:
        selected = [row for row in cost_rows if row["method"] == DISPLAY[method]]
        training_seconds = sum(float(row["training_seconds"]) for row in selected)
        evaluation_seconds = sum(
            float(row["full_evaluation_seconds"]) for row in selected
        )
        parameters = sorted({int(row["parameter_count"]) for row in selected})
        parameter_text = (
            f"{parameters[0]:,}"
            if len(parameters) == 1
            else f"{parameters[0]:,}–{parameters[-1]:,}"
        )
        peak_gib = max(int(row["peak_gpu_bytes"]) for row in selected) / 2**30
        summary_row = {
            "method": DISPLAY[method],
            "training_gpu_hours": training_seconds / 3600.0,
            "full_evaluation_gpu_hours": evaluation_seconds / 3600.0,
            "total_gpu_hours": (training_seconds + evaluation_seconds) / 3600.0,
            "parameter_count_range": parameter_text,
            "peak_gpu_memory_gib": peak_gib,
        }
        cost_summary_rows.append(summary_row)
        cost_markdown_rows.append(
            f"| {DISPLAY[method]} | {training_seconds / 3600.0:.4f} | "
            f"{evaluation_seconds / 3600.0:.4f} | "
            f"{(training_seconds + evaluation_seconds) / 3600.0:.4f} | "
            f"{parameter_text} | {peak_gib:.2f} |"
        )
    _write_csv(RESULTS_DIR / "main_cost_summary.csv", cost_summary_rows)
    actor_rows = _load_actor_rows()
    _write_csv(RESULTS_DIR / "actor_training.csv", actor_rows)

    method_configs = {
        "experiment_config": config,
        "method_names": DISPLAY,
        "state_convention": "z=sqrt(2)*grad(u)",
        "test_distribution": "t=0, x~N(0,I_d), 960 frozen points per cell",
        "external_source_commits": config["external_sources"],
        "final_dpi_profile": dpi_final_profile,
    }
    atomic_json(RESULTS_DIR / "method_configs.json", method_configs)

    _accuracy_figure(accuracy)
    _compute_figure(accuracy, costs)
    _gradient_figure(accuracy)

    total_seconds = sum(float(row["total_gpu_seconds"]) for row in costs)
    accounting: dict[str, Any] = {}
    if accounting_path is not None and accounting_path.is_file():
        with accounting_path.open() as handle:
            accounting = json.load(handle)
    total_actual_hours = accounting.get("total_actual_gpu_hours", total_seconds / 3600.0)

    gm_ratios = {
        method: float(np.exp(np.mean(np.log(values)))) for method, values in ratios.items()
    }
    answer = (
        f"No. Complete Actor--MLP is best in {ours_best}/15 cells; "
        if ours_best < 15
        else "Yes. Complete Actor--MLP is best in all 15 cells; "
    )
    answer += (
        f"it beats Vanilla in {pairwise_wins['vanilla_mlp']}/15 and Double in "
        f"{pairwise_wins['double_mlp']}/15. The best of the other five methods wins "
        f"{15 - ours_best}/15 cells. Total measured experiment use is {total_actual_hours:.3f} GH200 h."
    )

    table_header = [
        "| d | k | MLP | Double | DPI | Deep BSDE | PINN+SDGD | Ours |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    report = [
        "# Main B1 six-method benchmark report",
        "",
        f"**Numerical answer:** {answer}",
        "",
        "All entries below are actual held-out results at `t=0`, `x~N(0,I_d)`. "
        "MLP-family methods use `n=3, M=16` and 10 evaluation seeds. Neural baselines "
        "use five independent training seeds. Values are RMSE with a 95% normal-approximation "
        "interval across seed-level RMSEs.",
        "",
        "## Main Table 1 — accuracy",
        "",
        *table_header,
        *markdown_rows,
        "",
        "## Comparative result",
        "",
        f"Ours/Vanilla geometric-mean RMSE ratio: {gm_ratios['vanilla_mlp']:.4g}; "
        f"Ours/Double: {gm_ratios['double_mlp']:.4g}; Ours/DPI: {gm_ratios['dpi']:.4g}; "
        f"Ours/Deep BSDE: {gm_ratios['deep_bsde']:.4g}; Ours/PINN+SDGD: "
        f"{gm_ratios['pinn_sdgd']:.4g}.",
        "",
        "## Cost and protocol",
        "",
        "| Method | Training GPU h | Full-test GPU h | Total GPU h | Parameters | Peak GiB |",
        "|---|---:|---:|---:|---:|---:|",
        *cost_markdown_rows,
        "",
        f"Summed measured method-cell work is {total_seconds / 3600.0:.3f} GPU h. "
        f"Scheduler-accounted total experiment use is {total_actual_hours:.3f} GH200 h. "
        "`main_cost_table.csv` separates training and full-test inference. Actor label generation "
        "and training are included under Ours; independent recursive branches are included for Double.",
        "",
        "Deep BSDE is a global/amortized protocol trained over the same initial-state distribution; "
        "its `Z(t,x)` is learned independently. DPI uses control-variate BEL gradient labels and "
        "warm-started Picard rounds. PINN+SDGD uses two independent coordinate subsets, giving an "
        "unbiased gradient for the squared full residual, and hard-enforces the terminal condition.",
        "",
        "The final DPI pass retains the authors' HJB value-network architecture (four width-512 ELU "
        "hidden layers) and at most 5,120 optimizer steps. Under the fixed allocation it uses 8 rather "
        "than 40 Picard rounds and 32/128 rather than 4,096 training/validation Monte Carlo samples per "
        "label. This profile was frozen from the public configuration and Stage A timing, without using "
        "held-out Stage B error. The preliminary 4x128 tanh pass is archived under `superseded/` and is "
        "excluded from all final tables.",
        "",
        "The compute plot compares the measured cost and accuracy of the frozen primary configurations. "
        "The methods do not consume equal time, and no additional budget sweep was run, so this is not "
        "presented as an equal-spend frontier.",
        "",
        "## Correctness and calibration",
        "",
        "Stage 0 passed the terminal identity, analytic and autodiff PDE residual, gradient, active "
        "Hessian rank, P1 equivalence, diffusion scaling, actor identity, SDGD unbiasedness, BSDE "
        "sign, and CPU/GPU identical-draw recursion checks. Reference-only calibration froze the "
        "coefficient scale at 2.0 before solver evaluation; no cell was flagged as weakly nonlinear.",
        "",
        "Analytic values and gradients were used only for Stage 0 and held-out diagnostics/evaluation. "
        "They were not training targets, regularizers, or checkpoint-selection signals.",
        "",
        "## Artifacts",
        "",
        "- `main_accuracy_table.csv`, `main_cost_table.csv`, and `main_cost_summary.csv`",
        "- `actor_training.csv` and `pde_calibration.csv`",
        "- `figures/accuracy_vs_effective_dimension.*`",
        "- `figures/accuracy_vs_total_compute.*`",
        "- `figures/gradient_accuracy.*`",
        "",
        "## Method references",
        "",
        "- [Deep Picard Iteration paper](https://arxiv.org/abs/2409.08526) and "
        "[authors' implementation](https://github.com/DeepOptimalControl/DeepPicardIteration).",
        "- [Deep BSDE paper](https://doi.org/10.1073/pnas.1718942115) and "
        "[authors' implementation](https://github.com/frankhan91/DeepBSDE).",
        "- [SDGD paper](https://arxiv.org/abs/2307.12306) and "
        "[authors' implementation](https://github.com/zheyuanhu01/SDGD_PINN).",
        "",
        "No Stage C, component ablation, `(M,n)` sweep, or SCaSML run was performed.",
    ]
    report_path = PROJECT_DIR / "docs" / "MAIN_B1_SIXMETHOD_REPORT.md"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(report) + "\n")

    summary = {
        "ours_best_cells": ours_best,
        "pairwise_ours_wins": pairwise_wins,
        "geometric_mean_rmse_ratios_ours_over_baseline": gm_ratios,
        "summed_method_gpu_hours": total_seconds / 3600.0,
        "total_actual_gpu_hours": total_actual_hours,
        "report": str(report_path),
    }
    atomic_json(RESULTS_DIR / "final_summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--accounting", type=Path)
    args = parser.parse_args()
    print(json.dumps(make_report(args.accounting), indent=2))


if __name__ == "__main__":
    main()
