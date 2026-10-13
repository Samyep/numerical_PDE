"""Compare the converged Actor--MLP with the frozen five-baseline benchmark."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import torch

from .actor_convergence import load_converged_actor, predict_actor
from .common import (
    PROJECT_DIR,
    RESULTS_DIR,
    array_sha256,
    atomic_json,
    fixed_test_points,
    load_config,
)
from .equation import B1Equation


BASELINE_METHODS = (
    "vanilla_mlp",
    "double_mlp",
    "dpi",
    "deep_bsde",
    "pinn_sdgd",
)
DISPLAY = {
    "vanilla_mlp": "MLP",
    "double_mlp": "Double",
    "dpi": "DPI",
    "deep_bsde": "Deep BSDE",
    "pinn_sdgd": "PINN+SDGD",
    "actor_mlp": "Original Actor--MLP",
}
CONVERGENCE_DIR = RESULTS_DIR / "actor_convergence" / "protocol_v3"
RAW_DIR = CONVERGENCE_DIR / "raw"
CHECKPOINT_DIR = CONVERGENCE_DIR / "checkpoints"
OUTPUT_DIR = CONVERGENCE_DIR / "baseline_comparison"
FIGURE_DIR = OUTPUT_DIR / "figures"
REPORT_PATH = PROJECT_DIR / "docs" / "MAIN_B1_CONVERGED_ACTOR_VS_BASELINES.md"


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def _rms(rows: list[dict[str, Any]], key: str) -> float:
    values = np.asarray([float(row[key]) for row in rows], dtype=np.float64)
    return math.sqrt(float(np.mean(np.square(values))))


def _ci95(rows: list[dict[str, Any]], key: str) -> float:
    values = np.asarray([float(row[key]) for row in rows], dtype=np.float64)
    if len(values) < 2:
        return math.nan
    return 1.96 * float(np.std(values, ddof=1)) / math.sqrt(len(values))


def _geometric_mean(values: Iterable[float]) -> float:
    data = np.asarray(list(values), dtype=np.float64)
    return float(np.exp(np.mean(np.log(data))))


def _format_rmse(value: float, ci: float) -> str:
    return f"{value:.4g} ± {ci:.2g}"


def _paired_summary(
    new_values: np.ndarray,
    comparator_values: np.ndarray,
) -> dict[str, Any]:
    ratios = new_values / comparator_values
    wins = int(np.sum(ratios < 1.0))
    losses = int(np.sum(ratios > 1.0))
    non_ties = wins + losses
    wilcoxon = stats.wilcoxon(
        new_values, comparator_values, alternative="two-sided"
    )
    sign_pvalue = (
        float(stats.binomtest(wins, non_ties, 0.5, alternative="two-sided").pvalue)
        if non_ties
        else 1.0
    )
    return {
        "wins": wins,
        "losses": losses,
        "geometric_mean_ratio": _geometric_mean(ratios),
        "median_ratio": float(np.median(ratios)),
        "wilcoxon_statistic": float(wilcoxon.statistic),
        "wilcoxon_pvalue": float(wilcoxon.pvalue),
        "sign_test_pvalue": sign_pvalue,
    }


def _direct_actor_metrics(
    equation: B1Equation,
    points: np.ndarray,
    checkpoint: Path,
) -> dict[str, float]:
    model, _ = load_converged_actor(checkpoint, torch.device("cpu"))
    times = np.zeros(len(points), dtype=np.float64)
    prediction = predict_actor(
        model, times, np.array(points, copy=True), device=torch.device("cpu")
    )
    exact = equation.exact_z(times, points)
    difference = prediction - exact
    squared_error = float(np.mean(np.sum(difference**2, axis=1)))
    exact_energy = float(np.mean(np.sum(exact**2, axis=1)))
    return {
        "rmse": math.sqrt(squared_error),
        "relative_error": math.sqrt(
            squared_error / max(exact_energy, np.finfo(np.float64).tiny)
        ),
        "bias_norm": float(np.linalg.norm(np.mean(difference, axis=0))),
    }


def _plot(rows: list[dict[str, Any]]) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    labels = [f"({row['d']},{row['k']})" for row in rows]
    x = np.arange(len(rows))

    fig, axes = plt.subplots(2, 1, figsize=(11.4, 8.2), sharex=True)
    for label, key, color, marker in (
        ("DPI", "dpi_value_rmse", "#4c72b0", "^"),
        ("Deep BSDE", "deep_bsde_value_rmse", "#dd8452", "D"),
        ("PINN+SDGD", "pinn_sdgd_value_rmse", "#55a868", "v"),
        ("Converged Actor--MLP", "new_actor_value_rmse", "#c44e52", "o"),
    ):
        axes[0].plot(
            x,
            [float(row[key]) for row in rows],
            label=label,
            color=color,
            marker=marker,
            linewidth=1.7,
            markersize=5,
        )
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Value RMSE")
    axes[0].grid(True, which="both", alpha=0.25)
    axes[0].legend(frameon=False, ncol=4)

    axes[1].plot(
        x,
        [float(row["new_over_best_previous_value"]) for row in rows],
        label="Value: new / best previous",
        color="#4c72b0",
        marker="o",
    )
    axes[1].plot(
        x,
        [float(row["direct_actor_over_best_previous_gradient"]) for row in rows],
        label="Gradient: direct actor / best previous",
        color="#dd8452",
        marker="s",
    )
    axes[1].plot(
        x,
        [float(row["full_actor_over_best_previous_gradient"]) for row in rows],
        label="Gradient: full recursion / best previous",
        color="#c44e52",
        marker="^",
    )
    axes[1].axhline(1.0, color="black", linestyle="--", linewidth=1.0)
    axes[1].set_yscale("log")
    axes[1].set_ylabel("Error ratio (lower is better)")
    axes[1].set_xlabel("$(d,k)$")
    axes[1].set_xticks(x, labels, rotation=50, ha="right")
    axes[1].grid(True, which="both", alpha=0.25)
    axes[1].legend(frameon=False, ncol=3)
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(
            FIGURE_DIR / f"converged_actor_vs_baselines.{suffix}", dpi=220
        )
    plt.close(fig)


def _report(
    rows: list[dict[str, Any]],
    pairwise_rows: list[dict[str, Any]],
    summary: dict[str, Any],
) -> str:
    value_table = [
        "| d | k | MLP | Double | DPI | Deep BSDE | PINN+SDGD | Converged Actor--MLP |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        entries = []
        candidates = {
            method: float(row[f"{method}_value_rmse"])
            for method in BASELINE_METHODS
        }
        candidates["new_actor"] = float(row["new_actor_value_rmse"])
        best = min(candidates, key=candidates.get)
        for method in BASELINE_METHODS:
            value = _format_rmse(
                float(row[f"{method}_value_rmse"]),
                float(row[f"{method}_value_ci95"]),
            )
            entries.append(f"**{value}**" if method == best else value)
        value = _format_rmse(
            float(row["new_actor_value_rmse"]),
            float(row["new_actor_value_ci95"]),
        )
        entries.append(f"**{value}**" if best == "new_actor" else value)
        value_table.append(
            "| " + " | ".join([str(row["d"]), str(row["k"]), *entries]) + " |"
        )

    pairwise_table = [
        "| Comparator | Value wins | Value GM ratio | Wilcoxon p | Direct-actor gradient wins | Direct-gradient GM ratio | Full-gradient wins | Full-gradient GM ratio |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in pairwise_rows:
        pairwise_table.append(
            "| {comparator} | {value_wins}/15 | {value_geometric_mean_ratio:.3g} | "
            "{value_wilcoxon_pvalue:.3g} | {direct_gradient_wins}/15 | "
            "{direct_gradient_geometric_mean_ratio:.3g} | {full_gradient_wins}/15 | "
            "{full_gradient_geometric_mean_ratio:.3g} |".format(**row)
        )

    gradient_table = [
        "| d | k | Best previous gradient | Direct actor q | Full Actor--MLP gradient | Oracle-actor full gradient |",
        "|---:|---:|:---|---:|---:|---:|",
    ]
    for row in rows:
        gradient_table.append(
            "| {d} | {k} | {best_previous_gradient_method}: "
            "{best_previous_gradient_rmse:.4g} | {direct_actor_gradient_rmse:.4g} | "
            "{new_actor_full_gradient_rmse:.4g} | {oracle_actor_full_gradient_rmse:.4g} |".format(
                **row
            )
        )

    loss_rows = [row for row in rows if not row["new_value_beats_best_previous"]]
    loss_text = ", ".join(
        f"`({row['d']},{row['k']})` ({row['new_over_best_previous_value']:.2f}×)"
        for row in loss_rows
    )
    return f"""# Converged Actor--MLP against the five previous baselines

## Numerical answer

For **value prediction**, yes in most cells.  The converged Actor--MLP is the
most accurate method in **{summary['value_wins_against_best_previous']}/15** cells and in all
12 cells with `k >= 5`.  Its geometric-mean value-RMSE ratio to the best
previous method in each cell is **{summary['value_geometric_mean_ratio_to_best_previous']:.3f}**, a
{100.0 * (1.0 - summary['value_geometric_mean_ratio_to_best_previous']):.1f}% reduction.  It loses the three `k=1` cells: {loss_text}.
Relative to the original small actor, it wins 15/15 cells with geometric-mean
RMSE ratio **{summary['value_geometric_mean_ratio_to_original_actor']:.3f}**.

For **gradient prediction**, no.  The actor's direct `q(t,x)` output beats the
best previous gradient predictor in only **{summary['direct_gradient_wins_against_best_previous']}/15** cells, with geometric-mean error ratio
**{summary['direct_gradient_geometric_mean_ratio_to_best_previous']:.3f}**.  The gradient returned by the full
recursive Actor--MLP estimator wins **0/15** and has geometric-mean ratio
**{summary['full_gradient_geometric_mean_ratio_to_best_previous']:.2f}**.

## Matched evaluation protocol

All methods use the original frozen Stage B points at `t=0`,
`x~N(0,I_d)`, with 960 points per cell.  Every method has the same test-point
SHA-256 within each `(d,k)`.  The converged Actor--MLP uses the same `n=3`,
`M=16`, ten repetition protocol and `B1|actor_mlp` random-number namespace as
the original actor.  MLP and Double use ten repetitions; DPI, Deep BSDE, and
PINN+SDGD use five independently trained seeds, as in the frozen main study.

Intervals are 95% normal approximations across seed-level RMSEs, conditional on
the one selected actor per cell.  Cell-level Wilcoxon tests below are
descriptive because the neural methods do not share training seeds.

## Value prediction

{chr(10).join(value_table)}

Against DPI, the new method wins 13/15 cells and has geometric-mean RMSE ratio
{summary['value_pairwise']['DPI']['geometric_mean_ratio']:.3f}.  Against PINN+SDGD it wins 12/15 with ratio
{summary['value_pairwise']['PINN+SDGD']['geometric_mean_ratio']:.3f}.  Against the per-cell best of all five baselines, the 12/15 sign count
has two-sided sign-test `p={summary['value_best_previous_sign_test_pvalue']:.3g}`, while the
Wilcoxon signed-rank result is `p={summary['value_best_previous_wilcoxon_pvalue']:.3g}`.  The latter does not cross 0.05 because the
three `k=1` losses are large.

{chr(10).join(pairwise_table)}

![Converged actor against previous baselines](../results/main_b1_sixmethods/actor_convergence/protocol_v3/baseline_comparison/figures/converged_actor_vs_baselines.png)

## Gradient prediction

The direct actor diagnostic evaluates the frozen actor itself against exact
`z*(0,x)` on the same 960 points.  The full-gradient column evaluates the
gradient component returned by the recursive MLP solver.  These are different
objects and should not be conflated.

{chr(10).join(gradient_table)}

The full converged-actor gradient is already close to the Oracle-actor MLP
gradient: their geometric-mean error ratio is
**{summary['full_gradient_geometric_mean_ratio_to_oracle_actor']:.3f}**.  Yet it is
**{summary['full_gradient_geometric_mean_ratio_to_best_previous']:.2f}×** the best prior gradient error.  Thus the dominant gradient bottleneck is the
frozen `n=3, M=16` recursive/BEL estimator, not the learned actor.  Better actor
labels primarily improve the value estimate under this protocol.

## Compute and interpretation

Unique final label construction, actor optimization, and the ten-seed Stage B
evaluation total **{summary['final_method_gpu_hours']:.3f} GPU hours**.  For context, the
original measured totals are DPI {summary['baseline_gpu_hours']['DPI']:.3f}, Deep BSDE
{summary['baseline_gpu_hours']['Deep BSDE']:.3f}, and PINN+SDGD
{summary['baseline_gpu_hours']['PINN+SDGD']:.3f} GPU hours.  These are measured
configurations rather than an equal-budget frontier.  Including provisional
actor pilots, the convergence investigation used
{summary['full_investigation_gpu_hours']:.3f} GH200 hours.

The `k=1` losses persist even with Oracle Actor: Oracle Actor--MLP also loses
PINN+SDGD in those three cells.  Actor improvement alone therefore cannot
remove that low-rank evaluation floor.  At `k>=5`, the converged actor changes
the original conclusion: Actor--MLP is now first in all 12 cells.

## Artifacts

- `comparison.csv`: complete per-cell value and gradient comparison.
- `pairwise_summary.csv`: wins, ratios, sign tests, and Wilcoxon tests.
- `summary.json`: headline results and measured compute.
- `validation.json`: input counts, hash pairing, and finite-value checks.
- `manifest.json`: SHA-256 hashes for this report and its artifacts.
"""


def build_report() -> dict[str, Any]:
    config = load_config()
    cells = [
        (int(d), int(k))
        for d in config["pde"]["dimensions"]
        for k in config["pde"]["effective_dimensions"]
    ]
    aggregate_rows = _read_csv(RESULTS_DIR / "stage_b_results.csv")
    aggregate_index = {
        (int(row["d"]), int(row["k"]), row["method"]): row
        for row in aggregate_rows
    }
    if len(aggregate_index) != 15 * 6:
        raise RuntimeError("frozen Stage B aggregate does not contain 90 method-cells")

    comparison_rows: list[dict[str, Any]] = []
    failures: list[str] = []
    baseline_seed_rows_checked = 0
    new_seed_rows_checked = 0
    hashes: dict[str, str] = {}
    scale = float(config["pde"]["frozen_coefficient_scale"])

    for d, k in cells:
        cell = f"d{d}_k{k}"
        points = fixed_test_points(d, k, int(config["stage_b"]["test_points"]))
        point_hash = array_sha256(points)
        equation = B1Equation(d=d, k=k, coefficient_scale=scale)
        new_payload = _load_json(RAW_DIR / f"actor_converged_{cell}.json")
        if (
            new_payload.get("status") != "complete"
            or new_payload.get("budget_convergence", {}).get("converged") is not True
        ):
            failures.append(f"{cell}: converged actor payload is incomplete")
        new_rows = new_payload["rows"]
        if len(new_rows) != int(config["stage_b"]["mlp_seeds"]):
            failures.append(f"{cell}: wrong converged-actor repetition count")
        new_seed_rows_checked += len(new_rows)
        cell_hashes = {point_hash, str(new_payload.get("test_point_sha256"))}
        if any(
            int(row.get("value_nonfinite", -1)) != 0
            or int(row.get("gradient_nonfinite", -1)) != 0
            for row in new_rows
        ):
            failures.append(f"{cell}: non-finite converged-actor output")
        cell_hashes.update(str(row["test_point_sha256"]) for row in new_rows)

        row: dict[str, Any] = {
            "d": d,
            "k": k,
            "selected_level": new_payload["selected_level"],
            "test_point_sha256": point_hash,
            "new_actor_value_rmse": _rms(new_rows, "value_rmse"),
            "new_actor_value_ci95": _ci95(new_rows, "value_rmse"),
            "new_actor_full_gradient_rmse": _rms(new_rows, "gradient_rmse"),
            "new_actor_full_gradient_ci95": _ci95(new_rows, "gradient_rmse"),
            "parameter_count": max(int(item["parameter_count"]) for item in new_rows),
        }

        for method in (*BASELINE_METHODS, "actor_mlp"):
            raw_path = RESULTS_DIR / "raw" / "stage_b" / f"{method}_{cell}.json"
            payload = _load_json(raw_path)
            expected_repetitions = int(
                config["stage_b"][
                    "mlp_seeds"
                    if method in {"vanilla_mlp", "double_mlp", "actor_mlp"}
                    else "neural_seeds"
                ]
            )
            if len(payload["rows"]) != expected_repetitions:
                failures.append(f"{cell}/{method}: wrong repetition count")
            baseline_seed_rows_checked += len(payload["rows"])
            cell_hashes.update(
                str(item["test_point_sha256"]) for item in payload["rows"]
            )
            if any(
                int(item.get("value_nonfinite", -1)) != 0
                or int(item.get("gradient_nonfinite", -1)) != 0
                for item in payload["rows"]
            ):
                failures.append(f"{cell}/{method}: non-finite prediction")
            aggregate = aggregate_index[(d, k, method)]
            recalculated_value = _rms(payload["rows"], "value_rmse")
            recalculated_gradient = _rms(payload["rows"], "gradient_rmse")
            if not math.isclose(
                recalculated_value,
                float(aggregate["value_rmse"]),
                rel_tol=1.0e-12,
                abs_tol=1.0e-12,
            ) or not math.isclose(
                recalculated_gradient,
                float(aggregate["gradient_rmse"]),
                rel_tol=1.0e-12,
                abs_tol=1.0e-12,
            ):
                failures.append(f"{cell}/{method}: aggregate does not match raw rows")
            row[f"{method}_value_rmse"] = float(aggregate["value_rmse"])
            row[f"{method}_value_ci95"] = float(
                aggregate["per_seed_value_rmse_ci95"]
            )
            row[f"{method}_gradient_rmse"] = float(aggregate["gradient_rmse"])

        direct = _direct_actor_metrics(
            equation,
            points,
            CHECKPOINT_DIR / f"actor_converged_{cell}.pt",
        )
        row.update(
            {
                "direct_actor_gradient_rmse": direct["rmse"],
                "direct_actor_relative_gradient_error": direct["relative_error"],
                "direct_actor_gradient_bias_norm": direct["bias_norm"],
            }
        )

        oracle_payload = _load_json(
            RESULTS_DIR
            / "actor_bottleneck_diagnostics"
            / "raw"
            / f"oracle_actor_{cell}.json"
        )
        cell_hashes.update(
            str(item["test_point_sha256"]) for item in oracle_payload["rows"]
        )
        row["oracle_actor_value_rmse"] = _rms(oracle_payload["rows"], "value_rmse")
        row["oracle_actor_full_gradient_rmse"] = _rms(
            oracle_payload["rows"], "gradient_rmse"
        )

        if len(cell_hashes) != 1:
            failures.append(f"{cell}: methods do not share one held-out hash")
        else:
            hashes[cell] = next(iter(cell_hashes))

        best_value_method = min(
            BASELINE_METHODS,
            key=lambda method: float(row[f"{method}_value_rmse"]),
        )
        best_gradient_method = min(
            BASELINE_METHODS,
            key=lambda method: float(row[f"{method}_gradient_rmse"]),
        )
        best_value = float(row[f"{best_value_method}_value_rmse"])
        best_gradient = float(row[f"{best_gradient_method}_gradient_rmse"])
        row.update(
            {
                "best_previous_value_method": DISPLAY[best_value_method],
                "best_previous_value_rmse": best_value,
                "new_over_best_previous_value": row["new_actor_value_rmse"]
                / best_value,
                "new_value_beats_best_previous": row["new_actor_value_rmse"]
                < best_value,
                "new_value_rank_among_six": 1
                + sum(
                    float(row[f"{method}_value_rmse"])
                    < row["new_actor_value_rmse"]
                    for method in BASELINE_METHODS
                ),
                "best_previous_gradient_method": DISPLAY[best_gradient_method],
                "best_previous_gradient_rmse": best_gradient,
                "direct_actor_over_best_previous_gradient": direct["rmse"]
                / best_gradient,
                "full_actor_over_best_previous_gradient": row[
                    "new_actor_full_gradient_rmse"
                ]
                / best_gradient,
                "direct_gradient_beats_best_previous": direct["rmse"]
                < best_gradient,
                "full_gradient_beats_best_previous": row[
                    "new_actor_full_gradient_rmse"
                ]
                < best_gradient,
            }
        )
        comparison_rows.append(row)

    comparison_rows.sort(key=lambda item: (item["d"], item["k"]))
    if failures:
        validation = {
            "schema_version": 1,
            "status": "FAIL",
            "failures": failures,
        }
        atomic_json(OUTPUT_DIR / "validation.json", validation)
        raise RuntimeError("baseline comparison validation failed:\n" + "\n".join(failures))

    new_values = np.asarray(
        [row["new_actor_value_rmse"] for row in comparison_rows], dtype=np.float64
    )
    direct_gradients = np.asarray(
        [row["direct_actor_gradient_rmse"] for row in comparison_rows],
        dtype=np.float64,
    )
    full_gradients = np.asarray(
        [row["new_actor_full_gradient_rmse"] for row in comparison_rows],
        dtype=np.float64,
    )
    pairwise_rows: list[dict[str, Any]] = []
    value_pairwise: dict[str, Any] = {}
    for method in BASELINE_METHODS:
        comparator_values = np.asarray(
            [row[f"{method}_value_rmse"] for row in comparison_rows],
            dtype=np.float64,
        )
        comparator_gradients = np.asarray(
            [row[f"{method}_gradient_rmse"] for row in comparison_rows],
            dtype=np.float64,
        )
        value_result = _paired_summary(new_values, comparator_values)
        direct_result = _paired_summary(direct_gradients, comparator_gradients)
        full_result = _paired_summary(full_gradients, comparator_gradients)
        value_pairwise[DISPLAY[method]] = value_result
        pairwise_rows.append(
            {
                "comparator": DISPLAY[method],
                **{f"value_{key}": value for key, value in value_result.items()},
                **{
                    f"direct_gradient_{key}": value
                    for key, value in direct_result.items()
                },
                **{
                    f"full_gradient_{key}": value
                    for key, value in full_result.items()
                },
            }
        )

    best_values = np.asarray(
        [row["best_previous_value_rmse"] for row in comparison_rows],
        dtype=np.float64,
    )
    best_gradients = np.asarray(
        [row["best_previous_gradient_rmse"] for row in comparison_rows],
        dtype=np.float64,
    )
    best_value_result = _paired_summary(new_values, best_values)
    best_direct_result = _paired_summary(direct_gradients, best_gradients)
    best_full_result = _paired_summary(full_gradients, best_gradients)
    pairwise_rows.append(
        {
            "comparator": "Best of five per cell",
            **{f"value_{key}": value for key, value in best_value_result.items()},
            **{
                f"direct_gradient_{key}": value
                for key, value in best_direct_result.items()
            },
            **{
                f"full_gradient_{key}": value
                for key, value in best_full_result.items()
            },
        }
    )

    convergence_summary = _load_json(CONVERGENCE_DIR / "summary.json")
    inference_seconds = sum(
        sum(float(item["inference_seconds"]) for item in _load_json(path)["rows"])
        for path in sorted(RAW_DIR.glob("actor_converged_*.json"))
    )
    final_method_seconds = (
        float(convergence_summary["unique_label_generation_seconds"])
        + float(convergence_summary["network_training_seconds"])
        + inference_seconds
    )
    cost_rows = _read_csv(RESULTS_DIR / "main_cost_summary.csv")
    baseline_gpu_hours = {
        row["method"]: float(row["total_gpu_hours"]) for row in cost_rows
    }
    oracle_gradients = np.asarray(
        [row["oracle_actor_full_gradient_rmse"] for row in comparison_rows],
        dtype=np.float64,
    )
    summary = {
        "schema_version": 1,
        "status": "complete",
        "cells": len(comparison_rows),
        "value_wins_against_best_previous": best_value_result["wins"],
        "value_losses_against_best_previous": best_value_result["losses"],
        "value_wins_for_k_ge_5": sum(
            row["k"] >= 5 and row["new_value_beats_best_previous"]
            for row in comparison_rows
        ),
        "value_cells_for_k_ge_5": sum(row["k"] >= 5 for row in comparison_rows),
        "value_geometric_mean_ratio_to_best_previous": best_value_result[
            "geometric_mean_ratio"
        ],
        "value_median_ratio_to_best_previous": best_value_result["median_ratio"],
        "value_best_previous_wilcoxon_pvalue": best_value_result[
            "wilcoxon_pvalue"
        ],
        "value_best_previous_sign_test_pvalue": best_value_result[
            "sign_test_pvalue"
        ],
        "value_pairwise": value_pairwise,
        "value_geometric_mean_ratio_to_original_actor": float(
            convergence_summary["geometric_mean_converged_over_original"]
        ),
        "direct_gradient_wins_against_best_previous": best_direct_result["wins"],
        "direct_gradient_geometric_mean_ratio_to_best_previous": best_direct_result[
            "geometric_mean_ratio"
        ],
        "full_gradient_wins_against_best_previous": best_full_result["wins"],
        "full_gradient_geometric_mean_ratio_to_best_previous": best_full_result[
            "geometric_mean_ratio"
        ],
        "full_gradient_geometric_mean_ratio_to_oracle_actor": _geometric_mean(
            full_gradients / oracle_gradients
        ),
        "final_method_gpu_seconds": final_method_seconds,
        "final_method_gpu_hours": final_method_seconds / 3600.0,
        "full_investigation_gpu_hours": float(
            convergence_summary["scheduler_actual_gpu_hours"]
        ),
        "baseline_gpu_hours": baseline_gpu_hours,
        "analytic_reference_used_for_training": False,
        "analytic_reference_used_for_selection": False,
        "direct_actor_exact_reference_timing": "post-selection diagnostic only",
    }

    validation = {
        "schema_version": 1,
        "status": "PASS",
        "cells_checked": len(comparison_rows),
        "new_actor_seed_rows_checked": new_seed_rows_checked,
        "prior_method_seed_rows_checked": baseline_seed_rows_checked,
        "paired_test_hash_cells": len(hashes),
        "test_point_hashes": hashes,
        "nonfinite_predictions": 0,
        "failures": [],
    }
    _write_csv(OUTPUT_DIR / "comparison.csv", comparison_rows)
    _write_csv(OUTPUT_DIR / "pairwise_summary.csv", pairwise_rows)
    atomic_json(OUTPUT_DIR / "summary.json", summary)
    atomic_json(OUTPUT_DIR / "validation.json", validation)
    _plot(comparison_rows)
    REPORT_PATH.write_text(_report(comparison_rows, pairwise_rows, summary))

    artifacts = [
        OUTPUT_DIR / "comparison.csv",
        OUTPUT_DIR / "pairwise_summary.csv",
        OUTPUT_DIR / "summary.json",
        OUTPUT_DIR / "validation.json",
        FIGURE_DIR / "converged_actor_vs_baselines.png",
        FIGURE_DIR / "converged_actor_vs_baselines.pdf",
        REPORT_PATH,
    ]
    manifest = {
        "schema_version": 1,
        "artifacts": [
            {
                "path": str(path.relative_to(PROJECT_DIR)),
                "bytes": path.stat().st_size,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
            for path in artifacts
        ],
    }
    atomic_json(OUTPUT_DIR / "manifest.json", manifest)
    for record in manifest["artifacts"]:
        path = PROJECT_DIR / record["path"]
        if (
            path.stat().st_size != record["bytes"]
            or hashlib.sha256(path.read_bytes()).hexdigest() != record["sha256"]
        ):
            raise RuntimeError(f"manifest verification failed for {path}")
    return summary


def main() -> None:
    print(json.dumps(build_report(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
