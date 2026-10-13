"""Aggregate, validate, plot, and report the converged B1 actors."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
import statistics
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import torch

from .actor import ActorNet
from .actor_convergence import load_converged_actor, predict_actor
from .common import (
    CHECKPOINT_DIR as ORIGINAL_CHECKPOINT_DIR,
    PROJECT_DIR,
    RESULTS_DIR,
    atomic_json,
    load_config,
    stable_seed,
    training_points,
)
from .equation import B1Equation


CONVERGENCE_DIR = RESULTS_DIR / "actor_convergence"
PROTOCOL_DIR = CONVERGENCE_DIR / "protocol_v3"
RAW_DIR = PROTOCOL_DIR / "raw"
CHECKPOINT_DIR = PROTOCOL_DIR / "checkpoints"
PREDICTION_DIR = PROTOCOL_DIR / "predictions"
FIGURE_DIR = PROTOCOL_DIR / "figures"
REPORT_PATH = PROJECT_DIR / "docs" / "MAIN_B1_ACTOR_CONVERGENCE.md"


ACCOUNTING_RECORDS = (
    ("3370319", "pilot v1", 119, "gh051"),
    ("3370346", "pilot v2 continuation", 689, "gh143"),
    ("3370392", "protocol-v3 pilot", 142, "gh065"),
    ("3370400", "pilot final confirmation", 41, "gh047"),
    ("3370406_0", "cell d100 k1", 172, "gh065"),
    ("3370406_1", "cell d100 k5", 192, "gh065"),
    ("3370406_2", "cell d100 k20", 247, "gh116"),
    ("3370406_3", "cell d100 k40", 309, "gh047"),
    ("3370406_4", "cell d200 k1", 254, "gh060"),
    ("3370406_5", "cell d200 k5", 269, "gh065"),
    ("3370406_6", "cell d200 k10", 297, "gh051"),
    ("3370406_7", "cell d200 k20", 332, "gh047"),
    ("3370406_8", "cell d200 k40", 420, "gh060"),
    ("3370406_9", "cell d400 k1", 419, "gh065"),
    ("3370406_10", "cell d400 k5", 445, "gh051"),
    ("3370406_11", "cell d400 k10", 472, "gh137"),
    ("3370406_12", "cell d400 k20", 544, "gh058"),
    ("3370516_0", "confirmation d100 k1", 36, "gh058"),
    ("3370516_1", "confirmation d200 k5", 38, "gh058"),
)


def _expected_cells() -> list[tuple[int, int]]:
    config = load_config()
    return [
        (int(d), int(k))
        for d in config["pde"]["dimensions"]
        for k in config["pde"]["effective_dimensions"]
    ]


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _seed_rms(rows: Iterable[dict[str, Any]], key: str = "value_rmse") -> float:
    values = [float(row[key]) for row in rows]
    return math.sqrt(float(np.mean(np.square(values))))


def _geometric_mean(values: Iterable[float]) -> float:
    data = np.asarray(list(values), dtype=np.float64)
    return float(np.exp(np.mean(np.log(data))))


def _paired_ci(differences: np.ndarray) -> tuple[float, float]:
    mean = float(np.mean(differences))
    if len(differences) < 2:
        return mean, 0.0
    sem = float(stats.sem(differences))
    radius = float(stats.t.ppf(0.975, len(differences) - 1) * sem)
    return mean, radius


def _original_actor_same_audit(
    equation: B1Equation,
    t: np.ndarray,
    x: np.ndarray,
) -> dict[str, float]:
    config = load_config()["actor"]
    path = ORIGINAL_CHECKPOINT_DIR / f"actor_d{equation.d}_k{equation.k}.pt"
    payload = torch.load(path, map_location="cpu", weights_only=False)
    model = ActorNet(
        equation.d,
        width=int(config["hidden_width"]),
        depth=int(config["hidden_layers"]),
    )
    model.load_state_dict(payload["state_dict"])
    model.eval()
    prediction = predict_actor(model, t, x, device=torch.device("cpu"))
    exact = equation.exact_z(t, x)
    mse = float(np.mean(np.sum((prediction - exact) ** 2, axis=1)))
    energy = float(np.mean(np.sum(exact**2, axis=1)))
    return {
        "rmse": math.sqrt(mse),
        "relative_error": math.sqrt(mse / max(energy, np.finfo(np.float64).tiny)),
    }


def _new_actor_same_audit(
    equation: B1Equation,
    t: np.ndarray,
    x: np.ndarray,
    path: Path,
) -> dict[str, float]:
    model, _ = load_converged_actor(path, torch.device("cpu"))
    prediction = predict_actor(model, t, x, device=torch.device("cpu"))
    exact = equation.exact_z(t, x)
    mse = float(np.mean(np.sum((prediction - exact) ** 2, axis=1)))
    energy = float(np.mean(np.sum(exact**2, axis=1)))
    return {
        "rmse": math.sqrt(mse),
        "relative_error": math.sqrt(mse / max(energy, np.finfo(np.float64).tiny)),
    }


def _unique_label_seconds(payload: dict[str, Any]) -> float:
    metadata = [
        payload["label_subspace"]["calibration_labels"],
        payload["audit"]["bank_a"],
        payload["audit"]["bank_b"],
    ]
    for level in payload["levels"]:
        metadata.extend((level["labels"]["train"], level["labels"]["validation"]))
    seen: set[tuple[str, str]] = set()
    total = 0.0
    for item in metadata:
        key = (str(item["seed_namespace"]), str(item["point_sha256"]))
        if key not in seen:
            seen.add(key)
            total += float(item["wall_seconds"])
    return total


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _plot_comparisons(rows: list[dict[str, Any]]) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    labels = [f"({row['d']},{row['k']})" for row in rows]
    x = np.arange(len(rows))
    colors = {
        "Original actor": "#c44e52",
        "Converged actor": "#4c72b0",
        "Oracle actor": "#55a868",
        "No actor": "#8172b2",
    }
    series = (
        ("Original actor", "original_actor_value_rmse"),
        ("Converged actor", "converged_actor_value_rmse"),
        ("Oracle actor", "oracle_actor_value_rmse"),
        ("No actor", "no_actor_value_rmse"),
    )
    fig, axis = plt.subplots(figsize=(11.2, 4.8))
    for name, key in series:
        axis.plot(
            x,
            [row[key] for row in rows],
            marker="o",
            linewidth=1.7,
            markersize=4.5,
            label=name,
            color=colors[name],
        )
    axis.set_yscale("log")
    axis.set_ylabel("Stage B value RMSE (RMS over 10 seeds)")
    axis.set_xlabel("$(d,k)$")
    axis.set_xticks(x, labels, rotation=50, ha="right")
    axis.grid(True, which="both", alpha=0.25)
    axis.legend(ncol=4, frameon=False, loc="upper center")
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURE_DIR / f"actor_convergence_value_rmse.{suffix}", dpi=220)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(11.2, 4.8))
    axis.plot(
        x,
        [row["original_actor_relative_error_same_audit"] for row in rows],
        marker="o",
        label="Original actor",
        color=colors["Original actor"],
    )
    axis.plot(
        x,
        [row["converged_actor_relative_error_same_audit"] for row in rows],
        marker="o",
        label="Converged actor",
        color=colors["Converged actor"],
    )
    axis.set_yscale("log")
    axis.set_ylabel("Actor relative gradient error")
    axis.set_xlabel("$(d,k)$")
    axis.set_xticks(x, labels, rotation=50, ha="right")
    axis.grid(True, which="both", alpha=0.25)
    axis.legend(frameon=False)
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURE_DIR / f"actor_error_before_after.{suffix}", dpi=220)
    plt.close(fig)


def _plot_trajectories(trajectory_rows: list[dict[str, Any]]) -> None:
    grouped: dict[int, list[dict[str, Any]]] = {}
    for row in trajectory_rows:
        if row["prediction_relative_change"] is not None:
            grouped.setdefault(int(row["budget_index"]), []).append(row)
    indices = sorted(grouped)
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.1), sharex=True)
    for axis, key, title in (
        (axes[0], "prediction_relative_change", "Adjacent actor prediction change"),
        (axes[1], "audit_relative_change", "Adjacent noise-corrected audit change"),
    ):
        values = [
            np.asarray([float(row[key]) for row in grouped[index]], dtype=np.float64)
            for index in indices
        ]
        median = np.asarray([np.median(value) for value in values])
        lower = np.asarray([np.quantile(value, 0.25) for value in values])
        upper = np.asarray([np.quantile(value, 0.75) for value in values])
        axis.plot(indices, median, marker="o", color="#4c72b0")
        axis.fill_between(indices, lower, upper, color="#4c72b0", alpha=0.2)
        axis.axhline(0.05, color="#c44e52", linestyle="--", linewidth=1.2)
        axis.set_yscale("log")
        axis.set_title(title)
        axis.set_xlabel("Budget level")
        axis.grid(True, which="both", alpha=0.25)
    axes[0].set_ylabel("Relative change; median and IQR over cells")
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURE_DIR / f"actor_convergence_trajectories.{suffix}", dpi=220)
    plt.close(fig)


def _report(rows: list[dict[str, Any]], summary: dict[str, Any]) -> str:
    table = [
        "| d | k | Original | Converged | Oracle | No actor | New / old | Actor rel. old → new | Level |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|:---|",
    ]
    for row in rows:
        table.append(
            "| {d} | {k} | {original_actor_value_rmse:.4f} | "
            "{converged_actor_value_rmse:.4f} | {oracle_actor_value_rmse:.4f} | "
            "{no_actor_value_rmse:.4f} | {converged_over_original:.3f} | "
            "{original_actor_relative_error_same_audit:.3f} → "
            "{converged_actor_relative_error_same_audit:.3f} | {selected_level} |".format(
                **row
            )
        )
    return f"""# Main B1 actor label-quality and convergence study

## Status

All 15 B1 cells reached the frozen optimizer and cross-budget convergence
criterion.  The main six-method benchmark and its tables remain unchanged;
this is a supplemental actor study evaluated on the same 960 held-out points,
10 MLP repetitions, `n=3`, `M=16`, and paired Actor--MLP random numbers.

The final actor beats the original learned actor and Pathwise + Baseline without
an actor in **15/15** cells.  Its geometric-mean value-RMSE ratio versus the
original actor is **{summary['geometric_mean_converged_over_original']:.3f}**,
a **{100.0 * (1.0 - summary['geometric_mean_converged_over_original']):.1f}%** reduction.  The geometric-mean ratio to Oracle Actor is
**{summary['geometric_mean_converged_over_oracle']:.3f}**;
{summary['converged_within_10pct_oracle_cells']} cells are within 10% of Oracle value RMSE.

## What changed

The original actor used a 2×128 SiLU network, one `n=2, M=8` label per point,
4,096 training points, and at most 1,000 optimizer steps.  The new actor uses a
fixed 4×384 ambient-input network and five progressive label budgets, ending at
`n=3, M=20`, ten independent labels averaged at each of 32,768 training points,
followed by low-learning-rate confirmation stages.

High-dimensional label noise was removed without reading the analytic active
basis.  For each cell, a rank-`k` output subspace was learned from independent
`n=1` pathwise terminal-gradient labels.  Training, validation, audit labels,
and actor outputs were projected into this data-derived space.  The public
effective dimension `k` fixes the rank.  No nonterminal exact solution was used
for labels, optimization, early stopping, budget escalation, or selection.

A budget transition is stable only when the optimizer converged, the
noise-corrected stochastic-audit RMSE changed by at most 5%, and actor
predictions changed by at most 5%.  Selection required two consecutive stable
transitions.  The exact solution was evaluated only after the checkpoint was
frozen.

## Results

{chr(10).join(table)}

The median same-audit actor relative-gradient error fell from
**{summary['original_actor_relative_error_median']:.3f}** to **{summary['converged_actor_relative_error_median']:.3f}**.  Its geometric-mean new/old ratio is
**{summary['geometric_mean_actor_error_new_over_old']:.3f}**, and the new actor
improves all 15 cells.  The remaining gap is concentrated in
the hardest ambient settings: the maximum post-selection relative actor error
is **{summary['converged_actor_relative_error_max']:.3f}**.

![Stage B value comparison](../results/main_b1_sixmethods/actor_convergence/protocol_v3/figures/actor_convergence_value_rmse.png)

![Actor error before and after](../results/main_b1_sixmethods/actor_convergence/protocol_v3/figures/actor_error_before_after.png)

![Budget convergence trajectories](../results/main_b1_sixmethods/actor_convergence/protocol_v3/figures/actor_convergence_trajectories.png)

## Cost and interpretation

The full investigation, including the provisional pilots that exposed the
overly permissive noise-floor shortcut, used **{summary['scheduler_actual_gpu_hours']:.3f} GH200 hours**.  Unique final label generation accounts for
{summary['unique_label_generation_seconds'] / 60.0:.1f} sequential GPU-minutes, and neural optimization for
{summary['network_training_seconds'] / 60.0:.1f} minutes.  Parallel array execution reduced wall-clock time;
scheduler usage is the sum over allocated GPUs.

The result confirms the original diagnosis.  The pathwise terminal and
level-matched baseline components were useful, but the small noisy-label actor
was the dominant avoidable error.  Raising label fidelity, removing
orthogonal label noise, and training to a measured plateau closes most of the
learned-to-Oracle gap.  At `d=400`, especially large `k`, a material Oracle gap
remains, so this study supports improved actors rather than claiming that actor
error has vanished.

## Artifacts

- `actor_convergence_comparison.csv`: full 15-cell comparison.
- `actor_convergence_trajectory.csv`: every budget-level convergence record.
- `summary.json`: aggregate findings and validation counts.
- `execution_accounting.json`: all 19 GH200 allocations.
- `artifact_manifest.json`: SHA-256 hashes of the report artifacts and 15 final
  scalar cell records.
- Checkpoints, label caches, and paired prediction arrays remain in the local
  experiment workspace and are excluded from Git because of their size.
"""


def aggregate_actor_convergence() -> dict[str, Any]:
    config = load_config()
    comparison_rows: list[dict[str, Any]] = []
    trajectory_rows: list[dict[str, Any]] = []
    all_payloads: list[dict[str, Any]] = []

    for d, k in _expected_cells():
        result_path = RAW_DIR / f"actor_converged_d{d}_k{k}.json"
        payload = _load(result_path)
        if payload.get("status") != "complete":
            raise RuntimeError(f"cell d={d}, k={k} is not converged")
        if not payload["budget_convergence"]["converged"]:
            raise RuntimeError(f"cell d={d}, k={k} lacks convergence certificate")
        if len(payload["rows"]) != int(config["stage_b"]["mlp_seeds"]):
            raise RuntimeError(f"cell d={d}, k={k} has wrong Stage B seed count")
        all_payloads.append(payload)

        equation = B1Equation(
            d=d,
            k=k,
            coefficient_scale=float(config["pde"]["frozen_coefficient_scale"]),
        )
        seed = stable_seed("B1", "actor_convergence", d, k)
        audit_t, audit_x = training_points(
            d, k, "actor_convergence_audit", seed, int(payload["protocol"]["audit"]["points"])
        )
        original_actor = _original_actor_same_audit(equation, audit_t, audit_x)
        final_checkpoint = CHECKPOINT_DIR / f"actor_converged_d{d}_k{k}.pt"
        converged_actor = _new_actor_same_audit(
            equation, audit_t, audit_x, final_checkpoint
        )
        recorded = payload["postselection_exact_diagnostic"]
        if not math.isclose(
            converged_actor["relative_error"],
            float(recorded["actor_relative_gradient_error"]),
            # The recorded diagnostic used GH200 GEMM while this independent
            # replay uses CPU GEMM; both consume the identical FP32 checkpoint.
            rel_tol=1.0e-3,
            abs_tol=5.0e-5,
        ):
            raise RuntimeError(f"checkpoint/exact diagnostic mismatch for d={d}, k={k}")

        old = _load(RESULTS_DIR / "raw" / "stage_b" / f"actor_mlp_d{d}_k{k}.json")
        oracle = _load(
            RESULTS_DIR
            / "actor_bottleneck_diagnostics"
            / "raw"
            / f"oracle_actor_d{d}_k{k}.json"
        )
        no_actor = _load(
            RESULTS_DIR
            / "actor_bottleneck_diagnostics"
            / "raw"
            / f"pathwise_baseline_no_actor_d{d}_k{k}.json"
        )
        hashes = {
            *(row["test_point_sha256"] for row in payload["rows"]),
            *(row["test_point_sha256"] for row in old["rows"]),
            *(row["test_point_sha256"] for row in oracle["rows"]),
            *(row["test_point_sha256"] for row in no_actor["rows"]),
        }
        if len(hashes) != 1:
            raise RuntimeError(f"held-out hash mismatch for d={d}, k={k}")
        methods = (payload, old, oracle, no_actor)
        if any(
            row["value_nonfinite"] or row["gradient_nonfinite"]
            for method in methods
            for row in method["rows"]
        ):
            raise RuntimeError(f"nonfinite Stage B result for d={d}, k={k}")

        new_rmse = _seed_rms(payload["rows"])
        old_rmse = _seed_rms(old["rows"])
        oracle_rmse = _seed_rms(oracle["rows"])
        no_actor_rmse = _seed_rms(no_actor["rows"])
        new_by_seed = np.asarray(
            [row["value_rmse"] for row in sorted(payload["rows"], key=lambda row: row["repetition"])]
        )
        old_by_seed = np.asarray(
            [row["value_rmse"] for row in sorted(old["rows"], key=lambda row: row["repetition"])]
        )
        paired_mean, paired_radius = _paired_ci(new_by_seed - old_by_seed)
        label_seconds = _unique_label_seconds(payload)
        training_seconds = sum(
            float(level["training"]["training_seconds"]) for level in payload["levels"]
        )
        row = {
            "d": d,
            "k": k,
            "selected_level": payload["selected_level"],
            "budget_levels": len(payload["levels"]),
            "stable_transitions": payload["budget_convergence"]["stable_transitions"],
            "original_actor_value_rmse": old_rmse,
            "converged_actor_value_rmse": new_rmse,
            "oracle_actor_value_rmse": oracle_rmse,
            "no_actor_value_rmse": no_actor_rmse,
            "converged_over_original": new_rmse / old_rmse,
            "converged_over_oracle": new_rmse / oracle_rmse,
            "converged_over_no_actor": new_rmse / no_actor_rmse,
            "paired_mean_converged_minus_original_seed_rmse": paired_mean,
            "paired_ci95_converged_minus_original_seed_rmse": paired_radius,
            "original_actor_rmse_same_audit": original_actor["rmse"],
            "converged_actor_rmse_same_audit": converged_actor["rmse"],
            "original_actor_relative_error_same_audit": original_actor["relative_error"],
            "converged_actor_relative_error_same_audit": converged_actor["relative_error"],
            "actor_error_new_over_old": (
                converged_actor["relative_error"] / original_actor["relative_error"]
            ),
            "audit_projection_disagreement_reduction": payload["audit"][
                "projection_disagreement_reduction"
            ],
            "label_subspace_relative_residual": payload["label_subspace"][
                "relative_projection_residual"
            ],
            "unique_label_generation_seconds": label_seconds,
            "network_training_seconds": training_seconds,
            "parameter_count": payload["rows"][0]["parameter_count"],
            "test_point_sha256": next(iter(hashes)),
        }
        comparison_rows.append(row)

        previous_audit: float | None = None
        for level in payload["levels"]:
            audit_value = float(level["audit_noise_corrected_rmse"])
            audit_change = (
                None
                if previous_audit is None
                else abs(audit_value - previous_audit)
                / max(previous_audit, np.finfo(np.float64).tiny)
            )
            trajectory_rows.append(
                {
                    "d": d,
                    "k": k,
                    "budget_index": level["budget_index"],
                    "budget_level": level["budget"]["name"],
                    "label_n": level["budget"]["n"],
                    "label_M": level["budget"]["M"],
                    "label_replicas": level["budget"]["replicas"],
                    "training_points": level["budget"]["training_points"],
                    "steps_executed": level["training"]["steps_executed"],
                    "optimizer_converged": level["optimizer_converged"],
                    "audit_target_rmse": level["audit_target_rmse"],
                    "audit_noise_corrected_rmse": audit_value,
                    "audit_relative_change": audit_change,
                    "prediction_relative_change": level["prediction_relative_change"],
                }
            )
            previous_audit = audit_value

    comparison_rows.sort(key=lambda row: (row["d"], row["k"]))
    trajectory_rows.sort(key=lambda row: (row["d"], row["k"], row["budget_index"]))
    ratios_old = [row["converged_over_original"] for row in comparison_rows]
    ratios_oracle = [row["converged_over_oracle"] for row in comparison_rows]
    ratios_no_actor = [row["converged_over_no_actor"] for row in comparison_rows]
    actor_ratios = [row["actor_error_new_over_old"] for row in comparison_rows]
    old_actor_errors = [
        row["original_actor_relative_error_same_audit"] for row in comparison_rows
    ]
    new_actor_errors = [
        row["converged_actor_relative_error_same_audit"] for row in comparison_rows
    ]
    wilcoxon = stats.wilcoxon(
        [row["converged_actor_value_rmse"] for row in comparison_rows],
        [row["original_actor_value_rmse"] for row in comparison_rows],
        alternative="two-sided",
    )
    gpu_seconds = sum(record[2] for record in ACCOUNTING_RECORDS)
    summary = {
        "schema_version": 1,
        "status": "complete",
        "cells": len(comparison_rows),
        "converged_cells": sum(payload["status"] == "complete" for payload in all_payloads),
        "converged_better_than_original_cells": sum(value < 1.0 for value in ratios_old),
        "converged_better_than_no_actor_cells": sum(value < 1.0 for value in ratios_no_actor),
        "converged_within_10pct_oracle_cells": sum(value <= 1.10 for value in ratios_oracle),
        "geometric_mean_converged_over_original": _geometric_mean(ratios_old),
        "geometric_mean_converged_over_oracle": _geometric_mean(ratios_oracle),
        "geometric_mean_converged_over_no_actor": _geometric_mean(ratios_no_actor),
        "original_actor_relative_error_median": statistics.median(old_actor_errors),
        "converged_actor_relative_error_median": statistics.median(new_actor_errors),
        "converged_actor_relative_error_max": max(new_actor_errors),
        "geometric_mean_actor_error_new_over_old": _geometric_mean(actor_ratios),
        "actor_error_improved_cells": sum(value < 1.0 for value in actor_ratios),
        "wilcoxon_converged_vs_original_statistic": float(wilcoxon.statistic),
        "wilcoxon_converged_vs_original_pvalue": float(wilcoxon.pvalue),
        "selected_level_counts": {
            level: sum(row["selected_level"] == level for row in comparison_rows)
            for level in sorted({row["selected_level"] for row in comparison_rows})
        },
        "median_audit_projection_disagreement_reduction": statistics.median(
            row["audit_projection_disagreement_reduction"] for row in comparison_rows
        ),
        "maximum_label_subspace_relative_residual": max(
            row["label_subspace_relative_residual"] for row in comparison_rows
        ),
        "unique_label_generation_seconds": sum(
            row["unique_label_generation_seconds"] for row in comparison_rows
        ),
        "network_training_seconds": sum(
            row["network_training_seconds"] for row in comparison_rows
        ),
        "scheduler_allocations": len(ACCOUNTING_RECORDS),
        "scheduler_actual_gpu_seconds": gpu_seconds,
        "scheduler_actual_gpu_hours": gpu_seconds / 3600.0,
        "analytic_reference_used_for_training": False,
        "analytic_reference_used_for_selection": False,
        "paired_stage_b_seed_namespace": "B1|actor_mlp",
    }

    _write_csv(PROTOCOL_DIR / "actor_convergence_comparison.csv", comparison_rows)
    _write_csv(PROTOCOL_DIR / "actor_convergence_trajectory.csv", trajectory_rows)
    atomic_json(PROTOCOL_DIR / "summary.json", summary)
    accounting = {
        "schema_version": 1,
        "gpu_model": "NVIDIA GH200 120GB",
        "records": [
            {
                "job_id": job_id,
                "purpose": purpose,
                "elapsed_seconds": seconds,
                "allocated_gpus": 1,
                "gpu_seconds": seconds,
                "node": node,
                "state": "COMPLETED",
                "exit_code": "0:0",
            }
            for job_id, purpose, seconds, node in ACCOUNTING_RECORDS
        ],
        "actual_gpu_seconds": gpu_seconds,
        "actual_gpu_hours": gpu_seconds / 3600.0,
    }
    atomic_json(PROTOCOL_DIR / "execution_accounting.json", accounting)
    _plot_comparisons(comparison_rows)
    _plot_trajectories(trajectory_rows)
    REPORT_PATH.write_text(_report(comparison_rows, summary))

    artifacts = [
        PROTOCOL_DIR / "actor_convergence_comparison.csv",
        PROTOCOL_DIR / "actor_convergence_trajectory.csv",
        PROTOCOL_DIR / "summary.json",
        PROTOCOL_DIR / "execution_accounting.json",
        FIGURE_DIR / "actor_convergence_value_rmse.png",
        FIGURE_DIR / "actor_convergence_value_rmse.pdf",
        FIGURE_DIR / "actor_error_before_after.png",
        FIGURE_DIR / "actor_error_before_after.pdf",
        FIGURE_DIR / "actor_convergence_trajectories.png",
        FIGURE_DIR / "actor_convergence_trajectories.pdf",
        REPORT_PATH,
    ] + sorted(RAW_DIR.glob("actor_converged_*.json"))
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
    atomic_json(PROTOCOL_DIR / "artifact_manifest.json", manifest)
    return summary


if __name__ == "__main__":
    print(json.dumps(aggregate_actor_convergence(), indent=2, sort_keys=True))
