"""Publication figures from measured checkpoints only (no interpolation/smoothing)."""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .protocol import FIGURE_DIR, PDES, RESULTS_DIR


METHOD_LABELS = {
    "tangent_actor": "Tangent-Actor",
    "dpi": "DPI",
    "pinn_hte": "PINN-HTE",
}
METHOD_COLORS = {
    "tangent_actor": "#1677B8",
    "dpi": "#D55E00",
    "pinn_hte": "#2A9D55",
}
CONFIG_STYLES = {"C0": "-", "C1": "--", "C2": ":"}
PDE_LABELS = {
    "b6_fullrank": "B6 Full-Rank HJB (100D)",
    "dense_fourier": "Dense Fourier (100D)",
}


def _read() -> list[dict]:
    with (RESULTS_DIR / "checkpoint_metrics.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        for key in (
            "main_axis_seconds", "total_compute_seconds", "rmse", "rmse_sd",
            "training_seconds", "selection_seconds", "inference_seconds",
        ):
            row[key] = float(row[key]) if row.get(key) not in (None, "") else 0.0
        row["is_original_budget"] = str(row.get("is_original_budget", "")).lower() == "true"
    return rows


def _x_limit(rows: list[dict], *, total_cost: bool) -> float:
    key = "total_compute_seconds" if total_cost else "main_axis_seconds"
    largest = max((row[key] / 60.0 for row in rows), default=30.0)
    # Preserve the requested approximately 0--30 minute view while never
    # clipping a genuinely measured checkpoint that crossed a budget at an
    # algorithmic boundary.
    return max(30.0, 5.0 * math.ceil(largest / 5.0)) + 0.5


def _panel(
    ax, rows: list[dict], pde: str, *, total_cost: bool = False,
    legend: bool = False, x_max: float | None = None,
    restrict_to_x_max: bool = False,
) -> None:
    subset = [row for row in rows if row["pde"] == pde]
    x_key = "total_compute_seconds" if total_cost else "main_axis_seconds"
    for method in METHOD_LABELS:
        for configuration in CONFIG_STYLES:
            selected = sorted(
                [
                    row for row in subset
                    if row["method"] == method and row["configuration_id"] == configuration
                ],
                key=lambda row: row[x_key],
            )
            if restrict_to_x_max and x_max is not None:
                # The main paper panel is the preregistered 0--30 minute
                # window.  Do not draw a segment toward an algorithmic-boundary
                # checkpoint beyond 30 minutes, since that would visually
                # interpolate an unobserved value at the boundary.
                selected = [row for row in selected if row[x_key] / 60.0 <= x_max]
            if not selected:
                continue
            label = METHOD_LABELS[method] if configuration == "C0" else f"{METHOD_LABELS[method]} {configuration}"
            ax.plot(
                [row[x_key] / 60.0 for row in selected],
                [row["rmse"] for row in selected],
                linestyle=CONFIG_STYLES[configuration],
                color=METHOD_COLORS[method],
                linewidth=2.35 if configuration == "C0" else 1.25,
                marker="o", markersize=4.2 if configuration == "C0" else 3.2,
                alpha=1.0 if configuration == "C0" else 0.67,
                label=label,
            )
    ax.set_yscale("log")
    ax.set_xlim(0.0, x_max if x_max is not None else _x_limit(subset, total_cost=total_cost))
    ax.set_title(PDE_LABELS[pde])
    ax.set_xlabel("Cumulative training + selection (min)" if not total_cost else "Total compute incl. inference (min)")
    ax.set_ylabel("Whole-function RMSE")
    ax.grid(True, which="both", alpha=0.22, linewidth=0.6)
    if legend:
        ax.legend(frameon=False, fontsize=8, ncol=2)


def _single(rows: list[dict], pde: str, filename: str) -> None:
    fig, ax = plt.subplots(figsize=(5.15, 3.55), constrained_layout=True)
    _panel(ax, rows, pde, legend=True, x_max=30.0, restrict_to_x_max=True)
    fig.savefig(FIGURE_DIR / f"{filename}.pdf", bbox_inches="tight")
    plt.close(fig)


def _joint(rows: list[dict]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.15))
    x_max = 30.0
    for ax, pde in zip(axes, PDES):
        _panel(ax, rows, pde, x_max=x_max, restrict_to_x_max=True)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels,
               loc="lower center", bbox_to_anchor=(0.5, 0.015),
               ncol=5, frameon=False, fontsize=8)
    fig.subplots_adjust(left=0.075, right=0.99, top=0.90, bottom=0.26, wspace=0.20)
    fig.savefig(FIGURE_DIR / "two_pde_accuracy_compute.pdf", bbox_inches="tight", pad_inches=0.04)
    fig.savefig(FIGURE_DIR / "two_pde_accuracy_compute.png", dpi=260, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


def _configuration_ablation(rows: list[dict]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.15))
    x_max = _x_limit(rows, total_cost=False)
    for ax, pde in zip(axes, PDES):
        _panel(ax, rows, pde, legend=False, x_max=x_max)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.015),
               ncol=5, frameon=False, fontsize=8)
    fig.subplots_adjust(left=0.075, right=0.99, top=0.90, bottom=0.26, wspace=0.20)
    fig.savefig(FIGURE_DIR / "configuration_ablation.pdf", bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


def _total_cost(rows: list[dict]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.15))
    x_max = _x_limit(rows, total_cost=True)
    for ax, pde in zip(axes, PDES):
        _panel(ax, rows, pde, total_cost=True, x_max=x_max)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.015),
               ncol=5, frameon=False, fontsize=8)
    fig.subplots_adjust(left=0.075, right=0.99, top=0.90, bottom=0.26, wspace=0.20)
    fig.savefig(FIGURE_DIR / "accuracy_vs_total_cost.pdf", bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


def _coverage(rows: list[dict]) -> None:
    labels = [
        f"{METHOD_LABELS[method]} {configuration}"
        for method in METHOD_LABELS for configuration in CONFIG_STYLES
    ]
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.15), sharey=True, constrained_layout=True)
    for ax, pde in zip(axes, PDES):
        for y, label in enumerate(labels):
            method_label, configuration = label.rsplit(" ", 1)
            method = next(key for key, value in METHOD_LABELS.items() if value == method_label)
            times = [
                row["main_axis_seconds"] / 60.0 for row in rows
                if row["pde"] == pde and row["method"] == method
                and row["configuration_id"] == configuration
            ]
            ax.scatter(times, np.full(len(times), y), s=17, color=METHOD_COLORS[method])
        ax.set_xlim(0.0, 30.8); ax.set_title(PDE_LABELS[pde]); ax.set_xlabel("Measured checkpoint time (min)")
        ax.grid(True, axis="x", alpha=0.25)
    axes[0].set_yticks(range(len(labels)), labels, fontsize=8)
    fig.savefig(FIGURE_DIR / "checkpoint_time_coverage.pdf", bbox_inches="tight")
    plt.close(fig)


def _cost_breakdown(rows: list[dict]) -> None:
    selected = []
    for pde in PDES:
        for method in METHOD_LABELS:
            candidates = [
                row for row in rows
                if row["pde"] == pde and row["method"] == method
                and row["configuration_id"] == "C0" and row["is_original_budget"]
            ]
            if len(candidates) == 1:
                selected.append(candidates[0])
    labels = [f"{PDE_LABELS[row['pde']].split(' (')[0]}\n{METHOD_LABELS[row['method']]}" for row in selected]
    training = np.asarray([row["training_seconds"] / 60 for row in selected])
    selection = np.asarray([row["selection_seconds"] / 60 for row in selected])
    inference = np.asarray([row["inference_seconds"] / 60 for row in selected])
    x = np.arange(len(selected))
    fig, ax = plt.subplots(figsize=(8.2, 3.65), constrained_layout=True)
    ax.bar(x, training, label="Training", color="#6BAED6")
    ax.bar(x, selection, bottom=training, label="Selection", color="#FDAE6B")
    ax.bar(x, inference, bottom=training + selection, label="Inference", color="#74C476")
    ax.set_xticks(x, labels, fontsize=8); ax.set_ylabel("Measured minutes")
    ax.legend(frameon=False, ncol=3); ax.grid(True, axis="y", alpha=0.22)
    fig.savefig(FIGURE_DIR / "cost_breakdown.pdf", bbox_inches="tight")
    plt.close(fig)


def plot() -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    rows = _read()
    _single(rows, "b6_fullrank", "b6_100d_accuracy_compute")
    _single(rows, "dense_fourier", "dense_fourier_100d_accuracy_compute")
    _joint(rows)
    _configuration_ablation(rows)
    _total_cost(rows)
    _coverage(rows)
    _cost_breakdown(rows)


if __name__ == "__main__":
    plot()
