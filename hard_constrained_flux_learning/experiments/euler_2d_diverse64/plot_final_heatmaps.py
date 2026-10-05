"""Plot final-time density fields for deterministic held-out 2-D Euler cases.

The first held-out trajectory from each data family is used.  This avoids
selecting examples according to model error after looking at the test set.
All solution panels in a row share one density scale; all error panels in a
row share one absolute-error scale.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch

import run_diverse64 as R


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
DISPLAY_NAMES = {
    "oblique_riemann": "Oblique Riemann",
    "contact_shear": "Contact / shear",
    "oblique_quadrant": "Oblique quadrant",
    "radial_interface": "Radial interface",
    "colliding_waves": "Colliding waves",
    "smooth_packet": "Smooth packet",
}


def representative_indices(families: np.ndarray) -> list[int]:
    """Return the first deterministic held-out case from every family."""

    indices: list[int] = []
    for family in R.EXPECTED_FAMILIES:
        matches = np.flatnonzero(families == family)
        if matches.size == 0:
            raise RuntimeError(f"Missing held-out family: {family}")
        indices.append(int(matches[0]))
    return indices


@torch.no_grad()
def final_fields(seed: int, device: torch.device) -> dict[str, object]:
    archive = R.load_npz(R.DATA / "test.npz")
    families = archive["family_names"]
    indices = representative_indices(families)
    times = archive["times"]

    reference = torch.from_numpy(archive["q"][indices]).float()
    roe = torch.from_numpy(archive["native_roe_coarse"][indices]).float()
    initial = reference[:, 0].to(device)
    hllc, _ = R.B.rollout_hllc(initial, times)
    model = R.load_model(seed, device)
    hcfl, _ = R.B.rollout_hcfl(model, initial, times)

    return {
        "indices": indices,
        "families": [str(families[index]) for index in indices],
        "time": float(times[-1]),
        "Reference 512→64": reference[:, -1, ..., 0].numpy(),
        "Roe-64": roe[:, -1, ..., 0].numpy(),
        "HLLC-64": hllc[:, -1, ..., 0].numpy(),
        f"HCFL-64 (seed {seed})": hcfl[:, -1, ..., 0].numpy(),
    }


def configure_axis(axis: plt.Axes, row: int, rows: int, first: bool) -> None:
    axis.set_xlim(0.0, 1.0)
    axis.set_ylim(0.0, 1.0)
    axis.set_aspect("equal")
    axis.set_xticks((0.0, 0.5, 1.0) if row == rows - 1 else ())
    axis.set_yticks((0.0, 0.5, 1.0) if first else ())
    if row == rows - 1:
        axis.set_xlabel("x")
    if first:
        axis.set_ylabel("y")
    axis.tick_params(labelsize=7, length=2)


def plot_solution_fields(payload: dict[str, object], output: Path) -> None:
    families = payload["families"]
    method_names = [key for key in payload if key not in {"indices", "families", "time"}]
    rows = len(families)
    figure = plt.figure(figsize=(12.5, 17.0), constrained_layout=True)
    grid = figure.add_gridspec(rows, len(method_names) + 1, width_ratios=[1] * len(method_names) + [0.045])
    axes = np.empty((rows, len(method_names)), dtype=object)

    for row, family in enumerate(families):
        arrays = [np.asarray(payload[name])[row] for name in method_names]
        vmin = min(float(array.min()) for array in arrays)
        vmax = max(float(array.max()) for array in arrays)
        image = None
        for column, (name, array) in enumerate(zip(method_names, arrays)):
            axis = figure.add_subplot(grid[row, column])
            axes[row, column] = axis
            image = axis.imshow(
                array,
                origin="lower",
                extent=(0.0, 1.0, 0.0, 1.0),
                interpolation="nearest",
                cmap="viridis",
                vmin=vmin,
                vmax=vmax,
                rasterized=True,
            )
            configure_axis(axis, row, rows, column == 0)
            if row == 0:
                axis.set_title(name, fontsize=10)
            if column == 0:
                axis.text(
                    -0.34,
                    0.5,
                    DISPLAY_NAMES[str(family)],
                    transform=axis.transAxes,
                    ha="right",
                    va="center",
                    fontsize=9,
                    rotation=90,
                )
        color_axis = figure.add_subplot(grid[row, -1])
        colorbar = figure.colorbar(image, cax=color_axis)
        colorbar.ax.tick_params(labelsize=7, length=2)
        colorbar.set_label(r"density $\rho$", fontsize=8)

    figure.suptitle(
        f"Held-out 2-D Euler density at t = {payload['time']:.3f}\n"
        "first test case per family; nearest-cell rendering; shared scale within each row",
        fontsize=13,
    )
    figure.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(figure)


def plot_errors(payload: dict[str, object], output: Path) -> dict[str, object]:
    families = payload["families"]
    reference = np.asarray(payload["Reference 512→64"])
    method_names = [
        key
        for key in payload
        if key not in {"indices", "families", "time", "Reference 512→64"}
    ]
    rows = len(families)
    figure = plt.figure(figsize=(10.0, 17.0), constrained_layout=True)
    grid = figure.add_gridspec(rows, len(method_names) + 1, width_ratios=[1] * len(method_names) + [0.045])
    metrics: dict[str, object] = {
        "selection_rule": "first held-out case in each family",
        "time": payload["time"],
        "checkpoint_seed": 0,
        "cases": [],
    }

    for row, family in enumerate(families):
        errors = [np.abs(np.asarray(payload[name])[row] - reference[row]) for name in method_names]
        vmax = max(float(error.max()) for error in errors)
        image = None
        case_metrics: dict[str, object] = {
            "case_index": int(payload["indices"][row]),
            "family": str(family),
            "density_mae": {},
            "density_max_absolute_error": {},
        }
        for column, (name, error) in enumerate(zip(method_names, errors)):
            axis = figure.add_subplot(grid[row, column])
            image = axis.imshow(
                error,
                origin="lower",
                extent=(0.0, 1.0, 0.0, 1.0),
                interpolation="nearest",
                cmap="magma",
                vmin=0.0,
                vmax=max(vmax, np.finfo(np.float32).eps),
                rasterized=True,
            )
            configure_axis(axis, row, rows, column == 0)
            mae = float(error.mean())
            maximum = float(error.max())
            case_metrics["density_mae"][name] = mae
            case_metrics["density_max_absolute_error"][name] = maximum
            if row == 0:
                axis.set_title(name, fontsize=10)
            axis.text(
                0.02,
                0.98,
                f"MAE {mae:.3e}",
                transform=axis.transAxes,
                ha="left",
                va="top",
                fontsize=7,
                color="white",
                bbox={"facecolor": "black", "alpha": 0.55, "pad": 1.5, "edgecolor": "none"},
            )
            if column == 0:
                axis.text(
                    -0.34,
                    0.5,
                    DISPLAY_NAMES[str(family)],
                    transform=axis.transAxes,
                    ha="right",
                    va="center",
                    fontsize=9,
                    rotation=90,
                )
        color_axis = figure.add_subplot(grid[row, -1])
        colorbar = figure.colorbar(image, cax=color_axis)
        colorbar.ax.tick_params(labelsize=7, length=2)
        colorbar.set_label(r"$|\rho-\rho_{ref}|$", fontsize=8)
        metrics["cases"].append(case_metrics)

    figure.suptitle(
        f"Final-density absolute error against 512→64 reference at t = {payload['time']:.3f}\n"
        "shared error scale within each row",
        fontsize=13,
    )
    figure.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(figure)
    return metrics


def main() -> None:
    seed = 0
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print({"device": str(device), "checkpoint_seed": seed}, flush=True)
    payload = final_fields(seed, device)
    RESULTS.mkdir(parents=True, exist_ok=True)
    plot_solution_fields(payload, RESULTS / "final_time_density_heatmaps.png")
    metrics = plot_errors(payload, RESULTS / "final_time_density_error_heatmaps.png")
    (RESULTS / "final_time_heatmap_metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True), encoding="utf-8"
    )
    print(json.dumps(metrics, indent=2), flush=True)


if __name__ == "__main__":
    main()
