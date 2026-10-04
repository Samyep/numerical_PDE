"""Evaluate the authors' official data-free L-NN2 PINN checkpoints.

The PINN is evaluated only on the two native tasks for which the authors ship
weights (Sod and Lax).  Its continuous predictions are integrated over each
finite-volume cell with Gauss--Legendre quadrature before comparison.  This is
a native-task comparison, not part of the amortized OOD ranking.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch
from torch import nn


HERE = Path(__file__).resolve().parent
EXPERIMENTS = HERE.parent
PAPER_DIR = EXPERIMENTS / "paper_64_benchmark"
NONPERIODIC_DIR = EXPERIMENTS / "euler_1d_nonperiodic_audit"
for module_dir in (PAPER_DIR, NONPERIODIC_DIR):
    if str(module_dir) not in sys.path:
        sys.path.insert(0, str(module_dir))

import run_64_benchmark as paper  # noqa: E402
import run_nonperiodic_audit as nonperiodic  # noqa: E402


CELLS = 64
SNAPSHOTS = 64
DT = paper.euler.base.DT_SNAPSHOT
OFFICIAL_COMMIT = "cebd5f0062903ac971ff7e18063ee546352d7127"
CASES: dict[str, dict[str, Any]] = {
    "sod": {
        "display": "Sod",
        "left": (1.0, 0.0, 1.0),
        "right": (0.125, 0.0, 0.1),
        "subdirectory": "SST_Local",
        "width": 192,
    },
    "lax_published": {
        "display": "Lax (authors' state)",
        "left": (0.445, 0.689, 3.528),
        "right": (0.5, 0.0, 0.571),
        "subdirectory": "LST_Local",
        "width": 96,
    },
}
LABELS = {
    "hllc64": "HLLC-64",
    "hcfl64": "HCFL-64",
    "official_lnn2": "official data-free L-NN2 PINN",
}


class OfficialFNN(nn.Module):
    """Parameter-compatible implementation of the authors' FNN."""

    def __init__(self, width: int) -> None:
        super().__init__()
        sizes = [2] + 5 * [width] + [3]
        self.register_buffer("mean", torch.zeros(2))
        self.register_buffer("std", torch.ones(2))
        self.hidden_weights = nn.ParameterList(
            [nn.Parameter(torch.empty(sizes[i], sizes[i + 1])) for i in range(6)]
        )
        self.hidden_biases = nn.ParameterList(
            [nn.Parameter(torch.empty(1, sizes[i + 1])) for i in range(6)]
        )

    def forward(self, coordinates: torch.Tensor) -> torch.Tensor:
        values = (coordinates - self.mean) / self.std
        for index, (weight, bias) in enumerate(
            zip(self.hidden_weights, self.hidden_biases)
        ):
            values = values @ weight + bias
            if index < len(self.hidden_weights) - 1:
                values = torch.tanh(values)
        return torch.stack(
            (torch.exp(values[:, 0]), torch.exp(values[:, 1]), values[:, 2]),
            dim=-1,
        )


def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, np.generic):
        return json_ready(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def primitive_to_conservative(values: torch.Tensor) -> torch.Tensor:
    density, pressure, velocity = values.unbind(dim=-1)
    energy = pressure / (paper.euler.base.GAMMA - 1.0) + 0.5 * density * velocity**2
    return torch.stack((density, density * velocity, energy), dim=-1)


def conservative_to_primitive(values: torch.Tensor) -> torch.Tensor:
    density = values[..., 0].clamp_min(1.0e-12)
    velocity = values[..., 1] / density
    pressure = paper.euler.base.t_pressure(values)
    return torch.stack((density, velocity, pressure), dim=-1)


@torch.no_grad()
def pinn_cell_averages(
    model: OfficialFNN,
    device: torch.device,
    quadrature_order: int,
) -> tuple[torch.Tensor, torch.Tensor]:
    nodes, weights = np.polynomial.legendre.leggauss(quadrature_order)
    centers = (np.arange(CELLS) + 0.5) / CELLS
    x = centers[:, None] + nodes[None, :] / (2.0 * CELLS)
    times = np.arange(SNAPSHOTS) * DT
    x_grid = np.broadcast_to(x[None, :, :], (SNAPSHOTS, CELLS, quadrature_order))
    t_grid = np.broadcast_to(
        times[:, None, None], (SNAPSHOTS, CELLS, quadrature_order)
    )
    coordinates = torch.tensor(
        np.stack((x_grid, t_grid), axis=-1).reshape(-1, 2),
        dtype=torch.float32,
        device=device,
    )
    primitive = model(coordinates).reshape(SNAPSHOTS, CELLS, quadrature_order, 3)
    conservative = primitive_to_conservative(primitive)
    quadrature_weights = torch.tensor(
        weights / 2.0, dtype=conservative.dtype, device=device
    )
    averages = torch.einsum("scqv,q->scv", conservative, quadrature_weights)

    boundary_coordinates = torch.tensor(
        [[x_value, time] for time in times for x_value in (0.0, 1.0)],
        dtype=torch.float32,
        device=device,
    )
    boundary_primitive = model(boundary_coordinates).reshape(SNAPSHOTS, 2, 3)
    boundary = primitive_to_conservative(boundary_primitive)
    return averages[None].cpu(), boundary.cpu()


def numerical_boundary_states(trajectory: torch.Tensor) -> torch.Tensor:
    return torch.stack((trajectory[0, :, 0], trajectory[0, :, -1]), dim=1)


def entropy_flux(state: torch.Tensor) -> torch.Tensor:
    density = state[..., 0].clamp_min(1.0e-12)
    velocity = state[..., 1] / density
    return velocity * paper.euler.base.entropy(state.double())


def balance_metrics(
    trajectory: torch.Tensor, boundary: torch.Tensor
) -> dict[str, float]:
    dx = 1.0 / CELLS
    totals = trajectory[0].double().sum(dim=-2) * dx
    boundary_flux = paper.euler.base.t_flux(boundary.double())
    flux_difference = boundary_flux[:, 1] - boundary_flux[:, 0]
    cumulative = torch.zeros_like(totals)
    increments = 0.5 * DT * (flux_difference[1:] + flux_difference[:-1])
    cumulative[1:] = torch.cumsum(increments, dim=0)
    residual = totals - totals[:1] + cumulative
    component_scale = trajectory[0, 0].double().abs().sum(dim=-2) * dx
    mass_scale = trajectory[0, 0, :, 0].double().abs().sum() * dx
    scale = torch.maximum(component_scale, mass_scale.expand_as(component_scale))

    entropy = paper.euler.base.entropy(trajectory[0].double()).sum(dim=-1) * dx
    boundary_entropy_flux = entropy_flux(boundary)
    entropy_difference = boundary_entropy_flux[:, 1] - boundary_entropy_flux[:, 0]
    entropy_cumulative = torch.zeros_like(entropy)
    entropy_increments = 0.5 * DT * (
        entropy_difference[1:] + entropy_difference[:-1]
    )
    entropy_cumulative[1:] = torch.cumsum(entropy_increments, dim=0)
    entropy_balance = entropy - entropy[0] + entropy_cumulative
    entropy_step_change = entropy_balance[1:] - entropy_balance[:-1]
    return {
        "maximum_relative_boundary_aware_conservation_residual": float(
            (residual.abs() / scale).max()
        ),
        "maximum_boundary_aware_entropy_balance": float(entropy_balance.max()),
        "entropy_balance_step_violation_rate": float(
            (entropy_step_change > 1.0e-8).double().mean()
        ),
    }


def metrics(
    case: str,
    method: str,
    trajectory: torch.Tensor,
    boundary: torch.Tensor,
    reference: torch.Tensor,
    state_std: torch.Tensor,
    seconds: float,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    finite = bool(torch.isfinite(trajectory).all())
    safe = torch.where(torch.isfinite(trajectory), trajectory, torch.zeros_like(trajectory))
    density = safe[..., 0]
    pressure = paper.euler.base.t_pressure(safe)
    completed = (
        finite
        and bool((density >= paper.euler.shared.RHO_FLOOR).all())
        and bool((pressure >= paper.euler.shared.PRESSURE_FLOOR).all())
    )
    error = (safe - reference) / state_std.reshape(1, 1, 1, 3)
    error_value = float(torch.sqrt(error.double().square().mean())) if completed else None
    final_error = (
        float(torch.sqrt(error[:, -1].double().square().mean())) if completed else None
    )
    balances = (
        balance_metrics(safe, boundary)
        if completed
        else {
            "maximum_relative_boundary_aware_conservation_residual": None,
            "maximum_boundary_aware_entropy_balance": None,
            "entropy_balance_step_violation_rate": None,
        }
    )
    result: dict[str, Any] = {
        "case": case,
        "method": method,
        "label": LABELS[method],
        "completed": completed,
        "rollout_nrmse": error_value,
        "initial_nrmse": float(torch.sqrt(error[:, 0].double().square().mean())),
        "final_nrmse": final_error,
        "minimum_density": float(density.min()),
        "minimum_pressure": float(pressure.min()),
        "inference_or_rollout_seconds": seconds,
        **balances,
    }
    if extra:
        result.update(extra)
    return result


def configure_nonperiodic_cases() -> None:
    nonperiodic.SNAPSHOTS = SNAPSHOTS
    for name, spec in CASES.items():
        nonperiodic.CASES[name] = {
            "display": spec["display"],
            "left": spec["left"],
            "right": spec["right"],
            "cut": 0.5,
            "group": "pinn_native",
        }


def plot_profiles(
    trajectories: dict[str, dict[str, torch.Tensor]],
    output: Path,
) -> None:
    x = (np.arange(CELLS) + 0.5) / CELLS
    colors = {
        "reference": "#B8C2CC",
        "hllc64": "#0072B2",
        "hcfl64": "#D55E00",
        "official_lnn2": "#7B2CBF",
    }
    styles = {
        "reference": "-",
        "hllc64": "--",
        "hcfl64": "-",
        "official_lnn2": "-.",
    }
    figure, axes = plt.subplots(2, 3, figsize=(11.5, 5.4), constrained_layout=False)
    figure.subplots_adjust(
        left=0.07,
        right=0.985,
        bottom=0.10,
        top=0.82,
        wspace=0.15,
        hspace=0.22,
    )
    for row, (case, spec) in enumerate(CASES.items()):
        for column, variable in enumerate(("density", "velocity", "pressure")):
            axis = axes[row, column]
            for method in ("reference", "hllc64", "hcfl64", "official_lnn2"):
                primitive = conservative_to_primitive(trajectories[case][method])
                axis.plot(
                    x,
                    primitive[0, -1, :, column],
                    color=colors[method],
                    linestyle=styles[method],
                    linewidth=3.0 if method == "reference" else 1.5,
                    label="HLLC-2048 reference" if method == "reference" else LABELS[method],
                )
            if row == 0:
                axis.set_title(variable)
            if column == 0:
                axis.set_ylabel(spec["display"])
            if row == 1:
                axis.set_xlabel("x")
            axis.grid(alpha=0.17)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.92),
        ncol=4,
        fontsize=8.5,
        frameon=False,
    )
    figure.suptitle(
        f"Official data-free PINN native tasks, t={(SNAPSHOTS - 1) * DT:.4f}",
        y=0.985,
    )
    figure.savefig(output / "official_lnn2_native_profiles.png", dpi=240)
    plt.close(figure)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--official-root",
        type=Path,
        default=HERE.parents[4]
        / ".hcfl_third_party_audit"
        / "Datafree_PINN_Compressible",
    )
    parser.add_argument("--quadrature-order", type=int, default=8)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--results-dir", type=Path, default=HERE / "results" / "official_pinn")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable")
    output = args.results_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    configure_nonperiodic_cases()
    paper_results = PAPER_DIR / "results"
    _, _, state_std = paper.training_statistics("euler", 0, paper_results)
    hcfl = paper.load_hcfl("euler", 0, paper_results)

    rows: list[dict[str, Any]] = []
    trajectories: dict[str, dict[str, torch.Tensor]] = {}
    for case, spec in CASES.items():
        print(json.dumps({"stage": "reference", "case": case}), flush=True)
        fine_reference, reference_stats = nonperiodic.strict_hllc_rollout(case, 2048)
        reference = nonperiodic.restrict_reference(fine_reference, CELLS)
        hllc64, hllc_stats = nonperiodic.strict_hllc_rollout(case, CELLS)
        started = time.perf_counter()
        hcfl64, hcfl_stats = nonperiodic.learned_rollout(hcfl, case, CELLS)
        hcfl_seconds = time.perf_counter() - started

        checkpoint = (
            args.official_root
            / "1DRiemann"
            / spec["subdirectory"]
            / "results"
            / "model_lbfgs.pth"
        )
        model = OfficialFNN(int(spec["width"])).to(device)
        model.load_state_dict(
            torch.load(checkpoint, map_location=device, weights_only=True)
        )
        model.eval()
        started = time.perf_counter()
        pinn64, pinn_boundary = pinn_cell_averages(
            model, device, args.quadrature_order
        )
        pinn_seconds = time.perf_counter() - started

        methods = {
            "hllc64": (hllc64, numerical_boundary_states(hllc64), hllc_stats, 0.0),
            "hcfl64": (
                hcfl64,
                numerical_boundary_states(hcfl64),
                hcfl_stats,
                hcfl_seconds,
            ),
            "official_lnn2": (
                pinn64,
                pinn_boundary,
                {
                    "official_source_commit": OFFICIAL_COMMIT,
                    "official_checkpoint": str(checkpoint.relative_to(args.official_root)),
                    "cell_average_quadrature_order": args.quadrature_order,
                },
                pinn_seconds,
            ),
        }
        trajectories[case] = {"reference": reference}
        for method, (trajectory, boundary, extra, seconds) in methods.items():
            trajectories[case][method] = trajectory
            rows.append(
                metrics(
                    case,
                    method,
                    trajectory,
                    boundary,
                    reference,
                    state_std,
                    seconds,
                    extra,
                )
            )

    write_csv(output / "official_lnn2_native_metrics.csv", rows)
    payload = {
        "scope": "official authors' L-NN2 checkpoints on their native Sod/Lax tasks",
        "is_same_task_amortized_comparison": False,
        "official_source_commit": OFFICIAL_COMMIT,
        "cells": CELLS,
        "reference_cells": 2048,
        "final_time": (SNAPSHOTS - 1) * DT,
        "rows": rows,
    }
    (output / "official_lnn2_native_summary.json").write_text(
        json.dumps(json_ready(payload), indent=2), encoding="utf-8"
    )
    torch.save(trajectories, output / "official_lnn2_native_trajectories.pt")
    plot_profiles(trajectories, output)
    print(json.dumps(json_ready(rows), indent=2), flush=True)


if __name__ == "__main__":
    main()
