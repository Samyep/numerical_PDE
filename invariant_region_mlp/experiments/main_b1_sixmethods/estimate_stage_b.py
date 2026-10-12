"""Conservative Stage B wall/GPU-time estimate from measured Stage A costs."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from .common import RESULTS_DIR, atomic_json, load_config


def estimate() -> dict[str, Any]:
    config = load_config()
    with (RESULTS_DIR / "stage_a_costs.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    totals: dict[str, float] = {}
    for row in rows:
        totals[row["method"]] = totals.get(row["method"], 0.0) + float(
            row["total_gpu_seconds"]
        )

    cell_ratio = 15 / 3
    neural_seed_ratio = config["stage_b"]["neural_seeds"] / config["stage_a"]["neural_seeds"]
    # Measured only at d=100 and k={1,10,20}. This explicit blanket factor
    # conservatively covers wider d-dependent input/output layers and K=81.
    full_grid_dimension_factor = 1.75
    step_ratios = {
        "dpi": (
            config["dpi"]["stage_b_picard_rounds"]
            * config["dpi"]["stage_b_steps_per_round"]
            / (
                config["dpi"]["stage_a_picard_rounds"]
                * config["dpi"]["stage_a_steps_per_round"]
            )
        ),
        "pinn_sdgd": config["sdgd_pinn"]["stage_b_steps"]
        / config["sdgd_pinn"]["stage_a_steps"],
        "deep_bsde": config["deep_bsde"]["stage_b_steps"]
        / config["deep_bsde"]["stage_a_steps"],
    }
    method_seconds: dict[str, float] = {}
    for method in ("dpi", "pinn_sdgd", "deep_bsde"):
        method_seconds[method] = (
            totals[method]
            * cell_ratio
            * neural_seed_ratio
            * step_ratios[method]
            * full_grid_dimension_factor
        )

    # The dedicated GH200 n=3,M=16 timing benchmark is a tighter estimate for
    # all three MLP-family evaluations and actor preparation than extrapolating
    # from the tiny Stage A point set.
    timing_path = Path(
        "/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/results/estimate.json"
    )
    with timing_path.open() as handle:
        timing = json.load(handle)
    mlp_actor_seconds = float(timing["mlp_stage_b_plus_actor_preparation_seconds"])
    method_seconds["mlp_family_and_actor_preparation"] = mlp_actor_seconds
    total = sum(method_seconds.values())
    payload = {
        "source": "measured Stage A method costs plus dedicated GH200 Stage B MLP timing",
        "stage_a_measured_gpu_seconds": sum(totals.values()),
        "cell_ratio": cell_ratio,
        "neural_seed_ratio": neural_seed_ratio,
        "full_grid_dimension_and_K_safety_factor": full_grid_dimension_factor,
        "maximum_step_ratios": step_ratios,
        "method_estimated_seconds": method_seconds,
        "stage_b_estimated_gpu_seconds": total,
        "stage_b_estimated_gpu_hours": total / 3600.0,
        "recommended_wall_limit_hours": 10,
        "interpretation": "conservative upper planning estimate; early stopping may reduce actual time",
    }
    atomic_json(RESULTS_DIR / "stage_b_runtime_estimate.json", payload)
    return payload


def main() -> None:
    print(json.dumps(estimate(), indent=2))


if __name__ == "__main__":
    main()
