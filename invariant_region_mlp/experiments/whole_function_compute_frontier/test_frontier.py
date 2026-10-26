from __future__ import annotations

import numpy as np

from .common import plateau_reached, value_metrics
from .protocol import CONFIGURATION_GRID, GLOBAL_GPU_SECONDS_CAP, protocol_payload


def test_grid_has_exactly_eighteen_non_sdgd_trajectories() -> None:
    assert set(CONFIGURATION_GRID) == {"tangent_actor", "dpi", "pinn_hte"}
    assert all(set(rows) == {"C0", "C1", "C2"} for rows in CONFIGURATION_GRID.values())
    assert sum(len(rows) for rows in CONFIGURATION_GRID.values()) * 2 == 18
    assert GLOBAL_GPU_SECONDS_CAP == 18 * 3600
    assert all("sdgd" not in method for method in CONFIGURATION_GRID)
    assert protocol_payload()["sdgd_policy"].startswith("no SDGD-PINN job")


def test_primary_configuration_changes_are_frozen() -> None:
    assert CONFIGURATION_GRID["tangent_actor"]["C1"]["paths_per_label"] == 256
    assert CONFIGURATION_GRID["tangent_actor"]["C2"]["optimizer_steps"] == 10000
    assert CONFIGURATION_GRID["dpi"]["C1"]["epochs_per_picard"] == 32
    assert CONFIGURATION_GRID["dpi"]["C2"]["integral_mc_samples"] == 8192
    assert CONFIGURATION_GRID["pinn_hte"]["C1"]["hutchinson_vectors"] == 32
    assert CONFIGURATION_GRID["pinn_hte"]["C2"]["batch_size"] == 1024


def test_relative_mae_definition() -> None:
    truth = np.array([1.0, -2.0, 3.0])
    prediction = np.array([2.0, -1.0, 1.0])
    result = value_metrics(prediction, truth)
    assert result["relmae"] == 4.0 / 6.0


def test_plateau_rule_needs_fifteen_minutes_and_three_metrics() -> None:
    assert not plateau_reached(899.0, [1.0, 0.999, 0.998])
    assert not plateau_reached(900.0, [1.0, 0.8])
    assert plateau_reached(900.0, [1.0, 0.999, 0.998])
    assert not plateau_reached(900.0, [1.0, 0.9, 0.8])
