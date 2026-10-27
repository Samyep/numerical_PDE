"""CPU-fast correctness checks for protocol and variable-horizon adaptation."""

from __future__ import annotations

import math
import json
from pathlib import Path

import numpy as np
import torch

from invariant_region_mlp.experiments.original_pde_benchmark_comparison.deep_bsde_port import (
    DeepBSDEConfig,
)

from .protocol import DISPLAY_POINTS, TRAJECTORIES, TRAJECTORY_BY_ID
from .solver import VariableStartDeepBSDE


class _TimeDriftBenchmark:
    name = "time_drift_test"
    d = 2
    T = 1.0
    sigma = 0.0

    @staticmethod
    def driver(y: torch.Tensor, z: torch.Tensor) -> torch.Tensor:
        return torch.zeros_like(y)

    @staticmethod
    def terminal(x: torch.Tensor) -> torch.Tensor:
        return x.sum(-1)

    @staticmethod
    def base_drift(t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        return t[:, None].expand_as(x)


def _tiny_config() -> DeepBSDEConfig:
    return DeepBSDEConfig(
        benchmark="time_drift_test",
        total_time=0.75,
        intervals=3,
        hidden_widths=(3, 3),
        learning_rate=1e-2,
        iterations=1,
        batch_size=4,
        validation_size=4,
        logging_frequency=1,
        y_init_range=(0.0, 1.0),
        dtype="float64",
    )


def test_nonzero_start_uses_remaining_horizon_and_absolute_drift_time() -> None:
    benchmark = _TimeDriftBenchmark()
    model = VariableStartDeepBSDE(
        benchmark, _tiny_config(), torch.tensor([1.0, -2.0]), 0.25
    ).double().eval()
    dw = torch.zeros((4, 3, 2), dtype=torch.float64)
    _, terminal_x = model.rollout(dw)
    dt = 0.25
    displacement = dt * (0.25 + 0.50 + 0.75)
    expected = torch.tensor([1.0 + displacement, -2.0 + displacement], dtype=torch.float64)
    assert model.dt == 0.25
    assert torch.allclose(terminal_x, expected.expand_as(terminal_x), atol=1e-14, rtol=0.0)


def test_invalid_full_horizon_at_nonzero_start_is_rejected() -> None:
    config = _tiny_config()
    bad = DeepBSDEConfig(**{**config.__dict__, "total_time": 1.0})
    try:
        VariableStartDeepBSDE(_TimeDriftBenchmark(), bad, torch.zeros(2), 0.25)
    except ValueError:
        return
    raise AssertionError("nonzero start accepted a full T horizon")


def test_display_reuse_is_exact_and_reference_is_unique() -> None:
    reference_points = [
        row for row in DISPLAY_POINTS if row["trajectory_id"] == "reference"
    ]
    assert {(row["sweep"], row["level"], row["checkpoint_updates"]) for row in reference_points} == {
        ("C0", 1, 1000), ("C0", 2, 2000), ("C0", 3, 4000),
        ("C1", 2, 2000), ("C2", 2, 2000), ("J", 2, 2000),
    }
    assert len({spec.signature for spec in TRAJECTORIES}) == len(TRAJECTORIES)
    assert TRAJECTORY_BY_ID["reference"].capture_updates == (1000, 2000, 4000)


def test_sigma_convention_and_driver_formulas() -> None:
    from invariant_region_mlp.experiments.main_highdim_expansion.benchmarks import make_benchmark

    y = torch.tensor([0.4, -0.7], dtype=torch.float64)
    z = torch.tensor([[1.0, 2.0] + [0.0] * 98, [3.0, 4.0] + [0.0] * 98], dtype=torch.float64)
    b6 = make_benchmark("b6_fullrank", 100)
    dense = make_benchmark("dense_fourier", 100)
    assert math.isclose(float(b6.sigma), math.sqrt(2.0), rel_tol=0.0, abs_tol=1e-15)
    assert math.isclose(float(dense.sigma), math.sqrt(2.0), rel_tol=0.0, abs_tol=1e-15)
    assert torch.allclose(b6.driver(y, z), -0.5 * z.square().sum(-1))
    assert torch.allclose(
        dense.driver(y, z), 0.5 * torch.tanh(y) * z.square().sum(-1)
    )


def test_no_truth_dependency_in_solver_module() -> None:
    import inspect
    from . import solver

    source = inspect.getsource(solver.train_query)
    assert "_truth.npz" not in source
    assert "truth_u" not in source
    assert "exact_u" not in source


def test_frozen_prediction_discovery_uses_metadata_not_legacy_folder(
    tmp_path: Path,
) -> None:
    from .aggregate import _find_prediction_json

    run = (
        tmp_path / "runs"
        / "formal_dpi_b6_fullrank_c1_integral_mc_samples-512_s0"
        / "predictions" / "dpi" / "b6_fullrank" / "C0"
    )
    run.mkdir(parents=True)
    metadata_path = run / "picard_0080.json"
    metadata_path.write_text(json.dumps({
        "pde": "b6_fullrank",
        "method": "dpi",
        "sweep": "C1",
        "parameter_name": "integral_mc_samples",
        "parameter_value": 512,
        "training_seed": 0,
    }), encoding="utf-8")
    selected = {
        "pde": "b6_fullrank",
        "method": "dpi",
        "sweep": "C1",
        "parameter_name": "integral_mc_samples",
        "parameter_value": "512",
    }
    assert _find_prediction_json(tmp_path, selected, 0) == metadata_path

