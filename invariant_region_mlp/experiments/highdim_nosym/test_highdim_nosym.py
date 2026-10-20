"""Fast regression tests for the new equation-only adapters."""

from __future__ import annotations

import numpy as np
import torch

from .analyze import _curtailed_rows, _json_rows
from .equations import fixed_points, make_actor_benchmark, make_b6_equation, make_s05_equation
from .protocol import fh_tasks
from .solver import B6TorchEquation, load_b6_selection


def test_frozen_grid_has_expected_size() -> None:
    # Per problem: 2 dimensions x 3 repetitions x (8 raw + 8 path + 6 repaired).
    assert len(fh_tasks()) == 2 * 2 * 3 * 22


def test_compute_curtailment_is_not_treated_as_an_estimate(tmp_path) -> None:
    result = {"pde_id": "S05", "dimension": 400, "method": "raw", "n": 2, "M": 16}
    (tmp_path / "rep00.json").write_text(__import__("json").dumps(result), encoding="utf-8")
    marker = {
        "reason": "threshold",
        "rep0_seconds": 1801.0,
        "task": {"pde_id": "S05", "d": 400, "method": "raw", "n": 2, "M": 16, "repetition": 1},
    }
    (tmp_path / "rep01.curtailed.json").write_text(
        __import__("json").dumps(marker), encoding="utf-8"
    )
    assert len(_json_rows(tmp_path)) == 1
    curtailed = _curtailed_rows(tmp_path)
    assert len(curtailed) == 1
    assert curtailed[0]["status"] == "curtailed_after_rep0"
    assert curtailed[0]["rep0_seconds"] == 1801.0


def test_b6_construction_and_terminal_identity() -> None:
    first = make_b6_equation(100, 1.5, 1.0)
    second = make_b6_equation(100, 1.5, 1.0)
    assert first is second
    assert first.K == 200
    assert np.ptp(first.n2) > 0.0
    rng = np.random.default_rng(12)
    x = rng.uniform(-1.0, 1.0, (7, 100))
    state = first.exact_state(np.full(7, first.T), x)
    assert np.array_equal(state[:, 0], first.terminal(x))
    assert np.array_equal(state[:, 1:], first.terminal_gradient_z(x))


def test_actor_equation_adapters_match_numpy() -> None:
    rng = np.random.default_rng(13)
    for pde_id, equation in (
        ("S05", make_s05_equation(100)),
        ("B6", make_b6_equation(100, 1.5, 1.0)),
    ):
        benchmark = make_actor_benchmark(
            pde_id,
            100,
            beta=1.5 if pde_id == "B6" else None,
            gamma=1.0 if pde_id == "B6" else None,
        )
        x = rng.uniform(-1.0, 1.0, (8, 100))
        tensor = torch.tensor(x, dtype=torch.float64)
        observed_u = benchmark.terminal(tensor).numpy()
        observed_z = (benchmark.sigma * benchmark.terminal_gradient(tensor)).numpy()
        assert np.allclose(observed_u, equation.terminal(x), rtol=0.0, atol=2e-14)
        assert np.allclose(
            observed_z, equation.terminal_gradient_z(x), rtol=0.0, atol=2e-14
        )


def test_torch_b6_adapter_matches_numpy() -> None:
    beta, gamma, _payload = load_b6_selection(100)
    equation = make_b6_equation(100, beta, gamma)
    adapter = B6TorchEquation(equation, torch.device("cpu"))
    t, x, _ = fixed_points("B6", 100)
    observed = adapter.exact_state(
        torch.tensor(t[:9], dtype=torch.float64),
        torch.tensor(x[:9], dtype=torch.float64),
    ).numpy()
    assert np.allclose(
        observed, equation.exact_state(t[:9], x[:9]), rtol=2e-14, atol=2e-14
    )
