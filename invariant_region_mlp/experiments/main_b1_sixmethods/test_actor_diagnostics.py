"""Correctness checks for the supplemental actor-bottleneck diagnostics."""

from __future__ import annotations

import numpy as np
import torch

from invariant_region_mlp.experiments.gpu_mlp.solver import TorchMLP

from .actor_diagnostics import OracleActorMLP, evaluate_actor_diagnostic
from .common import stable_seed
from .equation import B1Equation


def test_oracle_driver_matches_original_at_exact_state() -> None:
    equation = B1Equation(d=12, k=5)
    primary, auxiliary = np.random.SeedSequence(17).spawn(2)
    solver = OracleActorMLP(
        equation,
        2,
        torch.device("cpu"),
        seed_sequence=primary,
        aux_seed_sequence=auxiliary,
    )
    generator = torch.Generator().manual_seed(19)
    t = torch.rand(7, dtype=torch.float64, generator=generator) * equation.T
    x = torch.randn((7, equation.d), dtype=torch.float64, generator=generator)
    exact_state = solver.equation.exact_state(t, x)
    observed = solver._actor_value(t, x, exact_state)
    expected = solver.equation.generator(exact_state)
    torch.testing.assert_close(observed, expected, rtol=2e-13, atol=2e-13)


def test_oracle_actor_is_evaluated_at_recursive_child_states() -> None:
    equation = B1Equation(d=8, k=3)
    rng = np.random.default_rng(23)
    output, cost = evaluate_actor_diagnostic(
        equation,
        np.zeros(4),
        rng.standard_normal((4, equation.d)),
        method="oracle_actor",
        n=2,
        M=2,
        repetition=0,
        device=torch.device("cpu"),
    )
    assert np.all(np.isfinite(output))
    assert cost["oracle_reference_used_inside_solver"] is True
    assert cost["actor_calls"] > 1
    assert cost["actor_points"] > 4


def test_no_actor_diagnostic_is_exact_path_base_with_paired_seed() -> None:
    equation = B1Equation(d=8, k=3)
    rng = np.random.default_rng(29)
    t = np.zeros(5)
    x = rng.standard_normal((5, equation.d))
    observed, cost = evaluate_actor_diagnostic(
        equation,
        t,
        x,
        method="pathwise_baseline_no_actor",
        n=2,
        M=2,
        repetition=3,
        device=torch.device("cpu"),
    )
    seed = stable_seed("B1", "actor_mlp", equation.d, equation.k, 2, 2, 3, 0)
    primary, auxiliary = np.random.SeedSequence(seed).spawn(2)
    solver = TorchMLP(
        equation,
        2,
        "path_base",
        torch.device("cpu"),
        seed_sequence=primary,
        aux_seed_sequence=auxiliary,
    )
    with torch.no_grad():
        expected = solver.solve(
            2,
            torch.as_tensor(t, dtype=torch.float64),
            torch.as_tensor(x, dtype=torch.float64),
        ).numpy()
    np.testing.assert_array_equal(observed, expected)
    assert cost["driver"] == "original quadratic driver"
    assert cost["actor_calls"] == 0

