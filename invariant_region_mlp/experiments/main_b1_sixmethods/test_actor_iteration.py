"""Correctness checks for teacher-free Actor Iteration."""

from __future__ import annotations

import numpy as np
import torch

from .actor import ActorNet
from .actor_iteration import (
    policy_value_and_z_from_normals,
    terminal_pathwise_labels,
    torch_terminal,
    torch_terminal_z,
)
from .equation import B1Equation


def _constant_actor(d: int, value: torch.Tensor) -> ActorNet:
    actor = ActorNet(d, width=8, depth=1)
    for parameter in actor.parameters():
        torch.nn.init.zeros_(parameter)
    final = [layer for layer in actor.model if isinstance(layer, torch.nn.Linear)][-1]
    final.bias.data.copy_(value)
    return actor


def test_zero_policy_is_heat_expectation_and_terminal_is_exact() -> None:
    equation = B1Equation(d=4, k=2)
    actor = _constant_actor(4, torch.zeros(4))
    t = torch.tensor([0.03, 0.11], dtype=torch.float32)
    x = torch.randn(2, 4, generator=torch.Generator().manual_seed(3))
    normals = torch.randn(
        2, 6, 4, 4, generator=torch.Generator().manual_seed(4)
    )
    value, _ = policy_value_and_z_from_normals(
        equation, actor, t, x, normals, compute_z=False
    )
    dt = (equation.T - t) / normals.shape[2]
    terminal_state = x[:, None, :] + equation.sigma * torch.sqrt(dt)[
        :, None, None
    ] * normals.sum(dim=2)
    expected = torch_terminal(equation, terminal_state).mean(dim=1)
    torch.testing.assert_close(value, expected)

    terminal_t = torch.full((2,), equation.T)
    terminal_value, terminal_z = policy_value_and_z_from_normals(
        equation, actor, terminal_t, x, normals, compute_z=True
    )
    torch.testing.assert_close(terminal_value, torch_terminal(equation, x))
    assert terminal_z is not None
    torch.testing.assert_close(terminal_z, torch_terminal_z(equation, x))


def test_controlled_drift_sign_and_running_cost() -> None:
    equation = B1Equation(d=3, k=2)
    policy = torch.tensor([0.2, -0.1, 0.3])
    actor = _constant_actor(3, policy)
    t = torch.tensor([0.05])
    x = torch.tensor([[0.3, -0.4, 0.2]])
    normals = torch.zeros((1, 4, 5, 3))
    value, _ = policy_value_and_z_from_normals(
        equation, actor, t, x, normals, compute_z=False
    )
    duration = equation.T - float(t.item())
    terminal_state = x - equation.sigma * policy * duration
    expected = torch_terminal(equation, terminal_state) + (
        0.5 * torch.sum(policy.square()) * duration
    )
    torch.testing.assert_close(value, expected, rtol=2.0e-6, atol=2.0e-6)


def test_pathwise_gradient_matches_common_random_finite_difference() -> None:
    equation = B1Equation(d=3, k=2)
    torch.manual_seed(11)
    actor = ActorNet(3, width=12, depth=2)
    t = torch.tensor([0.07])
    x = torch.tensor([[0.2, -0.3, 0.4]])
    normals = torch.randn(
        1, 16, 4, 3, generator=torch.Generator().manual_seed(12)
    )
    _, z = policy_value_and_z_from_normals(
        equation, actor, t, x, normals, compute_z=True
    )
    assert z is not None
    epsilon = 2.0e-3
    finite_difference = torch.empty(3)
    for coordinate in range(3):
        direction = torch.zeros_like(x)
        direction[:, coordinate] = epsilon
        plus, _ = policy_value_and_z_from_normals(
            equation, actor, t, x + direction, normals, compute_z=False
        )
        minus, _ = policy_value_and_z_from_normals(
            equation, actor, t, x - direction, normals, compute_z=False
        )
        finite_difference[coordinate] = (plus - minus) / (2.0 * epsilon)
    torch.testing.assert_close(
        z[0] / equation.sigma,
        finite_difference,
        rtol=8.0e-3,
        atol=2.0e-3,
    )
    assert all(parameter.requires_grad for parameter in actor.parameters())
    assert all(parameter.grad is None for parameter in actor.parameters())


def test_terminal_pathwise_labels_are_finite_and_correctly_shaped() -> None:
    equation = B1Equation(d=5, k=2)
    rng = np.random.default_rng(15)
    t = rng.uniform(0.0, equation.T, 7)
    x = rng.standard_normal((7, equation.d))
    labels, metadata = terminal_pathwise_labels(
        equation,
        t,
        x,
        paths=8,
        repetition=0,
        seed_namespace="actor_iteration_test",
        device=torch.device("cpu"),
        parent_batch_size=3,
    )
    assert labels.shape == (7, equation.d)
    assert np.all(np.isfinite(labels))
    assert metadata["terminal_samples"] == 56
    assert metadata["nonterminal_mlp_used"] is False
    assert metadata["nonterminal_exact_gradient_used"] is False
