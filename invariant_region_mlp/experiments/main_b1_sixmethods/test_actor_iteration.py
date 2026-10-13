"""Correctness checks for teacher-free Actor Iteration."""

from __future__ import annotations

import numpy as np
import torch

from .actor import ActorNet
from .actor_iteration import (
    coarsen_normals,
    evaluate_policy_mc,
    evaluate_policy_mc_with_analytic_gap,
    evaluate_policy_mc_with_martingale_cv,
    policy_payoffs_with_analytic_gap_from_normals,
    policy_payoffs_with_martingale_from_normals,
    policy_value_and_z_from_normals,
    refine_antithetic_normals,
    terminal_pathwise_labels,
    torch_terminal,
    torch_terminal_z,
    torch_exact_z,
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


def test_martingale_cv_has_zero_variance_for_linear_terminal_policy() -> None:
    class LinearTerminalEquation:
        d = 4
        T = 0.25
        sigma = float(np.sqrt(2.0))
        A = np.asarray([[0.4, -0.2, 0.1, 0.3]], dtype=np.float64)
        logc = np.zeros(1, dtype=np.float64)

    equation = LinearTerminalEquation()
    probe = torch.zeros((1, equation.d))
    exact_policy = torch_terminal_z(equation, probe)[0]
    actor = _constant_actor(equation.d, exact_policy)
    t = torch.tensor([0.03, 0.11], dtype=torch.float32)
    x = torch.randn(2, equation.d, generator=torch.Generator().manual_seed(21))
    normals = torch.randn(
        2,
        64,
        8,
        equation.d,
        generator=torch.Generator().manual_seed(22),
    )
    raw, corrected, martingale = policy_payoffs_with_martingale_from_normals(
        equation, actor, t, x, normals
    )
    assert torch.max(torch.var(raw, dim=1)) > 1.0e-4
    assert torch.max(torch.var(corrected, dim=1)) < 2.0e-12
    torch.testing.assert_close(raw - corrected, martingale)
    direction = torch.as_tensor(equation.A[0], dtype=torch.float32)
    expected = -x @ direction - (equation.T - t) * torch.sum(
        direction.square()
    )
    torch.testing.assert_close(
        corrected.mean(dim=1), expected, rtol=2.0e-6, atol=2.0e-6
    )
    assert all(parameter.grad is None for parameter in actor.parameters())


def test_analytic_policy_gap_identity_for_linear_terminal_problem() -> None:
    class LinearTerminalEquation:
        d = 3
        T = 0.25
        sigma = float(np.sqrt(2.0))
        A = np.asarray([[0.3, -0.4, 0.2]], dtype=np.float64)
        logc = np.zeros(1, dtype=np.float64)

    equation = LinearTerminalEquation()
    policy = torch.tensor([0.1, -0.2, 0.05])
    actor = _constant_actor(equation.d, policy)
    t = torch.tensor([0.04, 0.12], dtype=torch.float32)
    x = torch.randn(2, equation.d, generator=torch.Generator().manual_seed(25))
    half = torch.randn(
        2, 32, 8, equation.d, generator=torch.Generator().manual_seed(26)
    )
    normals = torch.cat((half, -half), dim=1)
    raw, corrected, _, policy_gap = policy_payoffs_with_analytic_gap_from_normals(
        equation, actor, t, x, normals
    )
    exact_policy = torch_exact_z(equation, t, x)
    expected_policy = torch_terminal_z(equation, x)
    torch.testing.assert_close(exact_policy, expected_policy)
    expected_gap = (
        0.5
        * torch.sum((policy[None, :] - exact_policy).square(), dim=-1)
        * (equation.T - t)
    )
    torch.testing.assert_close(
        policy_gap.mean(dim=1), expected_gap, rtol=2.0e-6, atol=2.0e-6
    )
    exact_value = torch_terminal(equation, x) - (
        equation.T - t
    ) * torch.sum(torch.as_tensor(equation.A[0], dtype=x.dtype).square())
    torch.testing.assert_close(
        corrected.mean(dim=1) - exact_value,
        expected_gap,
        rtol=3.0e-6,
        atol=3.0e-6,
    )
    assert torch.max(torch.var(raw, dim=1)) > 1.0e-4


def test_brownian_bridge_refinement_preserves_coarse_antithetic_paths() -> None:
    coarse_generator = torch.Generator().manual_seed(31)
    bridge_generator = torch.Generator().manual_seed(32)
    coarse = torch.randn(3, 4, 5, 2, generator=coarse_generator)
    coarse[:, 2:] = -coarse[:, :2]
    fine = refine_antithetic_normals(
        coarse, 4, generator=bridge_generator
    )
    recovered = coarsen_normals(fine, 5)
    torch.testing.assert_close(recovered, coarse, rtol=2.0e-6, atol=2.0e-6)
    torch.testing.assert_close(fine[:, 2:], -fine[:, :2])
    assert coarsen_normals(fine, 10).shape == (3, 4, 10, 2)


def test_paired_cv_evaluator_reuses_paths_and_zero_actor_is_unchanged() -> None:
    equation = B1Equation(d=4, k=2)
    actor = _constant_actor(equation.d, torch.zeros(equation.d))
    rng = np.random.default_rng(41)
    t = rng.uniform(0.0, equation.T, 5)
    x = rng.standard_normal((5, equation.d))
    raw, corrected, raw_variance, corrected_variance, metadata = (
        evaluate_policy_mc_with_martingale_cv(
            equation,
            actor,
            t,
            x,
            time_steps=16,
            paths=8,
            repetition=0,
            seed_namespace="actor_iteration_direct_test",
            device=torch.device("cpu"),
            parent_batch_size=3,
            path_chunk_size=4,
        )
    )
    np.testing.assert_array_equal(raw, corrected)
    np.testing.assert_array_equal(raw_variance, corrected_variance)
    assert metadata["raw_and_corrected_share_paths"] is True
    assert metadata["additional_actor_calls_for_control_variate"] == 0
    legacy, _, _ = evaluate_policy_mc(
        equation,
        actor,
        t,
        x,
        time_steps=16,
        paths=8,
        repetition=0,
        seed_namespace="actor_iteration_direct_test",
        device=torch.device("cpu"),
        parent_batch_size=3,
        path_chunk_size=4,
        compute_z=False,
    )
    np.testing.assert_array_equal(raw, legacy)

    (
        diagnostic_raw,
        diagnostic_corrected,
        policy_gap,
        diagnostic_raw_variance,
        diagnostic_corrected_variance,
        policy_gap_variance,
        diagnostic_metadata,
    ) = evaluate_policy_mc_with_analytic_gap(
        equation,
        actor,
        t,
        x,
        time_steps=16,
        paths=8,
        repetition=0,
        seed_namespace="actor_iteration_direct_test",
        device=torch.device("cpu"),
        parent_batch_size=3,
        path_chunk_size=4,
    )
    np.testing.assert_array_equal(diagnostic_raw, raw)
    np.testing.assert_array_equal(diagnostic_corrected, corrected)
    np.testing.assert_array_equal(diagnostic_raw_variance, raw_variance)
    np.testing.assert_array_equal(
        diagnostic_corrected_variance, corrected_variance
    )
    assert np.all(policy_gap >= 0.0)
    assert np.all(policy_gap_variance >= 0.0)
    assert diagnostic_metadata["analytic_reference_use"] == (
        "post-training diagnostic only"
    )
