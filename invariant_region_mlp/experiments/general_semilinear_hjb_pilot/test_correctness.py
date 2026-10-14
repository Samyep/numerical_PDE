"""Mandatory mathematical and implementation checks for the rapid pilot."""

from __future__ import annotations

import inspect

import numpy as np
import pytest
import torch
from scipy.optimize import minimize_scalar

from . import core
from .core import _pathwise_chunk, antithetic_normals, evaluate_policy_plain_mc, frozen
from .models import ActorNet, ZeroPolicy
from .problems import (
    BoundedControlProblem,
    LQGProblem,
    QuarticControlProblem,
    quartic_root_numpy,
)
from .reference import RiccatiReference, validate_scalar_reference


def test_lqg_sign_and_greedy_formula() -> None:
    problem = LQGProblem(d=3)
    reference = RiccatiReference(problem)
    t = np.linspace(0.05, 0.45, 9)
    h = 1e-5
    pp, _ = reference.coefficients(t + h)
    pm, _ = reference.coefficients(t - h)
    P, _ = reference.coefficients(t)
    derivative = (pp - pm) / (2.0 * h)
    residual = derivative - 2.0 * problem.gamma * P - P**2 / problem.R + problem.Q
    assert np.max(np.abs(residual)) < 2e-7

    gradient = torch.randn(11, problem.d)
    action = problem.greedy_action(torch.zeros(11), torch.zeros_like(gradient), gradient)
    assert torch.allclose(action, -gradient / problem.R)
    perturbation = torch.randn_like(action)
    base = (action * gradient).sum(-1) + 0.5 * problem.R * action.square().sum(-1)
    other = (perturbation * gradient).sum(-1) + 0.5 * problem.R * perturbation.square().sum(-1)
    assert torch.all(base <= other + 1e-6)


def test_quartic_solver_matches_foc_and_direct_minimization() -> None:
    problem = QuarticControlProblem(d=4)
    p = np.linspace(-5.0, 5.0, 41)
    action = quartic_root_numpy(p, problem.lam)
    assert np.max(np.abs(action + problem.lam * action**3 + p)) < 1e-11
    direct = np.array([
        minimize_scalar(
            lambda a: a * value + 0.5 * a**2 + 0.25 * problem.lam * a**4,
            bounds=(-8.0, 8.0),
            method="bounded",
            options={"xatol": 1e-12},
        ).x
        for value in p
    ])
    assert np.max(np.abs(action - direct)) < 2e-6


def test_bounded_actor_is_feasible_and_greedy_clips() -> None:
    problem = BoundedControlProblem(d=5, A=0.35)
    raw = torch.linspace(-100.0, 100.0, 1001).reshape(-1, 1)
    feasible = problem.feasible_action(raw)
    assert torch.max(feasible).item() <= problem.A
    assert torch.min(feasible).item() >= -problem.A
    gradient = torch.outer(torch.linspace(-2.0, 2.0, 101), torch.as_tensor(problem.beta).float())
    target = problem.greedy_action(torch.zeros(101), torch.zeros(101, 5), gradient)
    expected = torch.clamp(-torch.linspace(-2.0, 2.0, 101), -problem.A, problem.A)
    assert torch.allclose(target[:, 0], expected, atol=2e-6)


class _StateDependentPolicy(torch.nn.Module):
    def __init__(self, problem: QuarticControlProblem) -> None:
        super().__init__()
        self.problem = problem
        self.weight = torch.nn.Parameter(torch.tensor(0.23))

    def predict(self, t: torch.Tensor, x: torch.Tensor, T: float) -> torch.Tensor:
        del T
        return (self.weight * torch.tanh(x[:, :1] + 0.3 * t[:, None]))


def test_pathwise_vjp_matches_crn_central_difference_and_freezes_source() -> None:
    problem = QuarticControlProblem(d=3, T=0.3)
    policy = _StateDependentPolicy(problem)
    generator = torch.Generator().manual_seed(421)
    normals = antithetic_normals((1, 512, 12, problem.d), generator=generator, device=torch.device("cpu"))
    t = torch.tensor([0.07])
    x = torch.tensor([[0.2, -0.4, 0.1]])
    with frozen(policy):
        value, gradient = _pathwise_chunk(problem, policy, t, x, normals)
        assert not policy.weight.requires_grad
        assert policy.weight.grad is None
    assert policy.weight.requires_grad
    eps = 2e-3
    finite_difference = []
    with frozen(policy):
        for coordinate in range(problem.d):
            direction = torch.zeros_like(x)
            direction[0, coordinate] = eps
            plus, _ = _pathwise_chunk(problem, policy, t, x + direction, normals)
            minus, _ = _pathwise_chunk(problem, policy, t, x - direction, normals)
            finite_difference.append(((plus - minus) / (2.0 * eps)).item())
    assert np.max(np.abs(gradient.numpy()[0] - np.asarray(finite_difference))) < 3e-3
    assert np.linalg.norm(gradient.numpy()) > 1e-3
    assert np.isfinite(value.item())


@pytest.mark.parametrize("name", ["quartic", "bounded"])
def test_scalar_reference_refines_and_matches_terminal(name: str) -> None:
    problem = QuarticControlProblem(d=3) if name == "quartic" else BoundedControlProblem(d=3)
    reference, audit = validate_scalar_reference(problem)
    assert audit["terminal_max_abs"] < 1e-12
    assert audit["value_refinement_rmse"] < 2e-3
    assert audit["action_refinement_rmse"] < 2e-2
    values = reference.value(np.zeros(41), np.linspace(-4.0, 4.0, 41))
    assert np.all(np.isfinite(values))


@pytest.mark.parametrize("problem", [LQGProblem(d=4), QuarticControlProblem(d=4), BoundedControlProblem(d=4)])
def test_policy_evaluation_at_terminal_is_exact_and_finite(problem) -> None:
    rng = np.random.default_rng(99)
    x = rng.standard_normal((7, problem.d))
    t = np.full(7, problem.T)
    output, metadata = evaluate_policy_plain_mc(
        problem,
        ZeroPolicy(problem),
        t,
        x,
        paths=8,
        time_steps=4,
        seed=44,
        device=torch.device("cpu"),
        parent_batch=7,
        path_chunk=8,
    )
    assert np.max(np.abs(output["value"] - problem.terminal_numpy(x))) < 2e-6
    assert np.max(output["mc_variance_of_mean"]) < 1e-10
    assert np.all(np.isfinite(output["value"]))
    assert metadata["martingale_control_variate"] is False


def test_training_core_has_no_reference_leakage() -> None:
    source = inspect.getsource(core)
    forbidden = ("RiccatiReference", "ScalarReference", "solve_scalar_reference")
    assert all(token not in source for token in forbidden)
    assert "analytic_nonterminal_reference_used\": False" in source

