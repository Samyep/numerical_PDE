"""Small, independent checks of signs, terminal constraints, and gradients."""

from __future__ import annotations

import numpy as np
import torch

from .actor_method import ActorPhase, evaluate_policy_cv, policy_gradient_labels
from .baselines import (
    AmortizedDeepBSDE,
    HardTerminalPINN,
    TerminalAwareDPI,
    sdgd_objective,
)
from .equations import HJBLogEquation, torch_terminal, torch_terminal_z
from .models import ActorNet


def _zero_actor(d: int) -> ActorNet:
    actor = ActorNet(d, width=8, depth=1).to(dtype=torch.float32)
    with torch.no_grad():
        for parameter in actor.parameters():
            parameter.zero_()
    actor.eval()
    return actor


def test_actor_cv_terminal_and_source_freezing() -> None:
    equation = HJBLogEquation(d=4)
    actor = _zero_actor(4)
    x = np.arange(12, dtype=np.float64).reshape(3, 4) / 10.0
    t = np.full(3, equation.T)
    value, variance, _ = evaluate_policy_cv(
        equation,
        actor,
        t,
        x,
        paths=4,
        time_steps=2,
        seed=7,
        device=torch.device("cpu"),
        parent_batch=3,
        path_chunk=4,
    )
    np.testing.assert_allclose(value, equation.terminal(x), rtol=2e-7, atol=2e-7)
    np.testing.assert_allclose(variance, 0.0, atol=1e-14)
    _, labels, label_metadata = policy_gradient_labels(
        equation,
        actor,
        np.zeros(2),
        x[:2],
        paths=4,
        time_steps=2,
        seed=11,
        device=torch.device("cpu"),
        parent_batch=2,
        path_chunk=4,
    )
    assert labels.shape == (2, 4)
    assert np.all(np.isfinite(labels))
    assert label_metadata["path_chunk"] == 4
    assert all(parameter.grad is None for parameter in actor.parameters())


def test_small_path_budget_uses_matching_chunk() -> None:
    """Regression test for the 8-path smoke configuration."""

    equation = HJBLogEquation(d=2)
    actor = _zero_actor(2)
    value, variance, metadata = evaluate_policy_cv(
        equation,
        actor,
        np.zeros(3),
        np.zeros((3, 2)),
        paths=8,
        time_steps=2,
        seed=19,
        device=torch.device("cpu"),
        parent_batch=3,
        path_chunk=min(64, 8),
    )
    assert value.shape == variance.shape == (3,)
    assert metadata["paths"] == 8


def test_baseline_terminal_constraints_and_sdgd_backward() -> None:
    equation = HJBLogEquation(d=3)
    x = torch.randn(5, 3)
    terminal_t = torch.full((5,), equation.T)
    pinn = HardTerminalPINN(3, width=8, depth=1)
    dpi = TerminalAwareDPI(3, width=8, depth=1)
    np.testing.assert_allclose(
        pinn.value(equation, terminal_t, x).detach().numpy(),
        torch_terminal(equation, x).numpy(),
        atol=1e-7,
    )
    dpi_u, dpi_z = dpi.value_z(equation, terminal_t, x)
    np.testing.assert_allclose(dpi_u.detach().numpy(), torch_terminal(equation, x).numpy(), atol=1e-7)
    np.testing.assert_allclose(dpi_z.detach().numpy(), torch_terminal_z(equation, x).numpy(), atol=1e-7)

    generator = torch.Generator().manual_seed(3)
    objective, _ = sdgd_objective(
        pinn,
        equation,
        torch.rand(4),
        torch.randn(4, 3),
        coordinate_budget=2,
        generator=generator,
    )
    assert torch.isfinite(objective)
    objective.backward()
    assert any(parameter.grad is not None for parameter in pinn.parameters())


def test_deep_bsde_positive_quadratic_drift_sign() -> None:
    equation = HJBLogEquation(d=2)
    model = AmortizedDeepBSDE(2, width=4, depth=1)
    with torch.no_grad():
        for parameter in model.parameters():
            parameter.zero_()
        # A constant z=(1,1) gives +0.5*||z||^2*T = +T when dW=0.
        model.z.model[-1].bias.fill_(1.0)
    terminal_y, terminal_x = model.rollout(
        equation, torch.zeros(3, 2), torch.zeros(3, 4, 2)
    )
    np.testing.assert_allclose(terminal_y.detach().numpy(), equation.T, atol=1e-7)
    np.testing.assert_allclose(terminal_x.detach().numpy(), 0.0, atol=1e-7)
