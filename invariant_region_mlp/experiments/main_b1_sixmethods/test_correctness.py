"""Unit checks required before Stage A."""

from __future__ import annotations

import math

import numpy as np
import pytest
import torch

from invariant_region_mlp.experiments.main_b1_sixmethods.actor import ActorNet
from invariant_region_mlp.experiments.main_b1_sixmethods.actor_mlp import (
    actor_linearization_identity,
    evaluate_actor_mlp,
)
from invariant_region_mlp.experiments.main_b1_sixmethods.deep_bsde_adapter import GlobalDeepBSDE
from invariant_region_mlp.experiments.main_b1_sixmethods.dpi_adapter import (
    DPINet,
    generate_picard_labels,
)
from invariant_region_mlp.experiments.main_b1_sixmethods.equation import B1Equation
from invariant_region_mlp.experiments.mechanism_suite.equations import RidgeLSEHJB
from invariant_region_mlp.experiments.main_b1_sixmethods.sdgd_pinn_adapter import (
    HardTerminalPINN,
    _sampled_laplacian,
    sdgd_loss,
)


@pytest.mark.parametrize("d,k", [(20, 1), (20, 5), (30, 10)])
def test_terminal_and_analytic_residual(d: int, k: int) -> None:
    equation = B1Equation(d=d, k=k)
    rng = np.random.default_rng(7 + k)
    x = rng.standard_normal((16, d))
    t = rng.uniform(0.0, equation.T, 16)
    np.testing.assert_allclose(equation.exact_u(equation.T, x), equation.terminal(x), atol=2e-14)
    np.testing.assert_allclose(equation.analytic_residual(t, x), 0.0, atol=2e-13)


def test_k1_exactly_reproduces_p1() -> None:
    d = 100
    b1 = B1Equation(d=d, k=1)
    p1 = RidgeLSEHJB(d=d)
    rng = np.random.default_rng(123)
    x = rng.standard_normal((32, d))
    t = rng.uniform(0.0, b1.T, 32)
    np.testing.assert_array_equal(b1.B[:, 0], p1.w)
    np.testing.assert_array_equal(b1.b[:, 0], p1.lambdas)
    np.testing.assert_allclose(b1.exact_state(t, x), p1.exact_state(t, x), atol=2e-14)


def test_closed_form_gradient_and_pde_with_autodiff() -> None:
    equation = B1Equation(d=12, k=5)
    dtype = torch.float64
    x = torch.randn(7, equation.d, dtype=dtype, requires_grad=True)
    t = torch.rand(7, dtype=dtype, requires_grad=True) * equation.T
    A = torch.as_tensor(equation.A, dtype=dtype)
    n2 = torch.as_tensor(equation.n2, dtype=dtype)
    logits = x @ A.T + (equation.T - t[:, None]) * n2 - math.log(equation.K)
    u = -torch.logsumexp(logits, dim=-1)
    grad = torch.autograd.grad(u.sum(), x, create_graph=True)[0]
    u_t = torch.autograd.grad(u.sum(), t, create_graph=True)[0]
    laplacian = torch.zeros_like(u)
    for j in range(equation.d):
        second = torch.autograd.grad(grad[:, j].sum(), x, retain_graph=True)[0][:, j]
        laplacian = laplacian + second
    z = math.sqrt(2.0) * grad
    residual = u_t + laplacian - 0.5 * torch.sum(z.square(), dim=-1)
    expected_z = equation.exact_z(t.detach().numpy(), x.detach().numpy())
    np.testing.assert_allclose(z.detach().numpy(), expected_z, rtol=2e-12, atol=2e-12)
    np.testing.assert_allclose(residual.detach().numpy(), 0.0, atol=2e-12)


@pytest.mark.parametrize("k", [2, 5, 10])
def test_coefficients_and_hessian_have_full_active_rank(k: int) -> None:
    equation = B1Equation(d=20, k=k)
    assert equation.affine_rank == k
    hessian = equation.exact_active_hessian(0.11, np.zeros((1, 20)))[0]
    assert np.linalg.matrix_rank(hessian, tol=1e-11) == k
    assert np.count_nonzero(np.abs(equation.B) < 1e-14) == 0
    assert np.count_nonzero(np.abs(equation.b) < 1e-14) == 0


def test_terminal_gradient_and_sigma_scaling() -> None:
    equation = B1Equation(d=16, k=5)
    x = np.random.default_rng(9).standard_normal((13, 16))
    np.testing.assert_allclose(
        equation.exact_z(equation.T, x),
        math.sqrt(2.0) * equation.terminal_grad(x),
        atol=2e-14,
    )


def test_actor_linearization_identity() -> None:
    equation = B1Equation(d=16, k=5)
    actor = ActorNet(16)
    error = actor_linearization_identity(actor, equation, torch.device("cpu"))
    assert error <= 2e-5


def test_actor_mlp_smoke_and_child_actor_calls() -> None:
    equation = B1Equation(d=8, k=3)
    actor = ActorNet(8, width=16, depth=2)
    rng = np.random.default_rng(81)
    output, cost = evaluate_actor_mlp(
        equation,
        actor,
        np.zeros(4),
        rng.standard_normal((4, 8)),
        n=2,
        M=2,
        repetition=0,
        device=torch.device("cpu"),
    )
    assert np.all(np.isfinite(output))
    assert cost["actor_calls"] > 1
    assert cost["actor_points"] > 4


def test_dpi_label_smoke_without_exact_reference() -> None:
    equation = B1Equation(d=8, k=3)
    model = DPINet(8, width=16, depth=2, activation="tanh")
    generator = torch.Generator().manual_seed(92)
    t = torch.rand(5, generator=generator) * 0.2
    x = torch.randn((5, 8), generator=generator)
    value, z = generate_picard_labels(
        equation, model, t, x, mc_samples=4, generator=generator
    )
    assert value.shape == (5,)
    assert z.shape == (5, 8)
    assert torch.all(torch.isfinite(value)) and torch.all(torch.isfinite(z))


def test_deep_bsde_rollout_sign() -> None:
    equation = B1Equation(d=4, k=2)
    model = GlobalDeepBSDE(
        4, width=8, depth=1, activation="tanh", z_bound=8.0
    )
    for parameter in model.parameters():
        torch.nn.init.zeros_(parameter)
    # Make Z identically one through the final bias.
    final_linear = [layer for layer in model.z.model if isinstance(layer, torch.nn.Linear)][-1]
    final_linear.bias.data.fill_(1.0)
    x0 = torch.zeros((3, 4))
    normals = torch.zeros((3, 2, 4))
    y, x_terminal = model.rollout(equation, x0, normals)
    # Integral 0.5*|z|^2 over T = 0.5*4*0.25 = 0.5.
    np.testing.assert_allclose(y.detach().numpy(), 0.5, atol=1e-7)
    np.testing.assert_allclose(x_terminal.detach().numpy(), 0.0, atol=0.0)


def test_sdgd_laplacian_unbiased_and_terminal_exact() -> None:
    d = 20
    x = torch.randn(7, d, dtype=torch.float64, requires_grad=True)
    coefficients = torch.arange(1, d + 1, dtype=torch.float64)
    value = torch.sum(coefficients * x.square(), dim=-1)
    gradient = torch.autograd.grad(value.sum(), x, create_graph=True)[0]
    generator = torch.Generator().manual_seed(101)
    estimates = []
    for _ in range(2000):
        coordinates = torch.randperm(d, generator=generator)[:4]
        estimates.append(_sampled_laplacian(gradient, x, coordinates, d).detach())
    estimate = torch.stack(estimates).mean(dim=0)
    exact = 2.0 * coefficients.sum()
    assert torch.max(torch.abs(estimate - exact)) / exact < 0.03

    equation = B1Equation(d=d, k=5)
    model = HardTerminalPINN(d, width=16, depth=2, activation="tanh")
    terminal_t = torch.full((7,), equation.T, dtype=torch.float32)
    terminal_x = x.detach().float()
    terminal = model.value(equation, terminal_t, terminal_x)
    np.testing.assert_allclose(
        terminal.detach().numpy(), equation.terminal(terminal_x.numpy()), atol=2e-6
    )


def test_sdgd_loss_is_finite() -> None:
    equation = B1Equation(d=8, k=3)
    model = HardTerminalPINN(8, width=16, depth=2, activation="tanh")
    generator = torch.Generator().manual_seed(103)
    t = torch.rand(4, generator=generator) * 0.2
    x = torch.randn((4, 8), generator=generator)
    loss, diagnostics = sdgd_loss(
        model, equation, t, x, coordinate_budget=3, generator=generator
    )
    assert torch.isfinite(loss)
    assert all(math.isfinite(value) for value in diagnostics.values())
