"""Numerical correctness tests for the compute-frontier equations."""

from __future__ import annotations

import numpy as np
import torch

from .equations import HJBLogEquation, make_equation, torch_terminal, torch_terminal_z


def test_hjb_log_terminal_and_torch_derivative() -> None:
    equation = HJBLogEquation(d=7)
    rng = np.random.default_rng(31)
    x = rng.standard_normal((6, 7))
    np.testing.assert_allclose(equation.exact_u(equation.T, x), equation.terminal(x), rtol=0, atol=1e-13)
    np.testing.assert_allclose(equation.exact_z(equation.T, x), equation.terminal_z(x), rtol=0, atol=1e-13)
    tx = torch.tensor(x, dtype=torch.float64, requires_grad=True)
    terminal = torch_terminal(equation, tx)
    grad = torch.autograd.grad(terminal.sum(), tx)[0]
    np.testing.assert_allclose(
        np.sqrt(2.0) * grad.detach().numpy(),
        torch_terminal_z(equation, tx.detach()).numpy(),
        rtol=2e-13,
        atol=2e-13,
    )


def test_hjb_log_quadrature_matches_adaptive_reference() -> None:
    equation = HJBLogEquation(d=12, quadrature_nodes=320, quadrature_upper=90.0)
    rng = np.random.default_rng(91)
    x = rng.standard_normal((4, 12))
    times = np.array([0.0, 0.2, 0.77, 1.0])
    values = equation.exact_u(times, x)
    z = equation.exact_z(times, x)
    for index in range(len(times)):
        reference_u, reference_z, _ = equation.adaptive_reference(times[index], x[index])
        np.testing.assert_allclose(values[index], reference_u, rtol=2e-11, atol=2e-11)
        np.testing.assert_allclose(z[index], reference_z, rtol=1e-9, atol=1e-10)


def test_b1_factory_sign_convention() -> None:
    equation = make_equation("b1")
    rng = np.random.default_rng(5)
    t = rng.uniform(0.01, equation.T - 0.01, 8)
    x = rng.standard_normal((8, equation.d))
    np.testing.assert_allclose(equation.analytic_residual(t, x), 0.0, atol=2e-14)
