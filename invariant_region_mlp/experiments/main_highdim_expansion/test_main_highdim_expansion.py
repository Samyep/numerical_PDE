from __future__ import annotations

from dataclasses import asdict
import math

import numpy as np
import torch

from .benchmarks import (
    B6Equation,
    dense_parameters,
    fixed_b6_points,
    make_dense_fourier,
    terminal_gradient_torch,
)
from .protocol import (
    ACTOR_UPDATES,
    DIMENSIONS,
    GRADIENT_ACTOR_CONFIG,
    JOINT_ACTOR_CONFIG,
    METHODS,
    PDES,
    TEST_STATES,
    frozen_protocol,
)


def _laplacian(value: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
    gradient = torch.autograd.grad(value.sum(), x, create_graph=True)[0]
    result = torch.zeros_like(value)
    for index in range(x.shape[1]):
        result = result + torch.autograd.grad(
            gradient[:, index].sum(), x, create_graph=True
        )[0][:, index]
    return result


def test_protocol_is_the_requested_frozen_grid() -> None:
    assert DIMENSIONS == (20, 50, 100, 200, 400)
    assert PDES == ("b6_fullrank", "dense_fourier")
    assert METHODS == ("tangent_actor", "dpi", "pinn_hte", "sdgd_pinn")
    assert TEST_STATES == 960
    assert ACTOR_UPDATES == 4
    assert asdict(GRADIENT_ACTOR_CONFIG) == {
        "training_points": 8192,
        "validation_points": 512,
        "paths_per_label": 128,
        "label_time_steps": 32,
        "optimizer_steps": 5000,
        "width": 384,
        "depth": 4,
        "activation": "silu",
        "learning_rate": 0.00015,
        "batch_size": 512,
        "parent_batch": 8,
        "path_chunk": 32,
        "checkpoint_interval": 250,
        "gradient_clip": 10.0,
        "training_distribution": "function",
        "label_estimator": "mcv",
    }
    joint = asdict(JOINT_ACTOR_CONFIG)
    assert joint["training_points"] == 8192
    assert joint["paths_per_label"] == 128
    assert joint["optimizer_steps"] == 5000
    assert joint["width"] == 384 and joint["depth"] == 4
    assert len(frozen_protocol()["protocol_sha256"]) == 64


def test_dense_parameter_draw_matches_validated_100d_pilot() -> None:
    driver = make_dense_fourier(100)
    assert driver.metadata()["parameter_sha256"] == (
        "b96804e2f7c20cbd4e72c0f961a4ebfc21bf46cae120d198b6afe6ce2e2684b6"
    )
    for d in DIMENSIONS:
        omega, signs, phases = dense_parameters(d)
        assert omega.shape == (2 * d, d)
        assert np.max(np.abs(np.linalg.norm(omega, axis=1) - 1.0)) < 1e-14
        assert set(np.unique(signs)) <= {-1.0, 1.0}
        assert np.all((phases >= 0.0) & (phases <= 2.0 * np.pi))


def test_b6_banks_match_validated_100d_and_400d_banks() -> None:
    from .benchmarks import array_sha256

    expected = {
        100: "84b2132607bc6b17e44adbea838915f8ce657813f6a0987e34deca308660f1c1",
        400: "0df59abc155e90aa674532760578dc2563cbdb4aade97a849961c8a8edd9f2f9",
    }
    for d, digest in expected.items():
        t, x, mask = fixed_b6_points(d)
        assert len(t) == 1200 and np.count_nonzero(~mask) == 960
        assert array_sha256(t, x, mask) == digest


def test_b6_terminal_gradient_and_pde_identity() -> None:
    equation = B6Equation(6, beta=1.5, gamma=1.0)
    rng = np.random.default_rng(31)
    t_np = rng.uniform(0.0, 0.24, 4)
    x_np = rng.uniform(-1.0, 1.0, (4, 6))
    state = equation.exact_state(t_np, x_np)
    assert np.max(np.abs(equation.exact_u(np.full(4, equation.T), x_np) - equation.terminal(x_np))) < 1e-13
    h = 1e-5
    finite = np.empty_like(x_np)
    for j in range(equation.d):
        plus, minus = x_np.copy(), x_np.copy()
        plus[:, j] += h
        minus[:, j] -= h
        finite[:, j] = (
            equation.exact_u(t_np, plus) - equation.exact_u(t_np, minus)
        ) / (2.0 * h)
    assert np.max(np.abs(math.sqrt(2.0) * finite - state[:, 1:])) < 2e-9

    b = torch.tensor(equation.b.copy(), dtype=torch.float64)
    n2 = torch.tensor(equation.n2.copy(), dtype=torch.float64)
    t = torch.as_tensor(t_np, dtype=torch.float64).requires_grad_(True)
    x = torch.as_tensor(x_np, dtype=torch.float64).requires_grad_(True)
    logits = x @ b.T + (equation.T - t)[:, None] * n2
    u = -(torch.logsumexp(logits, dim=1) - math.log(equation.K))
    u_t = torch.autograd.grad(u.sum(), t, create_graph=True)[0]
    grad = torch.autograd.grad(u.sum(), x, create_graph=True)[0]
    residual = u_t + _laplacian(u, x) - torch.sum(grad.square(), dim=1)
    assert float(residual.detach().abs().max()) < 2e-11


def test_dense_terminal_gradient_transform_and_pde_identity() -> None:
    driver = make_dense_fourier(8)
    rng = np.random.default_rng(73)
    t_np = rng.uniform(0.0, 0.49, 4)
    x_np = rng.uniform(-1.0, 1.0, (4, 8))
    state = driver.exact_state(t_np, x_np)
    assert np.max(np.abs(driver.exact_u(np.full(4, driver.T), x_np) - driver.terminal_numpy(x_np))) < 1e-13
    h = 1e-5
    finite = np.empty_like(x_np)
    for j in range(driver.d):
        plus, minus = x_np.copy(), x_np.copy()
        plus[:, j] += h
        minus[:, j] -= h
        finite[:, j] = (
            driver.exact_u(t_np, plus) - driver.exact_u(t_np, minus)
        ) / (2.0 * h)
    assert np.max(np.abs(math.sqrt(2.0) * finite - state[:, 1:])) < 2e-9

    t = torch.as_tensor(t_np, dtype=torch.float64).requires_grad_(True)
    x = torch.as_tensor(x_np, dtype=torch.float64).requires_grad_(True)
    u = driver.exact_u_torch(t, x)
    u_t = torch.autograd.grad(u.sum(), t, create_graph=True)[0]
    grad = torch.autograd.grad(u.sum(), x, create_graph=True)[0]
    z = math.sqrt(2.0) * grad
    residual = u_t + _laplacian(u, x) + 0.5 * torch.tanh(u) * z.square().sum(1)
    assert float(residual.detach().abs().max()) < 2e-11
    omega = torch.as_tensor(driver.omega, dtype=torch.float64)
    signs = torch.as_tensor(driver.signs, dtype=torch.float64)
    phases = torch.as_tensor(driver.phases, dtype=torch.float64)
    attenuation = torch.exp(
        -(driver.T - t)[:, None]
        * torch.as_tensor(driver.omega_norm_squared, dtype=torch.float64)
    )
    w = 0.5 + driver.scale * torch.sum(
        attenuation * signs * torch.sin(x @ omega.T + phases), dim=1
    )
    assert float(torch.max(torch.abs(torch.sinh(u) - w)).detach()) < 2e-13
    w_t = torch.autograd.grad(w.sum(), t, create_graph=True)[0]
    heat_residual = w_t + _laplacian(w, x)
    assert float(heat_residual.detach().abs().max()) < 2e-11


def test_benchmark_terminal_gradient_adapter() -> None:
    dense = make_dense_fourier(8)
    state = torch.randn(5, 8, dtype=torch.float64)
    assert torch.allclose(
        terminal_gradient_torch(dense, state),
        dense.terminal_gradient_torch(state),
        rtol=0.0,
        atol=0.0,
    )


def test_dense_joint_tangent_coefficients_match_autodiff() -> None:
    driver = make_dense_fourier(7)
    generator = torch.Generator().manual_seed(88)
    r = torch.randn(5, dtype=torch.float64, generator=generator, requires_grad=True)
    q = torch.randn(5, 7, dtype=torch.float64, generator=generator, requires_grad=True)
    f = driver.driver(r, q)
    expected_c = torch.autograd.grad(f.sum(), r, create_graph=True)[0]
    expected_b = torch.autograd.grad(f.sum(), q, create_graph=True)[0]
    c, b, h = driver.tangent_coefficients(r, q)
    assert torch.max(torch.abs(c - expected_c)).item() < 2e-13
    assert torch.max(torch.abs(b - expected_b)).item() < 2e-13
    assert torch.max(torch.abs(c * r + (b * q).sum(1) + h - f)).item() < 2e-13
