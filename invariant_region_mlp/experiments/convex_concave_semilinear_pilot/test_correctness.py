"""CPU preflight tests for the generic weighted tangent evaluator."""

from __future__ import annotations

import inspect
import math

import numpy as np
import pytest
import torch
from torch import nn

from . import core
from .core import _weighted_rollout_chunk, frozen
from .drivers import (
    AllenCahnDriver,
    ConcaveSemilinearDriver,
    JointConcaveDriver,
    JointConvexOriginalDriver,
    SignTransformedDriver,
    make_driver,
)
from .models import ConstantTangentActor
from .reference import validate_periodic_reference


@pytest.mark.parametrize("name", ["allen_cahn", "joint_concave", "joint_convex"])
def test_tangent_coefficients_and_autodiff(name: str) -> None:
    torch.manual_seed(10)
    driver = make_driver(name, d=5)
    r = driver.r_low + (driver.r_high - driver.r_low) * torch.rand(17, dtype=torch.float64)
    q = torch.randn(17, driver.d, dtype=torch.float64) if driver.has_gradient else None
    analytic = driver.tangent_coefficients(r, q)
    automatic = driver.tangent_coefficients_autodiff(r, q)
    for left, right in zip(analytic, automatic, strict=True):
        torch.testing.assert_close(left, right, rtol=1e-12, atol=1e-12)
    c, b, h = analytic
    reconstructed = c * r + h
    if q is not None:
        reconstructed = reconstructed + (b * q).sum(dim=-1)
    torch.testing.assert_close(reconstructed, driver.driver(r, q), rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("name", ["allen_cahn", "joint_concave", "joint_convex"])
def test_tangent_is_upper_bound_on_valid_domain(name: str) -> None:
    torch.manual_seed(11)
    driver = make_driver(name, d=4)
    y = driver.r_low + (driver.r_high - driver.r_low) * torch.rand(101, dtype=torch.float64)
    r = driver.r_low + (driver.r_high - driver.r_low) * torch.rand(101, dtype=torch.float64)
    z = torch.randn(101, driver.d, dtype=torch.float64) if driver.has_gradient else None
    q = torch.randn(101, driver.d, dtype=torch.float64) if driver.has_gradient else None
    c, b, h = driver.tangent_coefficients(r, q)
    tangent = c * y + h
    if z is not None:
        tangent = tangent + (b * z).sum(dim=-1)
    assert float(torch.min(tangent - driver.driver(y, z))) >= -2e-12
    divergence = driver.bregman_divergence(y, z, r, q)
    torch.testing.assert_close(divergence, tangent - driver.driver(y, z))


def test_sign_transform_terminal_driver_and_reconstruction() -> None:
    torch.manual_seed(12)
    original = JointConvexOriginalDriver(d=7)
    transformed = SignTransformedDriver(original)
    x = torch.randn(13, 7, dtype=torch.float64)
    w = -torch.rand(13, dtype=torch.float64)
    z_w = torch.randn(13, 7, dtype=torch.float64)
    torch.testing.assert_close(transformed.terminal(x), -original.terminal(x))
    torch.testing.assert_close(
        transformed.driver(w, z_w), -original.driver(-w, -z_w)
    )
    c, b, h = transformed.tangent_coefficients(w, z_w)
    torch.testing.assert_close(c, 3.0 * w.square() - 1.0)
    torch.testing.assert_close(b, -z_w)
    torch.testing.assert_close(h, -2.0 * w.pow(3) + 0.5 * z_w.square().sum(dim=-1))
    values = np.linspace(-1.0, 0.0, 9)
    np.testing.assert_allclose(transformed.working_to_original(values), -values)


class _LinearDriver(ConcaveSemilinearDriver):
    name, d, T, has_gradient = "linear_test", 1, 0.7, False
    r_low, r_high = -2.0, 2.0
    c_value, h_value, g_value = -0.8, 0.35, 1.2

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        return torch.full(x.shape[:-1], self.g_value, dtype=x.dtype, device=x.device)

    def terminal_numpy(self, x: np.ndarray) -> np.ndarray:
        return np.full(x.shape[:-1], self.g_value)

    def driver(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        del z
        return self.c_value * y + self.h_value

    def driver_y(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        del z
        return torch.full_like(y, self.c_value)

    def driver_numpy(self, y: np.ndarray, z_scalar: np.ndarray) -> np.ndarray:
        del z_scalar
        return self.c_value * y + self.h_value


def test_weighted_feynman_kac_constant_coefficients() -> None:
    driver = _LinearDriver()
    actor = ConstantTangentActor(driver)
    paths, steps = 10, 35
    t = torch.tensor([0.0, 0.2], dtype=torch.float64)
    x = torch.zeros(2, 1, dtype=torch.float64)
    normals = torch.zeros(2, paths, steps, 1, dtype=torch.float64)
    payoff, _, weight, _ = _weighted_rollout_chunk(
        driver, actor, t, x, normals, need_z=False
    )
    durations = driver.T - t
    dt = durations / steps
    ratio = torch.exp(driver.c_value * dt)
    discrete_source = driver.h_value * dt * (1.0 - ratio.pow(steps)) / (1.0 - ratio)
    expected = torch.exp(driver.c_value * durations) * driver.g_value + discrete_source
    torch.testing.assert_close(payoff, expected[:, None].expand_as(payoff), rtol=2e-12, atol=2e-12)
    torch.testing.assert_close(
        weight, torch.exp(driver.c_value * durations)[:, None].expand_as(weight),
        rtol=2e-12, atol=2e-12,
    )
    continuous = (
        torch.exp(driver.c_value * durations) * driver.g_value
        + driver.h_value * torch.expm1(driver.c_value * durations) / driver.c_value
    )
    assert float(torch.max(torch.abs(expected - continuous))) < 0.003


class _StateActor(nn.Module):
    def __init__(self, driver: ConcaveSemilinearDriver) -> None:
        super().__init__()
        self.driver = driver
        self.linear = nn.Linear(driver.d + 1, driver.d + 1, bias=True, dtype=torch.float64)
        torch.manual_seed(19)
        with torch.no_grad():
            self.linear.weight.normal_(0.0, 0.08)
            self.linear.bias.zero_()

    def predict(self, t: torch.Tensor, x: torch.Tensor, T: float):
        raw = self.linear(torch.cat((t[:, None] / T, x), dim=-1))
        return self.driver.constrain_r(raw[:, 0]), raw[:, 1:]


def test_full_vjp_matches_common_random_number_finite_difference() -> None:
    driver = JointConcaveDriver(d=3, T=0.35)
    actor = _StateActor(driver)
    generator = torch.Generator().manual_seed(339)
    normals = torch.randn(2, 96, 10, 3, generator=generator, dtype=torch.float64)
    t = torch.tensor([0.03, 0.11], dtype=torch.float64)
    x = torch.tensor([[0.2, -0.1, 0.3], [-0.4, 0.2, 0.1]], dtype=torch.float64)
    with frozen(actor):
        payoff, z, _, _ = _weighted_rollout_chunk(
            driver, actor, t, x, normals, need_z=True
        )
    assert z is not None
    assert all(parameter.grad is None for parameter in actor.parameters())
    assert all(parameter.requires_grad for parameter in actor.parameters())
    epsilon = 2e-5
    finite = torch.empty_like(x)
    for coordinate in range(driver.d):
        shift = torch.zeros_like(x)
        shift[:, coordinate] = epsilon
        plus, _, _, _ = _weighted_rollout_chunk(
            driver, actor, t, x + shift, normals, need_z=False
        )
        minus, _, _, _ = _weighted_rollout_chunk(
            driver, actor, t, x - shift, normals, need_z=False
        )
        finite[:, coordinate] = driver.sigma * (
            plus.mean(dim=1) - minus.mean(dim=1)
        ) / (2.0 * epsilon)
    torch.testing.assert_close(z, finite, rtol=3e-4, atol=3e-4)


def test_terminal_condition_exactly_recovered() -> None:
    for name in ("allen_cahn", "joint_concave", "joint_convex"):
        driver = make_driver(name, d=6)
        actor = ConstantTangentActor(driver)
        x = torch.randn(4, 6)
        t = torch.full((4,), driver.T)
        normals = torch.randn(4, 8, 5, 6)
        payoff, _, weight, _ = _weighted_rollout_chunk(
            driver, actor, t, x, normals, need_z=False
        )
        torch.testing.assert_close(payoff, driver.terminal(x)[:, None].expand_as(payoff))
        torch.testing.assert_close(weight, torch.ones_like(weight))


class _QuadraticGradientDriver(ConcaveSemilinearDriver):
    name, d, T, has_gradient = "quadratic_hjb_test", 2, 0.4, True
    r_low, r_high = -5.0, 5.0

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        return 0.2 * x.square().sum(dim=-1)

    def terminal_numpy(self, x: np.ndarray) -> np.ndarray:
        return 0.2 * np.square(x).sum(axis=-1)

    def driver(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        del y
        assert z is not None
        return -0.5 * z.square().sum(dim=-1)

    def driver_y(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        del z
        return torch.zeros_like(y)

    def driver_z(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        del y
        assert z is not None
        return -z

    def driver_numpy(self, y: np.ndarray, z_scalar: np.ndarray) -> np.ndarray:
        del y
        return -0.5 * np.square(z_scalar)


class _ConstantJointActor(nn.Module):
    def __init__(self, q: torch.Tensor) -> None:
        super().__init__()
        self.register_buffer("q", q)

    def predict(self, t: torch.Tensor, x: torch.Tensor, T: float):
        del T
        r = torch.zeros_like(t)
        return r, self.q.to(x).expand_as(x)


def test_quadratic_hjb_reduces_to_established_controlled_policy_formula() -> None:
    driver = _QuadraticGradientDriver()
    q = torch.tensor([0.3, -0.2], dtype=torch.float64)
    actor = _ConstantJointActor(q)
    t = torch.tensor([0.07], dtype=torch.float64)
    x = torch.tensor([[0.1, 0.4]], dtype=torch.float64)
    normals = torch.randn(1, 14, 9, 2, generator=torch.Generator().manual_seed(40), dtype=torch.float64)
    payoff, _, weight, _ = _weighted_rollout_chunk(
        driver, actor, t, x, normals, need_z=False
    )
    dt = (driver.T - t[0]) / normals.shape[2]
    terminal_state = (
        x[:, None, :]
        - driver.sigma * q[None, None, :] * (driver.T - t[0])
        + driver.sigma * math.sqrt(float(dt)) * normals.sum(dim=2)
    )
    expected = driver.terminal(terminal_state) + 0.5 * q.square().sum() * (driver.T - t[0])
    torch.testing.assert_close(payoff, expected, rtol=2e-12, atol=2e-12)
    torch.testing.assert_close(weight, torch.ones_like(weight))


@pytest.mark.parametrize("name", ["allen_cahn", "joint_concave", "joint_convex"])
def test_periodic_reference_refinement_residual_and_domain(name: str) -> None:
    _, diagnostics = validate_periodic_reference(make_driver(name, d=5))
    # Linear interpolation on a 512 point periodic grid sets this floor.
    assert diagnostics["terminal_max_abs"] < 5e-6
    assert diagnostics["value_refinement_rmse"] < 2e-5
    assert diagnostics["gradient_refinement_rmse"] < 2e-4
    assert diagnostics["pde_residual_rmse"] < 2e-3
    assert diagnostics["inside_tangent_domain"]


def test_training_core_has_no_reference_dependency() -> None:
    source = inspect.getsource(core)
    assert "from .reference" not in source
    assert "import reference" not in source
