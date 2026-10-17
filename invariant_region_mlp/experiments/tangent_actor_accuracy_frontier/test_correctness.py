"""CPU correctness gates for the nonconvex Tangent-Actor experiment."""

from __future__ import annotations

import inspect
import math

import numpy as np
import pytest
import torch

from invariant_region_mlp.experiments.tangent_actor_mcv_training.labels import (
    _weighted_rollout_chunk_mcv_labels,
)
from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis.models import (
    TerminalAwareActor,
)

from .drivers import SineQuadraticDriver
from invariant_region_mlp.experiments.nonconvex_tangent_actor.full_history_baseline import (
    RawDriver,
)
from invariant_region_mlp.experiments.nonconvex_tangent_actor.pinn_baseline import (
    PlainPINN,
    pde_residual,
)
from .reference2d import solve_strang_reference


@pytest.mark.parametrize("eta", [0.0, 1.0])
def test_analytic_first_derivatives_and_tangent_match_autograd(eta: float) -> None:
    driver = SineQuadraticDriver(d=5, lam=0.2, eta=eta)
    generator = torch.Generator().manual_seed(19)
    y = torch.rand(7, generator=generator, dtype=torch.float64, requires_grad=True)
    z = torch.randn(7, 5, generator=generator, dtype=torch.float64, requires_grad=True)
    value = driver.driver(y, z)
    auto_y, auto_z = torch.autograd.grad(value.sum(), (y, z))
    torch.testing.assert_close(driver.driver_y(y, z), auto_y, rtol=1e-11, atol=1e-11)
    torch.testing.assert_close(driver.driver_z(y, z), auto_z, rtol=1e-11, atol=1e-11)
    analytic = driver.tangent_coefficients(y.detach(), z.detach())
    autodiff = driver.tangent_coefficients_autodiff(y.detach(), z.detach())
    for left, right in zip(analytic, autodiff, strict=True):
        torch.testing.assert_close(left, right, rtol=1e-11, atol=1e-11)


@pytest.mark.parametrize("eta", [0.0, 1.0])
def test_full_joint_hessian_matches_autograd(eta: float) -> None:
    driver = SineQuadraticDriver(d=4, lam=0.2, eta=eta)
    point = torch.tensor([0.63, 0.2, -0.1, 0.4, 0.7], dtype=torch.float64)

    def scalar(argument: torch.Tensor) -> torch.Tensor:
        return driver.driver(argument[0:1], argument[1:].reshape(1, -1))[0]

    actual = torch.autograd.functional.hessian(scalar, point)
    expected = torch.zeros_like(actual)
    expected[0, 0] = driver.driver_uu(point[0:1])[0]
    expected[1:, 1:] = -torch.eye(4, dtype=torch.float64)
    expected[0, 1:] = expected[1:, 0] = -eta * torch.as_tensor(
        driver.e, dtype=torch.float64
    )
    torch.testing.assert_close(actual, expected, rtol=1e-10, atol=1e-10)


def test_formal_curvature_changes_sign_and_coupling_is_nonzero() -> None:
    separable = SineQuadraticDriver(d=4, lam=0.2, eta=0.0)
    y = torch.tensor([0.25, 0.75], dtype=torch.float64)
    curvature = separable.driver_uu(y)
    assert curvature[0] < 0.0 < curvature[1]
    coupled = SineQuadraticDriver(d=4, lam=0.2, eta=1.0)
    eigenvalues = coupled.reduced_joint_hessian_eigenvalues(np.asarray([0.75]))[0]
    assert eigenvalues[0] < 0.0 < eigenvalues[-1]


@pytest.mark.parametrize("eta", [0.0, 1.0])
def test_zero_one_barriers_and_terminal_range(eta: float) -> None:
    driver = SineQuadraticDriver(d=8, lam=0.2, eta=eta)
    audit = driver.barrier_audit()
    assert audit["invariant_interval_verified_analytically"]
    rng = np.random.default_rng(31)
    terminal = driver.terminal_numpy(rng.standard_normal((4096, driver.d)))
    assert terminal.min() >= 0.2 - 1e-12
    assert terminal.max() <= 0.8 + 1e-12


def test_terminal_aware_actor_is_exact_and_uses_full_ambient_input() -> None:
    driver = SineQuadraticDriver(d=7, lam=0.2)
    actor = TerminalAwareActor(driver, width=32, depth=2)
    assert actor.network[0].in_features == driver.d + 1
    x = torch.randn(11, driver.d)
    t = torch.full((11,), driver.T)
    r, q = actor.predict(t, x, driver.T)
    torch.testing.assert_close(r, driver.terminal(x), atol=2e-7, rtol=2e-7)
    assert q is not None
    torch.testing.assert_close(
        q, math.sqrt(2.0) * driver.terminal_gradient_torch(x), atol=2e-7, rtol=2e-7
    )


def test_small_periodic_reference_is_finite_and_inside_domain() -> None:
    driver = SineQuadraticDriver(d=4, lam=0.2)
    reference = solve_strang_reference(driver, n=32, dt=0.001, store_dt=0.01)
    assert np.all(np.isfinite(reference.values))
    assert reference.values.min() >= -1e-8
    assert reference.values.max() <= 1.0 + 1e-8
    rng = np.random.default_rng(41)
    coordinates = rng.uniform(0.0, 2.0 * np.pi, (128, 2))
    terminal = reference.value(np.full(128, driver.T), coordinates)
    exact = driver.terminal_from_coordinates_numpy(coordinates)
    assert np.sqrt(np.mean(np.square(terminal - exact))) < 3e-3


def _corrected_mean(
    driver: SineQuadraticDriver,
    actor: TerminalAwareActor,
    t: torch.Tensor,
    x: torch.Tensor,
    normals: torch.Tensor,
    *,
    need_z: bool,
) -> tuple[torch.Tensor, torch.Tensor | None, torch.Tensor]:
    corrected, z, _, martingale, _ = _weighted_rollout_chunk_mcv_labels(
        driver, actor, t, x, normals, need_z=need_z
    )
    return corrected.mean(dim=1), z, martingale


def test_mcv_mean_identity_and_full_path_gradient_common_random_difference() -> None:
    torch.manual_seed(7)
    driver = SineQuadraticDriver(d=4, lam=0.2, eta=1.0)
    actor = TerminalAwareActor(driver, width=32, depth=2)
    actor.eval()
    paths, steps = 2048, 8
    generator = torch.Generator().manual_seed(811)
    half = torch.randn((1, paths // 2, steps, driver.d), generator=generator)
    normals = torch.cat((half, -half), dim=1)
    t = torch.tensor([0.17])
    x = torch.tensor([[0.11, -0.23, 0.07, 0.19]])
    corrected_mean, z, martingale = _corrected_mean(
        driver, actor, t, x, normals, need_z=True
    )
    assert z is not None
    pair = 0.5 * (martingale[:, : paths // 2] + martingale[:, paths // 2 :])
    pair_mean = float(pair.mean())
    pair_se = float(pair.std(unbiased=True) / math.sqrt(pair.numel()))
    assert abs(pair_mean) <= 5.0 * max(pair_se, 1e-8)

    epsilon = 2e-3
    plus, minus = x.clone(), x.clone()
    plus[0, 0] += epsilon
    minus[0, 0] -= epsilon
    plus_mean, _, _ = _corrected_mean(
        driver, actor, t, plus, normals, need_z=False
    )
    minus_mean, _, _ = _corrected_mean(
        driver, actor, t, minus, normals, need_z=False
    )
    finite_difference_z = math.sqrt(2.0) * (plus_mean - minus_mean) / (2.0 * epsilon)
    relative = float(
        torch.abs(finite_difference_z - z[:, 0])
        / torch.clamp(torch.abs(z[:, 0]), min=1e-3)
    )
    assert relative < 0.08
    assert torch.isfinite(corrected_mean).all()


def test_training_dependencies_do_not_import_reference_or_active_coordinates() -> None:
    from invariant_region_mlp.experiments.tangent_actor_mcv_training import labels
    from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis import training

    source = inspect.getsource(labels) + inspect.getsource(training)
    assert "reference2d" not in source
    assert "@ driver.B" not in source
    assert "@driver.B" not in source


def test_full_history_adapter_does_not_project_recursive_values() -> None:
    raw = RawDriver(SineQuadraticDriver(d=4, lam=0.2))
    assert raw.r_low == -float("inf")
    assert raw.r_high == float("inf")
    y = torch.tensor([1.2], dtype=torch.float64)
    z = torch.zeros(1, 4, dtype=torch.float64)
    torch.testing.assert_close(raw.driver(y, z), raw.original.driver(y, z))


def test_plain_pinn_residual_supports_the_generic_joint_driver() -> None:
    driver = SineQuadraticDriver(d=4, lam=0.2, eta=1.0)
    model = PlainPINN(driver.d, width=16, depth=2)
    generator = torch.Generator().manual_seed(97)
    residual = pde_residual(
        model,
        driver,
        torch.rand(6, generator=generator) * driver.T,
        torch.randn(6, driver.d, generator=generator),
        probes=2,
        generator=generator,
    )
    assert residual.shape == (6,)
    assert torch.isfinite(residual).all()
