"""Numerical and structural checks required before allocating the GPU pilot."""

from __future__ import annotations

import inspect

import numpy as np
import pytest
import torch

from . import core
from .core import evaluate_tangent_mc
from .drivers import NonseparableDriver, make_driver, terminal_signal_diagnostics
from .models import ConstantTangentActor, TangentActor
from .reference2d import validate_periodic_reference_2d


@pytest.mark.parametrize(
    "family",
    ["allen_cahn", "joint_concave", "joint_convex", "nonseparable"],
)
def test_generic_tangent_coefficients_match_autodiff_and_support(family: str) -> None:
    driver = make_driver(family, d=8, k=2)
    generator = torch.Generator().manual_seed(17)
    r = driver.r_low + (driver.r_high - driver.r_low) * torch.rand(
        19, generator=generator, dtype=torch.float64
    )
    q = (
        torch.randn(19, driver.d, generator=generator, dtype=torch.float64) * 0.3
        if driver.has_gradient
        else None
    )
    c, b, h = driver.tangent_coefficients(r, q)
    ac, ab, ah = driver.tangent_coefficients_autodiff(r, q)
    assert torch.allclose(c, ac, atol=2e-12, rtol=2e-12)
    assert torch.allclose(b, ab, atol=2e-12, rtol=2e-12)
    assert torch.allclose(h, ah, atol=2e-12, rtol=2e-12)
    assert torch.allclose(c * r + h + (0.0 if q is None else (b * q).sum(-1)), driver.driver(r, q))

    y = driver.r_low + (driver.r_high - driver.r_low) * torch.rand(
        19, generator=generator, dtype=torch.float64
    )
    z = (
        torch.randn(19, driver.d, generator=generator, dtype=torch.float64) * 0.3
        if driver.has_gradient
        else None
    )
    tangent = c * y + h + (0.0 if z is None else (b * z).sum(-1))
    assert float((driver.driver(y, z) - tangent).max()) < 2e-12


def test_nonseparable_mixed_derivative_and_joint_concavity() -> None:
    driver = NonseparableDriver(d=5)
    for value in (0.0, 0.2, 0.7, 1.0):
        point = torch.zeros(6, dtype=torch.float64, requires_grad=True)
        point.data[0] = value

        def scalar(argument: torch.Tensor) -> torch.Tensor:
            return driver.driver(argument[0:1], argument[1:].reshape(1, -1))[0]

        hessian = torch.autograd.functional.hessian(scalar, point).detach().numpy()
        assert np.linalg.eigvalsh(hessian).max() < 2e-10
    y = torch.tensor([0.43], dtype=torch.float64, requires_grad=True)
    z = torch.randn(1, 5, dtype=torch.float64, requires_grad=True)
    mixed = torch.autograd.grad(driver.driver_y(y, z).sum(), z)[0]
    assert torch.allclose(mixed[0], -torch.as_tensor(driver.v, dtype=torch.float64), atol=2e-12)


@pytest.mark.parametrize("k", [2, 5, 10, 20])
def test_dense_basis_and_terminal_gradient_span(k: int) -> None:
    driver = make_driver("joint_concave", d=100, k=k)
    assert np.max(np.abs(driver.B.T @ driver.B - np.eye(k))) < 2e-14
    assert np.median(np.abs(driver.B)) > 5e-3
    diagnostics = terminal_signal_diagnostics(driver, samples=2048)
    assert diagnostics["ambient_gradient_matrix_rank"] == k
    assert diagnostics["centered_ambient_gradient_matrix_rank"] == k
    assert diagnostics["centered_smallest_to_largest_active_singular_ratio"] > 0.8
    assert diagnostics["expected_rank"] == k


def test_terminal_amplitude_scaling_matches_controlled_formula() -> None:
    variances, energies = [], []
    for k in (2, 5, 10, 20):
        item = terminal_signal_diagnostics(make_driver("allen_cahn", d=100, k=k), samples=32768)
        variances.append(item["terminal_variance"])
        energies.append(item["gradient_energy"])
    ratios_v = np.asarray(variances[:-1]) / np.asarray(variances[1:])
    ratios_e = np.asarray(energies[:-1]) / np.asarray(energies[1:])
    expected = np.asarray([5 / 2, 10 / 5, 20 / 10])
    assert np.allclose(ratios_v, expected, rtol=0.12)
    assert np.allclose(ratios_e, expected, rtol=0.12)


def test_actor_uses_full_ambient_state() -> None:
    driver = make_driver("joint_concave", d=100, k=5)
    actor = TangentActor(driver, width=32, depth=2)
    first = next(layer for layer in actor.network if isinstance(layer, torch.nn.Linear))
    last = [layer for layer in actor.network if isinstance(layer, torch.nn.Linear)][-1]
    assert first.in_features == 101
    assert last.out_features == 101
    source = inspect.getsource(TangentActor.predict)
    assert "@" not in source and ".B" not in source and "basis" not in source


def test_terminal_recovery_of_generic_evaluator() -> None:
    driver = make_driver("nonseparable", d=8, k=2)
    rng = np.random.default_rng(42)
    x = rng.standard_normal((12, driver.d))
    output, metadata = evaluate_tangent_mc(
        driver,
        ConstantTangentActor(driver),
        np.full(12, driver.T),
        x,
        paths=8,
        time_steps=4,
        seed=3,
        device=torch.device("cpu"),
        parent_batch=4,
        path_chunk=8,
    )
    assert np.max(np.abs(output["working_value"] - driver.terminal_numpy(x))) < 2e-7
    assert metadata["nonfinite_count"] == 0


@pytest.mark.parametrize(
    "family", ["allen_cahn", "joint_concave", "joint_convex", "nonseparable"]
)
def test_two_dimensional_reference_refinement_residual_and_domain(family: str) -> None:
    _, diagnostics = validate_periodic_reference_2d(make_driver(family, d=5, k=2))
    assert diagnostics["terminal_interpolation_max_abs"] < 2e-4
    assert diagnostics["value_refinement_rmse"] < 2e-4
    assert diagnostics["gradient_refinement_vector_rmse"] < 2e-3
    assert diagnostics["pde_residual_rmse"] < 3e-3
    assert diagnostics["inside_tangent_domain"]


def test_training_core_has_no_reference_or_projection_dependency() -> None:
    source = inspect.getsource(core)
    assert "from .reference" not in source
    assert "import reference" not in source
    assert ".B" not in source and "basis" not in source.lower()
