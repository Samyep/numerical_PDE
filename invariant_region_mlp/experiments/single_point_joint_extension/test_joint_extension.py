"""Fast unit tests for the frozen joint-PDE extension."""

from __future__ import annotations

import numpy as np
import torch

from .benchmark import frozen_basis, make_benchmark
from .protocol import DIMENSIONS, QUERY_COUNT, frozen_protocol, load_query_bank


def test_basis_and_query_banks_are_frozen_and_well_formed() -> None:
    for dimension in DIMENSIONS:
        basis = frozen_basis(dimension)
        assert np.max(np.abs(basis.T @ basis - np.eye(2))) <= 1e-12
        t, x, truth, truth_z = load_query_bank(dimension)
        assert t.shape == (QUERY_COUNT,)
        assert x.shape == (QUERY_COUNT, dimension)
        assert truth.shape == (QUERY_COUNT,)
        assert truth_z.shape == (QUERY_COUNT, dimension)
        assert np.array_equal(x[0], np.zeros(dimension))
        assert np.all(np.isfinite(truth))


def test_terminal_and_tangent_coefficients_match_autograd() -> None:
    driver = make_benchmark(20)
    generator = torch.Generator().manual_seed(19)
    r = torch.rand(32, dtype=torch.float64, generator=generator)
    q = torch.randn((32, 20), dtype=torch.float64, generator=generator)
    analytic = driver.tangent_coefficients(r, q)
    autodiff = driver.tangent_coefficients_autodiff(r, q)
    for left, right in zip(analytic, autodiff):
        assert torch.max(torch.abs(left - right)).item() <= 1e-11


def test_genuinely_nonseparable_reduction() -> None:
    driver = make_benchmark(50)
    y = np.linspace(0.0, 1.0, 11)
    reduced_z = np.stack((np.linspace(-1.0, 1.0, 11), np.ones(11)), axis=1)
    ambient_z = reduced_z @ driver.B.T
    ambient = driver.driver(
        torch.as_tensor(y), torch.as_tensor(ambient_z)
    ).numpy()
    reduced = driver.driver_reduced_numpy(y, reduced_z)
    assert np.max(np.abs(ambient - reduced)) <= 1e-12
    assert np.linalg.norm(driver.v) > 0.999999999999


def test_protocol_identifies_frozen_preregistration() -> None:
    protocol = frozen_protocol()
    assert protocol["preregistration_commit"] == "0f1485703"
    assert protocol["query_count"] == 32
    assert protocol["pde"] == "joint_nonseparable"
