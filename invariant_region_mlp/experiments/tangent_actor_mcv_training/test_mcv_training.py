"""Correctness checks for differentiable MCV supervision."""

from __future__ import annotations

import math

import numpy as np
import torch
from torch import nn

from .labels import _weighted_rollout_chunk_mcv_labels, mcv_tangent_labels


class _ZeroDriver:
    d = 1
    T = 1.0
    sigma = math.sqrt(2.0)
    has_gradient = True
    r_low = -10.0
    r_high = 10.0

    @staticmethod
    def tangent_coefficients(r: torch.Tensor, q: torch.Tensor):
        return torch.zeros_like(r), torch.zeros_like(q), torch.zeros_like(r)

    @staticmethod
    def terminal(x: torch.Tensor) -> torch.Tensor:
        return x[..., 0].square()


class _LinearActor(nn.Module):
    def __init__(self, slope: float = 1.0) -> None:
        super().__init__()
        self.slope = nn.Parameter(torch.tensor(float(slope)))

    def predict(self, t: torch.Tensor, x: torch.Tensor, T: float):
        del t, T
        return torch.zeros(len(x), dtype=x.dtype, device=x.device), self.slope * x


class _ZeroActor(nn.Module):
    def predict(self, t: torch.Tensor, x: torch.Tensor, T: float):
        del t, T
        return (
            torch.zeros(len(x), dtype=x.dtype, device=x.device),
            torch.zeros_like(x),
        )


def test_left_endpoint_and_complete_correction_gradient() -> None:
    driver = _ZeroDriver()
    actor = _LinearActor()
    t = torch.tensor([0.0])
    x = torch.tensor([[2.0]])
    normals = torch.tensor([[[[0.25]], [[-0.50]]]])
    payoff, z, _, correction, metadata = _weighted_rollout_chunk_mcv_labels(
        driver, actor, t, x, normals, need_z=True
    )
    assert z is not None
    # Left-point q is x=2 and dt=1, so H is [0.5, -1].
    torch.testing.assert_close(correction, torch.tensor([[0.5, -1.0]]))
    # The differentiated value is mean[(x+sqrt(2)dW)^2-x*dW].
    expected_derivative = np.mean([
        2.0 * (2.0 + math.sqrt(2.0) * 0.25) - 0.25,
        2.0 * (2.0 - math.sqrt(2.0) * 0.50) + 0.50,
    ])
    torch.testing.assert_close(
        z / driver.sigma,
        torch.tensor([[expected_derivative]], dtype=z.dtype),
        rtol=2e-6,
        atol=2e-6,
    )
    assert metadata["nonfinite_count"] == 0.0
    assert torch.all(torch.isfinite(payoff))


def test_zero_q_matches_plain_value_and_gradient_exactly() -> None:
    driver = _ZeroDriver()
    actor = _ZeroActor()
    t = torch.tensor([0.0])
    x = torch.tensor([[1.25]])
    normals = torch.tensor([[[[0.2]], [[-0.2]]]])
    corrected, corrected_z, _, correction, _ = _weighted_rollout_chunk_mcv_labels(
        driver, actor, t, x, normals, need_z=True
    )
    assert corrected_z is not None
    torch.testing.assert_close(correction, torch.zeros_like(correction), rtol=0, atol=0)
    terminal = (x[:, None, 0] + driver.sigma * normals[..., 0, 0]).square()
    torch.testing.assert_close(corrected, terminal)
    expected_z = driver.sigma * (2.0 * x[:, 0])
    torch.testing.assert_close(corrected_z[:, 0], expected_z)


def test_source_parameters_frozen_restored_and_receive_no_grad() -> None:
    driver = _ZeroDriver()
    actor = _LinearActor()
    assert actor.slope.requires_grad
    values, z, metadata = mcv_tangent_labels(
        driver,
        actor,
        np.asarray([0.1, 0.3]),
        np.asarray([[0.4], [-0.2]]),
        paths=8,
        time_steps=3,
        seed=3,
        device=torch.device("cpu"),
        parent_batch=2,
        path_chunk=4,
    )
    assert actor.slope.requires_grad
    assert actor.slope.grad is None
    assert np.all(np.isfinite(values)) and z is not None and np.all(np.isfinite(z))
    assert metadata["source_actor_frozen"]
    assert metadata["state_dependence_preserved"]
    assert metadata["complete_corrected_payoff_differentiated"]
    assert not metadata["mcv_q_detached"]


def test_autodiff_matches_crn_directional_finite_difference() -> None:
    driver = _ZeroDriver()
    actor = _LinearActor(0.7)
    t = np.asarray([0.2, 0.4])
    x = np.asarray([[0.7], [-0.5]])
    common = dict(
        paths=64,
        time_steps=4,
        seed=123,
        device=torch.device("cpu"),
        parent_batch=2,
        path_chunk=16,
    )
    _, z, _ = mcv_tangent_labels(driver, actor, t, x, need_z=True, **common)
    epsilon = 1e-3
    plus, _, _ = mcv_tangent_labels(
        driver, actor, t, x + epsilon, need_z=False, **common
    )
    minus, _, _ = mcv_tangent_labels(
        driver, actor, t, x - epsilon, need_z=False, **common
    )
    assert z is not None
    fd = (plus - minus) / (2.0 * epsilon)
    np.testing.assert_allclose(z[:, 0] / driver.sigma, fd, rtol=2e-3, atol=2e-3)
