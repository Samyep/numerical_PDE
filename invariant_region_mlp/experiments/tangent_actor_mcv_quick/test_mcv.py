"""Minimal correctness tests for the inference-only martingale correction."""

from __future__ import annotations

import math

import numpy as np
import torch
from torch import nn

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.core import (
    evaluate_tangent_mc,
)
from invariant_region_mlp.experiments.convex_concave_semilinear_extended.drivers import (
    JointConcaveDriver,
)
from invariant_region_mlp.experiments.convex_concave_semilinear_extended.models import (
    ConstantTangentActor,
)

from .mcv import (
    _weighted_rollout_chunk_mcv,
    antithetic_pair_average,
    calibrated_alpha,
    variance_of_antithetic_mean,
)


class _ZeroDriver:
    d = 1
    T = 1.0
    sigma = math.sqrt(2.0)

    @staticmethod
    def tangent_coefficients(r: torch.Tensor, q: torch.Tensor):
        return torch.zeros_like(r), torch.zeros_like(q), torch.zeros_like(r)

    @staticmethod
    def terminal(x: torch.Tensor) -> torch.Tensor:
        return torch.zeros(x.shape[:-1], dtype=x.dtype, device=x.device)


class _StateGradientActor(nn.Module):
    def predict(self, t: torch.Tensor, x: torch.Tensor, T: float):
        del t, T
        return torch.zeros(len(x), dtype=x.dtype, device=x.device), x.clone()


class _ZeroGradientActor(nn.Module):
    def predict(self, t: torch.Tensor, x: torch.Tensor, T: float):
        del t, T
        return torch.zeros(len(x), dtype=x.dtype, device=x.device), torch.zeros_like(x)


def test_martingale_uses_left_endpoint_dw_without_extra_sqrt2() -> None:
    driver = _ZeroDriver()
    actor = _StateGradientActor()
    t = torch.tensor([0.0])
    x = torch.tensor([[2.0]])
    normals = torch.tensor([[[[0.25]], [[-0.50]]]])
    _, martingale, _, _ = _weighted_rollout_chunk_mcv(
        driver, actor, t, x, normals, collect_martingale=True
    )
    assert martingale is not None
    # dt=1, so dW=normal and q(X_0)=2. No sigma=sqrt(2) factor belongs here.
    torch.testing.assert_close(martingale, torch.tensor([[0.5, -1.0]]))


def test_zero_q_mcv_is_exactly_the_plain_estimator() -> None:
    driver = _ZeroDriver()
    actor = _ZeroGradientActor()
    t = torch.tensor([0.0])
    x = torch.tensor([[1.25]])
    normals = torch.tensor([[[[0.25]], [[-0.25]]]])
    payoff, martingale, _, _ = _weighted_rollout_chunk_mcv(
        driver, actor, t, x, normals, collect_martingale=True
    )
    assert martingale is not None
    torch.testing.assert_close(martingale, torch.zeros_like(martingale), rtol=0.0, atol=0.0)
    torch.testing.assert_close(payoff - martingale, payoff, rtol=0.0, atol=0.0)


def test_antithetic_variance_uses_pair_averages_as_independent_units() -> None:
    paths = torch.tensor([[1.0, 3.0, -1.0, 1.0]])
    pairs = antithetic_pair_average(paths)
    torch.testing.assert_close(pairs, torch.tensor([[0.0, 2.0]]))
    # Sample variance of [0,2] is 2; variance of their mean is 2/2=1.
    torch.testing.assert_close(
        variance_of_antithetic_mean(paths), torch.tensor([1.0], dtype=torch.float64)
    )


def test_calibration_centers_within_state_and_handles_zero_denominator() -> None:
    h = np.asarray([[1.0, -1.0, 2.0, -2.0], [3.0, -3.0, 4.0, -4.0]])
    state_means = np.asarray([[10.0], [-7.0]])
    payoff = state_means + 2.0 * h
    alpha, diagnostics = calibrated_alpha(payoff, h)
    assert abs(alpha - 2.0) < 1e-12
    np.testing.assert_allclose(
        payoff - alpha * h, np.broadcast_to(state_means, payoff.shape)
    )
    assert not diagnostics["denominator_was_zero"]
    zero_alpha, zero_diagnostics = calibrated_alpha(payoff, np.zeros_like(h))
    assert zero_alpha == 0.0
    assert zero_diagnostics["denominator_was_zero"]


def test_existing_fk_metadata_confirms_no_equivalent_mcv() -> None:
    driver = JointConcaveDriver(d=2, k=2, T=0.5)
    actor = ConstantTangentActor(driver)
    output, metadata = evaluate_tangent_mc(
        driver,
        actor,
        np.zeros(2),
        np.zeros((2, 2)),
        paths=4,
        time_steps=2,
        seed=123,
        device=torch.device("cpu"),
        parent_batch=2,
        path_chunk=4,
    )
    assert not metadata["martingale_control_variate"]
    assert metadata["antithetic"]
    assert np.all(np.isfinite(output["working_value"]))
