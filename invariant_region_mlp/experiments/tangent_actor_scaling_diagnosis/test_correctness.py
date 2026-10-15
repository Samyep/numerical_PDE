"""Correctness tests for the controlled scaling diagnosis."""

from __future__ import annotations

import inspect

import numpy as np
import pytest
import torch

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.core import (
    tangent_labels,
)
from invariant_region_mlp.experiments.convex_concave_semilinear_extended.drivers import (
    make_driver,
)

from .common import nested_points
from .models import StandardActor, TerminalAwareActor


@pytest.mark.parametrize(
    "family", ["allen_cahn", "joint_concave", "joint_convex", "nonseparable"]
)
def test_terminal_aware_actor_enforces_terminal_value_and_gradient(family: str) -> None:
    driver = make_driver(family, d=8, k=2)
    actor = TerminalAwareActor(driver, width=16, depth=2).double()
    generator = torch.Generator().manual_seed(113)
    x = torch.randn(19, driver.d, generator=generator, dtype=torch.float64)
    t = torch.full((len(x),), driver.T, dtype=torch.float64)
    r, q = actor.predict(t, x, driver.T)
    assert torch.allclose(r, driver.terminal(x), atol=2e-14, rtol=2e-14)
    assert float(r.detach().min()) > driver.r_low
    assert float(r.detach().max()) < driver.r_high
    if driver.has_gradient:
        assert q is not None
        expected = driver.sigma * driver.terminal_gradient_torch(x)
        assert torch.allclose(q, expected, atol=2e-14, rtol=2e-14)
    else:
        assert q is None


def test_transformed_convex_terminal_parameterization_is_negative_and_exact() -> None:
    driver = make_driver("joint_convex", d=8, k=2)
    actor = TerminalAwareActor(driver, width=12, depth=1).double()
    x = torch.randn(11, driver.d, dtype=torch.float64)
    t = torch.full((len(x),), driver.T, dtype=torch.float64)
    r, _ = actor.predict(t, x, driver.T)
    assert float(r.detach().min()) > -1.0
    assert float(r.detach().max()) < 0.0
    assert torch.allclose(r, driver.terminal(x), atol=2e-14, rtol=2e-14)


def test_terminal_gradient_contribution_remains_in_controlled_path_graph() -> None:
    driver = make_driver("joint_concave", d=8, k=2)
    actor = TerminalAwareActor(driver, width=12, depth=1).double()
    x = torch.randn(7, driver.d, dtype=torch.float64, requires_grad=True)
    t = torch.full((len(x),), 0.2, dtype=torch.float64)
    _, q = actor.predict(t, x, driver.T)
    assert q is not None and q.requires_grad
    derivative = torch.autograd.grad(q.square().sum(), x)[0]
    assert torch.isfinite(derivative).all()
    assert float(derivative.norm()) > 1e-5
    assert ".detach" not in inspect.getsource(TerminalAwareActor.predict)


def test_terminal_aware_actor_supports_full_vjp_label_rollout() -> None:
    driver = make_driver("joint_concave", d=8, k=2)
    actor = TerminalAwareActor(driver, width=12, depth=1)
    t = np.asarray([0.1, 0.2])
    x = np.random.default_rng(3).standard_normal((2, driver.d))
    value, q, metadata = tangent_labels(
        driver,
        actor,
        t,
        x,
        paths=4,
        time_steps=2,
        seed=9,
        device=torch.device("cpu"),
        parent_batch=2,
        path_chunk=4,
    )
    assert np.isfinite(value).all()
    assert q is not None and np.isfinite(q).all()
    assert metadata["state_dependence_preserved"]


def test_actors_receive_full_ambient_state() -> None:
    driver = make_driver("joint_concave", d=100, k=5)
    for actor_class in (StandardActor, TerminalAwareActor):
        actor = actor_class(driver, width=16, depth=2)
        first = next(layer for layer in actor.network if isinstance(layer, torch.nn.Linear))
        assert first.in_features == 101
        assert not hasattr(actor, "B")


def test_nested_training_banks_and_disjoint_splits() -> None:
    small_t, small_x = nested_points(
        "joint_concave_d400_k5", 12, "iteration_1_train", 32,
        T=0.5, include_time=True, training_seed=0
    )
    large_t, large_x = nested_points(
        "joint_concave_d400_k5", 12, "iteration_1_train", 64,
        T=0.5, include_time=True, training_seed=0
    )
    valid_t, valid_x = nested_points(
        "joint_concave_d400_k5", 12, "iteration_1_validation", 32,
        T=0.5, include_time=True, training_seed=0
    )
    assert np.array_equal(small_t, large_t[:32])
    assert np.array_equal(small_x, large_x[:32])
    assert not np.array_equal(small_x, valid_x)
    assert not np.array_equal(small_t, valid_t)

