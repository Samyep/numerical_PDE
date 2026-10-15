from __future__ import annotations

import math

import numpy as np
import torch

from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis.models import (
    make_actor,
)

from .common import driver_for, metrics
from .deep_picard import ZeroPicardActor, picard_labels
from .full_history_mlp import predict_full_history_mlp
from .models import AmortizedDeepBSDE, HardTerminalValue
from .neural_baselines import sdgd_objective
from .scasml import pinn_residual_hutchinson


def test_terminal_aware_actor_is_exact_terminal_tangent() -> None:
    driver = driver_for("nonseparable_d100_k2")
    actor = make_actor("terminal_aware", driver, width=16, depth=2)
    x = torch.randn(5, driver.d)
    t = torch.full((5,), driver.T)
    r, q = actor.predict(t, x, driver.T)
    assert torch.allclose(r, driver.terminal(x), atol=2e-7, rtol=2e-7)
    assert q is not None
    assert torch.allclose(
        q, driver.sigma * driver.terminal_gradient_torch(x), atol=2e-7, rtol=2e-7
    )


def test_hard_terminal_pinn_enforces_terminal_condition() -> None:
    driver = driver_for("joint_concave_d400_k2")
    model = HardTerminalValue(driver.d, width=8, depth=1)
    x = torch.randn(3, driver.d)
    t = torch.full((3,), driver.T)
    assert torch.equal(model.value(driver, t, x), driver.terminal(x))


def test_deep_bsde_uses_minus_driver_sign() -> None:
    driver = driver_for("nonseparable_d100_k2")
    model = AmortizedDeepBSDE(driver.d, width=8, depth=1)
    x = torch.randn(2, driver.d)
    normals = torch.zeros(2, 1, driver.d)
    y0 = model.initial_value(driver, x)
    t0 = torch.zeros(2)
    z0 = model.z_value(driver, t0, x)
    expected = y0 - driver.T * driver.driver(y0, z0)
    observed, final_x, nonfinite = model.rollout(driver, x, normals)
    assert torch.allclose(final_x, x)
    assert torch.allclose(observed, expected, atol=1e-6, rtol=1e-6)
    assert nonfinite == 0


def test_sdgd_objective_handles_generic_nonseparable_driver() -> None:
    driver = driver_for("nonseparable_d100_k2")
    model = HardTerminalValue(driver.d, width=8, depth=1)
    generator = torch.Generator().manual_seed(4)
    loss, diagnostics = sdgd_objective(
        model,
        driver,
        torch.rand(3) * driver.T,
        torch.randn(3, driver.d),
        coordinate_budget=2,
        generator=generator,
    )
    loss.backward()
    assert torch.isfinite(loss)
    assert math.isfinite(diagnostics["residual_product"])
    assert any(parameter.grad is not None for parameter in model.parameters())


class _QuadraticValue(torch.nn.Module):
    def value(self, driver, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        del driver
        return t + x.square().sum(dim=1)


def test_hutchinson_trace_is_exact_for_isotropic_quadratic() -> None:
    driver = driver_for("joint_concave_d400_k2")
    generator = torch.Generator().manual_seed(8)
    t = torch.rand(4) * driver.T
    x = torch.randn(4, driver.d)
    residual, _ = pinn_residual_hutchinson(
        _QuadraticValue(),
        "joint_concave_d400_k2",
        t,
        x,
        probes=2,
        generator=generator,
    )
    value = t + x.square().sum(dim=1)
    z = driver.sigma * 2.0 * x
    expected = 1.0 + 2.0 * driver.d + driver.driver(value, z)
    assert torch.allclose(residual, expected, atol=2e-4, rtol=2e-6)


def test_picard_zero_source_produces_finite_generic_labels() -> None:
    driver = driver_for("nonseparable_d100_k2")
    source = ZeroPicardActor(driver.d, driver.has_gradient)
    t = np.asarray([0.0, 0.2])
    x = np.zeros((2, driver.d))
    value, z, metadata = picard_labels(
        "nonseparable_d100_k2",
        source,
        t,
        x,
        paths=4,
        seed=9,
        device=torch.device("cpu"),
        parent_batch=2,
        path_chunk=2,
    )
    assert np.all(np.isfinite(value))
    assert z is not None and np.all(np.isfinite(z))
    assert metadata["nonfinite_count"] == 0
    # f(0,0)=0, so round-one values contain only the heat-semigroup terminal term.
    assert metadata["driver_std"] == 0.0


def test_full_history_level_one_is_finite() -> None:
    driver = driver_for("nonseparable_d100_k2")
    prediction, metadata = predict_full_history_mlp(
        "nonseparable_d100_k2",
        np.zeros((3, driver.d)),
        torch.device("cpu"),
        n=1,
        M=4,
        seed=11,
    )
    assert np.all(np.isfinite(prediction))
    assert metadata["stats"]["terminal_samples"] == 12
    assert metadata["training_seconds"] == 0.0


def test_metrics_match_independent_formula() -> None:
    prediction = np.asarray([1.0, 3.0, 5.0])
    reference = np.asarray([1.0, 2.0, 3.0])
    observed = metrics(prediction, reference)
    expected_rmse = math.sqrt(5.0 / 3.0)
    assert observed["rmse"] == expected_rmse
    assert observed["relative_l2"] == expected_rmse / math.sqrt(14.0 / 3.0)
    assert observed["bias"] == 1.0

