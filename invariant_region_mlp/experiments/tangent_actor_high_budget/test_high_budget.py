from __future__ import annotations

import numpy as np
import torch

from .actor_scaling import nested_bank_audit
from .common import driver_for
from .models import OriginalAmortizedDeepBSDE, RawTanhPINN
from .original_baselines import sdgd_objective
from .evaluate import _decomposition


def test_original_baselines_have_no_terminal_aware_ansatz() -> None:
    driver = driver_for("nonseparable_d100_k2")
    x = torch.randn(4, driver.d)
    t = torch.full((4,), driver.T)
    pinn = RawTanhPINN(driver.d, width=8, depth=1)
    deep_bsde = OriginalAmortizedDeepBSDE(driver.d, time_steps=2, width=8)
    assert not torch.equal(pinn.value(driver, t, x), driver.terminal(x))
    assert not torch.equal(deep_bsde.initial_value(driver, x), driver.terminal(x))
    assert not any("terminal" in name for name, _ in deep_bsde.named_modules())


def test_original_deep_bsde_rollout_has_correct_sign() -> None:
    driver = driver_for("nonseparable_d100_k2")
    model = OriginalAmortizedDeepBSDE(driver.d, time_steps=1, width=8)
    model.eval()
    x = torch.randn(3, driver.d)
    normal = torch.zeros(3, 1, driver.d)
    y0 = model.initial_value(driver, x)
    z0 = model.z_nets[0](x) / driver.d
    observed, final_x, nonfinite = model.rollout(driver, x, normal)
    expected = y0 - driver.T * driver.driver(y0, z0)
    assert torch.allclose(observed, expected, atol=1e-6, rtol=1e-6)
    assert torch.equal(final_x, x)
    assert nonfinite == 0


def test_sdgd_raw_pinn_objective_is_finite() -> None:
    driver = driver_for("nonseparable_d100_k2")
    model = RawTanhPINN(driver.d, width=8, depth=1)
    generator = torch.Generator().manual_seed(3)
    objective, diagnostics = sdgd_objective(
        model, driver, torch.rand(3) * driver.T, torch.randn(3, driver.d),
        torch.randn(3, driver.d), coordinate_budget=2, generator=generator,
    )
    objective.backward()
    assert torch.isfinite(objective)
    assert np.isfinite(diagnostics["terminal_penalty"])


def test_nested_point_banks_are_exact_prefixes() -> None:
    audit = nested_bank_audit("nonseparable_d100_k2", 0)
    assert audit["prefix_equality_verified"]
    assert len(set(audit["prefix_fingerprints"].values())) == 3


def test_systematic_decomposition_subtracts_variance_of_mean() -> None:
    prediction = np.asarray([1.2, 1.8])
    reference = np.asarray([1.0, 2.0])
    variance = np.asarray([0.01, 0.01])
    result = _decomposition(
        prediction, reference, variance, seed=7, bootstrap_repetitions=100
    )
    assert np.isclose(result["observed_mse"], 0.04)
    assert np.isclose(result["mc_variance_contribution"], 0.01)
    assert np.isclose(result["systematic_mse_proxy"], 0.03)
