from __future__ import annotations

import math

import numpy as np
import torch

from .benchmarks import AllenCahnBenchmark, DPIOUBenchmark, HJBLogBenchmark
from .deep_bsde_port import CONFIGS, OfficialDeepBSDE
from .diagnostics import exact_z_torch
from .evaluate_actor import _select_reference_z
from .full_history_mlp_original import OriginalPDEMLPDriver
from .models import TerminalAwareActor
from .references import dpi_exact_reference, hjb_log_quad_scalar, hjb_log_reference
from .tangent_solver import _rollout_chunk, antithetic_normals
from invariant_region_mlp.experiments.tangent_actor_strong_baselines.full_history_mlp import (
    GenericFullHistoryMLP,
)


def test_hjb_reference_matches_adaptive_integral() -> None:
    rng = np.random.default_rng(42)
    x = rng.normal(size=(4, 100))
    t = np.asarray([0.0, 0.3, 0.8, 0.99])
    value, z = hjb_log_reference(t, x, order=128, return_gradient=True)
    expected = np.asarray(
        [hjb_log_quad_scalar(float(ti), float(xi @ xi)) for ti, xi in zip(t, x)]
    )
    np.testing.assert_allclose(value, expected, rtol=0.0, atol=2e-12)
    assert z.shape == x.shape


def test_terminal_actor_is_exact_at_terminal_time() -> None:
    for benchmark in (HJBLogBenchmark(), AllenCahnBenchmark(), DPIOUBenchmark()):
        actor = TerminalAwareActor(
            benchmark, width=16, depth=2, activation="silu"
        )
        x = torch.randn(5, benchmark.d)
        t = torch.full((5,), benchmark.T)
        r, q = actor.predict(t, x)
        if benchmark.actor_kind == "gradient":
            assert q is not None
            torch.testing.assert_close(
                q, benchmark.sigma * benchmark.terminal_gradient(x)
            )
        else:
            assert r is not None
            torch.testing.assert_close(r, benchmark.terminal(x))


def test_tangent_identity_at_actor_point() -> None:
    generator = torch.Generator().manual_seed(9)
    hjb = HJBLogBenchmark()
    q = torch.randn(7, hjb.d, generator=generator)
    c, b, h = hjb.tangent_coefficients(None, q)
    tangent = c * 0.0 + (b * q).sum(dim=1) + h
    torch.testing.assert_close(tangent, -0.5 * q.square().sum(dim=1))

    allen = AllenCahnBenchmark()
    r = torch.rand(7, generator=generator)
    c, b, h = allen.tangent_coefficients(r, None)
    torch.testing.assert_close(c * r + h, r - r.pow(3))
    assert torch.count_nonzero(b) == 0


def test_dpi_numpy_and_torch_exact_solution_and_gradient() -> None:
    benchmark = DPIOUBenchmark()
    rng = np.random.default_rng(18)
    t = rng.uniform(0.0, 1.0, 6)
    x = rng.normal(size=(6, 100))
    value, z = dpi_exact_reference(benchmark, t, x, return_gradient=True)
    tt = torch.as_tensor(t, dtype=torch.float64)
    xx = torch.as_tensor(x, dtype=torch.float64).requires_grad_(True)
    torch_value = benchmark.exact_value(tt, xx)
    autograd = torch.autograd.grad(torch_value.sum(), xx)[0]
    np.testing.assert_allclose(value, torch_value.detach().numpy(), atol=2e-12)
    np.testing.assert_allclose(z, autograd.detach().numpy(), atol=2e-12)


def test_torch_exact_z_matches_independent_numpy_references() -> None:
    rng = np.random.default_rng(5103)
    for benchmark in (HJBLogBenchmark(), DPIOUBenchmark()):
        t = rng.uniform(0.0, benchmark.T, 7)
        x = rng.normal(size=(7, benchmark.d))
        if isinstance(benchmark, HJBLogBenchmark):
            _, expected = hjb_log_reference(t, x, return_gradient=True)
        else:
            _, expected = dpi_exact_reference(
                benchmark, t, x, return_gradient=True
            )
        actual = exact_z_torch(
            benchmark,
            torch.as_tensor(t, dtype=torch.float64),
            torch.as_tensor(x, dtype=torch.float64),
        ).numpy()
        np.testing.assert_allclose(actual, expected, rtol=2e-10, atol=2e-10)


def test_mcv_uses_unscaled_left_endpoint_dw() -> None:
    benchmark = HJBLogBenchmark()
    actor = TerminalAwareActor(
        benchmark, width=8, depth=1, activation="silu"
    )
    # At T the path has zero duration, so use a nonterminal state and compare
    # the implementation's P-H difference to the directly reconstructed H.
    t = torch.tensor([0.5])
    x = torch.zeros((1, 100))
    normals = torch.zeros((1, 2, 1, 100))
    normals[0, 0, 0, 0] = 0.75
    normals[0, 1, 0, 0] = -0.75
    plain, corrected, _, _, _ = _rollout_chunk(
        benchmark, actor, t, x, normals, mcv=True, need_gradient=False
    )
    _, q = actor.predict(t, x)
    assert q is not None
    dw = math.sqrt(0.5) * normals[:, :, 0]
    expected_h = (q[:, None, :] * dw).sum(dim=-1)
    torch.testing.assert_close(plain - corrected, expected_h)


def test_deep_bsde_native_scalar_architecture_and_rollout() -> None:
    benchmark = HJBLogBenchmark()
    config = CONFIGS[benchmark.name]
    model = OfficialDeepBSDE(benchmark, config).double()
    assert model.y_init.shape == (1,)
    assert model.z_init.shape == (1, 100)
    assert len(model.subnets) == 19
    dw = torch.zeros((64, 20, 100), dtype=torch.float64)
    model.train()
    y, x = model.rollout(dw)
    assert y.shape == (64,)
    assert x.shape == (64, 100)
    assert torch.isfinite(y).all()


def test_antithetic_normals_are_exact_pairs() -> None:
    normals = antithetic_normals(
        3,
        8,
        2,
        100,
        generator=torch.Generator().manual_seed(5),
        device=torch.device("cpu"),
    )
    torch.testing.assert_close(normals[:, :4], -normals[:, 4:])


def test_value_actor_evaluation_does_not_require_gradient_truth() -> None:
    indices = np.asarray([0, 2])
    value_bank = {"t": np.zeros(3), "z": np.empty((0, 100))}
    assert _select_reference_z(value_bank, indices, "value") is None

    gradient_bank = {
        "t": np.zeros(3),
        "z": np.arange(300, dtype=np.float64).reshape(3, 100),
    }
    np.testing.assert_array_equal(
        _select_reference_z(gradient_bank, indices, "gradient"),
        gradient_bank["z"][indices],
    )


def test_full_history_value_only_path_skips_unused_gradient() -> None:
    driver = OriginalPDEMLPDriver(AllenCahnBenchmark())

    def forbidden_gradient(_: torch.Tensor) -> torch.Tensor:
        raise AssertionError("a value-only recursion must not request terminal gradients")

    driver.terminal_gradient_torch = forbidden_gradient  # type: ignore[method-assign]
    solver = GenericFullHistoryMLP(
        driver, M=4, device=torch.device("cpu"), seed=81
    )
    state = solver.solve(2, torch.zeros(2), torch.zeros((2, driver.d)))
    assert torch.isfinite(state).all()
    assert state.shape == (2, 1)
    assert torch.count_nonzero(state[:, 1:]) == 0

