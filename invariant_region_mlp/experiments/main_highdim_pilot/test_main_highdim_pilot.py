from __future__ import annotations

import math

import numpy as np
import torch

from .actors import JointTerminalActor
from .b6_adapter import load_b6_reuse
from .common import load_dense, load_lqg
from .equations import CoupledLQGEquation, DenseFourierDriver
from .prepare import prepare
from .protocol import DIMENSION, JOINT_CONFIG, frozen_configuration


def _dense_small() -> DenseFourierDriver:
    rng = np.random.default_rng(18)
    d = 4
    omega = rng.standard_normal((2 * d, d))
    omega /= np.linalg.norm(omega, axis=1, keepdims=True)
    signs = (2 * rng.integers(0, 2, 2 * d) - 1).astype(float)
    phases = rng.uniform(0.0, 2.0 * np.pi, 2 * d)
    return DenseFourierDriver(omega, signs, phases, T=0.5, parameter_seed=(18,))


def test_dense_exact_terminal_gradient_and_tangent() -> None:
    driver = _dense_small()
    rng = np.random.default_rng(20)
    x = rng.normal(size=(3, driver.d))
    state = driver.exact_state(np.full(3, driver.T), x)
    np.testing.assert_allclose(state[:, 0], driver.terminal_numpy(x), atol=2e-14)
    np.testing.assert_allclose(
        state[:, 1:],
        math.sqrt(2.0) * driver.terminal_gradient_numpy(x),
        atol=2e-14,
    )
    r = torch.tensor(rng.normal(size=3), dtype=torch.float64)
    q = torch.tensor(rng.normal(size=(3, driver.d)), dtype=torch.float64)
    c, b, h = driver.tangent_coefficients(r, q)
    ac, ab, ah = driver.tangent_coefficients_autodiff(r, q)
    torch.testing.assert_close(c, ac, atol=2e-13, rtol=2e-13)
    torch.testing.assert_close(b, ab, atol=2e-13, rtol=2e-13)
    torch.testing.assert_close(h, ah, atol=2e-13, rtol=2e-13)
    torch.testing.assert_close(
        c * r + (b * q).sum(dim=-1) + h,
        driver.driver(r, q),
        atol=2e-13,
        rtol=2e-13,
    )


def test_joint_actor_is_exactly_terminal_anchored() -> None:
    driver = _dense_small()
    actor = JointTerminalActor(driver, width=16, depth=2)
    x = torch.randn(5, driver.d)
    terminal_time = torch.full((5,), driver.T)
    r, q = actor.predict(terminal_time, x, driver.T)
    torch.testing.assert_close(r, driver.terminal(x))
    torch.testing.assert_close(q, driver.sigma * driver.terminal_gradient_torch(x))
    assert actor.specification()["value_projection"] is False


def test_lqg_matrix_reference_and_gradient() -> None:
    rng = np.random.default_rng(25)
    c1 = rng.uniform(0.5, 1.5, 5)
    c2 = rng.uniform(0.5, 1.5, 5)
    equation = CoupledLQGEquation(c1, c2, T=0.5, coefficient_seed=(25,))
    x = rng.normal(size=(4, equation.d))
    direct = np.sum(
        c1 * (x[:, :-1] - x[:, 1:]) ** 2 + c2 * x[:, 1:] ** 2,
        axis=1,
    )
    np.testing.assert_allclose(equation.quadratic(x), direct, atol=2e-14)
    assert equation.eigenvalues[0] > 0.0
    t = np.array([0.1, 0.2, 0.3, 0.4])
    coarse = equation.exact_state(t, x, order=96)
    fine = equation.exact_state(t, x, order=192)
    np.testing.assert_allclose(coarse, fine, atol=2e-9, rtol=2e-9)
    epsilon = 2e-5
    plus, minus = x[0].copy(), x[0].copy()
    plus[2] += epsilon
    minus[2] -= epsilon
    fd = (
        equation.exact_u(t[0], plus[None, :], order=192)[0]
        - equation.exact_u(t[0], minus[None, :], order=192)[0]
    ) / (2 * epsilon)
    analytic = equation.exact_z(t[0], x[:1], order=192)[0, 2] / math.sqrt(2.0)
    assert abs(fd - analytic) < 2e-7


def test_frozen_inputs_are_idempotent_and_dimension_100() -> None:
    first = prepare()
    second = prepare()
    assert first == second
    assert frozen_configuration()["dimension"] == DIMENSION == 100
    dense, lqg = load_dense(), load_lqg()
    assert dense.d == lqg.d == 100
    assert first["coupled_lqg"]["benchmark_label"] == "RECREATED_PUBLISHED_BENCHMARK"
    assert JOINT_CONFIG.optimizer_steps == 1500


def test_b6_read_only_adapter_audits_both_frozen_dimensions() -> None:
    expected_banks = {
        100: "84b2132607bc6b17e44adbea838915f8ce657813f6a0987e34deca308660f1c1",
        400: "0df59abc155e90aa674532760578dc2563cbdb4aade97a849961c8a8edd9f2f9",
    }
    for dimension, bank_hash in expected_banks.items():
        reuse = load_b6_reuse(dimension)
        assert reuse["audit_status"] == "PASS"
        assert reuse["label"] == "REUSED_MATURE_RESULT"
        assert reuse["calibration"]["bank_sha256"] == bank_hash
        assert [row["actor_request"] for row in reuse["iterations"]] == [
            "a0",
            "a1",
            "a2",
            "a3",
            "a4",
        ]

