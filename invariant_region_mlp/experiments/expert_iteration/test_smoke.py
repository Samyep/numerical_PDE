"""Fast numerical invariants for the expert-iteration implementation."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import torch

from lqg import LQGEquation, uniform_unit_ball
from network import FrozenMLP, PINN
from solver import run_repetition


def _checkpoint(path: Path, d: int) -> Path:
    torch.manual_seed(17)
    model = PINN(d).float()
    checkpoint = path / "checkpoint.pt"
    torch.save(
        {"state_dict": model.state_dict(), "d": d, "network_seed": 0, "step": 1},
        checkpoint,
    )
    return checkpoint


def test_scaled_quadrature_is_exact_at_terminal() -> None:
    equation = LQGEquation.create(8)
    x = uniform_unit_ball(np.random.default_rng(7), 12, equation.d)
    t = np.full(len(x), equation.T)
    u, z = equation.hopf_cole(t, x)
    np.testing.assert_allclose(u, equation.terminal(x), rtol=0.0, atol=1e-14)
    np.testing.assert_allclose(z, equation.sigma_grad_terminal(x), rtol=0.0, atol=1e-14)


def test_manual_float64_derivatives_match_autograd(tmp_path: Path) -> None:
    d = 6
    checkpoint = _checkpoint(tmp_path, d)
    frozen = FrozenMLP.load(
        checkpoint, d=d, network_seed=0, device="cpu", exact_laplacian=True
    )
    rng = np.random.default_rng(9)
    x = rng.normal(size=(4, d))
    t = rng.uniform(0.0, 0.5, 4)
    u, z, u_t, lap = frozen.evaluate(t, x, need_laplacian=True)

    payload = torch.load(checkpoint, map_location="cpu", weights_only=True)
    model = PINN(d).double()
    model.load_state_dict(payload["state_dict"])
    inputs = torch.tensor(np.concatenate([x, t[:, None]], axis=1), requires_grad=True)
    value = model(inputs)[:, 0]
    gradient = torch.autograd.grad(value.sum(), inputs, create_graph=True)[0]
    exact_lap = torch.zeros(len(x), dtype=torch.float64)
    for coordinate in range(d):
        exact_lap += torch.autograd.grad(
            gradient[:, coordinate].sum(), inputs, retain_graph=True
        )[0][:, coordinate]
    np.testing.assert_allclose(u, value.detach().numpy(), rtol=0.0, atol=2e-14)
    np.testing.assert_allclose(z, np.sqrt(2.0) * gradient[:, :d].detach().numpy(), rtol=0.0, atol=2e-14)
    np.testing.assert_allclose(u_t, gradient[:, d].detach().numpy(), rtol=0.0, atol=2e-14)
    assert lap is not None
    np.testing.assert_allclose(lap, exact_lap.detach().numpy(), rtol=0.0, atol=2e-13)


def test_methods_share_identical_tree_draws(tmp_path: Path) -> None:
    d = 4
    equation = LQGEquation.create(d)
    checkpoint = _checkpoint(tmp_path, d)
    surrogate = FrozenMLP.load(
        checkpoint, d=d, network_seed=0, device="cpu", exact_laplacian=True
    )
    rng = np.random.default_rng(11)
    t = rng.uniform(0.0, equation.T, 4)
    x = uniform_unit_ball(rng, 4, d)
    is_validation = np.array([True, False, False, False])
    u, z = equation.hopf_cole(t, x)
    reference = np.concatenate([u[:, None], z], axis=1)
    fingerprints = []
    for method in ("path", "scasml_noclip"):
        row, _ = run_repetition(
            equation=equation,
            surrogate=surrogate,
            method_name=method,
            checkpoint=1,
            network_seed=0,
            n=2,
            M=2,
            repetition=0,
            t=t,
            x=x,
            is_validation=is_validation,
            reference_state=reference,
            chunk_size=2,
            trace_draws=True,
        )
        fingerprints.append(row["draw_fingerprint"])
    assert fingerprints[0] == fingerprints[1]

