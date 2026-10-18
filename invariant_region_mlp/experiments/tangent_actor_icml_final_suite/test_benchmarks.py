from __future__ import annotations

import numpy as np
import torch

from .benchmarks import AllenCahnBenchmark, HJBLogBenchmark, MLPFailureBenchmark


def test_terminal_gradients_match_autodiff() -> None:
    x = torch.randn((16, 10), dtype=torch.float64, requires_grad=True)
    for benchmark in (HJBLogBenchmark(10), AllenCahnBenchmark(10)):
        truth = torch.autograd.grad(benchmark.terminal(x).sum(), x, retain_graph=True)[0]
        assert torch.allclose(benchmark.terminal_gradient(x), truth, atol=1e-12, rtol=1e-12)


def test_mlp_failure_exact_native_value() -> None:
    benchmark = MLPFailureBenchmark(20)
    value = benchmark.exact_value_numpy(np.array([0.0]), np.zeros((1, 20)))[0]
    assert abs(value - np.sqrt(2.0 / np.pi)) < 1e-14


def test_mlp_failure_minimum_norm_subgradient() -> None:
    benchmark = MLPFailureBenchmark(20)
    q = torch.zeros((3, 20), dtype=torch.float64)
    _, b, h = benchmark.tangent_coefficients(None, q)
    assert torch.equal(b, torch.zeros_like(b))
    assert torch.equal(h, torch.zeros_like(h))
