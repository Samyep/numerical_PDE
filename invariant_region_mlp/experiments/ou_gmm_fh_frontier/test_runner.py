"""CPU-only mathematics/small recursion sanity tests. Never benchmark CPU speed."""

import math

import numpy as np
import torch

from invariant_region_mlp.experiments.tangent_actor_icml_final_suite.benchmarks import OUGMMBenchmark
from .runner import (
    ALPHA, RawOUFullHistory, ou_state, ou_variance,
    weighted_ou_bel,
)


def _benchmark(d=5):
    generator = np.random.default_rng(17)
    mean = generator.uniform(-1, 1, size=(5, d))
    pi = generator.dirichlet(np.ones(5))
    return OUGMMBenchmark(d, mean, pi)


def test_OU_variance_and_transition_moments():
    dt = torch.tensor(0.65, dtype=torch.float64)
    x = torch.tensor([0.9, -0.3], dtype=torch.float64)
    gen = torch.Generator().manual_seed(20261020)
    normal = torch.randn((80000, 2), generator=gen, dtype=torch.float64)
    state = ou_state(x, dt.expand(len(normal)), normal)
    target_mean = math.exp(dt.item()) * x
    target_variance = ou_variance(dt).item()
    assert torch.allclose(state.mean(0), target_mean, atol=0.015, rtol=0)
    assert torch.allclose(state.var(0), torch.full((2,),target_variance,dtype=torch.float64),
                          atol=0.025,rtol=0)


def test_OU_BEL_gradient_for_quadratic_payoff():
    dt = torch.tensor(0.4, dtype=torch.float64)
    x = torch.tensor([0.5,-0.2], dtype=torch.float64)
    rng = torch.Generator().manual_seed(1234)
    normal = torch.randn((200000, 2), generator=rng, dtype=torch.float64)
    state = ou_state(x, dt.expand(len(normal)), normal)
    payoff = state.square().sum(-1)
    exact = 2.0 * math.exp(2.0 * dt.item()) * x
    bel = torch.exp(dt) / torch.sqrt(ou_variance(dt)) * normal
    empirical = (payoff[:,None]*bel).mean(0)
    assert torch.allclose(empirical, exact, atol=0.09, rtol=0)


def test_weighted_BEL_matches_direct_expression():
    t_remaining = torch.tensor([0.3,0.7],dtype=torch.float64)
    fractions = torch.tensor([[1e-8,0.02,0.3],[0.05,0.1,0.9]],dtype=torch.float64)
    delta = t_remaining[:,None]*fractions
    direct = (
        t_remaining[:,None]*fractions.pow(1-ALPHA)/ALPHA
        * delta.exp()/torch.sqrt(ou_variance(delta))
    )
    stable = weighted_ou_bel(t_remaining,fractions,delta)
    assert torch.allclose(stable,direct,rtol=1e-10,atol=1e-10)


def test_GMM_terminal_and_gradient_match_definition():
    bench = _benchmark(7)
    model = RawOUFullHistory(bench, 2, torch.device("cpu"), seed=15, dtype=torch.float64)
    x = torch.randn((9,7),generator=torch.Generator().manual_seed(22),dtype=torch.float64)
    our_value, our_gradient = model.terminal_value_gradient(x)
    expected_value = bench.terminal(x)
    expected_gradient = bench.terminal_gradient(x)
    assert torch.allclose(our_value,expected_value,rtol=1e-12,atol=1e-12)
    assert torch.allclose(our_gradient,expected_gradient,rtol=1e-12,atol=1e-12)


def test_terminal_time_estimates_exact_terminal_and_gradient():
    bench = _benchmark(5)
    model = RawOUFullHistory(bench, 3, torch.device("cpu"), seed=4, dtype=torch.float64)
    x = torch.randn((3,5),generator=torch.Generator().manual_seed(6),dtype=torch.float64)
    t = torch.ones(3,dtype=torch.float64)*bench.T
    with torch.no_grad():
        output = model.solve(2, t, x)
    value, grad = model.terminal_value_gradient(x)
    assert torch.allclose(output[:,0],value,rtol=1e-12,atol=1e-12)
    assert torch.allclose(output[:,1:],grad,rtol=1e-12,atol=1e-12)


def test_exact_OU_GMM_PDE_gaussian_single_component_residual():
    # Every GMM with all means=0 reduces to this isotropic Gaussian case.
    # Check the sign of the OU drift x·∇ and deterministic source -d.
    d = 100
    for tau in (0.05,0.4,0.9):
        v = 2*math.exp(-2*tau) + .5*(1-math.exp(-2*tau))
        x_squared = 360.0
        v_t = 2*v-1
        u_t = (d/(2*v)-x_squared/(2*v*v))*v_t
        base_drift = x_squared/v
        laplace = d/v
        grad_squared = x_squared/(v*v)
        residual = u_t + base_drift + .5*laplace - .5*grad_squared - d
        assert abs(residual)<1e-10


def test_small_raw_OU_FH_recursion_finite_and_no_oracle():
    bench = _benchmark(3)
    model = RawOUFullHistory(bench, M=2, device=torch.device("cpu"), seed=19)
    tt = torch.tensor([0.0,0.3],dtype=torch.float32)
    xx = torch.randn(2,3,generator=torch.Generator().manual_seed(20))
    with torch.no_grad():
        output = model.solve(3, tt, xx)
    assert output.shape == (2,4)
    assert torch.isfinite(output).all()
    assert model.stats.terminal_samples > 0
    assert model.stats.driver_evaluations > 0


def test_full_five_component_exact_solution_has_correct_OU_drift_sign():
    bench = _benchmark(4)
    t = torch.tensor([0.2,0.65], dtype=torch.float64, requires_grad=True)
    x = torch.randn((2,4),generator=torch.Generator().manual_seed(99),
                    dtype=torch.float64,requires_grad=True)
    u = bench.exact_value(t,x)
    u_t = torch.autograd.grad(u.sum(),t,create_graph=True)[0]
    grad = torch.autograd.grad(u.sum(),x,create_graph=True)[0]
    laplace = torch.zeros(2,dtype=torch.float64)
    for j in range(4):
        component = torch.autograd.grad(
            grad[:,j].sum(),x,retain_graph=True,create_graph=True
        )[0][:,j]
        laplace += component
    residual = (
        u_t + .5*laplace + (x*grad).sum(dim=-1)
        - .5*grad.square().sum(dim=-1) - 4.0
    )
    assert torch.max(torch.abs(residual)).item()<1e-10
