from __future__ import annotations

import numpy as np
import pytest
import torch

from invariant_region_mlp.experiments.baseline_n3.solver import (
    CorrectionBaselineMLP,
    DrawTape,
    METHODS,
)
from invariant_region_mlp.experiments.double_estimator.equations import (
    BASE_SEED,
    make_equation,
)

from .protocol import g3_tasks, g5_tasks
from .solver import TapeTorchDraws, TorchEquation, TorchMLP, torch_seed


def _relative_linf(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.max(np.abs(a - b))) / max(
        float(np.max(np.abs(b))), np.finfo(np.float64).tiny
    )


@pytest.mark.parametrize("pde_id", ["P1", "MR"])
def test_torch_equations_match_numpy_on_cpu(pde_id: str) -> None:
    equation = make_equation(pde_id, 20)
    adapter = TorchEquation(equation, torch.device("cpu"))
    rng = np.random.default_rng(8)
    t = rng.uniform(0.0, equation.T, 10)
    x = rng.uniform(-1.0, 1.0, (10, equation.d))
    state = rng.normal(size=(10, equation.d + 1))
    tt = torch.tensor(t, dtype=torch.float64)
    xx = torch.tensor(x, dtype=torch.float64)
    ss = torch.tensor(state, dtype=torch.float64)
    assert _relative_linf(adapter.terminal(xx).numpy(), equation.terminal(x)) < 1e-12
    assert _relative_linf(adapter.exact_state(tt, xx).numpy(), equation.exact_state(t, x)) < 1e-12
    assert _relative_linf(
        adapter.generator(ss).numpy(), equation.generator(state[:, 0], state[:, 1:])
    ) < 1e-12


@pytest.mark.parametrize(
    "method", ["raw", "path", "double", "double_path", "path_base", "double_path_base"]
)
def test_draw_tape_replay_matches_numpy_on_cpu(method: str) -> None:
    equation = make_equation("MR", 20)
    n, M = 3, 4
    t = np.array([0.01, 0.08], dtype=np.float64)
    x = np.linspace(-1.0, 1.0, 40, dtype=np.float64).reshape(2, 20)
    parts = [BASE_SEED, 20, n, M, 3, 0]
    tape = DrawTape("record")
    numpy_solver = CorrectionBaselineMLP(
        equation,
        M,
        METHODS[method],
        np.random.default_rng(np.random.SeedSequence(parts)),
        aux_seed_sequence=np.random.SeedSequence(parts + [0xD0B1E]),
        dose_rng=np.random.default_rng(np.random.SeedSequence(parts + [0xD05E])),
        draw_tape=tape,
    )
    expected = numpy_solver.solve(n, t, x)
    replay = tape.replay()
    torch_solver = TorchMLP(
        equation,
        M,
        method,
        torch.device("cpu"),
        seed_sequence=np.random.SeedSequence(parts),
        aux_seed_sequence=np.random.SeedSequence(parts + [0xD0B1E]),
        draw_source=TapeTorchDraws(torch.device("cpu"), replay),
    )
    observed = torch_solver.solve(
        n, torch.tensor(t, dtype=torch.float64), torch.tensor(x, dtype=torch.float64)
    ).numpy()
    assert _relative_linf(observed[:, 0], expected[:, 0]) < 1e-10
    assert _relative_linf(observed[:, 1:], expected[:, 1:]) < 1e-10
    assert replay.exhausted


def test_torch_seed_is_stable_and_uses_two_words() -> None:
    seed = np.random.SeedSequence([20261201, 400, 3, 10, 100, 0])
    first = torch_seed(seed)
    second = torch_seed(np.random.SeedSequence([20261201, 400, 3, 10, 100, 0]))
    assert first == second
    assert 0 <= first < 2**63 - 1


def test_registered_g3_and_g5_counts() -> None:
    assert len(g3_tasks("numpy_cpu")) == 960
    assert len(g3_tasks("torch_gpu")) == 960
    assert len(g5_tasks()) == 360
