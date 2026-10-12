from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest
import torch
from scipy.fft import dct, idct

from invariant_region_mlp.experiments.baseline_n3.solver import DrawTape, METHODS
from invariant_region_mlp.experiments.double_estimator.equations import BASE_SEED
from invariant_region_mlp.experiments.gpu_mlp.solver import TapeTorchDraws

from .equations import SpectralSeparableHJB, make_equation
from .analyze import _evaluate_block
from .protocol import (
    all_tasks,
    stage_a_tasks,
    stage_b_extension_tasks,
    stage_b_n4_tasks,
    stage_b_s1_tasks,
)
from .solver import (
    NonRidgeCorrectionBaselineMLP,
    NonRidgeTorchMLP,
    SpectralTorchEquation,
    torch_dct_ortho,
    torch_idct_ortho,
)


def _relative_linf(observed: np.ndarray, expected: np.ndarray) -> float:
    return float(np.max(np.abs(observed - expected))) / max(
        float(np.max(np.abs(expected))), np.finfo(np.float64).tiny
    )


@pytest.mark.parametrize("d", [7, 20, 100])
def test_torch_dct_matches_scipy(d: int) -> None:
    rng = np.random.default_rng(d)
    values = rng.normal(size=(5, d))
    tensor = torch.tensor(values, dtype=torch.float64)
    observed = torch_dct_ortho(tensor).numpy()
    expected = dct(values, type=2, norm="ortho", axis=-1)
    assert _relative_linf(observed, expected) < 2e-15
    reconstructed = torch_idct_ortho(torch.tensor(observed)).numpy()
    scipy_reconstructed = idct(expected, type=2, norm="ortho", axis=-1)
    assert _relative_linf(reconstructed, scipy_reconstructed) < 3e-15


def test_special_case_is_closed_form_p1() -> None:
    equation = SpectralSeparableHJB(d=20, alpha=0.0, pde_id="S0")
    rng = np.random.default_rng(19)
    t = rng.uniform(0.0, equation.T, 12)
    x = rng.uniform(-1.0, 1.0, (12, equation.d))
    y = equation.rotate(x)[:, 0]
    c = np.zeros(equation.d)
    c[0] = 1.0
    value, derivative = equation.one_dimensional_state(
        c[None, :], (equation.T - t)[:, None], equation.rotate(x), order=96
    )
    logits = (
        y[:, None] * np.array([0.5, 1.5, 3.0])
        + (equation.T - t)[:, None] * np.array([0.5, 1.5, 3.0]) ** 2
    )
    maximum = np.max(logits, axis=1, keepdims=True)
    logsum = maximum[:, 0] + np.log(np.sum(np.exp(logits - maximum), axis=1))
    expected_u = -(logsum - math.log(3.0))
    weights = np.exp(logits - logsum[:, None])
    expected_dy = -np.sum(weights * np.array([0.5, 1.5, 3.0]), axis=1)
    assert _relative_linf(np.sum(value, axis=1), expected_u) < 1e-13
    assert _relative_linf(derivative[:, 0], expected_dy) < 1e-13
    assert np.max(np.abs(derivative[:, 1:])) == 0.0


@pytest.mark.parametrize("pde_id", ["S0", "S05", "S1"])
def test_torch_equation_matches_numpy(pde_id: str) -> None:
    equation = make_equation(pde_id, 20)
    adapter = SpectralTorchEquation(equation, torch.device("cpu"))
    rng = np.random.default_rng(4)
    x = rng.uniform(-1.0, 1.0, (8, 20))
    first = rng.normal(size=(8, 21))
    second = rng.normal(size=(8, 21))
    xx = torch.tensor(x, dtype=torch.float64)
    one = torch.tensor(first, dtype=torch.float64)
    two = torch.tensor(second, dtype=torch.float64)
    assert _relative_linf(adapter.terminal(xx).numpy(), equation.terminal(x)) < 1e-13
    assert _relative_linf(
        adapter.terminal_gradient_z(xx).numpy(), equation.terminal_gradient_z(x)
    ) < 1e-13
    assert _relative_linf(
        adapter.generator(one).numpy(), equation.generator(first[:, 0], first[:, 1:])
    ) < 1e-14
    assert _relative_linf(
        adapter.double_generator(one, two).numpy(),
        equation.double_generator(first[:, 1:], second[:, 1:]),
    ) < 1e-14


@pytest.mark.parametrize("method", ["path", "double_path_base"])
def test_draw_tape_replay_matches_cpu(method: str) -> None:
    equation = make_equation("S05", 20)
    n, M = 2, 3
    t = np.array([0.02, 0.17], dtype=np.float64)
    x = np.linspace(-1.0, 1.0, 40, dtype=np.float64).reshape(2, 20)
    seed_parts = [BASE_SEED, 20, n, M, 0, 0]
    tape = DrawTape("record")
    cpu = NonRidgeCorrectionBaselineMLP(
        equation,
        M,
        METHODS[method],
        np.random.default_rng(np.random.SeedSequence(seed_parts)),
        aux_seed_sequence=np.random.SeedSequence(seed_parts + [0xD0B1E]),
        dose_rng=np.random.default_rng(np.random.SeedSequence(seed_parts + [0xD05E])),
        draw_tape=tape,
    )
    expected = cpu.solve(n, t, x)
    replay = tape.replay()
    gpu = NonRidgeTorchMLP(
        equation,
        M,
        method,
        torch.device("cpu"),
        seed_sequence=np.random.SeedSequence(seed_parts),
        aux_seed_sequence=np.random.SeedSequence(seed_parts + [0xD0B1E]),
        draw_source=TapeTorchDraws(torch.device("cpu"), replay),
    )
    observed = gpu.solve(
        n, torch.tensor(t, dtype=torch.float64), torch.tensor(x, dtype=torch.float64)
    ).numpy()
    assert _relative_linf(observed[:, 0], expected[:, 0]) < 1e-10
    assert _relative_linf(observed[:, 1:], expected[:, 1:]) < 1e-10
    assert replay.exhausted


def test_registered_task_counts_and_uniqueness() -> None:
    assert len(stage_a_tasks()) == 440
    assert len(stage_b_s1_tasks()) == 220
    assert len(stage_b_extension_tasks()) == 180
    assert len(stage_b_n4_tasks()) == 60
    tasks = all_tasks()
    identities = {
        (t["pde_id"], t["d"], t["n"], t["M"], t["method"], t["repetition"])
        for t in tasks
    }
    assert len(tasks) == 900
    assert len(identities) == len(tasks)


def test_mechanical_stage_a_criteria() -> None:
    rows: list[dict[str, object]] = []
    for family in ("S0", "S05"):
        for d in (100, 400):
            entries = (
                (2, 16, "path", 0.5),
                (2, 64, "path", 0.4),
                (2, 64, "raw", 2.0),
                (2, 256, "path", 0.3),
                (3, 16, "raw", 2.0),
                (3, 16, "path", 0.5),
                (3, 16, "path_base", 0.1),
                (3, 16, "double_path", 0.2),
                (3, 16, "double_path_base", 0.1),
                (3, 24, "double_path", 0.15),
                (3, 24, "double_path_base", 0.08),
            )
            for n, M, method, skill in entries:
                rows.append(
                    {
                        "stage": "stage_a",
                        "pde_id": family,
                        "dimension": d,
                        "n": n,
                        "M": M,
                        "method": method,
                        "repetitions": 10,
                        "mean_skill": skill,
                    }
                )
    result = _evaluate_block(pd.DataFrame(rows), ("S0", "S05"), n3_required=7)
    for criterion in ("N-1", "N-2", "N-3", "N-4", "N-5"):
        assert result[criterion]["verdict"] == "PASS"
