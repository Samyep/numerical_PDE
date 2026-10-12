"""Adapters around the established float64 full-history GPU MLP recursion."""

from __future__ import annotations

from dataclasses import asdict
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.baseline_n3.solver import (
    CorrectionBaselineMLP,
    DrawTape,
    METHODS,
)
from invariant_region_mlp.experiments.gpu_mlp.solver import TapeTorchDraws, TorchMLP

from .common import stable_seed
from .equation import B1Equation


MLP_METHODS = {
    "vanilla_mlp": "raw",
    "double_mlp": "double",
    "actor_label": "path_base",
}


def evaluate_mlp(
    equation: B1Equation,
    t: np.ndarray,
    x: np.ndarray,
    *,
    n: int,
    M: int,
    method: str,
    repetition: int,
    device: torch.device,
    chunk_size: int | None = None,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Evaluate one established MLP method with independent chunk seeds."""

    internal_method = MLP_METHODS.get(method, method)
    if internal_method not in METHODS:
        raise ValueError(f"unknown MLP method {method!r}")
    times = np.asarray(t, dtype=np.float64).reshape(-1)
    states = np.asarray(x, dtype=np.float64).reshape(len(times), equation.d)
    chunk = len(times) if chunk_size is None else min(chunk_size, len(times))
    output = np.empty((len(times), equation.d + 1), dtype=np.float64)
    totals: dict[str, int] = {
        "terminal_samples": 0,
        "transition_samples": 0,
        "generator_calls": 0,
        "recursively_evaluated_states": 0,
        "recursive_calls": 0,
    }
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    for chunk_index, start in enumerate(range(0, len(times), chunk)):
        stop = min(start + chunk, len(times))
        seed = stable_seed(
            "B1", "mlp", method, equation.d, equation.k, n, M, repetition, chunk_index
        )
        bundle = np.random.SeedSequence(seed)
        primary, auxiliary = bundle.spawn(2)
        solver = TorchMLP(
            equation,
            M,
            internal_method,
            device,
            seed_sequence=primary,
            aux_seed_sequence=auxiliary,
        )
        with torch.no_grad():
            result = solver.solve(
                n,
                torch.tensor(times[start:stop], dtype=torch.float64, device=device),
                torch.tensor(states[start:stop], dtype=torch.float64, device=device),
            )
        output[start:stop] = result.cpu().numpy()
        for key in totals:
            totals[key] += int(getattr(solver.stats, key))
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        peak = int(torch.cuda.max_memory_allocated(device))
    else:
        peak = None
    elapsed = time.perf_counter() - started
    cost: dict[str, Any] = {
        "wall_seconds": elapsed,
        "peak_gpu_bytes": peak,
        "n": n,
        "M": M,
        "repetition": repetition,
        "dtype": "float64",
        **totals,
        "gaussian_coordinate_work": equation.d
        * (totals["terminal_samples"] + totals["transition_samples"]),
    }
    return output, cost


def cpu_gpu_tape_check(
    equation: B1Equation,
    *,
    method: str,
    n: int,
    M: int,
    device: torch.device,
) -> dict[str, Any]:
    """Replay identical draws through NumPy and Torch implementations."""

    internal_method = MLP_METHODS.get(method, method)
    seed = stable_seed("B1", "tape", equation.d, equation.k, internal_method, n, M)
    rng = np.random.default_rng(seed)
    t = rng.uniform(0.0, equation.T, 5).astype(np.float64)
    x = rng.standard_normal((5, equation.d), dtype=np.float64)
    tape = DrawTape("record")
    seed_parts = [seed, equation.d, equation.k, n, M]
    numpy_solver = CorrectionBaselineMLP(
        equation,
        M,
        METHODS[internal_method],
        np.random.default_rng(np.random.SeedSequence(seed_parts)),
        aux_seed_sequence=np.random.SeedSequence(seed_parts + [0xD0B1E]),
        dose_rng=np.random.default_rng(np.random.SeedSequence(seed_parts + [0xD05E])),
        draw_tape=tape,
    )
    expected = numpy_solver.solve(n, t, x)
    replay = tape.replay()
    torch_solver = TorchMLP(
        equation,
        M,
        internal_method,
        device,
        seed_sequence=np.random.SeedSequence(seed_parts),
        aux_seed_sequence=np.random.SeedSequence(seed_parts + [0xD0B1E]),
        draw_source=TapeTorchDraws(device, replay),
    )
    with torch.no_grad():
        observed = torch_solver.solve(
            n,
            torch.as_tensor(t, dtype=torch.float64, device=device),
            torch.as_tensor(x, dtype=torch.float64, device=device),
        ).cpu().numpy()
    absolute = float(np.max(np.abs(observed - expected)))
    relative = absolute / max(float(np.max(np.abs(expected))), np.finfo(np.float64).tiny)
    return {
        "method": method,
        "n": n,
        "M": M,
        "draw_entries": len(tape.entries),
        "tape_exhausted": replay.exhausted,
        "max_absolute_error": absolute,
        "max_relative_error": relative,
        "passed": bool(replay.exhausted and relative <= 1e-10),
    }
