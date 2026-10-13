"""Paired MLP diagnostics for isolating learned-actor error."""

from __future__ import annotations

import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.gpu_mlp.solver import TorchMLP

from .actor_mlp import ActorMLP
from .common import stable_seed
from .equation import B1Equation


DIAGNOSTIC_METHODS = (
    "oracle_actor",
    "pathwise_baseline_no_actor",
)
PAIRED_SEED_NAMESPACE = "B1|actor_mlp"


class OracleActorMLP(ActorMLP):
    """Actor--MLP with the diagnostic oracle ``q(t,x)=z*(t,x)``.

    The inherited recursion evaluates :meth:`_actor_value` at every sampled
    child state.  Calling :class:`TorchMLP` directly here avoids constructing a
    learned network while retaining the exact Actor--MLP recursion and
    level-matched actor-driver baseline.
    """

    def __init__(
        self,
        equation: B1Equation,
        M: int,
        device: torch.device,
        *,
        seed_sequence: np.random.SeedSequence,
        aux_seed_sequence: np.random.SeedSequence,
    ) -> None:
        TorchMLP.__init__(
            self,
            equation,
            M,
            "path_base",
            device,
            seed_sequence=seed_sequence,
            aux_seed_sequence=aux_seed_sequence,
        )
        self.actor_calls = 0
        self.actor_points = 0

    def _actor_value(
        self, t: torch.Tensor, x: torch.Tensor, state: torch.Tensor
    ) -> torch.Tensor:
        q = self.equation.exact_z(t, x)
        z = state[..., 1:]
        self.actor_calls += 1
        self.actor_points += len(t)
        return -torch.sum(q * z, dim=-1) + 0.5 * torch.sum(q.square(), dim=-1)


def evaluate_actor_diagnostic(
    equation: B1Equation,
    t: np.ndarray,
    x: np.ndarray,
    *,
    method: str,
    n: int,
    M: int,
    repetition: int,
    device: torch.device,
    chunk_size: int | None = None,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Evaluate one diagnostic with learned-Actor common random numbers."""

    if method not in DIAGNOSTIC_METHODS:
        raise ValueError(f"unknown actor diagnostic {method!r}")
    times = np.asarray(t, dtype=np.float64).reshape(-1)
    states = np.asarray(x, dtype=np.float64).reshape(len(times), equation.d)
    chunk = len(times) if chunk_size is None else min(chunk_size, len(times))
    output = np.empty((len(times), equation.d + 1), dtype=np.float64)
    total_stats = {
        "terminal_samples": 0,
        "transition_samples": 0,
        "generator_calls": 0,
        "recursively_evaluated_states": 0,
        "recursive_calls": 0,
        "actor_calls": 0,
        "actor_points": 0,
    }
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    for chunk_index, start in enumerate(range(0, len(times), chunk)):
        stop = min(start + chunk, len(times))
        # This is intentionally identical to evaluate_actor_mlp so the learned,
        # oracle, and no-actor variants use common random numbers.
        seed = stable_seed(
            "B1", "actor_mlp", equation.d, equation.k, n, M, repetition, chunk_index
        )
        primary, auxiliary = np.random.SeedSequence(seed).spawn(2)
        if method == "oracle_actor":
            solver: TorchMLP = OracleActorMLP(
                equation,
                M,
                device,
                seed_sequence=primary,
                aux_seed_sequence=auxiliary,
            )
        else:
            solver = TorchMLP(
                equation,
                M,
                "path_base",
                device,
                seed_sequence=primary,
                aux_seed_sequence=auxiliary,
            )
        with torch.no_grad():
            result = solver.solve(
                n,
                torch.as_tensor(
                    times[start:stop], dtype=torch.float64, device=device
                ),
                torch.as_tensor(
                    states[start:stop], dtype=torch.float64, device=device
                ),
            )
        output[start:stop] = result.cpu().numpy()
        for key in total_stats:
            if hasattr(solver.stats, key):
                total_stats[key] += int(getattr(solver.stats, key))
            elif hasattr(solver, key):
                total_stats[key] += int(getattr(solver, key))
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        peak = int(torch.cuda.max_memory_allocated(device))
    else:
        peak = None
    total_stats.update(
        {
            "wall_seconds": time.perf_counter() - started,
            "peak_gpu_bytes": peak,
            "n": n,
            "M": M,
            "repetition": repetition,
            "dtype": "float64",
            "gaussian_coordinate_work": equation.d
            * (total_stats["terminal_samples"] + total_stats["transition_samples"]),
            "paired_seed_namespace": PAIRED_SEED_NAMESPACE,
            "oracle_reference_used_inside_solver": method == "oracle_actor",
            "driver": (
                "oracle actor-linearized quadratic driver"
                if method == "oracle_actor"
                else "original quadratic driver"
            ),
            "terminal_gradient": "pathwise explicit terminal derivative",
            "level_one_gradient_baseline": True,
        }
    )
    return output, total_stats

