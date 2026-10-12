"""Complete Actor--MLP inside the established full-history recursion."""

from __future__ import annotations

import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.gpu_mlp.solver import TorchMLP

from .actor import ActorNet
from .common import stable_seed
from .equation import B1Equation


class ActorMLP(TorchMLP):
    """Pathwise terminal, level-1 baseline, and frozen actor driver."""

    def __init__(
        self,
        equation: B1Equation,
        M: int,
        actor: ActorNet,
        device: torch.device,
        *,
        seed_sequence: np.random.SeedSequence,
        aux_seed_sequence: np.random.SeedSequence,
    ) -> None:
        super().__init__(
            equation,
            M,
            "path_base",
            device,
            seed_sequence=seed_sequence,
            aux_seed_sequence=aux_seed_sequence,
        )
        self.actor = actor
        self.actor.eval()
        self.actor_calls = 0
        self.actor_points = 0

    def _actor_value(
        self, t: torch.Tensor, x: torch.Tensor, state: torch.Tensor
    ) -> torch.Tensor:
        actor_dtype = next(self.actor.parameters()).dtype
        q = self.actor.predict(t.to(actor_dtype), x.to(actor_dtype), self.equation.T)
        q = q.to(dtype=torch.float64)
        self.actor_calls += 1
        self.actor_points += len(t)
        z = state[..., 1:]
        return -torch.sum(q * z, dim=-1) + 0.5 * torch.sum(q.square(), dim=-1)

    def solve(self, n: int, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        t = t.to(dtype=torch.float64, device=self.device).reshape(-1)
        x = x.to(dtype=torch.float64, device=self.device).reshape(len(t), self.equation.d)
        batch = len(t)
        self.stats.recursive_calls += 1
        self.stats.recursively_evaluated_states += batch
        if n <= 0:
            return torch.zeros(
                (batch, self.equation.d + 1), dtype=torch.float64, device=self.device
            )

        output = self._terminal_estimate(n, t, x)
        baseline: torch.Tensor | None = None
        if n >= 2:
            # The terminal estimate uses randomness drawn before, and therefore
            # independently of every subsequently drawn sibling correction.
            baseline = self._actor_value(t, x, output)
            self.stats.generator_calls += batch

        remaining = self.equation.T - t
        alpha = self.time_beta_alpha
        tiny = torch.finfo(torch.float64).tiny
        for level in range(1, n):
            siblings = self.M ** (n - level)
            uniform = self.draws.uniform((batch, siblings))
            r = torch.clamp(uniform, min=tiny) ** (1.0 / alpha)
            dt = remaining[:, None] * r
            normal = self.draws.normal((batch, siblings, self.equation.d))
            self.stats.transition_samples += batch * siblings
            xr = (
                x[:, None, :]
                + self.equation.mu * dt[:, :, None]
                + self.equation.sigma * torch.sqrt(dt)[:, :, None] * normal
            )
            tr = t[:, None] + dt
            flat_t = tr.reshape(-1)
            flat_x = xr.reshape(-1, self.equation.d)
            high = self.solve(level, flat_t, flat_x).reshape(
                batch, siblings, self.equation.d + 1
            )
            difference = self._actor_value(
                flat_t, flat_x, high.reshape(-1, self.equation.d + 1)
            ).reshape(batch, siblings)
            self.stats.generator_calls += batch * siblings
            if level > 1:
                low = self.solve(level - 1, flat_t, flat_x).reshape(
                    batch, siblings, self.equation.d + 1
                )
                low_value = self._actor_value(
                    flat_t, flat_x, low.reshape(-1, self.equation.d + 1)
                ).reshape(batch, siblings)
                difference = difference - low_value
                self.stats.generator_calls += batch * siblings

            importance = remaining[:, None] * r ** (1.0 - alpha) / alpha
            output[:, 0] += torch.mean(importance * difference, dim=1)
            gradient_difference = difference
            if level == 1 and baseline is not None:
                gradient_difference = difference - baseline[:, None]
            brownian_weight = normal / torch.sqrt(
                torch.clamp(dt, min=tiny)
            )[:, :, None]
            output[:, 1:] += torch.mean(
                importance[:, :, None]
                * gradient_difference[:, :, None]
                * brownian_weight,
                dim=1,
            )
        return output


def evaluate_actor_mlp(
    equation: B1Equation,
    actor: ActorNet,
    t: np.ndarray,
    x: np.ndarray,
    *,
    n: int,
    M: int,
    repetition: int,
    device: torch.device,
    chunk_size: int | None = None,
) -> tuple[np.ndarray, dict[str, Any]]:
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
        seed = stable_seed(
            "B1", "actor_mlp", equation.d, equation.k, n, M, repetition, chunk_index
        )
        primary, auxiliary = np.random.SeedSequence(seed).spawn(2)
        solver = ActorMLP(
            equation,
            M,
            actor,
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
        for key in total_stats:
            total_stats[key] += (
                int(getattr(solver.stats, key))
                if hasattr(solver.stats, key)
                else int(getattr(solver, key))
            )
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
            "dtype": "float64 solver / float32 frozen actor",
            "gaussian_coordinate_work": equation.d
            * (total_stats["terminal_samples"] + total_stats["transition_samples"]),
        }
    )
    return output, total_stats


def actor_linearization_identity(
    actor: ActorNet, equation: B1Equation, device: torch.device
) -> float:
    seed = stable_seed("B1", "actor_identity", equation.d, equation.k)
    rng = np.random.default_rng(seed)
    t = torch.as_tensor(rng.uniform(0.0, equation.T, 32), dtype=torch.float32, device=device)
    x = torch.as_tensor(rng.standard_normal((32, equation.d)), dtype=torch.float32, device=device)
    z = torch.as_tensor(rng.standard_normal((32, equation.d)), dtype=torch.float64, device=device)
    with torch.no_grad():
        # The recursive solver converts actor output to float64 before its
        # dot products.  Check the identity in that same arithmetic; separate
        # float32 reductions can differ by several ulps when d=400.
        q = actor.predict(t, x, equation.T).to(torch.float64)
        actor_driver = -torch.sum(q * z, dim=-1) + 0.5 * torch.sum(q.square(), dim=-1)
        original = -0.5 * torch.sum(z.square(), dim=-1)
        expected = 0.5 * torch.sum((z - q).square(), dim=-1)
    return float(torch.max(torch.abs((actor_driver - original) - expected)).cpu())
