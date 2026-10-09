"""Batched float64 PyTorch implementation of P1/MR MLP recursion."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import math
import time
from typing import Any, Protocol

import numpy as np
import torch

from invariant_region_mlp.experiments.baseline_n3.solver import DrawTape, METHODS
from invariant_region_mlp.experiments.double_estimator.equations import (
    BASE_SEED,
    fixed_points,
    make_equation,
)
from invariant_region_mlp.experiments.double_estimator.solver import _metrics


IMPLEMENTATION_REVISION = "gpu_mlp_torch_v1"
GPU_METHODS = (
    "raw",
    "path",
    "double",
    "double_path",
    "path_base",
    "double_path_base",
)


def torch_seed(seed_sequence: np.random.SeedSequence) -> int:
    words = seed_sequence.generate_state(2, dtype=np.uint32)
    value = int(words[0]) | (int(words[1]) << 32)
    return value % (2**63 - 1)


class DrawSource(Protocol):
    def normal(self, shape: tuple[int, ...]) -> torch.Tensor: ...
    def uniform(self, shape: tuple[int, ...]) -> torch.Tensor: ...


class NativeTorchDraws:
    def __init__(self, device: torch.device, seed_sequence: np.random.SeedSequence) -> None:
        self.device = device
        self.generator = torch.Generator(device=device)
        self.generator.manual_seed(torch_seed(seed_sequence))

    def normal(self, shape: tuple[int, ...]) -> torch.Tensor:
        return torch.randn(
            shape,
            dtype=torch.float64,
            device=self.device,
            generator=self.generator,
        )

    def uniform(self, shape: tuple[int, ...]) -> torch.Tensor:
        return torch.rand(
            shape,
            dtype=torch.float64,
            device=self.device,
            generator=self.generator,
        )


class TapeTorchDraws:
    def __init__(self, device: torch.device, tape: DrawTape) -> None:
        if tape.mode != "replay":
            raise ValueError("Torch tape must be in replay mode")
        self.device = device
        self.tape = tape

    def normal(self, shape: tuple[int, ...]) -> torch.Tensor:
        return torch.as_tensor(
            self.tape.take("normal", shape), dtype=torch.float64, device=self.device
        )

    def uniform(self, shape: tuple[int, ...]) -> torch.Tensor:
        return torch.as_tensor(
            self.tape.take("uniform", shape), dtype=torch.float64, device=self.device
        )


@dataclass
class TorchStats:
    terminal_samples: int = 0
    transition_samples: int = 0
    generator_calls: int = 0
    recursively_evaluated_states: int = 0
    recursive_calls: int = 0

    def merge(self, other: "TorchStats") -> None:
        self.terminal_samples += other.terminal_samples
        self.transition_samples += other.transition_samples
        self.generator_calls += other.generator_calls
        self.recursively_evaluated_states += other.recursively_evaluated_states
        self.recursive_calls += other.recursive_calls


class TorchEquation:
    def __init__(self, equation: Any, device: torch.device) -> None:
        if equation.family not in {"ridge_lse", "multiridge"}:
            raise ValueError(f"unsupported GPU equation family {equation.family!r}")
        self.numpy = equation
        self.device = device
        self.family = equation.family
        self.d = int(equation.d)
        self.T = float(equation.T)
        self.sigma = float(equation.sigma)
        self.mu = torch.as_tensor(equation.mu, dtype=torch.float64, device=device)
        if self.family == "ridge_lse":
            self.w = torch.as_tensor(equation.w, dtype=torch.float64, device=device)
            self.lambdas = torch.as_tensor(
                equation.lambdas, dtype=torch.float64, device=device
            )
            self.log_count = math.log(len(equation.lambdas))
        else:
            self.A = torch.as_tensor(equation.A, dtype=torch.float64, device=device)
            self.logc = torch.as_tensor(equation.logc, dtype=torch.float64, device=device)
            self.n2 = torch.as_tensor(equation.n2, dtype=torch.float64, device=device)

    def _logits(self, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        h = self.T - t
        if self.family == "ridge_lse":
            s = torch.matmul(x, self.w)
            return s.unsqueeze(-1) * self.lambdas + h.unsqueeze(-1) * self.lambdas.square()
        return self.logc + torch.matmul(x, self.A.T) + h.unsqueeze(-1) * self.n2

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        if self.family == "ridge_lse":
            s = torch.matmul(x, self.w)
            logits = s.unsqueeze(-1) * self.lambdas
            return -(torch.logsumexp(logits, dim=-1) - self.log_count)
        logits = self.logc + torch.matmul(x, self.A.T)
        return -torch.logsumexp(logits, dim=-1)

    def exact_state(self, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        logits = self._logits(t, x)
        normalization = torch.logsumexp(logits, dim=-1, keepdim=True)
        weights = torch.exp(logits - normalization)
        if self.family == "ridge_lse":
            u = -(normalization[..., 0] - self.log_count)
            coefficient = torch.sum(weights * self.lambdas, dim=-1)
            z = -self.sigma * coefficient.unsqueeze(-1) * self.w
        else:
            u = -normalization[..., 0]
            z = -self.sigma * torch.matmul(weights, self.A)
        return torch.cat((u.unsqueeze(-1), z), dim=-1)

    def exact_z(self, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        return self.exact_state(t, x)[..., 1:]

    @staticmethod
    def generator(state: torch.Tensor) -> torch.Tensor:
        return -0.5 * torch.sum(state[..., 1:].square(), dim=-1)

    @staticmethod
    def double_generator(first: torch.Tensor, second: torch.Tensor) -> torch.Tensor:
        return -0.5 * torch.sum(first[..., 1:] * second[..., 1:], dim=-1)


class TorchMLP:
    def __init__(
        self,
        equation: Any,
        M: int,
        method_name: str,
        device: torch.device,
        *,
        seed_sequence: np.random.SeedSequence,
        aux_seed_sequence: np.random.SeedSequence,
        draw_source: DrawSource | None = None,
        time_beta_alpha: float = 0.5,
    ) -> None:
        if method_name not in GPU_METHODS:
            raise ValueError(method_name)
        self.equation = TorchEquation(equation, device)
        self.numpy_equation = equation
        self.M = int(M)
        self.method_name = method_name
        self.method = METHODS[method_name]
        self.device = device
        self.time_beta_alpha = float(time_beta_alpha)
        self.draws: DrawSource = (
            NativeTorchDraws(device, seed_sequence) if draw_source is None else draw_source
        )
        self._shared_draw_source = draw_source is not None
        self._aux_seed_sequence = aux_seed_sequence
        self.stats = TorchStats()

    def _terminal_estimate(
        self, n: int, t: torch.Tensor, x: torch.Tensor
    ) -> torch.Tensor:
        b, d = len(t), self.equation.d
        samples = self.M**n
        remaining = torch.clamp(self.equation.T - t, min=0.0)
        gx = self.equation.terminal(x)
        output = torch.zeros((b, d + 1), dtype=torch.float64, device=self.device)
        at_terminal = remaining <= 4.0 * torch.finfo(torch.float64).eps
        if torch.any(at_terminal):
            output[at_terminal] = self.equation.exact_state(t[at_terminal], x[at_terminal])
        active = ~at_terminal
        if not torch.any(active):
            return output
        xa = x[active]
        ha = remaining[active]
        ga = gx[active]
        normal = self.draws.normal((len(xa), samples, d))
        xt = (
            xa[:, None, :]
            + self.equation.mu * ha[:, None, None]
            + self.equation.sigma * torch.sqrt(ha)[:, None, None] * normal
        )
        gt = self.equation.terminal(xt)
        self.stats.terminal_samples += len(xa) * samples
        difference = gt - ga[:, None]
        output[active, 0] = ga + torch.mean(difference, dim=1)
        if self.method.pathwise_terminal:
            terminal_times = torch.full(
                terminal_points_shape := xt.shape[:-1],
                self.equation.T,
                dtype=torch.float64,
                device=self.device,
            )
            if terminal_times.shape != terminal_points_shape:  # pragma: no cover
                raise AssertionError("terminal time shape mismatch")
            terminal_z = self.equation.exact_z(terminal_times, xt)
            mean_z = torch.mean(terminal_z, dim=1)
            constant = torch.all(terminal_z == terminal_z[:, :1, :], dim=(1, 2))
            if torch.any(constant):
                mean_z[constant] = terminal_z[constant, 0, :]
            output[active, 1:] = mean_z
        else:
            weights = normal / torch.sqrt(ha)[:, None, None]
            output[active, 1:] = torch.mean(
                difference[:, :, None] * weights, dim=1
            )
        return output

    def _independent_solve(
        self, n: int, t: torch.Tensor, x: torch.Tensor
    ) -> torch.Tensor:
        bundle = self._aux_seed_sequence.spawn(1)[0]
        primary_seed, auxiliary_seed, _dose_seed = bundle.spawn(3)
        child = TorchMLP(
            self.numpy_equation,
            self.M,
            self.method_name,
            self.device,
            seed_sequence=primary_seed,
            aux_seed_sequence=auxiliary_seed,
            draw_source=self.draws if self._shared_draw_source else None,
            time_beta_alpha=self.time_beta_alpha,
        )
        result = child.solve(n, t, x)
        self.stats.merge(child.stats)
        return result

    def solve(self, n: int, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        t = t.to(dtype=torch.float64, device=self.device).reshape(-1)
        x = x.to(dtype=torch.float64, device=self.device).reshape(len(t), self.equation.d)
        b = len(t)
        self.stats.recursive_calls += 1
        self.stats.recursively_evaluated_states += b
        if n <= 0:
            return torch.zeros(
                (b, self.equation.d + 1), dtype=torch.float64, device=self.device
            )
        output = self._terminal_estimate(n, t, x)
        baseline: torch.Tensor | None = None
        if n >= 2 and self.method.baseline_kind != "none":
            state = (
                output
                if self.method.baseline_kind == "terminal"
                else self.equation.exact_state(t, x)
            )
            baseline = self.equation.generator(state)
            self.stats.generator_calls += b

        remaining = self.equation.T - t
        alpha = self.time_beta_alpha
        tiny = torch.finfo(torch.float64).tiny
        for level in range(1, n):
            siblings = self.M ** (n - level)
            uniform = self.draws.uniform((b, siblings))
            r = torch.clamp(uniform, min=tiny) ** (1.0 / alpha)
            dt = remaining[:, None] * r
            normal = self.draws.normal((b, siblings, self.equation.d))
            self.stats.transition_samples += b * siblings
            xr = (
                x[:, None, :]
                + self.equation.mu * dt[:, :, None]
                + self.equation.sigma * torch.sqrt(dt)[:, :, None] * normal
            )
            tr = t[:, None] + dt
            flat_t = tr.reshape(-1)
            flat_x = xr.reshape(-1, self.equation.d)
            high_first = self.solve(level, flat_t, flat_x).reshape(
                b, siblings, self.equation.d + 1
            )
            if self.method.double_generator:
                high_second = self._independent_solve(level, flat_t, flat_x).reshape(
                    b, siblings, self.equation.d + 1
                )
                difference = self.equation.double_generator(high_first, high_second)
            else:
                difference = self.equation.generator(high_first)
            self.stats.generator_calls += b * siblings

            if level > 1:
                low_first = self.solve(level - 1, flat_t, flat_x).reshape(
                    b, siblings, self.equation.d + 1
                )
                if self.method.double_generator:
                    low_second = self._independent_solve(
                        level - 1, flat_t, flat_x
                    ).reshape(b, siblings, self.equation.d + 1)
                    difference = difference - self.equation.double_generator(
                        low_first, low_second
                    )
                else:
                    difference = difference - self.equation.generator(low_first)
                self.stats.generator_calls += b * siblings

            importance = remaining[:, None] * r ** (1.0 - alpha) / alpha
            value_increment = torch.mean(importance * difference, dim=1)
            gradient_difference = difference
            if level == 1 and baseline is not None:
                gradient_difference = difference - baseline[:, None]
            ebl = normal / torch.sqrt(torch.clamp(dt, min=tiny))[:, :, None]
            gradient_increment = torch.mean(
                importance[:, :, None]
                * gradient_difference[:, :, None]
                * ebl,
                dim=1,
            )
            output[:, 0] += value_increment
            output[:, 1:] += gradient_increment
        return output


def gpu_chunk_size(d: int, n: int, M: int, *, target_coordinates: int = 20_000_000) -> int:
    root_coordinates = max(1, int(d) * int(M) ** int(n))
    return max(1, min(64, target_coordinates // root_coordinates))


def run_gpu_task(
    task: dict[str, Any],
    *,
    device: str | torch.device = "cuda",
    chunk_override: int | None = None,
) -> dict[str, Any]:
    device_object = torch.device(device)
    pde_id, d = str(task["pde_id"]), int(task["d"])
    n, M = int(task["n"]), int(task["M"])
    repetition, method = int(task["repetition"]), str(task["method"])
    equation = make_equation(pde_id, d)
    t, x, is_validation = fixed_points(pde_id, d)
    truth = equation.exact_state(t, x)
    prediction = np.empty_like(truth)
    chunk = chunk_override or gpu_chunk_size(d, n, M)
    totals = TorchStats()
    if device_object.type == "cuda":
        torch.cuda.synchronize(device_object)
        torch.cuda.reset_peak_memory_stats(device_object)
    started = time.perf_counter()
    with torch.no_grad():
        for chunk_index, begin in enumerate(range(0, len(t), chunk)):
            end = min(begin + chunk, len(t))
            seed_parts = [BASE_SEED, d, n, M, repetition, chunk_index]
            solver = TorchMLP(
                equation,
                M,
                method,
                device_object,
                seed_sequence=np.random.SeedSequence(seed_parts),
                aux_seed_sequence=np.random.SeedSequence(seed_parts + [0xD0B1E]),
            )
            tt = torch.tensor(t[begin:end], dtype=torch.float64, device=device_object)
            xx = torch.tensor(x[begin:end], dtype=torch.float64, device=device_object)
            result = solver.solve(n, tt, xx)
            prediction[begin:end] = result.cpu().numpy()
            totals.merge(solver.stats)
    if device_object.type == "cuda":
        torch.cuda.synchronize(device_object)
    elapsed = time.perf_counter() - started
    metrics = _metrics(prediction, truth, ~is_validation)
    validation = _metrics(prediction, truth, is_validation)
    point_hash = hashlib.sha256()
    point_hash.update(np.ascontiguousarray(t).tobytes())
    point_hash.update(np.ascontiguousarray(x).tobytes())
    point_hash.update(np.ascontiguousarray(is_validation).tobytes())
    vectors = totals.terminal_samples + totals.transition_samples
    peak = (
        int(torch.cuda.max_memory_allocated(device_object))
        if device_object.type == "cuda"
        else None
    )
    return {
        "implementation_revision": IMPLEMENTATION_REVISION,
        "backend": "torch_gpu" if device_object.type == "cuda" else "torch_cpu",
        "device": str(device_object),
        "pde_id": pde_id,
        "family": equation.family,
        "dimension": d,
        "n": n,
        "M": M,
        "method": method,
        "repetition": repetition,
        "base_seed": BASE_SEED,
        "seed_scheme": "torch.Generator seeded from SeedSequence([20261201,d,n,M,rep,chunk_index]).generate_state(2)",
        "paired_with_cpu": False,
        "dtype": "float64",
        "chunk_size": chunk,
        "point_count": len(t),
        "point_fingerprint": point_hash.hexdigest(),
        "skill": metrics["skill"],
        "value_rmse": metrics["value_rmse"],
        "root_gradient_skill": metrics["gradient_relative_l2"],
        "validation_skill": validation["skill"],
        "generator_calls": totals.generator_calls,
        "gaussian_vectors_drawn": vectors,
        "gaussian_coordinates_drawn": vectors * d,
        "terminal_samples": totals.terminal_samples,
        "transition_samples": totals.transition_samples,
        "recursive_calls": totals.recursive_calls,
        "recursively_evaluated_states": totals.recursively_evaluated_states,
        "wall_time_seconds": elapsed,
        "peak_gpu_memory_bytes": peak,
        "nonfinite_state_count": metrics["nonfinite_state_count"],
    }


__all__ = [
    "GPU_METHODS",
    "IMPLEMENTATION_REVISION",
    "NativeTorchDraws",
    "TapeTorchDraws",
    "TorchEquation",
    "TorchMLP",
    "gpu_chunk_size",
    "run_gpu_task",
    "torch_seed",
]
