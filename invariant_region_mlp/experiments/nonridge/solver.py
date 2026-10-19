"""CPU validation adapter and native-float64 GPU solver for non-ridge HJBs."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import hashlib
import math
from pathlib import Path
import subprocess
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.baseline_n3.solver import (
    CorrectionBaselineMLP,
    METHODS,
)
from invariant_region_mlp.experiments.double_estimator.equations import BASE_SEED
from invariant_region_mlp.experiments.double_estimator.solver import _metrics
from invariant_region_mlp.experiments.gpu_mlp.solver import (
    GPU_METHODS,
    DrawSource,
    NativeTorchDraws,
    TorchMLP,
    TorchStats,
)

from .equations import (
    PROFILE_SLOPES,
    SpectralSeparableHJB,
    fixed_points,
    hermite_rule,
    make_equation,
)


IMPLEMENTATION_REVISION = "nonridge_gpu_v1"


def _repository_root() -> Path:
    return Path(__file__).resolve().parents[3]


def current_commit() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=_repository_root(),
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def torch_dct_ortho(x: torch.Tensor) -> torch.Tensor:
    """Orthonormal type-II DCT along the last axis (Makhoul FFT form)."""

    n = x.shape[-1]
    reordered = torch.cat((x[..., ::2], torch.flip(x[..., 1::2], dims=(-1,))), dim=-1)
    spectrum = torch.fft.fft(reordered, dim=-1)
    theta = (
        torch.arange(n, dtype=x.dtype, device=x.device)
        * (math.pi / (2.0 * n))
    )
    output = spectrum.real * torch.cos(theta) + spectrum.imag * torch.sin(theta)
    output = output.clone()
    output[..., 0] /= math.sqrt(n)
    if n > 1:
        output[..., 1:] /= math.sqrt(n / 2.0)
    return output


def torch_idct_ortho(x: torch.Tensor) -> torch.Tensor:
    """Inverse of :func:`torch_dct_ortho` (orthonormal type-II IDCT)."""

    n = x.shape[-1]
    coefficients = x.clone()
    coefficients[..., 0] *= math.sqrt(n)
    if n > 1:
        coefficients[..., 1:] *= math.sqrt(n / 2.0)
    theta = (
        torch.arange(n, dtype=x.dtype, device=x.device)
        * (math.pi / (2.0 * n))
    )
    imaginary_seed = torch.cat(
        (
            torch.zeros_like(coefficients[..., :1]),
            -torch.flip(coefficients, dims=(-1,))[..., :-1],
        ),
        dim=-1,
    )
    real = coefficients * torch.cos(theta) - imaginary_seed * torch.sin(theta)
    imaginary = coefficients * torch.sin(theta) + imaginary_seed * torch.cos(theta)
    reordered = torch.fft.ifft(torch.complex(real, imaginary), dim=-1).real
    output = torch.empty_like(reordered)
    even_count = n - n // 2
    output[..., ::2] = reordered[..., :even_count]
    output[..., 1::2] = torch.flip(reordered, dims=(-1,))[..., : n // 2]
    return output


class SpectralTorchEquation:
    """Torch equation wrapper using an FFT DCT rather than a dense rotation."""

    def __init__(self, equation: SpectralSeparableHJB, device: torch.device) -> None:
        self.numpy = equation
        self.device = device
        self.family = equation.family
        self.d = equation.d
        self.T = equation.T
        self.sigma = equation.sigma
        self.mu = torch.as_tensor(equation.mu, dtype=torch.float64, device=device)
        self.c = torch.as_tensor(equation.c, dtype=torch.float64, device=device)
        self.signs = torch.as_tensor(equation.signs, dtype=torch.float64, device=device)
        self.slopes = torch.as_tensor(PROFILE_SLOPES, dtype=torch.float64, device=device)
        nodes, _weights, log_weights = hermite_rule(equation.quadrature_order)
        self.nodes = torch.tensor(nodes.copy(), dtype=torch.float64, device=device)
        self.log_weights = torch.tensor(
            log_weights.copy(), dtype=torch.float64, device=device
        )

    def rotate(self, x: torch.Tensor) -> torch.Tensor:
        return torch_dct_ortho(x * self.signs)

    def inverse_rotate(self, y: torch.Tensor) -> torch.Tensor:
        return self.signs * torch_idct_ortho(y)

    def profile(self, y: torch.Tensor) -> torch.Tensor:
        return math.log(3.0) - torch.logsumexp(
            y.unsqueeze(-1) * self.slopes, dim=-1
        )

    def profile_prime(self, y: torch.Tensor) -> torch.Tensor:
        logits = y.unsqueeze(-1) * self.slopes
        weights = torch.softmax(logits, dim=-1)
        return -torch.sum(weights * self.slopes, dim=-1)

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        y = self.rotate(x)
        return torch.sum(self.c * self.profile(y), dim=-1)

    def terminal_gradient_z(self, x: torch.Tensor) -> torch.Tensor:
        y = self.rotate(x)
        return self.sigma * self.inverse_rotate(self.c * self.profile_prime(y))

    def terminal_and_z(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        y = self.rotate(x)
        value = torch.sum(self.c * self.profile(y), dim=-1)
        z = self.sigma * self.inverse_rotate(self.c * self.profile_prime(y))
        return value, z

    def exact_state(self, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        leading = x.shape[:-1]
        time = torch.broadcast_to(t, leading)
        y = self.rotate(x)
        tau = torch.clamp(self.T - time, min=0.0)
        shifted = (
            y.unsqueeze(-1)
            + 2.0 * torch.sqrt(tau).unsqueeze(-1).unsqueeze(-1) * self.nodes
        )
        log_integrand = (
            -self.c.unsqueeze(-1) * self.profile(shifted) + self.log_weights
        )
        log_denominator = torch.logsumexp(log_integrand, dim=-1)
        tilted = torch.exp(log_integrand - log_denominator.unsqueeze(-1))
        value = -torch.sum(log_denominator, dim=-1)
        derivative = torch.sum(
            tilted
            * self.c.unsqueeze(-1)
            * self.profile_prime(shifted),
            dim=-1,
        )
        z = self.sigma * self.inverse_rotate(derivative)
        return torch.cat((value.unsqueeze(-1), z), dim=-1)

    def exact_z(self, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        if torch.all(
            torch.abs(t - self.T) <= 4.0 * torch.finfo(torch.float64).eps
        ):
            return self.terminal_gradient_z(x)
        return self.exact_state(t, x)[..., 1:]

    @staticmethod
    def generator(state: torch.Tensor) -> torch.Tensor:
        return -0.5 * torch.sum(state[..., 1:].square(), dim=-1)

    @staticmethod
    def double_generator(first: torch.Tensor, second: torch.Tensor) -> torch.Tensor:
        return -0.5 * torch.sum(first[..., 1:] * second[..., 1:], dim=-1)


class NonRidgeCorrectionBaselineMLP(CorrectionBaselineMLP):
    """Frozen CPU solver with only the registered quadratic double driver supplied."""

    def _double_value(self, first: np.ndarray, second: np.ndarray) -> np.ndarray:
        return -0.5 * np.sum(first[..., 1:] * second[..., 1:], axis=-1)


@dataclass
class NonRidgeTorchStats(TorchStats):
    nonfinite_generator_values: int = 0

    def merge(self, other: TorchStats) -> None:
        super().merge(other)
        self.nonfinite_generator_values += int(
            getattr(other, "nonfinite_generator_values", 0)
        )


class NonRidgeTorchMLP(TorchMLP):
    """Frozen GPU recursion extended to the spectral equation and V1 diagnostic."""

    def __init__(
        self,
        equation: SpectralSeparableHJB,
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
        self.equation = SpectralTorchEquation(equation, device)
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
        self.stats = NonRidgeTorchStats()
        self.root_level1_integrand_variance: torch.Tensor | None = None

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
            terminal_value, terminal_z = self.equation.terminal_and_z(x[at_terminal])
            output[at_terminal, 0] = terminal_value
            output[at_terminal, 1:] = terminal_z
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
        if self.method.pathwise_terminal:
            gt, terminal_z = self.equation.terminal_and_z(xt)
            output[active, 1:] = torch.mean(terminal_z, dim=1)
        else:
            gt = self.equation.terminal(xt)
            weights = normal / torch.sqrt(ha)[:, None, None]
            output[active, 1:] = torch.mean(
                (gt - ga[:, None])[:, :, None] * weights,
                dim=1,
            )
        self.stats.terminal_samples += len(xa) * samples
        output[active, 0] = ga + torch.mean(gt - ga[:, None], dim=1)
        return output

    def _independent_solve(
        self, n: int, t: torch.Tensor, x: torch.Tensor
    ) -> torch.Tensor:
        bundle = self._aux_seed_sequence.spawn(1)[0]
        primary_seed, auxiliary_seed, _dose_seed = bundle.spawn(3)
        child = type(self)(
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

    def _plain_generator(self, state: torch.Tensor) -> torch.Tensor:
        result = self.equation.generator(state)
        self.stats.generator_calls += result.numel()
        self.stats.nonfinite_generator_values += int(
            torch.count_nonzero(~torch.isfinite(result)).item()
        )
        return result

    def _double_generator(
        self, first: torch.Tensor, second: torch.Tensor
    ) -> torch.Tensor:
        result = self.equation.double_generator(first, second)
        self.stats.generator_calls += result.numel()
        self.stats.nonfinite_generator_values += int(
            torch.count_nonzero(~torch.isfinite(result)).item()
        )
        return result

    def solve(
        self,
        n: int,
        t: torch.Tensor,
        x: torch.Tensor,
        *,
        collect_root: bool = False,
    ) -> torch.Tensor:
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
        if collect_root:
            self.root_level1_integrand_variance = torch.full(
                (b,), float("nan"), dtype=torch.float64, device=self.device
            )
        baseline: torch.Tensor | None = None
        if n >= 2 and self.method.baseline_kind != "none":
            state = (
                output
                if self.method.baseline_kind == "terminal"
                else self.equation.exact_state(t, x)
            )
            baseline = self._plain_generator(state)

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
                difference = self._double_generator(high_first, high_second)
            else:
                difference = self._plain_generator(high_first)

            if level > 1:
                low_first = self.solve(level - 1, flat_t, flat_x).reshape(
                    b, siblings, self.equation.d + 1
                )
                if self.method.double_generator:
                    low_second = self._independent_solve(
                        level - 1, flat_t, flat_x
                    ).reshape(b, siblings, self.equation.d + 1)
                    difference = difference - self._double_generator(
                        low_first, low_second
                    )
                else:
                    difference = difference - self._plain_generator(low_first)

            importance = remaining[:, None] * r ** (1.0 - alpha) / alpha
            value_increment = torch.mean(importance * difference, dim=1)
            gradient_difference = difference
            if level == 1 and baseline is not None:
                gradient_difference = difference - baseline[:, None]
            ebl = normal / torch.sqrt(torch.clamp(dt, min=tiny))[:, :, None]
            gradient_integrand = (
                importance[:, :, None]
                * gradient_difference[:, :, None]
                * ebl
            )
            gradient_increment = torch.mean(gradient_integrand, dim=1)
            output[:, 0] += value_increment
            output[:, 1:] += gradient_increment
            if collect_root and level == 1:
                assert self.root_level1_integrand_variance is not None
                correction = 1 if siblings > 1 else 0
                self.root_level1_integrand_variance[:] = torch.mean(
                    torch.var(gradient_integrand, dim=1, correction=correction),
                    dim=1,
                )
        return output


def gpu_chunk_size(
    d: int,
    n: int,
    M: int,
    *,
    target_coordinates: int = 20_000_000,
) -> int:
    root_coordinates = max(1, int(d) * int(M) ** int(n))
    return max(1, min(64, target_coordinates // root_coordinates))


@lru_cache(maxsize=None)
def _reference_data(
    pde_id: str, d: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Cache immutable points and the quadrature truth within one run process."""

    equation = make_equation(pde_id, d)
    t, x, is_validation = fixed_points(pde_id, d)
    truth = equation.exact_state(t, x)
    truth.setflags(write=False)
    return t, x, is_validation, truth


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
    t, x, is_validation, truth = _reference_data(pde_id, d)
    prediction = np.empty_like(truth)
    root_variance = np.full(len(t), np.nan, dtype=np.float64)
    chunk = chunk_override or gpu_chunk_size(d, n, M)
    totals = NonRidgeTorchStats()
    if device_object.type == "cuda":
        torch.cuda.synchronize(device_object)
        torch.cuda.reset_peak_memory_stats(device_object)
    started = time.perf_counter()
    with torch.no_grad():
        for chunk_index, begin in enumerate(range(0, len(t), chunk)):
            end = min(begin + chunk, len(t))
            seed_parts = [BASE_SEED, d, n, M, repetition, chunk_index]
            solver = NonRidgeTorchMLP(
                equation,
                M,
                method,
                device_object,
                seed_sequence=np.random.SeedSequence(seed_parts),
                aux_seed_sequence=np.random.SeedSequence(seed_parts + [0xD0B1E]),
            )
            tt = torch.tensor(
                t[begin:end].copy(), dtype=torch.float64, device=device_object
            )
            xx = torch.tensor(
                x[begin:end].copy(), dtype=torch.float64, device=device_object
            )
            result = solver.solve(n, tt, xx, collect_root=True)
            prediction[begin:end] = result.cpu().numpy()
            if solver.root_level1_integrand_variance is not None:
                root_variance[begin:end] = (
                    solver.root_level1_integrand_variance.cpu().numpy()
                )
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
    finite_variance = root_variance[(~is_validation) & np.isfinite(root_variance)]
    peak = (
        int(torch.cuda.max_memory_allocated(device_object))
        if device_object.type == "cuda"
        else None
    )
    return {
        "schema_version": 1,
        "implementation_revision": IMPLEMENTATION_REVISION,
        "code_commit": current_commit(),
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
        "methods_paired_within_cell": True,
        "dtype": "float64",
        "chunk_size": chunk,
        "point_count": len(t),
        "point_fingerprint": point_hash.hexdigest(),
        "skill": metrics["skill"],
        "value_rmse": metrics["value_rmse"],
        "value_relative_l2": metrics["value_relative_l2"],
        "value_bias": metrics["value_bias"],
        "root_gradient_skill": metrics["gradient_relative_l2"],
        "validation_skill": validation["skill"],
        "generator_calls": totals.generator_calls,
        "gaussian_vectors_drawn": vectors,
        "gaussian_coordinates_drawn": vectors * d,
        "terminal_samples": totals.terminal_samples,
        "transition_samples": totals.transition_samples,
        "recursive_calls": totals.recursive_calls,
        "recursively_evaluated_states": totals.recursively_evaluated_states,
        "root_level1_gradient_integrand_variance": (
            float(np.mean(finite_variance)) if len(finite_variance) else None
        ),
        "wall_time_seconds": elapsed,
        "peak_gpu_memory_bytes": peak,
        "nonfinite_state_count": metrics["nonfinite_state_count"],
        "nonfinite_generator_count": totals.nonfinite_generator_values,
        "equation": equation.provenance(),
    }


__all__ = [
    "IMPLEMENTATION_REVISION",
    "NonRidgeCorrectionBaselineMLP",
    "NonRidgeTorchMLP",
    "SpectralTorchEquation",
    "gpu_chunk_size",
    "run_gpu_task",
    "torch_dct_ortho",
    "torch_idct_ortho",
]
