"""GPU Full-History adapters for S05 and B6, preserving frozen recursion."""

from __future__ import annotations

import hashlib
import math
from pathlib import Path
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.baseline_n3.solver import METHODS
from invariant_region_mlp.experiments.double_estimator.solver import _metrics
from invariant_region_mlp.experiments.gpu_mlp.solver import (
    DrawSource,
    GPU_METHODS,
    NativeTorchDraws,
)
from invariant_region_mlp.experiments.nonridge.solver import (
    IMPLEMENTATION_REVISION as NONRIDGE_IMPLEMENTATION_REVISION,
    NonRidgeCorrectionBaselineMLP,
    NonRidgeTorchMLP,
    NonRidgeTorchStats,
    current_commit,
    gpu_chunk_size,
)

from .equations import B6Equation, fixed_points, make_b6_equation, make_s05_equation
from .protocol import (
    BASE_SEED,
    CALIBRATION_FREEZE_PATH,
    PREDICTION_ROOT,
    REVISION,
)


IMPLEMENTATION_REVISION = "highdim_nosym_gpu_v1"


def load_b6_selection(d: int) -> tuple[float, float, dict[str, Any]]:
    import json

    if not CALIBRATION_FREEZE_PATH.is_file():
        raise RuntimeError(
            "B6 calibration is not frozen; run calibrate.py and commit the freeze first"
        )
    payload = json.loads(CALIBRATION_FREEZE_PATH.read_text(encoding="utf-8"))
    row = payload["selections"][str(int(d))]
    return float(row["selected_beta"]), float(row["selected_gamma"]), payload


def make_equation(pde_id: str, d: int):
    if pde_id == "S05":
        return make_s05_equation(d)
    if pde_id == "B6":
        beta, gamma, _ = load_b6_selection(d)
        return make_b6_equation(d, beta, gamma)
    raise KeyError(pde_id)


class B6TorchEquation:
    """Native float64 PyTorch equation adapter for the frozen B6 parameters."""

    def __init__(self, equation: B6Equation, device: torch.device) -> None:
        self.numpy = equation
        self.device = device
        self.family = equation.family
        self.d = equation.d
        self.T = equation.T
        self.sigma = equation.sigma
        self.mu = torch.as_tensor(equation.mu, dtype=torch.float64, device=device)
        self.b = torch.tensor(equation.b, dtype=torch.float64, device=device)
        self.n2 = torch.tensor(equation.n2, dtype=torch.float64, device=device)
        self.log_count = math.log(equation.K)

    def _logits(self, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        tau = torch.clamp(self.T - t, min=0.0)
        return torch.matmul(x, self.b.T) + tau.unsqueeze(-1) * self.n2

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        return -(torch.logsumexp(torch.matmul(x, self.b.T), dim=-1) - self.log_count)

    def terminal_gradient_z(self, x: torch.Tensor) -> torch.Tensor:
        weights = torch.softmax(torch.matmul(x, self.b.T), dim=-1)
        return -self.sigma * torch.matmul(weights, self.b)

    def terminal_and_z(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        logits = torch.matmul(x, self.b.T)
        weights = torch.softmax(logits, dim=-1)
        value = -(torch.logsumexp(logits, dim=-1) - self.log_count)
        return value, -self.sigma * torch.matmul(weights, self.b)

    def exact_state(self, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        logits = self._logits(t, x)
        weights = torch.softmax(logits, dim=-1)
        value = -(torch.logsumexp(logits, dim=-1) - self.log_count)
        z = -self.sigma * torch.matmul(weights, self.b)
        return torch.cat((value.unsqueeze(-1), z), dim=-1)

    def exact_z(self, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        return self.exact_state(t, x)[..., 1:]

    @staticmethod
    def generator(state: torch.Tensor) -> torch.Tensor:
        return -0.5 * torch.sum(state[..., 1:].square(), dim=-1)

    @staticmethod
    def double_generator(first: torch.Tensor, second: torch.Tensor) -> torch.Tensor:
        return -0.5 * torch.sum(first[..., 1:] * second[..., 1:], dim=-1)


class B6TorchMLP(NonRidgeTorchMLP):
    """The imported GPU recursion with only its equation adapter replaced."""

    def __init__(
        self,
        equation: B6Equation,
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
        self.equation = B6TorchEquation(equation, device)
        self.numpy_equation = equation
        self.M = int(M)
        self.method_name = method_name
        self.method = METHODS[method_name]
        self.device = device
        self.time_beta_alpha = float(time_beta_alpha)
        self.draws: DrawSource = (
            NativeTorchDraws(device, seed_sequence)
            if draw_source is None
            else draw_source
        )
        self._shared_draw_source = draw_source is not None
        self._aux_seed_sequence = aux_seed_sequence
        self.stats = NonRidgeTorchStats()
        self.root_level1_integrand_variance: torch.Tensor | None = None


def make_torch_solver(
    equation: Any,
    M: int,
    method: str,
    device: torch.device,
    *,
    seed_sequence: np.random.SeedSequence,
    aux_seed_sequence: np.random.SeedSequence,
    draw_source: DrawSource | None = None,
):
    cls = B6TorchMLP if isinstance(equation, B6Equation) else NonRidgeTorchMLP
    return cls(
        equation,
        M,
        method,
        device,
        seed_sequence=seed_sequence,
        aux_seed_sequence=aux_seed_sequence,
        draw_source=draw_source,
    )


def _point_fingerprint(
    t: np.ndarray, x: np.ndarray, is_validation: np.ndarray
) -> str:
    digest = hashlib.sha256()
    digest.update(np.ascontiguousarray(t).tobytes())
    digest.update(np.ascontiguousarray(x).tobytes())
    digest.update(np.ascontiguousarray(is_validation).tobytes())
    return digest.hexdigest()


def run_fh_task(
    task: dict[str, Any],
    *,
    device: str | torch.device = "cuda",
    chunk_override: int | None = None,
    save_prediction: bool = True,
) -> dict[str, Any]:
    device_object = torch.device(device)
    pde_id, d = str(task["pde_id"]), int(task["d"])
    n, M = int(task["n"]), int(task["M"])
    repetition, method = int(task["repetition"]), str(task["method"])
    equation = make_equation(pde_id, d)
    t, x, is_validation = fixed_points(pde_id, d)
    truth = equation.exact_state(t, x)
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
            solver = make_torch_solver(
                equation,
                M,
                method,
                device_object,
                seed_sequence=np.random.SeedSequence(seed_parts),
                aux_seed_sequence=np.random.SeedSequence(seed_parts + [0xD0B1E]),
            )
            result = solver.solve(
                n,
                torch.tensor(t[begin:end].copy(), dtype=torch.float64, device=device_object),
                torch.tensor(x[begin:end].copy(), dtype=torch.float64, device=device_object),
                collect_root=True,
            )
            prediction[begin:end] = result.cpu().numpy()
            if solver.root_level1_integrand_variance is not None:
                root_variance[begin:end] = (
                    solver.root_level1_integrand_variance.cpu().numpy()
                )
            totals.merge(solver.stats)
    if device_object.type == "cuda":
        torch.cuda.synchronize(device_object)
    elapsed = time.perf_counter() - started
    mask = ~is_validation
    metrics = _metrics(prediction, truth, mask)
    validation = _metrics(prediction, truth, is_validation)
    value_error = prediction[mask, 0] - truth[mask, 0]
    gradient_error = prediction[mask, 1:] - truth[mask, 1:]
    finite_variance = root_variance[mask & np.isfinite(root_variance)]
    vectors = totals.terminal_samples + totals.transition_samples
    peak = (
        int(torch.cuda.max_memory_allocated(device_object))
        if device_object.type == "cuda"
        else None
    )
    prediction_path: Path | None = None
    if save_prediction:
        prediction_path = (
            PREDICTION_ROOT
            / "fh"
            / pde_id
            / f"d{d}"
            / f"n{n}_M{M}_{method}_rep{repetition:02d}.npz"
        )
        prediction_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            prediction_path,
            prediction=prediction,
            truth=truth,
            t=t,
            x=x,
            is_validation=is_validation,
            root_level1_variance=root_variance,
        )
    return {
        "schema_version": 1,
        "study_revision": REVISION,
        "implementation_revision": IMPLEMENTATION_REVISION,
        "inherited_nonridge_revision": NONRIDGE_IMPLEMENTATION_REVISION,
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
        "seed_scheme": (
            "torch.Generator seeded from "
            "SeedSequence([20261201,d,n,M,rep,chunk_index]).generate_state(2)"
        ),
        "dtype": "float64",
        "chunk_size": chunk,
        "point_count": len(t),
        "test_point_count": int(np.count_nonzero(mask)),
        "point_fingerprint": _point_fingerprint(t, x, is_validation),
        "skill": metrics["skill"],
        "value_rmse": metrics["value_rmse"],
        "value_relative_l2": metrics["value_relative_l2"],
        "value_bias": metrics["value_bias"],
        "gradient_relative_l2": metrics["gradient_relative_l2"],
        "gradient_vector_rmse": float(
            np.sqrt(np.mean(np.sum(gradient_error * gradient_error, axis=1)))
        ),
        "value_mae": float(np.mean(np.abs(value_error))),
        "truth_std": float(np.std(truth[mask, 0])),
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
        "prediction_path": None if prediction_path is None else str(prediction_path),
        "reused": False,
        "equation": equation.provenance(),
    }


__all__ = [
    "B6TorchEquation",
    "B6TorchMLP",
    "IMPLEMENTATION_REVISION",
    "NonRidgeCorrectionBaselineMLP",
    "load_b6_selection",
    "make_equation",
    "make_torch_solver",
    "run_fh_task",
]
