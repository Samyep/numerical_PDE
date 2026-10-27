"""Source-faithful scalar-Y Deep BSDE adapted to arbitrary initial times."""

from __future__ import annotations

from dataclasses import asdict
import math
from pathlib import Path
import random
import tempfile
import time
from typing import Any, Iterable

import numpy as np
import torch

from invariant_region_mlp.experiments.main_highdim_expansion.benchmarks import (
    make_benchmark,
    metadata as benchmark_metadata,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.common import (
    seed_everything,
    synchronize,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.deep_bsde_port import (
    DeepBSDEConfig,
    OfficialDeepBSDE,
    _sample_dw,
    clipped_terminal_loss,
)

from .protocol import (
    CHECKPOINT_DIR,
    DIMENSION,
    HIDDEN_WIDTHS,
    INPUT_DIR,
    LEARNING_RATE,
    RAW_DIR,
    REVISION,
    TERMINAL_LOSS_CLIP,
    TRAJECTORY_BY_ID,
    VALIDATION_SIZE,
    Y_INIT_RANGE,
    TrajectorySpec,
    array_sha256,
    canonical_sha256,
    ensure_directories,
    file_sha256,
    stable_seed,
    write_json_atomic,
)


RESUME_INTERVAL = 250


class VariableStartDeepBSDE(OfficialDeepBSDE):
    """Official time-slice architecture on [t_start,T] from fixed x_start."""

    def __init__(
        self,
        benchmark: Any,
        config: DeepBSDEConfig,
        x_start: torch.Tensor,
        t_start: float,
    ) -> None:
        super().__init__(benchmark, config)
        fixed = torch.as_tensor(x_start, dtype=torch.float64).reshape(1, benchmark.d)
        self.register_buffer("fixed_x_start", fixed)
        self.t_start = float(t_start)
        expected_horizon = float(benchmark.T) - self.t_start
        if not 0.0 <= self.t_start < float(benchmark.T):
            raise ValueError((self.t_start, benchmark.T))
        if not math.isclose(config.total_time, expected_horizon, rel_tol=0.0, abs_tol=2e-15):
            raise ValueError((config.total_time, expected_horizon))

    @property
    def dt(self) -> float:
        return self.config.total_time / self.config.intervals

    def driver(self, y: torch.Tensor, z: torch.Tensor) -> torch.Tensor:
        return self.benchmark.driver(y[:, 0], z).reshape(-1, 1)

    def rollout(self, dw: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        batch = len(dw)
        dt = self.dt
        x = self.fixed_x_start.to(dw).expand(batch, self.benchmark.d).clone()
        y = self.y_init.reshape(1, 1).expand(batch, 1)
        z = self.z_init.expand(batch, self.benchmark.d)
        for index in range(self.config.intervals - 1):
            y = y - dt * self.driver(y, z) + (z * dw[:, index]).sum(1, keepdim=True)
            absolute_time = torch.full(
                (batch,), self.t_start + index * dt, dtype=dw.dtype, device=dw.device
            )
            x = x + self.benchmark.base_drift(absolute_time, x) * dt + self.benchmark.sigma * dw[:, index]
            z = self.subnets[index](x) / self.benchmark.d
        y = y - dt * self.driver(y, z) + (z * dw[:, -1]).sum(1, keepdim=True)
        absolute_time = torch.full(
            (batch,),
            self.t_start + (self.config.intervals - 1) * dt,
            dtype=dw.dtype,
            device=dw.device,
        )
        x = x + self.benchmark.base_drift(absolute_time, x) * dt + self.benchmark.sigma * dw[:, -1]
        return y[:, 0], x


def require_h100(device: torch.device) -> str:
    if device.type != "cuda" or not torch.cuda.is_available():
        raise RuntimeError("formal Deep BSDE computation requires CUDA")
    name = torch.cuda.get_device_name(device)
    if "H100" not in name.upper():
        raise RuntimeError(f"formal timing requires NVIDIA H100, found {name!r}")
    return name


def deep_config(benchmark: Any, t_start: float, spec: TrajectorySpec, *, updates: int | None = None) -> DeepBSDEConfig:
    horizon = float(benchmark.T) - float(t_start)
    return DeepBSDEConfig(
        benchmark=str(benchmark.name),
        total_time=horizon,
        intervals=spec.time_steps,
        hidden_widths=HIDDEN_WIDTHS,
        learning_rate=LEARNING_RATE,
        iterations=spec.updates if updates is None else int(updates),
        batch_size=spec.batch_size,
        validation_size=VALIDATION_SIZE,
        logging_frequency=RESUME_INTERVAL,
        y_init_range=Y_INIT_RANGE,
        dtype="float64",
    )


def _rng_state() -> dict[str, Any]:
    result: dict[str, Any] = {
        "python": random.getstate(),
        "numpy": np.random.get_state(),
        "torch_cpu": torch.random.get_rng_state(),
    }
    if torch.cuda.is_available():
        result["torch_cuda"] = torch.cuda.get_rng_state_all()
    return result


def _restore_rng(state: dict[str, Any]) -> None:
    random.setstate(state["python"])
    np.random.set_state(state["numpy"])
    torch.random.set_rng_state(state["torch_cpu"])
    if torch.cuda.is_available() and "torch_cuda" in state:
        torch.cuda.set_rng_state_all(state["torch_cuda"])


def _atomic_torch(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, suffix=".pt", delete=False) as stream:
        temporary = Path(stream.name)
    try:
        torch.save(payload, temporary)
        temporary.replace(path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _finite_parameters(model: torch.nn.Module) -> bool:
    return all(bool(torch.isfinite(parameter).all()) for parameter in model.parameters())


def _validate(model: VariableStartDeepBSDE, dw: torch.Tensor) -> tuple[float, int]:
    model.eval()
    with torch.no_grad():
        terminal_y, terminal_x = model.rollout(dw)
        target = model.benchmark.terminal(terminal_x)
        loss = clipped_terminal_loss(terminal_y - target, clip=TERMINAL_LOSS_CLIP)
        bad = int((~torch.isfinite(terminal_y)).sum())
        bad += int((~torch.isfinite(terminal_x)).sum())
        bad += int((~torch.isfinite(target)).sum())
        bad += int(not bool(torch.isfinite(loss)))
    model.train()
    return float(loss), bad


def _paths(pde: str, trajectory_id: str, seed: int, position: int, *, preflight: bool) -> tuple[Path, Path, Path]:
    prefix = "preflight" if preflight else "formal"
    root = RAW_DIR / prefix / pde / trajectory_id / f"s{seed}"
    checkpoint_root = CHECKPOINT_DIR / prefix / pde / trajectory_id / f"s{seed}" / f"q{position:02d}"
    return root / f"q{position:02d}.json", checkpoint_root / "resume.pt", checkpoint_root


def train_query(
    pde: str,
    trajectory: TrajectorySpec,
    seed: int,
    position: int,
    original_index: int,
    t_start: float,
    x_start: np.ndarray,
    device: torch.device,
    *,
    preflight_updates: int | None = None,
) -> dict[str, Any]:
    """Train one native scalar initial value without opening any truth file."""

    ensure_directories()
    preflight = preflight_updates is not None
    output, resume, checkpoint_root = _paths(
        pde, trajectory.trajectory_id, seed, position, preflight=preflight
    )
    if output.is_file() and not preflight:
        import json
        existing = json.loads(output.read_text(encoding="utf-8"))
        if existing.get("status") in {"complete", "numerical_failure"}:
            return existing
    benchmark = make_benchmark(pde, DIMENSION)
    wanted_updates = (
        (int(preflight_updates),)
        if preflight
        else tuple(int(value) for value in trajectory.capture_updates)
    )
    max_updates = int(preflight_updates) if preflight else trajectory.updates
    config = deep_config(benchmark, t_start, trajectory, updates=max_updates)
    query_hash = array_sha256(
        np.asarray(t_start, dtype=np.float64), np.asarray(x_start, dtype=np.float64)
    )
    full_seed = stable_seed(
        REVISION, "preflight" if preflight else "formal", pde,
        int(original_index), int(seed), trajectory.signature,
    )
    seed_everything(full_seed)
    synchronize(device)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    initialization_started = time.perf_counter()
    model = VariableStartDeepBSDE(
        benchmark, config, torch.as_tensor(x_start, dtype=torch.float64), t_start
    ).to(device=device, dtype=torch.float64)
    optimizer = torch.optim.Adam(model.parameters(), lr=config.learning_rate, eps=1e-8)
    train_generator = torch.Generator(device=device).manual_seed(full_seed + 1)
    validation_generator = torch.Generator(device=device).manual_seed(full_seed + 2)
    validation_dw = _sample_dw(
        config, DIMENSION, config.validation_size,
        generator=validation_generator, device=device,
    )
    synchronize(device)
    initialization_seconds = time.perf_counter() - initialization_started
    completed = 0
    elapsed_previous = 0.0
    captures: list[dict[str, Any]] = []
    nonfinite_events: list[dict[str, Any]] = []
    validation_hash = array_sha256(validation_dw.detach().cpu().numpy())
    if resume.is_file() and not preflight:
        state = torch.load(resume, map_location=device, weights_only=False)
        expected_resume = {
            "revision": REVISION,
            "trajectory_signature": trajectory.signature,
            "query_input_sha256": query_hash,
            "full_seed": full_seed,
        }
        if any(state.get(key) != value for key, value in expected_resume.items()):
            raise RuntimeError("resume metadata mismatch")
        model.load_state_dict(state["model"])
        optimizer.load_state_dict(state["optimizer"])
        train_generator.set_state(state["train_generator"])
        validation_dw = state["validation_dw"].to(device=device, dtype=torch.float64)
        _restore_rng(state["rng"])
        completed = int(state["completed_updates"])
        elapsed_previous = float(state["active_training_seconds"])
        captures = list(state["captures"])
        nonfinite_events = list(state["nonfinite_events"])
    capture_set = set(wanted_updates)
    synchronize(device)
    segment_started = time.perf_counter()
    failure_reason: str | None = None
    while completed < max_updates:
        model.train()
        dw = _sample_dw(
            config, DIMENSION, config.batch_size,
            generator=train_generator, device=device,
        )
        terminal_y, terminal_x = model.rollout(dw)
        target = benchmark.terminal(terminal_x)
        loss = clipped_terminal_loss(
            terminal_y - target, clip=TERMINAL_LOSS_CLIP
        )
        tensors_finite = all(bool(torch.isfinite(value).all()) for value in (terminal_y, terminal_x, target, loss))
        if not tensors_finite:
            failure_reason = "nonfinite_internal_trajectory_or_loss"
            nonfinite_events.append({"completed_updates": completed, "reason": failure_reason})
            break
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        try:
            grad_norm = torch.nn.utils.clip_grad_norm_(
                model.parameters(), max_norm=float("inf"), error_if_nonfinite=True
            )
        except RuntimeError:
            failure_reason = "nonfinite_gradients"
            nonfinite_events.append({"completed_updates": completed, "reason": failure_reason})
            break
        optimizer.step()
        completed += 1
        if not math.isfinite(float(model.y_init.detach())) or not math.isfinite(float(grad_norm)):
            failure_reason = "nonfinite_post_update_state"
            nonfinite_events.append({"completed_updates": completed, "reason": failure_reason})
            break
        capture_due = completed in capture_set
        resume_due = completed % RESUME_INTERVAL == 0 or capture_due or completed == max_updates
        if resume_due:
            synchronize(device)
            active_seconds = elapsed_previous + time.perf_counter() - segment_started
            if capture_due:
                validation_loss, bad = _validate(model, validation_dw)
                if bad or not _finite_parameters(model):
                    failure_reason = "nonfinite_validation_or_parameters"
                    nonfinite_events.append({
                        "completed_updates": completed,
                        "reason": failure_reason,
                        "validation_nonfinite_count": bad,
                    })
                checkpoint = checkpoint_root / f"update_{completed:04d}.pt"
                _atomic_torch(checkpoint, {
                    "revision": REVISION,
                    "trajectory_signature": trajectory.signature,
                    "query_input_sha256": query_hash,
                    "full_seed": full_seed,
                    "completed_updates": completed,
                    "state_dict": {key: value.detach().cpu() for key, value in model.state_dict().items()},
                    "configuration": asdict(config),
                    "t_start": float(t_start),
                    "original_index": int(original_index),
                })
                captures.append({
                    "updates": completed,
                    "prediction": float(model.y_init.detach()),
                    "validation_terminal_clipped_mse": validation_loss,
                    "active_training_seconds": active_seconds,
                    "checkpoint": str(checkpoint),
                    "checkpoint_sha256": file_sha256(checkpoint),
                })
            state = {
                "revision": REVISION,
                "trajectory_signature": trajectory.signature,
                "query_input_sha256": query_hash,
                "full_seed": full_seed,
                "model": model.state_dict(),
                "optimizer": optimizer.state_dict(),
                "train_generator": train_generator.get_state(),
                "validation_dw": validation_dw.detach().cpu(),
                "validation_dw_sha256": validation_hash,
                "rng": _rng_state(),
                "completed_updates": completed,
                "captures": captures,
                "nonfinite_events": nonfinite_events,
                "active_training_seconds": active_seconds,
            }
            if not preflight:
                _atomic_torch(resume, state)
            elapsed_previous = active_seconds
            segment_started = time.perf_counter()
        if failure_reason is not None:
            break
    synchronize(device)
    active_seconds = elapsed_previous + time.perf_counter() - segment_started
    expected_captures = set(wanted_updates)
    observed_captures = {int(row["updates"]) for row in captures}
    status = (
        "complete"
        if completed == max_updates and expected_captures == observed_captures and failure_reason is None
        else "numerical_failure"
    )
    record = {
        "revision": REVISION,
        "status": status,
        "preflight": preflight,
        "pde": pde,
        "dimension": DIMENSION,
        "trajectory_id": trajectory.trajectory_id,
        "trajectory_signature": trajectory.signature,
        "seed": int(seed),
        "full_seed": full_seed,
        "position": int(position),
        "original_index": int(original_index),
        "t_start": float(t_start),
        "remaining_horizon": float(benchmark.T - t_start),
        "euler_dt": float((benchmark.T - t_start) / trajectory.time_steps),
        "query_input_sha256": query_hash,
        "configuration": asdict(config),
        "benchmark_metadata": benchmark_metadata(benchmark),
        "captures": captures,
        "completed_updates": completed,
        "expected_updates": max_updates,
        "initialization_seconds": initialization_seconds,
        "active_training_seconds": active_seconds,
        "nonfinite_events": nonfinite_events,
        "failure_reason": failure_reason,
        "validation_dw_sha256": validation_hash,
        "parameter_count": sum(parameter.numel() for parameter in model.parameters()),
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device)) if device.type == "cuda" else None,
        "native_scalar_initial_value": True,
        "amortized_decoder": False,
        "truth_loaded": False,
        "z_convention": "z=sqrt(2)*grad(u)",
        "selection": "fixed checkpoint schedule; no reference used",
        "slurm": {
            "job_id": __import__("os").environ.get("SLURM_JOB_ID"),
            "array_job_id": __import__("os").environ.get("SLURM_ARRAY_JOB_ID"),
            "array_task_id": __import__("os").environ.get("SLURM_ARRAY_TASK_ID"),
            "node": __import__("os").environ.get("SLURMD_NODENAME"),
        },
    }
    write_json_atomic(output, record)
    return record


def load_input(pde: str, split: str, position: int) -> tuple[int, float, np.ndarray]:
    path = INPUT_DIR / f"{pde}_{split}_inputs.npz"
    with np.load(path, allow_pickle=False) as payload:
        indices = np.asarray(payload["original_index"], dtype=np.int64)
        t = np.asarray(payload["t"], dtype=np.float64)
        x = np.asarray(payload["x"], dtype=np.float64)
    return int(indices[position]), float(t[position]), np.asarray(x[position], dtype=np.float64)


__all__ = [
    "VariableStartDeepBSDE", "deep_config", "load_input", "require_h100",
    "train_query",
]

