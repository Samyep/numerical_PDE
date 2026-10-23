"""Audited single-initial-state Deep BSDE with exact checkpoint/resume."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import random
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.original_pde_benchmark_comparison.common import (
    seed_everything,
    synchronize,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.deep_bsde_port import (
    DeepBSDEConfig,
    OfficialDeepBSDE,
    _sample_dw,
    _validation,
    clipped_terminal_loss,
)

from .benchmarks import make_benchmark
from .protocol import (
    CHECKPOINT_DIR,
    PREDICTION_DIR,
    QUERY_COUNT,
    RAW_DIR,
    REVISION,
    TRAINING_SEEDS,
    array_sha256,
    deep_bsde_configuration,
    ensure_directories,
    file_sha256,
    query_path,
    require_stage_a,
    require_stage_b,
    stable_seed,
    write_json_atomic,
    write_npz_atomic,
)


class FixedInitialDeepBSDE(OfficialDeepBSDE):
    """Official architecture with an arbitrary, immutable fixed initial state."""

    def __init__(self, benchmark, config: DeepBSDEConfig, x0: torch.Tensor) -> None:
        super().__init__(benchmark, config)
        fixed = torch.as_tensor(x0, dtype=torch.float64).reshape(1, benchmark.d)
        self.register_buffer("fixed_x0", fixed)

    def driver(self, y: torch.Tensor, z: torch.Tensor) -> torch.Tensor:
        return self.benchmark.driver(y[:, 0], z).reshape(-1, 1)

    def rollout(self, dw: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        batch = len(dw)
        dt = self.config.total_time / self.config.intervals
        x = self.fixed_x0.to(dw).expand(batch, self.benchmark.d).clone()
        y = self.y_init.reshape(1, 1).expand(batch, 1)
        z = self.z_init.expand(batch, self.benchmark.d)
        for index in range(self.config.intervals - 1):
            y = y - dt * self.driver(y, z) + (z * dw[:, index]).sum(1, keepdim=True)
            time_value = torch.full(
                (batch,), index * dt, dtype=dw.dtype, device=dw.device
            )
            x = x + self.benchmark.base_drift(time_value, x) * dt + self.benchmark.sigma * dw[:, index]
            z = self.subnets[index](x) / self.benchmark.d
        y = y - dt * self.driver(y, z) + (z * dw[:, -1]).sum(1, keepdim=True)
        time_value = torch.full(
            (batch,), (self.config.intervals - 1) * dt,
            dtype=dw.dtype, device=dw.device,
        )
        x = x + self.benchmark.base_drift(time_value, x) * dt + self.benchmark.sigma * dw[:, -1]
        return y[:, 0], x


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
    temporary = path.with_suffix(path.suffix + ".tmp")
    torch.save(payload, temporary)
    temporary.replace(path)


def _config(pde: str, d: int, *, pilot: bool) -> DeepBSDEConfig:
    payload = deep_bsde_configuration(pde, d, pilot=pilot)
    payload["hidden_widths"] = tuple(payload["hidden_widths"])
    payload["y_init_range"] = tuple(payload["y_init_range"])
    return DeepBSDEConfig(**payload)


def _paths(pde: str, d: int, point: int, seed: int, *, pilot: bool) -> tuple[str, Path, Path, Path]:
    tag = "pilot" if pilot else "formal"
    run_id = f"deep_bsde__{pde}__d{d}__q{point:02d}__s{seed}__{tag}"
    output = RAW_DIR / tag / "deep_bsde" / f"{run_id}.json"
    resume = CHECKPOINT_DIR / tag / "deep_bsde" / f"{run_id}_resume.pt"
    final = CHECKPOINT_DIR / tag / "deep_bsde" / f"{run_id}.pt"
    return run_id, output, resume, final


def execute(
    pde: str,
    d: int,
    point: int,
    seed: int,
    device: torch.device,
    *,
    pilot: bool = False,
    stop_after_updates: int | None = None,
) -> dict[str, Any]:
    if not 0 <= point < QUERY_COUNT:
        raise ValueError(point)
    if not pilot and seed not in TRAINING_SEEDS:
        raise ValueError(seed)
    require_stage_a()
    if not pilot:
        require_stage_b()
    ensure_directories()
    run_id, output, resume, final = _paths(pde, d, point, seed, pilot=pilot)
    if output.is_file():
        existing = json.loads(output.read_text(encoding="utf-8"))
        if existing.get("status") in {"complete", "pilot_complete"}:
            return existing
    x_path = query_path(pde, d) / "x.npy"
    x_bank = np.load(x_path, allow_pickle=False)
    x0 = np.asarray(x_bank[point], dtype=np.float64)
    config = _config(pde, d, pilot=pilot)
    benchmark = make_benchmark(pde, d)
    full_seed = stable_seed("deep_bsde", pde, d, point, seed, "pilot" if pilot else "formal")
    seed_everything(full_seed)
    model = FixedInitialDeepBSDE(benchmark, config, torch.as_tensor(x0)).to(
        device=device, dtype=torch.float64
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=config.learning_rate, eps=1e-8)
    train_generator = torch.Generator(device=device).manual_seed(full_seed + 1)
    valid_generator = torch.Generator(device=device).manual_seed(full_seed + 2)
    validation_dw = _sample_dw(
        config, d, config.validation_size, generator=valid_generator, device=device
    )
    history: list[dict[str, Any]] = []
    completed_updates = 0
    nonfinite = 0
    elapsed_previous = 0.0
    if resume.is_file():
        state = torch.load(resume, map_location=device, weights_only=False)
        if state["revision"] != REVISION or state["config"] != asdict(config):
            raise RuntimeError("Deep BSDE resume metadata mismatch")
        if state["x0_sha256"] != array_sha256(x0):
            raise RuntimeError("Deep BSDE resume initial-state mismatch")
        model.load_state_dict(state["model"])
        optimizer.load_state_dict(state["optimizer"])
        train_generator.set_state(state["train_generator"])
        validation_dw = state["validation_dw"].to(device=device, dtype=torch.float64)
        _restore_rng(state["rng"])
        completed_updates = int(state["completed_updates"])
        history = list(state["history"])
        nonfinite = int(state["nonfinite_count"])
        elapsed_previous = float(state.get("elapsed_seconds", 0.0))
    validation_hash = array_sha256(validation_dw.detach().cpu().numpy())
    synchronize(device)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    invocation_started = time.perf_counter()
    target_updates = config.iterations
    if stop_after_updates is not None:
        target_updates = min(target_updates, completed_updates + int(stop_after_updates))
    while completed_updates < target_updates:
        model.train()
        dw = _sample_dw(
            config, d, config.batch_size, generator=train_generator, device=device
        )
        terminal_y, terminal_x = model.rollout(dw)
        loss = clipped_terminal_loss(terminal_y - benchmark.terminal(terminal_x))
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        optimizer.step()
        completed_updates += 1
        if not torch.isfinite(loss):
            nonfinite += 1
        checkpoint_due = (
            completed_updates % config.logging_frequency == 0
            or completed_updates == config.iterations
            or completed_updates == target_updates
        )
        if checkpoint_due:
            validation_loss, bad = _validation(model, validation_dw)
            nonfinite += bad
            elapsed = elapsed_previous + time.perf_counter() - invocation_started
            history.append({
                "completed_updates": completed_updates,
                "validation_terminal_clipped_mse": validation_loss,
                "y0": float(model.y_init.detach()),
                "nonfinite_count_cumulative": nonfinite,
                "elapsed_seconds": elapsed,
            })
            _atomic_torch(resume, {
                "revision": REVISION,
                "run_id": run_id,
                "config": asdict(config),
                "x0_sha256": array_sha256(x0),
                "model": model.state_dict(),
                "optimizer": optimizer.state_dict(),
                "train_generator": train_generator.get_state(),
                "validation_dw": validation_dw.detach().cpu(),
                "validation_dw_sha256": validation_hash,
                "rng": _rng_state(),
                "completed_updates": completed_updates,
                "history": history,
                "nonfinite_count": nonfinite,
                "elapsed_seconds": elapsed,
            })
            print(
                f"{run_id} updates={completed_updates}/{config.iterations} "
                f"validation={validation_loss:.6g} y0={float(model.y_init.detach()):.9g}",
                flush=True,
            )
    if completed_updates < config.iterations:
        return {
            "run_id": run_id,
            "status": "checkpointed",
            "completed_updates": completed_updates,
            "target_updates": config.iterations,
            "resume_checkpoint": str(resume),
        }
    synchronize(device)
    total_elapsed = elapsed_previous + time.perf_counter() - invocation_started
    estimate = float(model.y_init.detach())
    _atomic_torch(final, {
        "revision": REVISION,
        "run_id": run_id,
        "state_dict": {key: value.detach().cpu() for key, value in model.state_dict().items()},
        "config": asdict(config),
        "pde": pde,
        "dimension": d,
        "point": point,
        "seed": seed,
        "x0": x0,
    })
    # Ground truth is loaded only after the fixed schedule has completed.
    reference = float(np.load(query_path(pde, d) / "reference.npy", allow_pickle=False)[point])
    prediction_path = PREDICTION_DIR / ("pilot" if pilot else "formal") / "deep_bsde" / f"{run_id}.npz"
    write_npz_atomic(
        prediction_path,
        point=np.asarray(point, dtype=np.int64),
        x0=x0,
        prediction=np.asarray(estimate, dtype=np.float64),
        reference=np.asarray(reference, dtype=np.float64),
    )
    classification = (
        "native_source_family" if pde in {"hjb_log", "allen_cahn"}
        else "faithful_transfer"
    )
    record = {
        "revision": REVISION,
        "run_id": run_id,
        "method": "deep_bsde_official_pytorch_compatibility_port",
        "status": "pilot_complete" if pilot else "complete",
        "classification": classification,
        "pde": pde,
        "dimension": d,
        "point": point,
        "seed": seed,
        "x0_sha256": array_sha256(x0),
        "configuration": asdict(config),
        "reported_y0": estimate,
        "reference": reference,
        "error": estimate - reference,
        "absolute_error": abs(estimate - reference),
        "history": history,
        "completed_updates": completed_updates,
        "nonfinite_count": nonfinite,
        "training_seconds": total_elapsed,
        "parameter_count": sum(parameter.numel() for parameter in model.parameters()),
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device)) if device.type == "cuda" else None,
        "checkpoint": str(final),
        "checkpoint_sha256": file_sha256(final),
        "resume_checkpoint": str(resume),
        "validation_dw_sha256": validation_hash,
        "prediction_file": str(prediction_path),
        "selection": "fixed final scheduled iterate; no truth used",
        "native_scalar_y0": True,
        "amortized_y0": False,
        "terminal_aware_features": False,
        "mcv": False,
        "full_training_cost_charged_per_initial_state": True,
    }
    write_json_atomic(output, record)
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pde", required=True)
    parser.add_argument("--dimension", type=int, required=True)
    parser.add_argument("--point", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument("--stop-after-updates", type=int)
    args = parser.parse_args()
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA device requested but unavailable")
    result = execute(
        args.pde, args.dimension, args.point, args.seed, device,
        pilot=args.pilot, stop_after_updates=args.stop_after_updates,
    )
    print(json.dumps({key: result.get(key) for key in (
        "run_id", "status", "completed_updates", "reported_y0", "absolute_error"
    )}, sort_keys=True))


if __name__ == "__main__":
    main()

