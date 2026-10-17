"""Audited PyTorch compatibility port of the official DeepBSDE solver.

Only framework plumbing and deterministic resource accounting differ.  The
native scalar initial value, free initial z, nonshared time-slice networks,
batch normalization, Euler scheme, clipped terminal loss, precision, and
published optimization schedules are retained.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
import time
from typing import Any

import numpy as np
import torch
from torch import nn

from .benchmarks import AllenCahnBenchmark, HJBLogBenchmark, OriginalBenchmark
from .common import (
    CHECKPOINT_DIR,
    RAW_DIR,
    RESULTS_DIR,
    atomic_json,
    atomic_torch_save,
    ensure_dirs,
    seed_everything,
    stable_seed,
    synchronize,
)
from .references import hjb_log_reference


@dataclass(frozen=True)
class DeepBSDEConfig:
    benchmark: str
    total_time: float
    intervals: int
    hidden_widths: tuple[int, int]
    learning_rate: float
    iterations: int
    batch_size: int
    validation_size: int
    logging_frequency: int
    y_init_range: tuple[float, float]
    dtype: str = "float64"


CONFIGS = {
    "hjb_log_d100": DeepBSDEConfig(
        benchmark="hjb_log_d100", total_time=1.0, intervals=20,
        hidden_widths=(110, 110), learning_rate=1e-2, iterations=2000,
        batch_size=64, validation_size=256, logging_frequency=100,
        y_init_range=(0.0, 1.0),
    ),
    "allen_cahn_d100": DeepBSDEConfig(
        benchmark="allen_cahn_d100", total_time=0.3, intervals=20,
        hidden_widths=(110, 110), learning_rate=5e-4, iterations=4000,
        batch_size=64, validation_size=256, logging_frequency=100,
        y_init_range=(0.3, 0.5),
    ),
}


class OfficialZSubNet(nn.Module):
    def __init__(self, d: int, widths: tuple[int, int]) -> None:
        super().__init__()
        self.input_bn = nn.BatchNorm1d(d, eps=1e-6, momentum=0.01)
        self.hidden1 = nn.Linear(d, widths[0], bias=False)
        self.bn1 = nn.BatchNorm1d(widths[0], eps=1e-6, momentum=0.01)
        self.hidden2 = nn.Linear(widths[0], widths[1], bias=False)
        self.bn2 = nn.BatchNorm1d(widths[1], eps=1e-6, momentum=0.01)
        self.output = nn.Linear(widths[1], d, bias=True)
        self.output_bn = nn.BatchNorm1d(d, eps=1e-6, momentum=0.01)
        for layer in (self.hidden1, self.hidden2, self.output):
            nn.init.xavier_uniform_(layer.weight)
            if layer.bias is not None:
                nn.init.zeros_(layer.bias)
        for layer in (self.input_bn, self.bn1, self.bn2, self.output_bn):
            nn.init.normal_(layer.bias, mean=0.0, std=0.1)
            nn.init.uniform_(layer.weight, a=0.1, b=0.5)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        h = self.input_bn(x)
        h = torch.relu(self.bn1(self.hidden1(h)))
        h = torch.relu(self.bn2(self.hidden2(h)))
        return self.output_bn(self.output(h))


class OfficialDeepBSDE(nn.Module):
    def __init__(
        self, benchmark: OriginalBenchmark, config: DeepBSDEConfig
    ) -> None:
        super().__init__()
        self.benchmark = benchmark
        self.config = config
        self.y_init = nn.Parameter(torch.empty(1))
        nn.init.uniform_(self.y_init, *config.y_init_range)
        self.z_init = nn.Parameter(torch.empty(1, benchmark.d))
        nn.init.uniform_(self.z_init, -0.1, 0.1)
        self.subnets = nn.ModuleList(
            [
                OfficialZSubNet(benchmark.d, config.hidden_widths)
                for _ in range(config.intervals - 1)
            ]
        )

    def driver(self, y: torch.Tensor, z: torch.Tensor) -> torch.Tensor:
        if self.benchmark.name == "hjb_log_d100":
            return -0.5 * z.square().sum(dim=1, keepdim=True)
        if self.benchmark.name == "allen_cahn_d100":
            return y - y.pow(3)
        raise KeyError(self.benchmark.name)

    def rollout(self, dw: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        batch = len(dw)
        dt = self.config.total_time / self.config.intervals
        x = torch.zeros(
            (batch, self.benchmark.d), dtype=dw.dtype, device=dw.device
        )
        y = self.y_init.reshape(1, 1).expand(batch, 1)
        z = self.z_init.expand(batch, self.benchmark.d)
        for index in range(self.config.intervals - 1):
            y = y - dt * self.driver(y, z) + (z * dw[:, index]).sum(
                dim=1, keepdim=True
            )
            x = x + self.benchmark.sigma * dw[:, index]
            z = self.subnets[index](x) / self.benchmark.d
        y = y - dt * self.driver(y, z) + (z * dw[:, -1]).sum(
            dim=1, keepdim=True
        )
        x = x + self.benchmark.sigma * dw[:, -1]
        return y[:, 0], x


def _sample_dw(
    config: DeepBSDEConfig,
    d: int,
    count: int,
    *,
    generator: torch.Generator,
    device: torch.device,
) -> torch.Tensor:
    return torch.randn(
        (count, config.intervals, d),
        generator=generator,
        dtype=torch.float64,
        device=device,
    ) * math.sqrt(config.total_time / config.intervals)


def clipped_terminal_loss(delta: torch.Tensor, clip: float = 50.0) -> torch.Tensor:
    absolute = delta.abs()
    return torch.where(
        absolute < clip, delta.square(), 2.0 * clip * absolute - clip**2
    ).mean()


@torch.no_grad()
def _validation(
    model: OfficialDeepBSDE, dw: torch.Tensor
) -> tuple[float, int]:
    model.eval()
    terminal_y, terminal_x = model.rollout(dw)
    delta = terminal_y - model.benchmark.terminal(terminal_x)
    loss = clipped_terminal_loss(delta)
    nonfinite = int((~torch.isfinite(terminal_y)).sum().item())
    model.train()
    return float(loss.item()), nonfinite


def train_deep_bsde(
    benchmark_name: str, seed: int, device: torch.device
) -> dict[str, Any]:
    invocation_started = time.perf_counter()
    ensure_dirs()
    if benchmark_name not in CONFIGS:
        raise KeyError(benchmark_name)
    config = CONFIGS[benchmark_name]
    benchmark: OriginalBenchmark = (
        HJBLogBenchmark()
        if benchmark_name == "hjb_log_d100"
        else AllenCahnBenchmark()
    )
    full_seed = stable_seed("official_deep_bsde_port", benchmark_name, seed)
    seed_everything(full_seed)
    model = OfficialDeepBSDE(benchmark, config).to(device=device, dtype=torch.float64)
    optimizer = torch.optim.Adam(
        model.parameters(), lr=config.learning_rate, eps=1e-8
    )
    train_generator = torch.Generator(device=device).manual_seed(full_seed + 1)
    validation_generator = torch.Generator(device=device).manual_seed(full_seed + 2)
    validation_dw = _sample_dw(
        config,
        benchmark.d,
        config.validation_size,
        generator=validation_generator,
        device=device,
    )
    initialization_seconds = time.perf_counter() - invocation_started
    history: list[dict[str, Any]] = []
    nonfinite = 0
    synchronize(device)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    official_reported_y = float(model.y_init.detach().item())
    for step in range(config.iterations + 1):
        if step % config.logging_frequency == 0:
            validation_loss, validation_nonfinite = _validation(model, validation_dw)
            nonfinite += validation_nonfinite
            synchronize(device)
            official_reported_y = float(model.y_init.detach().item())
            history.append(
                {
                    "step": step,
                    "validation_terminal_clipped_mse": validation_loss,
                    "y_init": official_reported_y,
                    "elapsed_seconds": time.perf_counter() - started,
                    "nonfinite_count_cumulative": nonfinite,
                }
            )
            print(
                f"deep_bsde {benchmark_name} seed={seed} step={step} "
                f"loss={validation_loss:.6g} y0={official_reported_y:.9g}",
                flush=True,
            )
        model.train()
        dw = _sample_dw(
            config,
            benchmark.d,
            config.batch_size,
            generator=train_generator,
            device=device,
        )
        terminal_y, terminal_x = model.rollout(dw)
        loss = clipped_terminal_loss(
            terminal_y - benchmark.terminal(terminal_x)
        )
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        optimizer.step()
        if not torch.isfinite(loss):
            nonfinite += 1
    synchronize(device)
    elapsed = time.perf_counter() - started
    checkpoint = CHECKPOINT_DIR / f"deep_bsde_{benchmark_name}_s{seed}.pt"
    atomic_torch_save(
        checkpoint,
        {
            "state_dict": {
                key: value.detach().cpu().clone()
                for key, value in model.state_dict().items()
            },
            "benchmark": benchmark_name,
            "seed": seed,
            "config": config.__dict__,
            "official_reported_y": official_reported_y,
        },
    )
    if benchmark_name == "hjb_log_d100":
        native_reference = float(
            hjb_log_reference(
                np.asarray([0.0]), np.zeros((1, benchmark.d), dtype=np.float64)
            )[0]
        )
        reference_status = "analytic Hopf-Cole quadrature"
    else:
        reference_audit = json.loads(
            (RESULTS_DIR / "reference_accuracy_audit.json").read_text()
        )
        native_reference = float(
            reference_audit["allen_cahn"]["native_u_0_0_extrapolated"]
        )
        reference_status = "independent radial numerical reference"
    post_final_z = model.z_init.detach().double().cpu().numpy()[0]
    result = {
        "method": "deep_bsde_official_pytorch_compatibility_port",
        "benchmark": benchmark_name,
        "seed": seed,
        "configuration": config.__dict__,
        "official_source_commit": "e76fed80b1995daf4a0dd6d3e2cf64a0931dabda",
        "framework_change": "TensorFlow/Keras to PyTorch 2.10",
        "algorithmic_changes": [],
        "native_scalar_y0": True,
        "amortized_y0": False,
        "terminal_aware_features": False,
        "mcv": False,
        "reported_y0": official_reported_y,
        "post_final_update_y0": float(model.y_init.detach().item()),
        "post_final_update_z0_component_rmse": float(
            np.sqrt(np.mean(np.square(post_final_z)))
        ),
        "post_final_update_z0_vector_rmse": float(np.linalg.norm(post_final_z)),
        "native_exact_z0": "zero by radial symmetry",
        "native_reference": native_reference,
        "native_reference_status": reference_status,
        "native_absolute_error": abs(official_reported_y - native_reference),
        "native_relative_error": abs(official_reported_y - native_reference)
        / abs(native_reference),
        "history": history,
        "model_data_optimizer_initialization_seconds": initialization_seconds,
        "training_seconds": elapsed,
        "total_process_seconds": time.perf_counter() - invocation_started,
        "parameter_count": sum(parameter.numel() for parameter in model.parameters()),
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device))
        if device.type == "cuda"
        else None,
        "nonfinite_count": nonfinite,
        "checkpoint": str(checkpoint),
        "selection": "published fixed iteration count; no reference used",
    }
    atomic_json(RAW_DIR / f"deep_bsde_{benchmark_name}_s{seed}.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark", choices=tuple(CONFIGS), required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("a CUDA allocation is required")
    result = train_deep_bsde(args.benchmark, args.seed, device)
    print(json.dumps({"status": "complete", "y0": result["reported_y0"]}))


if __name__ == "__main__":
    main()

