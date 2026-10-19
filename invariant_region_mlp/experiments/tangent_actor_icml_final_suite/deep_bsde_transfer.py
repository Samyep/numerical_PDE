"""Native single-point Deep BSDE transfers for the ICML suite.

This keeps the audited official PyTorch compatibility port: scalar ``Y0``,
free ``Z0``, 19 nonshared time networks, batch normalization, 20 Euler steps,
float64, clipped terminal loss, and the published per-run optimizer budget.
Only dimension and equation coefficients/terminal data change.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import math
from pathlib import Path
import random
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.drivers import (
    JointConcaveDriver, JointConvexOriginalDriver,
)
from invariant_region_mlp.experiments.nonconvex_tangent_actor.drivers import SineQuadraticDriver
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.common import (
    seed_everything, synchronize,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.deep_bsde_port import (
    DeepBSDEConfig, OfficialDeepBSDE, _sample_dw, _validation, clipped_terminal_loss,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.references import (
    hjb_log_reference,
)

from .benchmarks import AllenCahnBenchmark, HJBLogBenchmark, MLPFailureBenchmark, OUGMMBenchmark
from .prepare_suite import write_json
from .protocol import RESULTS_DIR, SOURCE_COMMITS, stable_seed
from .run_joint_actor import reference as periodic_reference
from .run_standard_actor import _allen_reference


FAMILIES = (
    "hjb_log", "ou_gmm", "joint_concave", "joint_convex",
    "nonconvex_coupled", "mlp_failure", "allen_cahn",
)
INPUT_DIR = RESULTS_DIR / "inputs"
RAW_DIR = RESULTS_DIR / "raw"
CHECKPOINT_DIR = RESULTS_DIR / "official_baseline_checkpoints"


class DriverBenchmark:
    actor_kind = "gradient"
    input_scale = 1.0
    value_bounds = None

    def __init__(self, family: str, d: int):
        self.family, self.d = family, int(d)
        self.name = f"{family}_d{d}"
        if family == "joint_concave":
            self.impl = JointConcaveDriver(d=d, k=2, T=0.5)
        elif family == "joint_convex":
            self.impl = JointConvexOriginalDriver(d=d, k=2, T=0.5)
        elif family == "nonconvex_coupled":
            self.impl = SineQuadraticDriver(d=d, T=0.5, lam=0.2, eta=1.0)
        else:
            raise KeyError(family)
        self.T, self.sigma = self.impl.T, self.impl.sigma

    def terminal(self, x):
        return self.impl.terminal(x)

    def terminal_gradient(self, x):
        return self.impl.terminal_gradient_torch(x)

    def base_drift(self, t, x):
        del t
        return torch.zeros_like(x)

    def driver(self, r, q):
        return self.impl.driver(r, q)


def make_benchmark(family: str, d: int):
    if family == "hjb_log":
        return HJBLogBenchmark(d)
    if family == "ou_gmm":
        return OUGMMBenchmark.from_npz(INPUT_DIR / f"ou_gmm_d{d}.npz")
    if family == "allen_cahn":
        return AllenCahnBenchmark(d)
    if family == "mlp_failure":
        return MLPFailureBenchmark(d)
    return DriverBenchmark(family, d)


class SuiteDeepBSDE(OfficialDeepBSDE):
    def driver(self, y: torch.Tensor, z: torch.Tensor) -> torch.Tensor:
        value = self.benchmark.driver(y[:, 0], z)
        return value.reshape(-1, 1)

    def rollout(self, dw: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        batch = len(dw)
        dt = self.config.total_time / self.config.intervals
        x = torch.zeros((batch, self.benchmark.d), dtype=dw.dtype, device=dw.device)
        y = self.y_init.reshape(1, 1).expand(batch, 1)
        z = self.z_init.expand(batch, self.benchmark.d)
        for index in range(self.config.intervals - 1):
            y = y - dt * self.driver(y, z) + (z * dw[:, index]).sum(1, keepdim=True)
            t = torch.full((batch,), index * dt, dtype=dw.dtype, device=dw.device)
            drift = self.benchmark.base_drift(t, x)
            x = x + drift * dt + self.benchmark.sigma * dw[:, index]
            z = self.subnets[index](x) / self.benchmark.d
        y = y - dt * self.driver(y, z) + (z * dw[:, -1]).sum(1, keepdim=True)
        t = torch.full((batch,), (self.config.intervals - 1) * dt, dtype=dw.dtype, device=dw.device)
        x = x + self.benchmark.base_drift(t, x) * dt + self.benchmark.sigma * dw[:, -1]
        return y[:, 0], x


def _config(family: str, d: int, *, smoke: bool) -> DeepBSDEConfig:
    allen = family == "allen_cahn"
    return DeepBSDEConfig(
        benchmark=f"{family}_d{d}",
        total_time=0.3 if allen else (0.5 if family in {
            "joint_concave", "joint_convex", "nonconvex_coupled"
        } else 1.0),
        intervals=20, hidden_widths=(d + 10, d + 10),
        learning_rate=5e-4 if allen else 1e-2,
        iterations=(4 if smoke else (4000 if allen else 2000)),
        batch_size=(8 if smoke else 64), validation_size=(16 if smoke else 256),
        logging_frequency=(1 if smoke else 100),
        y_init_range=(0.3, 0.5) if allen else (0.0, 1.0), dtype="float64",
    )


def _reference(family: str, benchmark) -> tuple[float, dict[str, Any]]:
    d = benchmark.d
    if family == "hjb_log":
        value = float(hjb_log_reference(np.asarray([0.0]), np.zeros((1, d)))[0])
        return value, {"kind": "analytic Hopf-Cole quadrature"}
    if family == "ou_gmm":
        value = float(benchmark.exact_value(
            torch.zeros(1, dtype=torch.float64), torch.zeros((1, d), dtype=torch.float64)
        )[0])
        return value, {"kind": "exact evolving Gaussian-mixture density"}
    if family == "mlp_failure":
        return math.sqrt(2.0 / math.pi), {"kind": "exact Gaussian absolute moment"}
    if family == "allen_cahn":
        value, audit = _allen_reference(
            d, np.asarray([0.0]), np.zeros((1, d), dtype=np.float64), smoke=False,
        )
        return float(value[0]), audit
    ref, metadata = periodic_reference(family)
    driver = benchmark.impl
    working = float(ref.value(np.asarray([0.0]), np.zeros((1, 2)))[0])
    sign = -1.0 if family == "joint_convex" else 1.0
    return sign * working, {**metadata, "kind": "independent periodic 2D reference"}


def _rng_state() -> dict[str, Any]:
    result = {
        "python": random.getstate(), "numpy": np.random.get_state(),
        "torch": torch.random.get_rng_state(),
    }
    if torch.cuda.is_available():
        result["cuda"] = torch.cuda.get_rng_state_all()
    return result


def _restore_rng(state: dict[str, Any]) -> None:
    random.setstate(state["python"]); np.random.set_state(state["numpy"])
    torch.random.set_rng_state(state["torch"])
    if torch.cuda.is_available() and "cuda" in state:
        torch.cuda.set_rng_state_all(state["cuda"])


def execute(family: str, d: int, seed: int, device: torch.device, *, smoke: bool = False):
    if not smoke and seed != 0:
        raise ValueError("single-seed benchmark protocol permits only seed 0")
    RAW_DIR.mkdir(parents=True, exist_ok=True); CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)
    tag = "smoke" if smoke else "official"
    run_id = f"deep_bsde_{family}_d{d}_{tag}_s{seed}"
    output = RAW_DIR / f"{run_id}.json"
    if output.exists() and json.loads(output.read_text()).get("status") in {"complete", "smoke_nonresult"}:
        return json.loads(output.read_text())
    config = _config(family, d, smoke=smoke)
    benchmark = make_benchmark(family, d)
    full_seed = stable_seed("icml_deep_bsde_transfer", family, d, seed)
    seed_everything(full_seed)
    model = SuiteDeepBSDE(benchmark, config).to(device=device, dtype=torch.float64)
    optimizer = torch.optim.Adam(model.parameters(), lr=config.learning_rate, eps=1e-8)
    train_generator = torch.Generator(device=device).manual_seed(full_seed + 1)
    valid_generator = torch.Generator(device=device).manual_seed(full_seed + 2)
    validation_dw = _sample_dw(
        config, d, config.validation_size, generator=valid_generator, device=device,
    )
    checkpoint = CHECKPOINT_DIR / f"{run_id}.pt"
    resume = CHECKPOINT_DIR / f"{run_id}_resume.pt"
    history: list[dict[str, Any]] = []
    start_step, nonfinite = 0, 0
    if resume.exists():
        state = torch.load(resume, map_location=device, weights_only=False)
        if state["config"] != asdict(config):
            raise RuntimeError("resume configuration mismatch")
        model.load_state_dict(state["model"]); optimizer.load_state_dict(state["optimizer"])
        train_generator.set_state(state["train_generator"])
        _restore_rng(state["rng"])
        start_step = int(state["next_step"]); history = state["history"]
        nonfinite = int(state["nonfinite_count"])
    synchronize(device)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    official_y = float(model.y_init.detach())
    for step in range(start_step, config.iterations + 1):
        if step % config.logging_frequency == 0:
            loss_value, bad = _validation(model, validation_dw)
            nonfinite += bad; official_y = float(model.y_init.detach())
            history.append({
                "step": step, "validation_terminal_clipped_mse": loss_value,
                "y_init": official_y, "elapsed_seconds_this_invocation": time.perf_counter() - started,
                "nonfinite_count_cumulative": nonfinite,
            })
            print(f"{run_id} step={step} loss={loss_value:.6g} y0={official_y:.9g}", flush=True)
        model.train()
        dw = _sample_dw(config, d, config.batch_size, generator=train_generator, device=device)
        terminal_y, terminal_x = model.rollout(dw)
        loss = clipped_terminal_loss(terminal_y - benchmark.terminal(terminal_x))
        optimizer.zero_grad(set_to_none=True); loss.backward(); optimizer.step()
        if not torch.isfinite(loss):
            nonfinite += 1
        if (step + 1) % config.logging_frequency == 0:
            temporary = resume.with_suffix(".tmp")
            torch.save({
                "model": model.state_dict(), "optimizer": optimizer.state_dict(),
                "train_generator": train_generator.get_state(), "rng": _rng_state(),
                "next_step": step + 1, "history": history,
                "nonfinite_count": nonfinite, "config": asdict(config),
            }, temporary); temporary.replace(resume)
    synchronize(device)
    elapsed = time.perf_counter() - started
    torch.save({
        "state_dict": {key: value.detach().cpu() for key, value in model.state_dict().items()},
        "config": asdict(config), "family": family, "dimension": d, "seed": seed,
    }, checkpoint)
    truth, reference = _reference(family, benchmark)
    result = {
        "run_id": run_id, "method": "deep_bsde_official_pytorch_compatibility_port",
        "family": family, "dimension": d, "seed": seed,
        "status": "smoke_nonresult" if smoke else "complete",
        "classification": (
            "native_official" if (d == 100 and family in {"hjb_log", "allen_cahn"})
            else ("adapted" if family == "mlp_failure" else "faithful_transfer")
        ),
        "configuration": asdict(config), "source_configuration": (
            "allencahn_d100.json" if family == "allen_cahn" else "hjb_lq_d100.json"
        ),
        "official_source_commit": SOURCE_COMMITS["deep_bsde"],
        "framework_change": "TensorFlow/Keras to audited PyTorch 2.10 compatibility port",
        "algorithmic_changes": [], "native_scalar_y0": True, "amortized_y0": False,
        "terminal_aware_features": False, "mcv": False,
        "reported_y0": official_y, "post_final_update_y0": float(model.y_init.detach()),
        "native_reference": truth, "native_reference_metadata": reference,
        "native_absolute_error": abs(official_y - truth),
        "native_relative_error": abs(official_y - truth) / max(abs(truth), 1e-300),
        "history": history, "training_seconds_this_invocation": elapsed,
        "parameter_count": sum(parameter.numel() for parameter in model.parameters()),
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device)) if device.type == "cuda" else None,
        "nonfinite_count": nonfinite, "checkpoint": str(checkpoint),
        "selection": "published fixed iteration count; no reference used",
        "full_training_cost_charged_per_initial_state": True,
    }
    write_json(output, result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--family", choices=FAMILIES, required=True)
    parser.add_argument("--dimension", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA allocation required")
    result = execute(args.family, args.dimension, args.seed, device, smoke=args.smoke)
    print(json.dumps({
        "run_id": result["run_id"], "status": result["status"],
        "y0": result["reported_y0"], "absolute_error": result["native_absolute_error"],
    }, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
