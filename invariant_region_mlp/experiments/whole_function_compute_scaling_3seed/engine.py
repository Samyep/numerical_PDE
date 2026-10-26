"""Adapters around the audited TA/DPI/PINN implementations.

Each invocation runs in an isolated result root.  The mature numerical code is
imported unchanged; this module only supplies the frozen run-specific seed and
one-parameter configuration, plus an exact-checkpoint PINN loop.
"""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import asdict
import json
from pathlib import Path
import time
from typing import Any, Iterator

import numpy as np
import torch

from invariant_region_mlp.experiments.whole_function_compute_frontier import (
    baseline_frontier as base,
)
from invariant_region_mlp.experiments.whole_function_compute_frontier import (
    ta_frontier as ta,
)
from invariant_region_mlp.experiments.whole_function_compute_frontier.common import (
    file_sha256,
    synchronize,
    write_json_atomic,
)

from .protocol import (
    CALIBRATION_DIR,
    FORMAL_JOB_SAFETY_SECONDS,
    REVISION,
    RUNS_DIR,
    TA_EVALUATION_PATHS,
    TA_EVALUATION_STEPS,
    TA_INFERENCE_REPLICATES,
    RunSpec,
    ensure_directories,
)


def run_dir(spec: RunSpec) -> Path:
    parent = CALIBRATION_DIR / "runs" if spec.phase == "calibration" else RUNS_DIR
    return parent / spec.run_id


def trajectory_path(spec: RunSpec) -> Path:
    root = run_dir(spec) / "raw" / "trajectories"
    short = "ta" if spec.method == "tangent_actor" else spec.method
    return root / f"{short}_{spec.pde}_d100_c0_s{spec.seed}.json"


def _grid(spec: RunSpec) -> dict[str, dict[str, dict[str, int]]]:
    values = spec.settings
    return {
        "tangent_actor": {
            "C0": {
                "paths_per_label": int(values.get("paths_per_label", 128)),
                "optimizer_steps": int(values.get("optimizer_steps", 5000)),
                "actor_updates": int(values.get("actor_updates", 4)),
                "continuation_actor_updates_max": int(values.get("actor_updates", 4)),
            }
        },
        "dpi": {
            "C0": {
                "epochs_per_picard": int(values.get("epochs_per_picard", 16)),
                "integral_mc_samples": int(values.get("integral_mc_samples", 4096)),
                "original_picard_iterations": int(values.get("picard_iterations", 80)),
                "picard_iterations_max": int(values.get("picard_iterations", 80)),
            }
        },
        "pinn_hte": {
            "C0": {
                "hutchinson_vectors": int(values.get("hutchinson_vectors", 16)),
                "batch_size": int(values.get("batch_size", 512)),
                "original_epochs": int(values.get("training_epochs", 30000)),
            }
        },
    }


@contextmanager
def _patched(spec: RunSpec) -> Iterator[None]:
    """Patch run-local constants in the already-audited adapters."""

    root = run_dir(spec)
    roots = {
        "RAW_DIR": root / "raw",
        "CHECKPOINT_DIR": root / "checkpoints",
        "PREDICTION_DIR": root / "predictions",
        "STATUS_DIR": root / "status",
        "WORK_DIR": root / "work",
    }
    for path in roots.values():
        path.mkdir(parents=True, exist_ok=True)
    grid = _grid(spec)
    changes: list[tuple[Any, str, Any]] = []

    def set_value(module: Any, name: str, value: Any) -> None:
        changes.append((module, name, getattr(module, name)))
        setattr(module, name, value)

    for module in (ta, base):
        set_value(module, "REVISION", REVISION)
        set_value(module, "TRAINING_SEED", int(spec.seed))
        set_value(module, "CONFIGURATION_GRID", grid)
        for name, value in roots.items():
            if hasattr(module, name):
                set_value(module, name, value)
    set_value(base, "ACTIVE_TRAINING_CAP_SECONDS", FORMAL_JOB_SAFETY_SECONDS)
    set_value(base, "TARGET_MINUTES", tuple())
    set_value(ta, "TA_EVALUATION_PATHS", TA_EVALUATION_PATHS)
    set_value(ta, "TA_EVALUATION_STEPS", TA_EVALUATION_STEPS)
    set_value(ta, "TA_INFERENCE_REPLICATES", TA_INFERENCE_REPLICATES)

    # C0 normally reuses a historical seed-0 trajectory.  This study requires
    # fresh runs for all three seeds and configurations, so present an empty,
    # configuration-compatible history to the unchanged trainer.
    config = ta._configs(spec.pde, "C0") if spec.method == "tangent_actor" else None
    set_value(
        ta,
        "_historical_training",
        lambda _pde: {"status": "complete", "configuration": asdict(config), "phases": []},
    )
    try:
        yield
    finally:
        for module, name, old in reversed(changes):
            setattr(module, name, old)


def _annotate_record(spec: RunSpec, record: dict[str, Any]) -> dict[str, Any]:
    record.update({
        "study_revision": REVISION,
        "study_spec": spec.payload(),
        "run_id": spec.run_id,
        "truth_loaded": False,
        "terminal_checkpoint_only_for_sensitivity": spec.sweep != "C0",
    })
    write_json_atomic(trajectory_path(spec), record)
    return record


def _train_pinn(spec: RunSpec, device: torch.device) -> dict[str, Any]:
    """Train PINN to exact preregistered epoch boundaries.

    Microbatch accumulation forms the same mean Monte Carlo objective and
    performs one optimizer step per nominal full batch.  The upstream network
    has no batch-dependent layers.
    """

    output = trajectory_path(spec)
    if output.is_file():
        existing = json.loads(output.read_text(encoding="utf-8"))
        if existing.get("status") in {"TRAINED_NO_TRUTH", "COMPLETE"}:
            return existing

    base._seed_all(int(spec.seed))
    cfg, _source, _overrides, _profile = base._configuration(
        "pinn_hte", spec.pde, "C0"
    )
    directory = run_dir(spec) / "checkpoints" / "pinn_hte" / spec.pde / "C0"
    directory.mkdir(parents=True, exist_ok=True)
    model, equation, optimizer = base._new_pinn(
        cfg, None, directory / "unused_capture.pt", device
    )
    from picard.utils import hutchinson_trace_estimation_batch

    final_epoch = int(spec.settings["training_epochs"])
    wanted = (
        set(int(v) for v in spec.checkpoint_values)
        if spec.sweep == "C0" else {final_epoch}
    )
    wanted.add(final_epoch)
    batch_size = int(spec.settings["batch_size"])
    hte_vectors = int(spec.settings["hutchinson_vectors"])
    microbatch_size = min(batch_size, 512)
    checkpoints: list[dict[str, Any]] = []
    cumulative = 0.0
    last_loss: float | None = None
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    synchronize(device)
    segment_started = time.perf_counter()

    for epoch in range(1, final_epoch + 1):
        optimizer.zero_grad(set_to_none=True)
        weighted_loss_value = 0.0
        processed = 0
        while processed < batch_size:
            count = min(microbatch_size, batch_size - processed)
            t = equation.T * torch.rand(count, 1, device=device)
            x = equation.sample_x(t)
            t.requires_grad_()
            x.requires_grad_()
            output_value = model(torch.cat((t, x), dim=-1))[:, :1]
            u_t = torch.autograd.grad(
                output_value.sum(), t, create_graph=True, retain_graph=True
            )[0]
            u_x = torch.autograd.grad(
                output_value.sum(), x, create_graph=True, retain_graph=True
            )[0]
            laplacian = hutchinson_trace_estimation_batch(
                model, t, x, hte_vectors
            )
            loss_in = equation.pinn_function(
                t, x, output_value, u_t, u_x, laplacian
            ).square().mean()
            terminal_t = equation.T * torch.ones_like(t)
            terminal_x = equation.sample_x(terminal_t)
            terminal_value = model(torch.cat((terminal_t, terminal_x), dim=-1))[:, :1]
            loss_terminal = (
                terminal_value - equation.g(terminal_x)
            ).square().mean()
            loss = loss_in + 10.0 * loss_terminal
            if not torch.isfinite(loss):
                raise FloatingPointError(
                    f"nonfinite PINN loss at epoch={epoch}, microbatch={processed}"
                )
            weight = count / batch_size
            (weight * loss).backward()
            weighted_loss_value += weight * float(loss.detach().cpu())
            processed += count
        optimizer.step()
        last_loss = weighted_loss_value

        if epoch not in wanted:
            continue
        synchronize(device)
        cumulative += time.perf_counter() - segment_started
        checkpoint = directory / f"epoch_{epoch:07d}.pt"
        torch.save({
            "revision": REVISION,
            "state_dict": model.state_dict(),
            "optimizer_state_dict": optimizer.state_dict(),
            "epoch": epoch,
            "training_seconds": cumulative,
            "study_spec": spec.payload(),
        }, checkpoint)
        row = {
            "checkpoint_id": f"epoch_{epoch:07d}",
            "algorithm_iteration": epoch,
            "checkpoint": str(checkpoint),
            "checkpoint_sha256": file_sha256(checkpoint),
            "training_seconds": cumulative,
            "selection_seconds": 0.0,
            "main_axis_seconds": cumulative,
            "validation_metric": None,
            "training_loss": last_loss,
            "publish": True,
            "status": "COMPLETE",
            "microbatch_size": microbatch_size,
            "nominal_batch_size": batch_size,
        }
        checkpoints.append(row)
        interim = base._write_training_record(
            "pinn_hte", spec.pde, "C0", checkpoints,
            status="TRAINING_IN_PROGRESS",
        )
        _annotate_record(spec, interim)
        segment_started = time.perf_counter()

    synchronize(device)
    if checkpoints:
        checkpoints[-1]["is_final"] = True
    peak = (
        int(torch.cuda.max_memory_allocated(device))
        if device.type == "cuda" else None
    )
    result = base._write_training_record(
        "pinn_hte", spec.pde, "C0", checkpoints,
        status="TRAINED_NO_TRUTH", stop_reason="configured_training_epochs",
        peak_gpu_bytes=peak,
    )
    result["microbatch_size"] = microbatch_size
    result["nominal_batch_size"] = batch_size
    return _annotate_record(spec, result)


def train_spec(spec: RunSpec, device: torch.device) -> dict[str, Any]:
    ensure_directories()
    root = run_dir(spec)
    root.mkdir(parents=True, exist_ok=True)
    write_json_atomic(root / "spec.json", spec.payload())
    failure = root / "failure.json"
    try:
        with _patched(spec):
            if spec.method == "tangent_actor":
                result = ta.train(spec.pde, "C0", device)
            elif spec.method == "dpi":
                result = base.train_dpi(spec.pde, "C0", device)
            elif spec.method == "pinn_hte":
                result = _train_pinn(spec, device)
            else:
                raise KeyError(spec.method)
            result = _annotate_record(spec, result)
        if failure.exists():
            failure.unlink()
        return result
    except Exception as error:
        write_json_atomic(failure, {
            "run_id": spec.run_id,
            "spec": spec.payload(),
            "error_type": type(error).__name__,
            "error": repr(error),
            "status": "FAILED",
        })
        raise


def _selected_checkpoints(spec: RunSpec, record: dict[str, Any]) -> list[dict[str, Any]]:
    checkpoints = list(record.get("checkpoints", []))
    if spec.sweep == "C0":
        wanted = set(int(value) for value in spec.checkpoint_values)
        selected = [
            row for row in checkpoints
            if int(row["algorithm_iteration"]) in wanted
        ]
    else:
        selected = checkpoints[-1:] if checkpoints else []
    if len(selected) != len(spec.checkpoint_values):
        raise RuntimeError(
            f"{spec.run_id}: expected {len(spec.checkpoint_values)} selected "
            f"checkpoints, found {len(selected)}"
        )
    return selected


def evaluate_spec(spec: RunSpec, device: torch.device) -> list[dict[str, Any]]:
    if spec.phase != "formal":
        raise ValueError("calibration runs must never be evaluated on test truth")
    record_file = trajectory_path(spec)
    record = json.loads(record_file.read_text(encoding="utf-8"))
    if record.get("status") not in {"TRAINED_NO_TRUTH", "COMPLETE"}:
        raise RuntimeError(f"training is incomplete for {spec.run_id}")
    full_checkpoints = list(record["checkpoints"])
    selected = _selected_checkpoints(spec, record)
    filtered = dict(record)
    filtered["checkpoints"] = selected
    filtered["status"] = "TRAINED_NO_TRUTH"

    with _patched(spec):
        gate = run_dir(spec) / "status" / "all_training_complete.json"
        write_json_atomic(gate, {
            "status": "GLOBAL_FORMAL_TRAINING_COMPLETE_TRUTH_GATE_OPEN",
            "run_id": spec.run_id,
        })
        write_json_atomic(record_file, filtered)
        try:
            if spec.method == "tangent_actor":
                rows = ta.evaluate(spec.pde, "C0", device)
            else:
                rows = base.evaluate(spec.method, spec.pde, "C0", device)
        finally:
            restored = dict(record)
            restored["checkpoints"] = full_checkpoints
            restored["status"] = "COMPLETE"
            restored["truth_opened_only_after_global_training_stop"] = True
            write_json_atomic(record_file, restored)

    for row in rows:
        if spec.sweep == "C0":
            point_value = int(row["algorithm_iteration"])
        else:
            point_value = int(spec.parameter_value)
        row.update({
            "study_revision": REVISION,
            "run_id": spec.run_id,
            "sweep": spec.sweep,
            "parameter_name": spec.parameter_name,
            "parameter_value": point_value,
            "training_seed": int(spec.seed),
            "training_seconds": float(row["training_seconds"]),
            "main_axis_seconds": float(row["training_seconds"]),
        })
        metric_path = Path(row["prediction_file"]).with_suffix(".json")
        write_json_atomic(metric_path, row)
    return rows


__all__ = [
    "evaluate_spec", "run_dir", "train_spec", "trajectory_path",
]
