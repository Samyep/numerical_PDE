"""Anytime, resumable DPI and PINN-HTE trajectories.

The numerical objectives and networks come from the pinned DeepPicardIteration
checkout.  This module adds only boundary-safe continuation, synchronized time
accounting, truth-free validation, and post-stop evaluation on the frozen bank.
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import time
from types import SimpleNamespace
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.main_highdim_expansion import baseline_transfer as transfer
from invariant_region_mlp.experiments.main_highdim_expansion.benchmarks import array_sha256

from .common import (
    file_sha256,
    plateau_reached,
    restore_rng,
    rng_state,
    save_npz_atomic,
    synchronize,
    value_metrics,
    write_json_atomic,
)
from .prepare import historical_results_root
from .protocol import (
    ACTIVE_TRAINING_CAP_SECONDS,
    CHECKPOINT_DIR,
    CONFIGURATION_GRID,
    DIMENSION,
    PINN_CHECK_INTERVAL_EPOCHS,
    PINN_VALIDATION_HTE_VECTORS,
    PINN_VALIDATION_SEED,
    PINN_VALIDATION_STATES,
    PREDICTION_DIR,
    RAW_DIR,
    REVISION,
    STATUS_DIR,
    TARGET_MINUTES,
    TRAINING_SEED,
    UPSTREAM_DPI_REVISION,
    WORK_DIR,
    ensure_directories,
)


METHODS = ("dpi", "pinn_hte")


def _external_root() -> Path:
    configured = os.environ.get("WCF_DPI_ROOT")
    root = Path(configured).resolve() if configured else transfer.EXTERNAL_ROOT.resolve()
    if not (root / "picard" / "picard_iteration.py").is_file():
        raise FileNotFoundError("set WCF_DPI_ROOT to the pinned DeepPicardIteration checkout")
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=root, text=True
    ).strip()
    if revision != UPSTREAM_DPI_REVISION:
        raise RuntimeError(f"pinned DPI source mismatch: {revision}")
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    # baseline_transfer resolves these constants at import time.  Rebind them
    # for isolated worktrees whose sibling directory is intentionally empty.
    transfer.EXTERNAL_ROOT = root
    transfer.SEMILINEAR_DPI_CONFIG = root / "scripts/burgers/base_100d_T1.0_w1.0_0.yaml"
    transfer.SEMILINEAR_PINN_CONFIG = root / "scripts/burgers/pinn_100d_T1.0_v16_beta10.0.yaml"
    return root


def _seed_all(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def _trajectory_id(method: str, pde: str, configuration: str) -> str:
    return f"{method}_{pde}_d{DIMENSION}_{configuration.lower()}_s{TRAINING_SEED}"


def _trajectory_path(method: str, pde: str, configuration: str) -> Path:
    return RAW_DIR / "trajectories" / f"{_trajectory_id(method, pde, configuration)}.json"


def _checkpoint_dir(method: str, pde: str, configuration: str) -> Path:
    return CHECKPOINT_DIR / method / pde / configuration


def _state_path(method: str, pde: str, configuration: str) -> Path:
    return _checkpoint_dir(method, pde, configuration) / "resume_state.pt"


def _configuration(method: str, pde: str, configuration: str):
    if method not in METHODS or configuration not in ("C0", "C1", "C2"):
        raise KeyError((method, configuration))
    _external_root()
    names = transfer.register_equations()
    cfg, source_path, source_overrides, _run_id, _old_dir, profile = transfer._config(
        method, pde, DIMENSION, TRAINING_SEED, "official", names[pde], "auto"
    )
    cfg = cfg.clone()
    cfg.defrost()
    cfg.NAME = str(WORK_DIR / _trajectory_id(method, pde, configuration))
    cfg.FORCE = True
    # The upstream callback evaluates exact nonterminal truth and is diagnostic
    # only.  It is disabled so truth cannot perturb RNG or stopping decisions.
    cfg.EVAL.FREQ = None
    values = CONFIGURATION_GRID[method][configuration]
    if method == "dpi":
        cfg.PICARD.N = int(values["original_picard_iterations"])
        cfg.TRAIN.N_EPOCHS = int(values["epochs_per_picard"])
        cfg.DATA.kwargs.n_estimate_integral = int(values["integral_mc_samples"])
    else:
        cfg.PICARD.N = 1
        cfg.TRAIN.N_EPOCHS = int(values["original_epochs"])
        cfg.TRAIN.BATCH_SIZE = int(values["batch_size"])
        cfg.METHOD.num_v_samples = int(values["hutchinson_vectors"])
    cfg.freeze()
    return cfg, source_path, source_overrides, profile


def configuration_snapshot(method: str, pde: str, configuration: str) -> dict[str, Any]:
    cfg, source, overrides, profile = _configuration(method, pde, configuration)
    dumped = cfg.dump()
    return {
        "method": method,
        "pde": pde,
        "configuration_id": configuration,
        "expanded_config": dumped,
        "expanded_config_sha256": __import__("hashlib").sha256(dumped.encode()).hexdigest(),
        "source_config": str(source),
        "source_config_sha256": file_sha256(source),
        "source_overrides_before_frontier_changes": overrides,
        "profile": profile,
        "frontier_changes": CONFIGURATION_GRID[method][configuration],
    }


def _fresh_runner(cfg, capture: Path):
    transfer._patch_upstream(capture)
    from picard import PicardRunner

    Path(cfg.NAME).parent.mkdir(parents=True, exist_ok=True)
    old = Path.cwd()
    os.chdir(WORK_DIR)
    try:
        runner = PicardRunner(cfg)
    finally:
        os.chdir(old)
    return runner


def _resume_runner(cfg, state: dict[str, Any], capture: Path):
    transfer._patch_upstream(capture)
    from picard import PicardRunner

    bootstrap = cfg.clone()
    bootstrap.defrost()
    bootstrap.NAME = str(WORK_DIR / (Path(cfg.NAME).name + "_resume_bootstrap"))
    bootstrap.FORCE = True
    bootstrap.freeze()
    old = Path.cwd()
    os.chdir(WORK_DIR)
    try:
        runner = PicardRunner(bootstrap)
    finally:
        os.chdir(old)
    runner.cfg = cfg
    runner.exp_dir = Path(cfg.NAME)
    runner.i = int(state["iteration"])
    runner.previous_dataset_size_args = state.get("previous_dataset_size_args")
    max_iterations = int(CONFIGURATION_GRID["dpi"][state["configuration_id"]]["picard_iterations_max"])
    runner.N = max_iterations if runner.i >= int(cfg.PICARD.N) else int(cfg.PICARD.N)
    runner.n_total_iteration_digits = len(str(max_iterations))
    checkpoint = Path(state["checkpoint"])
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    runner.u_current = transfer._reconstruct_model(cfg, runner.equation, checkpoint, device)
    runner.u_history = [runner.u_current]
    restore_rng(state["rng"])
    return runner


def _write_training_record(
    method: str,
    pde: str,
    configuration: str,
    checkpoints: list[dict[str, Any]],
    *,
    status: str,
    stop_reason: str | None = None,
    peak_gpu_bytes: int | None = None,
) -> dict[str, Any]:
    config = configuration_snapshot(method, pde, configuration)
    result = {
        "revision": REVISION,
        "trajectory_id": _trajectory_id(method, pde, configuration),
        "method": method,
        "pde": pde,
        "dimension": DIMENSION,
        "seed": TRAINING_SEED,
        "configuration_id": configuration,
        "configuration": config,
        "checkpoints": checkpoints,
        "training_seconds": float(checkpoints[-1]["training_seconds"]) if checkpoints else 0.0,
        "selection_seconds": float(checkpoints[-1].get("selection_seconds", 0.0)) if checkpoints else 0.0,
        "stop_reason": stop_reason,
        "peak_gpu_bytes": peak_gpu_bytes,
        "status": status,
        "ground_truth_used_for_training_selection_or_stopping": False,
        "exact_solution_callback_disabled": True,
    }
    write_json_atomic(_trajectory_path(method, pde, configuration), result)
    return result


def train_dpi(pde: str, configuration: str, device: torch.device) -> dict[str, Any]:
    del device  # upstream configuration owns device placement
    ensure_directories()
    output = _trajectory_path("dpi", pde, configuration)
    if output.is_file():
        existing = json.loads(output.read_text(encoding="utf-8"))
        if existing.get("status") == "TRAINED_NO_TRUTH":
            return existing
    _seed_all(TRAINING_SEED)
    cfg, _source, _overrides, _profile = _configuration("dpi", pde, configuration)
    directory = _checkpoint_dir("dpi", pde, configuration)
    directory.mkdir(parents=True, exist_ok=True)
    state_file = _state_path("dpi", pde, configuration)
    capture = directory / "unused_pinn_capture.pt"
    if state_file.is_file():
        state = torch.load(state_file, map_location="cpu", weights_only=False)
        checkpoints = list(state["checkpoints"])
        runner = _resume_runner(cfg, state, capture)
        cumulative = float(state["training_seconds"])
    else:
        runner = _fresh_runner(cfg, capture)
        checkpoints = []
        cumulative = 0.0
    max_iterations = int(CONFIGURATION_GRID["dpi"][configuration]["picard_iterations_max"])
    original_iterations = int(CONFIGURATION_GRID["dpi"][configuration]["original_picard_iterations"])
    target_seconds = [60.0 * value for value in TARGET_MINUTES]
    if torch.cuda.is_available():
        torch.cuda.reset_peak_memory_stats()
    stop_reason = "maximum_predeclared_picard_iterations"
    old = Path.cwd()
    os.chdir(WORK_DIR)
    try:
        while runner.i < max_iterations:
            if runner.i >= original_iterations:
                runner.N = max_iterations
                runner.n_total_iteration_digits = len(str(max_iterations))
            synchronize(torch.device("cuda" if torch.cuda.is_available() else "cpu"))
            started = time.perf_counter()
            ok = runner.run_one()
            synchronize(torch.device("cuda" if torch.cuda.is_available() else "cpu"))
            elapsed = time.perf_counter() - started
            if not ok:
                stop_reason = "upstream_interrupted"
                break
            runner.u_history.append(runner.u_current)
            cumulative += elapsed
            checkpoint = runner.get_model_checkpoint_path(runner.exp_dir, runner.i)
            previous = float(checkpoints[-1]["training_seconds"]) if checkpoints else 0.0
            crossings = [
                target / 60.0 for target in target_seconds
                if previous < target <= cumulative
            ]
            row = {
                "checkpoint_id": f"picard_{runner.i:04d}",
                "algorithm_iteration": int(runner.i),
                "checkpoint": str(checkpoint),
                "checkpoint_sha256": file_sha256(checkpoint),
                "training_seconds": cumulative,
                "selection_seconds": 0.0,
                "main_axis_seconds": cumulative,
                "target_minutes_crossed": crossings,
                "is_original_budget": int(runner.i) == original_iterations,
                "publish": bool(crossings or int(runner.i) == original_iterations),
                "status": "COMPLETE",
            }
            checkpoints.append(row)
            state = {
                "revision": REVISION,
                "method": "dpi",
                "pde": pde,
                "configuration_id": configuration,
                "iteration": int(runner.i),
                "checkpoint": str(checkpoint),
                "previous_dataset_size_args": runner.previous_dataset_size_args,
                "training_seconds": cumulative,
                "checkpoints": checkpoints,
                "rng": rng_state(),
            }
            torch.save(state, state_file)
            _write_training_record(
                "dpi", pde, configuration, checkpoints,
                status="TRAINING_IN_PROGRESS",
            )
            if runner.i >= original_iterations and cumulative >= ACTIVE_TRAINING_CAP_SECONDS:
                stop_reason = "active_training_cap"
                break
    finally:
        os.chdir(old)
    if checkpoints:
        checkpoints[-1]["publish"] = True
        checkpoints[-1]["is_final"] = True
    peak = int(torch.cuda.max_memory_allocated()) if torch.cuda.is_available() else None
    return _write_training_record(
        "dpi", pde, configuration, checkpoints,
        status="TRAINED_NO_TRUTH", stop_reason=stop_reason, peak_gpu_bytes=peak,
    )


def _pinn_validation_path(pde: str) -> Path:
    return RAW_DIR.parent / "validation" / f"pinn_{pde}_d{DIMENSION}.npz"


def _pinn_validation_bank(equation, pde: str, device: torch.device) -> Path:
    path = _pinn_validation_path(pde)
    if path.is_file():
        return path
    saved = rng_state()
    offset = 0 if pde == "b6_fullrank" else 1
    _seed_all(PINN_VALIDATION_SEED + offset)
    with torch.no_grad():
        t = equation.T * torch.rand(PINN_VALIDATION_STATES, 1, device=device)
        x = equation.sample_x(t)
        terminal_t = equation.T * torch.ones(PINN_VALIDATION_STATES, 1, device=device)
        terminal_x = equation.sample_x(terminal_t)
    save_npz_atomic(
        path,
        t=t.cpu().numpy(), x=x.cpu().numpy(),
        terminal_t=terminal_t.cpu().numpy(), terminal_x=terminal_x.cpu().numpy(),
    )
    restore_rng(saved)
    return path


def _pinn_validation(model, equation, pde: str, device: torch.device) -> float:
    from picard.utils import hutchinson_trace_estimation_batch

    path = _pinn_validation_bank(equation, pde, device)
    with np.load(path, allow_pickle=False) as payload:
        t0 = np.asarray(payload["t"])
        x0 = np.asarray(payload["x"])
        terminal_t0 = np.asarray(payload["terminal_t"])
        terminal_x0 = np.asarray(payload["terminal_x"])
    saved = rng_state()
    offset = 0 if pde == "b6_fullrank" else 1
    _seed_all(PINN_VALIDATION_SEED + 100 + offset)
    t = torch.as_tensor(t0, dtype=torch.float64, device=device).requires_grad_(True)
    x = torch.as_tensor(x0, dtype=torch.float64, device=device).requires_grad_(True)
    output = model(torch.cat((t, x), dim=-1))[:, :1]
    u_t = torch.autograd.grad(output.sum(), t, create_graph=True, retain_graph=True)[0]
    u_x = torch.autograd.grad(output.sum(), x, create_graph=True, retain_graph=True)[0]
    laplacian = hutchinson_trace_estimation_batch(
        model, t, x, PINN_VALIDATION_HTE_VECTORS
    )
    residual = equation.pinn_function(t, x, output, u_t, u_x, laplacian)
    terminal_t = torch.as_tensor(terminal_t0, dtype=torch.float64, device=device)
    terminal_x = torch.as_tensor(terminal_x0, dtype=torch.float64, device=device)
    terminal_error = (
        model(torch.cat((terminal_t, terminal_x), dim=-1))[:, :1]
        - equation.g(terminal_x)
    )
    score = torch.sqrt(residual.square().mean() + 10.0 * terminal_error.square().mean())
    synchronize(device)
    value = float(score.detach().cpu())
    restore_rng(saved)
    return value


def _save_pinn_state(
    path: Path,
    model,
    optimizer,
    *,
    epoch: int,
    training_seconds: float,
    selection_seconds: float,
    checkpoints: list[dict[str, Any]],
    configuration: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    torch.save({
        "revision": REVISION,
        "state_dict": model.state_dict(),
        "optimizer_state_dict": optimizer.state_dict(),
        "epoch": epoch,
        "training_seconds": training_seconds,
        "selection_seconds": selection_seconds,
        "checkpoints": checkpoints,
        "configuration_id": configuration,
        "rng": rng_state(),
    }, path)


def _new_pinn(cfg, state: dict[str, Any] | None, capture: Path, device: torch.device):
    runner = _fresh_runner(cfg, capture) if state is None else _resume_pinn_runner(cfg, capture)
    model = runner.get_solution().to(device)
    equation = runner.equation.to(device=device)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    if state is not None:
        model.load_state_dict(state["state_dict"])
        optimizer.load_state_dict(state["optimizer_state_dict"])
        restore_rng(state["rng"])
    return model, equation, optimizer


def _resume_pinn_runner(cfg, capture: Path):
    transfer._patch_upstream(capture)
    from picard import PicardRunner

    bootstrap = cfg.clone()
    bootstrap.defrost()
    bootstrap.NAME = str(WORK_DIR / (Path(cfg.NAME).name + "_resume_bootstrap"))
    bootstrap.FORCE = True
    bootstrap.freeze()
    old = Path.cwd()
    os.chdir(WORK_DIR)
    try:
        runner = PicardRunner(bootstrap)
    finally:
        os.chdir(old)
    runner.cfg = cfg
    runner.exp_dir = Path(cfg.NAME)
    return runner


def train_pinn(pde: str, configuration: str, device: torch.device) -> dict[str, Any]:
    ensure_directories()
    output = _trajectory_path("pinn_hte", pde, configuration)
    if output.is_file():
        existing = json.loads(output.read_text(encoding="utf-8"))
        if existing.get("status") == "TRAINED_NO_TRUTH":
            return existing
    _seed_all(TRAINING_SEED)
    cfg, _source, _overrides, _profile = _configuration("pinn_hte", pde, configuration)
    directory = _checkpoint_dir("pinn_hte", pde, configuration)
    directory.mkdir(parents=True, exist_ok=True)
    resume = _state_path("pinn_hte", pde, configuration)
    state = torch.load(resume, map_location=device, weights_only=False) if resume.is_file() else None
    model, equation, optimizer = _new_pinn(cfg, state, directory / "unused_capture.pt", device)
    values = CONFIGURATION_GRID["pinn_hte"][configuration]
    batch_size = int(values["batch_size"])
    hte_vectors = int(values["hutchinson_vectors"])
    original_epochs = int(values["original_epochs"])
    epoch = int(state["epoch"]) if state else 0
    cumulative = float(state["training_seconds"]) if state else 0.0
    cumulative_selection = float(state["selection_seconds"]) if state else 0.0
    checkpoints = list(state["checkpoints"]) if state else []
    target_seconds = [60.0 * value for value in TARGET_MINUTES]
    max_epochs = 1_000_000
    from picard.utils import hutchinson_trace_estimation_batch

    if torch.cuda.is_available():
        torch.cuda.reset_peak_memory_stats()
    if not checkpoints:
        selection_started = time.perf_counter()
        metric = _pinn_validation(model, equation, pde, device)
        cumulative_selection += time.perf_counter() - selection_started
        checkpoint = directory / "epoch_0000000.pt"
        row = {
            "checkpoint_id": "epoch_0000000", "algorithm_iteration": 0,
            "checkpoint": str(checkpoint), "training_seconds": 0.0,
            "selection_seconds": cumulative_selection,
            "main_axis_seconds": cumulative_selection,
            "validation_metric": metric, "training_loss": None,
            "target_minutes_crossed": [], "is_original_budget": False,
            "publish": True, "status": "COMPLETE",
        }
        checkpoints.append(row)
        _save_pinn_state(
            checkpoint, model, optimizer, epoch=0, training_seconds=0.0,
            selection_seconds=cumulative_selection, checkpoints=checkpoints,
            configuration=configuration,
        )
        row["checkpoint_sha256"] = file_sha256(checkpoint)
        _save_pinn_state(
            resume, model, optimizer, epoch=0, training_seconds=0.0,
            selection_seconds=cumulative_selection, checkpoints=checkpoints,
            configuration=configuration,
        )
    stop_reason = "active_training_cap"
    last_loss: float | None = None
    synchronize(device)
    segment_started = time.perf_counter()
    while epoch < max_epochs:
        epoch += 1
        t = equation.T * torch.rand(batch_size, 1, device=device)
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
        laplacian = hutchinson_trace_estimation_batch(model, t, x, hte_vectors)
        loss_in = equation.pinn_function(
            t, x, output_value, u_t, u_x, laplacian
        ).square().mean()
        terminal_t = equation.T * torch.ones_like(t)
        terminal_x = equation.sample_x(terminal_t)
        terminal_value = model(torch.cat((terminal_t, terminal_x), dim=-1))[:, :1]
        loss_terminal = (terminal_value - equation.g(terminal_x)).square().mean()
        loss = loss_in + 10.0 * loss_terminal
        if not torch.isfinite(loss):
            stop_reason = "nonfinite_training_loss"
            break
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        last_loss = float(loss.detach().cpu())
        if epoch % PINN_CHECK_INTERVAL_EPOCHS != 0 and epoch != original_epochs:
            continue
        synchronize(device)
        cumulative += time.perf_counter() - segment_started
        previous = float(checkpoints[-1]["training_seconds"])
        crossings = [
            target / 60.0 for target in target_seconds
            if previous < target <= cumulative
        ]
        should_checkpoint = bool(
            crossings or epoch == original_epochs
            or cumulative >= ACTIVE_TRAINING_CAP_SECONDS
        )
        if should_checkpoint:
            selection_started = time.perf_counter()
            metric = _pinn_validation(model, equation, pde, device)
            cumulative_selection += time.perf_counter() - selection_started
            checkpoint = directory / f"epoch_{epoch:07d}.pt"
            row = {
                "checkpoint_id": f"epoch_{epoch:07d}",
                "algorithm_iteration": epoch,
                "checkpoint": str(checkpoint),
                "training_seconds": cumulative,
                "selection_seconds": cumulative_selection,
                "main_axis_seconds": cumulative + cumulative_selection,
                "validation_metric": metric,
                "training_loss": last_loss,
                "target_minutes_crossed": crossings,
                "is_original_budget": epoch == original_epochs,
                "publish": True,
                "status": "COMPLETE",
            }
            checkpoints.append(row)
            _save_pinn_state(
                checkpoint, model, optimizer, epoch=epoch,
                training_seconds=cumulative, selection_seconds=cumulative_selection,
                checkpoints=checkpoints, configuration=configuration,
            )
            row["checkpoint_sha256"] = file_sha256(checkpoint)
            _save_pinn_state(
                resume, model, optimizer, epoch=epoch,
                training_seconds=cumulative, selection_seconds=cumulative_selection,
                checkpoints=checkpoints, configuration=configuration,
            )
            _write_training_record(
                "pinn_hte", pde, configuration, checkpoints,
                status="TRAINING_IN_PROGRESS",
            )
            eligible = [
                float(item["validation_metric"]) for item in checkpoints
                if int(item["algorithm_iteration"]) > 0
            ]
            if epoch >= original_epochs and plateau_reached(cumulative, eligible):
                stop_reason = "predeclared_validation_plateau"
                break
            if epoch >= original_epochs and cumulative >= ACTIVE_TRAINING_CAP_SECONDS:
                stop_reason = "active_training_cap"
                break
        segment_started = time.perf_counter()
    synchronize(device)
    if checkpoints:
        checkpoints[-1]["is_final"] = True
    peak = int(torch.cuda.max_memory_allocated()) if torch.cuda.is_available() else None
    return _write_training_record(
        "pinn_hte", pde, configuration, checkpoints,
        status="TRAINED_NO_TRUTH", stop_reason=stop_reason, peak_gpu_bytes=peak,
    )


def train(method: str, pde: str, configuration: str, device: torch.device) -> dict[str, Any]:
    if method == "dpi":
        return train_dpi(pde, configuration, device)
    if method == "pinn_hte":
        return train_pinn(pde, configuration, device)
    raise KeyError(method)


def _bank(pde: str):
    path = historical_results_root() / "inputs" / f"{pde}_d{DIMENSION}_test_bank.npz"
    with np.load(path, allow_pickle=False) as payload:
        return (
            np.asarray(payload["t"]), np.asarray(payload["x"]),
            np.asarray(payload["truth_u"]), path,
        )


def _model_for_checkpoint(method: str, cfg, checkpoint: Path, device: torch.device):
    from picard import PicardRunner
    from picard.solution_jac import PicardSolutionGradientWrapper

    equation = PicardRunner.get_equation(cfg).to(device=device)
    if method == "dpi":
        return transfer._reconstruct_model(cfg, equation, checkpoint, device), equation

    class Runner(SimpleNamespace):
        def get_solution_plain(self):
            cls = PicardRunner.get_solution_cls(self.cfg)
            return cls.construct_solution(self)

    runner = Runner(equation=equation, cfg=cfg)
    model = PicardSolutionGradientWrapper.construct_solution(runner).to(device)
    payload = torch.load(checkpoint, map_location=device, weights_only=False)
    model.load_state_dict(payload["state_dict"])
    return model.eval(), equation


def evaluate(
    method: str, pde: str, configuration: str, device: torch.device
) -> list[dict[str, Any]]:
    gate = STATUS_DIR / "all_training_complete.json"
    if not gate.is_file():
        raise RuntimeError("global truth gate is closed: all 18 trajectories must stop first")
    record_path = _trajectory_path(method, pde, configuration)
    trajectory = json.loads(record_path.read_text(encoding="utf-8"))
    if trajectory.get("status") not in {"TRAINED_NO_TRUTH", "COMPLETE"}:
        raise RuntimeError("trajectory has not stopped")
    cfg, _source, _overrides, _profile = _configuration(method, pde, configuration)
    t, x, truth, bank_path = _bank(pde)
    rows: list[dict[str, Any]] = []
    for checkpoint in trajectory["checkpoints"]:
        if not checkpoint.get("publish", False):
            continue
        output_path = (
            PREDICTION_DIR / method / pde / configuration
            / f"{checkpoint['checkpoint_id']}.npz"
        )
        metric_path = output_path.with_suffix(".json")
        if metric_path.is_file():
            rows.append(json.loads(metric_path.read_text(encoding="utf-8")))
            continue
        model, _equation = _model_for_checkpoint(
            method, cfg, Path(checkpoint["checkpoint"]), device
        )
        predictions: list[np.ndarray] = []
        synchronize(device)
        started = time.perf_counter()
        with torch.no_grad():
            for begin in range(0, len(t), 128):
                end = min(begin + 128, len(t))
                tt = torch.as_tensor(t[begin:end, None], dtype=torch.float64, device=device)
                xx = torch.as_tensor(x[begin:end], dtype=torch.float64, device=device)
                predictions.append(
                    model(torch.cat((tt, xx), dim=-1))[:, 0].detach().cpu().numpy()
                )
        synchronize(device)
        inference_seconds = time.perf_counter() - started
        prediction = np.concatenate(predictions).astype(np.float64)
        metrics = value_metrics(prediction, truth)
        save_npz_atomic(
            output_path, t=t, x=x, truth=truth, prediction=prediction,
        )
        row = {
            "revision": REVISION,
            "trajectory_id": trajectory["trajectory_id"],
            "method": method, "pde": pde, "dimension": DIMENSION,
            "configuration_id": configuration, "seed": TRAINING_SEED,
            "checkpoint_id": checkpoint["checkpoint_id"],
            "algorithm_iteration": checkpoint["algorithm_iteration"],
            "training_seconds": checkpoint["training_seconds"],
            "selection_seconds": checkpoint.get("selection_seconds", 0.0),
            "main_axis_seconds": checkpoint["main_axis_seconds"],
            "inference_seconds": inference_seconds,
            "validation_metric": checkpoint.get("validation_metric"),
            "is_original_budget": checkpoint.get("is_original_budget", False),
            "is_final": checkpoint.get("is_final", False),
            "checkpoint": checkpoint["checkpoint"],
            "checkpoint_sha256": checkpoint["checkpoint_sha256"],
            "test_bank": str(bank_path),
            "test_bank_sha256": file_sha256(bank_path),
            "test_input_sha256": array_sha256(t, x),
            "truth_value_sha256": array_sha256(truth),
            "prediction_file": str(output_path),
            "prediction_sha256": file_sha256(output_path),
            "nonfinite_count": int((~np.isfinite(prediction)).sum()),
            "status": "COMPLETE",
            **metrics,
        }
        write_json_atomic(metric_path, row)
        rows.append(row)
    trajectory["status"] = "COMPLETE"
    trajectory["truth_opened_only_after_global_training_stop"] = True
    trajectory["evaluation_rows"] = len(rows)
    write_json_atomic(record_path, trajectory)
    return rows


__all__ = [
    "configuration_snapshot", "evaluate", "train", "train_dpi", "train_pinn",
]
