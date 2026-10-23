"""Isolated wrappers around the mature Tangent-Actor implementations."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import math
from pathlib import Path
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.original_pde_benchmark_comparison.tangent_solver import (
    ActorConfig,
    antithetic_normals,
    evaluate_policy as evaluate_gradient_or_value_policy,
    fit_actor,
    frozen,
    generate_labels as generate_original_labels,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.common import (
    synchronize,
)
from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis.training import (
    FitConfig,
)

from .protocol import (
    ACTOR_UPDATES,
    CACHE_DIR,
    CHECKPOINT_DIR,
    GRADIENT_ACTOR_CONFIG,
    JOINT_ACTOR_CONFIG,
    PREDICTION_DIR,
    QUERY_COUNT,
    RAW_DIR,
    REFERENCE_VALIDATION,
    REVISION,
    STAGE_A_AUDIT,
    TA_EVALUATION_PATHS,
    TA_INFERENCE_REPLICATES,
    TA_PRIMARY_EULER_STEPS,
    TRAINING_SEEDS,
    array_sha256,
    ensure_directories,
    file_sha256,
    query_path,
    require_stage_a,
    require_stage_b,
    stable_seed,
    write_json_atomic,
    write_npz_atomic,
)


SMOKE_GRADIENT = ActorConfig(
    training_points=32,
    validation_points=16,
    paths_per_label=8,
    label_time_steps=4,
    optimizer_steps=4,
    width=16,
    depth=2,
    activation="silu",
    learning_rate=2e-4,
    batch_size=16,
    parent_batch=4,
    path_chunk=4,
    checkpoint_interval=2,
    gradient_clip=10.0,
    training_distribution="function",
    label_estimator="mcv",
)

SMOKE_JOINT = FitConfig(
    training_points=32,
    validation_points=16,
    paths_per_label=8,
    label_time_steps=4,
    optimizer_steps=4,
    width=16,
    depth=2,
    learning_rate=2e-4,
    batch_size=16,
    parent_batch=4,
    path_chunk=4,
    checkpoint_interval=2,
    gradient_clip=10.0,
)


def _roots(pilot: bool) -> dict[str, Path]:
    tag = "pilot" if pilot else "formal"
    return {
        "raw": RAW_DIR / tag / "tangent_actor",
        "cache": CACHE_DIR / tag / "tangent_actor",
        "checkpoint": CHECKPOINT_DIR / tag / "tangent_actor",
        "prediction": PREDICTION_DIR / tag / "tangent_actor",
    }


def _value_mcv_rollout_chunk(
    benchmark,
    actor,
    t: torch.Tensor,
    x: torch.Tensor,
    normals: torch.Tensor,
    *,
    use_mcv: bool,
):
    """Linear-tangent payoff with a value-Actor autograd martingale CV.

    For a value Actor ``r_theta``, ``sigma*grad_x(r_theta)`` is predictable at
    each left endpoint.  Its stochastic integral has zero expectation and is
    therefore a valid control variate without changing the value-only Actor.
    """
    batch, paths, steps, dimension = normals.shape
    times = t.detach().reshape(-1)
    dt = torch.clamp(benchmark.T - times, min=0.0) / steps
    state = x.detach()[:, None, :].expand(batch, paths, dimension).clone()
    weight = torch.ones((batch, paths), dtype=x.dtype, device=x.device)
    integral = torch.zeros_like(weight)
    martingale = torch.zeros_like(weight)
    c_min, c_max = float("inf"), -float("inf")
    for step in range(steps):
        step_t = times + step * dt
        flat_t = step_t[:, None].expand(batch, paths).reshape(-1)
        flat_x = state.reshape(-1, dimension)
        if use_mcv:
            with torch.enable_grad():
                differentiable_x = flat_x.detach().requires_grad_(True)
                r_graph, _ = actor.predict(flat_t, differentiable_x)
                if r_graph is None:
                    raise RuntimeError("value Actor did not return r")
                gradient = torch.autograd.grad(r_graph.sum(), differentiable_x)[0]
            r = r_graph.detach()
            q_cv = benchmark.sigma * gradient.detach()
        else:
            with torch.no_grad():
                r, _ = actor.predict(flat_t, flat_x)
            if r is None:
                raise RuntimeError("value Actor did not return r")
            q_cv = None
        c, b, h = benchmark.tangent_coefficients(r, None)
        c = c.reshape(batch, paths)
        b = b.reshape(batch, paths, dimension)
        h = h.reshape(batch, paths)
        dw = torch.sqrt(dt)[:, None, None] * normals[:, :, step, :]
        if q_cv is not None:
            martingale = martingale + weight * (
                q_cv.reshape(batch, paths, dimension) * dw
            ).sum(dim=-1)
        integral = integral + weight * h * dt[:, None]
        weight = weight * torch.exp(c * dt[:, None])
        drift = benchmark.base_drift(flat_t, flat_x).reshape(batch, paths, dimension)
        state = state + (drift + benchmark.sigma * b) * dt[:, None, None]
        state = state + benchmark.sigma * dw
        c_min = min(c_min, float(c.min().item()))
        c_max = max(c_max, float(c.max().item()))
    with torch.no_grad():
        plain = weight * benchmark.terminal(state) + integral
        corrected = plain - martingale if use_mcv else plain
    nonfinite = int((~torch.isfinite(plain)).sum().item())
    nonfinite += int((~torch.isfinite(corrected)).sum().item())
    nonfinite += int((~torch.isfinite(weight)).sum().item())
    return plain, corrected, weight, martingale, {
        "c_min": c_min,
        "c_max": c_max,
        "nonfinite_count": nonfinite,
    }


def _value_aware_generate_labels(
    benchmark,
    source,
    t: np.ndarray,
    x: np.ndarray,
    *,
    paths: int,
    time_steps: int,
    seed: int,
    device: torch.device,
    parent_batch: int,
    path_chunk: int,
    estimator: str,
):
    if benchmark.actor_kind == "gradient":
        return generate_original_labels(
            benchmark, source, t, x, paths=paths, time_steps=time_steps,
            seed=seed, device=device, parent_batch=parent_batch,
            path_chunk=path_chunk, estimator=estimator,
        )
    if estimator not in {"plain", "mcv"}:
        raise KeyError(estimator)
    use_mcv = estimator == "mcv"
    times = np.asarray(t, dtype=np.float64).reshape(-1)
    states = np.asarray(x, dtype=np.float64).reshape(-1, benchmark.d)
    labels = np.empty(len(times), dtype=np.float64)
    generator = torch.Generator(device=device).manual_seed(seed)
    weight_sum = weight_square_sum = 0.0
    martingale_sum = martingale_square_sum = 0.0
    c_min, c_max = float("inf"), -float("inf")
    nonfinite = 0
    synchronize(device)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    with frozen(source):
        for begin in range(0, len(times), parent_batch):
            end = min(begin + parent_batch, len(times))
            bt = torch.as_tensor(times[begin:end], dtype=torch.float32, device=device)
            bx = torch.as_tensor(states[begin:end], dtype=torch.float32, device=device)
            value_sum = torch.zeros(end - begin, dtype=torch.float64, device=device)
            for _ in range(0, paths, path_chunk):
                normals = antithetic_normals(
                    end - begin, path_chunk, time_steps, benchmark.d,
                    generator=generator, device=device,
                )
                _, corrected, weight, martingale, diagnostics = _value_mcv_rollout_chunk(
                    benchmark, source, bt, bx, normals, use_mcv=use_mcv,
                )
                value_sum += corrected.double().sum(dim=1)
                weight_sum += float(weight.double().sum().item())
                weight_square_sum += float(weight.double().square().sum().item())
                martingale_sum += float(martingale.double().sum().item())
                martingale_square_sum += float(martingale.double().square().sum().item())
                c_min = min(c_min, diagnostics["c_min"])
                c_max = max(c_max, diagnostics["c_max"])
                nonfinite += diagnostics["nonfinite_count"]
            labels[begin:end] = (value_sum / paths).cpu().numpy()
    synchronize(device)
    elapsed = time.perf_counter() - started
    count = len(times) * paths
    if nonfinite or not np.all(np.isfinite(labels)):
        raise FloatingPointError("nonfinite value-Actor MCV labels")
    martingale_mean = martingale_sum / count
    martingale_variance = max(
        martingale_square_sum / count - martingale_mean**2, 0.0
    )
    return labels, {
        "kind": "differentiable_linear_tangent_feynman_kac",
        "actor_target": "value",
        "estimator": "value_actor_autograd_direct_alpha1_mcv" if use_mcv else "plain_fk",
        "alpha": 1.0 if use_mcv else 0.0,
        "points": len(times),
        "paths": paths,
        "time_steps": time_steps,
        "seed": seed,
        "source_actor_frozen": True,
        "control_integrand": "sqrt(2)*grad_x(value_actor) at predictable left endpoints",
        "left_endpoint_predictable_martingale": use_mcv,
        "analytic_nonterminal_reference_used": False,
        "ground_truth_used": False,
        "weight_ess_fraction": float(
            weight_sum**2 / max(count * weight_square_sum, np.finfo(float).tiny)
        ),
        "c_min": c_min,
        "c_max": c_max,
        "martingale_path_mean": float(martingale_mean),
        "martingale_path_std": float(math.sqrt(martingale_variance)),
        "label_rms": float(np.sqrt(np.mean(np.square(labels)))),
        "nonfinite_count": nonfinite,
        "wall_seconds": elapsed,
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device))
        if device.type == "cuda" else None,
    }


def _evaluate_value_mcv_policy(
    benchmark,
    actor,
    t: np.ndarray,
    x: np.ndarray,
    *,
    paths: int,
    time_steps: int,
    seed: int,
    device: torch.device,
    parent_batch: int = 8,
    path_chunk: int = 32,
):
    times = np.asarray(t, dtype=np.float64).reshape(-1)
    states = np.asarray(x, dtype=np.float64).reshape(-1, benchmark.d)
    plain_mean = np.empty(len(times), dtype=np.float64)
    mcv_mean = np.empty(len(times), dtype=np.float64)
    plain_variance = np.empty(len(times), dtype=np.float64)
    mcv_variance = np.empty(len(times), dtype=np.float64)
    generator = torch.Generator(device=device).manual_seed(seed)
    weight_sum = weight_square_sum = 0.0
    martingale_sum = martingale_square_sum = 0.0
    nonfinite = 0
    synchronize(device)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    with frozen(actor):
        for begin in range(0, len(times), parent_batch):
            end = min(begin + parent_batch, len(times))
            bt = torch.as_tensor(times[begin:end], dtype=torch.float32, device=device)
            bx = torch.as_tensor(states[begin:end], dtype=torch.float32, device=device)
            plain_pairs, mcv_pairs = [], []
            for _ in range(0, paths, path_chunk):
                normals = antithetic_normals(
                    end - begin, path_chunk, time_steps, benchmark.d,
                    generator=generator, device=device,
                )
                plain, corrected, weight, martingale, diagnostics = _value_mcv_rollout_chunk(
                    benchmark, actor, bt, bx, normals, use_mcv=True,
                )
                half = path_chunk // 2
                plain_pairs.append(0.5 * (plain[:, :half] + plain[:, half:]))
                mcv_pairs.append(0.5 * (corrected[:, :half] + corrected[:, half:]))
                weight_sum += float(weight.double().sum().item())
                weight_square_sum += float(weight.double().square().sum().item())
                martingale_sum += float(martingale.double().sum().item())
                martingale_square_sum += float(martingale.double().square().sum().item())
                nonfinite += diagnostics["nonfinite_count"]
            local_plain = torch.cat(plain_pairs, dim=1).double()
            local_mcv = torch.cat(mcv_pairs, dim=1).double()
            plain_mean[begin:end] = local_plain.mean(1).cpu().numpy()
            mcv_mean[begin:end] = local_mcv.mean(1).cpu().numpy()
            plain_variance[begin:end] = local_plain.var(1, unbiased=True).cpu().numpy()
            mcv_variance[begin:end] = local_mcv.var(1, unbiased=True).cpu().numpy()
    synchronize(device)
    elapsed = time.perf_counter() - started
    count = len(times) * paths
    martingale_mean = martingale_sum / count
    martingale_variance = max(martingale_square_sum / count - martingale_mean**2, 0.0)
    return {
        "plain": plain_mean,
        "mcv": mcv_mean,
        "plain_pair_variance": plain_variance,
        "mcv_pair_variance": mcv_variance,
    }, {
        "paths": paths,
        "antithetic_pairs": paths // 2,
        "time_steps": time_steps,
        "seed": seed,
        "mcv": True,
        "mcv_alpha": 1.0,
        "control_integrand": "sqrt(2)*grad_x(value_actor)",
        "weight_ess_fraction": float(
            weight_sum**2 / max(count * weight_square_sum, np.finfo(float).tiny)
        ),
        "martingale_path_mean": float(martingale_mean),
        "martingale_path_std": float(math.sqrt(martingale_variance)),
        "mean_plain_pair_variance": float(np.mean(plain_variance)),
        "mean_mcv_pair_variance": float(np.mean(mcv_variance)),
        "variance_reduction": float(
            np.mean(plain_variance) / max(np.mean(mcv_variance), 1e-300)
        ),
        "mean_plain_mc_se": float(np.mean(np.sqrt(plain_variance / (paths // 2)))),
        "mean_mcv_mc_se": float(np.mean(np.sqrt(mcv_variance / (paths // 2)))),
        "nonfinite_count": nonfinite,
        "wall_seconds": elapsed,
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device))
        if device.type == "cuda" else None,
    }


def _configure_standard(pilot: bool):
    from invariant_region_mlp.experiments.tangent_actor_icml_final_suite import (
        run_standard_actor as impl,
    )

    roots = _roots(pilot)
    impl.RAW_DIR = roots["raw"]
    impl.CACHE_DIR = roots["cache"]
    impl.CHECKPOINT_DIR = roots["checkpoint"]
    impl.PREDICTION_DIR = roots["prediction"]
    impl.REFERENCE_DIR = roots["raw"] / "references"
    impl.stable_seed = stable_seed
    impl.write_json = write_json_atomic
    impl.sha256_file = file_sha256
    impl.generate_labels = _value_aware_generate_labels
    for value in roots.values():
        value.mkdir(parents=True, exist_ok=True)
    return impl


def _train_allen_value_mcv(impl, d: int, seed: int, device: torch.device, pilot: bool):
    """Mature value-Actor loop with the preregistered MCV label estimator."""
    impl._dirs()
    base = impl.SMOKE if pilot else impl.MATURE
    config = ActorConfig(**{
        **asdict(base),
        "training_distribution": impl._purpose("allen_cahn"),
        "label_estimator": "mcv",
    })
    benchmark = impl.make_benchmark("allen_cahn", d)
    tag = "smoke" if pilot else "mature"
    run_id = impl._run_id("allen_cahn", d, seed, tag)
    output = impl.RAW_DIR / f"training_{run_id}.json"
    actor = impl._initial_actor(run_id, benchmark, config, device)
    phases: list[dict[str, Any]] = []
    if output.is_file():
        existing = json.loads(output.read_text(encoding="utf-8"))
        if existing["configuration"] != asdict(config):
            raise RuntimeError("existing Allen value-MCV configuration differs")
        phases = existing.get("phases", [])
        if phases:
            actor = impl._load_actor(impl._checkpoint(run_id, len(phases)), benchmark, device)
        if existing.get("status") == "complete":
            return existing
    started = time.perf_counter()
    for iteration in range(len(phases) + 1, (1 if pilot else ACTOR_UPDATES) + 1):
        phase_started = time.perf_counter()
        arrays, label_meta = impl._labels(
            run_id, "allen_cahn", benchmark, actor, iteration, config, seed, device
        )
        fitted, fit_meta, trace = fit_actor(
            benchmark,
            config,
            actor,
            arrays["train_t"],
            arrays["train_x"],
            arrays["train_y"],
            arrays["valid_t"],
            arrays["valid_x"],
            arrays["valid_y"],
            seed=stable_seed(run_id, iteration, "fit"),
            device=device,
        )
        checkpoint = impl._checkpoint(run_id, iteration)
        impl._save_actor(checkpoint, fitted, benchmark, config, iteration)
        phases.append({
            "iteration": iteration,
            "checkpoint": str(checkpoint),
            "checkpoint_sha256": file_sha256(checkpoint),
            "phase_seconds": time.perf_counter() - phase_started,
            "labels": label_meta,
            "training": fit_meta,
            "trace": trace,
        })
        actor = fitted
        write_json_atomic(output, {
            "run_id": run_id,
            "method": "terminal_aware_value_tangent_actor_autograd_mcv",
            "family": "allen_cahn",
            "dimension": d,
            "seed": seed,
            "configuration": asdict(config),
            "phases": phases,
            "status": "training_in_progress",
            "full_ambient_input": True,
            "oracle_projection": False,
            "ground_truth_used_for_training_or_selection": False,
        })
    actors = [impl._initial_actor(run_id, benchmark, config, device)] + [
        impl._load_actor(impl._checkpoint(run_id, iteration), benchmark, device)
        for iteration in range(1, (1 if pilot else ACTOR_UPDATES) + 1)
    ]
    fixed_points = [
        impl._fixed_point(
            run_id, "allen_cahn", benchmark, f"a{iteration}", candidate,
            seed, device, pilot,
        )
        for iteration, candidate in enumerate(actors)
    ]
    eligible = [item for item in fixed_points if item["nonfinite_count"] == 0]
    if not eligible:
        raise FloatingPointError("all Allen value Actors failed the reference-free selector")
    selected = min(eligible, key=lambda item: item["normalized_total_inconsistency"])
    selected_index = int(selected["actor"][1:])
    record = {
        "run_id": run_id,
        "method": "terminal_aware_value_tangent_actor_autograd_mcv",
        "family": "allen_cahn",
        "dimension": d,
        "seed": seed,
        "configuration": asdict(config),
        "phases": phases,
        "selection": {
            "criterion": "minimum independent noise-corrected fixed-point inconsistency",
            "rows": fixed_points,
            "selected_actor": selected["actor"],
            "selected_checkpoint": None if selected_index == 0 else str(
                impl._checkpoint(run_id, selected_index)
            ),
            "ground_truth_used": False,
        },
        "control_variate": "sqrt(2)*grad_x(value_actor) predictable martingale",
        "full_ambient_input": True,
        "oracle_projection": False,
        "ground_truth_used_for_training_or_selection": False,
        "invocation_seconds": time.perf_counter() - started,
        "status": "complete",
    }
    write_json_atomic(output, record)
    return record


def _configure_b6(pilot: bool):
    from invariant_region_mlp.experiments.main_highdim_expansion import (
        run_b6_actor as impl,
    )

    roots = _roots(pilot)
    impl.RAW_DIR = roots["raw"]
    impl.CACHE_DIR = roots["cache"]
    impl.CHECKPOINT_DIR = roots["checkpoint"]
    impl.PREDICTION_DIR = roots["prediction"]
    impl.REVISION = REVISION
    impl.VALIDATION_PATH = STAGE_A_AUDIT
    impl.stable_seed = stable_seed
    impl.write_json_atomic = write_json_atomic
    impl.file_sha256 = file_sha256
    impl.ensure_directories = ensure_directories
    impl.ACTOR_UPDATES = 1 if pilot else ACTOR_UPDATES
    impl.GRADIENT_ACTOR_CONFIG = SMOKE_GRADIENT if pilot else GRADIENT_ACTOR_CONFIG
    for value in roots.values():
        value.mkdir(parents=True, exist_ok=True)
    return impl


def _configure_dense(pilot: bool, seed: int):
    from invariant_region_mlp.experiments.main_highdim_expansion import (
        run_dense_actor as impl,
    )

    roots = _roots(pilot)
    impl.RAW_DIR = roots["raw"]
    impl.CACHE_DIR = roots["cache"]
    impl.CHECKPOINT_DIR = roots["checkpoint"]
    impl.PREDICTION_DIR = roots["prediction"]
    impl.REVISION = REVISION
    impl.VALIDATION_PATH = STAGE_A_AUDIT
    impl.stable_seed = stable_seed
    impl.write_json_atomic = write_json_atomic
    impl.file_sha256 = file_sha256
    impl.ensure_directories = ensure_directories
    impl.ACTOR_UPDATES = 1 if pilot else ACTOR_UPDATES
    impl.JOINT_ACTOR_CONFIG = SMOKE_JOINT if pilot else JOINT_ACTOR_CONFIG
    # The source suite had one formal seed; this experiment preregisters three.
    impl.TRAINING_SEED = int(seed)
    for value in roots.values():
        value.mkdir(parents=True, exist_ok=True)
    return impl


def train(
    pde: str, d: int, seed: int, device: torch.device, *, pilot: bool = False
) -> dict[str, Any]:
    if not pilot and seed not in TRAINING_SEEDS:
        raise ValueError(seed)
    require_stage_a()
    if not pilot:
        require_stage_b()
    ensure_directories()
    if pde == "hjb_log":
        impl = _configure_standard(pilot)
        record = impl.train(pde, d, seed, device, smoke=pilot)
    elif pde == "allen_cahn":
        impl = _configure_standard(pilot)
        record = _train_allen_value_mcv(impl, d, seed, device, pilot)
    elif pde == "b6_fullrank":
        impl = _configure_b6(pilot)
        record = impl.train(pde, d, seed, device)
    elif pde == "dense_fourier":
        impl = _configure_dense(pilot, seed)
        record = impl.train(d, seed, device)
    else:
        raise KeyError(pde)
    record = dict(record)
    record["study_revision"] = REVISION
    record["formal_protocol"] = not pilot
    record["status"] = "pilot_complete" if pilot else record.get("status", "complete")
    return record


def _selected_standard(impl, pde: str, d: int, seed: int, device: torch.device, pilot: bool):
    tag = "smoke" if pilot else "mature"
    run_id = impl._run_id(pde, d, seed, tag)
    training_path = impl.RAW_DIR / f"training_{run_id}.json"
    training = json.loads(training_path.read_text(encoding="utf-8"))
    benchmark = impl.make_benchmark(pde, d)
    config = ActorConfig(**training["configuration"])
    selected = training["selection"]["selected_actor"]
    if selected == "a0":
        actor = impl._initial_actor(run_id, benchmark, config, device)
    else:
        actor = impl._load_actor(impl._checkpoint(run_id, int(selected[1:])), benchmark, device)
    return run_id, training_path, selected, benchmark, actor, evaluate_gradient_or_value_policy


def _selected_b6(impl, pde: str, d: int, seed: int, device: torch.device):
    run_id = impl._run_id(pde, d, seed)
    training_path = impl._training_path(run_id)
    training = json.loads(training_path.read_text(encoding="utf-8"))
    benchmark = impl.make_benchmark(pde, d)
    config = ActorConfig(**training["configuration"])
    selected = training["selection"]["selected_actor"]
    if selected == "a0":
        actor = impl._initial_actor(run_id, benchmark, config, device)
    else:
        actor = impl._load_actor(impl._checkpoint(run_id, int(selected[1:])), benchmark, device)
    return run_id, training_path, selected, benchmark, actor, evaluate_gradient_or_value_policy


def _selected_dense(impl, d: int, seed: int, device: torch.device):
    run_id = impl._run_id(d, seed)
    training_path = impl._training_path(run_id)
    training = json.loads(training_path.read_text(encoding="utf-8"))
    benchmark = impl.make_dense_fourier(d)
    selected = training["selection"]["selected_actor"]
    actor = impl._load_actor(impl._checkpoint(run_id, int(selected[1:])), benchmark, device)
    return run_id, training_path, selected, benchmark, actor, impl.evaluate_policy


def evaluate(
    pde: str,
    d: int,
    seed: int,
    replicate: int,
    device: torch.device,
    *,
    pilot: bool = False,
    time_steps: int | None = None,
) -> dict[str, Any]:
    if not pilot and not 0 <= replicate < TA_INFERENCE_REPLICATES:
        raise ValueError(replicate)
    require_stage_a()
    if not pilot:
        require_stage_b()
    steps = int(time_steps or (4 if pilot else TA_PRIMARY_EULER_STEPS))
    paths = 8 if pilot else TA_EVALUATION_PATHS
    tag = "pilot" if pilot else "formal"
    result_id = f"ta__{pde}__d{d}__s{seed}__rep{replicate:02d}__steps{steps}__{tag}"
    output_path = RAW_DIR / tag / "tangent_actor_evaluations" / f"{result_id}.json"
    if output_path.is_file():
        existing = json.loads(output_path.read_text(encoding="utf-8"))
        prediction_file = Path(existing.get("prediction_file", ""))
        if existing.get("status") in {"complete", "pilot_complete"} and prediction_file.is_file():
            return existing
    if pde in {"hjb_log", "allen_cahn"}:
        impl = _configure_standard(pilot)
        selected_data = _selected_standard(impl, pde, d, seed, device, pilot)
    elif pde == "b6_fullrank":
        impl = _configure_b6(pilot)
        selected_data = _selected_b6(impl, pde, d, seed, device)
    elif pde == "dense_fourier":
        impl = _configure_dense(pilot, seed)
        selected_data = _selected_dense(impl, d, seed, device)
    else:
        raise KeyError(pde)
    run_id, training_path, selected, benchmark, actor, evaluator = selected_data
    t = np.load(query_path(pde, d) / "t.npy", allow_pickle=False)
    x = np.load(query_path(pde, d) / "x.npy", allow_pickle=False)
    if pilot:
        t, x = t[:4], x[:4]
    if pde == "allen_cahn":
        output, metadata = _evaluate_value_mcv_policy(
            benchmark,
            actor,
            t,
            x,
            paths=paths,
            time_steps=steps,
            seed=stable_seed("ta_evaluation", pde, d, seed, replicate, steps,
                             "pilot" if pilot else "formal"),
            device=device,
            parent_batch=4 if pilot else 8,
            path_chunk=4 if pilot else 32,
        )
    else:
        output, metadata = evaluator(
            benchmark,
            actor,
            t,
            x,
            paths=paths,
            time_steps=steps,
            seed=stable_seed("ta_evaluation", pde, d, seed, replicate, steps,
                             "pilot" if pilot else "formal"),
            device=device,
            parent_batch=4 if pilot else 8,
            path_chunk=4 if pilot else 32,
        )
    estimator = "mcv"
    prediction = np.asarray(output[estimator], dtype=np.float64)
    pair_variance = np.asarray(output[f"{estimator}_pair_variance"], dtype=np.float64)
    # Truth remains unavailable until Actor selection and inference are complete.
    reference = np.load(query_path(pde, d) / "reference.npy", allow_pickle=False)[:len(t)]
    prediction_path = PREDICTION_DIR / tag / "tangent_actor" / f"{result_id}.npz"
    write_npz_atomic(
        prediction_path,
        t=t,
        x=x,
        prediction=prediction,
        reference=reference,
        pair_variance=pair_variance,
    )
    error = prediction - reference
    record = {
        "revision": REVISION,
        "run_id": result_id,
        "training_run_id": run_id,
        "method": "tangent_actor",
        "status": "pilot_complete" if pilot else "complete",
        "pde": pde,
        "dimension": d,
        "training_seed": seed,
        "inference_replicate": replicate,
        "selected_actor": selected,
        "training_artifact": str(training_path),
        "training_artifact_sha256": file_sha256(training_path),
        "paths": paths,
        "time_steps": steps,
        "estimator": (
            "value_actor_autograd_direct_alpha1_mcv"
            if pde == "allen_cahn" else "direct_alpha1_mcv"
        ),
        "mcv_available": True,
        "prediction_file": str(prediction_path),
        "query_input_sha256": array_sha256(t, x),
        "prediction_sha256": array_sha256(prediction),
        "rmse": float(np.sqrt(np.mean(np.square(error)))),
        "nonfinite_count": int((~np.isfinite(prediction)).sum()),
        "inference_metadata": metadata,
        "ground_truth_used_for_training_or_selection": False,
    }
    write_json_atomic(output_path, record)
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("train", "evaluate"))
    parser.add_argument("--pde", required=True)
    parser.add_argument("--dimension", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--replicate", type=int, default=0)
    parser.add_argument("--time-steps", type=int)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--pilot", action="store_true")
    args = parser.parse_args()
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA device requested but unavailable")
    if args.action == "train":
        result = train(args.pde, args.dimension, args.seed, device, pilot=args.pilot)
    else:
        result = evaluate(
            args.pde, args.dimension, args.seed, args.replicate, device,
            pilot=args.pilot, time_steps=args.time_steps,
        )
    print(json.dumps({key: result.get(key) for key in (
        "run_id", "status", "selected_actor", "rmse", "nonfinite_count"
    )}, sort_keys=True))


if __name__ == "__main__":
    main()

