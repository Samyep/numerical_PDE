"""Mature Terminal-Aware Tangent-Actor runner for PDEs 3--5."""

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

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.core import (
    make_targets,
)
from invariant_region_mlp.experiments.convex_concave_semilinear_extended.drivers import (
    ConcaveSemilinearDriver, JointConcaveDriver, JointConvexOriginalDriver,
    SignTransformedDriver,
)
from invariant_region_mlp.experiments.convex_concave_semilinear_extended.reference2d import (
    PeriodicReference2D, validate_periodic_reference_2d,
)
from invariant_region_mlp.experiments.nonconvex_tangent_actor.drivers import (
    SineQuadraticDriver,
)
from invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.solver import (
    evaluate_policy,
)
from invariant_region_mlp.experiments.tangent_actor_mcv_training.labels import (
    mcv_tangent_labels,
)
from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis.models import (
    TerminalAwareActor,
)
from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis.training import (
    FitConfig, train_actor,
)

from .prepare_suite import sha256_file, write_json
from .protocol import RESULTS_DIR, stable_seed


CHECKPOINT_DIR = RESULTS_DIR / "checkpoints"
CACHE_DIR = RESULTS_DIR / "cache"
PREDICTION_DIR = RESULTS_DIR / "predictions"
REFERENCE_DIR = RESULTS_DIR / "references"
RAW_DIR = RESULTS_DIR / "raw"
INPUT_DIR = RESULTS_DIR / "inputs"


MATURE = FitConfig(
    training_points=8192, validation_points=512, paths_per_label=128,
    label_time_steps=32, optimizer_steps=5000, width=384, depth=4,
    learning_rate=1.5e-4, batch_size=512, parent_batch=8, path_chunk=32,
    checkpoint_interval=250, gradient_clip=10.0,
)

SMOKE = FitConfig(
    training_points=32, validation_points=16, paths_per_label=8,
    label_time_steps=4, optimizer_steps=4, width=16, depth=2,
    learning_rate=2e-4, batch_size=16, parent_batch=4, path_chunk=4,
    checkpoint_interval=2, gradient_clip=10.0,
)


def _dirs() -> None:
    for path in (CHECKPOINT_DIR, CACHE_DIR, PREDICTION_DIR, REFERENCE_DIR, RAW_DIR):
        path.mkdir(parents=True, exist_ok=True)


def make_driver(family: str, d: int) -> ConcaveSemilinearDriver:
    if family == "joint_concave":
        return JointConcaveDriver(d=d, k=2, T=0.5)
    if family == "joint_convex":
        return SignTransformedDriver(JointConvexOriginalDriver(d=d, k=2, T=0.5))
    if family == "nonconvex_coupled":
        return SineQuadraticDriver(d=d, T=0.5, lam=0.2, eta=1.0)
    raise KeyError(family)


def _run_id(family: str, d: int, seed: int, tag: str) -> str:
    return f"{family}_d{d}_k2_{tag}_s{seed}"


def _state_bank(
    family: str, d: int, split: str, count: int, *, seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(stable_seed("icml_joint_state_bank", family, d, split, count, seed))
    t = rng.uniform(0.0, np.nextafter(0.5, 0.0), count).astype(np.float64)
    x = rng.standard_normal((count, d)).astype(np.float64)
    return t, x


def _array_hash(*arrays: np.ndarray) -> str:
    import hashlib
    digest = hashlib.sha256()
    for array in arrays:
        value = np.ascontiguousarray(array)
        digest.update(str(value.dtype).encode())
        digest.update(str(value.shape).encode())
        digest.update(value.tobytes())
    return digest.hexdigest()


def _checkpoint_path(run_id: str, iteration: int) -> Path:
    return CHECKPOINT_DIR / f"{run_id}_a{iteration}.pt"


def _save_actor(
    path: Path, actor: TerminalAwareActor, *, driver: ConcaveSemilinearDriver,
    config: FitConfig, iteration: int,
) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    torch.save({
        "state_dict": {key: value.detach().cpu() for key, value in actor.state_dict().items()},
        "driver": driver.metadata(), "config": asdict(config), "iteration": iteration,
        "actor_kind": "terminal_aware", "oracle_projection": False,
    }, temporary)
    temporary.replace(path)


def _load_actor(
    path: Path, driver: ConcaveSemilinearDriver, config: FitConfig,
    device: torch.device,
) -> TerminalAwareActor:
    payload = torch.load(path, map_location=device, weights_only=False)
    actor = TerminalAwareActor(driver, width=config.width, depth=config.depth).to(device)
    actor.load_state_dict(payload["state_dict"])
    actor.eval()
    return actor


def _initial_actor(
    driver: ConcaveSemilinearDriver, config: FitConfig, *, seed: int,
    device: torch.device,
) -> TerminalAwareActor:
    torch.manual_seed(seed)
    if device.type == "cuda":
        torch.cuda.manual_seed_all(seed)
    actor = TerminalAwareActor(driver, width=config.width, depth=config.depth).to(device)
    actor.eval()
    return actor


def _labels(
    run_id: str, driver: ConcaveSemilinearDriver, source: TerminalAwareActor,
    iteration: int, config: FitConfig, *, training_seed: int,
    device: torch.device,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    cache = CACHE_DIR / f"{run_id}_a{iteration}_labels.npz"
    meta_path = CACHE_DIR / f"{run_id}_a{iteration}_labels.json"
    if cache.exists() and meta_path.exists():
        with np.load(cache) as data:
            arrays = {key: data[key] for key in data.files}
        return arrays, json.loads(meta_path.read_text(encoding="utf-8"))
    train_t, train_x = _state_bank(
        driver.family, driver.d, f"train_a{iteration}", config.training_points,
        seed=training_seed,
    )
    valid_t, valid_x = _state_bank(
        driver.family, driver.d, "fixed_label_validation", config.validation_points,
        seed=training_seed,
    )
    kwargs = dict(
        paths=config.paths_per_label, time_steps=config.label_time_steps,
        device=device, parent_batch=config.parent_batch, path_chunk=config.path_chunk,
    )
    train_value, train_q, train_meta = mcv_tangent_labels(
        driver, source, train_t, train_x,
        seed=stable_seed(run_id, iteration, "train_labels"), **kwargs,
    )
    valid_value, valid_q, valid_meta = mcv_tangent_labels(
        driver, source, valid_t, valid_x,
        seed=stable_seed(run_id, iteration, "valid_labels"), **kwargs,
    )
    replicate_value, replicate_q, replicate_meta = mcv_tangent_labels(
        driver, source, valid_t, valid_x,
        seed=stable_seed(run_id, iteration, "replicate_labels"), **kwargs,
    )
    if train_q is None or valid_q is None or replicate_q is None:
        raise RuntimeError("joint benchmark requires gradient labels")
    train_r, train_q, train_projection = make_targets(driver, train_value, train_q)
    valid_r, valid_q, valid_projection = make_targets(driver, valid_value, valid_q)
    arrays = {
        "train_t": train_t, "train_x": train_x, "train_r": train_r, "train_q": train_q,
        "valid_t": valid_t, "valid_x": valid_x, "valid_r": valid_r, "valid_q": valid_q,
    }
    np.savez_compressed(cache, **arrays)
    metadata = {
        "iteration": iteration, "source": f"a{iteration-1}",
        "train": train_meta, "validation": valid_meta, "replication": replicate_meta,
        "train_projection": train_projection, "validation_projection": valid_projection,
        "value_label_replication_rmse": float(
            np.sqrt(np.mean(np.square(valid_value - replicate_value)))
        ),
        "gradient_label_replication_component_rmse": float(
            np.sqrt(np.mean(np.square(valid_q - replicate_q)))
        ),
        "gradient_label_replication_vector_rmse": float(
            np.sqrt(np.mean(np.square(valid_q - replicate_q).sum(axis=1)))
        ),
        "train_bank_sha256": _array_hash(train_t, train_x),
        "validation_bank_sha256": _array_hash(valid_t, valid_x),
        "ground_truth_used": False, "oracle_projection_used": False,
    }
    write_json(meta_path, metadata)
    return arrays, metadata


def _predict(
    driver: ConcaveSemilinearDriver, actor: TerminalAwareActor,
    t: np.ndarray, x: np.ndarray, device: torch.device,
) -> tuple[np.ndarray, np.ndarray]:
    r_parts, q_parts = [], []
    actor.eval()
    with torch.no_grad():
        for begin in range(0, len(t), 1024):
            end = min(begin + 1024, len(t))
            r, q = actor.predict(
                torch.as_tensor(t[begin:end], dtype=torch.float32, device=device),
                torch.as_tensor(x[begin:end], dtype=torch.float32, device=device),
                driver.T,
            )
            if q is None:
                raise RuntimeError("joint Actor omitted q")
            r_parts.append(r.double().cpu().numpy())
            q_parts.append(q.double().cpu().numpy())
    return np.concatenate(r_parts), np.concatenate(q_parts)


def _fixed_point(
    run_id: str, actor_name: str, driver: ConcaveSemilinearDriver,
    actor: TerminalAwareActor, *, training_seed: int, device: torch.device,
    smoke: bool,
) -> dict[str, Any]:
    states, paths, steps = (16, 8, 4) if smoke else (128, 256, 32)
    t, x = _state_bank(driver.family, driver.d, "fixed_point", states, seed=training_seed)
    labels = []
    metadata = []
    for replicate in (0, 1):
        value, q, meta = mcv_tangent_labels(
            driver, actor, t, x, paths=paths, time_steps=steps,
            seed=stable_seed(run_id, "fixed_point", actor_name, replicate),
            device=device, parent_batch=4 if smoke else 8,
            path_chunk=4 if smoke else 64,
        )
        if q is None:
            raise RuntimeError("missing q fixed-point label")
        labels.append((value, q)); metadata.append(meta)
    actor_r, actor_q = _predict(driver, actor, t, x, device)
    target_r = 0.5 * (labels[0][0] + labels[1][0])
    target_q = 0.5 * (labels[0][1] + labels[1][1])
    raw_r = float(np.mean(np.square(actor_r - target_r)))
    raw_q = float(np.mean(np.square(actor_q - target_q)))
    noise_r = float(np.mean(np.square(labels[0][0] - labels[1][0])) / 4.0)
    noise_q = float(np.mean(np.square(labels[0][1] - labels[1][1])) / 4.0)
    corrected_r = max(raw_r - noise_r, 0.0)
    corrected_q = max(raw_q - noise_q, 0.0)
    r_scale = max(float(np.std(target_r)), 1e-2)
    q_scale = max(float(np.sqrt(np.mean(np.square(target_q)))), 1e-3)
    return {
        "actor": actor_name, "states": states, "paths_per_replication": paths,
        "time_steps": steps, "state_bank_sha256": _array_hash(t, x),
        "value_raw_rmse": math.sqrt(raw_r),
        "value_noise_corrected_rmse": math.sqrt(corrected_r),
        "gradient_raw_component_rmse": math.sqrt(raw_q),
        "gradient_raw_vector_rmse": math.sqrt(driver.d * raw_q),
        "gradient_noise_corrected_component_rmse": math.sqrt(corrected_q),
        "gradient_noise_corrected_vector_rmse": math.sqrt(driver.d * corrected_q),
        "normalized_value_inconsistency": corrected_r / r_scale**2,
        "normalized_gradient_inconsistency": corrected_q / q_scale**2,
        "normalized_total_inconsistency": corrected_r / r_scale**2 + corrected_q / q_scale**2,
        "minimum_weight_ess_fraction": min(m["weight_ess_fraction"] for m in metadata),
        "nonfinite_count": int(sum(m["nonfinite_count"] for m in metadata)),
        "ground_truth_used": False,
    }


def train(
    family: str, d: int, seed: int, device: torch.device, *, smoke: bool = False,
) -> dict[str, Any]:
    _dirs()
    tag = "smoke" if smoke else "mature"
    config = SMOKE if smoke else MATURE
    run_id = _run_id(family, d, seed, tag)
    output = RAW_DIR / f"training_{run_id}.json"
    driver = make_driver(family, d)
    actor = _initial_actor(
        driver, config, seed=stable_seed(run_id, "a0"), device=device,
    )
    phases: list[dict[str, Any]] = []
    if output.exists():
        existing = json.loads(output.read_text(encoding="utf-8"))
        if existing["configuration"] != asdict(config):
            raise RuntimeError("existing run has a different configuration")
        phases = existing.get("phases", [])
        if phases:
            actor = _load_actor(_checkpoint_path(run_id, len(phases)), driver, config, device)
        if existing.get("status") == "complete":
            return existing
    started = time.perf_counter()
    for iteration in range(len(phases) + 1, 5):
        phase_start = time.perf_counter()
        arrays, label_meta = _labels(
            run_id, driver, actor, iteration, config,
            training_seed=seed, device=device,
        )
        fitted, fit_meta, trace = train_actor(
            driver, config, "terminal_aware",
            arrays["train_t"], arrays["train_x"], arrays["train_r"], arrays["train_q"],
            arrays["valid_t"], arrays["valid_x"], arrays["valid_r"], arrays["valid_q"],
            seed=stable_seed(run_id, iteration, "fit"), device=device, initial=actor,
        )
        checkpoint = _checkpoint_path(run_id, iteration)
        _save_actor(checkpoint, fitted, driver=driver, config=config, iteration=iteration)
        phases.append({
            "iteration": iteration, "checkpoint": str(checkpoint),
            "phase_seconds": time.perf_counter() - phase_start,
            "labels": label_meta, "training": fit_meta, "trace": trace,
        })
        actor = fitted
        write_json(output, {
            "run_id": run_id, "family": family, "dimension": d, "seed": seed,
            "configuration": asdict(config), "driver": driver.metadata(),
            "phases": phases, "status": "training_in_progress",
            "ground_truth_used_for_training_or_selection": False,
        })
        print(
            f"{run_id} a{iteration} labels={label_meta['train']['wall_seconds']:.1f}s "
            f"fit={fit_meta['wall_seconds']:.1f}s val={fit_meta['best_validation_normalized_loss']:.6g}",
            flush=True,
        )
    actors = [
        _initial_actor(driver, config, seed=stable_seed(run_id, "a0"), device=device)
    ] + [
        _load_actor(_checkpoint_path(run_id, iteration), driver, config, device)
        for iteration in range(1, 5)
    ]
    fixed_points = [
        _fixed_point(
            run_id, f"a{iteration}", driver, candidate,
            training_seed=seed, device=device, smoke=smoke,
        )
        for iteration, candidate in enumerate(actors)
    ]
    eligible = [row for row in fixed_points if row["nonfinite_count"] == 0]
    selected = min(eligible, key=lambda row: row["normalized_total_inconsistency"])
    selected_iteration = int(selected["actor"][1:])
    record = {
        "run_id": run_id, "method": "terminal_aware_tangent_actor_mcv",
        "family": family, "dimension": d, "seed": seed,
        "configuration": asdict(config), "driver": driver.metadata(), "phases": phases,
        "selection": {
            "criterion": "minimum independent noise-corrected fixed-point inconsistency",
            "rows": fixed_points, "selected_actor": selected["actor"],
            "selected_checkpoint": None if selected_iteration == 0 else str(_checkpoint_path(run_id, selected_iteration)),
            "ground_truth_used": False,
        },
        "full_ambient_input": True, "oracle_projection": False,
        "ground_truth_used_for_training_or_selection": False,
        "invocation_seconds": time.perf_counter() - started, "status": "complete",
    }
    write_json(output, record)
    return record


def _save_reference(path: Path, reference: PeriodicReference2D) -> None:
    np.savez_compressed(
        path, name=np.asarray(reference.name), times=reference.times,
        space=reference.space, values=reference.values,
        derivative_1=reference.derivative_1, derivative_2=reference.derivative_2,
        dt=np.asarray(reference.dt), store_dt=np.asarray(reference.store_dt),
        n=np.asarray(reference.n),
    )


def _load_reference(path: Path) -> PeriodicReference2D:
    with np.load(path) as data:
        return PeriodicReference2D(
            name=str(data["name"]), times=data["times"], space=data["space"],
            values=data["values"], derivative_1=data["derivative_1"],
            derivative_2=data["derivative_2"], dt=float(data["dt"]),
            store_dt=float(data["store_dt"]), n=int(data["n"]),
        )


def reference(family: str) -> tuple[PeriodicReference2D, dict[str, Any]]:
    _dirs()
    if family == "nonconvex_coupled":
        source = (
            RESULTS_DIR.parent / "tangent_actor_accuracy_frontier/references/"
            "sine_quadratic_coupled_d100_k2_lam0p2_eta1p0_refined.npz"
        )
        audit_path = (
            RESULTS_DIR.parent / "tangent_actor_accuracy_frontier/references/"
            "coupled_refinement_audit.json"
        )
        return _load_reference(source), {
            "source": str(source), "sha256": sha256_file(source),
            "audit": str(audit_path),
            "reference_uncertainty_value_rmse": json.loads(audit_path.read_text())["reference_uncertainty_value_rmse"],
        }
    path = REFERENCE_DIR / f"{family}_k2_n128.npz"
    audit_path = REFERENCE_DIR / f"{family}_k2_audit.json"
    if not path.exists() or not audit_path.exists():
        driver = make_driver(family, 10)
        value, audit = validate_periodic_reference_2d(driver)
        _save_reference(path, value)
        write_json(audit_path, audit)
    return _load_reference(path), {
        "source": str(path), "sha256": sha256_file(path),
        "audit": str(audit_path),
        "reference_uncertainty_value_rmse": json.loads(audit_path.read_text())["value_refinement_rmse"],
    }


def _test_bank(family: str, d: int) -> tuple[np.ndarray, np.ndarray, str]:
    if family == "nonconvex_coupled" and d == 100:
        path = (
            RESULTS_DIR.parent / "tangent_actor_accuracy_frontier/predictions/"
            "post_freeze_coupled_d100_k2_lam0p2_eta1_main_mcv_s0.npz"
        )
        with np.load(path) as data:
            return data["test_t"], data["test_x"], str(path)
    path = INPUT_DIR / f"test_bank_{family}_d{d}.npz"
    with np.load(path) as data:
        return data["t"], data["x"], str(path)


def evaluate(
    family: str, d: int, seed: int, device: torch.device, *, smoke: bool = False,
    paths: int = 1024, steps: tuple[int, ...] = (64, 128),
) -> dict[str, Any]:
    tag = "smoke" if smoke else "mature"
    run_id = _run_id(family, d, seed, tag)
    training = json.loads((RAW_DIR / f"training_{run_id}.json").read_text())
    config = FitConfig(**training["configuration"])
    driver = make_driver(family, d)
    selected = training["selection"]["selected_actor"]
    if selected == "a0":
        actor = _initial_actor(driver, config, seed=stable_seed(run_id, "a0"), device=device)
    else:
        actor = _load_actor(Path(training["selection"]["selected_checkpoint"]), driver, config, device)
    test_t, test_x, bank_path = _test_bank(family, d)
    if smoke:
        test_t, test_x, paths, steps = test_t[:16], test_x[:16], 8, (4,)
    numerical_reference, reference_meta = reference(family)
    coordinates = test_x @ driver.B[:, :2]
    truth_working = numerical_reference.value(test_t, coordinates)
    truth_gradient_reduced = numerical_reference.gradient(test_t, coordinates)
    truth_q_working = driver.sigma * (truth_gradient_reduced @ driver.B[:, :2].T)
    truth = driver.working_to_original(truth_working)
    truth_q = driver.working_to_original(truth_q_working)
    actor_r_working, actor_q_working = _predict(driver, actor, test_t, test_x, device)
    actor_r = driver.working_to_original(actor_r_working)
    actor_q = driver.working_to_original(actor_q_working)
    rows, native_rows = [], []
    predictions: dict[str, np.ndarray] = {
        "t": test_t, "x": test_x, "truth": truth, "truth_q": truth_q,
        "actor_value": actor_r, "actor_q": actor_q,
    }
    native_t = np.zeros(1, dtype=np.float64)
    native_x = np.zeros((1, d), dtype=np.float64)
    native_working = float(numerical_reference.value(
        native_t, np.zeros((1, 2), dtype=np.float64),
    )[0])
    native_truth = float(driver.working_to_original(
        np.asarray([native_working], dtype=np.float64),
    )[0])
    for time_steps in steps:
        output, meta = evaluate_policy(
            driver, actor, test_t, test_x, paths=paths, time_steps=time_steps,
            seed=stable_seed(run_id, "postfreeze", paths, time_steps),
            device=device, parent_batch=4 if smoke else 8,
            path_chunk=4 if smoke else min(64, paths),
        )
        for estimator in ("plain", "mcv"):
            prediction = driver.working_to_original(output[estimator])
            error = prediction - truth
            pair_variance = output[f"{estimator}_pair_variance"]
            rows.append({
                "steps": time_steps, "paths": paths, "estimator": estimator,
                "value_rmse": float(np.sqrt(np.mean(np.square(error)))),
                "relative_l2": float(np.linalg.norm(error) / np.linalg.norm(truth)),
                "bias": float(np.mean(error)), "mae": float(np.mean(np.abs(error))),
                "mc_rmse_contribution": float(np.sqrt(np.mean(pair_variance / (paths // 2)))),
                "variance_reduction": meta["variance_reduction"],
                "weight_ess_fraction": meta["weight_ess_fraction"],
                "inference_seconds": meta["wall_seconds"],
                "peak_gpu_bytes": meta["peak_gpu_bytes"],
                "nonfinite_count": meta["nonfinite_count"],
            })
            predictions[f"{estimator}_n{time_steps}"] = prediction
            predictions[f"{estimator}_pair_variance_n{time_steps}"] = pair_variance
        native_output, native_meta = evaluate_policy(
            driver, actor, native_t, native_x, paths=paths,
            time_steps=time_steps,
            seed=stable_seed(run_id, "native_postfreeze", paths, time_steps),
            device=device, parent_batch=1,
            path_chunk=4 if smoke else min(64, paths),
        )
        for estimator in ("plain", "mcv"):
            native_prediction = float(driver.working_to_original(
                native_output[estimator],
            )[0])
            native_pair_variance = float(
                native_output[f"{estimator}_pair_variance"][0]
            )
            native_rows.append({
                "steps": time_steps, "paths": paths, "estimator": estimator,
                "prediction": native_prediction, "reference": native_truth,
                "absolute_error": abs(native_prediction - native_truth),
                "mc_standard_error": math.sqrt(native_pair_variance / (paths // 2)),
                "variance_reduction": native_meta["variance_reduction"],
                "weight_ess_fraction": native_meta["weight_ess_fraction"],
                "inference_seconds": native_meta["wall_seconds"],
                "nonfinite_count": native_meta["nonfinite_count"],
            })
    prediction_path = PREDICTION_DIR / f"{run_id}.npz"
    np.savez_compressed(prediction_path, **predictions)
    value_error = actor_r - truth
    q_error = actor_q - truth_q
    result = {
        "run_id": run_id, "method": "terminal_aware_tangent_actor_mcv",
        "family": family, "dimension": d, "seed": seed,
        "selected_actor": selected, "rows": rows, "native_rows": native_rows,
        "value_actor_rmse": float(np.sqrt(np.mean(np.square(value_error)))),
        "value_actor_relative_l2": float(np.linalg.norm(value_error) / np.linalg.norm(truth)),
        "gradient_actor_component_rmse": float(np.sqrt(np.mean(np.square(q_error)))),
        "gradient_actor_vector_rmse": float(np.sqrt(np.mean(np.square(q_error).sum(axis=1)))),
        "test_bank": bank_path, "states": len(test_t), "reference": reference_meta,
        "prediction_path": str(prediction_path),
        "ground_truth_loaded_only_after_checkpoint_selection": True,
        "status": "complete",
    }
    write_json(RAW_DIR / f"evaluation_{run_id}.json", result)
    return result


def evaluate_iterations(
    family: str, d: int, seed: int, device: torch.device, *, smoke: bool = False,
    paths: int = 512, time_steps: int = 64,
) -> dict[str, Any]:
    """Post-freeze update diagnostic on one common bank and MC stream.

    This diagnostic is deliberately separate from checkpoint selection.  It
    opens the independent numerical reference only after all actors and the
    reference-free selected checkpoint have been frozen.
    """
    tag = "smoke" if smoke else "mature"
    run_id = _run_id(family, d, seed, tag)
    training_path = RAW_DIR / f"training_{run_id}.json"
    training = json.loads(training_path.read_text())
    if training.get("status") != "complete":
        raise RuntimeError(f"training is not complete: {training_path}")
    config = FitConfig(**training["configuration"])
    driver = make_driver(family, d)
    test_t, test_x, bank_path = _test_bank(family, d)
    if smoke:
        test_t, test_x, paths, time_steps = test_t[:16], test_x[:16], 8, 4
    numerical_reference, reference_meta = reference(family)
    coordinates = test_x @ driver.B[:, :2]
    truth_working = numerical_reference.value(test_t, coordinates)
    truth_gradient_reduced = numerical_reference.gradient(test_t, coordinates)
    truth_q_working = driver.sigma * (truth_gradient_reduced @ driver.B[:, :2].T)
    truth = driver.working_to_original(truth_working)
    truth_q = driver.working_to_original(truth_q_working)
    consistency = {
        row["actor"]: row for row in training["selection"]["rows"]
    }
    common_mc_seed = stable_seed(
        run_id, "all_iterations_postfreeze", paths, time_steps,
    )
    rows: list[dict[str, Any]] = []
    for iteration in range(5):
        actor_name = f"a{iteration}"
        checkpoint = None if iteration == 0 else _checkpoint_path(run_id, iteration)
        actor = (
            _initial_actor(driver, config, seed=stable_seed(run_id, "a0"), device=device)
            if iteration == 0
            else _load_actor(checkpoint, driver, config, device)
        )
        actor_r_working, actor_q_working = _predict(
            driver, actor, test_t, test_x, device,
        )
        actor_r = driver.working_to_original(actor_r_working)
        actor_q = driver.working_to_original(actor_q_working)
        value_error = actor_r - truth
        q_error = actor_q - truth_q
        output, meta = evaluate_policy(
            driver, actor, test_t, test_x, paths=paths, time_steps=time_steps,
            seed=common_mc_seed, device=device,
            parent_batch=4 if smoke else 8,
            path_chunk=4 if smoke else min(64, paths),
        )
        selection_row = consistency[actor_name]
        for estimator in ("plain", "mcv"):
            prediction = driver.working_to_original(output[estimator])
            error = prediction - truth
            pair_variance = output[f"{estimator}_pair_variance"]
            rows.append({
                "actor": actor_name, "iteration": iteration,
                "selected_by_reference_free_rule": (
                    actor_name == training["selection"]["selected_actor"]
                ),
                "checkpoint": None if checkpoint is None else str(checkpoint),
                "steps": time_steps, "paths": paths, "estimator": estimator,
                "value_rmse": float(np.sqrt(np.mean(np.square(error)))),
                "relative_l2": float(np.linalg.norm(error) / np.linalg.norm(truth)),
                "bias": float(np.mean(error)),
                "mc_rmse_contribution": float(
                    np.sqrt(np.mean(pair_variance / (paths // 2)))
                ),
                "variance_reduction": meta["variance_reduction"],
                "weight_ess_fraction": meta["weight_ess_fraction"],
                "inference_seconds": meta["wall_seconds"],
                "peak_gpu_bytes": meta["peak_gpu_bytes"],
                "nonfinite_count": meta["nonfinite_count"],
                "actor_value_rmse": float(np.sqrt(np.mean(np.square(value_error)))),
                "actor_value_relative_l2": float(
                    np.linalg.norm(value_error) / np.linalg.norm(truth)
                ),
                "actor_gradient_component_rmse": float(
                    np.sqrt(np.mean(np.square(q_error)))
                ),
                "actor_gradient_vector_rmse": float(
                    np.sqrt(np.mean(np.square(q_error).sum(axis=1)))
                ),
                "reference_free_total_inconsistency": selection_row[
                    "normalized_total_inconsistency"
                ],
                "reference_free_value_inconsistency": selection_row[
                    "normalized_value_inconsistency"
                ],
                "reference_free_gradient_inconsistency": selection_row[
                    "normalized_gradient_inconsistency"
                ],
            })
    result = {
        "run_id": run_id, "method": "terminal_aware_tangent_actor_mcv",
        "family": family, "dimension": d, "seed": seed,
        "diagnostic_role": "post_freeze_all_actor_iterations",
        "selected_actor": training["selection"]["selected_actor"],
        "common_mc_seed": common_mc_seed, "rows": rows,
        "test_bank": bank_path, "states": len(test_t), "reference": reference_meta,
        "ground_truth_loaded_only_after_all_checkpoints_and_selection_frozen": True,
        "status": "complete",
    }
    write_json(RAW_DIR / f"iteration_evaluation_{run_id}.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--family", choices=("joint_concave", "joint_convex", "nonconvex_coupled"), required=True)
    parser.add_argument("--dimension", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--train-only", action="store_true")
    parser.add_argument("--evaluate-only", action="store_true")
    parser.add_argument("--iteration-diagnostic-only", action="store_true")
    parser.add_argument("--evaluation-paths", type=int, default=1024)
    parser.add_argument("--evaluation-steps", default="64,128")
    parser.add_argument("--iteration-paths", type=int, default=512)
    parser.add_argument("--iteration-steps", type=int, default=64)
    args = parser.parse_args()
    if (
        not args.smoke and not args.evaluate_only
        and not args.iteration_diagnostic_only and args.seed != 0
    ):
        parser.error("single-seed benchmark protocol permits only seed 0 for new training")
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA allocation required")
    if not args.evaluate_only and not args.iteration_diagnostic_only:
        record = train(args.family, args.dimension, args.seed, device, smoke=args.smoke)
        print(json.dumps({"training": record["status"], "selected": record["selection"]["selected_actor"]}), flush=True)
    if not args.train_only and not args.iteration_diagnostic_only:
        result = evaluate(
            args.family, args.dimension, args.seed, device, smoke=args.smoke,
            paths=args.evaluation_paths,
            steps=tuple(int(value) for value in args.evaluation_steps.split(",")),
        )
        print(json.dumps({"evaluation": result["status"], "rows": result["rows"]}), flush=True)
    if args.iteration_diagnostic_only:
        result = evaluate_iterations(
            args.family, args.dimension, args.seed, device, smoke=args.smoke,
            paths=args.iteration_paths, time_steps=args.iteration_steps,
        )
        print(json.dumps({
            "iteration_evaluation": result["status"],
            "selected": result["selected_actor"], "rows": len(result["rows"]),
        }), flush=True)


if __name__ == "__main__":
    main()
