"""Mature Terminal-Aware Tangent-Actor runner for HJB, OU/GMM and Allen--Cahn."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.original_pde_benchmark_comparison.common import (
    seed_everything, synchronize,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.models import (
    TerminalAwareActor,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.references import (
    AllenRadialConfig, AllenRadialReference, hjb_log_reference,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.tangent_solver import (
    ActorConfig, evaluate_policy, fit_actor, generate_labels,
)

from .benchmarks import AllenCahnBenchmark, HJBLogBenchmark, OUGMMBenchmark
from .prepare_suite import sha256_file, write_json
from .protocol import RESULTS_DIR, stable_seed


FAMILIES = ("hjb_log", "ou_gmm", "allen_cahn")
CHECKPOINT_DIR = RESULTS_DIR / "checkpoints"
CACHE_DIR = RESULTS_DIR / "cache"
PREDICTION_DIR = RESULTS_DIR / "predictions"
RAW_DIR = RESULTS_DIR / "raw"
INPUT_DIR = RESULTS_DIR / "inputs"
REFERENCE_DIR = RESULTS_DIR / "references"

MATURE = ActorConfig(
    training_points=8192, validation_points=512, paths_per_label=128,
    label_time_steps=32, optimizer_steps=5000, width=384, depth=4,
    activation="silu", learning_rate=1.5e-4, batch_size=512,
    parent_batch=8, path_chunk=32, checkpoint_interval=250,
    gradient_clip=10.0, training_distribution="function", label_estimator="mcv",
)
SMOKE = ActorConfig(
    training_points=32, validation_points=16, paths_per_label=8,
    label_time_steps=4, optimizer_steps=4, width=16, depth=2,
    activation="silu", learning_rate=2e-4, batch_size=16,
    parent_batch=4, path_chunk=4, checkpoint_interval=2,
    gradient_clip=10.0, training_distribution="function", label_estimator="mcv",
)


def _dirs() -> None:
    for path in (CHECKPOINT_DIR, CACHE_DIR, PREDICTION_DIR, RAW_DIR, REFERENCE_DIR):
        path.mkdir(parents=True, exist_ok=True)


def make_benchmark(family: str, d: int):
    if family == "hjb_log":
        return HJBLogBenchmark(d)
    if family == "ou_gmm":
        return OUGMMBenchmark.from_npz(INPUT_DIR / f"ou_gmm_d{d}.npz")
    if family == "allen_cahn":
        return AllenCahnBenchmark(d)
    raise KeyError(family)


def _purpose(family: str) -> str:
    return {"hjb_log": "function", "ou_gmm": "official", "allen_cahn": "native"}[family]


def _run_id(family: str, d: int, seed: int, tag: str) -> str:
    return f"{family}_d{d}_{tag}_s{seed}"


def _array_hash(*arrays: np.ndarray) -> str:
    digest = hashlib.sha256()
    for item in arrays:
        value = np.ascontiguousarray(item)
        digest.update(str(value.dtype).encode()); digest.update(str(value.shape).encode())
        digest.update(value.tobytes())
    return digest.hexdigest()


def _checkpoint(run_id: str, iteration: int) -> Path:
    return CHECKPOINT_DIR / f"{run_id}_a{iteration}.pt"


def _save_actor(path: Path, actor: TerminalAwareActor, benchmark, config, iteration: int) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    torch.save({
        "state_dict": {key: value.detach().cpu() for key, value in actor.state_dict().items()},
        "benchmark": benchmark.name, "dimension": benchmark.d, "iteration": iteration,
        "width": config.width, "depth": config.depth, "activation": config.activation,
        "config": asdict(config), "oracle_projection": False,
    }, temporary)
    temporary.replace(path)


def _load_actor(path: Path, benchmark, device: torch.device) -> TerminalAwareActor:
    payload = torch.load(path, map_location=device, weights_only=False)
    actor = TerminalAwareActor(
        benchmark, width=int(payload["width"]), depth=int(payload["depth"]),
        activation=str(payload["activation"]),
    ).to(device)
    actor.load_state_dict(payload["state_dict"]); actor.eval()
    return actor


def _initial_actor(run_id: str, benchmark, config, device: torch.device) -> TerminalAwareActor:
    seed_everything(stable_seed(run_id, "a0"))
    actor = TerminalAwareActor(
        benchmark, width=config.width, depth=config.depth, activation=config.activation,
    ).to(device)
    actor.eval()
    return actor


def _states(family: str, benchmark, count: int, *, split: str, seed: int):
    return benchmark.sample_states(
        count, seed=stable_seed("icml_standard_states", family, benchmark.d, split, seed),
        purpose=_purpose(family),
    )


def _labels(run_id: str, family: str, benchmark, source, iteration: int, config, seed: int, device):
    cache = CACHE_DIR / f"{run_id}_a{iteration}_labels.npz"
    metadata_path = CACHE_DIR / f"{run_id}_a{iteration}_labels.json"
    if cache.exists() and metadata_path.exists():
        with np.load(cache) as data:
            arrays = {key: data[key] for key in data.files}
        return arrays, json.loads(metadata_path.read_text())
    train_t, train_x = _states(
        family, benchmark, config.training_points, split=f"train_a{iteration}", seed=seed,
    )
    valid_t, valid_x = _states(
        family, benchmark, config.validation_points, split="fixed_validation", seed=seed,
    )
    kwargs = dict(
        paths=config.paths_per_label, time_steps=config.label_time_steps,
        device=device, parent_batch=config.parent_batch, path_chunk=config.path_chunk,
        estimator=config.label_estimator,
    )
    train_y, train_meta = generate_labels(
        benchmark, source, train_t, train_x,
        seed=stable_seed(run_id, iteration, "train_labels"), **kwargs,
    )
    valid_y, valid_meta = generate_labels(
        benchmark, source, valid_t, valid_x,
        seed=stable_seed(run_id, iteration, "validation_labels"), **kwargs,
    )
    replicate_y, replicate_meta = generate_labels(
        benchmark, source, valid_t, valid_x,
        seed=stable_seed(run_id, iteration, "replicate_labels"), **kwargs,
    )
    arrays = {
        "train_t": train_t, "train_x": train_x, "train_y": train_y,
        "valid_t": valid_t, "valid_x": valid_x, "valid_y": valid_y,
    }
    np.savez_compressed(cache, **arrays)
    difference = valid_y - replicate_y
    metadata = {
        "iteration": iteration, "source": f"a{iteration-1}",
        "train": train_meta, "validation": valid_meta, "replication": replicate_meta,
        "label_replication_component_rmse": float(np.sqrt(np.mean(difference**2))),
        "label_replication_vector_rmse": float(np.sqrt(np.mean(np.sum(difference**2, axis=1))))
        if difference.ndim == 2 else None,
        "train_bank_sha256": _array_hash(train_t, train_x),
        "validation_bank_sha256": _array_hash(valid_t, valid_x),
        "ground_truth_used": False,
    }
    write_json(metadata_path, metadata)
    return arrays, metadata


def _predict(actor, benchmark, t: np.ndarray, x: np.ndarray, device: torch.device) -> np.ndarray:
    result = []
    actor.eval()
    with torch.no_grad():
        for begin in range(0, len(t), 1024):
            end = min(begin + 1024, len(t))
            r, q = actor.predict(
                torch.as_tensor(t[begin:end], dtype=torch.float32, device=device),
                torch.as_tensor(x[begin:end], dtype=torch.float32, device=device),
            )
            target = q if benchmark.actor_kind == "gradient" else r
            if target is None:
                raise RuntimeError("Actor returned the wrong head")
            result.append(target.double().cpu().numpy())
    return np.concatenate(result)


def _fixed_point(run_id: str, family: str, benchmark, actor_name: str, actor, seed: int, device, smoke: bool):
    states, paths, steps = (16, 8, 4) if smoke else (128, 256, 32)
    t, x = _states(family, benchmark, states, split="fixed_point", seed=seed)
    labels, metadata = [], []
    for replication in (0, 1):
        value, meta = generate_labels(
            benchmark, actor, t, x, paths=paths, time_steps=steps,
            seed=stable_seed(run_id, "fixed_point", actor_name, replication),
            device=device, parent_batch=4 if smoke else 8,
            path_chunk=4 if smoke else 64, estimator="mcv",
        )
        labels.append(value); metadata.append(meta)
    prediction = _predict(actor, benchmark, t, x, device)
    target = 0.5 * (labels[0] + labels[1])
    raw = float(np.mean((prediction - target) ** 2))
    noise = float(np.mean((labels[0] - labels[1]) ** 2) / 4.0)
    corrected = max(raw - noise, 0.0)
    scale = max(float(np.sqrt(np.mean(target**2))), 1e-3)
    row = {
        "actor": actor_name, "states": states, "paths_per_replication": paths,
        "time_steps": steps, "state_bank_sha256": _array_hash(t, x),
        "raw_component_rmse": math.sqrt(raw),
        "noise_corrected_component_rmse": math.sqrt(corrected),
        "normalized_total_inconsistency": corrected / scale**2,
        "minimum_weight_ess_fraction": min(item["weight_ess_fraction"] for item in metadata),
        "nonfinite_count": int(sum(item["nonfinite_count"] for item in metadata)),
        "ground_truth_used": False,
    }
    if target.ndim == 2:
        row["raw_vector_rmse"] = math.sqrt(benchmark.d * raw)
        row["noise_corrected_vector_rmse"] = math.sqrt(benchmark.d * corrected)
    return row


def train(family: str, d: int, seed: int, device: torch.device, *, smoke: bool = False):
    _dirs()
    tag = "smoke" if smoke else "mature"
    base = SMOKE if smoke else MATURE
    config = ActorConfig(**{
        **asdict(base), "training_distribution": _purpose(family),
        "label_estimator": "plain" if family == "allen_cahn" else "mcv",
    })
    benchmark = make_benchmark(family, d)
    run_id = _run_id(family, d, seed, tag)
    output = RAW_DIR / f"training_{run_id}.json"
    actor = _initial_actor(run_id, benchmark, config, device)
    phases: list[dict[str, Any]] = []
    if output.exists():
        existing = json.loads(output.read_text())
        if existing["configuration"] != asdict(config):
            raise RuntimeError("existing configuration differs")
        phases = existing.get("phases", [])
        if phases:
            actor = _load_actor(_checkpoint(run_id, len(phases)), benchmark, device)
        if existing.get("status") == "complete":
            return existing
    started = time.perf_counter()
    for iteration in range(len(phases) + 1, 5):
        phase_started = time.perf_counter()
        arrays, label_meta = _labels(
            run_id, family, benchmark, actor, iteration, config, seed, device,
        )
        fitted, fit_meta, trace = fit_actor(
            benchmark, config, actor,
            arrays["train_t"], arrays["train_x"], arrays["train_y"],
            arrays["valid_t"], arrays["valid_x"], arrays["valid_y"],
            seed=stable_seed(run_id, iteration, "fit"), device=device,
        )
        checkpoint = _checkpoint(run_id, iteration)
        _save_actor(checkpoint, fitted, benchmark, config, iteration)
        phases.append({
            "iteration": iteration, "checkpoint": str(checkpoint),
            "phase_seconds": time.perf_counter() - phase_started,
            "labels": label_meta, "training": fit_meta, "trace": trace,
        })
        actor = fitted
        write_json(output, {
            "run_id": run_id, "method": "terminal_aware_tangent_actor_mcv",
            "family": family, "dimension": d, "seed": seed,
            "configuration": asdict(config), "phases": phases,
            "status": "training_in_progress", "full_ambient_input": True,
            "oracle_projection": False,
            "ground_truth_used_for_training_or_selection": False,
        })
        print(
            f"{run_id} a{iteration} label={label_meta['train']['wall_seconds']:.1f}s "
            f"fit={fit_meta['wall_seconds']:.1f}s val={fit_meta['selected_validation_normalized_mse']:.6g}",
            flush=True,
        )
    actors = [_initial_actor(run_id, benchmark, config, device)] + [
        _load_actor(_checkpoint(run_id, iteration), benchmark, device)
        for iteration in range(1, 5)
    ]
    fixed_points = [
        _fixed_point(run_id, family, benchmark, f"a{iteration}", candidate, seed, device, smoke)
        for iteration, candidate in enumerate(actors)
    ]
    eligible = [item for item in fixed_points if item["nonfinite_count"] == 0]
    selected = min(eligible, key=lambda item: item["normalized_total_inconsistency"])
    index = int(selected["actor"][1:])
    record = {
        "run_id": run_id, "method": "terminal_aware_tangent_actor_mcv",
        "family": family, "dimension": d, "seed": seed,
        "configuration": asdict(config), "phases": phases,
        "selection": {
            "criterion": "minimum independent noise-corrected fixed-point inconsistency",
            "rows": fixed_points, "selected_actor": selected["actor"],
            "selected_checkpoint": None if index == 0 else str(_checkpoint(run_id, index)),
            "ground_truth_used": False,
        },
        "full_ambient_input": True, "oracle_projection": False,
        "ground_truth_used_for_training_or_selection": False,
        "invocation_seconds": time.perf_counter() - started, "status": "complete",
    }
    write_json(output, record)
    return record


def _test_bank(family: str, d: int):
    audit = json.loads((RESULTS_DIR / "PREFLIGHT_AUDIT.json").read_text())
    row = audit["test_banks"][f"{family}_d{d}"]
    path = Path(__file__).resolve().parents[3] / row["path"]
    with np.load(path) as data:
        if "test_indices" in data:
            index = data["test_indices"]
            t, x = data["t"][index], data["x"][index]
        else:
            t, x = data["t"], data["x"]
    return np.asarray(t), np.asarray(x), path, row["sha256"]


def _allen_reference(d: int, t: np.ndarray, x: np.ndarray, *, smoke: bool):
    if smoke:
        config = AllenRadialConfig(s_max=100.0, ds=1.0, dt=0.005, dimension=d, horizon=0.3)
        value = AllenRadialReference(config).evaluate(t, np.square(x).sum(1))
        return value, {"smoke": True, "config": asdict(config)}
    bank_hash = _array_hash(t, x)
    path = REFERENCE_DIR / f"allen_cahn_d{d}_{bank_hash[:12]}.npz"
    audit_path = path.with_suffix(".json")
    if path.exists() and audit_path.exists():
        with np.load(path) as data:
            return data["extrapolated"], json.loads(audit_path.read_text())
    radius2 = np.square(x).sum(1)
    # Keep the far boundary beyond all frozen queries.  The two levels expose
    # the spatial/time floor rather than presenting a neural reference as exact.
    s_max = max(100.0, math.ceil((float(radius2.max()) + 20.0) / 10.0) * 10.0)
    coarse_cfg = AllenRadialConfig(
        s_max=s_max, ds=0.1, dt=0.0001, dimension=d, horizon=0.3,
    )
    fine_cfg = AllenRadialConfig(
        s_max=s_max, ds=0.05, dt=0.00005, dimension=d, horizon=0.3,
    )
    coarse = AllenRadialReference(coarse_cfg).evaluate(t, radius2)
    fine = AllenRadialReference(fine_cfg).evaluate(t, radius2)
    extrapolated = 2.0 * fine - coarse
    np.savez_compressed(path, t=t, radius2=radius2, coarse=coarse, fine=fine, extrapolated=extrapolated)
    audit = {
        "reference": "independent radial PDE Richardson extrapolation",
        "coarse": asdict(coarse_cfg), "fine": asdict(fine_cfg),
        "test_bank_array_sha256": bank_hash,
        "maximum_coarse_fine_difference": float(np.max(np.abs(fine - coarse))),
        "rmse_coarse_fine_difference": float(np.sqrt(np.mean((fine - coarse) ** 2))),
        "conservative_uncertainty": float(np.max(np.abs(extrapolated - fine))),
        "reference_file": str(path),
    }
    write_json(audit_path, audit)
    return extrapolated, audit


def _truth(family: str, benchmark, t: np.ndarray, x: np.ndarray, *, smoke: bool):
    if family == "hjb_log":
        return (*hjb_log_reference(t, x, order=64 if smoke else 128, return_gradient=True), {
            "kind": "independent Hopf-Cole Gauss-Laguerre quadrature", "order": 64 if smoke else 128,
        })
    if family == "ou_gmm":
        tt = torch.as_tensor(t, dtype=torch.float64)
        xx = torch.as_tensor(x, dtype=torch.float64)
        value = benchmark.exact_value(tt, xx).numpy()
        q = benchmark.exact_gradient(tt, xx).numpy()
        return value, q, {"kind": "exact evolving Gaussian-mixture density"}
    value, metadata = _allen_reference(benchmark.d, t, x, smoke=smoke)
    return value, None, metadata


def evaluate(family: str, d: int, seed: int, device: torch.device, *, smoke: bool = False,
             paths: int = 1024, steps: tuple[int, ...] = (64, 128)):
    tag = "smoke" if smoke else "mature"
    run_id = _run_id(family, d, seed, tag)
    training = json.loads((RAW_DIR / f"training_{run_id}.json").read_text())
    benchmark = make_benchmark(family, d)
    selected = training["selection"]["selected_actor"]
    if selected == "a0":
        actor = _initial_actor(run_id, benchmark, ActorConfig(**training["configuration"]), device)
    else:
        actor = _load_actor(Path(training["selection"]["selected_checkpoint"]), benchmark, device)
    t, x, bank_path, bank_hash = _test_bank(family, d)
    if smoke:
        t, x, paths, steps = t[:16], x[:16], 8, (4,)
    truth, truth_q, reference = _truth(family, benchmark, t, x, smoke=smoke)
    actor_prediction = _predict(actor, benchmark, t, x, device)
    rows, native_rows = [], []
    predictions = {"t": t, "x": x, "truth": truth, "actor": actor_prediction}
    native_t = np.zeros(1, dtype=np.float64)
    native_x = np.zeros((1, d), dtype=np.float64)
    native_truth = float(_truth(
        family, benchmark, native_t, native_x, smoke=smoke,
    )[0][0])
    for time_steps in steps:
        output, meta = evaluate_policy(
            benchmark, actor, t, x, paths=paths, time_steps=time_steps,
            seed=stable_seed(run_id, "postfreeze", paths, time_steps), device=device,
            parent_batch=4 if smoke else 8, path_chunk=4 if smoke else min(64, paths),
        )
        estimators = ("plain", "mcv") if benchmark.actor_kind == "gradient" else ("plain",)
        for estimator in estimators:
            prediction = output[estimator]
            error = prediction - truth
            variance = output[f"{estimator}_pair_variance"]
            rows.append({
                "steps": time_steps, "paths": paths, "estimator": estimator,
                "value_rmse": float(np.sqrt(np.mean(error**2))),
                "relative_l2": float(np.linalg.norm(error) / np.linalg.norm(truth)),
                "bias": float(np.mean(error)), "mae": float(np.mean(np.abs(error))),
                "mc_rmse_contribution": float(np.sqrt(np.mean(variance / (paths // 2)))),
                "variance_reduction": meta["variance_reduction"],
                "weight_ess_fraction": meta["weight_ess_fraction"],
                "inference_seconds": meta["wall_seconds"],
                "peak_gpu_bytes": meta["peak_gpu_bytes"],
                "nonfinite_count": meta["nonfinite_count"],
            })
            predictions[f"{estimator}_n{time_steps}"] = prediction
            predictions[f"{estimator}_pair_variance_n{time_steps}"] = variance
        native_output, native_meta = evaluate_policy(
            benchmark, actor, native_t, native_x, paths=paths,
            time_steps=time_steps,
            seed=stable_seed(run_id, "native_postfreeze", paths, time_steps),
            device=device, parent_batch=1,
            path_chunk=4 if smoke else min(64, paths),
        )
        native_estimators = (
            ("plain", "mcv") if benchmark.actor_kind == "gradient" else ("plain",)
        )
        for estimator in native_estimators:
            native_prediction = float(native_output[estimator][0])
            native_pair_variance = float(native_output[f"{estimator}_pair_variance"][0])
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
    result: dict[str, Any] = {
        "run_id": run_id, "method": "terminal_aware_tangent_actor_mcv",
        "family": family, "dimension": d, "seed": seed,
        "selected_actor": selected, "rows": rows, "native_rows": native_rows,
        "test_bank": str(bank_path), "test_bank_sha256": bank_hash,
        "states": len(t), "reference": reference,
        "prediction_path": str(prediction_path),
        "ground_truth_loaded_only_after_checkpoint_selection": True,
        "status": "complete",
    }
    difference = actor_prediction - (truth_q if benchmark.actor_kind == "gradient" else truth)
    if benchmark.actor_kind == "gradient":
        result.update({
            "gradient_actor_component_rmse": float(np.sqrt(np.mean(difference**2))),
            "gradient_actor_vector_rmse": float(np.sqrt(np.mean(np.sum(difference**2, axis=1)))),
        })
    else:
        result.update({
            "value_actor_rmse": float(np.sqrt(np.mean(difference**2))),
            "value_actor_relative_l2": float(np.linalg.norm(difference) / np.linalg.norm(truth)),
        })
    write_json(RAW_DIR / f"evaluation_{run_id}.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--family", choices=FAMILIES, required=True)
    parser.add_argument("--dimension", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--train-only", action="store_true")
    parser.add_argument("--evaluate-only", action="store_true")
    parser.add_argument("--evaluation-paths", type=int, default=1024)
    parser.add_argument("--evaluation-steps", default="64,128")
    args = parser.parse_args()
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA allocation required")
    if not args.evaluate_only:
        training = train(args.family, args.dimension, args.seed, device, smoke=args.smoke)
        print(json.dumps({"training": training["status"], "selected": training["selection"]["selected_actor"]}), flush=True)
    if not args.train_only:
        result = evaluate(
            args.family, args.dimension, args.seed, device, smoke=args.smoke,
            paths=args.evaluation_paths,
            steps=tuple(int(value) for value in args.evaluation_steps.split(",")),
        )
        print(json.dumps({"evaluation": result["status"], "rows": result["rows"]}), flush=True)


if __name__ == "__main__":
    main()
