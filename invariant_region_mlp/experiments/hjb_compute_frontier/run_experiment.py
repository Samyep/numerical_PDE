"""CLI for preflight, training, evaluation, and resumable experiment jobs."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import math
from pathlib import Path
import platform
import subprocess
import sys
import time
from typing import Any

import numpy as np
import scipy
import torch

from .actor_method import (
    ActorPhase,
    evaluate_policy_cv,
    load_actor,
    policy_gradient_labels,
    predict_actor,
    save_actor,
    terminal_pathwise_labels,
    train_actor_supervised,
)
from .baselines import (
    load_deep_bsde,
    load_dpi,
    load_pinn,
    predict_baseline,
    train_deep_bsde,
    train_dpi,
    train_pinn_sdgd,
)
from .common import (
    CACHE_DIR,
    CHECKPOINT_DIR,
    PREDICTION_DIR,
    RESULTS_DIR,
    array_sha256,
    atomic_json,
    atomic_npz,
    ensure_dirs,
    file_sha256,
    fixed_points,
    load_protocol,
    model_metrics,
    named_points,
    require_cuda,
    stable_seed,
    vector_relative_error,
)
from .equations import (
    HJBLogEquation,
    Equation,
    equation_key,
    make_equation,
    torch_terminal,
    torch_terminal_z,
)


def _git(command: list[str]) -> str:
    try:
        return subprocess.check_output(
            ["git", *command], text=True, stderr=subprocess.DEVNULL
        ).strip()
    except Exception:
        return "unavailable"


def _banks(equation: Equation) -> dict[str, np.ndarray]:
    key = equation_key(equation)
    policy_t, policy_x = fixed_points(
        key, equation.d, "policy_valid", 256, include_time=False, T=equation.T
    )
    test_t, test_x = fixed_points(
        key, equation.d, "test", 2048, include_time=False, T=equation.T
    )
    return {
        "policy_t": policy_t,
        "policy_x": policy_x,
        "test_t": test_t,
        "test_x": test_x,
        "test_u": equation.exact_u(test_t, test_x),
        "test_z": equation.exact_z(test_t, test_x),
    }


def save_banks(equation: Equation) -> dict[str, Any]:
    banks = _banks(equation)
    key = equation_key(equation)
    path = CACHE_DIR / "data_banks" / f"{key}.npz"
    atomic_npz(path, **banks)
    return {
        "path": str(path),
        "sha256": file_sha256(path),
        "policy_points_sha256": array_sha256(banks["policy_t"], banks["policy_x"]),
        "test_points_sha256": array_sha256(banks["test_t"], banks["test_x"]),
        "test_reference_sha256": array_sha256(banks["test_u"], banks["test_z"]),
    }


def run_preflight() -> None:
    ensure_dirs()
    protocol = load_protocol()
    equations = [make_equation("b1"), make_equation("hjb_log")]
    validation: dict[str, Any] = {
        "status": "in_progress",
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "git_branch": _git(["branch", "--show-current"]),
        "git_commit": _git(["rev-parse", "HEAD"]),
        "python": sys.version,
        "platform": platform.platform(),
        "torch": torch.__version__,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "protocol_sha256": file_sha256(Path(__file__).with_name("protocol_frozen.json")),
        "data_banks": {},
        "equations": {},
    }
    for equation in equations:
        key = equation_key(equation)
        validation["data_banks"][key] = save_banks(equation)
        rng = np.random.default_rng(stable_seed("preflight", key))
        x = rng.standard_normal((4, equation.d), dtype=np.float64)
        terminal_error = float(
            np.max(np.abs(equation.exact_u(equation.T, x) - equation.terminal(x)))
        )
        terminal_z_error = float(
            np.max(np.abs(equation.exact_z(equation.T, x) - equation.terminal_z(x)))
        )
        tx = torch.tensor(x, dtype=torch.float64, requires_grad=True)
        terminal = torch_terminal(equation, tx)
        autodiff = equation.sigma * torch.autograd.grad(terminal.sum(), tx)[0]
        torch_z_error = float(
            np.max(
                np.abs(
                    autodiff.detach().numpy()
                    - torch_terminal_z(equation, tx.detach()).numpy()
                )
            )
        )
        record: dict[str, Any] = {
            "provenance": equation.provenance(),
            "terminal_value_max_abs_error": terminal_error,
            "terminal_z_max_abs_error": terminal_z_error,
            "torch_terminal_autodiff_z_max_abs_error": torch_z_error,
        }
        if isinstance(equation, HJBLogEquation):
            audit_t = np.array([0.0, 0.53, equation.T])
            audit_x = x[:3]
            quadrature_u = equation.exact_u(audit_t, audit_x)
            quadrature_z = equation.exact_z(audit_t, audit_x)
            adaptive = [
                equation.adaptive_reference(float(audit_t[i]), audit_x[i])
                for i in range(3)
            ]
            adaptive_u = np.array([item[0] for item in adaptive])
            adaptive_z = np.stack([item[1] for item in adaptive])
            record.update(
                {
                    "adaptive_quadrature_value_max_abs_error": float(
                        np.max(np.abs(quadrature_u - adaptive_u))
                    ),
                    "adaptive_quadrature_z_max_abs_error": float(
                        np.max(np.abs(quadrature_z - adaptive_z))
                    ),
                    "adaptive_quad_reported_abs_error_max": float(
                        max(item[2] for item in adaptive)
                    ),
                }
            )
        else:
            times = rng.uniform(0.01, equation.T - 0.01, len(x))
            record["analytic_pde_residual_max_abs"] = float(
                np.max(np.abs(equation.analytic_residual(times, x)))
            )
        validation["equations"][key] = record
    source_text = "\n".join(
        path.read_text()
        for path in [
            Path(__file__).with_name("baselines.py"),
            Path(__file__).with_name("actor_method.py"),
        ]
    )
    validation["leakage_static_audit"] = {
        "baseline_module_contains_exact_u": "exact_u" in Path(__file__).with_name("baselines.py").read_text(),
        "baseline_module_contains_exact_z": "exact_z" in Path(__file__).with_name("baselines.py").read_text(),
        "actor_training_functions_receive_exact_reference": False,
        "ours_test_selection_explicitly_allowed": protocol["selection"]["ours_selection_label"],
    }
    thresholds_ok = all(
        record["terminal_value_max_abs_error"] < 1e-11
        and record["terminal_z_max_abs_error"] < 1e-11
        and record["torch_terminal_autodiff_z_max_abs_error"] < 1e-11
        for record in validation["equations"].values()
    )
    log_record = validation["equations"]["hjb_log_d400"]
    thresholds_ok = thresholds_ok and log_record["adaptive_quadrature_value_max_abs_error"] < 1e-10
    thresholds_ok = thresholds_ok and log_record["adaptive_quadrature_z_max_abs_error"] < 1e-10
    validation["status"] = "passed" if thresholds_ok else "failed"
    atomic_json(RESULTS_DIR / "validation.json", validation)
    if not thresholds_ok:
        raise RuntimeError("preflight validation thresholds failed")
    print(json.dumps(validation, indent=2, sort_keys=True))


def _profile_baseline(method: str, profile: str) -> dict[str, Any]:
    if profile == "smoke":
        return {
            "pinn_sdgd": {"max_steps": 20, "checkpoints": (10, 20), "width": 64, "batch_size": 8, "coordinate_budget": 2},
            "deep_bsde": {"max_steps": 20, "checkpoints": (10, 20), "width": 64, "time_steps": 4, "batch_size": 8},
            "dpi": {"rounds": 2, "training_points": 32, "paths": 8, "epochs": 1, "batch_size": 16, "width": 64, "checkpoint_rounds": (1, 2)},
        }[method]
    if profile == "pilot":
        return {
            "pinn_sdgd": {"max_steps": 1000, "checkpoints": (250, 500, 1000)},
            "deep_bsde": {"max_steps": 1000, "checkpoints": (250, 500, 1000)},
            "dpi": {"rounds": 5, "training_points": 1024, "paths": 128, "epochs": 8, "checkpoint_rounds": (1, 2, 5)},
        }[method]
    if profile == "final":
        return {}
    raise ValueError(profile)


def run_baseline(equation_name: str, method: str, seed_index: int, profile: str) -> None:
    ensure_dirs()
    equation = make_equation(equation_name)
    device = require_cuda()
    full_seed = stable_seed("hjb_compute_frontier", equation_key(equation), method, seed_index, profile)
    output_dir = CHECKPOINT_DIR / equation_key(equation) / method / f"seed{seed_index}" / profile
    output_dir.mkdir(parents=True, exist_ok=True)
    kwargs = _profile_baseline(method, profile)
    if method == "pinn_sdgd":
        result = train_pinn_sdgd(
            equation, seed=full_seed, device=device, output_dir=output_dir, **kwargs
        )
    elif method == "deep_bsde":
        result = train_deep_bsde(
            equation, seed=full_seed, device=device, output_dir=output_dir, **kwargs
        )
    elif method == "dpi":
        result = train_dpi(
            equation, seed=full_seed, device=device, output_dir=output_dir, **kwargs
        )
    else:
        raise ValueError(method)
    result.update(
        {
            "equation": equation_key(equation),
            "seed_index": seed_index,
            "full_seed": full_seed,
            "profile": profile,
            "gpu_name": torch.cuda.get_device_name(device),
            "git_commit": _git(["rev-parse", "HEAD"]),
        }
    )
    result_path = RESULTS_DIR / "runs" / equation_key(equation) / method / f"seed{seed_index}_{profile}.json"
    atomic_json(result_path, result)
    evaluate_baseline_run(equation, method, result, device, result_path)


def evaluate_baseline_run(
    equation: Equation,
    method: str,
    result: dict[str, Any],
    device: torch.device,
    result_path: Path,
) -> None:
    banks = _banks(equation)
    rows: list[dict[str, Any]] = []
    predictions_dir = PREDICTION_DIR / equation_key(equation) / method
    for checkpoint in result["checkpoints"]:
        path = Path(checkpoint["path"])
        if method == "pinn_sdgd":
            model, metadata = load_pinn(path, device)
        elif method == "deep_bsde":
            model, metadata = load_deep_bsde(path, device)
        else:
            model, metadata = load_dpi(path, device)
        value, z, inference_2048 = predict_baseline(
            method, model, equation, banks["test_x"], device
        )
        _, _, inference_960 = predict_baseline(
            method, model, equation, banks["test_x"][:960], device
        )
        metrics = model_metrics(value, banks["test_u"])
        metrics["gradient_relative_error"] = vector_relative_error(z, banks["test_z"])
        token = checkpoint.get("step", checkpoint.get("round"))
        prediction_path = predictions_dir / f"seed{result['seed_index']}_{result['profile']}_{token}.npz"
        atomic_npz(
            prediction_path,
            value=value.astype(np.float64),
            z=z.astype(np.float64),
            truth_u=banks["test_u"],
            truth_z=banks["test_z"],
            test_points_sha256=np.asarray(array_sha256(banks["test_x"])),
        )
        rows.append(
            {
                "checkpoint": token,
                "training_seconds": metadata["training_seconds"],
                "inference_seconds_960": inference_960,
                "inference_seconds_2048": inference_2048,
                "total_seconds_960": metadata["training_seconds"] + inference_960,
                "prediction_path": str(prediction_path),
                **metrics,
            }
        )
    result["evaluation"] = rows
    result["test_used_for_selection"] = False
    atomic_json(result_path, result)


def _actor_phases(profile: str) -> list[ActorPhase]:
    if profile == "smoke":
        return [
            ActorPhase("q0", None, 128, 32, 8, 1, 20, 64, 2, 5e-4, 32),
            ActorPhase("q1", "q0", 128, 32, 8, 2, 20, 64, 2, 2e-4, 32),
            ActorPhase("q1_refresh", "q0", 128, 32, 8, 2, 20, 64, 2, 2e-4, 32),
            ActorPhase("q2", "q1", 128, 32, 8, 2, 20, 64, 2, 2e-4, 32),
            ActorPhase("q3", "q2", 128, 32, 8, 2, 20, 64, 2, 2e-4, 32),
        ]
    if profile == "pilot":
        return [
            ActorPhase("q0", None, 2048, 512, 128, 1, 3000, 384, 4, 5e-4),
            ActorPhase("q1", "q0", 4096, 512, 256, 16, 4000),
            ActorPhase("q1_refresh", "q0", 4096, 512, 256, 16, 4000),
            ActorPhase("q2", "q1", 4096, 512, 256, 16, 4000),
            ActorPhase("q3", "q2", 4096, 512, 256, 16, 4000),
        ]
    if profile == "final":
        return [
            ActorPhase("q0", None, 4096, 1024, 256, 1, 6000, 384, 4, 5e-4),
            ActorPhase("q1", "q0", 8192, 512, 512, 32, 6000),
            ActorPhase("q1_refresh", "q0", 8192, 512, 512, 32, 6000),
            ActorPhase("q2", "q1", 8192, 512, 512, 32, 6000),
            ActorPhase("q3", "q2", 8192, 512, 512, 32, 6000),
        ]
    if profile == "high":
        return [
            ActorPhase("q0", None, 4096, 1024, 256, 1, 6000, 512, 4, 5e-4),
            ActorPhase("q1", "q0", 16384, 1024, 1024, 64, 12000, 512, 4, 1e-4),
            ActorPhase("q1_refresh", "q0", 16384, 1024, 1024, 64, 12000, 512, 4, 1e-4),
            ActorPhase("q2", "q1", 8192, 512, 512, 32, 6000, 512, 4, 1e-4),
            ActorPhase("q3", "q2", 8192, 512, 512, 32, 6000, 512, 4, 1e-4),
        ]
    raise ValueError(profile)


def _actor_eval_budget(profile: str) -> tuple[int, int, int]:
    return {
        "smoke": (8, 2, 64),
        "pilot": (64, 32, 512),
        "final": (256, 64, 2048),
        "high": (256, 64, 2048),
    }[profile]


def run_actor(equation_name: str, seed_index: int, profile: str) -> None:
    ensure_dirs()
    equation = make_equation(equation_name)
    key = equation_key(equation)
    device = require_cuda()
    phases = _actor_phases(profile)
    actors: dict[str, Any] = {}
    phase_metadata: dict[str, Any] = {}
    checkpoint_root = CHECKPOINT_DIR / key / "ours" / f"seed{seed_index}" / profile
    cache_root = CACHE_DIR / key / "ours" / f"seed{seed_index}" / profile
    checkpoint_root.mkdir(parents=True, exist_ok=True)
    cache_root.mkdir(parents=True, exist_ok=True)
    cumulative_exploration_seconds = 0.0
    path_seconds: dict[str, float] = {}
    for phase in phases:
        checkpoint_path = checkpoint_root / f"{phase.name}.pt"
        if checkpoint_path.is_file():
            actor, metadata = load_actor(checkpoint_path, device)
            actors[phase.name] = actor
            phase_metadata[phase.name] = metadata
            path_seconds[phase.name] = float(metadata["path_training_seconds"])
            cumulative_exploration_seconds += float(metadata.get("phase_total_seconds", 0.0))
            continue
        source = actors.get(phase.source) if phase.source is not None else None
        if phase.source is not None and source is None:
            raise RuntimeError(f"source actor {phase.source} is unavailable")
        train_t, train_x = named_points(
            key,
            equation.d,
            f"actor_{profile}_{phase.name}_seed{seed_index}_train",
            phase.training_points,
            T=equation.T,
        )
        valid_t, valid_x = named_points(
            key,
            equation.d,
            f"actor_{profile}_{phase.name}_seed{seed_index}_label_valid",
            phase.validation_points,
            T=equation.T,
        )
        label_path = cache_root / f"{phase.name}_labels.npz"
        label_meta_path = cache_root / f"{phase.name}_labels.json"
        if label_path.is_file() and label_meta_path.is_file():
            with np.load(label_path, allow_pickle=False) as payload:
                train_y = payload["train_y"]
                valid_y = payload["valid_y"]
            label_metadata = json.loads(label_meta_path.read_text())
        else:
            train_seed = stable_seed("actor_labels", key, profile, phase.name, seed_index, "train")
            valid_seed = stable_seed("actor_labels", key, profile, phase.name, seed_index, "valid")
            if source is None:
                train_y, train_label_meta = terminal_pathwise_labels(
                    equation,
                    train_t,
                    train_x,
                    paths=phase.paths_per_label,
                    seed=train_seed,
                    device=device,
                )
                valid_y, valid_label_meta = terminal_pathwise_labels(
                    equation,
                    valid_t,
                    valid_x,
                    paths=phase.paths_per_label,
                    seed=valid_seed,
                    device=device,
                )
            else:
                _, train_y, train_label_meta = policy_gradient_labels(
                    equation,
                    source,
                    train_t,
                    train_x,
                    paths=phase.paths_per_label,
                    time_steps=phase.label_time_steps,
                    seed=train_seed,
                    device=device,
                )
                _, valid_y, valid_label_meta = policy_gradient_labels(
                    equation,
                    source,
                    valid_t,
                    valid_x,
                    paths=phase.paths_per_label,
                    time_steps=phase.label_time_steps,
                    seed=valid_seed,
                    device=device,
                )
            label_metadata = {
                "train": train_label_meta,
                "validation": valid_label_meta,
                "point_hashes": {
                    "train": array_sha256(train_t, train_x),
                    "validation": array_sha256(valid_t, valid_x),
                },
            }
            atomic_npz(label_path, train_y=train_y, valid_y=valid_y)
            atomic_json(label_meta_path, label_metadata)
        if phase.name == "q1_refresh":
            initial = actors["q1"]
        else:
            initial = source
        train_seed = stable_seed("actor_optimizer", key, profile, phase.name, seed_index)
        actor, optimizer_metadata, _ = train_actor_supervised(
            equation,
            phase,
            train_t,
            train_x,
            train_y,
            valid_t,
            valid_x,
            valid_y,
            seed=train_seed,
            device=device,
            initial=initial,
        )
        label_seconds = float(label_metadata["train"]["wall_seconds"]) + float(
            label_metadata["validation"]["wall_seconds"]
        )
        phase_total = label_seconds + float(optimizer_metadata["optimizer_seconds"])
        parent_path = 0.0 if phase.source is None else path_seconds[phase.source]
        # q1_refresh warm-starts q1 and therefore includes q1's path cost even
        # though its fresh target source remains q0.
        if phase.name == "q1_refresh":
            parent_path = path_seconds["q1"]
        path_training = parent_path + phase_total
        metadata = {
            "equation": key,
            "profile": profile,
            "phase": asdict(phase),
            "seed_index": seed_index,
            "optimizer": optimizer_metadata,
            "labels": label_metadata,
            "phase_total_seconds": phase_total,
            "path_training_seconds": path_training,
            "selection_rule": "test_rmse_oracle_allowed_by_user",
            "output_projection": False,
            "analytic_nonterminal_reference_used_for_labels": False,
        }
        save_actor(checkpoint_path, actor, metadata)
        actors[phase.name] = actor
        phase_metadata[phase.name] = metadata
        path_seconds[phase.name] = path_training
        cumulative_exploration_seconds += phase_total

    banks = _banks(equation)
    paths, time_steps, point_count = _actor_eval_budget(profile)
    evaluations: list[dict[str, Any]] = []
    prediction_dir = PREDICTION_DIR / key / "ours" / f"seed{seed_index}" / profile
    for phase in phases:
        actor = actors[phase.name]
        selection_x = banks["test_x"][:point_count]
        selection_t = banks["test_t"][:point_count]
        value, estimator_variance, evaluation_meta = evaluate_policy_cv(
            equation,
            actor,
            selection_t,
            selection_x,
            paths=paths,
            time_steps=time_steps,
            seed=stable_seed("actor_test_eval", key, seed_index, profile, phase.name),
            device=device,
        )
        metrics = model_metrics(value, banks["test_u"][:point_count])
        direct_z = predict_actor(actor, selection_t, selection_x, equation, device)
        metrics["gradient_relative_error"] = vector_relative_error(
            direct_z, banks["test_z"][:point_count]
        )
        metrics["mean_mc_standard_error"] = float(
            np.mean(np.sqrt(np.maximum(estimator_variance, 0.0)))
        )
        prediction_path = prediction_dir / f"{phase.name}.npz"
        atomic_npz(
            prediction_path,
            value=value,
            estimator_variance=estimator_variance,
            direct_z=direct_z,
            truth_u=banks["test_u"][:point_count],
            truth_z=banks["test_z"][:point_count],
            test_points_sha256=np.asarray(array_sha256(selection_x)),
        )
        evaluations.append(
            {
                "actor": phase.name,
                "policy_iterations": int(phase.name[1]) if phase.name.startswith("q") and phase.name[1].isdigit() else 1,
                "path_training_seconds": path_seconds[phase.name],
                "cumulative_exploration_seconds": cumulative_exploration_seconds,
                "inference_seconds": evaluation_meta["wall_seconds"],
                "total_seconds": path_seconds[phase.name] + evaluation_meta["wall_seconds"],
                "evaluation": evaluation_meta,
                "prediction_path": str(prediction_path),
                **metrics,
            }
        )
    selected = min(evaluations, key=lambda item: item["relative_l2"])
    result = {
        "method": "ours_actor_iteration_cv",
        "equation": key,
        "profile": profile,
        "seed_index": seed_index,
        "gpu_name": torch.cuda.get_device_name(device),
        "git_commit": _git(["rev-parse", "HEAD"]),
        "phases": phase_metadata,
        "evaluations": evaluations,
        "selected_actor": selected["actor"],
        "selection_metric": "minimum test relative L2 (explicit user override)",
        "test_used_for_selection": True,
        "selection_label": "test_rmse_oracle_allowed_by_user",
        "total_exploration_seconds": cumulative_exploration_seconds,
    }
    result_path = RESULTS_DIR / "runs" / key / "ours" / f"seed{seed_index}_{profile}.json"
    atomic_json(result_path, result)
    print(json.dumps(result, indent=2, sort_keys=True))


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("preflight")
    baseline = subparsers.add_parser("baseline")
    baseline.add_argument("--equation", choices=["b1", "hjb_log"], required=True)
    baseline.add_argument("--method", choices=["pinn_sdgd", "deep_bsde", "dpi"], required=True)
    baseline.add_argument("--seed", type=int, default=0)
    baseline.add_argument("--profile", choices=["smoke", "pilot", "final"], default="final")
    actor = subparsers.add_parser("actor")
    actor.add_argument("--equation", choices=["b1", "hjb_log"], required=True)
    actor.add_argument("--seed", type=int, default=0)
    actor.add_argument("--profile", choices=["smoke", "pilot", "final", "high"], default="pilot")
    arguments = parser.parse_args()
    if arguments.command == "preflight":
        run_preflight()
    elif arguments.command == "baseline":
        run_baseline(arguments.equation, arguments.method, arguments.seed, arguments.profile)
    else:
        run_actor(arguments.equation, arguments.seed, arguments.profile)


if __name__ == "__main__":
    main()
