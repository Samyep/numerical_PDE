"""Train, select, and evaluate the frozen component-ablation cells."""

from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import hashlib
import json
import math
from pathlib import Path
import time
from typing import Any, Callable

import numpy as np
import torch
from torch import nn

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.core import (
    make_targets,
    tangent_labels,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.common import (
    seed_everything as seed_original,
    stable_seed as stable_original,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.models import (
    TerminalAwareActor as OriginalTerminalActor,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.references import (
    hjb_log_reference,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.run_actor import (
    _state_bank as original_state_bank,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.tangent_solver import (
    evaluate_policy as evaluate_original,
    generate_labels as original_labels,
)
from invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.common import (
    state_bank as nonconvex_state_bank,
)
from invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.solver import (
    evaluate_policy as evaluate_joint,
)
from invariant_region_mlp.experiments.tangent_actor_icml_final_suite.protocol import (
    stable_seed as stable_final,
)
from invariant_region_mlp.experiments.tangent_actor_icml_final_suite.run_joint_actor import (
    _state_bank as final_joint_state_bank,
    make_driver,
    reference as joint_reference,
)
from invariant_region_mlp.experiments.tangent_actor_icml_final_suite.run_standard_actor import (
    make_benchmark,
)
from invariant_region_mlp.experiments.tangent_actor_mcv_training.labels import (
    mcv_tangent_labels,
)
from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis.common import (
    seed_everything as seed_joint,
)
from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis.models import (
    TerminalAwareActor as JointTerminalActor,
)

from .fitters import fit_joint, fit_original
from .models import UnanchoredJointActor, UnanchoredOriginalActor
from .protocol import (
    CACHE_DIR,
    CHECKPOINT_DIR,
    DIMENSION,
    EVALUATION_PATHS,
    EVALUATION_REPLICATES,
    EVALUATION_STEPS,
    HJB_CONFIG,
    INPUT_DIR,
    JOINT_CONFIG,
    PDES,
    PREDICTION_DIR,
    RAW_DIR,
    RESULTS_DIR,
    TRAINING_SEEDS,
    UPDATES,
    VARIANTS,
    ensure_dirs,
    stable_seed,
)


ORIGINAL_RESULTS = RESULTS_DIR.parent / "original_pde_benchmark_comparison"
FINAL_RESULTS = RESULTS_DIR.parent / "tangent_actor_icml_final_suite"
FRONTIER_RESULTS = RESULTS_DIR.parent / "tangent_actor_accuracy_frontier"


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _array_hash(*arrays: np.ndarray) -> str:
    digest = hashlib.sha256()
    for array in arrays:
        value = np.ascontiguousarray(array)
        digest.update(str(value.dtype).encode())
        digest.update(str(value.shape).encode())
        digest.update(value.tobytes())
    return digest.hexdigest()


def _run_id(pde: str, variant: str, seed: int) -> str:
    return f"{pde}_d100_{variant}_s{seed}"


def _base_run_id(pde: str, seed: int) -> str:
    if pde in ("joint_concave", "joint_convex"):
        return f"{pde}_d100_k2_mature_s{seed}"
    if pde == "nonconvex_coupled":
        return f"coupled_d100_k2_lam0p2_eta1_main_mcv_s{seed}"
    return f"hjb_log_d100_focused_s{seed}"


def _is_hjb(pde: str) -> bool:
    return pde == "hjb_log"


def _actor_factory(pde: str, variant: str, problem: Any) -> Callable[[], nn.Module]:
    if _is_hjb(pde):
        if variant == "a2_no_terminal_anchor":
            return lambda: UnanchoredOriginalActor(
                problem,
                width=HJB_CONFIG.width,
                depth=HJB_CONFIG.depth,
                activation=HJB_CONFIG.activation,
            )
        return lambda: OriginalTerminalActor(
            problem,
            width=HJB_CONFIG.width,
            depth=HJB_CONFIG.depth,
            activation=HJB_CONFIG.activation,
        )
    if variant == "a2_no_terminal_anchor":
        return lambda: UnanchoredJointActor(
            problem, width=JOINT_CONFIG.width, depth=JOINT_CONFIG.depth
        )
    return lambda: JointTerminalActor(
        problem, width=JOINT_CONFIG.width, depth=JOINT_CONFIG.depth
    )


def _initial_actor(
    pde: str, variant: str, seed: int, problem: Any, device: torch.device
) -> nn.Module:
    base = _base_run_id(pde, seed)
    if pde == "hjb_log":
        seed_original(stable_original("final_component_ablation_initial", base))
    elif pde == "nonconvex_coupled":
        seed_joint(seed + 1701)
    else:
        seed_joint(stable_final(base, "a0"))
    actor = _actor_factory(pde, variant, problem)().to(device)
    actor.eval()
    return actor


def _local_checkpoint(pde: str, variant: str, seed: int, iteration: int) -> Path:
    return CHECKPOINT_DIR / f"{_run_id(pde, variant, seed)}_a{iteration}.pt"


def _save_local_actor(
    path: Path,
    actor: nn.Module,
    *,
    pde: str,
    variant: str,
    seed: int,
    iteration: int,
) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    torch.save(
        {
            "state_dict": {
                key: value.detach().cpu().clone()
                for key, value in actor.state_dict().items()
            },
            "pde": pde,
            "variant": variant,
            "seed": seed,
            "iteration": iteration,
            "configuration": asdict(HJB_CONFIG if _is_hjb(pde) else JOINT_CONFIG),
            "ground_truth_used": False,
        },
        temporary,
    )
    temporary.replace(path)


def _load_local_actor(
    path: Path, pde: str, variant: str, problem: Any, device: torch.device
) -> nn.Module:
    payload = torch.load(path, map_location=device, weights_only=False)
    actor = _actor_factory(pde, variant, problem)().to(device)
    actor.load_state_dict(payload["state_dict"])
    actor.eval()
    return actor


def _reuse_checkpoint(pde: str, seed: int, iteration: int) -> Path | None:
    if pde == "hjb_log":
        return ORIGINAL_RESULTS / "checkpoints" / (
            f"actor_hjb_log_d100_focused_s{seed}_a{iteration}.pt"
        )
    if pde == "nonconvex_coupled":
        return FRONTIER_RESULTS / "checkpoints" / (
            f"coupled_d100_k2_lam0p2_eta1_main_mcv_s{seed}_mcv_a{iteration}.pt"
        )
    if seed == 0:
        return FINAL_RESULTS / "checkpoints" / (
            f"{pde}_d100_k2_mature_s0_a{iteration}.pt"
        )
    return None


def _load_reused_actor(
    pde: str, seed: int, iteration: int, problem: Any, device: torch.device
) -> nn.Module:
    path = _reuse_checkpoint(pde, seed, iteration)
    if path is None or not path.exists():
        raise FileNotFoundError(f"no exact reused checkpoint for {pde} seed {seed} a{iteration}")
    payload = torch.load(path, map_location=device, weights_only=False)
    actor = _actor_factory(pde, "a0_full", problem)().to(device)
    actor.load_state_dict(payload["state_dict"])
    actor.eval()
    return actor


def _training_states(
    pde: str,
    problem: Any,
    seed: int,
    iteration: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if pde == "hjb_log":
        bank = original_state_bank(
            problem, "focused", seed, iteration, HJB_CONFIG
        )
        return bank["train_t"], bank["train_x"], bank["valid_t"], bank["valid_x"]
    if pde in ("joint_concave", "joint_convex"):
        train_t, train_x = final_joint_state_bank(
            pde,
            DIMENSION,
            f"train_a{iteration}",
            JOINT_CONFIG.training_points,
            seed=seed,
        )
        valid_t, valid_x = final_joint_state_bank(
            pde,
            DIMENSION,
            "fixed_label_validation",
            JOINT_CONFIG.validation_points,
            seed=seed,
        )
        return train_t, train_x, valid_t, valid_x
    base = _base_run_id(pde, seed)
    train_t, train_x = nonconvex_state_bank(
        problem.name,
        DIMENSION,
        f"{base}_train_a{iteration}",
        JOINT_CONFIG.training_points,
        T=problem.T,
        seed=seed,
    )
    valid_t, valid_x = nonconvex_state_bank(
        problem.name,
        DIMENSION,
        f"{base}_fixed_label_validation",
        JOINT_CONFIG.validation_points,
        T=problem.T,
        seed=seed,
    )
    return train_t, train_x, valid_t, valid_x


def _label_seed(pde: str, seed: int, iteration: int, split: str) -> int:
    base = _base_run_id(pde, seed)
    if pde == "hjb_log":
        return stable_original(
            "original_benchmark_labels",
            "hjb_log_d100",
            "focused",
            seed,
            iteration,
            split,
        )
    if pde in ("joint_concave", "joint_convex"):
        suffix = {"train": "train_labels", "valid": "valid_labels", "replicate": "replicate_labels"}[split]
        return stable_final(base, iteration, suffix)
    # The estimator token is deliberately held at mcv so plain/MCV use CRN.
    suffix = {"train": "train_labels", "valid": "valid_labels", "replicate": "replicate_labels"}[split]
    from invariant_region_mlp.experiments.nonconvex_tangent_actor.common import stable_seed as stable_nonconvex

    return stable_nonconvex(base, "mcv", iteration, suffix)


def _generate_labels(
    pde: str,
    variant: str,
    seed: int,
    iteration: int,
    problem: Any,
    source: nn.Module,
    device: torch.device,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    run_id = _run_id(pde, variant, seed)
    cache = CACHE_DIR / f"{run_id}_a{iteration}_labels.npz"
    metadata_path = cache.with_suffix(".json")
    if cache.exists() and metadata_path.exists():
        with np.load(cache) as archive:
            arrays = {key: archive[key] for key in archive.files}
        return arrays, json.loads(metadata_path.read_text(encoding="utf-8"))
    train_t, train_x, valid_t, valid_x = _training_states(
        pde, problem, seed, iteration
    )
    estimator = "plain" if variant == "a1_no_mcv" else "mcv"
    if pde == "hjb_log":
        config = replace(HJB_CONFIG, label_estimator=estimator)
        kwargs = {
            "paths": config.paths_per_label,
            "time_steps": config.label_time_steps,
            "device": device,
            "parent_batch": config.parent_batch,
            "path_chunk": config.path_chunk,
            "estimator": estimator,
        }
        train_y, train_meta = original_labels(
            problem,
            source,
            train_t,
            train_x,
            seed=_label_seed(pde, seed, iteration, "train"),
            **kwargs,
        )
        valid_y, valid_meta = original_labels(
            problem,
            source,
            valid_t,
            valid_x,
            seed=_label_seed(pde, seed, iteration, "valid"),
            **kwargs,
        )
        replicate_y, replicate_meta = original_labels(
            problem,
            source,
            valid_t,
            valid_x,
            seed=_label_seed(pde, seed, iteration, "replicate"),
            **kwargs,
        )
        arrays = {
            "train_t": train_t,
            "train_x": train_x,
            "train_y": train_y,
            "valid_t": valid_t,
            "valid_x": valid_x,
            "valid_y": valid_y,
        }
        difference = valid_y - replicate_y
        metadata = {
            "estimator": estimator,
            "train": train_meta,
            "validation": valid_meta,
            "replication": replicate_meta,
            "label_replication_component_rmse": float(
                np.sqrt(np.mean(np.square(difference)))
            ),
            "label_replication_vector_rmse": float(
                np.sqrt(np.mean(np.square(difference).sum(axis=1)))
            ),
        }
    else:
        labeler = tangent_labels if estimator == "plain" else mcv_tangent_labels
        kwargs = {
            "paths": JOINT_CONFIG.paths_per_label,
            "time_steps": JOINT_CONFIG.label_time_steps,
            "device": device,
            "parent_batch": JOINT_CONFIG.parent_batch,
            "path_chunk": JOINT_CONFIG.path_chunk,
        }
        train_value, train_q, train_meta = labeler(
            problem,
            source,
            train_t,
            train_x,
            seed=_label_seed(pde, seed, iteration, "train"),
            **kwargs,
        )
        valid_value, valid_q, valid_meta = labeler(
            problem,
            source,
            valid_t,
            valid_x,
            seed=_label_seed(pde, seed, iteration, "valid"),
            **kwargs,
        )
        replicate_value, replicate_q, replicate_meta = labeler(
            problem,
            source,
            valid_t,
            valid_x,
            seed=_label_seed(pde, seed, iteration, "replicate"),
            **kwargs,
        )
        if train_q is None or valid_q is None or replicate_q is None:
            raise RuntimeError("joint benchmark requires gradient labels")
        train_r, train_q, train_projection = make_targets(
            problem, train_value, train_q
        )
        valid_r, valid_q, valid_projection = make_targets(
            problem, valid_value, valid_q
        )
        arrays = {
            "train_t": train_t,
            "train_x": train_x,
            "train_r": train_r,
            "train_q": train_q,
            "valid_t": valid_t,
            "valid_x": valid_x,
            "valid_r": valid_r,
            "valid_q": valid_q,
        }
        metadata = {
            "estimator": estimator,
            "train": train_meta,
            "validation": valid_meta,
            "replication": replicate_meta,
            "train_projection": train_projection,
            "validation_projection": valid_projection,
            "value_label_replication_rmse": float(
                np.sqrt(np.mean(np.square(valid_value - replicate_value)))
            ),
            "gradient_label_replication_component_rmse": float(
                np.sqrt(np.mean(np.square(valid_q - replicate_q)))
            ),
            "gradient_label_replication_vector_rmse": float(
                np.sqrt(np.mean(np.square(valid_q - replicate_q).sum(axis=1)))
            ),
        }
    np.savez_compressed(cache, **arrays)
    metadata.update({
        "pde": pde,
        "variant": variant,
        "seed": seed,
        "iteration": iteration,
        "train_bank_sha256": _array_hash(train_t, train_x),
        "validation_bank_sha256": _array_hash(valid_t, valid_x),
        "common_state_banks_across_variants": True,
        "ground_truth_used": False,
    })
    _write_json(metadata_path, metadata)
    return arrays, metadata


def _fit_seed(pde: str, seed: int, iteration: int) -> int:
    base = _base_run_id(pde, seed)
    if pde == "hjb_log":
        return stable_original(
            "original_benchmark_fit", "hjb_log_d100", "focused", seed, iteration
        )
    if pde in ("joint_concave", "joint_convex"):
        return stable_final(base, iteration, "fit")
    from invariant_region_mlp.experiments.nonconvex_tangent_actor.common import stable_seed as stable_nonconvex

    return stable_nonconvex(base, "mcv", iteration, "actor_fit", seed)


def _predict_actor(
    pde: str,
    problem: Any,
    actor: nn.Module,
    t: np.ndarray,
    x: np.ndarray,
    device: torch.device,
) -> tuple[np.ndarray | None, np.ndarray | None]:
    r_parts: list[np.ndarray] = []
    q_parts: list[np.ndarray] = []
    actor.eval()
    with torch.no_grad():
        for begin in range(0, len(t), 1024):
            end = min(begin + 1024, len(t))
            tt = torch.as_tensor(t[begin:end], dtype=torch.float32, device=device)
            xx = torch.as_tensor(x[begin:end], dtype=torch.float32, device=device)
            if pde == "hjb_log":
                r, q = actor.predict(tt, xx)
            else:
                r, q = actor.predict(tt, xx, problem.T)
            if r is not None:
                r_parts.append(r.double().cpu().numpy())
            if q is not None:
                q_parts.append(q.double().cpu().numpy())
    return (
        np.concatenate(r_parts) if r_parts else None,
        np.concatenate(q_parts) if q_parts else None,
    )


def _selector_states(
    pde: str, problem: Any, seed: int
) -> tuple[np.ndarray, np.ndarray]:
    if pde == "hjb_log":
        return problem.sample_states(
            128,
            seed=stable_final(
                "icml_standard_states", pde, DIMENSION, "fixed_point", seed
            ),
            purpose="function",
        )
    if pde in ("joint_concave", "joint_convex"):
        return final_joint_state_bank(pde, DIMENSION, "fixed_point", 128, seed=seed)
    return nonconvex_state_bank(
        problem.name,
        DIMENSION,
        f"{_base_run_id(pde, seed)}_fixed_point_states",
        128,
        T=problem.T,
        seed=20261017,
    )


def _fixed_point(
    pde: str,
    variant: str,
    seed: int,
    actor_name: str,
    problem: Any,
    actor: nn.Module,
    device: torch.device,
) -> dict[str, Any]:
    estimator = "plain" if variant == "a1_no_mcv" else "mcv"
    t, x = _selector_states(pde, problem, seed)
    labels: list[tuple[np.ndarray, np.ndarray | None]] = []
    metadata: list[dict[str, Any]] = []
    for replicate in (0, 1):
        label_seed = stable_seed(
            "component_ablation_selector", pde, seed, actor_name, replicate
        )
        if pde == "hjb_log":
            q, meta = original_labels(
                problem,
                actor,
                t,
                x,
                paths=256,
                time_steps=32,
                seed=label_seed,
                device=device,
                parent_batch=8,
                path_chunk=64,
                estimator=estimator,
            )
            labels.append((q, None))
        else:
            labeler = tangent_labels if estimator == "plain" else mcv_tangent_labels
            value, q, meta = labeler(
                problem,
                actor,
                t,
                x,
                paths=256,
                time_steps=32,
                seed=label_seed,
                device=device,
                parent_batch=8,
                path_chunk=64,
            )
            labels.append((value, q))
        metadata.append(meta)
    actor_r, actor_q = _predict_actor(pde, problem, actor, t, x, device)
    if pde == "hjb_log":
        if actor_q is None:
            raise RuntimeError("HJB Actor omitted q")
        target = 0.5 * (labels[0][0] + labels[1][0])
        raw = float(np.mean(np.square(actor_q - target)))
        noise = float(np.mean(np.square(labels[0][0] - labels[1][0])) / 4.0)
        corrected = max(raw - noise, 0.0)
        scale = max(float(np.sqrt(np.mean(np.square(target)))), 1e-3)
        total = corrected / scale**2
        row = {
            "raw_component_rmse": math.sqrt(raw),
            "noise_corrected_component_rmse": math.sqrt(corrected),
            "normalized_total_inconsistency": total,
        }
    else:
        if actor_r is None or actor_q is None or labels[0][1] is None or labels[1][1] is None:
            raise RuntimeError("joint fixed point is incomplete")
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
        total = corrected_r / r_scale**2 + corrected_q / q_scale**2
        row = {
            "value_raw_rmse": math.sqrt(raw_r),
            "value_noise_corrected_rmse": math.sqrt(corrected_r),
            "gradient_raw_component_rmse": math.sqrt(raw_q),
            "gradient_noise_corrected_component_rmse": math.sqrt(corrected_q),
            "normalized_value_inconsistency": corrected_r / r_scale**2,
            "normalized_gradient_inconsistency": corrected_q / q_scale**2,
            "normalized_total_inconsistency": total,
        }
    return {
        "actor": actor_name,
        "estimator": estimator,
        "states": 128,
        "paths_per_replication": 256,
        "time_steps": 32,
        "state_bank_sha256": _array_hash(t, x),
        "minimum_weight_ess_fraction": min(
            float(meta["weight_ess_fraction"]) for meta in metadata
        ),
        "nonfinite_count": int(sum(int(meta["nonfinite_count"]) for meta in metadata)),
        "ground_truth_used": False,
        **row,
    }


def _source_selection(pde: str, seed: int) -> tuple[int, str]:
    if pde == "hjb_log":
        path = ORIGINAL_RESULTS / "raw" / f"actor_hjb_log_d100_focused_s{seed}.json"
        return 4, str(path)
    if pde == "nonconvex_coupled":
        path = FRONTIER_RESULTS / f"training_{_base_run_id(pde, seed)}.json"
        payload = json.loads(path.read_text(encoding="utf-8"))
        return int(payload["selection"]["selected_iteration"]), str(path)
    if seed == 0:
        path = FINAL_RESULTS / "raw" / f"training_{_base_run_id(pde, seed)}.json"
        payload = json.loads(path.read_text(encoding="utf-8"))
        return int(payload["selection"]["selected_actor"][1:]), str(path)
    raise FileNotFoundError("source selection is not reusable")


def train(
    pde: str, variant: str, seed: int, device: torch.device
) -> dict[str, Any]:
    ensure_dirs()
    if variant == "a3_no_updates":
        raise ValueError("a3_no_updates does not train")
    output = RAW_DIR / f"training_{_run_id(pde, variant, seed)}.json"
    if output.exists():
        existing = json.loads(output.read_text(encoding="utf-8"))
        if existing.get("status") == "complete":
            return existing
    problem = make_benchmark(pde, DIMENSION) if _is_hjb(pde) else make_driver(pde, DIMENSION)
    if variant == "a0_full":
        try:
            selected_iteration, selection_source = _source_selection(pde, seed)
            checkpoints = [_reuse_checkpoint(pde, seed, i) for i in range(1, UPDATES + 1)]
            if all(path is not None and path.exists() for path in checkpoints):
                record = {
                    "run_id": _run_id(pde, variant, seed),
                    "pde": pde,
                    "variant": variant,
                    "seed": seed,
                    "status": "complete",
                    "reuse": True,
                    "reused_checkpoints": [str(path) for path in checkpoints],
                    "selected_iteration": selected_iteration,
                    "selected_checkpoint": str(checkpoints[selected_iteration - 1])
                    if selected_iteration
                    else None,
                    "selection_source": selection_source,
                    "ground_truth_used_for_training_or_selection": False,
                    "training_seconds": sum(
                        float(
                            phase.get(
                                "phase_seconds",
                                phase.get("phase_wall_seconds_this_invocation", 0.0),
                            )
                        )
                        for phase in json.loads(Path(selection_source).read_text()).get("phases", [])
                    ),
                }
                _write_json(output, record)
                return record
        except FileNotFoundError:
            pass
    actor = _initial_actor(pde, variant, seed, problem, device)
    phases: list[dict[str, Any]] = []
    if output.exists():
        existing = json.loads(output.read_text(encoding="utf-8"))
        phases = existing.get("phases", [])
        if phases:
            actor = _load_local_actor(
                _local_checkpoint(pde, variant, seed, len(phases)),
                pde,
                variant,
                problem,
                device,
            )
    started = time.perf_counter()
    for iteration in range(len(phases) + 1, UPDATES + 1):
        phase_started = time.perf_counter()
        arrays, label_meta = _generate_labels(
            pde, variant, seed, iteration, problem, actor, device
        )
        factory = _actor_factory(pde, variant, problem)
        if pde == "hjb_log":
            fitted, fit_meta, trace = fit_original(
                problem,
                replace(
                    HJB_CONFIG,
                    label_estimator="plain" if variant == "a1_no_mcv" else "mcv",
                ),
                actor,
                factory,
                arrays["train_t"],
                arrays["train_x"],
                arrays["train_y"],
                arrays["valid_t"],
                arrays["valid_x"],
                arrays["valid_y"],
                seed=_fit_seed(pde, seed, iteration),
                device=device,
            )
        else:
            fitted, fit_meta, trace = fit_joint(
                problem,
                JOINT_CONFIG,
                actor,
                factory,
                arrays["train_t"],
                arrays["train_x"],
                arrays["train_r"],
                arrays["train_q"],
                arrays["valid_t"],
                arrays["valid_x"],
                arrays["valid_r"],
                arrays["valid_q"],
                seed=_fit_seed(pde, seed, iteration),
                device=device,
            )
        checkpoint = _local_checkpoint(pde, variant, seed, iteration)
        _save_local_actor(
            checkpoint,
            fitted,
            pde=pde,
            variant=variant,
            seed=seed,
            iteration=iteration,
        )
        phases.append({
            "iteration": iteration,
            "checkpoint": str(checkpoint),
            "phase_seconds": time.perf_counter() - phase_started,
            "labels": label_meta,
            "training": fit_meta,
            "trace": trace,
        })
        actor = fitted
        _write_json(output, {
            "run_id": _run_id(pde, variant, seed),
            "pde": pde,
            "variant": variant,
            "seed": seed,
            "phases": phases,
            "status": "training_in_progress",
            "ground_truth_used_for_training_or_selection": False,
        })
        print(
            f"{pde} {variant} s{seed} a{iteration} "
            f"labels={label_meta['train']['wall_seconds']:.1f}s "
            f"fit={fit_meta['wall_seconds']:.1f}s",
            flush=True,
        )
    actors = [_initial_actor(pde, variant, seed, problem, device)] + [
        _load_local_actor(
            _local_checkpoint(pde, variant, seed, iteration),
            pde,
            variant,
            problem,
            device,
        )
        for iteration in range(1, UPDATES + 1)
    ]
    selector_rows = [
        _fixed_point(pde, variant, seed, f"a{iteration}", problem, actor, device)
        for iteration, actor in enumerate(actors)
    ]
    eligible = [row for row in selector_rows if row["nonfinite_count"] == 0]
    selected = min(eligible, key=lambda row: row["normalized_total_inconsistency"])
    selected_iteration = int(selected["actor"][1:])
    record = {
        "run_id": _run_id(pde, variant, seed),
        "pde": pde,
        "variant": variant,
        "seed": seed,
        "configuration": asdict(HJB_CONFIG if pde == "hjb_log" else JOINT_CONFIG),
        "phases": phases,
        "selection": {
            "criterion": "minimum independent noise-corrected fixed-point inconsistency",
            "estimator": "plain" if variant == "a1_no_mcv" else "mcv",
            "rows": selector_rows,
            "selected_actor": selected["actor"],
            "selected_iteration": selected_iteration,
            "selected_checkpoint": str(
                _local_checkpoint(pde, variant, seed, selected_iteration)
            )
            if selected_iteration
            else None,
            "ground_truth_used": False,
        },
        "selected_iteration": selected_iteration,
        "selected_checkpoint": str(
            _local_checkpoint(pde, variant, seed, selected_iteration)
        )
        if selected_iteration
        else None,
        "reuse": False,
        "training_seconds": sum(float(row["phase_seconds"]) for row in phases),
        "invocation_seconds": time.perf_counter() - started,
        "ground_truth_used_for_training_or_selection": False,
        "status": "complete",
    }
    _write_json(output, record)
    return record


def _selected_actor(
    pde: str, variant: str, seed: int, problem: Any, device: torch.device
) -> tuple[nn.Module, dict[str, Any]]:
    if variant == "a3_no_updates":
        return _initial_actor(pde, "a0_full", seed, problem, device), {
            "selected_iteration": 0,
            "training_seconds": 0.0,
            "reuse": True,
            "selection": None,
        }
    training_path = RAW_DIR / f"training_{_run_id(pde, variant, seed)}.json"
    training = json.loads(training_path.read_text(encoding="utf-8"))
    iteration = int(training["selected_iteration"])
    if iteration == 0:
        actor = _initial_actor(pde, variant, seed, problem, device)
    elif training.get("reuse"):
        actor = _load_reused_actor(pde, seed, iteration, problem, device)
    else:
        actor = _load_local_actor(
            Path(training["selected_checkpoint"]), pde, variant, problem, device
        )
    return actor, training


def _test_bank(pde: str) -> tuple[np.ndarray, np.ndarray, Path]:
    if pde == "hjb_log":
        path = ORIGINAL_RESULTS / "evaluation_bank_hjb_log_d100.npz"
        with np.load(path) as archive:
            index = archive["test_indices"]
            return archive["t"][index], archive["x"][index], path
    if pde in ("joint_concave", "joint_convex"):
        path = FINAL_RESULTS / "inputs" / f"test_bank_{pde}_d100.npz"
        with np.load(path) as archive:
            return archive["t"], archive["x"], path
    path = INPUT_DIR / "test_bank_nonconvex_coupled_d100_960.npz"
    with np.load(path) as archive:
        return archive["t"], archive["x"], path


def _truth(
    pde: str, problem: Any, t: np.ndarray, x: np.ndarray
) -> tuple[np.ndarray, np.ndarray | None, dict[str, Any]]:
    if pde == "hjb_log":
        value, q = hjb_log_reference(t, x, order=128, return_gradient=True)
        return value, q, {
            "kind": "independent Hopf-Cole Gauss-Laguerre quadrature",
            "order": 128,
            "uncertainty": "order-128 source audit",
        }
    reference, metadata = joint_reference(pde)
    coordinates = x @ problem.B[:, :2]
    working_value = reference.value(t, coordinates)
    reduced_gradient = reference.gradient(t, coordinates)
    working_q = problem.sigma * (reduced_gradient @ problem.B[:, :2].T)
    return (
        problem.working_to_original(working_value),
        problem.working_to_original(working_q),
        metadata,
    )


def evaluate(
    pde: str,
    variant: str,
    seed: int,
    replicate: int,
    device: torch.device,
) -> dict[str, Any]:
    ensure_dirs()
    problem = make_benchmark(pde, DIMENSION) if _is_hjb(pde) else make_driver(pde, DIMENSION)
    actor, training = _selected_actor(pde, variant, seed, problem, device)
    t, x, bank_path = _test_bank(pde)
    if len(t) != 960:
        raise RuntimeError(f"{pde} test bank has {len(t)} states, expected 960")
    truth, truth_q, reference = _truth(pde, problem, t, x)
    evaluation_seed = stable_seed(
        "component_ablation_evaluation", pde, seed, replicate
    )
    if pde == "hjb_log":
        output, metadata = evaluate_original(
            problem,
            actor,
            t,
            x,
            paths=EVALUATION_PATHS,
            time_steps=EVALUATION_STEPS,
            seed=evaluation_seed,
            device=device,
            parent_batch=8,
            path_chunk=64,
        )
    else:
        output, metadata = evaluate_joint(
            problem,
            actor,
            t,
            x,
            paths=EVALUATION_PATHS,
            time_steps=EVALUATION_STEPS,
            seed=evaluation_seed,
            device=device,
            parent_batch=8,
            path_chunk=64,
        )
    estimator = "plain" if variant == "a1_no_mcv" else "mcv"
    prediction_working = output[estimator]
    plain_working = output["plain"]
    mcv_working = output["mcv"]
    if pde == "hjb_log":
        prediction = prediction_working
        plain = plain_working
        mcv = mcv_working
    else:
        prediction = problem.working_to_original(prediction_working)
        plain = problem.working_to_original(plain_working)
        mcv = problem.working_to_original(mcv_working)
    error = prediction - truth
    actor_r, actor_q = _predict_actor(pde, problem, actor, t, x, device)
    result: dict[str, Any] = {
        "run_id": _run_id(pde, variant, seed),
        "pde": pde,
        "variant": variant,
        "seed": seed,
        "inference_replicate": replicate,
        "primary_replicate": replicate == 0,
        "selected_iteration": int(training["selected_iteration"]),
        "estimator": estimator,
        "states": len(t),
        "steps": EVALUATION_STEPS,
        "paths": EVALUATION_PATHS,
        "value_rmse": float(np.sqrt(np.mean(np.square(error)))),
        "relative_l2": float(np.linalg.norm(error) / np.linalg.norm(truth)),
        "rmae": float(np.sum(np.abs(error)) / np.sum(np.abs(truth))),
        "bias": float(np.mean(error)),
        "mae": float(np.mean(np.abs(error))),
        "mc_rmse_contribution": float(
            np.sqrt(
                np.mean(
                    output[f"{estimator}_pair_variance"]
                    / (EVALUATION_PATHS // 2)
                )
            )
        ),
        "plain_rmse_same_paths": float(
            np.sqrt(np.mean(np.square(plain - truth)))
        ),
        "mcv_rmse_same_paths": float(np.sqrt(np.mean(np.square(mcv - truth)))),
        "variance_reduction": metadata["variance_reduction"],
        "weight_ess_fraction": metadata["weight_ess_fraction"],
        "inference_seconds": metadata["wall_seconds"],
        "peak_gpu_bytes": metadata["peak_gpu_bytes"],
        "nonfinite_count": metadata["nonfinite_count"],
        "training_seconds": float(training.get("training_seconds", 0.0)),
        "reused_training": bool(training.get("reuse", False)),
        "test_bank": str(bank_path),
        "test_bank_file_sha256": _file_hash(bank_path),
        "test_bank_array_sha256": _array_hash(t, x),
        "reference": reference,
        "reference_array_sha256": _array_hash(truth),
        "ground_truth_loaded_only_after_checkpoint_selection": True,
        "status": "complete",
    }
    if pde == "hjb_log":
        if actor_q is not None and truth_q is not None:
            q_error = actor_q - truth_q
            result.update({
                "actor_gradient_component_rmse": float(
                    np.sqrt(np.mean(np.square(q_error)))
                ),
                "actor_gradient_vector_rmse": float(
                    np.sqrt(np.mean(np.square(q_error).sum(axis=1)))
                ),
            })
    elif actor_r is not None and actor_q is not None and truth_q is not None:
        physical_r = problem.working_to_original(actor_r)
        physical_q = problem.working_to_original(actor_q)
        q_error = physical_q - truth_q
        result.update({
            "actor_value_rmse": float(
                np.sqrt(np.mean(np.square(physical_r - truth)))
            ),
            "actor_gradient_component_rmse": float(
                np.sqrt(np.mean(np.square(q_error)))
            ),
            "actor_gradient_vector_rmse": float(
                np.sqrt(np.mean(np.square(q_error).sum(axis=1)))
            ),
        })
    prediction_path = PREDICTION_DIR / (
        f"{_run_id(pde, variant, seed)}_rep{replicate}.npz"
    )
    np.savez_compressed(
        prediction_path,
        t=t,
        x=x,
        truth=truth,
        prediction=prediction,
        plain=plain,
        mcv=mcv,
        pair_variance=output[f"{estimator}_pair_variance"],
    )
    result["prediction_path"] = str(prediction_path)
    result["prediction_file_sha256"] = _file_hash(prediction_path)
    output_path = RAW_DIR / (
        f"evaluation_{_run_id(pde, variant, seed)}_rep{replicate}.json"
    )
    _write_json(output_path, result)
    print(
        json.dumps({
            "pde": pde,
            "variant": variant,
            "seed": seed,
            "replicate": replicate,
            "rmse": result["value_rmse"],
            "rmae": result["rmae"],
        }),
        flush=True,
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pde", choices=PDES, required=True)
    parser.add_argument("--variant", choices=VARIANTS, required=True)
    parser.add_argument("--seed", type=int, choices=TRAINING_SEEDS, required=True)
    parser.add_argument("--replicate", type=int, choices=EVALUATION_REPLICATES, default=0)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--train-only", action="store_true")
    parser.add_argument("--evaluate-only", action="store_true")
    args = parser.parse_args()
    if args.variant == "a3_no_updates" and args.seed != 0:
        parser.error("a3 is deterministic; use seed 0 and inference replicates")
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("formal execution requires a CUDA allocation")
    if not args.evaluate_only and args.variant != "a3_no_updates":
        train(args.pde, args.variant, args.seed, device)
    if not args.train_only:
        evaluate(args.pde, args.variant, args.seed, args.replicate, device)


if __name__ == "__main__":
    main()
