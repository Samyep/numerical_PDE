"""Run the paired continued-q1 versus true-q2 B1 experiment.

Both branches start from the same frozen high-budget q1 checkpoint.  The
continued-q1 control sees fresh linear-policy-evaluation labels from q0; the
true-q2 branch sees fresh labels from q1.  Points, optimizer initialization,
minibatch order, architecture, and optimization budget are paired.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import math
from pathlib import Path
import time
import traceback
from typing import Any

import numpy as np
import torch

from .actor_convergence import (
    ACTOR_DEPTH,
    ACTOR_WIDTH,
    BudgetLevel,
    checkpoint_payload,
    load_converged_actor,
    make_actor,
    predict_actor,
    project_labels,
    train_budget_level,
)
from .actor_iteration import evaluate_policy_mc
from .common import (
    RESULTS_DIR,
    array_sha256,
    atomic_json,
    fixed_test_points,
    load_config,
    stable_seed,
    training_points,
)
from .equation import B1Equation
from .run_actor_quality_ceiling import _evaluate_actor
from .run_actor_quality_upgrade import (
    _atomic_npz,
    _atomic_torch_save,
    _projection_metadata,
    _sha256,
    _state_sha256,
)


D = 400
K = 40
TRAIN_POINTS = 8192
VALID_POINTS = 512
LABEL_PATHS = 512
LABEL_TIME_STEPS = 32
LABEL_PARENT_BATCH_SIZE = 16
LABEL_PATH_CHUNK_SIZE = 128
OPTIMIZER_STEPS = 6000
LEARNING_RATE = 1.0e-4

OUTPUT_DIR = RESULTS_DIR / "actor_q2_control"
CACHE_DIR = OUTPUT_DIR / "cache"
CHECKPOINT_DIR = OUTPUT_DIR / "checkpoints"
PREDICTION_DIR = OUTPUT_DIR / "predictions"
Q0_CHECKPOINT = (
    RESULTS_DIR
    / "actor_iteration_pilot"
    / "checkpoints"
    / "actor_iteration_d400_k40_q0.pt"
)
Q1_CHECKPOINT = (
    RESULTS_DIR
    / "actor_quality_upgrade"
    / "upgrade"
    / "checkpoints"
    / "actor_quality_upgrade_d400_k40_q1.pt"
)

BRANCHES = {
    "continued_q1": {
        "label_policy": "q0",
        "source_checkpoint": Q0_CHECKPOINT,
        "actor_iteration": 1,
        "checkpoint": CHECKPOINT_DIR / "actor_q2_control_continued_q1.pt",
    },
    "true_q2": {
        "label_policy": "q1",
        "source_checkpoint": Q1_CHECKPOINT,
        "actor_iteration": 2,
        "checkpoint": CHECKPOINT_DIR / "actor_q2_control_true_q2.pt",
    },
}


def _budget(branch: str) -> BudgetLevel:
    return BudgetLevel(
        f"actor_q2_control_{branch}",
        LABEL_TIME_STEPS,
        LABEL_PATHS,
        1,
        TRAIN_POINTS,
        VALID_POINTS,
        OPTIMIZER_STEPS,
        OPTIMIZER_STEPS,
        LEARNING_RATE,
    )


def protocol_record() -> dict[str, Any]:
    """Return the frozen protocol without requiring a GPU."""

    point_seed = stable_seed("B1", "actor_q2_control_points", D, K)
    optimizer_seed = stable_seed("B1", "actor_q2_control_optimizer", D, K)
    return {
        "d": D,
        "k": K,
        "training_points": TRAIN_POINTS,
        "validation_points": VALID_POINTS,
        "paths_per_label": LABEL_PATHS,
        "sde_steps_per_label": LABEL_TIME_STEPS,
        "architecture": {
            "hidden_width": ACTOR_WIDTH,
            "hidden_layers": ACTOR_DEPTH,
            "activation": "silu",
        },
        "optimizer_steps": OPTIMIZER_STEPS,
        "learning_rate": LEARNING_RATE,
        "point_seed": point_seed,
        "optimizer_seed": optimizer_seed,
        "shared_point_realization": True,
        "shared_optimizer_minibatch_sequence": True,
        "independent_label_mc_streams": True,
        "initial_actor": "frozen high-budget q1",
        "branches": {
            name: {
                "label_policy": spec["label_policy"],
                "actor_iteration": spec["actor_iteration"],
                "label_seed_namespaces": {
                    split: f"actor_q2_control_{name}_{split}_labels"
                    for split in ("train", "valid")
                },
                "budget": asdict(_budget(name)),
            }
            for name, spec in BRANCHES.items()
        },
        "final_evaluation": (
            "64-step martingale CV, 256 paths, 3 seeds, shared Brownian banks"
        ),
        "final_evaluation_seed_namespace": "actor_iteration_direct_test",
        "analytic_reference_first_used_after_both_checkpoints_frozen": True,
        "analytic_reference_used_for_training": False,
        "analytic_reference_used_for_checkpoint_selection": False,
        "actor_structure_prior": (
            "terminal-gradient-derived rank-40 output subspace with supplied k"
        ),
        "training_seed_repetitions": 1,
        "maximum_actor_iteration": 2,
        "q3_or_q4_permitted": False,
    }


def _validate_label_cache(
    cache: Path,
    metadata_path: Path,
    *,
    point_hash: str,
    actor_hash: str,
    namespace: str,
    expected_points: int,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    with np.load(cache, allow_pickle=False) as payload:
        observed_point_hash = str(payload["point_hash"].item())
        observed_actor_hash = str(payload["actor_hash"].item())
        values = np.asarray(payload["values"], dtype=np.float64)
        labels = np.asarray(payload["labels"], dtype=np.float64)
    metadata = json.loads(metadata_path.read_text())
    if observed_point_hash != point_hash:
        raise RuntimeError(f"point hash mismatch in {cache}")
    if observed_actor_hash != actor_hash:
        raise RuntimeError(f"source actor hash mismatch in {cache}")
    if metadata.get("seed_namespace") != namespace:
        raise RuntimeError(f"seed namespace mismatch in {metadata_path}")
    if values.shape != (expected_points,) or labels.shape != (expected_points, D):
        raise RuntimeError(f"invalid label-cache shapes in {cache}")
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(labels)):
        raise RuntimeError(f"nonfinite values in {cache}")
    return values, labels, metadata


def _label_cache(
    equation: B1Equation,
    source_actor: torch.nn.Module,
    t: np.ndarray,
    x: np.ndarray,
    *,
    branch: str,
    split: str,
    device: torch.device,
    force: bool,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    cache = CACHE_DIR / f"{branch}_{split}.npz"
    metadata_path = cache.with_suffix(".json")
    point_hash = array_sha256(t, x)
    actor_hash = _state_sha256(source_actor)
    namespace = f"actor_q2_control_{branch}_{split}_labels"
    if cache.is_file() and metadata_path.is_file() and not force:
        return _validate_label_cache(
            cache,
            metadata_path,
            point_hash=point_hash,
            actor_hash=actor_hash,
            namespace=namespace,
            expected_points=len(t),
        )

    values, labels, cost = evaluate_policy_mc(
        equation,
        source_actor,
        t,
        x,
        time_steps=LABEL_TIME_STEPS,
        paths=LABEL_PATHS,
        repetition=0,
        seed_namespace=namespace,
        device=device,
        parent_batch_size=LABEL_PARENT_BATCH_SIZE,
        path_chunk_size=LABEL_PATH_CHUNK_SIZE,
        compute_z=True,
    )
    if labels is None:
        raise RuntimeError("policy evaluation did not return gradient labels")
    metadata = {
        "status": "complete",
        "branch": branch,
        "split": split,
        "points": len(t),
        "paths_per_label": LABEL_PATHS,
        "sde_steps_per_label": LABEL_TIME_STEPS,
        "point_sha256": point_hash,
        "source_actor_sha256": actor_hash,
        "source_policy": BRANCHES[branch]["label_policy"],
        "seed_namespace": namespace,
        "label_source": (
            "pathwise gradient through controlled-SDE linear policy evaluation"
        ),
        "cost": cost,
        "analytic_reference_used": False,
    }
    _atomic_npz(
        cache,
        values=values,
        labels=labels,
        point_hash=np.asarray(point_hash),
        actor_hash=np.asarray(actor_hash),
    )
    atomic_json(metadata_path, metadata)
    return values, labels, metadata


def _load_frozen_checkpoint(
    path: Path,
    device: torch.device,
    *,
    expected_branch: str,
    expected_initial_state_hash: str,
) -> tuple[torch.nn.Module, dict[str, Any]]:
    model, metadata = load_converged_actor(path, device)
    if metadata.get("experiment") != "actor_q2_control":
        raise RuntimeError(f"checkpoint is from the wrong experiment: {path}")
    if metadata.get("branch") != expected_branch:
        raise RuntimeError(f"checkpoint branch mismatch: {path}")
    if metadata.get("initial_actor_state_sha256") != expected_initial_state_hash:
        raise RuntimeError(f"checkpoint initial-actor hash mismatch: {path}")
    return model, metadata


def _train_branch(
    equation: B1Equation,
    initial_actor: torch.nn.Module,
    source_actor: torch.nn.Module,
    train_t: np.ndarray,
    train_x: np.ndarray,
    valid_t: np.ndarray,
    valid_x: np.ndarray,
    *,
    branch: str,
    device: torch.device,
    force: bool,
) -> tuple[torch.nn.Module, dict[str, Any]]:
    spec = BRANCHES[branch]
    checkpoint = Path(spec["checkpoint"])
    initial_state_hash = _state_sha256(initial_actor)
    if checkpoint.is_file() and not force:
        return _load_frozen_checkpoint(
            checkpoint,
            device,
            expected_branch=branch,
            expected_initial_state_hash=initial_state_hash,
        )

    _, train_raw, train_label_metadata = _label_cache(
        equation,
        source_actor,
        train_t,
        train_x,
        branch=branch,
        split="train",
        device=device,
        force=force,
    )
    _, valid_raw, valid_label_metadata = _label_cache(
        equation,
        source_actor,
        valid_t,
        valid_x,
        branch=branch,
        split="valid",
        device=device,
        force=force,
    )
    label_basis = initial_actor.label_basis.detach().cpu().numpy()
    train_y = project_labels(train_raw, label_basis)
    valid_y = project_labels(valid_raw, label_basis)
    initial_validation_prediction = predict_actor(
        initial_actor, valid_t, valid_x, device=device
    )

    model = make_actor(equation, device, label_basis)
    model.load_state_dict(initial_actor.state_dict())
    optimizer_seed = stable_seed("B1", "actor_q2_control_optimizer", D, K)
    training = train_budget_level(
        model,
        equation,
        _budget(branch),
        train_t,
        train_x,
        train_y,
        valid_t,
        valid_x,
        valid_y,
        device=device,
        seed=optimizer_seed,
    )
    if training["steps_executed"] != OPTIMIZER_STEPS:
        raise RuntimeError(f"{branch} did not execute exactly {OPTIMIZER_STEPS} steps")
    current_validation_prediction = predict_actor(
        model, valid_t, valid_x, device=device
    )
    prediction_difference = current_validation_prediction - initial_validation_prediction
    record = {
        "status": "complete",
        "experiment": "actor_q2_control",
        "branch": branch,
        "actor_iteration": spec["actor_iteration"],
        "label_policy": spec["label_policy"],
        "initial_actor": "frozen high-budget q1",
        "initial_actor_checkpoint": str(Q1_CHECKPOINT.relative_to(RESULTS_DIR.parent)),
        "initial_actor_checkpoint_sha256": _sha256(Q1_CHECKPOINT),
        "initial_actor_state_sha256": initial_state_hash,
        "label_source_checkpoint": str(
            Path(spec["source_checkpoint"]).relative_to(RESULTS_DIR.parent)
        ),
        "label_source_checkpoint_sha256": _sha256(Path(spec["source_checkpoint"])),
        "label_source_actor_state_sha256": _state_sha256(source_actor),
        "train_point_sha256": array_sha256(train_t, train_x),
        "validation_point_sha256": array_sha256(valid_t, valid_x),
        "optimizer_seed": optimizer_seed,
        "train_labels": train_label_metadata,
        "validation_labels": valid_label_metadata,
        "train_projection": _projection_metadata(train_raw, train_y),
        "validation_projection": _projection_metadata(valid_raw, valid_y),
        "label_generation_seconds": (
            train_label_metadata["cost"]["wall_seconds"]
            + valid_label_metadata["cost"]["wall_seconds"]
        ),
        "training": training,
        "actor_validation_prediction_change_rmse": math.sqrt(
            float(np.mean(np.sum(prediction_difference**2, axis=1)))
        ),
        "analytic_reference_used_for_training": False,
        "analytic_reference_used_for_checkpoint_selection": False,
    }
    _atomic_torch_save(checkpoint, checkpoint_payload(model, record))
    print(
        json.dumps(
            {
                "branch": branch,
                "checkpoint_frozen": True,
                "best_step": training["best_step"],
                "steps_executed": training["steps_executed"],
                "selected_validation_label_rmse": training[
                    "selected_validation_label_rmse"
                ],
                "label_generation_seconds": record["label_generation_seconds"],
                "training_seconds": training["training_seconds"],
            }
        ),
        flush=True,
    )
    return model, record


def _mean_policy_gap_rms(row: dict[str, Any]) -> float:
    return float(
        np.mean([seed["analytic_policy_gap_rms"] for seed in row["seed_rows"]])
    )


def _evaluation_seconds(row: dict[str, Any]) -> float:
    return float(sum(cost["wall_seconds"] for cost in row["costs"]))


def _label_comparison() -> dict[str, Any]:
    arrays: dict[str, dict[str, np.ndarray]] = {}
    for branch in BRANCHES:
        arrays[branch] = {}
        for split in ("train", "valid"):
            with np.load(CACHE_DIR / f"{branch}_{split}.npz", allow_pickle=False) as payload:
                arrays[branch][split] = np.asarray(payload["labels"], dtype=np.float64)
    result: dict[str, Any] = {}
    for split in ("train", "valid"):
        difference = arrays["true_q2"][split] - arrays["continued_q1"][split]
        result[split] = {
            "raw_label_difference_vector_rmse": math.sqrt(
                float(np.mean(np.sum(difference**2, axis=1)))
            ),
            "continued_q1_raw_label_sha256": array_sha256(
                arrays["continued_q1"][split]
            ),
            "true_q2_raw_label_sha256": array_sha256(arrays["true_q2"][split]),
        }
    return result


def _comparison(rows: dict[str, dict[str, Any]]) -> dict[str, Any]:
    baseline = rows["frozen_q1"]
    continued = rows["continued_q1"]
    q2 = rows["true_q2"]

    def improvement(new: dict[str, Any], old: dict[str, Any], key: str) -> float:
        return 1.0 - float(new[key]) / float(old[key])

    return {
        "continued_q1_relative_value_rmse_improvement_over_frozen_q1": improvement(
            continued, baseline, "corrected_value_rmse"
        ),
        "true_q2_relative_value_rmse_improvement_over_frozen_q1": improvement(
            q2, baseline, "corrected_value_rmse"
        ),
        "true_q2_relative_value_rmse_improvement_over_continued_q1": improvement(
            q2, continued, "corrected_value_rmse"
        ),
        "continued_q1_relative_mean_policy_gap_improvement_over_frozen_q1": improvement(
            continued, baseline, "mean_analytic_policy_gap"
        ),
        "true_q2_relative_mean_policy_gap_improvement_over_frozen_q1": improvement(
            q2, baseline, "mean_analytic_policy_gap"
        ),
        "true_q2_relative_mean_policy_gap_improvement_over_continued_q1": improvement(
            q2, continued, "mean_analytic_policy_gap"
        ),
        "true_q2_beats_frozen_q1_on_value_rmse": (
            q2["corrected_value_rmse"] < baseline["corrected_value_rmse"]
        ),
        "true_q2_beats_continued_q1_on_value_rmse": (
            q2["corrected_value_rmse"] < continued["corrected_value_rmse"]
        ),
        "true_q2_beats_frozen_q1_on_mean_policy_gap": (
            q2["mean_analytic_policy_gap"]
            < baseline["mean_analytic_policy_gap"]
        ),
        "true_q2_beats_continued_q1_on_mean_policy_gap": (
            q2["mean_analytic_policy_gap"]
            < continued["mean_analytic_policy_gap"]
        ),
        "interpretation_requires_independent_training_seed_confirmation": True,
    }


def run(*, force: bool = False) -> dict[str, Any]:
    result_path = OUTPUT_DIR / "result.json"
    if result_path.is_file() and not force:
        result = json.loads(result_path.read_text())
        if result.get("status") == "complete":
            print(json.dumps(result["headline"], indent=2), flush=True)
            return result
    missing = [str(path) for path in (Q0_CHECKPOINT, Q1_CHECKPOINT) if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing frozen actor checkpoints: {missing}")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda":
        raise RuntimeError("q2-control experiment requires an allocated GPU")
    torch.set_float32_matmul_precision("high")
    config = load_config()
    equation = B1Equation(
        d=D,
        k=K,
        coefficient_scale=float(config["pde"]["frozen_coefficient_scale"]),
    )
    started = time.perf_counter()

    initial_actor, _ = load_converged_actor(Q1_CHECKPOINT, device)
    q0_actor, _ = load_converged_actor(Q0_CHECKPOINT, device)
    q1_label_actor, _ = load_converged_actor(Q1_CHECKPOINT, device)
    initial_state_hash = _state_sha256(initial_actor)
    if _state_sha256(q1_label_actor) != initial_state_hash:
        raise RuntimeError("q1 label source does not match the initial q1 actor")

    protocol = protocol_record()
    point_seed = int(protocol["point_seed"])
    train_t, train_x = training_points(
        D, K, "actor_q2_control_train", point_seed, TRAIN_POINTS
    )
    valid_t, valid_x = training_points(
        D, K, "actor_q2_control_valid", point_seed, VALID_POINTS
    )
    point_hashes = {
        "train": array_sha256(train_t, train_x),
        "valid": array_sha256(valid_t, valid_x),
    }

    trained: dict[str, tuple[torch.nn.Module, dict[str, Any]]] = {}
    for branch, source_actor in (
        ("continued_q1", q0_actor),
        ("true_q2", q1_label_actor),
    ):
        trained[branch] = _train_branch(
            equation,
            initial_actor,
            source_actor,
            train_t,
            train_x,
            valid_t,
            valid_x,
            branch=branch,
            device=device,
            force=force,
        )

    # Exact references are constructed only after both learned checkpoints exist.
    for branch, spec in BRANCHES.items():
        if not Path(spec["checkpoint"]).is_file():
            raise RuntimeError(f"{branch} checkpoint was not frozen")
    points = fixed_test_points(D, K, int(config["stage_b"]["test_points"]))
    times = np.zeros(len(points), dtype=np.float64)
    exact_u = equation.exact_u(times, points)
    exact_z = equation.exact_z(times, points)

    evaluation_checkpoints = {
        "frozen_q1": Q1_CHECKPOINT,
        "continued_q1": Path(BRANCHES["continued_q1"]["checkpoint"]),
        "true_q2": Path(BRANCHES["true_q2"]["checkpoint"]),
    }
    rows: list[dict[str, Any]] = []
    prediction_arrays: dict[str, np.ndarray] = {}
    for name, checkpoint in evaluation_checkpoints.items():
        row, arrays = _evaluate_actor(
            name,
            checkpoint,
            equation,
            times,
            points,
            exact_u,
            exact_z,
            device=device,
        )
        row["analytic_policy_gap_rms"] = _mean_policy_gap_rms(row)
        row["total_evaluation_seconds"] = _evaluation_seconds(row)
        rows.append(row)
        prediction_arrays.update(arrays)
        print(
            json.dumps(
                {
                    "actor": name,
                    "direct_actor_gradient_rmse": row[
                        "direct_actor_gradient_rmse"
                    ],
                    "corrected_value_rmse": row["corrected_value_rmse"],
                    "corrected_value_ci95": row["corrected_value_ci95"],
                    "mean_analytic_policy_gap": row[
                        "mean_analytic_policy_gap"
                    ],
                    "analytic_policy_gap_rms": row[
                        "analytic_policy_gap_rms"
                    ],
                }
            ),
            flush=True,
        )
    _atomic_npz(PREDICTION_DIR / "paired_predictions.npz", **prediction_arrays)

    row_lookup = {row["actor"]: row for row in rows}
    comparison = _comparison(row_lookup)
    label_comparison = _label_comparison()
    cost_summary = {
        branch: {
            "label_generation_seconds": record["label_generation_seconds"],
            "optimizer_training_seconds": record["training"]["training_seconds"],
            "final_three_seed_evaluation_seconds": row_lookup[branch][
                "total_evaluation_seconds"
            ],
            "label_path_steps": (
                record["train_labels"]["cost"]["path_steps"]
                + record["validation_labels"]["cost"]["path_steps"]
            ),
            "optimizer_steps_executed": record["training"]["steps_executed"],
        }
        for branch, (_, record) in trained.items()
    }
    cost_summary["frozen_q1"] = {
        "incremental_label_generation_seconds": 0.0,
        "incremental_optimizer_training_seconds": 0.0,
        "final_three_seed_evaluation_seconds": row_lookup["frozen_q1"][
            "total_evaluation_seconds"
        ],
    }
    headline = {
        name: {
            "corrected_value_rmse": row_lookup[name]["corrected_value_rmse"],
            "corrected_value_ci95": row_lookup[name]["corrected_value_ci95"],
            "corrected_value_bias": row_lookup[name]["corrected_value_bias"],
            "direct_actor_gradient_rmse": row_lookup[name][
                "direct_actor_gradient_rmse"
            ],
            "mean_analytic_policy_gap": row_lookup[name][
                "mean_analytic_policy_gap"
            ],
            "analytic_policy_gap_rms": row_lookup[name][
                "analytic_policy_gap_rms"
            ],
        }
        for name in evaluation_checkpoints
    }
    headline["comparison"] = comparison
    payload = {
        "schema_version": 1,
        "status": "complete",
        "device": str(device),
        "protocol": protocol,
        "point_hashes": point_hashes,
        "initial_actor_checkpoint_sha256": _sha256(Q1_CHECKPOINT),
        "initial_actor_state_sha256": initial_state_hash,
        "branch_training": {
            branch: record for branch, (_, record) in trained.items()
        },
        "label_comparison": label_comparison,
        "rows": rows,
        "comparison": comparison,
        "cost_summary": cost_summary,
        "headline": headline,
        "elapsed_seconds": time.perf_counter() - started,
        "failed_runs": 0,
        "nan_failures": 0,
        "memory_failures": 0,
        "q3_or_q4_run": False,
        "limitations": [
            "The paired comparison uses one shared training seed.",
            (
                "Evaluation confidence intervals condition on the frozen trained "
                "actors and do not measure training-seed variability."
            ),
            "The rank-40 output subspace is a structural prior.",
        ],
    }
    atomic_json(result_path, payload)
    print(json.dumps(headline, indent=2), flush=True)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--describe-protocol", action="store_true")
    args = parser.parse_args()
    if args.describe_protocol:
        print(json.dumps(protocol_record(), indent=2, sort_keys=True))
        return
    try:
        run(force=args.force)
    except Exception as error:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        atomic_json(
            OUTPUT_DIR / "failure.json",
            {
                "status": "failed",
                "error": repr(error),
                "traceback": traceback.format_exc(),
            },
        )
        raise


if __name__ == "__main__":
    main()
