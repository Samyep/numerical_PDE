"""Train progressively higher-quality actors and run paired Stage B evaluation."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import time
import traceback
from typing import Any

import numpy as np
import torch

from .actor_convergence import (
    ACTOR_DEPTH,
    ACTOR_WIDTH,
    AUDIT_BANK_REPLICAS,
    AUDIT_M,
    AUDIT_N,
    AUDIT_POINTS,
    BUDGET_LEVELS,
    PROTOCOL_VERSION,
    audit_record,
    budget_protocol,
    checkpoint_payload,
    convergence_decision,
    generate_label_ensemble,
    load_converged_actor,
    make_actor,
    predict_actor,
    train_budget_level,
)
from .actor_mlp import actor_linearization_identity, evaluate_actor_mlp
from .common import (
    RESULTS_DIR,
    array_sha256,
    atomic_json,
    fixed_test_points,
    load_config,
    seed_everything,
    stable_seed,
    training_points,
)
from .equation import B1Equation
from .metrics import prediction_metrics
from .neural_common import count_parameters


CONVERGENCE_DIR = RESULTS_DIR / "actor_convergence"
CACHE_DIR = CONVERGENCE_DIR / "cache"
CHECKPOINT_DIR = CONVERGENCE_DIR / "checkpoints"
RAW_DIR = CONVERGENCE_DIR / "raw"
PREDICTION_DIR = CONVERGENCE_DIR / "predictions"


def _atomic_torch_save(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    torch.save(payload, temporary)
    temporary.replace(path)


def _atomic_predictions(path: Path, predictions: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    arrays: dict[str, np.ndarray] = {}
    for item in predictions:
        repetition = int(item["repetition"])
        arrays[f"value_seed{repetition}"] = item["value"]
        arrays[f"z_seed{repetition}"] = item["z"]
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}.npz")
    np.savez_compressed(temporary, **arrays)
    temporary.replace(path)


def _cell_paths(d: int, k: int) -> dict[str, Path]:
    stem = f"d{d}_k{k}"
    return {
        "result": RAW_DIR / f"actor_converged_{stem}.json",
        "progress": RAW_DIR / f"training_progress_{stem}.json",
        "final_checkpoint": CHECKPOINT_DIR / f"actor_converged_{stem}.pt",
        "predictions": PREDICTION_DIR / f"actor_converged_{stem}.npz",
    }


def _level_paths(d: int, k: int, level_name: str) -> tuple[Path, Path]:
    stem = f"actor_d{d}_k{k}_{level_name}"
    return CHECKPOINT_DIR / f"{stem}.pt", RAW_DIR / f"{stem}.json"


def _audit_data(
    equation: B1Equation,
    *,
    device: torch.device,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    seed = stable_seed("B1", "actor_convergence", equation.d, equation.k)
    audit_t, audit_x = training_points(
        equation.d,
        equation.k,
        "actor_convergence_audit",
        seed,
        AUDIT_POINTS,
    )
    cell_cache = CACHE_DIR / f"d{equation.d}_k{equation.k}"
    bank_a, metadata_a = generate_label_ensemble(
        equation,
        audit_t,
        audit_x,
        n=AUDIT_N,
        M=AUDIT_M,
        replicas=AUDIT_BANK_REPLICAS,
        repetition_offset=300_000,
        seed_namespace="actor_convergence_audit_bank_a",
        device=device,
        cache_path=cell_cache / "audit_bank_a.npz",
    )
    bank_b, metadata_b = generate_label_ensemble(
        equation,
        audit_t,
        audit_x,
        n=AUDIT_N,
        M=AUDIT_M,
        replicas=AUDIT_BANK_REPLICAS,
        repetition_offset=400_000,
        seed_namespace="actor_convergence_audit_bank_b",
        device=device,
        cache_path=cell_cache / "audit_bank_b.npz",
    )
    metadata = {
        "point_sha256": array_sha256(audit_t, audit_x),
        "bank_a": metadata_a,
        "bank_b": metadata_b,
        "nonterminal_oracle_reference_used": False,
    }
    return audit_t, audit_x, bank_a, bank_b, metadata


def _level_data(
    equation: B1Equation,
    level_index: int,
    *,
    device: torch.device,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    dict[str, Any],
]:
    level = BUDGET_LEVELS[level_index]
    seed = stable_seed("B1", "actor_convergence", equation.d, equation.k)
    train_t, train_x = training_points(
        equation.d,
        equation.k,
        f"actor_convergence_{level.name}_train",
        seed,
        level.training_points,
    )
    valid_t, valid_x = training_points(
        equation.d,
        equation.k,
        f"actor_convergence_{level.name}_valid",
        seed,
        level.validation_points,
    )
    cell_cache = CACHE_DIR / f"d{equation.d}_k{equation.k}"
    train_y, train_metadata = generate_label_ensemble(
        equation,
        train_t,
        train_x,
        n=level.n,
        M=level.M,
        replicas=level.replicas,
        repetition_offset=100_000 + 1000 * level_index,
        seed_namespace=f"actor_convergence_{level.name}_train",
        device=device,
        cache_path=cell_cache / f"{level.name}_train.npz",
    )
    valid_y, valid_metadata = generate_label_ensemble(
        equation,
        valid_t,
        valid_x,
        n=level.n,
        M=level.M,
        replicas=level.replicas,
        repetition_offset=200_000 + 1000 * level_index,
        seed_namespace=f"actor_convergence_{level.name}_valid",
        device=device,
        cache_path=cell_cache / f"{level.name}_valid.npz",
    )
    metadata = {
        "train": train_metadata,
        "validation": valid_metadata,
        "nonterminal_oracle_reference_used": False,
    }
    return train_t, train_x, train_y, valid_t, valid_x, valid_y, metadata


def _postselection_exact_diagnostic(
    equation: B1Equation,
    model: torch.nn.Module,
    audit_t: np.ndarray,
    audit_x: np.ndarray,
    *,
    device: torch.device,
) -> dict[str, Any]:
    """This function is called only after the checkpoint has been frozen."""

    prediction = predict_actor(model, audit_t, audit_x, device=device)
    exact = equation.exact_z(audit_t, audit_x)
    squared_error = float(np.mean(np.sum((prediction - exact) ** 2, axis=1)))
    exact_energy = float(np.mean(np.sum(exact**2, axis=1)))
    return {
        "actor_gradient_rmse": float(np.sqrt(squared_error)),
        "actor_relative_gradient_error": float(
            np.sqrt(squared_error / max(exact_energy, np.finfo(np.float64).tiny))
        ),
        "actor_bias_norm": float(np.linalg.norm(np.mean(prediction - exact, axis=0))),
        "mean_driver_gap_at_exact_z": 0.5 * squared_error,
        "reference_timing": "after final checkpoint selection",
        "used_for_training": False,
        "used_for_early_stopping": False,
        "used_for_budget_convergence": False,
        "used_for_checkpoint_selection": False,
    }


def _paired_stage_b_evaluation(
    equation: B1Equation,
    model: torch.nn.Module,
    *,
    device: torch.device,
    training_seconds: float,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], str]:
    config = load_config()
    stage = config["stage_b"]
    points = fixed_test_points(equation.d, equation.k, int(stage["test_points"]))
    t = np.zeros(len(points), dtype=np.float64)
    test_hash = array_sha256(points)
    reference_u = equation.exact_u(t, points)
    reference_z = equation.exact_z(t, points)
    rows: list[dict[str, Any]] = []
    predictions: list[dict[str, Any]] = []
    identity_error = actor_linearization_identity(model, equation, device)
    for repetition in range(int(stage["mlp_seeds"])):
        state, cost = evaluate_actor_mlp(
            equation,
            model,
            t,
            points,
            n=int(stage["mlp_n"]),
            M=int(stage["mlp_M"]),
            repetition=repetition,
            device=device,
        )
        row: dict[str, Any] = {
            "stage": "actor_convergence_evaluation",
            "method": "actor_mlp_converged",
            "d": equation.d,
            "k": equation.k,
            "K": equation.K,
            "coefficient_scale": (
                None if equation.k == 1 else equation.coefficient_scale
            ),
            "repetition": repetition,
            "test_points": len(points),
            "test_point_sha256": test_hash,
            "test_time": 0.0,
            "test_distribution": "x~N(0,I_d)",
            "training_seconds": training_seconds,
            "inference_seconds": cost["wall_seconds"],
            "parameter_count": count_parameters(model),
            "actor_identity_max_abs_error": identity_error,
        }
        row.update(
            prediction_metrics(
                state[:, 0], state[:, 1:], reference_u, reference_z
            )
        )
        row.update({f"cost_{key}": value for key, value in cost.items()})
        rows.append(row)
        predictions.append(
            {
                "repetition": repetition,
                "value": state[:, 0],
                "z": state[:, 1:],
            }
        )
    return rows, predictions, test_hash


def run_cell(d: int, k: int, *, force: bool = False) -> dict[str, Any]:
    paths = _cell_paths(d, k)
    if paths["result"].is_file() and not force:
        existing = json.loads(paths["result"].read_text())
        existing_protocol = existing.get("protocol", {})
        refreshed_decision = convergence_decision(existing.get("levels", []))
        if (
            existing.get("status") == "complete"
            and existing_protocol.get("protocol_version") == PROTOCOL_VERSION
            and refreshed_decision["converged"]
        ):
            print(json.dumps({"status": "skipped_complete", "d": d, "k": k}), flush=True)
            return existing

    config = load_config()
    equation = B1Equation(
        d=d,
        k=k,
        coefficient_scale=float(config["pde"]["frozen_coefficient_scale"]),
    )
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda":
        raise RuntimeError("actor convergence study requires an allocated GPU")
    torch.set_float32_matmul_precision("high")
    cell_started = time.perf_counter()
    audit_t, audit_x, audit_a, audit_b, audit_metadata = _audit_data(
        equation, device=device
    )

    initialization_seed = stable_seed(
        "B1", "actor_convergence_model", equation.d, equation.k
    )
    seed_everything(initialization_seed)
    model = make_actor(equation, device)
    records: list[dict[str, Any]] = []
    previous_prediction: np.ndarray | None = None
    selected_checkpoint: Path | None = None
    decision: dict[str, Any] = {
        "converged": False,
        "reason": "not_started",
        "stable_transitions": 0,
    }

    for level_index, level in enumerate(BUDGET_LEVELS):
        checkpoint_path, record_path = _level_paths(d, k, level.name)
        if checkpoint_path.is_file() and record_path.is_file() and not force:
            model, _ = load_converged_actor(checkpoint_path, device)
            record = json.loads(record_path.read_text())
            prediction = predict_actor(model, audit_t, audit_x, device=device)
        else:
            (
                train_t,
                train_x,
                train_y,
                valid_t,
                valid_x,
                valid_y,
                label_metadata,
            ) = _level_data(equation, level_index, device=device)
            training_metadata = train_budget_level(
                model,
                equation,
                level,
                train_t,
                train_x,
                train_y,
                valid_t,
                valid_x,
                valid_y,
                device=device,
                seed=stable_seed(
                    "B1", "actor_convergence_optimizer", d, k, level.name
                ),
            )
            prediction = predict_actor(model, audit_t, audit_x, device=device)
            record = {
                "schema_version": 1,
                "status": "complete",
                "d": d,
                "k": k,
                "budget_index": level_index + 1,
                "budget": {
                    "name": level.name,
                    "n": level.n,
                    "M": level.M,
                    "replicas": level.replicas,
                    "training_points": level.training_points,
                    "validation_points": level.validation_points,
                    "max_steps": level.max_steps,
                    "minimum_steps": level.minimum_steps,
                },
                "labels": label_metadata,
                "training": training_metadata,
                "optimizer_converged": training_metadata["optimizer_converged"],
                **audit_record(
                    prediction,
                    audit_a,
                    audit_b,
                    previous_prediction=previous_prediction,
                ),
                "nonterminal_oracle_reference_used_for_training": False,
                "nonterminal_oracle_reference_used_for_selection": False,
            }
            _atomic_torch_save(
                checkpoint_path, checkpoint_payload(model, record)
            )
            atomic_json(record_path, record)

        records.append(record)
        previous_prediction = prediction
        selected_checkpoint = checkpoint_path
        decision = convergence_decision(records)
        progress = {
            "schema_version": 1,
            "status": "converged" if decision["converged"] else "running",
            "d": d,
            "k": k,
            "protocol": budget_protocol(),
            "audit": audit_metadata,
            "levels": records,
            "decision": decision,
        }
        atomic_json(paths["progress"], progress)
        print(
            json.dumps(
                {
                    "d": d,
                    "k": k,
                    "level": level.name,
                    "steps": record["training"]["steps_executed"],
                    "optimizer_converged": record["optimizer_converged"],
                    "audit_target_rmse": record["audit_target_rmse"],
                    "audit_noise_corrected_rmse": record[
                        "audit_noise_corrected_rmse"
                    ],
                    "prediction_relative_change": record[
                        "prediction_relative_change"
                    ],
                    "decision": decision,
                }
            ),
            flush=True,
        )
        if decision["converged"]:
            break

    if selected_checkpoint is None:
        raise RuntimeError("no actor budget level was executed")
    model, selected_metadata = load_converged_actor(selected_checkpoint, device)
    selected_level = str(selected_metadata["budget"]["name"])
    final_metadata = {
        "status": "frozen",
        "d": d,
        "k": k,
        "selected_level": selected_level,
        "decision": decision,
        "architecture": {
            "width": ACTOR_WIDTH,
            "depth": ACTOR_DEPTH,
            "activation": "silu",
        },
        "nonterminal_oracle_reference_used_for_training": False,
        "nonterminal_oracle_reference_used_for_selection": False,
    }
    _atomic_torch_save(
        paths["final_checkpoint"], checkpoint_payload(model, final_metadata)
    )

    # The checkpoint is now frozen.  Analytic references first appear below.
    exact_diagnostic = _postselection_exact_diagnostic(
        equation, model, audit_t, audit_x, device=device
    )
    pre_evaluation_seconds = time.perf_counter() - cell_started
    rows, predictions, test_hash = _paired_stage_b_evaluation(
        equation,
        model,
        device=device,
        training_seconds=pre_evaluation_seconds,
    )
    _atomic_predictions(paths["predictions"], predictions)
    status = "complete" if decision["converged"] else "max_budget_reached"
    payload = {
        "schema_version": 1,
        "status": status,
        "d": d,
        "k": k,
        "device": str(device),
        "selected_level": selected_level,
        "budget_convergence": decision,
        "protocol": budget_protocol(),
        "audit": audit_metadata,
        "levels": records,
        "postselection_exact_diagnostic": exact_diagnostic,
        "paired_stage_b_seed_namespace": "B1|actor_mlp",
        "test_point_sha256": test_hash,
        "rows": rows,
        "elapsed_seconds": time.perf_counter() - cell_started,
        "nonterminal_oracle_reference_used_for_training": False,
        "nonterminal_oracle_reference_used_for_selection": False,
        "exact_reference_first_used_after_checkpoint_frozen": True,
    }
    atomic_json(paths["result"], payload)
    print(
        json.dumps(
            {
                "status": status,
                "d": d,
                "k": k,
                "selected_level": selected_level,
                "actor_relative_gradient_error": exact_diagnostic[
                    "actor_relative_gradient_error"
                ],
                "elapsed_seconds": payload["elapsed_seconds"],
            },
            indent=2,
        ),
        flush=True,
    )
    return payload


def _ordered_cells(pilot: bool) -> list[tuple[int, int]]:
    if pilot:
        return [(100, 10), (400, 40)]
    config = load_config()
    dimensions = [int(value) for value in config["pde"]["dimensions"]]
    effective = [int(value) for value in config["pde"]["effective_dimensions"]]
    pilot_cells = [(100, 10), (400, 40)]
    return pilot_cells + [
        (d, k)
        for d in dimensions
        for k in effective
        if (d, k) not in pilot_cells
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--d", type=int)
    parser.add_argument("--k", type=int)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if not (args.pilot or args.all or (args.d is not None and args.k is not None)):
        parser.error("use --pilot, --all, or both --d and --k")
    try:
        if args.d is not None and args.k is not None:
            cells = [(args.d, args.k)]
        else:
            cells = _ordered_cells(args.pilot)
        results = [run_cell(d, k, force=args.force) for d, k in cells]
        statuses = {f"d{item['d']}_k{item['k']}": item["status"] for item in results}
        atomic_json(
            CONVERGENCE_DIR / ("pilot_status.json" if args.pilot else "run_status.json"),
            {"schema_version": 1, "cells": statuses},
        )
        print(json.dumps({"cells": statuses}, indent=2), flush=True)
    except Exception as error:
        CONVERGENCE_DIR.mkdir(parents=True, exist_ok=True)
        atomic_json(
            CONVERGENCE_DIR / "failure.json",
            {
                "status": "failed",
                "error": repr(error),
                "traceback": traceback.format_exc(),
            },
        )
        raise


if __name__ == "__main__":
    main()

