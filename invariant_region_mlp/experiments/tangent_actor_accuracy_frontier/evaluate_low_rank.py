"""Post-freeze audit for the learned low-rank Gradient-Actor ablation."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import time

import numpy as np
import torch

from .common import (
    CHECKPOINT_DIR,
    OLD_REFERENCE_DIR,
    PREDICTION_DIR,
    RAW_DIR,
    RESULTS_DIR,
    array_sha256,
    atomic_json,
    atomic_npz,
    environment_payload,
    stable_seed,
    state_bank,
)
from .configuration import ExperimentSpec, fit_config, protocol, updates
from .evaluate_postfreeze import _actor_truth_metrics, _policy_row
from .low_rank_actor import load_low_rank_actor, make_low_rank_actor
from .reference2d import load_reference
from .solver import evaluate_policy


def _run_id(spec: ExperimentSpec, rank: int) -> str:
    return f"{spec.run_id}_lowrank_r{rank}"


def run(spec: ExperimentSpec, rank: int, device: torch.device) -> dict:
    run_id = _run_id(spec, rank)
    selection = json.loads(
        (RESULTS_DIR / f"selection_{run_id}.json").read_text(encoding="utf-8")
    )
    driver = spec.driver()
    config = fit_config(protocol(), "main")
    actors = [make_low_rank_actor(driver, config, rank, seed=spec.seed + 1701, device=device)]
    for iteration in range(1, updates(protocol(), "main") + 1):
        actors.append(load_low_rank_actor(
            CHECKPOINT_DIR / f"{run_id}_mcv_a{iteration}.pt", driver, config, device
        ))
    evaluation = protocol()["post_freeze_evaluation"]
    test_t, test_x = state_bank(
        driver.name, driver.d, "post_freeze_common_test", int(evaluation["states"]),
        T=driver.T, seed=20261017,
    )
    audit_t, audit_x = state_bank(
        driver.name, driver.d, "post_freeze_actor_audit", int(evaluation["actor_audit_states"]),
        T=driver.T, seed=20261017,
    )
    reference_driver = type(driver)(d=100, lam=spec.lam, eta=spec.eta if spec.coupled else 0.0)
    reference_path = OLD_REFERENCE_DIR / f"{reference_driver.name}.npz"
    reference = load_reference(reference_path)
    truth = reference.value(test_t, test_x @ driver.B)
    common_seed = stable_seed(run_id, "post_freeze_common_policy_paths")
    rows = []
    predictions = {"test_t": test_t, "test_x": test_x, "truth": truth}
    started = time.perf_counter()
    for index, actor in enumerate(actors):
        actor_name = f"a{index}"
        output, details = evaluate_policy(
            driver, actor, test_t, test_x,
            paths=int(evaluation["paths"]), time_steps=int(evaluation["steps"]),
            seed=common_seed, device=device,
            parent_batch=int(evaluation["parent_batch"]),
            path_chunk=int(evaluation["path_chunk"]),
        )
        row = _policy_row(spec, actor_name, output, details, truth)
        row.update(_actor_truth_metrics(driver, actor, reference, audit_t, audit_x, device))
        row.update({"run_id": run_id, "architecture": "learned_low_rank_gradient_actor", "rank": rank})
        rows.append(row)
        predictions[f"{actor_name}_mcv"] = output["mcv"]
        predictions[f"{actor_name}_plain"] = output["plain"]
        print(f"{run_id} {actor_name}: RMSE={row['mcv_value_rmse']:.6g} Z={row['actor_gradient_vector_rmse']:.6g}", flush=True)
    selected_index = int(selection["selected_iteration"])
    euler_rows = []
    for steps in evaluation["euler_steps"]:
        output, details = evaluate_policy(
            driver, actors[selected_index], test_t, test_x,
            paths=int(evaluation["paths"]), time_steps=int(steps),
            seed=stable_seed(run_id, "euler_common", steps), device=device,
            parent_batch=int(evaluation["parent_batch"]),
            path_chunk=int(evaluation["path_chunk"]),
        )
        row = _policy_row(spec, selection["selected_actor"], output, details, truth)
        row.update({"run_id": run_id, "architecture": "learned_low_rank_gradient_actor", "rank": rank})
        euler_rows.append(row)
    prediction_path = PREDICTION_DIR / f"post_freeze_{run_id}.npz"
    atomic_npz(prediction_path, **predictions)
    result = {
        "spec": asdict(spec),
        "run_id": run_id,
        "architecture": "learned_low_rank_gradient_actor",
        "rank": rank,
        "selection": selection,
        "iteration_rows": rows,
        "euler_rows": euler_rows,
        "selected_row": rows[selected_index],
        "reference_path": str(reference_path),
        "test_bank_sha256": array_sha256(test_t, test_x),
        "prediction_path": str(prediction_path),
        "ground_truth_loaded_only_after_selection_file": True,
        "ground_truth_used_for_training_or_selection": False,
        "wall_seconds": time.perf_counter() - started,
    }
    atomic_json(RESULTS_DIR / f"post_freeze_{run_id}.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--family", choices=("primary", "coupled"), default="primary")
    parser.add_argument("--d", type=int, default=400)
    parser.add_argument("--lam", type=float, default=0.2)
    parser.add_argument("--eta", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--rank", type=int, choices=(8, 16, 32), required=True)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    spec = ExperimentSpec(
        family=args.family, d=args.d, lam=args.lam,
        eta=args.eta if args.family == "coupled" else 0.0,
        seed=args.seed, budget="main", estimator="mcv",
    )
    run_id = _run_id(spec, args.rank)
    atomic_json(
        RAW_DIR / f"environment_evaluation_{run_id}.json",
        environment_payload("low_rank_postfreeze_evaluation"),
    )
    run(spec, args.rank, torch.device(args.device))


if __name__ == "__main__":
    main()
