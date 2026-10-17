"""Post-selection Ground Truth, Euler, and MCV evaluation."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import math
import time
from typing import Any

import numpy as np
import torch

from .common import (
    CHECKPOINT_DIR,
    PREDICTION_DIR,
    RAW_DIR,
    REFERENCE_DIR,
    RESULTS_DIR,
    array_sha256,
    atomic_json,
    atomic_npz,
    ensure_dirs,
    environment_payload,
    metrics,
    stable_seed,
    state_bank,
)
from .configuration import ExperimentSpec, fit_config, protocol, updates
from .drivers import SineQuadraticDriver
from .reference2d import load_reference
from .solver import actor_prediction, evaluate_policy, initial_actor, load_actor


def _reference_path(spec: ExperimentSpec):
    # The reduced reference is independent of ambient d. Store one d=100 copy
    # for each (lambda,eta) PDE and reuse it only after Actor selection.
    reference_driver = SineQuadraticDriver(
        d=100, lam=spec.lam, eta=spec.eta if spec.coupled else 0.0
    )
    return REFERENCE_DIR / f"{reference_driver.name}.npz"


def _load_actors(
    spec: ExperimentSpec, device: torch.device
) -> tuple[list[torch.nn.Module], Any, Any]:
    frozen_protocol = protocol()
    config = fit_config(frozen_protocol, spec.budget)
    driver = spec.driver()
    actors: list[torch.nn.Module] = [
        initial_actor(driver, config, seed=spec.seed + 1701, device=device)
    ]
    for iteration in range(1, updates(frozen_protocol, spec.budget) + 1):
        path = CHECKPOINT_DIR / f"{spec.run_id}_{spec.estimator}_a{iteration}.pt"
        if not path.exists():
            raise FileNotFoundError(path)
        actors.append(load_actor(path, driver, config, device))
    return actors, driver, config


def _actor_truth_metrics(
    driver: SineQuadraticDriver,
    actor: torch.nn.Module,
    reference,
    t: np.ndarray,
    x: np.ndarray,
    device: torch.device,
) -> dict[str, float]:
    predicted_r, predicted_q = actor_prediction(driver, actor, t, x, device)
    coordinates = x @ driver.B
    exact_r = reference.value(t, coordinates)
    exact_reduced_gradient = reference.gradient(t, coordinates)
    exact_q = math.sqrt(2.0) * (exact_reduced_gradient @ driver.B.T)
    r_error = predicted_r - exact_r
    q_error = predicted_q - exact_q
    return {
        "actor_value_rmse": float(np.sqrt(np.mean(np.square(r_error)))),
        "actor_value_relative_l2": float(
            np.linalg.norm(r_error)
            / max(np.linalg.norm(exact_r), np.finfo(float).tiny)
        ),
        "actor_gradient_component_rmse": float(
            np.sqrt(np.mean(np.square(q_error)))
        ),
        "actor_gradient_vector_rmse": float(
            np.sqrt(np.mean(np.square(q_error).sum(axis=1)))
        ),
        "actor_gradient_relative_l2": float(
            np.linalg.norm(q_error)
            / max(np.linalg.norm(exact_q), np.finfo(float).tiny)
        ),
        "actor_gradient_norm_mean": float(
            np.mean(np.linalg.norm(predicted_q, axis=1))
        ),
        "truth_gradient_norm_mean": float(np.mean(np.linalg.norm(exact_q, axis=1))),
        "actor_value_min": float(np.min(predicted_r)),
        "actor_value_max": float(np.max(predicted_r)),
    }


def _policy_row(
    spec: ExperimentSpec,
    actor_name: str,
    prediction: dict[str, np.ndarray],
    details: dict[str, Any],
    truth: np.ndarray,
) -> dict[str, Any]:
    plain = metrics(prediction["plain"], truth)
    mcv = metrics(prediction["mcv"], truth)
    pairs = int(details["antithetic_pairs"])
    plain_mc_mse = float(np.mean(prediction["plain_pair_variance"] / pairs))
    mcv_mc_mse = float(np.mean(prediction["mcv_pair_variance"] / pairs))
    return {
        "run_id": spec.run_id,
        "family": spec.family,
        "d": spec.d,
        "lambda": spec.lam,
        "eta": spec.eta,
        "seed": spec.seed,
        "budget": spec.budget,
        "training_estimator": spec.estimator,
        "iteration": actor_name,
        "truth_value_mean": float(np.mean(truth)),
        "plain_value_mean": float(np.mean(prediction["plain"])),
        "plain_value_rmse": plain["rmse"],
        "plain_relative_l2": plain["relative_l2"],
        "plain_bias": plain["bias"],
        "plain_systematic_rmse": math.sqrt(max(plain["rmse"] ** 2 - plain_mc_mse, 0.0)),
        "mcv_value_mean": float(np.mean(prediction["mcv"])),
        "mcv_value_rmse": mcv["rmse"],
        "mcv_relative_l2": mcv["relative_l2"],
        "mcv_bias": mcv["bias"],
        "mcv_systematic_rmse": math.sqrt(max(mcv["rmse"] ** 2 - mcv_mc_mse, 0.0)),
        "plain_mc_rmse_contribution": math.sqrt(plain_mc_mse),
        "mcv_mc_rmse_contribution": math.sqrt(mcv_mc_mse),
        "mcv_variance_reduction": details["variance_reduction"],
        "mean_plain_mc_se": details["mean_plain_mc_se"],
        "mean_mcv_mc_se": details["mean_mcv_mc_se"],
        "weight_ess_fraction": details["weight_ess_fraction"],
        "weight_min": details["weight_min"],
        "weight_max": details["weight_max"],
        "paths": details["paths"],
        "time_steps": details["time_steps"],
        "inference_seconds": details["wall_seconds"],
        "nonfinite_count": details["nonfinite_count"],
    }


def evaluate(spec: ExperimentSpec, device: torch.device) -> dict[str, Any]:
    ensure_dirs()
    selection_path = RESULTS_DIR / f"selection_{spec.run_id}.json"
    if not selection_path.exists():
        raise RuntimeError("reference-free selection must be frozen before Ground Truth")
    selection = json.loads(selection_path.read_text(encoding="utf-8"))
    if selection.get("ground_truth_loaded"):
        raise RuntimeError("selection boundary metadata was unexpectedly modified")
    reference_path = _reference_path(spec)
    if not reference_path.exists():
        raise FileNotFoundError(
            f"independent reference must be validated before training: {reference_path}"
        )
    # This is intentionally the first reference access in the post-freeze process.
    reference = load_reference(reference_path)
    actors, driver, _ = _load_actors(spec, device)
    frozen_protocol = protocol()
    evaluation = frozen_protocol["post_freeze_evaluation"]
    test_t, test_x = state_bank(
        driver.name,
        driver.d,
        "post_freeze_common_test",
        int(evaluation["states"]),
        T=driver.T,
        seed=20261017,
    )
    truth = reference.value(test_t, test_x @ driver.B)
    audit_t, audit_x = state_bank(
        driver.name,
        driver.d,
        "post_freeze_actor_audit",
        int(evaluation["actor_audit_states"]),
        T=driver.T,
        seed=20261017,
    )
    rows: list[dict[str, Any]] = []
    predictions: dict[str, np.ndarray] = {
        "test_t": test_t,
        "test_x": test_x,
        "truth": truth,
    }
    common_seed = stable_seed(spec.run_id, "post_freeze_common_policy_paths")
    started = time.perf_counter()
    for index, actor in enumerate(actors):
        actor_name = f"a{index}"
        output, details = evaluate_policy(
            driver,
            actor,
            test_t,
            test_x,
            paths=int(evaluation["paths"]),
            time_steps=int(evaluation["steps"]),
            seed=common_seed,
            device=device,
            parent_batch=int(evaluation["parent_batch"]),
            path_chunk=int(evaluation["path_chunk"]),
        )
        row = _policy_row(spec, actor_name, output, details, truth)
        row.update(
            _actor_truth_metrics(
                driver, actor, reference, audit_t, audit_x, device
            )
        )
        rows.append(row)
        predictions[f"{actor_name}_plain"] = output["plain"]
        predictions[f"{actor_name}_mcv"] = output["mcv"]
        predictions[f"{actor_name}_plain_pair_variance"] = output[
            "plain_pair_variance"
        ]
        predictions[f"{actor_name}_mcv_pair_variance"] = output[
            "mcv_pair_variance"
        ]
        print(
            f"{spec.run_id} {actor_name}: RMSE={row['mcv_value_rmse']:.6g} "
            f"ActorZ={row['actor_gradient_vector_rmse']:.6g} "
            f"VR={row['mcv_variance_reduction']:.2f}x",
            flush=True,
        )

    selected_index = int(selection["selected_iteration"])
    selected_actor = actors[selected_index]
    euler_rows: list[dict[str, Any]] = []
    for steps in evaluation["euler_steps"]:
        output, details = evaluate_policy(
            driver,
            selected_actor,
            test_t,
            test_x,
            paths=int(evaluation["paths"]),
            time_steps=int(steps),
            seed=stable_seed(spec.run_id, "euler_common", steps),
            device=device,
            parent_batch=int(evaluation["parent_batch"]),
            path_chunk=int(evaluation["path_chunk"]),
        )
        row = _policy_row(
            spec, selection["selected_actor"], output, details, truth
        )
        row["evaluation_kind"] = "selected_actor_euler_sensitivity"
        euler_rows.append(row)
    prediction_path = PREDICTION_DIR / f"post_freeze_{spec.run_id}.npz"
    atomic_npz(prediction_path, **predictions)
    preflight = json.loads(
        (RESULTS_DIR / "preflight_audit.json").read_text(encoding="utf-8")
    )
    reference_key = "coupled" if spec.coupled else "primary"
    reference_audit = preflight.get("drivers", {}).get(reference_key)
    selected_row = rows[selected_index]
    reference_proxy = (
        None
        if reference_audit is None
        else reference_audit["reference"]["reference_error_proxy_value_rmse"]
    )
    result = {
        "spec": asdict(spec),
        "selection": selection,
        "iteration_rows": rows,
        "euler_rows": euler_rows,
        "selected_row": selected_row,
        "reference_path": str(reference_path),
        "reference_error_proxy_value_rmse": reference_proxy,
        "reference_proxy_over_selected_rmse": None
        if reference_proxy is None
        else reference_proxy / max(selected_row["mcv_value_rmse"], 1e-300),
        "test_bank_sha256": array_sha256(test_t, test_x),
        "prediction_path": str(prediction_path),
        "ground_truth_loaded_only_after_selection_file": True,
        "ground_truth_used_for_training_or_selection": False,
        "wall_seconds": time.perf_counter() - started,
    }
    atomic_json(RESULTS_DIR / f"post_freeze_{spec.run_id}.json", result)
    return result


def _spec_from_args(args: argparse.Namespace) -> ExperimentSpec:
    return ExperimentSpec(
        family=args.family,
        d=args.d,
        lam=args.lam,
        eta=args.eta if args.family == "coupled" else 0.0,
        seed=args.seed,
        budget=args.budget,
        estimator=args.estimator,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--family", choices=("primary", "coupled"), default="primary")
    parser.add_argument("--d", type=int, default=100)
    parser.add_argument("--lam", type=float, default=0.2)
    parser.add_argument("--eta", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--budget", choices=("pilot", "main"), default="pilot")
    parser.add_argument("--estimator", choices=("mcv", "plain"), default="mcv")
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA allocation is required")
    spec = _spec_from_args(args)
    atomic_json(
        RAW_DIR / f"environment_evaluation_{spec.run_id}.json",
        environment_payload("post_freeze_evaluation"),
    )
    evaluate(spec, device)


if __name__ == "__main__":
    main()
