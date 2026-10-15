"""Evaluate one benchmark/seed on the common validation and 960-state test bank."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time
import traceback
from typing import Any, Callable

import numpy as np
import torch

from .common import (
    PREDICTION_DIR,
    RAW_DIR,
    atomic_json,
    driver_for,
    ensure_dirs,
    evaluation_data,
    metrics,
    stable_seed,
    synchronize,
)
from .deep_picard import load_deep_picard
from .full_history_mlp import DEFAULT_CONFIGURATIONS, predict_full_history_mlp
from .neural_baselines import (
    load_deep_bsde,
    load_sdgd_pinn,
    predict_deep_bsde,
    predict_sdgd_pinn,
)
from .scasml import predict_scasml_level1
from .tangent_baseline import (
    evaluate_tangent_actor,
    load_tangent_actor,
    terminal_tangent_actor,
)


def _load_json(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        return json.load(handle)


def _scasml_prediction_key(checkpoint: str, paths: int, purpose: str) -> str:
    if purpose not in {"training_frontier", "inference_sweep"}:
        raise ValueError(purpose)
    return f"scasml_l1_{checkpoint}_p{paths}_{purpose}"


def _unique_picard_checkpoints(
    checkpoints: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Return one row per real saved network, retaining crossed budgets."""
    unique: dict[tuple[str, int], dict[str, Any]] = {}
    for checkpoint in checkpoints:
        key = (str(checkpoint["checkpoint"]), int(checkpoint["round"]))
        target = checkpoint.get("target_budget_seconds")
        if key not in unique:
            unique[key] = dict(
                checkpoint,
                crossed_target_budget_seconds=[] if target is None else [target],
            )
        elif target is not None:
            unique[key]["crossed_target_budget_seconds"].append(target)
    return list(unique.values())


def _metric_row(
    *,
    method: str,
    benchmark: str,
    seed: int,
    checkpoint: str,
    training_seconds: float,
    inference_seconds: float,
    prediction: np.ndarray,
    reference: np.ndarray,
    validation_objective: float | None,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    result = {
        "method": method,
        "benchmark": benchmark,
        "seed": seed,
        "checkpoint": checkpoint,
        "training_seconds": float(training_seconds),
        "inference_seconds_N960": float(inference_seconds),
        "total_seconds_N960": float(training_seconds + inference_seconds),
        "validation_objective": validation_objective,
        "final_test_used_for_checkpoint_selection": False,
        **metrics(prediction, reference),
    }
    if extra:
        result.update(extra)
    return result


def _direct_actor_prediction(
    actor: torch.nn.Module,
    benchmark: str,
    x: np.ndarray,
    device: torch.device,
) -> tuple[np.ndarray, float]:
    driver = driver_for(benchmark)
    state = torch.as_tensor(x, dtype=torch.float32, device=device)
    t = torch.zeros(len(state), dtype=torch.float32, device=device)
    synchronize(device)
    started = time.perf_counter()
    with torch.no_grad():
        value, _ = actor.predict(t, state, driver.T)
    synchronize(device)
    elapsed = time.perf_counter() - started
    return driver.working_to_original(value.cpu().numpy()).astype(np.float64), elapsed


def _tangent_evaluation(
    method: str,
    benchmark: str,
    seed: int,
    checkpoint_label: str,
    training_seconds: float,
    actor: torch.nn.Module,
    test_x: np.ndarray,
    test_reference: np.ndarray,
    device: torch.device,
    *,
    paths: int,
    validation_objective: float | None,
) -> tuple[dict[str, Any], np.ndarray]:
    prediction, variance, metadata = evaluate_tangent_actor(
        benchmark,
        actor,
        test_x,
        device,
        paths=paths,
        seed=stable_seed(
            "strong_tangent_inference", benchmark, seed, method, checkpoint_label, paths, "test"
        ),
    )
    row = _metric_row(
        method=method,
        benchmark=benchmark,
        seed=seed,
        checkpoint=checkpoint_label,
        training_seconds=training_seconds,
        inference_seconds=metadata["wall_seconds"],
        prediction=prediction,
        reference=test_reference,
        validation_objective=validation_objective,
        extra={
            "inference_paths": paths,
            "inference_time_steps": 32,
            "mc_mean_se": metadata["mean_value_mc_se"],
            "mc_rmse_standard_error": float(np.sqrt(np.mean(variance))),
            "weight_ess_fraction": metadata["weight_ess_fraction"],
            "nonfinite_count_inference": metadata["nonfinite_count"],
            "inference_kind": "tangent Feynman-Kac Monte Carlo",
        },
    )
    return row, prediction


def evaluate(benchmark: str, seed: int, device: torch.device) -> dict[str, Any]:
    data = evaluation_data(benchmark)
    validation_indices = data["validation_indices"]
    test_indices = data["test_indices"]
    validation_x = data["x"][validation_indices]
    validation_reference = data["reference_original"][validation_indices]
    test_x = data["x"][test_indices]
    test_reference = data["reference_original"][test_indices]
    training_rows: list[dict[str, Any]] = []
    total_rows: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    predictions: dict[str, np.ndarray] = {"test_reference": test_reference}

    # Continuous neural baselines: every real elapsed-time checkpoint.
    for method in ("deep_bsde", "sdgd_pinn"):
        raw_path = RAW_DIR / f"{method}_{benchmark}_s{seed}_main.json"
        result = _load_json(raw_path)
        for index, checkpoint in enumerate(result["checkpoints"]):
            path = Path(checkpoint["checkpoint"])
            if method == "deep_bsde":
                model, _ = load_deep_bsde(path, device)
                prediction, inference = predict_deep_bsde(model, driver_for(benchmark), test_x, device)
            else:
                model, _ = load_sdgd_pinn(path, device)
                prediction, inference = predict_sdgd_pinn(model, driver_for(benchmark), test_x, device)
            row = _metric_row(
                method=method,
                benchmark=benchmark,
                seed=seed,
                checkpoint=f"B{index + 1}",
                training_seconds=checkpoint["training_seconds"],
                inference_seconds=inference,
                prediction=prediction,
                reference=test_reference,
                validation_objective=checkpoint["validation_objective"],
                extra={
                    "target_budget_seconds": checkpoint["target_budget_seconds"],
                    "optimizer_step": checkpoint["step"],
                    "inference_kind": "amortized deterministic network",
                    "nonfinite_count_training": checkpoint["nonfinite_count"],
                },
            )
            row["prediction_key"] = f"{method}_B{index + 1}"
            training_rows.append(row)
            total_rows.append(dict(row, total_variant="direct"))
            predictions[f"{method}_B{index + 1}"] = prediction

            # SCaSML uses exactly this PINN as its charged backbone.
            if method == "sdgd_pinn":
                try:
                    sca_prediction, sca_variance, sca_meta = predict_scasml_level1(
                        model,
                        benchmark,
                        test_x,
                        device,
                        paths=32,
                        probes=4,
                        seed=stable_seed("scasml", benchmark, seed, index, 32, "test"),
                    )
                    sca_row = _metric_row(
                        method="scasml_l1",
                        benchmark=benchmark,
                        seed=seed,
                        checkpoint=f"PINN_B{index + 1}",
                        training_seconds=checkpoint["training_seconds"],
                        inference_seconds=sca_meta["inference_seconds"],
                        prediction=sca_prediction,
                        reference=test_reference,
                        validation_objective=checkpoint["validation_objective"],
                        extra={
                            "target_budget_seconds": checkpoint["target_budget_seconds"],
                            "inference_paths": 32,
                            "hutchinson_probes": 4,
                            "mc_mean_se": sca_meta["mean_mc_se_across_workload"],
                            "mc_rmse_standard_error": float(np.sqrt(np.mean(sca_variance))),
                            "inference_kind": "SCaSML first defect correction",
                            "nonfinite_count_inference": sca_meta["nonfinite_count"],
                        },
                    )
                    # Keep this common training-frontier evaluation distinct
                    # from the final-checkpoint inference sweep below.  The two
                    # calls intentionally use independent Monte Carlo seeds, so
                    # sharing an archive key would invalidate metric replay.
                    sca_key = _scasml_prediction_key(
                        f"B{index + 1}", 32, "training_frontier"
                    )
                    sca_row["prediction_key"] = sca_key
                    training_rows.append(sca_row)
                    total_rows.append(dict(sca_row, total_variant="paths32"))
                    predictions[sca_key] = sca_prediction
                except Exception as error:
                    failures.append(
                        {
                            "method": "scasml_l1",
                            "checkpoint": f"PINN_B{index + 1}",
                            "error": repr(error),
                            "traceback": traceback.format_exc(),
                        }
                    )

    # Deep Picard's real complete-round checkpoints.
    dpi = _load_json(RAW_DIR / f"deep_picard_{benchmark}_s{seed}_main.json")
    dpi_checkpoints = _unique_picard_checkpoints(list(dpi["checkpoints"]))
    final_history = dpi["history"][-1]
    if final_history["checkpoint"] and all(
        item["checkpoint"] != final_history["checkpoint"] for item in dpi_checkpoints
    ):
        dpi_checkpoints.append(
            {
                "checkpoint": final_history["checkpoint"],
                "training_seconds": final_history["cumulative_training_seconds"],
                "round": final_history["round"],
                "validation_objective": final_history["validation_normalized_loss"],
                "target_budget_seconds": None,
            }
        )
    for index, checkpoint in enumerate(dpi_checkpoints):
        actor, _ = load_deep_picard(Path(checkpoint["checkpoint"]), device)
        prediction, inference = _direct_actor_prediction(actor, benchmark, test_x, device)
        row = _metric_row(
            method="deep_picard",
            benchmark=benchmark,
            seed=seed,
            checkpoint=f"round{checkpoint['round']}",
            training_seconds=checkpoint["training_seconds"],
            inference_seconds=inference,
            prediction=prediction,
            reference=test_reference,
            validation_objective=checkpoint["validation_objective"],
            extra={
                "target_budget_seconds": checkpoint.get("target_budget_seconds"),
                "picard_round": checkpoint["round"],
                "inference_kind": "amortized deterministic network",
            },
        )
        row["prediction_key"] = f"deep_picard_round{checkpoint['round']}"
        training_rows.append(row)
        total_rows.append(dict(row, total_variant="direct"))
        predictions[f"deep_picard_round{checkpoint['round']}"] = prediction

    # Primary and standard Tangent Actors. A common 512-path test evaluator is
    # used for the training frontier; the final primary checkpoint also gets
    # validation-audited 128/512/2048 inference variants below.
    for method in ("tangent_actor", "standard_tangent_actor"):
        result = _load_json(RAW_DIR / f"{method}_{benchmark}_s{seed}_main.json")
        for phase in result["phases"]:
            actor, _ = load_tangent_actor(Path(phase["checkpoint"]), device)
            label = f"a{phase['iteration']}"
            row, prediction = _tangent_evaluation(
                method,
                benchmark,
                seed,
                label,
                phase["cumulative_training_seconds"],
                actor,
                test_x,
                test_reference,
                device,
                paths=512,
                validation_objective=phase["training"]["best_validation_normalized_loss"],
            )
            row.update(
                actor_iteration=phase["iteration"],
                actor_value_validation_rmse=phase["training"]["selected_validation_value_rmse"],
                actor_gradient_validation_vector_rmse=phase["training"][
                    "selected_validation_gradient_vector_rmse"
                ],
                label_weight_ess_fraction=phase["labels"]["train"]["weight_ess_fraction"],
                value_label_replication_rmse=phase["labels"]["value_label_replication_rmse"],
                gradient_label_replication_vector_rmse=phase["labels"][
                    "gradient_label_replication_vector_rmse"
                ],
            )
            row["prediction_key"] = f"{method}_{label}_p512"
            training_rows.append(row)
            total_rows.append(dict(row, total_variant="paths512"))
            predictions[f"{method}_{label}_p512"] = prediction

    # Untrained A_g internal ablation; no optimizer or label cost.
    ag = terminal_tangent_actor(benchmark, device)
    ag_row, ag_prediction = _tangent_evaluation(
        "terminal_tangent_Ag",
        benchmark,
        seed,
        "untrained",
        0.0,
        ag,
        test_x,
        test_reference,
        device,
        paths=512,
        validation_objective=None,
    )
    training_rows.append(ag_row)
    ag_row["prediction_key"] = "terminal_tangent_Ag_p512"
    total_rows.append(dict(ag_row, total_variant="paths512"))
    predictions["terminal_tangent_Ag_p512"] = ag_prediction

    # Multiple inference budgets for the final predeclared primary Actor.
    primary = _load_json(RAW_DIR / f"tangent_actor_{benchmark}_s{seed}_main.json")
    last_phase = primary["phases"][-1]
    final_actor, _ = load_tangent_actor(Path(last_phase["checkpoint"]), device)
    for paths in (128, 512, 2048):
        validation_prediction, validation_variance, validation_meta = evaluate_tangent_actor(
            benchmark,
            final_actor,
            validation_x,
            device,
            paths=paths,
            seed=stable_seed("strong_tangent_inference", benchmark, seed, "a3", paths, "validation"),
        )
        if paths == 512:
            existing = next(
                row
                for row in total_rows
                if row["method"] == "tangent_actor" and row["checkpoint"] == "a3"
            )
            variant = dict(existing)
            test_prediction = predictions["tangent_actor_a3_p512"]
        else:
            variant, test_prediction = _tangent_evaluation(
                "tangent_actor",
                benchmark,
                seed,
                "a3",
                last_phase["cumulative_training_seconds"],
                final_actor,
                test_x,
                test_reference,
                device,
                paths=paths,
                validation_objective=last_phase["training"]["best_validation_normalized_loss"],
            )
        variant.update(
            total_variant=f"paths{paths}",
            validation_rmse=metrics(validation_prediction, validation_reference)["rmse"],
            validation_relative_l2=metrics(validation_prediction, validation_reference)["relative_l2"],
            validation_mc_mean_se=validation_meta["mean_value_mc_se"],
            validation_inference_seconds=validation_meta["wall_seconds"],
        )
        variant["prediction_key"] = f"tangent_actor_a3_p{paths}"
        # Replace the generic 512 row with its validation-audited version.
        total_rows = [
            row
            for row in total_rows
            if not (
                row["method"] == "tangent_actor"
                and row["checkpoint"] == "a3"
                and row.get("total_variant") == f"paths{paths}"
            )
        ]
        total_rows.append(variant)
        predictions[f"tangent_actor_a3_p{paths}"] = test_prediction

    # SCaSML final checkpoint inference sweep (validation and test are timed separately).
    pinn_result = _load_json(RAW_DIR / f"sdgd_pinn_{benchmark}_s{seed}_main.json")
    last_pinn = pinn_result["checkpoints"][-1]
    pinn_model, _ = load_sdgd_pinn(Path(last_pinn["checkpoint"]), device)
    for paths in (8, 32, 128):
        try:
            val_prediction, _, val_meta = predict_scasml_level1(
                pinn_model,
                benchmark,
                validation_x,
                device,
                paths=paths,
                probes=4,
                seed=stable_seed("scasml", benchmark, seed, "final", paths, "validation"),
            )
            test_prediction, test_variance, test_meta = predict_scasml_level1(
                pinn_model,
                benchmark,
                test_x,
                device,
                paths=paths,
                probes=4,
                seed=stable_seed("scasml", benchmark, seed, "final", paths, "test"),
            )
            row = _metric_row(
                method="scasml_l1",
                benchmark=benchmark,
                seed=seed,
                checkpoint="PINN_B5",
                training_seconds=last_pinn["training_seconds"],
                inference_seconds=test_meta["inference_seconds"],
                prediction=test_prediction,
                reference=test_reference,
                validation_objective=last_pinn["validation_objective"],
                extra={
                    "total_variant": f"paths{paths}",
                    "inference_paths": paths,
                    "hutchinson_probes": 4,
                    "mc_mean_se": test_meta["mean_mc_se_across_workload"],
                    "mc_rmse_standard_error": float(np.sqrt(np.mean(test_variance))),
                    "validation_rmse": metrics(val_prediction, validation_reference)["rmse"],
                    "validation_relative_l2": metrics(val_prediction, validation_reference)[
                        "relative_l2"
                    ],
                    "validation_mc_mean_se": val_meta["mean_mc_se_across_workload"],
                    "validation_inference_seconds": val_meta["inference_seconds"],
                    "inference_kind": "SCaSML first defect correction",
                },
            )
            sca_key = _scasml_prediction_key("B5", paths, "inference_sweep")
            row["prediction_key"] = sca_key
            total_rows = [
                current
                for current in total_rows
                if not (
                    current["method"] == "scasml_l1"
                    and current["checkpoint"] == "PINN_B5"
                    and current.get("total_variant") == f"paths{paths}"
                )
            ]
            total_rows.append(row)
            predictions[sca_key] = test_prediction
        except Exception as error:
            failures.append(
                {
                    "method": "scasml_l1",
                    "checkpoint": "PINN_B5",
                    "paths": paths,
                    "error": repr(error),
                    "traceback": traceback.format_exc(),
                }
            )

    # Inference-only Full-History MLP configurations. Stop before the scheduler
    # wall limit if cumulative configuration time becomes excessive.
    mlp_cumulative = 0.0
    for n, M in DEFAULT_CONFIGURATIONS:
        try:
            prediction, metadata = predict_full_history_mlp(
                benchmark, test_x, device, n=n, M=M, seed=seed
            )
            mlp_cumulative += metadata["inference_seconds"]
            row = _metric_row(
                method="full_history_mlp",
                benchmark=benchmark,
                seed=seed,
                checkpoint=f"n{n}_M{M}",
                training_seconds=0.0,
                inference_seconds=metadata["inference_seconds"],
                prediction=prediction,
                reference=test_reference,
                validation_objective=None,
                extra={
                    "total_variant": f"n{n}_M{M}",
                    "mlp_level": n,
                    "mlp_sample_base": M,
                    "terminal_samples": metadata["stats"]["terminal_samples"],
                    "transition_samples": metadata["stats"]["transition_samples"],
                    "value_projection_count": metadata["stats"]["value_projection_count"],
                    "inference_kind": "inference-only recursive multilevel Picard",
                },
            )
            row["prediction_key"] = f"full_history_n{n}_M{M}"
            total_rows.append(row)
            predictions[f"full_history_n{n}_M{M}"] = prediction
            if mlp_cumulative > 2400.0:
                failures.append(
                    {
                        "method": "full_history_mlp",
                        "checkpoint": "remaining configurations",
                        "error": "predeclared 2400-second cumulative evaluation guard reached",
                    }
                )
                break
        except Exception as error:
            failures.append(
                {
                    "method": "full_history_mlp",
                    "checkpoint": f"n{n}_M{M}",
                    "error": repr(error),
                    "traceback": traceback.format_exc(),
                }
            )

    prediction_path = PREDICTION_DIR / f"frontier_{benchmark}_s{seed}.npz"
    np.savez_compressed(prediction_path, **predictions)
    payload = {
        "benchmark": benchmark,
        "seed": seed,
        "training_frontier": training_rows,
        "total_compute_frontier": total_rows,
        "failures": failures,
        "prediction_path": str(prediction_path),
        "final_test_points": len(test_x),
        "validation_points": len(validation_x),
        "final_test_used_for_training_or_checkpoint_selection": False,
    }
    atomic_json(RAW_DIR / f"evaluation_{benchmark}_s{seed}.json", payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--benchmark",
        choices=("nonseparable_d100_k2", "joint_concave_d400_k2"),
        required=True,
    )
    parser.add_argument("--seed", type=int, choices=(0, 1), required=True)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    ensure_dirs()
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA allocation required")
    payload = evaluate(args.benchmark, args.seed, device)
    print(
        json.dumps(
            {
                "status": "complete",
                "benchmark": args.benchmark,
                "seed": args.seed,
                "training_rows": len(payload["training_frontier"]),
                "total_rows": len(payload["total_compute_frontier"]),
                "failures": len(payload["failures"]),
            },
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
