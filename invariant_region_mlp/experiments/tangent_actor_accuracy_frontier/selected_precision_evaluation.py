"""High-step CRN audit for Actors trained inside the frontier study."""

from __future__ import annotations

import argparse
import json
import time

import numpy as np
import torch

from .common import (
    CHECKPOINT_DIR,
    OLD_REFERENCE_DIR,
    OLD_RESULTS_DIR,
    PREDICTION_DIR,
    RAW_DIR,
    RESULTS_DIR,
    array_sha256,
    atomic_json,
    atomic_npz,
    ensure_dirs,
    environment_payload,
    stable_seed,
)
from .configuration import ExperimentSpec, fit_config, protocol
from .precision_evaluation import _evaluate_level, _row
from .reference2d import load_reference
from .solver import initial_actor, load_actor


def run(
    spec: ExperimentSpec,
    *,
    paths: int,
    levels: tuple[int, ...],
    fine_steps: int,
    parent_batch: int,
    path_chunk: int,
    bootstrap_replicates: int,
    device: torch.device,
) -> dict:
    ensure_dirs()
    selection_path = RESULTS_DIR / f"selection_{spec.run_id}.json"
    selection = json.loads(selection_path.read_text(encoding="utf-8"))
    selected_index = int(selection["selected_iteration"])
    driver = spec.driver()
    config = fit_config(protocol(), spec.budget)
    if selected_index == 0:
        actor = initial_actor(driver, config, seed=spec.seed + 1701, device=device)
        checkpoint = None
    else:
        checkpoint = CHECKPOINT_DIR / (
            f"{spec.run_id}_{spec.estimator}_a{selected_index}.pt"
        )
        actor = load_actor(checkpoint, driver, config, device)
    bank_path = PREDICTION_DIR / f"post_freeze_{spec.run_id}.npz"
    with np.load(bank_path) as bank:
        t = bank["test_t"].astype(np.float64)
        x = bank["test_x"].astype(np.float64)
    reference_driver = type(driver)(
        d=100, lam=spec.lam, eta=spec.eta if spec.coupled else 0.0
    )
    reference_path = OLD_REFERENCE_DIR / f"{reference_driver.name}.npz"
    reference = load_reference(reference_path)
    truth = reference.value(t, x @ driver.B)
    preflight = json.loads(
        (OLD_RESULTS_DIR / "preflight_audit.json").read_text(encoding="utf-8")
    )
    key = "coupled" if spec.coupled else "primary"
    reference_uncertainty = float(
        preflight["drivers"][key]["reference"]["reference_error_proxy_value_rmse"]
    )
    common_seed = stable_seed("frontier_selected_precision", spec.run_id, paths, fine_steps)
    rows, outputs = [], {}
    started = time.perf_counter()
    for steps in levels:
        output, details = _evaluate_level(
            driver, actor, t, x, paths=paths, steps=steps, fine_steps=fine_steps,
            seed=common_seed, parent_batch=parent_batch, path_chunk=path_chunk,
            device=device,
        )
        row = _row(
            output, truth, details,
            seed=stable_seed("selected_precision_bootstrap", spec.run_id, steps),
            reference_uncertainty=reference_uncertainty,
            bootstrap_replicates=bootstrap_replicates,
        )
        row.update({
            "run_id": spec.run_id,
            "training_seed": spec.seed,
            "selected_actor": selection["selected_actor"],
        })
        rows.append(row)
        outputs[steps] = output
        print(
            f"{spec.run_id} {selection['selected_actor']} steps={steps} "
            f"RMSE={row['value_rmse']:.7g} MC={row['mc_rmse_contribution']:.7g}",
            flush=True,
        )
    prediction_path = PREDICTION_DIR / (
        f"selected_precision_{spec.run_id}_paths{paths}_fine{fine_steps}.npz"
    )
    arrays = {"test_t": t, "test_x": x, "truth": truth}
    for steps, output in outputs.items():
        arrays[f"mcv_{steps}"] = output["mcv"]
        arrays[f"plain_{steps}"] = output["plain"]
        arrays[f"mcv_pair_variance_{steps}"] = output["mcv_pair_variance"]
    atomic_npz(prediction_path, **arrays)
    result = {
        "spec": spec.__dict__,
        "selected_actor": selection["selected_actor"],
        "selection_path": str(selection_path),
        "checkpoint": None if checkpoint is None else str(checkpoint),
        "reference_path": str(reference_path),
        "reference_uncertainty_value_rmse": reference_uncertainty,
        "test_bank_sha256": array_sha256(t, x),
        "common_random_numbers_across_levels": True,
        "rows": rows,
        "prediction_path": str(prediction_path),
        "total_wall_seconds": time.perf_counter() - started,
        "ground_truth_used_for_training_or_selection": False,
    }
    atomic_json(RAW_DIR / f"selected_precision_{spec.run_id}.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--family", choices=("primary", "coupled"), required=True)
    parser.add_argument("--d", type=int, required=True)
    parser.add_argument("--lam", type=float, default=0.2)
    parser.add_argument("--eta", type=float, default=1.0)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--paths", type=int, default=1024)
    parser.add_argument("--levels", type=int, nargs="+", default=[64, 128, 256])
    parser.add_argument("--fine-steps", type=int, default=256)
    parser.add_argument("--parent-batch", type=int, default=32)
    parser.add_argument("--path-chunk", type=int, default=128)
    parser.add_argument("--bootstrap-replicates", type=int, default=3000)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    spec = ExperimentSpec(
        family=args.family, d=args.d, lam=args.lam,
        eta=args.eta if args.family == "coupled" else 0.0,
        seed=args.seed, budget="main", estimator="mcv",
    )
    atomic_json(
        RAW_DIR / f"environment_selected_precision_{spec.run_id}.json",
        environment_payload("selected_precision_evaluation"),
    )
    run(
        spec, paths=args.paths, levels=tuple(args.levels), fine_steps=args.fine_steps,
        parent_batch=args.parent_batch, path_chunk=args.path_chunk,
        bootstrap_replicates=args.bootstrap_replicates,
        device=torch.device(args.device),
    )


if __name__ == "__main__":
    main()
