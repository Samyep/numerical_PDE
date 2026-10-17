"""Train and select the learned low-rank Gradient-Actor ablation."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import time

import torch

from .common import (
    CHECKPOINT_DIR,
    OLD_RESULTS_DIR,
    RAW_DIR,
    RESULTS_DIR,
    atomic_json,
    ensure_dirs,
    environment_payload,
    stable_seed,
)
from .configuration import ExperimentSpec, fit_config, protocol, updates
from .low_rank_actor import (
    load_low_rank_actor,
    make_low_rank_actor,
    save_low_rank_actor,
    train_low_rank_actor,
)
from .solver import fixed_point_diagnostic, generate_phase_labels


def _run_id(spec: ExperimentSpec, rank: int) -> str:
    return f"{spec.run_id}_lowrank_r{rank}"


def run(spec: ExperimentSpec, rank: int, device: torch.device) -> dict:
    ensure_dirs()
    preflight = json.loads((OLD_RESULTS_DIR / "preflight_audit.json").read_text())
    key = "coupled" if spec.coupled else "primary"
    if not preflight["drivers"][key]["passed"]:
        raise RuntimeError("parent preflight failed")
    frozen_protocol = protocol()
    config = fit_config(frozen_protocol, spec.budget)
    run_id = _run_id(spec, rank)
    source = make_low_rank_actor(
        spec.driver(), config, rank, seed=spec.seed + 1701, device=device
    )
    driver = spec.driver()
    actors = [source]
    phases = []
    started = time.perf_counter()
    for iteration in range(1, updates(frozen_protocol, spec.budget) + 1):
        arrays, label_metadata = generate_phase_labels(
            run_id, driver, source, iteration, config,
            training_seed=spec.seed, estimator=spec.estimator, device=device,
        )
        checkpoint = CHECKPOINT_DIR / f"{run_id}_mcv_a{iteration}.pt"
        metadata_path = CHECKPOINT_DIR / f"{run_id}_mcv_a{iteration}.json"
        if checkpoint.exists() and metadata_path.exists():
            actor = load_low_rank_actor(checkpoint, driver, config, device)
            phase = json.loads(metadata_path.read_text(encoding="utf-8"))
        else:
            actor, training, trace = train_low_rank_actor(
                driver, config, rank, arrays,
                seed=stable_seed(run_id, iteration, "actor_fit", spec.seed),
                device=device, initial=source,
            )
            phase = {
                "run_id": run_id,
                "iteration": iteration,
                "source": f"a{iteration - 1}",
                "estimator": spec.estimator,
                "labels": label_metadata,
                "training": training,
                "trace": trace,
                "checkpoint_selection": "minimum held-out noisy tangent-label loss",
                "ground_truth_used": False,
            }
            save_low_rank_actor(
                checkpoint, actor, driver=driver, config=config, metadata=phase
            )
            atomic_json(metadata_path, phase)
        actors.append(actor)
        phases.append(phase)
        source = actor
        print(
            f"{run_id} a{iteration}: fit={phase['training']['best_validation_normalized_loss']:.6g} "
            f"step={phase['training']['selected_step']}", flush=True,
        )

    fp = frozen_protocol["fixed_point"]
    diagnostics = {"primary": [], "ranking_audit": []}
    for tag, paths, steps in (
        ("primary", int(fp["main_paths"]), int(fp["main_steps"])),
        ("ranking_audit", int(fp["ranking_audit_paths"]), int(fp["ranking_audit_steps"])),
    ):
        for index, actor in enumerate(actors):
            path = RAW_DIR / f"fixed_point_{run_id}_a{index}_{tag}.json"
            if path.exists():
                row = json.loads(path.read_text(encoding="utf-8"))
            else:
                row = fixed_point_diagnostic(
                    run_id, f"a{index}", driver, actor,
                    states=int(fp["states"]), paths=paths, time_steps=steps,
                    training_seed=spec.seed, device=device,
                    value_scale_floor=float(fp["value_scale_floor"]),
                    gradient_component_scale_floor=float(fp["gradient_component_scale_floor"]),
                    audit_tag=tag,
                )
                atomic_json(path, row)
            diagnostics[tag].append(row)
    preliminary = min(
        (row for row in diagnostics["primary"] if row["selection_eligible"]),
        key=lambda row: row["normalized_total_inconsistency"],
    )
    selected = min(
        (row for row in diagnostics["ranking_audit"] if row["selection_eligible"]),
        key=lambda row: row["normalized_total_inconsistency"],
    )
    minimum_ess = min(row["minimum_weight_ess_fraction"] for row in diagnostics["ranking_audit"])
    a0_score = diagnostics["ranking_audit"][0]["normalized_total_inconsistency"]
    reduction = (a0_score - selected["normalized_total_inconsistency"]) / max(a0_score, 1e-300)
    selection = {
        "spec": asdict(spec),
        "architecture": "learned_low_rank_gradient_actor",
        "rank": rank,
        "run_id": run_id,
        "selected_actor": selected["actor"],
        "selected_iteration": int(selected["actor"][1:]),
        "selected_by": "minimum independently replicated noise-corrected normalized fixed-point inconsistency",
        "preliminary_actor": preliminary["actor"],
        "ranking_stable_when_paths_and_steps_doubled": selected["actor"] == preliminary["actor"],
        "fixed_point_reduction_vs_a0": reduction,
        "minimum_weight_ess_fraction": minimum_ess,
        "ground_truth_loaded": False,
        "primary_diagnostics": diagnostics["primary"],
        "ranking_audit_diagnostics": diagnostics["ranking_audit"],
        "oracle_subspace_used": False,
    }
    atomic_json(RESULTS_DIR / f"selection_{run_id}.json", selection)
    result = {
        "method": "terminal_aware_tangent_actor_mcv_learned_low_rank_ablation",
        "spec": asdict(spec),
        "run_id": run_id,
        "rank": rank,
        "fit_config": asdict(config),
        "phases": phases,
        "selection": selection,
        "full_ambient_actor_input": True,
        "projection_matrix_learned_end_to_end": True,
        "oracle_projection_used": False,
        "ground_truth_used_for_training_or_selection": False,
        "wall_seconds_this_invocation": time.perf_counter() - started,
    }
    atomic_json(RESULTS_DIR / f"training_{run_id}.json", result)
    print(f"SELECTION {run_id}: {selected['actor']}", flush=True)
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
        RAW_DIR / f"environment_training_{run_id}.json",
        environment_payload("low_rank_training"),
    )
    run(spec, args.rank, torch.device(args.device))


if __name__ == "__main__":
    main()
