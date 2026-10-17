"""Train undamped Tangent-Actor iterations and freeze reference-free selection."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import time
from typing import Any

import torch

from .common import RAW_DIR, RESULTS_DIR, atomic_json, ensure_dirs, environment_payload
from .configuration import ExperimentSpec, fit_config, protocol, updates
from .solver import (
    fit_phase,
    fixed_point_diagnostic,
    generate_phase_labels,
    initial_actor,
)


def _preflight_gate(spec: ExperimentSpec) -> dict[str, Any]:
    path = RESULTS_DIR / "preflight_audit.json"
    if not path.exists():
        raise RuntimeError("Stage A preflight audit is missing")
    payload = json.loads(path.read_text(encoding="utf-8"))
    name = "coupled" if spec.coupled else "primary"
    if spec.lam == 0.2 and name in payload["drivers"]:
        audit = payload["drivers"][name]
        if not audit["passed"]:
            raise RuntimeError(f"Stage A gate failed for {name}")
        return audit
    # Additional curvature-strength references are written by the expansion
    # preflight before this entry point is called.
    path = RESULTS_DIR / f"preflight_{spec.run_id.split('_' + spec.budget)[0]}.json"
    if not path.exists():
        raise RuntimeError(f"missing frozen expansion preflight: {path}")
    audit = json.loads(path.read_text(encoding="utf-8"))
    if not audit["passed"]:
        raise RuntimeError("expansion preflight failed")
    return audit


def _cached_diagnostic(
    spec: ExperimentSpec,
    actor_name: str,
    driver,
    actor,
    *,
    tag: str,
    states: int,
    paths: int,
    steps: int,
    device: torch.device,
    fp: dict[str, Any],
) -> dict[str, Any]:
    path = RAW_DIR / f"fixed_point_{spec.run_id}_{actor_name}_{tag}.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    result = fixed_point_diagnostic(
        spec.run_id,
        actor_name,
        driver,
        actor,
        states=states,
        paths=paths,
        time_steps=steps,
        training_seed=spec.seed,
        device=device,
        value_scale_floor=float(fp["value_scale_floor"]),
        gradient_component_scale_floor=float(
            fp["gradient_component_scale_floor"]
        ),
        audit_tag=tag,
    )
    atomic_json(path, result)
    return result


def train_and_select(spec: ExperimentSpec, device: torch.device) -> dict[str, Any]:
    ensure_dirs()
    frozen_protocol = protocol()
    preflight = _preflight_gate(spec)
    config = fit_config(frozen_protocol, spec.budget)
    number_updates = updates(frozen_protocol, spec.budget)
    driver = spec.driver()
    started = time.perf_counter()
    source = initial_actor(
        driver,
        config,
        seed=spec.seed + 1701,
        device=device,
    )
    actors = [source]
    phases: list[dict[str, Any]] = []
    for iteration in range(1, number_updates + 1):
        phase_started = time.perf_counter()
        arrays, label_metadata = generate_phase_labels(
            spec.run_id,
            driver,
            source,
            iteration,
            config,
            training_seed=spec.seed,
            estimator=spec.estimator,
            device=device,
        )
        actor, phase = fit_phase(
            spec.run_id,
            driver,
            source,
            iteration,
            config,
            arrays,
            label_metadata,
            training_seed=spec.seed,
            estimator=spec.estimator,
            device=device,
        )
        phase["phase_wall_seconds_this_invocation"] = time.perf_counter() - phase_started
        actors.append(actor)
        phases.append(phase)
        source = actor
        print(
            f"{spec.run_id} a{iteration}: "
            f"fit={phase['training']['best_validation_normalized_loss']:.6g} "
            f"step={phase['training']['selected_step']} "
            f"ESS={label_metadata['train']['weight_ess_fraction']:.4f}",
            flush=True,
        )

    fp = frozen_protocol["fixed_point"]
    diagnostic_paths = int(
        fp["pilot_paths"] if spec.budget == "pilot" else fp["main_paths"]
    )
    diagnostic_steps = int(
        fp["pilot_steps"] if spec.budget == "pilot" else fp["main_steps"]
    )
    primary_diagnostics = []
    for index, actor in enumerate(actors):
        result = _cached_diagnostic(
            spec,
            f"a{index}",
            driver,
            actor,
            tag="primary",
            states=int(fp["states"]),
            paths=diagnostic_paths,
            steps=diagnostic_steps,
            device=device,
            fp=fp,
        )
        primary_diagnostics.append(result)
        print(
            f"{spec.run_id} a{index} FP={result['normalized_total_inconsistency']:.6g} "
            f"ESS={result['minimum_weight_ess_fraction']:.4f}",
            flush=True,
        )
    eligible_primary = [
        row for row in primary_diagnostics if row["selection_eligible"]
    ]
    if not eligible_primary:
        raise RuntimeError("no finite Actor is eligible for reference-free selection")
    preliminary = min(
        eligible_primary, key=lambda row: row["normalized_total_inconsistency"]
    )

    ranking_diagnostics = []
    for index, actor in enumerate(actors):
        ranking_diagnostics.append(
            _cached_diagnostic(
                spec,
                f"a{index}",
                driver,
                actor,
                tag="ranking_audit",
                states=int(fp["states"]),
                paths=int(fp["ranking_audit_paths"]),
                steps=int(fp["ranking_audit_steps"]),
                device=device,
                fp=fp,
            )
        )
    eligible_audit = [
        row for row in ranking_diagnostics if row["selection_eligible"]
    ]
    selected = min(
        eligible_audit, key=lambda row: row["normalized_total_inconsistency"]
    )
    ranking_stable = selected["actor"] == preliminary["actor"]
    a0_score = ranking_diagnostics[0]["normalized_total_inconsistency"]
    reduction = (a0_score - selected["normalized_total_inconsistency"]) / max(
        a0_score, 1e-300
    )
    gate = frozen_protocol["pilot_gate"]
    minimum_ess = min(
        row["minimum_weight_ess_fraction"] for row in ranking_diagnostics
    )
    passed = bool(
        all(row["nonfinite_count"] == 0 for row in ranking_diagnostics)
        and minimum_ess >= float(gate["minimum_weight_ess_fraction"])
        and reduction >= float(gate["minimum_fixed_point_reduction_fraction"])
        and ranking_stable
    )
    # This file is the irreversible model-selection boundary. Post-freeze
    # Ground Truth evaluation refuses to run without it.
    selection_path = RESULTS_DIR / f"selection_{spec.run_id}.json"
    selection = {
        "spec": asdict(spec),
        "selected_actor": selected["actor"],
        "selected_iteration": int(selected["actor"][1:]),
        "selected_by": "minimum noise-corrected normalized fixed-point inconsistency",
        "preliminary_actor": preliminary["actor"],
        "ranking_stable_when_paths_and_steps_doubled": ranking_stable,
        "fixed_point_reduction_vs_a0": reduction,
        "minimum_weight_ess_fraction": minimum_ess,
        "passed_reference_free_gate": passed,
        "ground_truth_loaded": False,
        "primary_diagnostics": primary_diagnostics,
        "ranking_audit_diagnostics": ranking_diagnostics,
    }
    atomic_json(selection_path, selection)
    result = {
        "method": "terminal_aware_tangent_actor_mcv",
        "spec": asdict(spec),
        "driver": driver.metadata(),
        "fit_config": asdict(config),
        "updates": number_updates,
        "phases": phases,
        "selection": selection,
        "selection_path": str(selection_path),
        "preflight_passed": bool(preflight["passed"]),
        "full_ambient_actor_input": True,
        "oracle_projection_used": False,
        "ground_truth_used_for_training_or_selection": False,
        "undamped_updates": True,
        "wall_seconds_this_invocation": time.perf_counter() - started,
    }
    atomic_json(RESULTS_DIR / f"training_{spec.run_id}.json", result)
    if (
        spec.budget == "pilot"
        and spec.family == "primary"
        and spec.d == 100
        and spec.lam == 0.2
        and spec.estimator == "mcv"
        and spec.seed == 0
    ):
        atomic_json(
            RESULTS_DIR / "pilot_decision.json",
            {
                "passed": passed,
                "basis": "reference-free fixed-point consistency, ranking stability, ESS, and finiteness",
                "run_id": spec.run_id,
                "selected_actor": selected["actor"],
                "fixed_point_reduction_vs_a0": reduction,
                "minimum_weight_ess_fraction": minimum_ess,
                "ground_truth_used": False,
            },
        )
    print(
        f"SELECTION {spec.run_id}: {selected['actor']} gate={'PASS' if passed else 'FAIL'} "
        f"reduction={reduction:.3%} stable={ranking_stable}",
        flush=True,
    )
    return result


def _spec_from_args(args: argparse.Namespace) -> ExperimentSpec:
    eta = args.eta if args.family == "coupled" else 0.0
    return ExperimentSpec(
        family=args.family,
        d=args.d,
        lam=args.lam,
        eta=eta,
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
    ensure_dirs()
    atomic_json(
        RAW_DIR / f"environment_training_{_spec_from_args(args).run_id}.json",
        environment_payload("training"),
    )
    train_and_select(_spec_from_args(args), device)


if __name__ == "__main__":
    main()
