"""Audit saved results, fairness invariants, and reproducibility metadata."""

from __future__ import annotations

from collections import defaultdict
import json
import platform
import subprocess
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.drivers import (
    make_driver,
)

from .common import RESULTS_DIR, array_sha256, atomic_json
from .models import TerminalAwareActor
from .run_scaling import load_protocol, stage_specs


def _terminal_identity() -> dict[str, Any]:
    rows = []
    for family in ("allen_cahn", "joint_concave", "joint_convex", "nonseparable"):
        driver = make_driver(family, d=8, k=2)
        actor = TerminalAwareActor(driver, width=16, depth=2).double()
        generator = torch.Generator().manual_seed(991)
        x = torch.randn(64, driver.d, generator=generator, dtype=torch.float64)
        t = torch.full((len(x),), driver.T, dtype=torch.float64)
        r, q = actor.predict(t, x, driver.T)
        expected_q = driver.sigma * driver.terminal_gradient_torch(x)
        rows.append(
            {
                "family": family,
                "value_max_abs_error": float(
                    (r - driver.terminal(x)).detach().abs().max().item()
                ),
                "gradient_max_abs_error": None
                if q is None
                else float((q - expected_q).detach().abs().max().item()),
                "value_min": float(r.detach().min().item()),
                "value_max": float(r.detach().max().item()),
                "domain_low": driver.r_low,
                "domain_high": driver.r_high,
            }
        )
    return {
        "rows": rows,
        "all_exact_to_float64_tolerance": all(
            row["value_max_abs_error"] < 5e-14
            and (row["gradient_max_abs_error"] is None or row["gradient_max_abs_error"] < 5e-14)
            for row in rows
        ),
    }


def main() -> None:
    protocol = load_protocol()
    payloads = []
    for path in sorted(RESULTS_DIR.glob("run_*.json")):
        with path.open() as handle:
            payloads.append(json.load(handle))
    completed = {payload["iteration_rows"][0]["run_id"] for payload in payloads}
    expected_stage_a = {
        spec.run_id for spec in stage_specs("stage-a", protocol, selected_variant=None)
    }
    missing_stage_a = sorted(expected_stage_a - completed)
    selection_path = RESULTS_DIR / "stage_b_selection.json"
    if selection_path.exists():
        with selection_path.open() as handle:
            selected_variant = json.load(handle)["selected_variant"]
        expected_transfer = {
            spec.run_id
            for stage in ("stage-b", "stage-c")
            for spec in stage_specs(stage, protocol, selected_variant)
        }
    else:
        expected_transfer = set()
    missing_transfer = sorted(expected_transfer - completed)
    errors: list[str] = []
    fairness: dict[tuple[str, int], dict[str, set[Any]]] = defaultdict(
        lambda: defaultdict(set)
    )
    run_records = []
    for payload in payloads:
        run_id = payload["iteration_rows"][0]["run_id"]
        rows = payload["iteration_rows"]
        if [row["iteration"] for row in rows] != ["a0", "a1", "a2", "a3"]:
            errors.append(f"{run_id}: incomplete iteration sequence")
        if payload["selection_used_exact_test"]:
            errors.append(f"{run_id}: exact test used for selection")
        if payload["projected_actor_input_used"]:
            errors.append(f"{run_id}: projected Actor input")
        for row in rows:
            if row["evaluation_nonfinite_count"] != 0:
                errors.append(f"{run_id} {row['iteration']}: nonfinite rollout")
            ess = row["evaluation_weight_ess_fraction"]
            if not (0.0 < ess <= 1.000001):
                errors.append(f"{run_id} {row['iteration']}: invalid ESS {ess}")
        for phase in payload["phases"]:
            trace = phase["trace"]
            selected = min(trace, key=lambda item: item["validation_normalized_loss"])
            if selected["step"] != phase["training"]["selected_step"]:
                errors.append(f"{run_id} a{phase['iteration']}: checkpoint mismatch")
        initial = payload["phases"][0]["labels"]
        group = fairness[(payload["spec"]["problem_id"] if "problem_id" in payload["spec"] else f"{payload['spec']['family']}_d{payload['spec']['d']}_k{payload['spec']['k']}", payload["spec"]["training_seed"])]
        group["full_train_bank_sha256"].add(initial["train_bank_sha256"])
        group["validation_bank_sha256"].add(initial["validation_bank_sha256"])
        group["label_random_seed"].add(initial["label_random_seed"])
        run_records.append(
            {
                "run_id": run_id,
                "spec": payload["spec"],
                "variant": payload["variant"],
                "test_bank_sha256": payload["test_bank_sha256"],
                "phase_train_bank_sha256": [
                    phase["labels"]["used_train_bank_sha256"]
                    for phase in payload["phases"]
                ],
                "phase_validation_bank_sha256": [
                    phase["labels"]["validation_bank_sha256"]
                    for phase in payload["phases"]
                ],
                "phase_label_random_seeds": [
                    phase["labels"]["label_random_seed"]
                    for phase in payload["phases"]
                ],
                "checkpoint_training_seeds": [
                    # The deterministic seed is reconstructible from run_id and iteration;
                    # its role is stated explicitly in the runner and protocol.
                    f"stable_seed(scaling_actor,{run_id},{phase['iteration']})"
                    for phase in payload["phases"]
                ],
            }
        )
    fairness_rows = []
    for (problem, seed), values in sorted(fairness.items()):
        row = {
            "problem": problem,
            "training_seed": seed,
            **{key: sorted(value) for key, value in values.items()},
            "shared_initial_bank_and_mc_stream": all(len(value) == 1 for value in values.values()),
        }
        fairness_rows.append(row)
        if not row["shared_initial_bank_and_mc_stream"]:
            errors.append(f"{problem} seed {seed}: unmatched a1 labels")
    terminal = _terminal_identity()
    if not terminal["all_exact_to_float64_tolerance"]:
        errors.append("terminal identities failed")
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        commit = "unknown"
    reproducibility = {
        "protocol": protocol,
        "current_git_commit": commit,
        "python": platform.python_version(),
        "torch": torch.__version__,
        "cuda": torch.version.cuda,
        "completed_run_ids": sorted(completed),
        "run_records": run_records,
        "initial_label_fairness": fairness_rows,
        "terminal_identity": terminal,
        "sample_bank_scheme": "count-independent deterministic streams; smaller training banks are exact prefixes",
        "validation_and_policy_evaluation_splits": "independent named deterministic streams",
    }
    integrity = {
        "passed_for_completed_runs": not errors,
        "stage_a_complete": not missing_stage_a,
        "missing_stage_a_run_ids": missing_stage_a,
        "planned_transfer_complete": not missing_transfer,
        "missing_transfer_run_ids": missing_transfer,
        "errors": errors,
        "completed_run_count": len(payloads),
        "terminal_identity_passed": terminal["all_exact_to_float64_tolerance"],
        "all_initial_comparisons_matched": all(
            row["shared_initial_bank_and_mc_stream"] for row in fairness_rows
        ),
    }
    atomic_json(RESULTS_DIR / "reproducibility_metadata.json", reproducibility)
    atomic_json(RESULTS_DIR / "result_integrity.json", integrity)
    print(json.dumps(integrity, indent=2))
    if errors:
        raise RuntimeError(errors)


if __name__ == "__main__":
    main()

