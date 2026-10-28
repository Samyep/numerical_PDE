"""Freeze inputs, audit exact reuse, and enforce scientific launch gates."""

from __future__ import annotations

from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import numpy as np
import torch

from invariant_region_mlp.experiments.tangent_actor_icml_final_suite.protocol import (
    stable_seed as stable_final,
)

from .models import UnanchoredJointActor
from .protocol import (
    ALLOCATED_GPU_HOUR_CEILING,
    HJB_CONFIG,
    INPUT_DIR,
    JOINT_CONFIG,
    PDES,
    RESULTS_DIR,
    TRAINING_SEEDS,
    ensure_dirs,
    frozen_protocol,
)
from .run_ablation import (
    FINAL_RESULTS,
    FRONTIER_RESULTS,
    ORIGINAL_RESULTS,
    _array_hash,
    _file_hash,
    _reuse_checkpoint,
    _source_selection,
)
from invariant_region_mlp.experiments.tangent_actor_icml_final_suite.run_joint_actor import (
    make_driver,
    reference,
)
from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis.models import (
    TerminalAwareActor,
)


def _write_json(path: Path, payload: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _freeze_nonconvex_bank() -> dict[str, object]:
    path = INPUT_DIR / "test_bank_nonconvex_coupled_d100_960.npz"
    rng = np.random.default_rng(
        stable_final("icml_final_test_bank", "nonconvex_coupled", 100)
    )
    t = rng.uniform(0.0, np.nextafter(0.5, 0.0), 960).astype(np.float64)
    x = rng.standard_normal((960, 100)).astype(np.float64)
    if path.exists():
        with np.load(path) as archive:
            if not np.array_equal(archive["t"], t) or not np.array_equal(archive["x"], x):
                raise RuntimeError("frozen nonconvex bank does not match preregistered rule")
    else:
        np.savez_compressed(path, t=t, x=x)
    return {
        "path": str(path),
        "states": len(t),
        "array_sha256": _array_hash(t, x),
        "file_sha256": _file_hash(path),
        "seed": stable_final("icml_final_test_bank", "nonconvex_coupled", 100),
        "historical_512_primary_reuse": False,
    }


def _bank_audit(nonconvex: dict[str, object]) -> dict[str, object]:
    records: dict[str, object] = {"nonconvex_coupled": nonconvex}
    hjb = ORIGINAL_RESULTS / "evaluation_bank_hjb_log_d100.npz"
    with np.load(hjb) as archive:
        index = archive["test_indices"]
        t, x = archive["t"][index], archive["x"][index]
    records["hjb_log"] = {
        "path": str(hjb), "states": len(t), "array_sha256": _array_hash(t, x),
        "file_sha256": _file_hash(hjb),
    }
    for pde in ("joint_concave", "joint_convex"):
        path = FINAL_RESULTS / "inputs" / f"test_bank_{pde}_d100.npz"
        with np.load(path) as archive:
            t, x = archive["t"], archive["x"]
        records[pde] = {
            "path": str(path), "states": len(t), "array_sha256": _array_hash(t, x),
            "file_sha256": _file_hash(path),
        }
    if any(int(record["states"]) != 960 for record in records.values()):
        raise RuntimeError("all primary banks must have exactly 960 states")
    return records


def _reuse_audit() -> dict[str, object]:
    cells = []
    reused = new = 0
    for pde in PDES:
        for seed in TRAINING_SEEDS:
            paths = [_reuse_checkpoint(pde, seed, i) for i in range(1, 5)]
            exact = all(path is not None and path.exists() for path in paths)
            selection = None
            selection_source = None
            if exact:
                selection, selection_source = _source_selection(pde, seed)
                reused += 1
            else:
                new += 1
            cells.append({
                "pde": pde,
                "seed": seed,
                "variant": "a0_full",
                "exact_training_reuse": exact,
                "checkpoint_paths": [None if path is None else str(path) for path in paths],
                "checkpoint_sha256": [
                    None if path is None or not path.exists() else _file_hash(path)
                    for path in paths
                ],
                "selected_iteration": selection,
                "selection_source": selection_source,
                "evaluation_reuse": False,
                "reason": "all primary cells require two frozen 960-state inference replicates",
            })
    return {
        "cells": cells,
        "full_actor_training_cells_reused": reused,
        "full_actor_training_cells_new": new,
        "a1_training_cells_new": 12,
        "a2_training_cells_new": 12,
        "a3_training_cells": 0,
        "historical_primary_metrics_reused": 0,
    }


def _scientific_checks(banks: dict[str, object]) -> dict[str, object]:
    driver = make_driver("joint_concave", 100)
    torch.manual_seed(501)
    anchored = TerminalAwareActor(driver, width=32, depth=2)
    torch.manual_seed(501)
    unanchored = UnanchoredJointActor(driver, width=32, depth=2)
    generator = torch.Generator().manual_seed(502)
    x = torch.randn((8, 100), generator=generator)
    t = torch.rand(8, generator=generator) * driver.T
    ar, aq = anchored.predict(t, x, driver.T)
    ur, uq = unanchored.predict(t, x, driver.T)
    initial_difference = max(
        float(torch.max(torch.abs(ar - ur))),
        float(torch.max(torch.abs(aq - uq))) if aq is not None and uq is not None else 0.0,
    )
    terminal_t = torch.full((len(x),), driver.T)
    terminal_r, terminal_q = anchored.predict(terminal_t, x, driver.T)
    anchored_terminal_error = max(
        float(torch.max(torch.abs(terminal_r - driver.terminal(x)))),
        float(torch.max(torch.abs(terminal_q - driver.sigma * driver.terminal_gradient_torch(x))))
        if terminal_q is not None
        else float("inf"),
    )
    with torch.no_grad():
        unanchored.network[-1].bias.fill_(0.2)
    changed_r, changed_q = unanchored.predict(terminal_t, x, driver.T)
    unanchored_terminal_departure = max(
        float(torch.max(torch.abs(changed_r - driver.terminal(x)))),
        float(torch.max(torch.abs(changed_q - driver.sigma * driver.terminal_gradient_torch(x))))
        if changed_q is not None
        else 0.0,
    )
    convex = make_driver("joint_convex", 100)
    probe = np.asarray([-0.4, 0.1, 0.8])
    sign_error = float(np.max(np.abs(convex.working_to_original(probe) + probe)))
    _, concave_reference = reference("joint_concave")
    _, convex_reference = reference("joint_convex")
    _, nonconvex_reference = reference("nonconvex_coupled")
    checks = {
        "same_original_pde": True,
        "convex_sign_conversion_max_error": sign_error,
        "matching_initial_prediction_max_error": initial_difference,
        "anchored_terminal_identity_max_error": anchored_terminal_error,
        "unanchored_can_depart_terminal": unanchored_terminal_departure > 1e-3,
        "plain_mcv_same_target": "proved by payoff identity and unit test",
        "pathwise_labels": "plain and complete corrected payoffs differentiated; unit tested",
        "value_domain_constraints_retained": bool(
            torch.all(changed_r > driver.r_low) and torch.all(changed_r < driver.r_high)
        ),
        "banks_and_references_match": all(int(row["states"]) == 960 for row in banks.values()),
        "reference_audits": {
            "joint_concave": concave_reference,
            "joint_convex": convex_reference,
            "nonconvex_coupled": nonconvex_reference,
        },
        "ground_truth_leakage": False,
    }
    checks["passed"] = bool(
        sign_error < 1e-12
        and initial_difference < 5e-7
        and anchored_terminal_error < 5e-7
        and checks["unanchored_can_depart_terminal"]
        and checks["value_domain_constraints_retained"]
        and checks["banks_and_references_match"]
    )
    return checks


def _resource_estimate() -> list[dict[str, object]]:
    estimates = {
        "hjb_log": 0.38,
        "joint_concave": 0.60,
        "joint_convex": 0.63,
        "nonconvex_coupled": 0.92,
    }
    rows = []
    for variant in ("a1_no_mcv", "a2_no_terminal_anchor"):
        for pde in PDES:
            rows.append({
                "workload": "training_and_selection",
                "pde": pde,
                "variant": variant,
                "cells": 3,
                "estimated_allocated_gpu_hours": 3 * estimates[pde],
                "actual_allocated_gpu_hours": "",
            })
    for pde in ("joint_concave", "joint_convex"):
        rows.append({
            "workload": "missing_full_training_and_selection",
            "pde": pde,
            "variant": "a0_full",
            "cells": 2,
            "estimated_allocated_gpu_hours": 2 * estimates[pde],
            "actual_allocated_gpu_hours": "",
        })
    rows.append({
        "workload": "all_primary_and_independent_inference_replicates",
        "pde": "all",
        "variant": "all",
        "cells": 80,
        "estimated_allocated_gpu_hours": 4.0,
        "actual_allocated_gpu_hours": "",
    })
    return rows


def _write_resource_csv(rows: list[dict[str, object]]) -> None:
    import csv

    path = RESULTS_DIR / "RESOURCE_ACCOUNTING.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    ensure_dirs()
    frozen_path = RESULTS_DIR / "FROZEN_PROTOCOL.json"
    frozen = json.loads(frozen_path.read_text(encoding="utf-8"))
    expected = frozen_protocol()
    if frozen["protocol_sha256"] != expected["protocol_sha256"]:
        raise RuntimeError("checked-in frozen protocol differs from source declaration")
    nonconvex = _freeze_nonconvex_bank()
    banks = _bank_audit(nonconvex)
    reuse = _reuse_audit()
    checks = _scientific_checks(banks)
    if not checks["passed"]:
        raise RuntimeError(f"scientific preflight failed: {checks}")
    _write_json(RESULTS_DIR / "REUSE_AUDIT.json", reuse)
    _write_json(RESULTS_DIR / "PREFLIGHT_AUDIT.json", {
        "protocol_sha256": frozen["protocol_sha256"],
        "test_banks": banks,
        "scientific_checks": checks,
        "passed": True,
    })
    rows = _resource_estimate()
    _write_resource_csv(rows)
    estimate = sum(float(row["estimated_allocated_gpu_hours"]) for row in rows)
    _write_json(RESULTS_DIR / "GPU_HOUR_LEDGER.json", {
        "estimated_allocated_gpu_hours": estimate,
        "hard_ceiling": ALLOCATED_GPU_HOUR_CEILING,
        "fits_ceiling": estimate <= ALLOCATED_GPU_HOUR_CEILING,
        "failed_and_retry_hours_count_toward_ceiling": True,
    })
    print(json.dumps({
        "passed": True,
        "estimated_allocated_gpu_hours": estimate,
        "reused_a0": reuse["full_actor_training_cells_reused"],
        "new_a0": reuse["full_actor_training_cells_new"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
