"""Audit the completed Stage B grid before paper-facing aggregation."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from .common import RESULTS_DIR, atomic_json, load_config
from .runner import METHODS


FINAL_DPI_PROFILE = "published-hjb-architecture-resource-adapted"


def validate() -> dict[str, Any]:
    config = load_config()
    raw = RESULTS_DIR / "raw" / "stage_b"
    failures: list[str] = []
    cells_checked = 0
    rows_checked = 0
    hashes: dict[str, str] = {}
    for d in config["pde"]["dimensions"]:
        for k in config["pde"]["effective_dimensions"]:
            cell_hashes: set[str] = set()
            for method in METHODS:
                path = raw / f"{method}_d{d}_k{k}.json"
                if not path.is_file():
                    failures.append(f"missing {path.name}")
                    continue
                with path.open() as handle:
                    payload = json.load(handle)
                cells_checked += 1
                if payload.get("status") != "complete":
                    failures.append(f"{path.name}: status={payload.get('status')}")
                expected_rows = int(
                    config["stage_b"][
                        "mlp_seeds"
                        if method in {"vanilla_mlp", "double_mlp", "actor_mlp"}
                        else "neural_seeds"
                    ]
                )
                rows = payload.get("rows", [])
                if len(rows) != expected_rows:
                    failures.append(
                        f"{path.name}: rows={len(rows)}, expected={expected_rows}"
                    )
                for row in rows:
                    rows_checked += 1
                    cell_hashes.add(str(row.get("test_point_sha256")))
                    if int(row.get("test_points", -1)) != 960:
                        failures.append(f"{path.name}: test_points is not 960")
                    if row.get("test_distribution") != "x~N(0,I_d)":
                        failures.append(f"{path.name}: wrong test distribution")
                    if int(row.get("value_nonfinite", -1)) != 0:
                        failures.append(f"{path.name}: non-finite value prediction")
                    if int(row.get("gradient_nonfinite", -1)) != 0:
                        failures.append(f"{path.name}: non-finite gradient prediction")
                    for key in ("value_rmse", "value_bias", "gradient_rmse"):
                        if not math.isfinite(float(row[key])):
                            failures.append(f"{path.name}: non-finite {key}")

                shared = payload.get("shared", {})
                if method == "actor_mlp":
                    actor = shared.get("actor", {})
                    for flag in (
                        "analytic_reference_used_for_training",
                        "analytic_reference_used_for_checkpoint_selection",
                    ):
                        if actor.get(flag) is not False:
                            failures.append(f"{path.name}: actor {flag} is not false")
                elif method in {"dpi", "pinn_sdgd", "deep_bsde"}:
                    for seed in range(expected_rows):
                        metadata = shared.get(f"seed_{seed}", {})
                        for flag in (
                            "analytic_reference_used_for_training",
                            "analytic_reference_used_for_checkpoint_selection",
                        ):
                            if metadata.get(flag) is not False:
                                failures.append(
                                    f"{path.name}: seed {seed} {flag} is not false"
                                )
                        if method == "dpi" and metadata.get(
                            "training_config", {}
                        ).get("profile") != FINAL_DPI_PROFILE:
                            failures.append(
                                f"{path.name}: seed {seed} is not final DPI profile"
                            )

            if len(cell_hashes) != 1:
                failures.append(
                    f"d={d},k={k}: methods have {len(cell_hashes)} test hashes"
                )
            elif cell_hashes:
                hashes[f"d{d}_k{k}"] = next(iter(cell_hashes))

    expected_cells = 15 * len(METHODS)
    if cells_checked != expected_cells:
        failures.append(
            f"completed method-cells={cells_checked}, expected={expected_cells}"
        )
    payload = {
        "schema_version": 1,
        "status": "PASS" if not failures else "FAIL",
        "method_cells_checked": cells_checked,
        "expected_method_cells": expected_cells,
        "seed_rows_checked": rows_checked,
        "test_point_hashes": hashes,
        "failures": failures,
    }
    atomic_json(RESULTS_DIR / "stage_b_validation.json", payload)
    if failures:
        raise RuntimeError("Stage B validation failed:\n" + "\n".join(failures))
    return payload


def main() -> None:
    print(json.dumps(validate(), indent=2))


if __name__ == "__main__":
    main()
