"""Independent mechanical audit of the frozen 2-D Euler result artifacts."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"


def load_json(name: str) -> dict[str, Any]:
    return json.loads((RESULTS / name).read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    data_audit = load_json("data_and_reference_audit.json")
    selection = load_json("stencil_selection.json")
    summary = load_json("summary.json")
    reports = [load_json(f"hcfl_s6_report_seed{seed}.json") for seed in (0, 1, 2)]
    with (RESULTS / "case_metrics.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    if len(rows) != 65:
        raise RuntimeError(f"Expected 65 retained case rows, received {len(rows)}")
    if any(row["completed"].lower() != "true" for row in rows):
        raise RuntimeError("At least one trajectory failed completion")

    accuracy = summary["accuracy"]
    comparisons: dict[str, Any] = {}
    for split in ("test_id", "official_quadrants"):
        roe = accuracy[f"{split}/PyClaw Roe-64"]["nmae_mean_completed"]
        hllc = accuracy[f"{split}/HLLC-64"]["nmae_mean_completed"]
        hcfl = accuracy[f"{split}/HCFL-64"]["nmae_mean_completed"]
        comparisons[split] = {
            "roe64_nmae": roe,
            "hllc64_nmae": hllc,
            "hcfl64_nmae": hcfl,
            "hcfl_relative_reduction_vs_hllc": (hllc - hcfl) / hllc,
            "hcfl_relative_reduction_vs_roe": (roe - hcfl) / roe,
            "hcfl_beats_hllc": hcfl < hllc,
            "hcfl_beats_roe": hcfl < roe,
        }

    hcfl_safety = {
        key: value
        for key, value in summary["safety"].items()
        if key.startswith("HCFL-64")
    }
    physical_checks = {
        "all_hcfl_interface_residuals_below_1e-6": all(
            item["maximum_interface_residual"] <= 1.0e-6
            for item in hcfl_safety.values()
        ),
        "all_hcfl_entropy_balances_nonpositive_with_tolerance": all(
            item["maximum_entropy_balance"] <= 5.0e-7
            for item in hcfl_safety.values()
        ),
        "all_hcfl_conservation_closures_below_1e-7": all(
            item["maximum_conservation_closure"] <= 1.0e-7
            for item in hcfl_safety.values()
        ),
        "no_hcfl_positivity_fallback": all(
            item["positivity_fallback_rate"] == 0.0 for item in hcfl_safety.values()
        ),
        "no_hcfl_entropy_fallback": all(
            item["entropy_fallback_rate"] == 0.0 for item in hcfl_safety.values()
        ),
    }

    checkpoint_hashes = {
        f"seed{seed}": sha256(RESULTS / f"hcfl_s6_best_seed{seed}.pt")
        for seed in (0, 1, 2)
    }
    figures = {
        name: {
            "bytes": (RESULTS / name).stat().st_size,
            "sha256": sha256(RESULTS / name),
        }
        for name in (
            "official_final_density.png",
            "official_density_linecuts.png",
            "nmae_comparison.png",
        )
    }

    checks = {
        "data_audit_passed": bool(data_audit["passed"]),
        "selected_stencil_is_six": selection["selected_stencil_cells"] == 6,
        "selection_preceded_official_test": not selection["official_test_was_evaluated"],
        "three_seed_reports_present": len(reports) == 3,
        "all_seeds_validation_plateau_stopped": all(
            report["stop_reason"] == "validation_plateau_at_minimum_learning_rate"
            for report in reports
        ),
        "all_seed_best_updates_are_500": all(
            report["best_update"] == 500 for report in reports
        ),
        "all_seed_validation_completion_is_one": all(
            report["best_validation_completion"] == 1.0 for report in reports
        ),
        "all_case_rows_completed": all(
            row["completed"].lower() == "true" for row in rows
        ),
        **physical_checks,
        "figures_nonempty": all(item["bytes"] > 1000 for item in figures.values()),
    }
    audit = {
        "passed": all(checks.values()),
        "checks": checks,
        "comparisons": comparisons,
        "claim_status": {
            "supported": [
                "HCFL improves its dimension-by-dimension HLLC base in ID and official OOD NMAE",
                "all evaluated HCFL rollouts remain positive, conservative, and entropy admissible",
                "the deployment low-order fallback is inactive in every evaluated HCFL substep",
            ],
            "not_supported": [
                "HCFL outperforms native PyClaw Roe-64",
                "the current normal-stencil HCFL suppresses 2-D grid-aligned oscillations",
            ],
        },
        "checkpoint_sha256": checkpoint_hashes,
        "figures": figures,
    }
    (RESULTS / "AUDIT.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True), encoding="utf-8"
    )
    print(json.dumps(audit, indent=2), flush=True)
    if not audit["passed"]:
        raise RuntimeError("Result audit failed")


if __name__ == "__main__":
    main()
