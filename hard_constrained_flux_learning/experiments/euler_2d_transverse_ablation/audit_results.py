"""Machine-check the reported transverse-ablation artifacts and claim bounds."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
BASE_RESULTS = HERE.parent / "euler_2d_clawpack_quadrants" / "results"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    protocol = load_json(HERE / "protocol.json")
    summary = load_json(RESULTS / "summary.json")
    reports = {
        variant: load_json(RESULTS / f"{variant}_report_seed0.json")
        for variant in ("normal6wide", "flat18", "gated18")
    }
    flat_reports = [
        load_json(RESULTS / f"flat18_report_seed{seed}.json")
        for seed in (0, 1, 2)
    ]
    old_reports = [
        load_json(BASE_RESULTS / f"hcfl_s6_report_seed{seed}.json")
        for seed in (0, 1, 2)
    ]

    require(protocol["selected_architecture"] == "flat18 after seed-0 validation ablation", "selection record changed")
    require(reports["normal6wide"]["parameters"] == 10804, "normal capacity count")
    require(reports["flat18"]["parameters"] == 10804, "flat capacity count")
    require(reports["gated18"]["parameters"] == 10872, "gated parameter count")
    require(
        reports["flat18"]["best_validation_nmae"]
        < reports["normal6wide"]["best_validation_nmae"],
        "flat18 did not beat the equal-capacity control on selection data",
    )
    require(
        reports["flat18"]["best_validation_nmae"]
        < reports["gated18"]["best_validation_nmae"],
        "flat18 did not beat gated18 on selection data",
    )
    require(
        all(report["stop_reason"] == "validation_plateau_at_minimum_learning_rate" for report in flat_reports),
        "a flat18 seed did not reach the locked stopping condition",
    )

    flat_validation = np.asarray(
        [report["best_validation_nmae"] for report in flat_reports], dtype=np.float64
    )
    old_validation = np.asarray(
        [report["best_validation_nmae"] for report in old_reports], dtype=np.float64
    )
    require(float(flat_validation.mean()) < float(old_validation.mean()), "mean validation did not improve")

    between_seed = summary["flat18_between_seed_accuracy"]
    id_metrics = between_seed["test_id"]
    official_metrics = between_seed["official_quadrants"]
    require(id_metrics["seeds"] == [0, 1, 2], "ID seed set")
    require(official_metrics["seeds"] == [0, 1, 2], "official seed set")
    require(summary["accuracy"]["test_id/flat18"]["completion_rate"] == 1.0, "ID completion")
    require(summary["accuracy"]["official_quadrants/flat18"]["completion_rate"] == 1.0, "official completion")

    with (RESULTS / "case_metrics.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    flat_rows = [row for row in rows if row["method"] == "flat18"]
    require(len(flat_rows) == 39, "expected 36 ID and three official flat18 rows")
    require(all(row["completed"] == "True" for row in flat_rows), "failed rows were hidden")

    flat_safety = {
        key: value
        for key, value in summary["safety"].items()
        if key.startswith("flat18 seed")
    }
    require(len(flat_safety) == 6, "missing flat18 safety records")
    require(
        all(value["positivity_fallback_rate"] == 0.0 for value in flat_safety.values()),
        "positivity fallback was active",
    )
    require(
        all(value["entropy_fallback_rate"] == 0.0 for value in flat_safety.values()),
        "entropy fallback was active",
    )
    max_residual = max(
        value["maximum_interface_residual"] for value in flat_safety.values()
    )
    max_closure = max(
        value["maximum_conservation_closure"] for value in flat_safety.values()
    )
    require(max_residual <= 1.0e-6, "hard interface residual tolerance")
    require(max_closure <= 1.0e-7, "conservation closure tolerance")
    require(
        all(value["maximum_entropy_balance"] <= 5.0e-7 for value in flat_safety.values()),
        "fully-discrete entropy balance tolerance",
    )

    old_summary = load_json(BASE_RESULTS / "summary.json")
    old_id = old_summary["accuracy"]["test_id/HCFL-64"]
    old_official = old_summary["accuracy"]["official_quadrants/HCFL-64"]
    require(
        id_metrics["nmae_mean_between_seeds"] < old_id["nmae_mean_completed"],
        "reported ID mean did not improve",
    )
    require(
        id_metrics["density_tv_excess_mean_between_seeds"]
        < old_id["density_tv_excess_mean_completed"],
        "reported ID TV excess did not improve",
    )
    require(
        official_metrics["nmae_mean_between_seeds"]
        >= old_official["nmae_mean_completed"],
        "claim boundary changed: official NMAE unexpectedly marked as a win",
    )

    required_figures = (
        "validation_curves.png",
        "official_final_density.png",
        "official_density_linecuts.png",
    )
    require(all((RESULTS / name).is_file() for name in required_figures), "missing figure")

    audit = {
        "passed": True,
        "flat18_validation_nmae_mean": float(flat_validation.mean()),
        "flat18_validation_nmae_std": float(flat_validation.std()),
        "old_validation_nmae_mean": float(old_validation.mean()),
        "id_nmae_mean_between_seeds": id_metrics["nmae_mean_between_seeds"],
        "id_tv_excess_mean_between_seeds": id_metrics[
            "density_tv_excess_mean_between_seeds"
        ],
        "official_nmae_mean_between_seeds": official_metrics[
            "nmae_mean_between_seeds"
        ],
        "official_tv_excess_mean_between_seeds": official_metrics[
            "density_tv_excess_mean_between_seeds"
        ],
        "maximum_hard_interface_residual": max_residual,
        "maximum_relative_conservation_closure": max_closure,
        "flat18_completed_rows": len(flat_rows),
        "fallback_activations": 0,
        "claim_boundary": "ID improvement; no official long-OOD NMAE win; still behind PyClaw Roe-64",
    }
    (RESULTS / "AUDIT.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True), encoding="utf-8"
    )
    print(json.dumps(audit, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
