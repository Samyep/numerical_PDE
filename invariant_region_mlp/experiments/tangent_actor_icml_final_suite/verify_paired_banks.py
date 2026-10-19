#!/usr/bin/env python3
"""Audit whether frozen seed-0 whole-function methods share evaluation banks.

This is deliberately a *separate* publication-readiness guard. It does not
rewrite MAIN_RESULTS.csv or the original 16-check artifact audit.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path


def audit(main_results: Path, frozen_protocol: Path, dimension: int = 100) -> dict:
    protocol = json.loads(frozen_protocol.read_text(encoding="utf-8"))
    pdes = [item["key"] for item in protocol["pdes"]]
    methods = sorted(set(protocol["methods"]) - {"deep_bsde"})
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)

    with main_results.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if (
                int(row["dimension"]) == dimension
                and int(row["seed"]) == 0
                and row["method"] in methods
            ):
                groups[row["pde"]].append(row)

    results = {}
    paired_wins = []
    for pde in pdes:
        rows = groups[pde]
        present = [r["method"] for r in rows]
        hashes = {r["test_bank_sha256"].strip() for r in rows}
        counts = {int(float(r["states"])) for r in rows if r.get("states")}
        completeness = sorted(present) == methods
        same_bank = completeness and len(hashes) == 1 and "" not in hashes
        same_states = completeness and len(counts) == 1 and all(r.get("states") for r in rows)
        errors = {
            r["method"]: float(r["value_rmse"])
            for r in rows if r.get("value_rmse")
        }
        finite = len(errors) == len(rows) and all(
            math.isfinite(e) and e >= 0.0 for e in errors.values()
        )
        status = "PASS" if (completeness and same_bank and same_states and finite) else "FAIL"
        winner = min(errors, key=errors.get) if finite and errors else None
        if status == "PASS" and winner == "tangent_actor":
            paired_wins.append(pde)
        uncertainties = sorted({r["reference_uncertainty"] for r in rows})
        results[pde] = {
            "status": status,
            "methods_found": sorted(present),
            "bank_hash_by_method": {r["method"]: r["test_bank_sha256"] for r in rows},
            "states_by_method": {r["method"]: r.get("states", "") for r in rows},
            "same_bank": same_bank,
            "same_states": same_states,
            "winner_if_paired": winner if status == "PASS" else None,
            "reference_uncertainty_labels_vary": len(uncertainties) > 1,
        }

    all_matched = all(v["status"] == "PASS" for v in results.values())
    return {
        "status": "PASS" if all_matched else "FAIL",
        "dimension": dimension,
        "seed": 0,
        "pde_count": len(pdes),
        "paired_pde_count": sum(v["status"] == "PASS" for v in results.values()),
        "tangent_actor_strict_paired_wins": paired_wins,
        "results": results,
    }


def main() -> None:
    default_results = Path(__file__).resolve().parents[2] / "results" / "tangent_actor_icml_final_suite"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=default_results)
    parser.add_argument("--dimension", type=int, default=100)
    args = parser.parse_args()
    result = audit(
        args.results_dir / "MAIN_RESULTS.csv",
        args.results_dir / "FROZEN_PROTOCOL.json",
        args.dimension,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    if result["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
