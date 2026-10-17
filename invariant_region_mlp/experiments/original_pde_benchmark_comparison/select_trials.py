"""Reference-free successive-halving decisions for Actor trials."""

from __future__ import annotations

import argparse
import json

from .benchmarks import BENCHMARK_NAMES
from .common import RAW_DIR, RESULTS_DIR, atomic_json
from .run_actor import TRIALS


BROAD_TRIALS = ("pilot192", "wide384", "deep_elu", "quality128")


def select(benchmark: str, seed: int = 0) -> dict[str, object]:
    rows: list[dict[str, object]] = []
    for trial in BROAD_TRIALS:
        path = RAW_DIR / f"actor_{benchmark}_{trial}_s{seed}.json"
        record = json.loads(path.read_text())
        if len(record["phases"]) < 1:
            raise RuntimeError(f"incomplete broad trial: {trial}")
        phase = record["phases"][0]
        policy = phase["reference_free_policy_validation"]
        rows.append(
            {
                "trial": trial,
                "policy_value": float(policy["new_mean"]),
                "source_minus_new": float(policy["source_minus_new_mean"]),
                "policy_improvement_state_se": float(
                    policy["source_minus_new_state_se"]
                ),
                "heldout_label_rmse": float(
                    phase["training"]["selected_validation_rmse"]
                ),
                "label_replication_rmse": float(
                    phase["labels"]["validation_label_replication_rmse"]
                ),
                "training_seconds": float(phase["cumulative_training_seconds"]),
                "ground_truth_used": False,
            }
        )
    ranked = sorted(
        rows,
        key=lambda row: (float(row["policy_value"]), float(row["heldout_label_rmse"])),
    )
    decision: dict[str, object] = {
        "benchmark": benchmark,
        "seed": seed,
        "stage": "B broad one-update successive halving",
        "criterion": (
            "lowest common-bank reference-free tangent-policy value; held-out "
            "tangent-label RMSE breaks exact ties"
        ),
        "ranking": ranked,
        "selected_for_continuation": [row["trial"] for row in ranked[:2]],
        "eliminated": [row["trial"] for row in ranked[2:]],
        "ground_truth_used": False,
        "final_test_opened": False,
    }
    atomic_json(RESULTS_DIR / f"stage_b_selection_{benchmark}.json", decision)
    return decision


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark", choices=BENCHMARK_NAMES, required=True)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    print(json.dumps(select(args.benchmark, args.seed), indent=2))


if __name__ == "__main__":
    main()

