"""Train independent 64-cell HCFL replicates for Euler and SWE.

Only the paper-default method is repeated here: a symmetric four-cell
primitive stencil, central physical flux, nonnegative learned Roe
dissipation, proposal-feasibility loss, hard Tadmor projection, conservative
positivity limiting, and a fully-discrete entropy limiter.

Seed 0 is an immutable input from the earlier controlled stencil studies.
This script trains additional seeds with exactly the same convergence rule.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch


HERE = Path(__file__).resolve().parent
EXPERIMENTS = HERE.parent
EULER_DIR = EXPERIMENTS / "euler_1d_stencil_ablation"
SWE_DIR = EXPERIMENTS / "swe_1d_consistent_hcfl"
for module_dir in (EULER_DIR, SWE_DIR):
    if str(module_dir) not in sys.path:
        sys.path.insert(0, str(module_dir))

import run_stencil_ablation as euler  # noqa: E402
import run_swe_consistent as swe  # noqa: E402


EULER_ARM = euler.ARM_BY_NAME["central_nonnegative_feas_s4"]
SWE_ARM = next(
    arm for arm in swe.ARMS if arm.name == "central_nonnegative_feas_s4"
)
SEED0_EULER_RESULTS = EULER_DIR / "results"
SEED0_SWE_RESULTS = SWE_DIR / "results"


def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"No rows to write: {path}")
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def euler_cache(output: Path, seed: int) -> tuple[torch.Tensor, torch.Tensor]:
    path = output / f"data_cache_euler_seed{seed}.pt"
    if path.exists():
        payload = torch.load(path, map_location="cpu", weights_only=True)
        train = payload["train"]
        validation = payload["validation"]
        if tuple(train.shape) == (580, 16, 64, 3) and tuple(
            validation.shape
        ) == (136, 16, 64, 3):
            return train, validation
    train = euler.shared.make_baseline_training_data(seed)
    validation = euler.convergence.make_validation_data(seed)
    torch.save({"train": train, "validation": validation}, path)
    return train, validation


def common_args(seed: int, parsed: argparse.Namespace) -> SimpleNamespace:
    return SimpleNamespace(
        seed=seed,
        width=parsed.width,
        batch_size=parsed.batch_size,
        lr=parsed.lr,
        fixed_budget_updates=parsed.fixed_budget_updates,
        max_updates=parsed.max_updates,
        validation_interval=parsed.validation_interval,
        plateau_patience=parsed.plateau_patience,
        minimum_relative_improvement=parsed.minimum_relative_improvement,
        learning_rate_factor=parsed.learning_rate_factor,
        minimum_learning_rate=parsed.minimum_learning_rate,
    )


def train_euler(
    seed: int,
    args: argparse.Namespace,
    output: Path,
) -> dict[str, Any]:
    checkpoint = output / f"euler_hcfl64_s4_converged_best_seed{seed}.pt"
    report_path = output / f"euler_hcfl64_s4_report_seed{seed}.json"
    if args.resume and checkpoint.exists() and report_path.exists():
        return json.loads(report_path.read_text(encoding="utf-8"))

    train_data, validation_data = euler_cache(output, seed)
    primitive = euler.base.primitive(train_data)
    mean = primitive.mean(dim=(0, 1, 2)).numpy()
    std = primitive.std(dim=(0, 1, 2)).numpy()
    state_std = train_data.std(dim=(0, 1, 2))
    temp_arm = f"paper_euler_hcfl64_s4_seed{seed}"
    model, result, curve = euler.convergence.train_to_convergence(
        arm=temp_arm,
        model_name=EULER_ARM.model_name,
        train_data=train_data,
        validation_data=validation_data,
        mean=mean,
        std=std,
        state_std=state_std,
        seed=seed,
        width=args.width,
        batch_size=args.batch_size,
        learning_rate=args.lr,
        max_updates=args.max_updates,
        validation_interval=args.validation_interval,
        fixed_budget_updates=args.fixed_budget_updates,
        plateau_patience=args.plateau_patience,
        minimum_relative_improvement=args.minimum_relative_improvement,
        learning_rate_factor=args.learning_rate_factor,
        minimum_learning_rate=args.minimum_learning_rate,
        output=output,
        proposal_feasibility_weight=EULER_ARM.feasibility_weight,
    )
    temporary_checkpoint = output / f"{temp_arm}_converged_best_seed{seed}.pt"
    if not bool(result["converged"]):
        temporary_checkpoint.unlink(missing_ok=True)
        checkpoint.unlink(missing_ok=True)
        raise RuntimeError(
            f"Euler seed {seed} reached the cap without validation convergence"
        )
    temporary_checkpoint.replace(checkpoint)
    heldout = euler.convergence.evaluate_state(
        EULER_ARM.model_name,
        result["best"],
        mean,
        std,
        args.width,
        state_std,
        euler.shared.make_evaluation_suite(seed),
        seed,
        "paper_euler_hcfl64_s4",
        "converged_best",
    )
    summary = {
        key: value for key, value in result.items() if key not in ("fixed", "best")
    }
    report = {
        "system": "Euler 1D",
        "training_cells": 64,
        "training_data_shape": list(train_data.shape),
        "validation_data_shape": list(validation_data.shape),
        "method": "central + nonnegative Roe + feasibility + hard safety",
        "stencil_cells": 4,
        "stencil_shifts": list(model.flux_net.stencil_shifts),
        "checkpoint": checkpoint.name,
        "convergence": summary,
        "heldout_periodic_metrics": heldout,
    }
    write_csv(output / f"euler_hcfl64_s4_curve_seed{seed}.csv", curve)
    report_path.write_text(
        json.dumps(json_ready(report), indent=2), encoding="utf-8"
    )
    return report


def train_swe(
    seed: int,
    args: argparse.Namespace,
    output: Path,
) -> dict[str, Any]:
    checkpoint = output / f"swe_hcfl64_s4_converged_best_seed{seed}.pt"
    report_path = output / f"swe_hcfl64_s4_report_seed{seed}.json"
    if args.resume and checkpoint.exists() and report_path.exists():
        return json.loads(report_path.read_text(encoding="utf-8"))

    train_data, validation_data = swe.load_or_make_data(output, seed)
    mean, std, state_std = swe.prepare_statistics(train_data)
    local_args = common_args(seed, args)
    model, summary, curve = swe.train_to_convergence(
        SWE_ARM,
        train_data,
        validation_data,
        mean,
        std,
        state_std,
        local_args,
        output,
    )
    temporary_checkpoint = swe.checkpoint_path(output, SWE_ARM, seed)
    if not bool(summary["converged"]):
        temporary_checkpoint.unlink(missing_ok=True)
        checkpoint.unlink(missing_ok=True)
        raise RuntimeError(
            f"SWE seed {seed} reached the cap without validation convergence"
        )
    temporary_checkpoint.replace(checkpoint)
    validation_metrics = swe.evaluate_dataset(
        model, validation_data, state_std, "validation"
    )
    report = {
        "system": "SWE 1D",
        "training_cells": 64,
        "training_data_shape": list(train_data.shape),
        "validation_data_shape": list(validation_data.shape),
        "method": "central + nonnegative Roe + feasibility + hard safety",
        "stencil_cells": 4,
        "stencil_shifts": list(swe.stencil_shifts(4)),
        "checkpoint": checkpoint.name,
        "convergence": summary,
        "validation_metrics": validation_metrics,
    }
    write_csv(output / f"swe_hcfl64_s4_curve_seed{seed}.csv", curve)
    report_path.write_text(
        json.dumps(json_ready(report), indent=2), encoding="utf-8"
    )
    return report


def seed0_reports() -> tuple[dict[str, Any], dict[str, Any]]:
    euler_report = json.loads(
        (
            SEED0_EULER_RESULTS
            / "report_central_nonnegative_feas_s4_seed0.json"
        ).read_text(encoding="utf-8")
    )
    swe_report = json.loads(
        (
            SEED0_SWE_RESULTS
            / "report_central_nonnegative_feas_s4_seed0.json"
        ).read_text(encoding="utf-8")
    )
    return euler_report, swe_report


def aggregate_reports(
    euler_reports: list[dict[str, Any]],
    swe_reports: list[dict[str, Any]],
    output: Path,
) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for system, reports in (("Euler", euler_reports), ("SWE", swe_reports)):
        for seed, report in enumerate(reports):
            convergence = report["convergence"]
            rows.append(
                {
                    "system": system,
                    "seed": int(convergence.get("seed", seed)),
                    "best_update": convergence["best_update"],
                    "stop_update": convergence["stop_update"],
                    "validation_rollout_nrmse": convergence[
                        "best_validation_rollout_nrmse"
                    ],
                    "converged": convergence["converged"],
                    "stop_reason": convergence["stop_reason"],
                    "parameter_count": convergence["parameter_count"],
                }
            )
    write_csv(output / "hcfl64_replicate_convergence.csv", rows)

    aggregate: dict[str, Any] = {}
    for system in ("Euler", "SWE"):
        values = np.array(
            [
                row["validation_rollout_nrmse"]
                for row in rows
                if row["system"] == system
            ],
            dtype=np.float64,
        )
        aggregate[system] = {
            "seeds": [
                row["seed"] for row in rows if row["system"] == system
            ],
            "all_converged": all(
                row["converged"] for row in rows if row["system"] == system
            ),
            "validation_rollout_nrmse_mean": float(values.mean()),
            "validation_rollout_nrmse_sample_std": float(values.std(ddof=1)),
            "validation_rollout_nrmse_values": values.tolist(),
        }
    payload = {
        "scope": (
            "three independently initialized/data-seeded 64-cell training "
            "runs of the common four-cell HCFL method"
        ),
        "rows": rows,
        "aggregate": aggregate,
    }
    (output / "hcfl64_replicate_summary.json").write_text(
        json.dumps(json_ready(payload), indent=2), encoding="utf-8"
    )

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.7), constrained_layout=True)
    for axis, system in zip(axes, ("Euler", "SWE")):
        system_rows = [row for row in rows if row["system"] == system]
        seeds = [row["seed"] for row in system_rows]
        metrics = [row["validation_rollout_nrmse"] for row in system_rows]
        axis.bar(
            [str(seed) for seed in seeds],
            metrics,
            color="#0072B2" if system == "Euler" else "#D55E00",
            alpha=0.82,
        )
        axis.axhline(np.mean(metrics), color="black", linestyle="--", linewidth=1.2)
        axis.set_title(system)
        axis.set_xlabel("seed")
        axis.set_ylabel("validation rollout NRMSE")
        axis.grid(axis="y", alpha=0.2)
    fig.suptitle("64-cell HCFL validation-controlled replicates")
    fig.savefig(output / "hcfl64_replicate_convergence.png", dpi=220)
    plt.close(fig)
    return payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", nargs="+", type=int, default=[1, 2])
    parser.add_argument("--width", type=int, default=72)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--lr", type=float, default=3.0e-4)
    parser.add_argument("--fixed-budget-updates", type=int, default=1100)
    parser.add_argument("--max-updates", type=int, default=50000)
    parser.add_argument("--validation-interval", type=int, default=100)
    parser.add_argument("--plateau-patience", type=int, default=10)
    parser.add_argument("--minimum-relative-improvement", type=float, default=1.0e-3)
    parser.add_argument("--learning-rate-factor", type=float, default=0.3)
    parser.add_argument("--minimum-learning-rate", type=float, default=3.0e-6)
    parser.add_argument("--results-dir", type=Path, default=HERE / "results")
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output = args.results_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    seed0_euler, seed0_swe = seed0_reports()
    euler_reports = [seed0_euler]
    swe_reports = [seed0_swe]
    for seed in args.seeds:
        if seed == 0:
            continue
        print(json.dumps({"stage": "train", "system": "Euler", "seed": seed}), flush=True)
        euler_reports.append(train_euler(seed, args, output))
        print(json.dumps({"stage": "train", "system": "SWE", "seed": seed}), flush=True)
        swe_reports.append(train_swe(seed, args, output))
    summary = aggregate_reports(euler_reports, swe_reports, output)
    print(
        json.dumps(
            {
                "stage": "complete",
                "wall_seconds": time.perf_counter() - started,
                "aggregate": summary["aggregate"],
            },
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
