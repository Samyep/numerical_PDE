"""Controlled reference-grid precision ablation for direct-vector 1D Euler HCFL.

The only scientific factor changed is the fine-grid resolution of the
Rusanov + SSP-RK2 trajectory teacher (512 versus 2048 cells). Both fine
solutions are conservatively restricted to the same 64-cell learning grid.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import torch


HERE = Path(__file__).resolve().parent
HCFL_ROOT = HERE.parents[1]
COVERAGE_EXPERIMENT = HCFL_ROOT / "experiments" / "euler_1d_wave_coverage_ablation"
if str(COVERAGE_EXPERIMENT) not in sys.path:
    sys.path.insert(0, str(COVERAGE_EXPERIMENT))

import run_coverage_ablation as shared  # noqa: E402


base = shared.base
ARMS = ("reference_512", "reference_2048")
LOW_REFERENCE_CELLS = 512
HIGH_REFERENCE_CELLS = 2048
CANONICAL_CASES = shared.CANONICAL_CASES
RHO_FLOOR = shared.RHO_FLOOR
PRESSURE_FLOOR = shared.PRESSURE_FLOOR
FD_ENTROPY_TOLERANCE = shared.FD_ENTROPY_TOLERANCE


def conservative_restrict(values: np.ndarray, target_cells: int) -> np.ndarray:
    """Cell-average a periodic fine-grid state onto a nested target grid."""
    fine_cells = values.shape[-2]
    if fine_cells % target_cells:
        raise ValueError(f"{fine_cells} is not divisible by {target_cells}")
    factor = fine_cells // target_cells
    return values.reshape(*values.shape[:-2], target_cells, factor, 3).mean(axis=-2)


def strict_rollout_reference(
    initial: np.ndarray,
    fine_cells: int,
    nsnap: int = base.NSNAP,
) -> np.ndarray:
    """Roll out a strict fine-grid teacher and average it onto 64 cells."""
    if initial.shape[-2:] != (fine_cells, 3):
        raise ValueError(
            f"Expected initial shape (..., {fine_cells}, 3), got {initial.shape}"
        )
    if fine_cells % base.NCOARSE:
        raise ValueError(f"{fine_cells} must be divisible by {base.NCOARSE}")

    state = np.asarray(initial, dtype=np.float64).copy()
    cell_width = 1.0 / fine_cells

    def coarse_snapshot(values: np.ndarray) -> np.ndarray:
        return conservative_restrict(values, base.NCOARSE).astype(np.float32)

    snapshots = [coarse_snapshot(state)]
    for _ in range(nsnap - 1):
        remaining = base.DT_SNAPSHOT
        while remaining > 1.0e-14:
            rho = state[..., 0]
            pressure = base.np_pressure(state)
            velocity = state[..., 1] / rho
            sound_speed = np.sqrt(base.GAMMA * pressure / rho)
            max_speed = float(np.max(np.abs(velocity) + sound_speed))
            dt = min(remaining, 0.25 * cell_width / max(max_speed, 1.0e-12))
            state = base.np_ssprk2(state, dt, cell_width)
            if (
                not np.isfinite(state).all()
                or (state[..., 0] <= 0).any()
                or (base.np_pressure(state) <= 0).any()
            ):
                raise RuntimeError(
                    f"Strict {fine_cells}-cell reference lost admissibility; "
                    "no state repair was applied."
                )
            remaining -= dt
        snapshots.append(coarse_snapshot(state))
    return np.stack(snapshots, axis=1)


def paired_rollouts(
    high_initial: np.ndarray,
    nsnap: int = base.NSNAP,
) -> tuple[torch.Tensor, torch.Tensor]:
    """Generate nested 512/2048 references from one common initial state."""
    low_initial = conservative_restrict(high_initial, LOW_REFERENCE_CELLS)
    low = strict_rollout_reference(low_initial, LOW_REFERENCE_CELLS, nsnap=nsnap)
    high = strict_rollout_reference(high_initial, HIGH_REFERENCE_CELLS, nsnap=nsnap)
    if not np.array_equal(low[:, 0], high[:, 0]):
        max_difference = float(np.max(np.abs(low[:, 0] - high[:, 0])))
        if max_difference > 2.0e-6:
            raise RuntimeError(
                f"Paired coarse initial states differ by {max_difference:.3e}"
            )
    return torch.from_numpy(low), torch.from_numpy(high)


def make_paired_training_data(seed: int) -> tuple[torch.Tensor, torch.Tensor]:
    """Create paired 580-trajectory training tensors at both precisions."""
    specifications = (
        (
            "ordinary",
            base.generate_ic(220, HIGH_REFERENCE_CELLS, 6000 + seed, ood=False),
        ),
        (
            "broad_random",
            base.generate_ic(260, HIGH_REFERENCE_CELLS, 7000 + seed, ood=True),
        ),
        (
            "random_extreme",
            base.generate_extreme_ic(100, HIGH_REFERENCE_CELLS, 8000 + seed),
        ),
    )
    low_groups: list[torch.Tensor] = []
    high_groups: list[torch.Tensor] = []
    for name, initial in specifications:
        print(json.dumps({"stage": "training_reference", "group": name}))
        low, high = paired_rollouts(initial)
        low_groups.append(low)
        high_groups.append(high)
    return torch.cat(low_groups), torch.cat(high_groups)


def make_high_precision_evaluation_suite(seed: int) -> dict[str, torch.Tensor]:
    """Generate the shared evaluation suite exclusively at 2048 cells."""
    initial_conditions = {
        "ordinary_id": base.generate_ic(
            90, HIGH_REFERENCE_CELLS, 2000 + seed, ood=False
        ),
        "broad_random_in_support": base.generate_ic(
            90, HIGH_REFERENCE_CELLS, 3000 + seed, ood=True
        ),
        "moderate_ood_high_frequency": shared.generate_moderate_ood_ic(
            90, HIGH_REFERENCE_CELLS, 4000 + seed
        ),
    }
    suite: dict[str, torch.Tensor] = {}
    for name, initial in initial_conditions.items():
        print(json.dumps({"stage": "evaluation_reference", "split": name}))
        suite[name] = torch.from_numpy(
            strict_rollout_reference(initial, HIGH_REFERENCE_CELLS)
        )

    states = {
        "sod": ((1.0, 0.0, 1.0), (0.125, 0.0, 0.1)),
        "lax": ((0.445, 0.698, 3.528), (0.5, 0.0, 0.571)),
        "collision": ((1.0, 2.0, 1.0), (1.0, -2.0, 1.0)),
        "strong_pressure": ((1.0, 0.0, 5.0), (1.0, 0.0, 0.05)),
        "near_vacuum_expansion": ((1.0, -2.0, 0.4), (1.0, 2.0, 0.4)),
    }
    for name, (left, right) in states.items():
        print(json.dumps({"stage": "evaluation_reference", "split": name}))
        initial = shared._two_state_periodic(
            HIGH_REFERENCE_CELLS,
            HIGH_REFERENCE_CELLS // 2,
            left,
            right,
            0,
        )[None, ...]
        suite[name] = torch.from_numpy(
            strict_rollout_reference(
                initial,
                HIGH_REFERENCE_CELLS,
                nsnap=shared.CANONICAL_NSNAP,
            )
        )
    return suite


def reference_discrepancy_rows(
    low: torch.Tensor,
    high: torch.Tensor,
    state_std: torch.Tensor,
    seed: int,
) -> list[dict[str, Any]]:
    """Quantify the actual target change caused by reference refinement."""
    slices = {
        "ordinary": slice(0, 220),
        "broad_random": slice(220, 480),
        "random_extreme": slice(480, 580),
        "all": slice(0, 580),
    }
    rows: list[dict[str, Any]] = []
    for group, subset in slices.items():
        difference = (low[subset] - high[subset]) / state_std
        initial_difference = difference[:, 0]
        final_difference = difference[:, -1]
        rows.append(
            {
                "seed": seed,
                "group": group,
                "ntraj": int(low[subset].shape[0]),
                "initial_normalized_rmse": float(
                    torch.sqrt(torch.mean(initial_difference**2))
                ),
                "all_snapshot_normalized_rmse": float(
                    torch.sqrt(torch.mean(difference**2))
                ),
                "final_snapshot_normalized_rmse": float(
                    torch.sqrt(torch.mean(final_difference**2))
                ),
            }
        )
    return rows


def comparison_rows(metrics: list[dict[str, Any]]) -> list[dict[str, Any]]:
    indexed = {(row["arm"], row["split"]): row for row in metrics}
    split_order = (
        "ordinary_id",
        "broad_random_in_support",
        "moderate_ood_high_frequency",
        *CANONICAL_CASES,
    )
    rows: list[dict[str, Any]] = []
    for split in split_order:
        low = indexed[("reference_512", split)]
        high = indexed[("reference_2048", split)]
        old = float(low["rollout_nrmse"])
        new = float(high["rollout_nrmse"])
        rows.append(
            {
                "seed": low["seed"],
                "split": split,
                "reference_512_nrmse": old,
                "reference_2048_nrmse": new,
                "relative_nrmse_change_percent": 100.0 * (new / old - 1.0),
                "reference_512_local_limiter_rate": low[
                    "local_limiter_intervention_rate"
                ],
                "reference_2048_local_limiter_rate": high[
                    "local_limiter_intervention_rate"
                ],
                "reference_512_fd_entropy_rate": low[
                    "fd_entropy_intervention_rate"
                ],
                "reference_2048_fd_entropy_rate": high[
                    "fd_entropy_intervention_rate"
                ],
            }
        )
    return rows


def screening_decision(metrics: list[dict[str, Any]]) -> dict[str, Any]:
    indexed = {(row["arm"], row["split"]): row for row in metrics}

    def relative_change(split: str) -> float:
        old = float(indexed[("reference_512", split)]["rollout_nrmse"])
        new = float(indexed[("reference_2048", split)]["rollout_nrmse"])
        return new / old - 1.0

    ordinary_change = relative_change("ordinary_id")
    moderate_change = relative_change("moderate_ood_high_frequency")
    canonical_changes = {case: relative_change(case) for case in CANONICAL_CASES}
    old_mean = float(
        np.mean(
            [
                indexed[("reference_512", case)]["rollout_nrmse"]
                for case in CANONICAL_CASES
            ]
        )
    )
    new_mean = float(
        np.mean(
            [
                indexed[("reference_2048", case)]["rollout_nrmse"]
                for case in CANONICAL_CASES
            ]
        )
    )
    improved_count = sum(change < 0.0 for change in canonical_changes.values())
    safe = all(
        float(row["min_rho"]) >= RHO_FLOOR
        and float(row["min_pressure"]) >= PRESSURE_FLOOR
        and float(row["entropy_violation_rate"]) == 0.0
        and float(row["max_total_entropy_change"]) <= FD_ENTROPY_TOLERANCE
        for row in metrics
    )
    passed = (
        ordinary_change <= 0.05
        and moderate_change <= 0.05
        and improved_count >= 3
        and new_mean / old_mean - 1.0 <= -0.05
        and max(canonical_changes.values()) <= 0.10
        and safe
    )
    return {
        "passed": passed,
        "recommendation": "confirm_with_more_seeds" if passed else "stop_or_modify",
        "ordinary_id_relative_change": ordinary_change,
        "moderate_ood_relative_change": moderate_change,
        "canonical_relative_changes": canonical_changes,
        "canonical_improved_count": improved_count,
        "reference_512_canonical_mean_nrmse": old_mean,
        "reference_2048_canonical_mean_nrmse": new_mean,
        "canonical_mean_relative_change": new_mean / old_mean - 1.0,
        "all_rollouts_safe": safe,
    }


def self_test() -> None:
    high_initial = base.generate_ic(3, HIGH_REFERENCE_CELLS, 123, ood=False)
    low, high = paired_rollouts(high_initial, nsnap=3)
    assert low.shape == high.shape == (3, 3, base.NCOARSE, 3)
    assert torch.allclose(low[:, 0], high[:, 0], atol=2.0e-6, rtol=0.0)
    assert bool(shared.admissible(low).all())
    assert bool(shared.admissible(high).all())
    print("self-test passed")


def run(args: argparse.Namespace) -> None:
    output = Path(args.outdir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    total_started = time.perf_counter()

    print("Generating paired 512/2048 training references...")
    data_started = time.perf_counter()
    low_data, high_data = make_paired_training_data(args.seed)
    evaluation_suite = make_high_precision_evaluation_suite(args.seed)
    data_seconds = time.perf_counter() - data_started

    primitive = base.primitive(low_data)
    mean = primitive.mean(dim=(0, 1, 2)).numpy()
    std = primitive.std(dim=(0, 1, 2)).numpy()
    state_std = torch.tensor(
        [float(low_data[..., channel].std()) for channel in range(3)]
    )

    discrepancy = reference_discrepancy_rows(
        low_data, high_data, state_std, args.seed
    )
    shared.write_csv(
        output / f"reference_discrepancy_seed{args.seed}.csv", discrepancy
    )

    datasets = {"reference_512": low_data, "reference_2048": high_data}
    metrics: list[dict[str, Any]] = []
    curves: list[dict[str, Any]] = []
    runtime_rows: list[dict[str, Any]] = []

    for arm in ARMS:
        model, curve, train_seconds = shared.train_model(
            datasets[arm],
            mean,
            std,
            state_std,
            args.seed,
            args.iters,
            args.width,
            args.batch_size,
            args.lr,
            arm,
        )
        curves.extend(curve)
        parameter_count = sum(parameter.numel() for parameter in model.parameters())
        torch.save(model.state_dict(), output / f"{arm}_seed{args.seed}.pt")

        eval_started = time.perf_counter()
        for split, data in evaluation_suite.items():
            row = shared.evaluate_dataset(model, data, state_std, split)
            row.update(
                {
                    "seed": args.seed,
                    "arm": arm,
                    "parameter_count": parameter_count,
                    "train_seconds": train_seconds,
                }
            )
            metrics.append(row)
            print(
                json.dumps(
                    {
                        "seed": args.seed,
                        "arm": arm,
                        "split": split,
                        "nrmse": row["rollout_nrmse"],
                    }
                )
            )
        runtime_rows.append(
            {
                "seed": args.seed,
                "arm": arm,
                "parameter_count": parameter_count,
                "train_seconds": train_seconds,
                "evaluation_seconds": time.perf_counter() - eval_started,
                "final_logged_loss": curve[-1]["loss"],
            }
        )

    comparison = comparison_rows(metrics)
    decision = screening_decision(metrics)
    shared.write_csv(output / f"metrics_seed{args.seed}.csv", metrics)
    shared.write_csv(output / f"comparison_seed{args.seed}.csv", comparison)
    shared.write_csv(output / f"training_curve_seed{args.seed}.csv", curves)
    shared.write_csv(output / f"runtime_seed{args.seed}.csv", runtime_rows)

    metadata = {
        "seed": args.seed,
        "single_scientific_factor": "training_reference_grid_resolution",
        "training_composition_per_arm": {
            "ordinary": 220,
            "broad_random": 260,
            "random_extreme": 100,
        },
        "reference_512": "512-cell Rusanov + SSP-RK2, strict/no repair",
        "reference_2048": "2048-cell Rusanov + SSP-RK2, strict/no repair",
        "paired_initial_condition_method": (
            "generate at 2048 cells; conservatively average to 512 cells; "
            "identical 64-cell initial states"
        ),
        "evaluation_reference_cells": HIGH_REFERENCE_CELLS,
        "model": "direct-vector HLLC-HCFL",
        "iterations": args.iters,
        "width": args.width,
        "batch_size": args.batch_size,
        "learning_rate": args.lr,
        "ordinary_and_moderate_nsnap": base.NSNAP,
        "canonical_nsnap": shared.CANONICAL_NSNAP,
        "entropy_residual_tolerance": shared.ENTROPY_RESIDUAL_TOLERANCE,
        "fully_discrete_entropy_tolerance": FD_ENTROPY_TOLERANCE,
        "stored_dtype": "float32",
        "solver_compute_dtype": "float64",
        "pyclaw_available": False,
        "data_generation_seconds": data_seconds,
        "total_wall_seconds": time.perf_counter() - total_started,
        "screening_decision": decision,
    }
    (output / f"metadata_seed{args.seed}.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(json.dumps(decision, indent=2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--iters", type=int, default=1100)
    parser.add_argument("--width", type=int, default=72)
    parser.add_argument("--batch-size", type=int, default=56)
    parser.add_argument("--lr", type=float, default=3.0e-4)
    parser.add_argument("--outdir", default=str(HERE / "results"))
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    parsed = parse_args()
    if parsed.self_test:
        self_test()
    else:
        run(parsed)
