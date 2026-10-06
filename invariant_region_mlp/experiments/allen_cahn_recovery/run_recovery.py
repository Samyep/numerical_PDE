"""Staged runner for the Allen--Cahn truncated-MLP recovery study."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
import gc
import hashlib
import json
import math
import os
from pathlib import Path
import time
from typing import Any, Iterable

import numpy as np

from .allen_cahn_equation import (
    AllenCahnEquation,
    BECK_NUMERICAL_RADIUS,
    CERTIFIED_INTERVAL,
    KeyedRandomTree,
    PUBLISHED_DS_REFERENCES,
    PUBLISHED_MLP_REFERENCES,
    deterministic_test_point,
    theorem_example_radius,
)
from .beck_truncated_mlp import BeckTruncatedMLP, RawFullHistoryMLP
from .ir_interval_mlp import IntervalIRMLP


HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[1]
RESULT_ROOT = PROJECT_ROOT / "results" / "allen_cahn_recovery"
RAW_ROOT = RESULT_ROOT / "raw"


@dataclass(frozen=True)
class Task:
    stage: str
    dimension: int
    depth: int
    sample_size: int
    repetition: int
    seed: int
    test_point: str = "zero"
    radius: float = BECK_NUMERICAL_RADIUS
    radius_label: str = "beck_r4"

    @property
    def stem(self) -> str:
        radius_tag = self.radius_label.replace(".", "p")
        return (
            f"d{self.dimension:04d}_n{self.depth:02d}_M{self.sample_size:03d}_"
            f"{self.test_point}_{radius_tag}_rep{self.repetition:03d}"
        )


def stable_seed(stage: str, dimension: int, depth: int, sample_size: int, rep: int) -> int:
    payload = f"allen-cahn-recovery-v1|{stage}|{dimension}|{depth}|{sample_size}|{rep}"
    return int.from_bytes(hashlib.sha256(payload.encode("ascii")).digest()[:8], "little")


def source_tasks(repetitions: int = 5) -> list[Task]:
    tasks: list[Task] = []
    for dimension in (10, 100, 1000):
        for depth in range(1, 6):
            for repetition in range(repetitions):
                tasks.append(
                    Task(
                        "source",
                        dimension,
                        depth,
                        depth,
                        repetition,
                        stable_seed("source", dimension, depth, depth, repetition),
                    )
                )
    return tasks


def pilot_tasks(repetitions: int = 3) -> list[Task]:
    """Adaptive Stage B pilot within the unchanged companion PDE."""

    configs = {
        10: [
            (2, 2),
            (2, 4),
            (2, 8),
            (2, 16),
            (3, 2),
            (3, 3),
            (3, 4),
            (3, 6),
            (3, 8),
            (4, 2),
            (4, 3),
            (4, 4),
            (5, 2),
            (5, 3),
            (6, 2),
            (7, 2),
            (8, 2),
        ],
        100: [(3, 2), (4, 2), (5, 2), (6, 2), (4, 3), (5, 3)],
        1000: [(4, 2), (5, 2), (6, 2)],
    }
    tasks: list[Task] = []
    for dimension, cells in configs.items():
        for depth, sample_size in cells:
            for repetition in range(repetitions):
                tasks.append(
                    Task(
                        "pilot",
                        dimension,
                        depth,
                        sample_size,
                        repetition,
                        stable_seed(
                            "pilot", dimension, depth, sample_size, repetition
                        ),
                    )
                )
    return tasks


def final_tasks(repetitions: int = 30) -> list[Task]:
    """Representative inactive/deep cells, frozen after inspecting the pilot."""

    cells = {
        10: [(2, 16), (5, 3), (8, 2)],
        100: [(3, 3), (5, 3), (6, 2)],
        1000: [(3, 3), (5, 2), (6, 2)],
    }
    tasks: list[Task] = []
    for dimension, configs in cells.items():
        for depth, sample_size in configs:
            for repetition in range(repetitions):
                tasks.append(
                    Task(
                        "final",
                        dimension,
                        depth,
                        sample_size,
                        repetition,
                        stable_seed(
                            "final", dimension, depth, sample_size, repetition
                        ),
                    )
                )
    return tasks


def equivalence_tasks() -> list[Task]:
    tasks: list[Task] = []
    configs = ((0, 2), (1, 2), (2, 2), (3, 2), (3, 3))
    for dimension in (1, 10, 100, 1000):
        for depth, sample_size in configs:
            for repetition in range(3):
                for point in ("zero", "ramp"):
                    for radius_label, radius in (
                        ("beck_r4", BECK_NUMERICAL_RADIUS),
                        ("theorem_schedule", theorem_example_radius(sample_size)),
                    ):
                        stage = "equivalence"
                        seed = stable_seed(
                            f"{stage}-{point}-{radius_label}",
                            dimension,
                            depth,
                            sample_size,
                            repetition,
                        )
                        tasks.append(
                            Task(
                                stage,
                                dimension,
                                depth,
                                sample_size,
                                repetition,
                                seed,
                                point,
                                radius,
                                radius_label,
                            )
                        )
    return tasks


def _method_row(
    task: Task,
    method: str,
    result: Any,
    projection_interval: tuple[float, float] | None,
) -> dict[str, Any]:
    reference = PUBLISHED_MLP_REFERENCES.get(task.dimension)
    row: dict[str, Any] = {
        **asdict(task),
        "method": method,
        "value": result.value,
        "wall_clock_seconds": result.wall_clock_seconds,
        "draw_fingerprint": result.draw_fingerprint,
        "finite_output": math.isfinite(result.value),
        "projection_low": (
            projection_interval[0] if projection_interval is not None else None
        ),
        "projection_high": (
            projection_interval[1] if projection_interval is not None else None
        ),
        **result.work,
    }
    if reference is not None and task.test_point == "zero":
        row.update(
            {
                "published_mlp_reference": reference,
                "published_ds_reference": PUBLISHED_DS_REFERENCES[task.dimension],
                "absolute_error": abs(result.value - reference),
                "relative_error": abs(result.value - reference) / abs(reference),
            }
        )
    else:
        row.update(
            {
                "published_mlp_reference": None,
                "published_ds_reference": None,
                "absolute_error": None,
                "relative_error": None,
            }
        )
    return row


def execute_task(task: Task) -> dict[str, Any]:
    equation = AllenCahnEquation(task.dimension)
    point = deterministic_test_point(task.dimension, task.test_point)
    tree = KeyedRandomTree(task.seed, task.dimension)
    capture = task.stage == "equivalence"

    beck = BeckTruncatedMLP(
        equation,
        task.sample_size,
        task.radius,
        tree,
        capture_trace=capture,
    ).run(task.depth, point)
    interval_ir = IntervalIRMLP(
        equation,
        task.sample_size,
        (-task.radius, task.radius),
        tree,
        capture_trace=capture,
        method_name="interval_ir",
    ).run(task.depth, point)

    methods = [
        _method_row(
            task, "beck_truncated", beck, (-task.radius, task.radius)
        ),
        _method_row(
            task, "interval_ir", interval_ir, (-task.radius, task.radius)
        ),
    ]
    if task.stage != "equivalence":
        raw = RawFullHistoryMLP(
            equation, task.sample_size, task.radius, tree
        ).run(task.depth, point)
        certified = IntervalIRMLP(
            equation,
            task.sample_size,
            CERTIFIED_INTERVAL,
            tree,
            method_name="certified_ir_0_1",
        ).run(task.depth, point)
        methods.extend(
            [
                _method_row(task, "raw", raw, None),
                _method_row(task, "certified_ir_0_1", certified, CERTIFIED_INTERVAL),
            ]
        )

    if beck.correction_trace.shape != interval_ir.correction_trace.shape:
        correction_max = math.inf
        correction_mean = math.inf
        correction_exact = 0
        correction_count = max(
            beck.correction_trace.size, interval_ir.correction_trace.size
        )
    elif beck.correction_trace.size:
        correction_diff = np.abs(
            beck.correction_trace - interval_ir.correction_trace
        )
        correction_max = float(np.max(correction_diff))
        correction_mean = float(np.mean(correction_diff))
        correction_exact = int(np.count_nonzero(correction_diff == 0.0))
        correction_count = int(correction_diff.size)
    else:
        correction_max = 0.0
        correction_mean = 0.0
        correction_exact = 0
        correction_count = 0

    equivalence = {
        **asdict(task),
        "beck_value": beck.value,
        "interval_ir_value": interval_ir.value,
        "absolute_value_discrepancy": abs(beck.value - interval_ir.value),
        "exact_value_match": beck.value == interval_ir.value,
        "correction_count": correction_count,
        "exact_correction_matches": correction_exact,
        "max_correction_discrepancy": correction_max,
        "mean_correction_discrepancy": correction_mean,
        "draw_fingerprints_match": (
            beck.draw_fingerprint == interval_ir.draw_fingerprint
        ),
        "work_counters_match": beck.work == interval_ir.work,
    }
    return {"task": asdict(task), "methods": methods, "equivalence": equivalence}


def _atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + f".{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=True),
        encoding="utf-8",
    )
    temporary.replace(path)


def execute_and_save(task: Task) -> tuple[str, bool, float]:
    output = RAW_ROOT / task.stage / f"{task.stem}.json"
    if output.exists():
        return str(output), True, 0.0
    started = time.perf_counter()
    payload = execute_task(task)
    elapsed = time.perf_counter() - started
    payload["task_wall_clock_seconds"] = elapsed
    _atomic_json(output, payload)
    gc.collect()
    return str(output), False, elapsed


def run_tasks(tasks: Iterable[Task], stage: str, workers: int) -> None:
    task_list = list(tasks)
    started = time.perf_counter()
    completed: list[dict[str, Any]] = []
    if workers == 1:
        iterator = (execute_and_save(task) for task in task_list)
        for index, (path, skipped, elapsed) in enumerate(iterator, start=1):
            completed.append(
                {"path": str(Path(path).relative_to(RESULT_ROOT)), "skipped": skipped, "seconds": elapsed}
            )
            print(f"[{index}/{len(task_list)}] {Path(path).name} {'cached' if skipped else 'complete'}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(execute_and_save, task): task for task in task_list}
            for index, future in enumerate(as_completed(futures), start=1):
                path, skipped, elapsed = future.result()
                completed.append(
                    {"path": str(Path(path).relative_to(RESULT_ROOT)), "skipped": skipped, "seconds": elapsed}
                )
                print(f"[{index}/{len(task_list)}] {Path(path).name} {'cached' if skipped else 'complete'}", flush=True)
    manifest = {
        "schema_version": 1,
        "stage": stage,
        "status": "complete",
        "workers": workers,
        "tasks_total": len(task_list),
        "tasks_completed": len(completed),
        "elapsed_seconds": time.perf_counter() - started,
        "task_seconds_sum": sum(item["seconds"] for item in completed),
        "tasks": sorted(completed, key=lambda item: item["path"]),
    }
    _atomic_json(RESULT_ROOT / f"{stage}_manifest.json", manifest)
    print(json.dumps({key: value for key, value in manifest.items() if key != "tasks"}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "stage", choices=("source", "pilot", "final", "equivalence")
    )
    parser.add_argument("--workers", type=int, default=max(1, min(8, os.cpu_count() or 1)))
    parser.add_argument("--repetitions", type=int)
    args = parser.parse_args()
    if args.workers < 1:
        raise SystemExit("--workers must be positive")
    if args.stage == "source":
        tasks = source_tasks(args.repetitions or 5)
    elif args.stage == "pilot":
        tasks = pilot_tasks(args.repetitions or 3)
    elif args.stage == "final":
        tasks = final_tasks(args.repetitions or 30)
    else:
        if args.repetitions is not None:
            raise SystemExit("equivalence has a fixed multi-seed grid")
        tasks = equivalence_tasks()
    run_tasks(tasks, args.stage, args.workers)


if __name__ == "__main__":
    main()
