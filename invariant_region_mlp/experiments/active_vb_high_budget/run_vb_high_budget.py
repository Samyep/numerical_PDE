"""Resumable runner for the published-sigma high-budget VB experiment."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[1]
REPOSITORY_ROOT = PROJECT_ROOT.parent
RESULT_ROOT = PROJECT_ROOT / "results" / "active_vb_high_budget"
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from vb_equation import (  # noqa: E402
    PUBLISHED_SIGMA,
    SCASML_AUDITED_COMMIT,
    ViscousBurgersEquation,
    make_test_points,
    nonlinear_signal_summary,
)
from vb_mlp_methods import (  # noqa: E402
    BATCH_BOX,
    F_ZERO,
    RAW,
    SAMPLE_BALL,
    SAMPLE_BOX,
    Z_ONLY_BOX,
    Z_ZERO,
    MethodSpec,
    core_methods,
    parse_method,
    pilot_methods,
    run_single_repetition,
    save_repetition,
    tuning_methods,
)


TARGET_GRID = [
    (2, 2),
    (2, 4),
    (2, 8),
    (2, 16),
    (2, 24),
    (2, 32),
    (3, 2),
    (3, 3),
    (3, 4),
    (3, 6),
    (3, 8),
    (3, 10),
    (4, 2),
    (4, 3),
    (4, 4),
    (5, 2),
]

STAGE2_GRID = [
    (2, 2),
    (2, 8),
    (2, 24),
    (2, 32),
    (3, 3),
    (3, 6),
    (3, 10),
    (4, 2),
    (4, 3),
    (5, 2),
]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def git_output(*args: str) -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(REPOSITORY_ROOT), *args],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp-{os.getpid()}")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=True), encoding="utf-8")
    os.replace(temporary, path)


def run_signal_gate(*, n_samples: int = 100_000, resume: bool = True) -> Path:
    output = RESULT_ROOT / "signal_strength.json"
    if output.exists() and resume:
        print(f"signal gate already exists: {output}", flush=True)
        return output
    rows = []
    started = time.perf_counter()
    for d in (20, 40, 60, 80):
        row = nonlinear_signal_summary(d, n_samples=n_samples, sigma=PUBLISHED_SIGMA)
        rows.append(row)
        print(
            "signal",
            d,
            "mean|f|=",
            f"{row['abs_f']['mean']:.6g}",
            "f0 relL2=",
            f"{row['fzero_relative_l2']:.6g}",
            flush=True,
        )
    payload = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "sigma": PUBLISHED_SIGMA,
        "sample_geometry": "t~Uniform([0,T]), x~Uniform([-0.5,0.5]^d)",
        "rows": rows,
        "elapsed_seconds": time.perf_counter() - started,
        "provenance": {
            "experiment_base_commit": git_output("rev-parse", "HEAD"),
            "branch": git_output("branch", "--show-current"),
            "SCaSML_audited_commit": SCASML_AUDITED_COMMIT,
        },
    }
    atomic_json(output, payload)
    return output


def parse_csv_ints(value: str) -> list[int]:
    return [int(item) for item in value.split(",") if item]


def parse_configs(value: str) -> list[tuple[int, int]]:
    result = []
    for item in value.split(","):
        if not item:
            continue
        n_text, m_text = item.split(":", maxsplit=1)
        result.append((int(n_text), int(m_text)))
    return result


def _point_selection(points: dict[str, np.ndarray], subset: str) -> dict[str, np.ndarray]:
    if subset == "all":
        mask = np.ones(len(points["t"]), dtype=bool)
    elif subset == "validation":
        mask = points["is_validation"]
    elif subset == "test":
        mask = ~points["is_validation"]
    else:
        raise ValueError(f"unknown point subset {subset!r}")
    return {key: value[mask] for key, value in points.items()}


def _task_filename(task: dict[str, Any]) -> str:
    method = task["method"]["name"].replace(".", "p")
    return (
        f"d{task['d']:03d}_n{task['n']}_M{task['M']:03d}_"
        f"{method}_rep{task['repetition']:02d}.npz"
    )


def _execute_task(task: dict[str, Any]) -> dict[str, Any]:
    output = Path(task["output"])
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    points = make_test_points(
        task["d"],
        n_interior=task["n_interior"],
        n_boundary=task["n_boundary"],
        seed=task["point_seed"],
    )
    selected = _point_selection(points, task["point_subset"])
    equation = ViscousBurgersEquation(task["d"], sigma=task["sigma"])
    method = MethodSpec(**task["method"])
    result = run_single_repetition(
        equation=equation,
        method=method,
        n=task["n"],
        M=task["M"],
        repetition=task["repetition"],
        t=selected["t"],
        x=selected["x"],
        is_validation=selected["is_validation"],
        base_seed=task["base_seed"],
        chunk_size=task["chunk_size"],
        time_beta_alpha=task["time_beta_alpha"],
        trace_draws=task["trace_draws"],
    )
    result["metadata"]["point_subset"] = task["point_subset"]
    result["metadata"]["n_interior_source"] = task["n_interior"]
    result["metadata"]["n_boundary_source"] = task["n_boundary"]
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.stem + f".tmp-{os.getpid()}.npz")
    save_repetition(str(temporary), result)
    os.replace(temporary, output)
    metadata = result["metadata"]
    return {
        "path": str(output.relative_to(RESULT_ROOT)).replace("\\", "/"),
        "dimension": task["d"],
        "n": task["n"],
        "M": task["M"],
        "method": method.name,
        "repetition": task["repetition"],
        "wall_clock_seconds": metadata["wall_clock_seconds"],
        "test_value_relative_l2": metadata["metrics"]["test"]["value_relative_l2"],
        "all_value_relative_l2": metadata["metrics"]["all"]["value_relative_l2"],
        "f_evals": metadata["work"]["f_evals"],
        "total_stochastic_samples": metadata["work"]["total_stochastic_samples"],
    }


def _read_task_summary(path: Path) -> dict[str, Any]:
    with np.load(path, allow_pickle=False) as data:
        metadata = json.loads(str(data["metadata_json"]))
    return {
        "path": str(path.relative_to(RESULT_ROOT)).replace("\\", "/"),
        "dimension": metadata["dimension"],
        "n": metadata["n"],
        "M": metadata["M"],
        "method": metadata["method"]["name"],
        "repetition": metadata["repetition"],
        "wall_clock_seconds": metadata["wall_clock_seconds"],
        "test_value_relative_l2": metadata["metrics"]["test"]["value_relative_l2"],
        "all_value_relative_l2": metadata["metrics"]["all"]["value_relative_l2"],
        "f_evals": metadata["work"]["f_evals"],
        "total_stochastic_samples": metadata["work"]["total_stochastic_samples"],
    }


def _default_plan(args: argparse.Namespace) -> list[dict[str, Any]]:
    if args.stage == "pilot":
        dimensions = [20]
        configs = TARGET_GRID
        methods = pilot_methods()
        repetitions = 2
    elif args.stage == "tuning":
        dimensions = [20, 40, 60, 80]
        configs = STAGE2_GRID
        methods = [RAW, SAMPLE_BOX, Z_ZERO] + tuning_methods()
        repetitions = 3
    else:
        dimensions = [20, 40, 60, 80]
        configs = STAGE2_GRID
        methods = core_methods()
        repetitions = 5
    if args.dimensions:
        dimensions = parse_csv_ints(args.dimensions)
    if args.configs:
        configs = parse_configs(args.configs)
    if args.methods:
        methods = [parse_method(name) for name in args.methods.split(",") if name]
    if args.repetitions is not None:
        repetitions = args.repetitions
    return [
        {
            "d": d,
            "n": n,
            "M": M,
            "repetitions": repetitions,
            "methods": [method.to_dict() for method in methods],
        }
        for d in dimensions
        for n, M in configs
    ]


def _load_plan(args: argparse.Namespace) -> list[dict[str, Any]]:
    if args.plan is None:
        return _default_plan(args)
    payload = json.loads(Path(args.plan).read_text(encoding="utf-8"))
    plan = payload["plan"] if isinstance(payload, dict) else payload
    normalized = []
    for entry in plan:
        methods = []
        for method in entry["methods"]:
            methods.append(parse_method(method).to_dict() if isinstance(method, str) else method)
        normalized.append({**entry, "methods": methods})
    return normalized


def _save_points(stage_name: str, args: argparse.Namespace, dimensions: list[int]) -> None:
    point_dir = RESULT_ROOT / "test_points"
    point_dir.mkdir(parents=True, exist_ok=True)
    for d in dimensions:
        path = point_dir / f"{stage_name}_d{d:03d}.npz"
        if path.exists():
            continue
        points = make_test_points(
            d,
            n_interior=args.n_interior,
            n_boundary=args.n_boundary,
            seed=args.point_seed,
        )
        np.savez_compressed(path, **points)


def run_experiment(args: argparse.Namespace) -> Path:
    if args.sigma != PUBLISHED_SIGMA and not args.allow_nonpublished_sigma:
        raise ValueError("non-published sigma requires --allow-nonpublished-sigma")
    stage_name = args.output_stage or (f"{args.stage}_smoke" if args.quick else args.stage)
    plan = _load_plan(args)
    dimensions = sorted({int(entry["d"]) for entry in plan})
    _save_points(stage_name, args, dimensions)
    repetition_dir = RESULT_ROOT / "repetitions" / stage_name
    tasks: list[dict[str, Any]] = []
    for entry in plan:
        for repetition in range(int(entry["repetitions"])):
            for method in entry["methods"]:
                task = {
                    "d": int(entry["d"]),
                    "n": int(entry["n"]),
                    "M": int(entry["M"]),
                    "repetition": repetition,
                    "method": method,
                    "n_interior": args.n_interior,
                    "n_boundary": args.n_boundary,
                    "point_seed": args.point_seed,
                    "point_subset": entry.get("point_subset", args.point_subset),
                    "base_seed": args.base_seed,
                    "chunk_size": args.chunk_size,
                    "time_beta_alpha": args.time_beta_alpha,
                    "sigma": args.sigma,
                    "trace_draws": args.trace_draws,
                }
                task["output"] = str(repetition_dir / _task_filename(task))
                tasks.append(task)
    if args.limit is not None:
        tasks = tasks[: args.limit]

    manifest_name = "pilot.json" if stage_name == "pilot" else f"{stage_name}_manifest.json"
    manifest_path = RESULT_ROOT / manifest_name
    if manifest_path.exists() and not args.resume:
        raise FileExistsError(f"refusing to overwrite {manifest_path}; use --resume")
    completed: list[dict[str, Any]] = []
    pending: list[dict[str, Any]] = []
    for task in tasks:
        path = Path(task["output"])
        if path.exists():
            if not args.resume:
                raise FileExistsError(f"refusing to overwrite {path}")
            completed.append(_read_task_summary(path))
        else:
            pending.append(task)

    manifest: dict[str, Any] = {
        "schema_version": 1,
        "stage": stage_name,
        "status": "running",
        "created_or_resumed_utc": utc_now(),
        "branch": git_output("branch", "--show-current"),
        "experiment_base_commit": git_output("rev-parse", "HEAD"),
        "SCaSML_audited_commit": SCASML_AUDITED_COMMIT,
        "sigma": args.sigma,
        "time_beta_alpha": args.time_beta_alpha,
        "point_geometry": {
            "n_interior": args.n_interior,
            "n_boundary": args.n_boundary,
            "point_seed": args.point_seed,
            "default_subset": args.point_subset,
            "validation_fraction": 0.2,
        },
        "plan": plan,
        "tasks_total": len(tasks),
        "tasks_completed": len(completed),
        "results": sorted(completed, key=lambda row: row["path"]),
    }
    atomic_json(manifest_path, manifest)
    print(
        f"stage={stage_name} total={len(tasks)} existing={len(completed)} pending={len(pending)} "
        f"workers={args.workers}",
        flush=True,
    )
    started = time.perf_counter()
    if pending:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = {executor.submit(_execute_task, task): task for task in pending}
            progress_stride = max(1, len(pending) // 50)
            for index, future in enumerate(as_completed(futures), start=1):
                task = futures[future]
                try:
                    summary = future.result()
                except BaseException:
                    manifest["status"] = "failed"
                    manifest["failed_task"] = task
                    manifest["tasks_completed"] = len(completed)
                    manifest["results"] = sorted(completed, key=lambda row: row["path"])
                    atomic_json(manifest_path, manifest)
                    for item in futures:
                        item.cancel()
                    raise
                completed.append(summary)
                if index % progress_stride == 0 or index == len(pending):
                    print(
                        f"progress {index}/{len(pending)} latest="
                        f"d{summary['dimension']} n{summary['n']} M{summary['M']} "
                        f"{summary['method']} rep{summary['repetition']} "
                        f"relL2={summary['all_value_relative_l2']:.5g} "
                        f"sec={summary['wall_clock_seconds']:.3f}",
                        flush=True,
                    )
                    manifest["tasks_completed"] = len(completed)
                    manifest["results"] = sorted(completed, key=lambda row: row["path"])
                    manifest["elapsed_this_invocation_seconds"] = time.perf_counter() - started
                    atomic_json(manifest_path, manifest)
    manifest["status"] = "complete"
    manifest["completed_utc"] = utc_now()
    manifest["tasks_completed"] = len(completed)
    manifest["results"] = sorted(completed, key=lambda row: row["path"])
    manifest["elapsed_this_invocation_seconds"] = time.perf_counter() - started
    atomic_json(manifest_path, manifest)
    print(f"complete: {manifest_path}", flush=True)
    return manifest_path


def build_parser(*, forced_stage: str | None = None) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    if forced_stage is None:
        parser.add_argument("stage", choices=("pilot", "tuning", "main", "custom"))
    parser.add_argument("--plan", type=Path)
    parser.add_argument("--output-stage", help="result subdirectory/manifest stem (useful for refinements)")
    parser.add_argument("--dimensions", help="comma-separated dimensions")
    parser.add_argument("--configs", help="comma-separated n:M pairs")
    parser.add_argument("--methods", help="comma-separated method names")
    parser.add_argument("--repetitions", type=int)
    parser.add_argument("--n-interior", type=int)
    parser.add_argument("--n-boundary", type=int)
    parser.add_argument("--point-subset", choices=("all", "validation", "test"))
    parser.add_argument("--point-seed", type=int, default=20261006)
    parser.add_argument("--base-seed", type=int, default=20261006)
    parser.add_argument("--chunk-size", type=int, default=16)
    parser.add_argument("--time-beta-alpha", type=float, default=0.5)
    parser.add_argument("--sigma", type=float, default=PUBLISHED_SIGMA)
    parser.add_argument("--allow-nonpublished-sigma", action="store_true")
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--trace-draws", action="store_true")
    parser.add_argument("--skip-signal-gate", action="store_true")
    if forced_stage is not None:
        parser.set_defaults(stage=forced_stage)
    return parser


def cli_main(argv: list[str] | None = None, *, forced_stage: str | None = None) -> None:
    parser = build_parser(forced_stage=forced_stage)
    args = parser.parse_args(argv)
    if args.n_interior is None:
        args.n_interior = 96 if args.quick else (192 if args.stage == "pilot" else 1000)
    if args.n_boundary is None:
        args.n_boundary = 24 if args.quick else (48 if args.stage == "pilot" else 200)
    if args.point_subset is None:
        args.point_subset = "validation" if args.stage == "tuning" else "all"
    if args.quick and args.repetitions is None:
        args.repetitions = 1
    if not args.skip_signal_gate:
        run_signal_gate(n_samples=10_000 if args.quick else 100_000, resume=args.resume)
    run_experiment(args)


if __name__ == "__main__":
    cli_main()
