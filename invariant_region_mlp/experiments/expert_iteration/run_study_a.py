"""Resumable driver for the preregistered Study-A decision subset.

The default ``all`` command runs only d={20,100,160}, matching the user's
explicit decision-point instruction.  It does not silently treat omitted
dimensions, secondary MLP cells, or Study B as completed.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
from typing import Any, Iterable

import numpy as np
import torch

from lqg import (
    BASE_SEED,
    LQGEquation,
    antithetic_reference,
    make_evaluation_points,
    relative_l2,
)
from network import CHECKPOINTS, FrozenMLP, checkpoint_paths, train_surrogate
from solver import METHODS, baseline_row, run_repetition


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
RESULTS = PROJECT / "results" / "expert_iteration"
ARTIFACTS = RESULTS / "artifacts"
ROW_DIR = ARTIFACTS / "rows"
NETWORK_DIR = ARTIFACTS
REFERENCE_JSON = RESULTS / "A_reference.json"
DEFAULT_DIMS = (20, 100, 160)
PRIMARY_CELL = (2, 10)
SECONDARY_CELLS = ((2, 32), (3, 6))
FINAL_METHODS = (
    "path",
    "scasml",
    "scasml_noclip",
    "mlp",
    "mlp_clip",
    "path_clip",
    "oracle_state",
)


def _atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=True), encoding="utf-8")
    temporary.replace(path)


def _git_hash() -> str:
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=PROJECT, text=True
    ).strip()


def _reference_path(d: int) -> Path:
    return ARTIFACTS / "references" / f"lqg_d{d}.npz"


def run_references(
    dimensions: Iterable[int],
    *,
    tolerance: float = 1e-4,
    initial_pairs: int = 16_384,
    batch_pairs: int = 16_384,
    max_pairs: int = 16_777_216,
    point_chunk: int = 8,
) -> None:
    summaries: dict[str, Any] = {}
    if REFERENCE_JSON.exists():
        summaries = json.loads(REFERENCE_JSON.read_text(encoding="utf-8"))
    for d in dimensions:
        path = _reference_path(d)
        if path.exists() and str(d) in summaries:
            print(f"REFERENCE d={d} already complete", flush=True)
            continue
        eq = LQGEquation.create(d)
        t, x, is_validation = make_evaluation_points(eq)
        quad_u, quad_z = eq.hopf_cole(t, x)

        def progress(done: int, total: int, maximum_se: float, pairs: int) -> None:
            print(
                f"REFERENCE d={d} points={done}/{total} pairs={pairs} "
                f"chunk_max_se={maximum_se:.3e}",
                flush=True,
            )

        started = time.perf_counter()
        mc = antithetic_reference(
            eq,
            t,
            x,
            tolerance=tolerance,
            initial_pairs=initial_pairs,
            batch_pairs=batch_pairs,
            max_pairs=max_pairs,
            point_chunk=point_chunk,
            progress=progress,
        )
        elapsed = time.perf_counter() - started
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(path.stem + ".tmp.npz")
        np.savez_compressed(
            temporary,
            t=t,
            x=x,
            is_validation=is_validation,
            reference_u=mc["u"],
            reference_z=mc["z"],
            f_zero_u=mc["f_zero_u"],
            u_se=mc["u_se"],
            pairs=mc["pairs"],
            quadrature_u=quad_u,
            quadrature_z=quad_z,
            c1=eq.c1,
            c2=eq.c2,
        )
        temporary.replace(path)
        test = ~is_validation
        summaries[str(d)] = {
            "dimension": d,
            "gate_A_G1": bool(np.max(mc["u_se"][test]) <= tolerance),
            "tolerance": tolerance,
            "max_u_standard_error_test": float(np.max(mc["u_se"][test])),
            "max_u_standard_error_all": float(np.max(mc["u_se"])),
            "min_antithetic_pairs": int(np.min(mc["pairs"])),
            "max_antithetic_pairs": int(np.max(mc["pairs"])),
            "terminal_samples_total": int(np.sum(mc["terminal_samples"])),
            "mc_vs_quadrature_u_relative_l2_test": relative_l2(mc["u"][test], quad_u[test]),
            "mc_vs_quadrature_z_relative_l2_test": relative_l2(mc["z"][test], quad_z[test]),
            "f_zero_relative_l2_test": relative_l2(mc["f_zero_u"][test], mc["u"][test]),
            "wall_clock_seconds": elapsed,
            "reference_file": str(path.relative_to(PROJECT)).replace("\\", "/"),
            "quadrature_role": "independent audit and recursive-child truth; reported errors use MC",
            "seed_scheme": "SeedSequence([20261101,1,1,d,0,20,point_chunk])",
        }
        _atomic_json(REFERENCE_JSON, summaries)
        print(f"REFERENCE d={d} complete in {elapsed:.1f}s", flush=True)


def run_training(dimensions: Iterable[int], network_seeds: Iterable[int], *, device: str) -> None:
    for d in dimensions:
        eq = LQGEquation.create(d)
        for network_seed in network_seeds:
            directory = NETWORK_DIR / "networks" / f"d{d}" / f"seed{network_seed}"
            train_surrogate(eq, network_seed, directory, device=device)


def _load_reference(d: int) -> dict[str, np.ndarray]:
    path = _reference_path(d)
    if not path.exists():
        raise FileNotFoundError(f"reference is missing: {path}")
    with np.load(path, allow_pickle=False) as data:
        return {key: data[key] for key in data.files}


def _row_stem(
    d: int,
    network_seed: int,
    checkpoint: int,
    n: int,
    M: int,
    repetition: int,
    method: str,
    laplacian: str,
) -> str:
    return (
        f"d{d}_seed{network_seed}_ckpt{checkpoint}_n{n}_M{M}_rep{repetition}_"
        f"{method}_{laplacian}"
    )


def _save_row(stem: str, row: dict[str, Any], prediction: np.ndarray) -> None:
    ROW_DIR.mkdir(parents=True, exist_ok=True)
    row_path = ROW_DIR / f"{stem}.json"
    prediction_path = ROW_DIR / f"{stem}.npz"
    temporary = prediction_path.with_name(prediction_path.stem + ".tmp.npz")
    np.savez_compressed(temporary, prediction_u=np.asarray(prediction[:, 0], dtype=np.float64))
    temporary.replace(prediction_path)
    row = dict(row)
    row["prediction_file"] = str(prediction_path.relative_to(PROJECT)).replace("\\", "/")
    _atomic_json(row_path, row)


def _run_one(
    *,
    eq: LQGEquation,
    reference: dict[str, np.ndarray],
    surrogate: FrozenMLP | None,
    method: str,
    checkpoint: int,
    network_seed: int,
    n: int,
    M: int,
    repetition: int,
    laplacian_label: str,
    chunk_size: int,
) -> None:
    stem = _row_stem(
        eq.d, network_seed, checkpoint, n, M, repetition, method, laplacian_label
    )
    if (ROW_DIR / f"{stem}.json").exists():
        print(f"ROW {stem} already complete", flush=True)
        return
    truth = np.concatenate([reference["reference_u"][:, None], reference["reference_z"]], axis=1)
    if method in {"surrogate", "f_zero"}:
        row, prediction = baseline_row(
            equation=eq,
            surrogate=surrogate,
            method_name=method,
            checkpoint=checkpoint,
            network_seed=network_seed,
            n=n,
            M=M,
            repetition=repetition,
            t=reference["t"],
            x=reference["x"],
            is_validation=reference["is_validation"],
            reference_state=truth,
            f_zero_u=reference["f_zero_u"],
        )
    else:
        row, prediction = run_repetition(
            equation=eq,
            surrogate=surrogate,
            method_name=method,
            checkpoint=checkpoint,
            network_seed=network_seed,
            n=n,
            M=M,
            repetition=repetition,
            t=reference["t"],
            x=reference["x"],
            is_validation=reference["is_validation"],
            reference_state=truth,
            chunk_size=chunk_size,
            trace_draws=True,
        )
    row["code_commit_at_run"] = _git_hash()
    row["base_seed"] = BASE_SEED
    row["subset"] = "decision_d20_d100_d160"
    _save_row(stem, row, prediction)
    print(
        f"ROW {stem} relL2={row['value_relative_l2']:.6g} "
        f"seconds={row['wall_clock_seconds']:.1f}",
        flush=True,
    )


def run_rows(
    dimensions: Iterable[int],
    network_seeds: Iterable[int],
    *,
    device: str,
    phase: str,
    chunk_size: int,
) -> None:
    for d in dimensions:
        eq = LQGEquation.create(d)
        reference = _load_reference(d)
        for network_seed in network_seeds:
            tasks: list[tuple[int, int, int, int, str, bool]] = []
            # tuple: checkpoint,n,M,repetition,method,exact_laplacian
            if phase in {"critical", "decision", "all"}:
                for repetition in range(5):
                    for method in ("surrogate", "f_zero", *FINAL_METHODS):
                        tasks.append((2500, *PRIMARY_CELL, repetition, method, d <= 50))
            if phase in {"quality", "decision", "all"}:
                for checkpoint in (500, 1000):
                    for repetition in range(5):
                        for method in ("surrogate", "path", "scasml_noclip"):
                            tasks.append((checkpoint, *PRIMARY_CELL, repetition, method, d <= 50))
            if phase in {"secondary", "all"}:
                for n, M in SECONDARY_CELLS:
                    for repetition in range(5):
                        for method in ("surrogate", "f_zero", *FINAL_METHODS):
                            tasks.append((2500, n, M, repetition, method, d <= 50))
            if phase in {"exact-control", "decision", "all"} and d == 100 and network_seed == 0:
                for repetition in range(5):
                    for method in ("path", "scasml_noclip"):
                        tasks.append((2500, *PRIMARY_CELL, repetition, method, True))

            loaded: dict[tuple[int, bool], FrozenMLP] = {}
            for checkpoint, n, M, repetition, method, exact_lap in tasks:
                needs_surrogate = method not in {"mlp", "mlp_clip", "f_zero"}
                surrogate: FrozenMLP | None = None
                lap_label = "none"
                if needs_surrogate:
                    key = (checkpoint, exact_lap)
                    if key not in loaded:
                        checkpoint_file = (
                            NETWORK_DIR
                            / "networks"
                            / f"d{d}"
                            / f"seed{network_seed}"
                            / f"checkpoint_{checkpoint}.pt"
                        )
                        if not checkpoint_file.exists():
                            raise FileNotFoundError(f"network checkpoint is missing: {checkpoint_file}")
                        loaded[key] = FrozenMLP.load(
                            checkpoint_file,
                            d=d,
                            network_seed=network_seed,
                            device=device,
                            exact_laplacian=exact_lap,
                        )
                    surrogate = loaded[key]
                    lap_label = "exact" if exact_lap else f"hutch{len(surrogate.probes)}"
                _run_one(
                    eq=eq,
                    reference=reference,
                    surrogate=surrogate,
                    method=method,
                    checkpoint=checkpoint,
                    network_seed=network_seed,
                    n=n,
                    M=M,
                    repetition=repetition,
                    laplacian_label=lap_label,
                    chunk_size=chunk_size,
                )


def collect_rows() -> list[dict[str, Any]]:
    rows = []
    if ROW_DIR.exists():
        for path in sorted(ROW_DIR.glob("*.json")):
            rows.append(json.loads(path.read_text(encoding="utf-8")))
    return rows


def write_rows_csv() -> None:
    rows = collect_rows()
    path = RESULTS / "A_rows.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = sorted({key for row in rows for key in row})
    temporary = path.with_suffix(".tmp.csv")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def write_environment() -> None:
    payload = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "torch": torch.__version__,
        "cuda_available": torch.cuda.is_available(),
        "cuda_runtime": torch.version.cuda,
        "device": torch.cuda.get_device_name(0) if torch.cuda.is_available() else "cpu",
        "logical_cpus": os.cpu_count(),
        "code_commit": _git_hash(),
    }
    _atomic_json(RESULTS / "environment.json", payload)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command", choices=("reference", "train", "run", "collect", "all")
    )
    parser.add_argument("--dimensions", nargs="+", type=int, default=list(DEFAULT_DIMS))
    parser.add_argument("--network-seeds", nargs="+", type=int, default=[0, 1, 2])
    parser.add_argument("--device", default="cuda")
    parser.add_argument(
        "--phase",
        choices=("critical", "quality", "secondary", "exact-control", "decision", "all"),
        default="decision",
    )
    parser.add_argument("--chunk-size", type=int, default=8)
    parser.add_argument("--reference-tolerance", type=float, default=1e-4)
    parser.add_argument("--reference-initial-pairs", type=int, default=16_384)
    parser.add_argument("--reference-batch-pairs", type=int, default=16_384)
    parser.add_argument("--reference-max-pairs", type=int, default=16_777_216)
    parser.add_argument("--reference-point-chunk", type=int, default=8)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    RESULTS.mkdir(parents=True, exist_ok=True)
    write_environment()
    if args.command in {"reference", "all"}:
        run_references(
            args.dimensions,
            tolerance=args.reference_tolerance,
            initial_pairs=args.reference_initial_pairs,
            batch_pairs=args.reference_batch_pairs,
            max_pairs=args.reference_max_pairs,
            point_chunk=args.reference_point_chunk,
        )
    if args.command in {"train", "all"}:
        run_training(args.dimensions, args.network_seeds, device=args.device)
    if args.command in {"run", "all"}:
        run_rows(
            args.dimensions,
            args.network_seeds,
            device=args.device,
            phase=args.phase,
            chunk_size=args.chunk_size,
        )
    if args.command in {"collect", "all", "run"}:
        write_rows_csv()


if __name__ == "__main__":
    main()
