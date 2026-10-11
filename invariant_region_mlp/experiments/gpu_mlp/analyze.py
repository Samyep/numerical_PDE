"""Mechanical G-3 agreement gate, summaries, figures, and Study-G report."""

from __future__ import annotations

import json
import math
from pathlib import Path
import subprocess
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .protocol import (
    G3_CELLS,
    G3_METHODS,
    G3_REPETITIONS,
    RESULTS_ROOT,
    g3_tasks,
    g5_tasks,
    task_is_complete,
    task_path,
    write_json_atomic,
)


PACKAGE_ROOT = Path(__file__).resolve().parents[2]
REPORT_PATH = PACKAGE_ROOT / "docs" / "GPU_MLP_REPORT.md"
FIGURES_ROOT = RESULTS_ROOT / "figures"
G3_EXPECTED_CELLS = 32
CRITERIA = {
    "G-0": (
        "**G-0 (environment):** record GPU model, driver, CUDA, framework versions, FP64 throughput class. If no GPU is "
        "available, implement and validate on the framework's CPU backend and mark G-4 and G-5 \"not evaluated (no GPU)\"."
    ),
    "G-1": (
        "**G-1 (implementation):** PyTorch (or JAX), **float64 only**, breadth-first batched like the numpy code, supporting "
        "`raw`, `path`, `double`, `double_path`, `path_base`, `double_path_base` for P1 and MR; points processed in chunks to fit "
        "memory. Re-implement the equations; verify terminal, exact state and generator against the numpy equations to relative "
        "1e-12 on random inputs. Seeding: a device generator seeded from "
        "`SeedSequence([20261201, d, n, M, rep, chunk_index]).generate_state(2)`, identical for all methods of a cell, so GPU "
        "methods are paired with each other (they are not paired with CPU rows; say so)."
    ),
    "G-2": (
        "**G-2 (exactness, verdict):** with a draw-tape interface, record every primary and auxiliary draw of the numpy solver "
        "for one chunk and replay it in the GPU solver. For P1 and MR at d=100, cells (2,8), (3,4), (3,6), all six methods: "
        "every output (u and z) agrees with numpy to relative 1e-10. PASS only if all agree."
    ),
    "G-3": (
        "**G-3 (statistical agreement, verdict):** with native GPU randomness, P1 and MR, d in {100, 400}, cells (2,32) and "
        "(3,10), methods `path`, `double_path`, `path_base`, `double_path_base`: 30 GPU repetitions versus 30 new CPU repetitions "
        "(stored in `results/gpu_mlp/`). PASS if |difference of mean skill| <= 3 x Welch standard error in every cell."
    ),
    "G-4": (
        "**G-4 (speed, reported):** wall time per repetition (1,200 points) for `double_path_base`, P1, d=400, cells (3,16), "
        "(3,32), (4,6), (4,10), (4,16), on GPU and on CPU (numpy, same machine, report core count); peak GPU memory. CPU times "
        "that would exceed 2 hours may be extrapolated from one chunk; label them as extrapolated."
    ),
    "G-5": (
        "**G-5 (exploratory, only if G-2 and G-3 pass):** on GPU, n=4 with M in {6, 10, 16}, methods `path`, `double_path`, "
        "`double_path_base`; P1 and MR; d in {100, 400}; 10 repetitions. Report skill and all three costs next to the best n=2 "
        "and n=3 results from Study F. No verdict."
    ),
}


def _git(*args: str) -> str:
    root = Path(__file__).resolve().parents[3]
    return subprocess.run(
        ["git", *args], cwd=root, check=True, capture_output=True, text=True
    ).stdout.strip()


def _load_task_rows(tasks: list[dict[str, Any]]) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for task in tasks:
        if not task_is_complete(task):
            continue
        row = json.loads(task_path(task).read_text(encoding="utf-8"))
        records.append({key: value for key, value in row.items() if not isinstance(value, dict)})
    if not records:
        return pd.DataFrame()
    frame = pd.DataFrame(records)
    for column in (
        "dimension",
        "n",
        "M",
        "repetition",
        "skill",
        "root_gradient_skill",
        "generator_calls",
        "gaussian_vectors_drawn",
        "gaussian_coordinates_drawn",
        "wall_time_seconds",
        "peak_gpu_memory_bytes",
        "nonfinite_state_count",
    ):
        if column in frame:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return frame


def agreement() -> tuple[pd.DataFrame, dict[str, Any], pd.DataFrame]:
    rows = _load_task_rows(g3_tasks())
    records: list[dict[str, Any]] = []
    if not rows.empty:
        groups = rows.groupby(["pde_id", "dimension", "n", "M", "method"])
        for key, cell in groups:
            cpu = cell[cell.backend == "numpy_cpu"].skill.to_numpy(dtype=np.float64)
            gpu = cell[cell.backend == "torch_gpu"].skill.to_numpy(dtype=np.float64)
            complete = len(cpu) == len(G3_REPETITIONS) and len(gpu) == len(G3_REPETITIONS)
            cpu_mean = float(np.mean(cpu)) if len(cpu) else math.nan
            gpu_mean = float(np.mean(gpu)) if len(gpu) else math.nan
            cpu_var = float(np.var(cpu, ddof=1)) if len(cpu) > 1 else math.nan
            gpu_var = float(np.var(gpu, ddof=1)) if len(gpu) > 1 else math.nan
            se = (
                math.sqrt(cpu_var / len(cpu) + gpu_var / len(gpu))
                if len(cpu) > 1 and len(gpu) > 1
                else math.nan
            )
            difference = abs(gpu_mean - cpu_mean)
            passed = bool(
                complete
                and np.isfinite(se)
                and (difference <= 3.0 * se if se > 0.0 else difference == 0.0)
            )
            pde, d, n, M, method = key
            records.append(
                {
                    "pde_id": pde,
                    "dimension": int(d),
                    "n": int(n),
                    "M": int(M),
                    "method": method,
                    "cpu_repetitions": len(cpu),
                    "gpu_repetitions": len(gpu),
                    "cpu_mean_skill": cpu_mean,
                    "gpu_mean_skill": gpu_mean,
                    "absolute_mean_difference": difference,
                    "welch_standard_error": se,
                    "three_welch_standard_errors": 3.0 * se,
                    "standardized_difference": difference / se if se > 0.0 else math.inf,
                    "passed": passed,
                }
            )
    agreement_frame = pd.DataFrame(records)
    complete = len(agreement_frame) == G3_EXPECTED_CELLS and bool(
        (agreement_frame.cpu_repetitions == len(G3_REPETITIONS)).all()
        and (agreement_frame.gpu_repetitions == len(G3_REPETITIONS)).all()
    ) if not agreement_frame.empty else False
    verdict = (
        "PASS"
        if complete and bool(agreement_frame.passed.all())
        else ("FAIL" if complete else "NOT EVALUATED (COMPUTE)")
    )
    payload = {
        "verdict": verdict,
        "cells": len(agreement_frame),
        "expected_cells": G3_EXPECTED_CELLS,
        "passed_cells": int(agreement_frame.passed.sum()) if not agreement_frame.empty else 0,
        "failed_cells": int((~agreement_frame.passed).sum()) if complete else None,
        "cpu_rows": int((rows.backend == "numpy_cpu").sum()) if not rows.empty else 0,
        "gpu_rows": int((rows.backend == "torch_gpu").sum()) if not rows.empty else 0,
        "cpu_gpu_paired": False,
        "repetitions": len(G3_REPETITIONS),
    }
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    agreement_frame.to_csv(RESULTS_ROOT / "agreement.csv", index=False)
    write_json_atomic(RESULTS_ROOT / "agreement_summary.json", payload)
    return agreement_frame, payload, rows


def _n4_outputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = _load_task_rows(g5_tasks())
    rows.to_csv(RESULTS_ROOT / "n4_rows.csv", index=False)
    if rows.empty:
        summary = pd.DataFrame()
    else:
        summary = (
            rows.groupby(["pde_id", "dimension", "n", "M", "method"], as_index=False)
            .agg(
                repetitions=("repetition", "nunique"),
                mean_skill=("skill", "mean"),
                std_skill=("skill", "std"),
                mean_generator_calls=("generator_calls", "mean"),
                mean_gaussian_coordinates_drawn=("gaussian_coordinates_drawn", "mean"),
                mean_wall_time_seconds=("wall_time_seconds", "mean"),
                mean_peak_gpu_memory_bytes=("peak_gpu_memory_bytes", "mean"),
                nonfinite_state_count=("nonfinite_state_count", "sum"),
            )
            .sort_values(["pde_id", "dimension", "M", "method"])
        )
    summary.to_csv(RESULTS_ROOT / "n4_summary.csv", index=False)
    return rows, summary


def _study_f_best() -> pd.DataFrame:
    path = PACKAGE_ROOT / "results" / "baseline_n3" / "summary.csv"
    if not path.is_file():
        return pd.DataFrame()
    data = pd.read_csv(path)
    records: list[dict[str, Any]] = []
    for (pde, d, n), group in data.groupby(["pde_id", "dimension", "n"]):
        allowed = group[
            group.method.isin(["path", "path_base", "double_path", "double_path_base"])
        ]
        if allowed.empty:
            continue
        best = allowed.loc[allowed.mean_skill.idxmin()]
        records.append(
            {
                "pde_id": pde,
                "dimension": int(d),
                "n": int(n),
                "M": int(best.M),
                "method": best.method,
                "mean_skill": float(best.mean_skill),
                "mean_generator_calls": float(best.mean_generator_calls),
                "mean_gaussian_coordinates_drawn": float(best.mean_gaussian_coordinates_drawn),
                "mean_wall_time_seconds": float(best.mean_wall_time_seconds),
            }
        )
    return pd.DataFrame(records)


def _plot_speed(speed: pd.DataFrame) -> None:
    FIGURES_ROOT.mkdir(parents=True, exist_ok=True)
    fig, axis = plt.subplots(figsize=(8, 4.8))
    if speed.empty:
        axis.text(0.5, 0.5, "G-4 not evaluated", ha="center", va="center")
        axis.set_axis_off()
    else:
        labels = [f"({int(row.n)},{int(row.M)})" for _, row in speed.iterrows()]
        axis.plot(labels, speed.speedup_cpu_over_gpu, marker="o", color="#4c78a8")
        axis.axhline(1.0, color="black", linestyle="--", linewidth=1)
        axis.set_yscale("log")
        axis.set_xlabel("cell (n,M)")
        axis.set_ylabel("CPU/GPU wall-time speed-up")
        axis.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIGURES_ROOT / "speedup.png", dpi=180)
    plt.close(fig)


def _plot_n4(summary: pd.DataFrame) -> None:
    FIGURES_ROOT.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for row_index, pde in enumerate(("P1", "MR")):
        for column, d in enumerate((100, 400)):
            axis = axes[row_index, column]
            problem = summary[
                (summary.pde_id == pde) & (summary.dimension == d)
            ] if not summary.empty else summary
            if problem.empty:
                axis.text(0.5, 0.5, "G-5 not evaluated", ha="center", va="center")
            else:
                for method, color in (
                    ("path", "#4c78a8"),
                    ("double_path", "#e45756"),
                    ("double_path_base", "#f58518"),
                ):
                    data = problem[problem.method == method].sort_values("M")
                    axis.plot(data.M, data.mean_skill, marker="o", label=method, color=color)
                axis.set_yscale("log")
                axis.legend(frameon=False)
            axis.set_title(f"{pde}, d={d}")
            axis.set_xlabel("M (n=4)")
            axis.set_ylabel("mean test skill")
            axis.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIGURES_ROOT / "n4_skill.png", dpi=180)
    plt.close(fig)


def _number(value: Any) -> str:
    if value is None:
        return "NA"
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    return "NA" if not np.isfinite(number) else f"{number:.6g}"


def _report(
    validation: dict[str, Any],
    agreement_frame: pd.DataFrame,
    agreement_summary: dict[str, Any],
    speed: pd.DataFrame,
    n4_summary: pd.DataFrame,
    study_f: pd.DataFrame,
    environment: dict[str, Any],
    prereg_commit: str,
    analysis_commit: str,
) -> None:
    has_gpu = bool(environment.get("cuda_available"))
    g4 = "REPORTED" if len(speed) == 5 else ("NOT EVALUATED (NO GPU)" if not has_gpu else "NOT EVALUATED (COMPUTE)")
    gates_pass = validation.get("G-2", {}).get("verdict") == "PASS" and agreement_summary["verdict"] == "PASS"
    expected_n4 = 36
    g5 = (
        "REPORTED"
        if len(n4_summary) == expected_n4 and bool((n4_summary.repetitions == 10).all())
        else ("NOT RUN (G-2/G-3 GATE)" if not gates_pass else "NOT EVALUATED (COMPUTE)")
    )
    outcomes = {
        "G-0": "RECORDED",
        "G-1": validation.get("G-1", {}).get("verdict", "NOT EVALUATED"),
        "G-2": validation.get("G-2", {}).get("verdict", "NOT EVALUATED"),
        "G-3": agreement_summary["verdict"],
        "G-4": g4,
        "G-5": g5,
    }
    lines = [
        "# GPU MLP port: pre-registered report",
        "",
        "## Outcome table",
        "",
        "| Criterion | Outcome |",
        "|---|---|",
    ]
    for key in ("G-0", "G-1", "G-2", "G-3", "G-4", "G-5"):
        lines.append(f"| {key} | {outcomes[key]} |")
    lines.extend(["", "## Pre-registered criteria (verbatim)", ""])
    for key in ("G-0", "G-1", "G-2", "G-3", "G-4", "G-5"):
        lines.extend([f"### {key}", "", f"> {CRITERIA[key]}", "", f"**Outcome: {outcomes[key]}.**", ""])
        if key == "G-2" and validation.get("G-2"):
            item = validation["G-2"]
            lines.append(
                f"All {item['checks']}/{item['expected_checks']} checks were evaluated; maximum relative error {_number(item['max_relative_error'])} against threshold {_number(item['threshold'])}."
            )
            lines.append("")
        if key == "G-1" and validation.get("G-1"):
            item = validation["G-1"]
            lines.append(
                f"All {item['checks']} equation checks passed; maximum relative error {_number(item['max_relative_error'])} against threshold {_number(item['threshold'])}."
            )
            lines.append("")
        if key == "G-3" and not agreement_frame.empty:
            lines.extend(
                [
                    "| PDE | d | cell | method | CPU mean | GPU mean | |difference| / Welch SE | pass |",
                    "|---|---:|---|---|---:|---:|---:|---|",
                ]
            )
            for _, row in agreement_frame.iterrows():
                lines.append(
                    f"| {row.pde_id} | {int(row.dimension)} | ({int(row.n)},{int(row.M)}) | {row.method} | {_number(row.cpu_mean_skill)} | {_number(row.gpu_mean_skill)} | {_number(row.standardized_difference)} | {'yes' if row.passed else 'no'} |"
                )
            lines.append("")
        if key == "G-4" and not speed.empty:
            lines.extend(
                [
                    "| cell | GPU seconds | CPU seconds | CPU extrapolated | speed-up | peak GPU GiB |",
                    "|---|---:|---:|---|---:|---:|",
                ]
            )
            for _, row in speed.iterrows():
                lines.append(
                    f"| ({int(row.n)},{int(row.M)}) | {_number(row.gpu_wall_time_seconds)} | {_number(row.cpu_wall_time_seconds)} | {'yes' if row.cpu_extrapolated else 'no'} | {_number(row.speedup_cpu_over_gpu)} | {_number(row.peak_gpu_memory_bytes / 2**30)} |"
                )
            lines.append("")
        if key == "G-5" and not n4_summary.empty:
            lines.extend(
                [
                    "| PDE | d | M | method | skill | generator calls | Gaussian x d | seconds |",
                    "|---|---:|---:|---|---:|---:|---:|---:|",
                ]
            )
            for _, row in n4_summary.iterrows():
                lines.append(
                    f"| {row.pde_id} | {int(row.dimension)} | {int(row.M)} | {row.method} | {_number(row.mean_skill)} | {_number(row.mean_generator_calls)} | {_number(row.mean_gaussian_coordinates_drawn)} | {_number(row.mean_wall_time_seconds)} |"
                )
            lines.append("")

    if not study_f.empty:
        lines.extend(
            [
                "## Best Study-F n=2/n=3 references",
                "",
                "| PDE | d | n | cell/method | skill | generator calls | Gaussian x d | seconds |",
                "|---|---:|---:|---|---:|---:|---:|---:|",
            ]
        )
        for _, row in study_f.iterrows():
            lines.append(
                f"| {row.pde_id} | {int(row.dimension)} | {int(row.n)} | M={int(row.M)}, {row.method} | {_number(row.mean_skill)} | {_number(row.mean_generator_calls)} | {_number(row.mean_gaussian_coordinates_drawn)} | {_number(row.mean_wall_time_seconds)} |"
            )

    failures = [key for key in ("G-1", "G-2", "G-3") if outcomes[key] == "FAIL"]
    lines.extend(
        [
            "",
            "## Failures and deviations",
            "",
            f"- Mechanical failures: {', '.join(failures) if failures else 'none'}.",
            "- GPU and CPU G-3 rows use independent native random-number implementations and are not paired; agreement is assessed only by the preregistered Welch rule.",
            "- G-3 repetitions are numbered 100–129 so all CPU draws are new relative to Study F.",
            "- No post-hoc thresholds were introduced.",
            "",
            "## Provenance",
            "",
            f"- Frozen preregistration commit: `{prereg_commit}`.",
            f"- Analysis code commit: `{analysis_commit}`.",
            f"- Solver/run code commit: `{environment.get('code_commit', 'NA')}`.",
            f"- GPU: {environment.get('gpu_name', 'none')}; driver {environment.get('nvidia_driver', 'NA')}; CUDA {environment.get('cuda_runtime', 'NA')}; PyTorch {environment.get('torch', 'NA')}.",
            f"- FP64 throughput class: {environment.get('fp64_throughput_class', 'NA')}.",
            f"- CPU timing host: {int(speed.cpu_logical_cores.iloc[0]) if not speed.empty and 'cpu_logical_cores' in speed else 'NA'} logical cores; NumPy timings use the same machine as the GPU timings.",
            f"- G-3 rows: CPU {agreement_summary['cpu_rows']}/960, GPU {agreement_summary['gpu_rows']}/960.",
            f"- G-4 timing rows: {len(speed)}/5; G-5 rows: {int(n4_summary.repetitions.sum()) if not n4_summary.empty else 0}/360 across {len(n4_summary)}/36 summary cells.",
            "- Pairing: GPU methods within a cell share the preregistered device seed; GPU and CPU G-3 rows are intentionally unpaired. G-2 replayed every recorded primary and auxiliary NumPy draw.",
            "- Arithmetic: float64 only; point sets are the frozen 1,200 points and verdicts use the frozen 960-point test split.",
            "",
            "## Figures",
            "",
            "- `results/gpu_mlp/figures/speedup.png`",
            "- `results/gpu_mlp/figures/n4_skill.png`",
            "",
        ]
    )
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text("\n".join(lines), encoding="utf-8")


def analyze() -> dict[str, Any]:
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    validation_path = RESULTS_ROOT / "validation_summary.json"
    environment_path = RESULTS_ROOT / "environment.json"
    validation = json.loads(validation_path.read_text(encoding="utf-8")) if validation_path.is_file() else {}
    environment = json.loads(environment_path.read_text(encoding="utf-8")) if environment_path.is_file() else {}
    agreement_frame, agreement_summary, _rows = agreement()
    n4_rows, n4_summary = _n4_outputs()
    speed_path = RESULTS_ROOT / "speed.csv"
    speed = pd.read_csv(speed_path) if speed_path.is_file() else pd.DataFrame()
    study_f = _study_f_best()
    _plot_speed(speed)
    _plot_n4(n4_summary)
    prereg_commit = _git(
        "log",
        "--diff-filter=A",
        "-1",
        "--format=%H",
        "--",
        "invariant_region_mlp/experiments/baseline_n3/PREREGISTRATION.md",
    )
    analysis_commit = _git("rev-parse", "HEAD")
    payload = {
        "G-0": environment,
        "G-1": validation.get("G-1", {"verdict": "NOT EVALUATED"}),
        "G-2": validation.get("G-2", {"verdict": "NOT EVALUATED"}),
        "G-3": agreement_summary,
        "G-4": {"rows": len(speed), "expected_rows": 5},
        "G-5": {
            "rows": len(n4_rows),
            "expected_rows": len(g5_tasks()),
            "summary_cells": len(n4_summary),
        },
        "preregistration_commit": prereg_commit,
        "analysis_code_commit": analysis_commit,
    }
    write_json_atomic(RESULTS_ROOT / "analysis_summary.json", payload)
    _report(
        validation,
        agreement_frame,
        agreement_summary,
        speed,
        n4_summary,
        study_f,
        environment,
        prereg_commit,
        analysis_commit,
    )
    print(json.dumps(payload, indent=2, allow_nan=True), flush=True)
    return payload


if __name__ == "__main__":
    analyze()
