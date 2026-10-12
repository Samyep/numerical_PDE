"""Mechanical Stage-A/interim and final analysis for the non-ridge study."""

from __future__ import annotations

import json
import math
from pathlib import Path
import subprocess
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from invariant_region_mlp.experiments.double_estimator.equations import PACKAGE_ROOT

from .protocol import (
    RESULTS_ROOT,
    all_tasks,
    stage_a_tasks,
    stage_b_extension_tasks,
    stage_b_n4_tasks,
    stage_b_s1_tasks,
    task_is_complete,
    task_path,
)


REPORT_PATH = PACKAGE_ROOT / "docs" / "NONRIDGE_REPORT.md"
FIGURES_ROOT = RESULTS_ROOT / "figures"
BASELINE_SUMMARY = PACKAGE_ROOT / "results" / "baseline_n3" / "summary.csv"
CURTAIL_PATH = RESULTS_ROOT / "curtailed_cells.json"

CRITERIA = {
    "N-1": (
        "**N-1 (pathwise matters off-ridge):** `path` skill <= 0.5 x `raw` skill at (2,64) and at (3,16), in all 8 "
        "(family, d, cell) combinations. A non-finite `raw` mean counts as satisfied."
    ),
    "N-2": (
        "**N-2 (baseline off-ridge):** `path_base` skill <= 0.5 x `path` skill at (3,16), in all 4 (family, d) combinations."
    ),
    "N-3": (
        "**N-3 (double complements baseline off-ridge):** `double_path_base` skill <= `double_path` skill at (3,16) and (3,24), "
        "in at least 7 of the 8 (family, d, M) combinations."
    ),
    "N-4": (
        "**N-4 (n = 3 below the n = 2 floor at d = 400):** for each family at d = 400, the minimum over the n = 3 cells and "
        "methods {`double_path`, `double_path_base`} of mean skill is <= 0.9 x the minimum over the n = 2 `path` cells "
        "(M = 16, 64, 256). PASS if it holds for both families. Report the same quantity at d = 100 (no verdict)."
    ),
    "N-5": (
        "**N-5 (dimension robustness):** `double_path_base` skill ratio d = 400 / d = 100 <= 1.3 at (3,16) and (3,24), for both "
        "families (4 comparisons; PASS only if all hold)."
    ),
}


def _git(*args: str) -> str:
    root = Path(__file__).resolve().parents[3]
    return subprocess.run(
        ["git", *args], cwd=root, check=True, capture_output=True, text=True
    ).stdout.strip()


def _load_rows() -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for task in all_tasks():
        if not task_is_complete(task):
            continue
        row = json.loads(task_path(task).read_text(encoding="utf-8"))
        records.append({key: value for key, value in row.items() if not isinstance(value, (dict, list))})
    if not records:
        return pd.DataFrame()
    frame = pd.DataFrame(records)
    frame.to_csv(RESULTS_ROOT / "rows.csv", index=False)
    return frame


def _summary(rows: pd.DataFrame) -> pd.DataFrame:
    if rows.empty:
        return pd.DataFrame()
    result = (
        rows.groupby(["stage", "pde_id", "dimension", "n", "M", "method"], as_index=False)
        .agg(
            repetitions=("repetition", "nunique"),
            mean_skill=("skill", "mean"),
            std_skill=("skill", "std"),
            mean_root_gradient_skill=("root_gradient_skill", "mean"),
            mean_generator_calls=("generator_calls", "mean"),
            mean_gaussian_vectors_drawn=("gaussian_vectors_drawn", "mean"),
            mean_gaussian_coordinates_drawn=("gaussian_coordinates_drawn", "mean"),
            mean_wall_time_seconds=("wall_time_seconds", "mean"),
            mean_root_level1_variance=("root_level1_gradient_integrand_variance", "mean"),
            nonfinite_state_count=("nonfinite_state_count", "sum"),
            nonfinite_generator_count=("nonfinite_generator_count", "sum"),
        )
        .sort_values(["pde_id", "dimension", "n", "M", "method"])
    )
    result.to_csv(RESULTS_ROOT / "summary.csv", index=False)
    return result


def _find(
    summary: pd.DataFrame,
    family: str,
    d: int,
    n: int,
    M: int,
    method: str,
) -> pd.Series | None:
    selected = summary[
        (summary.pde_id == family)
        & (summary.dimension == d)
        & (summary.n == n)
        & (summary.M == M)
        & (summary.method == method)
    ]
    if selected.empty:
        return None
    row = selected.iloc[0]
    expected = 10
    if CURTAIL_PATH.is_file():
        curtailments = json.loads(CURTAIL_PATH.read_text(encoding="utf-8"))
        key = "|".join(
            str(value)
            for value in (row.stage, family, d, n, M, method)
        )
        expected = int(curtailments.get(key, {}).get("max_repetitions", 10))
    return row if int(row.repetitions) == expected else None


def _evaluate_block(
    summary: pd.DataFrame,
    families: tuple[str, ...],
    *,
    n3_required: int,
) -> dict[str, Any]:
    n1_rows: list[dict[str, Any]] = []
    n2_rows: list[dict[str, Any]] = []
    n3_rows: list[dict[str, Any]] = []
    n4_rows: list[dict[str, Any]] = []
    n5_rows: list[dict[str, Any]] = []
    for family in families:
        for d in (100, 400):
            for n, M in ((2, 64), (3, 16)):
                path = _find(summary, family, d, n, M, "path")
                raw = _find(summary, family, d, n, M, "raw")
                complete = path is not None and raw is not None
                ratio = (
                    float(path.mean_skill / raw.mean_skill)
                    if complete and np.isfinite(raw.mean_skill) and raw.mean_skill != 0.0
                    else (0.0 if complete and not np.isfinite(raw.mean_skill) else math.nan)
                )
                n1_rows.append(
                    {
                        "family": family,
                        "dimension": d,
                        "n": n,
                        "M": M,
                        "path_over_raw": ratio,
                        "met": bool(complete and (not np.isfinite(raw.mean_skill) or path.mean_skill <= 0.5 * raw.mean_skill)),
                        "complete": complete,
                    }
                )
            path = _find(summary, family, d, 3, 16, "path")
            baseline = _find(summary, family, d, 3, 16, "path_base")
            complete = path is not None and baseline is not None
            ratio = float(baseline.mean_skill / path.mean_skill) if complete else math.nan
            n2_rows.append(
                {
                    "family": family,
                    "dimension": d,
                    "baseline_over_path": ratio,
                    "met": bool(complete and baseline.mean_skill <= 0.5 * path.mean_skill),
                    "complete": complete,
                }
            )
            for M in (16, 24):
                double = _find(summary, family, d, 3, M, "double_path")
                baseline = _find(summary, family, d, 3, M, "double_path_base")
                complete = double is not None and baseline is not None
                ratio = float(baseline.mean_skill / double.mean_skill) if complete else math.nan
                n3_rows.append(
                    {
                        "family": family,
                        "dimension": d,
                        "M": M,
                        "baseline_over_double": ratio,
                        "met": bool(complete and baseline.mean_skill <= double.mean_skill),
                        "complete": complete,
                    }
                )
            floor_cells = [
                _find(summary, family, d, 2, M, "path") for M in (16, 64, 256)
            ]
            deep_cells = [
                _find(summary, family, d, 3, M, method)
                for M in (16, 24)
                for method in ("double_path", "double_path_base")
            ]
            complete = all(cell is not None for cell in floor_cells + deep_cells)
            if complete:
                floor = min(floor_cells, key=lambda row: row.mean_skill)
                deep = min(deep_cells, key=lambda row: row.mean_skill)
                ratio = float(deep.mean_skill / floor.mean_skill)
                n4_rows.append(
                    {
                        "family": family,
                        "dimension": d,
                        "floor_M": int(floor.M),
                        "floor_skill": float(floor.mean_skill),
                        "best_n3_M": int(deep.M),
                        "best_n3_method": str(deep.method),
                        "best_n3_skill": float(deep.mean_skill),
                        "best_n3_over_floor": ratio,
                        "met": bool(ratio <= 0.9),
                        "complete": True,
                    }
                )
            else:
                n4_rows.append({"family": family, "dimension": d, "complete": False, "met": False})
        for M in (16, 24):
            low = _find(summary, family, 100, 3, M, "double_path_base")
            high = _find(summary, family, 400, 3, M, "double_path_base")
            complete = low is not None and high is not None
            ratio = float(high.mean_skill / low.mean_skill) if complete else math.nan
            n5_rows.append(
                {
                    "family": family,
                    "M": M,
                    "d400_over_d100": ratio,
                    "met": bool(complete and ratio <= 1.3),
                    "complete": complete,
                }
            )

    def verdict(rows: list[dict[str, Any]], required: int, *, scope: Any = None) -> dict[str, Any]:
        selected = rows if scope is None else [row for row in rows if scope(row)]
        complete = bool(selected) and all(row["complete"] for row in selected)
        wins = sum(bool(row["met"]) for row in selected)
        return {
            "verdict": "PASS" if complete and wins >= required else ("FAIL" if complete else "NOT EVALUATED (COMPUTE)"),
            "wins": wins,
            "total": len(selected),
            "required": required,
            "rows": selected,
        }

    return {
        "N-1": verdict(n1_rows, len(families) * 4),
        "N-2": verdict(n2_rows, len(families) * 2),
        "N-3": verdict(n3_rows, n3_required),
        "N-4": verdict(n4_rows, len(families), scope=lambda row: row["dimension"] == 400),
        "N-4-all-dimensions": n4_rows,
        "N-5": verdict(n5_rows, len(families) * 2),
    }


def _extended_floors(summary: pd.DataFrame) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for family in ("S0", "S05", "S1"):
        for d in (100, 400):
            initial_floor = [
                _find(summary, family, d, 2, M, "path") for M in (16, 64, 256)
            ]
            extended_floor = initial_floor + [_find(summary, family, d, 2, 512, "path")]
            initial_deep = [
                _find(summary, family, d, 3, M, method)
                for M in (16, 24)
                for method in ("double_path", "double_path_base")
            ]
            extended_deep = initial_deep + [
                _find(summary, family, d, 3, 32, method)
                for method in ("double_path", "double_path_base")
            ]
            if not all(row is not None for row in extended_floor + extended_deep):
                continue
            old_floor = min(initial_floor, key=lambda row: row.mean_skill)
            old_deep = min(initial_deep, key=lambda row: row.mean_skill)
            new_floor = min(extended_floor, key=lambda row: row.mean_skill)
            new_deep = min(extended_deep, key=lambda row: row.mean_skill)
            records.append(
                {
                    "family": family,
                    "dimension": d,
                    "initial_ratio": float(old_deep.mean_skill / old_floor.mean_skill),
                    "extended_ratio": float(new_deep.mean_skill / new_floor.mean_skill),
                    "extended_floor_M": int(new_floor.M),
                    "extended_n3_M": int(new_deep.M),
                    "extended_n3_method": str(new_deep.method),
                    "extended_met": bool(new_deep.mean_skill <= 0.9 * new_floor.mean_skill),
                }
            )
    return records


def _frontier(summary: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    if summary.empty:
        return pd.DataFrame()
    for (family, d), problem in summary.groupby(["pde_id", "dimension"]):
        for axis, column in (
            ("generator_calls", "mean_generator_calls"),
            ("gaussian_coordinates", "mean_gaussian_coordinates_drawn"),
        ):
            ordered = problem.sort_values([column, "mean_skill"])
            best = math.inf
            for _, row in ordered.iterrows():
                is_frontier = bool(row.mean_skill < best)
                if is_frontier:
                    best = float(row.mean_skill)
                records.append(
                    {
                        "pde_id": family,
                        "dimension": int(d),
                        "cost_axis": axis,
                        "n": int(row.n),
                        "M": int(row.M),
                        "method": row.method,
                        "mean_cost": float(row[column]),
                        "mean_skill": float(row.mean_skill),
                        "is_frontier": is_frontier,
                    }
                )
    result = pd.DataFrame(records)
    result.to_csv(RESULTS_ROOT / "frontier.csv", index=False)
    return result


def _baseline_comparisons(summary: pd.DataFrame) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for family in ("S0", "S05", "S1"):
        for d in (100, 400):
            path = _find(summary, family, d, 3, 16, "path")
            path_base = _find(summary, family, d, 3, 16, "path_base")
            double = _find(summary, family, d, 3, 16, "double_path")
            double_base = _find(summary, family, d, 3, 16, "double_path_base")
            if all(row is not None for row in (path, path_base, double, double_base)):
                floors = [_find(summary, family, d, 2, M, "path") for M in (16, 64, 256)]
                deeps = [
                    _find(summary, family, d, 3, M, method)
                    for M in (16, 24)
                    for method in ("double_path", "double_path_base")
                ]
                floor = min(floors, key=lambda row: row.mean_skill)
                deep = min(deeps, key=lambda row: row.mean_skill)
                records.append(
                    {
                        "family": family,
                        "dimension": d,
                        "source": "nonridge GPU",
                        "path_base_over_path_M16": float(path_base.mean_skill / path.mean_skill),
                        "double_base_over_double_M16": float(double_base.mean_skill / double.mean_skill),
                        "best_n3_over_n2_floor": float(deep.mean_skill / floor.mean_skill),
                    }
                )
    if BASELINE_SUMMARY.is_file():
        baseline = pd.read_csv(BASELINE_SUMMARY)
        for family in ("P1", "MR"):
            for d in (100, 400):
                problem = baseline[(baseline.pde_id == family) & (baseline.dimension == d)]
                def bfind(n: int, M: int, method: str) -> pd.Series:
                    return problem[(problem.n == n) & (problem.M == M) & (problem.method == method)].iloc[0]
                path = bfind(3, 16, "path")
                path_base = bfind(3, 16, "path_base")
                double = bfind(3, 16, "double_path")
                double_base = bfind(3, 16, "double_path_base")
                floor = problem[(problem.n == 2) & (problem.method == "path")].sort_values("mean_skill").iloc[0]
                deep = problem[
                    (problem.n == 3)
                    & problem.method.isin(["double_path", "double_path_base"])
                ].sort_values("mean_skill").iloc[0]
                records.append(
                    {
                        "family": family,
                        "dimension": d,
                        "source": "Study F CPU; ratios only, unpaired",
                        "path_base_over_path_M16": float(path_base.mean_skill / path.mean_skill),
                        "double_base_over_double_M16": float(double_base.mean_skill / double.mean_skill),
                        "best_n3_over_n2_floor": float(deep.mean_skill / floor.mean_skill),
                    }
                )
    return records


def _v1_slopes(summary: pd.DataFrame) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    selected = summary[
        (summary.n == 3)
        & summary.method.isin(["path", "path_base"])
        & np.isfinite(summary.mean_root_level1_variance)
        & (summary.mean_root_level1_variance > 0.0)
        & (summary.mean_skill > 0.0)
    ].copy()
    if selected.empty:
        return rows
    selected["v1_d_over_M"] = (
        selected.mean_root_level1_variance * selected.dimension / selected.M
    )
    for family, group in selected.groupby("pde_id"):
        if len(group) >= 2:
            slope, intercept = np.polyfit(
                np.log(group.v1_d_over_M), np.log(group.mean_skill), 1
            )
            rows.append(
                {
                    "family": family,
                    "points": len(group),
                    "log_log_slope": float(slope),
                    "log_intercept": float(intercept),
                }
            )
    return rows


def _plots(summary: pd.DataFrame, frontier: pd.DataFrame) -> None:
    if summary.empty:
        return
    FIGURES_ROOT.mkdir(parents=True, exist_ok=True)
    families = [family for family in ("S0", "S05", "S1") if family in set(summary.pde_id)]
    methods = ("raw", "path", "path_base", "double_path", "double_path_base")
    colors = dict(zip(methods, ("#9d755d", "#4c78a8", "#72b7b2", "#e45756", "#f58518")))
    fig, axes = plt.subplots(len(families), 2, figsize=(11, 3.8 * len(families)), squeeze=False)
    for row_index, family in enumerate(families):
        for column, d in enumerate((100, 400)):
            axis = axes[row_index, column]
            problem = summary[(summary.pde_id == family) & (summary.dimension == d)]
            floor = problem[(problem.n == 2) & (problem.method == "path")]
            if not floor.empty:
                axis.axhline(float(floor.mean_skill.min()), color="black", linestyle="--", label="n=2 path floor")
            for method in methods:
                data = problem[(problem.n == 3) & (problem.method == method)].sort_values("M")
                if not data.empty:
                    axis.plot(data.M, data.mean_skill, marker="o", color=colors[method], label=method)
            axis.set_yscale("log")
            axis.set_title(f"{family}, d={d}")
            axis.set_xlabel("M (n=3)")
            axis.set_ylabel("mean test skill")
            axis.grid(True, which="both", alpha=0.25)
            axis.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES_ROOT / "skill_vs_M.png", dpi=180)
    plt.close(fig)

    collapse = summary[
        (summary.n == 3)
        & summary.method.isin(["path", "path_base"])
        & np.isfinite(summary.mean_root_level1_variance)
        & (summary.mean_root_level1_variance > 0.0)
    ].copy()
    fig, axis = plt.subplots(figsize=(7, 5))
    for family in families:
        for method, marker in (("path", "o"), ("path_base", "s")):
            data = collapse[(collapse.pde_id == family) & (collapse.method == method)]
            if not data.empty:
                x = data.mean_root_level1_variance * data.dimension / data.M
                axis.scatter(x, data.mean_skill, label=f"{family} {method}", marker=marker)
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("V1 x d / M")
    axis.set_ylabel("mean test skill")
    axis.grid(True, which="both", alpha=0.25)
    axis.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES_ROOT / "v1_collapse.png", dpi=180)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(7, 5))
    for family in families:
        for M, marker in ((16, "o"), (24, "s"), (32, "^")):
            data = summary[
                (summary.pde_id == family)
                & (summary.n == 3)
                & (summary.M == M)
                & (summary.method == "double_path_base")
            ].sort_values("dimension")
            if len(data) == 2:
                axis.plot(data.dimension, data.mean_skill, marker=marker, label=f"{family}, M={M}")
    axis.set_yscale("log")
    axis.set_xlabel("dimension")
    axis.set_ylabel("double_path_base mean skill")
    axis.grid(True, which="both", alpha=0.25)
    axis.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES_ROOT / "skill_vs_dimension.png", dpi=180)
    plt.close(fig)

    for cost_axis, xlabel in (("generator_calls", "mean generator calls"), ("gaussian_coordinates", "mean Gaussian coordinates")):
        fig, axes = plt.subplots(len(families), 2, figsize=(11, 3.8 * len(families)), squeeze=False)
        for row_index, family in enumerate(families):
            for column, d in enumerate((100, 400)):
                axis = axes[row_index, column]
                data = frontier[
                    (frontier.pde_id == family)
                    & (frontier.dimension == d)
                    & (frontier.cost_axis == cost_axis)
                    & frontier.is_frontier
                ].sort_values("mean_cost")
                if not data.empty:
                    axis.plot(data.mean_cost, data.mean_skill, marker="o")
                    for _, item in data.iterrows():
                        axis.annotate(f"n{int(item.n)} M{int(item.M)}\n{item.method}", (item.mean_cost, item.mean_skill), fontsize=6)
                axis.set_xscale("log")
                axis.set_yscale("log")
                axis.set_title(f"{family}, d={d}")
                axis.set_xlabel(xlabel)
                axis.set_ylabel("mean test skill")
                axis.grid(True, which="both", alpha=0.25)
        fig.tight_layout()
        fig.savefig(FIGURES_ROOT / f"frontier_{cost_axis}.png", dpi=180)
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
    rows: pd.DataFrame,
    summary: pd.DataFrame,
    stage_a: dict[str, Any],
    s1: dict[str, Any],
    extended: list[dict[str, Any]],
    comparisons: list[dict[str, Any]],
    slopes: list[dict[str, Any]],
    validation: dict[str, Any],
    validation_frame: pd.DataFrame,
    environment: dict[str, Any],
    calibration: pd.DataFrame,
    prereg_commit: str,
    analysis_commit: str,
) -> None:
    lines = [
        "# Non-ridge spectral HJB benchmark: pre-registered report",
        "",
        "## Outcome table",
        "",
        "| Block | N-1 | N-2 | N-3 | N-4 | N-5 |",
        "|---|---|---|---|---|---|",
        "| Stage A: S0 + S05 | " + " | ".join(stage_a[key]["verdict"] for key in ("N-1", "N-2", "N-3", "N-4", "N-5")) + " |",
        "| Stage B: S1 | " + " | ".join(s1[key]["verdict"] for key in ("N-1", "N-2", "N-3", "N-4", "N-5")) + " |",
        "",
        "## Stage A criteria (verbatim)",
        "",
    ]
    for key in ("N-1", "N-2", "N-3", "N-4", "N-5"):
        result = stage_a[key]
        lines.extend(
            [
                f"### {key}",
                "",
                f"> {CRITERIA[key]}",
                "",
                f"**Verdict: {result['verdict']} ({result['wins']}/{result['total']}; required {result['required']}).**",
                "",
            ]
        )

    lines.extend(["## Stage A mechanical details", ""])
    for key, ratio_key in (("N-1", "path_over_raw"), ("N-2", "baseline_over_path"), ("N-3", "baseline_over_double"), ("N-5", "d400_over_d100")):
        lines.extend([f"### {key} comparisons", "", "| family | d | cell/M | ratio | met |", "|---|---:|---|---:|---|"])
        for item in stage_a[key]["rows"]:
            d = item.get("dimension", "400/100")
            cell = f"({item.get('n', 3)},{item.get('M', 16)})"
            lines.append(f"| {item['family']} | {d} | {cell} | {_number(item.get(ratio_key))} | {'yes' if item['met'] else 'no'} |")
        lines.append("")
    lines.extend(["### N-4 floor comparison", "", "| family | d | n=2 floor | best n=3 | ratio | verdict scope | met |", "|---|---:|---|---|---:|---|---|"])
    for item in stage_a["N-4-all-dimensions"]:
        if not item["complete"]:
            lines.append(f"| {item['family']} | {item['dimension']} | NA | NA | NA | {'yes' if item['dimension'] == 400 else 'reported only'} | NA |")
            continue
        lines.append(
            f"| {item['family']} | {item['dimension']} | M={item['floor_M']}; {_number(item['floor_skill'])} | {item['best_n3_method']} M={item['best_n3_M']}; {_number(item['best_n3_skill'])} | {_number(item['best_n3_over_floor'])} | {'yes' if item['dimension'] == 400 else 'reported only'} | {'yes' if item['met'] else 'no'} |"
        )

    lines.extend(["", "## Stage B: S1 pre-registered block", "", "The same N-1 through N-5 inequalities and thresholds are applied to S1 alone; N-3 requires all 4 S1 comparisons.", ""])
    for key in ("N-1", "N-2", "N-3", "N-4", "N-5"):
        result = s1[key]
        lines.append(f"- **{key}: {result['verdict']}** ({result['wins']}/{result['total']}; required {result['required']}).")
    if extended:
        lines.extend(["", "## Extended M=32 / M=512 floor re-evaluation", "", "| family | d | initial ratio | extended ratio | extended attaining cells | <= 0.9 |", "|---|---:|---:|---:|---|---|"])
        for item in extended:
            lines.append(
                f"| {item['family']} | {item['dimension']} | {_number(item['initial_ratio'])} | {_number(item['extended_ratio'])} | n2 M={item['extended_floor_M']}; n3 {item['extended_n3_method']} M={item['extended_n3_M']} | {'yes' if item['extended_met'] else 'no'} |"
            )

    lines.extend(["", "## Calibration (reported, no verdict)", ""])
    if not calibration.empty:
        lines.extend([
            "| PDE | d | std u* | mean f(z*) | std f(z*) | RMS norm(z*) | norm(mean z*) | heat skill | z var outside top-1 | outside top-10 | Hessian Fro | Hessian eff. rank |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ])
        for _, item in calibration.iterrows():
            lines.append(
                f"| {item.pde_id} | {int(item.dimension)} | {_number(item.std_exact_u)} | {_number(item.mean_exact_generator)} | {_number(item.std_exact_generator)} | {_number(item.rms_z_norm)} | {_number(item.mean_z_norm)} | {_number(item.heat_iterate_skill)} | {_number(item.centred_z_variance_outside_top1)} | {_number(item.centred_z_variance_outside_top10)} | {_number(item.mean_terminal_hessian_frobenius)} | {_number(item.mean_terminal_hessian_effective_rank)} |"
            )

    lines.extend(["", "## Equation validation", ""])
    for key in ("V-1", "V-2", "V-3", "V-4", "V-5", "V-6"):
        item = validation.get(key, {})
        lines.append(
            f"- **{key}: {item.get('verdict', 'NOT EVALUATED')}**; checks {item.get('checks', 0)}, failures {item.get('failures', 'NA')}."
        )
    v5 = validation_frame[validation_frame.stage == "V-5"] if not validation_frame.empty else pd.DataFrame()
    if not v5.empty:
        lines.extend(["", "The GPU rotation uses the Makhoul FFT DCT. Timings are for one 1,024-point forward DCT:", "", "| family | d | milliseconds/call | max equation relative error |", "|---|---:|---:|---:|"])
        for _, item in v5.iterrows():
            lines.append(
                f"| {item.pde_id} | {int(item.dimension)} | {_number(item.dct_milliseconds_per_call)} | {_number(item.max_relative_error)} |"
            )

    lines.extend(["", "## Ratio-level comparison with Study F (reported, unpaired)", "", "| family | d | source | path_base/path M16 | double_base/double M16 | best n3/n2 floor |", "|---|---:|---|---:|---:|---:|"])
    for item in comparisons:
        lines.append(
            f"| {item['family']} | {item['dimension']} | {item['source']} | {_number(item['path_base_over_path_M16'])} | {_number(item['double_base_over_double_M16'])} | {_number(item['best_n3_over_n2_floor'])} |"
        )
    lines.extend(["", "## V1-collapse fit (reported, no verdict)", ""])
    for item in slopes:
        lines.append(f"- {item['family']}: slope {_number(item['log_log_slope'])} from {item['points']} available points.")

    n4 = summary[(summary.n == 4) & (summary.method == "double_path_base")]
    if not n4.empty:
        lines.extend(["", "## Exploratory n=4 results (no verdict)", "", "| family | d | M | skill | generator calls | Gaussian coordinates | GPU seconds |", "|---|---:|---:|---:|---:|---:|---:|"])
        for _, item in n4.sort_values(["pde_id", "dimension", "M"]).iterrows():
            lines.append(
                f"| {item.pde_id} | {int(item.dimension)} | {int(item.M)} | {_number(item.mean_skill)} | {_number(item.mean_generator_calls)} | {_number(item.mean_gaussian_coordinates_drawn)} | {_number(item.mean_wall_time_seconds)} |"
            )

    lines.extend(["", "## Equal-cost frontiers (reported, no verdict)", "", "Attaining cells for both generator-call and Gaussian-coordinate lower envelopes are recorded in `results/nonridge/frontier.csv` and annotated in the two frontier figures."])

    failures = [f"Stage A {key}" for key in ("N-1", "N-2", "N-3", "N-4", "N-5") if stage_a[key]["verdict"] == "FAIL"]
    failures.extend(f"S1 {key}" for key in ("N-1", "N-2", "N-3", "N-4", "N-5") if s1[key]["verdict"] == "FAIL")
    curtailments = json.loads(CURTAIL_PATH.read_text(encoding="utf-8")) if CURTAIL_PATH.is_file() else {}
    lines.extend([
        "",
        "## Failures and deviations",
        "",
        f"- Mechanical criterion failures: {', '.join(failures) if failures else 'none'}.",
        "- No post-hoc thresholds were introduced; non-finite raw means are handled only by the pre-registered N-1 rule.",
        "- The registered Stage-A grid contains `path` and `path_base` at M=16 but not M=24, although the reported V1-collapse sentence names both M=16 and M=24. No unregistered M=24 rows were added; the fit uses the available registered path/baseline points and this mismatch is retained explicitly.",
        f"- Stage-A >30-minute curtailments: {len(curtailments)} cells; details are in `results/nonridge/curtailed_cells.json` if nonzero.",
        "",
        "## Provenance",
        "",
        f"- Frozen preregistration commit: `{prereg_commit}`.",
        f"- Analysis code commit: `{analysis_commit}`.",
        f"- Solver/run commits represented in rows: {', '.join(f'`{value}`' for value in sorted(set(rows.code_commit)) if isinstance(value, str)) if not rows.empty else 'none'}.",
        f"- GPU: {environment.get('gpu_name', 'NA')}; driver {environment.get('nvidia_driver', 'NA')}; CUDA {environment.get('cuda_runtime', 'NA')}; PyTorch {environment.get('torch', 'NA')}; DCT: {environment.get('dct_implementation', 'NA')}.",
        f"- Completed rows: {len(rows)}; Stage A {int((rows.stage == 'stage_a').sum()) if not rows.empty else 0}/{len(stage_a_tasks())}; S1 {int((rows.stage == 'stage_b_s1').sum()) if not rows.empty else 0}/{len(stage_b_s1_tasks())}; extension {int((rows.stage == 'stage_b_extension').sum()) if not rows.empty else 0}/{len(stage_b_extension_tasks())}; n=4 {int((rows.stage == 'stage_b_n4').sum()) if not rows.empty else 0}/{len(stage_b_n4_tasks())}.",
        f"- Non-finite state values: {int(rows.nonfinite_state_count.sum()) if not rows.empty else 0}; non-finite generator values: {int(rows.nonfinite_generator_count.sum()) if not rows.empty else 0}.",
        f"- GPU wall time represented by rows: {_number(rows.wall_time_seconds.sum() if not rows.empty else 0.0)} seconds.",
        f"- Validation: {', '.join(f'{key}={validation.get(key, {}).get("verdict", "NA")}' for key in ('V-1','V-2','V-3','V-4','V-5','V-6'))}.",
        "- Arithmetic is float64; verdicts use the frozen 960-point test split and 10 repetitions unless the pre-registered >30-minute rule is triggered.",
        "",
        "## Figures",
        "",
        "- `results/nonridge/figures/skill_vs_M.png`",
        "- `results/nonridge/figures/v1_collapse.png`",
        "- `results/nonridge/figures/skill_vs_dimension.png`",
        "- `results/nonridge/figures/frontier_generator_calls.png`",
        "- `results/nonridge/figures/frontier_gaussian_coordinates.png`",
        "",
    ])
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text("\n".join(lines), encoding="utf-8")


def analyze() -> dict[str, Any]:
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    rows = _load_rows()
    summary = _summary(rows)
    stage_a = _evaluate_block(summary, ("S0", "S05"), n3_required=7)
    s1 = _evaluate_block(summary, ("S1",), n3_required=4)
    extended = _extended_floors(summary)
    frontier = _frontier(summary)
    comparisons = _baseline_comparisons(summary)
    slopes = _v1_slopes(summary)
    validation_path = RESULTS_ROOT / "validation_summary.json"
    validation = json.loads(validation_path.read_text(encoding="utf-8")) if validation_path.is_file() else {}
    validation_csv = RESULTS_ROOT / "validation.csv"
    validation_frame = pd.read_csv(validation_csv) if validation_csv.is_file() else pd.DataFrame()
    environment_path = RESULTS_ROOT / "environment.json"
    environment = json.loads(environment_path.read_text(encoding="utf-8")) if environment_path.is_file() else {}
    calibration_path = RESULTS_ROOT / "calibration.csv"
    calibration = pd.read_csv(calibration_path) if calibration_path.is_file() else pd.DataFrame()
    prereg_commit = _git(
        "log", "--diff-filter=A", "-1", "--format=%H", "--",
        "invariant_region_mlp/experiments/nonridge/PREREGISTRATION.md",
    )
    analysis_commit = _git("rev-parse", "HEAD")
    payload = {
        "stage_a": stage_a,
        "s1_block": s1,
        "extended_floors": extended,
        "ratio_comparisons": comparisons,
        "v1_slopes": slopes,
        "completed_rows": len(rows),
        "stage_counts": rows.stage.value_counts().to_dict() if not rows.empty else {},
        "preregistration_commit": prereg_commit,
        "analysis_code_commit": analysis_commit,
        "validation": validation,
    }
    write_path = RESULTS_ROOT / "analysis_summary.json"
    write_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=True) + "\n",
        encoding="utf-8",
    )
    _plots(summary, frontier)
    _report(
        rows, summary, stage_a, s1, extended, comparisons, slopes,
        validation, validation_frame, environment, calibration, prereg_commit, analysis_commit,
    )
    print(json.dumps(payload, indent=2, sort_keys=True, allow_nan=True), flush=True)
    return payload


if __name__ == "__main__":
    analyze()
