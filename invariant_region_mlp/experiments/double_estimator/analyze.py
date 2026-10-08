"""Mechanical aggregation, verdicts, frontiers, figures, and report."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import subprocess
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .protocol import CORE_CELLS, EC_CELLS, REPETITIONS, RESULTS_ROOT, RAW_ROOT


REPORT_PATH = RESULTS_ROOT.parents[1] / "docs" / "DOUBLE_ESTIMATOR_REPORT.md"
FIGURES = RESULTS_ROOT / "figures"
BOOTSTRAP_DRAWS = 1000
BOOTSTRAP_SEED = 20261201

CRITERIA = {
    "D-1": (
        "For P1 and MR, in cells (3,6) and (4,6), |mean generator bias| of "
        "`double` <= 0.05 x that of `raw`, at every d."
    ),
    "D-2": (
        "For P1 and MR, `double` skill <= 0.1 x `raw` skill in cells (3,6) "
        "and (4,6) at every d, and `double` skill at the largest d <= 1.5 x "
        "its value at the smallest d in each of those cells."
    ),
    "D-3": (
        "`double_path` skill <= `path` skill in cells (3,6) and (4,6) for at "
        "least 75% of the (PDE, d, cell) combinations over P1 and MR."
    ),
    "D-4": (
        "For P1 and MR at d=100 and d=400, the frontier of the better of "
        "`double`/`double_path` is <= 0.8 x the `raw` frontier at >= 8 of 10 "
        "cost levels."
    ),
    "D-5": (
        "On C2, |bias(`double`)| <= 0.1 x |bias(`raw`)| for convex and flip, "
        "at cells (2,32), (3,6), every d."
    ),
    "D-6": (
        "Exploratory P4: report bias and skill of the Double-Q norm form; no verdict."
    ),
}


def _git(*args: str) -> str:
    root = Path(__file__).resolve().parents[3]
    return subprocess.run(
        ["git", *args], cwd=root, check=True, capture_output=True, text=True
    ).stdout.strip()


def _load_rows() -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for path in RAW_ROOT.rglob("rep*.json") if RAW_ROOT.exists() else []:
        try:
            row = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        rows.append({key: value for key, value in row.items() if not isinstance(value, dict)})
    if not rows:
        return pd.DataFrame()
    frame = pd.DataFrame(rows)
    return frame.sort_values(
        ["pde_id", "dimension", "n", "M", "method", "repetition"]
    ).reset_index(drop=True)


def _summary(rows: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "pde_id",
        "dimension",
        "n",
        "M",
        "method",
        "repetitions",
        "mean_skill",
        "median_skill",
        "skill_std",
        "mean_generator_bias",
        "generator_bias_std",
        "mean_generator_rmse",
        "mean_generator_calls",
        "mean_recursive_calls",
        "mean_wall_time_seconds",
        "total_nonfinite_states",
        "total_nonfinite_generators",
    ]
    if rows.empty:
        return pd.DataFrame(columns=columns)
    grouped = rows.groupby(
        ["pde_id", "dimension", "n", "M", "method"], dropna=False
    )
    result = grouped.agg(
        repetitions=("repetition", "nunique"),
        mean_skill=("skill", "mean"),
        median_skill=("skill", "median"),
        skill_std=("skill", "std"),
        mean_generator_bias=("mean_generator_bias", "mean"),
        generator_bias_std=("mean_generator_bias", "std"),
        mean_generator_rmse=("generator_rmse", "mean"),
        mean_generator_calls=("generator_calls", "mean"),
        mean_recursive_calls=("recursive_calls", "mean"),
        mean_wall_time_seconds=("wall_time_seconds", "mean"),
        total_nonfinite_states=("nonfinite_state_count", "sum"),
        total_nonfinite_generators=("nonfinite_generator_count", "sum"),
    ).reset_index()
    return result[columns]


def _lookup(
    summary: pd.DataFrame,
    pde_id: str,
    d: int,
    cell: tuple[int, int],
    method: str,
) -> pd.Series | None:
    n, M = cell
    subset = summary[
        (summary.pde_id == pde_id)
        & (summary.dimension == d)
        & (summary.n == n)
        & (summary.M == M)
        & (summary.method == method)
    ]
    if len(subset) != 1 or int(subset.iloc[0].repetitions) != REPETITIONS:
        return None
    return subset.iloc[0]


def _finite_pair(first: float, second: float) -> bool:
    return bool(np.isfinite(first) and np.isfinite(second))


def _criterion_d1(summary: pd.DataFrame) -> dict[str, Any]:
    details = []
    missing = []
    for pde_id, dimensions in (("P1", (20, 100, 400)), ("MR", (100, 400))):
        for d in dimensions:
            for cell in ((3, 6), (4, 6)):
                raw = _lookup(summary, pde_id, d, cell, "raw")
                double = _lookup(summary, pde_id, d, cell, "double")
                if raw is None or double is None:
                    missing.append([pde_id, d, *cell])
                    continue
                raw_bias = abs(float(raw.mean_generator_bias))
                double_bias = abs(float(double.mean_generator_bias))
                passed = _finite_pair(raw_bias, double_bias) and double_bias <= 0.05 * raw_bias
                details.append(
                    {
                        "pde_id": pde_id,
                        "dimension": d,
                        "cell": list(cell),
                        "raw_abs_bias": raw_bias,
                        "double_abs_bias": double_bias,
                        "ratio": double_bias / raw_bias if raw_bias > 0 else math.inf,
                        "passed": passed,
                    }
                )
    verdict = "NOT EVALUATED (COMPUTE)" if missing else (
        "PASS" if all(item["passed"] for item in details) else "FAIL"
    )
    return {"verdict": verdict, "details": details, "missing": missing}


def _criterion_d2(summary: pd.DataFrame) -> dict[str, Any]:
    point_details = []
    growth_details = []
    missing = []
    dimensions_by_pde = {"P1": (20, 100, 400), "MR": (100, 400)}
    for pde_id, dimensions in dimensions_by_pde.items():
        for cell in ((3, 6), (4, 6)):
            double_by_dimension: dict[int, float] = {}
            for d in dimensions:
                raw = _lookup(summary, pde_id, d, cell, "raw")
                double = _lookup(summary, pde_id, d, cell, "double")
                if raw is None or double is None:
                    missing.append([pde_id, d, *cell])
                    continue
                raw_skill = float(raw.mean_skill)
                double_skill = float(double.mean_skill)
                double_by_dimension[d] = double_skill
                passed = _finite_pair(raw_skill, double_skill) and double_skill <= 0.1 * raw_skill
                point_details.append(
                    {
                        "pde_id": pde_id,
                        "dimension": d,
                        "cell": list(cell),
                        "raw_skill": raw_skill,
                        "double_skill": double_skill,
                        "ratio": double_skill / raw_skill if raw_skill > 0 else math.inf,
                        "passed": passed,
                    }
                )
            if len(double_by_dimension) == len(dimensions):
                low = double_by_dimension[min(dimensions)]
                high = double_by_dimension[max(dimensions)]
                growth_passed = _finite_pair(low, high) and high <= 1.5 * low
                growth_details.append(
                    {
                        "pde_id": pde_id,
                        "cell": list(cell),
                        "smallest_d_skill": low,
                        "largest_d_skill": high,
                        "growth": high / low if low > 0 else math.inf,
                        "passed": growth_passed,
                    }
                )
    verdict = "NOT EVALUATED (COMPUTE)" if missing else (
        "PASS"
        if all(item["passed"] for item in point_details + growth_details)
        else "FAIL"
    )
    return {
        "verdict": verdict,
        "pointwise": point_details,
        "dimension_growth": growth_details,
        "missing": missing,
    }


def _criterion_d3(summary: pd.DataFrame) -> dict[str, Any]:
    details = []
    missing = []
    for pde_id, dimensions in (("P1", (20, 100, 400)), ("MR", (100, 400))):
        for d in dimensions:
            for cell in ((3, 6), (4, 6)):
                path = _lookup(summary, pde_id, d, cell, "path")
                double_path = _lookup(summary, pde_id, d, cell, "double_path")
                if path is None or double_path is None:
                    missing.append([pde_id, d, *cell])
                    continue
                baseline = float(path.mean_skill)
                candidate = float(double_path.mean_skill)
                passed = _finite_pair(baseline, candidate) and candidate <= baseline
                details.append(
                    {
                        "pde_id": pde_id,
                        "dimension": d,
                        "cell": list(cell),
                        "path_skill": baseline,
                        "double_path_skill": candidate,
                        "passed": passed,
                    }
                )
    fraction = float(np.mean([item["passed"] for item in details])) if details else math.nan
    verdict = "NOT EVALUATED (COMPUTE)" if missing else (
        "PASS" if fraction >= 0.75 else "FAIL"
    )
    return {
        "verdict": verdict,
        "success_fraction": fraction,
        "successes": int(sum(item["passed"] for item in details)),
        "combinations": len(details),
        "details": details,
        "missing": missing,
    }


def _frontier_value(points: pd.DataFrame, budget: float) -> tuple[float, str | None]:
    eligible = points[points.mean_generator_calls <= budget * (1.0 + 1e-12)]
    eligible = eligible[np.isfinite(eligible.mean_skill)]
    if eligible.empty:
        return math.nan, None
    best = eligible.loc[eligible.mean_skill.idxmin()]
    return float(best.mean_skill), str(best.method)


def _frontier_tables(
    rows: pd.DataFrame, summary: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    frontier_rows: list[dict[str, Any]] = []
    level_rows: list[dict[str, Any]] = []
    criterion_details: dict[str, Any] = {}
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    for pde_id in ("P1", "MR"):
        for d in (100, 400):
            points = summary[
                (summary.pde_id == pde_id)
                & (summary.dimension == d)
                & summary.apply(lambda row: (int(row.n), int(row.M)) in EC_CELLS, axis=1)
                & summary.method.isin(["raw", "double", "double_path"])
            ].copy()
            expected = len(EC_CELLS) * 3
            complete = len(points) == expected and bool(
                np.all(points.repetitions.to_numpy() == REPETITIONS)
            )
            for method in ("raw", "double", "double_path"):
                method_points = points[points.method == method]
                for _, point in method_points.iterrows():
                    dominated = method_points[
                        (method_points.mean_generator_calls <= point.mean_generator_calls)
                        & (method_points.mean_skill <= point.mean_skill)
                        & (
                            (method_points.mean_generator_calls < point.mean_generator_calls)
                            | (method_points.mean_skill < point.mean_skill)
                        )
                    ]
                    frontier_rows.append(
                        {
                            "pde_id": pde_id,
                            "dimension": d,
                            "method": method,
                            "n": int(point.n),
                            "M": int(point.M),
                            "mean_generator_calls": float(point.mean_generator_calls),
                            "mean_wall_time_seconds": float(point.mean_wall_time_seconds),
                            "mean_skill": float(point.mean_skill),
                            "is_frontier": dominated.empty,
                        }
                    )
            key = f"{pde_id}_d{d}"
            if not complete:
                criterion_details[key] = {
                    "verdict": "NOT EVALUATED (COMPUTE)",
                    "available_points": len(points),
                    "expected_points": expected,
                }
                continue

            raw_points = points[points.method == "raw"]
            double_points = points[points.method.isin(["double", "double_path"])]
            lower = max(
                float(raw_points.mean_generator_calls.min()),
                float(double_points.mean_generator_calls.min()),
            )
            upper = min(
                float(raw_points.mean_generator_calls.max()),
                float(double_points.mean_generator_calls.max()),
            )
            if not (np.isfinite(lower) and np.isfinite(upper) and 0 < lower <= upper):
                criterion_details[key] = {
                    "verdict": "FAIL",
                    "reason": "frontiers have no finite overlapping cost interval",
                }
                continue
            levels = np.geomspace(lower, upper, 10)
            main: list[dict[str, Any]] = []
            for index, level in enumerate(levels):
                raw_skill, _ = _frontier_value(raw_points, level)
                double_skill, winner = _frontier_value(double_points, level)
                ratio = double_skill / raw_skill if raw_skill > 0 else math.inf
                main.append(
                    {
                        "level_index": index,
                        "cost_level": float(level),
                        "raw_skill": raw_skill,
                        "double_skill": double_skill,
                        "double_method": winner,
                        "ratio": ratio,
                        "passed": bool(np.isfinite(ratio) and ratio <= 0.8),
                    }
                )

            boot_ratios = np.full((BOOTSTRAP_DRAWS, 10), np.nan)
            problem_rows = rows[
                (rows.pde_id == pde_id)
                & (rows.dimension == d)
                & rows.apply(lambda row: (int(row.n), int(row.M)) in EC_CELLS, axis=1)
                & rows.method.isin(["raw", "double", "double_path"])
            ]
            for draw in range(BOOTSTRAP_DRAWS):
                sampled_reps = rng.integers(0, REPETITIONS, REPETITIONS)
                sampled = pd.concat(
                    [problem_rows[problem_rows.repetition == rep] for rep in sampled_reps],
                    ignore_index=True,
                )
                boot = sampled.groupby(["n", "M", "method"], as_index=False).agg(
                    mean_skill=("skill", "mean"),
                    mean_generator_calls=("generator_calls", "mean"),
                )
                raw_boot = boot[boot.method == "raw"]
                double_boot = boot[boot.method.isin(["double", "double_path"])]
                for index, level in enumerate(levels):
                    raw_skill, _ = _frontier_value(raw_boot, level)
                    double_skill, _ = _frontier_value(double_boot, level)
                    if np.isfinite(raw_skill) and np.isfinite(double_skill) and raw_skill > 0:
                        boot_ratios[draw, index] = double_skill / raw_skill
            for item in main:
                values = boot_ratios[:, item["level_index"]]
                finite = values[np.isfinite(values)]
                level_rows.append(
                    {
                        "pde_id": pde_id,
                        "dimension": d,
                        **item,
                        "bootstrap_ratio_low": (
                            float(np.quantile(finite, 0.025)) if len(finite) else math.nan
                        ),
                        "bootstrap_ratio_high": (
                            float(np.quantile(finite, 0.975)) if len(finite) else math.nan
                        ),
                        "bootstrap_win_probability": (
                            float(np.mean(finite <= 0.8)) if len(finite) else math.nan
                        ),
                    }
                )
            wins = int(sum(item["passed"] for item in main))
            criterion_details[key] = {
                "verdict": "PASS" if wins >= 8 else "FAIL",
                "wins": wins,
                "levels": main,
            }
    frontier = pd.DataFrame(
        frontier_rows,
        columns=[
            "pde_id",
            "dimension",
            "method",
            "n",
            "M",
            "mean_generator_calls",
            "mean_wall_time_seconds",
            "mean_skill",
            "is_frontier",
        ],
    )
    levels = pd.DataFrame(
        level_rows,
        columns=[
            "pde_id",
            "dimension",
            "level_index",
            "cost_level",
            "raw_skill",
            "double_skill",
            "double_method",
            "ratio",
            "passed",
            "bootstrap_ratio_low",
            "bootstrap_ratio_high",
            "bootstrap_win_probability",
        ],
    )
    all_complete = len(criterion_details) == 4 and all(
        detail["verdict"] != "NOT EVALUATED (COMPUTE)"
        for detail in criterion_details.values()
    )
    verdict = (
        "NOT EVALUATED (COMPUTE)"
        if not all_complete
        else (
            "PASS"
            if all(detail["verdict"] == "PASS" for detail in criterion_details.values())
            else "FAIL"
        )
    )
    return frontier, levels, {"verdict": verdict, "problems": criterion_details}


def _criterion_d5(summary: pd.DataFrame) -> dict[str, Any]:
    details = []
    missing = []
    for pde_id in ("C2-convex", "C2-flip"):
        for d in (20, 100, 400):
            for cell in ((2, 32), (3, 6)):
                raw = _lookup(summary, pde_id, d, cell, "raw")
                double = _lookup(summary, pde_id, d, cell, "double")
                if raw is None or double is None:
                    missing.append([pde_id, d, *cell])
                    continue
                raw_bias = abs(float(raw.mean_generator_bias))
                double_bias = abs(float(double.mean_generator_bias))
                passed = _finite_pair(raw_bias, double_bias) and double_bias <= 0.1 * raw_bias
                details.append(
                    {
                        "pde_id": pde_id,
                        "dimension": d,
                        "cell": list(cell),
                        "raw_abs_bias": raw_bias,
                        "double_abs_bias": double_bias,
                        "ratio": double_bias / raw_bias if raw_bias > 0 else math.inf,
                        "passed": passed,
                    }
                )
    verdict = "NOT EVALUATED (COMPUTE)" if missing else (
        "PASS" if all(item["passed"] for item in details) else "FAIL"
    )
    return {"verdict": verdict, "details": details, "missing": missing}


def _criterion_d6(summary: pd.DataFrame) -> dict[str, Any]:
    details = []
    missing = []
    for d in (20, 100):
        for cell in CORE_CELLS:
            for method in ("double", "double_path"):
                item = _lookup(summary, "P4", d, cell, method)
                if item is None:
                    missing.append([d, *cell, method])
                    continue
                details.append(
                    {
                        "dimension": d,
                        "cell": list(cell),
                        "method": method,
                        "skill": float(item.mean_skill),
                        "mean_generator_bias": float(item.mean_generator_bias),
                        "mean_generator_rmse": float(item.mean_generator_rmse),
                    }
                )
    return {
        "verdict": "NOT EVALUATED (COMPUTE)" if missing else "EXPLORATORY",
        "details": details,
        "missing": missing,
    }


def _plot_skill(summary: pd.DataFrame) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    methods = ["raw", "path", "double", "double_path", "box", "sub_box", "oracle_state"]
    colors = {
        "raw": "#444444",
        "path": "#4c78a8",
        "double": "#f58518",
        "double_path": "#e45756",
        "box": "#54a24b",
        "sub_box": "#54a24b",
        "oracle_state": "#b279a2",
    }
    fig, axes = plt.subplots(2, 4, figsize=(16, 7), sharey=False)
    for row_index, pde_id in enumerate(("P1", "MR")):
        for column, cell in enumerate(CORE_CELLS):
            axis = axes[row_index, column]
            for method in methods:
                data = summary[
                    (summary.pde_id == pde_id)
                    & (summary.n == cell[0])
                    & (summary.M == cell[1])
                    & (summary.method == method)
                ].sort_values("dimension")
                if data.empty:
                    continue
                axis.plot(
                    data.dimension,
                    data.mean_skill,
                    marker="o",
                    label=method,
                    color=colors[method],
                )
            axis.set_yscale("log")
            axis.set_title(f"{pde_id}, (n,M)={cell}")
            axis.set_xlabel("dimension")
            axis.grid(True, which="both", alpha=0.25)
            if column == 0:
                axis.set_ylabel("mean test skill")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=7, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(FIGURES / "skill_vs_d_per_cell.png", dpi=180)
    plt.close(fig)


def _plot_bias(summary: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    for row_index, pde_id in enumerate(("P1", "MR")):
        for column, cell in enumerate(((3, 6), (4, 6))):
            axis = axes[row_index, column]
            for method, color in (("raw", "#444444"), ("double", "#f58518")):
                data = summary[
                    (summary.pde_id == pde_id)
                    & (summary.n == cell[0])
                    & (summary.M == cell[1])
                    & (summary.method == method)
                ].sort_values("dimension")
                if not data.empty:
                    axis.plot(
                        data.dimension,
                        data.mean_generator_bias,
                        marker="o",
                        label=method,
                        color=color,
                    )
            axis.axhline(0.0, color="black", linewidth=0.8)
            axis.set_title(f"{pde_id}, (n,M)={cell}")
            axis.set_xlabel("dimension")
            axis.set_ylabel("mean generator bias")
            axis.grid(True, alpha=0.25)
            axis.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURES / "generator_bias_vs_d.png", dpi=180)
    plt.close(fig)


def _plot_frontiers(frontier: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for row_index, pde_id in enumerate(("P1", "MR")):
        for column, d in enumerate((100, 400)):
            axis = axes[row_index, column]
            for method, color in (
                ("raw", "#444444"),
                ("double", "#f58518"),
                ("double_path", "#e45756"),
            ):
                data = frontier[
                    (frontier.pde_id == pde_id)
                    & (frontier.dimension == d)
                    & (frontier.method == method)
                ].sort_values("mean_generator_calls")
                if data.empty:
                    continue
                axis.plot(
                    data.mean_generator_calls,
                    data.mean_skill,
                    marker="o",
                    linestyle="-",
                    alpha=0.35,
                    color=color,
                )
                selected = data[data.is_frontier]
                axis.scatter(
                    selected.mean_generator_calls,
                    selected.mean_skill,
                    label=method,
                    color=color,
                    s=35,
                )
            axis.set_xscale("log")
            axis.set_yscale("log")
            axis.set_title(f"{pde_id}, d={d}")
            axis.set_xlabel("mean generator calls")
            axis.set_ylabel("mean test skill")
            axis.grid(True, which="both", alpha=0.25)
            axis.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURES / "equal_cost_frontiers.png", dpi=180)
    plt.close(fig)


def _format_number(value: Any) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return "NA"
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.6g}"
    return str(value)


def _report(
    rows: pd.DataFrame,
    summary: pd.DataFrame,
    verdicts: dict[str, Any],
    prereg_commit: str,
    analysis_commit: str,
) -> None:
    lines = [
        "# Double estimator for the generator: pre-registered report",
        "",
        "## Outcome table",
        "",
        "| Criterion | Verdict |",
        "|---|---|",
    ]
    for criterion in CRITERIA:
        lines.append(f"| {criterion} | {verdicts[criterion]['verdict']} |")
    lines.extend(["", "## Criteria (verbatim)", ""])
    for criterion, text in CRITERIA.items():
        result = verdicts[criterion]
        lines.extend(
            [
                f"### {criterion}",
                "",
                f"> {text}",
                "",
                f"**Verdict: {result['verdict']}**",
                "",
                "```json",
                json.dumps(result, indent=2, sort_keys=True, allow_nan=True),
                "```",
                "",
            ]
        )

    deep = summary[
        summary.pde_id.isin(["P1", "MR"])
        & summary.apply(lambda row: (int(row.n), int(row.M)) in ((3, 6), (4, 6)), axis=1)
        & summary.method.isin(["raw", "path", "double", "double_path"])
    ]
    lines.extend(
        [
            "## Bias--variance separation",
            "",
            "The table reports the requested generator RMSE and across-repetition skill spread. Missing rows mean that the corresponding priority block was not completed.",
            "",
            "| PDE | d | cell | method | mean bias | generator RMSE | mean skill | skill SD |",
            "|---|---:|---|---|---:|---:|---:|---:|",
        ]
    )
    for _, row in deep.sort_values(["pde_id", "dimension", "n", "M", "method"]).iterrows():
        lines.append(
            "| {pde_id} | {dimension} | ({n},{M}) | {method} | {bias} | {rmse} | {skill} | {spread} |".format(
                pde_id=row.pde_id,
                dimension=int(row.dimension),
                n=int(row.n),
                M=int(row.M),
                method=row.method,
                bias=_format_number(row.mean_generator_bias),
                rmse=_format_number(row.mean_generator_rmse),
                skill=_format_number(row.mean_skill),
                spread=_format_number(row.skill_std),
            )
        )

    completed_pdes = sorted(rows.pde_id.unique().tolist()) if not rows.empty else []
    missing_criteria = [
        criterion
        for criterion, result in verdicts.items()
        if result["verdict"] == "NOT EVALUATED (COMPUTE)"
    ]
    code_commits = sorted(rows.code_commit.dropna().unique().tolist()) if not rows.empty else []
    lines.extend(
        [
            "",
            "## Failures and deviations",
            "",
            "- No failed criterion is suppressed. Criteria without every required cell and all 10 repetitions are marked `NOT EVALUATED (COMPUTE)`.",
            f"- Criteria currently not evaluated: {', '.join(missing_criteria) if missing_criteria else 'none'}.",
            "- The two recursive estimates are used separately only inside the generator. Their arithmetic average is the state used for diagnostics and for every non-generator use, as pre-registered.",
            "- The primary RNG uses the exact registered six-part SeedSequence. Auxiliary trees use a separately seeded recursive spawn tree, so their extra draws cannot shift the primary stream.",
            "- For P4, `lambda` in the registered Double-Q formula is implemented as the equation's effective z-coordinate coefficient `lambda_f/sigma`, because this codebase stores `z=sigma*grad(u)`.",
            "- The MR `centre` control retains the established ambient box-hull centre from `effdim_mlp.py`; the certified reference is `sub_box`.",
            "",
            "## Provenance and totals",
            "",
            f"- Frozen preregistration commit: `{prereg_commit}`.",
            f"- Analysis code commit: `{analysis_commit}`.",
            f"- Result-producing code commits: {', '.join(f'`{value}`' for value in code_commits) if code_commits else 'none yet'}.",
            f"- Completed PDEs: {', '.join(completed_pdes) if completed_pdes else 'none'}.",
            f"- Completed rows: {len(rows)} / 6080 unique registered tasks.",
            f"- Aggregate worker wall time: {_format_number(rows.wall_time_seconds.sum() if not rows.empty else 0.0)} seconds.",
            f"- Non-finite test-state values: {int(rows.nonfinite_state_count.sum()) if not rows.empty else 0}.",
            f"- Non-finite generator values: {int(rows.nonfinite_generator_count.sum()) if not rows.empty else 0}.",
            "- Test points: 1,200 fixed points per PDE/dimension with a fixed 20% validation split; verdicts use the 960-point test subset.",
            "- Arithmetic: float64 throughout.",
            "",
            "## Figures",
            "",
            "- `results/double_estimator/figures/skill_vs_d_per_cell.png`",
            "- `results/double_estimator/figures/generator_bias_vs_d.png`",
            "- `results/double_estimator/figures/equal_cost_frontiers.png`",
            "",
        ]
    )
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text("\n".join(lines), encoding="utf-8")


def analyze() -> dict[str, Any]:
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    rows = _load_rows()
    summary = _summary(rows)
    summary.to_csv(RESULTS_ROOT / "summary.csv", index=False)
    frontier, levels, d4 = _frontier_tables(rows, summary)
    frontier.to_csv(RESULTS_ROOT / "frontier.csv", index=False)
    levels.to_csv(RESULTS_ROOT / "frontier_levels.csv", index=False)
    verdicts = {
        "D-1": _criterion_d1(summary),
        "D-2": _criterion_d2(summary),
        "D-3": _criterion_d3(summary),
        "D-4": d4,
        "D-5": _criterion_d5(summary),
        "D-6": _criterion_d6(summary),
    }
    prereg_commit = _git(
        "log",
        "--diff-filter=A",
        "-1",
        "--format=%H",
        "--",
        "invariant_region_mlp/experiments/double_estimator/PREREGISTRATION.md",
    )
    analysis_commit = _git("rev-parse", "HEAD")
    payload = {
        "criteria": verdicts,
        "completed_rows": len(rows),
        "registered_rows": 6080,
        "bootstrap_draws": BOOTSTRAP_DRAWS,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "preregistration_commit": prereg_commit,
        "analysis_code_commit": analysis_commit,
    }
    (RESULTS_ROOT / "analysis_summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=True) + "\n",
        encoding="utf-8",
    )
    _plot_skill(summary)
    _plot_bias(summary)
    _plot_frontiers(frontier)
    _report(rows, summary, verdicts, prereg_commit, analysis_commit)
    print(json.dumps(payload, indent=2, sort_keys=True, allow_nan=True))
    return payload


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.parse_args()
    analyze()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
