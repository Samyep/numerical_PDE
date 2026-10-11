"""Mechanical Study-F criteria, three cost frontiers, figures, and report."""

from __future__ import annotations

import json
import math
from pathlib import Path
import subprocess
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .protocol import EXPECTED_ROWS, PROBLEMS, RESULTS_ROOT, REPETITIONS


PACKAGE_ROOT = Path(__file__).resolve().parents[2]
REPORT_PATH = PACKAGE_ROOT / "docs" / "BASELINE_N3_REPORT.md"
FIGURES_ROOT = RESULTS_ROOT / "figures"
BOOTSTRAP_DRAWS = 1000
BOOTSTRAP_SEED = 20261201
COST_AXES = {
    "generator_calls": "mean_generator_calls",
    "gaussian_coordinates": "mean_gaussian_coordinates_drawn",
    "wall_time": "mean_wall_time_seconds",
}
CRITERIA = {
    "F-1": (
        "**F-1 (baseline repairs `path` at n=3):** `path_base` mean skill <= 0.5 x `path` mean skill in at least 10 of the 12 "
        "(PDE, d, M) n=3 combinations."
    ),
    "F-2": (
        "**F-2 (baseline complements double):** `double_path_base` mean skill <= `double_path` mean skill in at least 10 of the "
        "12 combinations."
    ),
    "F-3": (
        "**F-3 (n=3 below the n=2 floor at d=400):** for each PDE at d=400, the minimum over n=3 cells and methods "
        "{`double_path`, `double_path_base`} of mean skill is <= 0.9 x the minimum over all n=2 `path` cells "
        "(M = 8 ... 512, including reused rows). PASS if it holds for both PDEs. Report the same quantity at d=100 (no verdict)."
    ),
    "F-4": (
        "**F-4 (reported, no verdict):** frontiers (lower envelope of mean skill vs cost) of the n=2-only `path` family and of "
        "the combined family {`path`, `path_base`, `double_path`, `double_path_base`} over n in {2,3}, on each of the three cost "
        "axes, at 10 log-spaced cost levels where both exist; the attaining cell at each level; 1,000-draw bootstrap over reps."
    ),
}


NUMERIC_COLUMNS = (
    "dimension",
    "n",
    "M",
    "repetition",
    "skill",
    "root_gradient_skill",
    "generator_calls",
    "gaussian_vectors_drawn",
    "gaussian_coordinates_drawn",
    "root_level1_gradient_integrand_variance",
    "wall_time_seconds",
    "nonfinite_state_count",
    "nonfinite_generator_count",
)


def _git(*args: str) -> str:
    root = Path(__file__).resolve().parents[3]
    return subprocess.run(
        ["git", *args], cwd=root, check=True, capture_output=True, text=True
    ).stdout.strip()


def _load_rows() -> pd.DataFrame:
    path = RESULTS_ROOT / "rows.csv"
    if not path.is_file():
        return pd.DataFrame()
    rows = pd.read_csv(path)
    for column in NUMERIC_COLUMNS:
        if column in rows:
            rows[column] = pd.to_numeric(rows[column], errors="coerce")
    return rows


def _summary(rows: pd.DataFrame) -> pd.DataFrame:
    if rows.empty:
        return pd.DataFrame()
    return (
        rows.groupby(["pde_id", "dimension", "n", "M", "method"], as_index=False)
        .agg(
            repetitions=("repetition", "nunique"),
            mean_skill=("skill", "mean"),
            std_skill=("skill", "std"),
            mean_root_gradient_skill=("root_gradient_skill", "mean"),
            mean_generator_calls=("generator_calls", "mean"),
            mean_gaussian_vectors_drawn=("gaussian_vectors_drawn", "mean"),
            mean_gaussian_coordinates_drawn=("gaussian_coordinates_drawn", "mean"),
            mean_root_level1_gradient_integrand_variance=(
                "root_level1_gradient_integrand_variance",
                "mean",
            ),
            mean_wall_time_seconds=("wall_time_seconds", "mean"),
            nonfinite_state_count=("nonfinite_state_count", "sum"),
            nonfinite_generator_count=("nonfinite_generator_count", "sum"),
        )
        .sort_values(["pde_id", "dimension", "n", "M", "method"])
        .reset_index(drop=True)
    )


def _lookup(summary: pd.DataFrame, pde: str, d: int, n: int, M: int, method: str) -> pd.Series | None:
    rows = summary[
        (summary.pde_id == pde)
        & (summary.dimension == d)
        & (summary.n == n)
        & (summary.M == M)
        & (summary.method == method)
    ]
    return None if rows.empty else rows.iloc[0]


def _criteria(summary: pd.DataFrame) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    detail: list[dict[str, Any]] = []
    f1_wins = 0
    f2_wins = 0
    complete = True
    for pde, d in PROBLEMS:
        for M in (16, 24, 32):
            path = _lookup(summary, pde, d, 3, M, "path")
            path_base = _lookup(summary, pde, d, 3, M, "path_base")
            double_path = _lookup(summary, pde, d, 3, M, "double_path")
            double_base = _lookup(summary, pde, d, 3, M, "double_path_base")
            cell_complete = all(
                row is not None and int(row.repetitions) == REPETITIONS
                for row in (path, path_base, double_path, double_base)
            )
            complete = complete and cell_complete
            path_ratio = (
                float(path_base.mean_skill / path.mean_skill)
                if cell_complete and float(path.mean_skill) > 0.0
                else math.nan
            )
            double_ratio = (
                float(double_base.mean_skill / double_path.mean_skill)
                if cell_complete and float(double_path.mean_skill) > 0.0
                else math.nan
            )
            f1 = bool(cell_complete and np.isfinite(path_ratio) and path_ratio <= 0.5)
            f2 = bool(cell_complete and np.isfinite(double_ratio) and double_ratio <= 1.0)
            f1_wins += int(f1)
            f2_wins += int(f2)
            detail.append(
                {
                    "pde_id": pde,
                    "dimension": d,
                    "M": M,
                    "complete": cell_complete,
                    "path_skill": float(path.mean_skill) if path is not None else math.nan,
                    "path_base_skill": float(path_base.mean_skill) if path_base is not None else math.nan,
                    "path_base_over_path": path_ratio,
                    "F1_met": f1,
                    "double_path_skill": float(double_path.mean_skill) if double_path is not None else math.nan,
                    "double_path_base_skill": float(double_base.mean_skill) if double_base is not None else math.nan,
                    "double_path_base_over_double_path": double_ratio,
                    "F2_met": f2,
                }
            )

    f3_detail: list[dict[str, Any]] = []
    f3_d400: list[bool] = []
    for pde, d in PROBLEMS:
        problem = summary[(summary.pde_id == pde) & (summary.dimension == d)]
        floor = problem[(problem.n == 2) & (problem.method == "path")]
        deep = problem[
            (problem.n == 3)
            & problem.method.isin(["double_path", "double_path_base"])
        ]
        cell_complete = (
            len(floor) == 7
            and len(deep) == 6
            and bool((floor.repetitions == REPETITIONS).all())
            and bool((deep.repetitions == REPETITIONS).all())
        )
        complete = complete and cell_complete
        floor_best = floor.loc[floor.mean_skill.idxmin()] if not floor.empty else None
        deep_best = deep.loc[deep.mean_skill.idxmin()] if not deep.empty else None
        ratio = (
            float(deep_best.mean_skill / floor_best.mean_skill)
            if cell_complete and floor_best is not None and float(floor_best.mean_skill) > 0.0
            else math.nan
        )
        met = bool(cell_complete and np.isfinite(ratio) and ratio <= 0.9)
        if d == 400:
            f3_d400.append(met)
        f3_detail.append(
            {
                "pde_id": pde,
                "dimension": d,
                "complete": cell_complete,
                "n2_floor_skill": float(floor_best.mean_skill) if floor_best is not None else math.nan,
                "n2_floor_M": int(floor_best.M) if floor_best is not None else None,
                "n3_best_skill": float(deep_best.mean_skill) if deep_best is not None else math.nan,
                "n3_best_method": str(deep_best.method) if deep_best is not None else None,
                "n3_best_M": int(deep_best.M) if deep_best is not None else None,
                "n3_over_n2_floor": ratio,
                "threshold": 0.9,
                "criterion_met": met if d == 400 else None,
                "verdict_scope": d == 400,
            }
        )

    if not complete:
        f1_verdict = f2_verdict = f3_verdict = "NOT EVALUATED (COMPUTE)"
    else:
        f1_verdict = "PASS" if f1_wins >= 10 else "FAIL"
        f2_verdict = "PASS" if f2_wins >= 10 else "FAIL"
        f3_verdict = "PASS" if len(f3_d400) == 2 and all(f3_d400) else "FAIL"
    criteria = {
        "F-1": {"verdict": f1_verdict, "wins": f1_wins, "required": 10, "total": 12},
        "F-2": {"verdict": f2_verdict, "wins": f2_wins, "required": 10, "total": 12},
        "F-3": {"verdict": f3_verdict, "problems": f3_detail},
        "F-4": {"verdict": "REPORTED" if complete else "NOT EVALUATED (COMPUTE)"},
    }
    return criteria, detail


def _family_rows(summary: pd.DataFrame, family: str) -> pd.DataFrame:
    if family == "n2_path":
        return summary[(summary.n == 2) & (summary.method == "path")]
    return summary[
        summary.method.isin(["path", "path_base", "double_path", "double_path_base"])
        & summary.n.isin([2, 3])
    ]


def _candidate(points: pd.DataFrame, cost_column: str, budget: float) -> pd.Series | None:
    eligible = points[
        np.isfinite(points[cost_column]) & (points[cost_column] <= budget)
    ]
    if eligible.empty:
        return None
    return eligible.sort_values(
        ["mean_skill", cost_column, "n", "M", "method"], kind="stable"
    ).iloc[0]


def _frontiers(
    rows: pd.DataFrame, summary: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    point_records: list[dict[str, Any]] = []
    level_records: list[dict[str, Any]] = []
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    for pde, d in PROBLEMS:
        problem_summary = summary[
            (summary.pde_id == pde) & (summary.dimension == d)
        ].copy()
        problem_rows = rows[(rows.pde_id == pde) & (rows.dimension == d)]
        cell_keys = [
            (int(row.n), int(row.M), str(row.method))
            for _, row in problem_summary.iterrows()
        ]
        key_to_index = {key: index for index, key in enumerate(cell_keys)}
        matrices = {
            "mean_skill": np.full((len(cell_keys), REPETITIONS), np.nan),
            "mean_generator_calls": np.full((len(cell_keys), REPETITIONS), np.nan),
            "mean_gaussian_coordinates_drawn": np.full(
                (len(cell_keys), REPETITIONS), np.nan
            ),
            "mean_wall_time_seconds": np.full((len(cell_keys), REPETITIONS), np.nan),
        }
        source_columns = {
            "mean_skill": "skill",
            "mean_generator_calls": "generator_calls",
            "mean_gaussian_coordinates_drawn": "gaussian_coordinates_drawn",
            "mean_wall_time_seconds": "wall_time_seconds",
        }
        for _, row in problem_rows.iterrows():
            key = (int(row.n), int(row.M), str(row.method))
            repetition = int(row.repetition)
            if key not in key_to_index or not (0 <= repetition < REPETITIONS):
                continue
            index = key_to_index[key]
            for destination, source in source_columns.items():
                matrices[destination][index, repetition] = float(row[source])
        bootstrap_ready = bool(cell_keys) and all(
            np.all(np.isfinite(matrix)) for matrix in matrices.values()
        )
        sampled_repetitions = rng.integers(
            0, REPETITIONS, (BOOTSTRAP_DRAWS, REPETITIONS)
        )
        boot_skill = (
            np.mean(matrices["mean_skill"][:, sampled_repetitions], axis=2).T
            if bootstrap_ready
            else np.empty((0, len(cell_keys)))
        )
        n2_mask = np.array(
            [n == 2 and method == "path" for n, _M, method in cell_keys],
            dtype=bool,
        )
        combined_mask = np.array(
            [
                n in (2, 3)
                and method in {"path", "path_base", "double_path", "double_path_base"}
                for n, _M, method in cell_keys
            ],
            dtype=bool,
        )
        for axis, cost_column in COST_AXES.items():
            family_points: dict[str, pd.DataFrame] = {}
            for family in ("n2_path", "combined"):
                points = _family_rows(problem_summary, family).copy()
                points = points[np.isfinite(points[cost_column]) & np.isfinite(points.mean_skill)]
                family_points[family] = points
                for _, row in points.iterrows():
                    dominated = bool(
                        np.any(
                            (points[cost_column] <= float(row[cost_column]))
                            & (points.mean_skill < float(row.mean_skill))
                        )
                    )
                    point_records.append(
                        {
                            "pde_id": pde,
                            "dimension": d,
                            "cost_axis": axis,
                            "family": family,
                            "n": int(row.n),
                            "M": int(row.M),
                            "method": str(row.method),
                            "mean_skill": float(row.mean_skill),
                            "mean_cost": float(row[cost_column]),
                            "is_frontier": not dominated,
                        }
                    )
            low = max(
                float(family_points[family][cost_column].min())
                for family in family_points
                if not family_points[family].empty
            )
            high = min(
                float(family_points[family][cost_column].max())
                for family in family_points
                if not family_points[family].empty
            )
            levels = np.geomspace(low, high, 10) if low > 0.0 and high >= low else np.full(10, np.nan)
            base_records: list[dict[str, Any]] = []
            for index, budget in enumerate(levels):
                selections = {
                    family: _candidate(points, cost_column, float(budget))
                    for family, points in family_points.items()
                }
                record: dict[str, Any] = {
                    "pde_id": pde,
                    "dimension": d,
                    "cost_axis": axis,
                    "level_index": index,
                    "cost_level": float(budget),
                }
                for family, selected in selections.items():
                    prefix = family
                    record[f"{prefix}_skill"] = (
                        float(selected.mean_skill) if selected is not None else math.nan
                    )
                    record[f"{prefix}_method"] = (
                        str(selected.method) if selected is not None else None
                    )
                    record[f"{prefix}_n"] = int(selected.n) if selected is not None else None
                    record[f"{prefix}_M"] = int(selected.M) if selected is not None else None
                    record[f"{prefix}_cost"] = (
                        float(selected[cost_column]) if selected is not None else math.nan
                    )
                n2_skill = record["n2_path_skill"]
                combined_skill = record["combined_skill"]
                record["combined_over_n2"] = (
                    combined_skill / n2_skill
                    if np.isfinite(n2_skill) and n2_skill > 0.0
                    else math.nan
                )
                base_records.append(record)

            bootstrap = np.full((BOOTSTRAP_DRAWS, 10), np.nan, dtype=np.float64)
            if bootstrap_ready:
                boot_cost = np.mean(
                    matrices[cost_column][:, sampled_repetitions], axis=2
                ).T
                for index, budget in enumerate(levels):
                    family_skills: dict[str, np.ndarray] = {}
                    for family, mask in (
                        ("n2_path", n2_mask),
                        ("combined", combined_mask),
                    ):
                        eligible = (boot_cost <= float(budget)) & mask[None, :]
                        candidates = np.where(eligible, boot_skill, np.inf)
                        selected = np.min(candidates, axis=1)
                        selected[~np.isfinite(selected)] = np.nan
                        family_skills[family] = selected
                    valid = (
                        np.isfinite(family_skills["n2_path"])
                        & (family_skills["n2_path"] > 0.0)
                        & np.isfinite(family_skills["combined"])
                    )
                    bootstrap[valid, index] = (
                        family_skills["combined"][valid]
                        / family_skills["n2_path"][valid]
                    )
            for record in base_records:
                finite = bootstrap[:, int(record["level_index"])]
                finite = finite[np.isfinite(finite)]
                record["combined_over_n2_bootstrap_low"] = (
                    float(np.quantile(finite, 0.025)) if len(finite) else math.nan
                )
                record["combined_over_n2_bootstrap_high"] = (
                    float(np.quantile(finite, 0.975)) if len(finite) else math.nan
                )
                level_records.append(record)
    return pd.DataFrame(point_records), pd.DataFrame(level_records)


def _predictions(summary: pd.DataFrame) -> dict[str, Any]:
    oracle_rows: list[dict[str, Any]] = []
    variance_rows: list[dict[str, Any]] = []
    scaling_rows: list[dict[str, Any]] = []
    for pde, d in PROBLEMS:
        for M in (16, 24, 32):
            oracle = _lookup(summary, pde, d, 3, M, "path_base_oracle")
            base = _lookup(summary, pde, d, 3, M, "path_base")
            path = _lookup(summary, pde, d, 3, M, "path")
            if oracle is not None and base is not None:
                oracle_rows.append(
                    {
                        "pde_id": pde,
                        "dimension": d,
                        "M": M,
                        "oracle_skill": float(oracle.mean_skill),
                        "baseline_skill": float(base.mean_skill),
                        "oracle_over_baseline": float(oracle.mean_skill / base.mean_skill),
                        "prediction_met": bool(oracle.mean_skill <= base.mean_skill),
                    }
                )
            if path is not None:
                scaling_rows.append(
                    {
                        "pde_id": pde,
                        "dimension": d,
                        "M": M,
                        "path_skill": float(path.mean_skill),
                        "M_times_skill_over_d": float(M * path.mean_skill / d),
                    }
                )
            if path is not None and base is not None:
                raw_variance = float(path.mean_root_level1_gradient_integrand_variance)
                base_variance = float(base.mean_root_level1_gradient_integrand_variance)
                if np.isfinite(raw_variance) and np.isfinite(base_variance):
                    variance_rows.append(
                        {
                            "pde_id": pde,
                            "dimension": d,
                            "M": M,
                            "path_variance": raw_variance,
                            "path_base_variance": base_variance,
                            "baseline_over_path": (
                                base_variance / raw_variance if raw_variance > 0.0 else math.inf
                            ),
                        }
                    )
    return {
        "oracle": oracle_rows,
        "root_level1_variance": variance_rows,
        "path_scaling": scaling_rows,
    }


def _plot_skill(summary: pd.DataFrame) -> None:
    FIGURES_ROOT.mkdir(parents=True, exist_ok=True)
    methods = ["path", "path_base", "path_base_oracle", "double_path", "double_path_base"]
    colors = {
        "path": "#4c78a8",
        "path_base": "#72b7b2",
        "path_base_oracle": "#54a24b",
        "double_path": "#e45756",
        "double_path_base": "#f58518",
    }
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for row_index, pde in enumerate(("P1", "MR")):
        for column, d in enumerate((100, 400)):
            axis = axes[row_index, column]
            problem = summary[(summary.pde_id == pde) & (summary.dimension == d)]
            floor = problem[(problem.n == 2) & (problem.method == "path")]
            floor_skill = float(floor.mean_skill.min()) if not floor.empty else math.nan
            if np.isfinite(floor_skill):
                axis.axhline(floor_skill, color="black", linestyle="--", linewidth=1.3, label="n=2 path floor")
            for method in methods:
                data = problem[(problem.n == 3) & (problem.method == method)].sort_values("M")
                if not data.empty:
                    axis.plot(data.M, data.mean_skill, marker="o", label=method, color=colors[method])
            axis.set_yscale("log")
            axis.set_title(f"{pde}, d={d}")
            axis.set_xlabel("M (n=3)")
            axis.set_ylabel("mean test skill")
            axis.grid(True, which="both", alpha=0.25)
            axis.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES_ROOT / "skill_n3.png", dpi=180)
    plt.close(fig)


def _plot_frontiers(frontier: pd.DataFrame) -> None:
    for axis, label in (
        ("generator_calls", "mean generator calls"),
        ("gaussian_coordinates", "mean Gaussian vectors x d"),
        ("wall_time", "mean wall time (seconds)"),
    ):
        fig, axes = plt.subplots(2, 2, figsize=(11, 8))
        for row_index, pde in enumerate(("P1", "MR")):
            for column, d in enumerate((100, 400)):
                panel = axes[row_index, column]
                for family, color in (("n2_path", "#4c78a8"), ("combined", "#e45756")):
                    data = frontier[
                        (frontier.pde_id == pde)
                        & (frontier.dimension == d)
                        & (frontier.cost_axis == axis)
                        & (frontier.family == family)
                        & frontier.is_frontier
                    ].sort_values("mean_cost")
                    if not data.empty:
                        panel.plot(data.mean_cost, data.mean_skill, marker="o", label=family, color=color)
                panel.set_xscale("log")
                panel.set_yscale("log")
                panel.set_title(f"{pde}, d={d}")
                panel.set_xlabel(label)
                panel.set_ylabel("mean test skill")
                panel.grid(True, which="both", alpha=0.25)
                panel.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(FIGURES_ROOT / f"frontier_{axis}.png", dpi=180)
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
    criteria: dict[str, Any],
    cell_detail: list[dict[str, Any]],
    levels: pd.DataFrame,
    predictions: dict[str, Any],
    fingerprint_audit: dict[str, Any],
    environment: dict[str, Any],
    prereg_commit: str,
    analysis_commit: str,
) -> None:
    lines = [
        "# Correction-level baseline at n=3: pre-registered report",
        "",
        "## Outcome table",
        "",
        "| Criterion | Verdict |",
        "|---|---|",
        f"| F-1 | {criteria['F-1']['verdict']} ({criteria['F-1']['wins']}/12) |",
        f"| F-2 | {criteria['F-2']['verdict']} ({criteria['F-2']['wins']}/12) |",
        f"| F-3 | {criteria['F-3']['verdict']} |",
        f"| F-4 | {criteria['F-4']['verdict']} (no verdict by design) |",
        "",
        "## Pre-registered criteria (verbatim)",
        "",
    ]
    for key in ("F-1", "F-2"):
        result = criteria[key]
        lines.extend(
            [
                f"### {key}",
                "",
                f"> {CRITERIA[key]}",
                "",
                f"**Verdict: {result['verdict']} ({result['wins']}/{result['total']}; required {result['required']}).**",
                "",
                "| PDE | d | M | ratio | criterion met |",
                "|---|---:|---:|---:|---|",
            ]
        )
        ratio_key = "path_base_over_path" if key == "F-1" else "double_path_base_over_double_path"
        met_key = "F1_met" if key == "F-1" else "F2_met"
        for item in cell_detail:
            lines.append(
                f"| {item['pde_id']} | {item['dimension']} | {item['M']} | {_number(item[ratio_key])} | {'yes' if item[met_key] else 'no'} |"
            )
        lines.append("")

    lines.extend(
        [
            "### F-3",
            "",
            f"> {CRITERIA['F-3']}",
            "",
            f"**Verdict: {criteria['F-3']['verdict']}.**",
            "",
            "| PDE | d | best n=2 path floor | best n=3 double family | ratio n3/n2 | verdict scope | met |",
            "|---|---:|---|---|---:|---|---|",
        ]
    )
    for item in criteria["F-3"]["problems"]:
        floor = f"M={item['n2_floor_M']}; {_number(item['n2_floor_skill'])}"
        deep = f"{item['n3_best_method']} M={item['n3_best_M']}; {_number(item['n3_best_skill'])}"
        lines.append(
            f"| {item['pde_id']} | {item['dimension']} | {floor} | {deep} | {_number(item['n3_over_n2_floor'])} | {'yes' if item['verdict_scope'] else 'reported only'} | {('yes' if item['criterion_met'] else 'no') if item['verdict_scope'] else 'NA'} |"
        )

    lines.extend(
        [
            "",
            "### F-4",
            "",
            f"> {CRITERIA['F-4']}",
            "",
            f"**Verdict: {criteria['F-4']['verdict']} (reported, no pass/fail criterion).**",
            "",
            "All 120 attaining-cell records and paired bootstrap intervals are in `results/baseline_n3/frontier_levels.csv`.",
            "",
            "| axis | PDE | d | levels | combined better levels |",
            "|---|---|---:|---:|---:|",
        ]
    )
    if not levels.empty:
        for (axis, pde, d), group in levels.groupby(["cost_axis", "pde_id", "dimension"]):
            lines.append(
                f"| {axis} | {pde} | {int(d)} | {len(group)} | {int((group.combined_over_n2 < 1.0).sum())} |"
            )

    lines.extend(["", "## Predictions (reported, no verdict)", ""])
    oracle = predictions["oracle"]
    oracle_ratios = [item["oracle_over_baseline"] for item in oracle]
    lines.append(
        f"- `path_base_oracle <= path_base`: {sum(item['prediction_met'] for item in oracle)}/{len(oracle)} cells; the oracle was instead {_number(np.min(oracle_ratios))}–{_number(np.max(oracle_ratios))}x the baseline skill (lower is better)."
    )
    variance = predictions["root_level1_variance"]
    if variance:
        ratios = [item["baseline_over_path"] for item in variance]
        lines.append(
            f"- Root level-1 integrand variance baseline/path: median {_number(np.median(ratios))}, range {_number(np.min(ratios))}–{_number(np.max(ratios))}; reused M=16 path rows lack this new diagnostic and were not rerun."
        )
    else:
        lines.append("- Root level-1 integrand variance comparison: not evaluated (no paired diagnostic rows).")
    scaling = predictions["path_scaling"]
    if scaling:
        values = [item["M_times_skill_over_d"] for item in scaling]
        lines.append(
            f"- For path, M x skill / d ranges from {_number(np.min(values))} to {_number(np.max(values))}; cell values are retained in `analysis_summary.json`."
        )

    failures: list[str] = []
    for key in ("F-1", "F-2", "F-3"):
        if criteria[key]["verdict"] == "FAIL":
            failures.append(key)
    lines.extend(
        [
            "",
            "## Failures and deviations",
            "",
            f"- Mechanical failures: {', '.join(failures) if failures else 'none'}.",
            "- No post-hoc thresholds were introduced.",
            "- Existing rows were reused without rerunning. Consequently, the newly requested root-integrand variance is unavailable for reused M=16 path/double-path and n=2 rows; this does not enter F-1 through F-4.",
            "",
            "## Provenance",
            "",
            f"- Frozen preregistration commit: `{prereg_commit}`.",
            f"- Analysis code commit: `{analysis_commit}`.",
            f"- Solver/run code commit: `{environment.get('code_commit', 'NA')}`.",
            "- Parent path-frontier commit: `a42e2c5cf3b3ad78ba702d42f3de6ab559150f0d`.",
            f"- Environment: {environment.get('platform', 'NA')}; Python {environment.get('python', 'NA')}; NumPy {environment.get('numpy', 'NA')}; {environment.get('logical_cpus', 'NA')} logical CPUs, {environment.get('workers', 'NA')} workers.",
            f"- Rows: {len(rows)} / {EXPECTED_ROWS}; repetitions: {REPETITIONS} per cell.",
            f"- Reused rows: {int((rows.source_result != 'baseline_n3').sum()) if not rows.empty else 0}; new rows: {int((rows.source_result == 'baseline_n3').sum()) if not rows.empty else 0}.",
            f"- Primary-draw fingerprint audit: {'PASS' if fingerprint_audit.get('passed') else 'FAIL'}; {int(fingerprint_audit.get('pairs', 0))} paired chunks, {int(fingerprint_audit.get('mismatches', 0))} mismatches; records SHA-256 `{fingerprint_audit.get('records_sha256', 'NA')}`.",
            f"- Non-finite state values: {int(rows.nonfinite_state_count.sum()) if not rows.empty else 0}; non-finite generator values: {int(rows.nonfinite_generator_count.sum()) if not rows.empty else 0}.",
            "- Test points: frozen 1,200 points with the frozen 20% validation split; all criteria use the 960 test points.",
            "- Arithmetic: float64; bootstrap: 1,000 paired repetition draws with seed 20261201.",
            "",
            "## Figures",
            "",
            "- `results/baseline_n3/figures/skill_n3.png`",
            "- `results/baseline_n3/figures/frontier_generator_calls.png`",
            "- `results/baseline_n3/figures/frontier_gaussian_coordinates.png`",
            "- `results/baseline_n3/figures/frontier_wall_time.png`",
            "",
        ]
    )
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text("\n".join(lines), encoding="utf-8")


def analyze() -> dict[str, Any]:
    rows = _load_rows()
    summary = _summary(rows)
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    summary.to_csv(RESULTS_ROOT / "summary.csv", index=False)
    criteria, cell_detail = _criteria(summary)
    frontier, levels = _frontiers(rows, summary) if not rows.empty else (pd.DataFrame(), pd.DataFrame())
    frontier.to_csv(RESULTS_ROOT / "frontier.csv", index=False)
    levels.to_csv(RESULTS_ROOT / "frontier_levels.csv", index=False)
    predictions = _predictions(summary)
    prereg_commit = _git(
        "log",
        "--diff-filter=A",
        "-1",
        "--format=%H",
        "--",
        "invariant_region_mlp/experiments/baseline_n3/PREREGISTRATION.md",
    )
    analysis_commit = _git("rev-parse", "HEAD")
    audit_path = RESULTS_ROOT / "fingerprint_audit.json"
    audit = json.loads(audit_path.read_text(encoding="utf-8")) if audit_path.is_file() else {"passed": False, "pairs": 0}
    environment_path = RESULTS_ROOT / "environment.json"
    environment = json.loads(environment_path.read_text(encoding="utf-8")) if environment_path.is_file() else {}
    payload = {
        "criteria": criteria,
        "cell_detail": cell_detail,
        "predictions": predictions,
        "completed_rows": len(rows),
        "registered_rows": EXPECTED_ROWS,
        "frontier_level_rows": len(levels),
        "bootstrap_draws": BOOTSTRAP_DRAWS,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "preregistration_commit": prereg_commit,
        "analysis_code_commit": analysis_commit,
        "fingerprint_audit": audit,
        "environment": environment,
    }
    (RESULTS_ROOT / "analysis_summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=True) + "\n",
        encoding="utf-8",
    )
    if not summary.empty:
        _plot_skill(summary)
        _plot_frontiers(frontier)
    _report(rows, criteria, cell_detail, levels, predictions, audit, environment, prereg_commit, analysis_commit)
    print(json.dumps(payload, indent=2, sort_keys=True, allow_nan=True))
    return payload


def main() -> int:
    analyze()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
