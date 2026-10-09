"""Mechanical frontier analysis, bootstrap, figure, and report."""

from __future__ import annotations

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

from .protocol import (
    DOUBLE_RESULTS_ROOT,
    EC_CELLS,
    EXPECTED_COMBINED_ROWS,
    EXPECTED_PATH_ROWS,
    METHODS,
    PROBLEMS,
    REPETITIONS,
    RESULTS_ROOT,
)


REPORT_PATH = RESULTS_ROOT.parents[1] / "docs" / "PATH_FRONTIER_REPORT.md"
FIGURE_PATH = RESULTS_ROOT / "figures" / "frontier.png"
BOOTSTRAP_DRAWS = 1000
BOOTSTRAP_SEED = 20261201

CRITERIA = {
    "E-1": (
        "Let `best_double` be the lower envelope over `double` and "
        "`double_path`. \"Double adds value at the frontier\" if `best_double "
        "<= 0.9 x path` at >= 8 of the 10 cost levels, for at least 3 of the "
        "4 (PDE, d) problems. Otherwise report \"pathwise alone suffices at "
        "the frontier\"."
    ),
    "E-2": (
        "the `path` frontier is <= 0.8 x the `raw` frontier at >= 8 of 10 "
        "cost levels, for each of the 4 problems."
    ),
}


def _repository_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _git(*args: str) -> str:
    return subprocess.run(
        ["git", *args],
        cwd=_repository_root(),
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def _load_rows() -> pd.DataFrame:
    path = RESULTS_ROOT / "rows.csv"
    if not path.is_file():
        return pd.DataFrame()
    frame = pd.read_csv(path)
    if frame.empty:
        return frame
    numeric = (
        "dimension",
        "n",
        "M",
        "repetition",
        "skill",
        "generator_calls",
        "wall_time_seconds",
        "nonfinite_state_count",
        "nonfinite_generator_count",
    )
    for column in numeric:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return frame.sort_values(
        ["pde_id", "dimension", "n", "M", "method", "repetition"]
    ).reset_index(drop=True)


def _summarize(rows: pd.DataFrame) -> pd.DataFrame:
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
        "mean_generator_calls",
        "mean_wall_time_seconds",
        "total_nonfinite_states",
        "total_nonfinite_generators",
    ]
    if rows.empty:
        return pd.DataFrame(columns=columns)
    result = (
        rows.groupby(["pde_id", "dimension", "n", "M", "method"], as_index=False)
        .agg(
            repetitions=("repetition", "nunique"),
            mean_skill=("skill", "mean"),
            median_skill=("skill", "median"),
            skill_std=("skill", "std"),
            mean_generator_calls=("generator_calls", "mean"),
            mean_wall_time_seconds=("wall_time_seconds", "mean"),
            total_nonfinite_states=("nonfinite_state_count", "sum"),
            total_nonfinite_generators=("nonfinite_generator_count", "sum"),
        )
    )
    return result[columns]


def _candidate(points: pd.DataFrame, budget: float) -> pd.Series | None:
    eligible = points[
        (points.mean_generator_calls <= budget * (1.0 + 1e-12))
        & np.isfinite(points.mean_skill)
    ]
    if eligible.empty:
        return None
    ordered = eligible.sort_values(
        ["mean_skill", "mean_generator_calls", "n", "M", "method"],
        kind="stable",
    )
    return ordered.iloc[0]


def _selection(point: pd.Series | None, prefix: str) -> dict[str, Any]:
    if point is None:
        return {
            f"{prefix}_skill": math.nan,
            f"{prefix}_method": None,
            f"{prefix}_n": None,
            f"{prefix}_M": None,
            f"{prefix}_generator_calls": math.nan,
            f"{prefix}_wall_time_seconds": math.nan,
        }
    return {
        f"{prefix}_skill": float(point.mean_skill),
        f"{prefix}_method": str(point.method),
        f"{prefix}_n": int(point.n),
        f"{prefix}_M": int(point.M),
        f"{prefix}_generator_calls": float(point.mean_generator_calls),
        f"{prefix}_wall_time_seconds": float(point.mean_wall_time_seconds),
    }


def _frontier_points(summary: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for (pde_id, d, method), points in summary.groupby(
        ["pde_id", "dimension", "method"], sort=True
    ):
        points = points.copy()
        for _, point in points.iterrows():
            dominated = points[
                (points.mean_generator_calls <= point.mean_generator_calls)
                & (points.mean_skill <= point.mean_skill)
                & (
                    (points.mean_generator_calls < point.mean_generator_calls)
                    | (points.mean_skill < point.mean_skill)
                )
            ]
            records.append(
                {
                    "pde_id": pde_id,
                    "dimension": int(d),
                    "method": method,
                    "n": int(point.n),
                    "M": int(point.M),
                    "mean_generator_calls": float(point.mean_generator_calls),
                    "mean_wall_time_seconds": float(point.mean_wall_time_seconds),
                    "mean_skill": float(point.mean_skill),
                    "skill_std": float(point.skill_std),
                    "is_frontier": dominated.empty,
                }
            )
    return pd.DataFrame(records).sort_values(
        ["pde_id", "dimension", "method", "mean_generator_calls", "mean_skill"]
    )


def _registered_cost_levels() -> dict[tuple[str, int], np.ndarray]:
    source = pd.read_csv(DOUBLE_RESULTS_ROOT / "frontier_levels.csv")
    output: dict[tuple[str, int], np.ndarray] = {}
    for pde_id, d in PROBLEMS:
        subset = source[
            (source.pde_id == pde_id) & (source.dimension == d)
        ].sort_values("level_index")
        levels = subset.cost_level.to_numpy(dtype=np.float64)
        if len(levels) != 10 or len(np.unique(levels)) != 10:
            raise ValueError(f"D-4 cost levels are incomplete for {pde_id}, d={d}")
        output[(pde_id, d)] = levels
    return output


def _complete_problem(summary: pd.DataFrame, pde_id: str, d: int) -> bool:
    points = summary[(summary.pde_id == pde_id) & (summary.dimension == d)]
    return len(points) == len(EC_CELLS) * len(METHODS) and bool(
        np.all(points.repetitions.to_numpy() == REPETITIONS)
    )


def _level_rows(
    rows: pd.DataFrame,
    summary: pd.DataFrame,
    cost_levels: dict[tuple[str, int], np.ndarray],
) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    for pde_id, d in PROBLEMS:
        points = summary[(summary.pde_id == pde_id) & (summary.dimension == d)]
        complete = _complete_problem(summary, pde_id, d)
        levels = cost_levels[(pde_id, d)]
        problem_records: list[dict[str, Any]] = []
        for index, budget in enumerate(levels):
            record: dict[str, Any] = {
                "pde_id": pde_id,
                "dimension": d,
                "level_index": index,
                "cost_level": float(budget),
                "complete": complete,
            }
            for method in METHODS:
                point = _candidate(points[points.method == method], float(budget))
                record.update(_selection(point, method))
            best = _candidate(
                points[points.method.isin(["double", "double_path"])],
                float(budget),
            )
            record.update(_selection(best, "best_double"))
            path_skill = float(record["path_skill"])
            raw_skill = float(record["raw_skill"])
            best_skill = float(record["best_double_skill"])
            double_ratio = best_skill / path_skill if path_skill > 0.0 else math.inf
            path_ratio = path_skill / raw_skill if raw_skill > 0.0 else math.inf
            record.update(
                {
                    "best_double_over_path": double_ratio,
                    "double_adds": bool(np.isfinite(double_ratio) and double_ratio <= 0.9),
                    "path_over_raw": path_ratio,
                    "path_beats_raw": bool(np.isfinite(path_ratio) and path_ratio <= 0.8),
                }
            )
            problem_records.append(record)

        boot_double = np.full((BOOTSTRAP_DRAWS, 10), np.nan, dtype=np.float64)
        boot_path = np.full((BOOTSTRAP_DRAWS, 10), np.nan, dtype=np.float64)
        problem_rows = rows[(rows.pde_id == pde_id) & (rows.dimension == d)]
        if complete:
            for draw in range(BOOTSTRAP_DRAWS):
                sampled_reps = rng.integers(0, REPETITIONS, REPETITIONS)
                sampled = pd.concat(
                    [problem_rows[problem_rows.repetition == rep] for rep in sampled_reps],
                    ignore_index=True,
                )
                boot = (
                    sampled.groupby(["n", "M", "method"], as_index=False)
                    .agg(
                        mean_skill=("skill", "mean"),
                        mean_generator_calls=("generator_calls", "mean"),
                        mean_wall_time_seconds=("wall_time_seconds", "mean"),
                    )
                )
                for index, budget in enumerate(levels):
                    raw = _candidate(boot[boot.method == "raw"], float(budget))
                    path = _candidate(boot[boot.method == "path"], float(budget))
                    best = _candidate(
                        boot[boot.method.isin(["double", "double_path"])],
                        float(budget),
                    )
                    if raw is not None and path is not None and float(raw.mean_skill) > 0:
                        boot_path[draw, index] = float(path.mean_skill / raw.mean_skill)
                    if path is not None and best is not None and float(path.mean_skill) > 0:
                        boot_double[draw, index] = float(best.mean_skill / path.mean_skill)

        for record in problem_records:
            index = int(record["level_index"])
            for values, prefix, threshold in (
                (boot_double[:, index], "best_double_over_path", 0.9),
                (boot_path[:, index], "path_over_raw", 0.8),
            ):
                finite = values[np.isfinite(values)]
                record[f"{prefix}_bootstrap_low"] = (
                    float(np.quantile(finite, 0.025)) if len(finite) else math.nan
                )
                record[f"{prefix}_bootstrap_high"] = (
                    float(np.quantile(finite, 0.975)) if len(finite) else math.nan
                )
                record[f"{prefix}_bootstrap_win_probability"] = (
                    float(np.mean(finite <= threshold)) if len(finite) else math.nan
                )
            records.append(record)
    return pd.DataFrame(records)


def _criteria(levels: pd.DataFrame) -> dict[str, Any]:
    e1_problems: dict[str, Any] = {}
    e2_problems: dict[str, Any] = {}
    all_complete = True
    for pde_id, d in PROBLEMS:
        subset = levels[(levels.pde_id == pde_id) & (levels.dimension == d)]
        key = f"{pde_id}_d{d}"
        complete = len(subset) == 10 and bool(subset.complete.all())
        all_complete = all_complete and complete
        e1_wins = int(subset.double_adds.sum()) if complete else 0
        e2_wins = int(subset.path_beats_raw.sum()) if complete else 0
        e1_problems[key] = {
            "complete": complete,
            "wins": e1_wins,
            "required_wins": 8,
            "criterion_met": bool(complete and e1_wins >= 8),
        }
        e2_problems[key] = {
            "complete": complete,
            "wins": e2_wins,
            "required_wins": 8,
            "criterion_met": bool(complete and e2_wins >= 8),
        }

    e1_problem_wins = int(sum(item["criterion_met"] for item in e1_problems.values()))
    if not all_complete:
        e1_verdict = "NOT EVALUATED (COMPUTE)"
        e2_verdict = "NOT EVALUATED (COMPUTE)"
    else:
        e1_verdict = (
            "Double adds value at the frontier"
            if e1_problem_wins >= 3
            else "pathwise alone suffices at the frontier"
        )
        e2_verdict = (
            "PASS"
            if all(item["criterion_met"] for item in e2_problems.values())
            else "FAIL"
        )
    return {
        "E-1": {
            "verdict": e1_verdict,
            "problems_meeting_criterion": e1_problem_wins,
            "required_problems": 3,
            "problems": e1_problems,
        },
        "E-2": {"verdict": e2_verdict, "problems": e2_problems},
    }


def _dimension_robustness(
    summary: pd.DataFrame, cost_levels: dict[tuple[str, int], np.ndarray]
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for pde_id in ("P1", "MR"):
        levels_100 = cost_levels[(pde_id, 100)]
        levels_400 = cost_levels[(pde_id, 400)]
        if not np.allclose(levels_100, levels_400, rtol=0.0, atol=0.0):
            raise ValueError(f"dimension cost levels differ for {pde_id}")
        # Conventional numeric median of the ten frozen D-4 cost levels.
        budget = float(np.median(levels_100))
        for label, methods in (
            ("path", ["path"]),
            ("best_double", ["double", "double_path"]),
            ("raw", ["raw"]),
        ):
            selections: dict[int, pd.Series | None] = {}
            for d in (100, 400):
                points = summary[
                    (summary.pde_id == pde_id)
                    & (summary.dimension == d)
                    & summary.method.isin(methods)
                ]
                selections[d] = _candidate(points, budget)
            low = selections[100]
            high = selections[400]
            low_skill = float(low.mean_skill) if low is not None else math.nan
            high_skill = float(high.mean_skill) if high is not None else math.nan
            ratio = high_skill / low_skill if low_skill > 0.0 else math.inf
            target = ratio >= 2.0 if label == "raw" else ratio <= 1.3
            records.append(
                {
                    "pde_id": pde_id,
                    "method": label,
                    "median_cost_level": budget,
                    "d100_skill": low_skill,
                    "d100_method": str(low.method) if low is not None else None,
                    "d100_cell": [int(low.n), int(low.M)] if low is not None else None,
                    "d400_skill": high_skill,
                    "d400_method": str(high.method) if high is not None else None,
                    "d400_cell": [int(high.n), int(high.M)] if high is not None else None,
                    "ratio_d400_over_d100": ratio,
                    "prediction": ">= 2" if label == "raw" else "<= 1.3",
                    "prediction_met": bool(np.isfinite(ratio) and target),
                }
            )
    return records


def _depth_comparison(summary: pd.DataFrame) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for pde_id, d in PROBLEMS:
        problem = summary[(summary.pde_id == pde_id) & (summary.dimension == d)]
        for label, methods in (
            ("path", ["path"]),
            ("best_double", ["double", "double_path"]),
        ):
            selected: dict[int, pd.Series | None] = {}
            for n in (2, 3):
                candidates = problem[(problem.n == n) & problem.method.isin(methods)]
                selected[n] = _candidate(candidates, float("inf"))
            shallow = selected[2]
            deep = selected[3]
            shallow_skill = float(shallow.mean_skill) if shallow is not None else math.nan
            deep_skill = float(deep.mean_skill) if deep is not None else math.nan
            records.append(
                {
                    "pde_id": pde_id,
                    "dimension": d,
                    "method": label,
                    "n2_skill": shallow_skill,
                    "n2_method": str(shallow.method) if shallow is not None else None,
                    "n2_M": int(shallow.M) if shallow is not None else None,
                    "n2_generator_calls": (
                        float(shallow.mean_generator_calls) if shallow is not None else math.nan
                    ),
                    "n3_skill": deep_skill,
                    "n3_method": str(deep.method) if deep is not None else None,
                    "n3_M": int(deep.M) if deep is not None else None,
                    "n3_generator_calls": (
                        float(deep.mean_generator_calls) if deep is not None else math.nan
                    ),
                    "n3_over_n2": (
                        deep_skill / shallow_skill if shallow_skill > 0.0 else math.inf
                    ),
                }
            )
    return records


def _plot(frontier: pd.DataFrame) -> None:
    FIGURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    colors = {
        "raw": "#444444",
        "path": "#4c78a8",
        "double": "#f58518",
        "double_path": "#e45756",
    }
    for row_index, pde_id in enumerate(("P1", "MR")):
        for column, d in enumerate((100, 400)):
            axis = axes[row_index, column]
            for method in METHODS:
                data = frontier[
                    (frontier.pde_id == pde_id)
                    & (frontier.dimension == d)
                    & (frontier.method == method)
                    & frontier.is_frontier
                ].sort_values("mean_generator_calls")
                if data.empty:
                    continue
                axis.plot(
                    data.mean_generator_calls,
                    data.mean_skill,
                    marker="o",
                    linewidth=1.7,
                    markersize=5,
                    label=method,
                    color=colors[method],
                )
            axis.set_xscale("log")
            axis.set_yscale("log")
            axis.set_title(f"{pde_id}, d={d}")
            axis.set_xlabel("mean generator calls")
            axis.set_ylabel("mean test skill")
            axis.grid(True, which="both", alpha=0.25)
            axis.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURE_PATH, dpi=180)
    plt.close(fig)


def _number(value: Any) -> str:
    if value is None:
        return "NA"
    if isinstance(value, (float, np.floating)):
        return "NA" if not np.isfinite(value) else f"{float(value):.6g}"
    return str(value)


def _cell(row: pd.Series, prefix: str) -> str:
    if pd.isna(row[f"{prefix}_n"]):
        return "NA"
    method = str(row[f"{prefix}_method"])
    return (
        f"{method} ({int(row[f'{prefix}_n'])},{int(row[f'{prefix}_M'])}); "
        f"{_number(row[f'{prefix}_skill'])}"
    )


def _report(
    rows: pd.DataFrame,
    levels: pd.DataFrame,
    criteria: dict[str, Any],
    e3: list[dict[str, Any]],
    depth: list[dict[str, Any]],
    prereg_commit: str,
    analysis_commit: str,
    audit: dict[str, Any],
) -> None:
    lines = [
        "# Pathwise-only equal-cost frontier: pre-registered report",
        "",
        "## Outcome table",
        "",
        "| Criterion | Verdict |",
        "|---|---|",
        f"| E-1 | {criteria['E-1']['verdict']} |",
        f"| E-2 | {criteria['E-2']['verdict']} |",
        "",
        "## Pre-registered criteria (verbatim)",
        "",
    ]
    for criterion in ("E-1", "E-2"):
        result = criteria[criterion]
        lines.extend(
            [
                f"### {criterion}",
                "",
                f"> {CRITERIA[criterion]}",
                "",
                f"**Verdict: {result['verdict']}**",
                "",
                "| Problem | qualifying levels | required | problem criterion |",
                "|---|---:|---:|---|",
            ]
        )
        for key, detail in result["problems"].items():
            lines.append(
                f"| {key} | {detail['wins']} | {detail['required_wins']} | "
                f"{'yes' if detail['criterion_met'] else 'no'} |"
            )
        lines.append("")

    lines.extend(
        [
            "## E-3: dimension robustness (reported, no verdict)",
            "",
            "The median budget is the conventional numeric median of the ten frozen D-4 cost levels. The prediction column is reported mechanically but is not a criterion.",
            "",
            "| PDE | method | median cost | d=100 skill (cell) | d=400 skill (cell) | ratio 400/100 | prediction | met |",
            "|---|---|---:|---|---|---:|---|---|",
        ]
    )
    for item in e3:
        lines.append(
            "| {pde} | {method} | {cost} | {s100} ({c100}) | {s400} ({c400}) | {ratio} | {prediction} | {met} |".format(
                pde=item["pde_id"],
                method=item["method"],
                cost=_number(item["median_cost_level"]),
                s100=_number(item["d100_skill"]),
                c100=item["d100_cell"],
                s400=_number(item["d400_skill"]),
                c400=item["d400_cell"],
                ratio=_number(item["ratio_d400_over_d100"]),
                prediction=item["prediction"],
                met="yes" if item["prediction_met"] else "no",
            )
        )

    lines.extend(
        [
            "",
            "## Frontier-attaining cells (reported, no verdict)",
            "",
            "Each entry is `method (n,M); mean skill`. Wall-time selections and bootstrap intervals are retained in `frontier_levels.csv`.",
            "",
            "| PDE | d | level | calls | raw | path | double | double_path | best_double |",
            "|---|---:|---:|---:|---|---|---|---|---|",
        ]
    )
    for _, row in levels.sort_values(["pde_id", "dimension", "level_index"]).iterrows():
        lines.append(
            f"| {row.pde_id} | {int(row.dimension)} | {int(row.level_index)} | "
            f"{_number(row.cost_level)} | {_cell(row, 'raw')} | {_cell(row, 'path')} | "
            f"{_cell(row, 'double')} | {_cell(row, 'double_path')} | "
            f"{_cell(row, 'best_double')} |"
        )

    lines.extend(
        [
            "",
            "## Best n=3 versus n=2 (reported, no verdict)",
            "",
            "| PDE | d | method family | best n=2 | best n=3 | n3/n2 |",
            "|---|---:|---|---|---|---:|",
        ]
    )
    for item in depth:
        n2 = f"{item['n2_method']} M={item['n2_M']}; {_number(item['n2_skill'])}"
        n3 = f"{item['n3_method']} M={item['n3_M']}; {_number(item['n3_skill'])}"
        lines.append(
            f"| {item['pde_id']} | {item['dimension']} | {item['method']} | "
            f"{n2} | {n3} | {_number(item['n3_over_n2'])} |"
        )

    path_rows = rows[rows.method == "path"] if not rows.empty else rows
    new_rows = path_rows[path_rows.source_result == "path_frontier"] if not rows.empty else rows
    reused_rows = path_rows[path_rows.source_result == "double_estimator"] if not rows.empty else rows
    code_commits = (
        sorted(new_rows.code_commit.dropna().unique().tolist()) if not new_rows.empty else []
    )
    missing = EXPECTED_COMBINED_ROWS - len(rows)
    lines.extend(
        [
            "",
            "## Provenance",
            "",
            f"- Frozen preregistration commit: `{prereg_commit}`.",
            f"- Analysis code commit: `{analysis_commit}`.",
            "- Frozen double-estimator source/results commit: `50c3c7d8c83db707cd6c6878835132547a7d2fab`.",
            f"- New-row code commits: {', '.join(f'`{value}`' for value in code_commits) if code_commits else 'none' }.",
            f"- Combined rows: {len(rows)} / {EXPECTED_COMBINED_ROWS}; missing: {missing}.",
            f"- Path rows: {len(path_rows)} / {EXPECTED_PATH_ROWS}; reused: {len(reused_rows)}; newly computed: {len(new_rows)}.",
            f"- New-row aggregate worker wall time: {_number(new_rows.wall_time_seconds.sum() if not new_rows.empty else 0.0)} seconds.",
            f"- Non-finite state values: {int(rows.nonfinite_state_count.sum()) if not rows.empty else 0}; non-finite generator values: {int(rows.nonfinite_generator_count.sum()) if not rows.empty else 0}.",
            f"- Paired primary-draw audit: {audit.get('pairs', 0)} pairs, {audit.get('fingerprint_mismatches', 'NA')} fingerprint mismatches, verdict {'PASS' if audit.get('passed') else 'FAIL/NOT RUN'}.",
            "- The paired audit traces one complete registered chunk for every task identity. It does not rerun or overwrite any existing 1,200-point result row.",
            "- Test points: the frozen 1,200 points per PDE/dimension and frozen 20% validation split; verdicts use the 960-point test subset.",
            "- Arithmetic: float64 throughout; 10 repetitions per cell; bootstrap: 1,000 paired repetition draws with seed 20261201.",
            "- Primary cost: mean generator calls. Mean wall time is carried as the secondary cost in both frontier output tables.",
            "",
            "## Figure",
            "",
            "- `results/path_frontier/figures/frontier.png`",
            "",
        ]
    )
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text("\n".join(lines), encoding="utf-8")


def analyze() -> dict[str, Any]:
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    rows = _load_rows()
    summary = _summarize(rows)
    summary.to_csv(RESULTS_ROOT / "summary.csv", index=False)
    frontier = _frontier_points(summary)
    frontier.to_csv(RESULTS_ROOT / "frontier.csv", index=False)
    cost_levels = _registered_cost_levels()
    levels = _level_rows(rows, summary, cost_levels)
    levels.to_csv(RESULTS_ROOT / "frontier_levels.csv", index=False)
    criteria = _criteria(levels)
    e3 = _dimension_robustness(summary, cost_levels)
    depth = _depth_comparison(summary)
    audit_path = RESULTS_ROOT / "fingerprint_audit.json"
    audit = (
        json.loads(audit_path.read_text(encoding="utf-8"))
        if audit_path.is_file()
        else {"passed": False, "pairs": 0, "fingerprint_mismatches": None}
    )
    prereg_commit = _git(
        "log",
        "--diff-filter=A",
        "-1",
        "--format=%H",
        "--",
        "invariant_region_mlp/experiments/path_frontier/PREREGISTRATION.md",
    )
    analysis_commit = _git("rev-parse", "HEAD")
    payload = {
        "criteria": criteria,
        "dimension_robustness": e3,
        "depth_comparison": depth,
        "completed_rows": len(rows),
        "registered_rows": EXPECTED_COMBINED_ROWS,
        "bootstrap_draws": BOOTSTRAP_DRAWS,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "preregistration_commit": prereg_commit,
        "analysis_code_commit": analysis_commit,
        "fingerprint_audit": audit,
    }
    (RESULTS_ROOT / "analysis_summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=True) + "\n",
        encoding="utf-8",
    )
    _plot(frontier)
    _report(rows, levels, criteria, e3, depth, prereg_commit, analysis_commit, audit)
    print(json.dumps(payload, indent=2, sort_keys=True, allow_nan=True))
    return payload


def main() -> int:
    analyze()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
