"""Mechanical analysis, figures, and report for mechanism-suite round 3."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
from typing import Any, Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy

from invariant_region_mlp.experiments.mechanism_suite.equations import stable_seed

from .equations import BASE_SEED, RESULTS_ROOT
from .mlp import IMPLEMENTATION_REVISION, PREREGISTRATION_COMMIT, load_repetition
from .protocol import (
    D400_PILOT_PATH,
    GATES_PATH,
    S1_PDES,
    primary_method,
    s1_tasks,
    s2_tasks,
    s3_tasks,
    s4_tasks,
    task_identity,
    task_path,
)


HERE = Path(__file__).resolve().parent
PACKAGE_ROOT = HERE.parents[1]
REPOSITORY_ROOT = PACKAGE_ROOT.parent
REPORT_PATH = PACKAGE_ROOT / "docs" / "MECHANISM_SUITE_R3_REPORT.md"
FIGURE_ROOT = RESULTS_ROOT / "figures"
S1_PATH = RESULTS_ROOT / "s1_frontier.csv"
LEVELS_PATH = RESULTS_ROOT / "s1_frontier_levels.csv"
S2_PATH = RESULTS_ROOT / "s2_theory.csv"
S3_PATH = RESULTS_ROOT / "s3_dimension.csv"
S4_PATH = RESULTS_ROOT / "s4_geometry.csv"
TUNING_PATH = RESULTS_ROOT / "tuning_choices.json"
SUMMARY_PATH = RESULTS_ROOT / "analysis_summary.json"
ONESTEP_PATH = RESULTS_ROOT / "s2_onestep_check.json"
FULL_HISTORY_SOURCE = HERE.parent / "active_vb_high_budget" / "vb_mlp_methods.py"

CRITERIA = {
    "PF-1": "**PF-1 (frontier):** for each PDE in S1, at d=100 and at d=400, the primary method's frontier has skill <= 0.8 x the raw frontier at >= 8 of the 10 cost levels.",
    "PF-2": "**PF-2 (dimension robustness at equal cost):** for P1, C1, C2-convex, C3: at the median cost level, raw-frontier skill at d=400 is >= 2 x its value at d=20, while the primary-method frontier at d=400 is <= 1.3 x its value at d=20.",
    "PF-3": "**PF-3 (data matters):** in the S1 cells, the primary method wins >= 7/10 paired reps against `centre` in >= 75% of cells, for P1, C1, C2-convex, C3. Prediction for P4 (tight certificate): fails PF-3 (round-2 replication).",
    "PF-4": "**PF-4 (beyond generic suppression):** same as PF-3 against the best tuned `shrink_c`.",
    "T-1": "**T-1 (sign):** in S2 at cells (2,8) and (2,32), mean raw generator bias is < 0 for C2-convex and > 0 for C2-flip at every d, and |bias(C2-cancel)| <= 0.25 x min(|bias(convex)|, |bias(flip)|) at the same d and cell.",
    "T-2": "**T-2 (where the curse appears):** best raw skill over the S2 cells grows by >= 2x from d=20 to d=400 for C2-convex and C2-flip, and by <= 1.3x for C2-cancel.",
    "T-3": "**T-3 (one-step quantitative):** in the extended terminal-block check, measured Jensen gap / predicted `(c/M)K` lies in [0.9, 1.1] for every C2 configuration with |K| >= 5.",
    "DIM": "**DIM (S3):** for P1, C1, C2-convex, C3, primary-method skill at d=400 / d=20 <= 1.3 at both cells.",
    "LQG-G": "**LQG-G (exploratory, no claim):** report the gate outcome against the stated prediction.",
}


def _git(*args: str) -> str | None:
    try:
        return subprocess.check_output(
            ["git", *args], cwd=REPOSITORY_ROOT, text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, (np.floating, float)):
        number = float(value)
        if math.isnan(number):
            return None
        if math.isinf(number):
            return "Infinity" if number > 0 else "-Infinity"
        return number
    return value


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(_json_safe(payload), indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _write_csv(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    frame.to_csv(temporary, index=False)
    temporary.replace(path)


def _admission() -> tuple[dict[str, bool], dict[str, Any]]:
    gates = json.loads(GATES_PATH.read_text(encoding="utf-8"))
    admission = {
        pde: bool(summary["admitted"])
        for pde, summary in gates["candidates"].items()
    }
    return admission, gates


def _dropped_configs() -> set[tuple[int, int]]:
    if not D400_PILOT_PATH.exists():
        return set()
    payload = json.loads(D400_PILOT_PATH.read_text(encoding="utf-8"))
    return {tuple(map(int, item)) for item in payload.get("dropped_configs", [])}


def _study_tasks(study: str, admission: dict[str, bool]) -> list[dict[str, Any]]:
    if study == "s1":
        tasks = s1_tasks(include_lqg=admission.get("LQG", False))
        dropped = _dropped_configs()
        tasks = [
            task for task in tasks
            if not (
                int(task["d"]) == 400
                and (int(task["n"]), int(task["M"])) in dropped
            )
        ]
    elif study == "s2":
        tasks = s2_tasks()
    elif study == "s3":
        tasks = s3_tasks()
    elif study == "s4":
        tasks = s4_tasks()
    else:
        raise ValueError(study)
    unique: dict[tuple[Any, ...], dict[str, Any]] = {}
    for task in tasks:
        if admission.get(str(task["pde_id"]), False):
            unique[task_identity(task)] = task
    return list(unique.values())


def _flatten(task: dict[str, Any], metadata: dict[str, Any]) -> dict[str, Any]:
    work = metadata["work"]
    generator = work.get("generator") or {}
    test = metadata["metrics"]["test"]
    validation = metadata["metrics"]["validation"]
    extra = metadata.get("extra_diagnostics") or {}
    return {
        "study": task["study"],
        "artifact": task_path(task).relative_to(RESULTS_ROOT).as_posix(),
        "pde": task["pde_id"],
        "d": int(task["d"]),
        "n": int(task["n"]),
        "M": int(task["M"]),
        "method": task["method"]["name"],
        "transform": task["method"]["transform"],
        "factor": float(task["method"]["factor"]),
        "rep": int(task["repetition"]),
        "chunk_size": int(task["chunk_size"]),
        "test_skill": float(test["skill"]),
        "test_value_rmse": float(test["value_rmse"]),
        "test_value_bias": float(test["value_bias"]),
        "test_value_mae": float(test["value_mae"]),
        "test_gradient_relative_l2": float(test["gradient_relative_l2"]),
        "validation_skill": float(validation["skill"]),
        "generator_calls": int(work.get("f_evals", 0)),
        "recursively_evaluated_states": int(work.get("recursively_evaluated_states", 0)),
        "stochastic_samples": int(work.get("total_stochastic_samples", 0)),
        "terminal_samples": int(work.get("terminal_samples", 0)),
        "transition_samples": int(work.get("transition_samples", 0)),
        "wall_time_seconds": float(metadata["wall_clock_seconds"]),
        "nonfinite_state_count": int(test.get("nonfinite_state_count", 0)),
        "nonfinite_generator_count": int(work.get("nonfinite_generator_values", 0)),
        "nonfinite_count": int(test.get("nonfinite_state_count", 0))
        + int(work.get("nonfinite_generator_values", 0)),
        "generator_bias": generator.get("bias"),
        "generator_mae": generator.get("mae"),
        "generator_mse": generator.get("mse"),
        "projection_activation_rate": work.get("projection_activation_rate"),
        "pre_box_violation_rate": work.get("pre_box_violation_rate"),
        "child_noise_count": (extra.get("child_u_z_noise") or {}).get("count"),
        "code_commit": metadata.get("code_commit_at_execution"),
        "exact_reuse": bool(metadata.get("exact_reuse")),
    }


def load_study(
    study: str, admission: dict[str, bool]
) -> tuple[pd.DataFrame, list[str]]:
    rows = []
    missing = []
    for task in _study_tasks(study, admission):
        path = task_path(task)
        if not path.exists():
            missing.append(path.relative_to(RESULTS_ROOT).as_posix())
            continue
        loaded = load_repetition(path)
        rows.append(_flatten(task, loaded["metadata"]))
    return pd.DataFrame(rows), missing


def tune_shrink(frames: Iterable[pd.DataFrame]) -> tuple[dict[str, Any], pd.DataFrame]:
    combined = pd.concat([frame for frame in frames if len(frame)], ignore_index=True)
    shrink = combined[combined["method"].str.startswith("shrink_c")].copy()
    choices = []
    virtual_rows = []
    keys = ["pde", "d", "n", "M"]
    for key, group in shrink.groupby(keys, sort=True):
        scores = group.groupby("method")["validation_skill"].mean().sort_values()
        selected = str(scores.index[0])
        choices.append({
            **dict(zip(keys, key)),
            "selected_method": selected,
            "mean_validation_skill_by_method": {
                str(name): float(value) for name, value in scores.items()
            },
        })
        chosen = group[group["method"] == selected].copy()
        chosen["source_method"] = selected
        chosen["method"] = "best_shrink"
        virtual_rows.append(chosen)
    payload = {
        "selection_rule": "minimum mean validation skill per PDE/d/n/M cell; ties lexicographic",
        "validation_fraction": 0.2,
        "choices": choices,
    }
    virtual = pd.concat(virtual_rows, ignore_index=True) if virtual_rows else pd.DataFrame()
    return payload, virtual


def _cell_means(frame: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "test_skill", "validation_skill", "generator_calls",
        "recursively_evaluated_states", "stochastic_samples", "wall_time_seconds",
    ]
    return (
        frame.groupby(["pde", "d", "n", "M", "method"], as_index=False)[columns]
        .mean()
    )


def _at_cost(cells: pd.DataFrame, method: str, cost: float) -> tuple[float, dict[str, Any] | None]:
    rows = cells[(cells["method"] == method) & (cells["generator_calls"] <= cost)]
    if not len(rows):
        return float("nan"), None
    best = rows.loc[rows["test_skill"].idxmin()]
    return float(best["test_skill"]), {
        "n": int(best["n"]), "M": int(best["M"]),
        "actual_generator_calls": float(best["generator_calls"]),
        "wall_time_seconds": float(best["wall_time_seconds"]),
        "stochastic_samples": float(best["stochastic_samples"]),
    }


def _bootstrap_frontiers(
    frame: pd.DataFrame,
    levels: np.ndarray,
    methods: list[str],
    *,
    pde: str,
    d: int,
    draws: int = 1000,
) -> dict[tuple[str, int], tuple[float, float]]:
    reps = sorted(frame["rep"].unique())
    if not reps:
        return {}
    rng = np.random.default_rng(
        np.random.SeedSequence([BASE_SEED, stable_seed(pde), d, 1000])
    )
    results = {(method, index): [] for method in methods for index in range(len(levels))}
    for _ in range(draws):
        sampled = rng.choice(reps, len(reps), replace=True)
        pieces = []
        for synthetic, source_rep in enumerate(sampled):
            part = frame[frame["rep"] == source_rep].copy()
            part["bootstrap_rep"] = synthetic
            pieces.append(part)
        sampled_frame = pd.concat(pieces, ignore_index=True)
        means = (
            sampled_frame.groupby(["n", "M", "method"], as_index=False)[
                ["test_skill", "generator_calls"]
            ].mean()
        )
        for method in methods:
            for index, level in enumerate(levels):
                value, _ = _at_cost(means, method, float(level))
                results[(method, index)].append(value)
    intervals = {}
    for key, values in results.items():
        finite = np.asarray(values, dtype=np.float64)
        finite = finite[np.isfinite(finite)]
        intervals[key] = (
            (float(np.quantile(finite, 0.025)), float(np.quantile(finite, 0.975)))
            if len(finite) else (float("nan"), float("nan"))
        )
    return intervals


def frontier_levels(s1: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (pde, d), group in s1.groupby(["pde", "d"], sort=True):
        cells = _cell_means(group)
        primary = primary_method(str(pde)).name
        bounds = []
        for method in ("raw", primary):
            costs = cells.loc[
                (cells["method"] == method) & (cells["generator_calls"] > 0),
                "generator_calls",
            ]
            if len(costs):
                bounds.append((float(costs.min()), float(costs.max())))
        if len(bounds) != 2:
            continue
        lower = max(value[0] for value in bounds)
        upper = min(value[1] for value in bounds)
        if lower > upper:
            continue
        levels = np.geomspace(lower, upper, 10)
        methods = sorted(cells["method"].unique())
        intervals = _bootstrap_frontiers(
            group, levels, methods, pde=str(pde), d=int(d), draws=1000
        )
        for index, level in enumerate(levels):
            for method in methods:
                value, selected = _at_cost(cells, method, float(level))
                low, high = intervals[(method, index)]
                rows.append({
                    "pde": pde, "d": int(d), "level_index": index + 1,
                    "cost_level": float(level), "method": method,
                    "frontier_skill": value, "bootstrap_low": low,
                    "bootstrap_high": high,
                    "selected_n": None if selected is None else selected["n"],
                    "selected_M": None if selected is None else selected["M"],
                    "selected_generator_calls": None if selected is None else selected["actual_generator_calls"],
                    "selected_wall_time_seconds": None if selected is None else selected["wall_time_seconds"],
                    "selected_stochastic_samples": None if selected is None else selected["stochastic_samples"],
                    "bootstrap_draws": 1000,
                })
    return pd.DataFrame(rows)


def _overall(details: list[dict[str, Any]]) -> str:
    evaluated = [item for item in details if item.get("passed") is not None]
    if any(item.get("passed") is False for item in evaluated):
        return "FAIL"
    if len(evaluated) != len(details) or not details:
        return "NOT EVALUATED"
    return "PASS"


def evaluate_criteria(
    s1: pd.DataFrame,
    levels: pd.DataFrame,
    s2: pd.DataFrame,
    s3: pd.DataFrame,
    gates: dict[str, Any],
) -> dict[str, Any]:
    criteria: dict[str, Any] = {}

    pf1 = []
    for pde in sorted(s1["pde"].unique()):
        primary = primary_method(str(pde)).name
        for d in (100, 400):
            subset = levels[(levels["pde"] == pde) & (levels["d"] == d)]
            pivot = subset.pivot(index="level_index", columns="method", values="frontier_skill")
            if "raw" not in pivot or primary not in pivot or len(pivot) != 10:
                pf1.append({"pde": pde, "d": d, "passed": None, "reason": "frontier incomplete"})
                continue
            count = int(np.count_nonzero(pivot[primary] <= 0.8 * pivot["raw"]))
            pf1.append({"pde": pde, "d": d, "levels_passing": count, "required": 8, "passed": count >= 8})
    criteria["PF-1"] = {"verdict": _overall(pf1), "details": pf1}

    cells = _cell_means(s1)
    pf2 = []
    for pde in ("P1", "C1", "C2-convex", "C3"):
        primary = primary_method(pde).name
        groups = {}
        positive_ranges = []
        for d in (20, 400):
            group = cells[(cells["pde"] == pde) & (cells["d"] == d)]
            groups[d] = group
            for method in ("raw", primary):
                costs = group.loc[(group["method"] == method) & (group["generator_calls"] > 0), "generator_calls"]
                if len(costs):
                    positive_ranges.append((float(costs.min()), float(costs.max())))
        if len(positive_ranges) != 4:
            pf2.append({"pde": pde, "passed": None, "reason": "d=20/400 frontier incomplete"})
            continue
        lower = max(item[0] for item in positive_ranges)
        upper = min(item[1] for item in positive_ranges)
        if lower > upper:
            pf2.append({"pde": pde, "passed": None, "reason": "no shared cost range"})
            continue
        median_cost = math.sqrt(lower * upper)
        raw20, _ = _at_cost(groups[20], "raw", median_cost)
        raw400, _ = _at_cost(groups[400], "raw", median_cost)
        primary20, _ = _at_cost(groups[20], primary, median_cost)
        primary400, _ = _at_cost(groups[400], primary, median_cost)
        raw_ratio = raw400 / raw20
        primary_ratio = primary400 / primary20
        pf2.append({
            "pde": pde, "median_shared_cost": median_cost,
            "raw_d400_over_d20": raw_ratio,
            "primary_d400_over_d20": primary_ratio,
            "passed": bool(raw_ratio >= 2.0 and primary_ratio <= 1.3),
        })
    criteria["PF-2"] = {"verdict": _overall(pf2), "details": pf2}

    def paired_cell_test(pde: str, comparator: str) -> dict[str, Any]:
        data = s1[s1["pde"] == pde]
        primary = primary_method(pde).name
        detail = []
        for key, group in data.groupby(["d", "n", "M"]):
            left = group[group["method"] == primary].set_index("rep")["test_skill"]
            right = group[group["method"] == comparator].set_index("rep")["test_skill"]
            common = left.index.intersection(right.index)
            if len(common) != 10:
                detail.append({"d": key[0], "n": key[1], "M": key[2], "wins": None, "passed": None})
                continue
            wins = int(np.count_nonzero(left.loc[common].to_numpy() < right.loc[common].to_numpy()))
            detail.append({"d": key[0], "n": key[1], "M": key[2], "wins": wins, "passed": wins >= 7})
        evaluated = [item for item in detail if item["passed"] is not None]
        fraction = float(np.mean([item["passed"] for item in evaluated])) if evaluated else None
        return {"pde": pde, "passing_cell_fraction": fraction, "required": 0.75, "cells": detail,
                "passed": None if fraction is None else fraction >= 0.75}

    pf3 = [paired_cell_test(pde, "centre") for pde in ("P1", "C1", "C2-convex", "C3")]
    p4 = paired_cell_test("P4", "centre")
    p4["expected_failure_replicated"] = None if p4["passed"] is None else not p4["passed"]
    pf3_combined = pf3 + [{"pde": "P4 expected failure", "passed": p4["expected_failure_replicated"]}]
    criteria["PF-3"] = {"verdict": _overall(pf3_combined), "details": pf3, "P4_replication": p4}

    pf4 = [paired_cell_test(pde, "best_shrink") for pde in ("P1", "C1", "C2-convex", "C3")]
    criteria["PF-4"] = {"verdict": _overall(pf4), "details": pf4}

    t1 = []
    raw2 = s2[s2["method"] == "raw"]
    bias_means = raw2.groupby(["pde", "d", "n", "M"])["generator_bias"].mean()
    for d in (20, 100, 400):
        for n, M in ((2, 8), (2, 32)):
            key_values = {}
            try:
                for pde in ("C2-convex", "C2-cancel", "C2-flip"):
                    key_values[pde] = float(bias_means.loc[(pde, d, n, M)])
            except KeyError:
                t1.append({"d": d, "n": n, "M": M, "passed": None, "reason": "missing S2 rows"})
                continue
            limit = 0.25 * min(abs(key_values["C2-convex"]), abs(key_values["C2-flip"]))
            passed = key_values["C2-convex"] < 0 and key_values["C2-flip"] > 0 and abs(key_values["C2-cancel"]) <= limit
            t1.append({"d": d, "n": n, "M": M, "bias": key_values, "cancel_limit": limit, "passed": bool(passed)})
    criteria["T-1"] = {"verdict": _overall(t1), "details": t1}

    t2 = []
    cell_skill = raw2.groupby(["pde", "d", "n", "M"])["test_skill"].mean()
    for pde in ("C2-convex", "C2-cancel", "C2-flip"):
        try:
            best20 = float(cell_skill.loc[(pde, 20)].min())
            best400 = float(cell_skill.loc[(pde, 400)].min())
            ratio = best400 / best20
            passed = ratio <= 1.3 if pde == "C2-cancel" else ratio >= 2.0
            t2.append({"pde": pde, "best_d20": best20, "best_d400": best400, "ratio": ratio, "passed": bool(passed)})
        except KeyError:
            t2.append({"pde": pde, "passed": None, "reason": "missing S2 rows"})
    criteria["T-2"] = {"verdict": _overall(t2), "details": t2}

    if ONESTEP_PATH.exists():
        onestep = json.loads(ONESTEP_PATH.read_text(encoding="utf-8"))
        t3_pass = bool(onestep["T3"]["passed"])
        criteria["T-3"] = {"verdict": "PASS" if t3_pass else "FAIL", "details": onestep["T3"]}
    else:
        criteria["T-3"] = {"verdict": "NOT EVALUATED", "details": {"reason": "s2_onestep_check.json missing"}}

    dim = []
    means3 = s3.groupby(["pde", "d", "n", "M", "method"])["test_skill"].mean()
    for pde in ("P1", "C1", "C2-convex", "C3"):
        primary = primary_method(pde).name
        for n, M in ((2, 32), (3, 10)):
            try:
                low = float(means3.loc[(pde, 20, n, M, primary)])
                high = float(means3.loc[(pde, 400, n, M, primary)])
                ratio = high / low
                dim.append({"pde": pde, "n": n, "M": M, "d400_over_d20": ratio, "passed": ratio <= 1.3})
            except KeyError:
                dim.append({"pde": pde, "n": n, "M": M, "passed": None, "reason": "missing S3 rows"})
    criteria["DIM"] = {"verdict": _overall(dim), "details": dim}

    lqg = gates["candidates"].get("LQG", {})
    criteria["LQG-G"] = {
        "verdict": "OBSERVED",
        "details": {
            "predicted": "may fail G2",
            "gate_status": lqg.get("gate_status"),
            "failed_gates": lqg.get("failed_gates"),
            "prediction_matched": "G2" in lqg.get("failed_gates", []),
        },
    }
    return criteria


def make_figures(s1: pd.DataFrame, levels: pd.DataFrame, s2: pd.DataFrame, s3: pd.DataFrame) -> list[str]:
    FIGURE_ROOT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"figure.dpi": 140, "savefig.dpi": 180, "axes.grid": True, "grid.alpha": 0.25})
    made = []
    colors = {"raw": "#444444", "centre": "#d95f02", "best_shrink": "#7570b3", "f_zero": "#1b9e77"}
    for pde in sorted(s1["pde"].unique()):
        primary = primary_method(str(pde)).name
        methods = ["raw", primary, "centre", "best_shrink", "f_zero"]
        fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0), sharey=False)
        for axis, d in zip(axes, (20, 100, 400)):
            panel = levels[(levels["pde"] == pde) & (levels["d"] == d)]
            for method in methods:
                data = panel[panel["method"] == method].sort_values("cost_level")
                if not len(data):
                    continue
                label = "primary: " + method if method == primary else method
                axis.plot(data["cost_level"], data["frontier_skill"], marker="o", ms=3, label=label,
                          color=colors.get(method))
                axis.fill_between(data["cost_level"], data["bootstrap_low"], data["bootstrap_high"], alpha=0.12,
                                  color=colors.get(method))
            axis.set_xscale("log")
            axis.set_yscale("log")
            axis.set_title(f"{pde}, d={d}")
            axis.set_xlabel("generator calls")
            axis.set_ylabel("test skill")
        handles, labels = axes[0].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, loc="upper center", ncol=min(5, len(handles)), frameon=False)
        fig.tight_layout(rect=(0, 0, 1, 0.90))
        path = FIGURE_ROOT / f"s1_frontier_{str(pde).lower().replace('-', '_')}.png"
        fig.savefig(path, bbox_inches="tight")
        plt.close(fig)
        made.append(path.relative_to(RESULTS_ROOT).as_posix())

    raw = s2[(s2["method"] == "raw") & (s2["n"] == 2) & (s2["M"].isin([8, 32]))]
    means = raw.groupby(["pde", "d", "M"])["generator_bias"].mean().reset_index()
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8), sharey=True)
    for axis, M in zip(axes, (8, 32)):
        for pde, group in means[means["M"] == M].groupby("pde"):
            axis.plot(group["d"], group["generator_bias"], marker="o", label=pde)
        axis.axhline(0.0, color="black", lw=0.8)
        axis.set_title(f"C2 raw bias, n=2, M={M}")
        axis.set_xlabel("dimension")
        axis.set_ylabel("mean generator bias")
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    bias_path = FIGURE_ROOT / "s2_c2_bias_sign.png"
    fig.savefig(bias_path, bbox_inches="tight")
    plt.close(fig)
    made.append(bias_path.relative_to(RESULTS_ROOT).as_posix())

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
    for axis, pde in zip(axes.flat, ("P1", "C1", "C2-convex", "C3")):
        primary = primary_method(pde).name
        data = s3[(s3["pde"] == pde) & (s3["method"] == primary)]
        means = data.groupby(["d", "n", "M"])["test_skill"].mean().reset_index()
        for (n, M), group in means.groupby(["n", "M"]):
            axis.plot(group["d"], group["test_skill"], marker="o", label=f"({n},{M})")
        axis.set_title(pde)
        axis.set_yscale("log")
        axis.set_xlabel("dimension")
        axis.set_ylabel("primary test skill")
        axis.legend(frameon=False)
    fig.tight_layout()
    dim_path = FIGURE_ROOT / "s3_dimension.png"
    fig.savefig(dim_path, bbox_inches="tight")
    plt.close(fig)
    made.append(dim_path.relative_to(RESULTS_ROOT).as_posix())
    return made


def _format(value: Any) -> str:
    if value is None:
        return "—"
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def _table(headers: list[str], rows: list[list[Any]]) -> str:
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    lines.extend("| " + " | ".join(_format(value) for value in row) + " |" for row in rows)
    return "\n".join(lines)


def write_report(summary: dict[str, Any], gates: dict[str, Any]) -> None:
    criteria = summary["criteria"]
    outcome_rows = [[key, item["verdict"]] for key, item in criteria.items()]
    failures = [key for key, item in criteria.items() if item["verdict"] == "FAIL"]
    not_evaluated = [key for key, item in criteria.items() if item["verdict"] == "NOT EVALUATED"]
    gate_rows = []
    for pde, item in gates["candidates"].items():
        status = item.get("gate_status", {})
        gate_rows.append(
            [pde]
            + [
                "PASS" if status.get(gate) else "FAIL"
                for gate in ("G1", "G2", "G3", "G4", "G5", "G6")
            ]
        )

    lines = [
        "# Mechanism suite round 3 report",
        "",
        "## Outcome table",
        "",
        _table(["criterion", "verdict"], outcome_rows),
        "",
        "Round 3 was pre-registered after inspecting rounds 1–2 and was evaluated on the fresh base seed `20261207` and fresh 1,200-point sets. Value skill on the held-out 80% test split determines every verdict; gradient errors are diagnostic only.",
        "",
        "## Pre-registered criteria (verbatim) and verdicts",
        "",
    ]
    for key, text in CRITERIA.items():
        lines.extend([f"### {key} — {criteria[key]['verdict']}", "", f"> {text}", "", "Mechanical details are recorded in [`analysis_summary.json`](../results/mechanism_suite_r3/analysis_summary.json).", ""])
    lines.extend([
        "## S0 gates",
        "",
        _table(["PDE", "G1", "G2", "G3", "G4", "G5", "G6"], gate_rows),
        "",
        "P1 and P4 are carried unchanged from the round-1 gate artifact, as pre-registered. LQG is exploratory; its outcome is compared only with the stated advance prediction.",
        "",
        "## Failures and adverse findings",
        "",
    ])
    if failures:
        lines.extend([f"- {key} failed its frozen mechanical threshold." for key in failures])
    else:
        lines.append("- No evaluated confirmatory criterion failed.")
    for key in not_evaluated:
        lines.append(f"- {key} was not evaluated; the reason is recorded in `analysis_summary.json`.")
    for pde, item in gates["candidates"].items():
        if not item.get("admitted", False):
            lines.append(f"- {pde} failed gate(s): {', '.join(item.get('failed_gates', [])) or 'not reported'}; downstream studies were stopped for that PDE.")
    dropped = summary.get("d400_dropped_configs", [])
    if dropped:
        lines.append(f"- The registered d=400 30-minute rule dropped these (n,M) cells: {dropped}.")
    lines.extend([
        "",
        "## Artifacts and figures",
        "",
        "- [S1 repetition rows](../results/mechanism_suite_r3/s1_frontier.csv) and [frontier levels](../results/mechanism_suite_r3/s1_frontier_levels.csv)",
        "- [S2 theory rows](../results/mechanism_suite_r3/s2_theory.csv) and [one-step check](../results/mechanism_suite_r3/s2_onestep_check.json)",
        "- [S3 dimension rows](../results/mechanism_suite_r3/s3_dimension.csv) and [S4 geometry rows](../results/mechanism_suite_r3/s4_geometry.csv)",
        "- [Validation-only tuning choices](../results/mechanism_suite_r3/tuning_choices.json)",
        "",
        "Generated figures: " + ", ".join(f"[`{Path(path).name}`](../results/mechanism_suite_r3/{path})" for path in summary.get("figures", [])) + ".",
        "",
        "## Provenance and work totals",
        "",
        f"- Frozen pre-registration commit: `{PREREGISTRATION_COMMIT}`.",
        f"- Experiment code commit(s) recorded by artifacts: {', '.join(f'`{value}`' for value in summary['code_commits']) or 'none'}.",
        f"- Analysis commit at generation: `{summary['analysis_environment'].get('git_commit')}` on `{summary['analysis_environment'].get('git_branch')}`.",
        f"- Unchanged `FullHistoryMLP` SHA-256: `{summary['full_history_sha256']}`.",
        f"- Environment: Python {summary['analysis_environment']['python']}; NumPy {summary['analysis_environment']['numpy']}; SciPy {summary['analysis_environment']['scipy']}; {summary['analysis_environment']['platform']}.",
        f"- Unique computed artifacts represented: {summary['work_totals']['unique_artifacts']}; generator calls: {summary['work_totals']['generator_calls']}; recursively evaluated states: {summary['work_totals']['recursively_evaluated_states']}; stochastic samples: {summary['work_totals']['stochastic_samples']}; summed worker time: {summary['work_totals']['wall_time_seconds']/3600:.2f} h; non-finite count: {summary['work_totals']['nonfinite_count']}.",
        "",
        "The manuscript, `FullHistoryMLP`, both earlier mechanism-suite modules, and all earlier result directories were left unchanged.",
    ])
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def analyze() -> dict[str, Any]:
    admission, gates = _admission()
    s1, missing1 = load_study("s1", admission)
    s2, missing2 = load_study("s2", admission)
    s3, missing3 = load_study("s3", admission)
    s4, missing4 = load_study("s4", admission)
    tuning, virtual = tune_shrink((s1, s2, s3))
    _write_json(TUNING_PATH, tuning)
    if len(virtual):
        s1_virtual = virtual[virtual["study"] == "s1"]
        s2_virtual = virtual[virtual["study"] == "s2"]
        s3_virtual = virtual[virtual["study"] == "s3"]
        s1 = pd.concat([s1, s1_virtual], ignore_index=True)
        s2 = pd.concat([s2, s2_virtual], ignore_index=True)
        s3 = pd.concat([s3, s3_virtual], ignore_index=True)
    levels = frontier_levels(s1)
    criteria = evaluate_criteria(s1, levels, s2, s3, gates)
    _write_csv(S1_PATH, s1)
    _write_csv(LEVELS_PATH, levels)
    _write_csv(S2_PATH, s2)
    _write_csv(S3_PATH, s3)
    _write_csv(S4_PATH, s4)
    figures = make_figures(s1, levels, s2, s3)

    all_frames = [frame for frame in (s1, s2, s3, s4) if len(frame)]
    combined = pd.concat(all_frames, ignore_index=True)
    unique = combined.drop_duplicates("artifact")
    code_commits = sorted(str(value) for value in unique["code_commit"].dropna().unique())
    pilot = json.loads(D400_PILOT_PATH.read_text(encoding="utf-8")) if D400_PILOT_PATH.exists() else {}
    summary = {
        "schema_version": 1,
        "criteria": criteria,
        "criteria_verbatim": CRITERIA,
        "missing_artifacts": {"s1": missing1, "s2": missing2, "s3": missing3, "s4": missing4},
        "d400_dropped_configs": pilot.get("dropped_configs", []),
        "figures": figures,
        "code_commits": code_commits,
        "full_history_sha256": hashlib.sha256(FULL_HISTORY_SOURCE.read_bytes()).hexdigest(),
        "work_totals": {
            "unique_artifacts": len(unique),
            "generator_calls": int(unique["generator_calls"].sum()),
            "recursively_evaluated_states": int(unique["recursively_evaluated_states"].sum()),
            "stochastic_samples": int(unique["stochastic_samples"].sum()),
            "wall_time_seconds": float(unique["wall_time_seconds"].sum()),
            "nonfinite_count": int(unique["nonfinite_count"].sum()),
        },
        "analysis_environment": {
            "python": platform.python_version(), "numpy": np.__version__,
            "scipy": scipy.__version__, "pandas": pd.__version__,
            "platform": platform.platform(), "cpu_count": os.cpu_count(),
            "git_branch": _git("branch", "--show-current"),
            "git_commit": _git("rev-parse", "HEAD"),
        },
        "preregistration_commit": PREREGISTRATION_COMMIT,
        "implementation_revision": IMPLEMENTATION_REVISION,
    }
    _write_json(SUMMARY_PATH, summary)
    write_report(summary, gates)
    print(json.dumps({key: value["verdict"] for key, value in criteria.items()}, indent=2))
    print(f"report: {REPORT_PATH}")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    analyze()


if __name__ == "__main__":
    main()
