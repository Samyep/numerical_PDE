"""Mechanical aggregation, verdicts, figures, and report generation."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess
import sys
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch

from .equations import fixed_points, make_s05_equation
from .protocol import (
    CALIBRATION_FREEZE_PATH,
    FIGURES_ROOT,
    RAW_ROOT,
    REPORT_PATH,
    RESULTS_ROOT,
    REVISION,
    ensure_directories,
    fh_tasks,
    write_json_atomic,
)


CRITERIA = {
    "H1": (
        "H1 (learning is decisive): RMSE(TA) <= 0.5 x RMSE(TA_a0) "
        "in all 4 (PDE, d) cells."
    ),
    "H2": (
        "H2 (more than one update matters): RMSE(TA) <= 0.8 x RMSE(a1), "
        "both under the same post-freeze common-random-number diagnostic, "
        "in at least 2 of the 4 cells. A failure means \"one policy-improvement "
        "step suffices here\" and must be stated plainly."
    ),
    "H3": (
        "H3 (beats raw Full-History MLP): RMSE(TA) <= 0.5 x the minimum "
        "RMSE over the entire frozen FH_raw grid, in all 4 cells."
    ),
    "H4": (
        "H4 (beats every Full-History variant): RMSE(TA) <= the minimum RMSE "
        "over all frozen FH_raw, FH_path and FH_rep cells (any cost), in at "
        "least 3 of the 4 cells."
    ),
    "H5": (
        "H5 (efficiency vs the best Full-History variant): in at least 2 of "
        "the 4 cells, the best Full-History family does not reach 2 x RMSE(TA) "
        "within its frozen grid, or needs more GPU seconds than TA's total to "
        "do so."
    ),
    "H6": (
        "H6 (no dependence on a fixed subspace): for B6, r_glob >= d/4 at "
        "the frozen calibration, and H1 holds on B6 at both d. Report r_glob "
        "for S05 alongside."
    ),
}

PROBLEMS = (("S05", 100), ("B6", 100), ("S05", 400), ("B6", 400))
METHOD_LABELS = {
    "raw": "FH_raw",
    "path": "FH_path",
    "double_path_base": "FH_rep",
}


def _repository_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _git(*args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=_repository_root(), check=True,
        capture_output=True, text=True,
    ).stdout.strip()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_rows(root: Path) -> list[dict[str, Any]]:
    rows = []
    if not root.is_dir():
        return rows
    for path in sorted(root.rglob("*.json")):
        if path.name.endswith((".failure.json", ".curtailed.json")):
            continue
        try:
            row = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        row["artifact_path"] = str(path)
        rows.append(row)
    return rows


def _curtailed_rows(root: Path) -> list[dict[str, Any]]:
    """Read compute-curtailment markers without treating them as estimates."""
    rows: list[dict[str, Any]] = []
    if not root.is_dir():
        return rows
    for path in sorted(root.rglob("*.curtailed.json")):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        task = payload.get("task", {})
        rows.append(
            {
                **task,
                "status": "curtailed_after_rep0",
                "reason": payload.get("reason"),
                "rep0_seconds": payload.get("rep0_seconds"),
                "artifact_path": str(path),
            }
        )
    return rows


def _actor_frames() -> tuple[pd.DataFrame, pd.DataFrame, dict[tuple[str, int, int], dict[str, Any]]]:
    evaluations = pd.DataFrame(_json_rows(RAW_ROOT / "actor_evaluations"))
    training: dict[tuple[str, int, int], dict[str, Any]] = {}
    actor_root = RAW_ROOT / "actor"
    if actor_root.is_dir():
        for path in sorted(actor_root.glob("training_*.json")):
            payload = json.loads(path.read_text(encoding="utf-8"))
            training[(payload["pde_id"], int(payload["dimension"]), int(payload["seed"]))] = payload
    if evaluations.empty:
        return evaluations, pd.DataFrame(), training
    numeric = [
        "value_rmse", "relative_l2", "skill", "bias", "mae",
        "gradient_vector_rmse", "gradient_relative_l2",
        "mc_rmse_contribution", "non_mc_rmse_proxy", "variance_reduction",
        "weight_ess_fraction", "inference_seconds", "peak_gpu_bytes",
        "nonfinite_count",
    ]
    group = [
        "pde_id", "dimension", "training_seed", "actor_request",
        "resolved_actor", "steps", "paths",
    ]
    aggregates: list[dict[str, Any]] = []
    for keys, data in evaluations.groupby(group, dropna=False):
        row = dict(zip(group, keys, strict=True))
        row["replicates"] = len(data)
        for column in numeric:
            values = pd.to_numeric(data[column], errors="coerce")
            row[f"mean_{column}"] = float(values.mean())
            row[f"std_{column}"] = float(values.std(ddof=1)) if len(values) > 1 else 0.0
        aggregates.append(row)
    return evaluations, pd.DataFrame(aggregates), training


def _fh_frames() -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = pd.DataFrame(_json_rows(RAW_ROOT / "fh"))
    if rows.empty:
        return rows, pd.DataFrame()
    groups = ["pde_id", "dimension", "method", "n", "M"]
    summaries: list[dict[str, Any]] = []
    for keys, data in rows.groupby(groups, dropna=False):
        row = dict(zip(groups, keys, strict=True))
        row["method_label"] = METHOD_LABELS[str(row["method"])]
        row["repetitions"] = len(data)
        row["expected_repetitions"] = 3
        row["complete"] = len(data) == 3
        reused = (
            data["reused"].eq(True)
            if "reused" in data.columns
            else pd.Series(False, index=data.index, dtype=bool)
        )
        row["reused_repetitions"] = int(reused.sum())
        for column in (
            "value_rmse", "skill", "value_relative_l2", "value_bias",
            "gradient_vector_rmse", "gradient_relative_l2", "wall_time_seconds",
            "peak_gpu_memory_bytes", "generator_calls", "gaussian_coordinates_drawn",
            "nonfinite_state_count", "nonfinite_generator_count",
        ):
            if column not in data:
                row[f"mean_{column}"] = float("nan")
                row[f"std_{column}"] = float("nan")
                continue
            values = pd.to_numeric(data[column], errors="coerce")
            row[f"mean_{column}"] = float(values.mean())
            row[f"std_{column}"] = float(values.std(ddof=1)) if values.notna().sum() > 1 else 0.0
        prediction_paths = [
            Path(value) for value in data.get("prediction_path", pd.Series(dtype=object)).dropna()
            if isinstance(value, str) and value
        ]
        row["mc_rmse_contribution"] = float("nan")
        row["non_mc_rmse_proxy"] = float("nan")
        if len(prediction_paths) >= 2 and all(path.is_file() for path in prediction_paths):
            predictions = []
            truth = None
            mask = None
            for path in prediction_paths:
                with np.load(path) as arrays:
                    predictions.append(arrays["prediction"][:, 0])
                    truth = arrays["truth"][:, 0]
                    mask = ~arrays["is_validation"]
            stack = np.stack(predictions)[:, mask]
            reference = np.asarray(truth)[mask]
            variance = np.var(stack, axis=0, ddof=1)
            mc = float(np.sqrt(np.mean(variance)))
            total = float(np.sqrt(np.mean((stack - reference[None, :]) ** 2)))
            row["mc_rmse_contribution"] = mc
            row["non_mc_rmse_proxy"] = math.sqrt(max(total * total - mc * mc, 0.0))
        summaries.append(row)
    return rows, pd.DataFrame(summaries)


def _lookup_actor(
    actor: pd.DataFrame, pde_id: str, d: int, name: str
) -> pd.Series | None:
    if actor.empty:
        return None
    selected = actor[
        (actor.pde_id == pde_id)
        & (actor.dimension == d)
        & (actor.training_seed == 0)
        & (actor.actor_request == name)
        & (actor.steps == 64)
    ]
    if len(selected) != 1 or int(selected.iloc[0].replicates) != 4:
        return None
    return selected.iloc[0]


def _ta_total_seconds(
    actor: pd.DataFrame,
    training: dict[tuple[str, int, int], dict[str, Any]],
    pde_id: str,
    d: int,
) -> float | None:
    row = _lookup_actor(actor, pde_id, d, "selected")
    record = training.get((pde_id, d, 0))
    if row is None or record is None or record.get("status") != "complete":
        return None
    timing = record["timing"]
    return float(
        timing["phase_total_seconds"]
        + timing["reference_free_selection_seconds"]
        + row["mean_inference_seconds"]
    )


def _graded_verdict(
    details: list[dict[str, Any]], required: int, total: int
) -> tuple[str, int, int]:
    evaluated = [row for row in details if row["met"] is not None]
    wins = sum(bool(row["met"]) for row in evaluated)
    missing = total - len(evaluated)
    if wins >= required:
        return "PASS", wins, len(evaluated)
    if wins + missing < required:
        return "FAIL", wins, len(evaluated)
    return "NOT EVALUATED (COMPUTE)", wins, len(evaluated)


def _all_verdict(details: list[dict[str, Any]], total: int) -> tuple[str, int, int]:
    evaluated = [row for row in details if row["met"] is not None]
    wins = sum(bool(row["met"]) for row in evaluated)
    if len(evaluated) < total:
        return "NOT EVALUATED (COMPUTE)", wins, len(evaluated)
    return ("PASS" if wins == total else "FAIL"), wins, len(evaluated)


def _s05_global_rank(d: int) -> float:
    equation = make_s05_equation(d)
    t, x, mask = fixed_points("S05", d)
    z = equation.exact_z(t[~mask], x[~mask])
    centered = z - np.mean(z, axis=0, keepdims=True)
    covariance = centered.T @ centered / len(centered)
    eigenvalues = np.maximum(np.linalg.eigvalsh(covariance), 0.0)
    return float(np.sum(eigenvalues) ** 2 / np.sum(eigenvalues**2))


def _criteria(
    actor: pd.DataFrame,
    fh: pd.DataFrame,
    training: dict[tuple[str, int, int], dict[str, Any]],
    calibration: pd.DataFrame,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    outcomes: dict[str, Any] = {}
    threshold_rows: list[dict[str, Any]] = []

    h1 = []
    h2 = []
    h3 = []
    h4 = []
    h5 = []
    for pde_id, d in PROBLEMS:
        ta = _lookup_actor(actor, pde_id, d, "selected")
        a0 = _lookup_actor(actor, pde_id, d, "a0")
        a1 = _lookup_actor(actor, pde_id, d, "a1")
        ta_rmse = None if ta is None else float(ta.mean_value_rmse)
        a0_rmse = None if a0 is None else float(a0.mean_value_rmse)
        a1_rmse = None if a1 is None else float(a1.mean_value_rmse)
        h1.append(
            {
                "pde_id": pde_id, "dimension": d, "ta_rmse": ta_rmse,
                "comparison_rmse": a0_rmse,
                "ratio": None if ta_rmse is None or a0_rmse is None else ta_rmse / a0_rmse,
                "met": None if ta_rmse is None or a0_rmse is None else ta_rmse <= 0.5 * a0_rmse,
            }
        )
        h2.append(
            {
                "pde_id": pde_id, "dimension": d, "ta_rmse": ta_rmse,
                "comparison_rmse": a1_rmse,
                "ratio": None if ta_rmse is None or a1_rmse is None else ta_rmse / a1_rmse,
                "met": None if ta_rmse is None or a1_rmse is None else ta_rmse <= 0.8 * a1_rmse,
            }
        )
        problem = fh[(fh.pde_id == pde_id) & (fh.dimension == d) & fh.complete] if not fh.empty else pd.DataFrame()
        raw = problem[problem.method == "raw"] if not problem.empty else pd.DataFrame()
        # A minimum over the *entire frozen grid* is not defined until all
        # registered cells have their three repetitions.
        raw_complete = len(raw) == 8
        whole_grid_complete = len(problem) == 22
        best_raw = float(raw.mean_value_rmse.min()) if raw_complete else None
        best_all = float(problem.mean_value_rmse.min()) if whole_grid_complete else None
        h3.append(
            {
                "pde_id": pde_id, "dimension": d, "ta_rmse": ta_rmse,
                "comparison_rmse": best_raw,
                "ratio": None if ta_rmse is None or best_raw is None else ta_rmse / best_raw,
                "met": None if ta_rmse is None or best_raw is None else ta_rmse <= 0.5 * best_raw,
            }
        )
        h4.append(
            {
                "pde_id": pde_id, "dimension": d, "ta_rmse": ta_rmse,
                "comparison_rmse": best_all,
                "ratio": None if ta_rmse is None or best_all is None else ta_rmse / best_all,
                "met": None if ta_rmse is None or best_all is None else ta_rmse <= best_all,
            }
        )
        ta_seconds = _ta_total_seconds(actor, training, pde_id, d)
        threshold = None if ta_rmse is None else 2.0 * ta_rmse
        reached_rows = (
            problem[problem.mean_value_rmse <= threshold]
            if threshold is not None and not problem.empty else pd.DataFrame()
        )
        fh_seconds = None if reached_rows.empty else float(reached_rows.mean_wall_time_seconds.min())
        met = None
        if threshold is not None and ta_seconds is not None and whole_grid_complete:
            met = reached_rows.empty or fh_seconds > ta_seconds
        detail = {
            "pde_id": pde_id, "dimension": d, "ta_rmse": ta_rmse,
            "threshold": threshold, "ta_total_seconds": ta_seconds,
            "fh_reached": not reached_rows.empty,
            "fh_min_seconds": fh_seconds, "met": met,
        }
        h5.append(detail)
        threshold_rows.append(detail)

    verdict, wins, evaluated = _all_verdict(h1, 4)
    outcomes["H1"] = {"verdict": verdict, "wins": wins, "evaluated": evaluated, "required": 4, "details": h1}
    verdict, wins, evaluated = _graded_verdict(h2, 2, 4)
    outcomes["H2"] = {"verdict": verdict, "wins": wins, "evaluated": evaluated, "required": 2, "details": h2}
    verdict, wins, evaluated = _all_verdict(h3, 4)
    outcomes["H3"] = {"verdict": verdict, "wins": wins, "evaluated": evaluated, "required": 4, "details": h3}
    verdict, wins, evaluated = _graded_verdict(h4, 3, 4)
    outcomes["H4"] = {"verdict": verdict, "wins": wins, "evaluated": evaluated, "required": 3, "details": h4}
    verdict, wins, evaluated = _graded_verdict(h5, 2, 4)
    outcomes["H5"] = {"verdict": verdict, "wins": wins, "evaluated": evaluated, "required": 2, "details": h5}

    h6_details = []
    for d in (100, 400):
        selected = calibration[(calibration.dimension == d) & calibration.qualifies]
        freeze = json.loads(CALIBRATION_FREEZE_PATH.read_text(encoding="utf-8"))["selections"][str(d)]
        beta, gamma = freeze["selected_beta"], freeze["selected_gamma"]
        selected = calibration[(calibration.dimension == d) & (calibration.beta == beta) & (calibration.gamma == gamma)]
        rank = None if selected.empty else float(selected.iloc[0].r_glob)
        h1_row = next(row for row in h1 if row["pde_id"] == "B6" and row["dimension"] == d)
        met = None if rank is None or h1_row["met"] is None else rank >= d / 4.0 and bool(h1_row["met"])
        h6_details.append(
            {
                "pde_id": "B6", "dimension": d, "r_glob": rank,
                "rank_threshold": d / 4.0, "h1_met": h1_row["met"], "met": met,
            }
        )
    verdict, wins, evaluated = _all_verdict(h6_details, 2)
    outcomes["H6"] = {"verdict": verdict, "wins": wins, "evaluated": evaluated, "required": 2, "details": h6_details}
    return outcomes, threshold_rows


def _frontier(actor: pd.DataFrame, fh: pd.DataFrame, training: dict[tuple[str, int, int], dict[str, Any]]) -> pd.DataFrame:
    candidates: list[dict[str, Any]] = []
    if not actor.empty:
        for _, row in actor[(actor.training_seed == 0) & actor.actor_request.str.match(r"a[0-4]")].iterrows():
            record = training.get((row.pde_id, int(row.dimension), 0))
            if record is None:
                continue
            iteration = int(str(row.actor_request)[1:])
            training_seconds = sum(
                phase["phase_seconds"] for phase in record.get("phases", [])[:iteration]
            )
            candidates.append(
                {
                    "pde_id": row.pde_id, "dimension": int(row.dimension),
                    "family": "TA", "cell": row.actor_request,
                    "cost_seconds": training_seconds + float(row.mean_inference_seconds),
                    "rmse": float(row.mean_value_rmse),
                }
            )
    if not fh.empty:
        for _, row in fh[fh.complete].iterrows():
            candidates.append(
                {
                    "pde_id": row.pde_id, "dimension": int(row.dimension),
                    "family": METHOD_LABELS[row.method],
                    "cell": f"n{int(row.n)} M{int(row.M)}",
                    "cost_seconds": float(row.mean_wall_time_seconds),
                    "rmse": float(row.mean_value_rmse),
                }
            )
    frame = pd.DataFrame(candidates)
    if frame.empty:
        return frame
    output = []
    for _, data in frame.groupby(["pde_id", "dimension", "family"]):
        best = float("inf")
        for _, row in data.sort_values(["cost_seconds", "rmse"]).iterrows():
            item = row.to_dict()
            item["is_frontier"] = bool(row.rmse < best)
            best = min(best, float(row.rmse))
            item["envelope_rmse"] = best
            output.append(item)
    return pd.DataFrame(output)


def _resource_rows(
    actor: pd.DataFrame,
    fh: pd.DataFrame,
    training: dict[tuple[str, int, int], dict[str, Any]],
) -> pd.DataFrame:
    rows = []
    for (pde_id, d, seed), record in training.items():
        timing = record.get("timing", {})
        rows.append(
            {
                "pde_id": pde_id, "dimension": d, "method": "TA-training",
                "seed_or_repetitions": seed,
                "label_seconds": timing.get("label_generation_seconds"),
                "fit_seconds": timing.get("fitting_seconds"),
                "selection_seconds": timing.get("reference_free_selection_seconds"),
                "inference_seconds": None,
                "total_seconds": (
                    timing.get("phase_total_seconds", 0.0)
                    + timing.get("reference_free_selection_seconds", 0.0)
                ),
                "peak_gpu_bytes": max(
                    [
                        phase["training"].get("peak_gpu_bytes") or 0
                        for phase in record.get("phases", [])
                    ]
                    or [0]
                ),
            }
        )
    if not actor.empty:
        for _, row in actor.iterrows():
            rows.append(
                {
                    "pde_id": row.pde_id, "dimension": int(row.dimension),
                    "method": f"TA-{row.actor_request}-steps{int(row.steps)}",
                    "seed_or_repetitions": int(row.replicates),
                    "label_seconds": None, "fit_seconds": None,
                    "selection_seconds": None,
                    "inference_seconds": row.mean_inference_seconds,
                    "total_seconds": row.mean_inference_seconds,
                    "peak_gpu_bytes": row.mean_peak_gpu_bytes,
                }
            )
    if not fh.empty:
        for _, row in fh.iterrows():
            rows.append(
                {
                    "pde_id": row.pde_id, "dimension": int(row.dimension),
                    "method": f"{METHOD_LABELS[row.method]}-n{int(row.n)}-M{int(row.M)}",
                    "seed_or_repetitions": int(row.repetitions),
                    "label_seconds": None, "fit_seconds": None,
                    "selection_seconds": None,
                    "inference_seconds": row.mean_wall_time_seconds,
                    "total_seconds": row.mean_wall_time_seconds,
                    "peak_gpu_bytes": row.mean_peak_gpu_memory_bytes,
                }
            )
    return pd.DataFrame(rows)


def _plots(actor: pd.DataFrame, fh: pd.DataFrame, frontier: pd.DataFrame, calibration: pd.DataFrame) -> None:
    FIGURES_ROOT.mkdir(parents=True, exist_ok=True)
    colors = {"TA": "#1f77b4", "FH_raw": "#d62728", "FH_path": "#2ca02c", "FH_rep": "#ff7f0e"}
    markers = {"TA": "o", "FH_raw": "x", "FH_path": "s", "FH_rep": "^"}
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.6), constrained_layout=True)
    for axis, (pde_id, d) in zip(axes.flat, PROBLEMS, strict=True):
        if not frontier.empty:
            data = frontier[(frontier.pde_id == pde_id) & (frontier.dimension == d)]
            for family in ("TA", "FH_raw", "FH_path", "FH_rep"):
                family_data = data[data.family == family].sort_values("cost_seconds")
                if not family_data.empty:
                    axis.scatter(
                        family_data.cost_seconds, family_data.rmse,
                        label=family, color=colors[family], marker=markers[family], s=38,
                    )
                    envelope = family_data[family_data.is_frontier]
                    axis.plot(envelope.cost_seconds, envelope.envelope_rmse, color=colors[family], alpha=0.75)
        axis.set_xscale("log"); axis.set_yscale("log")
        axis.set_title(f"{pde_id}, d={d}")
        axis.set_xlabel("GPU seconds"); axis.set_ylabel("Value RMSE")
        axis.grid(True, which="both", alpha=0.22)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=4, frameon=False)
    fig.savefig(FIGURES_ROOT / "rmse_vs_gpu_seconds.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.2), constrained_layout=True)
    for axis, (pde_id, d) in zip(axes.flat, PROBLEMS, strict=True):
        if not actor.empty:
            data = actor[
                (actor.pde_id == pde_id)
                & (actor.dimension == d)
                & (actor.training_seed == 0)
                & (actor.steps == 64)
                & actor.actor_request.str.match(r"a[0-4]")
            ].copy()
            if not data.empty:
                data["iteration"] = data.actor_request.str[1:].astype(int)
                data = data.sort_values("iteration")
                axis.errorbar(
                    data.iteration, data.mean_value_rmse,
                    yerr=data.std_value_rmse, marker="o", capsize=3, color=colors["TA"],
                )
        axis.set_yscale("log")
        axis.set_xticks(range(5)); axis.set_title(f"{pde_id}, d={d}")
        axis.set_xlabel("Actor iteration"); axis.set_ylabel("Value RMSE")
        axis.grid(True, which="both", alpha=0.22)
    fig.savefig(FIGURES_ROOT / "actor_iterations.png", dpi=220)
    plt.close(fig)

    if not calibration.empty:
        fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.9), constrained_layout=True)
        freeze = json.loads(CALIBRATION_FREEZE_PATH.read_text(encoding="utf-8"))["selections"]
        for axis, d in zip(axes, (100, 400), strict=True):
            data = calibration[calibration.dimension == d]
            matrix = data.pivot(index="gamma", columns="beta", values="heat_only_skill").sort_index(ascending=False)
            image = axis.imshow(matrix.values, aspect="auto", cmap="viridis")
            axis.set_xticks(range(len(matrix.columns)), labels=[f"{v:g}" for v in matrix.columns])
            axis.set_yticks(range(len(matrix.index)), labels=[f"{v:g}" for v in matrix.index])
            axis.set_xlabel(r"$\beta$"); axis.set_ylabel(r"$\gamma$")
            axis.set_title(f"B6 calibration, d={d}")
            chosen = freeze[str(d)]
            for yi, gamma in enumerate(matrix.index):
                for xi, beta in enumerate(matrix.columns):
                    row = data[(data.beta == beta) & (data.gamma == gamma)].iloc[0]
                    text = f"{row.heat_only_skill:.2f}\nr={row.r_glob:.0f}"
                    selected = beta == chosen["selected_beta"] and gamma == chosen["selected_gamma"]
                    axis.text(
                        xi, yi, text + (" *" if selected else ""), ha="center", va="center",
                        color="white" if row.heat_only_skill > data.heat_only_skill.median() else "black",
                        fontsize=8, fontweight="bold" if selected else "normal",
                    )
            fig.colorbar(image, ax=axis, label="Heat-only skill")
        fig.savefig(FIGURES_ROOT / "b6_calibration_heatmap.png", dpi=220)
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
    outcomes: dict[str, Any], actor: pd.DataFrame, fh: pd.DataFrame,
    training: dict[tuple[str, int, int], dict[str, Any]], calibration: pd.DataFrame,
    validation: dict[str, Any], s05_ranks: dict[int, float], provenance: dict[str, Any],
    curtailed: pd.DataFrame,
) -> None:
    lines = [
        "# High-effective-dimension, symmetry-free benchmark",
        "",
        "## Outcome table",
        "",
        "| Criterion | Verdict | Wins / evaluated | Required |",
        "|---|---|---:|---:|",
    ]
    for key in CRITERIA:
        item = outcomes[key]
        lines.append(f"| {key} | **{item['verdict']}** | {item['wins']} / {item['evaluated']} | {item['required']} |")
    lines.extend(["", "## Pre-registered criteria (verbatim)", ""])
    for key, text in CRITERIA.items():
        item = outcomes[key]
        lines.extend([f"### {key}", "", f"> {text}", "", f"**Verdict: {item['verdict']}.**", ""])
        if item.get("details"):
            if key in {"H1", "H2", "H3", "H4"}:
                lines.extend([
                    "| PDE | d | TA RMSE | comparator RMSE | TA / comparator | met |",
                    "|---|---:|---:|---:|---:|---|",
                ])
                for row in item["details"]:
                    met = "not evaluated" if row.get("met") is None else ("yes" if row["met"] else "no")
                    lines.append(
                        f"| {row['pde_id']} | {row['dimension']} | "
                        f"{_number(row.get('ta_rmse'))} | "
                        f"{_number(row.get('comparison_rmse'))} | "
                        f"{_number(row.get('ratio'))} | {met} |"
                    )
            elif key == "H5":
                lines.extend([
                    "| PDE | d | TA RMSE | 2 x TA RMSE | FH reached | FH GPU s | TA total GPU s | met |",
                    "|---|---:|---:|---:|---|---:|---:|---|",
                ])
                for row in item["details"]:
                    met = "not evaluated" if row.get("met") is None else ("yes" if row["met"] else "no")
                    lines.append(
                        f"| {row['pde_id']} | {row['dimension']} | "
                        f"{_number(row.get('ta_rmse'))} | {_number(row.get('threshold'))} | "
                        f"{'yes' if row.get('fh_reached') else 'no'} | "
                        f"{_number(row.get('fh_min_seconds'))} | "
                        f"{_number(row.get('ta_total_seconds'))} | {met} |"
                    )
            else:
                lines.extend([
                    "| PDE | d | r_glob | d / 4 | H1 met | met |",
                    "|---|---:|---:|---:|---|---|",
                ])
                for row in item["details"]:
                    met = "not evaluated" if row.get("met") is None else ("yes" if row["met"] else "no")
                    lines.append(
                        f"| {row['pde_id']} | {row['dimension']} | "
                        f"{_number(row.get('r_glob'))} | "
                        f"{_number(row.get('rank_threshold'))} | "
                        f"{'yes' if row.get('h1_met') else 'no'} | {met} |"
                    )
            lines.append("")

    lines.extend(["## B6 calibration (reference-only, frozen before solver runs)", "", "| d | beta | gamma | ESS med | z RMS | f CV | heat skill | r glob | outside top 10 | r loc | qualifies | selected |", "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|"])
    freeze = json.loads(CALIBRATION_FREEZE_PATH.read_text(encoding="utf-8"))["selections"]
    for _, row in calibration.sort_values(["dimension", "beta", "gamma"]).iterrows():
        chosen = freeze[str(int(row.dimension))]
        selected = row.beta == chosen["selected_beta"] and row.gamma == chosen["selected_gamma"]
        lines.append(
            f"| {int(row.dimension)} | {row.beta:g} | {row.gamma:g} | {_number(row.ESS_med)} | {_number(row.z_rms)} | {_number(row.f_cv)} | {_number(row.heat_only_skill)} | {_number(row.r_glob)} | {_number(row.z_variance_outside_top10)} | {_number(row.r_loc)} | {'yes' if row.qualifies else 'no'} | {'yes' if selected else ''} |"
        )
    lines.extend(["", "S05 global gradient-covariance effective ranks (same definition as `r_glob`): " + ", ".join(f"d={d}: {_number(rank)}" for d, rank in s05_ranks.items()) + ".", ""])

    lines.extend(["## Validation", ""])
    for key, item in validation.get("stages", {}).items():
        lines.append(f"- **{key}: {item['verdict']}** ({item['checks']} checks, {item['failures']} failures).")
    lines.extend(["", "## Actor selection and iteration diagnostics", ""])
    if actor.empty:
        lines.append("Not evaluated (compute).")
    else:
        lines.extend(["| PDE | d | seed | actor | resolved | steps | RMSE | MC contribution | GPU s |", "|---|---:|---:|---|---|---:|---:|---:|---:|"])
        for _, row in actor.sort_values(["pde_id", "dimension", "training_seed", "steps", "actor_request"]).iterrows():
            lines.append(f"| {row.pde_id} | {int(row.dimension)} | {int(row.training_seed)} | {row.actor_request} | {row.resolved_actor} | {int(row.steps)} | {_number(row.mean_value_rmse)} | {_number(row.mean_mc_rmse_contribution)} | {_number(row.mean_inference_seconds)} |")
    lines.extend(["", "The headline uses the arithmetic mean of the four independently replicated 256-path RMSE values; the repetitions are not pooled into a 1,024-path estimator. Truth was loaded only after reference-free checkpoint selection.", ""])

    lines.extend(["## Full-History grid", ""])
    if fh.empty:
        lines.append("Not evaluated (compute).")
    else:
        lines.extend(["| PDE | d | method | n | M | reps | reused | RMSE | skill | GPU s | nonfinite |", "|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|"])
        for _, row in fh.sort_values(["pde_id", "dimension", "method", "n", "M"]).iterrows():
            nonfinite = (row.mean_nonfinite_state_count if np.isfinite(row.mean_nonfinite_state_count) else 0) + (row.mean_nonfinite_generator_count if np.isfinite(row.mean_nonfinite_generator_count) else 0)
            lines.append(f"| {row.pde_id} | {int(row.dimension)} | {METHOD_LABELS[row.method]} | {int(row.n)} | {int(row.M)} | {int(row.repetitions)} | {int(row.reused_repetitions)} | {_number(row.mean_value_rmse)} | {_number(row.mean_skill)} | {_number(row.mean_wall_time_seconds)} | {_number(nonfinite)} |")

    missing = len(fh_tasks()) - (0 if fh.empty else int(fh.repetitions.sum()))
    failures = sorted(key for key, item in outcomes.items() if item["verdict"] == "FAIL")
    lines.extend([
        "", "## Optional Stage C (exploratory)", "",
        "- **Residual-Jensen noise-replicate test: NOT EVALUATED (OPTIONAL COMPUTE).** "
        "The registered primary study, iteration audits, and d=100 seed-robustness runs were prioritized; "
        "no post-hoc exploratory update was started.",
        "- **DPI neural baseline: NOT EVALUATED (OPTIONAL DEPENDENCY).** "
        "The official implementation requires the external `picard` package, which is absent from this environment; "
        "no substitute implementation was introduced.",
    ])
    lines.extend([
        "", "## Failures and deviations", "",
        f"- Mechanical criterion failures: {', '.join(failures) if failures else 'none among evaluated criteria'}.",
        f"- Frozen Full-History rows not yet available: {missing}; affected criteria are marked `NOT EVALUATED (COMPUTE)` rather than removed.",
        "- Runtime hardware deviation: the requested GH200 was unavailable; all new timings are from an NVIDIA GeForce RTX 5070 Ti. Historical reused S05 rows were generated on the same GPU model.",
        "- The two final-suite DPI unit tests that import the external `picard` package were not runnable because that optional package is absent; the 18 relevant highdim/nonridge tests passed.",
        "- Reused S05 artifacts do not contain per-point predictions, so their MC/non-MC decomposition is unavailable; no proxy was invented.",
        "- B6 d=400 selected r_glob is below d/4. The frozen selector was applied unchanged and the resulting H6 failure is retained.",
    ])
    if curtailed.empty:
        lines.append("- No Full-History repetitions were curtailed by the pre-registered 30-minute rule.")
    else:
        for _, row in curtailed.sort_values(["pde_id", "d", "method", "n", "M", "repetition"]).iterrows():
            lines.append(
                f"- Compute-curtailed: {row.pde_id}, d={int(row.d)}, "
                f"{METHOD_LABELS[row.method]}, n={int(row.n)}, M={int(row.M)}, "
                f"rep={int(row.repetition)} after rep0 took {_number(row.rep0_seconds)} GPU seconds."
            )
    lines.extend(["", "## Provenance", ""])
    lines.extend([
        f"- Preregistration commit: `{provenance['preregistration_commit']}`.",
        f"- Calibration freeze commit: `{provenance['calibration_freeze_commit']}`.",
        f"- Analysis commit/worktree head: `{provenance['analysis_commit']}`.",
        f"- Imported nonridge commit: `d384f4e8`; byte-identity audit: {provenance['nonridge_byte_identity']}.",
        f"- Frozen final-suite protocol SHA-256: `{provenance['frozen_protocol_sha256']}`.",
        f"- Frozen review-gate SHA-256: `{provenance['review_gate_sha256']}`.",
        f"- GPU/software: {provenance['gpu_name']}; Python {platform.python_version()}, NumPy {np.__version__}, PyTorch {torch.__version__}, CUDA {torch.version.cuda}.",
        f"- Validation: {'PASS' if validation.get('all_passed') else 'FAIL or incomplete'}.",
        "", "## Figures", "",
        "- `results/highdim_nosym/figures/rmse_vs_gpu_seconds.png`",
        "- `results/highdim_nosym/figures/actor_iterations.png`",
        "- `results/highdim_nosym/figures/b6_calibration_heatmap.png`",
        "",
    ])
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text("\n".join(lines), encoding="utf-8")


def _nonridge_identity() -> str:
    source = _git("ls-tree", "-r", "d384f4e8", "invariant_region_mlp/experiments/nonridge")
    current_lines = []
    for line in source.splitlines():
        _mode, _kind, blob_and_path = line.split(maxsplit=2)
        blob, path = blob_and_path.split("\t", 1)
        current_blob = _git("hash-object", path)
        current_lines.append(blob == current_blob)
    return f"PASS ({sum(current_lines)}/{len(current_lines)} blobs identical)" if all(current_lines) else "FAIL"


def analyze() -> dict[str, Any]:
    ensure_directories()
    actor_rows, actor_summary, training = _actor_frames()
    fh_rows, fh_summary = _fh_frames()
    curtailed_rows = pd.DataFrame(_curtailed_rows(RAW_ROOT / "fh"))
    calibration = pd.read_csv(RESULTS_ROOT / "b6_calibration.csv")
    validation = json.loads((RESULTS_ROOT / "validation_summary.json").read_text(encoding="utf-8"))
    outcomes, threshold_rows = _criteria(actor_summary, fh_summary, training, calibration)
    frontier = _frontier(actor_summary, fh_summary, training)
    resources = _resource_rows(actor_summary, fh_summary, training)
    actor_export = actor_rows.copy()
    if not actor_export.empty:
        actor_export["row_family"] = "tangent_actor"
    fh_export = fh_rows.copy()
    if not fh_export.empty:
        fh_export["row_family"] = "full_history"
        fh_export["method"] = fh_export["method"].map(METHOD_LABELS)
    curtailed_export = curtailed_rows.copy()
    if not curtailed_export.empty:
        curtailed_export["row_family"] = "full_history_curtailed"
        curtailed_export["method"] = curtailed_export["method"].map(METHOD_LABELS)
    pd.concat((actor_export, fh_export, curtailed_export), ignore_index=True, sort=False).to_csv(
        RESULTS_ROOT / "main_rows.csv", index=False
    )
    iteration = (
        actor_summary[actor_summary.actor_request.str.match(r"a[0-4]")]
        if not actor_summary.empty else actor_summary
    )
    iteration.to_csv(RESULTS_ROOT / "iteration_diagnostics.csv", index=False)
    fh_summary.to_csv(RESULTS_ROOT / "fh_grid.csv", index=False)
    frontier.to_csv(RESULTS_ROOT / "frontier.csv", index=False)
    pd.DataFrame(threshold_rows).to_csv(RESULTS_ROOT / "time_to_threshold.csv", index=False)
    resources.to_csv(RESULTS_ROOT / "resource_accounting.csv", index=False)
    s05_ranks = {d: _s05_global_rank(d) for d in (100, 400)}
    frozen_protocol = _repository_root() / "invariant_region_mlp/results/tangent_actor_icml_final_suite/FROZEN_PROTOCOL.json"
    review_gate = _repository_root() / "invariant_region_mlp/results/tangent_actor_icml_final_suite/REVIEW_GATE_STATUS.json"
    calibration_freeze = json.loads(CALIBRATION_FREEZE_PATH.read_text(encoding="utf-8"))
    device_name = calibration_freeze["environment"].get("gpu_name", "NA")
    provenance = {
        "preregistration_commit": _git(
            "log", "--diff-filter=A", "-1", "--format=%H", "--",
            "invariant_region_mlp/experiments/highdim_nosym/PREREGISTRATION.md",
        ),
        "calibration_freeze_commit": _git(
            "log", "--diff-filter=A", "-1", "--format=%H", "--",
            "invariant_region_mlp/results/highdim_nosym/B6_CALIBRATION_FREEZE.json",
        ),
        "analysis_commit": _git("rev-parse", "HEAD"),
        "nonridge_byte_identity": _nonridge_identity(),
        "frozen_protocol_sha256": _sha256(frozen_protocol),
        "review_gate_sha256": _sha256(review_gate),
        "gpu_name": device_name,
    }
    payload = {
        "study_revision": REVISION,
        "criteria": outcomes,
        "s05_r_glob": {str(key): value for key, value in s05_ranks.items()},
        "completed_actor_replicates": len(actor_rows),
        "completed_fh_rows": len(fh_rows),
        "curtailed_fh_rows": len(curtailed_rows),
        "expected_fh_rows": len(fh_tasks()),
        "provenance": provenance,
        "validation_all_passed": validation.get("all_passed", False),
    }
    write_json_atomic(RESULTS_ROOT / "analysis_summary.json", payload)
    _plots(actor_summary, fh_summary, frontier, calibration)
    _report(
        outcomes, actor_summary, fh_summary, training, calibration, validation,
        s05_ranks, provenance, curtailed_rows,
    )
    print(json.dumps(payload, indent=2, sort_keys=True, allow_nan=True), flush=True)
    return payload


if __name__ == "__main__":
    analyze()
