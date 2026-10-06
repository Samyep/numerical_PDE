"""Aggregate VB repetitions, tune controls, draw work curves, and write the report."""

from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import sys
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[1]
RESULT_ROOT = PROJECT_ROOT / "results" / "active_vb_high_budget"
FIGURE_ROOT = RESULT_ROOT / "figures"
REPORT_PATH = PROJECT_ROOT / "docs" / "VB_HIGH_BUDGET_REPORT.md"
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from vb_equation import PUBLISHED_SIGMA, SCASML_AUDITED_COMMIT  # noqa: E402
from vb_mlp_methods import core_methods, parse_method  # noqa: E402


CORE_NAMES = [method.name for method in core_methods()]
DISPLAY = {
    "raw": "Raw",
    "sample_box": "Certified box",
    "z_only_box": "Certified z-only box",
    "sample_ball": "Certified ball",
    "batch_box": "Batch-IR box",
    "z_zero": "z=0",
    "f_zero": "f=0",
    "best_shrink": "Validation-tuned shrink",
    "best_tighter": "Validation-tuned tighter box",
    "best_looser": "Validation-tuned looser box",
}
COLORS = {
    "raw": "#222222",
    "sample_box": "#0072B2",
    "z_only_box": "#56B4E9",
    "sample_ball": "#009E73",
    "batch_box": "#CC79A7",
    "z_zero": "#D55E00",
    "f_zero": "#E69F00",
    "best_shrink": "#7A3E00",
    "best_tighter": "#6A3D9A",
    "best_looser": "#A6761D",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=True), encoding="utf-8")
    temporary.replace(path)


def _nested(payload: dict[str, Any], *keys: str, default: Any = np.nan) -> Any:
    current: Any = payload
    for key in keys:
        if current is None or key not in current:
            return default
        current = current[key]
    return current


def repetition_row(path: Path) -> tuple[dict[str, Any], np.ndarray]:
    with np.load(path, allow_pickle=False) as data:
        metadata = json.loads(str(data["metadata_json"]))
        correction = np.array(data["nonlinear_u_correction"], copy=True)
    work = metadata["work"]
    metrics = metadata["metrics"]
    generator = work.get("generator")
    alpha = work.get("batch_alpha")
    subset = metadata.get("point_subset", "all")
    selection_metrics = metrics["all"] if subset == "validation" else metrics["validation"]
    row = {
        "path": str(path.relative_to(RESULT_ROOT)).replace("\\", "/"),
        "d": metadata["dimension"],
        "n": metadata["n"],
        "M": metadata["M"],
        "method": metadata["method"]["name"],
        "method_transform": metadata["method"]["transform"],
        "method_factor": metadata["method"]["factor"],
        "repetition": metadata["repetition"],
        "point_subset": subset,
        "wall_clock_seconds": metadata["wall_clock_seconds"],
        "selection_value_relative_l2": selection_metrics["value_relative_l2"],
    }
    for split in ("all", "validation", "test"):
        for metric, value in metrics[split].items():
            row[f"{split}_{metric}"] = value
    for key in (
        "terminal_g_evals",
        "terminal_samples",
        "f_evals",
        "recursively_evaluated_states",
        "transition_samples",
        "total_stochastic_samples",
        "standard_normal_variates",
        "time_uniform_variates",
        "correction_opportunities",
        "pre_box_violation_rate",
        "pre_ball_violation_rate",
        "method_geometry_violation_rate",
        "projection_activation_rate",
        "mean_box_overshoot_energy",
        "mean_method_overshoot_energy",
        "nonfinite_generator_values",
    ):
        row[key] = work.get(key, np.nan)
    for key in ("mse", "bias", "absolute_bias", "mae", "mean_estimate", "mean_truth"):
        row[f"generator_{key}"] = np.nan if generator is None else generator[key]
    for key in ("mean", "median", "p10", "p90", "activation_rate"):
        row[f"batch_alpha_{key}"] = np.nan if alpha is None else alpha[key]
    return row, correction


def collect_stage(stage: str) -> tuple[pd.DataFrame, dict[tuple[Any, ...], list[np.ndarray]]]:
    root = RESULT_ROOT / "repetitions" / stage
    paths = sorted(root.glob("*.npz"))
    if not paths:
        raise FileNotFoundError(f"no repetition artifacts under {root}")
    rows: list[dict[str, Any]] = []
    corrections: dict[tuple[Any, ...], list[np.ndarray]] = defaultdict(list)
    for path in paths:
        row, correction = repetition_row(path)
        rows.append(row)
        key = (row["d"], row["n"], row["M"], row["method"])
        corrections[key].append(correction)
    return pd.DataFrame(rows), corrections


def aggregate_repetitions(
    repetitions: pd.DataFrame,
    corrections: dict[tuple[Any, ...], list[np.ndarray]],
) -> pd.DataFrame:
    keys = ["d", "n", "M", "method", "method_transform", "method_factor", "point_subset"]
    numeric = [column for column in repetitions.columns if column not in keys + ["path"]]
    aggregations: dict[str, list[str] | str] = {}
    for column in numeric:
        if column == "repetition":
            continue
        if column in {
            "terminal_g_evals",
            "terminal_samples",
            "f_evals",
            "recursively_evaluated_states",
            "transition_samples",
            "total_stochastic_samples",
            "standard_normal_variates",
            "time_uniform_variates",
            "correction_opportunities",
        }:
            aggregations[column] = "mean"
        else:
            aggregations[column] = ["mean", "std"]
    grouped = repetitions.groupby(keys, dropna=False).agg(aggregations)
    grouped.columns = [
        column if isinstance(stat, str) and stat == "" else f"{column}_{stat}"
        for column, stat in grouped.columns
    ]
    grouped = grouped.reset_index()
    counts = repetitions.groupby(keys, dropna=False).size().rename("repetitions").reset_index()
    grouped = grouped.merge(counts, on=keys, how="left")

    variance_rows = []
    for key, arrays in corrections.items():
        if len(arrays) >= 2:
            stacked = np.stack(arrays, axis=0)
            variance = float(np.mean(np.var(stacked, axis=0, ddof=1)))
        else:
            variance = float("nan")
        variance_rows.append(
            {"d": key[0], "n": key[1], "M": key[2], "method": key[3], "nonlinear_correction_variance": variance}
        )
    grouped = grouped.merge(pd.DataFrame(variance_rows), on=["d", "n", "M", "method"], how="left")

    paired = []
    for (d, n, M), frame in repetitions.groupby(["d", "n", "M"]):
        pivot = frame.pivot_table(index="repetition", columns="method", values="test_value_relative_l2")
        complete = pivot.dropna(axis=0, how="any")
        fractional_wins = {method: 0.0 for method in pivot.columns}
        if len(complete):
            for _, values in complete.iterrows():
                best = float(values.min())
                tied = [method for method, value in values.items() if abs(float(value) - best) <= 1e-14]
                for method in tied:
                    fractional_wins[method] += 1.0 / len(tied)
        for method in pivot.columns:
            item: dict[str, Any] = {
                "d": d,
                "n": n,
                "M": M,
                "method": method,
                "method_win_fraction": (
                    fractional_wins[method] / len(complete) if len(complete) else np.nan
                ),
            }
            for reference in ("raw", "sample_box", "z_zero", "f_zero"):
                column = f"paired_win_fraction_vs_{reference}"
                if reference in pivot:
                    common = pivot[[method, reference]].dropna()
                    item[column] = (
                        float(np.mean(common[method] < common[reference])) if len(common) else np.nan
                    )
                else:
                    item[column] = np.nan
            paired.append(item)
    if paired:
        grouped = grouped.merge(pd.DataFrame(paired), on=["d", "n", "M", "method"], how="left")
    else:
        for column in (
            "method_win_fraction",
            "paired_win_fraction_vs_raw",
            "paired_win_fraction_vs_sample_box",
            "paired_win_fraction_vs_z_zero",
            "paired_win_fraction_vs_f_zero",
        ):
            grouped[column] = np.nan
    return grouped


def choose_controls(aggregate: pd.DataFrame) -> dict[str, Any]:
    choices: dict[str, Any] = {"schema_version": 1, "created_utc": utc_now(), "choices": {}}
    for (d, n, M), frame in aggregate.groupby(["d", "n", "M"]):
        errors = dict(zip(frame["method"], frame["selection_value_relative_l2_mean"]))
        shrink_candidates: dict[float, tuple[str, float]] = {}
        tighter_candidates: dict[float, tuple[str, float]] = {}
        looser_candidates: dict[float, tuple[str, float]] = {}
        if "z_zero" in errors:
            shrink_candidates[0.0] = ("z_zero", errors["z_zero"])
            tighter_candidates[0.0] = ("z_zero", errors["z_zero"])
        if "raw" in errors:
            shrink_candidates[1.0] = ("raw", errors["raw"])
        if "sample_box" in errors:
            # The exact tight certificate is the a=1 endpoint of the valid
            # relaxed-box family, but is not part of the invalid tighter set.
            looser_candidates[1.0] = ("sample_box", errors["sample_box"])
        for method, error in errors.items():
            if method.startswith("shrink_c"):
                shrink_candidates[float(method.removeprefix("shrink_c"))] = (method, error)
            if method.startswith("box_factor_a"):
                factor = float(method.removeprefix("box_factor_a"))
                if factor < 1.0:
                    tighter_candidates[factor] = (method, error)
                elif factor > 1.0:
                    looser_candidates[factor] = (method, error)
        key = f"d{int(d)}_n{int(n)}_M{int(M)}"
        entry: dict[str, Any] = {}
        if shrink_candidates:
            c, (method, error) = min(shrink_candidates.items(), key=lambda item: item[1][1])
            entry["best_shrink"] = {"c": c, "method": method, "validation_relative_l2": error}
        if tighter_candidates:
            a, (method, error) = min(tighter_candidates.items(), key=lambda item: item[1][1])
            entry["best_tighter"] = {"a": a, "method": method, "validation_relative_l2": error}
        if looser_candidates:
            a, (method, error) = min(looser_candidates.items(), key=lambda item: item[1][1])
            entry["best_looser"] = {"a": a, "method": method, "validation_relative_l2": error}
        entry["candidate_errors"] = {method: float(error) for method, error in sorted(errors.items())}
        choices["choices"][key] = entry
    return choices


def make_main_plan(
    choices: dict[str, Any],
    aggregate: pd.DataFrame,
    *,
    repetitions: int,
    headline_repetitions: int,
) -> dict[str, Any]:
    plan = []
    by_dimension: dict[int, list[tuple[int, int, float, float]]] = defaultdict(list)
    unique_configs = aggregate[
        ["d", "n", "M", "f_evals_mean", "total_stochastic_samples_mean"]
    ].drop_duplicates()
    for row in unique_configs.itertuples(index=False):
        by_dimension[int(row.d)].append(
            (int(row.n), int(row.M), float(row.total_stochastic_samples_mean), float(row.f_evals_mean))
        )
    headline: set[tuple[int, int, int]] = set()
    for d, configs in by_dimension.items():
        configs.sort(key=lambda item: item[2])
        if configs:
            indices = sorted(set([0, len(configs) // 2, len(configs) - 1]))
            headline.update((d, configs[index][0], configs[index][1]) for index in indices)
            deepest = max(configs, key=lambda item: item[3])
            headline.add((d, deepest[0], deepest[1]))
    for d, n, M in sorted(
        {(int(row.d), int(row.n), int(row.M)) for row in unique_configs.itertuples(index=False)}
    ):
        key = f"d{d}_n{n}_M{M}"
        entry = choices["choices"].get(key, {})
        methods = list(CORE_NAMES)
        for family in ("best_shrink", "best_tighter", "best_looser"):
            method = _nested(entry, family, "method", default=None)
            if method and method not in methods:
                methods.append(method)
        plan.append(
            {
                "d": d,
                "n": n,
                "M": M,
                "repetitions": headline_repetitions if (d, n, M) in headline else repetitions,
                "methods": methods,
                "point_subset": "all",
            }
        )
    return {
        "schema_version": 1,
        "created_utc": utc_now(),
        "selection_source": "held-out validation split only",
        "headline_settings": [list(item) for item in sorted(headline)],
        "plan": plan,
    }


def _choice_method(choices: dict[str, Any], d: int, n: int, M: int, family: str) -> str | None:
    return _nested(choices, "choices", f"d{d}_n{n}_M{M}", family, "method", default=None)


def plotting_rows(aggregate: pd.DataFrame, choices: dict[str, Any]) -> pd.DataFrame:
    pieces = []
    for _, row in aggregate.iterrows():
        method = row["method"]
        if method in CORE_NAMES:
            copy = row.copy()
            copy["plot_method"] = method
            pieces.append(copy)
        d, n, M = int(row["d"]), int(row["n"]), int(row["M"])
        for family in ("best_shrink", "best_tighter", "best_looser"):
            if method == _choice_method(choices, d, n, M, family):
                copy = row.copy()
                copy["plot_method"] = family
                pieces.append(copy)
    if not pieces:
        return pd.DataFrame()
    result = pd.DataFrame(pieces)
    return result.drop_duplicates(
        subset=["d", "n", "M", "plot_method"], keep="first"
    ).reset_index(drop=True)


def _panel_plot(
    data: pd.DataFrame,
    *,
    x: str,
    y: str,
    ylabel: str,
    output: Path,
    log_x: bool = True,
    log_y: bool | str = True,
) -> None:
    dimensions = sorted(data["d"].unique())
    fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.4), squeeze=False)
    for axis, d in zip(axes.ravel(), dimensions):
        frame = data[data["d"] == d]
        for method, group in frame.groupby("plot_method"):
            group = group.sort_values(x)
            axis.plot(
                group[x],
                group[y],
                marker="o",
                markersize=3.5,
                linewidth=1.25,
                label=DISPLAY.get(method, method),
                color=COLORS.get(method),
            )
        axis.set_title(f"d={d}")
        axis.grid(True, which="both", alpha=0.24)
        if log_x:
            if np.any(frame[x] == 0):
                axis.set_xscale("symlog", linthresh=1.0)
            else:
                axis.set_xscale("log")
        if log_y == "symlog":
            axis.set_yscale("symlog", linthresh=1.0)
        elif log_y and np.all(frame[y].dropna() > 0):
            axis.set_yscale("log")
        axis.set_xlabel(x.replace("_", " "))
        axis.set_ylabel(ylabel)
    for axis in axes.ravel()[len(dimensions) :]:
        axis.set_visible(False)
    handles, labels = axes.ravel()[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def make_figures(aggregate: pd.DataFrame, choices: dict[str, Any]) -> list[str]:
    data = plotting_rows(aggregate, choices)
    if data.empty:
        return []
    outputs = []
    specifications = [
        ("f_evals_mean", "test_value_relative_l2_mean", "Held-out value relative L2", "error_vs_f_calls.pdf", True, True),
        ("wall_clock_seconds_mean", "test_value_relative_l2_mean", "Held-out value relative L2", "error_vs_wallclock.pdf", True, True),
        ("total_stochastic_samples_mean", "test_value_relative_l2_mean", "Held-out value relative L2", "method_vs_budget_by_dimension.pdf", True, True),
        ("f_evals_mean", "generator_bias_mean", "Generator bias", "generator_bias_vs_budget.pdf", True, "symlog"),
        ("f_evals_mean", "pre_box_violation_rate_mean", "Pre-correction box violation rate", "violation_rate_vs_budget.pdf", True, False),
    ]
    for x, y, ylabel, filename, log_x, log_y in specifications:
        if y not in data or data[y].notna().sum() == 0:
            continue
        path = FIGURE_ROOT / filename
        _panel_plot(data, x=x, y=y, ylabel=ylabel, output=path, log_x=log_x, log_y=log_y)
        outputs.append(str(path.relative_to(RESULT_ROOT)).replace("\\", "/"))

    # Compact categorical crossover summary: winner at every actual-work point.
    winner_rows = []
    for (d, n, M), frame in data.groupby(["d", "n", "M"]):
        finite = frame[np.isfinite(frame["test_value_relative_l2_mean"])]
        if finite.empty:
            continue
        winner = finite.loc[finite["test_value_relative_l2_mean"].idxmin()]
        winner_rows.append(
            {
                "d": d,
                "n": n,
                "M": M,
                "work": winner["total_stochastic_samples_mean"],
                "winner": winner["plot_method"],
            }
        )
    if winner_rows:
        winners = pd.DataFrame(winner_rows).sort_values(["d", "work"])
        methods = list(DISPLAY)
        code = {method: index for index, method in enumerate(methods)}
        fig, axes = plt.subplots(2, 2, figsize=(10.5, 5.8), squeeze=False)
        for axis, (d, frame) in zip(axes.ravel(), winners.groupby("d")):
            values = [code[item] for item in frame["winner"]]
            colors = [COLORS[item] for item in frame["winner"]]
            axis.scatter(np.arange(len(frame)), np.zeros(len(frame)), c=colors, s=110, marker="s")
            for index, row in enumerate(frame.itertuples(index=False)):
                axis.text(index, 0.0, f"{row.n},{row.M}", ha="center", va="center", fontsize=6, color="white")
            axis.set_title(f"d={d}: winning method by increasing stochastic work")
            axis.set_yticks([])
            axis.set_xlabel("budget rank (label is n,M)")
            axis.set_xlim(-0.7, len(frame) - 0.3)
        for axis in axes.ravel()[winners["d"].nunique() :]:
            axis.set_visible(False)
        legend_handles = [
            plt.Line2D([0], [0], marker="s", linestyle="", color=COLORS[name], label=DISPLAY[name])
            for name in sorted(set(winners["winner"]))
        ]
        fig.legend(handles=legend_handles, loc="upper center", ncol=4, fontsize=8, frameon=False)
        fig.tight_layout(rect=(0, 0, 1, 0.9))
        path = FIGURE_ROOT / "crossover_summary.pdf"
        fig.savefig(path, bbox_inches="tight")
        plt.close(fig)
        outputs.append(str(path.relative_to(RESULT_ROOT)).replace("\\", "/"))
    return outputs


def _fmt(value: Any, digits: int = 4) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not math.isfinite(number):
        return "NA"
    return f"{number:.{digits}g}"


def _verdict(aggregate: pd.DataFrame, choices: dict[str, Any], forced: str) -> tuple[str, str]:
    data = plotting_rows(aggregate, choices)
    cert = {"sample_box", "z_only_box", "sample_ball", "batch_box"}
    control = {"z_zero", "f_zero", "best_shrink", "best_tighter"}
    comparisons = []
    for _, frame in data.groupby(["d", "n", "M"]):
        ce = frame[frame["plot_method"].isin(cert)]["test_value_relative_l2_mean"]
        co = frame[frame["plot_method"].isin(control)]["test_value_relative_l2_mean"]
        if len(ce) and len(co):
            comparisons.append((float(ce.min()), float(co.min())))
    certified_wins = sum(left < right for left, right in comparisons)
    control_wins = sum(right < left for left, right in comparisons)
    if forced == "auto":
        if comparisons and certified_wins >= max(2, math.ceil(0.25 * len(comparisons))):
            forced = "A"
        elif comparisons and control_wins == len(comparisons):
            forced = "B"
        else:
            forced = "C"
    justification = (
        f"Across {len(comparisons)} comparable dimension/budget cells, the best certified method "
        f"beat all tuned suppression controls in {certified_wins}; a suppression control won in {control_wins}."
    )
    text = {
        "A": "Certified geometry has an empirical regime of advantage beyond generic shrinkage.",
        "B": "Generic gradient suppression explains the observed finite-budget gains better than certified geometry.",
        "C": "Evidence is mixed / insufficient; specify exactly which additional experiment is needed.",
    }[forced]
    return forced, justification + "\n\nVERDICT " + forced + ":\n" + text


def write_report(
    aggregate: pd.DataFrame,
    choices: dict[str, Any],
    signal: dict[str, Any] | None,
    stage: str,
    verdict: str,
) -> None:
    dimensions = sorted(int(value) for value in aggregate["d"].unique())
    configs = sorted({(int(row.n), int(row.M)) for row in aggregate.itertuples(index=False)})
    methods = sorted(aggregate["method"].unique())
    lines = [
        "# Active-gradient viscous-Burgers high-budget MLP study",
        "",
        f"Generated from the `{stage}` repetition artifacts on {utc_now()}. No manuscript source was modified.",
        "",
        "## Benchmark and exact solution",
        "",
        "The main experiment uses the published SCaSML equation with `T=0.5`, "
        "`sigma=sqrt(2)`, coordinatewise drift `mu=-(1/d+sigma^2/2)`, generator "
        "`f(u,z)=sigma*u*sum_i z_i`, and state convention `z=sigma*grad(u)`. The exact "
        "solution is `u*(t,x)=sigmoid(t+sum_i x_i)` and `z_i*=sigma*u*(1-u)`.",
        "",
        "Writing `q=u(1-u)`, the analytic residual is",
        "",
        "`q + d*mu*q + (sigma^2/2)d*q(1-2u) + d*sigma^2*u*q = 0`.",
        "",
        "The unit test evaluates this identity in dimensions 2, 20, and 80 and requires maximum "
        "absolute residual below `2e-13`.",
        "",
        "## Exact certificates and corrections",
        "",
        "Because `0<u<1` and `0<=u(1-u)<=1/4`, every exact state lies in "
        "`[0,1] x [0,sigma/4]^d` and in the z-ball of radius `sigma*sqrt(d)/4`. "
        "Samplewise box, z-only box, and ball corrections are the corresponding Euclidean "
        "projections. Batch-IR first replaces negative z coordinates by zero, then applies one "
        "common sibling factor `alpha=min(1,(sigma/4)/max_{i,j}(z_ij)_+)`; u is clipped "
        "samplewise. Thus Batch-IR is box-feasible and identity on feasible batches, but it is "
        "explicitly not called a coordinatewise box projection.",
        "",
        "## Published/public-repository discrepancy and solver provenance",
        "",
        f"The published report specifies `sigma=sqrt(2)`. Public SCaSML commit "
        f"`{SCASML_AUDITED_COMMIT}` instead returns `sigma=0.25` for this equation. The public "
        "full-history terminal estimator also divides a standard normal by `T-t`; this study uses "
        "the corrected EBL normalization `xi/sqrt(T-t)`. All scientific runs use float64, no "
        "SCaSML surrogate clipping, and Beta(1/2,1) time importance sampling. The algebraically "
        "zero level-0 generator summand is elided and recorded as such.",
        "",
        "## Nonlinear-signal gate",
        "",
    ]
    if signal:
        lines.extend(
            [
                "| d | mean u* | median ||z*|| | mean |f*| | max |f*| | f=0 relative L2 |",
                "|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for row in signal["rows"]:
            lines.append(
                f"| {row['d']} | {_fmt(row['u']['mean'])} | {_fmt(row['z_l2']['median'])} | "
                f"{_fmt(row['abs_f']['mean'])} | {_fmt(row['abs_f']['max'])} | {_fmt(row['fzero_relative_l2'])} |"
            )
    else:
        lines.append("Signal-gate artifact was not available for this analysis.")
    lines.extend(
        [
            "",
            "## Grid actually completed",
            "",
            f"Dimensions: `{dimensions}`. Configurations `(n,M)`: `{configs}`. Methods present: `{methods}`. "
            f"The aggregate contains `{int(aggregate['repetitions'].sum())}` per-method repetitions.",
            "",
            "Fixed points are sampled on `t in [0,0.5]`, `x in [-0.5,0.5]^d`, with boundary "
            "points forcing one random coordinate to `+/-0.5`. Hyperparameters are selected only "
            "on the stratified 20% validation split; reported headline errors use the 80% test split.",
            "",
            "## Work accounting",
            "",
            "Every repetition records terminal g evaluations, nonlinear f evaluations, recursively "
            "evaluated states, transition and terminal stochastic samples, scalar normal draws, and "
            "wall time. `z=0` traverses the complete recursive tree and evaluates the zeroed generator; "
            "`f=0` deletes recursion and therefore has zero f calls. Total stochastic samples—not f calls "
            "alone—are used for comparisons involving `f=0`.",
            "",
            "## Validation-tuned controls",
            "",
            "The full machine-readable choices, candidate losses, and aliases (`c=0` is `z=0`, "
            "`c=1` is Raw, and tighter-box `a=0` is value-equivalent to `z=0`) are in "
            "`results/active_vb_high_budget/tuning_choices.json`.",
            "",
            "## Complete aggregate result table",
            "",
            "The exact full table is `results/active_vb_high_budget/work_summary.csv`. The compact "
            "table below gives held-out value error, gradient error, work, wall time, generator bias, "
            "and pre-correction box violations for every completed cell.",
            "",
            "| d | n | M | method | reps | value rel L2 | grad rel L2 | f calls | samples | sec | gen bias | box viol. |",
            "|---:|---:|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in aggregate.sort_values(["d", "f_evals_mean", "n", "M", "method"]).itertuples(index=False):
        lines.append(
            f"| {int(row.d)} | {int(row.n)} | {int(row.M)} | {row.method} | {int(row.repetitions)} | "
            f"{_fmt(row.test_value_relative_l2_mean)} | {_fmt(row.test_gradient_relative_l2_mean)} | "
            f"{_fmt(row.f_evals_mean, 6)} | {_fmt(row.total_stochastic_samples_mean, 6)} | "
            f"{_fmt(row.wall_clock_seconds_mean)} | {_fmt(row.generator_bias_mean)} | "
            f"{_fmt(row.pre_box_violation_rate_mean)} |"
        )
    verdict_code, verdict_text = _verdict(aggregate, choices, verdict)
    lines.extend(
        [
            "",
            "## Generator, projection, and crossover diagnostics",
            "",
            "Generator MSE/bias/MAE, nonlinear-correction variance across paired repetitions, box and "
            "ball violation rates, overshoot energies, activation rates, Batch-alpha quantiles, and "
            "paired win fractions are columns in the CSV/JSON summaries. The six PDFs in the result "
            "figure directory show error against f calls, wall time, and total stochastic work, plus "
            "generator bias, violation rate, and the winner sequence by dimension.",
            "",
            "## Life-or-death conclusion and paper recommendation",
            "",
            verdict_text.split("\n\nVERDICT")[0],
            "The paper should cite this outcome directly and distinguish correctness-preserving "
            "certification from validation-tuned accuracy regularization. No manuscript edit has been "
            "made in this branch.",
            "",
            f"VERDICT {verdict_code}:",
            verdict_text.split("\n")[-1],
        ]
    )
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", default="main")
    parser.add_argument(
        "--selection-stage",
        help="comma-separated stage(s) used to select c/a; defaults to --stage",
    )
    parser.add_argument("--main-repetitions", type=int, default=5)
    parser.add_argument("--headline-repetitions", type=int, default=10)
    parser.add_argument("--verdict", choices=("auto", "A", "B", "C"), default="auto")
    parser.add_argument("--plan-only", action="store_true")
    args = parser.parse_args()

    primary_stage_names = args.stage.split(",")
    primary_frames = []
    corrections: dict[tuple[Any, ...], list[np.ndarray]] = defaultdict(list)
    for primary_stage in primary_stage_names:
        stage_repetitions, stage_corrections = collect_stage(primary_stage)
        primary_frames.append(stage_repetitions)
        for key, values in stage_corrections.items():
            corrections[key].extend(values)
    repetitions = pd.concat(primary_frames, ignore_index=True)
    aggregate = aggregate_repetitions(repetitions, corrections)
    selection_stage_names = (args.selection_stage or args.stage).split(",")
    if selection_stage_names == primary_stage_names:
        selection_aggregate = aggregate
    else:
        selection_frames = []
        combined_corrections: dict[tuple[Any, ...], list[np.ndarray]] = defaultdict(list)
        for selection_stage in selection_stage_names:
            selection_repetitions, selection_corrections = collect_stage(selection_stage)
            selection_frames.append(selection_repetitions)
            for key, values in selection_corrections.items():
                combined_corrections[key].extend(values)
        selection_aggregate = aggregate_repetitions(
            pd.concat(selection_frames, ignore_index=True), combined_corrections
        )
    choices = choose_controls(selection_aggregate)
    atomic_json(RESULT_ROOT / "tuning_choices.json", choices)
    plan = make_main_plan(
        choices,
        selection_aggregate,
        repetitions=args.main_repetitions,
        headline_repetitions=args.headline_repetitions,
    )
    atomic_json(RESULT_ROOT / "main_plan.json", plan)
    RESULT_ROOT.mkdir(parents=True, exist_ok=True)
    safe_stage = "_".join(primary_stage_names).replace("/", "_").replace("\\", "_")
    aggregate.to_csv(RESULT_ROOT / f"work_summary_{safe_stage}.csv", index=False)
    repetitions.to_csv(RESULT_ROOT / f"repetition_metrics_{safe_stage}.csv", index=False)
    selection_aggregate.to_csv(RESULT_ROOT / "selection_candidate_summary.csv", index=False)
    if args.plan_only:
        print(RESULT_ROOT / "main_plan.json")
        return

    aggregate.to_csv(RESULT_ROOT / "work_summary.csv", index=False)
    repetitions.to_csv(RESULT_ROOT / "repetition_metrics.csv", index=False)
    figures = make_figures(aggregate, choices)
    signal_path = RESULT_ROOT / "signal_strength.json"
    signal = json.loads(signal_path.read_text(encoding="utf-8")) if signal_path.exists() else None
    summary = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "stage": primary_stage_names,
        "selection_stage": selection_stage_names,
        "rows": json.loads(aggregate.to_json(orient="records")),
        "figures": figures,
        "choices_path": "tuning_choices.json",
        "report_path": str(REPORT_PATH.relative_to(PROJECT_ROOT)).replace("\\", "/"),
    }
    atomic_json(RESULT_ROOT / "full_summary.json", summary)
    write_report(aggregate, choices, signal, ",".join(primary_stage_names), args.verdict)
    print(f"wrote {RESULT_ROOT / 'full_summary.json'}")
    print(f"wrote {REPORT_PATH}")


if __name__ == "__main__":
    main()
