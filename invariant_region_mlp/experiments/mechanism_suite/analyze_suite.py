"""Aggregate mechanism-suite repetitions, evaluate criteria, and make report."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from .mechanism_mlp import U_DOSE_SCALE, load_repetition
from .run_suite import (
    GATES_PATH,
    RESULTS_ROOT,
    e1_tasks,
    e2_tasks,
    e3_tasks,
    e5_tasks,
    e6_tasks,
    task_path,
)


FIGURE_ROOT = RESULTS_ROOT / "figures"
DOC_PATH = RESULTS_ROOT.parent.parent / "docs" / "MECHANISM_SUITE_REPORT.md"
TUNING_PATH = RESULTS_ROOT / "e3_tuning_choices.json"
SUMMARY_PATH = RESULTS_ROOT / "analysis_summary.json"

ACCEPTANCE_CRITERIA = (
    "**C1 supported** for a PDE if `D >= 0.5` in at least 3 of 4 headline "
    "configurations at each admitted d, and the E2 dose curve for raw input "
    "increases monotonically in s.",
    "**C2 supported** if `box/oracle <= 0.5 * raw/oracle` (equivalently box "
    "removes at least half of the excess error) in at least 3 of 4 configurations.",
    "**C3 supported** if box wins at least 7/10 paired reps against each of: "
    "best tuned ->0 control, best tuned illegal set, centre, best tuned "
    "shrink_centre, in at least 75% of E3 cells.",
    "Any PDE where centre or shrink_centre beats box in more than 25% of cells "
    "must be reported as \"certificate-prior dominated\" for that regime.",
    "Negative controls must behave as predicted; if not, report and explain.",
)


def _atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    frame.to_csv(temporary, index=False)
    temporary.replace(path)


def _json_scalar(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_scalar(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_scalar(item) for item in value]
    if isinstance(value, np.ndarray):
        return _json_scalar(value.tolist())
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


def _atomic_json(payload: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(_json_scalar(payload), indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _finite_skill(prediction: np.ndarray, truth: np.ndarray, mask: np.ndarray) -> float:
    pred = np.asarray(prediction, dtype=np.float64)[mask]
    ref = np.asarray(truth, dtype=np.float64)[mask]
    if not np.all(np.isfinite(pred)):
        return float("inf")
    return float(np.sqrt(np.mean((pred - ref) ** 2)) / np.std(ref, ddof=0))


def _flatten_artifact(path: Path, stage: str) -> dict[str, Any]:
    path = path.resolve()
    loaded = load_repetition(path)
    metadata = loaded["metadata"]
    method = metadata["method"]
    work = metadata["work"]
    test = metadata["metrics"]["test"]
    validation = metadata["metrics"]["validation"]
    noise = metadata["extra_diagnostics"]["child_u_z_noise"]
    generator = work["generator"]
    batch_alpha = work["batch_alpha"]
    nonfinite = int(test["nonfinite_state_count"]) + int(work["nonfinite_generator_values"])
    reused_from_stage: str | None = None
    if stage == "e3" and method["name"] in {
        "raw",
        "box",
        "segment",
        "oracle_state",
        "oracle_z",
        "oracle_u",
    }:
        reused_from_stage = "e1"
    elif stage == "e6" and method["name"] in {"raw", "box", "segment", "oracle_state"}:
        reused_from_stage = "e1"
    elif stage == "e6" and method["name"] == "centre":
        dimension = int(metadata["dimension"])
        if dimension in {20, 50}:
            reused_from_stage = "e3"
        elif int(metadata["repetition"]) < 3:
            reused_from_stage = "e0"
    elif (
        stage == "e5"
        and int(metadata["repetition"]) < 3
        and int(metadata["dimension"]) >= 50
        and method["name"] in {"raw", "box", "oracle_state", "centre"}
    ):
        reused_from_stage = "e0"
    return {
        "artifact": path.relative_to(RESULTS_ROOT).as_posix(),
        "stage": stage,
        "reused_from_stage": reused_from_stage,
        "pde": metadata["pde_id"],
        "equation_name": metadata["equation_name"],
        "family": metadata["family"],
        "d": int(metadata["dimension"]),
        "n": int(metadata["n"]),
        "M": int(metadata["M"]),
        "method": method["name"],
        "transform": method["transform"],
        "factor": float(method["factor"]),
        "rep": int(metadata["repetition"]),
        "chunk_size": int(metadata["chunk_size"]),
        "dtype": metadata["dtype"],
        "test_skill": float(test["skill"]),
        "test_value_rmse": float(test["value_rmse"]),
        "test_value_bias": float(test["value_bias"]),
        "test_value_mae": float(test["value_mae"]),
        "test_value_relative_l2": float(test["value_relative_l2"]),
        "test_gradient_relative_l2": float(test["gradient_relative_l2"]),
        "validation_skill": float(validation["skill"]),
        "validation_value_rmse": float(validation["value_rmse"]),
        "f_calls": int(work["f_evals"]),
        "recursively_evaluated_states": int(work["recursively_evaluated_states"]),
        "stochastic_samples": int(work["total_stochastic_samples"]),
        "standard_normal_variates": int(work["standard_normal_variates"]),
        "dose_normal_variates": int(
            metadata["extra_diagnostics"].get("dose_normal_variates", 0)
        ),
        "wall_time_seconds": float(metadata["wall_clock_seconds"]),
        "pre_correction_violation_rate": float(work["pre_box_violation_rate"]),
        "pre_ball_violation_rate": float(work["pre_ball_violation_rate"]),
        "activation_rate": float(work["projection_activation_rate"]),
        "mean_box_overshoot_energy": float(work["mean_box_overshoot_energy"]),
        "mean_method_overshoot_energy": float(work["mean_method_overshoot_energy"]),
        "batch_alpha_mean": (
            float(batch_alpha["mean"]) if batch_alpha is not None else np.nan
        ),
        "batch_alpha_p10": (
            float(batch_alpha["p10"]) if batch_alpha is not None else np.nan
        ),
        "batch_alpha_p90": (
            float(batch_alpha["p90"]) if batch_alpha is not None else np.nan
        ),
        "nonfinite_count": nonfinite,
        "nonfinite_state_count": int(test["nonfinite_state_count"]),
        "nonfinite_generator_count": int(work["nonfinite_generator_values"]),
        "generator_bias": float(generator["bias"]) if generator is not None else np.nan,
        "generator_mae": float(generator["mae"]) if generator is not None else np.nan,
        "generator_mse": float(generator["mse"]) if generator is not None else np.nan,
        "child_u_z_noise_correlation": noise["correlation"],
        "child_noise_count": int(noise["count"]),
        "generator_bias_by_s_bin_json": json.dumps(
            _json_scalar(metadata["extra_diagnostics"]["generator_bias_by_s_bin"]),
            sort_keys=True,
            allow_nan=False,
        ),
    }


def load_stage(stage: str, tasks: list[dict[str, Any]]) -> pd.DataFrame:
    expected = [task_path(task) for task in tasks]
    missing = [path for path in expected if not path.exists()]
    if missing:
        examples = "\n".join(str(path) for path in missing[:5])
        raise RuntimeError(
            f"{stage} is incomplete: {len(missing)}/{len(expected)} artifacts missing. "
            f"First paths:\n{examples}"
        )
    frame = pd.DataFrame(_flatten_artifact(path, stage) for path in expected)
    if len(frame) != len(expected):
        raise AssertionError("artifact count changed during analysis")
    return add_repetition_averages(frame)


def add_repetition_averages(frame: pd.DataFrame) -> pd.DataFrame:
    keys = ["pde", "d", "n", "M", "method"]
    summaries: dict[tuple[Any, ...], dict[str, float]] = {}
    for key, group in frame.groupby(keys, sort=False):
        predictions: list[np.ndarray] = []
        truth: np.ndarray | None = None
        validation: np.ndarray | None = None
        for artifact in group["artifact"]:
            loaded = load_repetition(RESULTS_ROOT / artifact)
            predictions.append(np.asarray(loaded["prediction_u"], dtype=np.float64))
            if truth is None:
                truth = np.asarray(loaded["truth_u"], dtype=np.float64)
                validation = np.asarray(loaded["is_validation"], dtype=bool)
            else:
                if not np.array_equal(truth, loaded["truth_u"]):
                    raise RuntimeError(f"truth mismatch in paired group {key}")
                if not np.array_equal(validation, loaded["is_validation"]):
                    raise RuntimeError(f"split mismatch in paired group {key}")
        assert truth is not None and validation is not None
        stack = np.stack(predictions)
        mean_prediction = np.mean(stack, axis=0)
        test = ~validation
        ref_std = float(np.std(truth[test], ddof=0))
        if np.all(np.isfinite(stack)):
            sampling_variance_skill = float(
                np.sqrt(np.mean((stack[:, test] - mean_prediction[test]) ** 2)) / ref_std
            )
        else:
            sampling_variance_skill = float("inf")
        summaries[key] = {
            "mean_rep_skill": float(np.mean(group["test_skill"])),
            "median_rep_skill": float(np.median(group["test_skill"])),
            "repetition_averaged_skill": _finite_skill(mean_prediction, truth, test),
            "sampling_variation_skill": sampling_variance_skill,
        }
    result = frame.copy()
    for column in next(iter(summaries.values())):
        result[column] = [summaries[tuple(row[key] for key in keys)][column] for _, row in result.iterrows()]
    return result


def _limit_ratio(numerator: float, denominator: float) -> float:
    if denominator == 0.0:
        return float("inf") if numerator > 0.0 else float("nan")
    return numerator / denominator


def add_e1_cell_metrics(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    metrics: dict[tuple[Any, ...], dict[str, float]] = {}
    for key, group in frame.groupby(["pde", "d", "n", "M"], sort=False):
        means = group.groupby("method")["test_skill"].mean().to_dict()
        raw = float(means["raw"])
        box = float(means["box"])
        oracle = float(means["oracle_state"])
        if math.isinf(raw) and math.isfinite(oracle):
            damage = 1.0
            gap_closed = 1.0 if math.isfinite(box) else float("nan")
        else:
            damage = (raw - oracle) / raw
            gap_closed = (raw - box) / (raw - oracle)
        metrics[key] = {
            "damage_share_D": damage,
            "gap_closed_Gc": gap_closed,
            "box_over_oracle": _limit_ratio(box, oracle),
            "raw_over_oracle": _limit_ratio(raw, oracle),
        }
    for column in next(iter(metrics.values())):
        result[column] = [
            metrics[(row.pde, row.d, row.n, row.M)][column]
            for row in result.itertuples(index=False)
        ]
    return result


def add_e6_dimension_metrics(frame: pd.DataFrame) -> pd.DataFrame:
    """Attach the pre-requested box/oracle dimension diagnostic to every row."""

    result = frame.copy()
    ratios: dict[int, float] = {}
    for dimension, group in frame.groupby("d", sort=True):
        means = group.groupby("method")["test_skill"].mean().to_dict()
        ratios[int(dimension)] = _limit_ratio(
            float(means["box"]), float(means["oracle_state"])
        )
    result["box_over_oracle"] = [ratios[int(value)] for value in result["d"]]
    return result


def e6_dimension_check(frame: pd.DataFrame) -> dict[str, Any]:
    """Summarize dimension stability without inventing a post-hoc pass threshold."""

    ratio = frame.groupby("d")["box_over_oracle"].first().sort_index()
    values = ratio.to_numpy(dtype=np.float64)
    mean = float(np.mean(values))
    minimum = float(np.min(values))
    maximum = float(np.max(values))
    return {
        "metric": "ratio of mean box skill to mean oracle_state skill",
        "pre_registered_threshold": None,
        "ratio_by_dimension": {str(int(d)): float(value) for d, value in ratio.items()},
        "mean": mean,
        "minimum": minimum,
        "maximum": maximum,
        "relative_range": (maximum - minimum) / mean if mean != 0.0 else float("inf"),
        "coefficient_of_variation": (
            float(np.std(values, ddof=0)) / mean if mean != 0.0 else float("inf")
        ),
        "max_over_min": _limit_ratio(maximum, minimum),
    }


def tune_e3(frame: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    result = frame.copy()
    result["selected_by_validation"] = False
    tuning_transforms = ("shrink", "illegal", "looser", "shrink_centre")
    choices: list[dict[str, Any]] = []
    cell_keys = ["pde", "d", "n", "M"]
    for cell, group in result.groupby(cell_keys, sort=False):
        for transform in tuning_transforms:
            candidates = group[group["transform"] == transform]
            if candidates.empty:
                continue
            validation = (
                candidates.groupby(["method", "factor"], as_index=False)["validation_skill"]
                .mean()
                .sort_values(["validation_skill", "factor", "method"])
            )
            chosen = validation.iloc[0]
            mask = np.ones(len(result), dtype=bool)
            for column, value in zip(cell_keys, cell):
                mask &= result[column].to_numpy() == value
            mask &= result["method"].to_numpy() == chosen["method"]
            result.loc[mask, "selected_by_validation"] = True
            test_rows = result.loc[mask]
            choices.append(
                {
                    "pde": cell[0],
                    "d": int(cell[1]),
                    "n": int(cell[2]),
                    "M": int(cell[3]),
                    "family": transform,
                    "selected_method": chosen["method"],
                    "selected_factor": float(chosen["factor"]),
                    "mean_validation_skill": float(chosen["validation_skill"]),
                    "mean_test_skill": float(test_rows["test_skill"].mean()),
                    "candidate_validation_skills": {
                        row.method: float(row.validation_skill)
                        for row in validation.itertuples(index=False)
                    },
                }
            )
    payload = {
        "selection_rule": "minimum mean skill on the fixed 20% validation split within each (PDE,d,n,M) cell; test values were not consulted",
        "choices": choices,
    }
    return result, payload


def e3_win_table(frame: pd.DataFrame, tuning: dict[str, Any]) -> pd.DataFrame:
    choice_lookup = {
        (item["pde"], item["d"], item["n"], item["M"], item["family"]): item[
            "selected_method"
        ]
        for item in tuning["choices"]
    }
    records: list[dict[str, Any]] = []
    for cell, group in frame.groupby(["pde", "d", "n", "M"], sort=False):
        indexed = group.set_index(["method", "rep"])["test_skill"]
        box = indexed.loc["box"]
        comparators = {
            "best_shrink": choice_lookup[(*cell, "shrink")],
            "best_illegal": choice_lookup[(*cell, "illegal")],
            "centre": "centre",
            "best_shrink_centre": choice_lookup[(*cell, "shrink_centre")],
            "raw": "raw",
        }
        row: dict[str, Any] = {"pde": cell[0], "d": cell[1], "n": cell[2], "M": cell[3]}
        for label, method in comparators.items():
            comparison = indexed.loc[method]
            common = box.index.intersection(comparison.index)
            row[f"box_wins_vs_{label}"] = int(
                np.count_nonzero(box.loc[common].to_numpy() < comparison.loc[common].to_numpy())
            )
            row[f"comparator_{label}"] = method
            row[f"box_mean_minus_{label}_mean"] = float(box.mean() - comparison.mean())
        records.append(row)
    return pd.DataFrame(records)


def make_e4(frame: pd.DataFrame, tuning: dict[str, Any]) -> pd.DataFrame:
    choice_lookup = {
        (item["pde"], item["d"], item["n"], item["M"], item["family"]): item[
            "selected_method"
        ]
        for item in tuning["choices"]
    }
    keep = np.zeros(len(frame), dtype=bool)
    for index, row in frame.iterrows():
        if row["family"] in {"ridge_lse", "norm_hjb"}:
            allowed = {"span_only", "box", "ball", "segment"}
        else:
            looser = choice_lookup[(row["pde"], row["d"], row["n"], row["M"], "looser")]
            allowed = {"sign_only", "box", looser}
        keep[index] = row["method"] in allowed
    result = frame.loc[keep].copy()
    result["ablation_role"] = np.where(
        result["transform"] == "looser", "best_validation_looser", result["method"]
    )
    return result


def evaluate_claims(
    e1: pd.DataFrame,
    e2: pd.DataFrame,
    wins: pd.DataFrame,
    e5: pd.DataFrame,
) -> dict[str, Any]:
    verdicts: dict[str, Any] = {}
    for pde in sorted(e1["pde"].unique()):
        pde_e1 = e1[e1["pde"] == pde]
        cells = pde_e1.groupby(["d", "n", "M"], as_index=False).first()
        dimension_c1: dict[str, Any] = {}
        dimension_c2: dict[str, Any] = {}
        for d, group in cells.groupby("d"):
            c1_count = int(np.count_nonzero(group["damage_share_D"] >= 0.5))
            c2_count = int(
                np.count_nonzero(
                    group["box_over_oracle"] <= 0.5 * group["raw_over_oracle"]
                )
            )
            dimension_c1[str(int(d))] = {"passing_configs": c1_count, "passed": c1_count >= 3}
            dimension_c2[str(int(d))] = {"passing_configs": c2_count, "passed": c2_count >= 3}

        dose_details: list[dict[str, Any]] = []
        pde_e2 = e2[(e2["pde"] == pde) & (e2["transform"] == "dose")]
        for cell, group in pde_e2.groupby(["d", "n", "M"]):
            curve = group.groupby("factor")["test_skill"].mean().sort_index()
            values = curve.to_numpy()
            monotone = bool(np.all(np.diff(values) >= -1e-12))
            dose_details.append(
                {
                    "d": int(cell[0]),
                    "n": int(cell[1]),
                    "M": int(cell[2]),
                    "levels": curve.index.to_list(),
                    "mean_skills": values.tolist(),
                    "monotone": monotone,
                }
            )
        c1 = all(item["passed"] for item in dimension_c1.values()) and all(
            item["monotone"] for item in dose_details
        )
        c2 = all(item["passed"] for item in dimension_c2.values())

        pde_wins = wins[wins["pde"] == pde]
        comparator_columns = {
            "best tuned ->0 control": "box_wins_vs_best_shrink",
            "best tuned illegal set": "box_wins_vs_best_illegal",
            "centre": "box_wins_vs_centre",
            "best tuned shrink_centre": "box_wins_vs_best_shrink_centre",
        }
        c3_fractions = {
            label: float(np.mean(pde_wins[column] >= 7))
            for label, column in comparator_columns.items()
        }
        # The registered sentence says "against each ... in at least 75% of
        # E3 cells": all four comparisons therefore have to clear 7/10 in
        # the same cell.  Marginal fractions remain useful diagnostics but do
        # not replace this joint cell-level criterion.
        joint_pass = np.logical_and.reduce(
            [pde_wins[column].to_numpy() >= 7 for column in comparator_columns.values()]
        )
        c3_joint_fraction = float(np.mean(joint_pass))
        c3 = c3_joint_fraction >= 0.75
        centre_better = float(np.mean(pde_wins["box_mean_minus_centre_mean"] > 0.0))
        shrink_centre_better = float(
            np.mean(pde_wins["box_mean_minus_best_shrink_centre_mean"] > 0.0)
        )
        verdicts[pde] = {
            "C1_supported": c1,
            "C1_by_dimension": dimension_c1,
            "E2_raw_dose_curves": dose_details,
            "C2_supported": c2,
            "C2_by_dimension": dimension_c2,
            "C3_supported": c3,
            "C3_fraction_cells_with_at_least_7_wins": c3_fractions,
            "C3_fraction_cells_meeting_all_four": c3_joint_fraction,
            "centre_beats_box_fraction_cells": centre_better,
            "shrink_centre_beats_box_fraction_cells": shrink_centre_better,
            "certificate_prior_dominated": centre_better > 0.25
            or shrink_centre_better > 0.25,
        }

    negative: dict[str, Any] = {}
    for pde in ("N3", "N4"):
        subset = e5[e5["pde"] == pde]
        means = subset.groupby(["d", "n", "M", "method"])["test_skill"].mean().unstack()
        if pde == "N3":
            # "centre >= box" is interpreted as performance: lower skill is better.
            passed = bool(np.all(means["centre"] <= means["box"]))
            detail = "centre skill <= box skill in every cell"
        else:
            passed = bool(np.all(means["oracle_state"] > 0.10))
            detail = "oracle_state skill >0.10 in every cell"
        negative[pde] = {
            "passed": passed,
            "interpretation": detail,
            "cells": means.reset_index().to_dict(orient="records"),
        }
    negative.update(
        {
            "N1": {
                "passed": True,
                "interpretation": "historical cited result: f=0 is best; no new run",
            },
            "N2": {
                "passed": True,
                "interpretation": "historical cited result: z=0 beats certified corrections; no new run",
            },
            "N5": {
                "passed": True,
                "interpretation": "historical cited result: corrections are identical in value for f=0; no new run",
            },
        }
    )
    return {"primary": verdicts, "negative_controls": negative}


def _style() -> None:
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 150,
            "savefig.dpi": 180,
        }
    )


def plot_e1(e1: pd.DataFrame) -> None:
    methods = ["raw", "box", "segment", "oracle_state", "oracle_z", "oracle_u"]
    pdes = sorted(e1["pde"].unique())
    summary = e1.groupby(["pde", "method"])["test_skill"].median()
    x = np.arange(len(pdes), dtype=float)
    width = 0.8 / len(methods)
    fig, ax = plt.subplots(figsize=(max(9.0, 1.35 * len(pdes)), 5.2))
    colors = plt.cm.tab10(np.linspace(0.0, 0.8, len(methods)))
    for index, (method, color) in enumerate(zip(methods, colors)):
        values = [summary.get((pde, method), np.nan) for pde in pdes]
        ax.bar(x - 0.4 + width / 2 + index * width, values, width, label=method, color=color)
    ax.set_yscale("log")
    ax.set_ylabel("test skill (median over cells and repetitions; log scale)")
    ax.set_xticks(x, pdes, rotation=25, ha="right")
    ax.set_title("E1 oracle decomposition")
    ax.axhline(1.0, color="0.35", lw=0.8, ls="--", label="constant predictor")
    ax.legend(ncol=4, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / "e1_oracle_decomposition.png")
    plt.close(fig)


def plot_e2(e2: pd.DataFrame) -> None:
    pdes = sorted(e2["pde"].unique())
    columns = 3
    rows = math.ceil(len(pdes) / columns)
    fig, axes = plt.subplots(rows, columns, figsize=(12.5, 3.7 * rows), squeeze=False)
    transforms = [
        ("dose", "raw z dose"),
        ("dose_box", "z dose + box"),
        ("dose_segment", "z dose + segment"),
        ("dose_u", "raw u dose"),
        ("dose_u_box", "u dose + clip"),
    ]
    colors = {
        "dose": "#d62728",
        "dose_box": "#1f77b4",
        "dose_segment": "#2ca02c",
        "dose_u": "#9467bd",
        "dose_u_box": "#ff7f0e",
    }
    for ax, pde in zip(axes.flat, pdes):
        subset = e2[e2["pde"] == pde]
        smallest = int(subset["d"].min())
        subset = subset[subset["d"] == smallest]
        for (n, M), config in subset.groupby(["n", "M"]):
            for transform, label in transforms:
                selected = config[config["transform"] == transform]
                if selected.empty:
                    continue
                curve = selected.groupby("factor")["test_skill"].mean().sort_index()
                ax.plot(
                    curve.index,
                    curve.values,
                    marker="o",
                    color=colors[transform],
                    ls="-" if (n, M) == (3, 6) else "--",
                    label=f"{label}, ({n},{M})",
                )
        ax.set_yscale("log")
        ax.set_title(f"{pde}, d={smallest}")
        ax.set_xlabel("dose s")
        ax.set_ylabel("mean test skill")
        ax.grid(alpha=0.2)
    for ax in axes.flat[len(pdes) :]:
        ax.axis("off")
    legend: dict[str, Any] = {}
    for ax in axes.flat[: len(pdes)]:
        handles, labels = ax.get_legend_handles_labels()
        legend.update(zip(labels, handles))
    fig.legend(
        list(legend.values()),
        list(legend.keys()),
        loc="upper center",
        ncol=min(6, len(legend)),
        fontsize=8,
    )
    fig.suptitle("E2 one-shot dose response (no recursive feedback)", y=1.01)
    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / "e2_dose_response.png", bbox_inches="tight")
    plt.close(fig)


def plot_e2_raw_all_dimensions(e2: pd.DataFrame) -> None:
    """Show every raw z-dose curve used by the C1 monotonicity check."""

    raw = e2[e2["transform"] == "dose"]
    pdes = sorted(raw["pde"].unique())
    dimensions = sorted(int(value) for value in raw["d"].unique())
    colors = {
        dimension: plt.cm.viridis(index / max(len(dimensions) - 1, 1))
        for index, dimension in enumerate(dimensions)
    }
    columns = 3
    rows = math.ceil(len(pdes) / columns)
    fig, axes = plt.subplots(rows, columns, figsize=(12.5, 3.7 * rows), squeeze=False)
    for ax, pde in zip(axes.flat, pdes):
        subset = raw[raw["pde"] == pde]
        for (dimension, n, M), group in subset.groupby(["d", "n", "M"]):
            curve = group.groupby("factor")["test_skill"].mean().sort_index()
            ax.plot(
                curve.index,
                curve.values,
                marker="o",
                color=colors[int(dimension)],
                ls="-" if (n, M) == (3, 6) else "--",
                label=f"d={int(dimension)}, ({int(n)},{int(M)})",
            )
        ax.set_yscale("log")
        ax.set_title(pde)
        ax.set_xlabel("dose s")
        ax.set_ylabel("mean test skill")
        ax.grid(alpha=0.2)
    for ax in axes.flat[len(pdes) :]:
        ax.axis("off")
    legend: dict[str, Any] = {}
    for ax in axes.flat[: len(pdes)]:
        handles, labels = ax.get_legend_handles_labels()
        legend.update(zip(labels, handles))
    fig.legend(
        list(legend.values()),
        list(legend.keys()),
        loc="upper center",
        ncol=min(4, len(legend)),
        fontsize=8,
    )
    fig.suptitle("E2 raw z-dose curves used for C1 monotonicity", y=1.04)
    fig.tight_layout()
    fig.savefig(
        FIGURE_ROOT / "e2_raw_dose_all_dimensions.png", bbox_inches="tight"
    )
    plt.close(fig)


def plot_e3_wins(wins: pd.DataFrame) -> None:
    columns = [
        "box_wins_vs_best_shrink",
        "box_wins_vs_best_illegal",
        "box_wins_vs_centre",
        "box_wins_vs_best_shrink_centre",
        "box_wins_vs_raw",
    ]
    labels = ["best shrink", "best illegal", "centre", "best shrink-centre", "raw"]
    ordered = wins.sort_values(["pde", "d", "n", "M"]).reset_index(drop=True)
    matrix = ordered[columns].to_numpy(dtype=float)
    height = max(6.0, 0.22 * len(ordered))
    fig, ax = plt.subplots(figsize=(8.5, height))
    image = ax.imshow(matrix, vmin=0, vmax=10, cmap="RdYlGn", aspect="auto")
    ax.set_xticks(np.arange(len(labels)), labels, rotation=25, ha="right")
    ylabels = [f"{r.pde} d{r.d} ({r.n},{r.M})" for r in ordered.itertuples(index=False)]
    ax.set_yticks(np.arange(len(ylabels)), ylabels, fontsize=7)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(j, i, str(int(matrix[i, j])), ha="center", va="center", fontsize=7)
    fig.colorbar(image, ax=ax, label="paired box wins out of 10")
    ax.set_title("E3 paired win counts")
    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / "e3_win_count_heatmap.png")
    plt.close(fig)


def plot_e4(e4: pd.DataFrame) -> None:
    summary = (
        e4.groupby(["pde", "ablation_role"])["test_skill"]
        .median()
        .rename("median_skill")
        .reset_index()
    )
    pdes = sorted(summary["pde"].unique())
    roles = sorted(summary["ablation_role"].unique())
    x = np.arange(len(pdes), dtype=float)
    width = 0.8 / max(len(roles), 1)
    fig, ax = plt.subplots(figsize=(max(9.0, 1.4 * len(pdes)), 5.0))
    for index, role in enumerate(roles):
        lookup = summary[summary["ablation_role"] == role].set_index("pde")["median_skill"]
        values = [lookup.get(pde, np.nan) for pde in pdes]
        ax.bar(x - 0.4 + width / 2 + index * width, values, width, label=role)
    ax.set_yscale("log")
    ax.set_xticks(x, pdes, rotation=25, ha="right")
    ax.set_ylabel("median test skill (log scale)")
    ax.set_title("E4 certificate ablation")
    ax.legend(ncol=3, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / "e4_certificate_ablation.png")
    plt.close(fig)


def plot_e6(e6: pd.DataFrame) -> None:
    fig, (ax, ratio_ax) = plt.subplots(
        1, 2, figsize=(11.2, 4.7), gridspec_kw={"width_ratios": [1.35, 1.0]}
    )
    for method, group in e6.groupby("method"):
        curve = group.groupby("d")["test_skill"].mean().sort_index()
        ax.plot(curve.index, curve.values, marker="o", label=method)
    ax.set_yscale("log")
    ax.set_xlabel("dimension d")
    ax.set_ylabel("mean test skill (log scale)")
    ax.set_title("E6 P1 dimension sweep, (n,M)=(4,6)")
    ax.legend(ncol=3, fontsize=8)
    ax.grid(alpha=0.2)
    ratio = e6.groupby("d")["box_over_oracle"].first().sort_index()
    ratio_ax.plot(ratio.index, ratio.values, marker="o", color="#1f77b4")
    ratio_ax.set_xlabel("dimension d")
    ratio_ax.set_ylabel("box / oracle_state (mean skill ratio)")
    ratio_ax.set_title("Requested dimension-independence check")
    ratio_ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / "e6_dimension.png")
    plt.close(fig)


def p4_s_bin_table(e1: pd.DataFrame) -> pd.DataFrame:
    subset = e1[(e1["pde"] == "P4") & (e1["method"] == "raw")]
    aggregates: dict[tuple[float, float], dict[str, float]] = {}
    for text in subset["generator_bias_by_s_bin_json"]:
        for item in json.loads(text):
            key = (float(item["low"]), float(item["high"]))
            aggregate = aggregates.setdefault(key, {"count": 0.0, "bias_sum": 0.0, "mae_sum": 0.0})
            count = int(item["count"])
            if count and item["bias"] is not None:
                aggregate["count"] += count
                aggregate["bias_sum"] += count * float(item["bias"])
                aggregate["mae_sum"] += count * float(item["mae"])
    rows = []
    for (low, high), item in sorted(aggregates.items()):
        count = item["count"]
        rows.append(
            {
                "s_low": low,
                "s_high": high,
                "count": int(count),
                "generator_bias": item["bias_sum"] / count if count else np.nan,
                "generator_mae": item["mae_sum"] / count if count else np.nan,
            }
        )
    return pd.DataFrame(rows)


def plot_p4_bins(table: pd.DataFrame) -> None:
    if table.empty:
        return
    labels = [f"[{r.s_low:g},{r.s_high:g})" for r in table.itertuples(index=False)]
    x = np.arange(len(table))
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.bar(x - 0.18, table["generator_bias"], width=0.36, label="bias")
    ax.bar(x + 0.18, table["generator_mae"], width=0.36, label="MAE")
    ax.set_xticks(x, labels, rotation=30, ha="right")
    ax.set_ylabel("raw generator error")
    ax.set_title("P4 raw generator error by ridge coordinate s")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / "p4_generator_bias_by_s_bin.png")
    plt.close(fig)


def make_figures(
    e1: pd.DataFrame,
    e2: pd.DataFrame,
    wins: pd.DataFrame,
    e4: pd.DataFrame,
    e6: pd.DataFrame,
    p4_bins: pd.DataFrame,
) -> None:
    FIGURE_ROOT.mkdir(parents=True, exist_ok=True)
    _style()
    plot_e1(e1)
    plot_e2(e2)
    plot_e2_raw_all_dimensions(e2)
    plot_e3_wins(wins)
    plot_e4(e4)
    plot_e6(e6)
    plot_p4_bins(p4_bins)


def _format(value: Any) -> str:
    if value is None:
        return "—"
    if isinstance(value, (bool, np.bool_)):
        return "PASS" if value else "FAIL"
    if isinstance(value, (float, np.floating)):
        if math.isnan(float(value)):
            return "—"
        if math.isinf(float(value)):
            return "∞" if value > 0 else "−∞"
        magnitude = abs(float(value))
        if magnitude != 0.0 and (magnitude < 1e-3 or magnitude >= 1e4):
            return f"{float(value):.3e}"
        return f"{float(value):.4f}"
    return str(value).replace("|", "\\|")


def _markdown_table(frame: pd.DataFrame, columns: list[str] | None = None) -> str:
    data = frame if columns is None else frame[columns]
    headers = [str(column) for column in data.columns]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in data.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(_format(value) for value in row) + " |")
    return "\n".join(lines)


def _e1_summary(e1: pd.DataFrame) -> pd.DataFrame:
    means = e1.groupby(["pde", "d", "n", "M", "method"])["test_skill"].mean().unstack()
    medians = (
        e1.groupby(["pde", "d", "n", "M", "method"])["test_skill"]
        .median()
        .unstack()
        .add_suffix("_median")
    )
    averaged = (
        e1.groupby(["pde", "d", "n", "M", "method"])["repetition_averaged_skill"]
        .first()
        .unstack()
        .add_suffix("_repavg")
    )
    gradients = (
        e1.groupby(["pde", "d", "n", "M", "method"])["test_gradient_relative_l2"]
        .mean()
        .unstack()
        .add_suffix("_gradrel")
    )
    base = e1.groupby(["pde", "d", "n", "M"], as_index=True).first()
    result = (
        means.join(medians)
        .join(averaged)
        .join(gradients)
        .join(base[["damage_share_D", "gap_closed_Gc", "box_over_oracle", "raw_over_oracle"]])
        .reset_index()
    )
    desired = [
        "pde", "d", "n", "M", "raw", "raw_median", "raw_repavg", "box",
        "box_repavg", "segment", "oracle_state", "oracle_state_repavg", "oracle_z",
        "oracle_u", "raw_gradrel", "box_gradrel", "oracle_state_gradrel",
        "damage_share_D", "gap_closed_Gc", "box_over_oracle",
    ]
    for column in desired:
        if column not in result:
            result[column] = np.nan
    return result[desired]


def write_report(
    gates: dict[str, Any],
    e1: pd.DataFrame,
    e2: pd.DataFrame,
    e3: pd.DataFrame,
    wins: pd.DataFrame,
    e4: pd.DataFrame,
    e5: pd.DataFrame,
    e6: pd.DataFrame,
    p4_bins: pd.DataFrame,
    verdicts: dict[str, Any],
) -> None:
    gate_rows = []
    for pde, item in gates["candidate_summary"].items():
        gate_rows.append(
            {
                "PDE": pde,
                **{gate: item["gate_status"][gate] for gate in ("G1", "G2", "G3", "G4", "G5", "G6")},
                "verdict": "admit" if item["admitted"] else "exclude",
                "failed": ", ".join(item["failed_gates"]) or "none",
            }
        )
    for pde, item in gates.get("historical_negative_controls", {}).items():
        row: dict[str, Any] = {
            "PDE": pde,
            **{gate: "historical—not rerun" for gate in ("G1", "G2", "G3", "G4", "G5", "G6")},
            "verdict": "historical negative control",
            "failed": item["failed_gate"],
        }
        row[item["failed_gate"]] = False
        gate_rows.append(row)
    gate_table = pd.DataFrame(gate_rows)
    e1_summary = _e1_summary(e1)
    e4_summary = (
        e4.groupby(["pde", "ablation_role"])["test_skill"]
        .agg(["mean", "median"])
        .reset_index()
    )
    e5_summary = (
        e5.groupby(["pde", "d", "n", "M", "method"])["test_skill"]
        .mean()
        .unstack()
        .reset_index()
    )
    e6_summary = (
        e6.groupby(["d", "method"])["test_skill"].mean().unstack().reset_index()
    )
    e6_ratio = e6.groupby("d", as_index=False)["box_over_oracle"].first()
    e6_summary = e6_summary.merge(e6_ratio, on="d", validate="one_to_one")
    e6_check = e6_dimension_check(e6)
    noise_summary = (
        e1.groupby(["pde", "method"])["child_u_z_noise_correlation"]
        .mean()
        .unstack()
        .reset_index()
    )
    p2_channel = (
        e1[
            e1["pde"].str.startswith("P2")
            & e1["method"].isin(["raw", "oracle_state", "oracle_z", "oracle_u"])
        ]
        .groupby(["pde", "d", "method"])[
            ["test_skill", "child_u_z_noise_correlation", "generator_bias", "generator_mae"]
        ]
        .mean()
        .reset_index()
    )
    work = pd.concat([e1, e2, e3, e5, e6], ignore_index=True)
    # E0 is outside this concatenated accounting table, so an E5 or E6 row
    # sourced from E0 represents real work not counted elsewhere here.  Reuses
    # sourced from E1/E3 are duplicates of rows already present and excluded.
    unique_work = work[
        work["reused_from_stage"].isna() | (work["reused_from_stage"] == "e0")
    ]
    nonfinite_rows = int(np.count_nonzero(work["nonfinite_count"] > 0))
    total_wall = float(unique_work["wall_time_seconds"].sum())
    total_f = int(unique_work["f_calls"].sum())
    total_states = int(unique_work["recursively_evaluated_states"].sum())
    total_samples = int(unique_work["stochastic_samples"].sum())
    total_dose_normals = int(unique_work["dose_normal_variates"].sum())
    reuse_rows = int(work["reused_from_stage"].notna().sum())
    superseded_count = sum(
        1 for _ in (RESULTS_ROOT / "raw" / "superseded").rglob("*.npz")
    )
    derivative_validation_count = sum(
        1 for _ in (RESULTS_ROOT / "raw" / "derivative_validation").rglob("*.npz")
    )

    primary_rows = []
    for pde, item in verdicts["primary"].items():
        primary_rows.append(
            {
                "PDE": pde,
                "C1": item["C1_supported"],
                "C2": item["C2_supported"],
                "C3": item["C3_supported"],
                "C3 joint-cell fraction": item["C3_fraction_cells_meeting_all_four"],
                "certificate-prior dominated": (
                    "yes" if item["certificate_prior_dominated"] else "no"
                ),
                "centre better fraction": item["centre_beats_box_fraction_cells"],
                "shrink-centre better fraction": item[
                    "shrink_centre_beats_box_fraction_cells"
                ],
            }
        )
    primary_count = len(primary_rows)
    claim_counts = {
        claim: sum(
            bool(item[f"{claim}_supported"])
            for item in verdicts["primary"].values()
        )
        for claim in ("C1", "C2", "C3")
    }

    failures: list[str] = []
    for row in gate_rows:
        if row["failed"] != "none":
            failures.append(f"{row['PDE']} failed gate(s) {row['failed']} and was excluded from C1–C3.")
    for pde, item in verdicts["primary"].items():
        if not item["C1_supported"]:
            weak_dimensions = [
                f"d={dimension}: {detail['passing_configs']}/4 D cells"
                for dimension, detail in item["C1_by_dimension"].items()
                if not detail["passed"]
            ]
            nonmonotone = [
                f"d={curve['d']} ({curve['n']},{curve['M']})"
                for curve in item["E2_raw_dose_curves"]
                if not curve["monotone"]
            ]
            parts = []
            if weak_dimensions:
                parts.append("; ".join(weak_dimensions))
            if nonmonotone:
                parts.append("non-monotone raw dose: " + ", ".join(nonmonotone))
            failures.append(
                f"{pde}: pre-registered C1 criterion failed ({'; '.join(parts)})."
            )
        if not item["C2_supported"]:
            weak_dimensions = [
                f"d={dimension}: {detail['passing_configs']}/4 cells"
                for dimension, detail in item["C2_by_dimension"].items()
                if not detail["passed"]
            ]
            failures.append(
                f"{pde}: pre-registered C2 criterion failed ({'; '.join(weak_dimensions)})."
            )
        if not item["C3_supported"]:
            failures.append(
                f"{pde}: pre-registered C3 criterion failed "
                f"(joint passing-cell fraction={_format(item['C3_fraction_cells_meeting_all_four'])}; "
                f"required >=0.75)."
            )
        if item["certificate_prior_dominated"]:
            failures.append(
                f"{pde}: certificate-prior dominated in the E3 regime "
                f"(centre fraction={_format(item['centre_beats_box_fraction_cells'])}, "
                f"shrink-centre fraction={_format(item['shrink_centre_beats_box_fraction_cells'])})."
            )
    for pde, item in verdicts["negative_controls"].items():
        if not item["passed"]:
            failures.append(
                f"{pde}: negative-control prediction failed ({item['interpretation']})."
            )
    if nonfinite_rows:
        failures.append(f"{nonfinite_rows} method-repetition rows contained non-finite states or generator values.")
    p4_bias_direction = "P4 was not admitted, so no generator-bias direction is available."
    p4_localization = "P4 was not admitted, so the rectification localization check is unavailable."
    if not p4_bins.empty:
        p4_bias_direction = (
            "The bias is negative in every populated bin, matching the direction "
            "expected when noisy gradients inflate the norm inside the negative norm driver."
            if bool(np.all(p4_bins["generator_bias"] < 0.0))
            else "The generator-bias sign is not uniform across the populated P4 bins."
        )
        central_mask = np.isclose(p4_bins["s_low"], -0.25) & np.isclose(
            p4_bins["s_high"], 0.25
        )
        if np.count_nonzero(central_mask) == 1:
            central = p4_bins.loc[central_mask].iloc[0]
            outer = p4_bins.loc[~central_mask]
            largest_outer = outer.loc[outer["generator_bias"].abs().idxmax()]
            localization_passed = bool(
                abs(central["generator_bias"]) >= abs(largest_outer["generator_bias"])
            )
            p4_localization = (
                f"The central-bin raw generator bias/MAE is "
                f"{_format(central['generator_bias'])}/{_format(central['generator_mae'])}; "
                f"the largest outer-bin bias magnitude is "
                f"{_format(abs(largest_outer['generator_bias']))} "
                f"(MAE {_format(largest_outer['generator_mae'])}) "
                f"on [{_format(largest_outer['s_low'])},{_format(largest_outer['s_high'])}). "
                f"Thus the pre-stated localization prediction that the Jensen damage is "
                f"largest near psi_s=0 {'is supported' if localization_passed else 'is not supported'} "
                f"by the raw recursive E1 diagnostic."
            )
            if not localization_passed:
                failures.append(
                    "P4: the pre-stated rectification-localization prediction was not "
                    f"observed (central-bin |generator bias|="
                    f"{_format(abs(central['generator_bias']))}; largest outer-bin "
                    f"|generator bias|={_format(abs(largest_outer['generator_bias']))})."
                )
    if not failures:
        failures.append("No gate, claim, negative-control, or non-finiteness failure was observed.")

    p4_screen = gates["p4_parameter_screen"]
    selected_p4 = p4_screen["selected"]
    rejected_g2 = p4_screen.get("screened_out", {}).get("gates", {}).get("G2", {}).get(
        "time_only_skill", np.nan
    )
    selected_gate = gates["candidate_summary"].get("P4", {}).get("instances", {}).get("20", {})
    selected_g2 = selected_gate.get("analytic", {}).get("G2", {}).get("time_only_skill", np.nan)
    selected_g1 = selected_gate.get("analytic", {}).get("G1", {})
    selected_ref = selected_g1.get("reference_refinement", {})
    problem_table = pd.DataFrame(
        [
            {
                "ID": "P1",
                "parameters": "ridge LSE; sigma=sqrt(2), T=0.25, kappa=1",
                "d": "20,50,100,200",
                "certificate": "PDE-derived segment, box, ball, span",
            },
            {
                "ID": "P2_a4",
                "parameters": "VB-a; a=4, sigma=0.5, T=0.5",
                "d": "20,50,100",
                "certificate": "u/sign PDE-derived; z cap solution-informed",
            },
            {
                "ID": "P2_a8",
                "parameters": "VB-a; a=8, sigma=0.5, T=0.5",
                "d": "20,50,100",
                "certificate": "u/sign PDE-derived; z cap solution-informed",
            },
            {
                "ID": "P3_rho1",
                "parameters": "Burgers-Fisher; a=4, rho=1, sigma=0.5, T=0.5",
                "d": "20,50",
                "certificate": "u/sign PDE-derived; z cap solution-informed",
            },
            {
                "ID": "P3_rho2",
                "parameters": "Burgers-Fisher; a=4, rho=2, sigma=0.5, T=0.5",
                "d": "20,50",
                "certificate": "u/sign PDE-derived; z cap solution-informed",
            },
            {
                "ID": "P4",
                "parameters": "norm HJB; beta=2, lambda_f=1, T=0.5",
                "d": "20,50,100",
                "certificate": "PDE-derived segment, box, ball, span",
            },
            {
                "ID": "N3",
                "parameters": "5-direction LSE; strength=0.5, T=0.5",
                "d": "100",
                "certificate": "PDE-derived convex hull, box, ball, span",
            },
            {
                "ID": "N4",
                "parameters": "published VB; a=1, sigma=sqrt(2), T=0.5",
                "d": "20,50,100",
                "certificate": "u/sign PDE-derived; z cap solution-informed",
            },
        ]
    )
    study_table = pd.DataFrame(
        [
            {
                "study": "E0",
                "design": "all candidates; (3,6),(4,6); 3 reps",
                "purpose": "hard admission gates",
            },
            {
                "study": "E1",
                "design": "4 (n,M) cells; all admitted d; 10 reps",
                "purpose": "oracle/channel decomposition",
            },
            {
                "study": "E2",
                "design": "(4,3),(3,6); six doses; 10 reps",
                "purpose": "one-shot causal dose response",
            },
            {
                "study": "E3",
                "design": "4 cells; two smallest d; 10 reps; validation-only tuning",
                "purpose": "geometry and suppression controls",
            },
            {
                "study": "E4",
                "design": "derived from untuned/tuned E3 geometry rows",
                "purpose": "certificate component ablation",
            },
            {
                "study": "E5",
                "design": "N3/N4; (3,6),(4,6); 10 reps",
                "purpose": "negative-control predictions",
            },
            {
                "study": "E6",
                "design": "P1; (4,6); four d; 10 reps",
                "purpose": "box/oracle dimension stability",
            },
        ]
    )

    lines = [
        "# Mechanism benchmark suite: noise entering a nonlinear generator",
        "",
        "## Outcome",
        "",
        _markdown_table(pd.DataFrame(primary_rows)),
        "",
        f"Across the {primary_count} admitted problem variants, C1 is supported for {claim_counts['C1']}/{primary_count}, C2 for {claim_counts['C2']}/{primary_count}, and C3 for {claim_counts['C3']}/{primary_count}. The suite therefore gives broad support to nonlinear-generator noise damage, but not to universal half-error removal by box retraction or to the claim that the box advantage generally requires the data rather than a tuned certificate prior.",
        "",
        "The table applies the pre-registered thresholds mechanically; a failed criterion is not softened by qualitative trends. Lower skill is better, and a constant predictor has skill 1.",
        "",
        "## Protocol and provenance",
        "",
        f"All computations used float64, 1,200 fixed points with a fixed 20% validation split, paired tree seeds `SeedSequence([base_seed,d,n,M,rep,chunk_index])`, and identical chunk sizes within each study. E0 uses 8 at d=20 and 4 otherwise; E1/E3/E5/E6 use 16 at d=20 and 4 otherwise; E2 uses 16/8/4 at d=20, d=50–100, and d=200 respectively. The unchanged `FullHistoryMLP` from `experiments/active_vb_high_budget/vb_mlp_methods.py` supplies the recursion; only the state immediately handed to f is transformed, and the root is never clipped. Base seed: {gates['protocol']['base_seed']}.",
        "",
        f"The gate snapshot records code commit `{gates['environment']['git_commit']}` on branch `{gates['environment']['git_branch']}`, Python {gates['environment']['python'].split()[0]}, NumPy {gates['environment']['numpy']}, SciPy {gates['environment'].get('scipy', 'not recorded')}, and {gates['environment']['cpu_count']} logical CPUs. The unchanged recursion source has SHA-256 `{gates.get('unchanged_recursion_source', {}).get('sha256', 'not recorded')}`.",
        "",
        "P2, P3, and N4 use the inherited VB test geometry with 1,000 interior and 200 face-boundary points; P1, P4, and N3 use t uniform on [0,T) and x uniform on [-1,1]^d. E0 containment independently samples 100,000 points from the corresponding distribution, including the VB boundary stratum.",
        "",
        "### Problem matrix",
        "",
        _markdown_table(problem_table),
        "",
        "### Study matrix",
        "",
        _markdown_table(study_table),
        "",
        "The supplied attachment did not contain the referenced `reference_code/` directory. P1–P3 and N4 were implemented directly from the equations in the protocol; N3 is explicitly marked as a canonical mathematical reconstruction. This missing bundle is a reproducibility limitation, not hidden.",
        f"During E0 optimization, {superseded_count} preliminary P4 artifacts using FITPACK's direct pointwise derivative evaluator were preserved under `raw/superseded/`, and {derivative_validation_count} paired comparison artifact under `raw/derivative_validation/`; all are excluded from every result table. The accepted artifacts use the audited cached derivative spline (maximum discrepancy from the direct spline derivative is recorded in G1 and the [paired runtime audit](../results/mechanism_suite/reference_cache/derivative_evaluator_audit.json)).",
        "",
        f"P4 used beta={selected_p4['beta']:g}, lambda_f={selected_p4['lambda_f']:g}, T={selected_p4['T']:g}. Lambda_f=1 was the stronger pre-listed suggestion and was fixed before inspecting any MLP result; the lambda_f=0.5 alternatives were not exhaustively screened. The suggested T=0.25 setting was screened out because G2={_format(rejected_g2)}<0.15; T=0.5 gave G2={_format(selected_g2)}. Its six-grid monotone upwind reference has common-node Richardson disagreement {_format(selected_ref.get('richardson_refinement_max_difference'))}, stored on {selected_ref.get('stored_time_points','—')}×{selected_ref.get('stored_space_points','—')} points. The unextrapolated finest-pair difference is {_format(selected_ref.get('raw_finest_pair_max_difference'))}; because the monotone scheme is first order, the formal reference and G1 comparison use successive common-node third-Richardson estimates, and both numbers are reported.",
        f"For P4 d=20, the spline-differential residual is {_format(selected_g1.get('analytic', {}).get('max_abs_residual'))}, the worst value in the three-step independent fourth-order FD audit is {_format(selected_g1.get('finite_difference', {}).get('max_abs_residual'))}, and the central gradient-check error is {_format(selected_g1.get('terminal_gradient', {}).get('gradient_max_abs_error'))}.",
        "",
        "P4 had no dimension list in the supplied protocol, so E0/E1/E2 use d={20,50,100}; the compute-limited E3 panel uses d={20,50}. This choice was fixed before MLP results were inspected.",
        "",
        "## E0: hard gates",
        "",
        _markdown_table(gate_table),
        "",
        "N3 was predicted in advance to fail only G2; N4 was predicted to fail only G3. Only primary PDEs passing all six gates enter C1–C3.",
        "[Full gate diagnostics, containment counts, certificate labels, and reference audits](../results/mechanism_suite/gates.json)",
        "",
        "## Pre-registered acceptance criteria (verbatim)",
        "",
    ]
    lines.extend(f"- {criterion}" for criterion in ACCEPTANCE_CRITERIA)
    lines.extend(
        [
            "",
            "For C2, the written ratio inequality is algebraically stricter than the parenthetical ‘half of excess error’ statement. Verdicts use the explicit pre-registered inequality `box/oracle <= 0.5 * raw/oracle`; `Gc` is reported separately so the alternative excess-gap reading remains visible.",
            "",
            "## E1: oracle decomposition",
            "",
            "Per-cell entries below are means over 10 paired repetitions. `D=(raw-oracle)/raw`; `Gc=(raw-box)/(raw-oracle)`. The CSV also contains every per-repetition skill, median skill, repetition-averaged skill, gradient relative L2, work counters, violations, activations, and non-finite counts.",
            "",
            _markdown_table(e1_summary),
            "",
            "[Full E1 rows](../results/mechanism_suite/e1_oracle_decomposition.csv) · [Oracle bars](../results/mechanism_suite/figures/e1_oracle_decomposition.png)",
            "",
            "### Child u/z noise correlation",
            "",
            "The correlation is reported rather than suppressing the known P2 `oracle_u` anomaly. The z scalar is the mean coordinate error, proportional to the sum-z error used by the product driver.",
            "",
            _markdown_table(noise_summary),
            "",
            "Replacing only u can remove covariance-driven cancellation while leaving the noisy sum-z channel intact, but the correlation is dimension-dependent and need not stay negative. The following dimension-resolved table exposes the sign and magnitude together with skill and generator bias. Where the correlation weakens or reverses, covariance cancellation cannot by itself explain the adverse `oracle_u` result.",
            "",
            _markdown_table(p2_channel),
            "",
            "## E2: causal one-shot dose response",
            "",
            f"At each generator evaluation the child estimate is replaced by exact state plus the registered perturbation, so a dose-induced state error never enters a later nonlinear generator (the resulting estimator contributions can still aggregate linearly). The P3 u-channel uses scale {U_DOSE_SCALE:g}, the half-width of its certified interval [0,1], because the protocol requested that channel but did not state a scale. Raw-dose monotonicity is evaluated separately for both registered configurations and every admitted dimension in `analysis_summary.json`.",
            "",
            "[Dose rows](../results/mechanism_suite/e2_dose.csv) · [Method-comparison curves](../results/mechanism_suite/figures/e2_dose_response.png) · [All raw-dose dimensions](../results/mechanism_suite/figures/e2_raw_dose_all_dimensions.png)",
            "",
            "## E3: geometry or suppression?",
            "",
            "Every tuned factor was chosen solely by mean validation skill within its (PDE,d,n,M) cell, then frozen for the test split. The table gives strict paired box wins out of 10; ties are not wins.",
            "C3 is evaluated jointly: a cell passes only when box records at least 7/10 wins against all four registered comparators in that same cell; at least 75% of a PDE's E3 cells must pass. Per-comparator marginal fractions are retained in `analysis_summary.json` as diagnostics.",
            "The certificate-prior-dominated flag uses cell-level mean test skill: a prior control ‘beats’ box when its mean is lower.",
            "",
            _markdown_table(wins),
            "",
            "[Control rows](../results/mechanism_suite/e3_controls.csv) · [Frozen choices](../results/mechanism_suite/e3_tuning_choices.json) · [Win heatmap](../results/mechanism_suite/figures/e3_win_count_heatmap.png)",
            "",
            "## E4: certificate ablation",
            "",
            _markdown_table(e4_summary),
            "",
            "[Ablation rows](../results/mechanism_suite/e4_ablation.csv) · [Ablation bars](../results/mechanism_suite/figures/e4_certificate_ablation.png)",
            "",
            "## E5: negative controls",
            "",
            _markdown_table(e5_summary),
            "",
            "N3's statement “centre >= box” is interpreted as predictive performance, hence centre skill <= box skill. N4 passes its prediction only when oracle-state skill remains above 0.10 in every tested cell.",
            "",
            f"N3 prediction: {_format(verdicts['negative_controls']['N3']['passed'])}. N4 prediction: {_format(verdicts['negative_controls']['N4']['passed'])}.",
            "",
            "Historical no-new-run controls: N1 is documented in [HJB_LIFE_OR_DEATH_ABLATION.md](HJB_LIFE_OR_DEATH_ABLATION.md) and [hjb_life_or_death_summary.json](../results/hjb_life_or_death_summary.json); N2 in [FUNDING_LIFE_OR_DEATH_ABLATION.md](FUNDING_LIFE_OR_DEATH_ABLATION.md) and [funding_life_or_death_summary.json](../results/funding_life_or_death_summary.json); N5 in [batchir_negative_controls.md](batchir_negative_controls.md). Their supplied gate diagnoses (G2 for N1/N2, G4 for N5) were cited, not rerun.",
            "",
            "[Negative-control rows](../results/mechanism_suite/e5_negative.csv)",
            "",
            "## E6: P1 dimension sweep",
            "",
            _markdown_table(e6_summary),
            "",
            f"No numerical threshold for dimension-independence was pre-registered. The observed mean-skill `box/oracle_state` ratio ranges from {_format(e6_check['minimum'])} to {_format(e6_check['maximum'])}; relative range={_format(e6_check['relative_range'])}, coefficient of variation={_format(e6_check['coefficient_of_variation'])}, and max/min={_format(e6_check['max_over_min'])}.",
            "",
            "[Dimension rows](../results/mechanism_suite/e6_dimension.csv) · [Dimension plot](../results/mechanism_suite/figures/e6_dimension.png)",
            "",
            "## P4 rectification mechanism",
            "",
            "The table pools raw E1 generator calls and weights each s-bin by its number of calls; the central bin [-0.25,0.25) is the pre-specified neighborhood of the rectification point psi_s=0.",
            "",
            _markdown_table(p4_bins) if not p4_bins.empty else "P4 was not admitted, so no s-bin result exists.",
            "",
            p4_bias_direction,
            p4_localization,
            "The raw recursive diagnostic combines Jensen rectification with state-dependent recursive variance, so a failed localization pattern does not contradict the controlled E2 dose response; it does reject the stronger claim that raw E1 damage is maximized at the rectification point.",
            "",
            "[P4 s-bin plot](../results/mechanism_suite/figures/p4_generator_bias_by_s_bin.png)",
            "",
            "## Failures and adverse findings",
            "",
        ]
    )
    lines.extend(f"- {item}" for item in failures)
    lines.extend(
        [
            "",
            "## Work accounting",
            "",
            f"Across E1/E2/E3/E5/E6: {len(work):,} method-repetition rows ({reuse_rows:,} exact artifact reuses), {total_f:,} non-duplicated generator calls, {total_states:,} recursively evaluated states, {total_samples:,} MLP stochastic samples, {total_dose_normals:,} additional registered dose-normal variates, and {total_wall/3600:.2f} non-duplicated summed worker-hours. E0-sourced E5/E6 artifacts count once here because E0 is outside this table. Rows with any non-finite state/generator count: {nonfinite_rows}.",
            "",
            "## Reproduction",
            "",
            "```powershell",
            "python -m invariant_region_mlp.experiments.mechanism_suite.run_suite --stage all --workers 8",
            "python -m invariant_region_mlp.experiments.mechanism_suite.analyze_suite",
            "```",
            "",
            "The runner is resumable and validates existing per-repetition artifacts before skipping them.",
        ]
    )
    DOC_PATH.parent.mkdir(parents=True, exist_ok=True)
    temporary = DOC_PATH.with_name(DOC_PATH.name + ".tmp")
    temporary.write_text("\n".join(lines) + "\n", encoding="utf-8")
    temporary.replace(DOC_PATH)


def analyze() -> None:
    if not GATES_PATH.exists():
        raise RuntimeError("gates.json is missing; E0 must finish before analysis")
    gates = json.loads(GATES_PATH.read_text(encoding="utf-8"))
    e1 = add_e1_cell_metrics(load_stage("e1", e1_tasks()))
    e2 = load_stage("e2", e2_tasks())
    e2["dose_level"] = e2["factor"]
    e2["dose_channel"] = np.where(e2["transform"].str.startswith("dose_u"), "u", "z")
    e3, tuning = tune_e3(load_stage("e3", e3_tasks()))
    wins = e3_win_table(e3, tuning)
    e4 = make_e4(e3, tuning)
    e5 = load_stage("e5", e5_tasks())
    e6 = add_e6_dimension_metrics(load_stage("e6", e6_tasks()))
    p4_bins = p4_s_bin_table(e1)
    verdicts = evaluate_claims(e1, e2, wins, e5)

    sort_columns = ["pde", "d", "n", "M", "method", "rep"]
    for frame in (e1, e2, e3, e4, e5, e6):
        frame.sort_values(sort_columns, inplace=True, ignore_index=True)
    _atomic_csv(e1, RESULTS_ROOT / "e1_oracle_decomposition.csv")
    _atomic_csv(e2, RESULTS_ROOT / "e2_dose.csv")
    _atomic_csv(e3, RESULTS_ROOT / "e3_controls.csv")
    _atomic_json(tuning, TUNING_PATH)
    _atomic_csv(e4, RESULTS_ROOT / "e4_ablation.csv")
    _atomic_csv(e5, RESULTS_ROOT / "e5_negative.csv")
    _atomic_csv(e6, RESULTS_ROOT / "e6_dimension.csv")
    _atomic_csv(p4_bins, RESULTS_ROOT / "p4_generator_bias_by_s_bin.csv")
    summary = {
        "criteria": list(ACCEPTANCE_CRITERIA),
        "verdicts": verdicts,
        "e3_win_counts": wins.to_dict(orient="records"),
        "e6_dimension_check": e6_dimension_check(e6),
        "row_counts": {
            "e1": len(e1),
            "e2": len(e2),
            "e3": len(e3),
            "e4": len(e4),
            "e5": len(e5),
            "e6": len(e6),
        },
    }
    _atomic_json(summary, SUMMARY_PATH)
    make_figures(e1, e2, wins, e4, e6, p4_bins)
    write_report(gates, e1, e2, e3, wins, e4, e5, e6, p4_bins, verdicts)
    print(f"wrote {DOC_PATH}")
    print(json.dumps(_json_scalar(verdicts), indent=2, sort_keys=True))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    return parser.parse_args()


def main() -> None:
    parse_args()
    analyze()


if __name__ == "__main__":
    main()
