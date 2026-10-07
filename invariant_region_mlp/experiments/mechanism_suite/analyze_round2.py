"""Aggregate, judge, visualize, and report mechanism-suite round 2."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

from .mechanism_mlp import load_repetition
from .run_round2 import (
    CONTAINMENT_PATH,
    FULL_HISTORY_SOURCE,
    ILLEGAL_R2_FACTORS,
    R1_DEEP_CONFIGS,
    R1_SHALLOW_CONFIGS,
    RESULTS_ROOT,
    ROUND2_BASE_SEED,
    ROUND2_PREREGISTRATION_COMMIT,
    _source_candidates,
    r1_tasks,
    r2_tasks,
    r3_tasks,
    r4_tasks,
    r5_tasks,
    r6_tasks,
    task_path,
)


FIGURE_ROOT = RESULTS_ROOT / "figures"
DOC_PATH = RESULTS_ROOT.parent.parent / "docs" / "MECHANISM_SUITE_R2_REPORT.md"
SUMMARY_PATH = RESULTS_ROOT / "analysis_summary.json"
TUNING_PATH = RESULTS_ROOT / "r2_tuning_choices.json"
ROUND1_RESULTS_ROOT = RESULTS_ROOT.parent / "mechanism_suite"
ROUND1_REPORT_PATH = DOC_PATH.parent / "MECHANISM_SUITE_REPORT.md"


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
            return "Infinity" if number > 0.0 else "-Infinity"
        return number
    return value


def _atomic_json(payload: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(
            _json_scalar(payload), indent=2, sort_keys=True, allow_nan=False
        )
        + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    frame.to_csv(temporary, index=False)
    temporary.replace(path)


def _finite_skill(
    prediction: np.ndarray, truth: np.ndarray, mask: np.ndarray
) -> float:
    pred = np.asarray(prediction, dtype=np.float64)[mask]
    ref = np.asarray(truth, dtype=np.float64)[mask]
    if not np.all(np.isfinite(pred)):
        return float("inf")
    reference_std = float(np.std(ref, ddof=0))
    return (
        float(np.sqrt(np.mean((pred - ref) ** 2)) / reference_std)
        if reference_std > 0.0
        else float("inf")
    )


def _reuse_source_stage(task: dict[str, Any], path: Path) -> str | None:
    for source in _source_candidates(task):
        source_path = task_path(source)
        if not source_path.exists():
            continue
        try:
            if os.path.samefile(path, source_path):
                return str(source["stage"])
        except OSError:
            continue
    return None


def _flatten_artifact(
    path: Path, task: dict[str, Any]
) -> dict[str, Any]:
    path = path.resolve()
    loaded = load_repetition(path)
    metadata = loaded["metadata"]
    method = metadata["method"]
    work = metadata["work"]
    test = metadata["metrics"]["test"]
    validation = metadata["metrics"]["validation"]
    generator = work["generator"]
    noise = metadata["extra_diagnostics"]["child_u_z_noise"]
    nonfinite = int(test["nonfinite_state_count"]) + int(
        work["nonfinite_generator_values"]
    )
    return {
        "artifact": path.relative_to(RESULTS_ROOT).as_posix(),
        "stage": str(task["stage"]),
        "reused_from_stage": _reuse_source_stage(task, path),
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
        "recursively_evaluated_states": int(
            work["recursively_evaluated_states"]
        ),
        "stochastic_samples": int(work["total_stochastic_samples"]),
        "standard_normal_variates": int(work["standard_normal_variates"]),
        "dose_normal_variates": int(
            metadata["extra_diagnostics"].get("dose_normal_variates", 0)
        ),
        "wall_time_seconds": float(metadata["wall_clock_seconds"]),
        "pre_correction_violation_rate": float(
            work["pre_box_violation_rate"]
        ),
        "activation_rate": float(work["projection_activation_rate"]),
        "mean_method_overshoot_energy": float(
            work["mean_method_overshoot_energy"]
        ),
        "nonfinite_count": nonfinite,
        "nonfinite_state_count": int(test["nonfinite_state_count"]),
        "nonfinite_generator_count": int(
            work["nonfinite_generator_values"]
        ),
        "generator_bias": (
            float(generator["bias"]) if generator is not None else np.nan
        ),
        "generator_mae": (
            float(generator["mae"]) if generator is not None else np.nan
        ),
        "generator_mse": (
            float(generator["mse"]) if generator is not None else np.nan
        ),
        "child_u_z_noise_correlation": noise["correlation"],
        "child_noise_count": int(noise["count"]),
        "generator_bias_by_s_bin_json": json.dumps(
            _json_scalar(
                metadata["extra_diagnostics"]["generator_bias_by_s_bin"]
            ),
            sort_keys=True,
            allow_nan=False,
        ),
    }


def load_stage(stage: str, tasks: list[dict[str, Any]]) -> pd.DataFrame:
    expected = [(task_path(task), task) for task in tasks]
    missing = [path for path, _ in expected if not path.exists()]
    if missing:
        examples = "\n".join(str(path) for path in missing[:5])
        raise RuntimeError(
            f"{stage} incomplete: {len(missing)}/{len(expected)} missing. "
            f"First paths:\n{examples}"
        )
    frame = pd.DataFrame(
        _flatten_artifact(path, task) for path, task in expected
    )
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
            predictions.append(
                np.asarray(loaded["prediction_u"], dtype=np.float64)
            )
            if truth is None:
                truth = np.asarray(loaded["truth_u"], dtype=np.float64)
                validation = np.asarray(
                    loaded["is_validation"], dtype=np.bool_
                )
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
        sampling_variation = (
            float(
                np.sqrt(
                    np.mean((stack[:, test] - mean_prediction[test]) ** 2)
                )
                / ref_std
            )
            if np.all(np.isfinite(stack)) and ref_std > 0.0
            else float("inf")
        )
        summaries[key] = {
            "mean_rep_skill": float(np.mean(group["test_skill"])),
            "median_rep_skill": float(np.median(group["test_skill"])),
            "repetition_averaged_skill": _finite_skill(
                mean_prediction, truth, test
            ),
            "sampling_variation_skill": sampling_variation,
        }
    result = frame.copy()
    for column in next(iter(summaries.values())):
        result[column] = [
            summaries[tuple(row[key] for key in keys)][column]
            for _, row in result.iterrows()
        ]
    return result


def _ratio(numerator: float, denominator: float) -> float:
    if denominator == 0.0:
        return float("inf") if numerator > 0.0 else float("nan")
    return numerator / denominator


def _primary_name(pde: str) -> str:
    return "tight_segment" if pde == "P4" else "box"


def analyze_r1(
    frame: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    cell_records: list[dict[str, Any]] = []
    for cell, group in frame.groupby(["pde", "d", "n", "M"], sort=False):
        means = group.groupby("method")["test_skill"].mean().to_dict()
        primary = _primary_name(str(cell[0]))
        raw = float(means["raw"])
        method = float(means[primary])
        oracle = float(means["oracle_state"])
        gap = raw - oracle
        gc = (raw - method) / gap if gap != 0.0 else float("nan")
        cell_records.append(
            {
                "pde": cell[0],
                "d": int(cell[1]),
                "n": int(cell[2]),
                "M": int(cell[3]),
                "stratum": "deep" if int(cell[2]) >= 4 else "shallow",
                "primary_method": primary,
                "raw_mean_skill": raw,
                "primary_mean_skill": method,
                "oracle_state_mean_skill": oracle,
                "Gc": gc,
                "raw_over_oracle_state": _ratio(raw, oracle),
                "log_raw_over_oracle_state": math.log(_ratio(raw, oracle)),
            }
        )
    cells = pd.DataFrame(cell_records)
    result = frame.merge(
        cells,
        on=["pde", "d", "n", "M"],
        how="left",
        validate="many_to_one",
    )

    by_pde: dict[str, Any] = {}
    for pde, group in cells.groupby("pde", sort=True):
        deep = group[group["stratum"] == "deep"]
        shallow = group[group["stratum"] == "shallow"]
        deep_fraction = float(np.mean(deep["Gc"] >= 0.5))
        if pde in {"P2_a4", "P2_a8", "P3_rho1", "P3_rho2"}:
            shallow_rule = "Gc < 0.5"
            shallow_fraction = float(np.mean(shallow["Gc"] < 0.5))
            shallow_threshold = 0.5
        else:
            shallow_rule = "Gc >= 0.5"
            shallow_fraction = float(np.mean(shallow["Gc"] >= 0.5))
            shallow_threshold = 0.75
        by_pde[str(pde)] = {
            "C2_deep_fraction": deep_fraction,
            "C2_deep_passed": deep_fraction >= 0.75,
            "C2_shallow_rule": shallow_rule,
            "C2_shallow_fraction": shallow_fraction,
            "C2_shallow_threshold": shallow_threshold,
            "C2_shallow_passed": shallow_fraction >= shallow_threshold,
            "deep_cell_count": len(deep),
            "shallow_cell_count": len(shallow),
        }

    p23 = cells[
        cells["pde"].isin(
            ["P2_a4", "P2_a8", "P3_rho1", "P3_rho2"]
        )
    ]
    correlation, pvalue = spearmanr(
        p23["Gc"].to_numpy(dtype=np.float64),
        p23["log_raw_over_oracle_state"].to_numpy(dtype=np.float64),
    )
    verdict = {
        "by_pde": by_pde,
        "C2_deep_all_pdes_passed": all(
            item["C2_deep_passed"] for item in by_pde.values()
        ),
        "C2_shallow_all_pdes_passed": all(
            item["C2_shallow_passed"] for item in by_pde.values()
        ),
        "C2_graded": {
            "pde_scope": "pooled P2_a4, P2_a8, P3_rho1, P3_rho2",
            "cell_count": len(p23),
            "spearman_correlation": float(correlation),
            "two_sided_pvalue_descriptive": float(pvalue),
            "threshold": 0.5,
            "passed": bool(correlation >= 0.5),
        },
    }
    return result, cells, verdict


def tune_r2_controls(
    frame: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, Any], pd.DataFrame, dict[str, Any]]:
    result = frame.copy()
    result["selected_by_validation"] = False
    choices: list[dict[str, Any]] = []
    choice_lookup: dict[tuple[Any, ...], str] = {}
    cell_keys = ["pde", "d", "n", "M"]
    for cell, group in result.groupby(cell_keys, sort=False):
        families = {
            "B": group[
                group["transform"].isin(["shrink", "z_zero", "f_zero"])
            ],
            "C": group[
                group["transform"].isin(["illegal", "tight_illegal"])
            ],
        }
        for comparator_class, candidates in families.items():
            validation = (
                candidates.groupby(["method", "factor"], as_index=False)[
                    "validation_skill"
                ]
                .mean()
                .sort_values(["validation_skill", "factor", "method"])
            )
            chosen = validation.iloc[0]
            mask = np.ones(len(result), dtype=bool)
            for column, value in zip(cell_keys, cell):
                mask &= result[column].to_numpy() == value
            mask &= result["method"].to_numpy() == chosen["method"]
            result.loc[mask, "selected_by_validation"] = True
            choice_lookup[(*cell, comparator_class)] = str(chosen["method"])
            choices.append(
                {
                    "pde": cell[0],
                    "d": int(cell[1]),
                    "n": int(cell[2]),
                    "M": int(cell[3]),
                    "comparator_class": comparator_class,
                    "selected_method": str(chosen["method"]),
                    "selected_factor": float(chosen["factor"]),
                    "mean_validation_skill": float(chosen["validation_skill"]),
                    "candidate_validation_skills": {
                        str(row.method): float(row.validation_skill)
                        for row in validation.itertuples(index=False)
                    },
                }
            )

    win_records: list[dict[str, Any]] = []
    for cell, group in result.groupby(cell_keys, sort=False):
        indexed = group.set_index(["method", "rep"])["test_skill"]
        primary_name = _primary_name(str(cell[0]))
        primary = indexed.loc[primary_name]
        comparators = {
            "A": "centre",
            "B": choice_lookup[(*cell, "B")],
            "C": choice_lookup[(*cell, "C")],
        }
        record: dict[str, Any] = {
            "pde": cell[0],
            "d": int(cell[1]),
            "n": int(cell[2]),
            "M": int(cell[3]),
            "primary_method": primary_name,
            "primary_mean_skill": float(primary.mean()),
        }
        for comparator_class, method in comparators.items():
            comparison = indexed.loc[method]
            common = primary.index.intersection(comparison.index)
            record[f"class_{comparator_class}_method"] = method
            record[f"wins_vs_class_{comparator_class}"] = int(
                np.count_nonzero(
                    primary.loc[common].to_numpy()
                    < comparison.loc[common].to_numpy()
                )
            )
            record[f"mean_minus_class_{comparator_class}"] = float(
                primary.mean() - comparison.mean()
            )
        if cell[0] == "P4":
            segment = indexed.loc["segment"]
            illegal = indexed.loc[comparators["C"]]
            common = segment.index.intersection(illegal.index)
            record["segment_wins_vs_class_C"] = int(
                np.count_nonzero(
                    segment.loc[common].to_numpy()
                    < illegal.loc[common].to_numpy()
                )
            )
        win_records.append(record)
    wins = pd.DataFrame(win_records)

    verdict_by_pde: dict[str, Any] = {}
    for pde, group in wins.groupby("pde", sort=True):
        class_result: dict[str, Any] = {}
        for comparator_class in ("A", "B", "C"):
            fraction = float(
                np.mean(group[f"wins_vs_class_{comparator_class}"] >= 7)
            )
            class_result[comparator_class] = {
                "fraction_cells_with_at_least_7_wins": fraction,
                "threshold": 0.75,
                "passed": fraction >= 0.75,
            }
        centre_better_fraction = float(
            np.mean(group["mean_minus_class_A"] > 0.0)
        )
        entry: dict[str, Any] = {
            "classes": class_result,
            "certificate_prior_dominated": centre_better_fraction > 0.25,
            "centre_better_fraction": centre_better_fraction,
        }
        if pde == "P4":
            segment_fraction = float(
                np.mean(group["segment_wins_vs_class_C"] >= 7)
            )
            entry["round1_segment_vs_class_C"] = {
                "fraction_cells_with_at_least_7_wins": segment_fraction,
                "predicted_to_fail": True,
                "prediction_confirmed": segment_fraction < 0.75,
            }
        verdict_by_pde[str(pde)] = entry
    verdict = {
        "by_pde": verdict_by_pde,
        "all_primary_predictions_passed": all(
            all(item["passed"] for item in entry["classes"].values())
            for entry in verdict_by_pde.values()
        ),
        "P4_round1_segment_failure_prediction_confirmed": verdict_by_pde[
            "P4"
        ]["round1_segment_vs_class_C"]["prediction_confirmed"],
    }
    tuning = {
        "selection_rule": (
            "minimum mean validation skill within each (PDE,d,n,M) cell; "
            "test values were not consulted"
        ),
        "class_B_candidates": [
            "shrink_c0.1",
            "shrink_c0.25",
            "shrink_c0.5",
            "shrink_c0.75",
            "z_zero",
            "f_zero",
        ],
        "class_C_factors": list(ILLEGAL_R2_FACTORS),
        "choices": choices,
    }
    return result, tuning, wins, verdict


def analyze_r3(
    frame: pd.DataFrame, r2_verdict: dict[str, Any]
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    result = frame.copy()
    result["theta"] = np.where(
        result["method"] == "segment",
        0.0,
        np.where(
            result["method"] == "tight_segment", 1.0, result["factor"]
        ),
    )
    records: list[dict[str, Any]] = []
    for cell, group in result.groupby(["d", "n", "M"], sort=True):
        curve = group.groupby("theta")["test_skill"].mean().sort_index()
        monotone = bool(np.all(np.diff(curve.to_numpy()) <= 0.0))
        record: dict[str, Any] = {
            "d": int(cell[0]),
            "n": int(cell[1]),
            "M": int(cell[2]),
            "monotone_nonincreasing": monotone,
        }
        for theta, value in curve.items():
            record[f"mean_skill_theta_{float(theta):g}"] = float(value)
        records.append(record)
    cells = pd.DataFrame(records)
    fraction = float(np.mean(cells["monotone_nonincreasing"]))
    p4_classes = r2_verdict["by_pde"]["P4"]["classes"]
    tight2 = all(item["passed"] for item in p4_classes.values())
    verdict = {
        "C_tight_1": {
            "fraction_monotone_cells": fraction,
            "threshold": 0.75,
            "passed": fraction >= 0.75,
        },
        "C_tight_2": {
            "requires": "P4 tight_segment passes C3-A, C3-B, and C3-C",
            "passed": tight2,
        },
    }
    result = result.merge(
        cells[["d", "n", "M", "monotone_nonincreasing"]],
        on=["d", "n", "M"],
        how="left",
        validate="many_to_one",
    )
    return result, cells, verdict


def _bin_label(low: float, high: float) -> str:
    left = "-inf" if math.isinf(low) and low < 0.0 else f"{low:g}"
    right = "inf" if math.isinf(high) and high > 0.0 else f"{high:g}"
    return f"[{left},{right})"


def _expand_s_bins(
    frame: pd.DataFrame, *, source: str
) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for row in frame.itertuples(index=False):
        bins = json.loads(row.generator_bias_by_s_bin_json)
        for index, item in enumerate(bins):
            low = float(item["low"])
            high = float(item["high"])
            count = int(item["count"])
            bias = item["bias"]
            mae = item["mae"]
            z_count = int(item.get("z_error_count", 0))
            z_squared_sum = float(item.get("z_squared_error_sum", 0.0))
            records.append(
                {
                    "source": source,
                    "pde": row.pde,
                    "d": int(row.d),
                    "n": int(row.n),
                    "M": int(row.M),
                    "method": row.method,
                    "dose": (
                        float(row.factor) if source == "dose" else np.nan
                    ),
                    "rep": int(row.rep),
                    "bin_index": index,
                    "bin_low": low,
                    "bin_high": high,
                    "bin_label": _bin_label(low, high),
                    "is_central_bin": low == -0.25 and high == 0.25,
                    "generator_count": count,
                    "generator_bias": (
                        float(bias) if bias is not None else np.nan
                    ),
                    "generator_mae": (
                        float(mae) if mae is not None else np.nan
                    ),
                    "generator_error_sum": (
                        float(bias) * count if bias is not None else 0.0
                    ),
                    "generator_abs_error_sum": (
                        float(mae) * count if mae is not None else 0.0
                    ),
                    "z_error_count": z_count,
                    "z_squared_error_sum": z_squared_sum,
                    "z_error_rms": (
                        math.sqrt(z_squared_sum / z_count)
                        if z_count
                        else np.nan
                    ),
                }
            )
    return pd.DataFrame(records)


def analyze_r4(
    dose_frame: pd.DataFrame, r1_frame: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    dose_bins = _expand_s_bins(dose_frame, source="dose")
    raw = r1_frame[
        (r1_frame["pde"] == "P4")
        & (r1_frame["method"] == "raw")
        & (r1_frame["d"].isin([20, 50]))
        & (
            r1_frame[["n", "M"]]
            .apply(tuple, axis=1)
            .isin([(4, 3), (3, 6)])
        )
    ]
    raw_bins = _expand_s_bins(raw, source="raw_r1")
    expanded = pd.concat((dose_bins, raw_bins), ignore_index=True)

    dose_records: list[dict[str, Any]] = []
    for (dose, index, label, central), group in dose_bins.groupby(
        ["dose", "bin_index", "bin_label", "is_central_bin"], sort=True
    ):
        count = int(group["generator_count"].sum())
        dose_records.append(
            {
                "dose": float(dose),
                "bin_index": int(index),
                "bin_label": label,
                "is_central_bin": bool(central),
                "generator_count": count,
                "mean_absolute_dose_generator_error": (
                    float(group["generator_abs_error_sum"].sum() / count)
                    if count
                    else np.nan
                ),
                "mean_dose_generator_error": (
                    float(group["generator_error_sum"].sum() / count)
                    if count
                    else np.nan
                ),
            }
        )
    dose_summary = pd.DataFrame(dose_records)
    dose_verdicts: dict[str, Any] = {}
    for dose, group in dose_summary.groupby("dose", sort=True):
        central = float(
            group.loc[
                group["is_central_bin"],
                "mean_absolute_dose_generator_error",
            ].iloc[0]
        )
        outer = float(
            group.loc[
                ~group["is_central_bin"],
                "mean_absolute_dose_generator_error",
            ].max()
        )
        dose_verdicts[f"{dose:g}"] = {
            "central_mean_absolute_error": central,
            "largest_outer_mean_absolute_error": outer,
            "central_is_largest": central >= outer,
            "criterion_in_scope": dose >= 1.0,
        }

    raw_records: list[dict[str, Any]] = []
    for (index, label, central), group in raw_bins.groupby(
        ["bin_index", "bin_label", "is_central_bin"], sort=True
    ):
        count = int(group["generator_count"].sum())
        z_count = int(group["z_error_count"].sum())
        bias = (
            float(group["generator_error_sum"].sum() / count)
            if count
            else np.nan
        )
        z_rms = (
            math.sqrt(group["z_squared_error_sum"].sum() / z_count)
            if z_count
            else np.nan
        )
        raw_records.append(
            {
                "bin_index": int(index),
                "bin_label": label,
                "is_central_bin": bool(central),
                "generator_count": count,
                "generator_bias": bias,
                "z_error_rms": z_rms,
                "normalized_absolute_bias": abs(bias) / z_rms,
            }
        )
    raw_summary = pd.DataFrame(raw_records)
    central_normalized = float(
        raw_summary.loc[
            raw_summary["is_central_bin"], "normalized_absolute_bias"
        ].iloc[0]
    )
    largest_outer_normalized = float(
        raw_summary.loc[
            ~raw_summary["is_central_bin"], "normalized_absolute_bias"
        ].max()
    )
    verdict = {
        "C_loc_1": {
            "by_dose": dose_verdicts,
            "passed": all(
                item["central_is_largest"]
                for item in dose_verdicts.values()
                if item["criterion_in_scope"]
            ),
        },
        "C_loc_2": {
            "central_normalized_absolute_bias": central_normalized,
            "largest_outer_normalized_absolute_bias": largest_outer_normalized,
            "passed": central_normalized >= largest_outer_normalized,
        },
    }
    return expanded, dose_summary, raw_summary, verdict


def _bootstrap_interval(
    differences: np.ndarray, key: tuple[Any, ...], draws: int = 10_000
) -> tuple[float, float]:
    digest = hashlib.sha256(
        ("|".join(str(item) for item in key) + "|r5-bootstrap").encode(
            "utf-8"
        )
    ).digest()
    seed = int.from_bytes(digest[:8], "little")
    rng = np.random.default_rng(
        np.random.SeedSequence([ROUND2_BASE_SEED, seed])
    )
    indices = rng.integers(
        0, len(differences), size=(draws, len(differences))
    )
    means = np.mean(differences[indices], axis=1)
    return float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))


def analyze_r5(
    soft_frame: pd.DataFrame, r1_frame: pd.DataFrame
) -> tuple[pd.DataFrame, dict[str, Any]]:
    hard = r1_frame[
        r1_frame["pde"].isin(
            ["P2_a4", "P2_a8", "P3_rho1", "P3_rho2"]
        )
        & (r1_frame["method"] == "box")
    ]
    records: list[dict[str, Any]] = []
    for cell, candidates in soft_frame.groupby(
        ["pde", "d", "n", "M"], sort=False
    ):
        validation = (
            candidates.groupby(["method", "factor"], as_index=False)[
                "validation_skill"
            ]
            .mean()
            .sort_values(["validation_skill", "factor", "method"])
        )
        selected = validation.iloc[0]
        soft = candidates[candidates["method"] == selected["method"]].set_index(
            "rep"
        )
        hard_cell = hard[
            (hard["pde"] == cell[0])
            & (hard["d"] == cell[1])
            & (hard["n"] == cell[2])
            & (hard["M"] == cell[3])
        ].set_index("rep")
        common = soft.index.intersection(hard_cell.index)
        differences = (
            soft.loc[common, "test_skill"].to_numpy(dtype=np.float64)
            - hard_cell.loc[common, "test_skill"].to_numpy(dtype=np.float64)
        )
        low, high = _bootstrap_interval(differences, cell)
        mean_difference = float(np.mean(differences))
        records.append(
            {
                "pde": cell[0],
                "d": int(cell[1]),
                "n": int(cell[2]),
                "M": int(cell[3]),
                "stratum": "deep" if int(cell[2]) >= 4 else "shallow",
                "selected_soft_method": str(selected["method"]),
                "selected_factor": float(selected["factor"]),
                "mean_validation_skill": float(selected["validation_skill"]),
                "hard_mean_test_skill": float(
                    hard_cell.loc[common, "test_skill"].mean()
                ),
                "soft_mean_test_skill": float(
                    soft.loc[common, "test_skill"].mean()
                ),
                "soft_minus_hard_mean_skill": mean_difference,
                "paired_bootstrap_95_low": low,
                "paired_bootstrap_95_high": high,
                "paired_repetitions": len(common),
                "directional_winner": (
                    "soft" if mean_difference < 0.0 else "hard"
                ),
            }
        )
    result = pd.DataFrame(records)
    shallow = result[result["stratum"] == "shallow"]
    deep = result[result["stratum"] == "deep"]
    shallow_fraction = float(
        np.mean(shallow["soft_minus_hard_mean_skill"] < 0.0)
    )
    deep_fraction = float(
        np.mean(deep["soft_minus_hard_mean_skill"] > 0.0)
    )
    by_pde: dict[str, Any] = {}
    for pde, group in result.groupby("pde", sort=True):
        local_shallow = group[group["stratum"] == "shallow"]
        local_deep = group[group["stratum"] == "deep"]
        by_pde[str(pde)] = {
            "soft_better_shallow_fraction": float(
                np.mean(local_shallow["soft_minus_hard_mean_skill"] < 0.0)
            ),
            "hard_better_deep_fraction": float(
                np.mean(local_deep["soft_minus_hard_mean_skill"] > 0.0)
            ),
        }
    verdict = {
        "soft_better_shallow_fraction": shallow_fraction,
        "soft_better_shallow_threshold": 0.5,
        "soft_better_shallow_passed": shallow_fraction >= 0.5,
        "hard_better_deep_fraction": deep_fraction,
        "hard_better_deep_threshold": 0.75,
        "hard_better_deep_passed": deep_fraction >= 0.75,
        "directional_prediction_passed": (
            shallow_fraction >= 0.5 and deep_fraction >= 0.75
        ),
        "by_pde_descriptive": by_pde,
        "bootstrap_draws_per_cell": 10_000,
        "difference_definition": "selected soft test skill minus hard-box test skill",
    }
    return result, verdict


def analyze_r6(
    frame: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for d, group in frame.groupby("d", sort=True):
        means = group.groupby("method")["test_skill"].mean().to_dict()
        ratio = _ratio(float(means["box"]), float(means["oracle_state"]))
        records.append(
            {
                "d": int(d),
                "box_mean_skill": float(means["box"]),
                "oracle_state_mean_skill": float(means["oracle_state"]),
                "box_over_oracle_state": ratio,
            }
        )
    dimensions = pd.DataFrame(records)
    ratios = dimensions["box_over_oracle_state"].to_numpy(dtype=np.float64)
    max_over_min = float(np.max(ratios) / np.min(ratios))
    verdict = {
        "ratio_by_dimension": {
            str(int(row.d)): float(row.box_over_oracle_state)
            for row in dimensions.itertuples(index=False)
        },
        "max_over_min": max_over_min,
        "threshold": 1.25,
        "passed": max_over_min <= 1.25,
    }
    result = frame.merge(
        dimensions[["d", "box_over_oracle_state"]],
        on="d",
        how="left",
        validate="many_to_one",
    )
    return result, dimensions, verdict


def analyze_r7(
    r1_frame: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for cell, group in r1_frame.groupby(
        ["pde", "d", "n", "M"], sort=False
    ):
        indexed = group.set_index(["method", "rep"])
        box = indexed.loc["box"]
        oracle = indexed.loc["oracle_state"]
        common = box.index.intersection(oracle.index)
        for repetition in common:
            records.append(
                {
                    "pde": cell[0],
                    "d": int(cell[1]),
                    "n": int(cell[2]),
                    "M": int(cell[3]),
                    "stratum": "deep" if int(cell[2]) >= 4 else "shallow",
                    "rep": int(repetition),
                    "box_test_skill": float(
                        box.loc[repetition, "test_skill"]
                    ),
                    "oracle_state_test_skill": float(
                        oracle.loc[repetition, "test_skill"]
                    ),
                    "per_rep_skill_gap": float(
                        box.loc[repetition, "test_skill"]
                        - oracle.loc[repetition, "test_skill"]
                    ),
                    "box_repetition_averaged_skill": float(
                        box.loc[repetition, "repetition_averaged_skill"]
                    ),
                    "oracle_state_repetition_averaged_skill": float(
                        oracle.loc[
                            repetition, "repetition_averaged_skill"
                        ]
                    ),
                    "repetition_averaged_skill_gap": float(
                        box.loc[repetition, "repetition_averaged_skill"]
                        - oracle.loc[
                            repetition, "repetition_averaged_skill"
                        ]
                    ),
                    "box_generator_signed_error": float(
                        box.loc[repetition, "generator_bias"]
                    ),
                    "oracle_state_generator_signed_error": float(
                        oracle.loc[repetition, "generator_bias"]
                    ),
                }
            )
    result = pd.DataFrame(records)
    cell = result.drop_duplicates(["pde", "d", "n", "M"])
    by_pde_records: list[dict[str, Any]] = []
    for pde, group in result.groupby("pde", sort=True):
        local_cells = cell[cell["pde"] == pde]
        by_pde_records.append(
            {
                "pde": pde,
                "mean_per_rep_skill_gap": float(
                    group["per_rep_skill_gap"].mean()
                ),
                "mean_repetition_averaged_skill_gap": float(
                    local_cells["repetition_averaged_skill_gap"].mean()
                ),
                "mean_box_generator_signed_error": float(
                    group["box_generator_signed_error"].mean()
                ),
                "mean_oracle_generator_signed_error": float(
                    group["oracle_state_generator_signed_error"].mean()
                ),
            }
        )
    by_pde = pd.DataFrame(by_pde_records)
    per_rep_gap = float(result["per_rep_skill_gap"].mean())
    averaged_gap = float(cell["repetition_averaged_skill_gap"].mean())
    persistence = _ratio(averaged_gap, per_rep_gap)
    verdict = {
        "descriptive_only": True,
        "mean_per_rep_box_minus_oracle_skill": per_rep_gap,
        "mean_repetition_averaged_box_minus_oracle_skill": averaged_gap,
        "repetition_averaged_fraction_of_per_rep_gap": persistence,
        "mean_box_generator_signed_error": float(
            result["box_generator_signed_error"].mean()
        ),
        "mean_oracle_generator_signed_error": float(
            result["oracle_state_generator_signed_error"].mean()
        ),
        "interpretation": (
            "The repetition-averaged gap quantifies the component that "
            "persists after Monte Carlo averaging; compare it directly with "
            "the per-repetition gap and signed generator errors."
        ),
    }
    return result, by_pde, verdict


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


def plot_r1(cells: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    markers = {"deep": "o", "shallow": "s"}
    pdes = sorted(cells["pde"].unique())
    colors = plt.cm.tab10(np.linspace(0.0, 1.0, len(pdes)))
    for color, pde in zip(colors, pdes):
        subset = cells[cells["pde"] == pde]
        for stratum in ("deep", "shallow"):
            local = subset[subset["stratum"] == stratum]
            ax.scatter(
                local["log_raw_over_oracle_state"],
                local["Gc"],
                color=color,
                marker=markers[stratum],
                alpha=0.8,
                label=f"{pde}, {stratum}",
            )
    ax.axhline(0.5, color="black", lw=1.0, ls="--")
    ax.set_xlabel("log(raw / oracle_state mean skill)")
    ax.set_ylabel("Gc(primary)")
    ax.set_title("R1 regime hypothesis on fresh cells and configurations")
    ax.grid(alpha=0.2)
    ax.legend(ncol=3, fontsize=7)
    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / "r1_gc_regime.png")
    plt.close(fig)


def plot_r2(verdict: dict[str, Any]) -> None:
    pdes = sorted(verdict["by_pde"])
    classes = ["A", "B", "C"]
    values = np.array(
        [
            [
                verdict["by_pde"][pde]["classes"][name][
                    "fraction_cells_with_at_least_7_wins"
                ]
                for name in classes
            ]
            for pde in pdes
        ],
        dtype=np.float64,
    )
    fig, ax = plt.subplots(figsize=(5.5, 4.6))
    image = ax.imshow(values, vmin=0.0, vmax=1.0, cmap="viridis")
    ax.set_xticks(range(len(classes)), [f"C3-{name}" for name in classes])
    ax.set_yticks(range(len(pdes)), pdes)
    for row in range(len(pdes)):
        for column in range(len(classes)):
            ax.text(
                column,
                row,
                f"{values[row, column]:.3f}",
                ha="center",
                va="center",
                color="white" if values[row, column] < 0.6 else "black",
            )
    ax.set_title("R2 fraction of cells with at least 7/10 primary wins")
    fig.colorbar(image, ax=ax, label="cell fraction")
    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / "r2_comparator_classes.png")
    plt.close(fig)


def plot_r3(frame: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(8.5, 6.5), squeeze=False)
    for ax, (cell, group) in zip(
        axes.flat, frame.groupby(["d", "n", "M"], sort=True)
    ):
        summary = group.groupby("theta")["test_skill"].agg(["mean", "std"])
        ax.errorbar(
            summary.index,
            summary["mean"],
            yerr=summary["std"],
            marker="o",
            capsize=3,
        )
        ax.set_title(f"d={cell[0]}, (n,M)=({cell[1]},{cell[2]})")
        ax.set_xlabel("tightness theta")
        ax.set_ylabel("mean test skill")
        ax.grid(alpha=0.2)
    fig.suptitle("R3 valid-certificate tightness family")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.96))
    fig.savefig(FIGURE_ROOT / "r3_tightness.png")
    plt.close(fig)


def plot_r4(
    dose_summary: pd.DataFrame, raw_summary: pd.DataFrame
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    for dose, group in dose_summary.groupby("dose", sort=True):
        axes[0].plot(
            group["bin_index"],
            group["mean_absolute_dose_generator_error"],
            marker="o",
            label=f"dose={dose:g}",
        )
    labels = (
        dose_summary.sort_values("bin_index")
        .drop_duplicates("bin_index")["bin_label"]
        .tolist()
    )
    axes[0].set_xticks(range(len(labels)), labels, rotation=35, ha="right")
    axes[0].set_ylabel("mean |dose-induced generator error|")
    axes[0].set_title("Homoscedastic dose localization")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.2)

    axes[1].bar(
        raw_summary["bin_index"], raw_summary["normalized_absolute_bias"]
    )
    axes[1].set_xticks(range(len(labels)), labels, rotation=35, ha="right")
    axes[1].set_ylabel("|generator bias| / RMS(z error)")
    axes[1].set_title("Raw R1 normalized bias")
    axes[1].grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / "r4_localization.png")
    plt.close(fig)


def plot_r5(frame: pd.DataFrame) -> None:
    ordered = frame.sort_values(["stratum", "pde", "d", "n", "M"]).reset_index(
        drop=True
    )
    colors = np.where(ordered["stratum"] == "deep", "tab:blue", "tab:orange")
    lower = ordered["soft_minus_hard_mean_skill"] - ordered[
        "paired_bootstrap_95_low"
    ]
    upper = ordered["paired_bootstrap_95_high"] - ordered[
        "soft_minus_hard_mean_skill"
    ]
    fig, ax = plt.subplots(figsize=(12.0, 5.0))
    for index, row in ordered.iterrows():
        ax.errorbar(
            index,
            row["soft_minus_hard_mean_skill"],
            yerr=np.array([[lower.iloc[index]], [upper.iloc[index]]]),
            fmt="o",
            color=colors[index],
            capsize=2,
            markersize=4,
        )
    ax.axhline(0.0, color="black", lw=1.0)
    ax.set_xlabel("P2/P3 R1 cell (sorted by stratum and PDE)")
    ax.set_ylabel("selected soft skill - hard-box skill")
    ax.set_title("R5 paired mean differences and bootstrap 95% intervals")
    ax.grid(axis="y", alpha=0.2)
    ax.scatter([], [], color="tab:blue", label="deep")
    ax.scatter([], [], color="tab:orange", label="shallow")
    ax.legend(title="stratum")
    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / "r5_soft_vs_hard.png")
    plt.close(fig)


def plot_r6(dimensions: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.plot(
        dimensions["d"],
        dimensions["box_over_oracle_state"],
        marker="o",
    )
    minimum = float(dimensions["box_over_oracle_state"].min())
    ax.axhline(1.25 * minimum, color="black", ls="--", lw=1.0, label="1.25 x min")
    ax.set_xlabel("dimension d")
    ax.set_ylabel("mean box skill / mean oracle_state skill")
    ax.set_title("R6 pre-registered dimension stability")
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / "r6_dimension.png")
    plt.close(fig)


def plot_r7(by_pde: pd.DataFrame) -> None:
    x = np.arange(len(by_pde))
    width = 0.36
    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    ax.bar(
        x - width / 2,
        by_pde["mean_per_rep_skill_gap"],
        width,
        label="per repetition",
    )
    ax.bar(
        x + width / 2,
        by_pde["mean_repetition_averaged_skill_gap"],
        width,
        label="repetition averaged",
    )
    ax.set_xticks(x, by_pde["pde"], rotation=20)
    ax.set_ylabel("box skill - oracle_state skill")
    ax.set_title("R7 residual gap before and after repetition averaging")
    ax.axhline(0.0, color="black", lw=0.8)
    ax.legend()
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / "r7_gap.png")
    plt.close(fig)


def make_figures(
    r1_cells: pd.DataFrame,
    r2_verdict: dict[str, Any],
    r3: pd.DataFrame,
    dose_summary: pd.DataFrame,
    raw_summary: pd.DataFrame,
    r5: pd.DataFrame,
    r6_dimensions: pd.DataFrame,
    r7_by_pde: pd.DataFrame,
) -> None:
    FIGURE_ROOT.mkdir(parents=True, exist_ok=True)
    _style()
    plot_r1(r1_cells)
    plot_r2(r2_verdict)
    plot_r3(r3)
    plot_r4(dose_summary, raw_summary)
    plot_r5(r5)
    plot_r6(r6_dimensions)
    plot_r7(r7_by_pde)


def _format(value: Any) -> str:
    if value is None:
        return "—"
    if isinstance(value, (np.bool_, bool)):
        return "PASS" if bool(value) else "FAIL"
    if isinstance(value, (np.integer, int)):
        return str(int(value))
    if isinstance(value, (np.floating, float)):
        number = float(value)
        if math.isnan(number):
            return "—"
        if math.isinf(number):
            return "inf" if number > 0.0 else "-inf"
        magnitude = abs(number)
        if magnitude != 0.0 and (magnitude >= 1e4 or magnitude < 1e-3):
            return f"{number:.3e}"
        return f"{number:.4f}"
    return str(value).replace("|", "\\|")


def _markdown_table(
    frame: pd.DataFrame, columns: list[str] | None = None
) -> str:
    selected = frame if columns is None else frame[columns]
    headers = [str(column) for column in selected.columns]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in selected.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(_format(item) for item in row) + " |")
    return "\n".join(lines)


def _round1_outcome_table() -> str:
    text = ROUND1_REPORT_PATH.read_text(encoding="utf-8")
    lines = text.splitlines()
    start = next(
        index
        for index, line in enumerate(lines)
        if line.startswith("| PDE | C1 | C2 | C3 |")
    )
    end = start
    while end < len(lines) and lines[end].startswith("|"):
        end += 1
    return "\n".join(lines[start:end])


def _pass(value: bool) -> str:
    return "PASS" if value else "FAIL"


def _round2_outcome_table(
    r1_verdict: dict[str, Any], r2_verdict: dict[str, Any]
) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for pde in sorted(r1_verdict["by_pde"]):
        r1 = r1_verdict["by_pde"][pde]
        r2 = r2_verdict["by_pde"][pde]
        records.append(
            {
                "PDE": pde,
                "C2-deep": _pass(r1["C2_deep_passed"]),
                "C2-shallow": _pass(r1["C2_shallow_passed"]),
                "C3-A": _pass(r2["classes"]["A"]["passed"]),
                "C3-B": _pass(r2["classes"]["B"]["passed"]),
                "C3-C": _pass(r2["classes"]["C"]["passed"]),
                "class-A dominated": (
                    "yes" if r2["certificate_prior_dominated"] else "no"
                ),
            }
        )
    return pd.DataFrame(records)


def _criteria_rows(
    containment: dict[str, Any],
    r1: dict[str, Any],
    r2: dict[str, Any],
    r3: dict[str, Any],
    r4: dict[str, Any],
    r5: dict[str, Any],
    r6: dict[str, Any],
) -> list[dict[str, Any]]:
    all_c3 = all(
        all(item["passed"] for item in entry["classes"].values())
        for entry in r2["by_pde"].values()
    )
    return [
        {
            "id": "G5-R3",
            "criterion": (
                "R3 proceeds only if the violation decreases under refinement "
                "and the finest-grid violation is below 1e-5."
            ),
            "passed": bool(containment["proceed"]),
            "detail": (
                f"decreasing={containment['strictly_decreasing_under_refinement']}; "
                f"finest={containment['successive_grid_results'][-1]['maximum_violation']:.6g}"
            ),
        },
        {
            "id": "C2-deep",
            "criterion": (
                "C2-deep: primary has `Gc >= 0.5` in >= 75% of deep cells, "
                "for every PDE."
            ),
            "passed": bool(r1["C2_deep_all_pdes_passed"]),
            "detail": "; ".join(
                f"{pde}={entry['C2_deep_fraction']:.3f}"
                for pde, entry in r1["by_pde"].items()
            ),
        },
        {
            "id": "C2-shallow",
            "criterion": (
                "C2-shallow (boundary prediction): for P2/P3, primary has "
                "`Gc < 0.5` in >= 50% of shallow cells. For P1 and P4, "
                "`Gc >= 0.5` in >= 75% of shallow cells."
            ),
            "passed": bool(r1["C2_shallow_all_pdes_passed"]),
            "detail": "; ".join(
                f"{pde}={entry['C2_shallow_fraction']:.3f} ({entry['C2_shallow_rule']})"
                for pde, entry in r1["by_pde"].items()
            ),
        },
        {
            "id": "C2-graded",
            "criterion": (
                "C2-graded: pooled over all P2/P3 cells, Spearman correlation "
                "between `Gc` and `log(raw/oracle_state)` is >= 0.5."
            ),
            "passed": bool(r1["C2_graded"]["passed"]),
            "detail": f"rho={r1['C2_graded']['spearman_correlation']:.6g}",
        },
        {
            "id": "C3-A",
            "criterion": (
                "C3-A: primary wins >= 7/10 vs centre in >= 75% of cells."
            ),
            "passed": all(
                entry["classes"]["A"]["passed"]
                for entry in r2["by_pde"].values()
            ),
            "detail": "; ".join(
                f"{pde}={entry['classes']['A']['fraction_cells_with_at_least_7_wins']:.3f}"
                for pde, entry in r2["by_pde"].items()
            ),
        },
        {
            "id": "C3-B",
            "criterion": (
                "C3-B: primary wins >= 7/10 vs best class-B control in >= 75% of cells."
            ),
            "passed": all(
                entry["classes"]["B"]["passed"]
                for entry in r2["by_pde"].values()
            ),
            "detail": "; ".join(
                f"{pde}={entry['classes']['B']['fraction_cells_with_at_least_7_wins']:.3f}"
                for pde, entry in r2["by_pde"].items()
            ),
        },
        {
            "id": "C3-C",
            "criterion": (
                "C3-C: primary wins >= 7/10 vs best class-C control in >= 75% of cells."
            ),
            "passed": all(
                entry["classes"]["C"]["passed"]
                for entry in r2["by_pde"].values()
            ),
            "detail": "; ".join(
                f"{pde}={entry['classes']['C']['fraction_cells_with_at_least_7_wins']:.3f}"
                for pde, entry in r2["by_pde"].items()
            ),
        },
        {
            "id": "P4-round1-segment-C",
            "criterion": (
                "P4 with the round-1 `segment` is predicted to fail C "
                "(round-1 replication)."
            ),
            "passed": bool(
                r2["P4_round1_segment_failure_prediction_confirmed"]
            ),
            "detail": (
                "prediction confirmed"
                if r2["P4_round1_segment_failure_prediction_confirmed"]
                else "prediction not confirmed"
            ),
        },
        {
            "id": "C-tight-1",
            "criterion": (
                "C-tight-1: test skill decreases monotonically in theta "
                "(mean over reps) in >= 75% of the tightness-family cells."
            ),
            "passed": bool(r3["C_tight_1"]["passed"]),
            "detail": f"fraction={r3['C_tight_1']['fraction_monotone_cells']:.3f}",
        },
        {
            "id": "C-tight-2",
            "criterion": (
                "C-tight-2: tight_segment passes C3-A, C3-B and C3-C (R2 thresholds)."
            ),
            "passed": bool(r3["C_tight_2"]["passed"]),
            "detail": f"all P4 classes pass={r3['C_tight_2']['passed']}",
        },
        {
            "id": "C-loc-1",
            "criterion": (
                "C-loc-1: for every dose >= 1, the central bin "
                "`[-0.25, 0.25)` has the largest mean |dose-induced generator error|."
            ),
            "passed": bool(r4["C_loc_1"]["passed"]),
            "detail": "; ".join(
                f"dose {dose}: central={item['central_mean_absolute_error']:.4g}, "
                f"outer={item['largest_outer_mean_absolute_error']:.4g}"
                for dose, item in r4["C_loc_1"]["by_dose"].items()
                if item["criterion_in_scope"]
            ),
        },
        {
            "id": "C-loc-2",
            "criterion": (
                "C-loc-2: in raw E1 rows, the normalised bias is largest in the central bin."
            ),
            "passed": bool(r4["C_loc_2"]["passed"]),
            "detail": (
                f"central={r4['C_loc_2']['central_normalized_absolute_bias']:.4g}, "
                f"outer={r4['C_loc_2']['largest_outer_normalized_absolute_bias']:.4g}"
            ),
        },
        {
            "id": "R5-direction",
            "criterion": (
                "Directional prediction: soft has lower mean skill than hard in >= 50% "
                "of shallow cells, and hard has lower mean skill than soft in >= 75% "
                "of deep cells."
            ),
            "passed": bool(r5["directional_prediction_passed"]),
            "detail": (
                f"shallow={r5['soft_better_shallow_fraction']:.3f}; "
                f"deep={r5['hard_better_deep_fraction']:.3f}"
            ),
        },
        {
            "id": "C-dim",
            "criterion": (
                "C-dim: `max_d(box/oracle) / min_d(box/oracle) <= 1.25`."
            ),
            "passed": bool(r6["passed"]),
            "detail": f"max/min={r6['max_over_min']:.6g}",
        },
        {
            "id": "all-C3-predictions",
            "criterion": (
                "Predictions: P1–P3 pass A, B, C. P4 with `tight_segment` passes A, B, C."
            ),
            "passed": bool(all_c3),
            "detail": "joint summary of the separately evaluated C3 classes",
        },
    ]


def _work_accounting(frames: list[pd.DataFrame]) -> dict[str, Any]:
    combined = pd.concat(frames, ignore_index=True)
    computed = combined[combined["reused_from_stage"].isna()]
    return {
        "formal_rows": len(combined),
        "exact_reuse_rows": int(combined["reused_from_stage"].notna().sum()),
        "computed_rows": len(computed),
        "generator_calls": int(computed["f_calls"].sum()),
        "recursively_evaluated_states": int(
            computed["recursively_evaluated_states"].sum()
        ),
        "stochastic_samples": int(computed["stochastic_samples"].sum()),
        "dose_normal_variates": int(computed["dose_normal_variates"].sum()),
        "summed_worker_hours": float(computed["wall_time_seconds"].sum() / 3600.0),
        "nonfinite_rows": int(np.count_nonzero(combined["nonfinite_count"])),
    }


def _failure_lines(
    criteria: list[dict[str, Any]],
    r1: dict[str, Any],
    r2: dict[str, Any],
    r3_cells: pd.DataFrame,
) -> list[str]:
    failures: list[str] = []
    for item in criteria:
        if not item["passed"]:
            failures.append(f"{item['id']}: {item['detail']}")
    for pde, entry in r1["by_pde"].items():
        if not entry["C2_deep_passed"]:
            failures.append(
                f"{pde} C2-deep fraction={entry['C2_deep_fraction']:.4f}."
            )
        if not entry["C2_shallow_passed"]:
            failures.append(
                f"{pde} C2-shallow fraction={entry['C2_shallow_fraction']:.4f} "
                f"under rule {entry['C2_shallow_rule']}."
            )
    for pde, entry in r2["by_pde"].items():
        for comparator_class, result in entry["classes"].items():
            if not result["passed"]:
                failures.append(
                    f"{pde} C3-{comparator_class} fraction="
                    f"{result['fraction_cells_with_at_least_7_wins']:.4f}."
                )
    for row in r3_cells.itertuples(index=False):
        if not row.monotone_nonincreasing:
            failures.append(
                f"R3 tightness was not monotone for d={row.d}, "
                f"(n,M)=({row.n},{row.M})."
            )
    # Keep order while removing repeated global/per-PDE descriptions.
    return list(dict.fromkeys(failures))


def write_report(
    *,
    containment: dict[str, Any],
    r1_cells: pd.DataFrame,
    r1_verdict: dict[str, Any],
    r2_wins: pd.DataFrame,
    r2_verdict: dict[str, Any],
    r3_cells: pd.DataFrame,
    r3_verdict: dict[str, Any],
    r4_dose: pd.DataFrame,
    r4_raw: pd.DataFrame,
    r4_verdict: dict[str, Any],
    r5: pd.DataFrame,
    r5_verdict: dict[str, Any],
    r6_dimensions: pd.DataFrame,
    r6_verdict: dict[str, Any],
    r7_by_pde: pd.DataFrame,
    r7_verdict: dict[str, Any],
    criteria: list[dict[str, Any]],
    work: dict[str, Any],
) -> None:
    round2_table = _round2_outcome_table(r1_verdict, r2_verdict)
    failures = _failure_lines(criteria, r1_verdict, r2_verdict, r3_cells)
    full_history_hash = hashlib.sha256(FULL_HISTORY_SOURCE.read_bytes()).hexdigest()

    r1_summary_records: list[dict[str, Any]] = []
    for pde, item in r1_verdict["by_pde"].items():
        r1_summary_records.append(
            {
                "PDE": pde,
                "deep Gc>=0.5 fraction": item["C2_deep_fraction"],
                "deep verdict": _pass(item["C2_deep_passed"]),
                "shallow rule": item["C2_shallow_rule"],
                "shallow fraction": item["C2_shallow_fraction"],
                "shallow verdict": _pass(item["C2_shallow_passed"]),
            }
        )
    r1_summary = pd.DataFrame(r1_summary_records)

    r2_summary_records: list[dict[str, Any]] = []
    for pde, item in r2_verdict["by_pde"].items():
        r2_summary_records.append(
            {
                "PDE": pde,
                "C3-A fraction": item["classes"]["A"][
                    "fraction_cells_with_at_least_7_wins"
                ],
                "C3-B fraction": item["classes"]["B"][
                    "fraction_cells_with_at_least_7_wins"
                ],
                "C3-C fraction": item["classes"]["C"][
                    "fraction_cells_with_at_least_7_wins"
                ],
                "centre better fraction": item["centre_better_fraction"],
                "class-A dominated": (
                    "yes" if item["certificate_prior_dominated"] else "no"
                ),
            }
        )
    r2_summary = pd.DataFrame(r2_summary_records)

    pde_r5_records: list[dict[str, Any]] = []
    for pde, item in r5_verdict["by_pde_descriptive"].items():
        pde_r5_records.append(
            {
                "PDE": pde,
                "soft better shallow fraction": item[
                    "soft_better_shallow_fraction"
                ],
                "hard better deep fraction": item[
                    "hard_better_deep_fraction"
                ],
            }
        )
    r5_by_pde = pd.DataFrame(pde_r5_records)

    persistence = r7_verdict[
        "repetition_averaged_fraction_of_per_rep_gap"
    ]
    if persistence <= 1.0:
        r7_averaging_interpretation = (
            f"Repetition averaging retains {persistence:.3f} of the aggregate "
            "gap. The persistent part is consistent with systematic child-state/"
            "generator error, while the part that closes under averaging is "
            "consistent with Monte Carlo variation."
        )
    else:
        r7_averaging_interpretation = (
            f"The repetition-averaged gap is {persistence:.3f} times the mean "
            "per-repetition gap, so averaging does not close the gap; it enlarges "
            f"it by {(persistence - 1.0) * 100.0:.1f}%. Because normalised RMSE is "
            "nonlinear, these two gaps are not an additive bias--variance "
            "decomposition. Their persistence under averaging nevertheless supports "
            "a mainly systematic, bias-like residual rather than one dominated by "
            "Monte Carlo variation."
        )
    content: list[str] = [
        "# Mechanism suite, round 2: confirmatory re-run",
        "",
        "Round 2 used criteria written after the round-1 results were seen. "
        "It is therefore a fresh-randomness confirmatory study, not a replacement "
        "for the original pre-registered verdicts. Round-1 results and files were "
        "left unchanged.",
        "",
        "## Round-1 outcome (unchanged)",
        "",
        _round1_outcome_table(),
        "",
        "## Round-2 outcome",
        "",
        _markdown_table(round2_table),
        "",
        f"The pooled C2-graded Spearman correlation is "
        f"{r1_verdict['C2_graded']['spearman_correlation']:.4f} "
        f"({_pass(r1_verdict['C2_graded']['passed'])}). "
        f"C-tight-1={_pass(r3_verdict['C_tight_1']['passed'])}, "
        f"C-tight-2={_pass(r3_verdict['C_tight_2']['passed'])}, "
        f"C-loc-1={_pass(r4_verdict['C_loc_1']['passed'])}, "
        f"C-loc-2={_pass(r4_verdict['C_loc_2']['passed'])}, and "
        f"C-dim={_pass(r6_verdict['passed'])}.",
        "",
        "## Pre-registered criteria and mechanical verdicts",
        "",
    ]
    for item in criteria:
        content.extend(
            [
                f"### {item['id']}: {_pass(item['passed'])}",
                "",
                f"> {item['criterion']}",
                "",
                item["detail"],
                "",
            ]
        )
    content.extend(
        [
            "## All failures",
            "",
        ]
    )
    if failures:
        content.extend(f"- {failure}" for failure in failures)
    else:
        content.append("No pre-registered criterion or cell-level prediction failed.")

    content.extend(
        [
            "",
            "## Protocol and provenance",
            "",
            f"Base seed: `{ROUND2_BASE_SEED}`. The fixed 1,200-point data sets, "
            "including the 20% validation split, were regenerated from this seed. "
            "Tree seeds are `SeedSequence([20261107,d,n,M,rep,chunk_index])`; "
            "methods within a study use identical chunk sizes and paired trees. "
            "Every numerical array is float64. Value accuracy alone determines "
            "verdicts; gradient diagnostics are retained but do not determine a claim.",
            "",
            f"The round-2 pre-registration was frozen at commit "
            f"`{ROUND2_PREREGISTRATION_COMMIT}` before implementation or computation. "
            f"The containment artifact records code commit `{containment['code_commit']}` "
            f"on branch `{containment['git_branch']}`. The unchanged FullHistoryMLP "
            f"source has SHA-256 `{full_history_hash}`.",
            "",
            "The P1/P4 directions and all equation parameters are inherited unchanged "
            "from round 1. No round-1 result is reused as a round-2 random draw; only "
            "mathematically identical rows within round 2 are hard-linked and counted once "
            "in work accounting.",
            "",
            "## R3: tight P4 certificate and containment gate",
            "",
            "For `v=psi_s`, differentiation gives `v_t + v_ss + b v_s = 0` "
            "with `b=-lambda_f sign(v)` and `|b|<=lambda_f`. Monotonicity of "
            "`tanh` and constant-drift comparison give",
            "",
            "```text",
            "v_minus(t,s) <= psi_s(t,s) <= v_plus(t,s)",
            "v_pm(t,s) = E tanh(beta*(s +/- lambda_f*(T-t) + sqrt(2*(T-t))*xi))",
            "xi ~ N(0,1).",
            "```",
            "",
            "The bounds were computed with 80-node Gauss-Hermite quadrature. The "
            "recursive solver uses an audited cubic cache of those quadrature values; "
            "it clips exactly to the interpolated endpoints and applies no inward "
            "safety margin. The sharp-bound G5 audit used 100,000 fresh test-distribution "
            "points plus 10,000 points with `T-t in [0,0.1]` and "
            "`s in [-0.25,0.25]`.",
            "",
            "The first implementation audited only the first four levels of the "
            "already-existing six-level P4 reference hierarchy and stopped at "
            "`3.023e-5`. That failed audit is preserved under `audit_history/`. "
            "Before any round-2 MLP run, the implementation was corrected to use all "
            "six pre-existing levels; the `1e-5` threshold and every scientific "
            "criterion remained unchanged.",
            "",
            _markdown_table(
                pd.DataFrame(containment["successive_grid_results"])[
                    [
                        "n_space",
                        "n_steps",
                        "ds",
                        "maximum_violation",
                        "positive_violation_count",
                    ]
                ]
            ),
            "",
            f"The accepted Richardson reference has maximum bound violation "
            f"{containment['accepted_reference_check']['maximum_violation']:.3e}. "
            f"The quadrature-cache interpolation audit maximum is "
            f"{containment['bound_cache']['audit_maximum_absolute_error']:.3e}.",
            "",
            "### Tightness family",
            "",
            _markdown_table(r3_cells),
            "",
            "[Per-repetition R3 rows](../results/mechanism_suite_r2/r3_tightness.csv) "
            "· [tightness figure](../results/mechanism_suite_r2/figures/r3_tightness.png)",
            "",
            "## R1: C2 regime hypothesis",
            "",
            _markdown_table(r1_summary),
            "",
            f"Pooled P2/P3 Spearman rho="
            f"{r1_verdict['C2_graded']['spearman_correlation']:.4f}; "
            f"the threshold was 0.5.",
            "",
            "[All R1 rows](../results/mechanism_suite_r2/r1_regime.csv) "
            "· [regime figure](../results/mechanism_suite_r2/figures/r1_gc_regime.png)",
            "",
            "## R2: corrected C3 comparator classes",
            "",
            _markdown_table(r2_summary),
            "",
            "Class B and C choices were selected only by mean validation skill and "
            "then frozen for the test split. The full choices, including every "
            "candidate validation score, are in `r2_tuning_choices.json`.",
            "",
            "[All R2 rows](../results/mechanism_suite_r2/r2_controls.csv) "
            "· [class heatmap](../results/mechanism_suite_r2/figures/r2_comparator_classes.png)",
            "",
            "## R4: homoscedastic rectification localization",
            "",
            "### Dose-induced generator error",
            "",
            _markdown_table(r4_dose),
            "",
            "### Raw R1 normalized bias",
            "",
            _markdown_table(r4_raw),
            "",
            "[Expanded bin rows](../results/mechanism_suite_r2/r4_localization.csv) "
            "· [localization figure](../results/mechanism_suite_r2/figures/r4_localization.png)",
            "",
            "## R5: hard projection versus soft contraction",
            "",
            f"Across shallow cells, soft contraction has lower mean skill in "
            f"{r5_verdict['soft_better_shallow_fraction']:.4f} of cells (threshold "
            f"0.5). Across deep cells, hard projection has lower mean skill in "
            f"{r5_verdict['hard_better_deep_fraction']:.4f} of cells (threshold 0.75).",
            "",
            _markdown_table(r5_by_pde),
            "",
            "Each cell's mean paired difference and deterministic 10,000-draw paired "
            "bootstrap interval is in the CSV. This study is descriptive with a "
            "pre-registered direction and does not determine a C-claim.",
            "",
            "[R5 cell rows](../results/mechanism_suite_r2/r5_soft_vs_hard.csv) "
            "· [paired intervals](../results/mechanism_suite_r2/figures/r5_soft_vs_hard.png)",
            "",
            "## R6: dimension stability",
            "",
            _markdown_table(r6_dimensions),
            "",
            f"`max/min={r6_verdict['max_over_min']:.4f}` against the pre-registered "
            f"threshold 1.25: **{_pass(r6_verdict['passed'])}**.",
            "",
            "[All R6 rows](../results/mechanism_suite_r2/r6_dimension.csv) "
            "· [dimension figure](../results/mechanism_suite_r2/figures/r6_dimension.png)",
            "",
            "## R7: what remains between box and oracle",
            "",
            _markdown_table(r7_by_pde),
            "",
            f"The mean per-repetition box-minus-oracle skill gap is "
            f"{r7_verdict['mean_per_rep_box_minus_oracle_skill']:.4g}; the mean "
            f"gap after averaging the ten predictions is "
            f"{r7_verdict['mean_repetition_averaged_box_minus_oracle_skill']:.4g}. "
            f"{r7_averaging_interpretation} The mean signed "
            f"child-generator errors are {r7_verdict['mean_box_generator_signed_error']:.4g} "
            f"for box and {r7_verdict['mean_oracle_generator_signed_error']:.4g} for "
            "oracle_state. These numbers answer the question descriptively without "
            "introducing a post-hoc pass threshold.",
            "",
            "[Paired R7 rows](../results/mechanism_suite_r2/r7_gap.csv) "
            "· [gap figure](../results/mechanism_suite_r2/figures/r7_gap.png)",
            "",
            "## Work and integrity accounting",
            "",
            f"The six computational outputs contain {work['formal_rows']:,} formal "
            f"method-repetition rows, including {work['exact_reuse_rows']:,} exact "
            f"within-round hard-link reuses and {work['computed_rows']:,} newly "
            f"computed rows. Non-duplicated work totals {work['generator_calls']:,} "
            f"generator calls, {work['recursively_evaluated_states']:,} recursive "
            f"states, {work['stochastic_samples']:,} stochastic samples, and "
            f"{work['summed_worker_hours']:.2f} summed worker-hours. Rows with any "
            f"non-finite state/generator count: {work['nonfinite_rows']}.",
            "",
            "The manuscript, FullHistoryMLP source, and all round-1 results were not modified.",
        ]
    )
    DOC_PATH.parent.mkdir(parents=True, exist_ok=True)
    temporary = DOC_PATH.with_name(DOC_PATH.name + ".tmp")
    temporary.write_text("\n".join(content) + "\n", encoding="utf-8")
    temporary.replace(DOC_PATH)


def write_containment_failure_report(containment: dict[str, Any]) -> None:
    grids = pd.DataFrame(containment["successive_grid_results"])
    content = [
        "# Mechanism suite, round 2: stopped at R3 containment gate",
        "",
        "Round-2 criteria were written after round 1 was seen. Round-1 verdicts "
        "remain unchanged:",
        "",
        _round1_outcome_table(),
        "",
        "## Round-2 outcome",
        "",
        "R3 was stopped before any confirmatory P4 MLP computation because its "
        "sharp-bound numerical gate did not pass. All downstream round-2 criteria "
        "are NOT RUN, not failures.",
        "",
        "## G5-R3: FAIL",
        "",
        "> R3 proceeds only if the violation decreases under refinement and the "
        "finest-grid violation is below 1e-5.",
        "",
        _markdown_table(grids),
        "",
        f"decreasing={containment['strictly_decreasing_under_refinement']}; "
        f"finest_below_1e-5={containment['finest_grid_below_threshold']}.",
        "",
        "No manuscript or round-1 result was modified.",
    ]
    DOC_PATH.parent.mkdir(parents=True, exist_ok=True)
    DOC_PATH.write_text("\n".join(content) + "\n", encoding="utf-8")


def analyze() -> dict[str, Any]:
    if not CONTAINMENT_PATH.exists():
        raise RuntimeError("r3_containment.json is missing; run R3 first")
    containment = json.loads(CONTAINMENT_PATH.read_text(encoding="utf-8"))
    if not containment["proceed"]:
        summary = {
            "round": 2,
            "stopped_at": "R3 containment gate",
            "containment": containment,
        }
        _atomic_json(summary, SUMMARY_PATH)
        write_containment_failure_report(containment)
        return summary

    r1 = load_stage("r1", r1_tasks())
    r2 = load_stage("r2", r2_tasks())
    r3 = load_stage("r3", r3_tasks())
    r4 = load_stage("r4", r4_tasks())
    r5_soft = load_stage("r5", r5_tasks())
    r6 = load_stage("r6", r6_tasks())

    r1_output, r1_cells, r1_verdict = analyze_r1(r1)
    r2_output, tuning, r2_wins, r2_verdict = tune_r2_controls(r2)
    r2_output = r2_output.merge(
        r2_wins,
        on=["pde", "d", "n", "M"],
        how="left",
        validate="many_to_one",
    )
    r3_output, r3_cells, r3_verdict = analyze_r3(r3, r2_verdict)
    r4_output, r4_dose, r4_raw, r4_verdict = analyze_r4(r4, r1)
    r5_output, r5_verdict = analyze_r5(r5_soft, r1)
    r6_output, r6_dimensions, r6_verdict = analyze_r6(r6)
    r7_output, r7_by_pde, r7_verdict = analyze_r7(r1)

    criteria = _criteria_rows(
        containment,
        r1_verdict,
        r2_verdict,
        r3_verdict,
        r4_verdict,
        r5_verdict,
        r6_verdict,
    )
    work = _work_accounting([r1, r2, r3, r4, r5_soft, r6])
    summary = {
        "round": 2,
        "base_seed": ROUND2_BASE_SEED,
        "round1_verdicts_preserved": True,
        "criteria_written_after_round1": True,
        "containment": containment,
        "r1": r1_verdict,
        "r1_cells": r1_cells.to_dict(orient="records"),
        "r2": r2_verdict,
        "r2_win_counts": r2_wins.to_dict(orient="records"),
        "r3": r3_verdict,
        "r3_cells": r3_cells.to_dict(orient="records"),
        "r4": r4_verdict,
        "r5": r5_verdict,
        "r6": r6_verdict,
        "r7": r7_verdict,
        "criteria": criteria,
        "all_failures": _failure_lines(
            criteria, r1_verdict, r2_verdict, r3_cells
        ),
        "row_counts": {
            "r1_regime": len(r1_output),
            "r2_controls": len(r2_output),
            "r3_tightness": len(r3_output),
            "r4_localization_bins": len(r4_output),
            "r5_soft_vs_hard_cells": len(r5_output),
            "r6_dimension": len(r6_output),
            "r7_gap": len(r7_output),
        },
        "work": work,
    }

    _atomic_csv(r1_output, RESULTS_ROOT / "r1_regime.csv")
    _atomic_csv(r2_output, RESULTS_ROOT / "r2_controls.csv")
    _atomic_json(tuning, TUNING_PATH)
    _atomic_csv(r3_output, RESULTS_ROOT / "r3_tightness.csv")
    _atomic_csv(r4_output, RESULTS_ROOT / "r4_localization.csv")
    _atomic_csv(r5_output, RESULTS_ROOT / "r5_soft_vs_hard.csv")
    _atomic_csv(r6_output, RESULTS_ROOT / "r6_dimension.csv")
    _atomic_csv(r7_output, RESULTS_ROOT / "r7_gap.csv")
    _atomic_json(summary, SUMMARY_PATH)
    make_figures(
        r1_cells,
        r2_verdict,
        r3_output,
        r4_dose,
        r4_raw,
        r5_output,
        r6_dimensions,
        r7_by_pde,
    )
    write_report(
        containment=containment,
        r1_cells=r1_cells,
        r1_verdict=r1_verdict,
        r2_wins=r2_wins,
        r2_verdict=r2_verdict,
        r3_cells=r3_cells,
        r3_verdict=r3_verdict,
        r4_dose=r4_dose,
        r4_raw=r4_raw,
        r4_verdict=r4_verdict,
        r5=r5_output,
        r5_verdict=r5_verdict,
        r6_dimensions=r6_dimensions,
        r6_verdict=r6_verdict,
        r7_by_pde=r7_by_pde,
        r7_verdict=r7_verdict,
        criteria=criteria,
        work=work,
    )
    print(f"wrote {DOC_PATH}")
    print(json.dumps(_json_scalar(summary), indent=2, allow_nan=False))
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    return parser.parse_args()


def main() -> None:
    parse_args()
    analyze()


if __name__ == "__main__":
    main()
