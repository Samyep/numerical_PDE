"""Aggregate rescue results, generate work plots, and write technical reports."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[1]
REPOSITORY_ROOT = PROJECT_ROOT.parent
RESULT_ROOT = PROJECT_ROOT / "results" / "rescue_high_budget"
FIGURE_ROOT = RESULT_ROOT / "figures"
DOC_ROOT = PROJECT_ROOT / "docs"
FUNDING_REFERENCE = 21.299

FUNDING_F1 = [
    (2, 10), (3, 8), (4, 3), (2, 12), (2, 16), (2, 24), (2, 32),
    (3, 10), (3, 12), (3, 16), (4, 4), (4, 5), (4, 6), (5, 2), (5, 3),
]
FUNDING_HIGH_PROBE = [(2, 48), (2, 64), (3, 20), (3, 24), (4, 8)]
FUNDING_F2 = [
    (2, 10), (2, 32), (2, 64), (3, 8), (3, 16), (3, 24),
    (4, 3), (4, 8), (5, 3), (3, 48),
]
FUNDING_ULTRA_PROBE = [(2, 96), (3, 32), (3, 40), (3, 48), (4, 10), (4, 12)]
FUNDING_DIAGNOSTIC = [(2, 10), (3, 16), (3, 24), (3, 48)]
HJB_H1 = [
    (2, 16), (2, 24), (2, 32), (2, 48), (2, 64),
    (3, 8), (3, 12), (3, 16), (3, 24),
    (4, 4), (4, 6), (4, 8),
]
HJB_H2 = [(2, 48), (2, 64), (2, 96)]

COLORS = {"raw": "#222222", "samplewise": "#0072B2", "z_zero": "#E69F00", "f_zero": "#CC79A7"}
LABELS = {"raw": "Raw", "samplewise": "Samplewise certified IR", "z_zero": "z=0", "f_zero": "f=0"}
MARKERS = {2: "o", 3: "s", 4: "^", 5: "D"}


def _metadata(data: Any) -> dict[str, Any]:
    return json.loads(str(data["metadata_json"]))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git(*args: str) -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(REPOSITORY_ROOT), *args], text=True
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def funding_repetitions(base_seed: int = 20261005) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    root = RESULT_ROOT / "raw" / "funding" / f"seed{base_seed}"
    for path in sorted(root.glob("n*_M*/*/block_*.npz")):
        with np.load(path, allow_pickle=False) as data:
            metadata = _metadata(data)
            prediction = data["prediction_state"][:, 0]
            root_terminal = data["root_terminal_state"][:, 0]
            correction = data["nonlinear_u_correction"]
        count = int(metadata["root_count"])
        work = metadata["work"]
        method = metadata["method"]["name"]
        for offset in range(count):
            value = float(prediction[offset])
            rows.append(
                {
                    "pde": "funding",
                    "n": metadata["n"],
                    "M": metadata["M"],
                    "method": method,
                    "repetition": metadata["root_start"] + offset,
                    "prediction": value,
                    "absolute_error": abs(value - FUNDING_REFERENCE),
                    "squared_error": (value - FUNDING_REFERENCE) ** 2,
                    "bias": value - FUNDING_REFERENCE,
                    "root_terminal_u": float(root_terminal[offset]),
                    "nonlinear_u_correction": float(correction[offset]),
                    "terminal_g_evals": work["terminal_g_evals"] / count,
                    "f_evals": work["f_evals"] / count,
                    "recursively_evaluated_states": work["recursively_evaluated_states"] / count,
                    "total_stochastic_samples": work["total_stochastic_samples"] / count,
                    "wall_clock_seconds": metadata["wall_clock_seconds"] / count,
                    "constraint_violation_rate": work["constraint_violation_rate"],
                    "projection_activation_rate": work["projection_activation_rate"],
                    "mean_overshoot_energy": work["mean_overshoot_energy"],
                    "nonfinite_states": work["nonfinite_states"],
                    "nonfinite_generators": work["nonfinite_generators"],
                    "draw_fingerprint": metadata["draw_fingerprint"],
                    "source_path": str(path.relative_to(RESULT_ROOT)).replace("\\", "/"),
                }
            )
    return pd.DataFrame(rows).sort_values(["n", "M", "method", "repetition"])


def summarize_funding(repetitions: pd.DataFrame) -> pd.DataFrame:
    rows = []
    raw_lookup = repetitions[repetitions.method == "raw"].set_index(
        ["n", "M", "repetition"]
    )
    for (n, M, method), frame in repetitions.groupby(["n", "M", "method"]):
        row: dict[str, Any] = {
            "n": n,
            "M": M,
            "method": method,
            "repetitions": len(frame),
            "prediction_mean": frame.prediction.mean(),
            "prediction_std": frame.prediction.std(ddof=1),
            "mae": frame.absolute_error.mean(),
            "mae_se": frame.absolute_error.std(ddof=1) / math.sqrt(len(frame)),
            "rmse": math.sqrt(frame.squared_error.mean()),
            "bias": frame.bias.mean(),
            "nonlinear_correction_variance": frame.nonlinear_u_correction.var(ddof=1),
            "terminal_g_evals": frame.terminal_g_evals.mean(),
            "f_evals": frame.f_evals.mean(),
            "recursively_evaluated_states": frame.recursively_evaluated_states.mean(),
            "total_stochastic_samples": frame.total_stochastic_samples.mean(),
            "wall_clock_seconds": frame.wall_clock_seconds.mean(),
            "constraint_violation_rate": frame.constraint_violation_rate.mean(),
            "projection_activation_rate": frame.projection_activation_rate.mean(),
            "mean_overshoot_energy": frame.mean_overshoot_energy.mean(),
            "paired_improvement_vs_raw": np.nan,
            "paired_win_fraction_vs_raw": np.nan,
        }
        if method != "raw":
            differences = []
            wins = []
            for item in frame.itertuples(index=False):
                key = (item.n, item.M, item.repetition)
                if key not in raw_lookup.index:
                    continue
                raw_error = float(raw_lookup.loc[key, "absolute_error"])
                differences.append(raw_error - item.absolute_error)
                wins.append(item.absolute_error < raw_error)
            if differences:
                row["paired_improvement_vs_raw"] = float(np.mean(differences))
                row["paired_win_fraction_vs_raw"] = float(np.mean(wins))
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["n", "M", "method"])


def hjb_repetitions() -> tuple[pd.DataFrame, dict[tuple[Any, ...], list[np.ndarray]]]:
    rows: list[dict[str, Any]] = []
    corrections: dict[tuple[Any, ...], list[np.ndarray]] = {}
    root = RESULT_ROOT / "raw" / "hjb"
    for stage in ("h1", "h2", "h2_fzero"):
        stage_root = root / stage
        if not stage_root.exists():
            continue
        for path in sorted(stage_root.glob("*.npz")):
            with np.load(path, allow_pickle=False) as data:
                metadata = _metadata(data)
                correction = data["nonlinear_u_correction"].copy()
            analysis_stage = "h2" if stage == "h2_fzero" else stage
            work = metadata["work"]
            metrics = metadata["metrics"]
            generator = metadata["generator_diagnostic"] or {}
            key = (
                analysis_stage,
                metadata["dimension"],
                metadata["n"],
                metadata["M"],
                metadata["method"]["name"],
            )
            corrections.setdefault(key, []).append(correction)
            rows.append(
                {
                    "pde": "hjb",
                    "stage": analysis_stage,
                    "dimension": metadata["dimension"],
                    "n": metadata["n"],
                    "M": metadata["M"],
                    "method": metadata["method"]["name"],
                    "repetition": metadata["repetition"],
                    "n_points": metadata["n_points"],
                    "value_relative_l2": metrics["value_relative_l2"],
                    "value_mae": metrics["value_mae"],
                    "value_bias": metrics["value_bias"],
                    "gradient_relative_l2": metrics["gradient_relative_l2"],
                    "gradient_mae": metrics["gradient_mae"],
                    "generator_mse": generator.get("mse", np.nan),
                    "generator_bias": generator.get("bias", np.nan),
                    "generator_absolute_bias": generator.get("absolute_bias", np.nan),
                    "generator_mae": generator.get("mae", np.nan),
                    "generator_mean_truth": generator.get("mean_truth", np.nan),
                    "generator_mean_truth_squared": generator.get("mean_truth_squared", np.nan),
                    "terminal_g_evals": work["terminal_g_evals"],
                    "f_evals": work["f_evals"],
                    "recursively_evaluated_states": work["recursively_evaluated_states"],
                    "total_stochastic_samples": work["total_stochastic_samples"],
                    "wall_clock_seconds": metadata["wall_clock_seconds"],
                    "constraint_violation_rate": work["constraint_violation_rate"],
                    "projection_activation_rate": work["projection_activation_rate"],
                    "mean_overshoot_energy": work["mean_overshoot_energy"],
                    "nonfinite_states": work["nonfinite_states"],
                    "nonfinite_generators": work["nonfinite_generators"],
                    "draw_fingerprint": metadata["draw_fingerprint"],
                    "heuristic_clipping": metadata["heuristic_clipping"],
                    "source_path": str(path.relative_to(RESULT_ROOT)).replace("\\", "/"),
                }
            )
    return (
        pd.DataFrame(rows).sort_values(
            ["stage", "dimension", "n", "M", "method", "repetition"]
        ),
        corrections,
    )


def summarize_hjb(
    repetitions: pd.DataFrame,
    corrections: dict[tuple[Any, ...], list[np.ndarray]],
) -> pd.DataFrame:
    rows = []
    raw_lookup = repetitions[repetitions.method == "raw"].set_index(
        ["stage", "dimension", "n", "M", "repetition"]
    )
    for key, frame in repetitions.groupby(["stage", "dimension", "n", "M", "method"]):
        stage, dimension, n, M, method = key
        arrays = corrections[key]
        correction_variance = np.nan
        if len(arrays) > 1:
            correction_variance = float(
                np.mean(np.var(np.stack(arrays), axis=0, ddof=1))
            )
        row: dict[str, Any] = {
            "stage": stage,
            "dimension": dimension,
            "n": n,
            "M": M,
            "method": method,
            "repetitions": len(frame),
            "n_points": frame.n_points.iloc[0],
            "value_relative_l2_mean": frame.value_relative_l2.mean(),
            "value_relative_l2_std": frame.value_relative_l2.std(ddof=1),
            "value_mae_mean": frame.value_mae.mean(),
            "value_bias_mean": frame.value_bias.mean(),
            "gradient_relative_l2_mean": frame.gradient_relative_l2.mean(),
            "gradient_mae_mean": frame.gradient_mae.mean(),
            "generator_mse_mean": frame.generator_mse.mean(),
            "generator_bias_mean": frame.generator_bias.mean(),
            "generator_absolute_bias_mean": frame.generator_absolute_bias.mean(),
            "nonlinear_correction_variance": correction_variance,
            "terminal_g_evals": frame.terminal_g_evals.mean(),
            "f_evals": frame.f_evals.mean(),
            "recursively_evaluated_states": frame.recursively_evaluated_states.mean(),
            "total_stochastic_samples": frame.total_stochastic_samples.mean(),
            "wall_clock_seconds": frame.wall_clock_seconds.mean(),
            "constraint_violation_rate": frame.constraint_violation_rate.mean(),
            "projection_activation_rate": frame.projection_activation_rate.mean(),
            "mean_overshoot_energy": frame.mean_overshoot_energy.mean(),
            "paired_improvement_vs_raw": np.nan,
            "paired_win_fraction_vs_raw": np.nan,
        }
        if method != "raw":
            differences = []
            wins = []
            for item in frame.itertuples(index=False):
                raw_key = (stage, dimension, n, M, item.repetition)
                if raw_key not in raw_lookup.index:
                    continue
                raw_error = float(raw_lookup.loc[raw_key, "value_relative_l2"])
                differences.append(raw_error - item.value_relative_l2)
                wins.append(item.value_relative_l2 < raw_error)
            if differences:
                row["paired_improvement_vs_raw"] = float(np.mean(differences))
                row["paired_win_fraction_vs_raw"] = float(np.mean(wins))
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["stage", "dimension", "n", "M", "method"])


def _main_funding(summary: pd.DataFrame) -> pd.DataFrame:
    return summary[
        summary.method.isin(["raw", "samplewise"])
        & summary.apply(lambda row: (row.n, row.M) in set(FUNDING_F1 + FUNDING_HIGH_PROBE + FUNDING_ULTRA_PROBE), axis=1)
    ].copy()


def _line_by_depth(axis: Any, data: pd.DataFrame, x: str, y: str, method: str) -> None:
    frame = data[data.method == method]
    for n, group in frame.groupby("n"):
        group = group.sort_values(x)
        axis.plot(
            group[x], group[y], marker=MARKERS.get(int(n), "o"), markersize=5,
            linewidth=1.5 if method in ("raw", "samplewise") else 1.0,
            color=COLORS.get(method), linestyle="-" if method in ("raw", "samplewise") else "--",
            alpha=1.0 if method in ("raw", "samplewise") else 0.8,
            label=f"{LABELS.get(method, method)}, n={n}",
        )


def make_figures(funding: pd.DataFrame, hjb: pd.DataFrame) -> list[Path]:
    FIGURE_ROOT.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []

    main = _main_funding(funding)
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.6))
    for axis, x, label in zip(
        axes,
        ("f_evals", "total_stochastic_samples", "wall_clock_seconds"),
        ("f evaluations / root", "stochastic samples / root", "wall seconds / root"),
    ):
        for method in ("raw", "samplewise"):
            _line_by_depth(axis, main, x, "mae", method)
        diag = funding[(funding.method == "z_zero") & funding.apply(lambda row: (row.n, row.M) in FUNDING_DIAGNOSTIC, axis=1)]
        if not diag.empty:
            axis.scatter(diag[x], diag.mae, color=COLORS["z_zero"], marker="x", s=45, label="z=0 diagnostics")
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlabel(label)
        axis.set_ylabel("root MAE")
        axis.grid(True, which="both", alpha=0.25)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.925),
        ncol=5,
        frameon=False,
        fontsize=8,
    )
    fig.suptitle("Funding rescue: Raw vs Samplewise certified IR", y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.76))
    path = FIGURE_ROOT / "funding_error_vs_work.pdf"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    outputs.append(path)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6))
    paired = main[main.method == "samplewise"].sort_values("total_stochastic_samples")
    for n, frame in paired.groupby("n"):
        axes[0].scatter(
            frame.total_stochastic_samples,
            frame.paired_improvement_vs_raw,
            marker=MARKERS.get(int(n), "o"),
            s=55,
            label=f"n={n}",
        )
    axes[0].axhline(0.0, color="black", linewidth=0.9)
    axes[0].set_xscale("log")
    axes[0].set_xlabel("stochastic samples / root")
    axes[0].set_ylabel("paired MAE improvement (Raw - IR)")
    axes[0].grid(True, which="both", alpha=0.25)
    axes[0].legend(frameon=False)
    diagnostics = funding[
        funding.apply(lambda row: (row.n, row.M) in FUNDING_DIAGNOSTIC, axis=1)
        & funding.method.isin(["raw", "samplewise", "z_zero", "f_zero"])
    ].copy()
    diagnostics["config"] = diagnostics.apply(lambda row: f"({int(row.n)},{int(row.M)})", axis=1)
    configs = [f"({n},{M})" for n, M in FUNDING_DIAGNOSTIC]
    width = 0.24
    position = np.arange(len(configs))
    for offset, method in enumerate(("raw", "samplewise", "z_zero")):
        values = []
        for config in configs:
            match = diagnostics[(diagnostics.config == config) & (diagnostics.method == method)]
            values.append(float(match.mae.iloc[0]) if len(match) else np.nan)
        axes[1].bar(position + (offset - 1) * width, values, width, label=LABELS[method], color=COLORS[method])
    axes[1].set_yscale("log")
    axes[1].set_xticks(position, configs)
    axes[1].set_xlabel("(n,M)")
    axes[1].set_ylabel("root MAE")
    axes[1].grid(True, axis="y", which="both", alpha=0.25)
    axes[1].legend(frameon=False)
    fig.suptitle("Funding crossover and suppression-floor diagnostics")
    fig.tight_layout()
    path = FIGURE_ROOT / "funding_budget_crossover.pdf"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    outputs.append(path)

    h1 = hjb[(hjb.stage == "h1") & (hjb.dimension == 100)]
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.6))
    for axis, x, label in zip(
        axes,
        ("f_evals", "total_stochastic_samples", "wall_clock_seconds"),
        ("f evaluations / repetition", "stochastic samples / repetition", "wall seconds / repetition"),
    ):
        for method in ("raw", "samplewise", "f_zero"):
            frame = h1[h1.method == method]
            for n, group in frame.groupby("n"):
                group = group.sort_values(x)
                xvalues = group[x].copy()
                if method == "f_zero" and x == "f_evals":
                    xvalues = np.full(len(group), 0.5)
                axis.plot(
                    xvalues,
                    group.value_relative_l2_mean,
                    marker=MARKERS.get(int(n), "o"),
                    markersize=5,
                    linewidth=1.5 if method != "f_zero" else 1.0,
                    color=COLORS[method],
                    linestyle="-" if method != "f_zero" else "--",
                    label=f"{LABELS[method]}, n={n}",
                )
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlabel(label)
        axis.set_ylabel("value relative L2")
        axis.grid(True, which="both", alpha=0.25)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.925),
        ncol=5,
        frameon=False,
        fontsize=8,
    )
    fig.suptitle("HJB Stage H1: high work is not always high sampling accuracy", y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.76))
    path = FIGURE_ROOT / "hjb_error_vs_work.pdf"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    outputs.append(path)

    h2 = hjb[hjb.stage == "h2"]
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.4), sharey=True)
    for axis, dimension in zip(axes, (100, 140, 160)):
        frame = h2[h2.dimension == dimension]
        for method in ("raw", "samplewise"):
            group = frame[frame.method == method].sort_values("M")
            axis.plot(group.M, group.value_relative_l2_mean, marker="o", linewidth=1.8, color=COLORS[method], label=LABELS[method])
        zero = frame[frame.method == "f_zero"]
        if not zero.empty:
            axis.axhline(float(zero.value_relative_l2_mean.mean()), color=COLORS["f_zero"], linestyle="--", linewidth=1.4, label="f=0 floor")
        axis.set_yscale("log")
        axis.set_xlabel("M (n=2)")
        axis.set_title(f"d={dimension}")
        axis.grid(True, which="both", alpha=0.25)
    axes[0].set_ylabel("value relative L2")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.91),
        ncol=3,
        frameon=False,
    )
    fig.suptitle("HJB high-sampling probe: trending, but crossover not reached", y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.80))
    path = FIGURE_ROOT / "hjb_crossover_probe.pdf"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    outputs.append(path)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6))
    funding_ratio = paired.copy()
    funding_ratio["ratio"] = funding_ratio.apply(
        lambda row: row.mae
        / float(main[(main.n == row.n) & (main.M == row.M) & (main.method == "raw")].mae.iloc[0]),
        axis=1,
    )
    for n, frame in funding_ratio.groupby("n"):
        axes[0].scatter(frame.total_stochastic_samples, frame.ratio, marker=MARKERS.get(int(n), "o"), s=55, label=f"n={n}")
    axes[0].axhline(1.0, color="black", linewidth=0.9)
    axes[0].set_xscale("log")
    axes[0].set_xlabel("Funding stochastic samples / root")
    axes[0].set_ylabel("IR error / Raw error")
    axes[0].set_title("Funding")
    axes[0].grid(True, which="both", alpha=0.25)
    axes[0].legend(frameon=False)
    h2_ir = h2[h2.method == "samplewise"].copy()
    ratios = []
    for row in h2_ir.itertuples(index=False):
        raw = h2[(h2.dimension == row.dimension) & (h2.n == row.n) & (h2.M == row.M) & (h2.method == "raw")]
        ratios.append(row.value_relative_l2_mean / float(raw.value_relative_l2_mean.iloc[0]))
    h2_ir["ratio"] = ratios
    for dimension, frame in h2_ir.groupby("dimension"):
        axes[1].plot(frame.M, frame.ratio, marker="o", label=f"d={dimension}")
    axes[1].axhline(1.0, color="black", linewidth=0.9)
    axes[1].set_xlabel("HJB M (n=2)")
    axes[1].set_ylabel("IR error / Raw error")
    axes[1].set_title("Rosenbrock HJB")
    axes[1].grid(True, alpha=0.25)
    axes[1].legend(frameon=False)
    fig.suptitle("Primary comparison: Samplewise certified IR relative to Raw")
    fig.tight_layout()
    path = FIGURE_ROOT / "raw_vs_ir_summary.pdf"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    outputs.append(path)
    return outputs


def _fmt(value: Any, digits: int = 5) -> str:
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return "-"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    return f"{float(value):.{digits}g}"


def _table(headers: list[str], rows: Iterable[Iterable[Any]]) -> str:
    result = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    result.extend("| " + " | ".join(str(item) for item in row) + " |" for row in rows)
    return "\n".join(result)


def _manifest_wall_seconds() -> tuple[float, list[dict[str, Any]]]:
    rows = []
    total = 0.0
    for path in sorted(RESULT_ROOT.glob("*_manifest.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        elapsed = float(payload.get("elapsed_this_invocation_seconds", 0.0))
        rows.append({"manifest": path.name, "elapsed_seconds": elapsed, "status": payload.get("status")})
        total += elapsed
    legacy = RESULT_ROOT / "hjb_legacy_reproduction.json"
    if legacy.exists():
        payload = json.loads(legacy.read_text(encoding="utf-8"))
        elapsed = float(payload.get("elapsed_seconds", 0.0))
        rows.append({"manifest": legacy.name, "elapsed_seconds": elapsed, "status": "complete"})
        total += elapsed
    return total, rows


def write_reports(funding: pd.DataFrame, hjb: pd.DataFrame, compute: dict[str, Any]) -> None:
    funding_primary = funding[funding.method.isin(["raw", "samplewise"])]
    configs = sorted(set(zip(funding_primary.n, funding_primary.M)), key=lambda item: (item[0], item[1]))
    primary_rows = []
    for n, M in configs:
        raw = funding_primary[(funding_primary.n == n) & (funding_primary.M == M) & (funding_primary.method == "raw")]
        ir = funding_primary[(funding_primary.n == n) & (funding_primary.M == M) & (funding_primary.method == "samplewise")]
        if raw.empty or ir.empty:
            continue
        r, p = raw.iloc[0], ir.iloc[0]
        primary_rows.append([
            n, M, int(min(r.repetitions, p.repetitions)), _fmt(r.mae), _fmt(p.mae),
            _fmt(p.paired_improvement_vs_raw), _fmt(p.paired_win_fraction_vs_raw, 3),
            _fmt(p.f_evals, 7), _fmt(p.total_stochastic_samples, 7), _fmt(p.wall_clock_seconds, 4),
        ])
    diagnostics_rows = []
    for n, M in FUNDING_DIAGNOSTIC:
        for method in ("raw", "samplewise", "z_zero", "f_zero", "shrink_c0.25", "tight_a0.75"):
            match = funding[(funding.n == n) & (funding.M == M) & (funding.method == method)]
            if match.empty:
                continue
            row = match.iloc[0]
            diagnostics_rows.append([n, M, method, int(row.repetitions), _fmt(row.mae), _fmt(row.bias), _fmt(row.paired_win_fraction_vs_raw, 3)])

    funding_report = f"""# Funding high-budget rescue

## Outcome

**Funding verdict F-B - PARTIALLY RESCUED.** Samplewise certified IR improves Raw full-history MLP throughout the completed primary grid, including every 100-root final configuration. The improvement is therefore a real stabilization effect, not a single favorable low-budget point. However, the correct nonlinear-gradient estimator never beats the `z=0` wrong-gradient diagnostic at the representative low, medium, or highest sampled settings. Funding is not rescued as evidence that the correct z-channel is statistically useful at attainable work.

## Protocol and implementation

The benchmark is exactly the validated 100D Funding problem: `T=0.5`, `sigma=0.2`, `mu=0.06`, `R_l=0.04`, `R_b=0.06`, root `x_i=100`, and reference value `21.299`. The payoff and driver are

`g(x)=(max_i x_i-120)_+ - 2(max_i x_i-150)_+`,

`f(y,z)=-R_l y-((mu-R_l)/sigma) sum_i z_i +(R_b-R_l)(sum_i z_i/sigma-y)_+`.

The solver uses float64 geometric-Brownian transitions, Beta(1/2,1) random-time importance sampling, and the corrected terminal EBL weight `xi/sqrt(T-t)`. Samplewise IR projects `Delta_i=z_i/(sigma*x_i)` onto the model-derived envelope `||Delta||_2 <= exp(sigma^2(T-t)/2)` immediately before every generator reuse. The returned root state is never clipped. This report does not strengthen the theoretical status of the Funding envelope beyond the earlier life-or-death study.

## Reproduction gate

The old `(n,M)=(2,10)`, 100-root seed schedule was reproduced before new sweeps. Raw MAE was exactly `1.5268272487666519`; Samplewise IR was `0.5214655268205621`, differing from the archived value by `3.3e-16`. Raw/IR root terminal blocks and random-tree fingerprints are bitwise paired.

## Completed grid and primary results

Stage F1 covered 15 prescribed settings with 30 paired roots. Five higher-sampling probes and six ultra-high probes were then added adaptively. Ten representative configurations were extended to 100 paired roots, including `(3,48)`. Actual work is reported per root; depth and sampling allocation are not treated as interchangeable.

{_table(["n", "M", "roots", "Raw MAE", "IR MAE", "Raw-IR paired gain", "IR win frac.", "f evals", "samples", "sec/root"], primary_rows)}

## Suppression and shrinkage diagnostics

The constant shrinkage coefficient was selected on a separate 30-root validation seed across `(2,32)`, `(3,16)`, and `(3,24)`. Among `c in {{0.25,0.5,0.75}}`, `c=0.25` was best at all three validation settings and was frozen before the diagnostic test roots were examined. The deliberately invalid tighter envelope used factor `0.75`.

{_table(["n", "M", "method", "roots", "MAE", "bias", "win frac. vs Raw"], diagnostics_rows)}

At `(3,48)`, Raw and IR have entered the expected high-sampling convergence regime: their 100-root MAEs are {_fmt(float(funding[(funding.n==3)&(funding.M==48)&(funding.method=='raw')].mae.iloc[0]))} and {_fmt(float(funding[(funding.n==3)&(funding.M==48)&(funding.method=='samplewise')].mae.iloc[0]))}, respectively. Yet `z=0` is {_fmt(float(funding[(funding.n==3)&(funding.M==48)&(funding.method=='z_zero')].mae.iloc[0]))}. Thus the biased approximation has reached a much lower practical error floor before the correct z-channel becomes worthwhile. `f=0` remains substantially worse where tested, confirming that the value-dependent driver channel is active.

## Interpretation and recommendation

Funding is a strong demonstration that certified projection stabilizes Raw recursive nonlinear reuse: IR wins the paired Raw comparison across the full work range and prevents deep/under-sampled degradation. It is not a clean positive active-gradient benchmark because `z=0` remains better even when Raw approaches IR. It should remain as a secondary mechanism/partial-rescue result, not be promoted beside active VB as evidence that certified geometry makes the correct gradient information useful.

No manuscript file was modified.

**FUNDING VERDICT F-B: PARTIALLY RESCUED.**
"""
    (DOC_ROOT / "FUNDING_HIGH_BUDGET_RESCUE.md").write_text(funding_report, encoding="utf-8")

    h1 = hjb[(hjb.stage == "h1") & (hjb.dimension == 100)]
    h1_rows = []
    for n, M in HJB_H1:
        for method in ("raw", "samplewise", "f_zero"):
            match = h1[(h1.n == n) & (h1.M == M) & (h1.method == method)]
            if match.empty:
                continue
            row = match.iloc[0]
            h1_rows.append([n, M, method, int(row.repetitions), _fmt(row.value_relative_l2_mean), _fmt(row.gradient_relative_l2_mean), _fmt(row.generator_mse_mean), _fmt(row.f_evals, 8), _fmt(row.total_stochastic_samples, 8), _fmt(row.wall_clock_seconds, 5)])
    h2_rows = []
    h2 = hjb[hjb.stage == "h2"]
    for dimension in (100, 140, 160):
        for n, M in HJB_H2:
            for method in ("raw", "samplewise", "f_zero"):
                match = h2[(h2.dimension == dimension) & (h2.n == n) & (h2.M == M) & (h2.method == method)]
                if match.empty:
                    continue
                row = match.iloc[0]
                h2_rows.append([dimension, n, M, method, int(row.repetitions), _fmt(row.value_relative_l2_mean), _fmt(row.value_relative_l2_std), _fmt(row.generator_mse_mean), _fmt(row.total_stochastic_samples, 9), _fmt(row.wall_clock_seconds, 5)])

    raw16 = float(h1[(h1.n==2)&(h1.M==16)&(h1.method=="raw")].generator_mse_mean.iloc[0])
    raw64 = float(h1[(h1.n==2)&(h1.M==64)&(h1.method=="raw")].generator_mse_mean.iloc[0])
    ir16 = float(h1[(h1.n==2)&(h1.M==16)&(h1.method=="samplewise")].generator_mse_mean.iloc[0])
    ir64 = float(h1[(h1.n==2)&(h1.M==64)&(h1.method=="samplewise")].generator_mse_mean.iloc[0])
    hjb_report = f"""# Rosenbrock HJB high-budget rescue

## Outcome

**HJB verdict H-B - TRENDING BUT NOT REACHED.** Along the statistically sensible `n=2` high-sampling path, Raw and Samplewise certified IR errors decrease monotonically, generator MSE falls sharply, and Raw approaches IR. Nevertheless, even at `M=96` the nonlinear methods remain roughly two orders of magnitude above the `f=0` value-error floor. No attainable crossover was observed.

## Validated protocol

The PDE is `u_t + Delta u - ||grad u||^2 = 0`, with `z=sqrt(2) grad(u)` and `f(z)=-0.5||z||^2`. The terminal condition is `log((1+x^T A x)/2)`. The Rosenbrock coefficients exactly reuse the old JAX `PRNGKey(0/1)` construction; at d=100 this gives `trace(A)=304.8856089115143` and `lambda_max=6.383491595137865`. The non-oracle projected radius is `sqrt(15)`. Samplewise correction acts only on z immediately before generator reuse; final roots and values are never clipped. References use the stable scaled Hopf-Cole/Gauss-Laguerre implementation.

The implementation is float64, uses corrected `xi/sqrt(T-t)` terminal EBL normalization and the old uniform random-time protocol, shares complete random trees between Raw and IR, elides the algebraically zero level-0 generator term, and contains no SCaSML heuristic clipping.

## Reproduction gate

The complete old d=100 `(n,M)=(2,10)`, 1200-point, 10-repetition headline was reproduced exactly. Mean relative L2 is `2.2799487622119137` for Raw and `0.6207713181675171` for Samplewise IR, with zero difference from the archived means.

## Stage H1 and stopping-rule decision

H1 used d=100, 300 interior plus 60 boundary points, and 3 paired repetitions.

{_table(["n", "M", "method", "reps", "value relL2", "gradient relL2", "generator MSE", "f evals", "samples", "sec"], h1_rows)}

Generator metrics use a fixed, deterministic cap of 512 child states per H1 method/repetition and 256 per H2 method/repetition; they diagnose the same early nonlinear-reuse locations under paired trees rather than pretending to enumerate every recursive child.

The continuation criterion was met only by the n=2 high-sampling path. From M=16 to M=64, Raw generator MSE fell from {_fmt(raw16)} to {_fmt(raw64)} ({_fmt(raw16/raw64, 4)}x) and IR generator MSE fell from {_fmt(ir16)} to {_fmt(ir64)} ({_fmt(ir16/ir64, 4)}x). Value error also fell by more than 20 percent at successive high-budget points. In contrast, n=3 Raw reached errors from hundreds to tens of thousands, n=4 Raw reached `10^11-10^15`, and IR remained finite but worse than the n=2 path. These are deep/under-sampled failures, not evidence against the high-sampling trend.

## Stage H2

H2 therefore retained only n=2 with M in `{{48,64,96}}`, dimensions `{{100,140,160}}`, 500 interior plus 100 boundary points, and 5 paired repetitions. A matched M=96 `f=0` reference was also run.

{_table(["d", "n", "M", "method", "reps", "value relL2", "std", "generator MSE", "samples", "sec"], h2_rows)}

At d=100, M=96, Raw and IR are approximately 0.218 and 0.190 while `f=0` is 0.00264. At d=160 they are approximately 0.331 and 0.241 while `f=0` is 0.00150. The paired IR improvement is highly stable, but it does not imply that the correct nonlinear correction is practically estimable. H3 was not run: scaling to 1200 points and 10 repetitions would confirm an already stable mean without plausibly closing a 60-160x error gap.

## Interpretation and recommendation

HJB should be presented as the weak-nonlinearity boundary/mechanism case. A valid certified radius materially stabilizes Raw MLP, and the high-sampling trend is real, but the nonlinear correction is so small relative to gradient-estimation noise that the biased `f=0` approximation remains vastly better. It should not remain a positive benchmark for practical certified-geometry accuracy.

No manuscript file was modified.

**HJB VERDICT H-B: TRENDING BUT NOT REACHED.**
"""
    (DOC_ROOT / "HJB_HIGH_BUDGET_RESCUE.md").write_text(hjb_report, encoding="utf-8")

    summary_report = f"""# Old benchmark high-budget rescue summary

## Executive result

The two older benchmarks do not join active VB as clean positive benchmarks, although both show that Samplewise certified IR can stabilize Raw full-history MLP.

| PDE | Verdict | Raw vs IR | Suppression floor crossed? | Paper role |
| --- | --- | --- | --- | --- |
| 100D Funding | F-B, partially rescued | IR consistently improves Raw and Raw approaches IR at high sampling | No; `z=0` remains much better | Secondary stabilization/mechanism result |
| 100--160D Rosenbrock HJB | H-B, trending but not reached | IR strongly stabilizes Raw; both improve along n=2 high sampling | No; `f=0` remains 60--160x better at M=96 | Weak-nonlinearity limitation/boundary case |

## Scientific recommendation

Active VB should remain the principal positive nonlinear benchmark. Funding may be retained as realistic evidence that certified projection regularizes Raw recursion, but not as proof that the correct gradient channel is useful at attainable budget. Rosenbrock HJB should be explicitly reframed as a mechanism/limitation example: correctness-preserving geometry can be statistically unhelpful when the true nonlinear correction is much smaller than gradient Monte Carlo noise.

The primary comparison was Raw versus Samplewise certified IR throughout. Suppression, fixed shrinkage, and an invalid tighter envelope were restricted to representative diagnostics. Batch-IR was omitted from the new study.

## Compute and reproducibility

The recorded parallel orchestration wall time is {_fmt(compute['orchestration_wall_seconds']/60, 5)} minutes; summed per-task wall time is {_fmt(compute['task_wall_seconds']/3600, 5)} CPU-task hours. The experiment contains {compute['funding_repetitions']} Funding root repetitions and {compute['hjb_repetitions']} HJB method-repetitions. All scientific arrays are float64. Legacy headline reproduction, paired terminal blocks, random-tree fingerprints, finite-state checks, work accounting, and manifest completeness are covered by the integrity audit.

Detailed reports:

- `docs/FUNDING_HIGH_BUDGET_RESCUE.md`
- `docs/HJB_HIGH_BUDGET_RESCUE.md`

Results and figures are under `results/rescue_high_budget/`. No manuscript source was modified.
"""
    (DOC_ROOT / "OLD_BENCHMARK_RESCUE_SUMMARY.md").write_text(summary_report, encoding="utf-8")


def main() -> None:
    RESULT_ROOT.mkdir(parents=True, exist_ok=True)
    funding_rep = funding_repetitions()
    funding_summary = summarize_funding(funding_rep)
    funding_validation_rep = funding_repetitions(base_seed=20261105)
    funding_validation_summary = summarize_funding(funding_validation_rep)
    hjb_rep, corrections = hjb_repetitions()
    hjb_summary = summarize_hjb(hjb_rep, corrections)
    funding_rep.to_csv(RESULT_ROOT / "funding_repetitions.csv", index=False)
    funding_summary.to_csv(RESULT_ROOT / "funding_summary.csv", index=False)
    funding_validation_summary.to_csv(
        RESULT_ROOT / "funding_shrink_validation.csv", index=False
    )
    hjb_rep.to_csv(RESULT_ROOT / "hjb_repetitions.csv", index=False)
    hjb_summary.to_csv(RESULT_ROOT / "hjb_summary.csv", index=False)
    work = pd.concat(
        [
            funding_summary.assign(
                pde="funding", stage="all", dimension=100,
                primary_error=funding_summary.mae,
            )[["pde", "stage", "dimension", "n", "M", "method", "repetitions", "primary_error", "f_evals", "total_stochastic_samples", "wall_clock_seconds"]],
            hjb_summary.assign(
                pde="hjb", primary_error=hjb_summary.value_relative_l2_mean,
            )[["pde", "stage", "dimension", "n", "M", "method", "repetitions", "primary_error", "f_evals", "total_stochastic_samples", "wall_clock_seconds"]],
        ],
        ignore_index=True,
    )
    work.to_csv(RESULT_ROOT / "work_summary.csv", index=False)
    figures = make_figures(funding_summary, hjb_summary)

    orchestration_wall, manifest_times = _manifest_wall_seconds()
    # Sum task wall directly from one representative row per source and restore
    # its block size for Funding; HJB rows already represent one task.
    funding_task_wall = 0.0
    for _, frame in funding_rep.groupby("source_path"):
        funding_task_wall += float(frame.wall_clock_seconds.iloc[0] * len(frame))
    hjb_task_wall = float(hjb_rep.wall_clock_seconds.sum())
    compute = {
        "orchestration_wall_seconds": orchestration_wall,
        "task_wall_seconds": funding_task_wall + hjb_task_wall,
        "funding_task_wall_seconds": funding_task_wall,
        "hjb_task_wall_seconds": hjb_task_wall,
        "funding_repetitions": len(funding_rep),
        "hjb_repetitions": len(hjb_rep),
        "manifest_times": manifest_times,
    }
    provenance = {
        "schema_version": 1,
        "branch": _git("branch", "--show-current"),
        "analysis_commit": _git("rev-parse", "HEAD"),
        "source_base_commits": {
            "life_or_death": "fd2eece949b3a72ab24986cb8cba00e4adcc863c",
            "active_vb_high_budget": "b333bf8796b789e922c0e6c1ed4ecdeaa6593957",
        },
        "legacy_sources": {
            str(path.relative_to(REPOSITORY_ROOT)).replace("\\", "/"): _sha256(path)
            for path in (
                PROJECT_ROOT / "experiments" / "funding_life_or_death.py",
                PROJECT_ROOT / "experiments" / "batch_contraction" / "hjb_batchwise_headline.py",
            )
        },
        "new_sources": {
            path.name: _sha256(path)
            for path in sorted(HERE.glob("*.py"))
        },
        "completed_grids": {
            "funding_f1_30_roots": FUNDING_F1,
            "funding_high_probe_30_roots": FUNDING_HIGH_PROBE,
            "funding_final_100_roots": FUNDING_F2,
            "funding_ultra_probe_30_roots": FUNDING_ULTRA_PROBE,
            "funding_diagnostic_configs": FUNDING_DIAGNOSTIC,
            "hjb_h1_d100_360_points_3_reps": HJB_H1,
            "hjb_h2_d100_140_160_600_points_5_reps": HJB_H2,
        },
        "verdicts": {"funding": "F-B", "hjb": "H-B"},
        "compute": compute,
        "figures": [str(path.relative_to(RESULT_ROOT)).replace("\\", "/") for path in figures],
    }
    (RESULT_ROOT / "provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True), encoding="utf-8"
    )
    full_summary = {
        "funding": json.loads(funding_summary.to_json(orient="records")),
        "hjb": json.loads(hjb_summary.to_json(orient="records")),
        "compute": compute,
        "verdicts": provenance["verdicts"],
    }
    (RESULT_ROOT / "full_summary.json").write_text(
        json.dumps(full_summary, indent=2, allow_nan=False), encoding="utf-8"
    )
    write_reports(funding_summary, hjb_summary, compute)
    print(f"Funding repetition rows: {len(funding_rep)}")
    print(f"HJB repetition rows: {len(hjb_rep)}")
    print(f"Figures: {len(figures)}")
    print(f"Recorded orchestration minutes: {orchestration_wall/60:.3f}")


if __name__ == "__main__":
    main()
