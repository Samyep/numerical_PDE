"""Aggregate and report the paired actor-bottleneck diagnostics."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr

from .common import PROJECT_DIR, RESULTS_DIR, atomic_json, load_config
from .run_actor_diagnostics import DIAGNOSTIC_RESULTS_DIR, RAW_DIR


REPORT_PATH = PROJECT_DIR / "docs" / "MAIN_B1_ACTOR_DIAGNOSTICS.md"
FIGURE_DIR = DIAGNOSTIC_RESULTS_DIR / "figures"
VARIANTS = (
    "learned_actor",
    "oracle_actor",
    "pathwise_baseline_no_actor",
)
DISPLAY = {
    "learned_actor": "Learned Actor",
    "oracle_actor": "Oracle Actor",
    "pathwise_baseline_no_actor": "Pathwise + Baseline, no Actor",
}
COLORS = {
    "learned_actor": "#d62728",
    "oracle_actor": "#111111",
    "pathwise_baseline_no_actor": "#1f77b4",
}
MARKERS = {
    "learned_actor": "*",
    "oracle_actor": "D",
    "pathwise_baseline_no_actor": "o",
}


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = list(rows[0]) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _load_payload(variant: str, d: int, k: int) -> dict[str, Any]:
    if variant == "learned_actor":
        path = RESULTS_DIR / "raw" / "stage_b" / f"actor_mlp_d{d}_k{k}.json"
    else:
        path = RAW_DIR / f"{variant}_d{d}_k{k}.json"
    if not path.is_file():
        raise FileNotFoundError(path)
    payload = json.loads(path.read_text())
    if payload.get("status") != "complete" or len(payload.get("rows", [])) != 10:
        raise RuntimeError(f"incomplete diagnostic input: {path}")
    return payload


def _pooled(rows: list[dict[str, Any]], key: str) -> float:
    values = np.asarray([float(row[key]) for row in rows], dtype=np.float64)
    return float(np.sqrt(np.mean(values**2)))


def _mean(rows: list[dict[str, Any]], key: str) -> float:
    return float(np.mean([float(row[key]) for row in rows]))


def _ci95(values: list[float]) -> float:
    if len(values) <= 1:
        return 0.0
    return 1.96 * float(np.std(values, ddof=1)) / math.sqrt(len(values))


def _aggregate_payload(
    variant: str, d: int, k: int, payload: dict[str, Any]
) -> dict[str, Any]:
    rows = sorted(payload["rows"], key=lambda row: int(row["repetition"]))
    value_rmse = [float(row["value_rmse"]) for row in rows]
    gradient_rmse = [float(row["gradient_rmse"]) for row in rows]
    actor_metadata = payload.get("shared", {}).get("actor", {})
    return {
        "d": d,
        "k": k,
        "variant": variant,
        "repetitions": len(rows),
        "value_rmse": float(np.sqrt(np.mean(np.square(value_rmse)))),
        "mean_per_seed_value_rmse": float(np.mean(value_rmse)),
        "per_seed_value_rmse_ci95": _ci95(value_rmse),
        "normalized_value_rmse": _mean(rows, "normalized_value_rmse"),
        "gradient_rmse": float(np.sqrt(np.mean(np.square(gradient_rmse)))),
        "relative_gradient_error": _mean(rows, "relative_gradient_error"),
        "value_bias": _mean(rows, "value_bias"),
        "nonfinite_count": sum(
            int(row["value_nonfinite"]) + int(row["gradient_nonfinite"])
            for row in rows
        ),
        "one_time_training_seconds": float(
            actor_metadata.get("total_actor_seconds", 0.0)
        ),
        "total_inference_seconds": sum(
            float(row["inference_seconds"]) for row in rows
        ),
        "peak_gpu_bytes": max(int(row["cost_peak_gpu_bytes"]) for row in rows),
        "oracle_reference_used_inside_solver": variant == "oracle_actor",
        "test_point_sha256": rows[0]["test_point_sha256"],
    }


def _accuracy_figure(results: list[dict[str, Any]]) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.2), sharey=True)
    for axis, d in zip(axes, (100, 200, 400)):
        for variant in VARIANTS:
            selected = sorted(
                (
                    row
                    for row in results
                    if int(row["d"]) == d and row["variant"] == variant
                ),
                key=lambda row: int(row["k"]),
            )
            axis.plot(
                [int(row["k"]) for row in selected],
                [float(row["value_rmse"]) for row in selected],
                label=DISPLAY[variant],
                color=COLORS[variant],
                marker=MARKERS[variant],
                linewidth=1.8,
                markersize=9 if variant == "learned_actor" else 6,
            )
        axis.set_title(f"d = {d}")
        axis.set_xlabel("effective dimension k")
        axis.set_xticks([1, 5, 10, 20, 40])
        axis.set_yscale("log")
        axis.grid(True, which="both", alpha=0.25)
    axes[0].set_ylabel("value RMSE")
    handles, labels = axes[-1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURE_DIR / f"actor_diagnostic_value_rmse.{suffix}", dpi=220)
    plt.close(fig)


def _quality_figure(comparisons: list[dict[str, Any]]) -> None:
    fig, axis = plt.subplots(figsize=(7.2, 5.0))
    for d, color in zip((100, 200, 400), ("#1f77b4", "#ff7f0e", "#2ca02c")):
        selected = [row for row in comparisons if int(row["d"]) == d]
        axis.scatter(
            [float(row["actor_relative_gradient_error"]) for row in selected],
            [float(row["learned_over_oracle_rmse"]) for row in selected],
            label=f"d = {d}",
            color=color,
            s=55,
        )
        for row in selected:
            axis.annotate(
                f"k={row['k']}",
                (
                    float(row["actor_relative_gradient_error"]),
                    float(row["learned_over_oracle_rmse"]),
                ),
                xytext=(4, 3),
                textcoords="offset points",
                fontsize=7,
            )
    axis.axhline(1.0, color="#666666", linestyle="--", linewidth=1.0)
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("learned actor relative gradient error")
    axis.set_ylabel("learned / oracle value RMSE")
    axis.grid(True, which="both", alpha=0.25)
    axis.legend(frameon=False)
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURE_DIR / f"actor_quality_vs_oracle_gap.{suffix}", dpi=220)
    plt.close(fig)


def _artifact_manifest() -> dict[str, Any]:
    paths = [REPORT_PATH]
    paths.extend(
        path
        for path in DIAGNOSTIC_RESULTS_DIR.rglob("*")
        if path.is_file()
        and path.name != "artifact_manifest.json"
        and path.suffix not in {".npz", ".pt"}
    )
    artifacts = []
    for path in sorted(set(paths)):
        data = path.read_bytes()
        artifacts.append(
            {
                "path": str(path.relative_to(PROJECT_DIR)),
                "bytes": len(data),
                "sha256": hashlib.sha256(data).hexdigest(),
            }
        )
    manifest = {
        "schema_version": 1,
        "artifact_count": len(artifacts),
        "artifacts": artifacts,
        "excluded_local_binary_suffixes": [".npz", ".pt"],
    }
    atomic_json(DIAGNOSTIC_RESULTS_DIR / "artifact_manifest.json", manifest)
    return manifest


def aggregate_actor_diagnostics() -> dict[str, Any]:
    config = load_config()
    dimensions = [int(value) for value in config["pde"]["dimensions"]]
    effective_dimensions = [
        int(value) for value in config["pde"]["effective_dimensions"]
    ]
    actor_training_rows = list(
        csv.DictReader((RESULTS_DIR / "actor_training.csv").open())
    )
    actor_training = {
        (int(row["d"]), int(row["k"])): row for row in actor_training_rows
    }
    results: list[dict[str, Any]] = []
    payloads: dict[tuple[int, int, str], dict[str, Any]] = {}
    for d in dimensions:
        for k in effective_dimensions:
            hashes: set[str] = set()
            for variant in VARIANTS:
                payload = _load_payload(variant, d, k)
                payloads[(d, k, variant)] = payload
                result = _aggregate_payload(variant, d, k, payload)
                results.append(result)
                hashes.add(str(result["test_point_sha256"]))
            if len(hashes) != 1:
                raise RuntimeError(f"test point mismatch for d={d}, k={k}")
    if len(results) != 45:
        raise RuntimeError(f"expected 45 method-cells, found {len(results)}")
    if any(int(row["nonfinite_count"]) for row in results):
        raise RuntimeError("non-finite diagnostic predictions found")
    _write_csv(DIAGNOSTIC_RESULTS_DIR / "actor_diagnostic_results.csv", results)

    index = {
        (int(row["d"]), int(row["k"]), str(row["variant"])): row
        for row in results
    }
    comparisons: list[dict[str, Any]] = []
    for d in dimensions:
        for k in effective_dimensions:
            learned = index[(d, k, "learned_actor")]
            oracle = index[(d, k, "oracle_actor")]
            no_actor = index[(d, k, "pathwise_baseline_no_actor")]
            learned_rows = sorted(
                payloads[(d, k, "learned_actor")]["rows"],
                key=lambda row: int(row["repetition"]),
            )
            oracle_rows = sorted(
                payloads[(d, k, "oracle_actor")]["rows"],
                key=lambda row: int(row["repetition"]),
            )
            no_actor_rows = sorted(
                payloads[(d, k, "pathwise_baseline_no_actor")]["rows"],
                key=lambda row: int(row["repetition"]),
            )
            if not all(
                int(a["repetition"]) == int(b["repetition"]) == int(c["repetition"])
                for a, b, c in zip(learned_rows, oracle_rows, no_actor_rows)
            ):
                raise RuntimeError(f"repetition mismatch for d={d}, k={k}")
            oracle_minus_learned = [
                float(o["value_rmse"]) - float(l["value_rmse"])
                for l, o in zip(learned_rows, oracle_rows)
            ]
            learned_minus_no_actor = [
                float(l["value_rmse"]) - float(n["value_rmse"])
                for l, n in zip(learned_rows, no_actor_rows)
            ]
            actor = actor_training[(d, k)]
            actor_rmse = float(actor["actor_gradient_rmse"])
            comparisons.append(
                {
                    "d": d,
                    "k": k,
                    "learned_actor_value_rmse": learned["value_rmse"],
                    "oracle_actor_value_rmse": oracle["value_rmse"],
                    "pathwise_baseline_no_actor_value_rmse": no_actor["value_rmse"],
                    "oracle_over_learned_rmse": float(oracle["value_rmse"])
                    / float(learned["value_rmse"]),
                    "learned_over_oracle_rmse": float(learned["value_rmse"])
                    / float(oracle["value_rmse"]),
                    "learned_over_no_actor_rmse": float(learned["value_rmse"])
                    / float(no_actor["value_rmse"]),
                    "oracle_over_no_actor_rmse": float(oracle["value_rmse"])
                    / float(no_actor["value_rmse"]),
                    "paired_mean_oracle_minus_learned_seed_rmse": float(
                        np.mean(oracle_minus_learned)
                    ),
                    "paired_ci95_oracle_minus_learned_seed_rmse": _ci95(
                        oracle_minus_learned
                    ),
                    "paired_mean_learned_minus_no_actor_seed_rmse": float(
                        np.mean(learned_minus_no_actor)
                    ),
                    "paired_ci95_learned_minus_no_actor_seed_rmse": _ci95(
                        learned_minus_no_actor
                    ),
                    "actor_gradient_rmse": actor_rmse,
                    "actor_relative_gradient_error": float(
                        actor["actor_relative_gradient_error"]
                    ),
                    "mean_learned_driver_gap_at_exact_z": 0.5 * actor_rmse**2,
                }
            )
    _write_csv(
        DIAGNOSTIC_RESULTS_DIR / "actor_diagnostic_comparison.csv", comparisons
    )

    def geometric_mean(key: str) -> float:
        return float(
            np.exp(np.mean(np.log([float(row[key]) for row in comparisons])))
        )

    oracle_better = sum(
        float(row["oracle_actor_value_rmse"])
        < float(row["learned_actor_value_rmse"])
        for row in comparisons
    )
    learned_better_no_actor = sum(
        float(row["learned_actor_value_rmse"])
        < float(row["pathwise_baseline_no_actor_value_rmse"])
        for row in comparisons
    )
    oracle_better_no_actor = sum(
        float(row["oracle_actor_value_rmse"])
        < float(row["pathwise_baseline_no_actor_value_rmse"])
        for row in comparisons
    )
    correlation = spearmanr(
        [float(row["actor_relative_gradient_error"]) for row in comparisons],
        [float(row["learned_over_oracle_rmse"]) for row in comparisons],
    )
    summary: dict[str, Any] = {
        "schema_version": 1,
        "cells": 15,
        "diagnostic_method_cells": 30,
        "diagnostic_repetitions": 300,
        "oracle_better_than_learned_cells": int(oracle_better),
        "learned_better_than_no_actor_cells": int(learned_better_no_actor),
        "oracle_better_than_no_actor_cells": int(oracle_better_no_actor),
        "geometric_mean_oracle_over_learned_rmse": geometric_mean(
            "oracle_over_learned_rmse"
        ),
        "geometric_mean_learned_over_no_actor_rmse": geometric_mean(
            "learned_over_no_actor_rmse"
        ),
        "geometric_mean_oracle_over_no_actor_rmse": geometric_mean(
            "oracle_over_no_actor_rmse"
        ),
        "spearman_actor_relative_error_vs_learned_over_oracle_ratio": float(
            correlation.statistic
        ),
        "spearman_pvalue": float(correlation.pvalue),
        "learned_actor_training_seconds": sum(
            float(row["training_seconds"]) + float(row["label_generation_seconds"])
            for row in actor_training_rows
        ),
        "learned_actor_inference_seconds": sum(
            float(row["total_inference_seconds"])
            for row in results
            if row["variant"] == "learned_actor"
        ),
        "oracle_actor_inference_seconds": sum(
            float(row["total_inference_seconds"])
            for row in results
            if row["variant"] == "oracle_actor"
        ),
        "no_actor_inference_seconds": sum(
            float(row["total_inference_seconds"])
            for row in results
            if row["variant"] == "pathwise_baseline_no_actor"
        ),
        "all_nonfinite_counts": 0,
        "paired_seed_namespace": "B1|actor_mlp",
        "oracle_is_diagnostic_only": True,
    }
    atomic_json(DIAGNOSTIC_RESULTS_DIR / "summary.json", summary)
    _accuracy_figure(results)
    _quality_figure(comparisons)

    rows_markdown: list[str] = []
    for row in comparisons:
        values = {
            "learned": float(row["learned_actor_value_rmse"]),
            "oracle": float(row["oracle_actor_value_rmse"]),
            "no_actor": float(row["pathwise_baseline_no_actor_value_rmse"]),
        }
        best = min(values, key=values.get)

        def formatted(name: str) -> str:
            text = f"{values[name]:.4g}"
            return f"**{text}**" if name == best else text

        rows_markdown.append(
            f"| {row['d']} | {row['k']} | {formatted('learned')} | "
            f"{formatted('oracle')} | {formatted('no_actor')} | "
            f"{float(row['oracle_over_learned_rmse']):.3f} | "
            f"{float(row['learned_over_no_actor_rmse']):.3f} | "
            f"{float(row['actor_relative_gradient_error']):.3f} |"
        )

    oracle_ratio = summary["geometric_mean_oracle_over_learned_rmse"]
    learned_no_actor_ratio = summary["geometric_mean_learned_over_no_actor_rmse"]
    oracle_no_actor_ratio = summary["geometric_mean_oracle_over_no_actor_rmse"]
    if oracle_ratio < 0.8 and oracle_better >= 10:
        bottleneck_text = (
            "The oracle actor produces a large and broad accuracy improvement, so learned-actor "
            "error is a major current bottleneck."
        )
    elif oracle_ratio < 1.0 and oracle_better >= 8:
        bottleneck_text = (
            "The oracle actor improves most cells, but the aggregate gain is moderate; actor "
            "error is one bottleneck rather than the sole bottleneck."
        )
    else:
        bottleneck_text = (
            "The oracle actor does not give a broad material improvement, so actor prediction "
            "error is not supported as the dominant bottleneck at this MLP budget."
        )
    if learned_no_actor_ratio < 1.0:
        actor_benefit_text = (
            "The learned actor improves geometric-mean RMSE over Pathwise + Baseline without an "
            "actor, although the cell count shows how consistent that benefit is."
        )
    else:
        actor_benefit_text = (
            "The learned actor worsens geometric-mean RMSE relative to Pathwise + Baseline without "
            "an actor at the current training budget."
        )

    report = [
        "# B1 actor bottleneck diagnostics",
        "",
        "## Numerical diagnosis",
        "",
        f"{bottleneck_text} Oracle Actor beats the learned actor in "
        f"{oracle_better}/15 cells, with geometric-mean Oracle/Learned RMSE ratio "
        f"{oracle_ratio:.4f}.",
        "",
        f"{actor_benefit_text} Learned Actor beats the no-actor variant in "
        f"{learned_better_no_actor}/15 cells, with geometric-mean Learned/No-Actor "
        f"RMSE ratio {learned_no_actor_ratio:.4f}. Oracle Actor beats the no-actor "
        f"variant in {oracle_better_no_actor}/15 cells, with geometric-mean "
        f"Oracle/No-Actor ratio {oracle_no_actor_ratio:.4f}.",
        "",
        "## Scope and pairing",
        "",
        "These are supplemental diagnostics and do not replace any primary six-method result. "
        "Every cell uses the frozen Stage B PDE, 960 held-out points, `n=3`, `M=16`, and 10 "
        "repetitions. Learned Actor, Oracle Actor, and Pathwise + Baseline without Actor use "
        "the same random seed for each `(d,k,repetition)`, enabling paired comparisons.",
        "",
        "Oracle Actor evaluates `q(t,x)=z*(t,x)` at every recursive child and in the matching "
        "level-1 baseline. It uses forbidden oracle information and is only a diagnostic upper "
        "bound. The no-actor variant uses the original quadratic driver with the pathwise "
        "terminal gradient and original-driver level-matched baseline.",
        "",
        "## Value RMSE",
        "",
        "| d | k | Learned Actor | Oracle Actor | Pathwise + Baseline, no Actor | Oracle/Learned | Learned/No Actor | Actor relative error |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
        *rows_markdown,
        "",
        "## Actor quality relation",
        "",
        f"Across the 15 cells, the Spearman correlation between held-out learned-actor relative "
        f"gradient error and the Learned/Oracle value-RMSE ratio is "
        f"{summary['spearman_actor_relative_error_vs_learned_over_oracle_ratio']:.3f} "
        f"(two-sided p={summary['spearman_pvalue']:.3g}). For each cell, "
        "`actor_diagnostic_comparison.csv` also reports the implied mean tangent-driver gap "
        "`0.5 * actor_gradient_rmse^2` at the exact gradient.",
        "",
        "## Measured cost",
        "",
        f"The existing 15 learned actors used {summary['learned_actor_training_seconds']:.3f} "
        f"seconds for label generation plus training and "
        f"{summary['learned_actor_inference_seconds']:.3f} seconds for the common MLP workload. "
        f"Oracle diagnostic inference used {summary['oracle_actor_inference_seconds']:.3f} "
        f"seconds; Pathwise + Baseline without Actor used "
        f"{summary['no_actor_inference_seconds']:.3f} seconds. Oracle cost is descriptive only "
        "because exact-gradient access is unavailable in a practical solver.",
        "",
        "## Artifacts",
        "",
        "- `actor_diagnostic_results.csv`: all 45 learned/oracle/no-actor method-cells.",
        "- `actor_diagnostic_comparison.csv`: paired ratios, differences, and actor-quality metrics.",
        "- `summary.json`: aggregate counts, geometric means, correlation, and runtime.",
        "- `figures/actor_diagnostic_value_rmse.*` and `figures/actor_quality_vs_oracle_gap.*`.",
        "",
    ]
    REPORT_PATH.write_text("\n".join(report))
    manifest = _artifact_manifest()
    summary["report"] = str(REPORT_PATH)
    summary["artifact_count"] = manifest["artifact_count"]
    return summary
