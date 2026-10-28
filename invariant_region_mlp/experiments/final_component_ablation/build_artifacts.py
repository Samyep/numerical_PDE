"""Aggregate completed cells and render the publication ablation artifacts."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import shutil
import subprocess
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .protocol import FIGURE_DIR, PDES, PROJECT_DIR, RESULTS_DIR, TABLE_DIR, VARIANTS


METHOD_NAMES = {
    "a0_full": "Tangent-Actor (Full)",
    "a1_no_mcv": "w/o MCV",
    "a2_no_terminal_anchor": "w/o Terminal Anchoring",
    "a3_no_updates": "w/o Tangent Updates",
}
PDE_NAMES = {
    "hjb_log": "Radial Log-Terminal HJB",
    "joint_concave": "Concave Cubic--Quadratic",
    "joint_convex": "Convex Cubic--Quadratic",
    "nonconvex_coupled": "Oscillatory Value--Gradient",
}
COLORS = {
    "a0_full": "#1f4e79",
    "a1_no_mcv": "#c55a11",
    "a2_no_terminal_anchor": "#70ad47",
    "a3_no_updates": "#7f6000",
}


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise RuntimeError(f"cannot write empty CSV: {path}")
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _evaluation_rows() -> list[dict[str, Any]]:
    rows = []
    for path in sorted((RESULTS_DIR / "raw").glob("evaluation_*_rep*.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("status") != "complete":
            continue
        rows.append({**payload, "source_json": str(path)})
    return rows


def _expected_count() -> int:
    return 3 * len(PDES) * 3 * 2 + len(PDES) * 2


def _mean(values: Iterable[float]) -> float:
    return float(np.mean(list(values)))


def _sample_std(values: Iterable[float]) -> float:
    array = np.asarray(list(values), dtype=float)
    return float(np.std(array, ddof=1)) if len(array) > 1 else 0.0


def _summary(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    result = []
    for variant in VARIANTS:
        for pde in PDES:
            local = [row for row in rows if row["variant"] == variant and row["pde"] == pde]
            primary = [row for row in local if row["primary_replicate"]]
            repeat = [row for row in local if not row["primary_replicate"]]
            if variant == "a3_no_updates":
                if len(primary) != 1 or len(repeat) != 1:
                    raise RuntimeError(f"incomplete deterministic A3 cell {pde}")
                both = primary + repeat
                center_rows = both
                variability = "two independent inference replicates; no training-seed variation"
            else:
                if len(primary) != 3 or len(repeat) != 3:
                    raise RuntimeError(f"incomplete three-seed cell {variant}/{pde}")
                center_rows = primary
                variability = "sample standard deviation across three training-seed primary RMSEs"
            result.append({
                "variant": variant,
                "method": METHOD_NAMES[variant],
                "pde": pde,
                "pde_name": PDE_NAMES[pde],
                "rmse_mean": _mean(float(row["value_rmse"]) for row in center_rows),
                "rmse_std": _sample_std(float(row["value_rmse"]) for row in center_rows),
                "rmae_mean": _mean(float(row["rmae"]) for row in center_rows),
                "rmae_std": _sample_std(float(row["rmae"]) for row in center_rows),
                "primary_seed_count": len(primary),
                "inference_replicates_per_actor": 2,
                "variability_definition": variability,
                "independent_replicate_rmse_mean": _mean(float(row["value_rmse"]) for row in repeat),
                "primary_minus_repeat_rmse_mean": _mean(
                    float(p["value_rmse"]) - float(r["value_rmse"])
                    for p, r in zip(sorted(primary, key=lambda x: x["seed"]), sorted(repeat, key=lambda x: x["seed"]), strict=True)
                ),
                "training_seconds_mean": _mean(float(row["training_seconds"]) for row in primary),
                "inference_seconds_mean": _mean(float(row["inference_seconds"]) for row in local),
                "mc_rmse_contribution_mean": _mean(float(row["mc_rmse_contribution"]) for row in local),
                "weight_ess_fraction_min": min(float(row["weight_ess_fraction"]) for row in local),
                "numerical_failures": sum(int(row["nonfinite_count"]) > 0 for row in local),
                "peak_gpu_bytes_max": max(int(row["peak_gpu_bytes"]) for row in local),
            })
    return result


def _paired(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    primary = {
        (row["pde"], row["variant"], int(row["seed"])): row
        for row in rows
        if row["primary_replicate"] and row["variant"] != "a3_no_updates"
    }
    output = []
    for pde in PDES:
        for variant in ("a1_no_mcv", "a2_no_terminal_anchor"):
            for seed in (0, 1, 2):
                full = primary[(pde, "a0_full", seed)]
                ablated = primary[(pde, variant, seed)]
                output.append({
                    "pde": pde,
                    "variant": variant,
                    "seed": seed,
                    "full_rmse": full["value_rmse"],
                    "ablated_rmse": ablated["value_rmse"],
                    "ablated_minus_full_rmse": float(ablated["value_rmse"]) - float(full["value_rmse"]),
                    "rmse_ratio_ablated_over_full": float(ablated["value_rmse"]) / float(full["value_rmse"]),
                    "common_test_bank": full["test_bank_array_sha256"] == ablated["test_bank_array_sha256"],
                })
    return output


def _mcv_inference_only(
    rows: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Summarize plain versus MCV inference for the same frozen Full actors."""
    seed_rows: list[dict[str, Any]] = []
    full_rows = [row for row in rows if row["variant"] == "a0_full"]
    for row in full_rows:
        mcv = float(row["mcv_rmse_same_paths"])
        plain = float(row["plain_rmse_same_paths"])
        seed_rows.append({
            "pde": row["pde"],
            "seed": int(row["seed"]),
            "inference_replicate": int(row["inference_replicate"]),
            "primary_replicate": bool(row["primary_replicate"]),
            "full_mcv_rmse": mcv,
            "full_without_mcv_inference_rmse": plain,
            "plain_minus_mcv_rmse": plain - mcv,
            "plain_over_mcv_rmse_ratio": plain / mcv,
            "same_actor": True,
            "same_brownian_paths": True,
            "training_and_selection_used_mcv": True,
            "variance_reduction": float(row["variance_reduction"]),
        })
    summary: list[dict[str, Any]] = []
    for pde in PDES:
        primary = [
            row
            for row in seed_rows
            if row["pde"] == pde and row["primary_replicate"]
        ]
        repeat = [
            row
            for row in seed_rows
            if row["pde"] == pde and not row["primary_replicate"]
        ]
        if len(primary) != 3 or len(repeat) != 3:
            raise RuntimeError(f"incomplete Full-actor inference-only MCV diagnostic: {pde}")
        summary.append({
            "pde": pde,
            "pde_name": PDE_NAMES[pde],
            "training_seed_count": 3,
            "inference_replicates_per_actor": 2,
            "full_mcv_rmse_mean": _mean(row["full_mcv_rmse"] for row in primary),
            "full_mcv_rmse_std": _sample_std(row["full_mcv_rmse"] for row in primary),
            "full_without_mcv_inference_rmse_mean": _mean(
                row["full_without_mcv_inference_rmse"] for row in primary
            ),
            "full_without_mcv_inference_rmse_std": _sample_std(
                row["full_without_mcv_inference_rmse"] for row in primary
            ),
            "paired_plain_minus_mcv_rmse_mean": _mean(
                row["plain_minus_mcv_rmse"] for row in primary
            ),
            "plain_over_mcv_rmse_ratio_mean": _mean(
                row["plain_over_mcv_rmse_ratio"] for row in primary
            ),
            "variance_reduction_mean": _mean(
                row["variance_reduction"] for row in primary
            ),
            "independent_plain_rmse_mean": _mean(
                row["full_without_mcv_inference_rmse"] for row in repeat
            ),
            "independent_mcv_rmse_mean": _mean(
                row["full_mcv_rmse"] for row in repeat
            ),
        })
    return seed_rows, summary


def _iteration_rows() -> list[dict[str, Any]]:
    source = RESULTS_DIR.parent / "tangent_actor_icml_final_suite" / "ABLATION_RESULTS.csv"
    rows = _read_csv(source)
    output = []
    seen: set[tuple[str, int, int]] = set()
    for row in rows:
        if row.get("component") != "learned_tangent_updates":
            continue
        pde = row.get("pde")
        if pde not in ("joint_concave", "joint_convex", "nonconvex_coupled"):
            continue
        if row.get("dimension") != "100":
            continue
        iteration_text = row.get("iteration", "")
        if not iteration_text.startswith("a"):
            continue
        seed = int(float(row["seed"]))
        iteration = int(iteration_text[1:])
        key = (pde, seed, iteration)
        if key in seen:
            continue
        seen.add(key)
        output.append({
            "pde": pde,
            "seed": seed,
            "iteration": iteration,
            "value_rmse": float(row["value_rmse"]),
            "steps": int(float(row["steps"])),
            "paths": int(float(row["paths"])),
            "source_artifact": row["source_artifact"],
            "evidence_role": "audited_existing_within_curve_diagnostic",
            "bank_limitation": "512-state historical bank"
            if pde == "nonconvex_coupled"
            else "960-state frozen bank",
        })
    return sorted(output, key=lambda row: (row["pde"], row["seed"], row["iteration"]))


def _format_value(mean: float, std: float) -> str:
    return rf"${mean:.2e} \pm {std:.2e}$"


def _table(summary: list[dict[str, Any]]) -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    lookup = {(row["variant"], row["pde"]): row for row in summary}
    best = {
        pde: min(float(lookup[(variant, pde)]["rmse_mean"]) for variant in VARIANTS)
        for pde in PDES
    }
    lines = [
        r"\begin{table*}[t]",
        r"\centering",
        r"\caption{Final component ablation at $d=100$. Entries are whole-function value RMSE (mean $\pm$ sample standard deviation across three training seeds). For the deterministic fixed-$A_0$ row, the two values are independent inference replicates rather than training seeds.}",
        r"\label{tab:final-component-ablation}",
        r"\resizebox{\textwidth}{!}{%",
        r"\begin{tabular}{lcccc}",
        r"\toprule",
        r"Method & Radial Log-Terminal HJB & Concave Cubic--Quadratic & Convex Cubic--Quadratic & Oscillatory Value--Gradient \\",
        r"\midrule",
    ]
    for variant in VARIANTS:
        values = []
        for pde in PDES:
            row = lookup[(variant, pde)]
            text = _format_value(float(row["rmse_mean"]), float(row["rmse_std"]))
            if math.isclose(float(row["rmse_mean"]), best[pde], rel_tol=1e-12, abs_tol=0.0):
                text = r"\mathbf{" + text[1:-1] + "}"
                text = "$" + text + "$"
            values.append(text)
        lines.append(METHOD_NAMES[variant] + " & " + " & ".join(values) + r" \\")
    lines.extend([
        r"\bottomrule",
        r"\end{tabular}%",
        r"}",
        r"\end{table*}",
    ])
    table_path = TABLE_DIR / "final_component_ablation.tex"
    table_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    preview = TABLE_DIR / "final_component_ablation_preview_source.tex"
    preview.write_text("\n".join([
        r"\documentclass{article}",
        r"\usepackage[margin=0.4in,landscape]{geometry}",
        r"\usepackage{booktabs,graphicx}",
        r"\begin{document}",
        r"\input{final_component_ablation.tex}",
        r"\end{document}",
    ]) + "\n", encoding="utf-8")
    if shutil.which("pdflatex"):
        subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", preview.name],
            cwd=TABLE_DIR,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        generated = preview.with_suffix(".pdf")
        generated.replace(TABLE_DIR / "final_component_ablation_preview.pdf")
    else:
        # The compute image intentionally has no TeX distribution. Render the
        # exact table values to a vector PDF preview without modifying the
        # publication-ready table* source.
        rows = []
        for variant in VARIANTS:
            rows.append([
                METHOD_NAMES[variant],
                *[
                    f"{lookup[(variant, pde)]['rmse_mean']:.2e} ± "
                    f"{lookup[(variant, pde)]['rmse_std']:.2e}"
                    for pde in PDES
                ],
            ])
        fig, ax = plt.subplots(figsize=(11.4, 2.35), constrained_layout=True)
        ax.axis("off")
        rendered = ax.table(
            cellText=rows,
            colLabels=["Method", *[PDE_NAMES[pde] for pde in PDES]],
            loc="center",
            cellLoc="center",
            colWidths=[0.22, 0.195, 0.195, 0.195, 0.195],
        )
        rendered.auto_set_font_size(False)
        rendered.set_fontsize(8)
        rendered.scale(1.0, 1.55)
        fig.savefig(
            TABLE_DIR / "final_component_ablation_preview.pdf",
            bbox_inches="tight",
        )
        plt.close(fig)


def _mcv_inference_table(summary: list[dict[str, Any]]) -> None:
    lookup = {row["pde"]: row for row in summary}
    lines = [
        r"\begin{table*}[t]",
        r"\centering",
        r"\caption{Inference-only martingale-control-variate diagnostic at $d=100$. Both rows use the same frozen Full Tangent-Actor and identical Brownian paths; only the final estimator changes. Entries are whole-function value RMSE (mean $\pm$ sample standard deviation across three training seeds).}",
        r"\label{tab:mcv-inference-only}",
        r"\resizebox{\textwidth}{!}{%",
        r"\begin{tabular}{lcccc}",
        r"\toprule",
        r"Method & Radial Log-Terminal HJB & Concave Cubic--Quadratic & Convex Cubic--Quadratic & Oscillatory Value--Gradient \\",
        r"\midrule",
    ]
    for label, mean_key, std_key in (
        ("Tangent-Actor (Full; FK--MCV)", "full_mcv_rmse_mean", "full_mcv_rmse_std"),
        (
            "Full TA w/o MCV Inference",
            "full_without_mcv_inference_rmse_mean",
            "full_without_mcv_inference_rmse_std",
        ),
    ):
        values = [
            _format_value(float(lookup[pde][mean_key]), float(lookup[pde][std_key]))
            for pde in PDES
        ]
        lines.append(label + " & " + " & ".join(values) + r" \\")
    lines.extend([
        r"\bottomrule",
        r"\end{tabular}%",
        r"}",
        r"\end{table*}",
    ])
    (TABLE_DIR / "mcv_inference_only.tex").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def _component_figure(summary: list[dict[str, Any]]) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    lookup = {(row["variant"], row["pde"]): row for row in summary}
    x = np.arange(len(PDES), dtype=float)
    width = 0.19
    fig, ax = plt.subplots(figsize=(7.15, 3.6), constrained_layout=True)
    for index, variant in enumerate(VARIANTS):
        means = np.asarray([lookup[(variant, pde)]["rmse_mean"] for pde in PDES])
        stds = np.asarray([lookup[(variant, pde)]["rmse_std"] for pde in PDES])
        ax.bar(
            x + (index - 1.5) * width,
            means,
            width,
            yerr=stds,
            capsize=2.2,
            color=COLORS[variant],
            edgecolor="white",
            linewidth=0.5,
            label=METHOD_NAMES[variant],
        )
    values = [float(row["rmse_mean"]) for row in summary]
    if max(values) / min(values) > 20:
        ax.set_yscale("log")
    ax.set_ylabel("Whole-function value RMSE")
    ax.set_xticks(x)
    ax.set_xticklabels(["Log HJB", "Concave", "Convex", "Oscillatory"])
    ax.grid(axis="y", which="both", alpha=0.22, linewidth=0.6)
    ax.legend(ncol=2, fontsize=8, frameon=False, loc="upper left")
    for suffix, dpi in (("pdf", None), ("png", 300)):
        fig.savefig(FIGURE_DIR / f"component_ablation.{suffix}", dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _iteration_figure(rows: list[dict[str, Any]]) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(7.15, 2.7), constrained_layout=True, sharex=True)
    for ax, pde in zip(axes, ("joint_concave", "joint_convex", "nonconvex_coupled"), strict=True):
        local = [row for row in rows if row["pde"] == pde]
        seeds = sorted({int(row["seed"]) for row in local})
        for seed in seeds:
            seed_rows = sorted(
                (row for row in local if int(row["seed"]) == seed),
                key=lambda row: int(row["iteration"]),
            )
            ax.plot(
                [row["iteration"] for row in seed_rows],
                [row["value_rmse"] for row in seed_rows],
                marker="o",
                markersize=3.2,
                linewidth=1.2,
                alpha=0.9,
                label=f"seed {seed}",
            )
        ax.set_title({
            "joint_concave": "Concave",
            "joint_convex": "Convex",
            "nonconvex_coupled": "Oscillatory",
        }[pde])
        ax.set_yscale("log")
        ax.set_xticks(range(5))
        ax.set_xlabel("Actor update")
        ax.grid(which="both", alpha=0.22, linewidth=0.6)
        if len(seeds) > 1:
            ax.legend(fontsize=7, frameon=False)
    axes[0].set_ylabel("Value RMSE")
    for suffix, dpi in (("pdf", None), ("png", 300)):
        fig.savefig(FIGURE_DIR / f"actor_iteration_ablation.{suffix}", dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _resource_accounting(job_ids: list[str]) -> tuple[list[dict[str, Any]], float]:
    rows = _read_csv(RESULTS_DIR / "RESOURCE_ACCOUNTING.csv")
    # Keep the builder idempotent: a render retry must replace, rather than
    # duplicate, accounting rows for the same formal allocations.
    requested = set(job_ids)
    rows = [
        row
        for row in rows
        if not (
            row.get("workload") == "slurm_allocation"
            and row.get("job_id") in requested
        )
    ]
    actual_total = 0.0
    if job_ids:
        output = subprocess.check_output(
            [
                "sacct",
                "-X",
                "-n",
                "-P",
                "-j",
                ",".join(job_ids),
                "--format=JobIDRaw,JobName,Partition,State,ElapsedRaw,AllocTRES",
            ],
            text=True,
        )
        for line in output.splitlines():
            fields = line.split("|")
            if len(fields) < 6 or fields[0] not in job_ids:
                continue
            elapsed = float(fields[4] or 0.0)
            tres = fields[5]
            gpu = 0
            for token in tres.split(","):
                if token.startswith("gres/gpu="):
                    gpu = int(token.split("=", 1)[1])
            hours = elapsed * gpu / 3600.0
            actual_total += hours
            rows.append({
                "workload": "slurm_allocation",
                "pde": "all",
                "variant": "all",
                "cells": "",
                "estimated_allocated_gpu_hours": "",
                "actual_allocated_gpu_hours": hours,
                "job_id": fields[0],
                "job_name": fields[1],
                "partition": fields[2],
                "state": fields[3],
                "elapsed_seconds": elapsed,
                "allocated_gpus": gpu,
            })
    _write_csv(RESULTS_DIR / "RESOURCE_ACCOUNTING.csv", rows)
    return rows, actual_total


def _report(
    summary: list[dict[str, Any]],
    paired: list[dict[str, Any]],
    rows: list[dict[str, Any]],
    mcv_inference: list[dict[str, Any]],
    actual_gpu_hours: float,
) -> None:
    lookup = {(row["variant"], row["pde"]): row for row in summary}
    reuse = json.loads((RESULTS_DIR / "REUSE_AUDIT.json").read_text(encoding="utf-8"))
    preflight = json.loads((RESULTS_DIR / "PREFLIGHT_AUDIT.json").read_text(encoding="utf-8"))
    lines = [
        "# Final Component Ablation Study for Tangent-Actor",
        "",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Experimental definition and provenance",
        "",
        "The preregistered 4 PDE × 4 method study used the mature 100D budget, four tangent updates, seeds 0/1/2, a frozen 960-state bank, 64 Euler steps, and 1,024 paths per state. Replicate 0 is primary and replicate 1 is an independent inference repeat. No test reference entered fitting, noisy-label checkpointing, or reference-free selection. The convex working solution and gradient were sign-converted back to the original physical variable before metrics.",
        "",
        f"Exact Full-TA actor-training reuse: {reuse['full_actor_training_cells_reused']} cells; newly trained Full-TA: {reuse['full_actor_training_cells_new']} cells. A1 and A2 each contain 12 newly trained cells; A3 has no training.",
        "",
        "The historical oscillatory bank had 512 states, so its old RMSE was not reused. The preregistered final-suite generation rule supplied one new frozen 960-state bank shared by all variants. Full actors were re-evaluated on it.",
        "",
        "## Primary RMSE table",
        "",
        "| Method | Radial Log-Terminal HJB | Concave Cubic–Quadratic | Convex Cubic–Quadratic | Oscillatory Value–Gradient |",
        "|---|---:|---:|---:|---:|",
    ]
    for variant in VARIANTS:
        cells = [lookup[(variant, pde)] for pde in PDES]
        lines.append(
            "| " + METHOD_NAMES[variant] + " | "
            + " | ".join(f"{cell['rmse_mean']:.3e} ± {cell['rmse_std']:.3e}" for cell in cells)
            + " |"
        )
    lines.extend([
        "",
        "For A0–A2, ± is the sample standard deviation over the three primary training-seed RMSEs. A3 is a deterministic terminal-extension actor: its center and ± summarize two independent inference-replicate RMSEs, not fabricated training seeds. The replicate predictions remain separate and were never averaged pathwise.",
        "",
        "## Component interpretation",
        "",
    ])
    for variant, label in (("a1_no_mcv", "MCV"), ("a2_no_terminal_anchor", "Terminal anchoring")):
        lines.append(f"### {label}")
        lines.append("")
        for pde in PDES:
            differences = [float(row["ablated_minus_full_rmse"]) for row in paired if row["variant"] == variant and row["pde"] == pde]
            ratios = [float(row["rmse_ratio_ablated_over_full"]) for row in paired if row["variant"] == variant and row["pde"] == pde]
            lines.append(
                f"- {PDE_NAMES[pde]}: paired ablated-minus-full RMSE {np.mean(differences):+.3e} on average; mean ratio {np.mean(ratios):.2f}×."
            )
        lines.append("")
    lines.extend([
        "### Iterative tangent learning",
        "",
    ])
    for pde in PDES:
        full = lookup[("a0_full", pde)]
        fixed = lookup[("a3_no_updates", pde)]
        lines.append(
            f"- {PDE_NAMES[pde]}: fixed-A0 / Full mean-RMSE ratio {fixed['rmse_mean'] / full['rmse_mean']:.2f}×. A3 remains a fixed linearization and is not FH-MLP."
        )
    lines.extend([
        "",
        "## Inference-only MCV diagnostic",
        "",
        "`Full TA w/o MCV Inference` is a supplemental diagnostic, not a fifth primary method. It freezes the Full TA Actor trained and selected with MCV, then replaces only the final FK–MCV estimator by plain FK on identical Brownian paths. In contrast, primary `w/o MCV` removes MCV from labels, reference-free selection, and inference.",
        "",
        "| PDE | Full FK–MCV RMSE | Full w/o MCV inference RMSE | Plain / MCV |",
        "|---|---:|---:|---:|",
    ])
    for row in mcv_inference:
        lines.append(
            f"| {row['pde_name']} | {row['full_mcv_rmse_mean']:.3e} ± {row['full_mcv_rmse_std']:.3e} | "
            f"{row['full_without_mcv_inference_rmse_mean']:.3e} ± {row['full_without_mcv_inference_rmse_std']:.3e} | "
            f"{row['plain_over_mcv_rmse_ratio_mean']:.2f}× |"
        )
    lines.extend([
        "",
        "Seed- and replicate-level paired values are in `MCV_INFERENCE_ONLY_SEED_METRICS.csv`; the separate summary is in `MCV_INFERENCE_ONLY_SUMMARY.csv`. This comparison isolates inference variance reduction without conflating it with Actor-learning changes.",
        "",
        "## Costs, failures, and numerical resolution",
        "",
        f"Actual Slurm allocated GPU-hours across the recorded jobs: {actual_gpu_hours:.3f}. The preregistered hard ceiling was 36. No variant-specific stabilization or post-test tuning was permitted.",
        "",
        f"Numerical failures recorded in the 16 summary cells: {sum(int(row['numerical_failures']) for row in summary)}. Minimum exponential-weight ESS fraction: {min(float(row['weight_ess_fraction_min']) for row in summary):.4f}.",
        "",
        "Reference uncertainty must qualify differences near the independent reference floor. The exact reference metadata for every seed is retained in `ABLATION_SEED_METRICS.csv`; the preflight reference audits are in `PREFLIGHT_AUDIT.json`. Small apparent rank changes comparable to those uncertainties or the independent inference-replicate fluctuation are treated as unresolved.",
        "",
        "## Limitations",
        "",
        "Three training seeds are sufficient for an ablation signal but not a universal variance claim. The oscillatory a0–a4 supporting curve uses the audited historical 512-state bank, while its primary 4×4 cells use the newly frozen 960-state bank; no values are mixed within a curve. Joint iteration curves have only the audited seed-0 evidence where additional seeds were unavailable. These ablations neither establish universally lower variance than FH-MLP nor identify fixed A0 with FH-MLP.",
        "",
        "## Artifacts",
        "",
        "- `results/final_component_ablation/ABLATION_SEED_METRICS.csv`",
        "- `results/final_component_ablation/ABLATION_SUMMARY.csv`",
        "- `results/final_component_ablation/MCV_INFERENCE_ONLY_SUMMARY.csv`",
        "- `results/final_component_ablation/ACTOR_ITERATION_SUMMARY.csv`",
        "- `results/final_component_ablation/tables/final_component_ablation.tex`",
        "- `results/final_component_ablation/tables/mcv_inference_only.tex`",
        "- `results/final_component_ablation/figures/component_ablation.pdf`",
        "- `results/final_component_ablation/figures/actor_iteration_ablation.pdf`",
    ])
    docs = PROJECT_DIR / "docs" / "FINAL_COMPONENT_ABLATION_REPORT.md"
    docs.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _artifact_audit(expected_evaluations: int) -> None:
    required = [
        RESULTS_DIR / "PREREGISTRATION.md",
        RESULTS_DIR / "FROZEN_PROTOCOL.json",
        RESULTS_DIR / "REUSE_AUDIT.json",
        RESULTS_DIR / "SOURCE_FIDELITY_AUDIT.json",
        RESULTS_DIR / "ABLATION_SEED_METRICS.csv",
        RESULTS_DIR / "ABLATION_SUMMARY.csv",
        RESULTS_DIR / "MCV_INFERENCE_ONLY_SEED_METRICS.csv",
        RESULTS_DIR / "MCV_INFERENCE_ONLY_SUMMARY.csv",
        RESULTS_DIR / "ACTOR_ITERATION_SUMMARY.csv",
        RESULTS_DIR / "RESOURCE_ACCOUNTING.csv",
        TABLE_DIR / "final_component_ablation.tex",
        TABLE_DIR / "mcv_inference_only.tex",
        TABLE_DIR / "final_component_ablation_preview.pdf",
        FIGURE_DIR / "component_ablation.pdf",
        FIGURE_DIR / "component_ablation.png",
        FIGURE_DIR / "actor_iteration_ablation.pdf",
        FIGURE_DIR / "actor_iteration_ablation.png",
        PROJECT_DIR / "docs" / "FINAL_COMPONENT_ABLATION_REPORT.md",
    ]
    rows = _evaluation_rows()
    audit = {
        "required_files": [
            {"path": str(path), "exists": path.exists(), "bytes": path.stat().st_size if path.exists() else None}
            for path in required
        ],
        "evaluation_records": len(rows),
        "expected_evaluation_records": expected_evaluations,
        "all_test_banks_960": all(int(row["states"]) == 960 for row in rows),
        "all_evaluations_finite": all(
            np.isfinite(float(row["value_rmse"])) and int(row["nonfinite_count"]) == 0
            for row in rows
        ),
        "ground_truth_post_selection_only": all(
            bool(row["ground_truth_loaded_only_after_checkpoint_selection"]) for row in rows
        ),
    }
    audit["passed"] = bool(
        all(item["exists"] and int(item["bytes"] or 0) > 0 for item in audit["required_files"])
        and audit["evaluation_records"] == expected_evaluations
        and audit["all_test_banks_960"]
        and audit["all_evaluations_finite"]
        and audit["ground_truth_post_selection_only"]
    )
    _write_json(RESULTS_DIR / "ARTIFACT_AUDIT.json", audit)
    if not audit["passed"]:
        raise RuntimeError(f"artifact audit failed: {audit}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--job-ids", default="")
    args = parser.parse_args()
    rows = _evaluation_rows()
    if len(rows) != _expected_count():
        raise RuntimeError(f"expected {_expected_count()} evaluations, found {len(rows)}")
    seed_rows = []
    for row in rows:
        seed_rows.append({
            key: value
            for key, value in row.items()
            if key not in ("reference",)
        } | {"reference_json": json.dumps(row["reference"], sort_keys=True)})
    _write_csv(RESULTS_DIR / "ABLATION_SEED_METRICS.csv", seed_rows)
    summary = _summary(rows)
    _write_csv(RESULTS_DIR / "ABLATION_SUMMARY.csv", summary)
    _write_json(RESULTS_DIR / "ABLATION_SUMMARY.json", {
        "rows": summary,
        "primary_metric": "whole-function value RMSE",
        "a3_variability": "two independent inference replicates, no training seeds",
    })
    paired = _paired(rows)
    _write_csv(RESULTS_DIR / "PAIRED_ABLATION_DIFFERENCES.csv", paired)
    mcv_seed_rows, mcv_inference = _mcv_inference_only(rows)
    _write_csv(RESULTS_DIR / "MCV_INFERENCE_ONLY_SEED_METRICS.csv", mcv_seed_rows)
    _write_csv(RESULTS_DIR / "MCV_INFERENCE_ONLY_SUMMARY.csv", mcv_inference)
    iterations = _iteration_rows()
    _write_csv(RESULTS_DIR / "ACTOR_ITERATION_SUMMARY.csv", iterations)
    _, actual = _resource_accounting([value for value in args.job_ids.split(",") if value])
    _table(summary)
    _mcv_inference_table(mcv_inference)
    _component_figure(summary)
    _iteration_figure(iterations)
    _report(summary, paired, rows, mcv_inference, actual)
    _artifact_audit(_expected_count())
    print(json.dumps({
        "evaluations": len(rows),
        "summary_cells": len(summary),
        "actual_allocated_gpu_hours": actual,
        "artifact_audit": "passed",
    }, sort_keys=True))


if __name__ == "__main__":
    main()
