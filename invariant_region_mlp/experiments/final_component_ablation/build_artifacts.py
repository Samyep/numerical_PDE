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
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, NullFormatter
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
    "a0_full": "#2166ac",
    "a1_no_mcv": "#ef8a62",
    "a2_no_terminal_anchor": "#67a9cf",
    "a3_no_updates": "#b2182b",
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
    """Cross training-label and inference MCV choices on paired actors/paths."""
    seed_rows: list[dict[str, Any]] = []
    full_rows = [row for row in rows if row["variant"] == "a0_full"]
    no_mcv_lookup = {
        (row["pde"], int(row["seed"]), int(row["inference_replicate"])): row
        for row in rows
        if row["variant"] == "a1_no_mcv"
    }
    for row in full_rows:
        key = (
            row["pde"],
            int(row["seed"]),
            int(row["inference_replicate"]),
        )
        no_mcv = no_mcv_lookup[key]
        mcv = float(row["mcv_rmse_same_paths"])
        plain = float(row["plain_rmse_same_paths"])
        without_train = float(no_mcv["mcv_rmse_same_paths"])
        without_all = float(no_mcv["plain_rmse_same_paths"])
        seed_rows.append({
            "pde": row["pde"],
            "seed": int(row["seed"]),
            "inference_replicate": int(row["inference_replicate"]),
            "primary_replicate": bool(row["primary_replicate"]),
            "full_mcv_rmse": mcv,
            "full_without_mcv_inference_rmse": plain,
            "without_mcv_training_rmse": without_train,
            "without_all_mcv_rmse": without_all,
            "plain_minus_mcv_rmse": plain - mcv,
            "plain_over_mcv_rmse_ratio": plain / mcv,
            "without_training_over_full_ratio": without_train / mcv,
            "without_all_over_full_ratio": without_all / mcv,
            "same_actor": True,
            "same_brownian_paths": True,
            "training_and_selection_used_mcv": True,
            "without_training_actor_used_plain_labels_and_selector": True,
            "without_training_inference_used_mcv": True,
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
            "without_mcv_training_rmse_mean": _mean(
                row["without_mcv_training_rmse"] for row in primary
            ),
            "without_mcv_training_rmse_std": _sample_std(
                row["without_mcv_training_rmse"] for row in primary
            ),
            "without_all_mcv_rmse_mean": _mean(
                row["without_all_mcv_rmse"] for row in primary
            ),
            "without_all_mcv_rmse_std": _sample_std(
                row["without_all_mcv_rmse"] for row in primary
            ),
            "paired_plain_minus_mcv_rmse_mean": _mean(
                row["plain_minus_mcv_rmse"] for row in primary
            ),
            "plain_over_mcv_rmse_ratio_mean": _mean(
                row["plain_over_mcv_rmse_ratio"] for row in primary
            ),
            "plain_over_mcv_rmse_ratio_std": _sample_std(
                row["plain_over_mcv_rmse_ratio"] for row in primary
            ),
            "without_training_over_full_ratio_mean": _mean(
                row["without_training_over_full_ratio"] for row in primary
            ),
            "without_training_over_full_ratio_std": _sample_std(
                row["without_training_over_full_ratio"] for row in primary
            ),
            "without_all_over_full_ratio_mean": _mean(
                row["without_all_over_full_ratio"] for row in primary
            ),
            "without_all_over_full_ratio_std": _sample_std(
                row["without_all_over_full_ratio"] for row in primary
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
            "independent_without_training_rmse_mean": _mean(
                row["without_mcv_training_rmse"] for row in repeat
            ),
            "independent_without_all_rmse_mean": _mean(
                row["without_all_mcv_rmse"] for row in repeat
            ),
        })
    return seed_rows, summary


def _iteration_rows() -> list[dict[str, Any]]:
    cpu_paths = sorted(
        (RESULTS_DIR / "raw").glob(
            "iteration_cpu_diagnostic_*_d100_a0_full_s*_a*_rep0.json"
        )
    )
    gpu_paths = sorted(
        (RESULTS_DIR / "raw").glob(
            "iteration_diagnostic_*_d100_a0_full_s*_a*_rep0.json"
        )
    )
    formal_paths = cpu_paths or gpu_paths
    formal_role = (
        "formal_paired_cpu_iteration_diagnostic"
        if cpu_paths
        else "formal_paired_gpu_iteration_diagnostic"
    )
    formal_rows: list[dict[str, Any]] = []
    for path in formal_paths:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("status") != "complete":
            continue
        formal_rows.append({
            "pde": payload["pde"],
            "seed": int(payload["seed"]),
            "iteration": int(payload["iteration"]),
            "value_rmse": float(payload["value_rmse"]),
            "steps": int(payload["steps"]),
            "paths": int(payload["paths"]),
            "states": int(payload["states"]),
            "selected_iteration": int(payload["selected_iteration"]),
            "is_selected_full_actor": bool(payload["is_selected_full_actor"]),
            "source_artifact": str(path),
            "evidence_role": formal_role,
            "device": payload.get(
                "device", "cuda"
            ),
            "bank_limitation": "960-state frozen bank",
            "common_brownian_seed": int(payload["common_brownian_seed"]),
            "reused_formal_primary_evaluation": bool(
                payload["reused_formal_primary_evaluation"]
            ),
        })
    if formal_rows:
        keys = {
            (row["pde"], row["seed"], row["iteration"])
            for row in formal_rows
        }
        expected = {
            (pde, seed, iteration)
            for pde in PDES
            for seed in (0, 1, 2)
            for iteration in range(5)
        }
        if keys != expected:
            missing = sorted(expected - keys)
            raise RuntimeError(
                "incomplete formal iteration diagnostics: "
                f"{len(keys)}/60 present; first missing={missing[:5]}"
            )
        return sorted(
            formal_rows,
            key=lambda row: (row["pde"], row["seed"], row["iteration"]),
        )

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
        r"\caption{Crossed martingale-control-variate diagnostic at $d=100$. Training MCV is crossed with inference MCV while keeping each frozen Actor and paired Brownian paths fixed within the relevant comparison. Entries are whole-function value RMSE (mean $\pm$ sample standard deviation across three training seeds).}",
        r"\label{tab:mcv-crossed}",
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
        (
            "Full TA w/o MCV Training",
            "without_mcv_training_rmse_mean",
            "without_mcv_training_rmse_std",
        ),
        (
            "Full TA w/o All MCV",
            "without_all_mcv_rmse_mean",
            "without_all_mcv_rmse_std",
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


def _component_forest_heatmap_figure(
    summary: list[dict[str, Any]], rows: list[dict[str, Any]]
) -> None:
    """Retain the compact forest-plus-heatmap layout as an alternate renderer."""
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    lookup = {(row["variant"], row["pde"]): row for row in summary}
    plot_order = (
        "a0_full",
        "a2_no_terminal_anchor",
        "a1_no_mcv",
        "a3_no_updates",
    )
    method_labels = {
        "a0_full": "Full TA",
        "a2_no_terminal_anchor": "w/o terminal anchor",
        "a1_no_mcv": "w/o MCV",
        "a3_no_updates": "w/o tangent updates",
    }
    short_pde_names = {
        "hjb_log": "Log HJB",
        "joint_concave": "Concave",
        "joint_convex": "Convex",
        "nonconvex_coupled": "Oscillatory",
    }
    markers = {
        "a0_full": "D",
        "a2_no_terminal_anchor": "o",
        "a1_no_mcv": "s",
        "a3_no_updates": "X",
    }
    positions = {variant: 3 - index for index, variant in enumerate(plot_order)}
    paper_style = {
        "font.family": "serif",
        "font.serif": ["DejaVu Serif"],
        "mathtext.fontset": "dejavuserif",
        "axes.linewidth": 0.7,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
    }
    with plt.rc_context(paper_style):
        fig = plt.figure(figsize=(7.25, 5.05), constrained_layout=True)
        grid = fig.add_gridspec(2, 4, height_ratios=(2.2, 1.25))
        forest_axes = [fig.add_subplot(grid[0, index]) for index in range(4)]

        for pde_index, (ax, pde) in enumerate(
            zip(forest_axes, PDES, strict=True)
        ):
            # Faint lines retain the paired-seed structure for trained variants.
            for seed, jitter in zip((0, 1, 2), (-0.09, 0.0, 0.09), strict=True):
                paired_values = []
                paired_positions = []
                for variant in ("a0_full", "a2_no_terminal_anchor", "a1_no_mcv"):
                    local = [
                        row
                        for row in rows
                        if row["pde"] == pde
                        and row["variant"] == variant
                        and int(row["seed"]) == seed
                        and row["primary_replicate"]
                    ]
                    if len(local) != 1:
                        raise RuntimeError(
                            f"missing paired seed point: {pde} {variant} seed={seed}"
                        )
                    paired_values.append(float(local[0]["value_rmse"]))
                    paired_positions.append(positions[variant] + jitter)
                ax.plot(
                    paired_values,
                    paired_positions,
                    color="#6b6b6b",
                    alpha=0.16,
                    linewidth=0.55,
                    zorder=1,
                )

            for variant in plot_order:
                center = lookup[(variant, pde)]
                local = [
                    row
                    for row in rows
                    if row["pde"] == pde
                    and row["variant"] == variant
                    and (
                        variant == "a3_no_updates" or row["primary_replicate"]
                    )
                ]
                local = sorted(
                    local,
                    key=lambda row: (
                        int(row["seed"]), int(row["inference_replicate"])
                    ),
                )
                jitter = np.linspace(-0.09, 0.09, len(local))
                ax.scatter(
                    [float(row["value_rmse"]) for row in local],
                    positions[variant] + jitter,
                    s=15,
                    color=COLORS[variant],
                    alpha=0.38,
                    edgecolors="none",
                    zorder=2,
                )
                ax.errorbar(
                    float(center["rmse_mean"]),
                    positions[variant],
                    xerr=float(center["rmse_std"]),
                    fmt=markers[variant],
                    markersize=5.0,
                    markerfacecolor=COLORS[variant],
                    markeredgecolor="white",
                    markeredgewidth=0.55,
                    color=COLORS[variant],
                    ecolor=COLORS[variant],
                    elinewidth=1.25,
                    capsize=2.0,
                    zorder=4,
                )

            ax.axvline(
                float(lookup[("a0_full", pde)]["rmse_mean"]),
                color=COLORS["a0_full"],
                linewidth=0.75,
                linestyle=(0, (2, 2)),
                alpha=0.28,
                zorder=0,
            )
            ax.set_xscale("log")
            displayed_means = np.asarray(
                [float(lookup[(variant, pde)]["rmse_mean"]) for variant in plot_order]
            )
            lower = float(np.min(displayed_means))
            upper = float(np.max(displayed_means))
            ax.set_xlim(lower / 1.22, upper * 1.22)
            ax.set_xticks((lower, upper))
            ax.set_xticklabels((f"{lower:.1e}", f"{upper:.1e}"), fontsize=6.7)
            ax.xaxis.set_minor_locator(
                LogLocator(base=10, subs=np.arange(2, 10) * 0.1, numticks=24)
            )
            ax.xaxis.set_minor_formatter(NullFormatter())
            ax.set_ylim(-0.55, 3.55)
            ax.set_title(short_pde_names[pde], fontsize=10.0, fontweight="bold")
            ax.grid(axis="x", which="both", color="#c7c7c7", alpha=0.28, linewidth=0.55)
            ax.tick_params(axis="x", which="major", length=3.0, pad=2.0)
            ax.tick_params(axis="x", which="minor", length=1.8)
            ax.tick_params(axis="y", length=0)
            ax.set_xlabel("value RMSE", fontsize=6.8, labelpad=3.0)
            if pde_index == 0:
                ax.set_yticks(
                    [positions[variant] for variant in plot_order],
                    [method_labels[variant] for variant in plot_order],
                    fontsize=7.2,
                )
            else:
                ax.set_yticks([positions[variant] for variant in plot_order])
                ax.set_yticklabels([])
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)

        forest_axes[0].set_ylabel("Component configuration", fontsize=8.0)
        fig.text(
            0.012,
            0.985,
            "a",
            ha="left",
            va="top",
            fontsize=11,
            fontweight="bold",
        )

        heat_ax = fig.add_subplot(grid[1, :])
        heat_variants = (
            "a2_no_terminal_anchor",
            "a1_no_mcv",
            "a3_no_updates",
        )
        ratios = np.asarray(
            [
                [
                    float(lookup[(variant, pde)]["rmse_mean"])
                    / float(lookup[("a0_full", pde)]["rmse_mean"])
                    for pde in PDES
                ]
                for variant in heat_variants
            ]
        )
        log_ratios = np.log10(ratios)
        span = max(float(np.max(np.abs(log_ratios))), 0.08)
        cmap = LinearSegmentedColormap.from_list(
            "tangent_actor_rdbu", ("#2166ac", "#f7f7f7", "#b2182b"), N=256
        )
        norm = TwoSlopeNorm(vmin=-span, vcenter=0.0, vmax=span)
        image = heat_ax.imshow(log_ratios, cmap=cmap, norm=norm, aspect="auto")
        heat_ax.set_xticks(
            range(len(PDES)),
            [short_pde_names[pde] for pde in PDES],
            fontsize=8.0,
        )
        heat_ax.set_yticks(
            range(len(heat_variants)),
            [
                "w/o terminal anchor",
                "w/o MCV",
                "w/o tangent updates",
            ],
            fontsize=7.6,
        )
        heat_ax.set_title(
            "Ablation penalty relative to Full TA",
            fontsize=9.3,
            fontweight="bold",
            pad=5,
        )
        heat_ax.set_xticks(np.arange(-0.5, len(PDES), 1), minor=True)
        heat_ax.set_yticks(np.arange(-0.5, len(heat_variants), 1), minor=True)
        heat_ax.grid(which="minor", color="white", linewidth=1.2)
        heat_ax.tick_params(which="minor", bottom=False, left=False)
        for row_index, variant in enumerate(heat_variants):
            for column_index, pde in enumerate(PDES):
                ratio = ratios[row_index, column_index]
                delta = (
                    float(lookup[(variant, pde)]["rmse_mean"])
                    - float(lookup[("a0_full", pde)]["rmse_mean"])
                )
                text_color = "white" if log_ratios[row_index, column_index] > 0.34 else "#222222"
                heat_ax.text(
                    column_index,
                    row_index,
                    f"{ratio:.2f}×\n({delta:+.1e})",
                    ha="center",
                    va="center",
                    fontsize=7.2,
                    fontweight="bold" if ratio >= 2 else "normal",
                    color=text_color,
                )
        colorbar = fig.colorbar(image, ax=heat_ax, fraction=0.025, pad=0.018, aspect=17)
        colorbar.set_ticks((-span, 0.0, span))
        colorbar.set_ticklabels(
            (f"{10 ** (-span):.2f}×", "1×", f"{10 ** span:.0f}×")
        )
        colorbar.set_label("RMSE / Full TA", fontsize=7.0)
        colorbar.ax.tick_params(labelsize=6.5, width=0.6)
        heat_ax.text(
            -0.085,
            1.11,
            "b",
            transform=heat_ax.transAxes,
            fontsize=11,
            fontweight="bold",
            va="top",
        )
        fig.suptitle(
            "Tangent-Actor component ablation",
            fontsize=11.5,
            fontweight="bold",
        )
        fig.text(
            0.5,
            -0.012,
            "Small points: training seeds (A0–A2); fixed-$A_0$ uses two inference replicates. "
            "Heatmap annotations: RMSE ratio and absolute change.",
            ha="center",
            va="top",
            fontsize=6.4,
            color="#4a4a4a",
        )
        for suffix, dpi in (("pdf", None), ("png", 300)):
            fig.savefig(
                FIGURE_DIR / f"component_ablation.{suffix}",
                dpi=dpi,
                bbox_inches="tight",
                facecolor="white",
            )
        plt.close(fig)


def _component_matrix_horizontal_figure(
    summary: list[dict[str, Any]],
    rows: list[dict[str, Any]],
    mcv_inference: list[dict[str, Any]],
) -> None:
    """Retain the horizontal diagnostic matrix as an alternate renderer."""
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    lookup = {(row["variant"], row["pde"]): row for row in summary}
    mcv_lookup = {row["pde"]: row for row in mcv_inference}
    plot_order = (
        "a0_full",
        "a2_no_terminal_anchor",
        "a1_no_mcv",
        "a3_no_updates",
    )
    markers = {
        "a0_full": "D",
        "a2_no_terminal_anchor": "o",
        "a1_no_mcv": "s",
        "a3_no_updates": "X",
    }
    short_method_names = {
        "a0_full": "Full TA",
        "a2_no_terminal_anchor": "w/o anchor",
        "a1_no_mcv": "w/o MCV",
        "a3_no_updates": "w/o updates",
    }
    row_names = {
        "hjb_log": "Log HJB",
        "joint_concave": "Joint concave",
        "joint_convex": "Joint convex",
        "nonconvex_coupled": "Nonconvex",
    }
    diverging = LinearSegmentedColormap.from_list(
        "tangent_actor_surface",
        ("#2166ac", "#d1e5f0", "#f7f7f7", "#fddbc7", "#b2182b"),
        N=256,
    )
    paper_style = {
        "font.family": "serif",
        "font.serif": ["DejaVu Serif"],
        "mathtext.fontset": "dejavuserif",
        "axes.linewidth": 0.55,
        "xtick.major.width": 0.55,
        "ytick.major.width": 0.55,
    }

    def local_rows(pde: str, variant: str) -> list[dict[str, Any]]:
        selected = [
            row
            for row in rows
            if row["pde"] == pde
            and row["variant"] == variant
            and (variant == "a3_no_updates" or row["primary_replicate"])
        ]
        return sorted(
            selected,
            key=lambda row: (int(row["seed"]), int(row["inference_replicate"])),
        )

    def add_absolute_gradient(ax: plt.Axes, lower: float, upper: float) -> None:
        edges = np.geomspace(lower, upper, 181)
        color_values = np.linspace(-1.0, 1.0, len(edges) - 1)[None, :]
        ax.pcolormesh(
            edges,
            (-0.34, 0.34),
            color_values,
            cmap=diverging,
            norm=Normalize(vmin=-1.0, vmax=1.0),
            shading="flat",
            alpha=0.16,
            zorder=0,
        )

    ratio_lower, ratio_upper = 0.82, 30.0
    ratio_norm = TwoSlopeNorm(
        vmin=math.log10(ratio_lower),
        vcenter=0.0,
        vmax=math.log10(ratio_upper),
    )

    def add_ratio_gradient(ax: plt.Axes) -> None:
        edges = np.geomspace(ratio_lower, ratio_upper, 181)
        centers = np.sqrt(edges[:-1] * edges[1:])
        ax.pcolormesh(
            edges,
            (-0.29, 0.29),
            np.log10(centers)[None, :],
            cmap=diverging,
            norm=ratio_norm,
            shading="flat",
            alpha=0.24,
            zorder=0,
        )
        ax.axvline(
            1.0,
            color=COLORS["a0_full"],
            linewidth=0.75,
            linestyle=(0, (2, 2)),
            alpha=0.72,
            zorder=1,
        )
        ax.scatter(
            [1.0],
            [0.0],
            s=17,
            marker="D",
            color=COLORS["a0_full"],
            edgecolor="white",
            linewidth=0.45,
            zorder=4,
        )

    def style_cell(ax: plt.Axes, *, bottom: bool, ratio_axis: bool) -> None:
        ax.set_xscale("log")
        ax.set_yticks([])
        ax.tick_params(axis="y", length=0)
        ax.tick_params(axis="x", which="major", labelsize=6.0, length=2.6, pad=1.5)
        ax.tick_params(axis="x", which="minor", length=1.3)
        ax.xaxis.set_minor_formatter(NullFormatter())
        if ratio_axis:
            ax.set_xlim(ratio_lower, ratio_upper)
            ax.set_xticks((1.0, 3.0, 10.0, 30.0))
            if bottom:
                ax.set_xticklabels(("1×", "3×", "10×", "30×"))
                ax.set_xlabel("RMSE / Full TA", fontsize=6.5, labelpad=2.0)
            else:
                ax.set_xticklabels([])
        for spine in ax.spines.values():
            spine.set_color("#5a5a5a")
            spine.set_linewidth(0.55)
        ax.grid(axis="x", which="major", color="white", alpha=0.70, linewidth=0.55)

    with plt.rc_context(paper_style):
        fig = plt.figure(figsize=(7.25, 5.28), facecolor="white")
        grid = fig.add_gridspec(
            4,
            4,
            left=0.155,
            right=0.985,
            bottom=0.12,
            top=0.765,
            hspace=0.27,
            wspace=0.23,
            width_ratios=(1.45, 1.0, 1.0, 1.0),
        )
        axes = np.empty((4, 4), dtype=object)
        column_titles = (
            "Absolute value error\nwhole-function RMSE",
            "Martingale CV\n○ inference · ■ train+infer.",
            "Terminal anchoring\nremoved",
            "Tangent updates\nfixed $A_0$",
        )

        for row_index, pde in enumerate(PDES):
            # Column 1: all four absolute errors with seed-level observations.
            absolute_ax = fig.add_subplot(grid[row_index, 0])
            axes[row_index, 0] = absolute_ax
            means = np.asarray(
                [float(lookup[(variant, pde)]["rmse_mean"]) for variant in plot_order]
            )
            lower, upper = float(np.min(means) / 1.30), float(np.max(means) * 1.30)
            add_absolute_gradient(absolute_ax, lower, upper)
            y_offsets = (0.19, 0.065, -0.065, -0.19)
            for variant, y_value in zip(plot_order, y_offsets, strict=True):
                center = lookup[(variant, pde)]
                observations = local_rows(pde, variant)
                observation_jitter = np.linspace(-0.022, 0.022, len(observations))
                absolute_ax.scatter(
                    [float(row["value_rmse"]) for row in observations],
                    y_value + observation_jitter,
                    s=11,
                    marker=markers[variant],
                    color=COLORS[variant],
                    alpha=0.32,
                    edgecolors="none",
                    zorder=3,
                )
                absolute_ax.errorbar(
                    float(center["rmse_mean"]),
                    y_value,
                    xerr=float(center["rmse_std"]),
                    fmt=markers[variant],
                    markersize=4.2,
                    color=COLORS[variant],
                    ecolor=COLORS[variant],
                    elinewidth=1.0,
                    capsize=1.6,
                    markerfacecolor=COLORS[variant],
                    markeredgecolor="white",
                    markeredgewidth=0.4,
                    zorder=5,
                )
            absolute_ax.set_xlim(lower, upper)
            absolute_ax.set_ylim(-0.34, 0.34)
            absolute_ax.set_xscale("log")
            absolute_ax.set_yticks([])
            absolute_ax.set_xticks((float(np.min(means)), float(np.max(means))))
            absolute_ax.set_xticklabels(
                (f"{np.min(means):.1e}", f"{np.max(means):.1e}"), fontsize=6.0
            )
            absolute_ax.xaxis.set_minor_locator(
                LogLocator(base=10, subs=np.arange(2, 10) * 0.1, numticks=24)
            )
            absolute_ax.xaxis.set_minor_formatter(NullFormatter())
            absolute_ax.tick_params(axis="x", which="major", length=2.6, pad=1.5)
            absolute_ax.tick_params(axis="x", which="minor", length=1.3)
            absolute_ax.grid(
                axis="x", which="major", color="white", alpha=0.72, linewidth=0.55
            )
            for spine in absolute_ax.spines.values():
                spine.set_color("#5a5a5a")
                spine.set_linewidth(0.55)
            if row_index == len(PDES) - 1:
                absolute_ax.set_xlabel("value RMSE", fontsize=6.5, labelpad=2.0)
            absolute_ax.text(
                -0.30,
                0.5,
                row_names[pde],
                transform=absolute_ax.transAxes,
                ha="right",
                va="center",
                fontsize=8.2,
                fontweight="bold",
            )

            full_mean = float(lookup[("a0_full", pde)]["rmse_mean"])
            full_by_seed = {
                int(row["seed"]): float(row["value_rmse"])
                for row in local_rows(pde, "a0_full")
            }

            # Column 2: inference-only versus full MCV removal.
            mcv_ax = fig.add_subplot(grid[row_index, 1])
            axes[row_index, 1] = mcv_ax
            add_ratio_gradient(mcv_ax)
            full_removal_ratio = float(lookup[("a1_no_mcv", pde)]["rmse_mean"]) / full_mean
            inference_ratio = float(mcv_lookup[pde]["plain_over_mcv_rmse_ratio_mean"])
            no_mcv_by_seed = {
                int(row["seed"]): float(row["value_rmse"])
                for row in local_rows(pde, "a1_no_mcv")
            }
            seed_full_removal = [
                no_mcv_by_seed[seed] / full_by_seed[seed] for seed in sorted(full_by_seed)
            ]
            full_rows = local_rows(pde, "a0_full")
            seed_inference = [
                float(row["plain_rmse_same_paths"]) / float(row["mcv_rmse_same_paths"])
                for row in full_rows
            ]
            for ratio, y_value, marker, filled in (
                (inference_ratio, -0.10, "o", False),
                (full_removal_ratio, 0.10, "s", True),
            ):
                mcv_ax.plot(
                    (1.0, ratio),
                    (y_value, y_value),
                    color=COLORS["a1_no_mcv"],
                    linewidth=1.05,
                    alpha=0.82,
                    zorder=2,
                )
                mcv_ax.scatter(
                    [ratio],
                    [y_value],
                    s=24,
                    marker=marker,
                    facecolor=COLORS["a1_no_mcv"] if filled else "white",
                    edgecolor=COLORS["a1_no_mcv"],
                    linewidth=0.9,
                    zorder=5,
                )
            mcv_ax.scatter(
                seed_full_removal,
                0.10 + np.linspace(-0.027, 0.027, len(seed_full_removal)),
                s=9,
                color=COLORS["a1_no_mcv"],
                alpha=0.28,
                zorder=3,
            )
            mcv_ax.scatter(
                seed_inference,
                -0.10 + np.linspace(-0.027, 0.027, len(seed_inference)),
                s=9,
                facecolor="white",
                edgecolor=COLORS["a1_no_mcv"],
                linewidth=0.55,
                alpha=0.55,
                zorder=3,
            )
            mcv_ax.text(
                0.96,
                0.78,
                f"all  {full_removal_ratio:.2f}×",
                transform=mcv_ax.transAxes,
                ha="right",
                va="center",
                fontsize=6.2,
                fontweight="bold",
            )
            mcv_ax.text(
                0.96,
                0.19,
                f"infer.  {inference_ratio:.2f}×",
                transform=mcv_ax.transAxes,
                ha="right",
                va="center",
                fontsize=6.0,
            )
            style_cell(mcv_ax, bottom=row_index == len(PDES) - 1, ratio_axis=True)

            # Columns 3 and 4: single-component penalties.
            for column_index, variant in (
                (2, "a2_no_terminal_anchor"),
                (3, "a3_no_updates"),
            ):
                effect_ax = fig.add_subplot(grid[row_index, column_index])
                axes[row_index, column_index] = effect_ax
                add_ratio_gradient(effect_ax)
                ratio = float(lookup[(variant, pde)]["rmse_mean"]) / full_mean
                delta = float(lookup[(variant, pde)]["rmse_mean"]) - full_mean
                effect_ax.plot(
                    (1.0, ratio),
                    (0.0, 0.0),
                    color=COLORS[variant],
                    linewidth=1.2,
                    alpha=0.88,
                    zorder=2,
                )
                effect_ax.scatter(
                    [ratio],
                    [0.0],
                    s=26,
                    marker=markers[variant],
                    facecolor=COLORS[variant],
                    edgecolor="white",
                    linewidth=0.55,
                    zorder=5,
                )
                if variant == "a2_no_terminal_anchor":
                    ablated_by_seed = {
                        int(row["seed"]): float(row["value_rmse"])
                        for row in local_rows(pde, variant)
                    }
                    observed_ratios = [
                        ablated_by_seed[seed] / full_by_seed[seed]
                        for seed in sorted(full_by_seed)
                    ]
                else:
                    observed_ratios = [
                        float(row["value_rmse"]) / full_mean
                        for row in local_rows(pde, variant)
                    ]
                effect_ax.scatter(
                    observed_ratios,
                    np.linspace(-0.035, 0.035, len(observed_ratios)),
                    s=10,
                    marker=markers[variant],
                    color=COLORS[variant],
                    alpha=0.30,
                    edgecolors="none",
                    zorder=3,
                )
                effect_ax.text(
                    0.96,
                    0.73,
                    f"{ratio:.2f}×",
                    transform=effect_ax.transAxes,
                    ha="right",
                    va="center",
                    fontsize=6.8,
                    fontweight="bold",
                )
                effect_ax.text(
                    0.96,
                    0.22,
                    f"Δ {delta:+.1e}",
                    transform=effect_ax.transAxes,
                    ha="right",
                    va="center",
                    fontsize=5.9,
                    color="#3f3f3f",
                )
                style_cell(
                    effect_ax,
                    bottom=row_index == len(PDES) - 1,
                    ratio_axis=True,
                )

        for column_index, title in enumerate(column_titles):
            axes[0, column_index].set_title(
                title,
                fontsize=7.6,
                fontweight="bold",
                pad=8.0,
                linespacing=1.15,
            )

        legend_handles = [
            Line2D(
                [0],
                [0],
                marker=markers[variant],
                linestyle="none",
                markerfacecolor=COLORS[variant],
                markeredgecolor="white",
                markeredgewidth=0.45,
                markersize=5.2,
                label=short_method_names[variant],
            )
            for variant in plot_order
        ]
        fig.legend(
            handles=legend_handles,
            loc="upper center",
            bbox_to_anchor=(0.57, 0.855),
            ncol=4,
            frameon=False,
            fontsize=6.5,
            columnspacing=1.35,
            handletextpad=0.35,
        )
        fig.suptitle(
            "Tangent-Actor across four high-dimensional PDEs",
            x=0.57,
            y=0.965,
            fontsize=11.2,
            fontweight="bold",
        )
        fig.text(
            0.57,
            0.027,
            "Large markers: reported center; faint markers: individual training seeds "
            "($A_0$ uses inference replicates). Background: lower error blue, larger penalty red.",
            ha="center",
            va="bottom",
            fontsize=6.0,
            color="#4a4a4a",
        )
        for suffix, dpi in (("pdf", None), ("png", 300)):
            fig.savefig(
                FIGURE_DIR / f"component_ablation.{suffix}",
                dpi=dpi,
                bbox_inches="tight",
                facecolor="white",
            )
        plt.close(fig)


def _component_vertical_matrix_figure(
    summary: list[dict[str, Any]],
    rows: list[dict[str, Any]],
    mcv_inference: list[dict[str, Any]],
    iteration_rows: list[dict[str, Any]],
) -> None:
    """Render a vertically encoded diagnostic matrix with gradient fields."""
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    lookup = {(row["variant"], row["pde"]): row for row in summary}
    mcv_lookup = {row["pde"]: row for row in mcv_inference}
    plot_order = (
        "a0_full",
        "a2_no_terminal_anchor",
        "a1_no_mcv",
        "a3_no_updates",
    )
    markers = {
        "a0_full": "D",
        "a2_no_terminal_anchor": "o",
        "a1_no_mcv": "s",
        "a3_no_updates": "X",
    }
    absolute_bar_colors = {
        "a0_full": "#8CB9D7",
        "a2_no_terminal_anchor": "#B7D8E8",
        "a1_no_mcv": "#F3B18F",
        "a3_no_updates": "#D97B86",
    }
    row_names = {
        "hjb_log": "Log HJB",
        "joint_concave": "Concave",
        "joint_convex": "Convex",
        "nonconvex_coupled": "Nonconvex",
    }
    diverging = LinearSegmentedColormap.from_list(
        "tangent_actor_vertical_surface",
        ("#2166ac", "#d1e5f0", "#f7f7f7", "#fddbc7", "#b2182b"),
        N=256,
    )
    paper_style = {
        "font.family": "serif",
        "font.serif": ["DejaVu Serif"],
        "mathtext.fontset": "dejavuserif",
        "axes.linewidth": 0.55,
        "xtick.major.width": 0.55,
        "ytick.major.width": 0.55,
    }

    def local_rows(pde: str, variant: str) -> list[dict[str, Any]]:
        selected = [
            row
            for row in rows
            if row["pde"] == pde
            and row["variant"] == variant
            and (variant == "a3_no_updates" or row["primary_replicate"])
        ]
        return sorted(
            selected,
            key=lambda row: (int(row["seed"]), int(row["inference_replicate"])),
        )

    def absolute_background(ax: plt.Axes, upper: float) -> None:
        edges = np.linspace(0.0, upper, 181)
        ax.pcolormesh(
            (-0.55, 3.55),
            edges,
            np.linspace(-1.0, 1.0, len(edges) - 1)[:, None],
            cmap=diverging,
            norm=Normalize(vmin=-1.0, vmax=1.0),
            shading="flat",
            alpha=0.16,
            zorder=0,
        )

    ratio_lower, ratio_upper = 0.82, 30.0
    ratio_norm = TwoSlopeNorm(
        vmin=math.log10(ratio_lower),
        vcenter=0.0,
        vmax=math.log10(ratio_upper),
    )

    def ratio_background(ax: plt.Axes, x_limits: tuple[float, float]) -> None:
        edges = np.geomspace(ratio_lower, ratio_upper, 181)
        centers = np.sqrt(edges[:-1] * edges[1:])
        ax.pcolormesh(
            x_limits,
            edges,
            np.log10(centers)[:, None],
            cmap=diverging,
            norm=ratio_norm,
            shading="flat",
            alpha=0.24,
            zorder=0,
        )
        ax.axhline(
            1.0,
            color=COLORS["a0_full"],
            linewidth=0.75,
            linestyle=(0, (2, 2)),
            alpha=0.72,
            zorder=2,
        )

    def style_absolute_axis(
        ax: plt.Axes, upper: float, *, bottom: bool
    ) -> None:
        ax.set_xlim(-0.55, 3.55)
        ax.set_ylim(0.0, upper)
        ticks = tuple(np.linspace(0.0, upper, 5))
        ax.set_yticks(ticks)
        ax.set_yticklabels(
            ("0", "", "", "", f"{ticks[4]:.1e}"), fontsize=5.7
        )
        ax.set_xticks(range(4))
        if bottom:
            ax.set_xticklabels(
                ("Full", "w/o\nanchor", "w/o\nMCV", "w/o\nupdates"),
                fontsize=5.2,
                linespacing=1.0,
            )
        else:
            ax.set_xticklabels([])
        ax.tick_params(axis="x", length=0)
        ax.tick_params(axis="y", length=2.2, pad=1.6)
        ax.grid(
            axis="y",
            which="major",
            color="white",
            alpha=0.92,
            linewidth=0.68,
            linestyle="-",
            zorder=2,
        )
        for spine in ax.spines.values():
            spine.set_color("#5a5a5a")
            spine.set_linewidth(0.55)

    def style_ratio_axis(
        ax: plt.Axes,
        *,
        bottom: bool,
        show_y_labels: bool,
        x_positions: tuple[float, ...],
        x_labels: tuple[str, ...],
        x_limits: tuple[float, float],
    ) -> None:
        ax.set_xlim(*x_limits)
        ax.set_yscale("log")
        ax.set_ylim(ratio_lower, ratio_upper)
        ax.set_yticks((1.0, 3.0, 10.0, 30.0))
        if show_y_labels:
            ax.set_yticklabels(("1×", "3×", "10×", "30×"), fontsize=5.8)
        else:
            ax.set_yticklabels([])
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.tick_params(axis="y", which="major", length=2.2, pad=1.6)
        ax.tick_params(axis="y", which="minor", length=1.1)
        ax.set_xticks(x_positions)
        if bottom:
            ax.set_xticklabels(x_labels, fontsize=5.6)
        else:
            ax.set_xticklabels([])
        ax.tick_params(axis="x", length=0)
        ax.grid(axis="y", which="major", color="white", alpha=0.72, linewidth=0.55)
        for spine in ax.spines.values():
            spine.set_color("#5a5a5a")
            spine.set_linewidth(0.55)

    with plt.rc_context(paper_style):
        fig = plt.figure(figsize=(7.15, 4.55), facecolor="white")
        grid = fig.add_gridspec(
            4,
            4,
            left=0.13,
            right=0.99,
            bottom=0.11,
            top=0.98,
            hspace=0.22,
            wspace=0.25,
            width_ratios=(1.42, 1.0, 1.0, 1.0),
        )
        axes = np.empty((4, 4), dtype=object)

        for row_index, pde in enumerate(PDES):
            bottom_row = row_index == len(PDES) - 1

            # Column 1: vertical bars with seed-level uncertainty whiskers.
            absolute_ax = fig.add_subplot(grid[row_index, 0])
            axes[row_index, 0] = absolute_ax
            means = np.asarray(
                [float(lookup[(variant, pde)]["rmse_mean"]) for variant in plot_order]
            )
            stds = np.asarray(
                [float(lookup[(variant, pde)]["rmse_std"]) for variant in plot_order]
            )
            upper = float(np.max(means + stds) * 1.22)
            absolute_background(absolute_ax, upper)
            x_positions = np.arange(len(plot_order), dtype=float)
            absolute_ax.bar(
                x_positions,
                means,
                width=0.62,
                color=[absolute_bar_colors[variant] for variant in plot_order],
                alpha=0.96,
                edgecolor="white",
                linewidth=0.55,
                yerr=stds,
                capsize=2.2,
                error_kw={
                    "ecolor": "#292929",
                    "elinewidth": 0.8,
                    "capthick": 0.8,
                    "zorder": 6,
                },
                zorder=3,
            )
            style_absolute_axis(absolute_ax, upper, bottom=bottom_row)
            absolute_box = absolute_ax.get_position()
            fig.text(
                0.035,
                (absolute_box.y0 + absolute_box.y1) / 2.0,
                row_names[pde],
                ha="left",
                va="center",
                fontsize=7.4,
                fontweight="bold",
            )

            full_mean = float(lookup[("a0_full", pde)]["rmse_mean"])
            full_by_seed = {
                int(row["seed"]): float(row["value_rmse"])
                for row in local_rows(pde, "a0_full")
            }

            # Column 2: cross MCV removal in inference, training, or both.
            mcv_ax = fig.add_subplot(grid[row_index, 1])
            axes[row_index, 1] = mcv_ax
            mcv_x_limits = (-0.55, 2.55)
            ratio_background(mcv_ax, mcv_x_limits)
            inference_ratio = float(mcv_lookup[pde]["plain_over_mcv_rmse_ratio_mean"])
            training_ratio = float(
                mcv_lookup[pde]["without_training_over_full_ratio_mean"]
            )
            full_removal_ratio = float(
                mcv_lookup[pde]["without_all_over_full_ratio_mean"]
            )
            no_mcv_mcv_by_seed = {
                int(row["seed"]): float(row["mcv_rmse_same_paths"])
                for row in local_rows(pde, "a1_no_mcv")
            }
            no_mcv_plain_by_seed = {
                int(row["seed"]): float(row["plain_rmse_same_paths"])
                for row in local_rows(pde, "a1_no_mcv")
            }
            inference_seed_ratios = [
                float(row["plain_rmse_same_paths"]) / float(row["mcv_rmse_same_paths"])
                for row in local_rows(pde, "a0_full")
            ]
            training_seed_ratios = [
                no_mcv_mcv_by_seed[seed] / full_by_seed[seed]
                for seed in sorted(full_by_seed)
            ]
            removal_seed_ratios = [
                no_mcv_plain_by_seed[seed] / full_by_seed[seed]
                for seed in sorted(full_by_seed)
            ]
            for x_value, ratio, observed, facecolor in (
                (0.0, inference_ratio, inference_seed_ratios, "white"),
                (1.0, training_ratio, training_seed_ratios, "#fddbc7"),
                (2.0, full_removal_ratio, removal_seed_ratios, COLORS["a1_no_mcv"]),
            ):
                spread = float(np.std(observed, ddof=1)) if len(observed) > 1 else 0.0
                mcv_ax.vlines(
                    x_value,
                    1.0,
                    ratio,
                    color=COLORS["a1_no_mcv"],
                    linewidth=1.65,
                    alpha=0.82,
                    zorder=3,
                )
                mcv_ax.scatter(
                    [x_value],
                    [ratio],
                    s=18,
                    marker="o",
                    facecolor=facecolor,
                    edgecolor=COLORS["a1_no_mcv"],
                    linewidth=0.72,
                    zorder=5,
                )
                mcv_ax.errorbar(
                    [x_value],
                    [ratio],
                    yerr=[spread],
                    fmt="none",
                    ecolor="#292929",
                    elinewidth=0.72,
                    capsize=2.0,
                    capthick=0.72,
                    linestyle="none",
                    zorder=8,
                )
                mcv_ax.annotate(
                    f"{ratio:.2f}×",
                    xy=(x_value, ratio + spread),
                    xytext=(0, 4),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=5.7,
                    color="#252525",
                )
            style_ratio_axis(
                mcv_ax,
                bottom=bottom_row,
                show_y_labels=True,
                x_positions=(0.0, 1.0, 2.0),
                x_labels=("w/o infer", "w/o train", "w/o all"),
                x_limits=mcv_x_limits,
            )

            # Column 3: remove hard terminal anchoring.
            anchor_ax = fig.add_subplot(grid[row_index, 2])
            axes[row_index, 2] = anchor_ax
            anchor_x_limits = (-0.62, 0.62)
            ratio_background(anchor_ax, anchor_x_limits)
            anchor_ratio = (
                float(lookup[("a2_no_terminal_anchor", pde)]["rmse_mean"])
                / full_mean
            )
            anchor_by_seed = {
                int(row["seed"]): float(row["value_rmse"])
                for row in local_rows(pde, "a2_no_terminal_anchor")
            }
            anchor_observed = [
                anchor_by_seed[seed] / full_by_seed[seed]
                for seed in sorted(full_by_seed)
            ]
            anchor_spread = float(np.std(anchor_observed, ddof=1))
            anchor_ax.vlines(
                0.0,
                1.0,
                anchor_ratio,
                color=COLORS["a2_no_terminal_anchor"],
                linewidth=1.75,
                alpha=0.84,
                zorder=3,
            )
            anchor_ax.errorbar(
                [0.0],
                [anchor_ratio],
                yerr=[anchor_spread],
                fmt="none",
                ecolor="#292929",
                elinewidth=0.72,
                capsize=2.0,
                capthick=0.72,
                zorder=6,
            )
            anchor_ax.scatter(
                [0.0],
                [anchor_ratio],
                s=20,
                marker=markers["a2_no_terminal_anchor"],
                facecolor=COLORS["a2_no_terminal_anchor"],
                edgecolor="white",
                linewidth=0.55,
                zorder=5,
            )
            anchor_ax.annotate(
                f"{anchor_ratio:.2f}×",
                xy=(0.0, anchor_ratio + anchor_spread),
                xytext=(0, 4),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=5.7,
                color="#252525",
            )
            style_ratio_axis(
                anchor_ax,
                bottom=bottom_row,
                show_y_labels=False,
                x_positions=(0.0,),
                x_labels=("",),
                x_limits=anchor_x_limits,
            )
            if bottom_row:
                anchor_ax.set_xlabel("w/o anchor", fontsize=5.9, labelpad=2.0)

            # Column 4: formal paired Actor improvement curve, A0 through A4.
            update_ax = fig.add_subplot(grid[row_index, 3])
            axes[row_index, 3] = update_ax
            update_x_limits = (-0.45, 4.45)
            ratio_background(update_ax, update_x_limits)
            pde_iterations = [
                row for row in iteration_rows if row["pde"] == pde
            ]
            if len(pde_iterations) != 15:
                raise RuntimeError(
                    f"formal iteration curve requires 15 rows for {pde}; "
                    f"found {len(pde_iterations)}"
                )
            selected_iteration_by_seed = {
                int(row["seed"]): float(row["value_rmse"])
                for row in pde_iterations
                if int(row["iteration"])
                == int(row["selected_iteration"])
            }
            if len(selected_iteration_by_seed) != 3:
                raise RuntimeError(
                    f"formal iteration curve requires three paired selected "
                    f"Actors for {pde}; found {len(selected_iteration_by_seed)}"
                )
            update_ratios: list[float] = []
            update_spreads: list[float] = []
            for iteration in range(5):
                by_seed = {
                    int(row["seed"]): float(row["value_rmse"])
                    for row in pde_iterations
                    if int(row["iteration"]) == iteration
                }
                observed = [
                    by_seed[seed] / selected_iteration_by_seed[seed]
                    for seed in sorted(selected_iteration_by_seed)
                ]
                update_ratios.append(float(np.mean(observed)))
                update_spreads.append(float(np.std(observed, ddof=1)))
            update_x = np.arange(5, dtype=float)
            update_curve_color = "#C44755"
            update_ax.plot(
                update_x,
                update_ratios,
                color=update_curve_color,
                linewidth=1.15,
                alpha=0.78,
                zorder=3,
            )
            update_ax.errorbar(
                update_x,
                update_ratios,
                yerr=update_spreads,
                fmt="none",
                ecolor="#292929",
                elinewidth=0.68,
                capsize=1.7,
                capthick=0.68,
                zorder=6,
            )
            update_ax.scatter(
                update_x,
                update_ratios,
                s=14,
                marker="o",
                facecolor=update_curve_color,
                edgecolor="white",
                linewidth=0.45,
                zorder=5,
            )
            for point_index, (x_value, ratio, spread) in enumerate(zip(
                update_x,
                update_ratios,
                update_spreads,
                strict=True,
            )):
                # Keep labels close to their points.  Only the A1 labels in the
                # lower three panels need a little extra clearance.
                offset = (2, 9) if row_index > 0 and point_index == 1 else (0, 3)
                update_ax.annotate(
                    f"{ratio:.2f}×",
                    xy=(x_value, ratio + spread),
                    xytext=offset,
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=4.7,
                    color="#252525",
                )
            style_ratio_axis(
                update_ax,
                bottom=bottom_row,
                show_y_labels=False,
                x_positions=(0.0, 1.0, 2.0, 3.0, 4.0),
                x_labels=("$A_0$", "$A_1$", "$A_2$", "$A_3$", "$A_4$"),
                x_limits=update_x_limits,
            )

        axes[0, 0].text(
            0.025,
            0.955,
            "RMSE",
            transform=axes[0, 0].transAxes,
            ha="left",
            va="top",
            fontsize=5.5,
            color="#555555",
        )
        axes[0, 1].text(
            0.025,
            0.955,
            r"$\mathrm{RMSE}/\mathrm{RMSE}_{\mathrm{Full}}$",
            transform=axes[0, 1].transAxes,
            ha="left",
            va="top",
            fontsize=5.5,
            color="#555555",
        )
        for suffix, dpi in (("pdf", None), ("png", 300)):
            fig.savefig(
                FIGURE_DIR / f"component_ablation.{suffix}",
                dpi=dpi,
                bbox_inches="tight",
                facecolor="white",
            )
        plt.close(fig)


def _iteration_figure(rows: list[dict[str, Any]]) -> None:
    fig, axes = plt.subplots(
        1,
        4,
        figsize=(7.15, 2.4),
        constrained_layout=True,
        sharex=True,
    )
    titles = {
        "hjb_log": "Log HJB",
        "joint_concave": "Concave",
        "joint_convex": "Convex",
        "nonconvex_coupled": "Oscillatory",
    }
    for ax, pde in zip(axes, PDES, strict=True):
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
        ax.set_title(titles[pde], fontsize=8.0)
        ax.set_yscale("log")
        ax.set_xticks(range(5))
        ax.set_xticklabels(("$A_0$", "$A_1$", "$A_2$", "$A_3$", "$A_4$"), fontsize=7)
        ax.set_xlabel("Actor iterate", fontsize=7.2)
        ax.tick_params(axis="y", labelsize=7)
        ax.grid(which="both", alpha=0.22, linewidth=0.6)
        if len(seeds) > 1:
            ax.legend(fontsize=5.8, frameon=False, handlelength=1.4)
    axes[0].set_ylabel("Value RMSE", fontsize=7.5)
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
    iteration_rows: list[dict[str, Any]],
    actual_gpu_hours: float,
) -> None:
    lookup = {(row["variant"], row["pde"]): row for row in summary}
    primary_rows = [row for row in rows if row["primary_replicate"]]
    mcv_inference_lookup = {row["pde"]: row for row in mcv_inference}

    def actor_metric_mean(pde: str, variant: str, metric: str) -> float:
        values = [
            float(row[metric])
            for row in primary_rows
            if row["pde"] == pde
            and row["variant"] == variant
            and row.get(metric) is not None
        ]
        if not values:
            return float("nan")
        return _mean(values)

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
        "The formal curve below evaluates every stored Full-TA Actor iterate with the same 960-state bank, 64 Euler steps, 1,024 paths, and paired Brownian seed. The iteration diagnostic was executed on "
        + ", ".join(sorted({str(row.get("device", "unknown")) for row in iteration_rows}))
        + "; all five iterates for a seed use the same device-specific random stream. Each learned update receives the mature 5,000-step optimization budget, so A1 is the controlled one-update comparison rather than a prematurely stopped checkpoint.",
        "",
        "| PDE | A0 | A1 | A2 | A3 | A4 | Selected Full vs A1 |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ])
    iteration_interpretation: list[str] = []
    for pde in PDES:
        cells: list[tuple[float, float]] = []
        for iteration in range(5):
            values = [
                float(row["value_rmse"])
                for row in iteration_rows
                if row["pde"] == pde and int(row["iteration"]) == iteration
            ]
            if len(values) != 3:
                raise RuntimeError(
                    f"formal report requires three seeds for {pde} A{iteration}"
                )
            cells.append((_mean(values), _sample_std(values)))
        selected_values = [
            float(row["value_rmse"])
            for row in iteration_rows
            if row["pde"] == pde
            and int(row["iteration"])
            == int(row["selected_iteration"])
        ]
        if len(selected_values) != 3:
            raise RuntimeError(
                f"formal report requires three selected Actors for {pde}"
            )
        selected_mean = _mean(selected_values)
        improvement = (cells[1][0] - selected_mean) / cells[1][0]
        lines.append(
            f"| {PDE_NAMES[pde]} | "
            + " | ".join(f"{mean:.3e} ± {std:.1e}" for mean, std in cells)
            + f" | {improvement:+.1%} |"
        )
        selected_iterations = sorted({
            int(row["selected_iteration"])
            for row in iteration_rows
            if row["pde"] == pde
        })
        iteration_interpretation.append(
            f"- {PDE_NAMES[pde]}: the reference-free selected multi-update Actor "
            f"reduces mean RMSE by {improvement:.1%} relative to A1; selected "
            f"iteration(s) across seeds: {', '.join(f'A{value}' for value in selected_iterations)}."
        )
    lines.extend([
        "",
        *iteration_interpretation,
        "",
        "These curves resolve whether the existing mature protocol benefits from later tangent updates. They do not establish that no alternative, longer single-update optimization schedule could attain the same accuracy; that stronger claim would require a separate compute-matched continuation control.",
    ])
    lines.extend([
        "",
        "## Crossed MCV diagnostic",
        "",
        "The existing `w/o MCV` Actors were trained and selected entirely with plain labels. Each frozen Actor was also evaluated with FK–MCV on exactly the same Brownian paths already stored by the formal evaluator. This yields the missing crossed cell—`w/o MCV training` with MCV inference—without retraining or test-based reselection. The comparison therefore separates MCV in Actor learning from MCV in final inference.",
        "",
        "| PDE | Full (train + infer MCV) | w/o infer | w/o train | w/o all MCV | Infer ratio | Train ratio |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ])
    for row in mcv_inference:
        lines.append(
            f"| {row['pde_name']} | {row['full_mcv_rmse_mean']:.3e} ± {row['full_mcv_rmse_std']:.3e} | "
            f"{row['full_without_mcv_inference_rmse_mean']:.3e} ± {row['full_without_mcv_inference_rmse_std']:.3e} | "
            f"{row['without_mcv_training_rmse_mean']:.3e} ± {row['without_mcv_training_rmse_std']:.3e} | "
            f"{row['without_all_mcv_rmse_mean']:.3e} ± {row['without_all_mcv_rmse_std']:.3e} | "
            f"{row['plain_over_mcv_rmse_ratio_mean']:.2f}× | "
            f"{row['without_training_over_full_ratio_mean']:.2f}× |"
        )
    lines.extend([
        "",
        "Seed- and replicate-level paired values are in `MCV_INFERENCE_ONLY_SEED_METRICS.csv`; the separate summary is in `MCV_INFERENCE_ONLY_SUMMARY.csv`. The filenames are retained for artifact compatibility, but both files now contain the full 2×2 training/inference crossing.",
        "",
        "### MCV training-label effect",
        "",
        "The table below supplements the direct crossed value comparison with Actor diagnostics. Gradient RMSE is post-selection only and never entered fitting or checkpoint selection. `MCV-evaluated ratio` compares the no-MCV-trained Actor against Full TA while holding the final MCV estimator fixed.",
        "",
        "| PDE | Full Actor gradient RMSE | w/o-MCV Actor gradient RMSE | MCV-evaluated ratio |",
        "|---|---:|---:|---:|",
    ])
    for pde in PDES:
        full_gradient = actor_metric_mean(pde, "a0_full", "actor_gradient_vector_rmse")
        no_mcv_gradient = actor_metric_mean(pde, "a1_no_mcv", "actor_gradient_vector_rmse")
        lines.append(
            f"| {PDE_NAMES[pde]} | {full_gradient:.3e} | {no_mcv_gradient:.3e} | "
            f"{mcv_inference_lookup[pde]['without_training_over_full_ratio_mean']:.2f}× |"
        )
    lines.extend([
        "",
        "MCV-trained Actors have lower gradient RMSE on every PDE (about 1% for Log HJB, 55% for both joint cubic--quadratic PDEs, and 28% for the oscillatory PDE). With the final estimator held to FK–MCV, removing training MCV raises mean value RMSE by only about 2%–5% across the four PDEs. This is directionally consistent but small relative to three-seed variability in some columns; it must not be described as a universally resolved gain. In contrast, removing inference MCV raises RMSE by about 2.2×–7.0× and is the dominant measured MCV contribution under this budget.",
        "",
        "## Costs, failures, and numerical resolution",
        "",
        f"Actual Slurm allocated GPU-hours across the recorded jobs: {actual_gpu_hours:.3f}. The preregistered hard ceiling was 36. No variant-specific stabilization or post-test tuning was permitted.",
        "",
        f"Numerical failures recorded in the 16 summary cells: {sum(int(row['numerical_failures']) for row in summary)}. Minimum exponential-weight ESS fraction: {min(float(row['weight_ess_fraction_min']) for row in summary):.4f}.",
        "",
        "| PDE | Reference uncertainty (value RMSE) | Terminal-anchoring mean delta |",
        "|---|---:|---:|",
    ])
    for pde in PDES:
        reference_row = next(row for row in primary_rows if row["pde"] == pde)
        reference_uncertainty = float(reference_row["reference"]["reference_uncertainty_value_rmse"])
        terminal_delta = float(lookup[("a2_no_terminal_anchor", pde)]["rmse_mean"]) - float(
            lookup[("a0_full", pde)]["rmse_mean"]
        )
        lines.append(
            f"| {PDE_NAMES[pde]} | {reference_uncertainty:.3e} | {terminal_delta:+.3e} |"
        )
    lines.extend([
        "",
        "The anchoring ablation is consistently worse in mean RMSE, but the concave and convex deltas are below their independent reference-resolution estimates and are therefore unresolved. The oscillatory delta exceeds its reference floor but is comparable to three-seed variability. Log HJB has negligible quadrature-reference uncertainty and shows the clearest anchoring effect, although it remains modest relative to the MCV and tangent-update effects. The MCV and fixed-A0 degradations are much larger than all listed reference floors.",
        "",
        "The exact reference metadata for every seed is retained in `ABLATION_SEED_METRICS.csv`; the preflight reference audits are in `PREFLIGHT_AUDIT.json`. Small apparent rank changes comparable to reference resolution or independent inference-replicate fluctuation are treated as unresolved.",
        "",
        "## Limitations",
        "",
        "Three training seeds are sufficient for an ablation signal but not a universal variance claim. The formal A0–A4 curves use the same frozen 960-state banks, 64 Euler steps, 1,024 paths, and paired evaluation seeds as the primary study for all four PDEs. They compare stored per-update Actors under a common inference protocol; they do not prove that an arbitrarily longer-trained single update could never match later updates. These ablations neither establish universally lower variance than FH-MLP nor identify fixed A0 with FH-MLP.",
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
    docs.parent.mkdir(parents=True, exist_ok=True)
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
    iteration_rows = _iteration_rows()
    iteration_keys = {
        (row["pde"], int(row["seed"]), int(row["iteration"]))
        for row in iteration_rows
    }
    expected_iteration_keys = {
        (pde, seed, iteration)
        for pde in PDES
        for seed in (0, 1, 2)
        for iteration in range(5)
    }
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
        "iteration_records": len(iteration_rows),
        "expected_iteration_records": 60,
        "iteration_key_coverage_complete": (
            iteration_keys == expected_iteration_keys
        ),
        "iteration_protocol_uniform": all(
            int(row["states"]) == 960
            and int(row["steps"]) == 64
            and int(row["paths"]) == 1024
            for row in iteration_rows
        ),
        "iteration_values_finite": all(
            np.isfinite(float(row["value_rmse"]))
            for row in iteration_rows
        ),
        "iteration_devices": sorted({
            str(row.get("device", "unknown")) for row in iteration_rows
        }),
    }
    audit["passed"] = bool(
        all(item["exists"] and int(item["bytes"] or 0) > 0 for item in audit["required_files"])
        and audit["evaluation_records"] == expected_evaluations
        and audit["all_test_banks_960"]
        and audit["all_evaluations_finite"]
        and audit["ground_truth_post_selection_only"]
        and audit["iteration_records"] == 60
        and audit["iteration_key_coverage_complete"]
        and audit["iteration_protocol_uniform"]
        and audit["iteration_values_finite"]
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
    _component_vertical_matrix_figure(summary, rows, mcv_inference, iterations)
    _iteration_figure(iterations)
    _report(summary, paired, rows, mcv_inference, iterations, actual)
    _artifact_audit(_expected_count())
    print(json.dumps({
        "evaluations": len(rows),
        "summary_cells": len(summary),
        "actual_allocated_gpu_hours": actual,
        "artifact_audit": "passed",
    }, sort_keys=True))


if __name__ == "__main__":
    main()
