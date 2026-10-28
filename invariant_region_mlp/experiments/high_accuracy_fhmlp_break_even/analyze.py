"""Measured break-even analysis and publication figure generation."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from .protocol import (
    FIGURE_DIR,
    RESULTS_DIR,
    REVISION,
    SOURCE_ROOT,
    TARGETS,
    require_frozen_selections,
    require_frozen_timing_grid,
    write_json_atomic,
)


COLORS = {
    "ta_direct": "#2166ac",
    "ta_fk_mcv": "#4f9ed8",
    "fh_mlp_standard": "#d95f02",
}
LABELS = {
    "ta_direct": "Tangent-Actor (direct)",
    "ta_fk_mcv": "Tangent-Actor (FK-MCV)",
    "fh_mlp_standard": "Full-History MLP",
}


def _read(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def _records(method: str, candidate: str) -> list[dict[str, Any]]:
    if method in {"ta_direct", "ta_fk_mcv"}:
        paths = [
            RESULTS_DIR / "raw" / "timing" / "tangent_actor"
            / f"ta__timing__{candidate}__r{replicate}.json"
            for replicate in range(3)
        ]
    else:
        n_text, m_text = candidate.split("_")
        n, M = int(n_text[1:]), int(m_text[1:])
        paths = [
            RESULTS_DIR / "raw" / "timing" / "fh_mlp"
            / f"fh__timing__n{n}_M{M}__r{replicate}.json"
            for replicate in range(3)
        ]
    if not all(path.is_file() for path in paths):
        raise RuntimeError(f"missing timing records for {method}/{candidate}: {paths}")
    return [_read(path) for path in paths]


def _mean_curve(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    mappings = [
        {int(row["queries"]): float(row["total_seconds"]) for row in record["prefix_rows"]}
        for record in records
    ]
    common = sorted(set.intersection(*(set(mapping) for mapping in mappings)))
    rows = []
    for queries in common:
        values = np.asarray([mapping[queries] for mapping in mappings])
        rows.append({
            "queries": queries,
            "mean_seconds": float(values.mean()),
            "std_seconds": float(values.std(ddof=1)),
            "min_seconds": float(values.min()),
            "max_seconds": float(values.max()),
            "replicate_seconds": values.tolist(),
        })
    return rows


def _crossing(ta: list[dict[str, Any]], fh: list[dict[str, Any]]) -> dict[str, Any]:
    left = {row["queries"]: row for row in ta}
    right = {row["queries"]: row for row in fh}
    common = sorted(set(left) & set(right))
    if not common:
        return {"status": "NO_COMMON_EXECUTED_PREFIX"}
    replicate_crossings = []
    for replicate in range(3):
        previous, observed = None, None
        for query in common:
            if right[query]["replicate_seconds"][replicate] >= left[query]["replicate_seconds"][replicate]:
                observed = query
                break
            previous = query
        replicate_crossings.append({
            "replicate": replicate,
            "previous_queries": previous if observed is not None else common[-1],
            "crossing_queries": observed,
            "status": "OBSERVED_BREAK_EVEN" if observed is not None else "NOT_REACHED_WITHIN_TESTED_BUDGET",
        })
    previous = None
    for index, query in enumerate(common):
        if right[query]["mean_seconds"] >= left[query]["mean_seconds"]:
            beyond = common[min(index + 1, len(common) - 1)]
            savings = 1.0 - left[beyond]["mean_seconds"] / right[beyond]["mean_seconds"]
            return {
                "status": "OBSERVED_BREAK_EVEN",
                "previous_queries": 0 if previous is None else previous,
                "crossing_queries": query,
                "ta_mean_seconds_at_crossing": left[query]["mean_seconds"],
                "ta_std_seconds_at_crossing": left[query]["std_seconds"],
                "fh_mean_seconds_at_crossing": right[query]["mean_seconds"],
                "fh_std_seconds_at_crossing": right[query]["std_seconds"],
                "representative_beyond_queries": beyond,
                "ta_mean_seconds_beyond": left[beyond]["mean_seconds"],
                "fh_mean_seconds_beyond": right[beyond]["mean_seconds"],
                "ta_relative_computation_saved_beyond": savings,
                "per_replicate_crossings": replicate_crossings,
            }
        previous = query
    query = common[-1]
    return {
        "status": "NOT_REACHED_WITHIN_TESTED_BUDGET",
        "largest_common_queries": query,
        "ta_mean_seconds": left[query]["mean_seconds"],
        "fh_mean_seconds": right[query]["mean_seconds"],
        "fh_over_ta_time": right[query]["mean_seconds"] / left[query]["mean_seconds"],
        "per_replicate_crossings": replicate_crossings,
    }


def _draw_curve(axis, curve, method, label=None):
    x = np.asarray([row["queries"] for row in curve])
    mean = np.asarray([row["mean_seconds"] for row in curve]) / 60.0
    std = np.asarray([row["std_seconds"] for row in curve]) / 60.0
    axis.plot(
        x, mean, marker="o", markersize=3.6, linewidth=2.0,
        color=COLORS[method], label=label or LABELS[method], zorder=3,
    )
    axis.fill_between(
        x, np.maximum(mean - std, 0.0), mean + std,
        color=COLORS[method], alpha=0.14, linewidth=0, zorder=1,
    )


def _plot_target(
    target: str,
    threshold: float,
    ta_curve: list[dict[str, Any]],
    fh_curve: list[dict[str, Any]],
    crossing: dict[str, Any],
) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig, axis = plt.subplots(figsize=(7.15, 4.45), constrained_layout=True)
    _draw_curve(axis, ta_curve, "ta_direct")
    _draw_curve(axis, fh_curve, "fh_mlp_standard")
    common_max = min(ta_curve[-1]["queries"], fh_curve[-1]["queries"])
    if crossing["status"] == "OBSERVED_BREAK_EVEN":
        previous = crossing["previous_queries"]
        observed = crossing["crossing_queries"]
        plot_previous = max(1, previous)
        axis.axvspan(plot_previous, observed, color="#777777", alpha=0.10, zorder=0)
        axis.axvline(observed, color="#555555", linewidth=1.1, linestyle=":", zorder=2)
        bracket_mid = float(np.sqrt(plot_previous * observed))
        axis.text(
            bracket_mid, 0.955, "Measured break-even bracket",
            transform=axis.get_xaxis_transform(), ha="center", va="top",
            fontsize=8.2, color="#444444",
        )
        axis.annotate(
            f"First measured crossing\n{observed:,} queries",
            xy=(observed, crossing["fh_mean_seconds_at_crossing"] / 60.0),
            xytext=(-12, -35), textcoords="offset points", fontsize=8.2,
            ha="right", va="top",
            arrowprops={"arrowstyle": "-", "color": "#555555", "lw": 0.8},
        )
        beyond = crossing["representative_beyond_queries"]
        saved = 100.0 * crossing["ta_relative_computation_saved_beyond"]
        axis.text(
            0.975, 0.055,
            f"At {beyond:,} queries: TA uses {saved:.1f}% less computation",
            transform=axis.transAxes, ha="right", va="bottom", fontsize=8.2,
            color="#333333",
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.90},
        )
    axis.set_xscale("log", base=2)
    powers = range(0, 21, 4)
    axis.set_xticks([2**power for power in powers], [rf"$2^{{{power}}}$" for power in powers])
    axis.set_xlim(0.85, max(2, common_max) * 1.10)
    axis.set_ylim(bottom=0.0)
    axis.set_xlabel("Number of query points")
    axis.set_ylabel("Measured cumulative H100 computation (minutes)")
    axis.set_title(
        rf"100D Joint Nonseparable PDE  |  common RMSE $\leq {threshold:.0e}$",
        pad=9,
    )
    axis.grid(True, alpha=0.20, linewidth=0.6)
    axis.legend(loc="center left", frameon=False)
    stem = f"high_accuracy_break_even_{target}"
    fig.savefig(FIGURE_DIR / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIGURE_DIR / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(FIGURE_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def _plot_mcv(curves: dict[tuple[str, str], list[dict[str, Any]]], final: dict[str, Any]) -> None:
    primary = final["targets"]["primary"]
    mcv = primary["ta_fk_mcv"]
    fh = primary["fh_mlp_standard"]
    direct = primary["ta_direct"]
    if not (
        mcv["status"] == "FINAL_ACCURACY_PASSED"
        and fh["status"] == "FINAL_ACCURACY_PASSED"
        and direct["status"] == "FINAL_ACCURACY_PASSED"
    ):
        return
    fig, axis = plt.subplots(figsize=(7.0, 4.2), constrained_layout=True)
    _draw_curve(axis, curves[("ta_direct", "direct")], "ta_direct")
    _draw_curve(axis, curves[("ta_fk_mcv", mcv["selected_candidate"])], "ta_fk_mcv")
    _draw_curve(axis, curves[("fh_mlp_standard", fh["selected_candidate"])], "fh_mlp_standard")
    axis.set_xscale("log", base=2)
    axis.set_xlabel("Number of query points")
    axis.set_ylabel("Measured cumulative H100 computation (minutes)")
    axis.set_title("Direct and FK-MCV Tangent-Actor Inference", pad=9)
    axis.grid(True, alpha=0.20, linewidth=0.6)
    axis.legend(frameon=False)
    for suffix in ("pdf", "svg", "png"):
        kwargs = {"dpi": 300} if suffix == "png" else {}
        fig.savefig(FIGURE_DIR / f"ta_inference_modes_timing.{suffix}", bbox_inches="tight", **kwargs)
    plt.close(fig)


def analyze() -> dict[str, Any]:
    frozen = require_frozen_selections()
    require_frozen_timing_grid()
    final = _read(RESULTS_DIR / "final_accuracy_summary.json")
    coverage = _read(RESULTS_DIR / "timing_coverage.json")
    curves: dict[tuple[str, str], list[dict[str, Any]]] = {}
    replicate_rows, mean_rows = [], []
    needed: set[tuple[str, str]] = set()
    for target in TARGETS:
        for method in ("ta_direct", "ta_fk_mcv", "fh_mlp_standard"):
            item = final["targets"][target][method]
            if coverage["targets"][target][method] == "RUN" and item.get("selected_candidate"):
                needed.add((method, item["selected_candidate"]))
    for key in sorted(needed):
        records = _records(*key)
        curves[key] = _mean_curve(records)
        for replicate, record in enumerate(records):
            for row in record["prefix_rows"]:
                replicate_rows.append({
                    "method": key[0], "candidate": key[1], "replicate": replicate,
                    "queries": row["queries"], "total_seconds": row["total_seconds"],
                    "inference_seconds": row.get("inference_seconds"),
                })
        for row in curves[key]:
            mean_rows.append({
                "method": key[0], "candidate": key[1], "queries": row["queries"],
                "mean_seconds": row["mean_seconds"], "std_seconds": row["std_seconds"],
                "min_seconds": row["min_seconds"], "max_seconds": row["max_seconds"],
            })
    _write_csv(RESULTS_DIR / "timing_replicates.csv", replicate_rows)
    _write_csv(RESULTS_DIR / "timing_rows.csv", mean_rows)

    break_even, flat = {}, []
    for target, threshold in TARGETS.items():
        local = final["targets"][target]
        if not (
            coverage["targets"][target]["ta_direct"] == "RUN"
            and coverage["targets"][target]["fh_mlp_standard"] == "RUN"
        ):
            result = {"status": "PAIR_NOT_ACCURACY_QUALIFIED"}
        else:
            result = _crossing(
                curves[("ta_direct", local["ta_direct"]["selected_candidate"])],
                curves[("fh_mlp_standard", local["fh_mlp_standard"]["selected_candidate"])],
            )
            _plot_target(
                target,
                threshold,
                curves[("ta_direct", local["ta_direct"]["selected_candidate"])],
                curves[("fh_mlp_standard", local["fh_mlp_standard"]["selected_candidate"])],
                result,
            )
        break_even[target] = result
        flat.append({"target": target, "threshold": threshold, **result})
    _write_csv(RESULTS_DIR / "break_even_summary.csv", flat)
    write_json_atomic(RESULTS_DIR / "break_even_summary.json", break_even)
    _plot_mcv(curves, final)

    primary = break_even["primary"]
    if primary["status"] == "OBSERVED_BREAK_EVEN":
        caption = (
            "Measured high-accuracy amortization on the 100-dimensional Joint "
            "Nonseparable PDE. Both methods satisfy the independently frozen "
            f"RMSE threshold of $10^{{-2}}$. Each marker is an actually executed "
            "query prefix on H100 hardware; ribbons show one standard deviation "
            "over three timing replicates. The shaded band is the measured "
            f"break-even bracket $({primary['previous_queries']:,}, "
            f"{primary['crossing_queries']:,}]$; no interpolation is used."
        )
    else:
        caption = (
            "Measured high-accuracy amortization on the 100-dimensional Joint "
            "Nonseparable PDE. Both methods satisfy the independently frozen "
            "RMSE threshold of $10^{-2}$. Each marker is an actually executed "
            "query prefix on H100 hardware; no crossing is extrapolated."
        )
    (RESULTS_DIR / "ICML_FIGURE_CAPTION.tex").write_text(caption + "\n", encoding="utf-8")

    paragraph = [
        "## Ready-to-insert results paragraph",
        "",
        "At the prespecified high-accuracy thresholds, we calibrated each method on a fresh bank and verified the frozen choices on a disjoint final bank.",
    ]
    for target in TARGETS:
        item = break_even[target]
        if item["status"] == "OBSERVED_BREAK_EVEN":
            paragraph.append(
                f"For {target} (RMSE <= {TARGETS[target]:.3g}), the observed break-even lies in ({item['previous_queries']:,}, {item['crossing_queries']:,}] actually executed queries."
            )
        else:
            paragraph.append(f"For {target}, the break-even status was `{item['status']}`.")
    paragraph.append(
        "The main curve uses the learned direct value head; FK-MCV is reported separately because the current manuscript describes it as the final evaluation estimator."
    )
    (RESULTS_DIR / "RESULTS_PARAGRAPH.md").write_text("\n".join(paragraph) + "\n", encoding="utf-8")

    commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=SOURCE_ROOT, text=True
    ).strip()
    payload = {
        "revision": REVISION,
        "status": "COMPLETE",
        "analysis_code_commit": commit,
        "frozen_selections_canonical_sha256": frozen["canonical_sha256"],
        "targets": final["targets"],
        "break_even": break_even,
        "timing_coverage": coverage,
        "curve_keys": [list(key) for key in sorted(curves)],
        "analyzed_utc": datetime.now(timezone.utc).isoformat(),
    }
    write_json_atomic(RESULTS_DIR / "analysis_summary.json", payload)
    return payload


if __name__ == "__main__":
    print(json.dumps(analyze(), indent=2, sort_keys=True))
