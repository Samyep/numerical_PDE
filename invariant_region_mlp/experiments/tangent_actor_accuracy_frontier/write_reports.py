"""Render the final scientific reports from frozen aggregate tables."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np

from .common import RESULTS_DIR


def _rows(name: str) -> list[dict[str, str]]:
    with (RESULTS_DIR / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _f(value) -> float | None:
    if value in (None, "", "None", "nan"):
        return None
    return float(value)


def _mean_sd(values) -> tuple[float, float]:
    array = np.asarray(list(values), dtype=np.float64)
    return float(array.mean()), float(array.std(ddof=1)) if len(array) > 1 else 0.0


def _write(name: str, lines: list[str]) -> None:
    (RESULTS_DIR / name).write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def _high_precision(euler: list[dict[str, str]], reference: dict) -> None:
    primary = [row for row in euler if row["family"] == "primary"]
    lines = [
        "# High-Precision Euler Audit",
        "",
        "This is a follow-up precision experiment on the three already selected frozen `a4` Actors. No Actor was retrained or reselected. Fine Brownian increments were regenerated on an independent stream and aggregated to construct common random numbers across 64, 128, 256, and 512 Euler steps.",
        "",
        f"The numerical reference was refined to `n={reference['refined_n']}`, `dt={reference['refined_dt']:.7g}`. The complete previous-to-current discrepancy, `{reference['reference_uncertainty_value_rmse']:.3e}`, is propagated as reference uncertainty.",
        "",
        "| Euler steps | seeds | Value RMSE mean ± sd | residual non-MC proxy | MC contribution | mean bias | MCV VR | ESS | observed weak order |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for steps in (64, 128, 256, 512):
        local = [row for row in primary if int(row["steps"]) == steps]
        if not local:
            continue
        rm, rs = _mean_sd(_f(row["value_rmse"]) for row in local)
        residual = np.mean([_f(row["residual_error_proxy"]) for row in local])
        mc = np.mean([_f(row["mc_rmse_contribution"]) for row in local])
        bias = np.mean([_f(row["bias"]) for row in local])
        vr = np.mean([_f(row["mcv_variance_reduction"]) for row in local])
        ess = np.mean([_f(row["ess"]) for row in local])
        orders = [_f(row.get("observed_weak_order")) for row in local]
        orders = [value for value in orders if value is not None]
        order_text = "—" if not orders else f"{np.mean(orders):.2f}"
        lines.append(f"| {steps} | {len(local)} | {rm:.6g} ± {rs:.2g} | {residual:.6g} | {mc:.6g} | {bias:+.3e} | {vr:.1f}× | {ess:.4f} | {order_text} |")
    lines.extend([
        "",
        "The reported raw RMSE includes finite-path Monte Carlo error. The residual proxy removes the estimated antithetic-pair sampling contribution and the reference discrepancy in quadrature; it is a diagnostic rather than an alternative headline metric. Confidence intervals use a state bootstrap with Gaussian perturbations calibrated to per-state MC variance and the reference proxy.",
        "",
        "The observed order is reported only when two successive common-random-number refinement differences remain positive after noise correction. It is not forced to first order.",
    ])
    _write("HIGH_PRECISION_EULER_AUDIT.md", lines)


def _full_history_report(full: list[dict[str, str]], summary: dict) -> None:
    lines = [
        "# Strong Full-History Comparison",
        "",
        "The baseline is the raw recursive Full-History MLP estimator. It has no Tangent-Actor, MCV, terminal-aware learned correction, neural surrogate, value projection, or oracle active subspace. Large terminal Monte Carlo tensors are accumulated in blocks to fit device memory; this changes storage and summation order only, not the estimator or sample count.",
        "",
        "| family | n | M | role | finite seeds | Value RMSE mean ± sd | wall seconds mean | peak GPU GB |",
        "|---|---:|---:|---|---:|---:|---:|---:|",
    ]
    groups = {}
    for row in full:
        if row["source"] != "accuracy_frontier":
            continue
        groups.setdefault((row["family"], int(row["n"]), int(row["M"]), row["role"]), []).append(row)
    for key, local in sorted(groups.items()):
        finite = [row for row in local if row["finite"] == "True"]
        if finite:
            rm, rs = _mean_sd(_f(row["value_rmse"]) for row in finite)
            wall = np.mean([_f(row["wall_seconds"]) for row in finite])
            peaks = [_f(row["peak_gpu_bytes"]) for row in finite]
            peaks = [value for value in peaks if value is not None]
            peak = max(peaks) / 2**30 if peaks else math.nan
            metric = f"{rm:.6g} ± {rs:.2g}"
        else:
            wall, peak, metric = math.nan, math.nan, "nonfinite"
        lines.append(f"| {key[0]} | {key[1]} | {key[2]} | {key[3]} | {len(finite)}/{len(local)} | {metric} | {wall:.1f} | {peak:.2f} |")
    strongest = summary.get("strongest_full_history_endpoint")
    lines.extend(["", "## Fair conclusion", ""])
    if strongest:
        lines.append(f"The strongest completed primary endpoint was `n={strongest['n']}, M={strongest['M']}` with mean RMSE `{strongest['rmse_mean']:.6g}` across {strongest['seeds']} seeds. Configuration grids were frozen from resource feasibility and algorithmic budget, and failed or worsening rows remain in the table.")
    lines.append("Runtime is kept visible, while the primary comparison concerns maximum attainable numerical accuracy. Historical rows were recomputed against the refined reference whenever stored predictions were available.")
    _write("STRONG_FULL_HISTORY_COMPARISON.md", lines)


def _nonseparable_report(rows: list[dict[str, str]], precision: list[dict[str, str]]) -> None:
    lines = [
        "# Nonseparable High-Budget Report",
        "",
        "The unchanged undamped generic tangent update was trained on `f(u,z)=0.2 sin(2πu)-0.5||z+u e||²` with full ambient inputs. Three seeds use the same mature budget and reference-free selection rule as the primary experiment.",
        "",
        "| iteration | seeds | 64-step Value RMSE mean ± sd | Gradient-vector RMSE | Value-Actor RMSE | fixed-point score | MCV VR | ESS |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for iteration in ("a0", "a1", "a2", "a3", "a4"):
        local = [row for row in rows if row["iteration"] == iteration]
        if not local:
            continue
        rm, rs = _mean_sd(_f(row["value_rmse_64"]) for row in local)
        lines.append(
            f"| {iteration} | {len(local)} | {rm:.6g} ± {rs:.2g} | "
            f"{np.mean([_f(r['gradient_vector_rmse']) for r in local]):.6g} | "
            f"{np.mean([_f(r['value_actor_rmse']) for r in local]):.6g} | "
            f"{np.mean([_f(r['fixed_point_score']) for r in local]):.5g} | "
            f"{np.mean([_f(r['mcv_variance_reduction']) for r in local]):.2f}× | "
            f"{np.mean([_f(r['ess']) for r in local]):.4f} |"
        )
    selected = [row for row in rows if row["selected"] == "True"]
    if selected:
        initial = {row["seed"]: row for row in rows if row["iteration"] == "a0"}
        reductions = [
            (_f(initial[row["seed"]]["value_rmse_64"]) - _f(row["value_rmse_64"]))
            / _f(initial[row["seed"]]["value_rmse_64"])
            for row in selected
        ]
        lines.extend(["", f"Reference-free selection produced {len(selected)} frozen seeds. Their mean observed 64-step improvement from `a0` is `{np.mean(reductions):.1%}`. This benchmark supplies local empirical evidence for genuinely mixed `(u,z)` curvature, not a universal nonconvex convergence theorem."])
    coupled_precision = [row for row in precision if row["family"] == "coupled"]
    if coupled_precision:
        lines.extend(["", "## Selected-Actor Euler audit", "", "| steps | RMSE mean ± sd | MC contribution |", "|---:|---:|---:|"])
        for steps in sorted({int(row["steps"]) for row in coupled_precision}):
            local = [row for row in coupled_precision if int(row["steps"]) == steps]
            mean, sd = _mean_sd(_f(row["value_rmse"]) for row in local)
            lines.append(f"| {steps} | {mean:.6g} ± {sd:.2g} | {np.mean([_f(r['mc_rmse_contribution']) for r in local]):.6g} |")
    _write("NONSEPARABLE_HIGH_BUDGET_REPORT.md", lines)


def _ambient_report(rows: list[dict[str, str]]) -> None:
    lines = [
        "# Ambient 400D Generalization Audit",
        "",
        "Known active directions are used only in this post-freeze decomposition. Every trained Actor receives the full ambient state and every learned low-rank matrix is randomly initialized and optimized end to end.",
        "",
        "| source | d | architecture | rank | budget | selected | a0 RMSE | selected RMSE | improvement | Gradient RMSE | active / inactive | effective rank | FP score |",
        "|---|---:|---|---:|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        rank = "—" if not row["rank"] else row["rank"]
        effective = "—" if not row["gradient_effective_rank"] else f"{_f(row['gradient_effective_rank']):.2f}"
        lines.append(
            f"| {row['source']} | {row['d']} | {row['architecture']} | {rank} | {row['budget']} | {row['selected_actor']} | "
            f"{_f(row['a0_value_rmse']):.6g} | {_f(row['selected_value_rmse']):.6g} | {_f(row['value_improvement_fraction']):.1%} | "
            f"{_f(row['gradient_vector_rmse']):.6g} | {_f(row['active_gradient_rmse']):.5g} / {_f(row['inactive_gradient_rms']):.5g} | "
            f"{effective} | {_f(row['fixed_point_score']):.5g} |"
        )
    lines.extend([
        "",
        "The diagnosis distinguishes supervision noise from statistical generalization by combining label-replication error, train/validation fit, active/inactive post-freeze errors, and the learned output spectrum. A low-rank row is an explicitly named ablation and does not replace the full-ambient primary method.",
    ])
    _write("AMBIENT_400D_GENERALIZATION_AUDIT.md", lines)


def _main_report(summary: dict, euler, full, nonsep, ambient) -> None:
    decision = summary["decision"]
    strongest = summary.get("strongest_full_history_endpoint")
    actor = summary.get("primary_actor_512_rmse_mean")
    lines = [
        "# Tangent-Actor High-Accuracy Frontier",
        "",
        f"## Decision: **{decision}**",
        "",
        "This follow-up tests the maximum numerical accuracy of frozen Tangent-Actor models, a stronger raw Full-History frontier, a three-seed nonseparable problem, and matched-budget 400D training. Ground Truth is absent from Actor training and checkpoint selection throughout.",
        "",
        "## Decisive findings",
        "",
    ]
    if actor is not None and not math.isnan(actor):
        lines.append(f"- The independent 512-step primary Tangent-Actor mean RMSE is `{actor:.6g}` across three frozen training seeds.")
    if strongest:
        relation = "lower" if actor is not None and actor < strongest["rmse_mean"] else "higher"
        lines.append(f"- The strongest completed Full-History endpoint is `n={strongest['n']}, M={strongest['M']}` at `{strongest['rmse_mean']:.6g}` mean RMSE; the Actor value is {relation}.")
    lines.extend([
        f"- A primary accuracy advantage that clears MC and reference uncertainty is `{summary['robust_primary_advantage']}`.",
        f"- The three-seed high-budget nonseparable support criterion is `{summary['nonseparable_support']}`.",
        f"- Matched-budget full-ambient training fixes the 400D degradation under the predeclared 10% criterion: `{summary['d400_fixed']}`.",
        "",
        "## Answers to the research questions",
        "",
        "1. **128-step replication.** See `HIGH_PRECISION_EULER_AUDIT.md`; all three frozen Actors use a new MC stream and CRN refinement.",
        "2. **256 and 512 steps.** The Euler table reports raw, MC, reference, and residual components and measured weak differences.",
        "3. **Full-History ceiling.** `STRONG_FULL_HISTORY_COMPARISON.md` retains every completed or failed depth/budget row.",
        f"4. **Strongest primary comparison.** Robust Tangent-Actor advantage: **{summary['robust_primary_advantage']}**.",
        f"5. **Nonseparable improvement.** Three-seed supporting result: **{summary['nonseparable_support']}**.",
        f"6. **Seed stability.** The primary 512-step seed standard deviation is `{summary.get('primary_actor_512_rmse_sd', math.nan):.3g}`; nonseparable trajectories are reported without selecting the best seed.",
        f"7. **400D matched budget.** Predeclared success criterion met: **{summary['d400_fixed']}**.",
        f"8. **Learned low rank.** `{summary['low_rank_runs']}` completed ablation rows are reported individually; no oracle subspace is used.",
        "9. **Remaining error.** The Euler and ambient audits separate time discretization, finite-path MC, reference accuracy, Actor error, and full-ambient generalization.",
        "10. **Defensible claim.** Tangent-Actor is an accurate learned linearization method for the tested smooth bounded-curvature, locally stable semilinear PDEs. Claims of universal nonconvex convergence, a global reference-free certificate, or general speed superiority are not supported.",
        "",
        "## Artifacts",
        "",
        "- `HIGH_PRECISION_EULER_AUDIT.md`",
        "- `STRONG_FULL_HISTORY_COMPARISON.md`",
        "- `NONSEPARABLE_HIGH_BUDGET_REPORT.md`",
        "- `AMBIENT_400D_GENERALIZATION_AUDIT.md`",
        "- `REFERENCE_FREE_ERROR_THEORY_NOTE.md`",
        "- `accuracy_frontier.csv`, `euler_convergence.csv`, `full_history_extended_grid.csv`",
        "- `nonseparable_three_seed_results.csv`, `ambient_scaling_comparison.csv`, `compute_accounting.csv`",
        "- `reproducibility_manifest.json`",
        "- `figures/primary_accuracy_frontier.{png,pdf}`",
        "- `figures/nonseparable_high_budget.{png,pdf}`",
        "- `figures/ambient_400d_generalization.{png,pdf}`",
    ])
    _write("ACCURACY_FRONTIER_REPORT.md", lines)


def run() -> None:
    summary = json.loads((RESULTS_DIR / "FINAL_SUMMARY.json").read_text(encoding="utf-8"))
    reference = json.loads((REFERENCE_DIR := RESULTS_DIR / "references" / "primary_refinement_audit.json").read_text(encoding="utf-8"))
    euler = _rows("euler_convergence.csv")
    full = _rows("full_history_extended_grid.csv")
    nonsep = _rows("nonseparable_three_seed_results.csv")
    ambient = _rows("ambient_scaling_comparison.csv")
    _high_precision(euler, reference)
    _full_history_report(full, summary)
    _nonseparable_report(nonsep, euler)
    _ambient_report(ambient)
    _main_report(summary, euler, full, nonsep, ambient)


if __name__ == "__main__":
    run()
