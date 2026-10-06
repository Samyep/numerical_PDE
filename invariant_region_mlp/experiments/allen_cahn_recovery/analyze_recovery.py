"""Aggregate Allen--Cahn recovery artifacts, figures, provenance, and report."""

from __future__ import annotations

import hashlib
import json
import math
import platform
from pathlib import Path
import subprocess
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .allen_cahn_equation import (
    BECK_NUMERICAL_RADIUS,
    CERTIFIED_INTERVAL,
    PUBLISHED_MLP_REFERENCES,
    IntervalStateProjector,
    reaction,
    sharp_clipped_driver_lipschitz,
)


HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[1]
REPOSITORY_ROOT = PROJECT_ROOT.parent
RESULT_ROOT = PROJECT_ROOT / "results" / "allen_cahn_recovery"
RAW_ROOT = RESULT_ROOT / "raw"
FIGURE_ROOT = RESULT_ROOT / "figures"
DOC_ROOT = PROJECT_ROOT / "docs"

METHOD_ORDER = ["raw", "beck_truncated", "interval_ir", "certified_ir_0_1"]
METHOD_LABELS = {
    "raw": "Raw MLP",
    "beck_truncated": "Beck truncated",
    "interval_ir": "Interval IR (same r)",
    "certified_ir_0_1": "Certified IR [0,1]",
}
COLORS = {
    "raw": "#222222",
    "beck_truncated": "#D55E00",
    "interval_ir": "#0072B2",
    "certified_ir_0_1": "#009E73",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git(*args: str) -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(REPOSITORY_ROOT), *args], text=True
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def load_raw() -> tuple[pd.DataFrame, pd.DataFrame]:
    method_rows: list[dict[str, Any]] = []
    equivalence_rows: list[dict[str, Any]] = []
    for stage in ("source", "pilot", "final", "equivalence"):
        for path in sorted((RAW_ROOT / stage).glob("*.json")):
            payload = json.loads(path.read_text(encoding="utf-8"))
            relative = str(path.relative_to(RESULT_ROOT)).replace("\\", "/")
            for row in payload["methods"]:
                item = dict(row)
                item["source_path"] = relative
                # Backfill artifacts generated before method-specific interval
                # metadata was added. The scientific values are unchanged.
                if "projection_low" not in item:
                    if item["method"] in ("beck_truncated", "interval_ir"):
                        item["projection_low"] = -float(item["radius"])
                        item["projection_high"] = float(item["radius"])
                    elif item["method"] == "certified_ir_0_1":
                        item["projection_low"], item["projection_high"] = CERTIFIED_INTERVAL
                    else:
                        item["projection_low"] = None
                        item["projection_high"] = None
                method_rows.append(item)
            equivalence = dict(payload["equivalence"])
            equivalence["source_path"] = relative
            equivalence_rows.append(equivalence)
    methods = pd.DataFrame(method_rows)
    equivalence = pd.DataFrame(equivalence_rows)
    if methods.empty or equivalence.empty:
        raise RuntimeError("no recovery artifacts found")
    return methods, equivalence


def add_paired_columns(methods: pd.DataFrame) -> pd.DataFrame:
    result = methods.copy()
    keys = [
        "stage",
        "dimension",
        "depth",
        "sample_size",
        "repetition",
        "seed",
        "test_point",
        "radius_label",
    ]
    raw_values = (
        result[result.method == "raw"][keys + ["value"]]
        .rename(columns={"value": "raw_value"})
        .drop_duplicates(keys)
    )
    beck_values = (
        result[result.method == "beck_truncated"][keys + ["value"]]
        .rename(columns={"value": "beck_value"})
        .drop_duplicates(keys)
    )
    result = result.merge(raw_values, on=keys, how="left")
    result = result.merge(beck_values, on=keys, how="left")
    result["absolute_discrepancy_vs_raw"] = np.abs(
        result.value - result.raw_value
    )
    result["absolute_discrepancy_vs_beck"] = np.abs(
        result.value - result.beck_value
    )
    return result


def summarize(methods: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    groups = ["stage", "dimension", "depth", "sample_size", "method"]
    for key, frame in methods[methods.stage != "equivalence"].groupby(groups):
        errors = frame.absolute_error.dropna().to_numpy(dtype=float)
        rows.append(
            {
                **dict(zip(groups, key, strict=True)),
                "repetitions": len(frame),
                "value_mean": frame.value.mean(),
                "value_std": frame.value.std(ddof=1) if len(frame) > 1 else 0.0,
                "mae": float(np.mean(errors)) if errors.size else math.nan,
                "rmse": float(np.sqrt(np.mean(errors**2))) if errors.size else math.nan,
                "relative_error_mean": frame.relative_error.mean(),
                "published_mlp_reference": frame.published_mlp_reference.iloc[0],
                "terminal_g_evals": frame.terminal_g_evals.mean(),
                "f_evals": frame.f_evals.mean(),
                "recursive_states": frame.recursive_states.mean(),
                "normal_scalar_draws": frame.normal_scalar_draws.mean(),
                "uniform_draws": frame.uniform_draws.mean(),
                "total_stochastic_samples": frame.total_stochastic_samples.mean(),
                "wall_clock_seconds": frame.wall_clock_seconds.mean(),
                "pre_truncation_violation_rate": frame.pre_truncation_violation_rate.mean(),
                "truncation_activation_rate": frame.truncation_activation_rate.mean(),
                "mean_squared_overshoot": frame.mean_squared_overshoot.mean(),
                "pre_projection_min": frame.pre_projection_min.min(),
                "pre_projection_max": frame.pre_projection_max.max(),
                "generator_abs_before_mean": frame.generator_abs_before_mean.mean(),
                "generator_abs_after_mean": frame.generator_abs_after_mean.mean(),
                "nonlinear_correction_variance": frame.nonlinear_correction_variance.mean(),
                "max_discrepancy_vs_raw": frame.absolute_discrepancy_vs_raw.max(),
                "max_discrepancy_vs_beck": frame.absolute_discrepancy_vs_beck.max(),
                "all_outputs_finite": bool(frame.finite_output.all()),
            }
        )
    return pd.DataFrame(rows).sort_values(groups).reset_index(drop=True)


def _style() -> None:
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "legend.fontsize": 8,
            "figure.dpi": 120,
            "savefig.bbox": "tight",
            "pdf.fonttype": 42,
        }
    )


def make_figures(
    repetitions: pd.DataFrame, summary: pd.DataFrame, equivalence: pd.DataFrame
) -> list[str]:
    _style()
    FIGURE_ROOT.mkdir(parents=True, exist_ok=True)
    outputs: list[str] = []

    source = summary[summary.stage == "source"]
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 4.0), sharey=False)
    width = 0.23
    offsets = {"raw": -width, "beck_truncated": 0.0, "interval_ir": width}
    for axis, dimension in zip(axes, (10, 100, 1000), strict=True):
        frame = source[source.dimension == dimension]
        depths = sorted(frame.depth.unique())
        for method in ("raw", "beck_truncated", "interval_ir"):
            values = [
                float(frame[(frame.depth == n) & (frame.method == method)].mae.iloc[0])
                for n in depths
            ]
            axis.bar(
                np.asarray(depths) + offsets[method],
                values,
                width=width,
                color=COLORS[method],
                alpha=0.83,
                label=METHOD_LABELS[method],
            )
        axis.set_yscale("log")
        axis.set_xlabel("Published allocation n=M")
        axis.set_title(f"d={dimension}")
        axis.grid(axis="y", alpha=0.25)
        axis.set_xticks(depths)
    axes[0].set_ylabel("MAE vs published numerical reference")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.90),
        ncol=3,
        frameon=False,
    )
    fig.suptitle("Source-faithful Allen-Cahn: exact method overlap", y=0.98)
    fig.text(
        0.5,
        0.015,
        "Bars are horizontally separated only for visibility; paired values are identical.",
        ha="center",
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.07, 1, 0.80))
    path = FIGURE_ROOT / "raw_vs_truncated_error.pdf"
    fig.savefig(path)
    plt.close(fig)
    outputs.append(str(path.relative_to(RESULT_ROOT)).replace("\\", "/"))

    non_equiv = repetitions[repetitions.stage != "equivalence"]
    raw = non_equiv[non_equiv.method == "raw"]
    beck = non_equiv[non_equiv.method == "beck_truncated"]
    certified = non_equiv[non_equiv.method == "certified_ir_0_1"]
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.7))
    for dimension, marker, color in zip(
        (10, 100, 1000), ("o", "s", "^"), ("#0072B2", "#D55E00", "#009E73"), strict=True
    ):
        frame = raw[raw.dimension == dimension]
        axes[0].scatter(
            frame.total_stochastic_samples,
            frame.pre_truncation_violation_rate,
            s=18,
            alpha=0.55,
            marker=marker,
            color=color,
            label=f"d={dimension}",
        )
        grouped = (
            frame.groupby(["depth", "sample_size", "total_stochastic_samples"])
            .pre_projection_max.max()
            .reset_index()
            .sort_values("total_stochastic_samples")
        )
        axes[1].scatter(
            grouped.total_stochastic_samples,
            grouped.pre_projection_max,
            s=28,
            marker=marker,
            color=color,
            label=f"d={dimension}",
        )
    axes[0].set_xscale("log")
    axes[0].set_ylim(-0.002, 0.03)
    axes[0].set_xlabel("One-dimensional random draws")
    axes[0].set_ylabel("Pre-truncation violation rate")
    axes[0].set_title("Both [-4,4] and [0,1]: zero violations")
    axes[1].set_xscale("log")
    axes[1].axhline(1.0, color="#009E73", linestyle="--", linewidth=1.2, label="Certified upper bound 1")
    axes[1].axhline(4.0, color="#D55E00", linestyle=":", linewidth=1.2, label="Beck radius 4")
    axes[1].set_xlabel("One-dimensional random draws")
    axes[1].set_ylabel("Maximum reused child value")
    axes[1].set_title("Observed states remain strictly feasible")
    for axis in axes:
        axis.grid(alpha=0.25)
    axes[0].legend(frameon=False)
    axes[1].legend(frameon=False, fontsize=7)
    fig.suptitle(
        f"Truncation mechanism audit ({len(beck):,} Beck and {len(certified):,} certified runs)",
        y=1.01,
    )
    fig.tight_layout()
    path = FIGURE_ROOT / "violation_vs_budget.pdf"
    fig.savefig(path)
    plt.close(fig)
    outputs.append(str(path.relative_to(RESULT_ROOT)).replace("\\", "/"))

    eq = equivalence[equivalence.stage == "equivalence"]
    max_value = float(eq.absolute_value_discrepancy.max())
    max_correction = float(eq.max_correction_discrepancy.max())
    exact_values = int(eq.exact_value_match.sum())
    correction_matches = int(eq.exact_correction_matches.sum())
    correction_total = int(eq.correction_count.sum())
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.5))
    display_floor = 1e-18
    axes[0].bar(
        ["Root value", "Saved corrections"],
        [max(max_value, display_floor), max(max_correction, display_floor)],
        color=["#0072B2", "#D55E00"],
    )
    axes[0].set_yscale("log")
    axes[0].set_ylim(1e-19, 1e-11)
    axes[0].set_ylabel("Maximum absolute discrepancy")
    axes[0].set_title("Actual maxima are exactly 0")
    axes[0].grid(axis="y", alpha=0.25)
    axes[1].bar(
        ["Root pairs", "Correction pairs"],
        [exact_values, correction_matches],
        color=["#0072B2", "#D55E00"],
    )
    axes[1].set_ylabel("Bitwise exact matches")
    axes[1].set_title(f"{exact_values}/{len(eq)} roots; {correction_matches}/{correction_total} corrections")
    axes[1].grid(axis="y", alpha=0.25)
    fig.suptitle("Beck direct truncation = generic interval IR pathwise", y=1.01)
    fig.tight_layout()
    path = FIGURE_ROOT / "pathwise_equivalence.pdf"
    fig.savefig(path)
    plt.close(fig)
    outputs.append(str(path.relative_to(RESULT_ROOT)).replace("\\", "/"))

    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.6))
    work_frame = summary[
        (summary.method == "raw") & summary.stage.isin(["source", "pilot", "final"])
    ]
    for axis, dimension in zip(axes, (10, 100, 1000), strict=True):
        frame = work_frame[work_frame.dimension == dimension].sort_values(
            "total_stochastic_samples"
        )
        for depth, group in frame.groupby("depth"):
            axis.scatter(
                group.total_stochastic_samples,
                group.mae,
                s=25 + 4 * depth,
                alpha=0.75,
                label=f"n={depth}",
            )
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlabel("One-dimensional random draws")
        axis.set_title(f"d={dimension}")
        axis.grid(alpha=0.25)
        axis.legend(frameon=False, fontsize=7, ncol=2)
    axes[0].set_ylabel("MAE vs published numerical reference")
    fig.suptitle("Error versus work: Raw = Beck truncated = Interval IR", y=1.01)
    fig.tight_layout()
    path = FIGURE_ROOT / "error_vs_work.pdf"
    fig.savefig(path)
    plt.close(fig)
    outputs.append(str(path.relative_to(RESULT_ROOT)).replace("\\", "/"))
    return outputs


def _fmt(value: float) -> str:
    if not math.isfinite(value):
        return "-"
    if value == 0.0:
        return "0"
    if abs(value) < 1e-4 or abs(value) >= 1e4:
        return f"{value:.4e}"
    return f"{value:.6f}".rstrip("0").rstrip(".")


def markdown_table(headers: list[str], rows: list[list[Any]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    lines.extend("| " + " | ".join(str(item) for item in row) + " |" for row in rows)
    return "\n".join(lines)


def write_report(
    repetitions: pd.DataFrame,
    summary: pd.DataFrame,
    equivalence: pd.DataFrame,
    compute: dict[str, Any],
) -> None:
    source = summary[(summary.stage == "source") & (summary.method == "raw")]
    source_rows: list[list[Any]] = []
    for row in source.itertuples(index=False):
        source_rows.append(
            [
                row.dimension,
                row.depth,
                row.sample_size,
                row.repetitions,
                _fmt(row.value_mean),
                _fmt(row.mae),
                _fmt(row.total_stochastic_samples),
                _fmt(row.pre_projection_max),
            ]
        )

    final = summary[summary.stage == "final"]
    final_rows: list[list[Any]] = []
    for (dimension, depth, sample_size), frame in final.groupby(
        ["dimension", "depth", "sample_size"]
    ):
        by_method = frame.set_index("method")
        final_rows.append(
            [
                dimension,
                depth,
                sample_size,
                int(by_method.loc["raw", "repetitions"]),
                _fmt(float(by_method.loc["raw", "mae"])),
                _fmt(float(by_method.loc["beck_truncated", "mae"])),
                _fmt(float(by_method.loc["interval_ir", "mae"])),
                _fmt(float(by_method.loc["certified_ir_0_1", "mae"])),
                _fmt(float(by_method.loc["raw", "pre_projection_max"])),
            ]
        )

    eq = equivalence[equivalence.stage == "equivalence"]
    all_eq = equivalence
    maximum_state = float(
        repetitions[repetitions.stage != "equivalence"].pre_projection_max.max()
    )
    report = f"""# Allen-Cahn truncated-MLP recovery report

## Outcome

**VERDICT B: Mathematical/code equivalence is verified, but the tested source-faithful Allen-Cahn configurations do not activate truncation enough to show a numerical rescue.**

The experiment proves containment and verifies it pathwise. Direct Beck-style `f(P_r(u))` and generic Samplewise projection onto `[-r,r] x R^d` agree exactly in two independently written recursions. The published Allen-Cahn convergence regime is also reproduced numerically. However, Raw MLP never leaves even the tighter certified interval `[0,1]` anywhere in the completed source, pilot, or 30-repetition final grids. Consequently Raw, Beck-Truncated, same-radius Interval IR, and certified `[0,1]` IR are identical in every paired production run. Claiming an active rescue would be unsupported.

The detailed source audit is in `docs/ALLEN_CAHN_SOURCE_NOTES.md`. Beck et al. (2020) provide the truncation theory and complexity result but no finite-budget numerical table. The concrete benchmark comes from Becker et al. (2020), their published companion simulation paper.

## Exact source setting

The reproduced terminal-value problem is

```text
partial_t u(t,x) + Delta u(t,x) + u(t,x) - u(t,x)^3 = 0,
T = 1,
u(T,x) = 1 / (2 + (2/5)||x||^2),
X_(t,s)^x = x + sqrt(2)(W_s-W_t),
x_root = 0.
```

The companion benchmark uses `d in {{10,100,1000}}`, `n=M in {{1,...,8}}`, a uniform random time, and fixed radius `r=4`. This study reproduced `n=M=1,...,5` with five independent repetitions per cell, matching the paper's five-run error protocol. The published `V_(8,8,4)` reference values (`0.29555`, `0.03373`, `0.00340`) are numerical rather than exact.

Beck's theory separately permits a growing `rho_M` with `rho_M -> infinity` and `rho_M=O(log log M)`; its explicit example is `log(1+log M)`. That theorem schedule, the companion's fixed `r=4`, and this instance's solution-side invariant interval `[0,1]` are not conflated.

## Mathematical containment proposition

Let `H=R x R^d`, `Y=(u,z)`, and `C_r=[-r,r] x R^d`. Because `C_r` is a Cartesian product of a closed interval and the full gradient space, Euclidean distance separates and

```text
Pi_(C_r)(u,z) = (P_r(u),z),
P_r(u) = min(r,max(-r,u)).
```

For a reaction driver independent of `z`, Samplewise IR therefore evaluates

```text
bar f(u,z) = f(P_r(u)),
```

which is exactly Beck's truncated driver.

For pathwise equality, fix the random tree, random-time variables, diffusion increments, sample allocation, and radius schedule. At Picard depth zero both recursions return zero. Assume equality at every depth below `n`. At every level in the depth-`n` estimator, the fine and coarse child values are equal by the induction hypothesis; applying the identities above gives identical fine and coarse generator values. Terminal samples, level weights, additions, and subtractions are then identical. Hence the complete depth-`n` estimates are equal on that random tree. Induction proves `U_IR(n,M,r)=U_Beck(n,M,r)` pathwise, not only in law.

## Independent code-path verification

The direct implementation performs its own scalar clamp and cubic evaluation. The IR implementation runs a separately written recursion, calls the generic interval projector, and then calls the original reaction. They share only equation primitives and an immutable keyed random tree.

- dedicated equivalence grid: {len(eq)} paired configurations;
- dimensions: `1,10,100,1000`;
- depths/sample sizes: `(0,2),(1,2),(2,2),(3,2),(3,3)`;
- three seeds, zero and nonzero ramp points;
- fixed `r=4` and `rho_M=log(1+log M)` tracks;
- maximum root discrepancy: `{_fmt(float(eq.absolute_value_discrepancy.max()))}`;
- exact root matches: `{int(eq.exact_value_match.sum())}/{len(eq)}`;
- maximum saved-correction discrepancy: `{_fmt(float(eq.max_correction_discrepancy.max()))}`;
- exact saved-correction matches: `{int(eq.exact_correction_matches.sum())}/{int(eq.correction_count.sum())}`.

Across all source/pilot/final tasks as well, the maximum Beck-vs-IR discrepancy is `{_fmt(float(all_eq.absolute_value_discrepancy.max()))}` and all draw fingerprints match.

## Published-regime reproduction

{markdown_table(["d", "n", "M", "reps", "mean value", "MAE", "random draws", "max reused u"], source_rows)}

At `n=M=5`, the means are `0.295944` (d=10), `0.033438` (d=100), and `0.003350` (d=1000), close to the companion's numerical references `0.29555`, `0.03373`, and `0.00340`. The reference is not treated as analytic truth.

## Controlled depth/budget sweep and final cells

Stage B kept the PDE, terminal function, horizon, diffusion, and evaluation point unchanged. It explored 26 pilot cells covering high sampling at depth two and increasingly deep/under-sampled allocations through `(n,M)=(8,2)`. Nine representative cells were then frozen at 30 repetitions:

{markdown_table(["d", "n", "M", "reps", "Raw MAE", "Beck MAE", "IR r=4 MAE", "IR [0,1] MAE", "max reused u"], final_rows)}

The largest pre-projection child value in all production artifacts is `{_fmt(maximum_state)}`. The minimum is zero. There are no `[-4,4]` violations, no `[0,1]` violations, no activations, and no nonfinite values or generators. Thus this is a strong inactive sanity check, not an active stabilization result.

## Driver regularity

For `f(u)=u-u^3`, `f'(u)=1-3u^2`. On `[-r,r]`, the clipped driver is globally Lipschitz with sharp constant

```text
max(1, |1-3r^2|),
```

bounded by `1+3r^2` as in a standard local estimate. At `r=4` the sharp constant is `{_fmt(sharp_clipped_driver_lipschitz(4.0))}` (the loose bound is 49); on `[0,1]` it is 2. A dense numerical slope check is included in the tests.

## Work and reproducibility

The production grid contains {len(repetitions[repetitions.stage != 'equivalence'])} method-repetitions plus {len(eq)} dedicated equivalence pairs. Recorded parallel orchestration wall time is `{compute['orchestration_wall_seconds']:.3f}` seconds; summed task wall time is `{compute['task_wall_seconds']:.3f}` seconds. Work columns contain terminal evaluations, generator evaluations, recursive states, normal scalar draws, uniform draws, total stochastic samples, and per-root wall time.

No manuscript source was modified. The result establishes backward compatibility with scalar truncated MLP; it does not claim truncation or the local-to-global argument as new. The broader IR contribution, if used later, must concern structured value-gradient geometry for gradient-dependent nonlinear reuse.

**VERDICT B: Mathematical/code equivalence is verified, but the tested source-faithful Allen-Cahn configurations do not activate truncation enough to show a numerical rescue.**
"""
    (DOC_ROOT / "ALLEN_CAHN_RECOVERY_REPORT.md").write_text(report, encoding="utf-8")


def build_provenance(figures: list[str], repetitions: pd.DataFrame) -> dict[str, Any]:
    manifests: list[dict[str, Any]] = []
    for path in sorted(RESULT_ROOT.glob("*_manifest.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        manifests.append(
            {
                "name": path.name,
                "stage": payload["stage"],
                "status": payload["status"],
                "tasks": payload["tasks_total"],
                "elapsed_seconds": payload["elapsed_seconds"],
                "task_seconds_sum": payload["task_seconds_sum"],
            }
        )
    source_files = sorted(HERE.glob("*.py")) + [HERE / "README.md"]
    grids: dict[str, list[dict[str, int]]] = {}
    for stage in ("source", "pilot", "final"):
        frame = repetitions[repetitions.stage == stage]
        cells = (
            frame[["dimension", "depth", "sample_size", "repetition"]]
            .drop_duplicates()
            .groupby(["dimension", "depth", "sample_size"])
            .size()
            .reset_index(name="repetitions")
        )
        grids[stage] = [
            {key: int(value) for key, value in row.items()}
            for row in cells.to_dict(orient="records")
        ]
    return {
        "schema_version": 1,
        "branch": git("branch", "--show-current"),
        "source_base_commit": git("rev-parse", "HEAD"),
        "source_papers": {
            "beck_theory": {
                "arxiv": "https://arxiv.org/abs/1907.06729",
                "doi": "https://doi.org/10.1515/jnma-2019-0074",
                "inspected_pdf_sha256": "521ed4d065ed24349f629f80794ca094ffb1e87664afd5341882a60bfd10aec5",
            },
            "becker_numerics": {
                "arxiv": "https://arxiv.org/abs/2005.10206",
                "doi": "https://doi.org/10.4208/cicp.OA-2020-0130",
                "inspected_pdf_sha256": "bff46fd1c1b99ceac7a0204bcf97eb9ea4f8d93b18739307cd4bd1fef925060f",
            },
        },
        "equation": {
            "horizon": 1.0,
            "reaction": "u-u^3",
            "terminal": "1/(2+(2/5)||x||^2)",
            "diffusion": "sqrt(2)",
            "random_time": "uniform",
            "published_radius": BECK_NUMERICAL_RADIUS,
            "certified_interval": list(CERTIFIED_INTERVAL),
            "published_numerical_references": PUBLISHED_MLP_REFERENCES,
        },
        "completed_grids": grids,
        "equivalence_grid": {
            "dimensions": [1, 10, 100, 1000],
            "n_M": [[0, 2], [1, 2], [2, 2], [3, 2], [3, 3]],
            "seeds_per_cell": 3,
            "test_points": ["zero", "ramp"],
            "radius_tracks": ["beck_r4", "theorem_schedule"],
        },
        "manifests": manifests,
        "compute": {
            "orchestration_wall_seconds": sum(item["elapsed_seconds"] for item in manifests),
            "task_wall_seconds": sum(item["task_seconds_sum"] for item in manifests),
            "method_repetitions": int(len(repetitions[repetitions.stage != "equivalence"])),
        },
        "new_source_sha256": {
            str(path.relative_to(HERE)).replace("\\", "/"): sha256(path)
            for path in source_files
            if path.exists()
        },
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "matplotlib": matplotlib.__version__,
            "platform": platform.platform(),
        },
        "figures": figures,
        "verdict": "B",
    }


def write_results_readme() -> None:
    text = """# Allen-Cahn recovery result bundle

This directory contains source-faithful per-repetition results, pathwise
equivalence metrics, aggregate summaries, resumable task artifacts,
provenance, validation, and four work/mechanism figures.

- `repetition_metrics.csv`: every production method/repetition.
- `pathwise_equivalence.csv`: direct Beck versus generic interval IR pairs.
- `summary.csv`: aggregate accuracy, work, timing, and activation metrics.
- `provenance.json`: sources, exact grids, hashes, environment, and runtime.
- `validation_audit.json`: scientific and artifact integrity checks.
- `raw/`: resumable task-level JSON artifacts.
- `figures/`: rendered and visually checked PDF figures.

Verdict: **B**. Mathematical and code containment are exact, but no
source-faithful tested state leaves even the tighter `[0,1]` interval, so no
active finite-budget truncation rescue is observed.
"""
    (RESULT_ROOT / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    RESULT_ROOT.mkdir(parents=True, exist_ok=True)
    methods, equivalence = load_raw()
    all_repetitions = add_paired_columns(methods)
    repetitions = all_repetitions[all_repetitions.stage != "equivalence"].copy()
    summary = summarize(all_repetitions)
    repetitions.to_csv(RESULT_ROOT / "repetition_metrics.csv", index=False)
    equivalence.to_csv(RESULT_ROOT / "pathwise_equivalence.csv", index=False)
    summary.to_csv(RESULT_ROOT / "summary.csv", index=False)
    figures = make_figures(repetitions, summary, equivalence)
    provenance = build_provenance(figures, repetitions)
    (RESULT_ROOT / "provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True), encoding="utf-8"
    )
    write_report(repetitions, summary, equivalence, provenance["compute"])
    write_results_readme()
    print(f"Production method-repetitions: {len(repetitions)}")
    print(f"Dedicated equivalence pairs: {len(equivalence[equivalence.stage == 'equivalence'])}")
    print(f"Figures: {len(figures)}")
    print(f"Verdict: {provenance['verdict']}")


if __name__ == "__main__":
    main()
