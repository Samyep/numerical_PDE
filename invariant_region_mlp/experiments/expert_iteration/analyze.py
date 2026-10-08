"""Mechanical analysis and report generation for the Study-A decision subset."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import subprocess
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
RESULTS = PROJECT / "results" / "expert_iteration"
DOC = PROJECT / "docs" / "EXPERT_ITERATION_REPORT.md"
FIGURES = RESULTS / "figures"
EXPECTED_ALL_DIMS = {20, 50, 100, 120, 140, 160}
EXPECTED_HIGH_DIMS = {100, 120, 140, 160}


CRITERIA = {
    "A-1": "`mlp` relative L2 >= 1 at d >= 100, and `scasml` / `surrogate` error ratio >= 0.5 at d = 160 in the primary cell. (If not reproduced, report and continue.)",
    "A-2": "In the primary cell, `path` <= 0.2 x `surrogate` and `path` <= 0.5 x `scasml` at every d in {100, ..., 160}.",
    "A-3": "The ratio `path / surrogate` at d = 160 is <= 1.3 x its value at d = 20, while the ratio `scasml_noclip / surrogate` at d = 160 is >= 1.5 x its value at d = 20.",
    "A-4": "Mean generator bias for `scasml_noclip` is negative at every d and its magnitude at d = 160 is >= 3 x that at d = 20; for `path` the magnitude at d = 160 is <= 1.5 x that at d = 20.",
    "A-5": "Plot `path` and `scasml_noclip` error against surrogate error across the three checkpoints; report whether `scasml_noclip` is ever worse than the surrogate alone (exploratory).",
    "B-pilot": "Proceed only if `EI-path` test error after round 3 is <= 0.5 x after round 0 on B2.",
    "B-1": "`EI-path` error after round 3 <= 0.2 x its error after round 0, on B1, B2, B3.",
    "B-2": "On B4, `EI-bismut` error after round 3 >= its error after round 1, while `EI-path` satisfies B-1 on B4.",
    "B-3": "At the wall-clock of `EI-path` round 6, `EI-path` error <= 0.5 x the best DPI error within the same wall-clock, on at least 2 of B1-B3.",
    "B-4": "Report whether DPI at 4x the wall-clock reaches `EI-path`'s round-6 error (no verdict).",
    "B-5": "For the final `EI-path` network, corrected error <= 0.3 x network error on B1-B3.",
}


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=PROJECT, text=True).strip()


def _atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=True), encoding="utf-8")
    temporary.replace(path)


def load_rows() -> pd.DataFrame:
    path = RESULTS / "A_rows.csv"
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    frame = pd.read_csv(path)
    for column in (
        "value_relative_l2",
        "generator_bias",
        "generator_calls",
        "wall_clock_seconds",
    ):
        if column in frame:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return frame


def standard_rows(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return frame
    # d=100 exact rows are the additional control; d=20 exact rows are the
    # preregistered standard because d<=50 uses the exact Laplacian.
    return frame[~((frame["dimension"] == 100) & (frame["laplacian"] == "exact"))].copy()


def build_summary(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame()
    grouping = ["dimension", "checkpoint", "n", "M", "method", "laplacian"]
    summary = (
        frame.groupby(grouping, dropna=False)
        .agg(
            rows=("value_relative_l2", "size"),
            network_seeds=("network_seed", "nunique"),
            repetitions=("repetition", "nunique"),
            median_value_relative_l2=("value_relative_l2", "median"),
            mean_value_relative_l2=("value_relative_l2", "mean"),
            median_skill=("skill", "median"),
            median_gradient_relative_l2=("gradient_relative_l2", "median"),
            mean_generator_bias=("generator_bias", "mean"),
            total_generator_calls=("generator_calls", "sum"),
            total_terminal_samples=("terminal_samples", "sum"),
            total_wall_clock_seconds=("wall_clock_seconds", "sum"),
            total_nonfinite_states=("nonfinite_state_count", "sum"),
            total_nonfinite_generators=("nonfinite_generator_values", "sum"),
        )
        .reset_index()
        .sort_values(grouping)
    )
    summary.to_csv(RESULTS / "A_summary.csv", index=False)
    return summary


def network_level(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return frame
    group = ["dimension", "network_seed", "checkpoint", "n", "M", "method", "laplacian"]
    return (
        frame.groupby(group, dropna=False)
        .agg(
            error=("value_relative_l2", "median"),
            generator_bias=(
                "generator_bias",
                lambda values: float(np.nanmean(values)) if np.any(np.isfinite(values)) else np.nan,
            ),
            rows=("value_relative_l2", "size"),
        )
        .reset_index()
    )


def _paired_ratios(network: pd.DataFrame, numerator: str, denominator: str) -> pd.DataFrame:
    keys = ["dimension", "network_seed", "checkpoint", "n", "M"]
    left = network[network.method == numerator][keys + ["error"]].rename(columns={"error": "num"})
    right = network[network.method == denominator][keys + ["error"]].rename(columns={"error": "den"})
    result = left.merge(right, on=keys, how="inner")
    result["ratio"] = result["num"] / result["den"]
    return result


def evaluate(frame: pd.DataFrame, references: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    standard = standard_rows(frame)
    net = network_level(standard)
    primary = net[(net.checkpoint == 2500) & (net.n == 2) & (net.M == 10)]
    available_dims = set(map(int, primary.dimension.unique())) if not primary.empty else set()
    outcome: dict[str, Any] = {}

    g1_values = {d: references.get(str(d), {}).get("gate_A_G1") for d in sorted(map(int, references))}
    outcome["A-G1"] = {
        "verdict": "PASS" if g1_values and all(value is True for value in g1_values.values()) else "FAIL",
        "details": g1_values,
    }
    outcome["A-G0"] = {
        "verdict": "REPORTED",
        "details": {
            d: references.get(str(d), {}).get("f_zero_relative_l2_test")
            for d in sorted(map(int, references))
        },
    }

    mlp = primary[primary.method == "mlp"]
    mlp_medians = mlp.groupby("dimension").error.median().to_dict()
    sc_ratio = _paired_ratios(primary, "scasml", "surrogate")
    sc160 = sc_ratio[sc_ratio.dimension == 160].ratio.median() if np.any(sc_ratio.dimension == 160) else np.nan
    a1_restricted = bool(
        all(mlp_medians.get(d, -np.inf) >= 1.0 for d in sorted(available_dims & {100, 120, 140, 160}))
        and np.isfinite(sc160)
        and sc160 >= 0.5
    )
    outcome["A-1"] = {
        "verdict": "PASS" if EXPECTED_HIGH_DIMS.issubset(available_dims) and a1_restricted else (
            "FAIL" if EXPECTED_HIGH_DIMS.issubset(available_dims) else "NOT EVALUATED (USER-DIRECTED SUBSET)"
        ),
        "restricted_subset_pass": a1_restricted,
        "mlp_median_errors": {str(k): float(v) for k, v in mlp_medians.items()},
        "scasml_over_surrogate_d160": None if not np.isfinite(sc160) else float(sc160),
        "missing_dimensions": sorted(EXPECTED_HIGH_DIMS - available_dims),
    }

    path_sur = _paired_ratios(primary, "path", "surrogate")
    path_sca = _paired_ratios(primary, "path", "scasml")
    high_path_sur = path_sur[path_sur.dimension.isin(EXPECTED_HIGH_DIMS)].groupby("dimension").ratio.median().to_dict()
    high_path_sca = path_sca[path_sca.dimension.isin(EXPECTED_HIGH_DIMS)].groupby("dimension").ratio.median().to_dict()
    a2_restricted = bool(
        high_path_sur
        and all(value <= 0.2 for value in high_path_sur.values())
        and all(value <= 0.5 for value in high_path_sca.values())
    )
    outcome["A-2"] = {
        "verdict": "PASS" if EXPECTED_HIGH_DIMS.issubset(available_dims) and a2_restricted else (
            "FAIL" if EXPECTED_HIGH_DIMS.issubset(available_dims) else "NOT EVALUATED (USER-DIRECTED SUBSET)"
        ),
        "restricted_subset_pass": a2_restricted,
        "path_over_surrogate": {str(k): float(v) for k, v in high_path_sur.items()},
        "path_over_scasml": {str(k): float(v) for k, v in high_path_sca.items()},
        "missing_dimensions": sorted(EXPECTED_HIGH_DIMS - available_dims),
    }

    noclip_sur = _paired_ratios(primary, "scasml_noclip", "surrogate")
    path_ratio_by_d = path_sur.groupby("dimension").ratio.median().to_dict()
    noclip_ratio_by_d = noclip_sur.groupby("dimension").ratio.median().to_dict()
    a3_complete = {20, 160}.issubset(available_dims)
    path_growth = path_ratio_by_d.get(160, np.nan) / path_ratio_by_d.get(20, np.nan)
    noclip_growth = noclip_ratio_by_d.get(160, np.nan) / noclip_ratio_by_d.get(20, np.nan)
    a3_pass = bool(a3_complete and path_growth <= 1.3 and noclip_growth >= 1.5)
    outcome["A-3"] = {
        "verdict": "PASS" if a3_pass else ("FAIL" if a3_complete else "NOT EVALUATED (MISSING DATA)"),
        "path_ratio_growth_d160_over_d20": None if not np.isfinite(path_growth) else float(path_growth),
        "scasml_noclip_ratio_growth_d160_over_d20": None if not np.isfinite(noclip_growth) else float(noclip_growth),
    }

    bias_rows = primary[primary.method.isin(["path", "scasml_noclip"])]
    bias = bias_rows.groupby(["dimension", "method"]).generator_bias.median().unstack()
    bias_payload = {
        str(int(d)): {method: float(value) for method, value in row.dropna().items()}
        for d, row in bias.iterrows()
    }
    present_bias_dims = set(map(int, bias.index))
    if {20, 160}.issubset(present_bias_dims) and {"path", "scasml_noclip"}.issubset(bias.columns):
        no_growth = abs(bias.loc[160, "scasml_noclip"]) / max(abs(bias.loc[20, "scasml_noclip"]), 1e-300)
        path_bias_growth = abs(bias.loc[160, "path"]) / max(abs(bias.loc[20, "path"]), 1e-300)
        restricted_a4 = bool(
            np.all(bias.loc[:, "scasml_noclip"] < 0)
            and no_growth >= 3.0
            and path_bias_growth <= 1.5
        )
    else:
        no_growth = path_bias_growth = np.nan
        restricted_a4 = False
    outcome["A-4"] = {
        "verdict": "PASS" if EXPECTED_ALL_DIMS.issubset(present_bias_dims) and restricted_a4 else (
            "FAIL" if EXPECTED_ALL_DIMS.issubset(present_bias_dims) else "NOT EVALUATED (USER-DIRECTED SUBSET)"
        ),
        "restricted_subset_pass": restricted_a4,
        "median_generator_bias": bias_payload,
        "scasml_noclip_bias_growth": None if not np.isfinite(no_growth) else float(no_growth),
        "path_bias_growth": None if not np.isfinite(path_bias_growth) else float(path_bias_growth),
        "missing_dimensions": sorted(EXPECTED_ALL_DIMS - present_bias_dims),
    }

    quality = net[
        (net.n == 2)
        & (net.M == 10)
        & net.method.isin(["surrogate", "path", "scasml_noclip"])
    ]
    wide = quality.pivot_table(
        index=["dimension", "network_seed", "checkpoint"], columns="method", values="error"
    ).reset_index()
    comparisons = int(np.count_nonzero(wide.get("scasml_noclip", np.nan) > wide.get("surrogate", np.nan)))
    outcome["A-5"] = {
        "verdict": "EXPLORATORY",
        "scasml_noclip_worse_count": comparisons,
        "comparisons": len(wide),
        "ever_worse": bool(comparisons > 0),
    }

    exact = frame[
        (frame.dimension == 100)
        & (frame.network_seed == 0)
        & (frame.checkpoint == 2500)
        & (frame.n == 2)
        & (frame.M == 10)
        & frame.method.isin(["path", "scasml_noclip"])
    ]
    exact_control: dict[str, Any] = {"status": "NOT EVALUATED"}
    if {"exact", "hutch25"}.issubset(set(exact.laplacian)):
        table = exact.groupby(["method", "laplacian"]).value_relative_l2.median().unstack()
        exact_control = {"status": "REPORTED"}
        exact_control.update({
            str(method): {
                "exact": float(row["exact"]),
                "hutchinson": float(row["hutch25"]),
                "relative_difference": float(abs(row["exact"] - row["hutch25"]) / row["exact"]),
            }
            for method, row in table.iterrows()
        })

    for key in ("B-pilot", "B-1", "B-2", "B-3", "B-4", "B-5"):
        outcome[key] = {
            "verdict": "NOT EVALUATED (USER-DIRECTED DECISION SUBSET)",
            "details": "The user explicitly requested stopping after the inexpensive Study-A decision point.",
        }
    return outcome, exact_control


def make_figures(frame: pd.DataFrame) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    standard = standard_rows(frame)
    net = network_level(standard)
    final = net[(net.checkpoint == 2500) & (net.n == 2) & (net.M == 10)]
    plt.figure(figsize=(8, 5))
    for method in ("surrogate", "f_zero", "mlp", "mlp_clip", "scasml", "scasml_noclip", "path", "path_clip", "oracle_state"):
        subset = final[final.method == method]
        if subset.empty:
            continue
        grouped = subset.groupby("dimension").error.median().sort_index()
        plt.plot(grouped.index, grouped.values, marker="o", label=method)
    plt.yscale("log")
    plt.xlabel("dimension d")
    plt.ylabel("test relative L2 error of u")
    plt.title("Study A: primary cell (n=2, M=10), checkpoint 2500")
    plt.grid(True, which="both", alpha=0.25)
    plt.legend(ncol=2, fontsize=8)
    plt.tight_layout()
    plt.savefig(FIGURES / "A_error_vs_d.png", dpi=180)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    for method in ("scasml_noclip", "path"):
        subset = final[final.method == method]
        if subset.empty:
            continue
        grouped = subset.groupby("dimension").generator_bias.median().sort_index()
        plt.plot(grouped.index, grouped.values, marker="o", label=method)
    plt.axhline(0.0, color="black", linewidth=0.8)
    plt.xlabel("dimension d")
    plt.ylabel("mean generator bias")
    plt.title("Noise entering the defect generator")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIGURES / "A_generator_bias_vs_d.png", dpi=180)
    plt.close()

    quality = net[
        (net.n == 2)
        & (net.M == 10)
        & net.method.isin(["surrogate", "path", "scasml_noclip"])
    ]
    wide = quality.pivot_table(
        index=["dimension", "network_seed", "checkpoint"], columns="method", values="error"
    ).reset_index()
    plt.figure(figsize=(6.5, 5))
    for method, marker in (("path", "o"), ("scasml_noclip", "s")):
        if method not in wide:
            continue
        plt.scatter(wide["surrogate"], wide[method], label=method, marker=marker, alpha=0.8)
    if len(wide):
        low = max(min(wide["surrogate"].min(), wide.get("path", wide["surrogate"]).min()), 1e-6)
        high = max(wide["surrogate"].max(), wide.get("scasml_noclip", wide["surrogate"]).max())
        plt.plot([low, high], [low, high], "k--", linewidth=0.8, label="no improvement")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("surrogate relative L2 error")
    plt.ylabel("corrected relative L2 error")
    plt.title("A-5: correction versus surrogate quality")
    plt.grid(True, which="both", alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIGURES / "A_error_vs_surrogate.png", dpi=180)
    plt.close()


def _write_empty_b_outputs() -> None:
    schemas = {
        "B_tuning.csv": ["status", "reason"],
        "B_pilot.csv": ["status", "reason"],
        "B_rounds.csv": ["status", "reason"],
        "B_summary.csv": ["criterion", "verdict", "reason"],
    }
    reason = "not evaluated: user-directed Study-A-only decision subset"
    for name, fields in schemas.items():
        path = RESULTS / name
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            if fields[0] == "criterion":
                for criterion in ("B-pilot", "B-1", "B-2", "B-3", "B-4", "B-5"):
                    writer.writerow({"criterion": criterion, "verdict": "NOT EVALUATED", "reason": reason})
            else:
                writer.writerow({"status": "NOT EVALUATED", "reason": reason})
    _atomic_json(
        RESULTS / "tuning_choices.json",
        {"status": "NOT EVALUATED", "reason": reason, "choices": []},
    )


def write_report(
    frame: pd.DataFrame,
    summary: pd.DataFrame,
    outcomes: dict[str, Any],
    exact_control: dict[str, Any],
    references: dict[str, Any],
) -> None:
    prereg_commit = _git("log", "--format=%H", "--reverse", "--", str(HERE / "PREREGISTRATION.md")).splitlines()[0]
    code_commit = _git("rev-parse", "HEAD")
    environment = json.loads((RESULTS / "environment.json").read_text(encoding="utf-8")) if (RESULTS / "environment.json").exists() else {}
    total_wall = float(frame.wall_clock_seconds.sum()) if not frame.empty else 0.0
    lines = [
        "# Noise-aware defect correction and expert iteration: decision-subset report",
        "",
        "> Status: Study A decision subset only. The user explicitly limited this run to d={20,100,160}; "
        "Study B and omitted Study-A dimensions/cells are not evaluated, not negative results.",
        "",
        "## Outcome table",
        "",
        "| Gate / criterion | Verdict |",
        "|---|---|",
    ]
    for key, value in outcomes.items():
        lines.append(f"| {key} | {value['verdict']} |")
    lines += ["", "## Gates and criteria (verbatim)", ""]
    lines.append(f"- **A-G0:** Report the relative L2 error of the f=0 solution. Result: `{outcomes['A-G0']['details']}`.")
    lines.append(f"- **A-G1:** MC estimated SE of u <=1e-4 at every test point. Verdict: **{outcomes['A-G1']['verdict']}**. Details: `{outcomes['A-G1']['details']}`.")
    for key in ("A-1", "A-2", "A-3", "A-4", "A-5", "B-pilot", "B-1", "B-2", "B-3", "B-4", "B-5"):
        lines += [
            "",
            f"### {key}",
            "",
            CRITERIA[key],
            "",
            f"Verdict: **{outcomes[key]['verdict']}**.",
            "",
            f"Mechanical details: `{json.dumps(outcomes[key], sort_keys=True, allow_nan=True)}`",
        ]
    lines += [
        "",
        "## Exact-Laplacian control at d=100",
        "",
        f"`{json.dumps(exact_control, sort_keys=True, allow_nan=True)}`",
        "",
        "## Failures, limitations, and deviations",
        "",
        "- The run is intentionally the user-requested decision subset d={20,100,160}. Dimensions 50, 120, and 140, secondary cells, and all of Study B are marked not evaluated.",
        "- The frozen prompt's interpretation is lambda=1, zero drift, T=0.5, and unit-ball test points. The audited public SCaSML LQG markdown says T=1; this run follows the frozen prompt.",
        "- The public SCaSML training code samples d/4 coordinate Hessian entries and also subsamples the gradient norm. The frozen prompt specifically says Hutchinson d/4 for the Laplacian; this implementation uses Rademacher Hutchinson for the Laplacian and the full gradient norm.",
        "- Reported test errors use the required adaptive antithetic MC reference. A scaled Gauss--Laguerre Hopf--Cole integral supplies truth only for recursive child-state diagnostics/oracle calls and is audited against the MC reference.",
        "- Antithetic reference points use common random numbers within each eight-point chunk. Each point retains the correct marginal iid Gaussian sample and its own SE; correlations across test-point errors do not enter any registered threshold.",
        "- The pathwise change is restricted to the terminal gradient estimator; level terms retain the unchanged Bismut estimator.",
        "",
        "## Provenance and totals",
        "",
        f"- Frozen preregistration commit: `{prereg_commit}`",
        f"- Code commit at analysis: `{code_commit}`",
        f"- Rows: {len(frame)}; aggregate inference wall-clock: {total_wall:.1f} s",
        f"- Environment: `{json.dumps(environment, sort_keys=True)}`",
        f"- Reference summary: `{json.dumps(references, sort_keys=True, allow_nan=True)}`",
        "",
        "## Figures",
        "",
        "- `results/expert_iteration/figures/A_error_vs_d.png`",
        "- `results/expert_iteration/figures/A_generator_bias_vs_d.png`",
        "- `results/expert_iteration/figures/A_error_vs_surrogate.png`",
    ]
    DOC.parent.mkdir(parents=True, exist_ok=True)
    temporary = DOC.with_suffix(".tmp.md")
    temporary.write_text("\n".join(lines) + "\n", encoding="utf-8")
    temporary.replace(DOC)


def main() -> None:
    frame = load_rows()
    summary = build_summary(frame)
    references = json.loads((RESULTS / "A_reference.json").read_text(encoding="utf-8")) if (RESULTS / "A_reference.json").exists() else {}
    outcomes, exact_control = evaluate(frame, references)
    make_figures(frame)
    _write_empty_b_outputs()
    analysis = {
        "scope": "Study A decision subset d={20,100,160}",
        "outcomes": outcomes,
        "exact_laplacian_control": exact_control,
        "rows": len(frame),
        "reference_dimensions": sorted(map(int, references)),
    }
    _atomic_json(RESULTS / "A_gates.json", outcomes)
    _atomic_json(RESULTS / "analysis_summary.json", analysis)
    write_report(frame, summary, outcomes, exact_control, references)
    print(json.dumps(analysis, indent=2, sort_keys=True, allow_nan=True))


if __name__ == "__main__":
    main()
