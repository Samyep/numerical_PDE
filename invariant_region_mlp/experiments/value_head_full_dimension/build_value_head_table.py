#!/usr/bin/env python3
"""Build the frozen direct-value-head comparison across dimensions.

This script performs no training and no checkpoint selection.  It reuses the
official, same-bank whole-function results and audits that every comparison in
a PDE/dimension cell uses the same frozen state bank.  The one value-head
metric absent from MAIN_RESULTS.csv (100D nonconvex_coupled) is recovered from
the selected-a4 post-freeze evaluation on the byte-identical canonical bank.
"""

from __future__ import print_function

import csv
import hashlib
import io
import json
import subprocess
from collections import Counter, OrderedDict
from pathlib import Path

import numpy as np


REPO = Path(__file__).resolve().parents[3]
RESULTS = REPO / "invariant_region_mlp/results/value_head_full_dimension"
TABLES = RESULTS / "tables"

FINAL_CSV = (
    REPO
    / "invariant_region_mlp/results/tangent_actor_icml_final_suite/MAIN_RESULTS.csv"
)
ACCURACY_ROOT = (
    REPO / "invariant_region_mlp/results/tangent_actor_accuracy_frontier"
)
SELECTION_JSON = (
    ACCURACY_ROOT / "selection_coupled_d100_k2_lam0p2_eta1_main_mcv_s0.json"
)
POSTFREEZE_JSON = (
    ACCURACY_ROOT / "post_freeze_coupled_d100_k2_lam0p2_eta1_main_mcv_s0.json"
)
POSTFREEZE_NPZ = (
    ACCURACY_ROOT
    / "predictions/post_freeze_coupled_d100_k2_lam0p2_eta1_main_mcv_s0.npz"
)
PRECISION_NPZ = (
    ACCURACY_ROOT
    / "predictions/selected_precision_coupled_d100_k2_lam0p2_eta1_main_mcv_s0_paths2048_fine512.npz"
)

EXPANSION_REF = "origin/ir-mlp-icml-two-highdim-expansion-login-ice-gnr-2"
EXPANSION_CSV = (
    "invariant_region_mlp/results/main_highdim_expansion/"
    "main_whole_function_results.csv"
)

DIMENSIONS = (20, 50, 100, 200, 400)
PDE_ORDER = (
    "allen_cahn",
    "joint_concave",
    "joint_convex",
    "nonconvex_coupled",
    "dense_fourier",
)
PDE_LABELS = OrderedDict(
    [
        ("allen_cahn", "Allen--Cahn"),
        ("joint_concave", "Concave Cubic--Quadratic"),
        ("joint_convex", "Convex Cubic--Quadratic"),
        ("nonconvex_coupled", "Oscillatory Value--Gradient"),
        ("dense_fourier", "Fourier--Tanh Coupled"),
    ]
)
METHODS = ("dpi", "pinn_hte", "sdgd_pinn")
METHOD_LABELS = OrderedDict(
    [
        ("direct_value_head", "Direct value head"),
        ("dpi", "DPI"),
        ("pinn_hte", "PINN-HTE"),
        ("sdgd_pinn", "SDGD-PINN"),
    ]
)


def sha256_file(path):
    h = hashlib.sha256()
    with open(str(path), "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_csv(path):
    with open(str(path), newline="") as handle:
        return list(csv.DictReader(handle))


def git_output(*args):
    return subprocess.check_output(
        ["git", "-C", str(REPO)] + list(args), universal_newlines=True
    ).strip()


def remote_expansion_rows():
    payload = git_output("show", EXPANSION_REF + ":" + EXPANSION_CSV)
    return list(csv.DictReader(io.StringIO(payload))), payload


def ensure_same_bank(rows, pde, dimension):
    bank_hashes = {r["test_bank_array_sha256"] for r in rows}
    state_counts = {int(float(r["states"])) for r in rows}
    if len(bank_hashes) != 1 or len(state_counts) != 1:
        raise RuntimeError(
            "Unpaired bank for {} d={}: hashes={}, states={}".format(
                pde, dimension, sorted(bank_hashes), sorted(state_counts)
            )
        )
    return next(iter(bank_hashes)), next(iter(state_counts))


def recover_nonconvex_d100_head():
    selection = json.load(open(str(SELECTION_JSON)))
    if selection["selected_actor"] != "a4" or selection["selected_iteration"] != 4:
        raise RuntimeError("Unexpected selected nonconvex d100 actor")

    postfreeze = json.load(open(str(POSTFREEZE_JSON)))
    selected_rows = [
        row for row in postfreeze["iteration_rows"] if row["iteration"] == "a4"
    ]
    if len(selected_rows) != 1:
        raise RuntimeError("Expected exactly one a4 post-freeze row")

    with np.load(str(POSTFREEZE_NPZ), allow_pickle=False) as post, np.load(
        str(PRECISION_NPZ), allow_pickle=False
    ) as precision:
        equality = {
            key: bool(np.array_equal(post[key], precision[key]))
            for key in ("test_t", "test_x", "truth")
        }
    if not all(equality.values()):
        raise RuntimeError("Post-freeze and canonical precision banks differ")

    return float(selected_rows[0]["actor_value_rmse"]), equality


def finite_float(value, field, pde, dimension):
    if value in (None, ""):
        raise RuntimeError("Missing {} for {} d={}".format(field, pde, dimension))
    result = float(value)
    if not np.isfinite(result):
        raise RuntimeError("Nonfinite {} for {} d={}".format(field, pde, dimension))
    return result


def fmt_sci(value):
    mantissa, exponent = "{:.2e}".format(float(value)).split("e")
    return "{}e{:+d}".format(mantissa, int(exponent))


def latex_sci(value):
    mantissa, exponent = "{:.2e}".format(float(value)).split("e")
    return r"{}\mathrm{{e}}{{{:+d}}}".format(mantissa, int(exponent))


def build_rows():
    final_rows = read_csv(FINAL_CSV)
    expansion_rows, expansion_payload = remote_expansion_rows()
    recovered_head, bank_equality = recover_nonconvex_d100_head()

    output = []
    for pde in PDE_ORDER:
        source_rows = expansion_rows if pde == "dense_fourier" else final_rows
        for dimension in DIMENSIONS:
            cell = [
                row
                for row in source_rows
                if row["pde"] == pde and int(row["dimension"]) == dimension
            ]
            by_method = {row["method"]: row for row in cell}
            expected = {"tangent_actor", "dpi", "pinn_hte", "sdgd_pinn"}
            if set(by_method) != expected:
                raise RuntimeError(
                    "Incomplete cell {} d={}: {}".format(
                        pde, dimension, sorted(by_method)
                    )
                )
            bank_hash, states = ensure_same_bank(cell, pde, dimension)
            ta = by_method["tangent_actor"]
            if pde == "dense_fourier":
                if ta["status"] != "COMPLETE":
                    raise RuntimeError("Incomplete dense Fourier TA row")
            elif ta["execution_status"] not in ("complete", "complete_reused"):
                raise RuntimeError("Incomplete final-suite TA row")

            head_source = "official actor_value_rmse field"
            head_value = ta.get("actor_value_rmse", "")
            if not head_value:
                if (pde, dimension) != ("nonconvex_coupled", 100):
                    raise RuntimeError(
                        "Unexpected missing value head for {} d={}".format(
                            pde, dimension
                        )
                    )
                head_value = recovered_head
                head_source = (
                    "selected-a4 post-freeze metric; canonical bank arrays byte-identical"
                )

            values = OrderedDict(
                [
                    (
                        "direct_value_head",
                        finite_float(head_value, "actor_value_rmse", pde, dimension),
                    ),
                    (
                        "dpi",
                        finite_float(
                            by_method["dpi"]["value_rmse"], "DPI RMSE", pde, dimension
                        ),
                    ),
                    (
                        "pinn_hte",
                        finite_float(
                            by_method["pinn_hte"]["value_rmse"],
                            "PINN-HTE RMSE",
                            pde,
                            dimension,
                        ),
                    ),
                    (
                        "sdgd_pinn",
                        finite_float(
                            by_method["sdgd_pinn"]["value_rmse"],
                            "SDGD-PINN RMSE",
                            pde,
                            dimension,
                        ),
                    ),
                ]
            )
            full_ta = finite_float(
                ta["value_rmse"], "Full TA FK-MCV RMSE", pde, dimension
            )
            head = values["direct_value_head"]
            wins = [method for method in METHODS if head < values[method]]
            rank = 1 + sum(values[method] < head for method in METHODS)
            best_baseline = min(values[method] for method in METHODS)
            output.append(
                OrderedDict(
                    [
                        ("pde", pde),
                        ("pde_label", PDE_LABELS[pde].replace("--", "–")),
                        ("dimension", dimension),
                        ("states", states),
                        ("test_bank_array_sha256", bank_hash),
                        ("full_ta_fk_mcv_rmse", full_ta),
                        ("direct_value_head_rmse", head),
                        ("dpi_rmse", values["dpi"]),
                        ("pinn_hte_rmse", values["pinn_hte"]),
                        ("sdgd_pinn_rmse", values["sdgd_pinn"]),
                        ("direct_head_rank_among_four", rank),
                        ("direct_head_pairwise_wins", len(wins)),
                        ("direct_head_beats_dpi", head < values["dpi"]),
                        ("direct_head_beats_pinn_hte", head < values["pinn_hte"]),
                        ("direct_head_beats_sdgd_pinn", head < values["sdgd_pinn"]),
                        ("direct_head_to_best_baseline_ratio", head / best_baseline),
                        ("head_metric_source", head_source),
                    ]
                )
            )

    provenance = {
        "final_suite_csv": {
            "path": str(FINAL_CSV),
            "sha256": sha256_file(FINAL_CSV),
            "git_blob": git_output("hash-object", str(FINAL_CSV)),
        },
        "highdim_expansion_csv": {
            "git_ref": EXPANSION_REF,
            "commit": git_output("rev-parse", EXPANSION_REF),
            "path": EXPANSION_CSV,
            "git_blob": git_output("rev-parse", EXPANSION_REF + ":" + EXPANSION_CSV),
            "sha256": hashlib.sha256(expansion_payload.encode("utf-8")).hexdigest(),
        },
        "nonconvex_d100_recovery": {
            "selection_json": str(SELECTION_JSON),
            "selection_json_sha256": sha256_file(SELECTION_JSON),
            "postfreeze_json": str(POSTFREEZE_JSON),
            "postfreeze_json_sha256": sha256_file(POSTFREEZE_JSON),
            "postfreeze_npz": str(POSTFREEZE_NPZ),
            "postfreeze_npz_sha256": sha256_file(POSTFREEZE_NPZ),
            "canonical_precision_npz": str(PRECISION_NPZ),
            "canonical_precision_npz_sha256": sha256_file(PRECISION_NPZ),
            "selected_actor": "a4",
            "bank_array_equality": bank_equality,
            "actor_value_rmse": recovered_head,
        },
    }
    return output, provenance


def summarize(rows):
    counts = Counter()
    per_pde = OrderedDict()
    for row in rows:
        counts["cells"] += 1
        wins = int(row["direct_head_pairwise_wins"])
        counts["pairwise_wins"] += wins
        counts["beats_any_baseline_cells"] += int(wins >= 1)
        counts["beats_majority_cells"] += int(wins >= 2)
        counts["best_among_four_cells"] += int(wins == 3)
        counts["beats_dpi"] += int(row["direct_head_beats_dpi"])
        counts["beats_pinn_hte"] += int(row["direct_head_beats_pinn_hte"])
        counts["beats_sdgd_pinn"] += int(row["direct_head_beats_sdgd_pinn"])
        counts["beats_full_ta"] += int(
            row["direct_value_head_rmse"] < row["full_ta_fk_mcv_rmse"]
        )
        bucket = per_pde.setdefault(
            row["pde"], {"label": row["pde_label"], "cells": 0, "wins": 0}
        )
        bucket["cells"] += 1
        bucket["wins"] += wins

    total_pairwise = counts["cells"] * len(METHODS)
    summary = OrderedDict(
        [
            ("pdes", len(PDE_ORDER)),
            ("dimensions", list(DIMENSIONS)),
            ("cells", counts["cells"]),
            ("pairwise_comparisons", total_pairwise),
            ("pairwise_wins", counts["pairwise_wins"]),
            ("pairwise_win_rate", counts["pairwise_wins"] / total_pairwise),
            (
                "beats_any_baseline_cells",
                counts["beats_any_baseline_cells"],
            ),
            (
                "beats_any_baseline_cell_rate",
                counts["beats_any_baseline_cells"] / counts["cells"],
            ),
            ("beats_majority_cells", counts["beats_majority_cells"]),
            (
                "beats_majority_cell_rate",
                counts["beats_majority_cells"] / counts["cells"],
            ),
            ("best_among_four_cells", counts["best_among_four_cells"]),
            ("beats_dpi", counts["beats_dpi"]),
            ("beats_pinn_hte", counts["beats_pinn_hte"]),
            ("beats_sdgd_pinn", counts["beats_sdgd_pinn"]),
            ("beats_full_ta_fk_mcv", counts["beats_full_ta"]),
            ("per_pde", per_pde),
            (
                "supported_claim",
                "The frozen direct value head is already a competitive solver: it "
                "beats at least one strong baseline in 24/25 cells and wins 47/75 "
                "pairwise comparisons.",
            ),
            (
                "unsupported_claim",
                "The data do not support saying that the direct value head "
                "outperforms almost all methods: it beats DPI in only 8/25 cells.",
            ),
        ]
    )
    return summary


def write_csv(rows):
    path = RESULTS / "VALUE_HEAD_FULL_DIMENSION.csv"
    with open(str(path), "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0].keys()),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    return path


def latex_cell(value, is_best):
    rendered = "$" + latex_sci(value) + "$"
    return r"\textbf{" + rendered + "}" if is_best else rendered


def build_table(rows):
    row_lookup = {(row["pde"], int(row["dimension"])): row for row in rows}
    lines = [
        r"\begin{table*}[t]",
        r"\centering",
        r"\caption{Whole-function value RMSE from the frozen direct Tangent-Actor value head, compared on the same state bank with three independently trained baselines. Lower is better; the smallest entry among the direct head and the three baselines is bold. The final FK--MCV estimator is intentionally excluded from the ranking because this table asks how accurate the learned head is by itself.}",
        r"\label{tab:value-head-full-dimension}",
        r"\small",
        r"\setlength{\tabcolsep}{5.2pt}",
        r"\begin{tabular}{llccccc}",
        r"\toprule",
        r"PDE & $d$ & Direct value head & DPI & PINN-HTE & SDGD-PINN & Head wins \\",
        r"\midrule",
    ]
    for pde_index, pde in enumerate(PDE_ORDER):
        for dimension_index, dimension in enumerate(DIMENSIONS):
            row = row_lookup[(pde, dimension)]
            values = OrderedDict(
                [
                    ("direct_value_head", row["direct_value_head_rmse"]),
                    ("dpi", row["dpi_rmse"]),
                    ("pinn_hte", row["pinn_hte_rmse"]),
                    ("sdgd_pinn", row["sdgd_pinn_rmse"]),
                ]
            )
            best = min(values.values())
            pde_label = PDE_LABELS[pde] if dimension_index == 0 else ""
            fields = [
                pde_label,
                str(dimension),
                latex_cell(values["direct_value_head"], values["direct_value_head"] == best),
                latex_cell(values["dpi"], values["dpi"] == best),
                latex_cell(values["pinn_hte"], values["pinn_hte"] == best),
                latex_cell(values["sdgd_pinn"], values["sdgd_pinn"] == best),
                "{}/3".format(row["direct_head_pairwise_wins"]),
            ]
            lines.append(" & ".join(fields) + r" \\")
        if pde_index != len(PDE_ORDER) - 1:
            lines.append(r"\addlinespace[2pt]")
    lines.extend(
        [
            r"\midrule",
            r"\multicolumn{2}{l}{Pairwise wins over all 25 cells} & \multicolumn{4}{c}{47/75 (62.7\%)} & \\",
            r"\multicolumn{2}{l}{Wins by baseline} & \multicolumn{4}{c}{DPI: 8/25; PINN-HTE: 15/25; SDGD-PINN: 24/25} & \\",
            r"\bottomrule",
            r"\end{tabular}",
            r"\end{table*}",
            "",
        ]
    )
    path = TABLES / "value_head_full_dimension.tex"
    path.write_text("\n".join(lines))
    return path


def build_preview(table_path):
    table = table_path.read_text()
    document = "\n".join(
        [
            r"\documentclass[10pt]{article}",
            r"\usepackage[margin=0.42in]{geometry}",
            r"\usepackage{booktabs}",
            r"\usepackage[T1]{fontenc}",
            r"\usepackage{lmodern}",
            r"\usepackage{microtype}",
            r"\pagestyle{empty}",
            r"\begin{document}",
            table,
            r"\end{document}",
            "",
        ]
    )
    path = TABLES / "value_head_full_dimension_preview.tex"
    path.write_text(document)
    return path


def manuscript_method_values(row):
    return OrderedDict(
        [
            ("direct_value_head", row["direct_value_head_rmse"]),
            ("dpi", row["dpi_rmse"]),
            ("pinn_hte", row["pinn_hte_rmse"]),
            ("sdgd_pinn", row["sdgd_pinn_rmse"]),
        ]
    )


def manuscript_number(value, is_best):
    rendered = fmt_sci(value)
    return r"\textbf{" + rendered + "}" if is_best else rendered


def manuscript_method_row(method, pdes, lookup):
    label = METHOD_LABELS[method]
    fields = [label]
    for pde in pdes:
        for dimension in DIMENSIONS:
            row = lookup[(pde, dimension)]
            values = manuscript_method_values(row)
            fields.append(
                manuscript_number(values[method], values[method] == min(values.values()))
            )
    prefix = r"\rowcolor{TAblue}" + "\n" if method == "direct_value_head" else ""
    return prefix + " & ".join(fields) + r" \\"


def build_paired_manuscript_panel(pde_left, pde_right, lookup):
    lines = [
        r"\toprule",
        (
            r"\textbf{Method} & \multicolumn{5}{c!{\color{PDEdivider}\vrule width 0.75pt}}"
            r"{\textbf{" + PDE_LABELS[pde_left] + r"}} & \multicolumn{5}{c}{\textbf{"
            + PDE_LABELS[pde_right]
            + r"}} \\"
        ),
        r"\cmidrule(lr){2-6}\cmidrule(lr){7-11}",
        r"\textbf{Dimension} & $20$ & $50$ & $100$ & $200$ & $400$ & $20$ & $50$ & $100$ & $200$ & $400$ \\",
        r"\midrule",
    ]
    methods = tuple(METHOD_LABELS)
    for index, method in enumerate(methods):
        lines.append(manuscript_method_row(method, (pde_left, pde_right), lookup))
        if index != len(methods) - 1:
            lines.append(r"\cmidrule(lr){1-11}")
    lines.append(r"\bottomrule")
    return lines


def build_single_manuscript_panel(pde, lookup):
    grouped_dimension_header = []
    for dimension_index, dimension in enumerate(DIMENSIONS):
        suffix = (
            r"!{\color{CellDivider}\vrule width 0.35pt}"
            if dimension_index != len(DIMENSIONS) - 1
            else ""
        )
        grouped_dimension_header.append(
            r"\multicolumn{2}{c" + suffix + "}{$" + str(dimension) + "$}"
        )

    lines = [
        r"\toprule",
        r"\textbf{Method} & \multicolumn{10}{c}{\textbf{" + PDE_LABELS[pde] + r"}} \\",
        r"\cmidrule(lr){2-11}",
        r"\textbf{Dimension} & " + " & ".join(grouped_dimension_header) + r" \\",
        r"\midrule",
    ]
    methods = tuple(METHOD_LABELS)
    for method_index, method in enumerate(methods):
        cells = []
        for dimension_index, dimension in enumerate(DIMENSIONS):
            row = lookup[(pde, dimension)]
            values = manuscript_method_values(row)
            suffix = (
                r"!{\color{CellDivider}\vrule width 0.35pt}"
                if dimension_index != len(DIMENSIONS) - 1
                else ""
            )
            value = manuscript_number(
                values[method], values[method] == min(values.values())
            )
            cells.append(r"\multicolumn{2}{c" + suffix + "}{" + value + "}")
        prefix = r"\rowcolor{TAblue}" + "\n" if method == "direct_value_head" else ""
        lines.append(prefix + METHOD_LABELS[method] + " & " + " & ".join(cells) + r" \\")
        if method_index != len(methods) - 1:
            lines.append(r"\cmidrule(lr){1-11}")
    lines.append(r"\bottomrule")
    return lines


def build_manuscript_table(rows):
    lookup = {(row["pde"], int(row["dimension"])): row for row in rows}
    column_spec = (
        r"@{}l!{\color{CellDivider}\vrule width 0.35pt}"
        r"c!{\color{CellDivider}\vrule width 0.35pt}"
        r"c!{\color{CellDivider}\vrule width 0.35pt}"
        r"c!{\color{CellDivider}\vrule width 0.35pt}"
        r"c!{\color{CellDivider}\vrule width 0.35pt}"
        r"c!{\hspace{5pt}\color{PDEdivider}\vrule width 0.75pt\hspace{5pt}}"
        r"c!{\color{CellDivider}\vrule width 0.35pt}"
        r"c!{\color{CellDivider}\vrule width 0.35pt}"
        r"c!{\color{CellDivider}\vrule width 0.35pt}"
        r"c!{\color{CellDivider}\vrule width 0.35pt}c@{}"
    )
    lines = [
        r"\begin{table}[!ht]",
        r"\centering",
        r"\caption{Whole-function value RMSE of the frozen direct Tangent-Actor value head across five PDEs and five dimensions (100 method cells). Comparisons within each PDE--dimension cell use the same frozen test bank; best values are bolded.}",
        r"\label{tab:direct_value_head_five_pdes}",
        r"\vspace{2pt}",
        r"{\fontsize{7.3}{8.6}\selectfont",
        r"\setlength{\tabcolsep}{1.1pt}",
        r"\renewcommand{\arraystretch}{1.10}",
        r"\resizebox{\columnwidth}{!}{%",
        r"\begin{tabular}{" + column_spec + "}",
    ]
    lines.extend(build_paired_manuscript_panel("allen_cahn", "joint_concave", lookup))
    lines.append(r"\addlinespace[7pt]")
    lines.extend(
        build_paired_manuscript_panel(
            "joint_convex", "nonconvex_coupled", lookup
        )
    )
    lines.append(r"\addlinespace[7pt]")
    lines.extend(build_single_manuscript_panel("dense_fourier", lookup))
    lines.extend(
        [
            r"\end{tabular}",
            r"}%",
            r"}",
            r"\end{table}",
            "",
        ]
    )
    path = TABLES / "value_head_full_dimension_icml.tex"
    path.write_text("\n".join(lines))
    return path


def build_report(rows, summary, provenance):
    lines = [
        "# Direct value-head full-dimensional audit",
        "",
        "This is a frozen-artifact analysis: no model was retrained, no checkpoint "
        "was reselected, and no test value entered training or selection. The scope "
        "is the five PDEs whose audited Tangent-Actor architecture has a value "
        "output, at dimensions 20, 50, 100, 200, and 400. Gradient-only HJB Actors "
        "are excluded because they do not have a value head.",
        "",
        "## Result",
        "",
        "The direct head beats at least one independently trained baseline in "
        "24/25 cells (96%), a majority of the three baselines in 17/25 cells "
        "(68%), and all three baselines in 6/25 cells (24%). Across all matched "
        "pairwise comparisons it wins 47/75 (62.7%): 8/25 against DPI, 15/25 "
        "against PINN-HTE, and 24/25 against SDGD-PINN.",
        "",
        "The strongest defensible wording is therefore that the auxiliary value "
        "head is already a competitive whole-function solver without FK inference. "
        "The stronger claim that it outperforms almost all methods is not supported, "
        "primarily because DPI is better in 17/25 cells. The formal Tangent-Actor "
        "FK--MCV estimator remains better than the direct head in all 25 cells.",
        "",
        "## Scope and comparability",
        "",
        "- Every comparison within a PDE/dimension cell uses the exact same frozen "
        "  state bank and numerical reference.",
        "- The analysis is seed 0, matching the official single-seed benchmark.",
        "- The 100D Oscillatory value-head metric is recovered from the selected-a4 "
        "  post-freeze record; its test time, states, and truth arrays are byte-identical "
        "  to the canonical 512-state benchmark artifact.",
        "- Fourier--Tanh rows come from the audited ICE expansion result at Git "
        "  commit `{}`.".format(provenance["highdim_expansion_csv"]["commit"]),
        "",
        "## Files",
        "",
        "- `VALUE_HEAD_FULL_DIMENSION.csv`: raw 25-cell comparison",
        "- `VALUE_HEAD_SUMMARY.json`: aggregate win counts",
        "- `SOURCE_AUDIT.json`: hashes and source provenance",
        "- `tables/value_head_full_dimension_icml.tex`: compact ICML single-column table",
        "- `tables/value_head_full_dimension.tex`: publication table",
        "- `tables/value_head_full_dimension_preview.pdf`: rendered preview",
        "",
    ]
    path = RESULTS / "VALUE_HEAD_FULL_DIMENSION_REPORT.md"
    path.write_text("\n".join(lines))
    return path


def main():
    RESULTS.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    rows, provenance = build_rows()
    if len(rows) != len(PDE_ORDER) * len(DIMENSIONS):
        raise RuntimeError("Expected 25 complete cells")
    summary = summarize(rows)
    csv_path = write_csv(rows)
    summary_path = RESULTS / "VALUE_HEAD_SUMMARY.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=False) + "\n")
    audit_path = RESULTS / "SOURCE_AUDIT.json"
    audit_path.write_text(json.dumps(provenance, indent=2, sort_keys=False) + "\n")
    table_path = build_table(rows)
    manuscript_table_path = build_manuscript_table(rows)
    preview_path = build_preview(table_path)
    report_path = build_report(rows, summary, provenance)
    print("wrote", csv_path)
    print("wrote", summary_path)
    print("wrote", audit_path)
    print("wrote", table_path)
    print("wrote", manuscript_table_path)
    print("wrote", preview_path)
    print("wrote", report_path)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
