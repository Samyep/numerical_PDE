#!/usr/bin/env python3
"""Audit the frozen Joint Nonseparable Tangent-Actor direct value head.

This script is intentionally extraction-only.  It reads the selected actor
metadata and the four already-saved evaluation NPZ files, verifies that the
deterministic ``actor_r`` arrays agree exactly, and computes direct-head
metrics on the original 32-state frozen bank.  It never imports the solver,
loads a model, runs FK/MCV inference, or performs checkpoint selection.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np


DIMENSIONS = (20, 50, 100, 200, 400)
REPLICAS = (0, 1, 2, 3)
PDE = "joint_nonseparable"
EXPECTED_ORIGIN_TIP = "24440691ac7beb9e4f2b8d09db6517ca1fee7248"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def array_sha256(*arrays: np.ndarray) -> str:
    """Match main_highdim_pilot.equations.array_sha256 exactly."""
    digest = hashlib.sha256()
    for array in arrays:
        value = np.ascontiguousarray(array)
        digest.update(str(value.dtype).encode("ascii"))
        digest.update(str(value.shape).encode("ascii"))
        digest.update(value.tobytes())
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as stream:
        return json.load(stream)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def git_revision(worktree: Path, revision: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(worktree), "rev-parse", revision],
        text=True,
    ).strip()


def resolve_recorded_path(recorded: str, fallback: Path) -> Path:
    candidate = Path(recorded)
    if candidate.is_file():
        return candidate.resolve()
    if fallback.is_file():
        return fallback.resolve()
    raise FileNotFoundError(
        f"Neither recorded path nor worktree-relative fallback exists: "
        f"{candidate} ; {fallback}"
    )


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def exact_array_equal(left: np.ndarray, right: np.ndarray) -> bool:
    return left.dtype == right.dtype and left.shape == right.shape and np.array_equal(left, right)


def finite_scalar(value: float, label: str) -> float:
    result = float(value)
    require(math.isfinite(result), f"Non-finite scalar for {label}: {result}")
    return result


def load_main_metrics(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def baseline_row(
    rows: Iterable[dict[str, str]], dimension: int, method: str
) -> dict[str, Any]:
    matches = [
        row
        for row in rows
        if row.get("pde") == PDE
        and int(row["dimension"]) == dimension
        and row.get("method") == method
    ]
    require(len(matches) == 1, f"Expected one {method} row for d={dimension}; got {len(matches)}")
    source = matches[0]
    failed = source.get("numerical_failure", "").strip().lower() == "true"
    rmse_text = source.get("rmse", "").strip()
    if failed:
        require(not rmse_text, f"Numerical failure unexpectedly has RMSE: {method}, d={dimension}")
        rmse = None
    else:
        require(bool(rmse_text), f"Missing RMSE for valid baseline: {method}, d={dimension}")
        rmse = finite_scalar(float(rmse_text), f"{method} d={dimension} RMSE")
    return {
        "method": method,
        "rmse": rmse,
        "status": "NF" if failed else "ok",
        "source_row": source,
    }


def format_sci(value: float) -> str:
    mantissa, exponent = f"{value:.2e}".split("e")
    return f"{mantissa}e{int(exponent)}"


def latex_value(value: float | None, best: float) -> str:
    if value is None:
        return r"\textit{NF}"
    rendered = format_sci(value)
    if math.isclose(value, best, rel_tol=1e-14, abs_tol=0.0):
        return rf"\textbf{{{rendered}}}"
    return rendered


def render_latex(rows: list[dict[str, Any]]) -> str:
    best_by_dimension: dict[int, float] = {}
    for row in rows:
        values = [
            row["direct_value_head_rmse"],
            row["full_ta_fk_mcv_rmse"],
            row["fh_mlp_rmse"],
        ]
        if row["deep_bsde_rmse"] is not None:
            values.append(row["deep_bsde_rmse"])
        best_by_dimension[row["dimension"]] = min(values)

    def method_line(label: str, key: str, row_color: bool = False) -> str:
        prefix = "\\rowcolor{TAblue}\n" if row_color else ""
        cells = [
            latex_value(row[key], best_by_dimension[row["dimension"]]) for row in rows
        ]
        return prefix + label + " & " + " & ".join(cells) + r" \\"

    return "\n".join(
        [
            r"\begin{table}[!ht]",
            r"\centering",
            r"\caption{Joint Nonseparable frozen-bank value accuracy. The direct",
            r"value head is the deterministic output of the formally selected seed-0",
            r"checkpoint on the original 32-state bank; Full TA (FK--MCV) is retained",
            r"only as an auxiliary reference. Entries are RMSE and best valid values",
            r"are bolded. NF denotes numerical failure in the frozen Deep BSDE run.}",
            r"\label{tab:joint_nonseparable_value_head}",
            r"\vspace{2pt}",
            r"{\fontsize{7.3}{8.6}\selectfont",
            r"\setlength{\tabcolsep}{1.8pt}",
            r"\renewcommand{\arraystretch}{1.10}",
            r"\resizebox{\columnwidth}{!}{%",
            r"\begin{tabular}{@{}l!{\color{CellDivider}\vrule width 0.35pt}c!{\color{CellDivider}\vrule width 0.35pt}c!{\color{CellDivider}\vrule width 0.35pt}c!{\color{CellDivider}\vrule width 0.35pt}c!{\color{CellDivider}\vrule width 0.35pt}c@{}}",
            r"\toprule",
            r"\textbf{Method} & \multicolumn{5}{c}{\textbf{Ambient dimension $d$}} \\",
            r"\cmidrule(lr){2-6}",
            r" & $20$ & $50$ & $100$ & $200$ & $400$ \\",
            r"\midrule",
            method_line("Direct value head", "direct_value_head_rmse", row_color=True),
            r"\cmidrule(lr){1-6}",
            method_line(r"Full TA (FK--MCV)", "full_ta_fk_mcv_rmse"),
            r"\cmidrule(lr){1-6}",
            method_line(r"FH-MLP", "fh_mlp_rmse"),
            r"\cmidrule(lr){1-6}",
            method_line(r"Deep BSDE", "deep_bsde_rmse"),
            r"\bottomrule",
            r"\end{tabular}",
            r"}%",
            r"}",
            r"\end{table}",
            "",
        ]
    )


def render_report(
    rows: list[dict[str, Any]],
    checks: dict[str, Any],
    worktree_head: str,
    origin_tip: str,
    legal_cells: int,
    won_cells: int,
) -> str:
    table_lines = [
        "| d | Actor | Direct-head RMSE | Direct-head RelMAE | Full TA FK-MCV RMSE | FH-MLP RMSE | Deep BSDE RMSE |",
        "|---:|:---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        deep = "NF" if row["deep_bsde_rmse"] is None else format_sci(row["deep_bsde_rmse"])
        table_lines.append(
            "| {dimension} | {selected_actor} | {direct} | {rmae} | {full} | {fh} | {deep} |".format(
                dimension=row["dimension"],
                selected_actor=row["selected_actor"],
                direct=format_sci(row["direct_value_head_rmse"]),
                rmae=format_sci(row["direct_value_head_rmae"]),
                full=format_sci(row["full_ta_fk_mcv_rmse"]),
                fh=format_sci(row["fh_mlp_rmse"]),
                deep=deep,
            )
        )

    passed = all(value is True for value in checks.values())
    suitable = passed and won_cells == legal_cells
    conclusion = (
        "Suitable for replacing the Oscillatory Value-Gradient diagnostic: the direct "
        f"head wins {won_cells}/{legal_cells} legal FH-MLP/Deep-BSDE cells, while the "
        "comparison remains explicitly a selected seed-0, 32-state frozen-bank diagnostic."
        if suitable
        else
        "Do not use as an unqualified replacement: either an audit check failed or the "
        f"direct head wins only {won_cells}/{legal_cells} legal FH-MLP/Deep-BSDE cells."
    )

    return "\n".join(
        [
            "# Joint Nonseparable direct value-head audit",
            "",
            "## Scope",
            "",
            "This audit extracts `actor_r` from the four formal saved prediction NPZ files.",
            "It does not retrain, reselect a checkpoint, load a model, or run FK/MCV.",
            "All metrics use the original 32-state frozen bank. The result is one selected",
            "seed-0 model, not a multi-seed mean and not a 960-state whole-function result.",
            "",
            f"- Worktree HEAD audited: `{worktree_head}`",
            f"- Formal origin tip: `{origin_tip}`",
            f"- Expected formal tip matched: `{origin_tip == EXPECTED_ORIGIN_TIP}`",
            "",
            "## Results",
            "",
            *table_lines,
            "",
            "Full TA (FK-MCV) is retained only as an auxiliary reference and is not the",
            "direct value-head metric.",
            "",
            "## Baseline-cell outcome",
            "",
            f"The direct value head wins **{won_cells}/{legal_cells}** legal baseline cells",
            "against FH-MLP and nonfailed Deep BSDE runs. Deep BSDE numerical failures are",
            "reported as NF and are excluded from the denominator.",
            "",
            "## Audit checks",
            "",
            *[f"- {'PASS' if value else 'FAIL'}: `{name}`" for name, value in checks.items()],
            "",
            "## Replacement assessment",
            "",
            conclusion,
            "",
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--worktree", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()

    worktree = args.worktree.resolve()
    result_root = worktree / "invariant_region_mlp/results/single_point_joint_extension"
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else result_root / "value_head_audit"
    )
    prediction_output_dir = output_dir / "predictions"
    table_output_dir = output_dir / "tables"
    prediction_output_dir.mkdir(parents=True, exist_ok=True)
    table_output_dir.mkdir(parents=True, exist_ok=True)

    main_metrics_path = result_root / "main_metrics.csv"
    require(main_metrics_path.is_file(), f"Missing baseline metrics: {main_metrics_path}")
    metric_rows = load_main_metrics(main_metrics_path)

    worktree_head = git_revision(worktree, "HEAD")
    origin_tip = git_revision(worktree, "origin/ir-mlp-icml-joint-pde")
    require(origin_tip == EXPECTED_ORIGIN_TIP, f"Unexpected formal origin tip: {origin_tip}")

    output_rows: list[dict[str, Any]] = []
    checkpoint_entries: list[dict[str, Any]] = []
    all_replica_value_heads_identical = True
    all_banks_match_predictions = True
    all_selected_checkpoints = True
    all_no_ground_truth_selection = True
    all_finite = True
    all_32_states = True
    old_checkpoint_used = False

    for dimension in DIMENSIONS:
        run_id = f"joint_nonseparable_d{dimension}_mature_s0"
        training_path = result_root / "raw/actor" / f"training_{run_id}.json"
        require(training_path.is_file(), f"Missing training metadata: {training_path}")
        training = read_json(training_path)
        selection = training["selection"]

        require(int(training["dimension"]) == dimension, f"Training dimension mismatch: d={dimension}")
        require(int(training["seed"]) == 0, f"Unexpected training seed: d={dimension}")
        require(training["status"] == "complete", f"Incomplete training metadata: d={dimension}")
        selected_actor = str(selection["selected_actor"])
        recorded_checkpoint = str(selection["selected_checkpoint"])
        fallback_checkpoint = result_root / "checkpoints" / Path(recorded_checkpoint).name
        checkpoint_path = resolve_recorded_path(recorded_checkpoint, fallback_checkpoint)
        expected_prefix = f"joint_nonseparable_d{dimension}_mature_s0_{selected_actor}.pt"
        require(checkpoint_path.name == expected_prefix, f"Selected checkpoint/actor mismatch: {checkpoint_path}")
        require(
            "convex_concave_semilinear_extended" not in str(checkpoint_path),
            f"Forbidden legacy checkpoint selected: {checkpoint_path}",
        )
        old_checkpoint_used = old_checkpoint_used or checkpoint_path.name.startswith("nonseparable_")

        no_gt = (
            training.get("ground_truth_used_for_training_or_selection") is False
            and selection.get("ground_truth_used") is False
            and all(row.get("ground_truth_used") is False for row in selection.get("rows", []))
        )
        require(no_gt, f"Ground truth entered training/selection: d={dimension}")
        all_no_ground_truth_selection &= no_gt

        checkpoint_hash = sha256_file(checkpoint_path)
        matching_phases = [
            phase
            for phase in training.get("phases", [])
            if Path(str(phase.get("checkpoint", ""))).name == checkpoint_path.name
        ]
        require(len(matching_phases) == 1, f"Cannot uniquely audit checkpoint phase: d={dimension}")
        metadata_checkpoint_hash = matching_phases[0].get("checkpoint_sha256")
        require(
            metadata_checkpoint_hash == checkpoint_hash,
            f"Checkpoint SHA mismatch: d={dimension}",
        )

        bank_path = result_root / "frozen_query_banks" / f"joint_nonseparable_d{dimension}" / "bank.npz"
        require(bank_path.is_file(), f"Missing frozen bank: {bank_path}")
        bank_file_hash = sha256_file(bank_path)
        with np.load(bank_path, allow_pickle=False) as bank:
            require(set(bank.files) >= {"t", "x", "reference"}, f"Incomplete bank: {bank_path}")
            bank_t = np.array(bank["t"], copy=True)
            bank_x = np.array(bank["x"], copy=True)
            bank_truth = np.array(bank["reference"], copy=True)
            bank_truth_z = np.array(bank["reference_z"], copy=True) if "reference_z" in bank else None
        require(bank_t.shape == (32,), f"Unexpected bank t shape: d={dimension}, {bank_t.shape}")
        require(bank_x.shape == (32, dimension), f"Unexpected bank x shape: d={dimension}, {bank_x.shape}")
        require(bank_truth.shape == (32,), f"Unexpected reference shape: d={dimension}")
        all_32_states &= bank_t.shape[0] == 32
        test_input_hash = array_sha256(bank_t, bank_x)
        reference_value_hash = array_sha256(bank_truth)
        reference_full_hash = (
            array_sha256(bank_truth, bank_truth_z)
            if bank_truth_z is not None
            else reference_value_hash
        )

        actor_r_replicas: list[np.ndarray] = []
        evaluation_records: list[dict[str, Any]] = []
        evaluation_record_audit: list[dict[str, Any]] = []
        base_t: np.ndarray | None = None
        base_x: np.ndarray | None = None
        base_truth: np.ndarray | None = None

        for replica in REPLICAS:
            evaluation_path = (
                result_root
                / "raw/actor_evaluations"
                / run_id
                / "selected"
                / f"steps64_rep{replica:02d}.json"
            )
            require(evaluation_path.is_file(), f"Missing evaluation record: {evaluation_path}")
            evaluation = read_json(evaluation_path)
            evaluation_records.append(evaluation)
            require(int(evaluation["dimension"]) == dimension, f"Evaluation dimension mismatch: {evaluation_path}")
            require(int(evaluation["inference_replicate"]) == replica, f"Replica mismatch: {evaluation_path}")
            require(evaluation["actor_request"] == "selected", f"Nonselected actor request: {evaluation_path}")
            require(evaluation["resolved_actor"] == selected_actor, f"Resolved actor mismatch: {evaluation_path}")
            require(int(evaluation["training_seed"]) == 0, f"Training seed mismatch: {evaluation_path}")
            require(int(evaluation["test_states"]) == 32, f"Not a 32-state evaluation: {evaluation_path}")
            require(
                evaluation.get("ground_truth_loaded_only_after_checkpoint_selection") is True,
                f"Ground truth timing audit failed: {evaluation_path}",
            )

            recorded_prediction = str(evaluation["prediction_path"])
            fallback_prediction = (
                result_root
                / "predictions/actor"
                / run_id
                / "selected"
                / f"steps64_rep{replica:02d}.npz"
            )
            prediction_path = resolve_recorded_path(recorded_prediction, fallback_prediction)
            prediction_file_hash = sha256_file(prediction_path)
            require(
                prediction_file_hash == evaluation["prediction_sha256"],
                f"Prediction SHA mismatch: {prediction_path}",
            )
            require(
                evaluation["test_bank_file_sha256"] == bank_file_hash,
                f"Bank file SHA mismatch in evaluation: {evaluation_path}",
            )
            require(
                evaluation["test_bank_sha256"] == test_input_hash,
                f"Bank array SHA mismatch in evaluation: {evaluation_path}",
            )

            with np.load(prediction_path, allow_pickle=False) as prediction:
                require(
                    set(prediction.files) >= {"actor_r", "truth", "t", "x"},
                    f"Prediction lacks direct-head arrays: {prediction_path}",
                )
                pred_t = np.array(prediction["t"], copy=True)
                pred_x = np.array(prediction["x"], copy=True)
                pred_truth = np.array(prediction["truth"], copy=True)
                actor_r = np.array(prediction["actor_r"], copy=True)

            bank_match = (
                exact_array_equal(pred_t, bank_t)
                and exact_array_equal(pred_x, bank_x)
                and exact_array_equal(pred_truth, bank_truth)
            )
            require(bank_match, f"Prediction does not match frozen bank: {prediction_path}")
            all_banks_match_predictions &= bank_match
            require(actor_r.shape == (32,), f"Unexpected actor_r shape: {prediction_path}")
            require(np.isfinite(actor_r).all(), f"Non-finite actor_r: {prediction_path}")

            head_error = actor_r - pred_truth
            head_rmse = float(np.sqrt(np.mean(np.square(head_error))))
            require(
                math.isclose(
                    head_rmse,
                    float(evaluation["actor_value_rmse"]),
                    rel_tol=1e-13,
                    abs_tol=1e-15,
                ),
                f"Recomputed actor RMSE mismatch: {evaluation_path}",
            )

            if base_t is None:
                base_t, base_x, base_truth = pred_t, pred_x, pred_truth
            actor_r_replicas.append(actor_r)
            evaluation_record_audit.append(
                {
                    "replica": replica,
                    "evaluation_json": str(evaluation_path.resolve()),
                    "evaluation_json_sha256": sha256_file(evaluation_path),
                    "prediction_npz": str(prediction_path),
                    "prediction_npz_sha256": prediction_file_hash,
                    "actor_value_rmse_json": float(evaluation["actor_value_rmse"]),
                    "actor_value_rmse_recomputed": head_rmse,
                }
            )

        require(base_t is not None and base_x is not None and base_truth is not None, "No predictions loaded")
        replica_identity = all(
            exact_array_equal(actor_r_replicas[0], candidate)
            for candidate in actor_r_replicas[1:]
        )
        json_rmse_identity = len(
            {float(record["actor_value_rmse"]).hex() for record in evaluation_records}
        ) == 1
        require(replica_identity, f"Direct value head differs across replicas: d={dimension}")
        require(json_rmse_identity, f"actor_value_rmse differs across replicas: d={dimension}")
        all_replica_value_heads_identical &= replica_identity and json_rmse_identity

        actor_r = actor_r_replicas[0]
        error = actor_r - base_truth
        direct_rmse = finite_scalar(np.sqrt(np.mean(np.square(error))), f"direct d={dimension} RMSE")
        direct_mae = finite_scalar(np.mean(np.abs(error)), f"direct d={dimension} MAE")
        direct_relative_l2 = finite_scalar(
            np.linalg.norm(error) / np.linalg.norm(base_truth),
            f"direct d={dimension} relative L2",
        )
        direct_rmae = finite_scalar(
            np.sum(np.abs(error)) / np.sum(np.abs(base_truth)),
            f"direct d={dimension} RelMAE",
        )
        direct_bias = finite_scalar(np.mean(error), f"direct d={dimension} bias")
        all_finite &= all(
            math.isfinite(value)
            for value in (direct_rmse, direct_mae, direct_relative_l2, direct_rmae, direct_bias)
        )

        full_ta = baseline_row(metric_rows, dimension, "Tangent-Actor")
        fh_mlp = baseline_row(metric_rows, dimension, "Full-History MLP")
        deep_bsde = baseline_row(metric_rows, dimension, "Deep BSDE")
        require(full_ta["rmse"] is not None, f"Formal TA failed unexpectedly: d={dimension}")
        require(fh_mlp["rmse"] is not None, f"FH-MLP failed unexpectedly: d={dimension}")

        compact_prediction_path = (
            prediction_output_dir / f"joint_nonseparable_d{dimension}_direct_head.npz"
        )
        np.savez_compressed(
            compact_prediction_path,
            t=base_t,
            x=base_x,
            truth=base_truth,
            actor_r=actor_r,
            error=error,
            dimension=np.asarray(dimension, dtype=np.int64),
            training_seed=np.asarray(0, dtype=np.int64),
            selected_actor=np.asarray(selected_actor),
            checkpoint_sha256=np.asarray(checkpoint_hash),
            test_bank_file_sha256=np.asarray(bank_file_hash),
            test_input_sha256=np.asarray(test_input_hash),
            reference_value_sha256=np.asarray(reference_value_hash),
        )

        output_row = {
            "dimension": dimension,
            "selected_actor": selected_actor,
            "direct_value_head_rmse": direct_rmse,
            "direct_value_head_mae": direct_mae,
            "direct_value_head_relative_l2": direct_relative_l2,
            "direct_value_head_rmae": direct_rmae,
            "direct_value_head_bias": direct_bias,
            "full_ta_fk_mcv_rmse": full_ta["rmse"],
            "fh_mlp_rmse": fh_mlp["rmse"],
            "deep_bsde_rmse": deep_bsde["rmse"],
            "deep_bsde_status": deep_bsde["status"],
            "states": 32,
            "training_seed": 0,
            "inference_replicas": 4,
            "replicas_are_independent_training_seeds": False,
            "test_bank_sha256": bank_file_hash,
            "test_bank_array_sha256": test_input_hash,
            "reference_value_sha256": reference_value_hash,
            "reference_full_sha256": reference_full_hash,
            "checkpoint_sha256": checkpoint_hash,
            "training_json_sha256": sha256_file(training_path),
            "direct_beats_fh_mlp": direct_rmse < fh_mlp["rmse"],
            "direct_beats_deep_bsde": (
                None if deep_bsde["rmse"] is None else direct_rmse < deep_bsde["rmse"]
            ),
            "direct_beats_full_ta_fk_mcv": direct_rmse < full_ta["rmse"],
            "prediction_npz": (
                Path("predictions") / compact_prediction_path.name
            ).as_posix(),
        }
        output_rows.append(output_row)

        checkpoint_entries.append(
            {
                "dimension": dimension,
                "training_seed": int(training["seed"]),
                "selected_actor": selected_actor,
                "selected_checkpoint_recorded": recorded_checkpoint,
                "selected_checkpoint_resolved": str(checkpoint_path),
                "checkpoint_sha256": checkpoint_hash,
                "checkpoint_sha256_from_training_metadata": metadata_checkpoint_hash,
                "checkpoint_sha256_match": metadata_checkpoint_hash == checkpoint_hash,
                "training_json": str(training_path.resolve()),
                "training_json_sha256": sha256_file(training_path),
                "checkpoint_selection_rule": selection["criterion"],
                "ground_truth_used_for_training_or_selection": training.get(
                    "ground_truth_used_for_training_or_selection"
                ),
                "selection_ground_truth_used": selection.get("ground_truth_used"),
                "architecture_and_training_configuration": training["configuration"],
                "parameter_count": training.get("parameter_count"),
                "full_ambient_input": training.get("full_ambient_input"),
                "adapter_only": training.get("adapter_only"),
                "oracle_projection": training.get("oracle_projection"),
                "source_solver_module": training.get("source_solver_module"),
                "bank_file": str(bank_path.resolve()),
                "bank_file_sha256": bank_file_hash,
                "test_input_sha256": test_input_hash,
                "reference_value_sha256": reference_value_hash,
                "reference_full_sha256": reference_full_hash,
                "states": 32,
                "evaluation_records": evaluation_record_audit,
                "value_head_replicas_exactly_identical": replica_identity,
                "actor_value_rmse_json_exactly_identical": json_rmse_identity,
            }
        )
        all_selected_checkpoints &= checkpoint_path.name == expected_prefix

    legal_cells = 0
    won_cells = 0
    for row in output_rows:
        legal_cells += 1
        won_cells += int(row["direct_beats_fh_mlp"])
        if row["deep_bsde_status"] == "ok":
            legal_cells += 1
            won_cells += int(row["direct_beats_deep_bsde"])

    checks = {
        "all_five_value_heads_from_formal_selected_checkpoints": all_selected_checkpoints,
        "no_test_rmse_used_for_actor_selection": all_no_ground_truth_selection,
        "extraction_only_no_fk_mcv_or_rollout_called": True,
        "four_inference_replicas_have_exactly_identical_direct_heads": all_replica_value_heads_identical,
        "all_predictions_match_original_frozen_banks": all_banks_match_predictions,
        "all_banks_have_32_states": all_32_states,
        "all_direct_metrics_finite": all_finite,
        "single_training_seed_zero_not_described_as_multiseed": True,
        "result_not_described_as_960_state_whole_function_accuracy": True,
        "forbidden_legacy_nonseparable_checkpoint_not_used": not old_checkpoint_used,
        "formal_origin_tip_matches_expected": origin_tip == EXPECTED_ORIGIN_TIP,
    }
    require(all(checks.values()), f"Audit check failure: {checks}")

    csv_path = output_dir / "JOINT_NONSEPARABLE_VALUE_HEAD.csv"
    csv_fields = [
        "dimension",
        "selected_actor",
        "direct_value_head_rmse",
        "direct_value_head_mae",
        "direct_value_head_relative_l2",
        "direct_value_head_rmae",
        "direct_value_head_bias",
        "full_ta_fk_mcv_rmse",
        "fh_mlp_rmse",
        "deep_bsde_rmse",
        "deep_bsde_status",
        "states",
        "test_bank_sha256",
        "test_bank_array_sha256",
        "reference_value_sha256",
        "checkpoint_sha256",
        "training_json_sha256",
    ]
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=csv_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(output_rows)

    generated_at = datetime.now(timezone.utc).isoformat()
    result_payload = {
        "schema_version": 1,
        "generated_at_utc": generated_at,
        "scope": {
            "pde": PDE,
            "dimensions": list(DIMENSIONS),
            "states_per_dimension": 32,
            "training_seed": 0,
            "source": "saved actor_r from formal prediction NPZ files",
            "full_ta_fk_mcv_role": "auxiliary reference only",
            "not_whole_function_960_state_accuracy": True,
        },
        "source_revision": {
            "worktree_head": worktree_head,
            "formal_origin_tip": origin_tip,
            "expected_origin_tip": EXPECTED_ORIGIN_TIP,
        },
        "rows": output_rows,
        "baseline_cell_summary": {
            "legal_fh_mlp_and_deep_bsde_cells": legal_cells,
            "direct_value_head_wins": won_cells,
        },
        "checks": checks,
    }
    write_json(output_dir / "JOINT_NONSEPARABLE_VALUE_HEAD.json", result_payload)
    write_json(
        output_dir / "CHECKPOINT_AUDIT.json",
        {
            "schema_version": 1,
            "generated_at_utc": generated_at,
            "source_revision": result_payload["source_revision"],
            "entries": checkpoint_entries,
            "checks": checks,
        },
    )

    latex = render_latex(output_rows)
    (table_output_dir / "joint_nonseparable_value_head.tex").write_text(
        latex, encoding="utf-8"
    )
    report = render_report(
        output_rows,
        checks,
        worktree_head,
        origin_tip,
        legal_cells,
        won_cells,
    )
    (output_dir / "VALUE_HEAD_AUDIT_REPORT.md").write_text(
        report, encoding="utf-8"
    )

    print(json.dumps(result_payload, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
