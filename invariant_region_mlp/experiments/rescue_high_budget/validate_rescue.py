"""Integrity audit for Funding/HJB high-budget rescue artifacts."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
from typing import Any

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[1]
REPOSITORY_ROOT = PROJECT_ROOT.parent
RESULT_ROOT = PROJECT_ROOT / "results" / "rescue_high_budget"


def metadata(path: Path) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    with np.load(path, allow_pickle=False) as data:
        payload = json.loads(str(data["metadata_json"]))
        arrays = {name: data[name].copy() for name in data.files if name != "metadata_json"}
    return payload, arrays


def main() -> None:
    checks: dict[str, Any] = {}

    manifests = sorted(RESULT_ROOT.glob("*_manifest.json"))
    manifest_payloads = [json.loads(path.read_text(encoding="utf-8")) for path in manifests]
    checks["manifest_count"] = len(manifests)
    checks["all_manifests_complete"] = all(
        item.get("status") == "complete"
        and item.get("tasks_completed") == item.get("tasks_total")
        for item in manifest_payloads
    )

    funding_root = RESULT_ROOT / "raw" / "funding" / "seed20261005"
    funding_pairs = 0
    funding_root_terminal_paired = True
    funding_draws_paired = True
    all_funding_finite = True
    all_funding_fzero_zero = True
    all_corrected_ebl = True
    funding_nonfinite_counters_zero = True
    for raw_path in sorted(funding_root.glob("n*_M*/raw/block_*.npz")):
        sample_path = raw_path.parent.parent / "samplewise" / raw_path.name
        if not sample_path.exists():
            continue
        raw_meta, raw_arrays = metadata(raw_path)
        ir_meta, ir_arrays = metadata(sample_path)
        funding_pairs += 1
        funding_root_terminal_paired &= np.array_equal(
            raw_arrays["root_terminal_state"], ir_arrays["root_terminal_state"]
        )
        funding_draws_paired &= (
            raw_meta["draw_fingerprint"] == ir_meta["draw_fingerprint"]
        )
        all_corrected_ebl &= (
            raw_meta["corrected_terminal_ebl"] == "standard_normal / sqrt(T-t)"
            and ir_meta["corrected_terminal_ebl"] == "standard_normal / sqrt(T-t)"
        )
        for arrays in (raw_arrays, ir_arrays):
            all_funding_finite &= all(np.all(np.isfinite(array)) for array in arrays.values())
        funding_nonfinite_counters_zero &= (
            raw_meta["work"]["nonfinite_states"] == 0
            and raw_meta["work"]["nonfinite_generators"] == 0
            and ir_meta["work"]["nonfinite_states"] == 0
            and ir_meta["work"]["nonfinite_generators"] == 0
        )
    for path in funding_root.glob("n*_M*/f_zero/block_*.npz"):
        meta, arrays = metadata(path)
        all_funding_fzero_zero &= (
            meta["work"]["f_evals"] == 0
            and np.array_equal(
                arrays["prediction_state"][:, 0], arrays["root_terminal_state"][:, 0]
            )
            and np.count_nonzero(arrays["nonlinear_u_correction"]) == 0
        )
    checks.update(
        {
            "funding_raw_ir_paired_blocks": funding_pairs,
            "funding_root_terminal_bitwise_paired": funding_root_terminal_paired,
            "funding_draw_fingerprints_paired": funding_draws_paired,
            "all_funding_arrays_finite": all_funding_finite,
            "funding_fzero_has_zero_correction_and_f_calls": all_funding_fzero_zero,
            "funding_nonfinite_counters_zero": funding_nonfinite_counters_zero,
        }
    )

    hjb_root = RESULT_ROOT / "raw" / "hjb"
    hjb_pairs = 0
    hjb_root_terminal_paired = True
    hjb_draws_paired = True
    all_hjb_finite = True
    hjb_nonfinite_counters_zero = True
    no_heuristic_clipping = True
    all_hjb_fzero_zero = True
    for stage in ("h1", "h2"):
        for raw_path in sorted((hjb_root / stage).glob("*_raw_*.npz")):
            sample_path = Path(str(raw_path).replace("_raw_", "_samplewise_"))
            if not sample_path.exists():
                continue
            raw_meta, raw_arrays = metadata(raw_path)
            ir_meta, ir_arrays = metadata(sample_path)
            hjb_pairs += 1
            hjb_root_terminal_paired &= np.array_equal(
                raw_arrays["root_terminal_u"], ir_arrays["root_terminal_u"]
            )
            hjb_draws_paired &= (
                raw_meta["draw_fingerprint"] == ir_meta["draw_fingerprint"]
            )
            no_heuristic_clipping &= (
                raw_meta["heuristic_clipping"] is False
                and ir_meta["heuristic_clipping"] is False
            )
            all_corrected_ebl &= (
                raw_meta["corrected_terminal_ebl"] == "standard_normal / sqrt(T-t)"
                and ir_meta["corrected_terminal_ebl"] == "standard_normal / sqrt(T-t)"
            )
            for arrays in (raw_arrays, ir_arrays):
                all_hjb_finite &= all(np.all(np.isfinite(array)) for array in arrays.values())
            hjb_nonfinite_counters_zero &= (
                raw_meta["work"]["nonfinite_states"] == 0
                and raw_meta["work"]["nonfinite_generators"] == 0
                and ir_meta["work"]["nonfinite_states"] == 0
                and ir_meta["work"]["nonfinite_generators"] == 0
            )
    for stage in ("h1", "h2_fzero"):
        for path in (hjb_root / stage).glob("*_f_zero_*.npz"):
            meta, arrays = metadata(path)
            all_hjb_fzero_zero &= (
                meta["work"]["f_evals"] == 0
                and np.array_equal(arrays["prediction_u"], arrays["root_terminal_u"])
                and np.count_nonzero(arrays["nonlinear_u_correction"]) == 0
            )
            all_hjb_finite &= all(np.all(np.isfinite(array)) for array in arrays.values())
    checks.update(
        {
            "hjb_raw_ir_paired_repetitions": hjb_pairs,
            "hjb_root_terminal_bitwise_paired": hjb_root_terminal_paired,
            "hjb_draw_fingerprints_paired": hjb_draws_paired,
            "all_hjb_arrays_finite": all_hjb_finite,
            "hjb_nonfinite_counters_zero": hjb_nonfinite_counters_zero,
            "hjb_fzero_has_zero_correction_and_f_calls": all_hjb_fzero_zero,
            "no_heuristic_clipping": no_heuristic_clipping,
            "all_corrected_ebl_metadata": all_corrected_ebl,
        }
    )

    legacy_hjb = json.loads(
        (RESULT_ROOT / "hjb_legacy_reproduction.json").read_text(encoding="utf-8")
    )
    checks["hjb_legacy_headline_exactly_reproduced"] = (
        legacy_hjb["max_absolute_difference"] == 0.0
    )
    funding_summary = pd.read_csv(RESULT_ROOT / "funding_summary.csv")
    old_raw = funding_summary[
        (funding_summary.n == 2)
        & (funding_summary.M == 10)
        & (funding_summary.method == "raw")
    ].iloc[0]
    old_ir = funding_summary[
        (funding_summary.n == 2)
        & (funding_summary.M == 10)
        & (funding_summary.method == "samplewise")
    ].iloc[0]
    checks["funding_legacy_headline_exactly_reproduced"] = (
        old_raw.repetitions == 100
        and old_ir.repetitions == 100
        and abs(old_raw.mae - 1.5268272487666519) < 5e-15
        and abs(old_ir.mae - 0.5214655268205618) < 5e-15
    )
    funding_validation = pd.read_csv(
        RESULT_ROOT / "funding_shrink_validation.csv"
    )
    validation_winners = (
        funding_validation.loc[
            funding_validation.groupby(["n", "M"])["mae"].idxmin(),
            ["n", "M", "method"],
        ]
        .sort_values(["n", "M"])
        .reset_index(drop=True)
    )
    checks["funding_shrink_validation_configs"] = len(validation_winners)
    checks["funding_validation_selects_c025_all_three"] = (
        len(validation_winners) == 3
        and set(validation_winners.method) == {"shrink_c0.25"}
    )

    funding_repetitions = pd.read_csv(RESULT_ROOT / "funding_repetitions.csv")
    hjb_repetitions = pd.read_csv(RESULT_ROOT / "hjb_repetitions.csv")
    checks["funding_repetition_rows"] = len(funding_repetitions)
    checks["hjb_repetition_rows"] = len(hjb_repetitions)
    checks["csv_metrics_all_finite"] = (
        np.all(np.isfinite(funding_repetitions["prediction"]))
        and np.all(np.isfinite(hjb_repetitions["value_relative_l2"]))
    )

    hjb_summary = pd.read_csv(RESULT_ROOT / "hjb_summary.csv")
    h2_m96 = hjb_summary[
        (hjb_summary.stage == "h2") & (hjb_summary.M == 96)
    ]
    checks["hjb_h2_dimensions_complete"] = set(h2_m96.dimension) == {100, 140, 160}
    for dimension in (100, 140, 160):
        frame = h2_m96[h2_m96.dimension == dimension]
        ir = float(frame[frame.method == "samplewise"].value_relative_l2_mean.iloc[0])
        zero = float(frame[frame.method == "f_zero"].value_relative_l2_mean.iloc[0])
        checks[f"hjb_d{dimension}_m96_ir_above_fzero_floor"] = ir > 50.0 * zero
    checks["hjb_h3_not_run_after_stopping_decision"] = not (hjb_root / "h3").exists()

    report_text = "\n".join(
        path.read_text(encoding="utf-8")
        for path in (
            PROJECT_ROOT / "docs" / "FUNDING_HIGH_BUDGET_RESCUE.md",
            PROJECT_ROOT / "docs" / "HJB_HIGH_BUDGET_RESCUE.md",
            PROJECT_ROOT / "docs" / "OLD_BENCHMARK_RESCUE_SUMMARY.md",
        )
    )
    checks["reports_contain_required_verdicts"] = (
        "FUNDING VERDICT F-B" in report_text and "HJB VERDICT H-B" in report_text
    )
    status = subprocess.check_output(
        ["git", "-C", str(REPOSITORY_ROOT), "status", "--porcelain"], text=True
    )
    checks["no_manuscript_path_changed"] = not any(
        "invariant_region_mlp/paper/" in line for line in status.splitlines()
    )

    checks = {
        key: (
            bool(value)
            if isinstance(value, np.bool_)
            else int(value)
            if isinstance(value, np.integer)
            else float(value)
            if isinstance(value, np.floating)
            else value
        )
        for key, value in checks.items()
    }
    boolean_checks = {key: value for key, value in checks.items() if isinstance(value, bool)}
    checks["all_boolean_checks_pass"] = all(boolean_checks.values())
    output = RESULT_ROOT / "validation_audit.json"
    output.write_text(json.dumps(checks, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(checks, indent=2, sort_keys=True))
    if not checks["all_boolean_checks_pass"]:
        raise SystemExit("rescue integrity audit failed")


if __name__ == "__main__":
    main()
