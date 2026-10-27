"""Bind the frozen banks and choose input-only stratified query subsets."""

from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

import numpy as np

from .protocol import (
    CALIBRATION_COUNT,
    DIMENSION,
    EXPECTED_INPUTS,
    FINAL_COUNT,
    INPUT_DIR,
    PDES,
    RESULTS_DIR,
    SOURCE_BANK_COUNT,
    SOURCE_ROOT,
    SUBSET_SEED,
    SUBSET_STRATA,
    array_sha256,
    ensure_directories,
    file_sha256,
    historical_results_root,
    protocol_payload,
    stable_seed,
    task_specs,
    write_json_atomic,
    write_npz_atomic,
)


def _stratified_indices(t: np.ndarray, pde: str) -> tuple[np.ndarray, np.ndarray]:
    """Take one calibration and two final rows from each time stratum."""

    if t.shape != (SOURCE_BANK_COUNT,):
        raise ValueError(t.shape)
    ordered = np.argsort(t, kind="stable")
    strata = np.array_split(ordered, SUBSET_STRATA)
    if any(len(group) < 3 for group in strata):
        raise RuntimeError("time strata are too small")
    rng = np.random.default_rng(stable_seed(SUBSET_SEED, pde, "input_only_time_strata"))
    calibration: list[int] = []
    final: list[int] = []
    for group in strata:
        chosen = rng.permutation(group)[:3]
        calibration.append(int(chosen[0]))
        final.extend((int(chosen[1]), int(chosen[2])))
    result_cal = np.asarray(calibration, dtype=np.int64)
    result_final = np.asarray(final, dtype=np.int64)
    if len(result_cal) != CALIBRATION_COUNT or len(result_final) != FINAL_COUNT:
        raise RuntimeError((len(result_cal), len(result_final)))
    if len(set(result_cal) & set(result_final)):
        raise RuntimeError("calibration and final subsets overlap")
    return result_cal, result_final


def _git_revision() -> str:
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=SOURCE_ROOT, text=True
    ).strip()


def prepare() -> dict:
    ensure_directories()
    source = historical_results_root()
    protocol = protocol_payload()
    write_json_atomic(RESULTS_DIR / "EXPERIMENT_PROTOCOL.json", protocol)
    manifest: dict = {
        "status": "FROZEN_INPUT_ONLY_QUERY_SELECTION",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "checked_out_revision": _git_revision(),
        "protocol_sha256": protocol["protocol_sha256"],
        "historical_results_root": str(source),
        "selection_used_reference_values": False,
        "selection_rule": (
            "stable-sort 960 inputs by t; split into 32 contiguous strata; "
            "input-seeded permutation; select one calibration and two final rows per stratum"
        ),
        "pdes": {},
    }
    for pde in PDES:
        bank_path = source / "inputs" / f"{pde}_d{DIMENSION}_test_bank.npz"
        parameter_path = source / "inputs" / f"{pde}_d{DIMENSION}_parameters.npz"
        expected = EXPECTED_INPUTS[pde]
        if file_sha256(bank_path) != expected["test_bank_file_sha256"]:
            raise RuntimeError(f"{pde}: frozen bank file hash mismatch")
        if file_sha256(parameter_path) != expected["parameter_file_sha256"]:
            raise RuntimeError(f"{pde}: frozen parameter file hash mismatch")
        # Point selection is completed before the truth array is touched.
        with np.load(bank_path, allow_pickle=False) as payload:
            t = np.asarray(payload["t"], dtype=np.float64)
            x = np.asarray(payload["x"], dtype=np.float64)
        if t.shape != (SOURCE_BANK_COUNT,) or x.shape != (SOURCE_BANK_COUNT, DIMENSION):
            raise RuntimeError(f"{pde}: unexpected bank shapes {(t.shape, x.shape)}")
        if array_sha256(t, x) != expected["test_input_sha256"]:
            raise RuntimeError(f"{pde}: frozen input array hash mismatch")
        calibration, final = _stratified_indices(t, pde)

        input_files: dict[str, str] = {}
        for label, indices in (("calibration", calibration), ("final", final)):
            output = INPUT_DIR / f"{pde}_{label}_inputs.npz"
            write_npz_atomic(
                output,
                original_index=indices,
                t=t[indices],
                x=x[indices],
            )
            input_files[label] = str(output)

        # Truth is copied only after the input-only indices are immutable and
        # is stored in a distinct file never opened by the trainer.
        with np.load(bank_path, allow_pickle=False) as payload:
            truth = np.asarray(payload["truth_u"], dtype=np.float64)
        if truth.shape != (SOURCE_BANK_COUNT,):
            raise RuntimeError(f"{pde}: unexpected truth shape {truth.shape}")
        if array_sha256(truth) != expected["truth_value_sha256"]:
            raise RuntimeError(f"{pde}: frozen truth hash mismatch")
        truth_files: dict[str, str] = {}
        for label, indices in (("calibration", calibration), ("final", final)):
            output = INPUT_DIR / f"{pde}_{label}_truth.npz"
            write_npz_atomic(
                output,
                original_index=indices,
                truth_u=truth[indices],
            )
            truth_files[label] = str(output)

        manifest["pdes"][pde] = {
            "source_bank": str(bank_path),
            "source_bank_file_sha256": file_sha256(bank_path),
            "parameter_file": str(parameter_path),
            "parameter_file_sha256": file_sha256(parameter_path),
            "source_input_sha256": array_sha256(t, x),
            "source_truth_sha256": array_sha256(truth),
            "calibration_indices": calibration.tolist(),
            "final_indices": final.tolist(),
            "calibration_input_sha256": array_sha256(t[calibration], x[calibration]),
            "final_input_sha256": array_sha256(t[final], x[final]),
            "calibration_truth_sha256": array_sha256(truth[calibration]),
            "final_truth_sha256": array_sha256(truth[final]),
            "calibration_input_file": input_files["calibration"],
            "final_input_file": input_files["final"],
            "calibration_truth_file": truth_files["calibration"],
            "final_truth_file": truth_files["final"],
            "calibration_input_file_sha256": file_sha256(Path(input_files["calibration"])),
            "final_input_file_sha256": file_sha256(Path(input_files["final"])),
            "calibration_truth_file_sha256": file_sha256(Path(truth_files["calibration"])),
            "final_truth_file_sha256": file_sha256(Path(truth_files["final"])),
            "calibration_t_range": [float(t[calibration].min()), float(t[calibration].max())],
            "final_t_range": [float(t[final].min()), float(t[final].max())],
        }
    write_json_atomic(RESULTS_DIR / "query_manifest.json", manifest)
    tasks = task_specs()
    task_manifest = {
        "status": "FORMAL_GRID_FROZEN_BEFORE_PREFLIGHT",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "query_manifest_sha256": file_sha256(RESULTS_DIR / "query_manifest.json"),
        "task_count": len(tasks),
        "tasks": tasks,
    }
    write_json_atomic(RESULTS_DIR / "task_manifest.json", task_manifest)
    return manifest


if __name__ == "__main__":
    print(json.dumps(prepare(), indent=2, sort_keys=True))

