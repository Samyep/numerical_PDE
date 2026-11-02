#!/usr/bin/env python3
"""Extract frozen Joint-Nonseparable direct-head metrics without inference.

The formal single-point Joint-Nonseparable suite already writes the direct
Actor value prediction into every selected-Actor prediction NPZ.  This tool
reads those frozen artifacts, verifies the four inference replicas agree on
the deterministic direct head, and exports one audited row per dimension.
It performs no training, checkpoint selection, or GPU work.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np


DIMENSIONS = (20, 50, 100, 200, 400)
RUN_TEMPLATE = "joint_nonseparable_d{dimension}_mature_s0"


def file_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path):
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def resolve_record_path(results_root, value):
    path = Path(value)
    if path.is_file():
        return path
    # Formal records contain absolute ICE paths.  The fallback makes the
    # extractor work after a result tree has been relocated.
    marker = "invariant_region_mlp/results/single_point_joint_extension/"
    if marker in value:
        relative = value.split(marker, 1)[1]
        relocated = results_root / relative
        if relocated.is_file():
            return relocated
    raise FileNotFoundError(value)


def extract_dimension(results_root, dimension):
    run_id = RUN_TEMPLATE.format(dimension=dimension)
    training_path = results_root / "raw" / "actor" / f"training_{run_id}.json"
    training = load_json(training_path)
    if training.get("status") != "complete":
        raise RuntimeError(f"incomplete training record: {training_path}")
    selection = training["selection"]
    selected_actor = selection["selected_actor"]
    checkpoint = resolve_record_path(results_root, selection["selected_checkpoint"])
    expected_checkpoint_sha = next(
        phase["checkpoint_sha256"]
        for phase in training["phases"]
        if f"a{phase['iteration']}" == selected_actor
    )
    observed_checkpoint_sha = file_sha256(checkpoint)
    if observed_checkpoint_sha != expected_checkpoint_sha:
        raise RuntimeError(f"checkpoint SHA mismatch for {run_id}")

    evaluation_dir = (
        results_root / "raw" / "actor_evaluations" / run_id / "selected"
    )
    records = sorted(evaluation_dir.glob("steps64_rep*.json"))
    if len(records) != 4:
        raise RuntimeError(f"expected four selected-Actor replicas for {run_id}")
    payloads = [load_json(path) for path in records]
    for payload in payloads:
        if payload["resolved_actor"] != selected_actor:
            raise RuntimeError(f"selected Actor mismatch for {run_id}")
        if int(payload["dimension"]) != dimension:
            raise RuntimeError(f"dimension mismatch for {run_id}")
        if int(payload["test_states"]) != 32:
            raise RuntimeError(f"unexpected query count for {run_id}")

    recorded_head = np.asarray(
        [float(payload["actor_value_rmse"]) for payload in payloads],
        dtype=np.float64,
    )
    if not np.all(recorded_head == recorded_head[0]):
        raise RuntimeError(f"direct head changed across inference replicas: {run_id}")
    bank_hashes = {payload["test_bank_sha256"] for payload in payloads}
    bank_file_hashes = {payload["test_bank_file_sha256"] for payload in payloads}
    if len(bank_hashes) != 1 or len(bank_file_hashes) != 1:
        raise RuntimeError(f"test bank changed across replicas: {run_id}")

    prediction_path = resolve_record_path(results_root, payloads[0]["prediction_path"])
    with np.load(prediction_path, allow_pickle=False) as arrays:
        actor_value = np.asarray(arrays["actor_r"], dtype=np.float64).reshape(-1)
        truth = np.asarray(arrays["truth"], dtype=np.float64).reshape(-1)
    if actor_value.shape != truth.shape or actor_value.size != 32:
        raise RuntimeError(f"unexpected direct-head array shape for {run_id}")
    error = actor_value - truth
    rmse = float(np.sqrt(np.mean(np.square(error))))
    relative_l2 = float(np.linalg.norm(error) / np.linalg.norm(truth))
    rmae = float(np.abs(error).sum() / np.abs(truth).sum())
    if not np.isclose(rmse, recorded_head[0], rtol=2e-7, atol=1e-10):
        raise RuntimeError(f"recomputed direct-head RMSE mismatch for {run_id}")

    return {
        "pde": "joint_nonseparable",
        "pde_label": "Joint Nonseparable",
        "dimension": dimension,
        "training_seed": 0,
        "selected_actor": selected_actor,
        "states": actor_value.size,
        "direct_value_head_rmse": rmse,
        "direct_value_head_relative_l2": relative_l2,
        "direct_value_head_rmae": rmae,
        "test_bank_sha256": next(iter(bank_hashes)),
        "test_bank_file_sha256": next(iter(bank_file_hashes)),
        "checkpoint": str(checkpoint),
        "checkpoint_sha256": observed_checkpoint_sha,
        "training_record": str(training_path),
        "training_record_sha256": file_sha256(training_path),
        "prediction_record": str(prediction_path),
        "prediction_record_sha256": file_sha256(prediction_path),
        "replicate_records": [str(path) for path in records],
        "replicate_record_sha256": [file_sha256(path) for path in records],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--results-root",
        type=Path,
        required=True,
        help="Path to invariant_region_mlp/results/single_point_joint_extension",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.results_root.resolve()
    rows = [extract_dimension(root, dimension) for dimension in DIMENSIONS]
    args.output_dir.mkdir(parents=True, exist_ok=True)

    json_path = args.output_dir / "joint_nonseparable_direct_value_head.json"
    csv_path = args.output_dir / "joint_nonseparable_direct_value_head.csv"
    json_path.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    csv_rows = [
        {key: value for key, value in row.items() if not isinstance(value, list)}
        for row in rows
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(csv_rows[0]))
        writer.writeheader()
        writer.writerows(csv_rows)
    print(json.dumps({"csv": str(csv_path), "json": str(json_path), "rows": rows}, indent=2))


if __name__ == "__main__":
    main()
