"""Audit immutable inputs inherited from the completed parent experiment."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import torch

from .common import (
    OLD_CHECKPOINT_DIR,
    OLD_PREDICTION_DIR,
    OLD_REFERENCE_DIR,
    OLD_RESULTS_DIR,
    RESULTS_DIR,
    array_sha256,
    atomic_json,
    ensure_dirs,
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run() -> dict:
    ensure_dirs()
    actors = []
    bank_hashes = []
    for seed in range(3):
        run_id = f"primary_d100_k2_lam0p2_eta0_main_mcv_s{seed}"
        checkpoint = OLD_CHECKPOINT_DIR / f"{run_id}_mcv_a4.pt"
        metadata = checkpoint.with_suffix(".json")
        selection_path = OLD_RESULTS_DIR / f"selection_{run_id}.json"
        prediction_path = OLD_PREDICTION_DIR / f"post_freeze_{run_id}.npz"
        for path in (checkpoint, metadata, selection_path, prediction_path):
            if not path.exists():
                raise FileNotFoundError(path)
        selection = json.loads(selection_path.read_text(encoding="utf-8"))
        if selection["selected_actor"] != "a4" or selection["ground_truth_loaded"]:
            raise RuntimeError(f"invalid frozen selection boundary for seed {seed}")
        payload = torch.load(checkpoint, map_location="cpu", weights_only=True)
        if not payload.get("state_dict"):
            raise RuntimeError(f"empty Actor checkpoint {checkpoint}")
        with np.load(prediction_path) as archive:
            bank_hash = array_sha256(archive["test_t"], archive["test_x"])
            shape = list(archive["test_x"].shape)
        bank_hashes.append(bank_hash)
        actors.append({
            "training_seed": seed,
            "selected_actor": selection["selected_actor"],
            "checkpoint": str(checkpoint),
            "checkpoint_sha256": _sha256(checkpoint),
            "metadata_sha256": _sha256(metadata),
            "selection_sha256": _sha256(selection_path),
            "test_bank_sha256": bank_hash,
            "test_bank_shape": shape,
            "state_dict_tensors": len(payload["state_dict"]),
        })
    if len(set(bank_hashes)) != 1:
        raise RuntimeError("prior training seeds do not share the test bank")
    references = []
    for name in (
        "sine_quadratic_d100_k2_lam0p2_eta0p0.npz",
        "sine_quadratic_coupled_d100_k2_lam0p2_eta1p0.npz",
    ):
        path = OLD_REFERENCE_DIR / name
        if not path.exists():
            raise FileNotFoundError(path)
        references.append({"path": str(path), "bytes": path.stat().st_size, "sha256": _sha256(path)})
    result = {
        "parent_commit": "ef2ec67af8da6791e2bbc6086b9b1a1c55227420",
        "actors": actors,
        "shared_test_bank": True,
        "shared_test_bank_sha256": bank_hashes[0],
        "references": references,
        "historical_results_modified": False,
        "passed": True,
    }
    atomic_json(RESULTS_DIR / "INPUT_AUDIT.json", result)
    return result


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))

