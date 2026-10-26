"""Freeze the configuration grid and bind it to the historical 100D inputs."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
from typing import Any

import numpy as np

from invariant_region_mlp.experiments.main_highdim_expansion.benchmarks import array_sha256

from .common import file_sha256, write_json_atomic
from .protocol import (
    BASE_REVISION,
    CONFIGURATION_GRID,
    DIMENSION,
    EXPECTED_INPUTS,
    PDES,
    RESULTS_DIR,
    SOURCE_ROOT,
    UPSTREAM_DPI_REVISION,
    ensure_directories,
    protocol_payload,
)


def historical_results_root() -> Path:
    configured = os.environ.get("WCF_HISTORICAL_RESULTS")
    candidates = [
        Path(configured) if configured else None,
        SOURCE_ROOT / "invariant_region_mlp" / "results" / "main_highdim_expansion",
    ]
    for candidate in candidates:
        if candidate is not None and (
            candidate / "inputs" / f"b6_fullrank_d{DIMENSION}_test_bank.npz"
        ).is_file():
            return candidate.resolve()
    raise FileNotFoundError(
        "set WCF_HISTORICAL_RESULTS to the audited main_highdim_expansion result root"
    )


def _git_revision(path: Path) -> str:
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=path, text=True
    ).strip()


def freeze() -> dict[str, Any]:
    ensure_directories()
    source = historical_results_root()
    manifest = json.loads(
        (source / "parameter_and_test_bank_manifest.json").read_text(encoding="utf-8")
    )
    cells = {
        (str(row["pde"]), int(row["dimension"])): row
        for row in manifest["cells"]
    }
    frozen_inputs: list[dict[str, Any]] = []
    for pde in PDES:
        row = cells[(pde, DIMENSION)]
        bank = source / "inputs" / f"{pde}_d{DIMENSION}_test_bank.npz"
        parameters = source / "inputs" / f"{pde}_d{DIMENSION}_parameters.npz"
        with np.load(bank, allow_pickle=False) as payload:
            t = np.asarray(payload["t"])
            x = np.asarray(payload["x"])
            truth = np.asarray(payload["truth_u"])
        observed = {
            "test_input_sha256": array_sha256(t, x),
            "test_bank_file_sha256": file_sha256(bank),
            "truth_value_sha256": array_sha256(truth),
            "parameter_file_sha256": file_sha256(parameters),
        }
        if observed != EXPECTED_INPUTS[pde]:
            raise RuntimeError(f"frozen input drift for {pde}: {observed}")
        if t.shape != (960,) or x.shape != (960, 100) or truth.shape != (960,):
            raise RuntimeError(f"unexpected bank shapes for {pde}")
        frozen_inputs.append({
            "pde": pde,
            "dimension": DIMENSION,
            "test_bank": str(bank),
            "parameter_file": str(parameters),
            **observed,
        })
    protocol = protocol_payload()
    write_json_atomic(RESULTS_DIR / "EXPERIMENT_PROTOCOL.json", protocol)
    payload = {
        "status": "FROZEN_BEFORE_FORMAL_TRAINING",
        "protocol_sha256": protocol["protocol_sha256"],
        "base_revision": BASE_REVISION,
        "checked_out_revision": _git_revision(SOURCE_ROOT),
        "upstream_dpi_revision": UPSTREAM_DPI_REVISION,
        "configuration_grid": CONFIGURATION_GRID,
        "frozen_inputs": frozen_inputs,
        "historical_results_root": str(source),
        "no_sdgd_jobs": True,
        "test_truth_permitted_only_after_all_18_trajectories_stop": True,
    }
    write_json_atomic(RESULTS_DIR / "configuration_manifest.json", payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    print(json.dumps(freeze(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

