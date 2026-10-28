"""Materialize the preregistered protocol before formal GPU execution."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess

from .protocol import PROJECT_DIR, RESULTS_DIR, ensure_dirs, frozen_protocol


SOURCES = (
    "experiments/tangent_actor_icml_final_suite/run_standard_actor.py",
    "experiments/tangent_actor_icml_final_suite/run_joint_actor.py",
    "experiments/original_pde_benchmark_comparison/run_actor.py",
    "experiments/original_pde_benchmark_comparison/tangent_solver.py",
    "experiments/original_pde_benchmark_comparison/models.py",
    "experiments/tangent_actor_scaling_diagnosis/models.py",
    "experiments/tangent_actor_scaling_diagnosis/training.py",
    "experiments/tangent_actor_mcv_training/labels.py",
    "experiments/tangent_actor_accuracy_frontier/solver.py",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def main() -> None:
    ensure_dirs()
    protocol = frozen_protocol()
    write_json(RESULTS_DIR / "FROZEN_PROTOCOL.json", protocol)
    head = subprocess.check_output(
        ["git", "-C", str(PROJECT_DIR.parent), "rev-parse", "HEAD"], text=True
    ).strip()
    source_rows = []
    for relative in SOURCES:
        path = PROJECT_DIR / relative
        source_rows.append({
            "path": str(path.relative_to(PROJECT_DIR)),
            "sha256": sha256(path),
            "exists": path.exists(),
        })
    write_json(RESULTS_DIR / "SOURCE_FIDELITY_AUDIT.json", {
        "base_commit": head,
        "protocol_sha256": protocol["protocol_sha256"],
        "source_files": source_rows,
        "frozen_sources_modified": False,
        "status": "preregistered_before_formal_gpu_results",
    })
    print(json.dumps({"protocol_sha256": protocol["protocol_sha256"], "base_commit": head}))


if __name__ == "__main__":
    main()
