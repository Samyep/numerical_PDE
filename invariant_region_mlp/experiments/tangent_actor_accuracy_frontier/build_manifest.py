"""Build reproducibility metadata without embedding large artifacts in Git."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import platform
import subprocess

import numpy as np
import torch

from .common import (
    CHECKPOINT_DIR,
    OLD_CHECKPOINT_DIR,
    OLD_REFERENCE_DIR,
    PREDICTION_DIR,
    REFERENCE_DIR,
    RESULTS_DIR,
    atomic_json,
)


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _entry(path: Path) -> dict:
    return {"path": str(path.resolve()), "bytes": path.stat().st_size, "sha256": _sha(path)}


def run() -> dict:
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    large = []
    patterns = [
        (OLD_CHECKPOINT_DIR, "primary_d100_k2_lam0p2_eta0_main_mcv_s*_mcv_a4.pt"),
        (CHECKPOINT_DIR, "*.pt"),
        (REFERENCE_DIR, "*.npz"),
        (OLD_REFERENCE_DIR, "sine_quadratic*_lam0p2_*.npz"),
        (PREDICTION_DIR, "*.npz"),
    ]
    seen = set()
    for directory, pattern in patterns:
        for path in sorted(directory.glob(pattern)):
            resolved = str(path.resolve())
            if resolved not in seen:
                seen.add(resolved)
                large.append(_entry(path))
    payload = {
        "branch": "ir-mlp-tangent-actor-accuracy-frontier",
        "git_commit": commit,
        "parent_commit": "ef2ec67af8da6791e2bbc6086b9b1a1c55227420",
        "python": platform.python_version(),
        "torch": torch.__version__,
        "numpy": np.__version__,
        "cuda_runtime": torch.version.cuda,
        "large_local_artifacts": large,
        "large_artifacts_committed": False,
        "historical_results_overwritten": False,
        "core_commands": [
            "python -m invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.reference_refinement --family primary --n 192 --dt 0.0000625 --query-points 32768",
            "python -m invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.reference_refinement --family primary --n 256 --dt 0.00003125 --query-points 32768",
            "python -m invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.reference_refinement --family coupled --n 192 --dt 0.0000625 --query-points 32768",
            "python -m invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.reference_refinement --family coupled --n 256 --dt 0.00003125 --query-points 32768",
            "sbatch invariant_region_mlp/experiments/tangent_actor_accuracy_frontier/slurm/run_stage_a_n4_and_upgrade.slurm",
            "sbatch invariant_region_mlp/experiments/tangent_actor_accuracy_frontier/slurm/run_stage_a_replicates.slurm",
            "sbatch invariant_region_mlp/experiments/tangent_actor_accuracy_frontier/slurm/run_nonseparable_seeds01.slurm",
            "sbatch --export=ALL,TRAINING_SEED=2 invariant_region_mlp/experiments/tangent_actor_accuracy_frontier/slurm/run_nonseparable_main.slurm",
            "sbatch invariant_region_mlp/experiments/tangent_actor_accuracy_frontier/slurm/run_coupled_full_history.slurm",
            "sbatch invariant_region_mlp/experiments/tangent_actor_accuracy_frontier/slurm/run_coupled_precision_512.slurm (conditional frozen gate)",
            "sbatch invariant_region_mlp/experiments/tangent_actor_accuracy_frontier/slurm/run_ambient_400_main.slurm",
            "sbatch invariant_region_mlp/experiments/tangent_actor_accuracy_frontier/slurm/run_low_rank_400.slurm (conditional RANK=8,16,32)",
            "sbatch --export=ALL,RANK=<selected-rank>,TRAINING_SEED=0 invariant_region_mlp/experiments/tangent_actor_accuracy_frontier/slurm/run_low_rank_coupled_transfer.slurm (conditional on the frozen helpfulness rule)",
            "python -m invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.analyze_results",
            "python -m invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.write_reports",
            "python -m invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.build_manifest",
            "python -m invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.audit_results",
        ],
        "selection_rule": "independently replicated noise-corrected normalized fixed-point consistency; no Ground Truth",
        "actor_input": "full ambient (t,x)",
        "oracle_subspace_used": False,
    }
    atomic_json(RESULTS_DIR / "reproducibility_manifest.json", payload)
    return payload


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
