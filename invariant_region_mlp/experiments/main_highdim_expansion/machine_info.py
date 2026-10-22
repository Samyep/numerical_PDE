"""Record the exact machine and software context used by Stage A."""

from __future__ import annotations

import importlib.metadata
import json
import os
import platform
from pathlib import Path
import subprocess
import sys
from typing import Any

import torch

from .protocol import RESULTS_DIR, SOURCE_ROOT, write_json_atomic


DEPENDENCIES = (
    "torch",
    "numpy",
    "scipy",
    "pandas",
    "matplotlib",
    "lightning",
    "h5py",
    "yacs",
    "typer",
    "wandb",
    "tensorboardX",
)


def _command(*args: str) -> str:
    try:
        return subprocess.run(
            list(args), check=True, capture_output=True, text=True
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError) as error:
        return f"UNAVAILABLE: {error}"


def _git(path: Path, *args: str) -> str:
    return _command("git", "-C", str(path), *args)


def collect() -> dict[str, Any]:
    dpi_root = Path(os.environ["DPI_ROOT"]).resolve()
    gpu: dict[str, Any] = {
        "available": torch.cuda.is_available(),
        "count": torch.cuda.device_count(),
        "torch_cuda": torch.version.cuda,
        "cudnn": torch.backends.cudnn.version(),
        "nvidia_smi": _command(
            "nvidia-smi",
            "--query-gpu=index,uuid,name,memory.total,driver_version",
            "--format=csv,noheader",
        ),
    }
    if torch.cuda.is_available():
        gpu["devices"] = [
            {
                "index": index,
                "name": torch.cuda.get_device_name(index),
                "capability": list(torch.cuda.get_device_capability(index)),
                "total_memory": torch.cuda.get_device_properties(index).total_memory,
            }
            for index in range(torch.cuda.device_count())
        ]
    dependencies = {}
    for name in DEPENDENCIES:
        try:
            dependencies[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            dependencies[name] = "MISSING"
    job_id = os.environ.get("SLURM_JOB_ID")
    payload = {
        "source_root": str(SOURCE_ROOT.resolve()),
        "git_head": _git(SOURCE_ROOT, "rev-parse", "HEAD"),
        "git_branch": _git(SOURCE_ROOT, "branch", "--show-current"),
        "git_status": _git(SOURCE_ROOT, "status", "--short"),
        "origin": _git(SOURCE_ROOT, "remote", "get-url", "origin"),
        "dpi_root": str(dpi_root),
        "dpi_head": _git(dpi_root, "rev-parse", "HEAD"),
        "hostname": platform.node(),
        "platform": platform.platform(),
        "python": sys.version,
        "python_executable": sys.executable,
        "torch": torch.__version__,
        "gpu": gpu,
        "cpu": _command("lscpu"),
        "memory": _command("free", "-h"),
        "slurm_job_id": job_id,
        "slurm_job": _command("scontrol", "show", "job", "-o", job_id)
        if job_id
        else "not in allocation",
        "slurm_environment": {
            key: value for key, value in sorted(os.environ.items()) if key.startswith("SLURM_")
        },
        "dependencies": dependencies,
    }
    return payload


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    payload = collect()
    write_json_atomic(RESULTS_DIR / "machine_environment.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
