"""PACE ICE accelerator preflight and auditable transient-failure handling."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import traceback
from typing import Any

from .protocol import RESULTS_DIR, write_json_atomic


TEMPORARY_FAILURE_EXIT_CODE = 75
EVENT_DIR = RESULTS_DIR / "infrastructure_events"

# Keep this deliberately narrow. Numerical exceptions, CUDA OOM, device-side
# assertions, and application errors are scientific/program failures and must
# remain terminal. These signatures describe an allocated accelerator that the
# CUDA runtime cannot use.
_TRANSIENT_ACCELERATOR_PATTERNS = (
    re.compile(r"CUDA-capable device\(s\) is/are busy or unavailable", re.I),
    re.compile(r"cudaErrorDevicesUnavailable", re.I),
    re.compile(r"CUDA driver initialization failed", re.I),
    re.compile(r"cudaErrorSystemNotReady", re.I),
)
_OUT_OF_MEMORY_PATTERNS = (
    re.compile(r"CUDA out of memory", re.I),
    re.compile(r"OutOfMemoryError", re.I),
    re.compile(r"oom[-_]kill", re.I),
)


def is_transient_accelerator_error(text: str) -> bool:
    """Return true only for retryable device/driver availability failures."""

    if any(pattern.search(text) for pattern in _OUT_OF_MEMORY_PATTERNS):
        return False
    return any(pattern.search(text) for pattern in _TRANSIENT_ACCELERATOR_PATTERNS)


def _command(*args: str) -> str:
    try:
        return subprocess.run(
            list(args), check=True, capture_output=True, text=True, timeout=20
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
        return f"UNAVAILABLE: {error}"


def slurm_node_name() -> str:
    """Return the Slurm node name rather than its DNS-qualified hostname."""

    explicit = os.environ.get("SLURMD_NODENAME")
    if explicit:
        return explicit
    return platform.node().split(".", 1)[0]


def _gpu_inventory() -> str:
    return _command(
        "nvidia-smi",
        "--query-gpu=index,uuid,name,memory.total,driver_version,compute_mode",
        "--format=csv,noheader",
    )


def _safe_token(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "-", value).strip("-") or "unknown"


def record_event(
    *,
    context: str,
    message: str,
    transient: bool,
    traceback_text: str | None = None,
    stderr_path: Path | None = None,
    node_override: str | None = None,
) -> tuple[Path, dict[str, Any]]:
    """Persist an immutable machine-failure event without changing cell data."""

    now = datetime.now(timezone.utc)
    job_id = os.environ.get("SLURM_JOB_ID", "local")
    node = node_override or slurm_node_name()
    payload: dict[str, Any] = {
        "schema_version": 1,
        "status": "TRANSIENT_INFRASTRUCTURE_FAILURE" if transient else "INFRASTRUCTURE_CHECK_FAILED",
        "transient": transient,
        "quarantine_node": transient,
        "context": context,
        "message": message,
        "traceback": traceback_text,
        "stderr": str(stderr_path) if stderr_path else None,
        "job_id": job_id,
        "node": node,
        "hostname": platform.node(),
        "gpu_inventory": _gpu_inventory(),
        "created_utc": now.isoformat(),
    }
    EVENT_DIR.mkdir(parents=True, exist_ok=True)
    stamp = now.strftime("%Y%m%dT%H%M%S.%fZ")
    path = EVENT_DIR / f"{stamp}__job{_safe_token(job_id)}__{_safe_token(context)}.json"
    write_json_atomic(path, payload)
    return path, payload


def quarantined_nodes(event_dir: Path = EVENT_DIR) -> tuple[str, ...]:
    """Return nodes quarantined by recorded transient accelerator failures."""

    nodes: set[str] = set()
    if not event_dir.is_dir():
        return ()
    for path in event_dir.glob("*.json"):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        node = str(payload.get("node") or "").strip()
        audited_transient = (
            payload.get("status") == "TRANSIENT_INFRASTRUCTURE_FAILURE"
            and payload.get("transient") is True
            and payload.get("quarantine_node") is True
        )
        if audited_transient and node:
            nodes.add(node)
    return tuple(sorted(nodes))


def merged_excluded_nodes(extra: str = "", event_dir: Path = EVENT_DIR) -> tuple[str, ...]:
    nodes = set(quarantined_nodes(event_dir))
    nodes.update(value.strip() for value in extra.split(",") if value.strip())
    return tuple(sorted(nodes))


def probe_cuda(context: str) -> int:
    """Exercise the allocated GPU, not merely CUDA discovery metadata."""

    try:
        import torch

        if torch.version.cuda is None:
            raise RuntimeError("installed PyTorch build has no CUDA runtime")
        if not torch.cuda.is_available():
            message = (
                f"CUDA {torch.version.cuda} build reports no available device "
                "inside a Slurm GPU allocation"
            )
            path, payload = record_event(context=context, message=message, transient=True)
            print(json.dumps({**payload, "event": str(path)}, sort_keys=True), file=sys.stderr)
            return TEMPORARY_FAILURE_EXIT_CODE

        torch.cuda.init()
        device = torch.device("cuda:0")
        values = torch.arange(16, device=device, dtype=torch.float32)
        observed = float(torch.sum(values * values).item())
        torch.cuda.synchronize(device)
        torch.cuda.reset_peak_memory_stats(device)
        if observed != 1240.0:
            raise RuntimeError(f"CUDA arithmetic probe returned {observed}, expected 1240.0")
        payload = {
            "status": "PASS",
            "context": context,
            "node": slurm_node_name(),
            "hostname": platform.node(),
            "job_id": os.environ.get("SLURM_JOB_ID"),
            "torch": torch.__version__,
            "torch_cuda": torch.version.cuda,
            "device_name": torch.cuda.get_device_name(device),
            "device_capability": list(torch.cuda.get_device_capability(device)),
            "gpu_inventory": _gpu_inventory(),
        }
        print(json.dumps(payload, sort_keys=True))
        return 0
    except Exception as error:  # the traceback is required for the immutable audit
        trace = traceback.format_exc()
        transient = is_transient_accelerator_error(trace)
        path, payload = record_event(
            context=context,
            message=str(error),
            transient=transient,
            traceback_text=trace,
        )
        print(json.dumps({**payload, "event": str(path)}, sort_keys=True), file=sys.stderr)
        return TEMPORARY_FAILURE_EXIT_CODE if transient else 1


def classify_runtime_failure(
    stderr_path: Path,
    context: str,
    record: bool,
    node_override: str | None = None,
) -> int:
    try:
        text = stderr_path.read_text(encoding="utf-8", errors="replace")
    except OSError as error:
        print(f"cannot read runner stderr {stderr_path}: {error}", file=sys.stderr)
        return 1
    transient = is_transient_accelerator_error(text)
    payload: dict[str, Any] = {
        "status": "TRANSIENT_INFRASTRUCTURE_FAILURE" if transient else "APPLICATION_FAILURE",
        "transient": transient,
        "stderr": str(stderr_path),
    }
    if transient and record:
        path, event = record_event(
            context=context,
            message="runner encountered a retryable CUDA device-availability failure",
            transient=True,
            stderr_path=stderr_path,
            node_override=node_override,
        )
        payload.update({"event": str(path), "node": event["node"]})
    print(json.dumps(payload, sort_keys=True))
    return 0 if transient else 1


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    probe = subparsers.add_parser("probe")
    probe.add_argument("--context", required=True)

    classify = subparsers.add_parser("classify")
    classify.add_argument("--stderr", required=True, type=Path)
    classify.add_argument("--context", required=True)
    classify.add_argument("--record", action="store_true")
    classify.add_argument("--node")

    excluded = subparsers.add_parser("excluded-nodes")
    excluded.add_argument("--extra", default="")

    args = parser.parse_args()
    if args.command == "probe":
        raise SystemExit(probe_cuda(args.context))
    if args.command == "classify":
        raise SystemExit(
            classify_runtime_failure(args.stderr, args.context, args.record, args.node)
        )
    print(",".join(merged_excluded_nodes(args.extra)))


if __name__ == "__main__":
    main()
