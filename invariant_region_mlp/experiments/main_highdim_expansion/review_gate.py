"""Fail-closed gate for every formal expansion cell."""

from __future__ import annotations

import json
import subprocess

from .prepare import file_sha256
from .protocol import MANIFEST_PATH, VALIDATION_PATH, frozen_protocol


def review() -> dict[str, object]:
    if not MANIFEST_PATH.is_file() or not VALIDATION_PATH.is_file():
        raise RuntimeError("frozen inputs and mathematical validation are required")
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    validation = json.loads(VALIDATION_PATH.read_text(encoding="utf-8"))
    checks = {
        "manifest_frozen_before_training": manifest.get("status") == "FROZEN_BEFORE_TRAINING",
        "protocol_exact": manifest.get("protocol") == frozen_protocol(),
        "mathematical_validation_pass": validation.get("status") == "VALIDATED" and validation.get("overall_pass") is True,
        "validation_binds_manifest": validation.get("manifest_sha256") == file_sha256(MANIFEST_PATH),
        "ten_parameter_bank_cells": len(manifest.get("cells", [])) == 10,
        "no_nonfinite_frozen_inputs": all(
            row.get("nonfinite_parameters") == 0 and row.get("nonfinite_bank") == 0
            for row in manifest.get("cells", [])
        ),
    }
    active_old = subprocess.run(
        ["squeue", "-h", "-u", __import__("getpass").getuser(), "-o", "%j"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.splitlines()
    checks["old_suite_not_active"] = not any(name.startswith("icml_whole") for name in active_old)
    result = {"status": "PASS" if all(checks.values()) else "FAIL", "checks": checks}
    if result["status"] != "PASS":
        raise RuntimeError(json.dumps(result, indent=2, sort_keys=True))
    return result


def main() -> None:
    print(json.dumps(review(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
