"""Hash compact final B1 artifacts while leaving local binary arrays out."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from .common import PROJECT_DIR, RESULTS_DIR, atomic_json


ALLOWED_SUFFIXES = {".csv", ".json", ".md", ".pdf", ".png"}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def make_manifest() -> dict[str, object]:
    output = RESULTS_DIR / "artifact_manifest.json"
    candidates = [
        path
        for path in RESULTS_DIR.rglob("*")
        if path.is_file()
        and path != output
        and path.suffix.lower() in ALLOWED_SUFFIXES
    ]
    report = PROJECT_DIR / "docs" / "MAIN_B1_SIXMETHOD_REPORT.md"
    if report.is_file():
        candidates.append(report)
    artifacts = [
        {
            "path": str(path.relative_to(PROJECT_DIR)),
            "bytes": path.stat().st_size,
            "sha256": _sha256(path),
        }
        for path in sorted(set(candidates))
    ]
    payload: dict[str, object] = {
        "schema_version": 1,
        "excluded_local_binary_suffixes": [".npz", ".pt"],
        "artifact_count": len(artifacts),
        "artifacts": artifacts,
    }
    atomic_json(output, payload)
    return payload


def main() -> None:
    payload = make_manifest()
    print(json.dumps({"artifact_count": payload["artifact_count"]}, indent=2))


if __name__ == "__main__":
    main()
