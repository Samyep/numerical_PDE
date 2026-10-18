"""Fail closed unless the independent Stage-B review explicitly approves scaling."""

from __future__ import annotations

import json

from .protocol import RESULTS_DIR


def require_approved() -> None:
    path = RESULTS_DIR / "REVIEW_GATE_STATUS.json"
    if not path.exists():
        raise RuntimeError(f"dimension sweep blocked: missing review gate {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("status") != "approved":
        raise RuntimeError(
            "dimension sweep blocked pending independent scientific review; "
            f"current gate status={payload.get('status')!r}"
        )


if __name__ == "__main__":
    require_approved()
    print("independent review gate approved")
