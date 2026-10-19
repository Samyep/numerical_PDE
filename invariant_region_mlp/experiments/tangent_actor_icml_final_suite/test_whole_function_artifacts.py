from __future__ import annotations

from pathlib import Path
import json

import pandas as pd


def test_empty_whole_function_grid_has_140_pending_cells(monkeypatch) -> None:
    from . import whole_function_artifacts as artifacts

    monkeypatch.setattr(artifacts, "_read_main", lambda: pd.DataFrame())
    monkeypatch.setattr(artifacts, "_attempt_groups", lambda: {})
    frame, attempts = artifacts.build_results_frame()

    assert len(frame) == 7 * 5 * 4 == 140
    assert set(frame["dimension"]) == {20, 50, 100, 200, 400}
    assert set(frame["method"]) == {"tangent_actor", "dpi", "pinn_hte", "sdgd_pinn"}
    assert set(frame["status"]) == {"PENDING"}
    assert frame[["value_rmse", "value_relative_l2"]].isna().all().all()
    assert attempts == []


def test_stage1_pilot_covers_both_dimensions_and_all_methods() -> None:
    from .whole_function_status import METHODS

    path = Path(__file__).resolve().parent / "slurm/whole_function_stage1_pilot.csv"
    rows = [line.strip().split(",") for line in path.read_text().splitlines() if line.strip()]
    assert len(rows) == len(set(tuple(row) for row in rows)) == 8
    assert {int(row[2]) for row in rows} == {20, 50}
    for dimension in (20, 50):
        assert {row[0] for row in rows if int(row[2]) == dimension} == set(METHODS)
    assert {row[1] for row in rows} == {"joint_concave"}


def test_checkpointed_allocation_is_pending_until_cumulative_cap(
    tmp_path: Path, monkeypatch,
) -> None:
    from . import whole_function_artifacts as artifacts
    from . import whole_function_status as status

    monkeypatch.setattr(status, "STATUS_DIR", tmp_path)
    payload = {
        "method": "tangent_actor", "pde": "joint_concave", "dimension": 20,
        "status": "CHECKPOINTED", "actual_gpu_seconds": 6900,
        "finished_utc": "2026-10-19T00:00:00Z",
    }
    (tmp_path / "attempt.json").write_text(json.dumps(payload))

    budget = status.cell_budget("tangent_actor", "joint_concave", 20)
    assert budget["status"] == "CHECKPOINTED"
    assert budget["used_gpu_seconds"] == 6900
    assert budget["remaining_gpu_seconds"] == 21600 - 6900
    assert artifacts._terminal_status([payload]) == ("PENDING", False)
