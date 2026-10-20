from __future__ import annotations

from pathlib import Path
import json
from types import SimpleNamespace

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


def test_complete_ledger_waits_for_canonical_aggregation(monkeypatch) -> None:
    from . import whole_function_artifacts as artifacts

    attempt = {
        "method": "dpi", "pde": "joint_concave", "dimension": 20,
        "status": "COMPLETE", "actual_gpu_seconds": 10,
        "finished_utc": "2026-10-20T00:00:00Z",
    }
    monkeypatch.setattr(artifacts, "_read_main", lambda: pd.DataFrame())
    monkeypatch.setattr(
        artifacts, "_attempt_groups",
        lambda: {"dpi__joint_concave__d20": [attempt]},
    )
    frame, _ = artifacts.build_results_frame()
    row = frame[
        (frame["method"] == "dpi")
        & (frame["pde"] == "joint_concave")
        & (frame["dimension"] == 20)
    ].iloc[0]
    assert row["status"] == "PENDING"
    assert pd.isna(row["value_rmse"])


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


def test_attempt_keeps_git_head_from_start(tmp_path: Path, monkeypatch) -> None:
    from . import whole_function_status as status

    monkeypatch.setattr(status, "STATUS_DIR", tmp_path)
    heads = iter(("start-head\n", "finish-head\n"))
    monkeypatch.setattr(status.subprocess, "check_output", lambda *args, **kwargs: next(heads))
    values = {
        "method": "dpi", "pde": "joint_concave", "dimension": 20,
        "job_id": "1", "task_id": "0", "started_utc": "start",
        "finished_utc": "", "actual_gpu_seconds": 0.0, "exit_code": 0,
        "stdout": "out", "stderr": "err", "notes": "running",
    }
    path = status.record(SimpleNamespace(status="RUNNING", **values))
    values.update(finished_utc="finish", actual_gpu_seconds=10.0, notes="done")
    status.record(SimpleNamespace(status="COMPLETE", **values))

    assert json.loads(path.read_text())["git_head"] == "start-head"


def test_balanced_queue_is_frozen_without_accuracy(monkeypatch) -> None:
    from . import whole_function_queue as queue

    rows = []
    for dimension in (20, 50, 200):
        for method in ("tangent_actor", "dpi", "pinn_hte", "sdgd_pinn"):
            rows.append({
                "pde": "joint_concave", "dimension": dimension, "method": method,
                "status": "PENDING", "active_slurm_attempt": False,
            })
    monkeypatch.setattr(queue, "_active_suite_jobs", lambda: [])
    monkeypatch.setattr(queue, "build_results_frame", lambda: (pd.DataFrame(rows), []))
    estimates = pd.DataFrame([
        {"pde": row["pde"], "dimension": row["dimension"], "method": row["method"],
         "estimated_gpu_hours": 1 + index / 10}
        for index, row in enumerate(rows)
    ])
    monkeypatch.setattr(queue.pd, "read_csv", lambda _: estimates)
    first, meta = queue.generate("stage2", slots=4)
    second, _ = queue.generate("stage2", slots=4)
    assert first == second
    assert len(first) == 12
    assert {row["slot"] for row in first} == {0, 1, 2, 3}
    assert all(row["high_gpu_hours"] >= row["estimated_gpu_hours"] for row in first)
    assert meta["ordering"].endswith("no accuracy or test metric used")
    for slot in {row["slot"] for row in first}:
        methods = [row["method"] for row in first if row["slot"] == slot]
        first_resumable = next(
            (index for index, method in enumerate(methods)
             if method not in {"dpi", "pinn_hte"}),
            len(methods),
        )
        assert all(method in {"dpi", "pinn_hte"} for method in methods[:first_resumable])
        assert all(method not in {"dpi", "pinn_hte"} for method in methods[first_resumable:])


def test_slurm_duration_parser_for_dispatch_reservation() -> None:
    from .whole_function_job_accounting import _billing_units, _duration_seconds

    assert _duration_seconds("02:00:00") == 7200
    assert _duration_seconds("1-02:03:04") == 93784
    assert _billing_units("billing=4000,cpu=24,gres/gpu=2,mem=220G") == 4000


def test_historical_compute_falls_back_from_nan_allocation() -> None:
    from .whole_function_artifacts import _historical_gpu_seconds

    assert _historical_gpu_seconds({
        "allocated_gpu_seconds": float("nan"), "total_seconds": 3600.0,
    }) == 3600.0
    assert _historical_gpu_seconds({
        "allocated_gpu_seconds": 1800.0, "total_seconds": 3600.0,
    }) == 1800.0
