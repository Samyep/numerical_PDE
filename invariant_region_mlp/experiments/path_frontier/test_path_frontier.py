from __future__ import annotations

import pandas as pd

from .analyze import _criteria, _registered_cost_levels
from .protocol import (
    EXPECTED_COMBINED_ROWS,
    EXPECTED_PATH_ROWS,
    PROBLEMS,
    has_reusable_path,
    new_task_path,
    path_tasks,
    source_path,
)
from .run import _fingerprint_worker


def test_registered_matrix_and_frozen_reuse_partition() -> None:
    tasks = path_tasks()
    assert len(tasks) == EXPECTED_PATH_ROWS == 480
    identities = {
        (task["pde_id"], task["d"], task["n"], task["M"], task["repetition"])
        for task in tasks
    }
    assert len(identities) == len(tasks)
    reused = [task for task in tasks if has_reusable_path(task)]
    assert len(reused) == 160
    assert EXPECTED_COMBINED_ROWS == 1920


def test_new_artifacts_cannot_overlap_frozen_results() -> None:
    missing = next(task for task in path_tasks() if not has_reusable_path(task))
    frozen, source = source_path(missing, "raw")
    assert source == "double_estimator"
    assert frozen != new_task_path(missing)
    assert "results\\double_estimator" in str(frozen) or "results/double_estimator" in str(frozen)
    assert "results\\path_frontier" in str(new_task_path(missing)) or "results/path_frontier" in str(new_task_path(missing))


def test_cost_levels_are_exactly_the_frozen_d4_levels() -> None:
    levels = _registered_cost_levels()
    assert set(levels) == set(PROBLEMS)
    for values in levels.values():
        assert len(values) == 10
        assert values[0] == 9600.0
        assert values[-1] == 1267200.0


def test_full_recursive_path_and_raw_chunk_fingerprints_match() -> None:
    task = {
        "stage": "path_frontier",
        "pde_id": "P1",
        "d": 100,
        "n": 2,
        "M": 8,
        "method": "path",
        "repetition": 0,
        "chunk_size": 4,
    }
    audit = _fingerprint_worker(task)
    assert audit["match"]
    assert audit["raw_fingerprint"]
    assert audit["raw_terminal_samples"] == audit["path_terminal_samples"]
    assert audit["raw_transition_samples"] == audit["path_transition_samples"]


def test_mechanical_criteria_thresholds() -> None:
    records = []
    for problem_index, (pde_id, d) in enumerate(PROBLEMS):
        for level in range(10):
            records.append(
                {
                    "pde_id": pde_id,
                    "dimension": d,
                    "complete": True,
                    "double_adds": level < (8 if problem_index < 3 else 7),
                    "path_beats_raw": level < 8,
                }
            )
    result = _criteria(pd.DataFrame(records))
    assert result["E-1"]["verdict"] == "Double adds value at the frontier"
    assert result["E-1"]["problems_meeting_criterion"] == 3
    assert result["E-2"]["verdict"] == "PASS"

    frame = pd.DataFrame(records)
    frame["double_adds"] = False
    result = _criteria(frame)
    assert result["E-1"]["verdict"] == "pathwise alone suffices at the frontier"
