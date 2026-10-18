from __future__ import annotations

import pandas as pd
import pytest

from .build_publication_artifacts import _select_native_primary_rows


def test_native_table_uses_the_single_frozen_primary_seed() -> None:
    rows = []
    for seed in (0, 1):
        for steps in (64, 128):
            for estimator in ("plain", "mcv"):
                rows.append({
                    "pde": "joint_convex", "dimension": 100,
                    "method": "tangent_actor", "method_status": "adapted",
                    "seed": seed, "steps": steps, "estimator": estimator,
                    "absolute_error": seed + steps / 1000,
                })
    selected = _select_native_primary_rows(pd.DataFrame(rows).query("seed == 0"))
    assert len(selected) == 1
    assert set(selected["seed"]) == {0}
    assert set(selected["steps"]) == {64}
    assert set(selected["estimator"]) == {"mcv"}


def test_native_table_rejects_ambiguous_duplicate_baseline_rows() -> None:
    duplicate = pd.DataFrame([
        {"pde": "hjb_log", "dimension": 100, "method": "deep_bsde",
         "method_status": "native_official", "seed": 0, "steps": 20,
         "estimator": "native", "absolute_error": 0.1},
        {"pde": "hjb_log", "dimension": 100, "method": "deep_bsde",
         "method_status": "native_official", "seed": 0, "steps": 20,
         "estimator": "native", "absolute_error": 0.2},
    ])
    with pytest.raises(RuntimeError, match="exactly one primary estimate"):
        _select_native_primary_rows(duplicate)
