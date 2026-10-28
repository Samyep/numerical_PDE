from __future__ import annotations

import inspect
from pathlib import Path

import numpy as np
import torch

from . import accounting, audit, dispatch, final_accuracy, freeze, prepare
from .benchmark import make_benchmark
from .protocol import (
    BANK_SIZES,
    EXPECTED_TA_ARTIFACTS,
    FH_CONFIGS,
    GEOMETRIC_TIMING_COUNTS,
    PREREGISTRATION_COMMIT,
    RESULTS_DIR,
    TA_PATHS,
    TARGETS,
    refined_timing_counts,
    stable_seed,
)


def test_preregistration_and_thresholds_are_frozen() -> None:
    assert PREREGISTRATION_COMMIT == "60c04015682a5f77ff03fd829f575d059620a1bd"
    assert TARGETS == {"primary": 1.0e-2, "secondary": 5.0e-3}
    assert BANK_SIZES == {"calibration": 64, "final": 128, "timing": 2**20}


def test_full_history_candidate_grid_is_exact() -> None:
    assert FH_CONFIGS == (
        (2, 32), (2, 48), (2, 64), (2, 96), (2, 128),
        (2, 192), (2, 256), (2, 384), (2, 512),
        (3, 24), (3, 32), (3, 48), (3, 64), (3, 96),
    )


def test_ta_modes_and_charged_training_are_separate() -> None:
    assert TA_PATHS == (32, 128, 512)
    assert all(row["charged_training_seconds"] > 1300 for row in EXPECTED_TA_ARTIFACTS.values())
    source = inspect.getsource(dispatch._run)
    assert 'candidate"] == "direct"' in source
    assert "MCV_TIMING_TASK_CAP_SECONDS" in source


def test_refinement_preserves_geometric_grid_and_adds_actual_counts() -> None:
    refined = refined_timing_counts(200_000)
    assert set(GEOMETRIC_TIMING_COUNTS).issubset(refined)
    assert len(refined) > len(GEOMETRIC_TIMING_COUNTS)
    assert refined == tuple(sorted(set(refined)))
    assert refined[-1] == 2**20


def test_fresh_bank_seeds_are_split_specific() -> None:
    seeds = {stable_seed("query_bank", name) for name in BANK_SIZES}
    assert len(seeds) == 3


def test_bank_separation_ignores_only_prior_prior_overlap(monkeypatch) -> None:
    prior_shared = np.asarray([[9.0, 9.0]], dtype=np.float64)
    monkeypatch.setattr(
        prepare,
        "_load_previous_small_banks",
        lambda: {"old_a": prior_shared, "old_b": prior_shared.copy()},
    )
    states = {
        "calibration": np.asarray([[1.0, 1.0]], dtype=np.float64),
        "final": np.asarray([[2.0, 2.0]], dtype=np.float64),
        "timing": np.asarray([[3.0, 3.0], [4.0, 4.0]], dtype=np.float64),
    }
    audit = prepare._separation_audit(states)
    assert audit["status"] == "VALIDATED"
    assert not any(
        {row["left"], row["right"]} == {"old_a", "old_b"}
        for row in audit["comparisons"]
    )


def test_joint_driver_formula() -> None:
    benchmark = make_benchmark(100)
    rng = np.random.default_rng(821)
    y_np = rng.normal(size=4)
    z_np = rng.normal(size=(4, 100))
    y = torch.as_tensor(y_np, dtype=torch.float64)
    z = torch.as_tensor(z_np, dtype=torch.float64)
    expected = y_np - y_np**3 - 0.5 * np.square(z_np + y_np[:, None] * benchmark.v).sum(axis=1)
    np.testing.assert_allclose(benchmark.driver(y, z).numpy(), expected, rtol=0, atol=1e-12)


def test_final_truth_is_gated_by_frozen_selection() -> None:
    assert "require_frozen_selections" in inspect.getsource(final_accuracy.aggregate)
    source = inspect.getsource(freeze.freeze)
    assert source.index("_required_paths") < source.index('load_bank("calibration"')


def test_new_results_are_isolated() -> None:
    assert RESULTS_DIR.name == "high_accuracy_fhmlp_break_even"
    assert "single_point_break_even_joint" not in str(RESULTS_DIR)


def test_no_out_of_scope_solver_is_dispatched() -> None:
    source = inspect.getsource(dispatch._run).lower()
    for token in ("deep_bsde", "pinn", "sdgd", "dpi"):
        assert token not in source


def test_accounting_counts_each_slurm_array_allocation(monkeypatch) -> None:
    output = "\n".join([
        "77_0|101|COMPLETED|10|gres/gpu:h100=1|node-a|0:0",
        "77_0.batch|101.batch|COMPLETED|10|gres/gpu:h100=1|node-a|0:0",
        "77_1|102|COMPLETED|12|gres/gpu:h100=1|node-b|0:0",
        "77_1.extern|102.extern|COMPLETED|12|gres/gpu:h100=1|node-b|0:0",
    ])
    monkeypatch.setattr(accounting.subprocess, "check_output", lambda *args, **kwargs: output)
    rows = accounting._sacct("77")
    assert [row["job_id"] for row in rows] == ["77_0", "77_1"]
    assert [row["job_id_raw"] for row in rows] == ["101", "102"]
    assert sum(row["elapsed_raw_seconds"] for row in rows) == 22


def test_artifact_audit_has_source_root_for_relative_hash_manifest() -> None:
    assert audit.SOURCE_ROOT.is_dir()
