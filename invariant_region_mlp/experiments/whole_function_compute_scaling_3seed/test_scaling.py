from __future__ import annotations

import pytest

from .calibration import (
    _monotone_anchors,
    _nearest_distinct,
    _predicted_candidates,
    _predict_piecewise_log_linear,
    _supplemental_spec,
)
from .protocol import (
    CALIBRATION_ANCHORS,
    METHODS,
    PDES,
    SWEEPS,
    TRAINING_SEEDS,
    calibration_specs,
    settings_for,
)


def test_calibration_covers_every_line() -> None:
    specs = calibration_specs()
    observed = {(spec.pde, spec.method, spec.sweep) for spec in specs}
    expected = {(pde, method, sweep) for pde in PDES for method in METHODS for sweep in SWEEPS}
    assert observed == expected
    assert all(spec.seed == 314159 for spec in specs)
    assert len(specs) == 30


def test_exactly_one_hyperparameter_changes() -> None:
    assert settings_for("tangent_actor", "C1", 256) == {
        "actor_updates": 4, "paths_per_label": 256, "optimizer_steps": 5000,
    }
    assert settings_for("dpi", "C2", 32) == {
        "picard_iterations": 80, "integral_mc_samples": 4096,
        "epochs_per_picard": 32,
    }
    assert settings_for("pinn_hte", "C2", 2048) == {
        "training_epochs": 30000, "hutchinson_vectors": 16,
        "batch_size": 2048,
    }


def test_formal_seeds_are_frozen() -> None:
    assert TRAINING_SEEDS == (0, 1, 2)
    assert set(CALIBRATION_ANCHORS) == {
        (method, sweep) for method in METHODS for sweep in SWEEPS
    }


def test_calibration_model_is_monotone_piecewise_log_linear() -> None:
    anchors = _monotone_anchors([(8, 80.0), (2, 20.0), (4, 15.0)])
    assert anchors == [(2, 20.0), (4, 20.0), (8, 80.0)]
    assert _predict_piecewise_log_linear(anchors, 4) == pytest.approx(20.0)
    assert _predict_piecewise_log_linear(anchors, 6) == pytest.approx(45.0)
    assert _predict_piecewise_log_linear([(2, 20.0)], 8) == pytest.approx(80.0)


def test_budget_reduction_selects_distinct_targets() -> None:
    values = [1, 2, 3, 4, 5]
    predictions = {1: 300.0, 2: 600.0, 3: 900.0, 4: 1200.0, 5: 1800.0}
    assert _nearest_distinct(values, predictions, (300.0, 900.0, 1800.0)) == [1, 3, 5]


def test_slow_ladder_requests_one_measured_doubling() -> None:
    candidates, _predictions, _anchors, inserted, reason = _predicted_candidates(
        "dpi", "C1", [(1024, 350.0), (8192, 400.0)]
    )
    assert inserted == 131072
    assert candidates[-1] == inserted
    assert reason == "ladder_below_25_double"


def test_c0_supplement_executes_the_target_work() -> None:
    spec = _supplemental_spec("b6_fullrank", "dpi", "C0", 450)
    assert spec.settings["picard_iterations"] == 450
    assert spec.checkpoint_values == (450,)
    assert spec.calibration_scale == 1.0
