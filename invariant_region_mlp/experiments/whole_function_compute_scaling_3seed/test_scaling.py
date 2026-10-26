from __future__ import annotations

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
