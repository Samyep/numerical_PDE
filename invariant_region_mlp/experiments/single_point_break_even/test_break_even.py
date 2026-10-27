"""Fast structural tests; formal CUDA checks live in :mod:`preflight`."""

from __future__ import annotations

import json

import numpy as np

from . import deep_bsde, tangent_actor
from .protocol import (
    DB_CANDIDATE_UPDATES,
    DIMENSION,
    FH_CONFIGS,
    REVISION,
    TA_PATHS,
    TARGETS,
    canonical_sha256,
    stable_seed,
)


def test_frozen_constants() -> None:
    assert REVISION == "single_point_break_even_v1"
    assert DIMENSION == 100
    assert TARGETS == {"epsilon_1": 2.5e-4, "epsilon_2": 1.0e-4}
    assert TA_PATHS == (32, 64, 128, 256, 512, 1024)
    assert DB_CANDIDATE_UPDATES == (2000, 4000, 8000)
    assert FH_CONFIGS[0] == (3, 12) and FH_CONFIGS[-1] == (4, 16)


def test_seed_sequence_is_stable_and_namespaced() -> None:
    first = stable_seed("bank", "calibration")
    assert first == stable_seed("bank", "calibration")
    assert first != stable_seed("bank", "final")
    assert 0 <= first < 2**32


def test_canonical_hash_ignores_mapping_order() -> None:
    assert canonical_sha256({"a": 1, "b": 2}) == canonical_sha256({"b": 2, "a": 1})


def test_deep_bsde_configuration_is_exact() -> None:
    config = deep_bsde.CONFIG
    assert config.intervals == 20
    assert config.hidden_widths == (110, 110)
    assert config.batch_size == 64
    assert config.validation_size == 256
    assert config.learning_rate == 5e-4
    assert config.dtype == "float64"


def test_json_protocol_values_are_finite() -> None:
    payload = {"targets": TARGETS, "paths": TA_PATHS, "configs": FH_CONFIGS}
    encoded = json.dumps(payload, allow_nan=False)
    assert np.isfinite(len(encoded))


def test_tangent_actor_timing_cap_includes_training() -> None:
    assert tangent_actor._next_prefix_fits(
        fixed_seconds=800.0,
        inference_seconds=100.0,
        previous_queries=100,
        next_queries=200,
        max_active_seconds=3600.0,
    )
    assert not tangent_actor._next_prefix_fits(
        fixed_seconds=800.0,
        inference_seconds=1500.0,
        previous_queries=100,
        next_queries=200,
        max_active_seconds=3600.0,
    )
    assert not tangent_actor._next_prefix_fits(
        fixed_seconds=3600.0,
        inference_seconds=0.0,
        previous_queries=0,
        next_queries=1,
        max_active_seconds=3600.0,
    )

