"""Fast structural tests for the frozen joint addendum."""

from dataclasses import asdict
import json

import numpy as np

from . import deep_bsde
from .benchmark import frozen_basis, make_benchmark
from .protocol import (
    BANK_SIZES, DB_CANDIDATE_UPDATES, DIMENSION, FH_CONFIGS,
    JOINT_ACTOR_CONFIG, REVISION, TA_PATHS, TARGETS, canonical_sha256,
    stable_seed,
)


def test_frozen_protocol() -> None:
    assert REVISION == "single_point_break_even_joint_v1"
    assert DIMENSION == 100
    assert TARGETS == {"epsilon_1": 6e-2, "epsilon_2": 2e-2}
    assert BANK_SIZES["timing"] == 1_048_576
    assert TA_PATHS == (32, 128, 512)
    assert FH_CONFIGS[0] == (2, 8) and FH_CONFIGS[-1] == (3, 24)
    assert DB_CANDIDATE_UPDATES == (500, 1000, 2000, 4000)
    assert asdict(JOINT_ACTOR_CONFIG)["optimizer_steps"] == 5000


def test_seed_namespace() -> None:
    assert stable_seed("bank", "calibration") == stable_seed("bank", "calibration")
    assert stable_seed("bank", "calibration") != stable_seed("bank", "final")


def test_benchmark_identity() -> None:
    benchmark = make_benchmark(100)
    basis = frozen_basis(100)
    assert np.max(np.abs(basis.T @ basis - np.eye(2))) < 1e-12
    x = np.zeros((3, 100))
    assert np.allclose(benchmark.terminal_numpy(x), 0.65)
    assert benchmark.has_gradient
    assert np.linalg.norm(benchmark.v) > 0.999999


def test_deep_configuration() -> None:
    config = deep_bsde.CONFIG
    assert config.total_time == 0.5
    assert config.intervals == 20
    assert config.hidden_widths == (110, 110)
    assert config.learning_rate == 1e-2
    assert config.batch_size == 64
    assert config.dtype == "float64"


def test_canonical_hash() -> None:
    assert canonical_sha256({"a": 1, "b": 2}) == canonical_sha256({"b": 2, "a": 1})
    assert json.dumps(TARGETS, allow_nan=False)

