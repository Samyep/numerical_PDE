"""Tests for the post-training B6/Dense relative-MAE derivation."""

from __future__ import annotations

import numpy as np
import pytest

from .derive_b6_dense_rmae import relative_mae


def test_relative_mae_is_relative_l1_error() -> None:
    prediction = np.asarray([2.0, -1.0])
    reference = np.asarray([1.0, -3.0])
    assert relative_mae(prediction, reference) == 0.75


def test_relative_mae_rejects_zero_reference_scale() -> None:
    with pytest.raises(ValueError, match="zero L1 reference"):
        relative_mae(np.ones(2), np.zeros(2))
