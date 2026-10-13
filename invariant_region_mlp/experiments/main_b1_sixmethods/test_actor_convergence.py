"""Protocol checks for the progressive actor convergence study."""

from __future__ import annotations

import numpy as np

from .actor_convergence import (
    BUDGET_LEVELS,
    audit_record,
    budget_protocol,
    convergence_decision,
)


def test_budget_label_work_and_data_are_strictly_increasing() -> None:
    work = [
        (level.M**level.n) * level.replicas * (
            level.training_points + level.validation_points
        )
        for level in BUDGET_LEVELS
    ]
    assert all(
        right > left for left, right in zip(work[:-1], work[1:], strict=True)
    )
    assert all(
        right.training_points > left.training_points
        for left, right in zip(
            BUDGET_LEVELS[:-1], BUDGET_LEVELS[1:], strict=True
        )
    )
    protocol = budget_protocol()
    assert protocol["nonterminal_oracle_reference_used_for_training"] is False
    assert protocol["nonterminal_oracle_reference_used_for_selection"] is False


def test_audit_noise_floor_formula_and_prediction_change() -> None:
    bank_a = np.array([[2.0, 0.0], [0.0, 2.0]])
    bank_b = np.zeros_like(bank_a)
    consensus = 0.5 * (bank_a + bank_b)
    first = audit_record(consensus, bank_a, bank_b, previous_prediction=None)
    assert first["audit_target_rmse"] == 0.0
    assert first["audit_consensus_noise_floor_rmse"] == 1.0
    assert first["audit_noise_corrected_rmse"] == 0.0
    assert first["prediction_relative_change"] is None
    second = audit_record(
        consensus,
        bank_a,
        bank_b,
        previous_prediction=np.zeros_like(consensus),
    )
    assert second["prediction_relative_change"] == 1.0


def test_convergence_requires_two_label_only_stable_transitions() -> None:
    records = [
        {
            "optimizer_converged": True,
            "audit_target_rmse": 2.0,
            "audit_consensus_noise_floor_rmse": 0.1,
            "audit_noise_corrected_rmse": 1.00,
            "prediction_relative_change": None,
            "oracle_error_that_must_be_ignored": 1000.0,
        },
        {
            "optimizer_converged": True,
            "audit_target_rmse": 1.9,
            "audit_consensus_noise_floor_rmse": 0.1,
            "audit_noise_corrected_rmse": 0.98,
            "prediction_relative_change": 0.02,
            "oracle_error_that_must_be_ignored": 0.0,
        },
        {
            "optimizer_converged": True,
            "audit_target_rmse": 1.8,
            "audit_consensus_noise_floor_rmse": 0.1,
            "audit_noise_corrected_rmse": 0.97,
            "prediction_relative_change": 0.03,
            "oracle_error_that_must_be_ignored": 500.0,
        },
    ]
    assert convergence_decision(records[:2])["converged"] is False
    decision = convergence_decision(records)
    assert decision["converged"] is True
    assert decision["reason"] == "two_consecutive_budget_plateaus"

