from __future__ import annotations

from .protocol import METHODS, PDES, compatibility, dimensions_for, frozen_protocol


def test_grid_has_42_pde_dimension_configurations() -> None:
    assert sum(len(dimensions_for(pde.key)) for pde in PDES) == 42


def test_grid_has_210_method_cells() -> None:
    cells = [
        compatibility(pde.key, method, d)
        for pde in PDES for d in dimensions_for(pde.key) for method in METHODS
    ]
    assert len(cells) == 210


def test_native_anchors_are_frozen() -> None:
    assert compatibility("hjb_log", "sdgd_pinn", 100)[0] == "native_official"
    assert compatibility("ou_gmm", "dpi", 100)[0] == "native_official"
    assert compatibility("ou_gmm", "pinn_hte", 100)[0] == "native_official"
    assert compatibility("allen_cahn", "deep_bsde", 100)[0] == "native_official"


def test_nonsmooth_problem_is_never_silently_called_faithful() -> None:
    for method in METHODS:
        assert compatibility("mlp_failure", method, 100)[0] == "adapted"


def test_protocol_hash_is_stable_length() -> None:
    assert len(frozen_protocol()["protocol_sha256"]) == 64
