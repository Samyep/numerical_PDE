from __future__ import annotations

from dataclasses import replace
import hashlib
import math
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from invariant_region_mlp.experiments.original_pde_benchmark_comparison.deep_bsde_port import (
    DeepBSDEConfig,
    OfficialDeepBSDE,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.models import (
    TerminalAwareActor,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.tangent_solver import (
    antithetic_normals,
)

from .benchmarks import make_benchmark
from .aggregate import aggregate_metrics
from .deep_bsde import FixedInitialDeepBSDE
from .fh_mlp import AuditedFullHistoryMLP
from .orchestration import build_manifest, load_manifest, manifest_path, write_manifest
from .prepare import _query_bank
from .protocol import DIMENSIONS, PDES, QUERY_SEEDS, REVISION, frozen_protocol
from .tangent_actor import _value_mcv_rollout_chunk


def test_frozen_query_seeds_match_sha256_derivation() -> None:
    for pde in PDES:
        for d in DIMENSIONS:
            payload = f"{REVISION}::query_bank::{pde}::{d}".encode("utf-8")
            expected = int.from_bytes(hashlib.sha256(payload).digest()[:4], "little")
            assert QUERY_SEEDS[pde][d] == expected


def test_query_banks_are_reproducible_and_begin_at_origin() -> None:
    for pde in PDES:
        t1, x1 = _query_bank(pde, 20)
        t2, x2 = _query_bank(pde, 20)
        assert t1.shape == (32,)
        assert x1.shape == (32, 20)
        assert t1.dtype == x1.dtype == np.float64
        assert np.array_equal(t1, t2)
        assert np.array_equal(x1, x2)
        assert np.array_equal(x1[0], np.zeros(20))
        assert np.all(t1 == 0.0)


def test_frozen_protocol_is_json_native_and_stable() -> None:
    first = frozen_protocol()
    second = frozen_protocol()
    assert first == second
    assert len(first["protocol_sha256"]) == 64
    assert first["query_count"] == 32


def test_drivers_have_the_preregistered_signs() -> None:
    y = torch.tensor([0.4], dtype=torch.float64)
    z = torch.tensor([[1.0, 2.0]], dtype=torch.float64)
    hjb = make_benchmark("hjb_log", 2)
    allen = make_benchmark("allen_cahn", 2)
    dense = make_benchmark("dense_fourier", 2)
    assert torch.allclose(hjb.driver(y, z), torch.tensor([-2.5], dtype=torch.float64))
    assert torch.allclose(allen.driver(y, z), y - y**3)
    assert torch.allclose(
        dense.driver(y, z), 0.5 * torch.tanh(y) * z.square().sum(1)
    )


def test_fixed_initial_deep_bsde_uses_arbitrary_x0() -> None:
    benchmark = make_benchmark("hjb_log", 2)
    config = DeepBSDEConfig(
        benchmark="hjb_log_d2",
        total_time=1.0,
        intervals=2,
        hidden_widths=(4, 4),
        learning_rate=1e-2,
        iterations=2,
        batch_size=4,
        validation_size=4,
        logging_frequency=1,
        y_init_range=(0.0, 1.0),
    )
    x0 = torch.tensor([1.25, -0.75], dtype=torch.float64)
    model = FixedInitialDeepBSDE(benchmark, config, x0).to(dtype=torch.float64)
    model.eval()
    dw = torch.zeros((4, 2, 2), dtype=torch.float64)
    _, terminal_x = model.rollout(dw)
    assert torch.equal(terminal_x, x0.expand_as(terminal_x))


def test_zero_initial_state_matches_native_hjb_rollout() -> None:
    benchmark = make_benchmark("hjb_log", 100)
    config = DeepBSDEConfig(
        benchmark="hjb_log_d100",
        total_time=1.0,
        intervals=20,
        hidden_widths=(110, 110),
        learning_rate=1e-2,
        iterations=2,
        batch_size=4,
        validation_size=4,
        logging_frequency=1,
        y_init_range=(0.0, 1.0),
    )
    torch.manual_seed(123)
    native = OfficialDeepBSDE(benchmark, config).to(dtype=torch.float64)
    fixed = FixedInitialDeepBSDE(benchmark, config, torch.zeros(100)).to(dtype=torch.float64)
    fixed.load_state_dict(native.state_dict(), strict=False)
    native.eval()
    fixed.eval()
    generator = torch.Generator().manual_seed(456)
    dw = torch.randn((4, 20, 100), generator=generator, dtype=torch.float64) / math.sqrt(20)
    native_y, native_x = native.rollout(dw)
    fixed_y, fixed_x = fixed.rollout(dw)
    assert torch.equal(native_x, fixed_x)
    assert torch.equal(native_y, fixed_y)


class _ConstantJointDriver:
    d = 3
    T = 0.5
    sigma = math.sqrt(2.0)
    has_gradient = True
    r_low = -math.inf
    r_high = math.inf

    @staticmethod
    def terminal(x: torch.Tensor) -> torch.Tensor:
        return torch.full(x.shape[:-1], 2.5, dtype=x.dtype, device=x.device)

    @staticmethod
    def terminal_gradient_torch(x: torch.Tensor) -> torch.Tensor:
        return torch.zeros_like(x)

    @staticmethod
    def driver(value: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        assert z is not None
        return torch.zeros_like(value)


def test_fh_mlp_constant_joint_solution_is_exact() -> None:
    solver = AuditedFullHistoryMLP(
        _ConstantJointDriver(), 8, torch.device("cpu"), seed=7, dtype=torch.float64
    )
    state = solver.solve(
        3,
        torch.zeros(2, dtype=torch.float64),
        torch.tensor([[1.0, 2.0, 3.0], [-2.0, 0.5, 1.0]], dtype=torch.float64),
    )
    assert torch.equal(state[:, 0], torch.full((2,), 2.5, dtype=torch.float64))
    assert torch.equal(state[:, 1:], torch.zeros((2, 3), dtype=torch.float64))


def test_dense_fourier_joint_driver_depends_on_value_and_gradient() -> None:
    driver = make_benchmark("dense_fourier", 20)
    z1 = torch.ones((2, 20), dtype=torch.float64)
    z2 = 2.0 * z1
    y = torch.tensor([0.2, 0.7], dtype=torch.float64)
    base = driver.driver(y, z1)
    changed_y = driver.driver(y + 0.3, z1)
    changed_z = driver.driver(y, z2)
    assert not torch.equal(base, changed_y)
    assert torch.allclose(changed_z, 4.0 * base)


def test_value_actor_autograd_mcv_is_a_predictable_zero_mean_control() -> None:
    benchmark = make_benchmark("allen_cahn", 2).impl
    actor = TerminalAwareActor(benchmark, width=8, depth=2, activation="silu")
    actor.eval()
    generator = torch.Generator().manual_seed(2345)
    normals = antithetic_normals(
        1, 512, 8, 2, generator=generator, device=torch.device("cpu")
    )
    plain, corrected, _, martingale, diagnostics = _value_mcv_rollout_chunk(
        benchmark,
        actor,
        torch.zeros(1),
        torch.zeros((1, 2)),
        normals,
        use_mcv=True,
    )
    assert torch.allclose(corrected, plain - martingale)
    values = martingale.double().reshape(-1)
    standard_error = values.std(unbiased=True) / math.sqrt(len(values))
    assert abs(float(values.mean())) <= 6.0 * float(standard_error) + 1e-12
    assert diagnostics["nonfinite_count"] == 0


def test_frozen_four_gpu_task_manifests_are_complete_and_unique(
    tmp_path: Path,
) -> None:
    pilot = build_manifest("B", 100)
    assert pilot["task_count"] == 44
    assert sum(pilot["slot_task_counts"]) == 44
    assert len({row["task_id"] for row in pilot["tasks"]}) == 44

    formal = build_manifest("C", 20)
    assert formal["task_count"] == 636
    assert sum(formal["slot_task_counts"]) == 636
    assert len({row["task_id"] for row in formal["tasks"]}) == 636
    assert sum(row["method"] == "tangent_actor" for row in formal["tasks"]) == 12
    assert sum(row["method"] == "deep_bsde" for row in formal["tasks"]) == 384
    assert sum(row["method"] == "fh_mlp_standard" for row in formal["tasks"]) == 240
    assert max(formal["slot_weights"]) / min(formal["slot_weights"]) < 1.01

    path = tmp_path / "manifest.json"
    write_manifest(path, formal)
    assert load_manifest(path) == formal
    assert load_manifest(manifest_path("B", 100))["task_count"] == 44
    for dimension in DIMENSIONS:
        assert load_manifest(manifest_path("C", dimension))["task_count"] == 636


def test_rmse_32_r_uses_every_point_and_repeat() -> None:
    rows = []
    for seed in range(3):
        for point in range(32):
            rows.append({
                "method": "Deep BSDE",
                "pde": "hjb_log",
                "dimension": 20,
                "point": point,
                "repeat_id": f"s{seed}",
                "prediction": 1.0,
                "reference": 0.0,
                "error": 1.0,
            })
    scaling, pointwise = aggregate_metrics(pd.DataFrame(rows))
    assert len(pointwise) == 32
    assert scaling.loc[0, "status"] == "COMPLETE"
    assert scaling.loc[0, "observations"] == 96
    assert scaling.loc[0, "rmse_32_R"] == 1.0
    assert scaling.loc[0, "rmse_random_31_R"] == 1.0
    assert scaling.loc[0, "origin_error_rms_over_repeats"] == 1.0

