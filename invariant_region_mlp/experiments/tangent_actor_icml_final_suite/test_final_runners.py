from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import sys

import numpy as np
import torch


ROOT = Path(__file__).resolve().parents[3]
DPI = ROOT.parent / "external" / "DeepPicardIteration"
DPI_DEPS = Path("/work/nvme/bgvf/ywang151/original_pde_source_audit/dpi_deps")
for path in (DPI_DEPS, DPI):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))


def test_mature_actor_budget_is_frozen() -> None:
    from .run_joint_actor import MATURE as joint
    from .run_mlp_failure_actor import MATURE as nonsmooth
    from .run_standard_actor import MATURE as standard

    required = {
        "training_points": 8192, "paths_per_label": 128,
        "label_time_steps": 32, "optimizer_steps": 5000,
        "width": 384, "depth": 4,
    }
    for config in (joint, nonsmooth, standard):
        record = asdict(config)
        assert all(record[key] == value for key, value in required.items())


def test_official_dpi_transfer_equations_match_frozen_pdes() -> None:
    from .official_dpi_transfer import register_equations
    import picard.equations as eqs

    names = register_equations()
    generator = torch.Generator().manual_seed(17)
    x = torch.randn((7, 10), generator=generator, dtype=torch.float64)
    y = torch.rand(7, generator=generator, dtype=torch.float64)
    z = torch.randn((7, 10), generator=generator, dtype=torch.float64)
    t = torch.rand((7, 1), generator=generator, dtype=torch.float64)

    hjb = getattr(eqs, names["hjb_log"])(nx=10, T=1.0, alpha=2.0)
    assert torch.allclose(hjb.fff(t, x, y[:, None], z)[:, 0], -0.5 * z.square().sum(1))
    assert hjb.g(x[0]).shape == (1,)  # PISGradNet differentiates one sample at a time.

    concave = getattr(eqs, names["joint_concave"])(nx=10, T=0.5, alpha=2.0)
    expected = y - y.pow(3) - 0.5 * z.square().sum(1)
    assert torch.allclose(concave.fff(t, x, y[:, None], z)[:, 0], expected)

    convex = getattr(eqs, names["joint_convex"])(nx=10, T=0.5, alpha=2.0)
    expected = y.pow(3) - y + 0.5 * z.square().sum(1)
    assert torch.allclose(convex.fff(t, x, y[:, None], z)[:, 0], expected)

    nonsmooth = getattr(eqs, names["mlp_failure"])(nx=10, T=1.0, alpha=1.0)
    assert torch.allclose(
        nonsmooth.fff(t, x, y[:, None], z)[:, 0],
        torch.linalg.vector_norm(z[:, 1:], dim=1),
    )

    allen = getattr(eqs, names["allen_cahn"])(nx=10, T=0.3, alpha=2.0)
    assert torch.allclose(allen.f(t, x, y[:, None])[:, 0], y - y.pow(3))


def test_external_baseline_budgets_and_terminal_embedding() -> None:
    from .deep_bsde_transfer import _config as deep_config, make_benchmark
    from .sdgd_transfer import Config, HardTerminalNet

    deep = deep_config("joint_concave", 50, smoke=False)
    assert deep.intervals == 20 and deep.iterations == 2000
    assert deep.hidden_widths == (60, 60) and deep.batch_size == 64
    deep_allen = deep_config("allen_cahn", 100, smoke=False)
    assert deep_allen.iterations == 4000
    assert deep_allen.learning_rate == 5e-4
    assert deep_allen.y_init_range == (0.3, 0.5)

    config = Config()
    assert config.epochs == 10_000
    assert config.width == 1024 and config.hidden_layers == 3
    assert config.pgd_steps == 20
    assert config.training_dimension_samples == 100
    assert config.adversarial_dimension_samples == 10
    benchmark = make_benchmark("joint_concave", 10)
    model = HardTerminalNet(benchmark, Config(width=16, hidden_layers=2))
    x = torch.randn(5, 10)
    terminal_t = torch.full((5, 1), benchmark.T)
    assert torch.allclose(
        model(torch.cat((x, terminal_t), 1))[:, 0], benchmark.terminal(x), atol=1e-7,
    )


def test_dpi_and_pinn_profiles_are_selected_before_metrics() -> None:
    from .official_dpi_transfer import _config, _resolve_profile, register_equations

    assert _resolve_profile("ou_gmm", "auto") == "hjb_pis"
    for family in (
        "hjb_log", "joint_concave", "joint_convex",
        "nonconvex_coupled", "mlp_failure", "allen_cahn",
    ):
        assert _resolve_profile(family, "auto") == "semilinear_paper"

    names = register_equations()
    dpi, _, _, _, _, profile = _config(
        "dpi", "joint_concave", 100, 0, "official",
        names["joint_concave"], "auto",
    )
    assert profile == "semilinear_paper"
    assert dpi.PICARD.N == 80
    assert dpi.DATA.DATA_SIZE == 4096
    assert dpi.DATA.kwargs.n_estimate_terminal == 4096
    assert dpi.DATA.kwargs.n_estimate_integral == 4096
    assert dpi.TRAIN.N_EPOCHS == 16 and dpi.TRAIN.BATCH_SIZE == 512
    assert list(dpi.NETWORK.NEURONS) == [128, 128, 128, 128]
    assert list(dpi.NETWORK.ACTIVATIONS) == ["ELU"] * 4

    pinn, _, _, _, _, profile = _config(
        "pinn_hte", "joint_concave", 100, 0, "official",
        names["joint_concave"], "auto",
    )
    assert profile == "semilinear_paper"
    assert pinn.METHOD.num_v_samples == 16
    assert pinn.TRAIN.N_EPOCHS == 30_000
    assert pinn.TRAIN.LOSS.beta == 10.0
    assert list(pinn.NETWORK.NEURONS) == [128, 128, 128, 128]

    native, _, _, _, _, profile = _config(
        "dpi", "ou_gmm", 100, 0, "official", names["ou_gmm"], "auto",
    )
    assert profile == "hjb_pis"
    assert native.PICARD.N == 40
    assert list(native.NETWORK.NEURONS) == [512, 512, 512, 512]
    assert native.NETWORK.PISGRADNET is True


def test_sdgd_profiles_match_the_two_published_training_protocols() -> None:
    from .sdgd_transfer import _configuration, resolve_profile

    assert resolve_profile("hjb_log", "auto") == "hjb_adversarial"
    assert resolve_profile("ou_gmm", "auto") == "hjb_adversarial"
    for family in (
        "joint_concave", "joint_convex", "nonconvex_coupled",
        "mlp_failure", "allen_cahn",
    ):
        assert resolve_profile(family, "auto") == "semilinear_sde"

    hjb = _configuration("hjb_adversarial", smoke=False)
    assert hjb.epochs == 10_000 and hjb.collocation_points == 100
    assert hjb.pgd_steps == 20 and hjb.pgd_step_size == 0.05
    assert hjb.training_dimension_samples == 100
    assert hjb.adversarial_dimension_samples == 10

    semilinear = _configuration("semilinear_sde", smoke=False)
    assert semilinear.epochs == 10_000
    assert semilinear.trajectories == 100 and semilinear.time_steps == 20
    assert semilinear.training_dimension_samples == 10
    assert semilinear.terminal_value_weight == 20.0
    assert semilinear.terminal_gradient_weight == 1.0
    assert semilinear.exponential_lr_gamma == 0.9995


def test_sdgd_semilinear_resamples_coordinates_at_each_time_slice(monkeypatch) -> None:
    from .sdgd_transfer import (
        SemilinearSolver, _configuration, make_benchmark,
    )

    config = _configuration("semilinear_sde", smoke=True)
    benchmark = make_benchmark("joint_concave", 10)
    calls = 0
    original_choice = np.random.choice

    def counted_choice(*args, **kwargs):
        nonlocal calls
        calls += 1
        return original_choice(*args, **kwargs)

    monkeypatch.setattr(np.random, "choice", counted_choice)
    solver = SemilinearSolver(
        "joint_concave", benchmark, config, torch.device("cpu"), seed=0, smoke=True,
    )
    loss = solver.step()
    assert np.isfinite(loss)
    assert calls == 2 * config.time_steps
