"""Mathematical unit tests for the frozen component ablations."""

from __future__ import annotations

import numpy as np
import torch

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.core import (
    _weighted_rollout_chunk,
    antithetic_normals,
    tangent_labels,
)
from invariant_region_mlp.experiments.tangent_actor_icml_final_suite.run_joint_actor import (
    make_driver,
)
from invariant_region_mlp.experiments.tangent_actor_mcv_training.labels import (
    _weighted_rollout_chunk_mcv_labels,
    mcv_tangent_labels,
)
from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis.models import (
    TerminalAwareActor,
)

from .models import UnanchoredJointActor


def _actors(family: str = "joint_concave", d: int = 6):
    driver = make_driver(family, d)
    torch.manual_seed(1123)
    anchored = TerminalAwareActor(driver, width=16, depth=2)
    torch.manual_seed(1123)
    unanchored = UnanchoredJointActor(driver, width=16, depth=2)
    return driver, anchored, unanchored


def test_initial_actors_match_everywhere_and_have_same_parameters() -> None:
    driver, anchored, unanchored = _actors()
    assert anchored.state_dict().keys() == unanchored.state_dict().keys()
    for key in anchored.state_dict():
        torch.testing.assert_close(anchored.state_dict()[key], unanchored.state_dict()[key])
    generator = torch.Generator().manual_seed(77)
    t = torch.rand(19, generator=generator) * driver.T
    x = torch.randn((19, driver.d), generator=generator)
    anchored_r, anchored_q = anchored.predict(t, x, driver.T)
    unanchored_r, unanchored_q = unanchored.predict(t, x, driver.T)
    torch.testing.assert_close(anchored_r, unanchored_r, rtol=0.0, atol=2e-7)
    assert anchored_q is not None and unanchored_q is not None
    torch.testing.assert_close(anchored_q, unanchored_q, rtol=0.0, atol=2e-7)


def test_terminal_identity_is_hard_only_for_anchored_actor() -> None:
    driver, anchored, unanchored = _actors()
    with torch.no_grad():
        anchored.network[-1].bias.fill_(0.25)
        unanchored.network[-1].bias.fill_(0.25)
    x = torch.randn((17, driver.d), generator=torch.Generator().manual_seed(41))
    t = torch.full((len(x),), driver.T)
    terminal = driver.terminal(x)
    terminal_q = driver.sigma * driver.terminal_gradient_torch(x)
    anchored_r, anchored_q = anchored.predict(t, x, driver.T)
    unanchored_r, unanchored_q = unanchored.predict(t, x, driver.T)
    torch.testing.assert_close(anchored_r, terminal, rtol=0.0, atol=2e-7)
    assert anchored_q is not None and unanchored_q is not None
    torch.testing.assert_close(anchored_q, terminal_q, rtol=0.0, atol=2e-7)
    assert float(torch.max(torch.abs(unanchored_r - terminal))) > 1e-3
    assert float(torch.max(torch.abs(unanchored_q - terminal_q))) > 1e-3
    assert torch.all(unanchored_r > driver.r_low)
    assert torch.all(unanchored_r < driver.r_high)


def test_plain_and_mcv_pathwise_payoffs_have_exact_correction_identity() -> None:
    driver, actor, _ = _actors(d=4)
    t = torch.tensor([0.1, 0.3], dtype=torch.float32)
    x = torch.randn((2, driver.d), generator=torch.Generator().manual_seed(5))
    generator = torch.Generator().manual_seed(6)
    normals = antithetic_normals(
        (2, 8, 4, driver.d), generator=generator, device=torch.device("cpu")
    )
    plain, plain_z, plain_weight, _ = _weighted_rollout_chunk(
        driver, actor, t, x, normals, need_z=True
    )
    corrected, corrected_z, corrected_weight, martingale, _ = (
        _weighted_rollout_chunk_mcv_labels(
            driver, actor, t, x, normals, need_z=True
        )
    )
    torch.testing.assert_close(plain, corrected + martingale, rtol=2e-6, atol=2e-6)
    torch.testing.assert_close(plain_weight, corrected_weight, rtol=2e-6, atol=2e-6)
    assert plain_z is not None and corrected_z is not None
    assert torch.isfinite(plain_z).all()
    assert torch.isfinite(corrected_z).all()
    assert plain_z.shape == corrected_z.shape == x.shape


def test_plain_and_mcv_label_interfaces_return_full_pathwise_gradients() -> None:
    driver, actor, _ = _actors(d=4)
    rng = np.random.default_rng(19)
    t = rng.uniform(0.0, driver.T, 3)
    x = rng.standard_normal((3, driver.d))
    common = dict(
        paths=8,
        time_steps=2,
        seed=29,
        device=torch.device("cpu"),
        parent_batch=3,
        path_chunk=4,
    )
    plain_r, plain_q, plain_meta = tangent_labels(driver, actor, t, x, **common)
    mcv_r, mcv_q, mcv_meta = mcv_tangent_labels(driver, actor, t, x, **common)
    assert plain_q is not None and mcv_q is not None
    assert plain_r.shape == mcv_r.shape == (3,)
    assert plain_q.shape == mcv_q.shape == (3, driver.d)
    assert np.isfinite(plain_q).all() and np.isfinite(mcv_q).all()
    assert plain_meta["state_dependence_preserved"]
    assert mcv_meta["complete_corrected_payoff_differentiated"]


def test_convex_sign_transform_is_applied_to_value_and_gradient() -> None:
    driver = make_driver("joint_convex", 6)
    values = np.asarray([-0.2, 0.4, 1.3])
    gradients = np.arange(18, dtype=np.float64).reshape(3, 6)
    np.testing.assert_allclose(driver.working_to_original(values), -values)
    np.testing.assert_allclose(driver.working_to_original(gradients), -gradients)
