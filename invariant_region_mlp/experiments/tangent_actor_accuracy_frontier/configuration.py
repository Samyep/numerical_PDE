"""Frozen experiment specifications shared by training and post-freeze evaluation."""

from __future__ import annotations

from dataclasses import dataclass, fields
import json
from typing import Any

from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis.training import (
    FitConfig,
)

from .common import PACKAGE_DIR
from .drivers import SineQuadraticDriver


@dataclass(frozen=True)
class ExperimentSpec:
    family: str
    d: int
    lam: float
    eta: float
    seed: int
    budget: str
    estimator: str = "mcv"

    @property
    def coupled(self) -> bool:
        return self.family == "coupled"

    @property
    def run_id(self) -> str:
        lam_tag = f"{self.lam:g}".replace(".", "p")
        eta_tag = f"{self.eta:g}".replace(".", "p")
        return (
            f"{self.family}_d{self.d}_k2_lam{lam_tag}_eta{eta_tag}_"
            f"{self.budget}_{self.estimator}_s{self.seed}"
        )

    def driver(self) -> SineQuadraticDriver:
        return SineQuadraticDriver(
            d=self.d, lam=self.lam, eta=self.eta if self.coupled else 0.0
        )


def protocol() -> dict[str, Any]:
    return json.loads((PACKAGE_DIR / "protocol.json").read_text(encoding="utf-8"))


def fit_config(payload: dict[str, Any], budget: str) -> FitConfig:
    names = {field.name for field in fields(FitConfig)}
    return FitConfig(
        **{key: value for key, value in payload[budget].items() if key in names}
    )


def updates(payload: dict[str, Any], budget: str) -> int:
    return int(payload[budget]["updates"])


__all__ = ["ExperimentSpec", "fit_config", "protocol", "updates"]
