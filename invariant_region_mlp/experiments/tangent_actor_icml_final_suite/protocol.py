"""Single source of truth for the frozen ICML final-suite protocol.

This module deliberately contains definitions and compatibility decisions, not
measured results.  Result aggregation must never mutate these declarations.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path
from typing import Final


BASE_COMMIT: Final = "d36ab9e79d0e6ffe49c43d60cdb28a064dc9ac09"
BRANCH: Final = "ir-mlp-tangent-actor-icml-final-suite"
DIMENSIONS: Final = (10, 20, 50, 100, 200, 400)
MLP_FAILURE_DIMENSIONS: Final = (20, 50, 100, 200, 400, 1000)
METHODS: Final = ("tangent_actor", "pinn_hte", "sdgd_pinn", "dpi", "deep_bsde")
TRAINING_SEEDS: Final = (0,)
TEST_STATES: Final = 960

# Protocol amendment requested after Stage-B execution had started. Completed
# multi-seed artifacts remain immutable historical evidence, but every new
# PDE/dimension/method cell uses seed 0 only. The large dimension sweep is
# blocked until the independent baseline-fidelity review is approved.
PROTOCOL_AMENDMENT: Final = {
    "effective_date": "2026-10-18",
    "name": "single_seed_and_mandatory_fidelity_gate",
    "training_seeds": [0],
    "historical_multiseed_policy": "preserve, report separately, never select best seed",
    "dimension_sweep_gate": "independent_review_required_after_100d_stage_b",
}

SOURCE_COMMITS: Final = {
    "deep_bsde": "e76fed80b1995daf4a0dd6d3e2cf64a0931dabda",
    "dpi": "e8616c2f531c9f95299b246a13c1180b8f02822a",
    "pinn_hte": "e8616c2f531c9f95299b246a13c1180b8f02822a",
    "sdgd_pinn": "c597cc4fe9caac0b3ca118ea84453405f25ff833",
    "tangent_actor": BASE_COMMIT,
}


@dataclass(frozen=True)
class PDEDefinition:
    key: str
    title: str
    equation: str
    terminal: str
    horizon: float
    sigma: str
    reference: str
    test_distribution: str
    effective_dimension: int | None = None


PDES: Final = (
    PDEDefinition(
        "hjb_log", "HJB-Log",
        "u_t + Delta u - ||grad u||^2 = 0",
        "log((1+||x||^2)/2)", 1.0, "sqrt(2) I",
        "independent Hopf-Cole one-dimensional quadrature",
        "t~Uniform[0,T), x~N(0,I_d); d=100 reuses the frozen common bank",
    ),
    PDEDefinition(
        "ou_gmm", "OU/GMM HJB",
        "u_t + 0.5 Delta u - (mu-x).grad u - 0.5||grad u||^2 - d = 0",
        "negative log density of the frozen five-component GMM", 1.0, "I",
        "exact evolving Gaussian-mixture density",
        "t~Uniform[0,T), x|t~N(0,(4+t)I_d); d=100 uses official tensors/bank",
    ),
    PDEDefinition(
        "joint_concave", "Jointly concave value-gradient PDE",
        "u_t + Delta u + u-u^3-0.5||sqrt(2)grad u||^2 = 0",
        "0.5+0.15(sin(b1.x)+sin(b2.x))", 0.5, "sqrt(2) I",
        "independently refined two-dimensional periodic spectral solve",
        "t~Uniform[0,T), x~N(0,I_d)", 2,
    ),
    PDEDefinition(
        "joint_convex", "Jointly convex value-gradient PDE",
        "u_t + Delta u + u^3-u+0.5||sqrt(2)grad u||^2 = 0",
        "0.5+0.15(sin(b1.x)+sin(b2.x))", 0.5, "sqrt(2) I",
        "independently refined two-dimensional periodic spectral solve",
        "t~Uniform[0,T), x~N(0,I_d)", 2,
    ),
    PDEDefinition(
        "nonconvex_coupled", "Nonconvex-nonconcave coupled PDE",
        "u_t + Delta u + 0.2 sin(2pi u)-0.5||sqrt(2)grad u+u e||^2 = 0",
        "0.5+0.3 sin(b1.x) cos(b2.x)", 0.5, "sqrt(2) I",
        "independently refined two-dimensional Fourier solve",
        "t~Uniform[0,T), x~N(0,I_d)", 2,
    ),
    PDEDefinition(
        "mlp_failure", "Hutzenthaler-Nguyen MLP failure example",
        "u_t + 0.5 Delta u + ||grad(u)_{2:d}|| = 0",
        "|x_1|", 1.0, "I",
        "exact E|x_1+sqrt(1-t)Z|",
        "t~Uniform[0,T), x~N(0,I_d); native point x=0 reported separately", 1,
    ),
    PDEDefinition(
        "allen_cahn", "Allen-Cahn",
        "u_t + Delta u + u-u^3 = 0",
        "1/(2+0.4||x||^2)", 0.3, "sqrt(2) I",
        "dimension-specific radial PDE with mesh/time refinement and Richardson audit",
        "t~Uniform[0,T), x|t~N(0,2t I_d); d=100 reuses the frozen common bank",
    ),
)


NATIVE = {
    ("hjb_log", "deep_bsde", 100),
    ("ou_gmm", "dpi", 100),
    ("ou_gmm", "pinn_hte", 100),
    ("allen_cahn", "deep_bsde", 100),
}


def dimensions_for(pde: str) -> tuple[int, ...]:
    return MLP_FAILURE_DIMENSIONS if pde == "mlp_failure" else DIMENSIONS


def compatibility(pde: str, method: str, dimension: int) -> tuple[str, str, str]:
    """Return (classification, execution_state, reason).

    Classification follows the user-specified native/transfer/adapted/
    unsupported vocabulary.  Execution state is updated only by the manifest,
    never by changing the classification after seeing accuracy.
    """
    if (pde, method, dimension) in NATIVE:
        return "native_official", "pending_or_reuse_audit", "published PDE/dimension/configuration"
    if pde == "hjb_log" and method == "sdgd_pinn":
        return (
            "faithful_transfer", "pending",
            "paper HJB equation but 100D is below the published 250D minimum; native 250D fidelity anchor is separate",
        )
    if method == "deep_bsde":
        if pde == "mlp_failure":
            return "adapted", "pending", "native single-point BSDE with explicit minimum-norm subgradient convention"
        return "faithful_transfer", "pending", "native scalar-Y0 algorithm; reported only on the single-point task"
    if pde == "mlp_failure":
        if method == "tangent_actor":
            return "adapted", "pending", "minimum-norm subgradient b=0 at zero; no global smooth tangent theorem"
        return (
            "adapted",
            "pending",
            "source-compatible minimum-norm subgradient preflight passed; no smooth-driver theory claim",
        )
    if method == "tangent_actor":
        if pde == "joint_convex":
            return "adapted", "pending", "mathematically equivalent sign transformation w=-u, used only by Tangent-Actor"
        return "faithful_transfer", "pending_or_reuse_audit", "frozen final Tangent-Actor configuration"
    return "faithful_transfer", "pending", "original method/budget with only PDE and dimension tensors changed"


def stable_seed(*parts: object) -> int:
    digest = sha256("::".join(map(str, parts)).encode()).digest()
    return int.from_bytes(digest[:4], "little")


def frozen_protocol() -> dict[str, object]:
    payload = {
        "base_commit": BASE_COMMIT,
        "branch": BRANCH,
        "dimensions": list(DIMENSIONS),
        "mlp_failure_dimensions": list(MLP_FAILURE_DIMENSIONS),
        "methods": list(METHODS),
        "training_seeds": list(TRAINING_SEEDS),
        "protocol_amendment": PROTOCOL_AMENDMENT,
        "test_states_default": TEST_STATES,
        "pdes": [asdict(item) for item in PDES],
        "source_commits": SOURCE_COMMITS,
        "selection": {
            "tangent_actor": "reference-free held-out tangent labels and paired policy-value consistency",
            "external_baselines": "published fixed budget; no final-test checkpoint selection",
            "euler_grid": [64, 128, 256, 512],
            "euler_role": "predeclared post-freeze numerical refinement audit",
        },
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    payload["protocol_sha256"] = sha256(canonical).hexdigest()
    return payload


PACKAGE_DIR = Path(__file__).resolve().parent
SOURCE_ROOT = PACKAGE_DIR.parents[2]
RESULTS_DIR = SOURCE_ROOT / "invariant_region_mlp" / "results" / "tangent_actor_icml_final_suite"

