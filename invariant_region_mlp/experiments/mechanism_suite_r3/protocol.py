"""Frozen round-3 task matrix and method registry."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Any

from invariant_region_mlp.experiments.mechanism_suite.mechanism_mlp import (
    BALL,
    BOX,
    CENTRE,
    F_ZERO,
    ORACLE_STATE,
    ORACLE_U,
    ORACLE_Z,
    RAW,
    SEGMENT,
    MethodSpec,
    SHRINK_FACTORS,
    shrink_method,
)
from invariant_region_mlp.experiments.mechanism_suite.round2_mlp import (
    TIGHT_CENTRE,
    TIGHT_SEGMENT,
)

from .equations import BASE_SEED, POINT_COUNT, RESULTS_ROOT, SPECS


RAW_ROOT = RESULTS_ROOT / "raw"
GATES_PATH = RESULTS_ROOT / "gates.json"
D400_PILOT_PATH = RESULTS_ROOT / "s1_d400_pilots.json"
SPAN = MethodSpec("span", "span_only", 1.0)
RAW_ANTITHETIC = MethodSpec("raw_antithetic", "raw", 1.0)

S0_PDES = ("C1", "C2-convex", "C2-cancel", "C2-flip", "C3", "LQG")
S1_PDES = ("P1", "P4", "C1", "C2-convex", "C2-flip", "C3")
S2_PDES = ("C2-convex", "C2-cancel", "C2-flip")
S3_PDES = ("P1", "C1", "C2-convex", "C3")
S4_PDES = ("C1", "C3")

S0_CONFIGS = ((3, 6), (4, 6))
S1_CONFIGS = (
    (1, 16), (1, 64), (1, 256), (1, 1024),
    (2, 8), (2, 16), (2, 32), (2, 64), (2, 128), (2, 256),
    (3, 4), (3, 6), (3, 10), (3, 16), (3, 24),
    (4, 3), (4, 4), (4, 6), (4, 8),
)
S2_CONFIGS = ((2, 8), (2, 32), (3, 6), (4, 6))
S3_CONFIGS = ((2, 32), (3, 10))
S4_CONFIGS = ((2, 32), (3, 6), (4, 6))


def primary_method(pde_id: str) -> MethodSpec:
    name = str(SPECS[pde_id]["primary"])
    return {
        "box": BOX,
        "segment": SEGMENT,
        "tight_segment": TIGHT_SEGMENT,
    }[name]


def centre_method(pde_id: str) -> MethodSpec:
    return TIGHT_CENTRE if pde_id == "P4" else CENTRE


def geometry_methods(pde_id: str) -> list[MethodSpec]:
    methods = [SEGMENT, BOX, BALL, SPAN]
    if pde_id == "P4":
        methods.insert(0, TIGHT_SEGMENT)
    return methods


def _unique(methods: Iterable[MethodSpec]) -> list[MethodSpec]:
    by_name: dict[str, MethodSpec] = {}
    for method in methods:
        prior = by_name.get(method.name)
        if prior is not None and prior.to_dict() != method.to_dict():
            raise RuntimeError(f"conflicting method definitions for {method.name}")
        by_name[method.name] = method
    return list(by_name.values())


def regular_methods(pde_id: str) -> list[MethodSpec]:
    return _unique(
        [RAW]
        + geometry_methods(pde_id)
        + [ORACLE_STATE, centre_method(pde_id), F_ZERO]
        + [shrink_method(value) for value in SHRINK_FACTORS]
    )


def s0_methods(pde_id: str) -> list[MethodSpec]:
    methods = [RAW, primary_method(pde_id), ORACLE_STATE, centre_method(pde_id)]
    if pde_id == "LQG":
        methods.append(RAW_ANTITHETIC)
    return _unique(methods)


def s2_methods(pde_id: str) -> list[MethodSpec]:
    return _unique(regular_methods(pde_id) + [ORACLE_Z, ORACLE_U])


def chunk_size(d: int) -> int:
    if d <= 20:
        return 8
    if d <= 100:
        return 4
    return 1


def make_task(
    study: str,
    pde_id: str,
    d: int,
    n: int,
    M: int,
    method: MethodSpec,
    repetition: int,
) -> dict[str, Any]:
    return {
        "study": study,
        "pde_id": pde_id,
        "d": int(d),
        "n": int(n),
        "M": int(M),
        "method": method.to_dict(),
        "repetition": int(repetition),
        "chunk_size": chunk_size(int(d)),
        "base_seed": BASE_SEED,
        "point_count": POINT_COUNT,
    }


def cross_tasks(
    study: str,
    pdes: Iterable[str],
    dimensions: Iterable[int],
    configs: Iterable[tuple[int, int]],
    method_factory: Any,
    repetitions: int,
) -> list[dict[str, Any]]:
    return [
        make_task(study, pde, d, n, M, method, repetition)
        for pde in pdes
        for d in dimensions
        for n, M in configs
        for method in method_factory(pde)
        for repetition in range(repetitions)
    ]


def s0_tasks() -> list[dict[str, Any]]:
    return cross_tasks("s0", S0_PDES, (20, 100), S0_CONFIGS, s0_methods, 3)


def s2_tasks() -> list[dict[str, Any]]:
    return cross_tasks("s2", S2_PDES, (20, 100, 400), S2_CONFIGS, s2_methods, 10)


def s1_tasks(*, include_lqg: bool = False) -> list[dict[str, Any]]:
    pdes = list(S1_PDES) + (["LQG"] if include_lqg else [])
    return cross_tasks("s1", pdes, (20, 100, 400), S1_CONFIGS, regular_methods, 10)


def s3_tasks() -> list[dict[str, Any]]:
    return cross_tasks(
        "s3", S3_PDES, (20, 50, 100, 200, 400), S3_CONFIGS,
        regular_methods, 10,
    )


def s4_tasks() -> list[dict[str, Any]]:
    return cross_tasks(
        "s4", S4_PDES, (20, 100), S4_CONFIGS, geometry_methods, 10
    )


def safe_name(value: str) -> str:
    return value.replace("-", "m").replace(".", "p").replace("+", "")


def task_path(task: dict[str, Any]) -> Path:
    """Canonical path: exact overlaps across studies have one artifact."""

    return (
        RAW_ROOT
        / str(task["pde_id"])
        / f"d{int(task['d'])}"
        / f"n{int(task['n'])}_M{int(task['M'])}"
        / safe_name(str(task["method"]["name"]))
        / f"rep{int(task['repetition']):02d}.npz"
    )


def task_identity(task: dict[str, Any]) -> tuple[Any, ...]:
    return (
        task["pde_id"], int(task["d"]), int(task["n"]), int(task["M"]),
        task["method"]["name"], int(task["repetition"]),
    )

