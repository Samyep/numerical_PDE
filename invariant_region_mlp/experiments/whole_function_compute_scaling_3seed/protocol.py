"""Frozen declarations and run-spec construction for the three-seed study."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import json
from typing import Any, Final


REVISION: Final = "whole_function_compute_scaling_3seed_v1"
BRANCH: Final = "ir-mlp-compute-scaling-3seed"
PREREGISTRATION_COMMIT: Final = "51e0c5262d09fccf5bb379cc74eaf3fa025e225e"
BASE_REVISION: Final = "1245eb98edefd94f5bd0c55b886df229e2257011"
UPSTREAM_DPI_REVISION: Final = "e8616c2f531c9f95299b246a13c1180b8f02822a"

PDES: Final = ("b6_fullrank", "dense_fourier")
METHODS: Final = ("tangent_actor", "dpi", "pinn_hte")
SWEEPS: Final = ("C0", "C1", "C2")
TRAINING_SEEDS: Final = (0, 1, 2)
DIMENSION: Final = 100
TEST_STATES: Final = 960
CALIBRATION_SEED: Final = 314159
FORMAL_PROJECTED_CAP_SECONDS: Final = 82 * 3600
GLOBAL_GPU_SECONDS_CAP: Final = 90 * 3600
FORMAL_JOB_SAFETY_SECONDS: Final = 45 * 60

TA_EVALUATION_PATHS: Final = 1024
TA_EVALUATION_STEPS: Final = 64
TA_INFERENCE_REPLICATES: Final = 3

EXPECTED_INPUTS: Final = {
    "b6_fullrank": {
        "test_input_sha256": "0ee05fff7ce16f30d0d31563f60df7db3eb23a9c46bbe7bc9be8eb3dbee1311b",
        "test_bank_file_sha256": "6ee6c56ca3d4f413764f28865cb8c00823ea2881ebc90f2643ff4a36fe878730",
        "truth_value_sha256": "6262db7f6b7ee2c41cb6011dadb7288226157533488f708be7d6ea6773f2e13a",
        "parameter_file_sha256": "21b508c30ad63de912979cc9c85ecda0d81aedab1357c88a40e04692495972a4",
    },
    "dense_fourier": {
        "test_input_sha256": "d409296128f48c5ae00c5de2c6ddd998768c01dff18a89c385245a50ff581b70",
        "test_bank_file_sha256": "998fde8b2c0556de4b8fa9a55dec5e489224828cf74e5f3b36f77cb1a8e2c011",
        "truth_value_sha256": "3b33ea36e21595b5c20e5d91097d1fadfcbf3f7884c0328f20579e69995fc816",
        "parameter_file_sha256": "042588d7c3cc8175cc7b9344284725443650837f72008806c668d35b9a005307",
    },
}

PARAMETER_NAMES: Final = {
    ("tangent_actor", "C0"): "actor_updates",
    ("tangent_actor", "C1"): "paths_per_label",
    ("tangent_actor", "C2"): "optimizer_steps",
    ("dpi", "C0"): "picard_iterations",
    ("dpi", "C1"): "integral_mc_samples",
    ("dpi", "C2"): "epochs_per_picard",
    ("pinn_hte", "C0"): "training_epochs",
    ("pinn_hte", "C1"): "hutchinson_vectors",
    ("pinn_hte", "C2"): "batch_size",
}

CANDIDATE_LADDERS: Final = {
    ("tangent_actor", "C0"): [1, 2, 3, 4, 5, 6, 8, 10, 12, 16],
    ("tangent_actor", "C1"): [32, 64, 128, 192, 256, 384, 512, 768, 1024],
    ("tangent_actor", "C2"): [1250, 2500, 5000, 7500, 10000, 15000, 20000, 30000, 40000],
    ("dpi", "C0"): [20, 40, 80, 160, 320, 640, 1280],
    ("dpi", "C1"): [512, 1024, 2048, 4096, 8192, 16384, 32768, 65536],
    ("dpi", "C2"): [4, 8, 16, 32, 64, 128, 256],
    ("pinn_hte", "C0"): [10000, 30000, 60000, 120000, 240000, 480000],
    ("pinn_hte", "C1"): [4, 8, 16, 32, 64, 128, 256],
    ("pinn_hte", "C2"): [128, 256, 512, 1024, 2048, 4096, 8192, 16384],
}

BASELINE_SETTINGS: Final = {
    "tangent_actor": {
        "actor_updates": 4,
        "paths_per_label": 128,
        "optimizer_steps": 5000,
    },
    "dpi": {
        "picard_iterations": 80,
        "integral_mc_samples": 4096,
        "epochs_per_picard": 16,
    },
    "pinn_hte": {
        "training_epochs": 30000,
        "hutchinson_vectors": 16,
        "batch_size": 512,
    },
}

# Each tuple is (varied parameter, shortened calibration work amount,
# multiplier used to predict the fixed formal work amount).
CALIBRATION_ANCHORS: Final = {
    ("tangent_actor", "C0"): [(2, 2, 1.0)],
    ("tangent_actor", "C1"): [(64, 1, 4.0), (256, 1, 4.0)],
    ("tangent_actor", "C2"): [(2500, 1, 4.0), (10000, 1, 4.0)],
    ("dpi", "C0"): [(20, 20, 1.0)],
    ("dpi", "C1"): [(1024, 10, 8.0), (8192, 10, 8.0)],
    ("dpi", "C2"): [(8, 10, 8.0), (32, 10, 8.0)],
    ("pinn_hte", "C0"): [(10000, 10000, 1.0)],
    ("pinn_hte", "C1"): [(8, 5000, 6.0), (32, 5000, 6.0)],
    ("pinn_hte", "C2"): [(256, 5000, 6.0), (1024, 5000, 6.0)],
}

PACKAGE_DIR = Path(__file__).resolve().parent
SOURCE_ROOT = PACKAGE_DIR.parents[2]
RESULTS_DIR = SOURCE_ROOT / "invariant_region_mlp" / "results" / "whole_function_compute_scaling_3seed"
CALIBRATION_DIR = RESULTS_DIR / "calibration"
RUNS_DIR = RESULTS_DIR / "runs"
STATUS_DIR = RESULTS_DIR / "status"
FIGURE_DIR = RESULTS_DIR / "figures"
CALIBRATION_MANIFEST = RESULTS_DIR / "CALIBRATION_MANIFEST.json"
RUN_MANIFEST = RESULTS_DIR / "run_manifest.json"
REPORT_PATH = SOURCE_ROOT / "invariant_region_mlp" / "docs" / "WHOLE_FUNCTION_COMPUTE_SCALING_3SEED_REPORT.md"


@dataclass(frozen=True)
class RunSpec:
    phase: str
    pde: str
    method: str
    sweep: str
    parameter_name: str
    parameter_value: int
    seed: int
    settings: dict[str, int]
    checkpoint_values: tuple[int, ...]
    calibration_scale: float = 1.0

    @property
    def run_id(self) -> str:
        prefix = "cal" if self.phase == "calibration" else "formal"
        method = "ta" if self.method == "tangent_actor" else self.method
        return (
            f"{prefix}_{method}_{self.pde}_{self.sweep.lower()}_"
            f"{self.parameter_name}-{self.parameter_value}_s{self.seed}"
        )

    def payload(self) -> dict[str, Any]:
        return {**asdict(self), "run_id": self.run_id}


def settings_for(method: str, sweep: str, value: int) -> dict[str, int]:
    settings = dict(BASELINE_SETTINGS[method])
    settings[PARAMETER_NAMES[(method, sweep)]] = int(value)
    return settings


def calibration_specs() -> list[RunSpec]:
    specs: list[RunSpec] = []
    for pde in PDES:
        for method in METHODS:
            for sweep in SWEEPS:
                for value, shortened_work, scale in CALIBRATION_ANCHORS[(method, sweep)]:
                    settings = settings_for(method, sweep, value)
                    if method == "tangent_actor":
                        settings["actor_updates"] = shortened_work
                        checkpoints = (shortened_work,)
                    elif method == "dpi":
                        settings["picard_iterations"] = shortened_work
                        checkpoints = (shortened_work,)
                    else:
                        settings["training_epochs"] = shortened_work
                        checkpoints = (shortened_work,)
                    specs.append(RunSpec(
                        phase="calibration", pde=pde, method=method, sweep=sweep,
                        parameter_name=PARAMETER_NAMES[(method, sweep)],
                        parameter_value=value, seed=CALIBRATION_SEED,
                        settings=settings, checkpoint_values=checkpoints,
                        calibration_scale=scale,
                    ))
    return specs


def load_calibration_manifest() -> dict[str, Any]:
    if not CALIBRATION_MANIFEST.is_file():
        raise FileNotFoundError(f"calibration grid is not frozen: {CALIBRATION_MANIFEST}")
    return json.loads(CALIBRATION_MANIFEST.read_text(encoding="utf-8"))


def formal_specs() -> list[RunSpec]:
    manifest = load_calibration_manifest()
    specs: list[RunSpec] = []
    for line in manifest["selected_lines"]:
        pde, method, sweep = line["pde"], line["method"], line["sweep"]
        values = tuple(int(v) for v in line["parameter_values"])
        if sweep == "C0":
            final = max(values)
            for seed in TRAINING_SEEDS:
                specs.append(RunSpec(
                    phase="formal", pde=pde, method=method, sweep=sweep,
                    parameter_name=PARAMETER_NAMES[(method, sweep)],
                    parameter_value=final, seed=seed,
                    settings=settings_for(method, sweep, final),
                    checkpoint_values=values,
                ))
        else:
            for value in values:
                for seed in TRAINING_SEEDS:
                    specs.append(RunSpec(
                        phase="formal", pde=pde, method=method, sweep=sweep,
                        parameter_name=PARAMETER_NAMES[(method, sweep)],
                        parameter_value=value, seed=seed,
                        settings=settings_for(method, sweep, value),
                        checkpoint_values=(value,),
                    ))
    return specs


def ensure_directories() -> None:
    for path in (RESULTS_DIR, CALIBRATION_DIR, RUNS_DIR, STATUS_DIR, FIGURE_DIR):
        path.mkdir(parents=True, exist_ok=True)


__all__ = [
    name for name in globals()
    if name.isupper()
] + [
    "RunSpec", "calibration_specs", "ensure_directories", "formal_specs",
    "load_calibration_manifest", "settings_for",
]
