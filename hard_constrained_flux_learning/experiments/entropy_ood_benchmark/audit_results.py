"""Machine-checkable consistency audit for the entropy/OOD artifacts."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess

import numpy as np
import pandas as pd
import torch

import swe_radial_common as C


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def git_head(path: Path) -> str:
    return subprocess.check_output(
        ["git", "-C", str(path), "rev-parse", "HEAD"], text=True
    ).strip()


def algebra_tests() -> dict[str, float]:
    state = torch.zeros(2, C.N_COARSE, C.N_COARSE, 3)
    state[..., 0] = 1.3
    state[..., 1] = 0.2
    state[..., 2] = -0.1
    fluxes = C.hll_fluxes(state)
    uniform_error = float((C.flux_divergence(state, *fluxes, 0.01) - state).abs().max())
    require(uniform_error == 0.0, "Uniform SWE state is not preserved")

    generator = torch.Generator().manual_seed(991)
    h_left = 0.5 + torch.rand(3, 7, generator=generator)
    h_right = 0.5 + torch.rand(3, 7, generator=generator)
    left = torch.stack(
        (
            h_left,
            h_left * (torch.rand(3, 7, generator=generator) - 0.5),
            h_left * (torch.rand(3, 7, generator=generator) - 0.5),
        ),
        dim=-1,
    )
    right = torch.stack(
        (
            h_right,
            h_right * (torch.rand(3, 7, generator=generator) - 0.5),
            h_right * (torch.rand(3, 7, generator=generator) - 0.5),
        ),
        dim=-1,
    )
    matrix, alpha, _ = C.roe_waves_oriented(left, right)
    reconstructed = torch.einsum("...ij,...j->...i", matrix, alpha)
    roe_error = float((reconstructed - (right - left)).abs().max())
    require(roe_error < 2.0e-6, "Roe characteristic reconstruction failed")

    cells = 12
    random_state = torch.zeros(2, 4, cells, 3)
    random_state[..., 0] = 0.7 + torch.rand(2, 4, cells, generator=generator)
    random_state[..., 1:] = 0.4 * (
        torch.rand(2, 4, cells, 2, generator=generator) - 0.5
    )
    proposal = C.hll_faces_oriented(random_state)
    proposal[..., 1:-1, :] += 0.5 * torch.randn(
        proposal[..., 1:-1, :].shape, generator=generator
    )
    projected = C.strict_entropy_projection_oriented(proposal, random_state)
    maximum_residual = float(
        C.interface_entropy_residual_oriented(
            projected.double(), random_state.double()
        ).max()
    )
    require(maximum_residual <= 1.0e-11, "Strict Tadmor projection is infeasible")
    return {
        "uniform_state_max_error": uniform_error,
        "roe_reconstruction_max_error": roe_error,
        "projection_max_residual": maximum_residual,
    }


def main() -> None:
    protocol = read_json(HERE / "protocol.json")
    require(
        protocol["status"] == "exploratory_design_record"
        and not protocol["immutable_preregistration"],
        "Protocol scope/chronology label changed",
    )
    require(protocol["failed_trajectories_are_retained"], "Failure accounting changed")

    external = HERE.parents[4] / ".hcfl_third_party_audit"
    commits = {
        "clawNO": git_head(external / "clawNO"),
        "RoeNet": git_head(external / "RoeNet"),
        "PINN": git_head(external / "Datafree_PINN_Compressible"),
        "PDEBench": git_head(external / "PDEBench"),
    }
    require(
        commits["clawNO"] == protocol["swe_2d"]["official_clawfno_commit"],
        "clawNO revision drift",
    )
    require(commits["RoeNet"] == protocol["official_roenet_commit"], "RoeNet revision drift")

    euler = read_json(RESULTS / "frozen64" / "frozen64_summary.json")
    require(len(euler["rows"]) == 5 * 7, "Euler row count changed")
    require(
        all(euler["aggregate"][f"hcfl64_seed{seed}"]["completed_cases"] == 5 for seed in (0, 1, 2)),
        "An HCFL Euler failure was hidden",
    )
    require(euler["aggregate"]["fno64"]["completed_cases"] == 1, "Euler FNO failure count changed")
    require(euler["aggregate"]["roenet64"]["completed_cases"] == 3, "RoeNet failure count changed")
    for row in euler["rows"]:
        if not row["completed"]:
            require(
                row["maximum_saved_snapshot_entropy_increase"] is None,
                "A failed Euler state was assigned thermodynamic entropy",
            )
    roenet_report = read_json(RESULTS / "roenet" / "roenet_adaptation64_report_seed0.json")
    require(roenet_report["best_update"] == 2000, "Corrected RoeNet winner changed")
    require(
        roenet_report["checkpoint_selection_rule"]
        == "lexicographic(completed_trajectories, -completed_only_nrmse)",
        "RoeNet checkpoint selection is not strict completion-first",
    )

    pinn = read_json(RESULTS / "official_pinn" / "official_lnn2_native_summary.json")
    require(len(pinn["rows"]) == 6, "PINN native table must contain exactly six rows")
    require(not pinn["is_same_task_amortized_comparison"], "PINN scope label became misleading")

    dataset = torch.load(RESULTS / "swe_radial" / "swe_radial_dataset.pt", weights_only=False)
    for name, spec in C.SPLIT_SPECS.items():
        require(
            dataset[name]["trajectory"].shape[0] == spec["count"],
            f"Wrong SWE count for {name}",
        )
        require(
            dataset[name]["trajectory"].shape[1] == spec["states"],
            f"Wrong SWE horizon for {name}",
        )

    swe_frame = pd.read_csv(RESULTS / "swe_radial" / "swe_radial_case_metrics.csv")
    require(len(swe_frame) == 4 * (20 + 20 + 20 + 10), "SWE failure denominator changed")
    for split, total in (("test_id", 20), ("test_radius_ood", 20), ("test_height_ood", 20), ("test_long", 10)):
        subset = swe_frame[swe_frame.split == split]
        for method in ("HLL-32", "HCFL-s6", "FNO", "clawFNO"):
            rows = subset[subset.method == method]
            require(len(rows) == total, f"Missing {split}/{method} rows")
        require(subset[subset.method == "HCFL-s6"].completed.all(), f"HCFL failed {split}")
        require(subset[subset.method == "FNO"].completed.all(), f"FNO completion changed on {split}")
        claw = subset[subset.method == "clawFNO"]
        require(not claw.completed.any(), f"clawFNO failure count changed on {split}")
        require(claw.rollout_nrmse.isna().all(), "Failed clawFNO rows were assigned finite error")
        require(
            claw.maximum_entropy_balance.isna().all(),
            "Failed clawFNO rows were assigned thermodynamic entropy",
        )

    swe = read_json(RESULTS / "swe_radial" / "swe_radial_summary.json")
    require(
        swe["method_qualification"]["FNO"]["status"] == "not_physics_qualified"
        and swe["method_qualification"]["FNO"]["lower_nrmse_is_not_solver_success"],
        "FNO raw error was incorrectly promoted to physical solver success",
    )
    require(
        swe["aggregate"]["test_height_ood"]["FNO"][
            "mean_centerline_curvature_ratio_completed"
        ]
        > 2.0,
        "FNO strong-height ringing diagnostic changed unexpectedly",
    )
    safety_maxima = []
    for split, values in swe["exact_flux_method_safety"].items():
        hcfl = values["HCFL-s6"]
        require(hcfl["maximum_interface_residual"] <= 1.0e-11, f"Interface entropy failure on {split}")
        require(hcfl["maximum_entropy_balance"] <= C.ENTROPY_TOL + 1.0e-12, f"Fully discrete entropy failure on {split}")
        require(hcfl["maximum_conservation_closure"] <= 1.0e-6, f"Conservation failure on {split}")
        require(hcfl["depth_limiter_active"] == 0.0, f"Unexpected depth fallback count on {split}")
        safety_maxima.append(hcfl["maximum_entropy_balance"])

    ledger = read_json(RESULTS / "swe_radial" / "TUNING_LEDGER.json")
    require(ledger["selection_used_validation_only"], "Validation-only selection flag is false")
    require(not ledger["test_sets_used_for_tuning"], "Test-set tuning flag is true")
    require(ledger["selected_hcfl_stencil"] == 6, "Selected stencil differs from validation winner")
    require(
        ledger["operator_reports"]["fno"]["best_validation_completion"] == 1.0,
        "FNO validation completion changed",
    )
    require(
        ledger["operator_reports"]["claw"]["best_validation_completion"] == 0.0
        and not ledger["operator_reports"]["claw"][
            "has_fully_admissible_validation_checkpoint"
        ],
        "clawFNO must remain explicitly marked as having no admissible checkpoint",
    )

    refinement = read_json(RESULTS / "swe_radial" / "reference_grid_audit.json")
    errors = refinement["primitive_nrmse"]
    require(errors["128_vs_256_rollout"] < errors["64_vs_128_rollout"], "Reference is not grid-converging")

    report = {
        "status": "PASS",
        "algebra": algebra_tests(),
        "official_commits": commits,
        "euler_rows": len(euler["rows"]),
        "swe_rows": len(swe_frame),
        "fno_physics_qualification": swe["method_qualification"]["FNO"]["status"],
        "fno_lower_nrmse_is_not_solver_success": True,
        "maximum_exact_hcfl_entropy_balance": max(safety_maxima),
        "fully_discrete_tolerance": C.ENTROPY_TOL,
        "reference_grid_audit": errors,
    }
    output = RESULTS / "AUDIT.json"
    output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
