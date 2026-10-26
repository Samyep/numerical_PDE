"""Truth-free, explicitly non-result preflight for all three method paths."""

from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
import time

import torch

from .baseline_frontier import _configuration, _fresh_runner, _seed_all
from .common import file_sha256, synchronize, write_json_atomic
from .protocol import CONFIGURATIONS, PDES, RAW_DIR, REVISION, WORK_DIR, ensure_directories
from .ta_frontier import _configs as ta_configs, _initial_and_benchmark


def _ta_smoke(pde: str, device: torch.device) -> list[dict]:
    rows = []
    for configuration in CONFIGURATIONS:
        config = ta_configs(pde, configuration)
        benchmark, actor = _initial_and_benchmark(pde, config, device)
        t = torch.zeros(2, dtype=torch.float64, device=device)
        x = torch.zeros(2, 100, dtype=torch.float64, device=device)
        with torch.no_grad():
            output = actor.predict(t, x)
        tensors = [item for item in output if item is not None]
        finite = all(bool(torch.isfinite(item).all()) for item in tensors)
        rows.append({
            "method": "tangent_actor", "pde": pde,
            "configuration_id": configuration, "finite": finite,
            "output_shapes": [list(item.shape) for item in tensors], "nonresult": True,
        })
        if not finite:
            raise FloatingPointError(rows[-1])
    return rows


def _dpi_smoke(pde: str, device: torch.device) -> dict:
    del device
    cfg, source, _overrides, _profile = _configuration("dpi", pde, "C2")
    cfg = cfg.clone()
    cfg.defrost()
    smoke_dir = WORK_DIR / "smoke" / f"dpi_c2_{pde}_{time.time_ns()}"
    cfg.NAME = str(smoke_dir)
    cfg.PICARD.N = 1
    cfg.DATA.DATA_SIZE = 32
    cfg.DATA.N_BUFFER = 1
    cfg.DATA.NEW_SAMPLING = False
    cfg.DATA.kwargs.n_estimate_terminal = 8
    # Keep C2's full 8192 integral paths: this is the proposed change under test.
    cfg.TRAIN.N_EPOCHS = 1
    cfg.TRAIN.BATCH_SIZE = 32
    cfg.freeze()
    runner = _fresh_runner(cfg, smoke_dir / "unused.pt")
    started = time.perf_counter()
    ok = runner.run_one()
    synchronize(torch.device("cuda" if torch.cuda.is_available() else "cpu"))
    checkpoint = runner.get_model_checkpoint_path(runner.exp_dir, 1)
    return {
        "method": "dpi", "pde": pde, "configuration_id": "C2",
        "status": "PASS" if ok and checkpoint.is_file() else "FAIL",
        "seconds": time.perf_counter() - started,
        "integral_mc_samples": int(cfg.DATA.kwargs.n_estimate_integral),
        "checkpoint_sha256": file_sha256(checkpoint) if checkpoint.is_file() else None,
        "source_config_sha256": file_sha256(source), "nonresult": True,
    }


def _pinn_smoke(pde: str, device: torch.device) -> dict:
    from picard.utils import hutchinson_trace_estimation_batch

    cfg, source, _overrides, _profile = _configuration("pinn_hte", pde, "C2")
    cfg = cfg.clone()
    cfg.defrost()
    smoke_dir = WORK_DIR / "smoke" / f"pinn_c2_{pde}_{time.time_ns()}"
    cfg.NAME = str(smoke_dir)
    cfg.freeze()
    runner = _fresh_runner(cfg, smoke_dir / "unused.pt")
    model = runner.get_solution().to(device)
    equation = runner.equation.to(device=device)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    batch = int(cfg.TRAIN.BATCH_SIZE)
    t = equation.T * torch.rand(batch, 1, device=device)
    x = equation.sample_x(t)
    t.requires_grad_(); x.requires_grad_()
    u = model(torch.cat((t, x), dim=-1))[:, :1]
    u_t = torch.autograd.grad(u.sum(), t, create_graph=True, retain_graph=True)[0]
    u_x = torch.autograd.grad(u.sum(), x, create_graph=True, retain_graph=True)[0]
    laplacian = hutchinson_trace_estimation_batch(
        model, t, x, int(cfg.METHOD.num_v_samples)
    )
    loss_in = equation.pinn_function(t, x, u, u_t, u_x, laplacian).square().mean()
    terminal_t = equation.T * torch.ones_like(t)
    terminal_x = equation.sample_x(terminal_t)
    terminal_u = model(torch.cat((terminal_t, terminal_x), dim=-1))[:, :1]
    loss = loss_in + 10.0 * (terminal_u - equation.g(terminal_x)).square().mean()
    optimizer.zero_grad(); loss.backward(); optimizer.step(); synchronize(device)
    return {
        "method": "pinn_hte", "pde": pde, "configuration_id": "C2",
        "status": "PASS" if bool(torch.isfinite(loss)) else "FAIL",
        "batch_size": batch, "hutchinson_vectors": int(cfg.METHOD.num_v_samples),
        "loss": float(loss.detach().cpu()),
        "source_config_sha256": file_sha256(source), "nonresult": True,
    }


def run(device: torch.device) -> dict:
    ensure_directories()
    _seed_all(0)
    rows: list[dict] = []
    for pde in PDES:
        rows.extend(_ta_smoke(pde, device))
        rows.append(_dpi_smoke(pde, device))
        rows.append(_pinn_smoke(pde, device))
    if any(row.get("status") == "FAIL" or row.get("finite") is False for row in rows):
        status = "FAIL"
    else:
        status = "PASS"
    payload = {
        "revision": REVISION, "status": status,
        "classification": "NONRESULT_PREFLIGHT",
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "cuda_device": torch.cuda.get_device_name(0) if device.type == "cuda" else None,
        "rows": rows,
        "test_truth_loaded": False,
    }
    write_json_atomic(RAW_DIR / "smoke_nonresults.json", payload)
    if status != "PASS":
        raise RuntimeError(payload)
    return payload


def main() -> None:
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(json.dumps(run(device), indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
