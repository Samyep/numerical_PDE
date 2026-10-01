"""Offline PyClaw/SharpClaw reference-data generator for HCFL.

This script is intentionally independent from the differentiable PyTorch HCFL
rollout. It generates high-resolution reference trajectories and conservatively
restricts them to a coarse grid.

Supported benchmarks:
  - sod
  - shu_osher
  - woodward_colella
  - quadrants (2D Euler)

Recommended installation:
  pip install clawpack==v5.10.0
"""

import argparse
import json
from pathlib import Path
import numpy as np

GAMMA = 1.4


def require_clawpack():
    try:
        from clawpack import pyclaw, riemann
    except Exception as exc:
        raise RuntimeError(
            "PyClaw is not installed. Recommended: pip install clawpack==v5.10.0"
        ) from exc
    return pyclaw, riemann


def prim_to_cons_1d(rho, u, p):
    E = p / (GAMMA - 1.0) + 0.5 * rho * u * u
    return np.stack([rho, rho * u, E], axis=-1)


def prim_to_cons_2d(rho, u, v, p):
    E = p / (GAMMA - 1.0) + 0.5 * rho * (u*u + v*v)
    return np.stack([rho, rho*u, rho*v, E], axis=-1)


def pressure_1d(q):
    rho, m, E = q[...,0], q[...,1], q[...,2]
    return (GAMMA - 1.0) * (E - 0.5*m*m/rho)


def pressure_2d(q):
    rho, mx, my, E = q[...,0], q[...,1], q[...,2], q[...,3]
    return (GAMMA - 1.0) * (E - 0.5*(mx*mx+my*my)/rho)


def coarsen_1d(q, coarse_nx):
    fine_nx = q.shape[-2]
    if fine_nx % coarse_nx:
        raise ValueError("fine_nx must be divisible by coarse_nx")
    r = fine_nx // coarse_nx
    return q.reshape(*q.shape[:-2], coarse_nx, r, q.shape[-1]).mean(axis=-2)


def coarsen_2d(q, coarse_ny, coarse_nx):
    fine_ny, fine_nx, neq = q.shape[-3:]
    if fine_ny % coarse_ny or fine_nx % coarse_nx:
        raise ValueError("fine shape must be divisible by coarse shape")
    ry, rx = fine_ny // coarse_ny, fine_nx // coarse_nx
    return q.reshape(
        *q.shape[:-3], coarse_ny, ry, coarse_nx, rx, neq
    ).mean(axis=(-4, -2))


def ic_sod(x):
    left = x < 0.5
    rho = np.where(left, 1.0, 0.125)
    u = np.zeros_like(x)
    p = np.where(left, 1.0, 0.1)
    return prim_to_cons_1d(rho,u,p)


def ic_shu_osher(x):
    left = x < 1.0
    rho = np.where(left, 3.857143, 1.0 + 0.2*np.sin(5.0*x))
    u = np.where(left, 2.629369, 0.0)
    p = np.where(left, 10.33333, 1.0)
    return prim_to_cons_1d(rho,u,p)


def ic_woodward_colella(x):
    rho = np.ones_like(x)
    u = np.zeros_like(x)
    p = np.where(x < 0.1, 1000.0, np.where(x < 0.9, 0.01, 100.0))
    return prim_to_cons_1d(rho,u,p)


def ic_quadrants(xc, yc):
    X, Y = np.meshgrid(xc, yc, indexing="xy")
    l, b = X < 0.8, Y < 0.8
    r, t = ~l, ~b

    rho = (
        1.5*r*t
        + 0.532258064516129*l*t
        + 0.137992831541219*l*b
        + 0.532258064516129*r*b
    )
    u = (
        0.0*r*t
        + 1.206045378311055*l*t
        + 1.206045378311055*l*b
        + 0.0*r*b
    )
    v = (
        0.0*r*t
        + 0.0*l*t
        + 1.206045378311055*l*b
        + 1.206045378311055*r*b
    )
    p = (
        1.5*r*t
        + 0.3*l*t
        + 0.029032258064516*l*b
        + 0.3*r*b
    )
    return prim_to_cons_2d(rho,u,v,p)


def run_1d(problem, fine_nx, coarse_nx, tfinal, nout, solver_type):
    pyclaw, riemann = require_clawpack()

    if solver_type == "sharpclaw":
        solver = pyclaw.SharpClawSolver1D(riemann.euler_with_efix_1D)
        solver.weno_order = 5
        solver.time_integrator = "SSP33"
        solver.cfl_desired = 0.6
        solver.cfl_max = 0.7
    else:
        solver = pyclaw.ClawSolver1D(riemann.euler_with_efix_1D)
        solver.cfl_desired = 0.8
        solver.cfl_max = 0.9

    if problem == "sod":
        xlower, xupper, ic = 0.0, 1.0, ic_sod
        defaults = 0.20
        solver.bc_lower[0] = pyclaw.BC.extrap
        solver.bc_upper[0] = pyclaw.BC.extrap
    elif problem == "shu_osher":
        xlower, xupper, ic = 0.0, 10.0, ic_shu_osher
        defaults = 1.80
        solver.bc_lower[0] = pyclaw.BC.extrap
        solver.bc_upper[0] = pyclaw.BC.extrap
    elif problem == "woodward_colella":
        xlower, xupper, ic = 0.0, 1.0, ic_woodward_colella
        defaults = 0.038
        solver.bc_lower[0] = pyclaw.BC.wall
        solver.bc_upper[0] = pyclaw.BC.wall
    else:
        raise ValueError(problem)

    if tfinal is None:
        tfinal = defaults

    xdim = pyclaw.Dimension(xlower, xupper, fine_nx, name="x")
    domain = pyclaw.Domain(xdim)
    state = pyclaw.State(domain, 3)
    state.problem_data["gamma"] = GAMMA
    state.problem_data["gamma1"] = GAMMA - 1.0
    state.q[...] = np.moveaxis(ic(domain.grid.x.centers), -1, 0)

    claw = pyclaw.Controller()
    claw.solution = pyclaw.Solution(state, domain)
    claw.solver = solver
    claw.tfinal = tfinal
    claw.num_output_times = nout
    claw.keep_copy = True
    claw.output_format = None
    claw.verbosity = 0
    claw.run()

    qfine = np.stack(
        [np.moveaxis(np.asarray(frame.q), 0, -1) for frame in claw.frames],
        axis=0,
    )
    times = np.asarray([float(frame.t) for frame in claw.frames])
    q = coarsen_1d(qfine, coarse_nx).astype(np.float32)

    if (q[...,0] <= 0).any() or (pressure_1d(q) <= 0).any():
        raise RuntimeError("Generated 1D dataset contains non-admissible states")

    return q, times, {
        "problem": problem,
        "dimension": 1,
        "system": "euler",
        "gamma": GAMMA,
        "solver_type": solver_type,
        "fine_nx": fine_nx,
        "coarse_nx": coarse_nx,
        "tfinal": tfinal,
    }


def run_quadrants(fine_nx, fine_ny, coarse_nx, coarse_ny, tfinal, nout, solver_type, riemann_name):
    pyclaw, riemann = require_clawpack()

    if solver_type == "sharpclaw":
        solver = pyclaw.SharpClawSolver2D(riemann.euler_5wave_2D)
        solver.weno_order = 5
        solver.time_integrator = "SSP33"
    else:
        if riemann_name == "hlle":
            solver = pyclaw.ClawSolver2D(riemann.euler_hlle_2D)
            solver.transverse_waves = 0
        elif riemann_name == "5wave":
            solver = pyclaw.ClawSolver2D(riemann.euler_5wave_2D)
        else:
            solver = pyclaw.ClawSolver2D(riemann.euler_4wave_2D)
            solver.transverse_waves = 2

    solver.all_bcs = pyclaw.BC.extrap

    xdim = pyclaw.Dimension(0.0, 1.0, fine_nx, name="x")
    ydim = pyclaw.Dimension(0.0, 1.0, fine_ny, name="y")
    domain = pyclaw.Domain([xdim, ydim])
    state = pyclaw.State(domain, 4)
    state.problem_data["gamma"] = GAMMA
    state.problem_data["gamma1"] = GAMMA - 1.0

    q0 = ic_quadrants(domain.grid.x.centers, domain.grid.y.centers)
    state.q[...] = np.moveaxis(np.swapaxes(q0,0,1), -1, 0)

    claw = pyclaw.Controller()
    claw.solution = pyclaw.Solution(state, domain)
    claw.solver = solver
    claw.tfinal = tfinal
    claw.num_output_times = nout
    claw.keep_copy = True
    claw.output_format = None
    claw.verbosity = 0
    claw.run()

    frames = []
    for frame in claw.frames:
        arr = np.moveaxis(np.asarray(frame.q), 0, -1)  # (mx,my,neq)
        arr = np.swapaxes(arr, 0, 1)                   # (my,mx,neq)
        frames.append(arr)
    qfine = np.stack(frames, axis=0)
    times = np.asarray([float(frame.t) for frame in claw.frames])
    q = coarsen_2d(qfine, coarse_ny, coarse_nx).astype(np.float32)

    if (q[...,0] <= 0).any() or (pressure_2d(q) <= 0).any():
        raise RuntimeError("Generated 2D dataset contains non-admissible states")

    return q, times, {
        "problem": "quadrants",
        "dimension": 2,
        "system": "euler",
        "gamma": GAMMA,
        "solver_type": solver_type,
        "riemann_solver": riemann_name,
        "fine_nx": fine_nx,
        "fine_ny": fine_ny,
        "coarse_nx": coarse_nx,
        "coarse_ny": coarse_ny,
        "tfinal": tfinal,
    }


def save_npz(output, q, times, metadata):
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        q=q,
        times=times,
        metadata_json=np.asarray(json.dumps(metadata, sort_keys=True)),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--problem", choices=["sod","shu_osher","woodward_colella","quadrants"], required=True)
    p.add_argument("--solver-type", choices=["classic","sharpclaw"], default="sharpclaw")
    p.add_argument("--riemann-solver", choices=["roe","hlle","5wave"], default="roe")
    p.add_argument("--fine-nx", type=int, default=1024)
    p.add_argument("--fine-ny", type=int, default=None)
    p.add_argument("--coarse-nx", type=int, default=64)
    p.add_argument("--coarse-ny", type=int, default=64)
    p.add_argument("--tfinal", type=float, default=None)
    p.add_argument("--num-output-times", type=int, default=16)
    p.add_argument("--output", required=True)
    a = p.parse_args()

    if a.problem == "quadrants":
        q, times, meta = run_quadrants(
            fine_nx=a.fine_nx,
            fine_ny=a.fine_ny or a.fine_nx,
            coarse_nx=a.coarse_nx,
            coarse_ny=a.coarse_ny,
            tfinal=0.8 if a.tfinal is None else a.tfinal,
            nout=a.num_output_times,
            solver_type=a.solver_type,
            riemann_name=a.riemann_solver,
        )
    else:
        q, times, meta = run_1d(
            problem=a.problem,
            fine_nx=a.fine_nx,
            coarse_nx=a.coarse_nx,
            tfinal=a.tfinal,
            nout=a.num_output_times,
            solver_type=a.solver_type,
        )

    save_npz(a.output, q, times, meta)
    print(json.dumps({"shape": list(q.shape), "metadata": meta}, indent=2))
