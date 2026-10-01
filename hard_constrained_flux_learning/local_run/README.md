# Local HCFL runs

This folder is the clean entry point for running the recent HCFL ML ablations on a local machine.

## 1. Environment

Recommended:

```bash
conda create -n hcfl python=3.11 -y
conda activate hcfl
pip install numpy torch
```

The scripts use CPU by default, but PyTorch will use CUDA automatically if you modify the device placement.

## 2. 1D Euler model ablations

`euler_ablation_runner.py` is self-contained and does **not** depend on the older `/mnt/data` experiment files.

Run the current direct-vector HLLC-HCFL baseline:

```bash
python euler_ablation_runner.py --model direct --seed 0 --iters 1100
```

Other architectures:

```bash
python euler_ablation_runner.py --model invariant --seed 0 --iters 1100
python euler_ablation_runner.py --model characteristic --seed 0 --iters 1100
python euler_ablation_runner.py --model dissipation --seed 0 --iters 1100
python euler_ablation_runner.py --model conv --seed 0 --iters 1100
```

For several seeds:

```bash
for s in 0 1 2 3 4; do
  python euler_ablation_runner.py --model direct --seed $s --iters 1100
done
```

Outputs are written under `outputs/euler_ablation/`.

### Model meanings

- `direct`: HLLC + learned 3-vector flux correction.
- `invariant`: dimensionless / Galilean-aware local inputs, still a direct flux correction.
- `characteristic`: correction predicted in a normalized Roe characteristic basis.
- `dissipation`: learns characteristic-wise changes in numerical viscosity.
- `conv`: same direct-vector idea, but the 5-cell receptive field is implemented with Conv1d.

All five use the same hard Tadmor projection after the learned proposal.

## 3. Reference data with PyClaw

See `pyclaw_reference.py`.

Recommended install:

```bash
pip install clawpack==v5.10.0
```

A Fortran compiler is normally required.

The PyClaw file is intended only for **offline reference data / classical benchmarks**. The HCFL training rollout remains in PyTorch.

## 4. Important status

The older files under `hard_constrained_flux_learning/experiments/` record the exploratory research history. Some of them were originally written in the ChatGPT runtime and contain historical `/mnt/data` defaults. For new local experiments, use this `local_run/` folder first.

The current reference-data generator inside `euler_ablation_runner.py` is still the exploratory high-resolution Rusanov + SSP-RK2 teacher. The planned publication runs should replace that teacher with converged PyClaw / SharpClaw references.
