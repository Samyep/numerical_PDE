# Float64 GPU MLP port

This directory contains the isolated PyTorch implementation used in Study G.
It supports `raw`, `path`, `double`, `double_path`, `path_base`, and
`double_path_base` for P1 and MR.  A shared draw-tape interface provides
sample-for-sample replay of the NumPy implementation for G-2.

```bash
python -m invariant_region_mlp.experiments.gpu_mlp.validate
python -m invariant_region_mlp.experiments.gpu_mlp.run g3 --workers 8
python -m invariant_region_mlp.experiments.gpu_mlp.run g4
python -m invariant_region_mlp.experiments.gpu_mlp.run g5
python -m invariant_region_mlp.experiments.gpu_mlp.analyze
```
