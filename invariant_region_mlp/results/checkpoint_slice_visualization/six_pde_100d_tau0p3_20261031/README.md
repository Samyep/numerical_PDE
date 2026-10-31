# Six-PDE 100D frozen slice export at `t/T = 0.3`

This directory publishes the complete portable plotting data for six 100D
PDEs evaluated with Tangent-Actor, DPI, SDGD-PINN, and PINN-HTE.

Download `six_pde_100d_tau0p3_20261031.zip`, then run:

```bash
unzip six_pde_100d_tau0p3_20261031.zip
cd six_pde_100d_tau0p3_20261031
python -m pip install numpy matplotlib
python scripts/render_long_pde_figures.py
```

The archive contains 6 reference NPZ files, 24 method-prediction NPZ files,
all metrics and provenance audits, PDF/PNG/SVG figures, and per-file checksums.

Archive SHA-256:

```text
4e367b24929214829b38e64fafcee313401af30b46f82057bbdce00731ce5ed0
```
