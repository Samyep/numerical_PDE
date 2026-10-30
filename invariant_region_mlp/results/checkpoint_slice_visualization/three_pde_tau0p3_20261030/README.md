# Three-PDE frozen slice export at `t/T = 0.3`

This directory publishes the complete portable data bundle used to draw the
long frozen-checkpoint comparison figures for:

- 400D Jointly Concave;
- 100D Nonconvex--Nonconcave;
- 100D OU/GMM HJB.

Each PDE compares Tangent-Actor, DPI, SDGD-PINN, and PINN-HTE on the same
81-by-81 physical slice. No training or checkpoint selection is performed by
the plotting code.

## Download and reproduce

Download `three_pde_tau0p3_20261030.zip`, then run:

```bash
unzip three_pde_tau0p3_20261030.zip
cd three_pde_tau0p3_20261030
python -m pip install numpy matplotlib
python scripts/render_long_pde_figures.py
```

The archive contains three reference NPZ files, twelve method-prediction NPZ
files, all value/gradient/error arrays, metrics and provenance audits, plotting
scripts, vector PDFs, PNG previews, and per-file checksums. See
`BUNDLE_README.md` and `DATA_MANIFEST.json` for the exact schema.

Archive SHA-256:

```text
a816300222749fb185ae58754850774af7d96687a5f0a29aed43535ce5465e65
```

The archive passed `zip -T`; all internal payloads passed `sha256sum -c`.
