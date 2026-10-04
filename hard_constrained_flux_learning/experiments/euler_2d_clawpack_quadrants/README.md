# 2-D Euler Clawpack quadrant audit

This directory contains the prospective 64-grid HCFL experiment based on the
official PyClaw quadrant gallery problem.  Read `PROTOCOL.md` before running it.

The Clawpack dependency is isolated in a Linux container:

```text
docker build -f Dockerfile.clawpack -t hcfl-clawpack:5.10.0 ../../..
```

Generated datasets and checkpoints are intentionally not versioned.  Scripts,
metadata, scalar results, plots, and audit reports are versioned.
