# Six-cell normal HCFL with fixed transverse transport

This experiment keeps the learned 2-D flux as close as possible to the retained
1-D method:

- one shared network for x and y faces;
- six primitive-variable cells only along the face normal;
- central physical flux plus nonnegative entropy-fixed Roe dissipation;
- proposal feasibility loss and final hard interface entropy projection;
- a fixed Roe transverse increment-wave split with zero trainable parameters.

The transverse term implements the increment-wave portion of Clawpack's classic
2-D wave-propagation algorithm (`transverse_waves=1`).  It is deliberately not
labelled `transverse_waves=2`, because level 2 additionally transports an
explicit second-order normal correction wave that is not part of the 1-D HCFL
method.  This choice preserves an exact reduction to the 1-D update whenever
the state is constant in the transverse direction.

The implementation was checked against the Clawpack 5.10 reference routines
[`flux2.f90`](https://github.com/clawpack/pyclaw/blob/v5.10.0/src/pyclaw/classic/flux2.f90)
and [`rpt2_euler.f90`](https://github.com/clawpack/riemann/blob/v5.10.0/src/rpt2_euler.f90).

Training and evaluation use the disjoint diverse64 data splits generated on a
fine grid and conservatively restricted to 64 by 64.  Pass their location with
`--data-dir` or the `HCFL_DIVERSE64_DATA` environment variable.

The complete three-seed outcome, same-checkpoint transverse ablation, safety
audit, and figure index are recorded in [RESULTS.md](RESULTS.md).
