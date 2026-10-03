# Nonperiodic 1D Euler audit

This experiment applies the periodic-trained checkpoints zero-shot to a
nonperiodic finite-volume problem with transmissive (constant-extrapolation)
boundaries.

- Neural fluxes are evaluated only at the `N-1` interior interfaces.
- The two boundary fluxes are the physical Euler fluxes of the endpoint
  states; the network never predicts a boundary flux.
- Replicated ghost cells provide the five-cell stencil near each boundary,
  without circular wrapping.
- Tadmor projection is applied only to interior learned proposals.
- The fully discrete entropy check includes the physical boundary entropy
  flux, rather than incorrectly requiring total entropy to be nonincreasing.
- The checkpoints are not retrained; this is a zero-shot boundary-condition
  transfer test.

The test suite contains the five centered canonical shock-tube problems plus
three cases whose contact or pressure wave reaches a boundary during the
rollout.  The common metric reference is strict nonperiodic HLLC-2048,
conservatively restricted to 512 cells only for scoring.

Run:

```text
python run_nonperiodic_audit.py --seed 0
python audit_results.py --seed 0
```

Outputs are written to `results/nonperiodic_transmissive512_seed0.{png,json,csv}`.
See `RESULTS.md` for the measured outcome and limitations.
