# 100D paired-bank remediation — 2026-10-19

## Disposition

The repaired seed-0 100D whole-function table passes a strict, artifact-level
same-bank audit for all seven PDEs and all four amortized methods.  No model was
retrained, no checkpoint was changed, and no test error was used for model
selection.  The original generated tables and report were copied to files with
the suffix `PRE_PAIRING_REPAIR_20261019` before canonical regeneration.

## Nonconvex coupled finding

The apparent 960-versus-512 mismatch was a metadata defect, rather than a
different test population.  Two incompatible quantities had both been called
`test_bank_sha256`:

- `85a9430c6f60762f5964a9de46122e8010cd1d96b67b286a6997a824d20c2180`
  is the content hash of the numerical `(t,x)` arrays.
- `3aef1dc9855657a66f1744f0cdfe3d12e80b009f52b158154a1da799e42044de`
  is the byte-level SHA-256 of the NPZ container holding those arrays.

The historical Tangent-Actor aggregator also hard-coded 960 states.  Opening
the frozen prediction artifacts shows that Tangent-Actor, DPI, PINN-HTE and
SDGD-PINN all use exactly the same 512 `t` values, the same `512 x 100` state
array and the same reference array.  Array equality is exact, with zero maximum
absolute difference.  Canonical metadata now reports 512 states, stores file
and array hashes in separate fields, and recomputes every metric from the
frozen prediction and canonical reference.  All four nonconvex RMSEs are
bit-for-bit unchanged.

## Allen–Cahn reference repair

All four Allen–Cahn methods already used the same 410 test states.  The strict
canonical state-array hash is
`8c1a0e9f69953a06660f1ef10db8ceac261f2d06318d8ad5d1fa5892f08ce1f2`.
The baseline prediction NPZs embedded a slightly coarser radial reference,
while the original high-accuracy evaluation bank contains the Richardson
reference with hash
`d377fba6f91460e87e072c6960e246244230af69d5946a718434123d1687e986`.
The two truth arrays differ by at most `9.159e-7` and have RMS difference
`6.301e-8`.  Every method is now evaluated against the latter reference, whose
recorded uncertainty is `9.645e-6`.

The resulting RMSE changes are:

| Method | Previous RMSE | Canonical RMSE | Change |
|---|---:|---:|---:|
| Tangent-Actor | 0.0001733811866 | 0.0001733811866 | 0 |
| DPI | 0.001066548999 | 0.001066535425 | -1.36e-8 |
| PINN-HTE | 0.0009479992857 | 0.0009479920998 | -7.19e-9 |
| SDGD-PINN | 0.002692391823 | 0.002692403597 | +1.18e-8 |

These changes are far below the reference uncertainty and do not alter the
ranking.

## New guards

`verify_paired_banks.py` now opens every canonical prediction artifact rather
than trusting CSV metadata.  For each PDE and method it independently checks:

1. numerical `(t,x)` array hash and state count;
2. numerical reference-array hash;
3. canonical prediction-file SHA-256;
4. frozen source-prediction SHA-256;
5. finite arrays and exact shape compatibility;
6. RMSE, relative L2, bias and MAE recomputed from the arrays; and
7. equality of the actual bank and reference hashes across all four methods;
8. exact equality between each canonical prediction and the prediction
   extracted from its immutable source artifact; and
9. exact equality between canonical `t`, `x`, and truth arrays and the
   authoritative frozen bank/reference arrays.

The guard writes `PAIRWISE_BANK_VERIFICATION_100D.json` and exits nonzero on
any mismatch.  Its current status is `PASS`, and Tangent-Actor has the smallest
frozen RMSE in 7/7 strictly paired 100D banks.  A unit test modifies an actual
state array while leaving every CSV string unchanged and confirms that the
audit fails.  Future Tangent-Actor, DPI, PINN-HTE and SDGD runners now record
array hashes and NPZ file hashes separately.

After the independent review, the recurring verifier was hardened with items
8--9 above.  The aggregator now also copies shared row dictionaries before
repair, so every generated `main_*.npz` and `frontier_*.npz` is referenced by
its own table.  The Allen--Cahn free-form reference label is regenerated from
the canonical high-accuracy source.  The resulting package test status is
`24 passed`.

The broader scientific caveats from the independent review remain in force:
these are single-seed results, most external baselines are method transfers,
several PDEs have low effective dimension, and the comparison does not imply
uniform native-task or compute superiority.
