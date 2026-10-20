# High-dimensional closeout and OU/GMM FH-MLP frontier

## Outcome

The S05/B6 study is complete and pushed. The OU/GMM FH-MLP timing sweep is
**BLOCKED_MISSING_ARTIFACTS** and was not started. No timing point, RMSE point,
crossing, PNG, or PDF was fabricated.

## Stage 0: S05/B6 closeout

| Criterion | Verdict | Result |
|---|---|---|
| H1 | FAIL | 0/4; learned TA did not reach 0.5 times the terminal-tangent RMSE |
| H2 | FAIL | 0/4; one policy-improvement step suffices here under the registered threshold |
| H3 | PASS | 4/4; TA/raw-FH RMSE ratios were 0.00662, 0.0154, 0.000887, and 0.00468 |
| H4 | NOT EVALUATED (COMPUTE) | 2/2 evaluated wins, but 3 were required and two d=400 cells were curtailed |
| H5 | PASS | 2/2 required wins, both at d=100 |
| H6 | FAIL | H1 failed on B6 and selected B6 d=400 had r_glob=77.14<100 |

The registered artifacts contain 288 Actor evaluation replicates and 260/264
Full-History repetitions. Four FH repetitions were curtailed under the frozen
30-minute rule. The recorded new GPU work is approximately 55,682.44 seconds
(15.47 GPU-hours): 12,637.05 seconds of Actor training/selection, 20,978.77
seconds of Actor inference audits, and 22,066.62 seconds of newly run FH cells.
An additional 3,859.95 seconds belongs to legally reused historical S05 FH rows
and is not counted as new work.

- Branch: `ir-mlp-highdim-nosym`
- Final commit: `281977883e0a5ce4f05895878680b582620461be`
- Remote push: successful to `origin/ir-mlp-highdim-nosym`
- Full report: `invariant_region_mlp/docs/HIGHDIM_NOSYM_REPORT.md`

## Stage 1: OU/GMM pre-run audit

The RTX 5070 Ti was idle and available. The mathematical implementation passed:

- 9/9 dedicated OU/FH tests;
- maximum relative analytic PDE residual on 200 states: `4.631436141997569e-15`;
- terminal Value and Gradient identities: maximum error `0`;
- stabilized OU BEL/direct formula relative error: `2.026980648162402e-16`;
- correct OU transition variance and positive drift kernel;
- no nonterminal exact solution, oracle Actor, or reference gradient in the FH recurrence;
- the known spatially constant `-d` source is integrated deterministically.

Across the broader relevant suite, 48 tests passed. Two final-suite DPI tests
could not import the optional external `picard` package; these failures do not
exercise the OU FH implementation.

The formal artifact gate failed because these immutable files are absent locally:

1. `invariant_region_mlp/results/tangent_actor_icml_final_suite/inputs/ou_gmm_d100.npz`
2. `invariant_region_mlp/results/tangent_actor_icml_final_suite/paired_predictions/main_ou_gmm_d100_tangent_actor_s0_n64p0_p1024_mcv.npz`

The search covered all local worktrees and user storage, local ZIP archives,
all fetched Git objects/remote branches, and accessible Anvil home/scratch
storage. The documented DeltaAI source could not be read because the current
session has no noninteractive NCSA authentication.

As an audit only, the tracked official tensors and frozen seed reproduced the
960-point bank array hash exactly (`632c00f9...`) and reproduced the parameter
array hash exactly (`2437348e...`). Those reconstructed arrays were not written
or substituted for the missing originals. In particular, the frozen
Tangent-Actor RMSE could not be recomputed from its saved per-point prediction,
so the explicit fail-closed rule applies.

## OU/GMM frontier result

- Status: `BLOCKED_MISSING_ARTIFACTS`
- Completed FH configurations: 0
- RTX 5070 Ti numerical GPU time: 0 seconds
- Target reached: not evaluated
- Best FH RMSE: not available
- Curve PNG/PDF: intentionally not generated
- Target scalar retained only as registered provenance: `0.2197627563608737`
- Historical TA time retained only as GH200 provenance: `1675.9737 s`

Therefore the paper cannot presently claim that FH reaches, fails to reach, or
is faster/slower than Tangent-Actor on this frontier. The honest high-dimensional
result remains narrower: TA strongly beats raw FH on all four S05/B6 cells, but
the preregistered claims that learning is decisive and that later Actor updates
matter were falsified.

## Artifacts

- `OU_GMM_FH_TIME_ACCURACY_AUDIT.md`
- `ou_gmm_fh_time_accuracy_audit.json`
- `ou_gmm_fh_time_accuracy_summary.json`

