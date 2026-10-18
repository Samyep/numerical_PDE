# Ambient 400D Generalization Audit

Known active directions are used only in this post-freeze decomposition. Every trained Actor receives the full ambient state and every learned low-rank matrix is randomly initialized and optimized end to end.

| source | d | architecture | rank | data | selected | a0 RMSE | selected RMSE | improvement | Gradient RMSE | active / inactive | effective rank | train / held-out fit | label replication | step | gap phases | FP score |
|---|---:|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| parent | 100 | full_ambient | — | 8192 | a4 | 0.0060305 | 0.000684943 | 88.6% | 0.007139 | 0.0071331 / 0.00028946 | — | 0.00221 / 0.007259 | 0.003359 | 2750 | 4/4 | 0.0035376 |
| parent | 100 | full_ambient | — | 8192 | a4 | 0.00590316 | 0.000617239 | 89.5% | 0.00684501 | 0.0068405 / 0.00024943 | — | 0.002142 / 0.006632 | 0.003171 | 2750 | 4/4 | 0.0030008 |
| parent | 100 | full_ambient | — | 8192 | a4 | 0.00603108 | 0.000666913 | 88.9% | 0.00680163 | 0.0067961 / 0.0002736 | — | 0.002284 / 0.006764 | 0.003321 | 2500 | 4/4 | 0.0031085 |
| parent | 100 | full_ambient | — | 2048 | a2 | 0.00577928 | 0.00150872 | 73.9% | 0.0607912 | 0.060729 / 0.0027515 | — | 0.009382 / 0.05346 | 0.01382 | 2000 | 2/2 | 0.14507 |
| parent | 400 | full_ambient | — | 2048 | a2 | 0.00641646 | 0.00629103 | 2.0% | 0.153572 | 0.15318 / 0.010965 | — | 0.1678 / 0.1747 | 0.01702 | 1 | 0/2 | 2.4101 |
| frontier | 400 | full_ambient | — | 8192 | a4 | 0.00639172 | 0.00122297 | 80.9% | 0.0353662 | 0.035331 / 0.0015857 | 2.00 | 0.009599 / 0.03596 | 0.005661 | 750 | 4/4 | 0.18055 |
| frontier | 400 | full_ambient | — | 8192 | a4 | 0.00663148 | 0.000920406 | 86.1% | 0.0346877 | 0.034673 / 0.0010174 | 1.99 | 0.004265 / 0.03712 | 0.005924 | 1750 | 4/4 | 0.12422 |
| frontier | 400 | full_ambient | — | 8192 | a4 | 0.00662054 | 0.00107583 | 83.8% | 0.0340846 | 0.034071 / 0.00097226 | 1.99 | 0.004194 / 0.0365 | 0.005348 | 1750 | 4/4 | 0.090316 |

The diagnosis distinguishes supervision noise from statistical generalization by combining label-replication error, train/held-out fit, checkpoint step, active/inactive post-freeze errors, and the learned output spectrum. The frozen persistent-gap trigger requires at least three update phases with held-out Gradient RMSE at least 1.5× the train RMSE and an excess of at least half the independent label-replication RMSE. A low-rank row is an explicitly named ablation and does not replace the full-ambient primary method.

## Measured cost

Across three matched-budget runs, production training plus reference-free selection costs `2871.4` GPU-seconds per seed on average. The post-freeze Actor/policy audit costs `540.8` GPU-seconds per seed. The training total sums disjoint synchronized label-generation, Actor-fit, and fixed-point-diagnostic timers across resumed Slurm invocations; queue time is excluded.

The seed-0 full-ambient run exceeded the predeclared 10% Value-improvement gate and the result was replicated with seeds 1 and 2. The conditional learned low-rank branch was therefore not triggered; its absence is not a negative low-rank result.
