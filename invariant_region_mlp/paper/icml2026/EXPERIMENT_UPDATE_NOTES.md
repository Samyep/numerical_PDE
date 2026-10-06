# 本轮实验并入说明

## 版本与范围

以 `Constraint_Aware_MLP_ICML2026_related_work.zip` 为底稿。直接核对 GitHub
`Samyep/numerical_PDE` 的 `771f26594aa7ee364d426ad29744b7716d76cdc3` 版本，读取用户指定的
四个实验脚本、两个 JSON 结果文件及两份报告。没有重新运行 PDE 实验。
当前是正文到第 7 页、全文 18 页的匿名 ICML 样式稿。

## 实质修改

1. Abstract 和第三项 contribution 更新为 headline HJB + deeper funding + matched diagnostics。
   不增加新的 convergence theorem；Related Work 保留。
2. Section 6.2 将 100–160D HJB、n=3/4 funding 提升为主实验，新增 Figure 3 两个面板。
3. Section 6.3 新增 Figure 4 两个面板，分别展示同批数据上的 generator MSE 和
   **实际平均 level correction 的跨重复方差**。
4. Section 6.4 加入 M=2 Neufeld–Wu 平局和 batchwise structural sanity，并保留 negative controls。
5. Appendix E 提供完整 protocol、7 张新增结果表、指标公式及假设边界。
6. 旧 15 个图形面板全部保留；其中 6 个旧正文面板移至匹配的附录小节，未覆盖任何旧图文件。
   当前 19 个面板；32 个原图/数据/归档文件 SHA-256 检查均未变。

## 严格区分的事项

- 新 HJB Raw = **无修正**，不是旧版 heuristic clipping；不把新 Batch 曲线混入旧 paired 曲线。
- 新 HJB 在每个维度 10/10 paired repetitions 更低，不能说有 40 个独立 PDE 实例。
- HJB 新图误差线是每个方法跨 10 次重复的 **sample SD**，不是 SE 或 CI。
- HJB generator 诊断有 6,400 个 child states（64 parents × 10 children × 10 repetitions）；
  correction variance 则在 1,200 个固定 parents 上各自跨 10 次重复计算，再平均。
- Individual-increment pooled variance 在 d=120/140/160 **更高**，仍全部写入 Table 8；
  “所有 variance 都下降”的说法不成立。
- Funding 改变 n 的同时也改变 M，不能把三个点当作 fixed-work depth trend。
- Funding p-values 来自仓库 aggregate summary，没有重新计算或捏造未提供的 replica-level CI。
- Neufeld–Wu 两方法是平局，不计作 Batch win；代码使用 12-step Euler，误差包含离散化影响。
- Counterexample batch sanity 的 **value outputs** 逐次相同；batch root gradient 没有投影，
  因此不是完整 (u,z) MSE 相等的验证，也没有替代原来的 65,536-replica 独立 MSE 检查。
- Finance envelope 的现有解析证明限制保持不变；新实验没有自动补齐该证明。
- Batch-IR sharp consistency / cost-to-accuracy 理论仍为 future work，不声称 universal dominance。

## 数值复核

对 HJB JSON 内 120 个 error 值（4 dimensions × 3 methods × 10 repetitions）复算均值、
样本标准差、paired gain、win count 和 paired t-test，均与记录一致。
Funding reductions、Neufeld 平局、counterexample 最大 paired value discrepancy 也作了数值一致性检查。
可运行 `python scripts/audit_extended_results.py`。

## 仍待补的旧图

此前 31 页稿的 overshoot-energy / final-gain 原图及完整 36 点数据仍不在输入包中。
本次保留其文字统计，不猜测或 digitize 新数据点。所有**输入版本已经包含**的图均保留。

## GitHub 状态

本轮仅交付更新后的稿件 ZIP / PDF，未写入 GitHub。上述 SHA 是实验输入版本。
