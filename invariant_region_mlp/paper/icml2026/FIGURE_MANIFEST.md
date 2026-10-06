# 全图版逐图目录 / Complete evidence manifest

最新成稿：正文 7 页，全文 27 页；保留输入版全部 19 个面板，新增 15 个独立图。共 34 个独立面板、35 处展示（Neufeld–Wu 两预算图同时用于主文和附录）。

| 图号 | PDF页 | 面板数 | 状态 | 内容 |
|---|---:|---:|---|---|
| 1 | 2 | 1 | Original | Recursive interface / method overview |
| 2 | 5 | 2 | Original | Nine-dimension full-state rescue and value-only ablation |
| 3 | 5 | 2 | Original extended suite | Headline HJB and deeper funding with Batch-IR |
| 4 | 6 | 2 | Original extended suite | Matched generator MSE and averaged-correction variance |
| 5 | 7 | 2 | NEW reproduction + old panel repeated | Reproduced old Figure-7 diagnostic; independent two-budget Neufeld–Wu returned to main |
| A1 | 10 | 2 | Original | Independent 65,536-replica MSE and separate 512-replica primary-sweep RMSE |
| A2 | 11 | 1 | Original | Complete four-depth counterexample sweep |
| A3 | 12 | 2 | Original | Corrected HJB geometry comparison and Neufeld–Wu budgets |
| A4 | 13 | 2 | Original | Earlier two-method generator MSE and correction variance |
| A5 | 13 | 1 | Original PNG | Original cumulative funding work envelope |
| A6 | 14 | 2 | Original | Shallow HJB and funding mean-preserving tradeoff |
| A7 | 15 | 1 | Original | First-level HJB activity; not the 36-setting overshoot sweep |
| A8 | 15 | 1 | Original | Earlier Allen–Cahn negative control |
| A9 | 20 | 1 | NEW reproduction | Companion violation-frequency / gain curve, all 36 settings |
| A10 | 22 | 2 | Recovered archive data | All seven budget points plus n=4/n=5 joint-versus-box ablation |
| A11 | 23 | 2 | Git report values | Mean-preserving feasibility and n=M=3 Neufeld–Wu ties |
| A12 | 24 | 2 | Recovered 40-run archive | Controlled analytic surrogate at M=1 and M=2 |
| A13 | 24 | 1 | Recovered 40-run archive | M=4 surrogate boundary case, including IR losses |
| A14 | 26 | 2 | Git JSON metric subset | Three-method Allen–Cahn and credit-risk overlaps |
| A15 | 26 | 2 | Git JSON metric subset | Linear value overlap despite active gradient correction |
| A16 | 27 | 2 | Archived rounded report | Outside-MLP neural-BSDE diagnostics, with unfavorable results retained |

## 原图 7：恢复的是可复现的实验，而不是猜测图片点位

原来的 `overshoot_energy_gain.png` 仍只能预览，未取得原始 PNG 字节。正文 Figure 5(a) 是同类 36 设置诊断的独立复现，全部点来自本次实际计算，不是从原图读点，也不是原始随机样本。

- d=100，n=2；M=2,3,4,6,8,10；半径倍率=1,1.25,1.5,2,3,4；每个预算十次配对重复。
- 仅约束 child gradient；不裁剪最终 u；使用修正的 scaled Hopf–Cole 参考积分。
- 本次 Spearman(overshoot,gain)=0.886229；Spearman(violation,gain)=0.994595。旧 corrected-reference 报告对应 0.889 和 0.996，文中分别标明。
- 原始评估点、矩阵、预测数组、分重复统计、36行JSON及脚本都在包内。

## 不能混用的结果

旧 unscaled-reference HJB error 图和旧 sweep 数据保留在 `archive_superseded/`，但不恢复为当前有效证据。当前 HJB geometry 图使用修正后的参考值。

Neufeld–Wu 的“相对 Raw 有优势”和“Samplewise/Batchwise 打平”是两个不同比较，均保留。Funding 原有 work-envelope 图、n=4/n=5 的 joint-versus-box 消融，与新 Batch-IR 三方法实验分别展示，不拼接成一条假配对曲线。

Controlled surrogate 只使用 q*u*，不是 trained SCaSML；elliptic 是 neural BSDE，不是 MLP。本版保留其中 IR 不占优的结果。

## 机器核对

44 个输入版图像、向量图、数据和归档图资产逐字节校验：缺失 0，修改 0。15 个新增图全部被 LaTeX 引用。36 设置统计由原始数组重新计算，最大差异 0。详见 `data/restored_evidence/restoration_audit.json`。

本次主仓库实验读取版本：`771f26594aa7ee364d426ad29744b7716d76cdc3`。恢复稿尚未写回 GitHub。
