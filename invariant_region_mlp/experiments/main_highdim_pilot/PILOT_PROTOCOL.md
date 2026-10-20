# Main Experiment Expansion — Three High-Dimensional PDE Pilots on RTX 5070 Ti

你是负责 `Samyep/numerical_PDE` 中 `invariant_region_mlp` 研究的 coding/research agent。

**任务：先在本地 NVIDIA RTX 5070 Ti 上，对三个主实验候选 PDE 做小规模、严格可复现的 pilot，判断是否值得进入正式 ICML 主实验。完成验证和小实验后自动分析、画图、提交代码及结果。不要启动完整大规模实验。**

## 0. 当前研究背景与已有结果

我们的方法是 Tangent-Actor：用 learned Actor 对半线性 PDE 的 nonlinear driver 做局部线性化，通过线性 Feynman–Kac 和 martingale control variates (MCV) 求解，并迭代更新 Actor。

已完成的 `ir-mlp-highdim-nosym` 实验包括 B6、S05 的 100D/400D 结果。其中 B6 是我们这次三个主实验之一。

已知 B6 结果：

| d | a0 RMSE | a4 RMSE |
|---|---:|---:|
| 100 | 0.00164973 | 0.00133404 |
| 400 | 0.00387309 | 0.00288083 |

B6 的现有结果、参数选择、checkpoint 和测试集已经冻结。**优先复用并审计，禁止未经必要性证明而重新训练 B6。**

原来七个 PDE 不需要重新运行。

## 1. Git 与实验安全

先检查当前 Git branch、工作目录、已有 GPU 进程及最近实验状态。

- 创建独立分支 `ir-mlp-main-highdim-pilot-5070` 和新的工作目录。
- 不覆盖或删除 `highdim_nosym`、ICML final suite、OU/GMM 和其他实验的任何结果。
- 新代码放入 `invariant_region_mlp/experiments/main_highdim_pilot/`。
- 新结果放入 `invariant_region_mlp/results/main_highdim_pilot/`。
- 在运行前提交 `PILOT_PROTOCOL.md`，冻结公式、随机种子、样本分布和预算。
- 只使用当前 5070 Ti 一张 GPU；如果已有其他 GPU 训练任务，先排队而不是抢占。
- 支持 checkpoint/resume、原子写入和异常恢复。

禁止根据测试 Ground Truth 反复调超参数。

## 2. 三个 PDE

### PDE A — B6 Full-Rank Log-Sum-Exp HJB

复用 `highdim_nosym` 中已经验证的 B6。

\[
u_t+\Delta u-\frac12\|z\|^2=0,
\qquad z=\sqrt2\nabla u.
\]

\[
g(x)=-\log\left(\frac1K\sum_{m=1}^K e^{b_m^\top x}\right).
\]

精确解：

\[
u^*(t,x)=
-\log\left(
\frac1K\sum_{m=1}^K
e^{b_m^\top x+(T-t)\|b_m\|^2}
\right).
\]

任务：

1. 核对原始 B6 参数冻结文件、公式、checkpoint、bank hash、seed、训练和预测记录。
2. 读取已有 100D/400D a0、a1、a4 的 RMSE、梯度有效秩、MC 分量和训练时间。
3. 不重新校准 beta/gamma，不重新挑选更有利的测试分布。
4. 形成一个 B6 reusable benchmark adapter，供后面的正式 main suite 使用。
5. 将 B6 现有结果标注为 `REUSED_MATURE_RESULT`，不得与新 pilot 的缩减训练预算伪装成同等预算比较。

B6 作为已有 full-rank quadratic HJB anchor。

### PDE B — Dense Fourier Joint Value–Gradient Nonlinear PDE

新增：

\[
u_t+\Delta u+\frac12\tanh(u)\|z\|^2=0,
\qquad z=\sqrt2\nabla u,
\quad T=0.5.
\]

固定

\[
g(x)=\operatorname{arsinh}(h(x)),
\]

\[
h(x)=0.5+\frac2{\sqrt{2d}}
\sum_{j=1}^{2d}a_j\sin(\omega_j^\top x+\phi_j).
\]

使用冻结随机种子；\(\omega_j\) 为稠密随机单位向量，\(a_j\in\{-1,1\}\)，\(\phi_j\sim U[0,2\pi]\)。参数只生成一次，保存完整数组及 hash。

令 \(v=\sinh u\)，则 \(v_t+\Delta v=0\)，由此得到精确解：

\[
u^*(t,x)=
\operatorname{arsinh}\left[
0.5+\frac2{\sqrt{2d}}
\sum_{j=1}^{2d}
a_j e^{-(T-t)\|\omega_j\|^2}
\sin(\omega_j^\top x+\phi_j)
\right].
\]

解析计算 \(u^*\)、\(\nabla u^*\) 和 \(z^*=\sqrt2\nabla u^*\)。

**特别注意：这个 PDE 同时依赖 \(u,z\)，不能套用 gradient-only Actor。**

应使用 value-gradient 双头：

\[
r_\theta(t,x)=g(x)+(T-t)R_\theta(t,x),
\]

\[
q_\theta(t,x)=\sqrt2\nabla g(x)+(T-t)Q_\theta(t,x).
\]

在 Actor \((r,q)\) 处线性化：

\[
f(u,z)=\frac12\tanh(u)\|z\|^2,
\]

\[
c=f_u(r,q)
=\frac12\operatorname{sech}^2(r)\|q\|^2,
\]

\[
b=f_z(r,q)=\tanh(r)q,
\]

\[
h=f(r,q)-cr-b^\top q.
\]

线性 PDE：

\[
v_t+\Delta v+c\,v+b^\top(\sqrt2\nabla v)+h=0.
\]

必须正确处理 FK 中的 drift、reaction exponential weight、source integral 以及 MCV martingale。

优先复用并扩展仓库已验证的 joint value-gradient FK 代码，不要自己重新发明一套未经验证的 solver。

若现有 MCV 不支持 joint driver，增加独立适配和单元测试。不能通过忽略 \(u\)-dependence 或强制 \(c=0\) 来让代码运行。

进行验证：

- 终端条件恒等式；
- 精确解与梯度对照 autodiff；
- PDE residual；
- \(c,b,h\) 与自动微分一致；
- tangent reconstruction identity；
- 零 Actor 更新和终端极限；
- discount weight 的符号与梯度正确性；
- plain FK / MCV 一致性及 MC 方差检查。

该 PDE 可以通过 \(\sinh\) 变量变换线性化，应明确标注为 **analytically integrable high-dimensional stress test**，不能称为任意不可积非线性 PDE。

### PDE C — Published Coupled LQG HJB

参考 ICLR 2026：

*Physics-Informed Inference Time Scaling for Solving High-Dimensional Partial Differential Equations*（SCaSML）。

优先核对论文 Section 3.3 和官方代码：

https://proceedings.iclr.cc/paper_files/paper/2026/file/cf8353a0efc3644454a5591a7cb55c17-Paper-Conference.pdf

https://github.com/Francis-Fan-create/SCaSML

其 HJB 是：

\[
u_t+\Delta u-\|\nabla u\|^2=0.
\]

终端：

\[
g(x)=\log\left[
\frac{
1+\sum_{i=1}^{d-1}
\big(c_{1,i}(x_i-x_{i+1})^2+c_{2,i}x_{i+1}^2\big)
}{2}
\right],
\]

其中 \(c_{1,i},c_{2,i}\in[0.5,1.5]\)。

**以论文和官方实现为最终依据。** 优先复用官方 coefficients、时间区间和测试分布；如果只能按照论文重新生成随机参数，清楚标注 `RECREATED_PUBLISHED_BENCHMARK`，不能称为 official exact reproduction。

记

\[
g(x)=\log\frac{1+x^\top Ax}{2},
\]

其中 \(A\) 是由相邻坐标耦合构造的对称正定三对角矩阵。

通过 Hopf–Cole：

\[
u^*(t,x)=
-\log\mathbb E\left[
\frac2{1+X_T^\top AX_T}
\right],
\]

\[
X_T=x+\sqrt{2(T-t)}Z.
\]

进一步实现独立的一维积分 reference：

\[
e^{-u^*(t,x)}
=
2\int_0^\infty
\frac{
e^{-s-sx^\top A(I+4(T-t)sA)^{-1}x}
}{
\sqrt{\det(I+4(T-t)sA)}
}\,ds.
\]

用 \(A\) 的 eigen-decomposition 高效计算 determinant 和 quadratic form，不能对每个求积节点显式求大矩阵逆。

同时推导、实现相应解析积分形式的 \(z^*=\sqrt2\nabla u^*\)。

严格验证：

- 矩阵 \(A\) 的构造；
- Hopf–Cole 公式；
- 求积阶数加密收敛；
- 与独立高样本 MC reference 一致；
- PDE residual；
- terminal condition 与 gradient。

若数值 reference 不稳定，必须停在 `REFERENCE_NOT_VALIDATED`，不能将不可靠的 MC 估计冒充 Ground Truth。

## 3. 第一阶段：只做数学和问题难度预检

先不训练网络。

维度只取 `d=100`。

- B6：读取已有冻结 calibration。
- Dense Fourier：使用冻结参数生成新的测试集。
- LQG：使用论文兼容参数和数据分布。

Dense Fourier 使用 \(t\sim U[0,T)\)、\(x\sim U[-1,1]^{100}\) 作为 pilot 分布。

LQG 优先使用论文中的时间区间和 unit-ball 状态分布，不要为了得到更好的数字偷偷换成 Gaussian。

每个新 PDE 固定 512 个测试点和独立验证点。

报告：

1. \(u^*\) 的标准差及取值范围；
2. \(z^*\) RMS；
3. centred gradient covariance effective rank

\[
r_{\rm eff}=\frac{(\operatorname{tr}C_z)^2}
{\operatorname{tr}(C_z^2)};
\]

4. 前十个方向之外的梯度方差比例；
5. 非线性 driver 的 RMS 和空间变化；
6. heat-only approximation 相对精确解的误差

\[
S_{\rm nonlinear}=
\frac{\operatorname{RMSE}(u^*-u_{\rm heat})}
{\operatorname{Std}(u^*)}.
\]

其中 \(u_{\rm heat}=E[g(x+\sqrt{2(T-t)}Z)]\)，对没有闭式表达的情况使用独立高精度 MC，并报告不确定性。

这些指标要真实重新计算，不能直接复制此前讨论中未经正式复核的 Dense Fourier 有效秩数字。

如果发现 LQG 在论文分布下非线性修正极弱，也要保留观察结果，将其标记为外部 benchmark 而不是强非线性成功案例。

所有测试集、参数和预检结果冻结后再进入训练阶段。

## 4. 第二阶段：5070 Ti 小规模 Tangent-Actor pilot

只对 **Dense Fourier 100D 和 LQG 100D** 进行新训练。

B6 使用已有成熟结果作为 anchor，避免重复花算力。

对两个新 PDE 统一使用固定轻量配置：

| 参数 | Pilot 配置 |
|---|---|
| Dimension | 100 |
| Training states | 2,048 |
| Validation states | 256 |
| Paths per training label | 64 |
| Label time steps | 16 |
| Neural network | 3 hidden layers × width 256, SiLU |
| Optimizer | AdamW |
| Optimizer steps / update | 1,500 |
| Learning rate | 1.5e-4 |
| Batch size | 256 |
| Actor updates | 2 |
| Seed | 0 |
| Gradient clipping | 10 |
| Inference | 512 states, 1,024 paths, 64 steps |
| MC evaluation replicates | 2 independent replicates |

不要进行网格搜索。

保存：

- a0：untrained terminal tangent；
- a1：一次 learned update；
- a2：两次 learned update。

使用 reference-free validation/selector 选择 checkpoint，不能依据精确解挑选。

训练及推理全部计时，包括 label generation、optimizer、selector 和 FK evaluation。

评价 a0、a1、a2 时尽量使用配对 test bank 与 common random numbers，另保留独立 MC replication，以估计数值误差。

如果总 RMSE 已接近 MC noise floor，不要直接判定学习无效。可以分析两份已保存独立预测的 pooling 效果，但标记为诊断分析，不擅自增加训练预算。

同时报告 direct Actor head errors 和经过 linear FK 后的 solution errors，避免混淆。

## 5. GPU 时间、安全和自动停止

5070 Ti 显存约 16 GiB。运行时自动检查显存和可用 GPU；不得依赖假定的 GH200 环境。

建议预算：

- 每个新 PDE pilot 最多 60 分钟 GPU 运行时间；
- 全部新 pilot 最多 3 GPU 小时，包括测试、失败重试和推理；
- 同时只运行一个 GPU-heavy process。

支持小 batch / path chunk，遇到 CUDA OOM 可以按预先声明的规则缩小 chunk（不改变总 MC 样本数），并完整记录。

如果基础数学验证失败，不进入训练。

如果两个 PDE 的数据和数学验证都通过，但 GPU 预算不足，优先 Dense Fourier 100D，然后 LQG 100D。

当全部 pilot 完成或者达到预算上限时停止。不自动扩展到 200D、400D，也不自行增加训练 seeds。

## 6. 暂时不做的事情

本阶段禁止：

- FH-MLP 的大网格计算；
- 重新训练 B6；
- 重训旧七个 PDE；
- 新的超参数搜索；
- 400D 大规模训练；
- 长时间 DPI / SDGD / PINN-HTE 训练；
- 未经授权的额外研究方向。

可以检查 DPI、PINN-HTE、SDGD 的现有接口是否支持三个新 PDE，做 CPU/small-batch compatibility smoke test。正式强基线留到下一阶段。

FH-MLP 将来放在单独的 Recursive Monte Carlo Ablation / Accuracy–Compute 小节，不纳入本次新增主实验 pilot。

## 7. 结果与图表

生成：

`invariant_region_mlp/results/main_highdim_pilot/`

至少包含：

- `validation.json`
- `benchmark_calibration.csv`
- `pilot_metrics.csv`
- `actor_iteration_metrics.csv`
- `resource_accounting.csv`
- `test_bank_manifest.json`
- `figures/actor_iteration_rmse.png`
- `figures/effective_dimension.png`
- `figures/nonlinearity_strength.png`
- `figures/pilot_accuracy_compute.png`

报告：

`invariant_region_mlp/docs/MAIN_HIGHDIM_PILOT_5070_REPORT.md`

需要一张主摘要表：

| PDE | d | Effective rank | Nonlinear skill | a0 RMSE | a1 RMSE | a2 RMSE | MC contribution | Total GPU time |
|---|---|---|---|---|---|---|---|---|

B6 的 a4 可另列，但明确说明使用 mature training budget，不与新 pilot 的两轮缩减预算直接比较时间优劣。

## 8. 最终科学判断

每个 PDE 最后给出：

- `GO_MAIN`: 数学 reference 可信、有效维度足够高、非线性项影响非平凡、算法能稳定运行。
- `GO_SUPPORTING`: 数学与算法通过，但有效维度或非线性强度不足以担任主要高维困难案例。
- `NEEDS_METHOD_FIX`: 算法存在实现或数值稳定问题，需要修复后再决定。
- `NO_GO`: 数学基准无法可靠验证，或有明确不适合本研究目标的问题。

这些是 pilot 的科学适用性判断，不能修改已有 H1–H6 预注册结果。

重点回答：

1. Dense Fourier 是否真的体现 joint \(u,z\) 非线性的困难？
2. 它是否足够高有效维度，是否存在明显的 Actor 学习收益？
3. Published LQG 在论文分布下是否真正困难，还是更接近 heat-only？
4. B6、Dense Fourier、LQG 是否可以互补，构成比旧七 PDE 更强的 ICML 主实验？
5. 下一阶段应该在哪些维度、用哪些强基线比较？给出具体预估 GPU 小时，但不要自动启动。

## 9. 交付要求

完成后：

1. 运行所有相关 correctness tests。
2. 核对结果和原始预测数组，报告 numerical reference uncertainty。
3. 保存完整代码、PDE 参数、随机种子、checkpoint、输入 bank hash 和 GPU 环境。
4. Git commit 并尝试 push 到独立分支。
5. 将真正完成的实验与失败、未运行、预算截断的实验严格区分。
6. 用一页研究级总结向用户报告结果和下一阶段建议。

**不要只交一个 plan。请实际实现、验证和运行这个 5070 Ti pilot。完成后停止，等待用户决定是否扩展为正式主实验。**

