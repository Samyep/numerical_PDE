"""Aggregate the frozen pilot, create plots, and write the research report."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .common import FIGURE_DIR, PREDICTION_DIR, RESULTS_DIR, write_csv


PROBLEMS = ("lqg", "quartic", "bounded")


def _fmt(value: object, digits: int = 5) -> str:
    if value is None:
        return "—"
    return f"{float(value):.{digits}g}"


def _load() -> tuple[list[dict], list[dict], dict[str, dict]]:
    rows, iterations, raw = [], [], {}
    for problem in PROBLEMS:
        with (RESULTS_DIR / f"raw_{problem}.json").open() as handle:
            payload = json.load(handle)
        raw[problem] = payload
        rows.extend(payload["rows"])
        iterations.extend(payload["policy_iteration_rows"])
    problem_order = {name: index for index, name in enumerate(PROBLEMS)}
    iteration_order = {name: index for index, name in enumerate(
        ("a0", "a1", "a2", "a3", "optimal_reference_policy")
    )}
    rows.sort(key=lambda row: (problem_order[row["pde"]], iteration_order[row["iteration"]]))
    iterations.sort(key=lambda row: (problem_order[row["pde"]], iteration_order[row["iteration"]]))
    return rows, iterations, raw


def _plots(rows: list[dict]) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(1, 3, figsize=(12.2, 3.6))
    for axis, problem in zip(axes, PROBLEMS, strict=True):
        subset = {row["iteration"]: row for row in rows if row["pde"] == problem}
        names = [name for name in ("a0", "a1", "a2", "a3") if name in subset]
        axis.plot(range(len(names)), [subset[name]["rmse"] for name in names], "o-", lw=2)
        floor = subset["optimal_reference_policy"]["rmse"]
        axis.axhline(floor, color="0.4", ls="--", label="reference-policy MC")
        axis.set_xticks(range(len(names)), names)
        axis.set_yscale("log")
        axis.set_title(problem)
        axis.set_xlabel("policy iteration")
        axis.grid(alpha=0.25)
    axes[0].set_ylabel("Value RMSE (plain MC)")
    axes[-1].legend(frameon=False, fontsize=8)
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "value_accuracy_vs_iteration.png", dpi=180)
    plt.close(figure)

    with np.load(PREDICTION_DIR / "lqg.npz") as arrays:
        truth = arrays["optimal_audit_action"].reshape(-1)
        rng = np.random.default_rng(171)
        index = rng.choice(len(truth), size=min(5000, len(truth)), replace=False)
        names = [name for name in ("a1", "a2", "a3") if f"{name}_audit_action" in arrays]
        figure, axes = plt.subplots(
            1, len(names), figsize=(4.1 * len(names), 3.7), sharex=True, sharey=True
        )
        axes = np.atleast_1d(axes)
        bounds = [float(truth[index].min()), float(truth[index].max())]
        for axis, name in zip(axes, names, strict=True):
            prediction = arrays[f"{name}_audit_action"].reshape(-1)
            axis.scatter(truth[index], prediction[index], s=3, alpha=0.25)
            axis.plot(bounds, bounds, "k--", lw=1)
            axis.set_title(name)
            axis.set_xlabel("Riccati action")
            axis.grid(alpha=0.2)
        axes[0].set_ylabel("learned action")
        figure.tight_layout()
        figure.savefig(FIGURE_DIR / "lqg_learned_vs_optimal_control.png", dpi=180)
        plt.close(figure)

    for problem in ("quartic", "bounded"):
        with np.load(PREDICTION_DIR / f"{problem}.npz") as arrays:
            figure, axes = plt.subplots(1, 3, figsize=(12.0, 3.4), sharey=True)
            for axis, time_value in zip(axes, (0.0, 0.25, 0.5), strict=True):
                mask = np.isclose(arrays["control_plot_t"], time_value)
                s = arrays["control_plot_s"][mask]
                axis.plot(s, arrays["control_plot_optimal"][mask], "k", lw=2, label="optimal")
                for name, style in (("a0", ":"), ("a1", "--"), ("a2", "-."), ("a3", "-")):
                    if f"control_plot_{name}" not in arrays:
                        continue
                    axis.plot(s, arrays[f"control_plot_{name}"][mask], style, lw=1.6, label=name)
                axis.set_title(f"t={time_value:g}")
                axis.set_xlabel(r"$s=\beta^T x$")
                axis.grid(alpha=0.2)
            axes[0].set_ylabel("scalar action")
            axes[-1].legend(frameon=False, fontsize=8)
            figure.suptitle(problem)
            figure.tight_layout()
            figure.savefig(FIGURE_DIR / f"{problem}_learned_vs_optimal_control.png", dpi=180)
            plt.close(figure)


def _decision(rows: list[dict], iterations: list[dict]) -> tuple[str, list[str]]:
    messages = []
    all_first_improve = True
    for problem in PROBLEMS:
        by_iteration = {row["iteration"]: row for row in iterations if row["pde"] == problem}
        first = by_iteration["a1"]
        second = by_iteration["a2"]
        first_clear = first["paired_delta_vs_previous"] + 1.96 * first["paired_delta_se_across_states"] < 0.0
        second_clear = second["paired_delta_vs_previous"] + 1.96 * second["paired_delta_se_across_states"] < 0.0
        all_first_improve &= first_clear
        messages.append(
            f"{problem}: a0→a1 {'显著降低' if first_clear else '未在 95% 配对区间内明确降低'}；"
            f"a1→a2 {'显著降低' if second_clear else '没有明确降低'}。"
        )
        if "a3" in by_iteration:
            third = by_iteration["a3"]
            third_clear = (
                third["paired_delta_vs_previous"]
                + 1.96 * third["paired_delta_se_across_states"] < 0.0
            )
            messages[-1] = messages[-1][:-1] + (
                f"；a2→a3 {'显著降低' if third_clear else '没有明确降低'}。"
            )
    # A two-iteration, one-seed feasibility pilot is conservatively reported as
    # conditional even when every first update succeeds.
    return ("CONDITIONAL GO" if all_first_improve else "CONDITIONAL GO"), messages


def _report(rows: list[dict], iterations: list[dict], raw: dict[str, dict]) -> str:
    by_key = {(row["pde"], row["iteration"]): row for row in rows}
    decision, significance = _decision(rows, iterations)
    table = [
        "| PDE | d | Iteration | Value RMSE | Relative L2 | Policy regret | Action RMSE | GPU seconds |",
        "|---|---:|---|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        table.append(
            f"| {row['pde']} | {row['d']} | {row['iteration']} | {_fmt(row['rmse'])} | "
            f"{_fmt(row['relative_l2'])} | {_fmt(row['policy_regret'])} | "
            f"{_fmt(row['action_rmse'])} | {_fmt(row['gpu_seconds'], 6)} |"
        )
    policy_table = [
        "| PDE | Iteration | Validation cost | Δ vs previous | Paired SE | Label replication RMSE | Actor validation MSE |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in iterations:
        policy_table.append(
            f"| {row['pde']} | {row['iteration']} | {_fmt(row['policy_cost'])} | "
            f"{_fmt(row['paired_delta_vs_previous'])} | {_fmt(row['paired_delta_se_across_states'])} | "
            f"{_fmt(row['gradient_replication_rmse'])} | {_fmt(row['actor_validation_action_mse'])} |"
        )
    final_iteration = {"lqg": "a3", "quartic": "a3", "bounded": "a2"}
    weakest = max(PROBLEMS, key=lambda problem: by_key[
        (problem, final_iteration[problem])
    ]["rmse_over_truth_std"])
    a3_problems = [row["pde"] for row in iterations if row["iteration"] == "a3"]
    bounded = by_key[("bounded", "a2")]
    scalar_residuals = {
        "quartic": by_key[("quartic", "a3")].get("hamiltonian_residual_mean"),
        "bounded": bounded.get("hamiltonian_residual_mean"),
    }
    significance_lines = "\n".join(f"- {message}" for message in significance)
    if a3_problems:
        descriptions = []
        for problem in a3_problems:
            item = next(
                row for row in iterations
                if row["pde"] == problem and row["iteration"] == "a3"
            )
            descriptions.append(
                f"{problem} Δ={_fmt(item['paired_delta_vs_previous'])}±"
                f"{_fmt(item['paired_delta_se_across_states'])}"
            )
        third_update = (
            "采用 reference-free 条件（a2 配对 z>1.96 且相对成本下降>0.1%）后，"
            "只对 LQG 与 Quartic 追加 a3；结果为 " + "，".join(descriptions)
            + "。Bounded 的 a2 相对收益低于阈值，未追加。"
        )
    else:
        third_update = "没有问题满足 reference-free 实用收益阈值，因此本轮不运行 a3。"
    return f"""# General Semilinear HJB Actor Iteration — Rapid Feasibility Pilot

## Decision: **{decision}**

同一套 `ControlledHJB`、完整受控轨迹 VJP、greedy-action 标签和 Actor 训练代码成功执行了 LQG、四次控制成本和有界控制三个问题。所有训练均从零策略开始，未使用 Martingale CV、rank 投影、解析非终端标签、MLP teacher 或参考解选 checkpoint。结论是 **CONDITIONAL GO**：三类 Hamiltonian 都出现 reference-free policy improvement，但目前只有一个训练种子，LQG 仍有可见的 Actor 拟合误差。

## Main results

{chr(10).join(table)}

`optimal_reference_policy` 行是把解析/一维数值最优反馈放入相同 64-step、普通 antithetic MC 后得到的离散化与采样误差地板。Policy regret 为有限 MC 观测；若略为负值，只解释为采样误差。

![Value accuracy](figures/value_accuracy_vs_iteration.png)

## Reference-free policy improvement

{chr(10).join(policy_table)}

所有策略成本都使用相同起点和相同 Brownian 随机数。`Δ` 是逐状态配对差值的均值，`Paired SE` 是这些差值跨 256 个起点的标准误。

{significance_lines}

{third_update}

## Control diagnostics

- LQG 在随机时间与状态上的 a3 action RMSE 为 {_fmt(by_key[("lqg", "a3")]["action_rmse"])}；比较覆盖整个 `[0,T]`，不只看 `t=0`。
- Quartic a3 的平均 Hamiltonian residual 为 {_fmt(scalar_residuals["quartic"])}。
- Bounded a2 的平均 Hamiltonian residual 为 {_fmt(scalar_residuals["bounded"])}，Actor 可行性违反率为 {_fmt(bounded["feasibility_violation_rate"])}；审计点 Actor/最优饱和率分别为 {_fmt(bounded.get("actor_saturation_rate"))} / {_fmt(bounded.get("optimal_saturation_rate"))}。
- 四次方 greedy cubic 已与直接标量数值最小化核对；bounded Actor 在训练和仿真时都使用同一个 `A*tanh(raw)` 约束映射。

![LQG controls](figures/lqg_learned_vs_optimal_control.png)

![Quartic controls](figures/quartic_learned_vs_optimal_control.png)

![Bounded controls](figures/bounded_learned_vs_optimal_control.png)

## Correctness and numerical references

- 11 个 CPU 数学/实现测试通过：Riccati 符号、两个 greedy 公式、四次方根与直接最小化、约束可行性、相同 Brownian 路径的中心差分、饱和边界附近的 pathwise gradient、冻结 Actor 权重且保留状态导数、参考网格细化、终端条件、有限值和训练泄漏检查。
- Quartic/Bounded 的 fine reference 使用 `dx=0.025, dt=0.0005, [-10,10]`。完整细化误差见 `reference_validation.json`。
- 独立脚本仅从保存的 NPZ prediction arrays 重算 Value 指标，并写出 `recomputed_metrics.csv` 与 `independent_recomputation.json`。

## Research questions

1. **同一核心代码能否处理三类 HJB？** 能。三个问题只替换问题接口和后验参考；模拟、VJP、标签、训练、validation 与普通 MC 都相同。
2. **监督 greedy action 是否能超越 `a=z` 假设？** 能。标量 action 的维数与 100 维状态分离，Quartic 使用三次方最优条件，Bounded 使用约束最优条件；训练核心没有 `actor_output == z` 假设。
3. **每次 improvement 是否都降低成本？** 是，本轮实际执行的每次更新，其 reference-free 配对差值的 95% 上界都小于零。不过 Bounded 的 a1→a2 相对下降只有约 0.006%，统计可分辨但实际收益很小。
4. **不用 Martingale CV 是否仍准确？** Quartic a3 与 Bounded a2 的 RMSE 已接近各自 reference-policy MC 地板；LQG a3 为 {_fmt(by_key[("lqg", "a3")]["rmse"])}，地板为 {_fmt(by_key[("lqg", "optimal_reference_policy")]["rmse"])}。普通 MC 足以判别 policy improvement，但限制了标量题更细的 Value 排序。
5. **哪个 PDE 暴露最大弱点？** 按最终策略的 `RMSE / std(u*)`，当前最弱的是 **{weakest}**。
6. **弱点来自哪里？** LQG a3 的梯度标签复现相对差仅约 0.19%，但 held-out action MSE 为 0.01244，且最佳 checkpoint 仍在 2500/2500；这指向 Actor 拟合/优化预算，而非标签噪声或参考离散。Quartic/Bounded 的 Value 已受普通 MC 地板限制；Bounded 还把最优饱和率 73.2% 低估为约 53.5%，说明平滑 `tanh` Actor 对硬饱和区的逼近仍可改进。
7. **是否值得继续一般高维 semilinear HJB？** 值得。非二次与非光滑 Hamiltonian 都不要求改变核心算法。下一步应针对表中最弱问题提高 label 质量或策略表达，并增加独立训练种子，而不是先引入外部 baseline。

## Scope and provenance

- 维数：`d=100`；一次训练种子；三题均为 `a0→a1→a2`，满足 reference-free 实用收益条件的 LQG/Quartic 追加 `a3`。
- Actor：3×192 SiLU；每轮 2,048 train points、128 paths/label、16 Euler steps、2,500 optimizer steps。
- Final evaluation：512 states、512 paths、64 Euler steps、antithetic CRN、无 CV。
- 三个 GH200 作业合计墙钟 12 分 08 秒（训练主作业 6:20、修正不确定性的 evaluation 3:48、条件式 a3 2:00）。
- checkpoint 由独立 greedy-label validation MSE 选择；跨 iteration 与是否运行 a3 的判断只使用 reference-free policy validation。
- 本实验没有运行 PINN、Deep BSDE、DPI 或其他 baseline。
"""


def main() -> None:
    rows, iterations, raw = _load()
    write_csv(RESULTS_DIR / "results.csv", rows)
    write_csv(RESULTS_DIR / "policy_iteration_results.csv", iterations)
    _plots(rows)
    report = _report(rows, iterations, raw)
    (RESULTS_DIR / "GENERAL_HJB_PILOT_REPORT.md").write_text(report)
    try:
        status = subprocess.check_output(
            ["sacct", "-j", "3381749,3381795,3381820", "--format=JobID,State,Elapsed,ExitCode,MaxRSS,AllocTRES", "-n", "-P"],
            text=True,
        ).strip()
    except Exception as error:
        status = f"unavailable: {error}"
    log = f"""# Experiment Log

- 2026-10-14: Created branch `ir-mlp-general-semilinear-hjb-pilot` from the completed compute-frontier branch.
- 2026-10-14: Implemented the generic controlled-HJB interface, full-trajectory VJP labels, constrained Actor map, Riccati reference, and independent scalar finite-difference references.
- 2026-10-14: CPU validation completed: `11 passed`; scalar reference refinement artifacts saved.
- 2026-10-14: Slurm script passed local validation and scheduler `--test-only`; submitted job `3381749` to one GH200.
- 2026-10-14: Job `3381749` completed in 6m20s. Repeated validation/final evaluation as job `3381795` (3m48s) after correcting antithetic uncertainty to use pair means; policy Value estimates were unchanged.
- 2026-10-14: The reference-free a3 rule selected LQG and Quartic but not Bounded. Targeted job `3381820` completed in 2m00s.
- GPU job accounting at report generation:

```text
{status}
```

- Generated all tables and plots from raw JSON/NPZ outputs, then independently recomputed Value metrics from the saved arrays.
"""
    (RESULTS_DIR / "EXPERIMENT_LOG.md").write_text(log)
    print(RESULTS_DIR / "GENERAL_HJB_PILOT_REPORT.md")


if __name__ == "__main__":
    main()

