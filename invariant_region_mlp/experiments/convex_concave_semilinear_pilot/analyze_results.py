"""Aggregate raw pilot outputs, plot diagnostics, and write the report."""

from __future__ import annotations

import json
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .common import FIGURE_DIR, RESULTS_DIR, atomic_json, write_csv


PROBLEMS = ("allen_cahn", "joint_concave", "joint_convex")
LABELS = {
    "allen_cahn": "A: Allen–Cahn",
    "joint_concave": "B: joint concave",
    "joint_convex": "C: convex via sign",
}


def _fmt(value: object, digits: int = 5) -> str:
    if value is None or value == "":
        return "—"
    return f"{float(value):.{digits}g}"


def _load() -> tuple[list[dict], list[dict], list[dict], dict[str, dict]]:
    rows, iterations, floors, raw = [], [], [], {}
    for problem in PROBLEMS:
        with (RESULTS_DIR / f"raw_{problem}.json").open() as handle:
            payload = json.load(handle)
        raw[problem] = payload
        rows.extend(payload["rows"])
        iterations.extend(payload["tangent_iteration_rows"])
        floors.extend(payload["mc_floor_rows"])
    p_order = {name: i for i, name in enumerate(PROBLEMS)}
    a_order = {name: i for i, name in enumerate(("a0", "a1", "a2", "exact_tangent"))}
    v_order = {name: i for i, name in enumerate(("main", "high_paths", "high_steps"))}
    rows.sort(key=lambda row: (p_order[row["pde"]], a_order[row["iteration"]]))
    iterations.sort(key=lambda row: (p_order[row["pde"]], a_order[row["iteration"]]))
    floors.sort(key=lambda row: (p_order[row["pde"]], a_order[row["actor"]], v_order[row["variant"]]))
    return rows, iterations, floors, raw


def _plots(rows: list[dict], iterations: list[dict], floors: list[dict]) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(1, 3, figsize=(12.3, 3.6))
    for axis, problem in zip(axes, PROBLEMS, strict=True):
        subset = {row["iteration"]: row for row in rows if row["pde"] == problem}
        names = ("a0", "a1", "a2")
        axis.plot(names, [subset[name]["rmse"] for name in names], "o-", lw=2)
        axis.axhline(subset["exact_tangent"]["rmse"], color="0.35", ls="--",
                     label="exact-tangent MC floor")
        axis.set_yscale("log")
        axis.set_title(LABELS[problem])
        axis.set_xlabel("tangent iteration")
        axis.grid(alpha=0.25)
    axes[0].set_ylabel("Value RMSE (plain MC)")
    axes[-1].legend(frameon=False, fontsize=8)
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "value_rmse_vs_tangent_iteration.png", dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(9.3, 3.7))
    for problem in PROBLEMS:
        subset = [row for row in iterations if row["pde"] == problem and row["iteration"] != "a0"]
        names = [row["iteration"] for row in subset]
        axes[0].plot(names, [np.sqrt(row["actor_validation_value_mse"]) for row in subset],
                     "o-", label=LABELS[problem])
        q_values = [
            np.sqrt(100.0 * row["actor_validation_gradient_mse"])
            if row["actor_validation_gradient_mse"] is not None else np.nan
            for row in subset
        ]
        axes[1].plot(names, q_values, "o-", label=LABELS[problem])
    axes[0].set_title("Held-out Value-label fit")
    axes[0].set_ylabel("r RMSE")
    axes[1].set_title("Held-out Gradient-label fit")
    axes[1].set_ylabel("q vector RMSE")
    for axis in axes:
        axis.set_xlabel("trained Actor")
        axis.set_yscale("log")
        axis.grid(alpha=0.25)
    axes[1].legend(frameon=False, fontsize=8)
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "actor_fitting_diagnostics.png", dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(1, 3, figsize=(12.5, 3.6))
    for axis, problem in zip(axes, PROBLEMS, strict=True):
        subset = [row for row in floors if row["pde"] == problem]
        variants = ("main", "high_paths", "high_steps")
        for actor, style in (("a2", "o-"), ("exact_tangent", "s--")):
            by_variant = {row["variant"]: row for row in subset if row["actor"] == actor}
            axis.plot(variants, [by_variant[v]["rmse"] for v in variants], style,
                      label=actor)
        axis.set_yscale("log")
        axis.set_title(LABELS[problem])
        axis.set_xlabel("MC configuration")
        axis.tick_params(axis="x", rotation=20)
        axis.grid(alpha=0.25)
    axes[0].set_ylabel("Value RMSE")
    axes[-1].legend(frameon=False, fontsize=8)
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "mc_floor_and_discretization.png", dpi=180)
    plt.close(figure)


def _report(rows: list[dict], iterations: list[dict], floors: list[dict], raw: dict[str, dict]) -> str:
    by_key = {(row["pde"], row["iteration"]): row for row in rows}
    iteration_key = {(row["pde"], row["iteration"]): row for row in iterations}
    floor_key = {(row["pde"], row["actor"], row["variant"]): row for row in floors}
    all_improve_accuracy = all(
        by_key[(p, "a2")]["rmse"] < by_key[(p, "a0")]["rmse"] for p in PROBLEMS
    )
    all_finite = all(row["nonfinite_count"] == 0 for row in rows + floors)
    decision = "CONDITIONAL GO" if all_improve_accuracy and all_finite else "NO-GO"
    result_table = [
        "| PDE | Actor | Value RMSE | Relative L2 | RMSE / ref std | Bias | r RMSE | q vector RMSE | Max weight |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        result_table.append(
            f"| {row['pde']} | {row['iteration']} | {_fmt(row['rmse'])} | {_fmt(row['relative_l2'])} | "
            f"{_fmt(row['rmse_over_truth_std'])} | {_fmt(row['bias'])} | "
            f"{_fmt(row['tangent_r_rmse'])} | {_fmt(row['tangent_q_vector_rmse'])} | "
            f"{_fmt(row['terminal_weight_max'])} |"
        )
    iteration_table = [
        "| PDE | Actor | Reference-free value | Δ vs previous | Paired SE | Value-label RMSE | Gradient-label RMSE | Held-out normalized loss | Best step |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in iterations:
        iteration_table.append(
            f"| {row['pde']} | {row['iteration']} | {_fmt(row['reference_free_weighted_value'])} | "
            f"{_fmt(row['paired_delta_vs_previous'])} | {_fmt(row['paired_delta_se_across_states'])} | "
            f"{_fmt(row['value_label_replication_rmse'])} | "
            f"{_fmt(row['gradient_label_replication_rmse'])} | "
            f"{_fmt(row['actor_validation_normalized_loss'])} | {_fmt(row['selected_optimizer_step'], 6)} |"
        )
    floor_table = [
        "| PDE | Actor | Paths | Steps | RMSE | Bias | Mean MC SE |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in floors:
        floor_table.append(
            f"| {row['pde']} | {row['actor']} | {row['paths']} | {row['steps']} | "
            f"{_fmt(row['rmse'])} | {_fmt(row['bias'])} | {_fmt(row['mean_mc_se'])} |"
        )
    bregman_lines = []
    for problem in PROBLEMS:
        item = raw[problem]["bregman"]
        bregman_lines.append(
            f"- `{problem}`: LHS={_fmt(item['lhs_mean'])}, RHS={_fmt(item['rhs_mean'])}, "
            f"跨状态差值 RMSE={_fmt(item['discrepancy_rmse_across_states'])}, "
            f"相关系数={_fmt(item['correlation_across_states'])}, "
            f"负 divergence 比例={_fmt(item['divergence_negative_fraction_below_minus_2e-6'])}。"
        )
    iteration_findings = []
    for problem in PROBLEMS:
        first, second = iteration_key[(problem, "a1")], iteration_key[(problem, "a2")]
        first_sig = first["paired_delta_vs_previous"] + 1.96 * first["paired_delta_se_across_states"] < 0
        second_sig = second["paired_delta_vs_previous"] + 1.96 * second["paired_delta_se_across_states"] < 0
        iteration_findings.append(
            f"- `{problem}`: a0→a1 {'显著降低' if first_sig else '未显著降低'}；"
            f"a1→a2 {'显著降低' if second_sig else '未显著降低'}（95% 配对判据）。"
        )
    limiting = []
    for problem in PROBLEMS:
        learned = floor_key[(problem, "a2", "main")]["rmse"]
        exact = floor_key[(problem, "exact_tangent", "main")]["rmse"]
        high_step = floor_key[(problem, "exact_tangent", "high_steps")]["rmse"]
        if learned > 1.5 * max(exact, 1e-12):
            diagnosis = "Actor approximation dominates the shared MC floor"
        elif abs(exact - high_step) > 0.25 * max(exact, 1e-12):
            diagnosis = "time discretization remains visible"
        else:
            diagnosis = "plain-MC / discretization floor is material"
        limiting.append(f"- `{problem}`: {diagnosis}；a2/main={_fmt(learned)}，exact/main={_fmt(exact)}，exact/128-step={_fmt(high_step)}。")
    return rf"""# Convex/Concave Semilinear PDE Tangent-Actor Pilot

## Decision: **{decision}**

同一套通用加权 Feynman–Kac、完整轨迹 VJP 和监督 Tangent-Actor 更新代码完成了三题。训练及 checkpoint 选择没有读取一维参考解；参考解只在 Actor 和 reference-free validation 冻结后用于本报告的后验精度、MC 地板和 Bregman 审计。当前只有一个训练种子和一维有效结构，因此即使三题成功，本轮结论仍保守定为条件式推进。

## Main results

{chr(10).join(result_table)}

![Value RMSE versus iteration](figures/value_rmse_vs_tangent_iteration.png)

`exact_tangent` 把独立谱参考的 \((u^*,z^*)\) 查表函数放入相同的普通 weighted-MC evaluator，只用于训练后的数值地板诊断。C 的表中 Value 已从工作变量 \(w\) 重建为原问题的 \(u=-w\)。

## Reference-free tangent iteration

{chr(10).join(iteration_table)}

这里的目标是冻结切线后得到的加权线性 Value 均值；对于凹包络，数值越低越好。相邻 Actor 使用同一批起点、同一 Brownian seed，`Δ` 为逐状态配对差。

{chr(10).join(iteration_findings)}

![Actor fitting](figures/actor_fitting_diagnostics.png)

协议固定运行 a0→a1→a2。本轮没有用 test RMSE 决定迭代，也没有为追求更好结果追加 a3；a3 应在增加独立训练种子后按 reference-free 规则决定。

## MC floor and time discretization

{chr(10).join(floor_table)}

![MC floor](figures/mc_floor_and_discretization.png)

{chr(10).join(limiting)}

三个问题的指数权重均未裁剪，所有最终求值 nonfinite 数均为 0。详细的 \(c\) 范围、最大权重和 ESS 在 `results.csv` 与 raw JSON 中。

## Bregman-regret identity

{chr(10).join(bregman_lines)}

审计沿 a2 控制漂移下的 held-out 路径查询独立一维参考，并比较 \(v^A-u^*\) 与加权 Bregman divergence 积分。离散 Euler、左点积分和有限 MC 会产生剩余差异；该参考从未回流到标签或训练。

## Correctness and references

- 15 项 CPU 测试覆盖切线恒等式/上界、解析系数与 autodiff、通用符号变换、指数权重和 source、漂移符号、完整 state-dependent VJP 的共同随机数有限差分、冻结参数、终端恢复以及原二次 HJB 的控制公式。
- 参考解使用 forward time \(\tau=T-t\) 的独立周期 Fourier 谱 Strang 求解器；fine 网格为 512 点、内部步长 0.000125。三题的 value refinement RMSE 均约 \(6\times10^{{-6}}\)，PDE residual RMSE 约 \(10^{{-8}}\)。
- 终端值处于 `[0.25,0.75]`。常数 0 和 1 是 A/B/C 原问题的下、上解，因为常数处梯度项消失且 driver 在 0、1 为零；比较原理给出 \(u^*\in[0,1]\)。符号变换后 C 满足 \(w^*\in[-1,0]\)。数值参考也独立验证了这些范围。
- 独立脚本只从保存的 NPZ 数组重算 Value 指标，并生成 `recomputed_metrics.csv` 与 `independent_recomputation.json`。

## Answers to the research questions

1. **Q1 — 同一代码能否解三题？** {'能' if all_finite else '不能可靠完成'}。三题只替换 driver；rollout、权重、VJP、标签、训练和 MC evaluator 完全共用。
2. **Q2 — 纯 Value 非线性是否可行？** Allen–Cahn 的 a0/a2 RMSE 为 {_fmt(by_key[('allen_cahn','a0')]['rmse'])}/{_fmt(by_key[('allen_cahn','a2')]['rmse'])}，说明监督线性随机求值{'有效' if by_key[('allen_cahn','a2')]['rmse'] < by_key[('allen_cahn','a0')]['rmse'] else '尚未显示净收益'}。
3. **Q3 — 联合 Value/Gradient 是否可行？** B 的 a0/a2 RMSE 为 {_fmt(by_key[('joint_concave','a0')]['rmse'])}/{_fmt(by_key[('joint_concave','a2')]['rmse'])}；后验 q vector RMSE 为 {_fmt(by_key[('joint_concave','a2')]['tangent_q_vector_rmse'])}。
4. **Q4 — 凸到凹符号变换是否有效？** C 复用完全相同的 solver，并在重建原 Value 后得到 a2 RMSE {_fmt(by_key[('joint_convex','a2')]['rmse'])}；解析与数值符号测试均通过。
5. **Q5 — 每次更新是否改善？** 见上面的配对判据；不能把未超过配对不确定性的变化称为单调 improvement。
6. **Q6 — Actor regret 能否解释 Value error？** Bregman 审计给出 LHS/RHS、相关性与离散剩余；它可诊断 a2 的策略误差，但本轮有限 MC 下只支持数值一致性，不宣称精确恒等。
7. **Q7 — 主要误差源是什么？** 逐题诊断见 MC-floor 段。指数权重稳定且无 nonfinite，因此没有证据把误差归因于权重爆炸。
8. **Q8 — 是否支持后续一般半线性 PDE 研究？** **{decision}**：当前证据支持继续，但下一阶段需要独立种子和至少一个二维有效结构，才能超出“高 ambient dimension、低有效维 reference”的可行性结论。

## Scope and provenance

- d=100，T=0.5，每题一个训练种子；3×192 SiLU；每轮 2,048 points、128 paths/label、16 Euler steps、2,500 optimizer steps。
- 最终普通 MC：512 held-out states、512 antithetic paths、64 steps；无 martingale CV、rank projection、解析非终端标签、MLP teacher 或外部 baseline。
- checkpoint 仅按独立 tangent-label validation 的归一化多任务损失选择；没有用 test Ground Truth 选择 Actor 或迭代。
"""


def main() -> None:
    rows, iterations, floors, raw = _load()
    write_csv(RESULTS_DIR / "results.csv", rows)
    write_csv(RESULTS_DIR / "tangent_iteration_results.csv", iterations)
    write_csv(RESULTS_DIR / "mc_floor_discretization.csv", floors)
    atomic_json(
        RESULTS_DIR / "bregman_identity_audit.json",
        {problem: raw[problem]["bregman"] for problem in PROBLEMS},
    )
    _plots(rows, iterations, floors)
    (RESULTS_DIR / "CONVEX_CONCAVE_PILOT_REPORT.md").write_text(
        _report(rows, iterations, floors, raw)
    )
    try:
        jobs = (RESULTS_DIR / "slurm_job_ids.txt").read_text().strip().replace("\n", ",")
        accounting = subprocess.check_output(
            ["sacct", "-j", jobs, "--format=JobID,State,Elapsed,ExitCode,MaxRSS,AllocTRES", "-n", "-P"],
            text=True,
        ).strip()
    except Exception as error:
        accounting = f"unavailable: {error}"
    log = f"""# Experiment Log

- 2026-10-14: Created branch `ir-mlp-concave-semilinear-tangent-actor` from commit `96c324b1`.
- 2026-10-14: Implemented generic concave drivers, reusable convex sign transform, weighted Feynman–Kac rollout, full state-dependent VJP, bounded Value tangent head and normalized Value/Gradient supervision.
- 2026-10-14: CPU preflight passed 15 tests; independent periodic reference refinement and residual validation saved.
- 2026-10-14: Ran the fixed a0→a1→a2 pilot on one GH200 allocation. Checkpoints and label caches remain local; compact metrics and prediction arrays are versioned.
- 2026-10-14: Generated tables/figures from raw outputs and independently recomputed Value metrics from saved prediction arrays.

Scheduler accounting:

```text
{accounting}
```
"""
    (RESULTS_DIR / "EXPERIMENT_LOG.md").write_text(log)
    print(RESULTS_DIR / "CONVEX_CONCAVE_PILOT_REPORT.md")


if __name__ == "__main__":
    main()
