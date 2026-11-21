"""Replay and plot every 725 speed-grid condition with fixed source models."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from scipy.io import loadmat

from audit_subjects import ROOT
from fit_subject import many, record, rollout, score_rows, boundary_torque, SHARED
from prepare_exp_speed_grid import BASE, SUBJECT, LOADS, SPEEDS, excluded_panels

DOCS = ROOT / "docs/exp_speed_grid"
OLD = ROOT / "results/exp_three_loads" / SUBJECT / "neural_models.mat"
COLORS = ("#1469b0", "#8450ad")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replay(path, prepared, generator=None):
    cache = loadmat(path, simplify_cells=True)
    models = many(cache["models"])
    if generator is None:
        generator = prepared["generator"]
    records = [record(s, generator) for s in many(prepared["sequences"])]
    values = [
        [
            rollout(
                m["parameters"],
                m["neuralParameters"],
                cache["controller"],
                r,
                source,
                10,
            )
            for r in records
        ]
        for source, m in enumerate(models)
    ]
    return cache, models, records, values


def summarize(rows):
    result = {}
    for source in ("Experimental activation", "MAP activation"):
        result[source] = {}
        for label, speeds in (
            ("priority_35_55", (35, 45, 55)),
            ("all_25_65", (25, 35, 45, 55, 65)),
        ):
            group = [
                r for r in rows if r["source"] == source and r["speedCode"] in speeds
            ]
            result[source][label] = dict(
                conditions=len(group),
                meanShapeCorrelation=float(
                    np.mean([r["shapeCorrelation"] for r in group])
                ),
                minimumShapeCorrelation=float(
                    min(r["shapeCorrelation"] for r in group)
                ),
                shapeCorrelationAtLeast0p7=sum(
                    r["shapeCorrelation"] >= 0.7 for r in group
                ),
                twoMainValleyConditions=sum(r["predictedValleys"] == 2 for r in group),
                meanThetaRMSEDegrees=float(
                    np.mean([r["thetaRMSEDegrees"] for r in group])
                ),
                meanAvailableTorqueRMSENm=float(
                    np.mean(
                        [
                            r["torqueRMSENm"]
                            for r in group
                            if r["torqueRMSENm"] is not None
                        ]
                    )
                ),
                boundedAngleConditions=sum(
                    r["minAngleRad"] >= -0.2 and r["maxAngleRad"] <= 2.1 for r in group
                ),
                maximumActionBoundaryTorqueNm=max(r["boundaryPeakNm"] for r in group),
                maximumBoundaryTorqueNm=max(
                    r["wholeSequenceBoundaryPeakNm"] for r in group
                ),
                maximumAbdomenBackMeanRatio=max(r["abdomenBackRatio"] for r in group),
                failed=sum(r["failed"] for r in group),
            )
    return result


def plots(data, models, records, values, old_values, rows, output, decisions):
    figures = []
    for source, label in enumerate(
        ("Experimental activation", "Trained MAP activation")
    ):
        for kind in ("theta", "torque", "muscles"):
            fig, axes = plt.subplots(3, 5, figsize=(18, 9), sharey=True)
            excluded_panels(axes, decisions)
            for index, rec in enumerate(records):
                seg = rec["segments"][0]
                ax = axes[LOADS.index(seg["loadKg"]), SPEEDS.index(seg["speedCode"])]
                ids = seg["ids"]
                t = rec["seq"]["t"]
                o = values[source][index][0]
                old = old_values[source][index][0]
                ref = rec["seq"]
                measured = ~np.asarray(ref["waiting"], dtype=bool)
                metric = rows[source * len(records) + index]
                if kind == "theta":
                    ax.plot(
                        t,
                        np.where(measured, ref["q"], np.nan),
                        "--",
                        color="#c13b31",
                        lw=1.7,
                        label="Original-scale IK reference",
                    )
                    ax.plot(
                        t,
                        old[:, 0],
                        color=".65",
                        lw=1.0,
                        label="Previous three-anchor BSE",
                    )
                    ax.plot(
                        t,
                        o[:, 0],
                        color=COLORS[source],
                        lw=1.7,
                        label="Grid-fitted BSE",
                    )
                    note = f"r={metric['shapeCorrelation']:.2f}; {o[:,0].min():.2f}..{o[:,0].max():.2f} rad"
                elif kind == "torque":
                    ax.plot(
                        t,
                        np.where(
                            np.asarray(ref["torqueReferenceValid"], dtype=bool),
                            ref["torque"],
                            np.nan,
                        ),
                        "--",
                        color="#c13b31",
                        lw=1.5,
                        label="CEINMS estimate (uncalibrated)",
                    )
                    ax.plot(
                        t,
                        old[:, 2],
                        color=".65",
                        lw=1.0,
                        label="Previous three-anchor BSE",
                    )
                    ax.plot(
                        t,
                        o[:, 2],
                        color=COLORS[source],
                        lw=1.6,
                        label="Grid-fitted net muscle torque",
                    )
                    valid = np.asarray(ref["torqueReferenceValid"]).astype(bool)
                    ax.fill_between(
                        t,
                        0,
                        1,
                        where=~valid,
                        transform=ax.get_xaxis_transform(),
                        color=".9",
                        zorder=-2,
                    )
                    note = f"Action reference coverage {valid[ids].mean():.0%}"
                else:
                    ax.plot(
                        t, o[:, 3], color="#1469b0", lw=1.6, label="Back muscle torque"
                    )
                    ax.plot(
                        t,
                        o[:, 4],
                        color="#8450ad",
                        lw=1.5,
                        label="Abdominal antagonist torque magnitude",
                    )
                    note = f"Mean abdomen / back = {metric['abdomenBackRatio']:.1%}"
                ax.set_title(
                    f"{seg['loadKg']:g} kg / {seg['speedCode']}\n{note}", fontsize=10
                )
                ax.set_xlim(0, t[-1])
                ax.grid(alpha=0.2)
                ax.axvline(t[ids[0]], color=".6", ls=":", lw=0.7)
                ax.axvline(t[ids[-1]], color=".6", ls=":", lw=0.7)
                ax.set_xlabel("Time with standing waits (s)")
                if 35 <= seg["speedCode"] <= 55:
                    ax.set_facecolor("#fff8e8")
            for ax in axes[:, 0]:
                ax.set_ylabel("Theta (rad)" if kind == "theta" else "Torque (N m)")
            handles, labels = axes.flat[0].get_legend_handles_labels()
            fig.legend(
                handles, labels, loc="lower center", ncol=3, frameon=False, fontsize=10
            )
            fig.suptitle(
                f"725 | {label} | {kind}\nOne fixed parameter vector across {len(records)} retained conditions; yellow = priority 35--55; grey = excluded"
            )
            fig.tight_layout(rect=(0, 0.045, 1, 1))
            filename = f"{('experimental', 'map')[source]}_{kind}.png"
            fig.savefig(output / filename, dpi=150)
            figures.append(fig)
    fig, axes = plt.subplots(3, 5, figsize=(18, 8), sharey=True)
    excluded_panels(axes, decisions)
    for index, rec in enumerate(records):
        seg = rec["segments"][0]
        ax = axes[LOADS.index(seg["loadKg"]), SPEEDS.index(seg["speedCode"])]
        ids = seg["ids"]
        t = rec["seq"]["t"]
        peaks = []
        for source in (0, 1):
            o = values[source][index][0]
            boundary = boundary_torque(models[source]["parameters"], o[:, 0])
            peaks.append(float(np.max(np.abs(boundary))))
            ax.plot(
                t, boundary, color=COLORS[source], label=("Experimental", "MAP")[source]
            )
        ax.axhline(0, color=".55", ls="--", lw=0.6)
        ax.set_title(
            f"{seg['loadKg']:g} kg / {seg['speedCode']}\nPeak E / MAP: {peaks[0]:.2f} / {peaks[1]:.2f} N m",
            fontsize=10,
        )
        ax.axvline(t[ids[0]], color=".6", ls=":", lw=0.7)
        ax.axvline(t[ids[-1]], color=".6", ls=":", lw=0.7)
        ax.set_xlabel("Time with standing waits (s)")
        ax.grid(alpha=0.2)
    for ax in axes[:, 0]:
        ax.set_ylabel("Signed soft-boundary torque (N m)")
    fig.legend(
        *axes.flat[0].get_legend_handles_labels(),
        loc="lower center",
        ncol=2,
        frameon=False,
    )
    fig.suptitle(
        "725 | soft-boundary torque from the same Hill simulation; no MPC assistance"
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(output / "soft_boundary.png", dpi=150)
    figures.append(fig)
    with PdfPages(output / "725_bse_speed_grid.pdf") as pdf:
        for fig in figures:
            pdf.savefig(fig)
            plt.close(fig)


def validate(cache, models, records):
    assert np.array_equal(
        models[0]["parameters"][SHARED], models[1]["parameters"][SHARED]
    )
    checks, standing = [], []
    for source, model in enumerate(models):
        parameter_hash = hashlib.sha256(
            np.r_[model["parameters"], model["neuralParameters"]].tobytes()
        ).hexdigest()
        for rec in records:
            a, failed, _ = rollout(
                model["parameters"],
                model["neuralParameters"],
                cache["controller"],
                rec,
                source,
                10,
            )
            fine, fine_failed, _ = rollout(
                model["parameters"],
                model["neuralParameters"],
                cache["controller"],
                rec,
                source,
                20,
            )
            changed = copy.deepcopy(rec)
            changed["seq"]["q"][1:] += 0.27
            changed["seq"]["torque"] *= 1.7
            independent, _, _ = rollout(
                model["parameters"],
                model["neuralParameters"],
                cache["controller"],
                changed,
                source,
                10,
            )
            checks.append(
                dict(
                    source=source,
                    id=rec["segments"][0]["recordId"],
                    parameterSHA256=parameter_hash,
                    failed=failed or fine_failed,
                    halfStepAngleDifferenceRad=float(np.max(abs(a[:, 0] - fine[:, 0]))),
                    halfStepTorqueDifferenceNm=float(np.max(abs(a[:, 2] - fine[:, 2]))),
                    referenceChangeMaximumOutputDifference=float(
                        np.max(abs(a - independent))
                    ),
                )
            )
        for initial, velocity in ((1.7, 0.6), (1.25, -0.3)):
            n = 3001
            zero = np.zeros(n)
            q = np.full(n, 1.5)
            q[0] = initial
            stress = dict(
                q=q,
                torque=zero,
                a=np.zeros((n, 2)),
                front=zero,
                m=zero,
                loadKg=zero,
                speedCode=zero,
                stageId=zero,
                waiting=np.ones(n),
                referenceIsFilled=np.ones(n),
                t=np.arange(n) * 0.01,
                dt=0.01,
                standingAngle=1.5,
                v0=velocity,
                segments=[],
            )
            stress_rec = record(stress, model["generator"])
            out, bad, _ = rollout(
                model["parameters"],
                model["neuralParameters"],
                cache["controller"],
                stress_rec,
                source,
                10,
            )
            error = float(abs(out[-1, 0] - 1.5))
            standing.append(
                dict(
                    source=source,
                    initialAngleRad=initial,
                    initialVelocity=velocity,
                    finalErrorRad=error,
                    minimumAngleRad=float(out[:, 0].min()),
                    maximumAngleRad=float(out[:, 0].max()),
                    passed=bool(
                        not bad
                        and error < 0.01
                        and out[:, 0].min() > -0.2
                        and out[:, 0].max() < 2.1
                    ),
                )
            )
    result = dict(
        sharedBodyGeometry=True,
        fixedVectorPerSource=True,
        checks=checks,
        standingRecovery=standing,
        passed=all(
            not c["failed"]
            and c["halfStepAngleDifferenceRad"] < 0.02
            and c["halfStepTorqueDifferenceNm"] < 3
            and c["referenceChangeMaximumOutputDifference"] == 0
            for c in checks
        )
        and all(c["passed"] for c in standing),
    )
    return result


def overview(old_rows, rows, output):
    fig, axes = plt.subplots(2, 2, figsize=(12, 6))
    for source, label in enumerate(("Experimental activation", "MAP activation")):
        for col, (values, title) in enumerate(
            ((old_rows, "Previous BSE"), (rows, "Updated BSE"))
        ):
            grid = np.full((3, 5), np.nan)
            for row in values:
                if row["source"] == label:
                    grid[
                        LOADS.index(int(row["loadKg"])), SPEEDS.index(row["speedCode"])
                    ] = row["shapeCorrelation"]
            ax = axes[source, col]
            palette = plt.get_cmap("RdYlGn").copy()
            palette.set_bad(".88")
            plotted = ax.imshow(grid, cmap=palette, vmin=-0.2, vmax=1, aspect="auto")
            ax.set_xticks(range(5), SPEEDS)
            ax.set_yticks(range(3), [f"{load} kg" for load in LOADS])
            ax.set_xlabel("Speed code")
            ax.set_title(f"{label} | {title}")
            for i in range(3):
                for j in range(5):
                    ax.text(
                        j,
                        i,
                        f"{grid[i,j]:.2f}" if np.isfinite(grid[i, j]) else "Excluded",
                        ha="center",
                        va="center",
                        fontsize=10 if np.isfinite(grid[i, j]) else 8,
                    )
    fig.colorbar(
        plotted, ax=axes, label="Theta shape correlation", fraction=0.025, pad=0.03
    )
    fig.suptitle(
        "725 | same 10 retained conditions before / after\nCalibration trend metric only; 5 cases without CEINMS support excluded"
    )
    fig.subplots_adjust(top=0.83, bottom=0.1, hspace=0.6, wspace=0.3, right=0.86)
    fig.savefig(output / "trend_summary.png", dpi=170)
    plt.close(fig)


def write_summary(report, models, output):
    lines = [
        "# 725 BSE 跨速度标定结果",
        "",
        "保留 10 个工况，其中优先速度 35–55 为 6 个；5 个无 torque 支持的工况完全排除。全部结果是标定回放，不是独立泛化验证。",
        "",
        "| 激活来源 | 范围 | 旧平均相关性 | 新平均相关性 | 新最差相关性 | 新 r≥0.7 工况数 |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for source in ("Experimental activation", "MAP activation"):
        for region, label in (
            ("priority_35_55", "35–55"),
            ("all_25_65", "全部保留工况"),
        ):
            old = report["oldSummary"][source][region]
            new = report["summary"][source][region]
            lines.append(
                f"| {source} | {label} | {old['meanShapeCorrelation']:.3f} | {new['meanShapeCorrelation']:.3f} | {new['minimumShapeCorrelation']:.3f} | {new['shapeCorrelationAtLeast0p7']}/{new['conditions']} |"
            )
    lines += [
        "",
        "r≥0.7 是趋势浏览指标，不等于物理验证或完整 W 形验收。角度、肌肉扭矩和软限位曲线需同时检查。",
        "",
        "| 激活来源 | 新角度总范围(rad) | 最大平均腹/背扭矩比 | soft-boundary 峰值(N·m) |",
        "|---|---:|---:|---:|",
    ]
    for source in ("Experimental activation", "MAP activation"):
        rows = [r for r in report["metrics"] if r["source"] == source]
        lines.append(
            f"| {source} | {min(r['sequenceMinAngleRad'] for r in rows):.3f}–{max(r['sequenceMaxAngleRad'] for r in rows):.3f} | {max(r['abdomenBackRatio'] for r in rows):.1%} | {max(r['wholeSequenceBoundaryPeakNm'] for r in rows):.3f} |"
        )
    p = models[0]["parameters"]
    lines += [
        "",
        "| 激活来源 | 旧→新角度 RMSE(°) | 旧→新可用 torque RMSE(N·m) |",
        "|---|---:|---:|",
    ]
    for source in ("Experimental activation", "MAP activation"):
        old = report["oldSummary"][source]["all_25_65"]
        new = report["summary"][source]["all_25_65"]
        lines.append(
            f"| {source} | {old['meanThetaRMSEDegrees']:.2f} → {new['meanThetaRMSEDegrees']:.2f} | {old['meanAvailableTorqueRMSENm']:.2f} → {new['meanAvailableTorqueRMSENm']:.2f} |"
        )
    lines += [
        "",
        "torque RMSE 是逐工况有效区间 RMSE 的等权平均，不包含任何缺失参考样本。腹/背比例仅统计动作段，边界峰值和角度范围统计包含等待的完整序列。",
        "",
    ]
    lines += [
        "",
        f"共享参数：有效躯干 M={p[18]:.3f} kg，身体力臂 L={p[19]:.4f} m，负载力臂={p[12]:.4f} m，负载偏置角={p[13]:.4f} rad。soft boundary 区间 {p[30]:.4f}–{p[31]:.4f} rad，刚度 {p[28]:.3f} N·m/rad，饱和参数 {p[29]:.3f} N·m。这些是有界标定值，不是独立测量的解剖参数。",
        "",
        "[全部逐工况指标](evaluation.json) · [数值验证](validation.json) · [完整图集](725_bse_speed_grid.pdf)",
        "",
        "本轮没有重调 NMPC，不能将旧 NMPC 结果当作这个新 BSE 的控制结果。",
        "",
    ]
    (output / "summary.zh-CN.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=BASE)
    parser.add_argument("--output", type=Path, default=DOCS)
    parser.add_argument("--baseline-only", action="store_true")
    args = parser.parse_args()
    folder = args.base / SUBJECT
    data = loadmat(folder / "prepared.mat", simplify_cells=True)
    gate = json.loads((args.base / "quality_gate.json").read_text())
    assert digest(args.base / "quality_gate.json") == data["qualityGateSHA256"]
    included = {d["id"] for d in gate["trials"] if d["included"]}
    cycles = many(data["cycles"])
    assert {c["id"] for c in cycles} == included
    assert all(np.asarray(c["torqueReferenceValid"], dtype=bool).any() for c in cycles)
    old_cache = loadmat(OLD, simplify_cells=True)
    old_generator = many(old_cache["models"])[0]["generator"]
    old_cache, old_models, old_records, old_values = replay(OLD, data, old_generator)
    old_rows = score_rows(
        old_models, old_cache["controller"], old_records, "previous_three_anchor"
    )
    (folder / "baseline_metrics.json").write_text(json.dumps(old_rows, indent=2) + "\n")
    print("PREVIOUS", json.dumps(summarize(old_rows)), flush=True)
    if args.baseline_only:
        return
    path = folder / "neural_models.mat"
    cache, models, records, values = replay(path, data)
    rows = score_rows(models, cache["controller"], records, "speed_grid")
    args.output.mkdir(parents=True, exist_ok=True)
    plots(data, models, records, values, old_values, rows, args.output, gate["trials"])
    overview(old_rows, rows, args.output)
    validation = validate(cache, models, records)
    (args.output / "validation.json").write_text(
        json.dumps(validation, indent=2) + "\n"
    )
    report = dict(
        subject=SUBJECT,
        modelSHA256=digest(path),
        preparedSHA256=digest(folder / "prepared.mat"),
        baselineModelSHA256=digest(OLD),
        oldSummary=summarize(old_rows),
        summary=summarize(rows),
        conditions=len(records),
        priorityConditions=sum(
            35 <= r["segments"][0]["speedCode"] <= 55 for r in records
        ),
        excludedConditions=[d for d in gate["trials"] if not d["included"]],
        allConditionsAreCalibration=True,
        zeroTorqueSupportCasesEnterFitOrMetrics=False,
        independentValidation=False,
        MAPRefitted=False,
        MPCRefitted=False,
        criteriaNote="Correlation >= 0.7 is an engineering trend indicator, not physiological validation.",
        oldMetrics=old_rows,
        metrics=rows,
        numericalValidationPassed=validation["passed"],
    )
    for dest in (folder, args.output):
        (dest / "evaluation.json").write_text(json.dumps(report, indent=2) + "\n")
    write_summary(report, models, args.output)
    print("NEW", json.dumps(report["summary"]), flush=True)
    print("NUMERICS", validation["passed"], flush=True)


if __name__ == "__main__":
    main()
