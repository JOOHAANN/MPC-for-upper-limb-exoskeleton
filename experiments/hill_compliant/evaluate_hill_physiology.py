"""Compare a fixed parameter pair over all nine original conditions.

Use the unchanged forward Hill solver, exact references and continuous waits.
Screening criteria describe this demo; they are not universal physiological
limits and do not establish an individual subject's true muscle torques.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat, savemat

import fit_hill_physiology as fit
from fit_hill_compliant import ROOT, PARAMETER_NAMES, boundary_torque


def shade(ax, trial):
    t = trial["t"]
    waiting = np.asarray(trial["waiting"], dtype=bool)
    edges = np.diff(np.r_[False, waiting, False].astype(int))
    for start, end in zip(np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)):
        ax.axvspan(t[start], t[min(end, len(t) - 1)], color=".95", zorder=0)
    ax.set_xlim(t[0], t[-1])
    ax.grid(alpha=.25)


def label_loads(ax, trial):
    """Retain original speed codes rather than inventing a BPM conversion."""
    for segment in trial["segments"]:
        first, last = int(segment["startIndex"]) - 1, int(segment["endIndex"]) - 1
        middle = (trial["t"][first] + trial["t"][last]) / 2
        ax.text(middle, .97, f'{int(segment["loadKg"])} kg / code {int(segment["speedCode"])}',
                transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=8,
                bbox={"facecolor": "white", "alpha": .8, "edgecolor": "none", "pad": 2})


def savefig(fig, out, stem):
    for suffix in ("png", "pdf"):
        fig.savefig(out / f"{stem}.{suffix}", dpi=180)
    plt.close(fig)


def plots(out, record, previous, candidate, pair, previous_label="Previous calibration"):
    trial, group = record["trial"], record["group"]
    fig, axes = plt.subplots(2, 2, figsize=(13, 7.5), sharex=True, sharey="row", layout="constrained")
    for source, title in enumerate(("Experimental activation", "MAP activation")):
        old, new = previous[source], candidate[source]
        for row, col, ref, ylabel in ((0, 0, trial["q"], "Angle (rad)"),
                                     (1, 2, trial["torque"], "Net Hill torque (N·m)")):
            ax = axes[row, source]
            ax.plot(trial["t"], ref, color="#d02929", ls="--", label="Original reference")
            ax.plot(trial["t"], old[:, col], color=".5", alpha=.7, label=previous_label)
            ax.plot(trial["t"], new[:, col], color="#0877b7", label="Depth-priority calibration")
            ax.set_ylabel(ylabel)
            shade(ax, trial)
        axes[0, source].set_title(title)
        label_loads(axes[0, source], trial)
        axes[1, source].set_xlabel("Time (s)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=3, fontsize=8)
    fig.suptitle(f"Group {group} | Original Scale retained | One fixed parameter set per source")
    savefig(fig, out, f"torque_theta_group{group}")
    fig, axes = plt.subplots(3, 2, figsize=(13, 11), sharex=True, sharey="row", layout="constrained")
    for source, title in enumerate(("Experimental activation", "MAP activation")):
        old, new = previous[source], candidate[source]
        for row, pred, label in ((0, old, "Previous"), (1, new, "New")):
            ax = axes[row, source]
            ax.plot(trial["t"], pred[:, 3], color="#0877b7", label="Back extensor")
            ax.plot(trial["t"], pred[:, 4], color="#db6500", label="Abdominal flexor")
            ax.plot(trial["t"], pred[:, 2], "k--", label="Net = back - abdomen")
            ax.set_ylabel(f"{label} torque (N·m)")
            ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(0, .88))
            shade(ax, trial)
            label_loads(ax, trial)
        ax = axes[2, source]
        ax.plot(trial["t"], record["base"][source][:, 1], color=".4", label="Front base input")
        ax.plot(trial["t"], old[:, 14], color=".7", label="Previous abdominal Hill activation")
        ax.plot(trial["t"], new[:, 14], color="#db6500", label="New abdominal Hill activation")
        ax.set_ylabel("Activation (not verified %MVC)")
        ax.set_ylim(0, 1)
        ax.legend(fontsize=8)
        ax.set_xlabel("Time (s)")
        shade(ax, trial)
        axes[0, source].set_title(title)
    fig.suptitle(f"Group {group} | Opposing muscle torque magnitudes | Predictions, not independent measurements")
    savefig(fig, out, f"muscle_recruitment_group{group}")
    fig, axes = plt.subplots(2, 2, figsize=(13, 7.5), sharex=True, sharey="row", layout="constrained")
    for source, title in enumerate(("Experimental activation", "MAP activation")):
        new = candidate[source]
        for col, label in ((4, "Total abdominal"), (8, "Active CE"), (9, "Passive PE"), (10, "Fiber damping")):
            axes[0, source].plot(trial["t"], new[:, col], label=label)
        axes[0, source].set_title(title)
        axes[0, source].set_ylabel("Abdominal torque (N·m)")
        axes[0, source].legend(fontsize=8, loc="upper left", bbox_to_anchor=(0, .88))
        label_loads(axes[0, source], trial)
        axes[1, source].plot(trial["t"], new[:, 13], color="#0877b7", label="Back Hill activation")
        axes[1, source].plot(trial["t"], new[:, 14], color="#db6500", label="Abdominal Hill activation")
        axes[1, source].set_ylabel("Model activation")
        axes[1, source].set_xlabel("Time (s)")
        axes[1, source].legend(fontsize=8)
        for row in (0, 1):
            shade(axes[row, source], trial)
    fig.suptitle(f"Group {group} | New abdominal components and paired Hill activations")
    savefig(fig, out, f"abdominal_components_group{group}")
    fig, axes = plt.subplots(2, 2, figsize=(13, 7.5), sharex=True, sharey="row", layout="constrained")
    for source, title in enumerate(("Experimental activation", "MAP activation")):
        new, par = candidate[source], pair[source]
        axes[0, source].set_title(title)
        axes[0, source].plot(trial["t"], new[:, 0], label="Actual Hill angle")
        axes[0, source].axhline(par[30], color=".5", ls=":", label="Lower soft onset")
        axes[0, source].axhline(par[31], color=".6", ls="--", label="Upper soft onset")
        axes[0, source].set_ylabel("Angle (rad)")
        axes[0, source].legend(fontsize=8, loc="upper left", bbox_to_anchor=(0, .88))
        label_loads(axes[0, source], trial)
        axes[1, source].plot(trial["t"], boundary_torque(par, new[:, 0]), label="Soft-limit torque")
        axes[1, source].plot(trial["t"], -par[3] * new[:, 1], label="Joint damping torque")
        axes[1, source].set_ylabel("Other torque (N·m)")
        axes[1, source].set_xlabel("Time (s)")
        axes[1, source].legend(fontsize=8)
        for row in (0, 1):
            shade(axes[row, source], trial)
    fig.suptitle(f"Group {group} | Elastic soft boundary is separate from muscle torque | No angle clipping")
    savefig(fig, out, f"soft_boundary_group{group}")


def table(out, name, records):
    with (out / name).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def depth_summary_plot(out, rows):
    """Compare both flexion valleys; lower theta means deeper trunk flexion."""
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharey=True, layout="constrained")
    for source, title in enumerate(("Experimental activation", "MAP activation"), start=1):
        new = [r for r in rows if r["source"] == source and r["model"] == "candidate"]
        old = [r for r in rows if r["source"] == source and r["model"] == "previous"]
        old_by_case = {(r["loadKg"], r["speedCode"]): r for r in old}
        labels = [f'{int(r["loadKg"])} kg\ncode {int(r["speedCode"])}' for r in new]
        x = np.arange(len(new))
        for valley in (1, 2):
            ax = axes[valley - 1, source - 1]
            field = f"predictedValley{valley}Rad"
            reference = [r[f"referenceValley{valley}Rad"] for r in new]
            previous = [old_by_case[(r["loadKg"], r["speedCode"])][field] for r in new]
            current = [r[field] for r in new]
            ax.bar(x - .24, reference, .24, color="#d02929", alpha=.8, label="Original reference")
            ax.bar(x, previous, .24, color=".65", label="Previous calibration")
            ax.bar(x + .24, current, .24, color="#0877b7", label="Deeper-flexion calibration")
            ax.set_xticks(x, labels, fontsize=8)
            ax.set_ylabel(f"Valley {valley} angle (rad)")
            ax.set_title(title if valley == 1 else "Second flexion valley")
            ax.grid(axis="y", alpha=.2)
            ax.set_axisbelow(True)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=3, fontsize=8)
    fig.suptitle("All nine cases | Smaller angle = deeper flexion | Same fixed parameters within each source")
    savefig(fig, out, "valley_depths_all_conditions")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=ROOT / "data/processed/hill_compliant/neural_models.mat")
    parser.add_argument("--output", default="hill_neural_selected_evaluation")
    parser.add_argument("--previous-model", type=Path, default=ROOT / "data/processed/hill_compliant/neural_models_shape_baseline.mat")
    parser.add_argument("--angle-profile", choices=("shape", "depth"), default="depth")
    args = parser.parse_args()
    package = loadmat(args.model, simplify_cells=True)
    pair = np.asarray(package["parameters"]) if "parameters" in package else np.vstack([m["parameters"] for m in package["models"]])
    generators = [m["generator"] for m in package["models"]]
    assert pair.shape == (2, 32) and np.all(np.isfinite(pair))
    assert np.array_equal(pair[0, fit.SHARED], pair[1, fit.SHARED])
    records = fit.prepared([1, 2, 3], generators)
    fit.ANGLE_PROFILE = args.angle_profile
    angle_policy = fit.ANGLE_PROFILES[args.angle_profile]
    original_rollout = fit.rollout
    previous_package = loadmat(args.previous_model, simplify_cells=True) if args.previous_model.exists() else fit.original.PACKAGE
    old = np.vstack([m["parameters"] for m in previous_package["models"]])
    out = ROOT / "results" / args.output
    out.mkdir(parents=True, exist_ok=True)
    if "controller" in previous_package:
        fit.CONTROLLER[1:] = np.asarray(previous_package["controller"]).reshape(6)
    previous_neural = all("neuralParameters" in m for m in previous_package["models"])
    if previous_neural:
        import fit_hill_neural as neural
        neural.NEURAL[:] = neural.parameter_rows(np.vstack([m["neuralParameters"] for m in previous_package["models"]]))
        fit.rollout = neural.rollout
    rows = fit.metrics(old, records, label="previous")
    previous_outputs = {(record["group"], source): fit.rollout(old[source], record, source, steps=10)[0]
                        for record in records for source in (0, 1)}
    has_neural = "neuralParameters" in package or all("neuralParameters" in m for m in package["models"])
    if "controller" in package:
        fit.CONTROLLER[1:] = np.asarray(package["controller"]).reshape(6)
    if has_neural:
        import fit_hill_neural as neural
        neural.NEURAL[:] = neural.parameter_rows(package["neuralParameters"] if "neuralParameters" in package else np.vstack([m["neuralParameters"] for m in package["models"]]))
        fit.rollout = neural.rollout
    else:
        fit.rollout = original_rollout
    rows += fit.metrics(pair, records)
    fit.write_metrics(out, rows)
    depth_summary_plot(out, rows)
    numerical, waiting, bundles, phases = [], [], [], []
    for record in records:
        t, group = record["trial"], record["group"]
        previous, candidate = [], []
        for source in (0, 1):
            o, failed, a, inferred = fit.rollout(pair[source], record, source, steps=10)
            coarse, coarse_fail, _, _ = fit.rollout(pair[source], record, source, steps=20)
            assert not failed and not coarse_fail
            previous.append(previous_outputs[(group, source)])
            candidate.append(o)
            residual = np.max(np.abs(o[:, 4] - o[:, 8] - o[:, 9] - o[:, 10]))
            assert residual < 1e-6
            numerical.append({"group": group, "source": source + 1, "frontComponentResidualNm": float(residual),
                              "fiberEquilibriumResidualMax": float(np.max(o[:, 15])),
                              "momentumResidualMaxNm": float(np.max(np.abs(o[:, 16]))),
                              "halfStepAngleDifferenceRad": float(np.max(np.abs(o[:, 0] - coarse[:, 0]))),
                              "halfStepTorqueDifferenceNm": float(np.max(np.abs(o[:, 2] - coarse[:, 2])))})
            settled = np.array(t["waiting"], dtype=bool)
            for s in t["segments"]:
                stop = int(s["endIndex"])
                settled[stop: min(len(settled), stop + 200)] = False
            waiting.append({"group": group, "source": source + 1,
                            "maxSettledAngleErrorRad": float(np.max(np.abs(o[settled, 0] - t["standingAngle"]))),
                            "minimumSequenceAngleRad": float(np.min(o[:, 0])), "maximumSequenceAngleRad": float(np.max(o[:, 0])),
                            "softLimitPeakNm": float(np.max(np.abs(boundary_torque(pair[source], o[:, 0])))),
                            "frontFullSequencePeakNm": float(np.max(o[:, 4])),
                            "frontWaitingPeakNm": float(np.max(o[np.asarray(t["waiting"], dtype=bool), 4])),
                            "frontFullActivationPeak": float(np.max(o[:, 14]))})
            for segment in record["segments"]:
                ids = segment["indices"]
                for phase, mask in (
                    ("unloaded_flexion", (t["m"][ids] <= .5) & (o[ids, 1] < -.05)),
                    ("loaded_flexion", (t["m"][ids] > .5) & (o[ids, 1] < -.05)),
                    ("loaded_extension", (t["m"][ids] > .5) & (o[ids, 1] > .05)),
                    ("unloaded_extension", (t["m"][ids] <= .5) & (o[ids, 1] > .05)),
                ):
                    selected = ids[mask]
                    if len(selected):
                        phases.append({"group": group, "source": source + 1, "loadKg": segment["loadKg"],
                                       "speedCode": segment["speed"], "phase": phase, "samples": len(selected),
                                       "backMeanTorqueNm": float(np.mean(o[selected, 3])),
                                       "frontMeanTorqueNm": float(np.mean(o[selected, 4])),
                                       "frontToBackTorqueRatio": float(np.sum(o[selected, 4]) / max(np.sum(o[selected, 3]), 1e-12)),
                                       "backHillActivationMean": float(np.mean(o[selected, 13])),
                                       "frontHillActivationMean": float(np.mean(o[selected, 14]))})
            np.savez_compressed(out / f"group{group}_source{source + 1}.npz", output=o, inputs=a, inferred=inferred,
                                parameters=pair[source], previous=previous[-1], qReference=t["q"], torqueReference=t["torque"])
        plots(out, record, previous, candidate, pair,
              previous_label="Previous shallow calibration" if previous_neural else "Previous calibration")
        bundles.append({"group": group, "previous": previous, "candidate": candidate})
    table(out, "numerical_checks.csv", numerical)
    table(out, "waiting_checks.csv", waiting)
    table(out, "phase_recruitment.csv", phases)
    table(out, "parameters.csv", [{"parameter": name, "experimental": pair[0, i], "map": pair[1, i],
                                   "lower": fit.LOW[i], "upper": fit.HIGH[i]} for i, name in enumerate(PARAMETER_NAMES)])
    chosen = [r for r in rows if r["model"] == "candidate"]
    checks = {
        "sameParametersAcrossNineConditionsPerSource": True,
        "sameModelStructureAcrossSources": True,
        "sharedBodyLoadAndSoftLimitParameters": True,
        "allMovementAnglesWithinBroadRange": all(r["minimumAngleRad"] >= angle_policy["lower"] and r["maximumAngleRad"] <= angle_policy["upper"] for r in chosen),
        "entireSequenceAnglesWithinAllowedRange": all(r["minimumSequenceAngleRad"] >= angle_policy["sequenceLower"] and r["maximumSequenceAngleRad"] <= angle_policy["sequenceUpper"] for r in waiting),
        "allValleysWithinReferenceDepthTolerance": args.angle_profile != "depth" or all(r["maximumValleyTooShallowRad"] <= angle_policy["valleyTolerance"] for r in chosen),
        "allExcursionsAtLeast075Reference": args.angle_profile != "depth" or all(r["excursionRatio"] >= angle_policy["excursionRatio"] for r in chosen),
        "allMiddleReboundsAtLeast065Reference": args.angle_profile != "depth" or all(r["middleReboundRatio"] >= .65 for r in chosen),
        "allConditionsHaveTwoMajorValleys": all(r["predictedValleys"] == 2 for r in chosen),
        "minimumShapeCorrelationAbove06": min(r["shapeCorrelation"] for r in chosen) >= .6,
        "allMeanAbdominalToBackTorqueRatiosBelow03": max(r["frontBackMeanTorqueRatio"] for r in chosen) <= .3,
        "allAbdominalPeakTorquesBelow100Nm": max(r["frontTorquePeakNm"] for r in chosen) <= 100,
        "fullSequenceAbdominalPeakBelow50Nm": max(r["frontFullSequencePeakNm"] for r in waiting) <= 50,
        "softBoundaryPeakBelow50Nm": max(r["softLimitPeakNm"] for r in waiting) <= 50,
        "settledWaitingErrorBelow001Rad": max(r["maxSettledAngleErrorRad"] for r in waiting) <= .01,
        "halfStepAngleDifferenceBelow002Rad": max(r["halfStepAngleDifferenceRad"] for r in numerical) <= .02,
        "halfStepTorqueDifferenceBelow3Nm": max(r["halfStepTorqueDifferenceNm"] for r in numerical) <= 3,
    }
    status = {"checks": checks, "screenPassed": all(checks.values()),
              "criteriaStatus": "Demo engineering screening, not universal physiological ratios or individual validation",
              "modelSHA256": hashlib.sha256(args.model.read_bytes()).hexdigest(),
              "angleProfile": args.angle_profile, "anglePolicy": angle_policy,
              "previousModel": (args.previous_model.relative_to(ROOT).as_posix() if args.previous_model.is_relative_to(ROOT)
                                else args.previous_model.as_posix()) if args.previous_model.exists() else "Original Hill calibration",
              "referenceStatus": "Exact original project Scale",
              "activationStatus": "Not verified %MVC"}
    status["modelIncludesNeuralControl"] = has_neural
    status["motionIntentIsNotOutputAngle"] = True
    if has_neural:
        table(out, "neural_parameters.csv", [{"parameter": name, "experimental": neural.NEURAL[0, i],
                                               "map": neural.NEURAL[1, i], "lower": neural.NLOW[i], "upper": neural.NHIGH[i]}
                                              for i, name in enumerate(neural.NEURAL_NAMES)])
    (out / "screening.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
    for source in (1, 2):
        selected = [r for r in chosen if r["source"] == source]
        print("SOURCE", source, "MEAN CORR", np.mean([r["shapeCorrelation"] for r in selected]),
              "W", sum(r["predictedValleys"] == 2 for r in selected),
              "SHALLOW ERROR MAX", max(r["maximumValleyTooShallowRad"] for r in selected),
              "ABD MEAN", np.mean([r["frontTorqueMeanNm"] for r in selected]),
              "ABD PEAK", max(r["frontTorquePeakNm"] for r in selected))
    print(json.dumps(status, indent=2))
    print("Evaluation:", out)


if __name__ == "__main__":
    main()
