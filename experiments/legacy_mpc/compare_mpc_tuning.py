"""Compare saved MPC settings on identical Hill plants, inputs and references."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from scipy.signal import find_peaks
from matplotlib.backends.backend_pdf import PdfPages
from plot_hill_mpc_results import (ROOT, OUT, LABELS, NOMINAL, ASSISTED, REFERENCE,
                                   plt, segment_slice, shade_movements,
                                   back_threshold_lines, BACK_THRESHOLDS_NM, publish, sha256)

BEFORE = ROOT / "results/hill_mpc_before_predictive_tuning"


def residual_peak_figure(old_groups, groups, options):
    """Expose the remaining peak and actuator saturation without clipping data."""
    worst = None
    for group in groups:
        for source in range(2):
            p = group["predictions"][source]
            for segment in group["sequence"]["segments"]:
                selected = segment_slice(segment)
                threshold = BACK_THRESHOLDS_NM[int(segment["loadKg"])]
                index = selected.start+int(np.argmax(p["backTorque"][selected]))
                excess = float(p["backTorque"][index]-threshold)
                if worst is None or excess > worst[0]:
                    worst = (excess, group, source, segment, index, threshold)
    if worst[0] <= 0:
        return None, {"peakExcessNm": 0.}
    excess, group, source, segment, index, threshold = worst
    seq = group["sequence"]
    p0 = old_groups[int(group["group"])]["predictions"][source]
    p1 = group["predictions"][source]
    half_width = round(.35/seq["dt"])
    chosen = slice(max(0, index-half_width), min(len(seq["t"]), index+half_width+1))
    t = seq["t"][chosen]
    fig, axes = plt.subplots(2, 1, figsize=(11.5, 6.5), sharex=True)
    for ax, key, ylabel in ((axes[0], "backTorque", "Back muscle torque (N·m)"),
                             (axes[1], "assistanceTorque", "Applied assistance (N·m)")):
        ax.plot(t, p0[key][chosen], **NOMINAL, label="Previous MPC")
        ax.plot(t, p1[key][chosen], **ASSISTED, label="MPC with headroom + velocity term")
        ax.axvline(seq["t"][index], color="#777777", ls=":", lw=.9)
        ax.set_ylabel(ylabel); ax.grid(alpha=.2)
    axes[0].axhline(threshold, color="#333333", ls="--", lw=1.2, label=f"Back threshold {threshold:g} N·m")
    axes[0].fill_between(t, threshold, p1["backTorque"][chosen], where=p1["backTorque"][chosen]>threshold,
                         color="#cc3434", alpha=.25, interpolate=True, label="Remaining exceedance")
    axes[0].set_title(f"Remaining peak: +{excess:.2f} N·m at t = {seq['t'][index]:.2f} s", fontsize=11)
    axes[1].axhline(options["capNm"], color="#333333", ls="--", lw=1.2, label=f"Actuator cap {options['capNm']:g} N·m")
    axes[1].annotate(f"Actuator cap {options['capNm']:g} N·m", (.02, options["capNm"]),
                     xycoords=("axes fraction", "data"), xytext=(0, 4), textcoords="offset points",
                     fontsize=9, va="bottom", bbox={"facecolor": "white", "edgecolor": "none", "alpha": .9, "pad": 1})
    axes[1].set_xlabel("Time (s)")
    for ax in axes:
        ax.set_xlim(t[0], t[-1])
    axes[0].legend(frameon=False, fontsize=9, loc="upper right")
    fig.suptitle(f"Residual peak detail — {int(segment['loadKg'])} kg, speed code {int(segment['speedCode'])}, {LABELS[source]}", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, .955), h_pad=1.8)
    record = {"group": int(group["group"]), "source": LABELS[source], "loadKg": int(segment["loadKg"]),
              "speedCode": int(segment["speedCode"]), "timeSeconds": float(seq["t"][index]),
              "peakExcessNm": excess, "backTorqueNm": float(p1["backTorque"][index]),
              "appliedAssistanceNm": float(p1["assistanceTorque"][index]), "actuatorCapNm": options["capNm"]}
    return fig, record


def compare(publish_docs=False, before_directory=BEFORE, prefix="mpc_headroom"):
    before_directory = Path(before_directory)
    before = loadmat(before_directory / "hill_mpc.mat", simplify_cells=True)["result"]
    after = loadmat(OUT / "hill_mpc.mat", simplify_cells=True)["result"]
    old_meta = json.loads((before_directory / "evaluation.json").read_text())
    metadata = json.loads((OUT / "evaluation.json").read_text())
    assert old_meta["modelSHA256"] == metadata["modelSHA256"]
    for source in range(2):
        for key in ("parameters", "neuralParameters"):
            assert np.array_equal(before["models"][source][key], after["models"][source][key])
    old_groups = {int(g["group"]): g for g in before["groups"]}
    rows, samples, figure_names = [], [[], []], []
    old_weight = old_meta["options"].get("backTorqueWeight", 2.)
    new_weight = metadata["options"]["backTorqueWeight"]
    plt.rcParams.update({"font.size": 10, "pdf.fonttype": 42})
    bundle = f"{prefix}_before_after.pdf"
    with PdfPages(OUT / bundle) as pdf:
        for group in after["groups"]:
            old = old_groups[int(group["group"])]
            seq, t = group["sequence"], group["sequence"]["t"]
            for key in ("t", "q", "torque", "waiting", "loadKg"):
                assert np.array_equal(old["sequence"][key], seq[key])
            fig, axes = plt.subplots(3, 2, figsize=(14.2, 9.6), sharex="col")
            for source, label in enumerate(LABELS):
                assert np.array_equal(old["inputs"][source], group["inputs"][source])
                p0, p1 = old["predictions"][source], group["predictions"][source]
                for row, key, ylabel in ((0, "backTorque", "Back muscle torque (N·m)"),
                                          (1, "q", "Trunk angle (rad)"),
                                          (2, "assistanceTorque", "Assistance (N·m)")):
                    ax = axes[row, source]
                    ax.plot(t, p0[key], **NOMINAL, label=f"Previous MPC (back weight {old_weight:g})")
                    ax.plot(t, p1[key], **ASSISTED, label=f"Tuned MPC (back weight {new_weight:g})")
                    ax.set_ylabel(ylabel)
                axes[0, source].set_title(label, pad=12)
                axes[0, source].margins(y=.22)
                back_threshold_lines(axes[0, source], seq)
                waiting = np.asarray(seq["waiting"], dtype=bool)
                axes[1, source].plot(t, np.where(waiting, np.nan, seq["q"]), **REFERENCE,
                                     label="Experimental angle reference")
                axes[2, source].axhline(metadata["options"]["capNm"], color="#333333", ls="--", lw=1.)
                axes[2, source].set_ylim(-3, 89)
                axes[2, source].set_xlabel("Time (s)")
                shade_movements(axes[:, source], seq)
                for segment in seq["segments"]:
                    selected = segment_slice(segment)
                    threshold = BACK_THRESHOLDS_NM[int(segment["loadKg"])]
                    e0 = np.maximum(p0["backTorque"][selected]-threshold, 0.)
                    e1 = np.maximum(p1["backTorque"][selected]-threshold, 0.)
                    q0, q1, ref = p0["q"][selected], p1["q"][selected], seq["q"][selected]
                    valleys = lambda q: len(find_peaks(-q, prominence=.12, distance=50)[0])
                    entry = {
                        "group": int(group["group"]), "source": label,
                        "loadKg": int(segment["loadKg"]), "speedCode": int(segment["speedCode"]),
                        "thresholdNm": threshold,
                        "previousBackPeakNm": float(p0["backTorque"][selected].max()),
                        "tunedBackPeakNm": float(p1["backTorque"][selected].max()),
                        "previousPeakExcessNm": float(e0.max()), "tunedPeakExcessNm": float(e1.max()),
                        "previousMeanPositiveExcessNm": float(e0.mean()), "tunedMeanPositiveExcessNm": float(e1.mean()),
                        "previousExcessIntegralNmS": float(e0.sum()*seq["dt"]),
                        "tunedExcessIntegralNmS": float(e1.sum()*seq["dt"]),
                        "previousSecondsAboveThresholdPlus0_5Nm": float(np.sum(e0>.5)*seq["dt"]),
                        "tunedSecondsAboveThresholdPlus0_5Nm": float(np.sum(e1>.5)*seq["dt"]),
                        "previousThetaRMSEDegrees": float(np.sqrt(np.mean((q0-ref)**2))*180/np.pi),
                        "tunedThetaRMSEDegrees": float(np.sqrt(np.mean((q1-ref)**2))*180/np.pi),
                        "maximumAngleChangeRad": float(np.max(np.abs(q1-q0))),
                        "referenceValleys": valleys(ref), "previousValleys": valleys(q0), "tunedValleys": valleys(q1),
                        "tunedAssistancePeakNm": float(p1["assistanceTorque"][selected].max()),
                        "tunedAbdomenBackMeanRatio": float(p1["frontTorque"][selected].mean()/p1["backTorque"][selected].mean()),
                        "previousSecondsAboveThresholdPlus1Nm": float(np.sum(e0>1.)*seq["dt"]),
                        "tunedSecondsAboveThresholdPlus1Nm": float(np.sum(e1>1.)*seq["dt"]),
                    }
                    assert entry["tunedValleys"] == entry["referenceValleys"] == 2
                    rows.append(entry)
                    samples[source].append((e0, e1))
            handles, labels = axes[0, 0].get_legend_handles_labels()
            handles.append(axes[1, 0].lines[-1]); labels.append("Experimental angle reference")
            fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(.5, .027))
            fig.suptitle(f"MPC tuning before / after — sequence {int(group['group'])}", fontsize=15, y=.994)
            opts = metadata["options"]
            fig.text(.5, .009,
                     f"Back weight {new_weight:g}; headroom {opts.get('backTorqueMarginNm',0):g} N·m; velocity weight {opts.get('velocityWeight',0):g}; "
                     f"update {1000*opts['controlPeriod']:g} ms; actuator 80 N·m / 200 N·m/s. Same Hill parameters and inputs.",
                     fontsize=8.5, ha="center")
            fig.tight_layout(rect=(0, .108, 1, .96), h_pad=1.8)
            name = f"{prefix}_comparison_group{int(group['group'])}"
            for suffix in ("png", "pdf"):
                filename = f"{name}.{suffix}"
                fig.savefig(OUT / filename, dpi=180, bbox_inches="tight", facecolor="white")
                figure_names.append(filename)
            pdf.savefig(fig, bbox_inches="tight", facecolor="white")
            plt.close(fig)
        fig, worst_peak = residual_peak_figure(old_groups, after["groups"], metadata["options"])
        if fig is not None:
            for suffix in ("png", "pdf"):
                filename = f"{prefix}_residual_peak.{suffix}"
                fig.savefig(OUT / filename, dpi=180, bbox_inches="tight", facecolor="white")
                figure_names.append(filename)
            pdf.savefig(fig, bbox_inches="tight", facecolor="white")
            plt.close(fig)
    figure_names.append(bundle)
    assert len(rows) == 18
    csv_name = f"{prefix}_comparison.csv"
    with (OUT / csv_name).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    figure_names.append(csv_name)
    summary = []
    for source, label in enumerate(LABELS):
        selected = [row for row in rows if row["source"] == label]
        e0 = np.concatenate([e[0] for e in samples[source]])
        e1 = np.concatenate([e[1] for e in samples[source]])
        summary.append({"source": label, "previousMeanPositiveExcessNm": float(e0.mean()),
                        "tunedMeanPositiveExcessNm": float(e1.mean()),
                        "meanPositiveExcessReductionPercent": float(100*(1-e1.mean()/e0.mean())),
                        "previousWorstPeakExcessNm": float(e0.max()), "tunedWorstPeakExcessNm": float(e1.max()),
                        "previousMeanCaseThetaRMSEDegrees": float(np.mean([row["previousThetaRMSEDegrees"] for row in selected])),
                        "tunedMeanCaseThetaRMSEDegrees": float(np.mean([row["tunedThetaRMSEDegrees"] for row in selected])),
                        "maximumAngleChangeRad": max(row["maximumAngleChangeRad"] for row in selected),
                        "allNineWShapes": all(row["tunedValleys"] == 2 for row in selected),
                        "allNineMeanExcessReduced": all(row["tunedMeanPositiveExcessNm"] < row["previousMeanPositiveExcessNm"] for row in selected),
                        "allNinePeakExcessReduced": all(row["tunedPeakExcessNm"] < row["previousPeakExcessNm"] for row in selected)})
        summary[-1].update({"previousSecondsAboveThresholdPlus1Nm": float(sum(row["previousSecondsAboveThresholdPlus1Nm"] for row in selected)),
                            "tunedSecondsAboveThresholdPlus1Nm": float(sum(row["tunedSecondsAboveThresholdPlus1Nm"] for row in selected)),
                            "casesWithoutExceedance": sum(row["tunedPeakExcessNm"] == 0 for row in selected),
                            "casesWithin1Nm": sum(row["tunedPeakExcessNm"] <= 1. for row in selected)})
    report = {"previousResultSHA256": sha256(before_directory / "hill_mpc.mat"),
              "tunedResultSHA256": sha256(OUT / "hill_mpc.mat"), "modelSHA256": metadata["modelSHA256"],
              "previousBackWeight": old_weight, "tunedBackWeight": new_weight,
              "previousOptions": old_meta["options"], "tunedOptions": metadata["options"],
              "worstResidualPeak": worst_peak,
              "metricDefinition": "Mean positive excess is mean(max(backTorque - loadThreshold, 0)) over all movement samples, including zero-excess samples; waiting is excluded.",
              "summary": summary}
    json_name = f"{prefix}_summary.json"
    (OUT / json_name).write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    figure_names.append(json_name)
    if publish_docs:
        publish(figure_names, metadata)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--publish", action="store_true")
    parser.add_argument("--before", type=Path, default=BEFORE)
    parser.add_argument("--prefix", default="mpc_headroom")
    args = parser.parse_args()
    compare(publish_docs=args.publish, before_directory=args.before, prefix=args.prefix)
