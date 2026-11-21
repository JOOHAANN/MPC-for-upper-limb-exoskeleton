"""Diagnose post-movement transients with matched-state feedback ablations.

The selected controller and published trajectories are never modified.
Counterfactuals start at each original movement end and replay the same
subsequent inputs, mass and recorded actuator commands for at most one second.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from matplotlib.backends.backend_pdf import PdfPages
from run_hill_mpc import Plant, ROOT, LABELS, prediction_dict
from plot_hill_mpc_results import plt, sha256, RED, BLUE, AMBER, GREEN

OUT = ROOT / "results/hill_mpc"
TARGET = ROOT / "docs/hill_mpc_evaluation"


def modified_plant(model, controller, angle, variant):
    modified = dict(model)
    modified["neuralParameters"] = np.array(model["neuralParameters"], copy=True)
    cfg = np.array(controller, copy=True)
    if variant in ("noStandingP", "noStandingPD", "noBothFeedbacks"):
        cfg[0] = 0.
    if variant in ("noStandingD", "noStandingPD", "noBothFeedbacks"):
        cfg[1] = 0.
    if variant in ("noMovementFeedback", "noBothFeedbacks"):
        modified["neuralParameters"][3:5] = 0.
    return Plant(modified, cfg, angle)


def write_csv(name, records):
    with (OUT / name).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader(); writer.writerows(records)


def main(publish=False):
    saved = loadmat(OUT / "hill_mpc.mat", simplify_cells=True)["result"]
    cache_path = ROOT / "data/processed/hill_compliant/neural_models.mat"
    cache = loadmat(cache_path, simplify_cells=True)
    controller = np.asarray(cache["controller"])
    rows, windows, diagnostics, replay_checks = [], [], [], []
    variants = ("noStandingP", "noStandingD", "noStandingPD", "noMovementFeedback", "noBothFeedbacks")
    for group in saved["groups"]:
        seq, dt = group["sequence"], group["sequence"]["dt"]
        ends = [int(segment["endIndex"])-1 for segment in seq["segments"]]
        for source, label in enumerate(LABELS):
            model = cache["models"][source]
            neural = model["neuralParameters"]
            inputs = group["inputs"][source]
            actual, nominal = group["predictions"][source], group["nominal"][source]
            plant = Plant(model, controller, seq["standingAngle"])
            state, first = plant.init(inputs[0], seq["m"][0], seq["waiting"][0], np.array([seq["q"][0], seq["v0"]]))
            replay = np.zeros((len(inputs), 27)); replay[0] = first
            snapshots = {}
            detail = np.zeros((len(inputs), 5))
            for k in range(1, len(inputs)):
                # Capture the exact q,v used by the last internal feedback step.
                fraction = .9
                mid_input = inputs[k-1]+fraction*(inputs[k]-inputs[k-1])
                mid_mass = seq["m"][k-1]+fraction*(seq["m"][k]-seq["m"][k-1])
                mid_wait = float(seq["waiting"][k-1])+fraction*(float(seq["waiting"][k])-float(seq["waiting"][k-1]))
                previous_u, final_u = actual["assistanceTorque"][k-1:k+1]
                mid_u = previous_u+fraction*(final_u-previous_u)
                plant.step(state, mid_input, mid_mass, mid_wait, mid_u, dt*.9, 9)
                q_for_control, v_for_control = state[:2].copy()
                replay[k] = plant.step(state, inputs[k], seq["m"][k], seq["waiting"][k], final_u, dt*.1, 1)
                p_term = controller[0]*(seq["standingAngle"]-q_for_control)
                d_term = -controller[1]*v_for_control
                stance = controller[3]*np.tanh((p_term+d_term+controller[2]*state[7])/controller[3])
                moving = neural[3]*(state[10]-q_for_control)-neural[4]*v_for_control
                mixed = state[6]*stance+(1-state[6])*moving
                detail[k] = [p_term, d_term, stance, moving, mixed]
                if k in ends:
                    snapshots[k] = state.copy()
            assert np.max(np.abs(replay[:, 0]-actual["q"])) < 1e-7
            assert np.max(np.abs(replay[:, 3]-actual["backTorque"])) < 1e-4
            assert np.max(np.abs(detail[1:, 4]-actual["postureCorrectionRequest"][1:])) < 1e-5
            replay_checks.append({"group": int(group["group"]), "source": label,
                                  "maximumAngleDifferenceRad": float(np.max(np.abs(replay[:, 0]-actual["q"]))),
                                  "maximumBackTorqueDifferenceNm": float(np.max(np.abs(replay[:, 3]-actual["backTorque"]))),
                                  "maximumFeedbackReconstructionDifferenceNm": float(np.max(np.abs(detail[1:, 4]-actual["postureCorrectionRequest"][1:])))})
            for segment, end in zip(seq["segments"], ends):
                stop = min(len(inputs), end+1+round(1./dt))
                post = slice(end+1, stop)
                selected = slice(max(1, end-round(.45/dt)), stop)
                index = post.start+int(np.argmax(actual["backTorque"][post]))
                counterfactuals = {}
                for variant in variants:
                    changed = modified_plant(model, controller, seq["standingAngle"], variant)
                    copied = snapshots[end].copy()
                    outputs = []
                    for k in range(end+1, stop):
                        outputs.append(changed.step(copied, inputs[k], seq["m"][k], seq["waiting"][k], actual["assistanceTorque"][k], dt, 10))
                    counterfactuals[variant] = prediction_dict(np.asarray(outputs), model["parameters"])
                row = {"group": int(group["group"]), "source": label, "loadKg": int(segment["loadKg"]),
                       "speedCode": int(segment["speedCode"]), "movementEndSeconds": float(seq["t"][end]),
                       "peakTimeSeconds": float(seq["t"][index]), "peakDelaySeconds": float(seq["t"][index]-seq["t"][end]),
                       "postBackPeakNm": float(actual["backTorque"][index]),
                       "postNetTorqueAtBackPeakNm": float(actual["torque"][index]),
                       "unassistedPostBackPeakNm": float(nominal["backTorque"][post].max()),
                       "angleAtPeakRad": float(actual["q"][index]), "velocityAtPeakRadPerSec": float(actual["v"][index]),
                       "standingPAtPeakNm": float(detail[index, 0]), "standingDAtPeakNm": float(detail[index, 1]),
                       "standingSaturatedPDAtPeakNm": float(detail[index, 2]), "movementFeedbackAtPeakNm": float(detail[index, 3]),
                       "mixedFeedbackAtPeakNm": float(detail[index, 4]), "standingWeightAtPeak": float(actual["balanceWeight"][index]),
                       "gravityDemandAtPeakNm": float(-actual["gravityTorque"][index]),
                       "assistanceAtPeakNm": float(actual["assistanceTorque"][index]),
                       "maximumPostAssistanceNm": float(actual["assistanceTorque"][post].max()),
                       "boundaryAtPeakNm": float(actual["boundaryTorque"][index]),
                       "backActiveAtPeakNm": float(actual["backActive"][index]), "backPassiveAtPeakNm": float(actual["backPassive"][index]),
                       "backFiberDampingAtPeakNm": float(actual["backFiberDamping"][index]),
                       "motionIntentAtPeakRad": float(actual["motionIntent"][index]),
                       "sourceBackActivationAtPeak": float(inputs[index, 0]),
                       "feedbackRecruitmentFractionAtPeak": float(neural[6]+(1-neural[6])*actual["balanceWeight"][index]),
                       "effectiveLoadAtPeakKg": float(seq["m"][index])}
                for variant in variants:
                    row[f"{variant}PostBackPeakNm"] = float(counterfactuals[variant]["backTorque"].max())
                    row[f"{variant}AtOriginalPeakNm"] = float(counterfactuals[variant]["backTorque"][index-post.start])
                rows.append(row)
                windows.append({"row": row, "seq": seq, "source": source, "actual": actual, "nominal": nominal,
                                "detail": detail, "inputs": inputs, "selected": selected, "post": post,
                                "counterfactuals": counterfactuals})
                for k in range(selected.start, selected.stop):
                    diagnostics.append({"group": int(group["group"]), "source": label, "loadKg": int(segment["loadKg"]),
                                        "timeSeconds": float(seq["t"][k]), "relativeToEndSeconds": float(seq["t"][k]-seq["t"][end]),
                                        "standingPRequestNm": float(detail[k, 0]), "standingDRequestNm": float(detail[k, 1]),
                                        "standingPDRequestNm": float(detail[k, 2]), "movementFeedbackRequestNm": float(detail[k, 3]),
                                        "mixedFeedbackRequestNm": float(detail[k, 4])})
    write_csv("post_movement_peak_audit.csv", rows)
    write_csv("post_movement_feedback_traces.csv", diagnostics)
    plt.rcParams.update({"font.size": 10, "pdf.fonttype": 42})
    files = ["post_movement_peak_audit.csv", "post_movement_feedback_traces.csv"]
    bundle = "post_movement_spike_diagnosis.pdf"
    with PdfPages(OUT / bundle) as pdf:
        for group_id in (1, 2, 3):
            fig, axes = plt.subplots(3, 2, figsize=(13.8, 9.6), sharex=True)
            for source, label in enumerate(LABELS):
                chosen = [w for w in windows if w["row"]["group"] == group_id and w["source"] == source]
                for ax, window in zip(axes[:, source], chosen):
                    row, selected = window["row"], window["selected"]
                    t = window["seq"]["t"][selected]-row["movementEndSeconds"]
                    post_t = window["seq"]["t"][window["post"]]-row["movementEndSeconds"]
                    ax.plot(t, window["actual"]["backTorque"][selected], color=BLUE, lw=1.6, label="Current back torque")
                    ax.plot(t, window["nominal"]["backTorque"][selected], color=AMBER, ls="--", lw=1.4, label="Unassisted model (different history)")
                    ax.plot(post_t, window["counterfactuals"]["noStandingPD"]["backTorque"], color=GREEN, ls=":", lw=1.7,
                            label="Same end state; standing PD disabled")
                    ax.plot(post_t, window["counterfactuals"]["noBothFeedbacks"]["backTorque"], color="#8054a6", ls="-.", lw=1.2,
                            label="Same end state; both feedbacks disabled")
                    ax.axvline(0, color="#444444", ls="--", lw=.8); ax.axvspan(0, t[-1], color="#eeeeee", alpha=.35)
                    ax.set_xlim(t[0], t[-1]); ax.grid(alpha=.2); ax.set_ylabel("Back torque (N·m)")
                    ax.set_title(f"{row['loadKg']} kg / code {row['speedCode']} — {label}", fontsize=11)
                    if ax in axes[-1, :]: ax.set_xlabel("Time relative to movement end (s)")
            handles, labels = axes[0, 0].get_legend_handles_labels()
            fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(.5, .008), fontsize=9)
            fig.suptitle(f"Post-movement feedback ablation — sequence {group_id}", fontsize=15, y=.995)
            fig.tight_layout(rect=(0, .085, 1, .965), h_pad=1.8)
            for suffix in ("png", "pdf"):
                name = f"post_movement_ablation_group{group_id}.{suffix}"
                fig.savefig(OUT / name, dpi=180, bbox_inches="tight", facecolor="white"); files.append(name)
            pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)
        # The largest tail peak is a useful example of both feedback layers.
        example = max(windows, key=lambda w: w["row"]["postBackPeakNm"])
        row, selected = example["row"], example["selected"]
        seq, p, detail = example["seq"], example["actual"], example["detail"]
        t = seq["t"][selected]-row["movementEndSeconds"]
        fig, axes = plt.subplots(4, 1, figsize=(11.7, 10.5), sharex=True)
        for key, label, color in (("backTorque", "Back total", BLUE), ("backActive", "Back active", GREEN),
                                  ("frontTorque", "Abdominal total", AMBER), ("assistanceTorque", "MPC actuator", "#8054a6"),
                                  ("boundaryTorque", "Soft boundary", "#555555")):
            axes[0].plot(t, p[key][selected], label=label, color=color, ls="--" if key == "backActive" else "-")
        axes[0].set_ylabel("Actual torques (N·m)"); axes[0].legend(frameon=False, ncol=3, fontsize=9)
        axes[1].plot(t, detail[selected, 0], label="Standing P request", color=BLUE)
        axes[1].plot(t, detail[selected, 1], label="Standing D request", color=AMBER)
        axes[1].plot(t, detail[selected, 3], label="Movement feedback request", color="#8054a6", ls="--")
        axes[1].plot(t, p["postureCorrectionRequest"][selected], label="Recorded blended feedback request", color="#222222", lw=1.6)
        axes[1].set_ylabel("Requests (N·m)"); axes[1].legend(frameon=False, ncol=2, fontsize=9)
        axes[2].plot(t, p["q"][selected], label="Model angle", color=BLUE)
        axes[2].axhline(seq["standingAngle"], color="#444444", ls="--", label="Standing target")
        post_t = seq["t"][example["post"]]-row["movementEndSeconds"]
        axes[2].plot(post_t, example["counterfactuals"]["noStandingPD"]["q"], color=GREEN, ls=":", label="Standing PD disabled")
        axes[2].plot(post_t, example["counterfactuals"]["noBothFeedbacks"]["q"], color="#8054a6", ls="-.", label="Both feedbacks disabled")
        axes[2].set_ylabel("Angle (rad)")
        twin = axes[2].twinx(); twin.plot(t, p["v"][selected], color=AMBER, ls=":", label="Angular velocity")
        twin.set_ylabel("Angular velocity (rad/s)", color=AMBER)
        axes[2].legend(frameon=True, facecolor="white", edgecolor="none", framealpha=.95,
                       fontsize=8, ncol=2, loc="upper left")
        axes[3].plot(t, p["balanceWeight"][selected], label="Standing blend weight", color="#222222")
        axes[3].plot(t, p["backCommand"][selected], label="Back excitation command", color=GREEN)
        axes[3].plot(t, p["backActivation"][selected], label="Back activation state", color=BLUE)
        axes[3].plot(t, example["inputs"][selected, 0], label="Source back activation", color=AMBER, ls="--")
        axes[3].set_ylabel("Weight / activation"); axes[3].set_xlabel("Time relative to movement end (s)")
        axes[3].legend(frameon=False, ncol=2, fontsize=9)
        for ax in axes:
            ax.axvline(0, color="#444444", ls="--", lw=.8); ax.grid(alpha=.2); ax.set_xlim(t[0], t[-1])
        fig.suptitle(f"Largest post-movement peak: sequence {row['group']}, {row['loadKg']} kg, {row['source']}", fontsize=13, y=.994)
        fig.text(.5, .009, "Feedback requests act through muscle excitation; they are not additional joint torques. The peak is active Hill muscle torque.", ha="center", fontsize=8.5)
        fig.tight_layout(rect=(0, .027, 1, .965), h_pad=1.6)
        for suffix in ("png", "pdf"):
            name = f"post_movement_component_detail.{suffix}"
            fig.savefig(OUT / name, dpi=180, bbox_inches="tight", facecolor="white"); files.append(name)
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)
    files.append(bundle)
    report = {"sourceResultSHA256": sha256(OUT / "hill_mpc.mat"), "modelSHA256": sha256(cache_path),
              "standingController": {"kp": controller[0], "kd": controller[1], "ki": controller[2],
                                     "requestCapNm": controller[3], "coContractionNm": controller[4], "blendTimeSeconds": controller[5]},
              "counterfactualPolicy": "Same state at original movement end, same future inputs/load/recorded actuator command for up to 1 s; only named feedback gains changed. Diagnostic only, not a replacement controller.",
              "replayChecks": replay_checks,
              "movementFeedback": [{"source": LABELS[k], "kp": float(model["neuralParameters"][3]),
                                     "kd": float(model["neuralParameters"][4]), "feedbackRecruitmentFraction": float(model["neuralParameters"][6])}
                                    for k, model in enumerate(cache["models"])],
              "largestPostPeak": row, "transitionCount": len(rows)}
    name = "post_movement_diagnosis.json"
    (OUT / name).write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8"); files.append(name)
    if publish:
        import shutil
        for name in files: shutil.copy2(OUT / name, TARGET / name)
        manifest_path = TARGET / "export_metadata.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["files"] = [{"file": p.name, "sha256": sha256(p)} for p in sorted(TARGET.iterdir()) if p.is_file() and p.name != "export_metadata.json"]
        manifest_path.write_text(json.dumps(manifest, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--publish", action="store_true")
    main(publish=parser.parse_args().publish)
