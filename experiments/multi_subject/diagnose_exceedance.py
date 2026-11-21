"""Matched-state MPC forensics; future-input branches are OFFLINE diagnostics.

An assistance-rate/cap counterfactual is not a proof of global infeasibility.
Do not use diagnostic future samples in controller tuning or online prediction.
"""

import argparse
import json
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from scipy.io import loadmat

from audit_subjects import ROOT
from fit_subject import many

sys.path.insert(0, str(ROOT / "experiments/legacy_mpc"))
from run_hill_mpc import Plant


def diagnose(folder, output, quality):
    cache = loadmat(folder / "neural_models.mat", simplify_cells=True)
    result = loadmat(folder / "hill_mpc.mat", simplify_cells=True)["result"]
    worst = None
    for group in many(result["groups"])[:-1]:
        for source, pred in enumerate(many(group["predictions"])):
            for seg in many(group["sequence"]["segments"]):
                ids = np.arange(int(seg["startIndex"]) - 1, int(seg["endIndex"]))
                limit = float(np.interp(seg["loadKg"], [5, 10, 15], [100, 140, 180]))
                k = int(ids[np.argmax(pred["backTorque"][ids])])
                score = pred["backTorque"][k] - limit
                if worst is None or score > worst[0]:
                    worst = (score, group, source, seg, k, limit)
    score, group, source, seg, peak, limit = worst
    seq = group["sequence"]
    saved = many(group["predictions"])[source]
    inputs = group["inputs"][source]
    model = many(cache["models"])[source]
    options = result["options"]
    plant = Plant(model, cache["controller"], seq["standingAngle"])
    state, first = plant.init(
        inputs[0], seq["m"][0], seq["waiting"][0], np.array([seq["q"][0], seq["v0"]])
    )
    states = [state.copy()]
    parity = 0.0
    for k in range(1, len(inputs)):
        out = plant.step(
            state,
            inputs[k],
            seq["m"][k],
            seq["waiting"][k],
            saved["assistanceTorque"][k],
            seq["dt"],
            10,
        )
        parity = max(parity, abs(out[3] - saved["backTorque"][k]))
        states.append(state.copy())
    assert parity < 1e-6
    left = max(0, peak - 30)
    right = min(len(inputs) - 1, peak + 20)
    start = states[left]
    slope = (inputs[left] - inputs[max(0, left - 5)]) / max(
        0.01, (left - max(0, left - 5)) * 0.01
    )
    traces = {}
    for branch in (
        "recorded_inputs_and_controls",
        "held_EMG_and_load",
        "held_EMG_actual_load",
        "causal_EMG_trend",
        "maximum_rate_limited_assistance",
    ):
        state = start.copy()
        trace = []
        for k in range(left + 1, right + 1):
            a = inputs[k]
            mass = seq["m"][k]
            wait = seq["waiting"][k]
            u = saved["assistanceTorque"][k]
            if branch.startswith("held_EMG"):
                a = inputs[left]
                if branch == "held_EMG_and_load":
                    mass = seq["m"][left]
            if branch == "causal_EMG_trend":
                a = np.clip(
                    inputs[left]
                    + np.clip(
                        slope, -options["inputSlopeLimit"], options["inputSlopeLimit"]
                    )
                    * 0.1
                    * (1 - np.exp(-(k - left) * 0.01 / 0.1)),
                    0,
                    1,
                )
                mass = seq["m"][left]
            if branch == "maximum_rate_limited_assistance":
                u = min(options["capNm"], state[28] + options["rateNmPerSec"] * 0.01)
            trace.append(plant.step(state, a, mass, wait, u, 0.01, 10))
        traces[branch] = np.asarray(trace)
    idx = peak - left - 1
    row = dict(
        subject=folder.name,
        recordId=seg["recordId"],
        source=source,
        loadKg=float(seg["loadKg"]),
        split=seg["split"],
        peakTimeSeconds=float(seq["t"][peak]),
        thresholdNm=limit,
        peakBackTorqueNm=float(saved["backTorque"][peak]),
        excessPeakNm=float(max(0, score)),
        commandAtPeakNm=float(saved["assistanceTorque"][peak]),
        commandCapNm=float(options["capNm"]),
        capActiveAtPeak=bool(
            saved["assistanceTorque"][peak] >= options["capNm"] - 0.01
        ),
        sourceActivationAtPeak=float(inputs[peak, 0]),
        backCommandAtPeak=float(saved["backCommand"][peak]),
        backEffectiveActivationAtPeak=float(saved["backActivation"][peak]),
        activeBackAtPeakNm=float(saved["backActive"][peak]),
        passiveBackAtPeakNm=float(saved["backPassive"][peak]),
        fiberDampingAtPeakNm=float(saved["backFiberDamping"][peak]),
        boundaryAtPeakNm=float(saved["boundaryTorque"][peak]),
        replayMaximumTorqueDifferenceNm=parity,
        checkpointLeadSeconds=(peak - left) * 0.01,
        diagnosticFutureSamplesUsedOnlyInOfflineBranches=True,
        counterfactualPeakTimeTorqueNm={
            name: float(trace[idx, 3]) for name, trace in traces.items()
        },
        heldForecastTorqueErrorAtPeakNm=float(
            traces["recorded_inputs_and_controls"][idx, 3]
            - traces["held_EMG_and_load"][idx, 3]
        ),
        controllerTorqueCeilingWasHardConstraint=False,
    )
    matching = [
        d
        for d in quality["cycles"]
        if d["subject"] == folder.name and d["id"] == seg["recordId"]
    ]
    row["qualityGateDecision"] = (
        matching[0] if matching else "No matching complete filtered action"
    )
    with PdfPages(output / (folder.name + "_exceedance_forensics.pdf")) as pdf:
        fig, axes = plt.subplots(4, 1, figsize=(11, 9), sharex=True)
        ids = np.arange(max(0, peak - 100), min(len(inputs), peak + 101))
        t = seq["t"][ids]
        axes[0].plot(t, inputs[ids, 0], label="Source back activation")
        axes[0].plot(t, saved["backCommand"][ids], label="Neural back excitation")
        axes[0].plot(t, saved["backActivation"][ids], label="Hill effective activation")
        axes[0].set_ylabel("Activation")
        axes[1].plot(
            t,
            saved["backTorque"][ids],
            color="#1469b0",
            lw=2,
            label="Applied MPC result",
        )
        branch_time = seq["t"][left + 1 : right + 1]
        styles = {
            "held_EMG_and_load": ("--", "#e69935"),
            "causal_EMG_trend": (":", "#8054a6"),
            "maximum_rate_limited_assistance": ("--", "#23844a"),
        }
        for name, (ls, color) in styles.items():
            axes[1].plot(
                branch_time,
                traces[name][:, 3],
                ls=ls,
                color=color,
                label=name.replace("_", " "),
            )
        axes[1].axhline(limit, color="red", ls="--", label="Back torque threshold")
        axes[1].set_ylabel("Back torque (N m)")
        for name in ("backActive", "backPassive", "backFiberDamping", "boundaryTorque"):
            axes[2].plot(t, saved[name][ids], label=name)
        axes[2].set_ylabel("Torque components (N m)")
        axes[3].plot(
            t,
            saved["assistanceTorque"][ids],
            color="#23844a",
            label="Applied assistance",
        )
        axes[3].axhline(options["capNm"], color="red", ls="--", label="Actuator cap")
        axes[3].set(ylabel="Assistance (N m)", xlabel="Replay time (s)")
        for ax in axes:
            ax.axvline(seq["t"][peak], color=".4", ls=":")
            ax.grid(alpha=0.2)
            ax.legend(fontsize=7, loc="upper right")
        fig.suptitle(
            folder.name
            + f" | worst case +{max(0,score):.2f} N m | matched-state offline diagnosis"
        )
        fig.tight_layout()
        pdf.savefig(fig)
        fig.savefig(output / (folder.name + "_exceedance_forensics.png"), dpi=140)
        plt.close(fig)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=ROOT / "results/multi_subject")
    parser.add_argument(
        "--subjects", nargs="+", default=["240617Pilot", "240627Pilot", "240719Pilot"]
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "results/multi_subject_quality_v2/forensics_before",
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    quality = json.loads(
        (ROOT / "results/multi_subject_quality_v2/quality_gate.json").read_text()
    )
    rows = [
        diagnose(args.results / subject, args.output, quality)
        for subject in args.subjects
    ]
    (args.output / "exceedance_diagnosis.json").write_text(
        json.dumps(rows, indent=2) + "\n"
    )
    for row in rows:
        print(
            row["subject"],
            row["recordId"],
            "excess",
            row["excessPeakNm"],
            "cap",
            row["capActiveAtPeak"],
            "hold error",
            row["heldForecastTorqueErrorAtPeakNm"],
            "branches",
            row["counterfactualPeakTimeTorqueNm"],
            flush=True,
        )


if __name__ == "__main__":
    main()
