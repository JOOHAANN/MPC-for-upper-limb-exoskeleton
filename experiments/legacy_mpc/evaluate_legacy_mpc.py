"""Plot actual archived MPC replays, keeping their model identity explicit."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "legacy_mpc"
SOURCES = ("control_experimental", "control_map")
LABELS = ("Experimental activation", "MAP activation")
# These are the movement windows in the archived plotting code, in samples.
WINDOWS = ((7757, 13853, 5), (20956, 27639, 10), (36908, 41734, 15))
BLUE, RED, GREY, GREEN = "#1768ac", "#d62728", "#858585", "#26824c"


def write_csv(path: Path, records: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def read_replay(source: str) -> dict:
    status = json.loads((OUT / f"{source}_checked_status.json").read_text(encoding="utf-8"))
    assert status["completed"], f"{source} did not complete; do not plot a zero-filled tail."
    work = Path(status["workDirectory"])
    data = loadmat(work / "legacy_execution.mat", simplify_cells=True)
    count = int(data["NumberOfData"])
    assert status["lastStep"] == count - 1
    data["status"] = status
    data["count"] = count
    data["time"] = np.arange(count) * float(data["Ts"])
    data["backTorque"] = np.asarray(data["MTF"])[:count] * data["coeffs"][0]
    data["assistanceTorque"] = np.asarray(data["Tac"])[:count]
    data["netMuscleTorque"] = data["backTorque"] - data["T_front"]
    # The stored Torque_Tac0 is the archived back-tendon output, not a fresh run.
    data["inputHash"] = hashlib.sha256((work / "Combination.mat").read_bytes()).hexdigest()
    manifest = json.loads((work / "workflow_manifest.json").read_text(encoding="utf-8"))
    for name in ("Combination.mat", "Theta_Tac0.mat", "Torque_Tac0.mat"):
        canonical = ROOT / next(e["source"] for e in manifest["files"] if e["target"] == name)
        assert canonical.read_bytes() == (work / name).read_bytes(), f"Archived input changed: {name}"
    return data


def style(ax, data: dict) -> None:
    ax.grid(alpha=.22)
    ax.set_xlim(0, data["time"][-1])
    for start, end, _ in WINDOWS:
        ax.axvspan(start * data["Ts"], end * data["Ts"], color="#f3b66d", alpha=.09, zorder=0)


def save(fig, name: str) -> None:
    fig.savefig(OUT / f"{name}.png", dpi=170, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    data = [read_replay(source) for source in SOURCES]
    summaries, movement, solvers, provenance = [], [], [], []
    for source, d in zip(SOURCES, data):
        t, q, back, assist = d["time"], d["theta"], d["backTorque"], d["assistanceTorque"]
        used = np.isfinite(d["legacyExitFlag"])
        flags = d["legacyExitFlag"][used]
        latency = 1000 * d["legacySolveSeconds"][used]
        finite = all(np.isfinite(d[name]).all() for name in
                     ("theta", "dtheta", "backTorque", "assistanceTorque", "T_front", "F_k_up", "T_balance"))
        summaries.append({
            "source": source, "durationSeconds": t[-1], "wallSeconds": d["status"]["elapsedSeconds"],
            "executionComplete": True, "allFinite": bool(finite),
            "minimumAngleRad": float(q.min()), "maximumAngleRad": float(q.max()),
            "thetaRMSEReferenceDeg": float(np.sqrt(np.mean((q - d["theta_Ceinms"]) ** 2)) * 180 / np.pi),
            "peakBackTorqueNm": float(back.max()), "peakAssistanceNm": float(assist.max()),
            "assistanceCapNm": float(np.max(d["u_max"])),
            "assistanceOverCapSamples": int(np.sum(assist > d["u_max"] + 1e-6)),
            "softBoundaryAbsPeakNm": float(np.max(np.abs(d["F_k_up"]))),
            "externalBalanceAbsPeakNm": float(np.max(np.abs(d["T_balance"]))),
            "minimumFrontTorqueNm": float(d["T_front"].min()),
            "maximumFrontTorqueNm": float(d["T_front"].max()),
            "solverCalls": int(used.sum()), "solverNonpositiveExits": int(np.sum(flags <= 0)),
            "solveMedianMs": float(np.median(latency)), "solveP95Ms": float(np.percentile(latency, 95)),
            "solveMaximumMs": float(latency.max()),
        })
        for start, end, load in WINDOWS:
            selected = (np.arange(1, len(t) + 1) > start) & (np.arange(1, len(t) + 1) < end)
            baseline = d["Torque_real"][selected]
            current = back[selected]
            movement.append({
                "source": source, "loadKg": load, "startSeconds": start * d["Ts"], "endSeconds": end * d["Ts"],
                "archivedUnassistedBackPeakNm": float(baseline.max()), "assistedBackPeakNm": float(current.max()),
                "backPeakReductionPercent": float(100 * (1 - current.max() / baseline.max())),
                "assistedBackMeanNm": float(current.mean()), "archivedUnassistedBackMeanNm": float(baseline.mean()),
                "referenceMinimumAngleRad": float(d["theta_Ceinms"][selected].min()),
                "assistedMinimumAngleRad": float(q[selected].min()), "assistedMaximumAngleRad": float(q[selected].max()),
                "thetaRMSEReferenceDeg": float(np.sqrt(np.mean((q[selected] - d["theta_Ceinms"][selected]) ** 2)) * 180 / np.pi),
                "peakAssistanceNm": float(assist[selected].max()),
                "backOverThresholdSamples": int(np.sum(current > d["Threshold"][selected] + 1e-6)),
            })
        for i in np.flatnonzero(used):
            solvers.append({"source": source, "sample": int(i + 1), "timeSeconds": float(t[i]),
                            "exitFlag": int(d["legacyExitFlag"][i]),
                            "solveSeconds": float(d["legacySolveSeconds"][i]),
                            "iterations": int(d["legacyIterations"][i]),
                            "functionEvaluations": int(d["legacyEvaluations"][i]),
                            "objective": float(d["legacyObjective"][i])})
        work = Path(d["status"]["workDirectory"])
        provenance.append({"source": source, "originalSimulator": d["status"]["simulator"],
                           "originalSimulatorSHA256": hashlib.sha256((work / (d["status"]["simulator"] + ".m")).read_bytes()).hexdigest(),
                           "executedCopySHA256": hashlib.sha256((work / (d["status"]["executedEntry"] + ".m")).read_bytes()).hexdigest(),
                           "predictorSHA256": hashlib.sha256((work / "nonlinear_system_period.m").read_bytes()).hexdigest(),
                           "inputSHA256": d["inputHash"],
                           "replaySHA256": hashlib.sha256((work / "legacy_execution.mat").read_bytes()).hexdigest(),
                           "parameters": d["coeffs"].tolist(), "executionChanges": d["status"]["executionChanges"],
                           "modelPolicy": d["status"]["modelPolicy"]})
    write_csv(OUT / "summary.csv", summaries)
    write_csv(OUT / "movement_metrics.csv", movement)
    write_csv(OUT / "solver_diagnostics.csv", solvers)
    (OUT / "evaluation.json").write_text(json.dumps({"summary": summaries, "provenance": provenance,
        "scope": "Each archived simulator's original 45.4-second, three-load recording",
        "controller": {"sampleSeconds": .001, "predictionSamples": 100, "controlBlocks": 5,
                       "blockSamples": 10, "actuatorCapNm": 80, "incrementBoundsNm": [0, .5],
                       "Q": 12, "R": 3, "W": 3, "experimentalThresholdsNm": [100, 140, 180],
                       "MAPThresholdsNm": [90, 120, 150]},
        "baselinePolicy": "Archived Theta_Tac0/Torque_Tac0; not an independently regenerated baseline",
        "torquePolicy": "Legacy plotted output is R*MTF (back tendon torque); abdomen is a signed polynomial"},
        indent=2) + "\n", encoding="utf-8")

    plt.rcParams.update({"font.size": 10, "axes.titlesize": 12, "axes.labelsize": 11})
    fig, axes = plt.subplots(3, 2, figsize=(14, 10.2), sharex="col")
    for col, (label, d, summary) in enumerate(zip(LABELS, data, summaries)):
        t = d["time"]
        ax = axes[0, col]
        ax.plot(t, d["Torque_Ceinms"], color=RED, ls="--", lw=1.3, label="L5/S1 reference")
        ax.plot(t, d["Torque_real"], color=GREY, lw=1.1, label="Archived unassisted output")
        ax.plot(t, d["backTorque"], color=BLUE, lw=1.6, label="Legacy MPC assisted output")
        ax.plot(t, d["Threshold"], color="#333333", ls=":", lw=1.1, label="Back torque threshold")
        ax.set_title(label)
        ax.set_ylabel("Legacy back torque (N·m)")
        if col == 0:
            handles, labels = ax.get_legend_handles_labels()
        ax = axes[1, col]
        ax.plot(t, d["theta_Ceinms"], color=RED, ls="--", lw=1.3)
        ax.plot(t, d["theta_real"], color=GREY, lw=1.1)
        ax.plot(t, d["theta"], color=BLUE, lw=1.6)
        ax.set_ylabel("Trunk angle (rad)")
        ax.set_title(f"Actual assisted range: {summary['minimumAngleRad']:.2f} to {summary['maximumAngleRad']:.2f} rad", fontsize=11)
        ax = axes[2, col]
        ax.plot(t, d["assistanceTorque"], color=GREEN, lw=1.6, label="Applied actuator torque")
        ax.plot(t, d["u_max"], color="#a83232", ls="--", lw=1.0, label="Actuator cap")
        above = d["assistanceTorque"] > d["u_max"] + 1e-6
        ax.fill_between(t, d["u_max"], d["assistanceTorque"], where=above, color="#e88a8a", alpha=.5)
        ax.set_ylabel("Applied assistance (N·m)")
        ax.set_title(f"Peak: {summary['peakAssistanceNm']:.2f} N·m; {summary['solverNonpositiveExits']} solver exits <= 0", fontsize=11)
        ax.set_xlabel("Time (s)")
        ax.legend(loc="upper right", frameon=False, fontsize=9)
        for row in range(3):
            style(axes[row, col], d)
    fig.suptitle("Archived MPC simulator test: original BSE and controller settings", y=.993, fontsize=15)
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(.5, .01))
    fig.tight_layout(rect=(0, .047, 1, .97), h_pad=1.8)
    save(fig, "legacy_mpc_torque_theta")

    fig, axes = plt.subplots(3, 2, figsize=(14, 9.8), sharex="col")
    for col, (label, d) in enumerate(zip(LABELS, data)):
        t = d["time"]
        ax = axes[0, col]
        ax.plot(t, d["F_k_up"], color="#904cad", label="Old soft-boundary torque")
        ax.plot(t, d["T_balance"], color="#d48714", label="External low-activation balance torque", lw=1.1)
        ax.set_title(label)
        ax.set_ylabel("External terms (N·m)")
        ax.legend(frameon=False, loc="best", fontsize=9)
        ax = axes[1, col]
        ax.plot(t, d["T_front"], color="#a83232")
        ax.axhline(0, color="#555555", lw=.7)
        ax.set_ylabel("Old abdominal term (N·m)")
        ax.set_title("Signed angle polynomial in the legacy model; not a Hill muscle", fontsize=10)
        ax = axes[2, col]
        used = np.isfinite(d["legacyExitFlag"])
        ax.scatter(t[used], d["legacyExitFlag"][used], s=3, color=BLUE, alpha=.4)
        ax.axhline(0, color="#a83232", ls="--", lw=.9)
        ax.set_yticks([-2, 0, 1, 2])
        ax.set_ylabel("fmincon exit flag")
        ax.set_xlabel("Time (s)")
        ax.set_title("Positive: solver termination met a convergence criterion", fontsize=10)
        for row in range(3):
            style(axes[row, col], d)
    fig.suptitle("Legacy MPC: external torques, abdominal term and solver outcomes", y=.995, fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, .97), h_pad=2)
    save(fig, "legacy_mpc_components")
    for summary in summaries:
        print(json.dumps(summary))


if __name__ == "__main__":
    main()
