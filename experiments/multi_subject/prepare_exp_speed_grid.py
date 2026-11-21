"""Audit and prepare the 725 5/10/15 kg by 25--65 BSE calibration grid.

The existing MAP and historical IK scale are frozen. Missing CEINMS torque
samples remain missing. Central speeds receive priority without case-specific
model parameters. All grid conditions are calibration, not unseen validation.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat, savemat
from scipy.signal import find_peaks

from audit_subjects import DATA, ROOT, REQUIRED
from prepare_exp_three_loads import SUBJECT, angle_reference, supported
from quality_gate import native_tables, sha, smoothness
from prepare_subjects import combine, structs, map_predict
from fit_hill_compliant import input_generator

BASE = ROOT / "results/exp_speed_grid_torque"
VERSION = "excitation-ik-speed-grid-torque-v1"
LOADS = (5, 10, 15)
SPEEDS = (25, 35, 45, 55, 65)


def collect():
    review = json.loads((ROOT / "docs/exp_one_cycle_review/review.json").read_text())
    previous = json.loads(
        (ROOT / "results/exp_subjects_A01/quality_gate.json").read_text()
    )
    map_path = ROOT / "results/exp_subjects_A01/map" / (SUBJECT + "_map.mat")
    cache = loadmat(map_path, simplify_cells=True)
    model = cache["mapModel"]
    records, decisions = [], []
    for load in LOADS:
        for speed in SPEEDS:
            row = next(
                r
                for r in review["conditions"]
                if r["subject"] == "240725Exp"
                and r["loadKg"] == load
                and r["speedCode"] == speed
            )
            source = DATA / SUBJECT / "selected_windows" / f"{load}kg_{speed}.mat"
            item = loadmat(source, simplify_cells=True)
            epoch = item["timeEpochSeconds"]
            path = DATA / SUBJECT / row["trial"]
            q, q_audit = angle_reference(path, epoch)
            tables, clocks, identity, exact, issues, repairs = native_tables(path)
            if issues:
                raise ValueError(f"{load}/{speed}: native pipeline issues {issues}")
            six = item["activation"][:, :6]
            back = six @ model["lumpWeights"]
            front = item["activation"][:, 6]
            if (
                not np.isfinite(np.c_[six, front]).all()
                or min(six.min(), front.min()) < -1e-6
                or max(six.max(), front.max()) > 1 + 1e-6
            ):
                raise ValueError(f"{load}/{speed}: invalid activation")
            diagnostic = smoothness(back)
            rough = [
                k
                for k, v in previous["smoothnessUpperLimits"].items()
                if diagnostic[k] > v
            ]
            if rough:
                raise ValueError(
                    f"{load}/{speed}: activation smoothness review {rough}"
                )
            table = tables["Torque.sto"][1]
            valid = supported(table[:, 0], epoch)
            if not valid.any():
                decisions.append(
                    dict(
                        id=f"{load}kg_{speed}_complete_M",
                        trialPath=row["trial"],
                        loadKg=load,
                        speedCode=speed,
                        included=False,
                        exclusionReason="No CEINMS torque support in the selected action; excluded by user instruction.",
                        torqueReferenceCoverageFraction=0.0,
                        sourceFiles={
                            n: sha(path / n)
                            for n in (*REQUIRED, "emg.sto", "emgFilt.sto")
                        },
                        activationWindowSHA256=sha(source),
                    )
                )
                continue
            torque = np.interp(epoch, table[:, 0], table[:, 1])
            torque[~valid] = np.nan
            correlation = (
                float(np.corrcoef(q[valid], torque[valid])[0, 1])
                if valid.sum() > 2
                else None
            )
            duration = float(epoch[-1] - epoch[0])
            r = dict(
                id=f"{load}kg_{speed}_complete_M",
                trialPath=row["trial"],
                loadKg=load,
                speedCode=speed,
                cadenceBPM=60 / duration,
                duration=duration,
                sourceStartSeconds=float(epoch[0] - row["originEpochSeconds"]),
                sourceEndSeconds=float(epoch[-1] - row["originEpochSeconds"]),
                q=q,
                torque=torque,
                torqueReferenceValid=valid,
                back=back,
                front=front,
                six=six,
                weights=np.tile(model["lumpWeights"], (len(q), 1)),
                angleTorqueCorrelation=correlation if correlation is not None else 0.0,
                angleUsable=True,
                rightCensored=False,
                split="train",
                fitWeight=2.0 if 35 <= speed <= 55 else 1.0,
                boundaryResidualWeight=2.0,
            )
            r["mappedBack"] = map_predict(model, r)
            records.append(r)
            decisions.append(
                dict(
                    id=r["id"],
                    trialPath=r["trialPath"],
                    loadKg=load,
                    speedCode=speed,
                    split="train",
                    included=True,
                    fitWeight=r["fitWeight"],
                    sourceFiles={
                        n: sha(path / n) for n in (*REQUIRED, "emg.sto", "emgFilt.sto")
                    },
                    activationWindowSHA256=sha(source),
                    theta=q_audit,
                    referenceValleys=len(
                        find_peaks(-q, prominence=0.12, distance=50)[0]
                    ),
                    torqueReferenceCoverageFraction=float(valid.mean()),
                    angleTorqueCorrelation=correlation,
                    torqueIdentity=identity,
                    smoothness=diagnostic,
                    structuralRepairs=repairs,
                    CEINMSTorqueUsedAsGroundTruth=False,
                )
            )
    return records, decisions, model, map_path, previous


def prepare(output=BASE, audit_only=False):
    records, decisions, model, map_path, previous = collect()
    output.mkdir(parents=True, exist_ok=True)
    gate = dict(
        version=VERSION,
        subjects=[SUBJECT],
        BSELoadsKg=list(LOADS),
        speedCodes=list(SPEEDS),
        prioritySpeedCodes=[35, 45, 55],
        allConditionsAreCalibration=True,
        excludeConditionsWithZeroTorqueSupport=True,
        independentTorqueReferenceAvailable=False,
        missingTorqueIsMasked=True,
        referenceAnglePolarityChanged=False,
        historicalScaleRetained=True,
        MAPAlreadyTrained=True,
        MAPFileSHA256=sha(map_path),
        trials=decisions,
        filterPolicy=previous["filterPolicy"],
    )
    gate_path = output / "quality_gate.json"
    if not audit_only and (output / SUBJECT / "prepared.mat").exists():
        old = json.loads(gate_path.read_text())
        if old != gate:
            raise RuntimeError(
                "Prepared input provenance changed; use a new output folder."
            )
    gate_path.write_text(json.dumps(gate, indent=2) + "\n")
    for d in decisions:
        if not d["included"]:
            print(
                f"{d['loadKg']:2g} kg / {d['speedCode']} | EXCLUDED: no torque support",
                flush=True,
            )
            continue
        print(
            f"{d['loadKg']:2g} kg / {d['speedCode']} | theta {d['theta']['referenceMin']:.3f}..{d['theta']['referenceMax']:.3f} | valleys {d['referenceValleys']} | torque {d['torqueReferenceCoverageFraction']:.0%}",
            flush=True,
        )
    fig, axes = plt.subplots(3, 5, figsize=(17, 9), sharey=True)
    excluded_panels(axes, decisions)
    for r in records:
        ax = axes[LOADS.index(r["loadKg"]), SPEEDS.index(r["speedCode"])]
        d = next(d for d in decisions if d["id"] == r["id"])
        t = np.arange(len(r["q"])) / 100
        ax.plot(t, r["q"], color="#bb3d2d", label="Original-scale IK theta")
        ax.set_title(
            f"{r['loadKg']} kg / {r['speedCode']} | torque {d['torqueReferenceCoverageFraction']:.0%}",
            fontsize=10,
        )
        ax.set_xlabel("Action time (s)")
        ax.grid(alpha=0.2)
        if 35 <= r["speedCode"] <= 55:
            ax.set_facecolor("#fff8e8")
    for ax in axes[:, 0]:
        ax.set_ylabel("Reference theta (rad)")
    fig.suptitle(
        "725 | recorded reference audit, historical Scale retained\nYellow panels: priority speeds; torque coverage is availability, not accuracy"
    )
    fig.tight_layout()
    fig.savefig(output / "reference_grid.png", dpi=150)
    plt.close(fig)
    fig, axes = plt.subplots(3, 5, figsize=(17, 9), sharey=True)
    excluded_panels(axes, decisions)
    for r in records:
        ax = axes[LOADS.index(r["loadKg"]), SPEEDS.index(r["speedCode"])]
        t = np.arange(len(r["q"])) / 100
        flexion = 1.5 - r["q"]
        flexion = (flexion - flexion.mean()) / max(flexion.std(), 0.08)
        ax.plot(
            t, flexion, "--", color="#c13b31", label="Recorded flexion (standardized)"
        )
        for key, color, label in (
            ("back", "#1469b0", "Experimental activation"),
            ("mappedBack", "#8450ad", "Trained MAP activation"),
        ):
            signal = r[key]
            ax.plot(
                t,
                (signal - signal.mean()) / max(signal.std(), 0.02),
                color=color,
                label=label,
            )
        ax.set_title(f"{r['loadKg']} kg / {r['speedCode']}", fontsize=10)
        ax.set_xlabel("Action time (s)")
        ax.grid(alpha=0.2)
    for ax in axes[:, 0]:
        ax.set_ylabel("Standardized signal (diagnostic only)")
    fig.legend(
        *axes.flat[0].get_legend_handles_labels(),
        loc="lower center",
        ncol=3,
        frameon=False,
    )
    fig.suptitle(
        "725 | measured flexion vs input activation timing\nStandardization is for this diagnostic plot only; fitting retains original angles and activation"
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(output / "activation_timing_grid.png", dpi=150)
    plt.close(fig)
    if audit_only:
        return
    folder = output / SUBJECT
    folder.mkdir(exist_ok=True)
    # The auxiliary abdomen predictor uses activation labels only.
    generator = input_generator(
        [
            dict(
                a=np.c_[r["back"], r["mappedBack"]],
                front=r["front"],
                m=np.full(len(r["q"]), r["loadKg"]),
                loadKg=r["loadKg"],
                q=r["q"],
            )
            for r in records
        ]
    )
    sequences = [combine([r], SUBJECT, r["id"]) for r in records]
    if not (folder / "prepared.mat").exists():
        savemat(
            folder / "prepared.mat",
            dict(
                sequences=structs(sequences),
                trainSequences=structs(sequences),
                cycles=structs(
                    [
                        {k: v for k, v in r.items() if k not in ("six", "weights")}
                        for r in records
                    ]
                ),
                mapModel=model,
                generator=generator,
                qualityGateVersion=VERSION,
                qualityGateSHA256=sha(gate_path),
            ),
            long_field_names=True,
            do_compression=True,
        )
    (folder / "preparation.json").write_text(
        json.dumps(
            dict(
                subject=SUBJECT,
                trainingCycles=len(records),
                priorityCycles=sum(35 <= r["speedCode"] <= 55 for r in records),
                excludedCycles=sum(not d["included"] for d in decisions),
                independentTestCycles=0,
                allConditionsAreCalibration=True,
                MAPParametersRefitted=False,
                qualityGateVersion=VERSION,
                qualityGateSHA256=sha(gate_path),
                conditionPriorityWeights={"25": 1, "35": 2, "45": 2, "55": 2, "65": 1},
            ),
            indent=2,
        )
        + "\n"
    )


def excluded_panels(axes, decisions):
    for d in decisions:
        if not d["included"]:
            ax = axes[LOADS.index(d["loadKg"]), SPEEDS.index(d["speedCode"])]
            ax.set_facecolor(".92")
            ax.set_title(f"{d['loadKg']} kg / {d['speedCode']}", fontsize=10)
            ax.text(
                0.5,
                0.5,
                "Excluded\nNo CEINMS torque\nin selected action",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.set_xticks([])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=BASE)
    parser.add_argument("--audit-only", action="store_true")
    args = parser.parse_args()
    prepare(args.output, args.audit_only)
