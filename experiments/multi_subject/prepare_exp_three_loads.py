"""Prepare three BSE anchor loads from complete excitation/IK action bouts.

The user accepts only 5/10/15 kg BSE scenarios. Missing CEINMS torque samples
are masked, never manufactured. Angle scaling retains its historical affine
coefficients, and all already-supported historical samples are retained.
MAP coefficients remain frozen: the existing MAP already used full EMG bouts.
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

from audit_subjects import DATA, ROOT, extract, read_sto, REQUIRED
from quality_gate import gaussian_matlab, native_tables, sha, smoothness
from prepare_subjects import combine, structs, map_predict
from fit_hill_compliant import input_generator

BASE = ROOT / "results/exp_three_loads"
VERSION = "excitation-ik-three-loads-v1"
SUBJECT = "240725Exp_A01"
# Select complete, already-reviewed sensor bouts. Extra scenarios are holdouts,
# not an admission requirement. The selection precedes BSE/MPC fitting.
TRAIN = {5: 45, 10: 45, 15: 55}
TEST = {5: 35, 10: 25, 15: 65}


def supported(clock, query, maximum_gap=0.25):
    """No endpoint extrapolation and no interpolation across a large gap."""
    index = np.searchsorted(clock, query)
    left = np.maximum(index - 1, 0)
    right = np.minimum(index, len(clock) - 1)
    exact = np.minimum(abs(query - clock[left]), abs(query - clock[right])) < 1e-6
    return (
        (query >= clock[0])
        & (query <= clock[-1])
        & (exact | (clock[right] - clock[left] <= maximum_gap))
    )


def angle_reference(path, query):
    original = extract(path)
    names, ik = read_sto(path / "ik.sto")
    unique, first = np.unique(ik[:, 0], return_index=True)
    ik = ik[np.sort(first)]
    if np.any(np.diff(ik[:, 0]) <= 0) or not np.isfinite(ik).all():
        raise ValueError("Invalid IK timestamps/values")
    raw = original["qRaw"]
    filtered = gaussian_matlab(raw, 5)
    if np.ptp(filtered) <= 0.015:
        raise ValueError("Insufficient raw angle excursion")
    gain = 0.9 / np.ptp(filtered)
    offset = 0.5 - gain * filtered.min()
    native = original["absoluteStart"] + original["t"]
    # Use the same five-sample operation at the native LMT median cadence for
    # IK-only extensions. No second min/max scale is fitted to the chosen bout.
    delta = float(np.median(np.diff(native)))
    start = max(ik[0, 0], min(native[0], query[0]) - 3 * delta)
    finish = min(ik[-1, 0], max(native[-1], query[-1]) + 3 * delta)
    grid = start + np.arange(int((finish - start) / delta) + 1) * delta
    qraw = np.interp(grid, ik[:, 0], ik[:, names.index("L5_S1_Flex_Ext")])
    q = gain * np.interp(query, grid, gaussian_matlab(qraw, 5)) + offset
    retained = supported(native, query)
    q[retained] = np.interp(query[retained], native, original["q"])
    parity = (
        float(
            np.max(abs(q[retained] - np.interp(query[retained], native, original["q"])))
        )
        if retained.any()
        else 0.0
    )
    if not supported(ik[:, 0], query).all():
        raise ValueError("Chosen action lacks continuous measured IK support")
    return q, dict(
        gain=float(gain),
        offset=float(offset),
        historicalSupportedFraction=float(retained.mean()),
        historicalReferenceMaximumDifference=parity,
        referenceMin=float(q.min()),
        referenceMax=float(q.max()),
        ikSHA256=sha(path / "ik.sto"),
    )


def prepare(output=BASE):
    review = json.loads((ROOT / "docs/exp_one_cycle_review/review.json").read_text())
    previous_gate = json.loads(
        (ROOT / "results/exp_subjects_A01/quality_gate.json").read_text()
    )
    model_file = ROOT / "results/exp_subjects_A01/map" / (SUBJECT + "_map.mat")
    map_cache = loadmat(model_file, simplify_cells=True)
    model = map_cache["mapModel"]
    map_train_paths = set(np.atleast_1d(map_cache["trainingTrialPaths"]).tolist())
    records, decisions = [], []
    for split, choices in (("train", TRAIN), ("test_condition", TEST)):
        for load, speed in choices.items():
            row = next(
                r
                for r in review["conditions"]
                if r["subject"] == "240725Exp"
                and r["loadKg"] == load
                and r["speedCode"] == speed
            )
            source = DATA / SUBJECT / "selected_windows" / f"{load}kg_{speed}.mat"
            item = loadmat(source, simplify_cells=True)
            t = item["timeEpochSeconds"]
            path = DATA / SUBJECT / row["trial"]
            q, q_audit = angle_reference(path, t)
            tables, clocks, identity, exact, issues, repairs = native_tables(path)
            if issues:
                raise ValueError("Native pipeline issue: " + str(issues))
            six = item["activation"][:, :6]
            back = six @ model["lumpWeights"]
            front = item["activation"][:, 6]
            if (
                min(six.min(), front.min()) < -1e-6
                or max(six.max(), front.max()) > 1 + 1e-6
            ):
                raise ValueError("Nonphysiological activation input")
            diagnostic = smoothness(back)
            rough = [
                k
                for k, v in previous_gate["smoothnessUpperLimits"].items()
                if diagnostic[k] > v
            ]
            if rough:
                raise ValueError("Activation needs smoothing review: " + str(rough))
            torque_table = tables["Torque.sto"][1]
            valid = supported(torque_table[:, 0], t)
            torque = np.interp(t, torque_table[:, 0], torque_table[:, 1])
            torque[~valid] = np.nan
            correlation = (
                float(np.corrcoef(q[valid], torque[valid])[0, 1])
                if valid.sum() > 2
                else 0.0
            )
            duration = float(t[-1] - t[0])
            r = dict(
                id=f"{load}kg_{speed}_complete_M",
                trialPath=row["trial"],
                loadKg=load,
                speedCode=speed,
                cadenceBPM=60 / duration,
                duration=duration,
                sourceStartSeconds=float(t[0] - row["originEpochSeconds"]),
                sourceEndSeconds=float(t[-1] - row["originEpochSeconds"]),
                q=q,
                torque=torque,
                torqueReferenceValid=valid,
                back=back,
                front=front,
                six=six,
                weights=np.tile(model["lumpWeights"], (len(t), 1)),
                angleTorqueCorrelation=correlation,
                angleUsable=True,
                rightCensored=False,
                split=split,
            )
            r["mappedBack"] = map_predict(model, r)
            records.append(r)
            decisions.append(
                dict(
                    id=r["id"],
                    trialPath=r["trialPath"],
                    loadKg=load,
                    speedCode=speed,
                    split=split,
                    sourceFiles={
                        n: sha(path / n) for n in (*REQUIRED, "emg.sto", "emgFilt.sto")
                    },
                    activationWindowSHA256=sha(source),
                    theta=q_audit,
                    torqueReferenceCoverageFraction=float(valid.mean()),
                    angleTorqueCorrelation=correlation,
                    torquePhaseWarning=correlation >= -0.25,
                    torqueIdentity=identity,
                    smoothness=diagnostic,
                    MAPWindowWasInTraining=row["trial"] in map_train_paths,
                    CEINMSTorqueUsedAsGroundTruth=False,
                )
            )
    train = [r for r in records if r["split"] == "train"]
    generator = input_generator(
        [
            dict(
                a=np.c_[r["back"], r["mappedBack"]],
                front=r["front"],
                m=np.full(len(r["q"]), r["loadKg"]),
                loadKg=r["loadKg"],
                q=r["q"],
            )
            for r in train
        ]
    )
    output.mkdir(parents=True, exist_ok=True)
    gate = dict(
        version=VERSION,
        subjects=[SUBJECT],
        BSELoadsKg=[5, 10, 15],
        scenarioRule="One complete excitation/IK action per 5/10/15 kg anchor suffices for calibration; extra complete bouts are BSE holdouts.",
        MAPAlreadyTrained=True,
        MAPFileSHA256=sha(model_file),
        minimumNumberOfHoldoutsIsAdmissionRequirement=False,
        referenceAnglePolarityChanged=False,
        independentTorqueReferenceAvailable=False,
        missingTorqueIsMasked=True,
        CEINMSPhaseCorrelationIsWarningNotAdmissionVeto=True,
        trials=decisions,
        filterPolicy=previous_gate["filterPolicy"],
    )
    gate_file = output / "quality_gate.json"
    gate_file.write_text(json.dumps(gate, indent=2) + "\n")
    folder = output / SUBJECT
    folder.mkdir(exist_ok=True)
    if not (folder / "prepared.mat").exists():
        savemat(
            folder / "prepared.mat",
            dict(
                sequences=structs([combine([r], SUBJECT, r["id"]) for r in records]),
                trainSequences=structs([combine([r], SUBJECT, r["id"]) for r in train]),
                cycles=structs(
                    [
                        {k: v for k, v in r.items() if k not in ("six", "weights")}
                        for r in records
                    ]
                ),
                mapModel=model,
                generator=generator,
                qualityGateVersion=VERSION,
                qualityGateSHA256=sha(gate_file),
            ),
            long_field_names=True,
            do_compression=True,
        )
    metrics = [
        {
            k: r[k]
            for k in (
                "id",
                "trialPath",
                "loadKg",
                "speedCode",
                "cadenceBPM",
                "duration",
                "split",
                "angleUsable",
                "angleTorqueCorrelation",
                "rightCensored",
                "mapOutOfDomain",
            )
        }
        | dict(
            mapRMSE=float(np.sqrt(np.mean((r["mappedBack"] - r["back"]) ** 2))),
            mapCorrelation=float(np.corrcoef(r["mappedBack"], r["back"])[0, 1]),
        )
        for r in records
    ]
    report = dict(
        subject=SUBJECT,
        trainingCycles=3,
        testCycles=3,
        knownConditions=[[r["loadKg"], r["speedCode"]] for r in records],
        angleUsableCycles=6,
        qualityCheckedBeforeTraining=True,
        qualityGateVersion=VERSION,
        qualityGateSHA256=sha(gate_file),
        mapArchitecture=model["architecture"],
        cycles=metrics,
        MAPParametersRefitted=False,
        MAPHoldoutIndependent=False,
        BSEHoldoutIndependent=True,
        limitation="The separate BSE holdouts may have contributed to the previously trained MAP; these are not end-to-end pipeline holdouts.",
    )
    (folder / "preparation.json").write_text(json.dumps(report, indent=2) + "\n")
    (output / "data_quality.json").write_text(json.dumps(decisions, indent=2) + "\n")
    fig, axes = plt.subplots(3, 3, figsize=(13, 9))
    for i, load in enumerate((5, 10, 15)):
        r = next(r for r in train if r["loadKg"] == load)
        time = np.arange(len(r["q"])) / 100
        axes[i, 0].plot(time, r["back"], label="Experimental activation")
        axes[i, 0].plot(time, r["mappedBack"], label="Frozen trained MAP")
        axes[i, 1].plot(time, r["q"], color="#1469b0")
        axes[i, 2].plot(time, r["torque"], color="#c8841c")
        valid = r["torqueReferenceValid"]
        axes[i, 2].set_title(
            f"CEINMS estimate: {valid.mean():.0%} available\nMissing samples omitted",
            fontsize=10,
        )
        axes[i, 2].set_xlim(0, time[-1])
        if not valid.all():
            axes[i, 2].fill_between(
                time,
                0,
                1,
                where=~valid,
                transform=axes[i, 2].get_xaxis_transform(),
                color=".92",
                zorder=-2,
            )
        axes[i, 0].set_ylabel(f"{load} kg\nActivation")
        axes[i, 1].set_ylabel("Original-scale theta (rad)")
        axes[i, 2].set_ylabel("Torque (N m)")
        for ax in axes[i]:
            ax.grid(alpha=0.2)
            ax.set_xlabel("Complete action time (s)")
    axes[0, 0].legend(fontsize=8)
    axes[0, 1].set_title("Recorded IK\nHistorical Scale retained", fontsize=10)
    fig.suptitle("240725Exp_A01 | only 5/10/15 kg BSE calibration anchors")
    fig.tight_layout()
    fig.savefig(output / "three_anchor_inputs.png", dpi=160)
    plt.close(fig)
    print(
        json.dumps(
            [
                dict(
                    id=r["id"],
                    split=r["split"],
                    qRange=[float(r["q"].min()), float(r["q"].max())],
                    valleys=len(find_peaks(-r["q"], prominence=0.12, distance=50)[0]),
                    torqueCoverage=float(r["torqueReferenceValid"].mean()),
                    MAPWindowWasInTraining=r["trialPath"] in map_train_paths,
                )
                for r in records
            ],
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=BASE)
    prepare(parser.parse_args().output)
